"""Resolve SRA/ENA run and experiment accessions through ENA file reports."""

from __future__ import annotations

import csv
import io
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Iterable


ENA_FILE_REPORT_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = (
    "run_accession",
    "experiment_accession",
    "library_layout",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
)
RUN_ACCESSION_RE = re.compile(r"^(?:SRR|ERR)\d+$")
EXPERIMENT_ACCESSION_RE = re.compile(r"^(?:SRX|ERX)\d+$")


class AcquisitionError(RuntimeError):
    """Raised when accessions cannot be resolved or downloaded safely."""


@dataclass
class FilePlan:
    """A remote FASTQ and its intended local destination."""

    url: str
    md5: str
    size_bytes: int | None
    path: Path
    mate: str = "extra"


@dataclass
class RunPlan:
    """Resolved download information for one sequencing run."""

    requested_accession: str
    experiment_accession: str
    run_accession: str
    library_layout: str
    backend: str
    run_dir: Path
    files: list[FilePlan] = field(default_factory=list)
    status: str = "planned"


def normalize_accession(value: str) -> str:
    accession = value.strip().upper()
    if not (RUN_ACCESSION_RE.fullmatch(accession) or EXPERIMENT_ACCESSION_RE.fullmatch(accession)):
        raise AcquisitionError(
            f"Unsupported accession {value!r}; expected SRR, SRX, ERR, or ERX followed by digits"
        )
    return accession


def accession_kind(accession: str) -> str:
    accession = normalize_accession(accession)
    return "run" if RUN_ACCESSION_RE.fullmatch(accession) else "experiment"


def _split_semicolon(value: str) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def _https_url(ena_path: str) -> str:
    if ena_path.startswith("https://"):
        return ena_path
    if ena_path.startswith("http://"):
        return "https://" + ena_path.removeprefix("http://")
    if ena_path.startswith("ftp://"):
        return "https://" + ena_path.removeprefix("ftp://")
    return "https://" + ena_path.lstrip("/")


def _mate_from_name(name: str) -> str:
    lowered = name.lower()
    if lowered.endswith("_1.fastq.gz") or lowered.endswith("_r1.fastq.gz"):
        return "r1"
    if lowered.endswith("_2.fastq.gz") or lowered.endswith("_r2.fastq.gz"):
        return "r2"
    return "extra"


def classify_files(files: list[FilePlan], layout: str) -> list[FilePlan]:
    """Assign mate labels and return files in R1, R2, extra order."""

    for item in files:
        item.mate = _mate_from_name(item.path.name)

    if layout == "PAIRED":
        r1 = [item for item in files if item.mate == "r1"]
        r2 = [item for item in files if item.mate == "r2"]
        extras = [item for item in files if item.mate == "extra"]
        if not r1 and not r2 and len(files) == 2:
            files[0].mate = "r1"
            files[1].mate = "r2"
            r1, r2, extras = [files[0]], [files[1]], []
        if len(r1) != 1 or len(r2) != 1:
            raise AcquisitionError(
                "Paired run did not resolve to exactly one R1 and one R2 FASTQ: "
                + ", ".join(item.path.name for item in files)
            )
        return r1 + r2 + extras

    if not files:
        return files
    files[0].mate = "r1"
    for item in files[1:]:
        item.mate = "extra"
    return files


def build_file_report_url(accession: str) -> str:
    query = urllib.parse.urlencode(
        {
            "accession": normalize_accession(accession),
            "result": "read_run",
            "fields": ",".join(ENA_FIELDS),
            "format": "tsv",
        }
    )
    return f"{ENA_FILE_REPORT_URL}?{query}"


def fetch_ena_report(accession: str, *, timeout: int = 60, retries: int = 4) -> str:
    """Fetch a small ENA run file report with bounded retries."""

    request = urllib.request.Request(
        build_file_report_url(accession),
        headers={"User-Agent": "short-read-processing/0.1"},
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8")
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504}:
                break
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
        if attempt + 1 < retries:
            time.sleep(min(2**attempt, 8))
    raise AcquisitionError(f"ENA lookup failed for {accession}: {last_error}")


def parse_ena_report(
    text: str,
    *,
    requested_accession: str,
    output_dir: Path,
    backend: str = "auto",
) -> list[RunPlan]:
    """Parse ENA TSV into run-level download plans."""

    requested_accession = normalize_accession(requested_accession)
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames or "run_accession" not in reader.fieldnames:
        raise AcquisitionError(f"ENA returned an invalid file report for {requested_accession}")

    plans: list[RunPlan] = []
    for row in reader:
        run = row.get("run_accession", "").strip().upper()
        if not RUN_ACCESSION_RE.fullmatch(run):
            continue
        layout = row.get("library_layout", "").strip().upper()
        if layout not in {"SINGLE", "PAIRED"}:
            raise AcquisitionError(f"Unsupported library layout {layout!r} for {run}")

        urls = [_https_url(item) for item in _split_semicolon(row.get("fastq_ftp", ""))]
        md5s = _split_semicolon(row.get("fastq_md5", ""))
        sizes = _split_semicolon(row.get("fastq_bytes", ""))
        if md5s and len(md5s) != len(urls):
            raise AcquisitionError(f"ENA returned mismatched FASTQ and MD5 counts for {run}")
        if sizes and len(sizes) != len(urls):
            raise AcquisitionError(f"ENA returned mismatched FASTQ and byte counts for {run}")

        run_dir = (output_dir / run).resolve()
        files: list[FilePlan] = []
        for index, url in enumerate(urls):
            filename = PurePosixPath(urllib.parse.urlparse(url).path).name
            if not filename or filename in {".", ".."} or os.sep in filename:
                raise AcquisitionError(f"Unsafe FASTQ filename returned for {run}: {filename!r}")
            size = int(sizes[index]) if sizes and sizes[index] else None
            files.append(
                FilePlan(
                    url=url,
                    md5=md5s[index] if md5s else "",
                    size_bytes=size,
                    path=run_dir / filename,
                )
            )

        direct_usable = bool(files)
        if direct_usable:
            try:
                files = classify_files(files, layout)
            except AcquisitionError:
                direct_usable = False

        selected_backend = backend
        if backend == "auto":
            selected_backend = "ena" if direct_usable else "sra"
        elif backend == "ena" and not direct_usable:
            raise AcquisitionError(f"ENA has no complete compressed FASTQ set for {run}")
        elif backend not in {"ena", "sra"}:
            raise AcquisitionError(f"Unsupported backend: {backend}")

        plans.append(
            RunPlan(
                requested_accession=requested_accession,
                experiment_accession=row.get("experiment_accession", "").strip().upper(),
                run_accession=run,
                library_layout=layout,
                backend=selected_backend,
                run_dir=run_dir,
                files=files if selected_backend == "ena" else [],
            )
        )

    if not plans:
        raise AcquisitionError(f"No sequencing runs found for {requested_accession}")
    return sorted(plans, key=lambda item: item.run_accession)


def resolve_accession(accession: str, *, output_dir: Path, backend: str = "auto") -> list[RunPlan]:
    accession = normalize_accession(accession)
    return parse_ena_report(
        fetch_ena_report(accession),
        requested_accession=accession,
        output_dir=output_dir,
        backend=backend,
    )


def resolve_accessions(
    accessions: Iterable[str],
    *,
    output_dir: Path,
    backend: str = "auto",
    workers: int = 8,
) -> list[RunPlan]:
    """Resolve accessions concurrently while returning deterministic ordering."""

    normalized = list(dict.fromkeys(normalize_accession(item) for item in accessions))
    if not normalized:
        raise AcquisitionError("No accessions were provided")

    by_accession: dict[str, list[RunPlan]] = {}
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=max(1, min(workers, len(normalized)))) as executor:
        futures = {
            executor.submit(resolve_accession, item, output_dir=output_dir, backend=backend): item
            for item in normalized
        }
        for future in as_completed(futures):
            accession = futures[future]
            try:
                by_accession[accession] = future.result()
            except Exception as exc:  # preserve all resolution failures for batch users
                errors.append(f"{accession}: {exc}")
    if errors:
        raise AcquisitionError("Could not resolve all accessions:\n  " + "\n  ".join(sorted(errors)))

    return [plan for accession in normalized for plan in by_accession[accession]]

