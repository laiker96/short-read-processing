"""High-throughput ENA FASTQ downloads with an SRA Toolkit fallback."""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .accessions import AcquisitionError, FilePlan, RunPlan, classify_files


@dataclass(frozen=True)
class DownloadOptions:
    file_jobs: int
    connections: int
    sra_jobs: int
    threads: int
    keep_sra_cache: bool = False


def _require_executable(name: str) -> str:
    executable = shutil.which(name)
    if not executable:
        raise AcquisitionError(
            f"Required executable {name!r} was not found on PATH; recreate or update the repo-local .venv"
        )
    return executable


def _run(command: list[str], *, label: str) -> None:
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        raise AcquisitionError(f"{label} failed with exit code {exc.returncode}") from exc


def _aria2_input(files: Iterable[FilePlan]) -> str:
    blocks: list[str] = []
    seen: dict[Path, FilePlan] = {}
    for item in files:
        previous = seen.get(item.path)
        if previous:
            if previous.url != item.url or previous.md5 != item.md5:
                raise AcquisitionError(f"Conflicting sources for destination {item.path}")
            continue
        seen[item.path] = item
        item.path.parent.mkdir(parents=True, exist_ok=True)
        blocks.append(item.url)
        blocks.append(f"  dir={item.path.parent}")
        blocks.append(f"  out={item.path.name}")
        if item.md5:
            blocks.append(f"  checksum=md5={item.md5}")
    return "\n".join(blocks) + "\n"


def _discard_untracked_partial_files(files: Iterable[FilePlan]) -> list[Path]:
    """Remove size-mismatched files that aria2 cannot safely resume.

    Aria2 writes a sibling ``.aria2`` control file for its own partial
    downloads. A partial file copied by another tool (for example, an
    interrupted rsync) has no segment map and must be restarted from byte zero.
    Complete-size files are left for aria2's checksum verification.
    """

    discarded: list[Path] = []
    for item in files:
        control = item.path.with_name(item.path.name + ".aria2")
        if (
            item.path.is_file()
            and not control.exists()
            and item.size_bytes is not None
            and item.path.stat().st_size != item.size_bytes
        ):
            size = item.path.stat().st_size
            item.path.unlink()
            discarded.append(item.path)
            print(
                f"Restarting untracked partial FASTQ {item.path} "
                f"({size} of {item.size_bytes} bytes)",
                flush=True,
            )
    return discarded


def download_ena(plans: list[RunPlan], options: DownloadOptions) -> None:
    if not plans:
        return
    aria2c = _require_executable("aria2c")
    files = [item for plan in plans for item in plan.files]
    if not files:
        raise AcquisitionError("ENA download was selected but no FASTQ files were resolved")

    _discard_untracked_partial_files(files)
    input_text = _aria2_input(files)
    output_root = Path(os.path.commonpath([str(plan.run_dir.parent) for plan in plans]))
    output_root.mkdir(parents=True, exist_ok=True)
    input_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".aria2-fastq-",
            suffix=".txt",
            dir=output_root,
            delete=False,
        ) as handle:
            handle.write(input_text)
            input_path = Path(handle.name)

        command = [
            aria2c,
            f"--input-file={input_path}",
            f"--max-concurrent-downloads={max(1, options.file_jobs)}",
            f"--max-connection-per-server={max(1, min(options.connections, 16))}",
            f"--split={max(1, min(options.connections, 16))}",
            "--min-split-size=1M",
            "--continue=true",
            "--check-integrity=true",
            "--file-allocation=none",
            "--auto-file-renaming=false",
            "--allow-overwrite=false",
            "--max-tries=0",
            "--retry-wait=5",
            "--connect-timeout=30",
            "--timeout=60",
            "--disk-cache=64M",
            "--summary-interval=10",
            "--console-log-level=notice",
        ]
        _run(command, label="aria2c FASTQ download")
    finally:
        if input_path:
            input_path.unlink(missing_ok=True)

    missing = [str(item.path) for item in files if not item.path.is_file() or item.path.stat().st_size == 0]
    if missing:
        raise AcquisitionError("Download completed without expected FASTQ files:\n  " + "\n  ".join(missing))
    for plan in plans:
        plan.status = "downloaded"


def _gzip_fastq(path: Path, *, threads: int) -> Path:
    pigz = shutil.which("pigz")
    if pigz:
        _run([pigz, "--processes", str(max(1, threads)), str(path)], label=f"pigz compression of {path.name}")
        return path.with_suffix(path.suffix + ".gz")

    output = path.with_suffix(path.suffix + ".gz")
    with path.open("rb") as source, gzip.open(output, "wb", compresslevel=6) as destination:
        shutil.copyfileobj(source, destination, length=1024 * 1024)
    path.unlink()
    return output


def _existing_sra_fastqs(plan: RunPlan) -> list[Path]:
    marker = plan.run_dir / ".download-complete"
    files = sorted(plan.run_dir.glob(f"{plan.run_accession}*.fastq.gz"))
    return files if marker.is_file() and files else []


def _assign_local_sra_files(plan: RunPlan, paths: list[Path]) -> None:
    files = [
        FilePlan(url="", md5="", size_bytes=path.stat().st_size, path=path.resolve())
        for path in sorted(paths)
    ]
    plan.files = classify_files(files, plan.library_layout)


def _download_one_sra(plan: RunPlan, *, threads: int, keep_cache: bool) -> None:
    existing = _existing_sra_fastqs(plan)
    if existing:
        _assign_local_sra_files(plan, existing)
        plan.status = "existing"
        return

    prefetch = _require_executable("prefetch")
    fasterq_dump = _require_executable("fasterq-dump")
    plan.run_dir.mkdir(parents=True, exist_ok=True)
    output_root = plan.run_dir.parent
    cache_dir = output_root / ".sra-cache" / plan.run_accession
    temporary_dir = output_root / ".fasterq-tmp" / plan.run_accession
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir.mkdir(parents=True, exist_ok=True)

    _run(
        [prefetch, plan.run_accession, "--max-size", "u", "--output-directory", str(cache_dir)],
        label=f"prefetch {plan.run_accession}",
    )
    _run(
        [
            fasterq_dump,
            str(cache_dir),
            "--split-files",
            "--threads",
            str(max(1, threads)),
            "--temp",
            str(temporary_dir),
            "--outdir",
            str(plan.run_dir),
        ],
        label=f"fasterq-dump {plan.run_accession}",
    )

    uncompressed = sorted(plan.run_dir.glob(f"{plan.run_accession}*.fastq"))
    if not uncompressed:
        raise AcquisitionError(f"fasterq-dump produced no FASTQ files for {plan.run_accession}")
    compressed = [_gzip_fastq(path, threads=threads) for path in uncompressed]
    _assign_local_sra_files(plan, compressed)
    (plan.run_dir / ".download-complete").write_text("complete\n", encoding="utf-8")
    plan.status = "downloaded"

    shutil.rmtree(temporary_dir, ignore_errors=True)
    if not keep_cache:
        shutil.rmtree(cache_dir, ignore_errors=True)


def download_sra(plans: list[RunPlan], options: DownloadOptions) -> None:
    if not plans:
        return
    jobs = max(1, min(options.sra_jobs, len(plans)))
    threads_per_job = max(1, options.threads // jobs)
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(
                _download_one_sra,
                plan,
                threads=threads_per_job,
                keep_cache=options.keep_sra_cache,
            ): plan.run_accession
            for plan in plans
        }
        for future in as_completed(futures):
            run = futures[future]
            try:
                future.result()
            except Exception as exc:
                errors.append(f"{run}: {exc}")
    if errors:
        raise AcquisitionError("SRA fallback failed:\n  " + "\n  ".join(sorted(errors)))


def download_plans(plans: list[RunPlan], options: DownloadOptions) -> None:
    """Download ENA plans in one aria2 queue, then concurrent SRA fallbacks."""

    download_ena([plan for plan in plans if plan.backend == "ena"], options)
    download_sra([plan for plan in plans if plan.backend == "sra"], options)
