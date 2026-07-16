"""Stable TSV manifest I/O shared by download and config commands."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path
from typing import Iterable

from .accessions import FilePlan, RunPlan


MANIFEST_FIELDS = (
    "requested_accession",
    "experiment_accession",
    "run_accession",
    "library_layout",
    "backend",
    "status",
    "fastq_1",
    "fastq_2",
    "extra_fastqs",
    "md5_1",
    "md5_2",
    "extra_md5s",
)


def _one(files: list[FilePlan], mate: str) -> FilePlan | None:
    return next((item for item in files if item.mate == mate), None)


def plan_to_row(plan: RunPlan) -> dict[str, str]:
    r1 = _one(plan.files, "r1")
    r2 = _one(plan.files, "r2")
    extras = [item for item in plan.files if item.mate == "extra"]
    return {
        "requested_accession": plan.requested_accession,
        "experiment_accession": plan.experiment_accession,
        "run_accession": plan.run_accession,
        "library_layout": plan.library_layout,
        "backend": plan.backend,
        "status": plan.status,
        "fastq_1": str(r1.path) if r1 else "",
        "fastq_2": str(r2.path) if r2 else "",
        "extra_fastqs": ";".join(str(item.path) for item in extras),
        "md5_1": r1.md5 if r1 else "",
        "md5_2": r2.md5 if r2 else "",
        "extra_md5s": ";".join(item.md5 for item in extras),
    }


def read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Manifest not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(MANIFEST_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest {path} is missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def write_manifest(path: Path, plans: Iterable[RunPlan], *, merge: bool = True) -> None:
    """Atomically write plans, updating matching requested-accession/run rows."""

    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_by_key: dict[tuple[str, str], dict[str, str]] = {}
    if merge and path.exists():
        for row in read_manifest(path):
            rows_by_key[(row["requested_accession"], row["run_accession"])] = row
    for plan in plans:
        row = plan_to_row(plan)
        key = (row["requested_accession"], row["run_accession"])
        previous = rows_by_key.get(key)
        if row["status"] == "planned" and previous and previous["status"] in {"downloaded", "existing"}:
            continue
        rows_by_key[key] = row

    rows = sorted(rows_by_key.values(), key=lambda row: (row["requested_accession"], row["run_accession"]))
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            encoding="utf-8",
            prefix=f".{path.name}.",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_name = handle.name
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    finally:
        if temporary_name and os.path.exists(temporary_name):
            os.unlink(temporary_name)
