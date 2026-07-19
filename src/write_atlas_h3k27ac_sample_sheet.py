#!/usr/bin/env python3
"""Convert selected atlas H3K27ac IP metadata rows to the canonical pipeline TSV."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


OUTPUT_FIELDS = (
    "accession",
    "sample_id",
    "assay",
    "genome",
    "role",
    "control_id",
    "replicate",
    "peak_caller",
    "notes",
)
REPLICATE_RE = re.compile(r"rep(?:licate)?[_ -]?([1-9][0-9]*)", re.IGNORECASE)
SAFE_CONTEXT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _replicate_number(value: str, *, accession: str) -> int:
    value = value.strip()
    match = REPLICATE_RE.fullmatch(value)
    if match:
        return int(match.group(1))
    if value.lower() == "single":
        return 1
    raise ValueError(f"Selected run {accession} has invalid biological_replicate: {value!r}")


def selected_sample_rows(metadata: Path) -> list[dict[str, str | int]]:
    with metadata.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "context_id",
            "context",
            "atlas_use",
            "selection_status",
            "assay_role",
            "biological_replicate",
            "run_accession",
            "sample_title",
            "library_layout",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing atlas metadata columns: " + ", ".join(sorted(missing)))
        selected = [
            row
            for row in reader
            if row["atlas_use"].strip().lower() == "yes"
            and row["selection_status"].strip().lower() == "selected"
        ]

    output: list[dict[str, str | int]] = []
    accessions: set[str] = set()
    for row in selected:
        accession = row["run_accession"].strip().upper()
        if accession in accessions:
            raise ValueError(f"Duplicate selected run accession: {accession}")
        accessions.add(accession)
        if row["assay_role"].strip().lower() != "h3k27ac_ip":
            raise ValueError(f"Selected row {accession} is not an H3K27ac IP")

        layout = row["library_layout"].strip().upper()
        if layout not in {"PAIRED", "SINGLE"}:
            raise ValueError(f"Selected run {accession} has unsupported layout: {layout!r}")
        replicate = _replicate_number(row["biological_replicate"], accession=accession)

        context_id = row["context_id"].strip()
        if not SAFE_CONTEXT_RE.fullmatch(context_id):
            raise ValueError(f"Selected run {accession} has invalid context_id: {context_id!r}")
        output.append(
            {
                "accession": accession,
                "sample_id": f"{context_id.lower()}_h3k27ac_rep{replicate}",
                "assay": "chip_histone",
                "genome": "dm6",
                "role": "treatment",
                "control_id": "",
                "replicate": replicate,
                "peak_caller": "callpeak",
                "notes": (
                    f"atlas selected IP-only: {row['context'].strip()}; "
                    f"{row['sample_title'].strip()}; layout={layout}"
                ),
            }
        )
    if not output:
        raise ValueError("Atlas metadata contains no selected H3K27ac IP rows")
    return output


def write_sample_sheet(metadata: Path, output: Path | None) -> int:
    rows = selected_sample_rows(metadata)
    handle = output.open("w", encoding="utf-8", newline="") if output else sys.stdout
    try:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if output:
            handle.close()
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("--output", type=Path, help="Output TSV; omit to write to stdout")
    args = parser.parse_args()
    count = write_sample_sheet(args.metadata, args.output)
    if args.output:
        print(f"Wrote {count} selected H3K27ac IP runs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
