#!/usr/bin/env python3
"""Convert selected atlas ATAC metadata rows to the canonical pipeline TSV."""

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


def selected_sample_rows(metadata: Path) -> list[dict[str, str]]:
    with metadata.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "context_id",
            "context",
            "atlas_use",
            "selection_status",
            "biological_replicate",
            "technical_run",
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

    output = []
    accessions = set()
    for row in selected:
        accession = row["run_accession"].strip().upper()
        if accession in accessions:
            raise ValueError(f"Duplicate selected run accession: {accession}")
        accessions.add(accession)
        if row["library_layout"].strip().upper() != "PAIRED":
            raise ValueError(f"Selected run {accession} is not paired-end; HMMRATAC is invalid")

        replicate = row["biological_replicate"].strip()
        if not replicate:
            match = REPLICATE_RE.search(row["sample_title"])
            if not match:
                raise ValueError(
                    f"Selected run {accession} lacks a biological replicate number "
                    "and none can be inferred from sample_title"
                )
            replicate = match.group(1)

        context_id = row["context_id"].strip().lower()
        technical_run = row["technical_run"].strip()
        output.append(
            {
                "accession": accession,
                "sample_id": f"{context_id}_atac_rep{replicate}",
                "assay": "atac",
                "genome": "dm6",
                "role": "treatment",
                "control_id": "",
                "replicate": replicate,
                "peak_caller": "hmmratac",
                "notes": (
                    f"atlas selected: {row['context'].strip()}; "
                    f"{row['sample_title'].strip()}; technical_run={technical_run}"
                ),
            }
        )
    if not output:
        raise ValueError("Atlas metadata contains no selected ATAC rows")
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
        print(f"Wrote {count} selected ATAC runs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
