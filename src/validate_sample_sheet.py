#!/usr/bin/env python3
"""Validate and summarize the canonical accession sample sheet."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from short_read_processing.cli import cli_main
from short_read_processing.sample_sheet import DEFAULT_SCHEMA, read_sample_sheet


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_sheet", type=Path)
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    args = parser.parse_args()

    rows = read_sample_sheet(args.sample_sheet, schema_path=args.schema)
    samples = {str(row["library_id"]) for row in rows}
    assays = Counter(str(row["assay"]) for row in rows)
    callers = Counter(
        str(row["peak_caller"]) for row in rows if str(row["role"]) == "treatment"
    )
    print(f"valid: {len(rows)} accession row(s), {len(samples)} biological library/libraries")
    print("assays: " + ", ".join(f"{key}={value}" for key, value in sorted(assays.items())))
    print("treatment peak callers: " + ", ".join(
        f"{key}={value}" for key, value in sorted(callers.items())
    ))
    return 0


if __name__ == "__main__":
    cli_main(main)
