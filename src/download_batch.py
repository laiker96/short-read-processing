#!/usr/bin/env python3
"""Validate an accession sample sheet and download its FASTQ files."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.cli import (
    add_download_arguments,
    cli_main,
    execute_download,
)
from short_read_processing.sample_sheet import DEFAULT_SCHEMA, sample_sheet_accessions


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate a canonical CSV/TSV sample sheet and download all FASTQs concurrently."
    )
    parser.add_argument("sample_sheet", type=Path, help="Canonical comma- or tab-separated file")
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA)
    add_download_arguments(parser)
    args = parser.parse_args()
    accessions = sample_sheet_accessions(args.sample_sheet, schema_path=args.schema)
    execute_download(accessions, args)
    return 0


if __name__ == "__main__":
    cli_main(main)
