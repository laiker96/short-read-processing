#!/usr/bin/env python3
"""Download FASTQ files for one SRR/SRX/ERR/ERX accession."""

from __future__ import annotations

import argparse

from short_read_processing.cli import add_download_arguments, cli_main, execute_download


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve one run or experiment accession and download compressed FASTQs. "
            "ENA HTTPS plus aria2c is preferred; SRA Toolkit is the fallback."
        )
    )
    parser.add_argument("accession", help="SRR, SRX, ERR, or ERX accession")
    add_download_arguments(parser)
    args = parser.parse_args()
    execute_download([args.accession], args)
    return 0


if __name__ == "__main__":
    cli_main(main)

