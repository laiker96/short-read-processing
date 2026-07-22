#!/usr/bin/env python3
"""Refine lenient ATAC candidates with MACS3 qpois signal."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.qpois_refinement import run_refinement


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qpois", type=Path, required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--name-prefix", required=True)
    parser.add_argument("--minimum-exponent", type=int, default=2)
    parser.add_argument("--maximum-exponent", type=int, default=325)
    parser.add_argument("--minimum-length", type=int, default=50)
    parser.add_argument("--maximum-length", type=int, default=400)
    parser.add_argument("--merge-gap", type=int, default=1)
    args = parser.parse_args()
    metrics = run_refinement(
        qpois_bedgraph=args.qpois,
        candidate_peaks=args.peaks,
        output_bed=args.output,
        excluded_bed=args.excluded,
        stats_json=args.stats,
        name_prefix=args.name_prefix,
        minimum_exponent=args.minimum_exponent,
        maximum_exponent=args.maximum_exponent,
        minimum_length=args.minimum_length,
        maximum_length=args.maximum_length,
        merge_gap=args.merge_gap,
    )
    print(
        f"refined={metrics['refined_peaks']} excluded={metrics['excluded_peaks']} "
        f"candidates={metrics['candidate_peaks']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
