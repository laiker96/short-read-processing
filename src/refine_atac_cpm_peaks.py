#!/usr/bin/env python3
"""Progressively refine MACS3 peaks using a CPM BigWig."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.cpm_refinement import refine_cpm_bigwig


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--signal-bigwig", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--merge-gap-bp", type=int, default=1)
    parser.add_argument("--minimum-length", type=int, default=50)
    parser.add_argument("--maximum-length", type=int, default=400)
    parser.add_argument("--minimum-mean-cpm", type=float, default=0.0)
    parser.add_argument("--minimum-mode-prominence", type=float, default=0.25)
    args = parser.parse_args()
    refine_cpm_bigwig(
        peaks=args.peaks,
        signal_bigwig=args.signal_bigwig,
        output=args.output,
        excluded=args.excluded,
        stats=args.stats,
        merge_gap_bp=args.merge_gap_bp,
        minimum_length=args.minimum_length,
        maximum_length=args.maximum_length,
        minimum_mean_cpm=args.minimum_mean_cpm,
        minimum_mode_prominence=args.minimum_mode_prominence,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
