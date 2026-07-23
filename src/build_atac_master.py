#!/usr/bin/env python3
"""Build a summit-aware master DHS registry across final ATAC contexts."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.master_dhs import build_master_registry


def named_path(value: str) -> tuple[str, Path]:
    name, separator, path = value.partition("=")
    if not separator or not name or not path:
        raise argparse.ArgumentTypeError("use CONTEXT=PATH")
    return name, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-peaks", action="append", type=named_path, required=True)
    parser.add_argument("--context-signal", action="append", type=named_path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--summit-max-distance", type=int, default=150)
    parser.add_argument("--minimum-summit-separation", type=int, default=50)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output-bed", type=Path, required=True)
    parser.add_argument("--summit-bed", type=Path, required=True)
    parser.add_argument("--membership-tsv", type=Path, required=True)
    parser.add_argument("--context-matrix-tsv", type=Path, required=True)
    parser.add_argument("--stats-json", type=Path, required=True)
    args = parser.parse_args()
    context_peaks = dict(args.context_peaks)
    context_signals = dict(args.context_signal)
    if len(context_peaks) != len(args.context_peaks):
        parser.error("--context-peaks context IDs must be unique")
    if len(context_signals) != len(args.context_signal):
        parser.error("--context-signal context IDs must be unique")
    metrics = build_master_registry(
        context_peaks=context_peaks,
        context_signals=context_signals,
        chrom_sizes_path=args.chrom_sizes,
        output_bed=args.output_bed,
        summit_bed=args.summit_bed,
        membership_tsv=args.membership_tsv,
        context_matrix_tsv=args.context_matrix_tsv,
        stats_json=args.stats_json,
        summit_max_distance=args.summit_max_distance,
        minimum_summit_separation=args.minimum_summit_separation,
        workers=args.workers,
    )
    print(
        f"source_peaks={metrics['source_peak_count']} "
        f"master_dhs={metrics['master_dhs_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
