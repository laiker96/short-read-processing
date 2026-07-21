#!/usr/bin/env python3
"""Build within-condition consensus, cross-condition atlas, or signal shapes."""

from __future__ import annotations

import argparse
from pathlib import Path

from short_read_processing.atlas import (
    build_condition_consensus,
    build_dhs_center_mode_half_prominence,
    build_dhs_support_fwhm,
    build_global_atlas,
    build_signal_shaped_atlas,
)


def _replicate(value: str) -> tuple[str, Path]:
    sample, separator, path = value.partition("=")
    if not separator or not sample or not path:
        raise argparse.ArgumentTypeError("use SAMPLE_ID=REFINED_BED")
    return sample, Path(path)


def _condition_path(value: str) -> tuple[str, Path]:
    condition, separator, path = value.partition("=")
    if not separator or not condition or not path:
        raise argparse.ArgumentTypeError("use CONDITION_ID=PATH")
    return condition, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    condition = subparsers.add_parser("condition")
    condition.add_argument("--condition-id", required=True)
    condition.add_argument("--pooled-peaks", type=Path, required=True)
    condition.add_argument("--replicate", action="append", type=_replicate, required=True)
    condition.add_argument("--minimum-replicates", type=int, default=2)
    condition.add_argument("--overlap-fraction", type=float, default=0.5)
    condition.add_argument("--output-bed", type=Path, required=True)
    condition.add_argument("--support-tsv", type=Path, required=True)
    condition.add_argument("--stats-json", type=Path, required=True)

    atlas = subparsers.add_parser("atlas")
    atlas.add_argument(
        "--condition",
        action="append",
        nargs=3,
        metavar=("CONDITION_ID", "CONSENSUS_BED", "CPM_BIGWIG"),
        required=True,
    )
    atlas.add_argument("--peak-width", type=int, default=250)
    atlas.add_argument(
        "--grouping-method",
        choices=("fixed_window", "fixed_window_narrow_first", "dhs_seed"),
        default="fixed_window",
    )
    atlas.add_argument("--output-bed", type=Path, required=True)
    atlas.add_argument("--variable-bed", type=Path, required=True)
    atlas.add_argument("--membership-tsv", type=Path, required=True)
    atlas.add_argument("--presence-tsv", type=Path, required=True)
    atlas.add_argument("--coverage-tsv", type=Path, required=True)
    atlas.add_argument("--mean-cpm-tsv", type=Path, required=True)
    atlas.add_argument("--maximum-cpm-tsv", type=Path, required=True)
    atlas.add_argument("--stats-json", type=Path, required=True)

    shape = subparsers.add_parser("shape")
    shape.add_argument("--membership-tsv", type=Path, required=True)
    shape.add_argument(
        "--condition-bigwig",
        action="append",
        type=_condition_path,
        required=True,
    )
    shape.add_argument("--output-bed", type=Path, required=True)
    shape.add_argument("--aggregate-bigwig", type=Path, required=True)
    shape.add_argument("--diagnostics-tsv", type=Path, required=True)
    shape.add_argument("--stats-json", type=Path, required=True)
    shape.add_argument("--window-size", type=int, default=1000)
    shape.add_argument("--bin-size", type=int, default=10)
    shape.add_argument("--smoothing-bins", type=int, default=3)
    shape.add_argument("--relative-threshold", type=float, default=0.2)
    shape.add_argument("--background-mad-multiplier", type=float, default=3.0)
    shape.add_argument("--minimum-length", type=int, default=50)
    shape.add_argument("--maximum-length", type=int, default=400)

    fwhm = subparsers.add_parser("fwhm")
    fwhm.add_argument("--anchors-bed", type=Path, required=True)
    fwhm.add_argument(
        "--condition-dhs",
        action="append",
        type=_condition_path,
        required=True,
    )
    fwhm.add_argument("--chrom-sizes", type=Path, required=True)
    fwhm.add_argument("--support-bigwig", type=Path, required=True)
    fwhm.add_argument("--output-bed", type=Path, required=True)
    fwhm.add_argument("--diagnostics-tsv", type=Path, required=True)
    fwhm.add_argument("--stats-json", type=Path, required=True)

    center_mode = subparsers.add_parser("center-mode")
    center_mode.add_argument("--anchors-bed", type=Path, required=True)
    center_mode.add_argument(
        "--condition-dhs",
        action="append",
        type=_condition_path,
        required=True,
    )
    center_mode.add_argument("--chrom-sizes", type=Path, required=True)
    center_mode.add_argument("--output-bed", type=Path, required=True)
    center_mode.add_argument("--diagnostics-tsv", type=Path, required=True)
    center_mode.add_argument("--stats-json", type=Path, required=True)
    args = parser.parse_args()

    if args.command == "condition":
        replicates = dict(args.replicate)
        if len(replicates) != len(args.replicate):
            parser.error("--replicate sample IDs must be unique")
        build_condition_consensus(
            condition_id=args.condition_id,
            pooled_peaks=args.pooled_peaks,
            replicate_peaks=replicates,
            output_bed=args.output_bed,
            support_tsv=args.support_tsv,
            stats_json=args.stats_json,
            minimum_replicates=args.minimum_replicates,
            overlap_fraction=args.overlap_fraction,
        )
    elif args.command == "atlas":
        build_global_atlas(
            conditions=[
                (condition_id, Path(bed), Path(bigwig))
                for condition_id, bed, bigwig in args.condition
            ],
            output_bed=args.output_bed,
            variable_bed=args.variable_bed,
            membership_tsv=args.membership_tsv,
            presence_tsv=args.presence_tsv,
            coverage_tsv=args.coverage_tsv,
            mean_cpm_tsv=args.mean_cpm_tsv,
            maximum_cpm_tsv=args.maximum_cpm_tsv,
            stats_json=args.stats_json,
            peak_width=args.peak_width,
            grouping_method=args.grouping_method,
        )
    elif args.command == "shape":
        condition_bigwigs = dict(args.condition_bigwig)
        if len(condition_bigwigs) != len(args.condition_bigwig):
            parser.error("--condition-bigwig condition IDs must be unique")
        build_signal_shaped_atlas(
            membership_tsv=args.membership_tsv,
            condition_bigwigs=condition_bigwigs,
            output_bed=args.output_bed,
            aggregate_bigwig=args.aggregate_bigwig,
            diagnostics_tsv=args.diagnostics_tsv,
            stats_json=args.stats_json,
            window_size=args.window_size,
            bin_size=args.bin_size,
            smoothing_bins=args.smoothing_bins,
            relative_threshold=args.relative_threshold,
            background_mad_multiplier=args.background_mad_multiplier,
            minimum_length=args.minimum_length,
            maximum_length=args.maximum_length,
        )
    elif args.command == "fwhm":
        condition_dhs = dict(args.condition_dhs)
        if len(condition_dhs) != len(args.condition_dhs):
            parser.error("--condition-dhs condition IDs must be unique")
        build_dhs_support_fwhm(
            anchors_bed=args.anchors_bed,
            condition_dhs=condition_dhs,
            chromosome_sizes_path=args.chrom_sizes,
            support_bigwig=args.support_bigwig,
            output_bed=args.output_bed,
            diagnostics_tsv=args.diagnostics_tsv,
            stats_json=args.stats_json,
        )
    else:
        condition_dhs = dict(args.condition_dhs)
        if len(condition_dhs) != len(args.condition_dhs):
            parser.error("--condition-dhs condition IDs must be unique")
        build_dhs_center_mode_half_prominence(
            anchors_bed=args.anchors_bed,
            condition_dhs=condition_dhs,
            chromosome_sizes_path=args.chrom_sizes,
            output_bed=args.output_bed,
            diagnostics_tsv=args.diagnostics_tsv,
            stats_json=args.stats_json,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
