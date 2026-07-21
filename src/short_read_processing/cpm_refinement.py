"""Refine MACS3 candidates from high to low CPM BigWig signal."""

from __future__ import annotations

import bisect
import json
import math
import os
from dataclasses import dataclass, replace
from pathlib import Path
from statistics import median


REFINEMENT_ALGORITHM = "descending_cpm_prominence_watershed_v1"


@dataclass(frozen=True)
class SignalInterval:
    chrom: str
    start: int
    end: int
    maximum_cpm: float
    weighted_signal: float
    signal_bases: int
    selection_cutoff_cpm: float = 0.0

    @classmethod
    def from_bigwig(
        cls, chrom: str, start: int, end: int, signal: float
    ) -> "SignalInterval":
        length = end - start
        return cls(chrom, start, end, signal, signal * length, length)

    @property
    def mean_cpm(self) -> float:
        return self.weighted_signal / self.signal_bases


class PeakContainmentIndex:
    """Answer whether an interval is fully contained in any candidate peak."""

    def __init__(self, peaks: dict[str, list[tuple[int, int]]]):
        self.starts: dict[str, list[int]] = {}
        self.maximum_ends: dict[str, list[int]] = {}
        for chrom, intervals in peaks.items():
            intervals.sort()
            maximum = -1
            starts: list[int] = []
            maximum_ends: list[int] = []
            for start, end in intervals:
                starts.append(start)
                maximum = max(maximum, end)
                maximum_ends.append(maximum)
            self.starts[chrom] = starts
            self.maximum_ends[chrom] = maximum_ends

    def contains(self, chrom: str, start: int, end: int) -> bool:
        starts = self.starts.get(chrom)
        if not starts:
            return False
        index = bisect.bisect_right(starts, start) - 1
        return index >= 0 and self.maximum_ends[chrom][index] >= end


def _read_peaks(
    path: Path,
) -> tuple[PeakContainmentIndex, dict[str, list[tuple[int, int]]], int]:
    peaks: dict[str, list[tuple[int, int]]] = {}
    count = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 columns")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            peaks.setdefault(fields[0], []).append((start, end))
            count += 1
    return PeakContainmentIndex(peaks), peaks, count


def progressively_refine_cpm(
    intervals: list[SignalInterval],
    *,
    merge_gap_bp: int = 1,
    minimum_length: int = 50,
    maximum_length: int = 400,
    minimum_mean_cpm: float = 0.0,
    minimum_mode_prominence: float = 0.25,
) -> tuple[list[SignalInterval], list[SignalInterval], list[float]]:
    """Grow high-CPM modes and retain sufficiently prominent mode mergers."""

    if (
        merge_gap_bp < 0
        or minimum_length < 1
        or maximum_length < minimum_length
        or not math.isfinite(minimum_mean_cpm)
        or minimum_mean_cpm < 0
        or not math.isfinite(minimum_mode_prominence)
        or not 0 <= minimum_mode_prominence <= 1
    ):
        raise ValueError("Invalid refinement parameter")
    if not intervals:
        return [], [], []
    intervals = sorted(
        intervals,
        key=lambda item: (item.chrom, item.start, item.end),
    )
    thresholds = sorted(
        {interval.maximum_cpm for interval in intervals}, reverse=True
    )
    activation_order = sorted(
        range(len(intervals)),
        key=lambda index: intervals[index].maximum_cpm,
        reverse=True,
    )
    parent = list(range(len(intervals)))
    active = [False] * len(intervals)
    starts = [item.start for item in intervals]
    ends = [item.end for item in intervals]
    maxima = [item.maximum_cpm for item in intervals]
    weighted_signals = [item.weighted_signal for item in intervals]
    signal_bases = [item.signal_bases for item in intervals]
    component_mode_ids: list[list[int] | None] = [None] * len(intervals)
    mode_snapshots: list[SignalInterval] = []
    retained_modes: list[bool] = []

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int, saddle_cpm: float) -> int:
        left_root, right_root = find(left), find(right)
        if left_root == right_root:
            return left_root
        if left_root > right_root:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        starts[left_root] = min(starts[left_root], starts[right_root])
        ends[left_root] = max(ends[left_root], ends[right_root])
        maxima[left_root] = max(maxima[left_root], maxima[right_root])
        weighted_signals[left_root] += weighted_signals[right_root]
        signal_bases[left_root] += signal_bases[right_root]
        left_modes = component_mode_ids[left_root]
        right_modes = component_mode_ids[right_root]
        if left_modes and right_modes:
            left_maximum = max(
                mode_snapshots[mode_id].maximum_cpm for mode_id in left_modes
            )
            right_maximum = max(
                mode_snapshots[mode_id].maximum_cpm for mode_id in right_modes
            )
            lower_maximum = min(left_maximum, right_maximum)
            prominence = (lower_maximum - saddle_cpm) / lower_maximum
            if prominence >= minimum_mode_prominence:
                left_modes.extend(right_modes)
            else:
                if right_maximum > left_maximum:
                    retained, discarded = right_modes, left_modes
                else:
                    retained, discarded = left_modes, right_modes
                for mode_id in discarded:
                    retained_modes[mode_id] = False
                component_mode_ids[left_root] = retained
        elif right_modes:
            component_mode_ids[left_root] = right_modes
        component_mode_ids[right_root] = None
        return left_root

    def component(root: int) -> SignalInterval:
        root = find(root)
        return SignalInterval(
            chrom=intervals[root].chrom,
            start=starts[root],
            end=ends[root],
            maximum_cpm=maxima[root],
            weighted_signal=weighted_signals[root],
            signal_bases=signal_bases[root],
        )

    activation_index = 0

    for threshold in thresholds:
        activated: list[int] = []
        while (
            activation_index < len(activation_order)
            and intervals[activation_order[activation_index]].maximum_cpm == threshold
        ):
            index = activation_order[activation_index]
            active[index] = True
            activated.append(index)
            activation_index += 1
        for index in activated:
            if (
                index > 0
                and active[index - 1]
                and intervals[index - 1].chrom == intervals[index].chrom
                and intervals[index].start <= intervals[index - 1].end + merge_gap_bp
            ):
                union(index - 1, index, threshold)
            if (
                index + 1 < len(intervals)
                and active[index + 1]
                and intervals[index + 1].chrom == intervals[index].chrom
                and intervals[index + 1].start <= intervals[index].end + merge_gap_bp
            ):
                union(index, index + 1, threshold)

        changed = sorted(
            {find(index) for index in activated},
            key=lambda root: (
                intervals[root].chrom,
                starts[root],
                ends[root],
            ),
        )
        for root in changed:
            interval = component(root)
            qualifies = (
                minimum_length <= interval.end - interval.start <= maximum_length
                and interval.mean_cpm >= minimum_mean_cpm
            )
            modes = component_mode_ids[root]
            if not qualifies or (modes and len(modes) > 1):
                continue
            snapshot = replace(interval, selection_cutoff_cpm=threshold)
            if modes:
                mode_snapshots[modes[0]] = snapshot
            else:
                component_mode_ids[root] = [len(mode_snapshots)]
                mode_snapshots.append(snapshot)
                retained_modes.append(True)

    final_roots = [
        root
        for root in sorted(
            {find(index) for index in range(len(intervals))},
            key=lambda root: (
                intervals[root].chrom,
                starts[root],
                ends[root],
            ),
        )
        if ends[root] - starts[root] >= minimum_length
    ]
    final_threshold_intervals = [component(root) for root in final_roots]
    excluded = [
        interval
        for root, interval in zip(final_roots, final_threshold_intervals)
        if not component_mode_ids[root]
    ]
    return [
        interval
        for interval, retained in zip(mode_snapshots, retained_modes)
        if retained
    ], excluded, thresholds


def _write_intervals(
    path: Path,
    intervals: list[SignalInterval],
    *,
    maximum_signal: float,
    chromosome_order: dict[str, int],
) -> None:
    intervals.sort(
        key=lambda item: (chromosome_order[item.chrom], item.start, item.end)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for index, interval in enumerate(intervals, start=1):
            score = round(1000 * interval.maximum_cpm / maximum_signal)
            handle.write(
                f"{interval.chrom}\t{interval.start}\t{interval.end}\t"
                f"CPM_refined_{index}\t{score}\t.\t{interval.mean_cpm:.6g}\t"
                f"{interval.maximum_cpm:.6g}\t{interval.selection_cutoff_cpm:.6g}\n"
            )
    os.replace(temporary, path)


def refine_cpm_bigwig(
    *,
    peaks: Path,
    signal_bigwig: Path,
    output: Path,
    excluded: Path,
    stats: Path,
    merge_gap_bp: int = 1,
    minimum_length: int = 50,
    maximum_length: int = 400,
    minimum_mean_cpm: float = 0.0,
    minimum_mode_prominence: float = 0.25,
) -> dict[str, int | float | str]:
    """Refine MACS3 peaks with positive CPM bins fully contained by the peaks."""

    import pyBigWig

    peak_index, candidate_intervals, peak_count = _read_peaks(peaks)
    bigwig = pyBigWig.open(str(signal_bigwig))
    if bigwig is None:
        raise ValueError(f"Could not open BigWig: {signal_bigwig}")
    contained: dict[tuple[str, int, int], SignalInterval] = {}
    try:
        chromosome_sizes = bigwig.chroms()
        chromosome_order = {chrom: index for index, chrom in enumerate(chromosome_sizes)}
        for chrom, chrom_peaks in candidate_intervals.items():
            if chrom not in chromosome_sizes:
                raise ValueError(f"Peak chromosome is absent from BigWig: {chrom}")
            for peak_start, peak_end in chrom_peaks:
                for start, end, value in bigwig.intervals(
                    chrom, peak_start, peak_end
                ) or ():
                    signal = float(value)
                    if signal <= 0 or not math.isfinite(signal):
                        continue
                    if peak_index.contains(chrom, start, end):
                        contained[(chrom, start, end)] = SignalInterval.from_bigwig(
                            chrom, start, end, signal
                        )
    finally:
        bigwig.close()

    intervals = sorted(
        contained.values(),
        key=lambda item: (chromosome_order[item.chrom], item.start, item.end),
    )
    refined, excluded_intervals, thresholds = progressively_refine_cpm(
        intervals,
        merge_gap_bp=merge_gap_bp,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        minimum_mean_cpm=minimum_mean_cpm,
        minimum_mode_prominence=minimum_mode_prominence,
    )
    maximum_signal = max((item.maximum_cpm for item in intervals), default=1.0)
    _write_intervals(
        output,
        refined,
        maximum_signal=maximum_signal,
        chromosome_order=chromosome_order,
    )
    _write_intervals(
        excluded,
        excluded_intervals,
        maximum_signal=maximum_signal,
        chromosome_order=chromosome_order,
    )

    widths = [item.end - item.start for item in refined]
    if peak_count == 0:
        status = "no_candidate_peaks"
    elif not intervals:
        status = "no_contained_positive_signal"
    else:
        status = "ok"
    metrics: dict[str, int | float | str] = {
        "status": status,
        "candidate_macs3_peaks": peak_count,
        "refinement_algorithm": REFINEMENT_ALGORITHM,
        "contained_positive_signal_intervals": len(intervals),
        "observed_positive_cpm_thresholds": len(thresholds),
        "minimum_positive_cpm": min(thresholds, default=0.0),
        "maximum_cpm": max(thresholds, default=0.0),
        "refined_intervals": len(refined),
        "excluded_intervals": len(excluded_intervals),
        "minimum_width": min(widths) if widths else 0,
        "median_width": float(median(widths)) if widths else 0.0,
        "maximum_width": max(widths) if widths else 0,
        "minimum_mean_cpm": minimum_mean_cpm,
        "minimum_mode_prominence": minimum_mode_prominence,
        "merge_gap_bp": merge_gap_bp,
        "minimum_length": minimum_length,
        "maximum_length_parameter": maximum_length,
        "threshold_rule": (
            "all observed positive CPM levels, descending; signal >= cutoff; "
            "merge modes below the configured relative saddle prominence; "
            "freeze boundaries before retained prominent mode mergers"
        ),
        "output_columns": "BED6,mean_cpm,max_cpm,selection_cutoff_cpm",
    }
    stats.parent.mkdir(parents=True, exist_ok=True)
    temporary_stats = stats.with_suffix(stats.suffix + ".tmp")
    temporary_stats.write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary_stats, stats)
    return metrics
