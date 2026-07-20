"""Progressively refine MACS3 peaks using a CPM BigWig signal."""

from __future__ import annotations

import bisect
import json
import math
import os
from dataclasses import dataclass, replace
from pathlib import Path
from statistics import median


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


class NonOverlappingIndex:
    def __init__(self):
        self.intervals: dict[str, list[tuple[int, int]]] = {}

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        intervals = self.intervals.setdefault(chrom, [])
        index = bisect.bisect_left(intervals, (start, -1))
        if index and intervals[index - 1][1] > start:
            return True
        return index < len(intervals) and intervals[index][0] < end

    def add(self, chrom: str, start: int, end: int) -> None:
        bisect.insort(self.intervals.setdefault(chrom, []), (start, end))


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


def _merge(
    intervals: list[SignalInterval], *, merge_gap_bp: int, minimum_length: int
) -> list[SignalInterval]:
    merged: list[SignalInterval] = []
    for interval in intervals:
        if (
            not merged
            or interval.chrom != merged[-1].chrom
            or interval.start > merged[-1].end + merge_gap_bp
        ):
            merged.append(interval)
            continue
        previous = merged[-1]
        merged[-1] = SignalInterval(
            chrom=previous.chrom,
            start=previous.start,
            end=max(previous.end, interval.end),
            maximum_cpm=max(previous.maximum_cpm, interval.maximum_cpm),
            weighted_signal=previous.weighted_signal + interval.weighted_signal,
            signal_bases=previous.signal_bases + interval.signal_bases,
        )
    return [item for item in merged if item.end - item.start >= minimum_length]


def progressively_refine_cpm(
    intervals: list[SignalInterval],
    *,
    merge_gap_bp: int = 1,
    minimum_length: int = 50,
    maximum_length: int = 400,
    minimum_mean_cpm: float = 0.0,
) -> tuple[list[SignalInterval], list[SignalInterval], list[float]]:
    """Apply the rDHS geometric selection over all observed CPM cutoffs."""

    if (
        merge_gap_bp < 0
        or minimum_length < 1
        or maximum_length < minimum_length
        or not math.isfinite(minimum_mean_cpm)
        or minimum_mean_cpm < 0
    ):
        raise ValueError("Invalid merge-gap or interval-length parameters")
    if not intervals:
        return [], [], []
    thresholds = sorted({interval.maximum_cpm for interval in intervals})
    active = intervals
    selected: list[SignalInterval] = []
    selected_index = NonOverlappingIndex()
    final_threshold_intervals: list[SignalInterval] = []

    for threshold_index, threshold in enumerate(thresholds):
        if threshold_index:
            active = [item for item in active if item.maximum_cpm >= threshold]
        threshold_intervals = _merge(
            active,
            merge_gap_bp=merge_gap_bp,
            minimum_length=minimum_length,
        )
        if threshold_index == len(thresholds) - 1:
            final_threshold_intervals = threshold_intervals
        for interval in threshold_intervals:
            if interval.mean_cpm < minimum_mean_cpm:
                continue
            if interval.end - interval.start > maximum_length:
                continue
            if selected_index.overlaps(interval.chrom, interval.start, interval.end):
                continue
            selected.append(replace(interval, selection_cutoff_cpm=threshold))
            selected_index.add(interval.chrom, interval.start, interval.end)

    excluded = [
        interval
        for interval in final_threshold_intervals
        if not selected_index.overlaps(interval.chrom, interval.start, interval.end)
    ]
    return selected, excluded, thresholds


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
        "contained_positive_signal_intervals": len(intervals),
        "observed_positive_cpm_thresholds": len(thresholds),
        "minimum_positive_cpm": thresholds[0] if thresholds else 0.0,
        "maximum_cpm": thresholds[-1] if thresholds else 0.0,
        "refined_intervals": len(refined),
        "excluded_intervals": len(excluded_intervals),
        "minimum_width": min(widths) if widths else 0,
        "median_width": float(median(widths)) if widths else 0.0,
        "maximum_width": max(widths) if widths else 0,
        "minimum_mean_cpm": minimum_mean_cpm,
        "merge_gap_bp": merge_gap_bp,
        "minimum_length": minimum_length,
        "maximum_length_parameter": maximum_length,
        "threshold_rule": "all observed positive CPM levels, ascending; signal >= cutoff",
        "output_columns": "BED6,mean_cpm,max_cpm,selection_cutoff_cpm",
    }
    stats.parent.mkdir(parents=True, exist_ok=True)
    temporary_stats = stats.with_suffix(stats.suffix + ".tmp")
    temporary_stats.write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary_stats, stats)
    return metrics
