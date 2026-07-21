"""Build reproducible condition ATAC peaks and a bounded cross-condition atlas."""

from __future__ import annotations

import bisect
import csv
import json
import math
import os
import re
import statistics
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, TextIO

from .accessions import AcquisitionError


SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
CONDITION_MAP_COLUMNS = {"condition_id", "condition_label", "sample_id"}


def _read_delimited_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the small atlas CSV/TSV map without workflow-only dependencies."""

    text = path.read_text(encoding="utf-8-sig")
    if path.suffix.lower() == ".csv":
        delimiter = ","
    elif path.suffix.lower() in {".tsv", ".tab"}:
        delimiter = "\t"
    else:
        try:
            delimiter = csv.Sniffer().sniff(text[:8192], delimiters=",\t").delimiter
        except csv.Error as exc:
            raise AcquisitionError(
                f"Could not determine whether {path} is CSV or TSV"
            ) from exc
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    fieldnames = [field.strip() for field in (reader.fieldnames or []) if field]
    if not fieldnames:
        raise AcquisitionError(f"ATAC atlas condition map {path} has no header")
    if len(fieldnames) != len(set(fieldnames)):
        raise AcquisitionError(
            f"ATAC atlas condition map {path} has duplicate column names"
        )
    rows = [
        {
            str(key).strip(): (value or "").strip()
            for key, value in row.items()
            if key is not None
        }
        for row in reader
    ]
    return fieldnames, [row for row in rows if any(row.values())]


@dataclass(frozen=True)
class ConditionSpec:
    condition_id: str
    label: str
    samples: tuple[str, ...]


@dataclass(frozen=True)
class RefinedPeak:
    chrom: str
    start: int
    end: int
    name: str
    score: int
    strand: str
    mean_cpm: float
    maximum_cpm: float
    selection_cutoff_cpm: float
    condition_id: str = ""
    support_n: int = 0
    replicate_n: int = 0
    support_fraction: float = 0.0
    supporting_samples: str = ""


@dataclass(frozen=True)
class AtlasCandidate:
    peak: RefinedPeak
    summit: int
    start: int
    end: int
    priority: float


@dataclass
class VariableBoundary:
    atlas_id: str
    chrom: str
    start: int
    end: int
    condition_votes: int
    collision_trimmed: bool = False


@dataclass(frozen=True)
class AtlasMembershipSource:
    atlas_id: str
    condition_id: str
    peak_id: str
    chrom: str
    start: int
    end: int
    summit: int
    atlas_center: int
    support_fraction: float
    priority: float


@dataclass(frozen=True)
class SupportSegment:
    chrom: str
    start: int
    end: int
    support: int


@dataclass(frozen=True)
class FwhmBoundary:
    atlas_id: str
    chrom: str
    anchor_start: int
    anchor_end: int
    start: int
    end: int
    summit: int
    maximum_support: int
    half_maximum_support: int
    component_maximum_support: int
    maximum_plateaus: int
    bridged_higher_peak: bool
    status: str


@dataclass(frozen=True)
class LocalModeBoundary:
    atlas_id: str
    chrom: str
    anchor_start: int
    anchor_end: int
    start: int
    end: int
    summit: int
    peak_support: int
    left_base_support: int
    right_base_support: int
    prominence: int
    half_prominence_support: int
    local_mode_n: int
    fallback_nonmaximum: bool
    bridged_higher_peak: bool
    status: str


class IntervalUnion:
    """Merged interval index used for exact covered-base calculations."""

    def __init__(self, intervals: Iterable[tuple[str, int, int]]):
        grouped: dict[str, list[tuple[int, int]]] = {}
        for chrom, start, end in intervals:
            grouped.setdefault(chrom, []).append((start, end))
        self.intervals: dict[str, list[tuple[int, int]]] = {}
        self.starts: dict[str, list[int]] = {}
        for chrom, values in grouped.items():
            merged: list[tuple[int, int]] = []
            for start, end in sorted(values):
                if not merged or start > merged[-1][1]:
                    merged.append((start, end))
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            self.intervals[chrom] = merged
            self.starts[chrom] = [start for start, _end in merged]

    def covered_bases(self, chrom: str, start: int, end: int) -> int:
        intervals = self.intervals.get(chrom, [])
        if not intervals:
            return 0
        index = max(0, bisect.bisect_right(self.starts[chrom], start) - 1)
        covered = 0
        while index < len(intervals) and intervals[index][0] < end:
            interval_start, interval_end = intervals[index]
            covered += max(0, min(end, interval_end) - max(start, interval_start))
            index += 1
        return covered


def read_condition_map(
    path: Path,
    *,
    sample_ids: Iterable[str],
    minimum_replicates: int,
) -> list[ConditionSpec]:
    """Read and validate an explicit sample-to-condition atlas map."""

    fieldnames, rows = _read_delimited_rows(path)
    missing_columns = sorted(CONDITION_MAP_COLUMNS - set(fieldnames))
    if missing_columns:
        raise AcquisitionError(
            f"ATAC atlas condition map {path} is missing columns: "
            + ", ".join(missing_columns)
        )
    if minimum_replicates < 2:
        raise AcquisitionError("ATAC atlas minimum_replicates must be at least 2")

    expected = set(sample_ids)
    assigned: dict[str, str] = {}
    labels: dict[str, str] = {}
    samples_by_condition: dict[str, list[str]] = {}
    for line, row in enumerate(rows, start=2):
        condition = row["condition_id"]
        label = row["condition_label"]
        sample = row["sample_id"]
        if not SAFE_ID_RE.fullmatch(condition):
            raise AcquisitionError(f"{path}:{line}: invalid condition_id {condition!r}")
        if not label:
            raise AcquisitionError(f"{path}:{line}: condition_label is blank")
        if sample not in expected:
            raise AcquisitionError(f"{path}:{line}: unknown ATAC sample_id {sample!r}")
        if sample in assigned:
            raise AcquisitionError(
                f"{path}:{line}: sample_id {sample!r} is assigned more than once"
            )
        if condition in labels and labels[condition] != label:
            raise AcquisitionError(
                f"{path}:{line}: condition {condition!r} has inconsistent labels"
            )
        assigned[sample] = condition
        labels[condition] = label
        samples_by_condition.setdefault(condition, []).append(sample)

    missing_samples = sorted(expected - set(assigned))
    if missing_samples:
        raise AcquisitionError(
            f"ATAC atlas condition map {path} does not assign: "
            + ", ".join(missing_samples)
        )
    undersized = [
        condition
        for condition, samples in samples_by_condition.items()
        if len(samples) < minimum_replicates
    ]
    if undersized:
        raise AcquisitionError(
            "ATAC atlas conditions below minimum_replicates: "
            + ", ".join(undersized)
        )
    return [
        ConditionSpec(condition, labels[condition], tuple(samples))
        for condition, samples in samples_by_condition.items()
    ]


def _read_refined_peaks(path: Path, *, condition_bed: bool = False) -> list[RefinedPeak]:
    peaks: list[RefinedPeak] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            minimum_fields = 14 if condition_bed else 9
            if len(fields) < minimum_fields:
                raise ValueError(
                    f"{path}:{line_number}: expected at least {minimum_fields} columns"
                )
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            numeric = [float(fields[index]) for index in (6, 7, 8)]
            if any(not math.isfinite(value) or value < 0 for value in numeric):
                raise ValueError(f"{path}:{line_number}: invalid CPM value")
            peaks.append(
                RefinedPeak(
                    chrom=fields[0],
                    start=start,
                    end=end,
                    name=fields[3],
                    score=int(fields[4]),
                    strand=fields[5],
                    mean_cpm=numeric[0],
                    maximum_cpm=numeric[1],
                    selection_cutoff_cpm=numeric[2],
                    condition_id=fields[9] if condition_bed else "",
                    support_n=int(fields[10]) if condition_bed else 0,
                    replicate_n=int(fields[11]) if condition_bed else 0,
                    support_fraction=float(fields[12]) if condition_bed else 0.0,
                    supporting_samples=fields[13] if condition_bed else "",
                )
            )
    return peaks


def _atomic_writer(path: Path) -> tuple[TextIO, Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix=f".{path.name}.",
        dir=path.parent,
        delete=False,
    )
    return handle, Path(handle.name)


def _finish_atomic(handle: TextIO, temporary: Path, destination: Path) -> None:
    handle.close()
    os.replace(temporary, destination)


def _write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    handle, temporary = _atomic_writer(path)
    try:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        _finish_atomic(handle, temporary, path)
    finally:
        if not handle.closed:
            handle.close()
        temporary.unlink(missing_ok=True)


def build_condition_consensus(
    *,
    condition_id: str,
    pooled_peaks: Path,
    replicate_peaks: dict[str, Path],
    output_bed: Path,
    support_tsv: Path,
    stats_json: Path,
    minimum_replicates: int = 2,
    overlap_fraction: float = 0.5,
) -> dict[str, Any]:
    """Retain pooled refined peaks covered by enough biological replicates."""

    if not SAFE_ID_RE.fullmatch(condition_id):
        raise ValueError(f"Invalid condition ID: {condition_id!r}")
    if minimum_replicates < 2 or minimum_replicates > len(replicate_peaks):
        raise ValueError("minimum_replicates is inconsistent with replicate inputs")
    if not 0 < overlap_fraction <= 1:
        raise ValueError("overlap_fraction must be in (0, 1]")

    pooled = _read_refined_peaks(pooled_peaks)
    replicate_names = list(replicate_peaks)
    indexes = {
        sample: IntervalUnion(
            (peak.chrom, peak.start, peak.end)
            for peak in _read_refined_peaks(path)
        )
        for sample, path in replicate_peaks.items()
    }

    bed_handle, temporary_bed = _atomic_writer(output_bed)
    support_handle, temporary_support = _atomic_writer(support_tsv)
    retained_count = 0
    try:
        support_handle.write(
            "condition_id\tpooled_peak_id\tretained\tsupport_n\treplicate_n"
            "\tsupport_fraction\t"
            + "\t".join(f"{sample}_coverage_fraction" for sample in replicate_names)
            + "\n"
        )
        for peak in pooled:
            length = peak.end - peak.start
            fractions = {
                sample: indexes[sample].covered_bases(peak.chrom, peak.start, peak.end)
                / length
                for sample in replicate_names
            }
            supporting = [
                sample for sample in replicate_names if fractions[sample] >= overlap_fraction
            ]
            support_n = len(supporting)
            retained = support_n >= minimum_replicates
            support_fraction = support_n / len(replicate_names)
            support_handle.write(
                f"{condition_id}\t{peak.name}\t{int(retained)}\t{support_n}\t"
                f"{len(replicate_names)}\t{support_fraction:.6g}\t"
                + "\t".join(f"{fractions[sample]:.6g}" for sample in replicate_names)
                + "\n"
            )
            if not retained:
                continue
            retained_count += 1
            bed_handle.write(
                f"{peak.chrom}\t{peak.start}\t{peak.end}\t{peak.name}\t"
                f"{peak.score}\t{peak.strand}\t{peak.mean_cpm:.6g}\t"
                f"{peak.maximum_cpm:.6g}\t{peak.selection_cutoff_cpm:.6g}\t"
                f"{condition_id}\t{support_n}\t{len(replicate_names)}\t"
                f"{support_fraction:.6g}\t{','.join(supporting)}\n"
            )
        _finish_atomic(bed_handle, temporary_bed, output_bed)
        _finish_atomic(support_handle, temporary_support, support_tsv)
    finally:
        if not bed_handle.closed:
            bed_handle.close()
        if not support_handle.closed:
            support_handle.close()
        temporary_bed.unlink(missing_ok=True)
        temporary_support.unlink(missing_ok=True)

    metrics: dict[str, Any] = {
        "status": "ok" if retained_count else "no_reproducible_condition_peaks",
        "condition_id": condition_id,
        "pooled_refined_peaks": len(pooled),
        "retained_consensus_peaks": retained_count,
        "replicate_n": len(replicate_names),
        "minimum_replicates": minimum_replicates,
        "minimum_pooled_peak_coverage_fraction": overlap_fraction,
        "replicates": replicate_names,
    }
    _write_json_atomic(stats_json, metrics)
    return metrics


def _signal_summary(bigwig, chrom: str, start: int, end: int) -> tuple[float, float, int]:
    weighted = 0.0
    maximum = 0.0
    maximum_start = start
    maximum_end = end
    for interval_start, interval_end, value in bigwig.intervals(chrom, start, end) or ():
        signal = float(value)
        if not math.isfinite(signal) or signal < 0:
            continue
        overlap_start = max(start, int(interval_start))
        overlap_end = min(end, int(interval_end))
        if overlap_end <= overlap_start:
            continue
        weighted += signal * (overlap_end - overlap_start)
        if signal > maximum:
            maximum = signal
            maximum_start = overlap_start
            maximum_end = overlap_end
    summit = (maximum_start + maximum_end) // 2 if maximum > 0 else (start + end) // 2
    return weighted / (end - start), maximum, summit


def _fixed_window(center: int, width: int, chromosome_size: int) -> tuple[int, int]:
    if chromosome_size <= width:
        return 0, chromosome_size
    start = max(0, center - width // 2)
    end = start + width
    if end > chromosome_size:
        end = chromosome_size
        start = end - width
    return start, end


def _percentiles(values: list[float]) -> list[float]:
    ordered = sorted(values)
    return [bisect.bisect_right(ordered, value) / len(ordered) for value in values]


def _overlap(first: AtlasCandidate, second: AtlasCandidate) -> int:
    if first.peak.chrom != second.peak.chrom:
        return 0
    return max(0, min(first.end, second.end) - max(first.start, second.start))


def _candidate_sort_key(
    item: AtlasCandidate, chromosome_order: dict[str, int]
) -> tuple[Any, ...]:
    return (
        -item.priority,
        -item.peak.support_fraction,
        -item.peak.maximum_cpm,
        chromosome_order[item.peak.chrom],
        item.start,
        item.peak.condition_id,
        item.peak.name,
    )


def _select_nonoverlapping(
    candidates: list[AtlasCandidate], chromosome_order: dict[str, int]
) -> tuple[list[AtlasCandidate], list[tuple[AtlasCandidate, int]]]:
    ordered = sorted(
        candidates,
        key=lambda item: _candidate_sort_key(item, chromosome_order),
    )
    selected: list[AtlasCandidate] = []
    selected_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    assignments: list[tuple[AtlasCandidate, int]] = []
    for candidate in ordered:
        intervals = selected_by_chrom.setdefault(candidate.peak.chrom, [])
        position = bisect.bisect_left(intervals, (candidate.start, -1, -1))
        possible = []
        if position:
            possible.append(intervals[position - 1][2])
        if position < len(intervals):
            possible.append(intervals[position][2])
        overlapping = [index for index in possible if _overlap(candidate, selected[index]) > 0]
        if overlapping:
            selected_index = min(
                overlapping,
                key=lambda index: (
                    -_overlap(candidate, selected[index]),
                    -selected[index].priority,
                    selected[index].start,
                ),
            )
        else:
            selected_index = len(selected)
            selected.append(candidate)
            bisect.insort(intervals, (candidate.start, candidate.end, selected_index))
        assignments.append((candidate, selected_index))
    return selected, assignments


def _select_narrow_source_first(
    candidates: list[AtlasCandidate], chromosome_order: dict[str, int]
) -> tuple[list[AtlasCandidate], list[tuple[AtlasCandidate, int]]]:
    """Select narrow source DHSs first and annotate bridges to every anchor."""

    ordered = sorted(
        candidates,
        key=lambda item: (
            item.peak.end - item.peak.start,
            _candidate_sort_key(item, chromosome_order),
        ),
    )
    selected: list[AtlasCandidate] = []
    selected_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for candidate in ordered:
        intervals = selected_by_chrom.setdefault(candidate.peak.chrom, [])
        position = bisect.bisect_left(intervals, (candidate.start, -1, -1))
        overlaps_previous = (
            position > 0 and intervals[position - 1][1] > candidate.start
        )
        overlaps_following = (
            position < len(intervals) and intervals[position][0] < candidate.end
        )
        if overlaps_previous or overlaps_following:
            continue
        selected_index = len(selected)
        selected.append(candidate)
        bisect.insort(intervals, (candidate.start, candidate.end, selected_index))

    assignments: list[tuple[AtlasCandidate, int]] = []
    for candidate in ordered:
        intervals = selected_by_chrom[candidate.peak.chrom]
        position = bisect.bisect_left(intervals, (candidate.start, -1, -1))
        matching: list[int] = []
        if position > 0:
            start, end, selected_index = intervals[position - 1]
            if end > candidate.start and start < candidate.end:
                matching.append(selected_index)
        while position < len(intervals) and intervals[position][0] < candidate.end:
            start, end, selected_index = intervals[position]
            if end > candidate.start:
                matching.append(selected_index)
            position += 1
        if not matching:
            raise ValueError("Narrow-first selection left a candidate unassigned")
        assignments.extend((candidate, selected_index) for selected_index in matching)
    return selected, assignments


def _dhs_seed_match(candidate: AtlasCandidate, seed: AtlasCandidate) -> bool:
    if candidate.peak.chrom != seed.peak.chrom:
        return False
    overlap = max(
        0,
        min(candidate.peak.end, seed.peak.end)
        - max(candidate.peak.start, seed.peak.start),
    )
    if overlap == 0:
        return False
    return (
        seed.peak.start <= candidate.summit < seed.peak.end
        or candidate.peak.start <= seed.summit < candidate.peak.end
    )


def _select_dhs_seeds(
    candidates: list[AtlasCandidate], chromosome_order: dict[str, int]
) -> tuple[list[AtlasCandidate], list[tuple[AtlasCandidate, int]]]:
    """Group only DHSs that directly match a stronger selected seed."""

    ordered = sorted(
        candidates,
        key=lambda item: _candidate_sort_key(item, chromosome_order),
    )
    maximum_width = max(
        (candidate.peak.end - candidate.peak.start for candidate in candidates),
        default=0,
    )
    selected: list[AtlasCandidate] = []
    selected_by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    assignments: list[tuple[AtlasCandidate, int]] = []
    for candidate in ordered:
        intervals = selected_by_chrom.setdefault(candidate.peak.chrom, [])
        lower = bisect.bisect_left(
            intervals,
            (candidate.peak.start - maximum_width, -1, -1),
        )
        upper = bisect.bisect_left(intervals, (candidate.peak.end, -1, -1))
        matching = [
            selected_index
            for _start, _end, selected_index in intervals[lower:upper]
            if _dhs_seed_match(candidate, selected[selected_index])
        ]
        if matching:
            selected_index = min(
                matching,
                key=lambda index: (
                    -max(
                        0,
                        min(candidate.peak.end, selected[index].peak.end)
                        - max(candidate.peak.start, selected[index].peak.start),
                    ),
                    abs(candidate.summit - selected[index].summit),
                    _candidate_sort_key(selected[index], chromosome_order),
                ),
            )
        else:
            selected_index = len(selected)
            selected.append(candidate)
            bisect.insort(
                intervals,
                (candidate.peak.start, candidate.peak.end, selected_index),
            )
        assignments.append((candidate, selected_index))
    return selected, assignments


def _integer_median(values: list[int], *, round_up: bool = False) -> int:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    total = ordered[middle - 1] + ordered[middle]
    return (total + int(round_up)) // 2


def _build_variable_boundaries(
    selected_by_atlas_id: dict[str, AtlasCandidate],
    assignments_by_atlas: list[tuple[str, AtlasCandidate]],
) -> tuple[list[VariableBoundary], int, int]:
    """Take one unweighted boundary vote per condition and trim neighbor collisions."""

    grouped: dict[str, list[AtlasCandidate]] = {
        atlas_id: [] for atlas_id in selected_by_atlas_id
    }
    for atlas_id, candidate in assignments_by_atlas:
        grouped[atlas_id].append(candidate)

    boundaries: list[VariableBoundary] = []
    for atlas_id, representative in selected_by_atlas_id.items():
        candidates_by_condition: dict[str, list[AtlasCandidate]] = {}
        for candidate in grouped[atlas_id]:
            candidates_by_condition.setdefault(candidate.peak.condition_id, []).append(
                candidate
            )
        votes = [
            min(
                candidates,
                key=lambda item: (
                    -item.priority,
                    -item.peak.support_fraction,
                    -item.peak.maximum_cpm,
                    item.peak.start,
                    item.peak.end,
                    item.peak.name,
                ),
            )
            for candidates in candidates_by_condition.values()
        ]
        boundaries.append(
            VariableBoundary(
                atlas_id=atlas_id,
                chrom=representative.peak.chrom,
                start=_integer_median([vote.peak.start for vote in votes]),
                end=_integer_median(
                    [vote.peak.end for vote in votes], round_up=True
                ),
                condition_votes=len(votes),
            )
        )

    trimmed_ids: set[str] = set()
    for previous, current in zip(boundaries, boundaries[1:]):
        if previous.chrom != current.chrom or previous.end <= current.start:
            continue
        split = (previous.end + current.start) // 2
        if split - previous.start < 50 or current.end - split < 50:
            continue
        previous.end = split
        current.start = split
        previous.collision_trimmed = True
        current.collision_trimmed = True
        trimmed_ids.update((previous.atlas_id, current.atlas_id))
    if any(boundary.end <= boundary.start for boundary in boundaries):
        raise ValueError("Variable-boundary collision resolution produced an empty peak")
    overlapping_pairs = sum(
        previous.chrom == current.chrom and previous.end > current.start
        for previous, current in zip(boundaries, boundaries[1:])
    )
    return boundaries, len(trimmed_ids), overlapping_pairs


def build_global_atlas(
    *,
    conditions: list[tuple[str, Path, Path]],
    output_bed: Path,
    variable_bed: Path,
    membership_tsv: Path,
    presence_tsv: Path,
    coverage_tsv: Path,
    mean_cpm_tsv: Path,
    maximum_cpm_tsv: Path,
    stats_json: Path,
    peak_width: int = 250,
    grouping_method: str = "fixed_window",
) -> dict[str, Any]:
    """Iteratively select bounded peaks and retain condition coverage matrices."""

    import pyBigWig

    if peak_width < 1:
        raise ValueError("peak_width must be positive")
    if grouping_method not in {
        "fixed_window",
        "fixed_window_narrow_first",
        "dhs_seed",
    }:
        raise ValueError(f"Unsupported grouping method: {grouping_method}")
    if not conditions:
        raise ValueError("At least one condition is required")
    condition_ids = [condition for condition, _bed, _bigwig in conditions]
    if len(condition_ids) != len(set(condition_ids)):
        raise ValueError("Condition IDs must be unique")

    peaks_by_condition: dict[str, list[RefinedPeak]] = {}
    bigwigs: dict[str, Any] = {}
    chromosome_sizes: dict[str, int] | None = None
    try:
        for condition, bed_path, bigwig_path in conditions:
            peaks = _read_refined_peaks(bed_path, condition_bed=True)
            if any(peak.condition_id != condition for peak in peaks):
                raise ValueError(f"Condition BED does not match {condition!r}: {bed_path}")
            peaks_by_condition[condition] = peaks
            bigwig = pyBigWig.open(str(bigwig_path))
            if bigwig is None:
                raise ValueError(f"Could not open BigWig: {bigwig_path}")
            sizes = {str(chrom): int(size) for chrom, size in bigwig.chroms().items()}
            if chromosome_sizes is None:
                chromosome_sizes = sizes
            elif sizes != chromosome_sizes:
                raise ValueError("Condition BigWigs have different chromosome headers")
            bigwigs[condition] = bigwig

        assert chromosome_sizes is not None
        chromosome_order = {chrom: index for index, chrom in enumerate(chromosome_sizes)}
        candidates: list[AtlasCandidate] = []
        for condition in condition_ids:
            peaks = peaks_by_condition[condition]
            summaries = []
            for peak in peaks:
                if peak.chrom not in chromosome_sizes:
                    raise ValueError(f"Peak chromosome is absent from BigWig: {peak.chrom}")
                summaries.append(
                    _signal_summary(bigwigs[condition], peak.chrom, peak.start, peak.end)
                )
            priorities = _percentiles([maximum for _mean, maximum, _summit in summaries])
            for peak, (_mean, _maximum, summit), priority in zip(
                peaks, summaries, priorities, strict=True
            ):
                start, end = _fixed_window(
                    summit,
                    peak_width,
                    chromosome_sizes[peak.chrom],
                )
                candidates.append(AtlasCandidate(peak, summit, start, end, priority))

        if grouping_method == "fixed_window":
            selected, assignments = _select_nonoverlapping(
                candidates, chromosome_order
            )
            method = "fixed-width iterative overlap removal"
            grouping_rule = "direct overlap between summit-centered fixed windows"
        elif grouping_method == "fixed_window_narrow_first":
            selected, assignments = _select_narrow_source_first(
                candidates, chromosome_order
            )
            method = "narrow-source-first fixed-width overlap removal"
            grouping_rule = (
                "select by original source DHS width before signal priority; annotate a bridging source to every retained fixed window it overlaps"
            )
        else:
            selected, assignments = _select_dhs_seeds(candidates, chromosome_order)
            method = "direct seed DHS grouping"
            grouping_rule = (
                "source DHS intervals overlap and either source summit lies inside the other interval"
            )
        genomic_selected = sorted(
            enumerate(selected),
            key=lambda item: (
                chromosome_order[item[1].peak.chrom],
                item[1].start,
                item[1].end,
            ),
        )
        atlas_id_by_selected = {
            original_index: f"ATAC_ATLAS_{index:06d}"
            for index, (original_index, _candidate) in enumerate(genomic_selected, start=1)
        }
        selected_by_atlas_id = {
            atlas_id_by_selected[original_index]: candidate
            for original_index, candidate in genomic_selected
        }
        assigned_conditions: dict[str, set[str]] = {
            atlas_id: set() for atlas_id in selected_by_atlas_id
        }
        assignments_by_atlas: list[tuple[str, AtlasCandidate]] = []
        assignment_counts: dict[tuple[str, str, str, int, int], int] = {}
        for candidate, selected_index in assignments:
            atlas_id = atlas_id_by_selected[selected_index]
            assigned_conditions[atlas_id].add(candidate.peak.condition_id)
            assignments_by_atlas.append((atlas_id, candidate))
            source_key = (
                candidate.peak.condition_id,
                candidate.peak.name,
                candidate.peak.chrom,
                candidate.peak.start,
                candidate.peak.end,
            )
            assignment_counts[source_key] = assignment_counts.get(source_key, 0) + 1

        (
            variable_boundaries,
            collision_trimmed_peaks,
            variable_overlapping_pairs,
        ) = _build_variable_boundaries(selected_by_atlas_id, assignments_by_atlas)

        condition_indexes = {
            condition: IntervalUnion(
                (peak.chrom, peak.start, peak.end)
                for peak in peaks_by_condition[condition]
            )
            for condition in condition_ids
        }

        destinations = [
            output_bed,
            variable_bed,
            membership_tsv,
            presence_tsv,
            coverage_tsv,
            mean_cpm_tsv,
            maximum_cpm_tsv,
        ]
        opened = [_atomic_writer(path) for path in destinations]
        handles = [item[0] for item in opened]
        try:
            (
                bed_handle,
                variable_handle,
                membership_handle,
                presence_handle,
                coverage_handle,
                mean_handle,
                max_handle,
            ) = handles
            membership_handle.write(
                "atlas_peak_id\tcondition_id\tsource_peak_id\tsource_chrom\t"
                "source_start\tsource_end\tsource_summit\tatlas_center\t"
                "summit_distance\tsource_support_n\tsource_replicate_n\t"
                "source_support_fraction\tsource_priority_percentile\t"
                "source_coverage_fraction\n"
            )
            matrix_header = "atlas_peak_id\t" + "\t".join(condition_ids) + "\n"
            presence_handle.write(matrix_header)
            coverage_handle.write(matrix_header)
            mean_handle.write(matrix_header)
            max_handle.write(matrix_header)

            for atlas_id, candidate in selected_by_atlas_id.items():
                score = round(1000 * candidate.priority)
                bed_handle.write(
                    f"{candidate.peak.chrom}\t{candidate.start}\t{candidate.end}\t"
                    f"{atlas_id}\t{score}\t.\n"
                )
                presence_handle.write(
                    atlas_id
                    + "\t"
                    + "\t".join(
                        "1" if condition in assigned_conditions[atlas_id] else "0"
                        for condition in condition_ids
                    )
                    + "\n"
                )
                length = candidate.end - candidate.start
                coverage_values = []
                mean_values = []
                maximum_values = []
                for condition in condition_ids:
                    coverage_values.append(
                        condition_indexes[condition].covered_bases(
                            candidate.peak.chrom, candidate.start, candidate.end
                        )
                        / length
                    )
                    mean_signal, maximum_signal, _summit = _signal_summary(
                        bigwigs[condition],
                        candidate.peak.chrom,
                        candidate.start,
                        candidate.end,
                    )
                    mean_values.append(mean_signal)
                    maximum_values.append(maximum_signal)
                coverage_handle.write(
                    atlas_id + "\t" + "\t".join(f"{value:.6g}" for value in coverage_values) + "\n"
                )
                mean_handle.write(
                    atlas_id + "\t" + "\t".join(f"{value:.6g}" for value in mean_values) + "\n"
                )
                max_handle.write(
                    atlas_id + "\t" + "\t".join(f"{value:.6g}" for value in maximum_values) + "\n"
                )

            for boundary in variable_boundaries:
                candidate = selected_by_atlas_id[boundary.atlas_id]
                score = round(1000 * candidate.priority)
                variable_handle.write(
                    f"{boundary.chrom}\t{boundary.start}\t{boundary.end}\t"
                    f"{boundary.atlas_id}\t{score}\t.\n"
                )

            for atlas_id, source in sorted(
                assignments_by_atlas,
                key=lambda item: (
                    item[0],
                    item[1].peak.condition_id,
                    item[1].peak.name,
                ),
            ):
                atlas = selected_by_atlas_id[atlas_id]
                atlas_length = atlas.end - atlas.start
                source_coverage = max(
                    0,
                    min(atlas.end, source.peak.end) - max(atlas.start, source.peak.start),
                ) / atlas_length
                atlas_center = (atlas.start + atlas.end) // 2
                membership_handle.write(
                    f"{atlas_id}\t{source.peak.condition_id}\t{source.peak.name}\t"
                    f"{source.peak.chrom}\t{source.peak.start}\t{source.peak.end}\t"
                    f"{source.summit}\t{atlas_center}\t"
                    f"{abs(source.summit - atlas_center)}\t{source.peak.support_n}\t"
                    f"{source.peak.replicate_n}\t{source.peak.support_fraction:.6g}\t"
                    f"{source.priority:.6g}\t{source_coverage:.6g}\n"
                )

            for handle, temporary, destination in zip(
                handles,
                [item[1] for item in opened],
                destinations,
                strict=True,
            ):
                _finish_atomic(handle, temporary, destination)
        finally:
            for handle, temporary in opened:
                if not handle.closed:
                    handle.close()
                temporary.unlink(missing_ok=True)

        metrics: dict[str, Any] = {
            "status": "ok" if selected else "no_condition_consensus_peaks",
            "method": method,
            "grouping_method": grouping_method,
            "grouping_rule": grouping_rule,
            "condition_n": len(condition_ids),
            "conditions": condition_ids,
            "source_condition_peaks": len(candidates),
            "source_membership_assignments": len(assignments_by_atlas),
            "multi_anchor_source_peaks": sum(
                count > 1 for count in assignment_counts.values()
            ),
            "atlas_peaks": len(selected),
            "peak_width": peak_width,
            "variable_boundary_method": (
                "unweighted median start/end of the strongest assigned peak per condition"
            ),
            "variable_boundary_collision_rule": (
                "neighboring overlaps split midway only when both intervals remain at least 50 bp"
            ),
            "variable_boundary_collision_trimmed_peaks": collision_trimmed_peaks,
            "variable_boundary_overlapping_pairs": variable_overlapping_pairs,
            "variable_boundary_min_width": min(
                (boundary.end - boundary.start for boundary in variable_boundaries),
                default=0,
            ),
            "variable_boundary_max_width": max(
                (boundary.end - boundary.start for boundary in variable_boundaries),
                default=0,
            ),
            "priority": (
                "original source DHS width ascending, then within-condition maximum-CPM percentile"
                if grouping_method == "fixed_window_narrow_first"
                else "within-condition percentile of maximum CPM"
            ),
            "presence_rule": (
                "source condition peak assigned during iterative overlap removal"
                if grouping_method == "fixed_window"
                else (
                    "source condition peak annotated to every overlapping retained fixed window"
                    if grouping_method == "fixed_window_narrow_first"
                    else "source condition peak assigned by direct seed DHS grouping"
                )
            ),
            "coverage_rule": "covered bases by union of condition consensus peaks / atlas width",
        }
        _write_json_atomic(stats_json, metrics)
        return metrics
    finally:
        for bigwig in bigwigs.values():
            bigwig.close()


def _read_atlas_membership(
    path: Path,
) -> dict[str, list[AtlasMembershipSource]]:
    required = {
        "atlas_peak_id",
        "condition_id",
        "source_peak_id",
        "source_chrom",
        "source_start",
        "source_end",
        "source_summit",
        "atlas_center",
        "source_support_fraction",
        "source_priority_percentile",
    }
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(
                f"Atlas membership {path} is missing columns: {', '.join(missing)}"
            )
        grouped: dict[str, list[AtlasMembershipSource]] = {}
        seen: set[tuple[str, str, str, str, int, int]] = set()
        for line_number, row in enumerate(reader, start=2):
            start, end = int(row["source_start"]), int(row["source_end"])
            summit = int(row["source_summit"])
            center = int(row["atlas_center"])
            support = float(row["source_support_fraction"])
            priority = float(row["source_priority_percentile"])
            if start < 0 or end <= start or not start <= summit < end:
                raise ValueError(f"{path}:{line_number}: invalid source interval")
            if not 0 <= support <= 1 or not 0 <= priority <= 1:
                raise ValueError(f"{path}:{line_number}: invalid source weight")
            key = (
                row["atlas_peak_id"],
                row["condition_id"],
                row["source_peak_id"],
                row["source_chrom"],
                start,
                end,
            )
            if key in seen:
                raise ValueError(f"{path}:{line_number}: duplicate source peak")
            seen.add(key)
            source = AtlasMembershipSource(
                atlas_id=row["atlas_peak_id"],
                condition_id=row["condition_id"],
                peak_id=row["source_peak_id"],
                chrom=row["source_chrom"],
                start=start,
                end=end,
                summit=summit,
                atlas_center=center,
                support_fraction=support,
                priority=priority,
            )
            grouped.setdefault(source.atlas_id, []).append(source)
    return grouped


def _membership_condition_votes(
    sources: list[AtlasMembershipSource],
) -> list[AtlasMembershipSource]:
    grouped: dict[str, list[AtlasMembershipSource]] = {}
    for source in sources:
        grouped.setdefault(source.condition_id, []).append(source)
    return [
        min(
            condition_sources,
            key=lambda source: (
                -source.priority,
                -source.support_fraction,
                source.start,
                source.end,
                source.peak_id,
            ),
        )
        for condition_sources in grouped.values()
    ]


def _binned_signal(
    bigwig: Any,
    chrom: str,
    start: int,
    end: int,
    bin_size: int,
) -> list[float]:
    bin_count = math.ceil((end - start) / bin_size)
    weighted = [0.0] * bin_count
    for interval_start, interval_end, raw_value in bigwig.intervals(
        chrom, start, end
    ) or ():
        value = float(raw_value)
        if not math.isfinite(value) or value < 0:
            continue
        overlap_start = max(start, int(interval_start))
        overlap_end = min(end, int(interval_end))
        if overlap_end <= overlap_start:
            continue
        first_bin = (overlap_start - start) // bin_size
        last_bin = (overlap_end - 1 - start) // bin_size
        for index in range(first_bin, last_bin + 1):
            bin_start = start + index * bin_size
            bin_end = min(end, bin_start + bin_size)
            overlap = max(
                0,
                min(overlap_end, bin_end) - max(overlap_start, bin_start),
            )
            weighted[index] += value * overlap
    return [
        value / (min(end, start + (index + 1) * bin_size) - (start + index * bin_size))
        for index, value in enumerate(weighted)
    ]


def _smooth_signal(values: list[float], bins: int) -> list[float]:
    radius = bins // 2
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    smoothed = []
    for index in range(len(values)):
        lower = max(0, index - radius)
        upper = min(len(values), index + radius + 1)
        smoothed.append((prefix[upper] - prefix[lower]) / (upper - lower))
    return smoothed


def _signal_components(values: list[float], threshold: float) -> list[tuple[int, int]]:
    components: list[tuple[int, int]] = []
    start: int | None = None
    for index, value in enumerate(values):
        if value >= threshold and start is None:
            start = index
        elif value < threshold and start is not None:
            components.append((start, index))
            start = None
    if start is not None:
        components.append((start, len(values)))
    return components


def _write_bigwig(
    path: Path,
    chromosome_sizes: dict[str, int],
    values: dict[str, dict[int, float]],
    bin_size: int,
) -> None:
    import pyBigWig

    writer = pyBigWig.open(str(path), "w")
    if writer is None:
        raise ValueError(f"Could not create BigWig: {path}")
    try:
        writer.addHeader(list(chromosome_sizes.items()))
        for chrom in chromosome_sizes:
            chrom_values = values.get(chrom, {})
            if not chrom_values:
                continue
            starts = sorted(chrom_values)
            writer.addEntries(
                [chrom] * len(starts),
                starts,
                ends=[min(start + bin_size, chromosome_sizes[chrom]) for start in starts],
                values=[chrom_values[start] for start in starts],
            )
    finally:
        writer.close()


def build_signal_shaped_atlas(
    *,
    membership_tsv: Path,
    condition_bigwigs: dict[str, Path],
    output_bed: Path,
    aggregate_bigwig: Path,
    diagnostics_tsv: Path,
    stats_json: Path,
    window_size: int = 1000,
    bin_size: int = 10,
    smoothing_bins: int = 3,
    relative_threshold: float = 0.2,
    background_mad_multiplier: float = 3.0,
    minimum_length: int = 50,
    maximum_length: int = 400,
) -> dict[str, Any]:
    """Call signal-shaped boundaries from contributor-aware aggregate profiles."""

    import pyBigWig

    if window_size < maximum_length or window_size % 2:
        raise ValueError("window_size must be even and at least maximum_length")
    if bin_size < 1 or smoothing_bins < 1 or smoothing_bins % 2 == 0:
        raise ValueError("bin_size must be positive and smoothing_bins must be odd")
    if not 0 < relative_threshold <= 1:
        raise ValueError("relative_threshold must be in (0, 1]")
    if background_mad_multiplier < 0:
        raise ValueError("background_mad_multiplier must be non-negative")
    if minimum_length < 1 or maximum_length < minimum_length:
        raise ValueError("invalid signal-shaped peak length constraints")
    if not condition_bigwigs:
        raise ValueError("At least one condition BigWig is required")

    grouped = _read_atlas_membership(membership_tsv)
    bigwigs: dict[str, Any] = {}
    chromosome_sizes: dict[str, int] | None = None
    bed_handle, temporary_bed = _atomic_writer(output_bed)
    diagnostics_handle, temporary_diagnostics = _atomic_writer(diagnostics_tsv)
    aggregate_bigwig.parent.mkdir(parents=True, exist_ok=True)
    temporary_bigwig_handle = tempfile.NamedTemporaryFile(
        prefix=f".{aggregate_bigwig.name}.",
        dir=aggregate_bigwig.parent,
        delete=False,
    )
    temporary_bigwig = Path(temporary_bigwig_handle.name)
    temporary_bigwig_handle.close()
    temporary_bigwig.unlink()
    aggregate_values: dict[str, dict[int, float]] = {}
    widths: list[int] = []
    multimodal_count = 0
    no_signal_count = 0
    contributor_profiles = 0
    shaped_intervals: list[tuple[str, int, int]] = []
    try:
        for condition, path in condition_bigwigs.items():
            bigwig = pyBigWig.open(str(path))
            if bigwig is None:
                raise ValueError(f"Could not open BigWig: {path}")
            sizes = {str(chrom): int(size) for chrom, size in bigwig.chroms().items()}
            if chromosome_sizes is None:
                chromosome_sizes = sizes
            elif sizes != chromosome_sizes:
                raise ValueError("Condition BigWigs have different chromosome headers")
            bigwigs[condition] = bigwig
        assert chromosome_sizes is not None

        diagnostics_handle.write(
            "atlas_peak_id\tchrom\tstart\tend\twidth\tconsensus_summit\t"
            "contributing_conditions\tprofile_n\taggregate_max\tthreshold\t"
            "background_median\tbackground_mad\tstrong_mode_n\tmultimodal\n"
        )
        for atlas_id, sources in grouped.items():
            votes = _membership_condition_votes(sources)
            conditions = [vote.condition_id for vote in votes]
            missing_conditions = sorted(set(conditions) - set(bigwigs))
            if missing_conditions:
                raise ValueError(
                    f"Missing BigWigs for {atlas_id}: {', '.join(missing_conditions)}"
                )
            chroms = {vote.chrom for vote in votes}
            centers = {vote.atlas_center for vote in sources}
            if len(chroms) != 1 or len(centers) != 1:
                raise ValueError(f"Inconsistent membership coordinates for {atlas_id}")
            chrom = next(iter(chroms))
            if chrom not in chromosome_sizes:
                raise ValueError(f"Membership chromosome is absent from BigWigs: {chrom}")
            center = next(iter(centers))
            raw_start = max(0, center - window_size // 2)
            window_start = (raw_start // bin_size) * bin_size
            window_end = min(
                chromosome_sizes[chrom],
                math.ceil((center + window_size // 2) / bin_size) * bin_size,
            )
            profiles: list[list[float]] = []
            for vote in votes:
                profile = _smooth_signal(
                    _binned_signal(
                        bigwigs[vote.condition_id],
                        chrom,
                        window_start,
                        window_end,
                        bin_size,
                    ),
                    smoothing_bins,
                )
                source_bins = [
                    index
                    for index in range(len(profile))
                    if window_start + index * bin_size < vote.end
                    and min(window_end, window_start + (index + 1) * bin_size)
                    > vote.start
                ]
                maximum = max((profile[index] for index in source_bins), default=0.0)
                if maximum > 0:
                    profiles.append([min(1.0, value / maximum) for value in profile])
            contributor_profiles += len(profiles)
            if profiles:
                aggregate = [
                    statistics.median(values)
                    for values in zip(*profiles, strict=True)
                ]
            else:
                no_signal_count += 1
                aggregate = [0.0] * math.ceil(
                    (window_end - window_start) / bin_size
                )
            contributor_bins = [
                index
                for index in range(len(aggregate))
                if any(
                    window_start + index * bin_size < vote.end
                    and min(window_end, window_start + (index + 1) * bin_size)
                    > vote.start
                    for vote in votes
                )
            ]
            aggregate_maximum = max(
                (aggregate[index] for index in contributor_bins),
                default=0.0,
            )
            summit_index = max(
                contributor_bins,
                key=lambda index: (aggregate[index], -abs(
                    window_start + index * bin_size + bin_size // 2 - center
                )),
            )
            summit = min(
                chromosome_sizes[chrom] - 1,
                window_start + summit_index * bin_size + bin_size // 2,
            )
            background = []
            for index, value in enumerate(aggregate):
                bin_start = window_start + index * bin_size
                bin_end = min(window_end, bin_start + bin_size)
                if all(
                    bin_end <= vote.start or bin_start >= vote.end
                    for vote in votes
                ):
                    background.append(value)
            if not background:
                flank = min(10, len(aggregate) // 2)
                background = aggregate[:flank] + aggregate[-flank:]
            background_median = statistics.median(background) if background else 0.0
            background_mad = (
                statistics.median(
                    abs(value - background_median) for value in background
                )
                if background
                else 0.0
            )
            threshold = min(
                aggregate_maximum,
                max(
                    relative_threshold * aggregate_maximum,
                    background_median
                    + background_mad_multiplier * 1.4826 * background_mad,
                ),
            )
            components = (
                _signal_components(aggregate, threshold)
                if aggregate_maximum > 0
                else [(summit_index, summit_index + 1)]
            )
            component = next(
                (
                    (lower, upper)
                    for lower, upper in components
                    if lower <= summit_index < upper
                ),
                (summit_index, summit_index + 1),
            )
            strong_modes = (
                sum(
                    max(aggregate[lower:upper], default=0.0)
                    >= 0.5 * aggregate_maximum
                    for lower, upper in components
                )
                if aggregate_maximum > 0
                else 0
            )
            multimodal = strong_modes > 1
            multimodal_count += int(multimodal)
            start = window_start + component[0] * bin_size
            end = min(window_end, window_start + component[1] * bin_size)
            length = end - start
            if length < minimum_length:
                start, end = _fixed_window(
                    summit, minimum_length, chromosome_sizes[chrom]
                )
            elif length > maximum_length:
                start, end = _fixed_window(
                    summit, maximum_length, chromosome_sizes[chrom]
                )
            width = end - start
            widths.append(width)
            shaped_intervals.append((chrom, start, end))
            bed_handle.write(
                f"{chrom}\t{start}\t{end}\t{atlas_id}\t"
                f"{round(1000 * aggregate_maximum)}\t.\n"
            )
            diagnostics_handle.write(
                f"{atlas_id}\t{chrom}\t{start}\t{end}\t{width}\t{summit}\t"
                f"{','.join(conditions)}\t{len(profiles)}\t"
                f"{aggregate_maximum:.6g}\t{threshold:.6g}\t"
                f"{background_median:.6g}\t{background_mad:.6g}\t"
                f"{strong_modes}\t{int(multimodal)}\n"
            )
            chrom_values = aggregate_values.setdefault(chrom, {})
            for index, value in enumerate(aggregate):
                bin_start = window_start + index * bin_size
                bin_end = min(window_end, bin_start + bin_size)
                if bin_start < end and bin_end > start:
                    chrom_values[bin_start] = max(
                        chrom_values.get(bin_start, 0.0), value
                    )

        _write_bigwig(
            temporary_bigwig,
            chromosome_sizes,
            aggregate_values,
            bin_size,
        )
        bed_handle.close()
        diagnostics_handle.close()
        os.replace(temporary_bed, output_bed)
        os.replace(temporary_diagnostics, diagnostics_tsv)
        os.replace(temporary_bigwig, aggregate_bigwig)
    finally:
        for bigwig in bigwigs.values():
            bigwig.close()
        if not bed_handle.closed:
            bed_handle.close()
        if not diagnostics_handle.closed:
            diagnostics_handle.close()
        temporary_bed.unlink(missing_ok=True)
        temporary_diagnostics.unlink(missing_ok=True)
        temporary_bigwig.unlink(missing_ok=True)

    overlapping_pairs = sum(
        first_chrom == second_chrom and first_end > second_start
        for (first_chrom, _first_start, first_end), (
            second_chrom,
            second_start,
            _second_end,
        ) in zip(shaped_intervals, shaped_intervals[1:])
    )
    metrics: dict[str, Any] = {
        "status": "ok" if widths else "no_atlas_elements",
        "method": "contributor-aware median of locally normalized pooled-context signal",
        "atlas_peaks": len(widths),
        "condition_n": len(condition_bigwigs),
        "contributor_profiles": contributor_profiles,
        "no_signal_peaks": no_signal_count,
        "multimodal_peaks": multimodal_count,
        "overlapping_neighbor_pairs": overlapping_pairs,
        "window_size": window_size,
        "bin_size": bin_size,
        "smoothing_bins": smoothing_bins,
        "relative_threshold": relative_threshold,
        "background_mad_multiplier": background_mad_multiplier,
        "minimum_length": minimum_length,
        "maximum_length": maximum_length,
        "width_min": min(widths, default=0),
        "width_median": statistics.median(widths) if widths else 0,
        "width_mean": statistics.mean(widths) if widths else 0,
        "width_max": max(widths, default=0),
        "aggregation": "equal-weight median across contributing conditions",
        "per_condition_normalization": (
            "smoothed local profile divided by its maximum inside the assigned source DHS, capped at 1"
        ),
        "summit_constraint": "aggregate maximum must overlap a contributing source DHS",
        "boundary_rule": (
            "contiguous component around aggregate maximum above the configured relative and robust-background thresholds"
        ),
    }
    _write_json_atomic(stats_json, metrics)
    return metrics


def _read_chromosome_sizes(path: Path) -> dict[str, int]:
    chromosome_sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected chromosome and size")
            chrom, size = fields[0], int(fields[1])
            if chrom in chromosome_sizes or size < 1:
                raise ValueError(f"{path}:{line_number}: invalid chromosome size")
            chromosome_sizes[chrom] = size
    if not chromosome_sizes:
        raise ValueError(f"Chromosome sizes file is empty: {path}")
    return chromosome_sizes


def _read_merged_bed_intervals(
    path: Path, chromosome_sizes: dict[str, int]
) -> dict[str, list[tuple[int, int]]]:
    grouped: dict[str, list[tuple[int, int]]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 columns")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in chromosome_sizes:
                raise ValueError(f"{path}:{line_number}: unknown chromosome {chrom}")
            if start < 0 or end <= start or end > chromosome_sizes[chrom]:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            grouped.setdefault(chrom, []).append((start, end))

    merged_by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in grouped.items():
        merged: list[tuple[int, int]] = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        merged_by_chrom[chrom] = merged
    return merged_by_chrom


def _build_dhs_support_segments(
    condition_dhs: dict[str, Path], chromosome_sizes: dict[str, int]
) -> tuple[dict[str, list[SupportSegment]], int]:
    events: dict[str, dict[int, int]] = {}
    merged_interval_count = 0
    for path in condition_dhs.values():
        intervals_by_chrom = _read_merged_bed_intervals(path, chromosome_sizes)
        for chrom, intervals in intervals_by_chrom.items():
            chrom_events = events.setdefault(chrom, {})
            for start, end in intervals:
                chrom_events[start] = chrom_events.get(start, 0) + 1
                chrom_events[end] = chrom_events.get(end, 0) - 1
                merged_interval_count += 1

    segments_by_chrom: dict[str, list[SupportSegment]] = {}
    for chrom in chromosome_sizes:
        chrom_events = events.get(chrom)
        if not chrom_events:
            continue
        support = 0
        previous: int | None = None
        segments: list[SupportSegment] = []
        for position, delta in sorted(chrom_events.items()):
            if previous is not None and previous < position and support > 0:
                if segments and segments[-1].end == previous and segments[-1].support == support:
                    prior = segments[-1]
                    segments[-1] = SupportSegment(chrom, prior.start, position, support)
                else:
                    segments.append(SupportSegment(chrom, previous, position, support))
            support += delta
            if support < 0:
                raise ValueError(f"Negative DHS support on {chrom}:{position}")
            previous = position
        if support != 0:
            raise ValueError(f"Unbalanced DHS intervals on {chrom}")
        segments_by_chrom[chrom] = segments
    return segments_by_chrom, merged_interval_count


def _write_support_bigwig(
    path: Path,
    chromosome_sizes: dict[str, int],
    segments_by_chrom: dict[str, list[SupportSegment]],
    condition_n: int,
) -> None:
    import pyBigWig

    writer = pyBigWig.open(str(path), "w")
    if writer is None:
        raise ValueError(f"Could not create BigWig: {path}")
    try:
        writer.addHeader(list(chromosome_sizes.items()))
        for chrom in chromosome_sizes:
            segments = segments_by_chrom.get(chrom, [])
            if not segments:
                continue
            writer.addEntries(
                [chrom] * len(segments),
                [segment.start for segment in segments],
                ends=[segment.end for segment in segments],
                values=[segment.support / condition_n for segment in segments],
            )
    finally:
        writer.close()


def _read_atlas_anchors(
    path: Path, chromosome_sizes: dict[str, int]
) -> list[tuple[str, int, int, str]]:
    chromosome_order = {chrom: index for index, chrom in enumerate(chromosome_sizes)}
    anchors: list[tuple[str, int, int, str]] = []
    previous_key: tuple[int, int, int] | None = None
    seen_ids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least 4 columns")
            chrom, start, end, atlas_id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if chrom not in chromosome_sizes:
                raise ValueError(f"{path}:{line_number}: unknown chromosome {chrom}")
            if start < 0 or end <= start or end > chromosome_sizes[chrom]:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            if atlas_id in seen_ids:
                raise ValueError(f"{path}:{line_number}: duplicate atlas ID {atlas_id}")
            key = (chromosome_order[chrom], start, end)
            if previous_key is not None and key < previous_key:
                raise ValueError(f"{path}:{line_number}: anchors are not coordinate sorted")
            anchors.append((chrom, start, end, atlas_id))
            seen_ids.add(atlas_id)
            previous_key = key
    return anchors


def build_dhs_support_fwhm(
    *,
    anchors_bed: Path,
    condition_dhs: dict[str, Path],
    chromosome_sizes_path: Path,
    support_bigwig: Path,
    output_bed: Path,
    diagnostics_tsv: Path,
    stats_json: Path,
) -> dict[str, Any]:
    """Build a tissue-balanced DHS track and anchor-centered half-maximum widths."""

    if not condition_dhs:
        raise ValueError("At least one condition DHS BED is required")
    chromosome_sizes = _read_chromosome_sizes(chromosome_sizes_path)
    anchors = _read_atlas_anchors(anchors_bed, chromosome_sizes)
    segments_by_chrom, merged_interval_count = _build_dhs_support_segments(
        condition_dhs, chromosome_sizes
    )
    starts_by_chrom = {
        chrom: [segment.start for segment in segments]
        for chrom, segments in segments_by_chrom.items()
    }

    boundaries: list[FwhmBoundary] = []
    for chrom, anchor_start, anchor_end, atlas_id in anchors:
        segments = segments_by_chrom.get(chrom, [])
        starts = starts_by_chrom.get(chrom, [])
        index = max(0, bisect.bisect_right(starts, anchor_start) - 1)
        overlapping_indices: list[int] = []
        while index < len(segments) and segments[index].start < anchor_end:
            if segments[index].end > anchor_start:
                overlapping_indices.append(index)
            index += 1
        anchor_center = (anchor_start + anchor_end) // 2
        if not overlapping_indices:
            boundaries.append(
                FwhmBoundary(
                    atlas_id,
                    chrom,
                    anchor_start,
                    anchor_end,
                    anchor_start,
                    anchor_end,
                    anchor_center,
                    0,
                    0,
                    0,
                    0,
                    False,
                    "no_dhs_support",
                )
            )
            continue

        def maximum_key(segment_index: int) -> tuple[int, int, int]:
            segment = segments[segment_index]
            overlap_start = max(anchor_start, segment.start)
            overlap_end = min(anchor_end, segment.end)
            nearest = min(max(anchor_center, overlap_start), overlap_end - 1)
            return (-segment.support, abs(nearest - anchor_center), segment.start)

        selected_index = min(overlapping_indices, key=maximum_key)
        selected = segments[selected_index]
        maximum_support = selected.support
        half_maximum_support = math.ceil(maximum_support / 2)
        overlap_start = max(anchor_start, selected.start)
        overlap_end = min(anchor_end, selected.end)
        summit = min(max(anchor_center, overlap_start), overlap_end - 1)

        left = selected_index
        while (
            left > 0
            and segments[left - 1].end == segments[left].start
            and segments[left - 1].support >= half_maximum_support
        ):
            left -= 1
        right = selected_index
        while (
            right + 1 < len(segments)
            and segments[right].end == segments[right + 1].start
            and segments[right + 1].support >= half_maximum_support
        ):
            right += 1
        component = segments[left : right + 1]
        component_maximum = max(segment.support for segment in component)
        maximum_plateaus = sum(
            segment.support == maximum_support for segment in component
        )
        boundaries.append(
            FwhmBoundary(
                atlas_id,
                chrom,
                anchor_start,
                anchor_end,
                component[0].start,
                component[-1].end,
                summit,
                maximum_support,
                half_maximum_support,
                component_maximum,
                maximum_plateaus,
                component_maximum > maximum_support,
                "ok",
            )
        )

    support_bigwig.parent.mkdir(parents=True, exist_ok=True)
    temporary_bigwig_handle = tempfile.NamedTemporaryFile(
        prefix=f".{support_bigwig.name}.",
        dir=support_bigwig.parent,
        delete=False,
    )
    temporary_bigwig = Path(temporary_bigwig_handle.name)
    temporary_bigwig_handle.close()
    temporary_bigwig.unlink()
    bed_handle, temporary_bed = _atomic_writer(output_bed)
    diagnostics_handle, temporary_diagnostics = _atomic_writer(diagnostics_tsv)
    neighbor_contact_count = 0
    try:
        _write_support_bigwig(
            temporary_bigwig,
            chromosome_sizes,
            segments_by_chrom,
            len(condition_dhs),
        )
        diagnostics_handle.write(
            "atlas_peak_id\tchrom\tanchor_start\tanchor_end\tanchor_center\t"
            "start\tend\twidth\tsummit\tmaximum_support_n\tcondition_n\t"
            "maximum_support_fraction\thalf_maximum_support_n\t"
            "component_maximum_support_n\tmaximum_plateau_n\t"
            "bridged_higher_peak\tneighbor_contact\tstatus\n"
        )
        for index, boundary in enumerate(boundaries):
            previous = boundaries[index - 1] if index else None
            following = boundaries[index + 1] if index + 1 < len(boundaries) else None
            neighbor_contact = (
                previous is not None
                and previous.chrom == boundary.chrom
                and boundary.start < previous.anchor_end
            ) or (
                following is not None
                and following.chrom == boundary.chrom
                and boundary.end > following.anchor_start
            )
            neighbor_contact_count += int(neighbor_contact)
            score = round(
                1000 * boundary.maximum_support / len(condition_dhs)
            )
            bed_handle.write(
                f"{boundary.chrom}\t{boundary.start}\t{boundary.end}\t"
                f"{boundary.atlas_id}\t{score}\t.\n"
            )
            diagnostics_handle.write(
                f"{boundary.atlas_id}\t{boundary.chrom}\t{boundary.anchor_start}\t"
                f"{boundary.anchor_end}\t"
                f"{(boundary.anchor_start + boundary.anchor_end) // 2}\t"
                f"{boundary.start}\t{boundary.end}\t"
                f"{boundary.end - boundary.start}\t{boundary.summit}\t"
                f"{boundary.maximum_support}\t{len(condition_dhs)}\t"
                f"{boundary.maximum_support / len(condition_dhs):.6g}\t"
                f"{boundary.half_maximum_support}\t"
                f"{boundary.component_maximum_support}\t"
                f"{boundary.maximum_plateaus}\t"
                f"{int(boundary.bridged_higher_peak)}\t{int(neighbor_contact)}\t"
                f"{boundary.status}\n"
            )
        _finish_atomic(bed_handle, temporary_bed, output_bed)
        _finish_atomic(
            diagnostics_handle,
            temporary_diagnostics,
            diagnostics_tsv,
        )
        os.replace(temporary_bigwig, support_bigwig)
    finally:
        if not bed_handle.closed:
            bed_handle.close()
        if not diagnostics_handle.closed:
            diagnostics_handle.close()
        temporary_bed.unlink(missing_ok=True)
        temporary_diagnostics.unlink(missing_ok=True)
        temporary_bigwig.unlink(missing_ok=True)

    widths = [boundary.end - boundary.start for boundary in boundaries]
    maximum_support_values = [boundary.maximum_support for boundary in boundaries]
    overlapping_pairs = sum(
        first.chrom == second.chrom and first.end > second.start
        for first, second in zip(boundaries, boundaries[1:])
    )
    metrics: dict[str, Any] = {
        "status": "ok" if boundaries else "no_atlas_elements",
        "method": "anchor-centered half-maximum width of tissue-balanced DHS support",
        "atlas_peaks": len(boundaries),
        "condition_n": len(condition_dhs),
        "conditions": list(condition_dhs),
        "merged_condition_dhs_intervals": merged_interval_count,
        "track_value": "fraction of conditions covered by at least one consensus DHS",
        "local_maximum_rule": "maximum support inside the fixed anchor, nearest its center",
        "half_maximum_rule": "support >= ceil(local maximum support / 2)",
        "component_rule": "connected half-maximum component containing the selected local maximum",
        "no_support_peaks": sum(
            boundary.status == "no_dhs_support" for boundary in boundaries
        ),
        "bridged_higher_peak_count": sum(
            boundary.bridged_higher_peak for boundary in boundaries
        ),
        "neighbor_contact_peaks": neighbor_contact_count,
        "overlapping_neighbor_pairs": overlapping_pairs,
        "width_min": min(widths, default=0),
        "width_median": statistics.median(widths) if widths else 0,
        "width_mean": statistics.mean(widths) if widths else 0,
        "width_max": max(widths, default=0),
        "maximum_support_min": min(maximum_support_values, default=0),
        "maximum_support_median": (
            statistics.median(maximum_support_values)
            if maximum_support_values
            else 0
        ),
        "maximum_support_max": max(maximum_support_values, default=0),
    }
    _write_json_atomic(stats_json, metrics)
    return metrics


def _is_local_support_maximum(
    segments: list[SupportSegment], index: int
) -> bool:
    segment = segments[index]
    left_support = (
        segments[index - 1].support
        if index > 0 and segments[index - 1].end == segment.start
        else 0
    )
    right_support = (
        segments[index + 1].support
        if index + 1 < len(segments) and segment.end == segments[index + 1].start
        else 0
    )
    return segment.support >= max(left_support, right_support) and (
        segment.support > left_support or segment.support > right_support
    )


def _support_prominence_base(
    segments: list[SupportSegment], index: int, step: int
) -> int:
    peak_support = segments[index].support
    minimum_support = peak_support
    current = index
    while True:
        following = current + step
        if following < 0 or following >= len(segments):
            return 0
        if step < 0:
            contiguous = segments[following].end == segments[current].start
        else:
            contiguous = segments[current].end == segments[following].start
        if not contiguous:
            return 0
        if segments[following].support > peak_support:
            return minimum_support
        minimum_support = min(minimum_support, segments[following].support)
        current = following


def build_dhs_center_mode_half_prominence(
    *,
    anchors_bed: Path,
    condition_dhs: dict[str, Path],
    chromosome_sizes_path: Path,
    output_bed: Path,
    diagnostics_tsv: Path,
    stats_json: Path,
) -> dict[str, Any]:
    """Measure the center-associated DHS-support mode at half prominence."""

    if not condition_dhs:
        raise ValueError("At least one condition DHS BED is required")
    chromosome_sizes = _read_chromosome_sizes(chromosome_sizes_path)
    anchors = _read_atlas_anchors(anchors_bed, chromosome_sizes)
    segments_by_chrom, merged_interval_count = _build_dhs_support_segments(
        condition_dhs, chromosome_sizes
    )
    starts_by_chrom = {
        chrom: [segment.start for segment in segments]
        for chrom, segments in segments_by_chrom.items()
    }

    boundaries: list[LocalModeBoundary] = []
    for chrom, anchor_start, anchor_end, atlas_id in anchors:
        segments = segments_by_chrom.get(chrom, [])
        starts = starts_by_chrom.get(chrom, [])
        index = max(0, bisect.bisect_right(starts, anchor_start) - 1)
        overlapping_indices: list[int] = []
        while index < len(segments) and segments[index].start < anchor_end:
            if segments[index].end > anchor_start:
                overlapping_indices.append(index)
            index += 1
        anchor_center = (anchor_start + anchor_end) // 2
        if not overlapping_indices:
            boundaries.append(
                LocalModeBoundary(
                    atlas_id,
                    chrom,
                    anchor_start,
                    anchor_end,
                    anchor_start,
                    anchor_end,
                    anchor_center,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    True,
                    False,
                    "no_dhs_support",
                )
            )
            continue

        local_modes = [
            segment_index
            for segment_index in overlapping_indices
            if _is_local_support_maximum(segments, segment_index)
        ]
        fallback = not local_modes
        candidates = local_modes or overlapping_indices

        def center_key(segment_index: int) -> tuple[int, int, int, int]:
            segment = segments[segment_index]
            overlap_start = max(anchor_start, segment.start)
            overlap_end = min(anchor_end, segment.end)
            nearest = min(max(anchor_center, overlap_start), overlap_end - 1)
            return (
                abs(nearest - anchor_center),
                segment.end - segment.start,
                -segment.support,
                segment.start,
            )

        selected_index = min(candidates, key=center_key)
        selected = segments[selected_index]
        overlap_start = max(anchor_start, selected.start)
        overlap_end = min(anchor_end, selected.end)
        summit = min(max(anchor_center, overlap_start), overlap_end - 1)
        left_base = _support_prominence_base(segments, selected_index, -1)
        right_base = _support_prominence_base(segments, selected_index, 1)
        contour_base = max(left_base, right_base)
        prominence = max(0, selected.support - contour_base)
        half_prominence_support = math.ceil(
            contour_base + prominence / 2
        )

        left = selected_index
        while (
            left > 0
            and segments[left - 1].end == segments[left].start
            and segments[left - 1].support >= half_prominence_support
        ):
            left -= 1
        right = selected_index
        while (
            right + 1 < len(segments)
            and segments[right].end == segments[right + 1].start
            and segments[right + 1].support >= half_prominence_support
        ):
            right += 1
        component = segments[left : right + 1]
        boundaries.append(
            LocalModeBoundary(
                atlas_id,
                chrom,
                anchor_start,
                anchor_end,
                component[0].start,
                component[-1].end,
                summit,
                selected.support,
                left_base,
                right_base,
                prominence,
                half_prominence_support,
                len(local_modes),
                fallback,
                max(segment.support for segment in component) > selected.support,
                "ok" if prominence > 0 else "zero_prominence",
            )
        )

    bed_handle, temporary_bed = _atomic_writer(output_bed)
    diagnostics_handle, temporary_diagnostics = _atomic_writer(diagnostics_tsv)
    neighbor_contact_count = 0
    try:
        diagnostics_handle.write(
            "atlas_peak_id\tchrom\tanchor_start\tanchor_end\tanchor_center\t"
            "start\tend\twidth\tsummit\tpeak_support_n\tcondition_n\t"
            "peak_support_fraction\tleft_base_support_n\tright_base_support_n\t"
            "prominence_n\thalf_prominence_support_n\tlocal_mode_n\t"
            "fallback_nonmaximum\tbridged_higher_peak\tneighbor_contact\tstatus\n"
        )
        for index, boundary in enumerate(boundaries):
            previous = boundaries[index - 1] if index else None
            following = boundaries[index + 1] if index + 1 < len(boundaries) else None
            neighbor_contact = (
                previous is not None
                and previous.chrom == boundary.chrom
                and boundary.start < previous.anchor_end
            ) or (
                following is not None
                and following.chrom == boundary.chrom
                and boundary.end > following.anchor_start
            )
            neighbor_contact_count += int(neighbor_contact)
            score = round(1000 * boundary.peak_support / len(condition_dhs))
            bed_handle.write(
                f"{boundary.chrom}\t{boundary.start}\t{boundary.end}\t"
                f"{boundary.atlas_id}\t{score}\t.\n"
            )
            diagnostics_handle.write(
                f"{boundary.atlas_id}\t{boundary.chrom}\t{boundary.anchor_start}\t"
                f"{boundary.anchor_end}\t"
                f"{(boundary.anchor_start + boundary.anchor_end) // 2}\t"
                f"{boundary.start}\t{boundary.end}\t"
                f"{boundary.end - boundary.start}\t{boundary.summit}\t"
                f"{boundary.peak_support}\t{len(condition_dhs)}\t"
                f"{boundary.peak_support / len(condition_dhs):.6g}\t"
                f"{boundary.left_base_support}\t{boundary.right_base_support}\t"
                f"{boundary.prominence}\t{boundary.half_prominence_support}\t"
                f"{boundary.local_mode_n}\t{int(boundary.fallback_nonmaximum)}\t"
                f"{int(boundary.bridged_higher_peak)}\t{int(neighbor_contact)}\t"
                f"{boundary.status}\n"
            )
        _finish_atomic(bed_handle, temporary_bed, output_bed)
        _finish_atomic(
            diagnostics_handle,
            temporary_diagnostics,
            diagnostics_tsv,
        )
    finally:
        if not bed_handle.closed:
            bed_handle.close()
        if not diagnostics_handle.closed:
            diagnostics_handle.close()
        temporary_bed.unlink(missing_ok=True)
        temporary_diagnostics.unlink(missing_ok=True)

    widths = [boundary.end - boundary.start for boundary in boundaries]
    overlapping_pairs = sum(
        first.chrom == second.chrom and first.end > second.start
        for first, second in zip(boundaries, boundaries[1:])
    )
    metrics: dict[str, Any] = {
        "status": "ok" if boundaries else "no_atlas_elements",
        "method": "half-prominence width of the DHS-support local mode nearest the fixed-anchor center",
        "atlas_peaks": len(boundaries),
        "condition_n": len(condition_dhs),
        "conditions": list(condition_dhs),
        "merged_condition_dhs_intervals": merged_interval_count,
        "track_value": "fraction of conditions covered by at least one consensus DHS",
        "mode_selection_rule": "local support maximum nearest the fixed-anchor center; narrower plateau breaks distance ties",
        "prominence_rule": "peak support minus the higher of the nearest left and right valley bases before a higher peak or zero-support gap",
        "width_rule": "connected component at support >= ceil(base + prominence / 2)",
        "no_support_peaks": sum(
            boundary.status == "no_dhs_support" for boundary in boundaries
        ),
        "zero_prominence_peaks": sum(
            boundary.status == "zero_prominence" for boundary in boundaries
        ),
        "fallback_nonmaximum_peaks": sum(
            boundary.fallback_nonmaximum for boundary in boundaries
        ),
        "multiple_local_mode_peaks": sum(
            boundary.local_mode_n > 1 for boundary in boundaries
        ),
        "bridged_higher_peak_count": sum(
            boundary.bridged_higher_peak for boundary in boundaries
        ),
        "neighbor_contact_peaks": neighbor_contact_count,
        "overlapping_neighbor_pairs": overlapping_pairs,
        "width_min": min(widths, default=0),
        "width_median": statistics.median(widths) if widths else 0,
        "width_mean": statistics.mean(widths) if widths else 0,
        "width_max": max(widths, default=0),
    }
    _write_json_atomic(stats_json, metrics)
    return metrics
