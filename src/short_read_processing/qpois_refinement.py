"""Refine lenient MACS3 candidates with the corresponding qpois signal."""

from __future__ import annotations

import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


ALGORITHM = "event_driven_progressive_qpois_v1"


@dataclass(frozen=True)
class Candidate:
    chrom: str
    start: int
    end: int


@dataclass(frozen=True)
class Segment:
    start: int
    end: int
    score: float
    last_exponent: int


@dataclass(frozen=True)
class RefinedPeak:
    chrom: str
    start: int
    end: int
    maximum_qscore: float
    selection_exponent: int


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


def read_candidates(path: Path) -> tuple[dict[str, list[Candidate]], dict[str, int]]:
    candidates: dict[str, list[Candidate]] = {}
    chromosome_order: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least three columns")
            candidate = Candidate(fields[0], int(fields[1]), int(fields[2]))
            if candidate.start < 0 or candidate.end <= candidate.start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            chromosome_order.setdefault(candidate.chrom, len(chromosome_order))
            values = candidates.setdefault(candidate.chrom, [])
            if values and candidate.start < values[-1].end:
                raise ValueError(
                    f"{path}:{line_number}: candidates overlap or are not coordinate sorted"
                )
            values.append(candidate)
    return candidates, chromosome_order


def split_components(segments: Iterable[Segment], gap: int) -> list[list[Segment]]:
    components: list[list[Segment]] = []
    current: list[Segment] = []
    current_end = -1
    for segment in segments:
        if current and segment.start > current_end + gap:
            components.append(current)
            current = []
        current.append(segment)
        current_end = max(current_end, segment.end) if len(current) > 1 else segment.end
    if current:
        components.append(current)
    return components


def _overlaps(interval: RefinedPeak, others: list[RefinedPeak]) -> bool:
    return any(other.start < interval.end and interval.start < other.end for other in others)


def refine_candidate(
    candidate: Candidate,
    segments: list[Segment],
    *,
    minimum_exponent: int,
    maximum_exponent: int,
    minimum_length: int,
    maximum_length: int,
    merge_gap: int,
) -> tuple[list[RefinedPeak], list[RefinedPeak]]:
    selected: list[RefinedPeak] = []

    def visit(component: list[Segment], exponent: int) -> None:
        start = component[0].start
        end = max(segment.end for segment in component)
        width = end - start
        if width < minimum_length:
            return
        if width <= maximum_length:
            selected.append(
                RefinedPeak(
                    candidate.chrom,
                    start,
                    end,
                    max(segment.score for segment in component),
                    exponent,
                )
            )
            return

        next_exponent = min(segment.last_exponent for segment in component) + 1
        if next_exponent > maximum_exponent:
            return
        surviving = [
            segment for segment in component if segment.last_exponent >= next_exponent
        ]
        for child in split_components(surviving, merge_gap):
            visit(child, next_exponent)

    for component in split_components(segments, merge_gap):
        visit(component, minimum_exponent)

    excluded: list[RefinedPeak] = []
    final_segments = [segment for segment in segments if segment.score > maximum_exponent]
    for component in split_components(final_segments, merge_gap):
        interval = RefinedPeak(
            candidate.chrom,
            component[0].start,
            max(segment.end for segment in component),
            max(segment.score for segment in component),
            maximum_exponent,
        )
        if interval.end - interval.start >= minimum_length and not _overlaps(
            interval, selected
        ):
            excluded.append(interval)
    return selected, excluded


def _last_active_exponent(score: float, maximum_exponent: int) -> int:
    if not math.isfinite(score):
        return maximum_exponent
    return min(maximum_exponent, math.ceil(score) - 1)


def refine_qpois(
    *,
    qpois_bedgraph: Path,
    candidate_peaks: Path,
    minimum_exponent: int = 2,
    maximum_exponent: int = 325,
    minimum_length: int = 50,
    maximum_length: int = 400,
    merge_gap: int = 1,
) -> tuple[list[RefinedPeak], list[RefinedPeak], dict[str, int]]:
    if minimum_exponent < 0 or maximum_exponent < minimum_exponent:
        raise ValueError("qpois exponent range is invalid")
    if minimum_length < 1 or maximum_length < minimum_length or merge_gap < 0:
        raise ValueError("qpois refinement geometry is invalid")

    candidates, chromosome_order = read_candidates(candidate_peaks)
    selected: list[RefinedPeak] = []
    excluded: list[RefinedPeak] = []
    contained_segments = 0
    current_chrom = ""
    current_candidates: list[Candidate] = []
    candidate_index = 0
    candidate_segments: list[Segment] = []

    def flush_candidate() -> None:
        nonlocal candidate_segments
        if candidate_segments and candidate_index < len(current_candidates):
            refined, rejected = refine_candidate(
                current_candidates[candidate_index],
                candidate_segments,
                minimum_exponent=minimum_exponent,
                maximum_exponent=maximum_exponent,
                minimum_length=minimum_length,
                maximum_length=maximum_length,
                merge_gap=merge_gap,
            )
            selected.extend(refined)
            excluded.extend(rejected)
        candidate_segments = []

    with qpois_bedgraph.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{qpois_bedgraph}:{line_number}: expected four columns"
                )
            chrom, start, end, score = (
                fields[0],
                int(fields[1]),
                int(fields[2]),
                float(fields[3]),
            )
            if end <= start or not math.isfinite(score):
                raise ValueError(f"{qpois_bedgraph}:{line_number}: invalid segment")
            if chrom != current_chrom:
                flush_candidate()
                current_chrom = chrom
                current_candidates = candidates.get(chrom, [])
                candidate_index = 0
            while (
                candidate_index < len(current_candidates)
                and current_candidates[candidate_index].end <= start
            ):
                flush_candidate()
                candidate_index += 1
            if candidate_index >= len(current_candidates):
                continue
            candidate = current_candidates[candidate_index]
            if score > minimum_exponent and candidate.start <= start and end <= candidate.end:
                candidate_segments.append(
                    Segment(start, end, score, _last_active_exponent(score, maximum_exponent))
                )
                contained_segments += 1
    flush_candidate()

    sort_key = lambda peak: (
        chromosome_order.get(peak.chrom, len(chromosome_order)),
        peak.start,
        peak.end,
    )
    selected.sort(key=sort_key)
    excluded.sort(key=sort_key)
    return selected, excluded, {
        "candidate_peaks": sum(map(len, candidates.values())),
        "contained_qpois_segments": contained_segments,
    }


def _write_peaks(path: Path, peaks: list[RefinedPeak], prefix: str) -> None:
    handle, temporary = _atomic_writer(path)
    try:
        for index, peak in enumerate(peaks, start=1):
            bed_score = min(1000, max(0, int(round(peak.maximum_qscore))))
            handle.write(
                f"{peak.chrom}\t{peak.start}\t{peak.end}\t{prefix}_{index:06d}\t"
                f"{bed_score}\t.\t{peak.maximum_qscore:.10g}\t"
                f"{peak.selection_exponent}\n"
            )
        handle.close()
        os.replace(temporary, path)
    finally:
        if not handle.closed:
            handle.close()
        temporary.unlink(missing_ok=True)


def run_refinement(
    *,
    qpois_bedgraph: Path,
    candidate_peaks: Path,
    output_bed: Path,
    excluded_bed: Path,
    stats_json: Path,
    name_prefix: str,
    minimum_exponent: int = 2,
    maximum_exponent: int = 325,
    minimum_length: int = 50,
    maximum_length: int = 400,
    merge_gap: int = 1,
) -> dict[str, object]:
    selected, excluded, counts = refine_qpois(
        qpois_bedgraph=qpois_bedgraph,
        candidate_peaks=candidate_peaks,
        minimum_exponent=minimum_exponent,
        maximum_exponent=maximum_exponent,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        merge_gap=merge_gap,
    )
    _write_peaks(output_bed, selected, f"{name_prefix}_QPOIS")
    _write_peaks(excluded_bed, excluded, f"{name_prefix}_EXCLUDED")
    metrics: dict[str, object] = {
        "algorithm": ALGORITHM,
        "status": "ok" if selected else "no_refined_peaks",
        **counts,
        "refined_peaks": len(selected),
        "excluded_peaks": len(excluded),
        "parameters": {
            "minimum_exponent": minimum_exponent,
            "maximum_exponent": maximum_exponent,
            "minimum_length": minimum_length,
            "maximum_length": maximum_length,
            "merge_gap": merge_gap,
        },
    }
    handle, temporary = _atomic_writer(stats_json)
    try:
        json.dump(metrics, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.close()
        os.replace(temporary, stats_json)
    finally:
        if not handle.closed:
            handle.close()
        temporary.unlink(missing_ok=True)
    return metrics
