"""Build a cross-context master DHS registry from final ATAC peak sets."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
import json
import math
import os
from pathlib import Path
import statistics
import tempfile
from typing import Any, Iterable, TextIO


@dataclass(frozen=True)
class SourcePeak:
    context: str
    chrom: str
    start: int
    end: int
    name: str
    score: int = 0
    strand: str = "."
    support_n: int = 0
    replicate_n: int = 0
    support_fraction: float = 0.0
    supporting_libraries: str = ""
    peak_method: str = ""
    summit: int | None = None
    summit_signal: float = 0.0
    summit_fallback: bool = False
    input_start: int | None = None
    input_end: int | None = None
    coordinate_clipped: bool = False

    @property
    def width(self) -> int:
        return self.end - self.start


@dataclass
class MasterDHS:
    chrom: str
    start: int
    end: int
    summit: int
    members: list[SourcePeak]
    master_id: str = ""


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected two columns")
            chrom, size_text = fields
            size = int(size_text)
            if not chrom or size < 1 or chrom in sizes:
                raise ValueError(f"{path}:{line_number}: invalid chromosome entry")
            sizes[chrom] = size
    if not sizes:
        raise ValueError(f"Chromosome sizes file is empty: {path}")
    return sizes


def read_context_peaks(
    path: Path,
    *,
    context: str,
    chrom_sizes: dict[str, int],
) -> list[SourcePeak]:
    peaks: list[SourcePeak] = []
    name_counts: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least three columns")
            chrom, input_start, input_end = fields[0], int(fields[1]), int(fields[2])
            if chrom not in chrom_sizes:
                raise ValueError(f"{path}:{line_number}: unknown chromosome {chrom!r}")
            if input_start < 0 or input_end <= input_start or input_start >= chrom_sizes[chrom]:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            start = input_start
            end = min(input_end, chrom_sizes[chrom])
            raw_name = fields[3] if len(fields) > 3 and fields[3] else f"peak_{len(peaks) + 1}"
            occurrence = name_counts.get(raw_name, 0) + 1
            name_counts[raw_name] = occurrence
            name = raw_name if occurrence == 1 else f"{raw_name}__duplicate_{occurrence}"
            try:
                score = int(float(fields[4])) if len(fields) > 4 else 0
            except ValueError:
                score = 0
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-", "."} else "."
            if len(fields) >= 12 and fields[6] != context:
                raise ValueError(
                    f"{path}:{line_number}: BED context {fields[6]!r} does not match {context!r}"
                )
            peaks.append(
                SourcePeak(
                    context=context,
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=min(1000, max(0, score)),
                    strand=strand,
                    support_n=int(fields[7]) if len(fields) >= 12 else 0,
                    replicate_n=int(fields[8]) if len(fields) >= 12 else 0,
                    support_fraction=float(fields[9]) if len(fields) >= 12 else 0.0,
                    supporting_libraries=fields[10] if len(fields) >= 12 else "",
                    peak_method=fields[11] if len(fields) >= 12 else "",
                    input_start=input_start,
                    input_end=input_end,
                    coordinate_clipped=input_end != end,
                )
            )
    return peaks


def _summit_from_values(peak: SourcePeak, values: Iterable[float]) -> SourcePeak:
    finite = [float(value) if value is not None else math.nan for value in values]
    valid = [value for value in finite if math.isfinite(value)]
    if not valid:
        return replace(
            peak,
            summit=(peak.start + peak.end - 1) // 2,
            summit_signal=0.0,
            summit_fallback=True,
        )

    maximum = max(valid)
    maximal_indices = [
        index for index, value in enumerate(finite) if math.isfinite(value) and value == maximum
    ]
    runs: list[tuple[int, int]] = []
    run_start = maximal_indices[0]
    previous = run_start
    for index in maximal_indices[1:]:
        if index != previous + 1:
            runs.append((run_start, previous))
            run_start = index
        previous = index
    runs.append((run_start, previous))

    peak_midpoint = (peak.width - 1) / 2
    selected = min(
        runs,
        key=lambda run: (
            abs(((run[0] + run[1]) / 2) - peak_midpoint),
            run[0],
        ),
    )
    offset = (selected[0] + selected[1]) // 2
    return replace(
        peak,
        summit=peak.start + offset,
        summit_signal=maximum,
        summit_fallback=False,
    )


def load_context_peaks(
    context: str,
    peaks_path: Path,
    signal_path: Path,
    chrom_sizes: dict[str, int],
) -> list[SourcePeak]:
    """Read one context and find each peak's maximum in its pooled signal track."""

    try:
        import pyBigWig  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - exercised in the rule environment
        raise RuntimeError("pyBigWig is required to build master DHS summits") from exc

    peaks = read_context_peaks(peaks_path, context=context, chrom_sizes=chrom_sizes)
    bigwig = pyBigWig.open(str(signal_path))
    if bigwig is None:
        raise ValueError(f"Could not open BigWig: {signal_path}")
    try:
        bigwig_chroms = bigwig.chroms()
        missing = sorted({peak.chrom for peak in peaks} - set(bigwig_chroms))
        if missing:
            raise ValueError(
                f"BigWig {signal_path} is missing peak chromosomes: " + ", ".join(missing)
            )
        return [
            _summit_from_values(
                peak,
                bigwig.values(peak.chrom, peak.start, peak.end),
            )
            for peak in peaks
        ]
    finally:
        bigwig.close()


def _load_context_job(
    arguments: tuple[str, Path, Path, dict[str, int]],
) -> list[SourcePeak]:
    return load_context_peaks(*arguments)


def _reciprocally_contains_summits(left: SourcePeak, right: SourcePeak) -> bool:
    return (
        left.chrom == right.chrom
        and left.start <= int(right.summit) < left.end
        and right.start <= int(left.summit) < right.end
    )


def _representative_summit(members: list[SourcePeak]) -> int:
    summits = [peak.summit for peak in members]
    if any(summit is None for summit in summits):
        raise ValueError("Every source peak must have a summit before clustering")
    median = statistics.median(int(summit) for summit in summits if summit is not None)
    selected = min(
        members,
        key=lambda peak: (
            abs(int(peak.summit) - median),
            peak.width,
            peak.context,
            peak.name,
        ),
    )
    return int(selected.summit)


def _merge_close_context_shifted_clusters(
    clusters: list[list[SourcePeak]],
    *,
    summit_max_distance: int,
    minimum_summit_separation: int,
) -> tuple[list[list[SourcePeak]], int]:
    """Merge close clusters unless one context independently resolves both."""

    if minimum_summit_separation < 0:
        raise ValueError("minimum_summit_separation must be non-negative")
    ordered = sorted(clusters, key=lambda members: _representative_summit(members))
    merges = 0
    while True:
        eligible_pairs: list[tuple[int, int, int]] = []
        for index, (left, right) in enumerate(zip(ordered, ordered[1:])):
            left_summit = _representative_summit(left)
            right_summit = _representative_summit(right)
            context_overlap = {peak.context for peak in left} & {
                peak.context for peak in right
            }
            combined_summits = [int(peak.summit) for peak in left + right]
            gap = right_summit - left_summit
            if (
                gap < minimum_summit_separation
                and not context_overlap
                and max(combined_summits) - min(combined_summits)
                <= summit_max_distance
            ):
                eligible_pairs.append((gap, left_summit, index))
        if not eligible_pairs:
            break
        index = min(eligible_pairs)[2]
        ordered[index] = ordered[index] + ordered[index + 1]
        del ordered[index + 1]
        merges += 1
    return ordered, merges


def assemble_master_dhs(
    source_peaks: Iterable[SourcePeak],
    *,
    chrom_sizes: dict[str, int],
    summit_max_distance: int = 150,
    minimum_summit_separation: int = 50,
) -> tuple[list[MasterDHS], int, int]:
    """Cluster compatible source peaks and return non-overlapping master DHSs."""

    if summit_max_distance < 0:
        raise ValueError("summit_max_distance must be non-negative")
    grouped: dict[str, list[SourcePeak]] = {chrom: [] for chrom in chrom_sizes}
    for peak in source_peaks:
        if peak.chrom not in chrom_sizes:
            raise ValueError(f"Unknown chromosome: {peak.chrom}")
        if peak.summit is None or not peak.start <= peak.summit < peak.end:
            raise ValueError(f"Peak {peak.context}:{peak.name} has an invalid summit")
        grouped[peak.chrom].append(peak)

    masters: list[MasterDHS] = []
    clipped_boundaries = 0
    close_cluster_merges = 0
    for chrom in chrom_sizes:
        clusters: list[list[SourcePeak]] = []
        bin_width = max(1, summit_max_distance + 1)
        summit_bins: dict[int, set[int]] = {}
        for peak in sorted(
            grouped[chrom],
            key=lambda item: (item.width, int(item.summit), item.start, item.end, item.context, item.name),
        ):
            eligible: list[tuple[float, int, int]] = []
            peak_summit = int(peak.summit)
            first_bin = (peak_summit - summit_max_distance) // bin_width
            last_bin = (peak_summit + summit_max_distance) // bin_width
            candidate_indices = set().union(
                *(summit_bins.get(bin_id, set()) for bin_id in range(first_bin, last_bin + 1))
            )
            for index in sorted(candidate_indices):
                members = clusters[index]
                if any(member.context == peak.context for member in members):
                    continue
                summits = [int(member.summit) for member in members] + [peak_summit]
                if max(summits) - min(summits) > summit_max_distance:
                    continue
                if not any(
                    _reciprocally_contains_summits(member, peak) for member in members
                ):
                    continue
                cluster_median = statistics.median(int(member.summit) for member in members)
                envelope_width = max(member.end for member in members) - min(
                    member.start for member in members
                )
                eligible.append((abs(peak_summit - cluster_median), envelope_width, index))
            if eligible:
                selected_index = min(eligible)[2]
                clusters[selected_index].append(peak)
            else:
                clusters.append([peak])
                selected_index = len(clusters) - 1
            summit_bins.setdefault(peak_summit // bin_width, set()).add(selected_index)

        clusters, chromosome_close_merges = _merge_close_context_shifted_clusters(
            clusters,
            summit_max_distance=summit_max_distance,
            minimum_summit_separation=minimum_summit_separation,
        )
        close_cluster_merges += chromosome_close_merges

        chromosome_masters = [
            MasterDHS(
                chrom=chrom,
                start=min(peak.start for peak in members),
                end=max(peak.end for peak in members),
                summit=_representative_summit(members),
                members=sorted(members, key=lambda peak: (peak.context, peak.name)),
            )
            for members in clusters
        ]
        chromosome_masters.sort(key=lambda item: (item.summit, item.start, item.end))
        for left, right in zip(chromosome_masters, chromosome_masters[1:]):
            if left.summit >= right.summit:
                raise ValueError(
                    f"Distinct master DHSs on {chrom} have non-increasing summits at "
                    f"{left.summit} and {right.summit}; members are "
                    f"{[(peak.context, peak.name, peak.summit) for peak in left.members]} and "
                    f"{[(peak.context, peak.name, peak.summit) for peak in right.members]}"
                )
            if left.end > right.start:
                boundary = (left.summit + right.summit + 1) // 2
                left.end = min(left.end, boundary)
                right.start = max(right.start, boundary)
                clipped_boundaries += 1
        masters.extend(chromosome_masters)

    for index, master in enumerate(masters, start=1):
        if not master.start <= master.summit < master.end:
            raise ValueError(f"Boundary clipping excluded summit at {master.chrom}:{master.summit}")
        master.master_id = f"DHS{index:07d}"
    return masters, clipped_boundaries, close_cluster_merges


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


def _median_width(masters: list[MasterDHS]) -> float:
    return float(statistics.median(master.end - master.start for master in masters)) if masters else 0.0


def build_master_registry(
    *,
    context_peaks: dict[str, Path],
    context_signals: dict[str, Path],
    chrom_sizes_path: Path,
    output_bed: Path,
    summit_bed: Path,
    membership_tsv: Path,
    context_matrix_tsv: Path,
    stats_json: Path,
    summit_max_distance: int = 150,
    minimum_summit_separation: int = 50,
    workers: int = 1,
) -> dict[str, Any]:
    """Build and atomically write the complete cross-context DHS registry."""

    if not context_peaks or set(context_peaks) != set(context_signals):
        raise ValueError("Peak and signal inputs must name the same non-empty contexts")
    if workers < 1:
        raise ValueError("workers must be positive")
    chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    contexts = list(context_peaks)

    if workers == 1:
        by_context = {
            context: load_context_peaks(
                context,
                context_peaks[context],
                context_signals[context],
                chrom_sizes,
            )
            for context in contexts
        }
    else:
        jobs = [
            (context, context_peaks[context], context_signals[context], chrom_sizes)
            for context in contexts
        ]
        with ProcessPoolExecutor(max_workers=min(workers, len(contexts))) as executor:
            loaded = executor.map(_load_context_job, jobs)
            by_context = dict(zip(contexts, loaded))
    all_peaks = [peak for context in contexts for peak in by_context[context]]
    masters, clipped_boundaries, close_cluster_merges = assemble_master_dhs(
        all_peaks,
        chrom_sizes=chrom_sizes,
        summit_max_distance=summit_max_distance,
        minimum_summit_separation=minimum_summit_separation,
    )

    outputs = [output_bed, summit_bed, membership_tsv, context_matrix_tsv, stats_json]
    writers = [_atomic_writer(path) for path in outputs]
    handles = [item[0] for item in writers]
    temporary_paths = [item[1] for item in writers]
    bed_handle, summit_handle, membership_handle, matrix_handle, stats_handle = handles
    try:
        membership_handle.write(
            "master_dhs_id\tcontext_id\tsource_peak_id\tchrom\tinput_start\tinput_end\t"
            "source_start\tsource_end\tcoordinate_clipped\tsource_summit\t"
            "summit_signal\tsummit_fallback\tsource_score\tpeak_method\t"
            "support_n\treplicate_n\tsupport_fraction\tsupporting_libraries\n"
        )
        matrix_handle.write(
            "master_dhs_id\tchrom\tstart\tend\tsummit\tcontext_n\t"
            + "\t".join(contexts)
            + "\n"
        )
        for master in masters:
            member_contexts = sorted({peak.context for peak in master.members})
            bed_handle.write(
                f"{master.chrom}\t{master.start}\t{master.end}\t"
                f"{master.master_id}\t0\t.\n"
            )
            summit_handle.write(
                f"{master.chrom}\t{master.summit}\t{master.summit + 1}\t"
                f"{master.master_id}\t0\t.\n"
            )
            presence = {context: 0 for context in contexts}
            for peak in master.members:
                presence[peak.context] += 1
                membership_handle.write(
                    f"{master.master_id}\t{peak.context}\t{peak.name}\t{peak.chrom}\t"
                    f"{peak.input_start if peak.input_start is not None else peak.start}\t"
                    f"{peak.input_end if peak.input_end is not None else peak.end}\t"
                    f"{peak.start}\t{peak.end}\t{int(peak.coordinate_clipped)}\t"
                    f"{peak.summit}\t{peak.summit_signal:.12g}\t"
                    f"{int(peak.summit_fallback)}\t{peak.score}\t{peak.peak_method}\t"
                    f"{peak.support_n}\t{peak.replicate_n}\t{peak.support_fraction:.12g}\t"
                    f"{peak.supporting_libraries}\n"
                )
            matrix_handle.write(
                f"{master.master_id}\t{master.chrom}\t{master.start}\t{master.end}\t"
                f"{master.summit}\t{len(member_contexts)}\t"
                + "\t".join(str(presence[context]) for context in contexts)
                + "\n"
            )

        widths = [master.end - master.start for master in masters]
        metrics: dict[str, Any] = {
            "status": "ok" if masters else "no_master_dhs",
            "method": "reciprocal_summit_complete_linkage_v2",
            "contexts": contexts,
            "context_peak_counts": {
                context: len(by_context[context]) for context in contexts
            },
            "source_peak_count": len(all_peaks),
            "master_dhs_count": len(masters),
            "multi_context_master_dhs_count": sum(
                len({peak.context for peak in master.members}) > 1 for master in masters
            ),
            "summit_max_distance": summit_max_distance,
            "minimum_summit_separation": minimum_summit_separation,
            "close_context_shifted_cluster_merges": close_cluster_merges,
            "summit_fallback_count": sum(peak.summit_fallback for peak in all_peaks),
            "coordinate_clipped_source_peak_count": sum(
                peak.coordinate_clipped for peak in all_peaks
            ),
            "clipped_boundary_pairs": clipped_boundaries,
            "master_width_min": min(widths) if widths else 0,
            "master_width_median": _median_width(masters),
            "master_width_max": max(widths) if widths else 0,
        }
        json.dump(metrics, stats_handle, indent=2, sort_keys=True)
        stats_handle.write("\n")
        for handle in handles:
            handle.close()
        for temporary, output in zip(temporary_paths, outputs):
            os.replace(temporary, output)
    finally:
        for handle in handles:
            if not handle.closed:
                handle.close()
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)
    return metrics
