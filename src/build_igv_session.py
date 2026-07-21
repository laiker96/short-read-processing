#!/usr/bin/env python3
"""Build a portable IGV session for ATAC peaks and optional H3K27ac signal."""

from __future__ import annotations

import argparse
import json
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path


ATAC_BIGWIG_SUFFIX = ".fragments-lt150.Tn5-shifted.CPM.bw"
H3K27AC_BIGWIG_SUFFIX = ".CPM.bw"
ATAC_SAMPLE_RE = re.compile(
    r"^(?P<context>.+)_atac_rep(?P<replicate>[1-9][0-9]*)$"
)
H3K27AC_SAMPLE_RE = re.compile(
    r"^(?P<context>.+)_h3k27ac_rep(?P<replicate>[1-9][0-9]*)$"
)


def sample_parts(sample: str, pattern: re.Pattern[str]) -> tuple[str, str]:
    match = pattern.fullmatch(sample)
    if not match:
        raise ValueError(f"Unrecognized sample name: {sample}")
    return match.group("context"), match.group("replicate")


def relative_path(path: Path, session_path: Path) -> str:
    return Path(os.path.relpath(path, session_path.parent)).as_posix()


def add_track(
    resources: ET.Element,
    panel: ET.Element,
    *,
    path: Path,
    output: Path,
    name: str,
    color: str,
    renderer: str,
    height: str,
    window_function: str,
) -> None:
    track_id = relative_path(path, output)
    ET.SubElement(resources, "Resource", path=track_id)
    ET.SubElement(
        panel,
        "Track",
        color=color,
        expand="false",
        height=height,
        id=track_id,
        name=name,
        renderer=renderer,
        visible="true",
        windowFunction=window_function,
    )


def build_session(
    results: Path,
    output: Path,
    genome: str,
    locus: str,
    h3k27ac_tracks: Path | None = None,
) -> tuple[int, int]:
    bigwigs = sorted((results / "tracks").glob(f"*{ATAC_BIGWIG_SUFFIX}"))
    if not bigwigs:
        raise ValueError(f"No short-fragment CPM BigWigs found under {results / 'tracks'}")

    atac_by_context: dict[str, list[tuple[str, str, Path, Path, Path]]] = {}
    missing: list[Path] = []
    for bigwig in bigwigs:
        sample = bigwig.name.removesuffix(ATAC_BIGWIG_SUFFIX)
        context, replicate = sample_parts(sample, ATAC_SAMPLE_RE)
        narrowpeak = results / "macs3" / sample / f"{sample}_peaks.narrowPeak"
        refined = results / "refined" / f"{sample}.CPM-refined.bed"
        if not refined.is_file():
            refined = (
                results
                / "cpm_refinement"
                / "refined"
                / f"{sample}.CPM-refined.bed"
            )
        missing.extend(path for path in (narrowpeak, refined) if not path.is_file())
        atac_by_context.setdefault(context, []).append(
            (sample, replicate, bigwig, narrowpeak, refined)
        )
    if missing:
        paths = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Missing IGV input files:\n{paths}")

    h3k27ac_by_context: dict[str, list[tuple[str, Path]]] = {}
    if h3k27ac_tracks is not None:
        h3k27ac_bigwigs = sorted(
            h3k27ac_tracks.glob(f"*{H3K27AC_BIGWIG_SUFFIX}")
        )
        if not h3k27ac_bigwigs:
            raise ValueError(f"No H3K27ac CPM BigWigs found under {h3k27ac_tracks}")
        for bigwig in h3k27ac_bigwigs:
            sample = bigwig.name.removesuffix(H3K27AC_BIGWIG_SUFFIX)
            context, replicate = sample_parts(sample, H3K27AC_SAMPLE_RE)
            h3k27ac_by_context.setdefault(context, []).append((replicate, bigwig))

    output.parent.mkdir(parents=True, exist_ok=True)
    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="true",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name="ATAC and H3K27ac")

    h3k27ac_count = 0
    for context in sorted(atac_by_context):
        context_label = context.upper()
        for sample, replicate, bigwig, narrowpeak, refined in atac_by_context[context]:
            label = f"{context_label} ATAC rep {replicate}"
            add_track(
                resources,
                panel,
                path=bigwig,
                output=output,
                name=f"{label} | CPM signal (<150 bp)",
                color="49,104,174",
                renderer="BAR_CHART",
                height="42",
                window_function="mean",
            )
            add_track(
                resources,
                panel,
                path=narrowpeak,
                output=output,
                name=f"{label} | MACS3 narrowPeak",
                color="110,110,110",
                renderer="BASIC_FEATURE",
                height="20",
                window_function="count",
            )
            add_track(
                resources,
                panel,
                path=refined,
                output=output,
                name=f"{label} | CPM-refined peaks",
                color="0,145,130",
                renderer="BASIC_FEATURE",
                height="24",
                window_function="count",
            )
        for replicate, bigwig in h3k27ac_by_context.get(context, []):
            add_track(
                resources,
                panel,
                path=bigwig,
                output=output,
                name=f"{context_label} H3K27ac rep {replicate} | CPM signal",
                color="221,126,32",
                renderer="BAR_CHART",
                height="42",
                window_function="mean",
            )
            h3k27ac_count += 1

    ET.indent(session, space="  ")
    ET.ElementTree(session).write(output, encoding="utf-8", xml_declaration=True)
    return sum(map(len, atac_by_context.values())), h3k27ac_count


def build_condition_atlas_session(
    results: Path,
    output: Path,
    genome: str,
    locus: str,
    include_dhs_driven: bool = False,
) -> tuple[int, int]:
    """Build an IGV session for pooled condition signal and atlas boundaries."""

    conditions_root = results / "conditions"
    condition_ids = (
        sorted(path.name for path in conditions_root.iterdir() if path.is_dir())
        if conditions_root.is_dir()
        else []
    )
    stats_path = results / "atlas.stats.json"
    if stats_path.is_file():
        configured_order = json.loads(stats_path.read_text(encoding="utf-8")).get(
            "conditions", []
        )
        if set(configured_order) == set(condition_ids):
            condition_ids = configured_order
    if not condition_ids:
        raise ValueError(f"No condition directories found under {conditions_root}")

    condition_tracks: list[tuple[str, Path, Path]] = []
    for condition in condition_ids:
        condition_root = conditions_root / condition
        bigwig = condition_root / "tracks" / f"{condition}{ATAC_BIGWIG_SUFFIX}"
        consensus = condition_root / f"{condition}.consensus.bed"
        condition_tracks.append((condition, bigwig, consensus))
    fixed = results / "atlas.peaks.bed"
    variable = results / "atlas.variable.peaks.bed"
    narrow_first = results / "atlas.narrow-first.anchors250.bed"
    narrow_first_variable = results / "atlas.narrow-first.variable.peaks.bed"
    dhs_support = results / "atlas.dhs-support-fraction.bw"
    fwhm_boundaries = results / "atlas.fwhm-boundaries.bed"
    center_mode = results / "atlas.center-mode-half-prominence-boundaries.bed"
    dhs_anchors = results / "atlas.dhs-driven.anchors250.bed"
    dhs_peaks = results / "atlas.dhs-driven.peaks.bed"
    dhs_aggregate = results / "atlas.dhs-driven.aggregate-shape.bw"
    dhs_shaped = results / "atlas.dhs-driven.signal-shaped.peaks.bed"
    missing = [
        path
        for _condition, bigwig, consensus in condition_tracks
        for path in (bigwig, consensus)
        if not path.is_file()
    ]
    missing.extend(path for path in (fixed, variable) if not path.is_file())
    if include_dhs_driven:
        missing.extend(
            path
            for path in (
                dhs_support,
                fwhm_boundaries,
                center_mode,
                narrow_first,
                narrow_first_variable,
                dhs_anchors,
                dhs_peaks,
                dhs_aggregate,
                dhs_shaped,
            )
            if not path.is_file()
        )
    if missing:
        paths = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"Missing IGV input files:\n{paths}")

    output.parent.mkdir(parents=True, exist_ok=True)
    session = ET.Element(
        "Session",
        genome=genome,
        hasGeneTrack="true",
        hasSequenceTrack="true",
        locus=locus,
        version="3",
    )
    resources = ET.SubElement(session, "Resources")
    panel = ET.SubElement(session, "Panel", name="ATAC condition consensus atlas")

    for condition, bigwig, consensus in condition_tracks:
        label = condition.upper()
        add_track(
            resources,
            panel,
            path=bigwig,
            output=output,
            name=f"{label} | pooled ATAC CPM signal (<150 bp)",
            color="49,104,174",
            renderer="BAR_CHART",
            height="48",
            window_function="mean",
        )
        add_track(
            resources,
            panel,
            path=consensus,
            output=output,
            name=f"{label} | replicate-supported consensus peaks",
            color="0,145,130",
            renderer="BASIC_FEATURE",
            height="24",
            window_function="count",
        )

    add_track(
        resources,
        panel,
        path=fixed,
        output=output,
        name="Global ATAC atlas | fixed 250-bp anchors",
        color="202,61,52",
        renderer="BASIC_FEATURE",
        height="28",
        window_function="count",
    )

    add_track(
        resources,
        panel,
        path=variable,
        output=output,
        name="Global ATAC atlas | variable median boundaries",
        color="117,74,147",
        renderer="BASIC_FEATURE",
        height="28",
        window_function="count",
    )
    if include_dhs_driven:
        add_track(
            resources,
            panel,
            path=narrow_first,
            output=output,
            name="Global ATAC atlas | narrow-source-first 250-bp anchors",
            color="166,86,40",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )
        add_track(
            resources,
            panel,
            path=narrow_first_variable,
            output=output,
            name="Global ATAC atlas | narrow-source-first median boundaries",
            color="102,194,165",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )
        add_track(
            resources,
            panel,
            path=dhs_support,
            output=output,
            name="Global ATAC atlas | condition DHS support fraction",
            color="0,158,115",
            renderer="BAR_CHART",
            height="48",
            window_function="mean",
        )
        add_track(
            resources,
            panel,
            path=fwhm_boundaries,
            output=output,
            name="Global ATAC atlas | DHS-support FWHM boundaries",
            color="204,121,167",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )
        add_track(
            resources,
            panel,
            path=center_mode,
            output=output,
            name="Global ATAC atlas | center-mode half-prominence boundaries",
            color="128,64,0",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )
        add_track(
            resources,
            panel,
            path=dhs_anchors,
            output=output,
            name="DHS-driven atlas | fixed 250-bp measurement anchors",
            color="230,159,0",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )
        add_track(
            resources,
            panel,
            path=dhs_peaks,
            output=output,
            name="DHS-driven atlas | variable DHS consensus boundaries",
            color="86,180,233",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )
        add_track(
            resources,
            panel,
            path=dhs_aggregate,
            output=output,
            name="DHS-driven atlas | contributor-normalized aggregate shape",
            color="0,114,178",
            renderer="BAR_CHART",
            height="48",
            window_function="mean",
        )
        add_track(
            resources,
            panel,
            path=dhs_shaped,
            output=output,
            name="DHS-driven atlas | signal-shaped boundaries",
            color="213,94,0",
            renderer="BASIC_FEATURE",
            height="28",
            window_function="count",
        )

    ET.indent(session, space="  ")
    ET.ElementTree(session).write(output, encoding="utf-8", xml_declaration=True)
    return len(condition_tracks), len(condition_tracks) * 2 + 2 + 9 * int(
        include_dhs_driven
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="Export directory containing tracks/ and peaks")
    parser.add_argument("--output", type=Path, help="Output XML; defaults inside the export directory")
    parser.add_argument(
        "--h3k27ac-tracks",
        type=Path,
        help="Optional directory containing <context>_h3k27ac_rep<N>.CPM.bw",
    )
    parser.add_argument("--genome", default="dm6")
    parser.add_argument("--locus", default="All")
    parser.add_argument(
        "--condition-atlas",
        action="store_true",
        help="Use pooled condition tracks plus fixed and variable global atlas peaks",
    )
    parser.add_argument(
        "--include-dhs-driven",
        action="store_true",
        help="Add DHS-support, FWHM, and DHS-driven tracks to a condition-atlas session",
    )
    args = parser.parse_args()

    results = args.results.resolve()
    output = (
        args.output.resolve()
        if args.output
        else results / (
            (
                "atlas-dhs-driven-comparison.igv.xml"
                if args.include_dhs_driven
                else "atlas-condition-consensus.igv.xml"
            )
            if args.condition_atlas
            else "atlas-atac-short-fragments.igv.xml"
        )
    )
    if args.condition_atlas:
        if args.h3k27ac_tracks:
            parser.error("--h3k27ac-tracks is not supported with --condition-atlas")
        condition_count, track_count = build_condition_atlas_session(
            results,
            output,
            args.genome,
            args.locus,
            args.include_dhs_driven,
        )
        print(
            f"Wrote {condition_count} ATAC conditions and {track_count} tracks "
            f"to {output}"
        )
        return 0

    if args.include_dhs_driven:
        parser.error("--include-dhs-driven requires --condition-atlas")

    h3k27ac_tracks = args.h3k27ac_tracks.resolve() if args.h3k27ac_tracks else None
    atac_count, h3k27ac_count = build_session(
        results,
        output,
        args.genome,
        args.locus,
        h3k27ac_tracks,
    )
    print(
        f"Wrote {atac_count} ATAC samples, {h3k27ac_count} H3K27ac samples, "
        f"and {atac_count * 3 + h3k27ac_count} tracks to {output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
