#!/usr/bin/env python3
"""Build a portable IGV session for ATAC peaks and optional H3K27ac signal."""

from __future__ import annotations

import argparse
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
        refined = results / "cpm_refinement" / "refined" / f"{sample}.CPM-refined.bed"
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
    args = parser.parse_args()

    results = args.results.resolve()
    output = (
        args.output.resolve()
        if args.output
        else results / "atlas-atac-short-fragments.igv.xml"
    )
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
