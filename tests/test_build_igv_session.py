import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from build_igv_session import build_condition_atlas_session, build_session


@pytest.mark.parametrize("refined_directory", ["refined", "cpm_refinement/refined"])
def test_build_session_groups_atac_and_h3k27ac_by_context(
    tmp_path: Path, refined_directory: str
):
    results = tmp_path / "atac"
    h3k27ac = results / "h3k27ac_tracks"
    for sample in ("e5_atac_rep1", "e5_atac_rep2"):
        bigwig = results / "tracks" / f"{sample}.fragments-lt150.Tn5-shifted.CPM.bw"
        narrowpeak = results / "macs3" / sample / f"{sample}_peaks.narrowPeak"
        refined = results / refined_directory / f"{sample}.CPM-refined.bed"
        for path in (bigwig, narrowpeak, refined):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    h3k27ac.mkdir()
    (h3k27ac / "e5_h3k27ac_rep1.CPM.bw").touch()
    output = results / "session.xml"

    counts = build_session(results, output, "dm6", "All", h3k27ac)

    root = ET.parse(output).getroot()
    resources = root.findall("./Resources/Resource")
    tracks = root.findall("./Panel/Track")
    assert counts == (2, 1)
    assert len(resources) == len(tracks) == 7
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks][-1] == (
        "E5 H3K27ac rep 1 | CPM signal"
    )
    assert all("qpois" not in item.attrib["id"].lower() for item in tracks)


@pytest.mark.parametrize("include_dhs_driven,track_count", [(False, 6), (True, 15)])
def test_build_condition_atlas_session_adds_tissues_and_global_boundaries(
    tmp_path: Path,
    include_dhs_driven: bool,
    track_count: int,
):
    results = tmp_path / "atlas"
    for condition in ("e5", "wid"):
        bigwig = (
            results
            / "conditions"
            / condition
            / "tracks"
            / f"{condition}.fragments-lt150.Tn5-shifted.CPM.bw"
        )
        consensus = results / "conditions" / condition / f"{condition}.consensus.bed"
        for path in (bigwig, consensus):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
    (results / "atlas.peaks.bed").touch()
    (results / "atlas.variable.peaks.bed").touch()
    (results / "atlas.stats.json").write_text('{"conditions": ["wid", "e5"]}\n')
    if include_dhs_driven:
        (results / "atlas.narrow-first.anchors250.bed").touch()
        (results / "atlas.narrow-first.variable.peaks.bed").touch()
        (results / "atlas.dhs-support-fraction.bw").touch()
        (results / "atlas.fwhm-boundaries.bed").touch()
        (results / "atlas.center-mode-half-prominence-boundaries.bed").touch()
        (results / "atlas.dhs-driven.anchors250.bed").touch()
        (results / "atlas.dhs-driven.peaks.bed").touch()
        (results / "atlas.dhs-driven.aggregate-shape.bw").touch()
        (results / "atlas.dhs-driven.signal-shaped.peaks.bed").touch()
    output = results / "atlas.igv.xml"

    counts = build_condition_atlas_session(
        results,
        output,
        "dm6",
        "All",
        include_dhs_driven,
    )

    root = ET.parse(output).getroot()
    resources = root.findall("./Resources/Resource")
    tracks = root.findall("./Panel/Track")
    assert counts == (2, track_count)
    assert len(resources) == len(tracks) == track_count
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert tracks[0].attrib["name"].startswith("WID |")
    expected_global = [
        "Global ATAC atlas | fixed 250-bp anchors",
        "Global ATAC atlas | variable median boundaries",
    ]
    if include_dhs_driven:
        expected_global.extend(
            [
                "Global ATAC atlas | narrow-source-first 250-bp anchors",
                "Global ATAC atlas | narrow-source-first median boundaries",
                "Global ATAC atlas | condition DHS support fraction",
                "Global ATAC atlas | DHS-support FWHM boundaries",
                "Global ATAC atlas | center-mode half-prominence boundaries",
                "DHS-driven atlas | fixed 250-bp measurement anchors",
                "DHS-driven atlas | variable DHS consensus boundaries",
                "DHS-driven atlas | contributor-normalized aggregate shape",
                "DHS-driven atlas | signal-shaped boundaries",
            ]
        )
    assert [track.attrib["name"] for track in tracks][-len(expected_global):] == (
        expected_global
    )
