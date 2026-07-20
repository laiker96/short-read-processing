import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from build_igv_session import build_session


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
