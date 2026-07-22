import xml.etree.ElementTree as ET
from pathlib import Path

from build_igv_session import build_session


def touch_all(paths):
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()


def test_session_contains_final_qpois_consensus_and_chip_tracks(tmp_path: Path):
    atac = tmp_path / "atac"
    condition = atac / "conditions" / "e5"
    touch_all(
        [
            condition / "tracks/e5.MACS3-pileup.unscaled.bw",
            condition / "tracks/e5.qpois.bw",
            condition / "peaks/e5.candidates.narrowPeak",
            condition / "peaks/e5.qpois-refined.bed",
            condition / "peaks/e5.replicate-supported.bed",
        ]
    )
    chip = tmp_path / "chip"
    touch_all(
        [
            chip / "tracks/e5_h3k27ac_rep1.CPM.bw",
            chip / "peaks/e5_h3k27ac_rep1/e5_h3k27ac_rep1_peaks.broadPeak",
        ]
    )
    output = tmp_path / "session.xml"

    counts = build_session(atac, output, "dm6", "All", chip)

    root = ET.parse(output).getroot()
    resources = root.findall("./Resources/Resource")
    tracks = root.findall("./Panel/Track")
    assert counts == (1, 1, 7)
    assert len(resources) == len(tracks) == 7
    assert all((output.parent / item.attrib["path"]).is_file() for item in resources)
    assert [track.attrib["name"] for track in tracks][-2:] == [
        "e5_h3k27ac_rep1 | ChIP CPM signal",
        "e5_h3k27ac_rep1 | MACS3 peaks",
    ]
