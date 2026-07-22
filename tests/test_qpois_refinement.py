import json
from pathlib import Path

from short_read_processing.qpois_refinement import ALGORITHM, run_refinement


def test_qpois_refinement_splits_a_broad_candidate_at_higher_signal(tmp_path: Path):
    candidates = tmp_path / "candidates.narrowPeak"
    candidates.write_text(
        "chr1\t0\t250\tfirst\t10\t.\n"
        "chr1\t300\t900\tsecond\t10\t.\n"
    )
    qpois = tmp_path / "qpois.bdg"
    qpois.write_text(
        "chr1\t0\t100\t3\n"
        "chr1\t100\t200\t3\n"
        "chr1\t300\t400\t5\n"
        "chr1\t400\t700\t3\n"
        "chr1\t700\t800\t5\n"
    )
    output = tmp_path / "refined.bed"
    excluded = tmp_path / "excluded.bed"
    stats = tmp_path / "stats.json"

    metrics = run_refinement(
        qpois_bedgraph=qpois,
        candidate_peaks=candidates,
        output_bed=output,
        excluded_bed=excluded,
        stats_json=stats,
        name_prefix="test",
        maximum_length=400,
    )

    rows = [line.split("\t") for line in output.read_text().splitlines()]
    assert [(int(row[1]), int(row[2])) for row in rows] == [
        (0, 200),
        (300, 400),
        (700, 800),
    ]
    assert [int(row[7]) for row in rows] == [2, 3, 3]
    assert all(len(row) == 8 for row in rows)
    assert excluded.read_text() == ""
    assert metrics["algorithm"] == ALGORITHM
    assert metrics["refined_peaks"] == 3
    assert json.loads(stats.read_text()) == metrics


def test_qpois_refinement_writes_valid_empty_outputs(tmp_path: Path):
    candidates = tmp_path / "candidates.narrowPeak"
    candidates.write_text("")
    qpois = tmp_path / "qpois.bdg"
    qpois.write_text("")
    output = tmp_path / "refined.bed"
    excluded = tmp_path / "excluded.bed"
    stats = tmp_path / "stats.json"

    metrics = run_refinement(
        qpois_bedgraph=qpois,
        candidate_peaks=candidates,
        output_bed=output,
        excluded_bed=excluded,
        stats_json=stats,
        name_prefix="empty",
    )

    assert output.read_text() == excluded.read_text() == ""
    assert metrics["status"] == "no_refined_peaks"
    assert metrics["candidate_peaks"] == 0
