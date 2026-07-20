import json
from pathlib import Path

import pytest

from short_read_processing.cpm_refinement import (
    SignalInterval,
    progressively_refine_cpm,
)


def test_progressive_cpm_refinement_splits_a_long_low_signal_region():
    intervals = [
        SignalInterval.from_bigwig("chr1", 0, 100, 1.0),
        SignalInterval.from_bigwig("chr1", 100, 500, 2.0),
        SignalInterval.from_bigwig("chr1", 500, 600, 1.0),
        SignalInterval.from_bigwig("chr1", 700, 760, 1.0),
    ]

    refined, excluded, thresholds = progressively_refine_cpm(intervals)

    assert [(item.start, item.end, item.selection_cutoff_cpm) for item in refined] == [
        (700, 760, 1.0),
        (100, 500, 2.0),
    ]
    assert excluded == []
    assert thresholds == [1.0, 2.0]


def test_progressive_cpm_refinement_enforces_minimum_mean_cpm():
    intervals = [
        SignalInterval.from_bigwig("chr1", 0, 100, 1.0),
        SignalInterval.from_bigwig("chr1", 100, 500, 2.0),
        SignalInterval.from_bigwig("chr1", 500, 600, 1.0),
        SignalInterval.from_bigwig("chr1", 700, 760, 1.0),
    ]

    refined, excluded, thresholds = progressively_refine_cpm(
        intervals,
        minimum_mean_cpm=2.0,
    )

    assert [(item.start, item.end, item.mean_cpm) for item in refined] == [
        (100, 500, 2.0),
    ]
    assert excluded == []
    assert thresholds == [1.0, 2.0]


@pytest.mark.parametrize("minimum_mean_cpm", [-1.0, float("nan"), float("inf")])
def test_progressive_cpm_refinement_rejects_invalid_minimum_mean_cpm(
    minimum_mean_cpm,
):
    with pytest.raises(ValueError):
        progressively_refine_cpm([], minimum_mean_cpm=minimum_mean_cpm)


def test_cpm_bigwig_refinement_writes_signal_scores(tmp_path: Path):
    pybigwig = pytest.importorskip("pyBigWig")
    from short_read_processing.cpm_refinement import refine_cpm_bigwig

    signal = tmp_path / "signal.bw"
    bigwig = pybigwig.open(str(signal), "w")
    bigwig.addHeader([("chr1", 1000)])
    bigwig.addEntries(
        ["chr1", "chr1", "chr1", "chr1"],
        [0, 100, 500, 700],
        ends=[100, 500, 600, 760],
        values=[1.0, 2.0, 1.0, 1.0],
    )
    bigwig.close()
    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text("chr1\t0\t600\tpeak1\nchr1\t700\t760\tpeak2\n")
    output = tmp_path / "refined.bed"
    excluded = tmp_path / "excluded.bed"
    stats = tmp_path / "stats.json"

    metrics = refine_cpm_bigwig(
        peaks=peaks,
        signal_bigwig=signal,
        output=output,
        excluded=excluded,
        stats=stats,
        minimum_mean_cpm=1.5,
    )

    assert output.read_text().splitlines() == [
        "chr1\t100\t500\tCPM_refined_1\t1000\t.\t2\t2\t2",
    ]
    assert excluded.read_text() == ""
    assert metrics["candidate_macs3_peaks"] == 2
    assert metrics["refined_intervals"] == 1
    assert metrics["minimum_mean_cpm"] == 1.5
    assert json.loads(stats.read_text()) == metrics
