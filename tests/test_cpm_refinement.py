import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from short_read_processing.cpm_refinement import (
    SignalInterval,
    progressively_refine_cpm,
    refine_cpm_bigwig,
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
        (100, 500, 2.0),
        (700, 760, 1.0),
    ]
    assert excluded == []
    assert thresholds == [2.0, 1.0]


def test_descending_refinement_preserves_modes_before_low_signal_bridge():
    intervals = [
        SignalInterval.from_bigwig("chr1", 0, 100, 3.0),
        SignalInterval.from_bigwig("chr1", 100, 200, 1.0),
        SignalInterval.from_bigwig("chr1", 200, 300, 3.0),
    ]

    refined, excluded, thresholds = progressively_refine_cpm(
        intervals,
        minimum_mean_cpm=2.0,
    )

    assert [(item.start, item.end, item.selection_cutoff_cpm) for item in refined] == [
        (0, 100, 3.0),
        (200, 300, 3.0),
    ]
    assert excluded == []
    assert thresholds == [3.0, 1.0]

    reversed_refined, reversed_excluded, reversed_thresholds = progressively_refine_cpm(
        list(reversed(intervals)),
        minimum_mean_cpm=2.0,
    )
    assert reversed_refined == refined
    assert reversed_excluded == excluded
    assert reversed_thresholds == thresholds


def test_descending_refinement_expands_one_mode_to_its_last_valid_boundary():
    intervals = [
        SignalInterval.from_bigwig("chr1", 0, 100, 3.0),
        SignalInterval.from_bigwig("chr1", 100, 200, 2.0),
        SignalInterval.from_bigwig("chr1", 200, 300, 1.0),
    ]

    refined, excluded, thresholds = progressively_refine_cpm(
        intervals,
        minimum_mean_cpm=2.1,
    )

    assert [(item.start, item.end, item.selection_cutoff_cpm) for item in refined] == [
        (0, 200, 2.0),
    ]
    assert excluded == []
    assert thresholds == [3.0, 2.0, 1.0]


def test_prominence_filter_merges_modes_across_a_shallow_saddle():
    intervals = [
        SignalInterval.from_bigwig("chr1", 0, 100, 3.0),
        SignalInterval.from_bigwig("chr1", 100, 200, 2.4),
        SignalInterval.from_bigwig("chr1", 200, 300, 2.8),
    ]

    refined, excluded, thresholds = progressively_refine_cpm(intervals)

    assert [(item.start, item.end, item.selection_cutoff_cpm) for item in refined] == [
        (0, 300, 2.4),
    ]
    assert excluded == []
    assert thresholds == [3.0, 2.8, 2.4]

    unfiltered, _, _ = progressively_refine_cpm(
        intervals,
        minimum_mode_prominence=0.0,
    )
    assert [(item.start, item.end) for item in unfiltered] == [
        (0, 100),
        (200, 300),
    ]


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
    assert [(item.start, item.end, item.mean_cpm) for item in excluded] == [
        (700, 760, 1.0),
    ]
    assert thresholds == [2.0, 1.0]


@pytest.mark.parametrize("minimum_mean_cpm", [-1.0, float("nan"), float("inf")])
def test_progressive_cpm_refinement_rejects_invalid_minimum_mean_cpm(
    minimum_mean_cpm,
):
    with pytest.raises(ValueError):
        progressively_refine_cpm([], minimum_mean_cpm=minimum_mean_cpm)


@pytest.mark.parametrize(
    "minimum_mode_prominence",
    [-0.01, 1.01, float("nan"), float("inf")],
)
def test_progressive_cpm_refinement_rejects_invalid_mode_prominence(
    minimum_mode_prominence,
):
    with pytest.raises(ValueError):
        progressively_refine_cpm(
            [], minimum_mode_prominence=minimum_mode_prominence
        )


def test_empty_candidate_peaks_write_empty_outputs_and_status(
    tmp_path: Path, monkeypatch
):
    class FakeBigWig:
        def chroms(self):
            return {"chr1": 1000}

        def close(self):
            pass

    monkeypatch.setitem(
        sys.modules,
        "pyBigWig",
        SimpleNamespace(open=lambda _: FakeBigWig()),
    )
    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text("")
    signal = tmp_path / "signal.bw"
    signal.touch()
    output = tmp_path / "refined.bed"
    excluded = tmp_path / "excluded.bed"
    stats = tmp_path / "stats.json"

    metrics = refine_cpm_bigwig(
        peaks=peaks,
        signal_bigwig=signal,
        output=output,
        excluded=excluded,
        stats=stats,
    )

    assert output.read_text() == ""
    assert excluded.read_text() == ""
    assert metrics["status"] == "no_candidate_peaks"
    assert metrics["candidate_macs3_peaks"] == 0
    assert metrics["refined_intervals"] == 0
    assert metrics["minimum_positive_cpm"] == 0.0
    assert metrics["maximum_cpm"] == 0.0
    assert json.loads(stats.read_text()) == metrics


def test_cpm_bigwig_refinement_writes_signal_scores(tmp_path: Path):
    pybigwig = pytest.importorskip("pyBigWig")

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
    assert metrics["status"] == "ok"
    assert metrics["refinement_algorithm"] == "descending_cpm_prominence_watershed_v1"
    assert metrics["minimum_mean_cpm"] == 1.5
    assert metrics["minimum_mode_prominence"] == 0.25
    assert metrics["threshold_rule"].startswith(
        "all observed positive CPM levels, descending"
    )
    assert json.loads(stats.read_text()) == metrics
