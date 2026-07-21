import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.atlas import (
    AtlasCandidate,
    RefinedPeak,
    _build_variable_boundaries,
    _select_dhs_seeds,
    _select_narrow_source_first,
    build_condition_consensus,
    build_dhs_center_mode_half_prominence,
    build_dhs_support_fwhm,
    build_global_atlas,
    build_signal_shaped_atlas,
    read_condition_map,
)


def _refined(chrom, start, end, name, maximum=5):
    return f"{chrom}\t{start}\t{end}\t{name}\t1000\t.\t3\t{maximum}\t2\n"


def _condition_peak(condition, start, end, name, maximum=5):
    return (
        _refined("chr1", start, end, name, maximum).rstrip("\n")
        + f"\t{condition}\t2\t2\t1\t{condition}_rep1,{condition}_rep2\n"
    )


def test_condition_map_requires_every_sample_once(tmp_path: Path):
    condition_map = tmp_path / "conditions.tsv"
    condition_map.write_text(
        "condition_id\tcondition_label\tsample_id\n"
        "eye\tEye disc\teye_rep1\n"
        "eye\tEye disc\teye_rep2\n"
    )

    conditions = read_condition_map(
        condition_map,
        sample_ids=["eye_rep1", "eye_rep2"],
        minimum_replicates=2,
    )

    assert conditions[0].condition_id == "eye"
    assert conditions[0].samples == ("eye_rep1", "eye_rep2")

    with pytest.raises(AcquisitionError, match="does not assign"):
        read_condition_map(
            condition_map,
            sample_ids=["eye_rep1", "eye_rep2", "eye_rep3"],
            minimum_replicates=2,
        )


def test_condition_consensus_requires_half_peak_coverage_in_two_replicates(
    tmp_path: Path,
):
    pooled = tmp_path / "pooled.bed"
    pooled.write_text(
        _refined("chr1", 0, 100, "pooled_1")
        + _refined("chr1", 200, 300, "pooled_2")
    )
    replicate_1 = tmp_path / "rep1.bed"
    replicate_1.write_text(
        _refined("chr1", 0, 100, "rep1_1")
        + _refined("chr1", 200, 240, "rep1_2")
    )
    replicate_2 = tmp_path / "rep2.bed"
    replicate_2.write_text(
        _refined("chr1", 0, 60, "rep2_1")
        + _refined("chr1", 200, 300, "rep2_2")
    )
    output = tmp_path / "consensus.bed"
    support = tmp_path / "support.tsv"
    stats = tmp_path / "stats.json"

    metrics = build_condition_consensus(
        condition_id="eye",
        pooled_peaks=pooled,
        replicate_peaks={"eye_rep1": replicate_1, "eye_rep2": replicate_2},
        output_bed=output,
        support_tsv=support,
        stats_json=stats,
    )

    assert [line.split("\t")[3] for line in output.read_text().splitlines()] == [
        "pooled_1"
    ]
    assert support.read_text().splitlines()[1].split("\t")[2:4] == ["1", "2"]
    assert support.read_text().splitlines()[2].split("\t")[2:4] == ["0", "1"]
    assert metrics["retained_consensus_peaks"] == 1
    assert json.loads(stats.read_text()) == metrics


def test_iterative_atlas_avoids_daisy_chain_and_keeps_coverage(tmp_path, monkeypatch):
    signals = {
        "a.bw": [("chr1", 190, 210, 10.0)],
        "b.bw": [("chr1", 260, 280, 9.0)],
        "c.bw": [("chr1", 450, 470, 8.0)],
    }

    class FakeBigWig:
        def __init__(self, path):
            self.entries = signals[Path(path).name]

        def chroms(self):
            return {"chr1": 1000}

        def intervals(self, chrom, start, end):
            return [
                (interval_start, interval_end, value)
                for entry_chrom, interval_start, interval_end, value in self.entries
                if entry_chrom == chrom and interval_start < end and interval_end > start
            ]

        def close(self):
            pass

    monkeypatch.setitem(
        sys.modules,
        "pyBigWig",
        SimpleNamespace(open=lambda path: FakeBigWig(path)),
    )
    conditions = []
    for condition, start, end, maximum in (
        ("a", 100, 300, 10),
        ("b", 180, 380, 9),
        ("c", 360, 560, 8),
    ):
        bed = tmp_path / f"{condition}.bed"
        bed.write_text(_condition_peak(condition, start, end, f"{condition}_peak", maximum))
        bigwig = tmp_path / f"{condition}.bw"
        bigwig.touch()
        conditions.append((condition, bed, bigwig))

    outputs = {
        name: tmp_path / name
        for name in (
            "atlas.bed",
            "variable.bed",
            "membership.tsv",
            "presence.tsv",
            "coverage.tsv",
            "mean.tsv",
            "maximum.tsv",
            "stats.json",
        )
    }
    metrics = build_global_atlas(
        conditions=conditions,
        output_bed=outputs["atlas.bed"],
        variable_bed=outputs["variable.bed"],
        membership_tsv=outputs["membership.tsv"],
        presence_tsv=outputs["presence.tsv"],
        coverage_tsv=outputs["coverage.tsv"],
        mean_cpm_tsv=outputs["mean.tsv"],
        maximum_cpm_tsv=outputs["maximum.tsv"],
        stats_json=outputs["stats.json"],
        peak_width=250,
    )

    atlas = [line.split("\t") for line in outputs["atlas.bed"].read_text().splitlines()]
    assert [(int(row[1]), int(row[2])) for row in atlas] == [(75, 325), (335, 585)]
    assert all(int(row[2]) - int(row[1]) == 250 for row in atlas)
    assert int(atlas[0][2]) <= int(atlas[1][1])
    variable = [
        line.split("\t") for line in outputs["variable.bed"].read_text().splitlines()
    ]
    assert [(int(row[1]), int(row[2])) for row in variable] == [
        (140, 340),
        (360, 560),
    ]
    assert outputs["presence.tsv"].read_text().splitlines() == [
        "atlas_peak_id\ta\tb\tc",
        "ATAC_ATLAS_000001\t1\t1\t0",
        "ATAC_ATLAS_000002\t0\t0\t1",
    ]
    first_coverage = outputs["coverage.tsv"].read_text().splitlines()[1].split("\t")
    assert [float(value) for value in first_coverage[1:]] == pytest.approx([0.8, 0.58, 0])
    memberships = outputs["membership.tsv"].read_text().splitlines()
    assert any(line.startswith("ATAC_ATLAS_000001\tb\tb_peak\t") for line in memberships)
    assert metrics["source_condition_peaks"] == 3
    assert metrics["atlas_peaks"] == 2
    assert json.loads(outputs["stats.json"].read_text()) == metrics

    dhs_outputs = {
        name: tmp_path / f"dhs-{name}"
        for name in (
            "anchors.bed",
            "variable.bed",
            "membership.tsv",
            "presence.tsv",
            "coverage.tsv",
            "mean.tsv",
            "maximum.tsv",
            "stats.json",
        )
    }
    dhs_metrics = build_global_atlas(
        conditions=conditions,
        output_bed=dhs_outputs["anchors.bed"],
        variable_bed=dhs_outputs["variable.bed"],
        membership_tsv=dhs_outputs["membership.tsv"],
        presence_tsv=dhs_outputs["presence.tsv"],
        coverage_tsv=dhs_outputs["coverage.tsv"],
        mean_cpm_tsv=dhs_outputs["mean.tsv"],
        maximum_cpm_tsv=dhs_outputs["maximum.tsv"],
        stats_json=dhs_outputs["stats.json"],
        peak_width=250,
        grouping_method="dhs_seed",
    )
    assert dhs_metrics["grouping_method"] == "dhs_seed"
    assert dhs_metrics["atlas_peaks"] == 2
    assert dhs_metrics["presence_rule"] == (
        "source condition peak assigned by direct seed DHS grouping"
    )


def test_variable_boundaries_use_one_vote_per_condition_and_split_collisions():
    def candidate(condition, start, end, name, priority, fixed_start, fixed_end):
        peak = RefinedPeak(
            chrom="chr1",
            start=start,
            end=end,
            name=name,
            score=1000,
            strand=".",
            mean_cpm=3,
            maximum_cpm=10,
            selection_cutoff_cpm=2,
            condition_id=condition,
            support_n=2,
            replicate_n=2,
            support_fraction=1,
        )
        return AtlasCandidate(
            peak=peak,
            summit=(start + end) // 2,
            start=fixed_start,
            end=fixed_end,
            priority=priority,
        )

    strong = candidate("a", 100, 500, "a_strong", 1, 150, 400)
    weak_same_condition = candidate("a", 200, 300, "a_weak", 0.5, 175, 425)
    right = candidate("b", 400, 800, "b_peak", 1, 500, 750)

    boundaries, trimmed, overlapping = _build_variable_boundaries(
        {"ATAC_ATLAS_000001": strong, "ATAC_ATLAS_000002": right},
        [
            ("ATAC_ATLAS_000001", strong),
            ("ATAC_ATLAS_000001", weak_same_condition),
            ("ATAC_ATLAS_000002", right),
        ],
    )

    assert [(boundary.start, boundary.end) for boundary in boundaries] == [
        (100, 450),
        (450, 800),
    ]
    assert [boundary.condition_votes for boundary in boundaries] == [1, 1]
    assert trimmed == 2
    assert overlapping == 0


def test_variable_boundaries_keep_minimum_width_when_neighbors_overlap():
    def candidate(condition, start, end, fixed_start, fixed_end):
        return AtlasCandidate(
            peak=RefinedPeak(
                chrom="chr1",
                start=start,
                end=end,
                name=f"{condition}_peak",
                score=1000,
                strand=".",
                mean_cpm=3,
                maximum_cpm=10,
                selection_cutoff_cpm=2,
                condition_id=condition,
                support_n=2,
                replicate_n=2,
                support_fraction=1,
            ),
            summit=(start + end) // 2,
            start=fixed_start,
            end=fixed_end,
            priority=1,
        )

    left = candidate("a", 100, 150, 25, 275)
    right = candidate("b", 120, 170, 300, 550)
    boundaries, trimmed, overlapping = _build_variable_boundaries(
        {"ATAC_ATLAS_000001": left, "ATAC_ATLAS_000002": right},
        [("ATAC_ATLAS_000001", left), ("ATAC_ATLAS_000002", right)],
    )

    assert [(boundary.start, boundary.end) for boundary in boundaries] == [
        (100, 150),
        (120, 170),
    ]
    assert trimmed == 0
    assert overlapping == 1


def test_dhs_seed_grouping_is_direct_and_does_not_use_fixed_window_overlap():
    def candidate(condition, start, end, summit, priority):
        peak = RefinedPeak(
            chrom="chr1",
            start=start,
            end=end,
            name=f"{condition}_peak",
            score=1000,
            strand=".",
            mean_cpm=3,
            maximum_cpm=10,
            selection_cutoff_cpm=2,
            condition_id=condition,
            support_n=2,
            replicate_n=2,
            support_fraction=1,
        )
        return AtlasCandidate(
            peak=peak,
            summit=summit,
            start=summit - 125,
            end=summit + 125,
            priority=priority,
        )

    seed = candidate("a", 100, 260, 200, 1)
    direct_match = candidate("b", 220, 380, 240, 0.9)
    chain_only = candidate("c", 350, 510, 370, 0.8)
    fixed_window_only = candidate("d", 270, 320, 295, 0.7)

    selected, assignments = _select_dhs_seeds(
        [seed, direct_match, chain_only, fixed_window_only],
        {"chr1": 0},
    )

    assert selected == [seed, chain_only, fixed_window_only]
    assert [(candidate.peak.condition_id, index) for candidate, index in assignments] == [
        ("a", 0),
        ("b", 0),
        ("c", 1),
        ("d", 2),
    ]


def test_narrow_source_first_preserves_two_modes_around_a_broad_bridge():
    def candidate(condition, source_start, source_end, summit, maximum):
        peak = RefinedPeak(
            chrom="chr1",
            start=source_start,
            end=source_end,
            name=f"{condition}_peak",
            score=1000,
            strand=".",
            mean_cpm=3,
            maximum_cpm=maximum,
            selection_cutoff_cpm=2,
            condition_id=condition,
            support_n=2,
            replicate_n=2,
            support_fraction=1,
        )
        return AtlasCandidate(
            peak=peak,
            summit=summit,
            start=summit - 125,
            end=summit + 125,
            priority=1,
        )

    broad = candidate("broad", 100, 500, 300, 20)
    left = candidate("left", 100, 180, 140, 10)
    right = candidate("right", 420, 500, 460, 9)

    selected, assignments = _select_narrow_source_first(
        [broad, left, right], {"chr1": 0}
    )

    assert selected == [left, right]
    assigned = [
        (candidate.peak.condition_id, selected_index)
        for candidate, selected_index in assignments
    ]
    assert assigned == [
        ("left", 0),
        ("right", 1),
        ("broad", 0),
        ("broad", 1),
    ]


def test_signal_shaped_atlas_uses_only_equal_weight_contributors(
    tmp_path: Path, monkeypatch
):
    membership = tmp_path / "membership.tsv"
    membership.write_text(
        "atlas_peak_id\tcondition_id\tsource_peak_id\tsource_chrom\t"
        "source_start\tsource_end\tsource_summit\tatlas_center\t"
        "summit_offset\tsource_support_n\tsource_replicate_n\t"
        "source_support_fraction\tsource_priority_percentile\t"
        "source_atlas_coverage_fraction\n"
        "ATAC_ATLAS_000001\ta\ta1\tchr1\t150\t250\t200\t200\t0\t2\t2\t1\t1\t1\n"
        "ATAC_ATLAS_000002\ta\ta2\tchr1\t650\t750\t700\t700\t0\t2\t2\t1\t1\t1\n"
        "ATAC_ATLAS_000002\tb\tb2\tchr1\t660\t760\t710\t700\t10\t2\t2\t1\t1\t1\n"
    )
    signals = {
        "a.bw": [
            ("chr1", 170, 230, 20.0),
            ("chr1", 300, 340, 1000.0),
            ("chr1", 680, 720, 20.0),
        ],
        "b.bw": [("chr1", 690, 730, 2.0)],
    }
    written: dict[str, object] = {"entries": []}

    class FakeReader:
        def __init__(self, path):
            self.entries = signals[Path(path).name]

        def chroms(self):
            return {"chr1": 1000}

        def intervals(self, chrom, start, end):
            return [
                (interval_start, interval_end, value)
                for entry_chrom, interval_start, interval_end, value in self.entries
                if entry_chrom == chrom and interval_start < end and interval_end > start
            ]

        def close(self):
            pass

    class FakeWriter:
        def __init__(self, path):
            self.path = Path(path)

        def addHeader(self, header):
            written["header"] = header

        def addEntries(self, chroms, starts, *, ends, values):
            written["entries"].extend(zip(chroms, starts, ends, values, strict=True))

        def close(self):
            self.path.touch()

    def fake_open(path, mode=None):
        return FakeWriter(path) if mode == "w" else FakeReader(path)

    monkeypatch.setitem(sys.modules, "pyBigWig", SimpleNamespace(open=fake_open))
    bigwigs = {condition: tmp_path / f"{condition}.bw" for condition in ("a", "b")}
    for path in bigwigs.values():
        path.touch()
    output = tmp_path / "signal-shaped.bed"
    aggregate = tmp_path / "aggregate.bw"
    diagnostics = tmp_path / "signal-shape.tsv"
    stats = tmp_path / "signal-shape.json"

    metrics = build_signal_shaped_atlas(
        membership_tsv=membership,
        condition_bigwigs=bigwigs,
        output_bed=output,
        aggregate_bigwig=aggregate,
        diagnostics_tsv=diagnostics,
        stats_json=stats,
        window_size=400,
        bin_size=10,
        smoothing_bins=3,
        minimum_length=50,
        maximum_length=200,
    )

    peaks = [line.split("\t") for line in output.read_text().splitlines()]
    assert [row[3] for row in peaks] == ["ATAC_ATLAS_000001", "ATAC_ATLAS_000002"]
    assert all(50 <= int(row[2]) - int(row[1]) <= 200 for row in peaks)
    diagnostic_rows = diagnostics.read_text().splitlines()
    assert diagnostic_rows[1].split("\t")[7] == "1"
    assert diagnostic_rows[2].split("\t")[7] == "2"
    assert 150 <= int(diagnostic_rows[1].split("\t")[5]) < 250
    assert metrics["atlas_peaks"] == 2
    assert metrics["contributor_profiles"] == 3
    assert metrics["no_signal_peaks"] == 0
    assert json.loads(stats.read_text()) == metrics
    assert aggregate.is_file()
    assert written["header"] == [("chr1", 1000)]
    assert written["entries"]
    assert all(0 <= entry[3] <= 1 for entry in written["entries"])


def test_dhs_support_fwhm_counts_each_condition_once_per_base(
    tmp_path: Path, monkeypatch
):
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr1\t1000\n")
    anchors = tmp_path / "anchors.bed"
    anchors.write_text(
        "chr1\t100\t350\tATAC_ATLAS_000001\t1000\t.\n"
        "chr1\t400\t650\tATAC_ATLAS_000002\t1000\t.\n"
    )
    condition_dhs = {
        "a": tmp_path / "a.bed",
        "b": tmp_path / "b.bed",
        "c": tmp_path / "c.bed",
    }
    condition_dhs["a"].write_text(
        "chr1\t120\t300\ta1\n"
        "chr1\t180\t320\ta1_overlap\n"
        "chr1\t470\t600\ta2\n"
    )
    condition_dhs["b"].write_text("chr1\t150\t280\tb1\n")
    condition_dhs["c"].write_text("chr1\t160\t290\tc1\n")
    written: dict[str, object] = {"entries": []}

    class FakeWriter:
        def __init__(self, path):
            self.path = Path(path)

        def addHeader(self, header):
            written["header"] = header

        def addEntries(self, chroms, starts, *, ends, values):
            written["entries"].extend(zip(chroms, starts, ends, values, strict=True))

        def close(self):
            self.path.touch()

    monkeypatch.setitem(
        sys.modules,
        "pyBigWig",
        SimpleNamespace(open=lambda path, mode=None: FakeWriter(path)),
    )
    support_bigwig = tmp_path / "support.bw"
    output = tmp_path / "fwhm.bed"
    diagnostics = tmp_path / "fwhm.tsv"
    stats = tmp_path / "fwhm.json"

    metrics = build_dhs_support_fwhm(
        anchors_bed=anchors,
        condition_dhs=condition_dhs,
        chromosome_sizes_path=chrom_sizes,
        support_bigwig=support_bigwig,
        output_bed=output,
        diagnostics_tsv=diagnostics,
        stats_json=stats,
    )

    peaks = [line.split("\t") for line in output.read_text().splitlines()]
    assert [(int(row[1]), int(row[2])) for row in peaks] == [
        (150, 290),
        (470, 600),
    ]
    rows = [line.split("\t") for line in diagnostics.read_text().splitlines()]
    assert rows[1][9:13] == ["3", "3", "1", "2"]
    assert rows[2][9:13] == ["1", "3", "0.333333", "1"]
    assert metrics["atlas_peaks"] == 2
    assert metrics["merged_condition_dhs_intervals"] == 4
    assert metrics["no_support_peaks"] == 0
    assert json.loads(stats.read_text()) == metrics
    assert support_bigwig.is_file()
    assert written["header"] == [("chr1", 1000)]
    assert max(entry[3] for entry in written["entries"]) == 1


def test_center_mode_half_prominence_separates_lower_peak_from_taller_neighbor(
    tmp_path: Path,
):
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr1\t1000\n")
    anchors = tmp_path / "anchors.bed"
    anchors.write_text("chr1\t100\t350\tATAC_ATLAS_000001\t1000\t.\n")
    condition_dhs = {name: tmp_path / f"{name}.bed" for name in "abcde"}
    condition_dhs["a"].write_text("chr1\t180\t320\ta1\n")
    condition_dhs["b"].write_text("chr1\t180\t320\tb1\n")
    condition_dhs["c"].write_text(
        "chr1\t180\t240\tc1\n"
        "chr1\t260\t320\tc2\n"
    )
    condition_dhs["d"].write_text("chr1\t260\t320\td1\n")
    condition_dhs["e"].write_text("")
    output = tmp_path / "center-mode.bed"
    diagnostics = tmp_path / "center-mode.tsv"
    stats = tmp_path / "center-mode.json"

    metrics = build_dhs_center_mode_half_prominence(
        anchors_bed=anchors,
        condition_dhs=condition_dhs,
        chromosome_sizes_path=chrom_sizes,
        output_bed=output,
        diagnostics_tsv=diagnostics,
        stats_json=stats,
    )

    peak = output.read_text().split("\t")
    assert (int(peak[1]), int(peak[2])) == (180, 240)
    row = diagnostics.read_text().splitlines()[1].split("\t")
    assert row[9:16] == ["3", "5", "0.6", "0", "2", "1", "3"]
    assert row[16:19] == ["2", "0", "0"]
    assert metrics["atlas_peaks"] == 1
    assert metrics["multiple_local_mode_peaks"] == 1
    assert metrics["bridged_higher_peak_count"] == 0
    assert json.loads(stats.read_text()) == metrics
