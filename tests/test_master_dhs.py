import json
import math
from pathlib import Path

from short_read_processing import master_dhs
from short_read_processing.master_dhs import SourcePeak, assemble_master_dhs


def source(
    context: str,
    start: int,
    end: int,
    summit: int,
    name: str,
) -> SourcePeak:
    return SourcePeak(
        context=context,
        chrom="chr1",
        start=start,
        end=end,
        name=name,
        summit=summit,
        summit_signal=10.0,
        support_n=2,
        replicate_n=2,
        support_fraction=1.0,
        peak_method="qpois",
    )


def test_summit_uses_center_of_maximum_plateau_and_midpoint_fallback():
    peak = source("A", 100, 106, 102, "A1")
    peak = SourcePeak(**{**peak.__dict__, "summit": None})

    measured = master_dhs._summit_from_values(peak, [0, 2, 2, 2, 0, 0])
    fallback = master_dhs._summit_from_values(peak, [math.nan] * 6)

    assert measured.summit == 102
    assert measured.summit_signal == 2
    assert measured.summit_fallback is False
    assert fallback.summit == 102
    assert fallback.summit_fallback is True


def test_reader_disambiguates_names_and_records_contig_clipping(tmp_path: Path):
    bed = tmp_path / "context.bed"
    bed.write_text(
        "chr1\t10\t20\trepeated\t0\t.\n"
        "chr1\t90\t105\trepeated\t0\t.\n"
    )

    peaks = master_dhs.read_context_peaks(
        bed,
        context="A",
        chrom_sizes={"chr1": 100},
    )

    assert [peak.name for peak in peaks] == ["repeated", "repeated__duplicate_2"]
    assert peaks[1].end == 100
    assert peaks[1].input_end == 105
    assert peaks[1].coordinate_clipped is True


def test_narrow_same_context_peaks_anchor_two_dhs_despite_broad_bridge():
    peaks = [
        source("resolved", 100, 160, 130, "left"),
        source("resolved", 220, 280, 250, "right"),
        source("broad", 110, 270, 150, "bridge"),
    ]

    masters, clipped, close_merges = assemble_master_dhs(
        peaks,
        chrom_sizes={"chr1": 1000},
        summit_max_distance=150,
    )

    assert len(masters) == 2
    assert [master.summit for master in masters] == [130, 250]
    assert {peak.name for peak in masters[0].members} == {"left", "bridge"}
    assert {peak.name for peak in masters[1].members} == {"right"}
    assert masters[0].end <= masters[1].start
    assert clipped == 1
    assert close_merges == 0


def test_complete_linkage_prevents_transitive_summit_chaining():
    peaks = [
        source("A", 0, 201, 100, "A1"),
        source("B", 50, 301, 200, "B1"),
        source("C", 150, 401, 300, "C1"),
    ]

    masters, _clipped, _close_merges = assemble_master_dhs(
        peaks,
        chrom_sizes={"chr1": 1000},
        summit_max_distance=150,
    )

    assert len(masters) == 2
    assert max(int(peak.summit) for peak in masters[0].members) - min(
        int(peak.summit) for peak in masters[0].members
    ) <= 150


def test_overlap_without_reciprocal_summit_containment_stays_separate():
    peaks = [
        source("narrow", 100, 160, 130, "left"),
        source("broad", 120, 300, 250, "right-shifted"),
    ]

    masters, _clipped, _close_merges = assemble_master_dhs(
        peaks,
        chrom_sizes={"chr1": 1000},
        summit_max_distance=150,
    )

    assert len(masters) == 2


def test_close_shifted_peaks_merge_when_contexts_are_disjoint():
    peaks = [
        source("EAD", 221, 274, 222, "ead"),
        source("LB", 232, 315, 246, "lb"),
    ]

    masters, clipped, close_merges = assemble_master_dhs(
        peaks,
        chrom_sizes={"chr1": 1000},
        summit_max_distance=150,
        minimum_summit_separation=50,
    )

    assert len(masters) == 1
    assert (masters[0].start, masters[0].end) == (221, 315)
    assert clipped == 0
    assert close_merges == 1


def test_close_peaks_remain_separate_when_one_context_resolves_both():
    peaks = [
        source("EAD", 100, 150, 149, "ead_left"),
        source("EAD", 150, 200, 150, "ead_right"),
    ]

    masters, clipped, close_merges = assemble_master_dhs(
        peaks,
        chrom_sizes={"chr1": 1000},
        summit_max_distance=150,
        minimum_summit_separation=50,
    )

    assert len(masters) == 2
    assert clipped == 0
    assert close_merges == 0


def test_registry_outputs_membership_and_context_matrix(tmp_path: Path, monkeypatch):
    chrom_sizes = tmp_path / "chrom.sizes"
    chrom_sizes.write_text("chr1\t1000\n")
    by_context = {
        "A": [source("A", 100, 180, 140, "A1")],
        "B": [source("B", 110, 190, 145, "B1")],
    }

    def fake_load(context, _peaks, _signal, _sizes):
        return by_context[context]

    monkeypatch.setattr(master_dhs, "load_context_peaks", fake_load)
    outputs = {
        "output_bed": tmp_path / "master.bed",
        "summit_bed": tmp_path / "summits.bed",
        "membership_tsv": tmp_path / "membership.tsv",
        "context_matrix_tsv": tmp_path / "matrix.tsv",
        "stats_json": tmp_path / "stats.json",
    }
    metrics = master_dhs.build_master_registry(
        context_peaks={"A": Path("A.bed"), "B": Path("B.bed")},
        context_signals={"A": Path("A.bw"), "B": Path("B.bw")},
        chrom_sizes_path=chrom_sizes,
        **outputs,
    )

    bed_fields = outputs["output_bed"].read_text().rstrip("\n").split("\t")
    assert bed_fields == ["chr1", "100", "190", "DHS0000001", "0", "."]
    assert len(outputs["membership_tsv"].read_text().splitlines()) == 3
    assert outputs["context_matrix_tsv"].read_text().splitlines()[1].endswith("\t1\t1")
    assert json.loads(outputs["stats_json"].read_text()) == metrics
    assert metrics["master_dhs_count"] == 1
    assert metrics["minimum_summit_separation"] == 50
