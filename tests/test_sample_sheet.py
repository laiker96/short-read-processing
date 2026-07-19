from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.sample_sheet import read_sample_sheet
from write_atlas_atac_sample_sheet import selected_sample_rows as selected_atac_rows
from write_atlas_h3k27ac_sample_sheet import selected_sample_rows as selected_h3k27ac_rows


def test_csv_is_accepted_and_defaults_are_resolved(tmp_path: Path):
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "accession,sample_id,assay,genome,role,control_id,replicate,peak_caller\n"
        "SRR123,atac_rep1,atac,dm6,treatment,,1,\n"
    )
    row = read_sample_sheet(sheet)[0]
    assert row["accession"] == "SRR123"
    assert row["peak_caller"] == "hmmratac"
    assert row["adapter_preset"] == "nextera"
    assert row["mapq_minimum"] == 30


def test_missing_required_column_fails(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("accession\tsample_id\nSRR123\tatac_rep1\n")
    with pytest.raises(AcquisitionError, match="missing required columns"):
        read_sample_sheet(sheet)


def test_technical_runs_must_have_identical_parameters(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tsample_id\tassay\tgenome\trole\tcontrol_id\treplicate\tpeak_caller\n"
        "SRR123\tatac_rep1\tatac\tdm6\ttreatment\t\t1\thmmratac\n"
        "SRR124\tatac_rep1\tatac\tdm6\ttreatment\t\t2\thmmratac\n"
    )
    with pytest.raises(AcquisitionError, match="disagree on: replicate"):
        read_sample_sheet(sheet)


def test_callpeak_shift_requires_bam_format(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tsample_id\tassay\tgenome\trole\tcontrol_id\treplicate\tpeak_caller"
        "\tmacs3_shift\tmacs3_extsize\n"
        "SRR123\tatac_rep1\tatac\tdm6\ttreatment\t\t1\tcallpeak\t-75\t150\n"
    )
    with pytest.raises(AcquisitionError, match="requires macs3_format=BAM"):
        read_sample_sheet(sheet)


def test_ip_only_histone_defaults_are_valid(tmp_path: Path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "accession\tsample_id\tassay\tgenome\trole\tcontrol_id\treplicate\tpeak_caller\n"
        "SRR123\th3k27ac_rep1\tchip_histone\tdm6\ttreatment\t\t1\tcallpeak\n"
    )
    row = read_sample_sheet(sheet)[0]

    assert row["control_id"] is None
    assert row["macs3_broad"] is True
    assert row["macs3_broad_cutoff"] == 0.1
    assert row["adapter_preset"] == "truseq"


def test_atlas_selected_rows_preserve_biological_and_technical_replicates():
    metadata = Path(__file__).parents[1] / "resources/atlas_atac_seq_metadata.tsv"
    rows = selected_atac_rows(metadata)

    assert len(rows) == 23
    assert len({row["accession"] for row in rows}) == 23
    assert len({row["sample_id"] for row in rows}) == 22
    assert [row["sample_id"] for row in rows].count("e5_atac_rep1") == 2
    assert {row["sample_id"] for row in rows if row["sample_id"].startswith("e11_")} == {
        "e11_atac_rep1",
        "e11_atac_rep2",
    }
    assert {row["sample_id"] for row in rows if row["sample_id"].startswith("wid_")} == {
        "wid_atac_rep1",
        "wid_atac_rep2",
    }


def test_atlas_h3k27ac_selected_rows_are_ip_only():
    metadata = Path(__file__).parents[1] / "resources/atlas_h3k27ac_metadata_ip_only.tsv"
    rows = selected_h3k27ac_rows(metadata)

    assert len(rows) == 15
    assert len({row["accession"] for row in rows}) == 15
    assert len({row["sample_id"] for row in rows}) == 15
    assert all(row["assay"] == "chip_histone" for row in rows)
    assert all(row["control_id"] == "" for row in rows)
    assert {row["sample_id"] for row in rows if str(row["sample_id"]).startswith("e5_")} == {
        "e5_h3k27ac_rep1",
        "e5_h3k27ac_rep2",
    }
