from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.sample_sheet import read_sample_sheet


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
