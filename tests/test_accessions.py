from pathlib import Path

import pytest

from short_read_processing.accessions import (
    AcquisitionError,
    normalize_accession,
    parse_ena_report,
)
from short_read_processing.downloader import _aria2_input


FIXTURES = Path(__file__).parent / "fixtures"


def test_normalize_supported_accessions():
    assert normalize_accession(" srr123 ") == "SRR123"
    assert normalize_accession("ERX42") == "ERX42"


@pytest.mark.parametrize("value", ["SRP123", "SRR", "SRX12x", "", "../SRR123"])
def test_reject_unsupported_accessions(value):
    with pytest.raises(AcquisitionError):
        normalize_accession(value)


def test_parse_paired_ena_report(tmp_path):
    plans = parse_ena_report(
        (FIXTURES / "ena_report_paired.tsv").read_text(),
        requested_accession="SRX123456",
        output_dir=tmp_path,
    )
    assert len(plans) == 1
    plan = plans[0]
    assert plan.backend == "ena"
    assert plan.library_layout == "PAIRED"
    assert [item.mate for item in plan.files] == ["r1", "r2"]
    assert plan.files[0].url.startswith("https://ftp.sra.ebi.ac.uk/")
    assert plan.files[0].path == tmp_path / "SRR123456" / "SRR123456_1.fastq.gz"


def test_aria2_input_has_parallel_download_metadata(tmp_path):
    plan = parse_ena_report(
        (FIXTURES / "ena_report_paired.tsv").read_text(),
        requested_accession="SRX123456",
        output_dir=tmp_path,
    )[0]
    text = _aria2_input(plan.files)
    assert "  out=SRR123456_1.fastq.gz" in text
    assert "  checksum=md5=11111111111111111111111111111111" in text
    assert f"  dir={tmp_path / 'SRR123456'}" in text

