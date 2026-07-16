from pathlib import Path

import pytest

from short_read_processing.accessions import AcquisitionError
from short_read_processing.cli import read_accession_column


def test_read_configurable_accession_column(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("run_id\tlabel\nSRR123\ta\nERX456\tb\nSRR123\tduplicate\n")
    assert read_accession_column(metadata, "run_id") == ["SRR123", "ERX456"]


def test_read_accession_column_from_csv(tmp_path):
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("run_id,label\nSRR123,a\nERR456,b\n")
    assert read_accession_column(metadata, "run_id") == ["SRR123", "ERR456"]


def test_missing_accession_column_is_clear(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("wrong\nSRR123\n")
    with pytest.raises(AcquisitionError, match="available columns: wrong"):
        read_accession_column(metadata, "run_id")
