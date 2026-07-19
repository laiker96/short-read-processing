from pathlib import Path

from short_read_processing.accessions import FilePlan
from short_read_processing.downloader import _discard_untracked_partial_files


def _file_plan(path: Path, size: int) -> FilePlan:
    return FilePlan(url="https://example.org/reads.fastq.gz", md5="", size_bytes=size, path=path)


def test_discards_size_mismatched_file_without_aria2_state(tmp_path):
    partial = tmp_path / "reads.fastq.gz"
    partial.write_bytes(b"partial")

    discarded = _discard_untracked_partial_files([_file_plan(partial, 100)])

    assert discarded == [partial]
    assert not partial.exists()


def test_keeps_aria2_managed_partial_file(tmp_path):
    partial = tmp_path / "reads.fastq.gz"
    partial.write_bytes(b"partial")
    partial.with_name(partial.name + ".aria2").touch()

    discarded = _discard_untracked_partial_files([_file_plan(partial, 100)])

    assert discarded == []
    assert partial.exists()
