from pathlib import Path

from short_read_processing.accessions import FilePlan, RunPlan
from short_read_processing.downloader import (
    _discard_untracked_partial_files,
    _download_one_sra,
)


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


def test_sra_download_stages_fastqs_before_completion(tmp_path, monkeypatch):
    run_dir = tmp_path / "raw" / "SRR123456"
    run_dir.mkdir(parents=True)
    (run_dir / "SRR123456_1.fastq.gz").write_bytes(b"partial")
    plan = RunPlan(
        requested_accession="SRR123456",
        experiment_accession="SRX999999",
        run_accession="SRR123456",
        library_layout="PAIRED",
        backend="sra",
        run_dir=run_dir,
    )

    monkeypatch.setattr(
        "short_read_processing.downloader._require_executable", lambda name: name
    )

    def fake_run(command, *, label):
        if label.startswith("fasterq-dump"):
            output = Path(command[command.index("--outdir") + 1])
            (output / "SRR123456_1.fastq").write_bytes(b"r1")
            (output / "SRR123456_2.fastq").write_bytes(b"r2")

    def fake_gzip(path, *, threads):
        compressed = path.with_suffix(path.suffix + ".gz")
        compressed.write_bytes(path.read_bytes())
        path.unlink()
        return compressed

    monkeypatch.setattr("short_read_processing.downloader._run", fake_run)
    monkeypatch.setattr("short_read_processing.downloader._gzip_fastq", fake_gzip)

    _download_one_sra(plan, threads=2, keep_cache=False)

    assert (run_dir / ".download-complete").read_text() == "complete\n"
    assert (run_dir / "SRR123456_1.fastq.gz").read_bytes() == b"r1"
    assert (run_dir / "SRR123456_2.fastq.gz").read_bytes() == b"r2"
    assert [item.mate for item in plan.files] == ["r1", "r2"]
