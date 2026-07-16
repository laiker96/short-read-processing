from short_read_processing.accessions import FilePlan, RunPlan
from short_read_processing.manifest import read_manifest, write_manifest


def test_dry_run_does_not_downgrade_completed_manifest(tmp_path):
    fastq = tmp_path / "SRR123" / "SRR123.fastq.gz"
    fastq.parent.mkdir()
    fastq.write_bytes(b"fastq")
    completed = RunPlan(
        requested_accession="SRR123",
        experiment_accession="SRX123",
        run_accession="SRR123",
        library_layout="SINGLE",
        backend="ena",
        run_dir=fastq.parent,
        files=[FilePlan("https://example/fastq", "abc", 5, fastq, "r1")],
        status="downloaded",
    )
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, [completed])

    planned = RunPlan(
        requested_accession="SRR123",
        experiment_accession="SRX123",
        run_accession="SRR123",
        library_layout="SINGLE",
        backend="ena",
        run_dir=fastq.parent,
        files=[FilePlan("https://example/fastq", "abc", 5, fastq, "r1")],
        status="planned",
    )
    write_manifest(manifest, [planned])
    assert read_manifest(manifest)[0]["status"] == "downloaded"
