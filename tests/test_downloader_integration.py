import hashlib
import shutil
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

import pytest

from short_read_processing.accessions import FilePlan, RunPlan
from short_read_processing.downloader import DownloadOptions, download_ena


class _QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        return


@pytest.mark.skipif(shutil.which("aria2c") is None, reason="aria2c is not installed")
def test_real_aria2_download_and_checksum(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source = source_dir / "SRR123.fastq.gz"
    payload = b"deterministic-fastq-block\n" * 100_000
    source.write_bytes(payload)
    checksum = hashlib.md5(payload).hexdigest()

    handler = partial(_QuietHandler, directory=str(source_dir))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        destination = tmp_path / "output" / "SRR123" / source.name
        plan = RunPlan(
            requested_accession="SRR123",
            experiment_accession="SRX123",
            run_accession="SRR123",
            library_layout="SINGLE",
            backend="ena",
            run_dir=destination.parent,
            files=[
                FilePlan(
                    url=f"http://127.0.0.1:{server.server_port}/{source.name}",
                    md5=checksum,
                    size_bytes=len(payload),
                    path=destination,
                    mate="r1",
                )
            ],
        )
        download_ena(
            [plan],
            DownloadOptions(file_jobs=2, connections=4, sra_jobs=1, threads=2),
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    assert destination.read_bytes() == payload
    assert plan.status == "downloaded"
