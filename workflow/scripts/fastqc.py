"""Run FastQC in an isolated directory and move its deterministic outputs."""

import shutil
import subprocess
import tempfile
from pathlib import Path


input_path = Path(str(snakemake.input[0]))
html_path = Path(str(snakemake.output.html))
zip_path = Path(str(snakemake.output.zip))
log_path = Path(str(snakemake.log[0]))
html_path.parent.mkdir(parents=True, exist_ok=True)
log_path.parent.mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(dir=html_path.parent) as temporary:
    with log_path.open("w", encoding="utf-8") as log_handle:
        subprocess.run(
            ["fastqc", "--threads", str(snakemake.threads), "--outdir", temporary, str(input_path)],
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    generated_html = list(Path(temporary).glob("*_fastqc.html"))
    generated_zip = list(Path(temporary).glob("*_fastqc.zip"))
    if len(generated_html) != 1 or len(generated_zip) != 1:
        raise RuntimeError(f"FastQC did not produce exactly one HTML/ZIP pair for {input_path}")
    shutil.move(str(generated_html[0]), html_path)
    shutil.move(str(generated_zip[0]), zip_path)
