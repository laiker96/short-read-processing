"""Calculate fragment-level FRiP for paired reads and read-level FRiP for single reads."""

import json
import subprocess
import tempfile
from pathlib import Path


bam = Path(str(snakemake.input.bam))
peaks = Path(str(snakemake.input.peaks))
output_tsv = Path(str(snakemake.output.tsv))
output_json = Path(str(snakemake.output.json))
log_path = Path(str(snakemake.log[0]))
layout = str(snakemake.params.layout)
output_tsv.parent.mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(dir=output_tsv.parent) as temporary:
    fragments = Path(temporary) / "fragments.bed"
    if layout == "paired":
        collate = subprocess.Popen(
            ["samtools", "collate", "-@", str(snakemake.threads), "-O", str(bam)],
            stdout=subprocess.PIPE,
        )
        bed_process = subprocess.Popen(
            ["bedtools", "bamtobed", "-bedpe", "-i", "stdin"],
            stdin=collate.stdout,
            stdout=subprocess.PIPE,
            text=True,
        )
        assert collate.stdout is not None
        collate.stdout.close()
        assert bed_process.stdout is not None
        denominator = 0
        with fragments.open("w", encoding="utf-8") as handle:
            for line in bed_process.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[0] != fields[3]:
                    continue
                start = min(int(fields[1]), int(fields[4]))
                end = max(int(fields[2]), int(fields[5]))
                if end <= start:
                    continue
                handle.write(f"{fields[0]}\t{start}\t{end}\n")
                denominator += 1
        if bed_process.wait() != 0 or collate.wait() != 0:
            raise RuntimeError("Failed to convert paired BAM to fragment BED")
        unit = "fragments"
    else:
        with fragments.open("w", encoding="utf-8") as handle:
            subprocess.run(
                ["bedtools", "bamtobed", "-i", str(bam)],
                stdout=handle,
                check=True,
                text=True,
            )
        denominator = sum(1 for line in fragments.open(encoding="utf-8") if line.strip())
        unit = "reads"

    overlap = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(fragments), "-b", str(peaks)],
        capture_output=True,
        text=True,
        check=True,
    )
    numerator = sum(1 for line in overlap.stdout.splitlines() if line.strip())

ratio = numerator / denominator if denominator else 0.0
metrics = {
    "sample": str(snakemake.wildcards.sample),
    "unit": unit,
    "numerator": numerator,
    "denominator": denominator,
    "frip": ratio,
}
output_tsv.write_text(
    "sample\tunit\tnumerator\tdenominator\tfrip\n"
    f"{metrics['sample']}\t{unit}\t{numerator}\t{denominator}\t{ratio:.8f}\n",
    encoding="utf-8",
)
output_json.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
log_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text(
    f"Calculated {unit}-level FRiP: {numerator}/{denominator}={ratio:.8f}\n",
    encoding="utf-8",
)
