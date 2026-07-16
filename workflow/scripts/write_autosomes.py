"""Write configured autosomes after checking that they exist in the FASTA index."""

from pathlib import Path


fai = Path(snakemake.input.fai)
output = Path(snakemake.output.contigs)
log = Path(snakemake.log[0])
available = {line.split("\t", 1)[0] for line in fai.read_text().splitlines() if line}
autosomes = list(snakemake.params.autosomes)
missing = [contig for contig in autosomes if contig not in available]
if missing:
    raise ValueError("Configured autosomes missing from FASTA: " + ", ".join(missing))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
temporary = Path(str(output) + ".tmp")
temporary.write_text("".join(f"{contig}\n" for contig in autosomes), encoding="utf-8")
temporary.replace(output)
log.write_text(f"Wrote {len(autosomes)} autosomes from {fai}\n", encoding="utf-8")
