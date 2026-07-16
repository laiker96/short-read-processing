"""Create a sorted, unique, zero-based TSS BED from a GTF annotation."""

import gzip
from pathlib import Path
import re


ATTRIBUTE_RE = re.compile(r'(?:^|;\s*)(gene_id|transcript_id)\s+"([^"]+)"')


def attributes(value: str) -> dict[str, str]:
    return {key: item for key, item in ATTRIBUTE_RE.findall(value)}


annotation = Path(snakemake.input.annotation)
fai = Path(snakemake.input.fai)
output = Path(snakemake.output.bed)
log = Path(snakemake.log[0])
contigs = {line.split("\t", 1)[0] for line in fai.read_text().splitlines() if line}
tss_records: set[tuple[str, int, int, str, str]] = set()
opener = gzip.open if annotation.suffix == ".gz" else open
with opener(annotation, "rt", encoding="utf-8") as handle:
    for line in handle:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "transcript" or fields[0] not in contigs:
            continue
        start, end = int(fields[3]), int(fields[4])
        strand = fields[6]
        position = start - 1 if strand == "+" else end - 1
        identifiers = attributes(fields[8])
        name = identifiers.get("transcript_id") or identifiers.get("gene_id") or "TSS"
        tss_records.add((fields[0], position, position + 1, name, strand))

if not tss_records:
    raise ValueError(f"No transcript TSS records found in {annotation}")
ordered = sorted(tss_records, key=lambda item: (item[0], item[1], item[2], item[3], item[4]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
temporary = Path(str(output) + ".tmp")
with temporary.open("w", encoding="utf-8") as handle:
    for chrom, start, end, name, strand in ordered:
        handle.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
temporary.replace(output)
log.write_text(f"Wrote {len(ordered)} transcript TSS records from {annotation}\n")
