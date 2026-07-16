"""Combine stable per-sample alignment, peak, and FRiP metrics."""

import csv
import json
import re
from pathlib import Path


manifest = json.loads(str(snakemake.params.manifest))
rows = []
for sample, paths in manifest.items():
    flagstat_text = Path(paths["flagstat"]).read_text(encoding="utf-8")

    def count(pattern: str):
        match = re.search(pattern, flagstat_text, flags=re.MULTILINE)
        return int(match.group(1)) if match else None

    row = {
        "sample": sample,
        "assay": paths["assay"],
        "layout": paths["layout"],
        "role": paths["role"],
        "total_alignments": count(r"^(\d+) \+ \d+ in total"),
        "mapped_alignments": count(r"^(\d+) \+ \d+ mapped \("),
        "properly_paired": count(r"^(\d+) \+ \d+ properly paired \("),
        "frip_numerator": None,
        "frip_denominator": None,
        "frip": None,
        "peak_count": None,
    }
    if paths["frip"]:
        frip = json.loads(Path(paths["frip"]).read_text(encoding="utf-8"))
        row["frip_numerator"] = frip["numerator"]
        row["frip_denominator"] = frip["denominator"]
        row["frip"] = frip["frip"]
    if paths["peaks"]:
        row["peak_count"] = sum(
            1
            for line in Path(paths["peaks"]).open(encoding="utf-8")
            if line.strip() and not line.startswith(("#", "track", "browser"))
        )
    rows.append(row)

fields = list(rows[0])
output_tsv = Path(str(snakemake.output.tsv))
output_json = Path(str(snakemake.output.json))
output_tsv.parent.mkdir(parents=True, exist_ok=True)
with output_tsv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
output_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
log_path = Path(str(snakemake.log[0]))
log_path.parent.mkdir(parents=True, exist_ok=True)
log_path.write_text(f"Aggregated metrics for {len(rows)} sample(s)\n", encoding="utf-8")
