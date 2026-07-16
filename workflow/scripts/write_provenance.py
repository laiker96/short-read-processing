"""Persist the fully resolved workflow configuration in the result namespace."""

import json
from pathlib import Path


output = Path(str(snakemake.output.config))
log = Path(str(snakemake.log[0]))
output.parent.mkdir(parents=True, exist_ok=True)
log.parent.mkdir(parents=True, exist_ok=True)
resolved = json.loads(str(snakemake.params.config))
output.write_text(json.dumps(resolved, indent=2, sort_keys=True) + "\n", encoding="utf-8")
log.write_text("Wrote fully resolved workflow configuration\n", encoding="utf-8")
