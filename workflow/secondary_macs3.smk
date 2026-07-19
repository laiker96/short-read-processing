"""Run MACS3 callpeak from BAMs produced by an existing primary workflow run."""

import json
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(workflow.basedir).parent.resolve()
sys.path.insert(0, str(REPO_ROOT / "src"))

from short_read_processing.macs3 import callpeak_arguments
from short_read_processing.workflow_config import (
    SAFE_ID_RE,
    resolve_input_paths,
    validate_workflow_config,
    wildcard_regex,
)


if not config:
    raise ValueError("Pass the primary resolved YAML with --configfile")
validate_workflow_config(config)
resolve_input_paths(config, Path.cwd())
if config["assay"] != "atac":
    raise ValueError("The secondary shifted MACS3 workflow currently requires assay=atac")

PROJECT = str(config["project"])
SOURCE_RUN_ID = str(config["run_id"])
OUTPUT_RUN_ID = str(
    config.get("secondary_macs3_run_id", "macs3-shift-neg75-extsize-150")
)
if not SAFE_ID_RE.fullmatch(OUTPUT_RUN_ID):
    raise ValueError(f"Invalid secondary_macs3_run_id: {OUTPUT_RUN_ID!r}")

SHIFT = int(config.get("secondary_macs3_shift", -75))
EXTSIZE = int(config.get("secondary_macs3_extsize", 150))
QVALUE = float(config.get("secondary_macs3_qvalue", 0.01))
if EXTSIZE < 1:
    raise ValueError("secondary_macs3_extsize must be positive")
if not 0 < QVALUE <= 1:
    raise ValueError("secondary_macs3_qvalue must be in (0, 1]")

RESULTS_ROOT = Path(str(config["output_dir"]))
SOURCE_ROOT = RESULTS_ROOT / PROJECT / SOURCE_RUN_ID
OUTPUT_ROOT = RESULTS_ROOT / PROJECT / OUTPUT_RUN_ID
SAMPLES = {
    str(sample["id"]): sample
    for sample in config["samples"]
    if sample["role"] == "treatment"
}
SAMPLE_IDS = list(SAMPLES)
SAMPLE_RE = wildcard_regex(SAMPLE_IDS)
SOURCE_BAMS = {
    sample: str(SOURCE_ROOT / "bam" / f"{sample}.final.bam") for sample in SAMPLE_IDS
}
PEAKS = {
    sample: str(OUTPUT_ROOT / "peaks" / sample / f"{sample}_peaks.narrowPeak")
    for sample in SAMPLE_IDS
}
BEDGRAPHS = [
    str(OUTPUT_ROOT / "peaks" / sample / f"{sample}_{suffix}.bdg")
    for sample in SAMPLE_IDS
    for suffix in ("treat_pileup", "control_lambda")
]
PROVENANCE = str(OUTPUT_ROOT / "provenance" / "secondary_macs3.json")
PROVENANCE_PAYLOAD = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "project": PROJECT,
    "source_run_id": SOURCE_RUN_ID,
    "output_run_id": OUTPUT_RUN_ID,
    "source_bams": SOURCE_BAMS,
    "macs3": {
        "command": "callpeak",
        "format": "BAM",
        "qvalue": QVALUE,
        "nomodel": True,
        "shift": SHIFT,
        "extsize": EXTSIZE,
        "write_bedgraph": True,
        "spmr": True,
    },
}


def secondary_callpeak_argv(wildcards):
    peak_config = {
        "command": "callpeak",
        "format": "BAM",
        "qvalue": QVALUE,
        "broad": False,
        "nomodel": True,
        "shift": SHIFT,
        "extsize": EXTSIZE,
        "write_bedgraph": True,
        "spmr": True,
    }
    arguments = callpeak_arguments(
        peak_config,
        treatment_bam=Path(SOURCE_BAMS[wildcards.sample]),
        control_bam=None,
        name=wildcards.sample,
        genome_size=config["reference"]["macs3_genome_size"],
        output_dir=OUTPUT_ROOT / "peaks" / wildcards.sample,
    )
    return shlex.join(arguments)


rule all:
    input:
        list(PEAKS.values()) + BEDGRAPHS + [PROVENANCE]


rule secondary_callpeak:
    input:
        bam=lambda wc: SOURCE_BAMS[wc.sample],
        bai=lambda wc: SOURCE_BAMS[wc.sample] + ".bai"
    output:
        peaks=str(OUTPUT_ROOT / "peaks" / "{sample}" / "{sample}_peaks.narrowPeak"),
        summits=str(OUTPUT_ROOT / "peaks" / "{sample}" / "{sample}_summits.bed"),
        xls=str(OUTPUT_ROOT / "peaks" / "{sample}" / "{sample}_peaks.xls"),
        treat_bdg=str(OUTPUT_ROOT / "peaks" / "{sample}" / "{sample}_treat_pileup.bdg"),
        control_bdg=str(
            OUTPUT_ROOT / "peaks" / "{sample}" / "{sample}_control_lambda.bdg"
        )
    params:
        command=secondary_callpeak_argv
    wildcard_constraints:
        sample=SAMPLE_RE
    conda:
        "envs/peaks.yaml"
    resources:
        mem_mb=8000
    log:
        str(OUTPUT_ROOT / "logs" / "peaks" / "{sample}.callpeak.log")
    shell:
        "mkdir -p $(dirname {output.peaks:q}) $(dirname {log:q}) && "
        "{params.command} > {log:q} 2>&1"


rule secondary_macs3_provenance:
    input:
        peaks=list(PEAKS.values()),
        bedgraphs=BEDGRAPHS
    output:
        config=PROVENANCE
    params:
        config=json.dumps(PROVENANCE_PAYLOAD)
    conda:
        "envs/reporting.yaml"
    log:
        str(OUTPUT_ROOT / "logs" / "provenance" / "secondary_macs3.log")
    script:
        "scripts/write_provenance.py"
