import copy
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError
from short_read_processing.configuration import REFERENCE_SOURCES
from short_read_processing.workflow_config import validate_workflow_config


REPO_ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = yaml.safe_load((REPO_ROOT / "tests/fixtures/workflow_config.yaml").read_text())


def _callpeak(*, broad=False):
    config = {
        "command": "callpeak",
        "format": "BAMPE",
        "qvalue": 0.01,
        "broad": broad,
        "nomodel": False,
        "shift": None,
        "extsize": None,
        "write_bedgraph": True,
        "spmr": True,
    }
    if broad:
        config["broad_cutoff"] = 0.1
    return config


@pytest.mark.parametrize("branch", ["atac_hmmratac", "atac_se", "chip_tf", "chip_histone"])
def test_workflow_branches_dry_run(tmp_path, branch):
    config = copy.deepcopy(BASE_CONFIG)
    config["output_dir"] = str(tmp_path / "results")
    treatment = config["samples"][0]
    if branch == "atac_hmmratac":
        treatment["peak_caller"] = {
            "command": "hmmratac",
            "lower": 10,
            "upper": 20,
            "prescan_cutoff": 1.2,
        }
    elif branch == "atac_se":
        treatment["layout"] = "single"
        treatment.pop("r2")
        treatment["peak_caller"] = _callpeak()
        treatment["peak_caller"]["format"] = "BAM"
    else:
        config["assay"] = branch
        treatment["peak_caller"] = _callpeak(broad=branch == "chip_histone")
        treatment["control"] = "input_rep1"
        control = copy.deepcopy(treatment)
        control["id"] = "input_rep1"
        control["accessions"] = ["SRR123457"]
        control["role"] = "control"
        control.pop("control")
        control.pop("peak_caller")
        config["samples"].append(control)

    config_path = tmp_path / f"{branch}.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    snakemake = Path(sys.executable).with_name("snakemake")
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    result = subprocess.run(
        [
            str(snakemake),
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "2",
            "--dry-run",
            "--quiet",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_technical_lanes_align_separately_then_merge(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["output_dir"] = str(tmp_path / "results")
    sample = config["samples"][0]
    sample["accessions"].append("SRR123457")
    sample["r1"].append(sample["r1"][0])
    sample["r2"].append(sample["r2"][0])

    config_path = tmp_path / "two-lanes.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    snakemake = Path(sys.executable).with_name("snakemake")
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    result = subprocess.run(
        [
            str(snakemake),
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "8",
            "--dry-run",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert re.search(r"align_lane\s+2", output)
    assert re.search(r"merge_and_mark_duplicates\s+1", output)


def test_workflow_config_rejects_callpeak_without_bedgraphs():
    config = copy.deepcopy(BASE_CONFIG)
    config["samples"][0]["peak_caller"]["write_bedgraph"] = False
    with pytest.raises(AcquisitionError, match="must write -B --SPMR"):
        validate_workflow_config(config)


def test_auto_reference_preparation_dry_run(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    reference_root = tmp_path / "references" / "dm6"
    config["output_dir"] = str(tmp_path / "results")
    config["reference"].update(
        {
            "fasta": str(reference_root / "dm6.fa"),
            "bowtie2_index": str(reference_root / "bowtie2" / "dm6"),
            "chrom_sizes": str(reference_root / "dm6.chrom.sizes"),
            "blacklist_bed": str(reference_root / "dm6.blacklist.bed"),
            "tss_bed": str(reference_root / "dm6.tss.bed"),
            "autosomes_file": str(reference_root / "dm6.autosomes.txt"),
            "preparation": {"mode": "download", **REFERENCE_SOURCES["dm6"]},
        }
    )
    config_path = tmp_path / "auto-reference.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    snakemake = Path(sys.executable).with_name("snakemake")
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    result = subprocess.run(
        [
            str(snakemake),
            "--snakefile",
            "workflow/Snakefile",
            "--configfile",
            str(config_path),
            "--cores",
            "2",
            "--dry-run",
            "--quiet",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_workflow_config_rejects_unchecked_reference_source():
    config = copy.deepcopy(BASE_CONFIG)
    config["reference"]["preparation"] = {"mode": "download", **REFERENCE_SOURCES["dm6"]}
    config["reference"]["preparation"]["fasta"] = {
        "url": "https://example.org/dm6.fa.gz",
        "checksum": "",
    }
    with pytest.raises(AcquisitionError, match="checksum"):
        validate_workflow_config(config)
