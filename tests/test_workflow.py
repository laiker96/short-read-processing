import copy
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError
from short_read_processing.configuration import ATAC_QPOIS_DEFAULTS, REFERENCE_SOURCES
from short_read_processing.workflow_config import validate_workflow_config


REPO_ROOT = Path(__file__).resolve().parents[1]
BASE_CONFIG = yaml.safe_load((REPO_ROOT / "tests/fixtures/workflow_config.yaml").read_text())


def _chip_callpeak(*, broad=False):
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


def _hmmratac():
    return {
        "command": "hmmratac",
        "lower": 10,
        "upper": 20,
        "prescan_cutoff": 1.2,
    }


def _dry_run(tmp_path: Path, config: dict, name: str = "workflow") -> str:
    config["output_dir"] = str(tmp_path / "results")
    config_path = tmp_path / f"{name}.yaml"
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
            "--printshellcmds",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    return output


@pytest.mark.parametrize(
    "branch",
    ["atac_qpois", "atac_hmmratac", "atac_se", "chip_tf", "chip_histone", "chip_histone_ip_only"],
)
def test_workflow_branches_dry_run(tmp_path, branch):
    config = copy.deepcopy(BASE_CONFIG)
    treatment = config["samples"][0]
    if branch == "atac_hmmratac":
        treatment["peak_caller"] = _hmmratac()
    elif branch == "atac_se":
        treatment["layout"] = "single"
        treatment.pop("r2")
    elif branch.startswith("chip"):
        assay = "chip_histone" if branch == "chip_histone_ip_only" else branch
        config["assay"] = assay
        config.pop("atac_qpois")
        treatment["peak_caller"] = _chip_callpeak(broad=assay == "chip_histone")
        if branch != "chip_histone_ip_only":
            treatment["control"] = "input_rep1"
            control = copy.deepcopy(treatment)
            control["id"] = "input_rep1"
            control["accessions"] = ["SRR123457"]
            control["role"] = "control"
            control.pop("control")
            control.pop("peak_caller")
            config["samples"].append(control)

    output = _dry_run(tmp_path, config, branch)
    if branch in {"atac_qpois", "atac_se"}:
        assert "prepare_atac_tn5_insertions" in output
        assert "call_atac_replicate_qpois" in output
        assert "refine_atac_replicate_qpois" in output
    if branch == "atac_se":
        assert "filter_atac_short_fragments" not in output
    if branch == "atac_hmmratac":
        assert "hmmratac_replicate" in output
        assert ".hmmratac.narrowPeak" in output
        assert "call_atac_replicate_qpois" not in output
    if branch == "chip_histone_ip_only":
        assert "callpeak_broad" in output
        assert "chip_fingerprint" not in output


def test_technical_lanes_align_separately_then_merge(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    sample = config["samples"][0]
    sample["accessions"].append("SRR123457")
    sample["r1"].append(sample["r1"][0])
    sample["r2"].append(sample["r2"][0])

    output = _dry_run(tmp_path, config, "two-lanes")

    assert re.search(r"align_lane\s+2", output)
    assert re.search(r"merge_and_mark_duplicates\s+1", output)


def test_workflow_config_rejects_scaled_atac_qpois_signal():
    config = copy.deepcopy(BASE_CONFIG)
    config["samples"][0]["peak_caller"]["spmr"] = True
    with pytest.raises(AcquisitionError, match="invalid two-ended Tn5 qpois"):
        validate_workflow_config(config)


@pytest.mark.parametrize(
    ("field", "value"),
    [("minimum_exponent", -1), ("maximum_exponent", 1), ("minimum_length", 0)],
)
def test_workflow_config_rejects_invalid_qpois_parameters(field, value):
    config = copy.deepcopy(BASE_CONFIG)
    config["atac_qpois"] = dict(ATAC_QPOIS_DEFAULTS)
    config["atac_qpois"][field] = value
    with pytest.raises(AcquisitionError, match="ATAC qpois parameters are invalid"):
        validate_workflow_config(config)


def _add_second_replicate(config: dict) -> None:
    second = copy.deepcopy(config["samples"][0])
    second["id"] = "atac_rep2"
    second["accessions"] = ["SRR123457"]
    config["samples"].append(second)


def _enable_consensus(config: dict) -> None:
    config["atac_consensus"] = {
        "enabled": True,
        "conditions": [
            {
                "id": "example",
                "label": "example",
                "samples": ["atac_rep1", "atac_rep2"],
            }
        ],
        "minimum_replicates": 2,
        "replicate_overlap_fraction": 0.5,
    }


def test_qpois_condition_consensus_builds_master_dhs_as_final_atac_step(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)

    output = _dry_run(tmp_path, config, "qpois-consensus")

    for rule in (
        "pool_atac_condition_insertions",
        "call_atac_condition_qpois",
        "refine_atac_condition_qpois",
        "atac_condition_pileup_bigwig",
        "atac_condition_qpois_bigwig",
        "filter_atac_qpois_replicate_support",
        "build_atac_master_dhs",
    ):
        assert rule in output
    assert "replicate-supported.bed" in output
    assert "master_dhs.bed" in output


def test_hmmratac_condition_consensus_dry_run(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["samples"][0]["peak_caller"] = _hmmratac()
    _add_second_replicate(config)
    _enable_consensus(config)

    output = _dry_run(tmp_path, config, "hmmratac-consensus")

    assert "merge_atac_condition_bams" in output
    assert "hmmratac_condition" in output
    assert "filter_atac_hmmratac_replicate_support" in output
    assert "build_atac_master_dhs" in output
    assert "call_atac_condition_qpois" not in output


def test_workflow_config_rejects_negative_master_summit_distance():
    config = copy.deepcopy(BASE_CONFIG)
    _add_second_replicate(config)
    _enable_consensus(config)
    config["atac_master"] = {
        "summit_max_distance": -1,
        "minimum_summit_separation": 50,
    }

    with pytest.raises(AcquisitionError, match="master DHS parameters"):
        validate_workflow_config(config)


def test_auto_reference_preparation_dry_run(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    reference_root = tmp_path / "references" / "dm6"
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
    output = _dry_run(tmp_path, config, "auto-reference")
    assert "prepare_reference_fasta" in output


def test_workflow_config_rejects_unchecked_reference_source():
    config = copy.deepcopy(BASE_CONFIG)
    config["reference"]["preparation"] = {
        "mode": "download",
        **REFERENCE_SOURCES["dm6"],
    }
    config["reference"]["preparation"]["fasta"] = {
        "url": "https://example.org/dm6.fa.gz",
        "checksum": "",
    }
    with pytest.raises(AcquisitionError, match="checksum"):
        validate_workflow_config(config)
