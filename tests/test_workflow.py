import copy
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError
from short_read_processing.configuration import ATAC_REFINEMENT_DEFAULTS, REFERENCE_SOURCES
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


@pytest.mark.parametrize(
    "branch",
    ["atac_hmmratac", "atac_se", "chip_tf", "chip_histone", "chip_histone_ip_only"],
)
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
        assay = "chip_histone" if branch == "chip_histone_ip_only" else branch
        config["assay"] = assay
        treatment["peak_caller"] = _callpeak(broad=assay == "chip_histone")
        if branch != "chip_histone_ip_only":
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
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    if branch == "atac_hmmratac":
        assert "atac_rep1_accessible_regions.narrowPeak" in output
    if branch == "chip_histone_ip_only":
        assert "callpeak_broad" in output
        assert "chip_fingerprint" not in output


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


def test_workflow_config_rejects_atlas_when_refinement_is_disabled():
    config = copy.deepcopy(BASE_CONFIG)
    config["atac_refinement"] = {
        "enabled": False,
        "fragment_maximum": 150,
        "macs3_qvalue": 0.1,
        "macs3_shift": -75,
        "macs3_extsize": 150,
        "bigwig_bin_size": 10,
        "minimum_mean_cpm": 2.0,
        "minimum_mode_prominence": 0.25,
        "merge_gap_bp": 1,
        "minimum_length": 50,
        "maximum_length": 400,
    }
    config["atac_atlas"] = {
        "enabled": True,
        "condition_map": "conditions.tsv",
        "peak_width": 250,
        "minimum_replicates": 2,
        "replicate_overlap_fraction": 0.5,
    }

    with pytest.raises(AcquisitionError, match="requires atac_refinement.enabled=true"):
        validate_workflow_config(config)


@pytest.mark.parametrize("prominence", [-0.01, 1.01])
def test_workflow_config_rejects_invalid_mode_prominence(prominence):
    config = copy.deepcopy(BASE_CONFIG)
    config["atac_refinement"] = dict(ATAC_REFINEMENT_DEFAULTS)
    config["atac_refinement"]["minimum_mode_prominence"] = prominence

    with pytest.raises(AcquisitionError, match="refinement parameters are invalid"):
        validate_workflow_config(config)


def test_secondary_macs3_dry_run_reuses_primary_bam(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["output_dir"] = str(tmp_path / "results")
    source_bam = (
        Path(config["output_dir"])
        / config["project"]
        / config["run_id"]
        / "bam"
        / "atac_rep1.final.bam"
    )
    source_bam.parent.mkdir(parents=True)
    source_bam.touch()
    source_bam.with_suffix(source_bam.suffix + ".bai").touch()
    config_path = tmp_path / "secondary-macs3.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))

    snakemake = Path(sys.executable).with_name("snakemake")
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(tmp_path / "cache")
    result = subprocess.run(
        [
            str(snakemake),
            "--snakefile",
            "workflow/secondary_macs3.smk",
            "--configfile",
            str(config_path),
            "--cores",
            "2",
            "--dry-run",
            "--config",
            "secondary_macs3_run_id=macs3-test",
            "secondary_macs3_shift=-75",
            "secondary_macs3_extsize=150",
        ],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert re.search(r"secondary_callpeak\s+1", output)


def test_primary_atac_workflow_includes_lenient_cpm_refinement(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["output_dir"] = str(tmp_path / "results")
    config["samples"][0]["peak_caller"] = {
        "command": "hmmratac",
        "lower": 10,
        "upper": 20,
        "prescan_cutoff": 1.2,
    }
    config["atac_refinement"] = {
        "enabled": True,
        "fragment_maximum": 150,
        "macs3_qvalue": 0.1,
        "macs3_shift": -75,
        "macs3_extsize": 150,
        "bigwig_bin_size": 10,
        "minimum_mean_cpm": 2.0,
        "minimum_mode_prominence": 0.25,
        "merge_gap_bp": 1,
        "minimum_length": 50,
        "maximum_length": 400,
    }
    config_path = tmp_path / "atac.yaml"
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
    assert re.search(r"hmmratac\s+1", output)
    assert re.search(r"filter_atac_short_fragments\s+1", output)
    assert re.search(r"atac_short_fragment_bigwig\s+1", output)
    assert re.search(r"call_lenient_atac_short_fragment_peaks\s+1", output)
    assert re.search(r"refine_atac_short_fragment_cpm\s+1", output)
    assert "tlen > -150 && tlen < 150" in output
    assert "-q 0.1" in output
    assert "--shift -75 --extsize 150" in output
    assert "--keep-dup all" in output
    assert "--minimum-mean-cpm 2.0" in output
    assert "--minimum-mode-prominence 0.25" in output


def test_optional_atac_atlas_branch_dry_run(tmp_path):
    config = copy.deepcopy(BASE_CONFIG)
    config["output_dir"] = str(tmp_path / "results")
    second = copy.deepcopy(config["samples"][0])
    second["id"] = "atac_rep2"
    second["accessions"] = ["SRR123457"]
    second["replicate"] = 2
    config["samples"].append(second)
    condition_map = tmp_path / "conditions.tsv"
    condition_map.write_text(
        "condition_id\tcondition_label\tsample_id\n"
        "embryo\tEmbryo\tatac_rep1\n"
        "embryo\tEmbryo\tatac_rep2\n"
    )
    config["atac_atlas"] = {
        "enabled": True,
        "condition_map": str(condition_map),
        "peak_width": 250,
        "minimum_replicates": 2,
        "replicate_overlap_fraction": 0.5,
    }
    config_path = tmp_path / "atac-atlas.yaml"
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
    assert re.search(r"merge_atac_condition_short_fragments\s+1", output)
    assert re.search(r"call_atac_condition_candidates\s+1", output)
    assert re.search(r"refine_atac_condition_cpm\s+1", output)
    assert re.search(r"filter_atac_condition_replicate_support\s+1", output)
    assert re.search(r"build_cross_condition_atac_atlas\s+1", output)
    assert re.search(
        r"build_narrow_first_cross_condition_atac_atlas\s+1", output
    )
    assert re.search(r"build_atac_dhs_support_fwhm\s+1", output)
    assert re.search(r"build_atac_dhs_center_mode_width\s+1", output)
    assert re.search(r"build_dhs_driven_atac_atlas\s+1", output)
    assert re.search(r"shape_dhs_driven_atac_atlas\s+1", output)
    assert "--replicate atac_rep1=" in output
    assert "--replicate atac_rep2=" in output
    assert "--peak-width 250" in output
    assert "--variable-bed" in output
    assert "--grouping-method fixed_window" in output
    assert "--grouping-method fixed_window_narrow_first" in output
    assert "--grouping-method dhs_seed" in output
    assert "--condition-bigwig embryo=" in output
    assert "--condition-dhs embryo=" in output
    assert "center-mode --anchors-bed" in output
    assert "--relative-threshold 0.2" in output


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
