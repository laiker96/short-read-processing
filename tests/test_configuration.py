from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError, FilePlan, RunPlan
from short_read_processing.configuration import generate_configs
from short_read_processing.manifest import write_manifest


HEADER = "accession\tsample_id\tassay\tgenome\trole\tcontrol_id\treplicate\tpeak_caller"


def _run_plan(root: Path, requested: str, run: str, layout: str = "PAIRED") -> RunPlan:
    run_dir = root / run
    run_dir.mkdir(parents=True)
    r1 = run_dir / f"{run}_1.fastq.gz"
    r1.write_bytes(b"r1")
    files = [FilePlan("https://example/r1", "", 2, r1, "r1")]
    if layout == "PAIRED":
        r2 = run_dir / f"{run}_2.fastq.gz"
        r2.write_bytes(b"r2")
        files.append(FilePlan("https://example/r2", "", 2, r2, "r2"))
    return RunPlan(
        requested_accession=requested,
        experiment_accession=requested if requested.startswith(("SRX", "ERX")) else "SRX999999",
        run_accession=run,
        library_layout=layout,
        backend="ena",
        run_dir=run_dir,
        files=files,
        status="downloaded",
    )


def _generate(tmp_path: Path, plans: list[RunPlan], sheet_text: str):
    manifest = tmp_path / "manifest.tsv"
    write_manifest(manifest, plans)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(sheet_text)
    return generate_configs(
        manifest_path=manifest,
        sample_sheet_path=sheet,
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="baseline",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
    )


def test_atac_defaults_to_hmmratac_and_groups_technical_runs(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR123456", "SRR123456"),
        _run_plan(tmp_path / "raw", "SRR123457", "SRR123457"),
    ]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tdm6\ttreatment\t\t1\t"
        + "\nSRR123457\tatac_rep1\tatac\tdm6\ttreatment\t\t1\t\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    config = yaml.safe_load(output.read_text())
    sample = config["samples"][0]

    assert config["assay"] == "atac"
    assert sample["accessions"] == ["SRR123456", "SRR123457"]
    assert sample["layout"] == "paired"
    assert len(sample["r1"]) == 2
    assert len(sample["r2"]) == 2
    assert sample["peak_caller"]["command"] == "hmmratac"
    assert sample["peak_caller"]["lower"] == 10
    assert sample["parameters"]["trimming"]["adapter_preset"] == "nextera"
    preparation = config["reference"]["preparation"]
    assert preparation["mode"] == "download"
    assert preparation["fasta"]["checksum"].startswith("md5:")
    assert preparation["annotation"]["url"].endswith("dm6.ncbiRefSeq.gtf.gz")
    assert preparation["autosomes"] == ["chr2L", "chr2R", "chr3L", "chr3R", "chr4"]
    assert config["atac_refinement"] == {
        "enabled": True,
        "fragment_maximum": 150,
        "macs3_qvalue": 0.1,
        "macs3_shift": -75,
        "macs3_extsize": 150,
        "bigwig_bin_size": 10,
        "minimum_mean_cpm": 2.0,
        "merge_gap_bp": 1,
        "minimum_length": 50,
        "maximum_length": 400,
    }


def test_identical_config_generation_does_not_replace_file(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR123456", "SRR123456")
    sheet = HEADER + "\nSRR123456\tatac_rep1\tatac\tdm6\ttreatment\t\t1\t\n"
    output = _generate(tmp_path, [plan], sheet)[0]
    original = output.read_bytes()
    output.touch()
    timestamp = output.stat().st_mtime_ns

    regenerated = generate_configs(
        manifest_path=tmp_path / "manifest.tsv",
        sample_sheet_path=tmp_path / "samples.tsv",
        output_dir=tmp_path / "configs",
        project="test-project",
        run_id="baseline",
        reference_root=tmp_path / "references",
        path_base=tmp_path,
        require_fastq_files=True,
    )[0]

    assert regenerated.read_bytes() == original
    assert regenerated.stat().st_mtime_ns == timestamp


def test_atac_callpeak_shift_override_writes_bedgraphs(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR123456", "SRR123456")
    header = HEADER + "\tmacs3_format\tmacs3_qvalue\tmacs3_nomodel\tmacs3_shift\tmacs3_extsize"
    sheet = (
        header
        + "\nSRR123456\tatac_rep1\tatac\tdm6\ttreatment\t\t1\tcallpeak"
        + "\tBAM\t0.01\ttrue\t-75\t150\n"
    )
    output = _generate(tmp_path, [plan], sheet)[0]
    peak = yaml.safe_load(output.read_text())["samples"][0]["peak_caller"]

    assert peak["command"] == "callpeak"
    assert peak["format"] == "BAM"
    assert peak["qvalue"] == 0.01
    assert peak["nomodel"] is True
    assert peak["shift"] == -75
    assert peak["extsize"] == 150
    assert peak["write_bedgraph"] is True
    assert peak["spmr"] is True
    assert peak["bedgraph_outputs"]["treatment_suffix"] == "_treat_pileup.bdg"
    assert peak["bedgraph_outputs"]["control_suffix"] == "_control_lambda.bdg"


def test_histone_callpeak_is_broad_and_resolves_control(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR100001", "SRR100001"),
        _run_plan(tmp_path / "raw", "SRR100002", "SRR100002"),
    ]
    sheet = (
        HEADER
        + "\nSRR100001\tH3K27ac_rep1\tchip_histone\thg38\ttreatment\tinput_rep1\t1\t"
        + "\nSRR100002\tinput_rep1\tchip_histone\thg38\tcontrol\t\t1\t\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    config = yaml.safe_load(output.read_text())
    treatment = config["samples"][0]

    assert treatment["control"] == "input_rep1"
    assert treatment["peak_caller"]["command"] == "callpeak"
    assert treatment["peak_caller"]["broad"] is True
    assert treatment["peak_caller"]["broad_cutoff"] == 0.1
    assert treatment["peak_caller"]["write_bedgraph"] is True
    assert "peak_caller" not in config["samples"][1]


def test_histone_callpeak_can_run_without_control(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001")
    sheet = (
        HEADER
        + "\nSRR100001\tH3K27ac_rep1\tchip_histone\tdm6\ttreatment\t\t1\tcallpeak\n"
    )
    output = _generate(tmp_path, [plan], sheet)[0]
    treatment = yaml.safe_load(output.read_text())["samples"][0]

    assert "control" not in treatment
    assert treatment["peak_caller"]["format"] == "BAMPE"
    assert treatment["peak_caller"]["broad"] is True
    assert treatment["peak_caller"]["broad_cutoff"] == 0.1
    assert treatment["parameters"]["trimming"]["adapter_preset"] == "truseq"


def test_hmmratac_default_rejects_single_end_atac(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001", layout="SINGLE")
    sheet = HEADER + "\nSRR100001\tatac_rep1\tatac\tdm6\ttreatment\t\t1\t\n"
    with pytest.raises(AcquisitionError, match="requires paired-end"):
        _generate(tmp_path, [plan], sheet)
