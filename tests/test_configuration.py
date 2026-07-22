from pathlib import Path

import pytest
import yaml

from short_read_processing.accessions import AcquisitionError, FilePlan, RunPlan
from short_read_processing.configuration import generate_configs
from short_read_processing.manifest import write_manifest


HEADER = "accession\tlibrary_id\tassay\tcontext\trole\tcontrol_library\tpeak_caller"


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


def _generate(tmp_path: Path, plans: list[RunPlan], sheet_text: str, **kwargs):
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
        **kwargs,
    )


def test_atac_defaults_and_contexts_group_biological_and_technical_runs(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", accession, accession)
        for accession in ("SRR123456", "SRR123457", "SRR123458")
    ]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123458\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    config = yaml.safe_load(output.read_text())
    sample = config["samples"][0]

    assert config["assay"] == "atac"
    assert config["reference"]["name"] == "dm6"
    assert sample["accessions"] == ["SRR123456", "SRR123457"]
    assert sample["context"] == "embryo"
    assert sample["layout"] == "paired"
    assert len(sample["r1"]) == 2
    assert len(sample["r2"]) == 2
    assert sample["peak_caller"] == {
        "command": "callpeak",
        "mode": "tn5_qpois",
        "format": "BED",
        "qvalue": 0.1,
        "broad": False,
        "nomodel": True,
        "shift": -75,
        "extsize": 150,
        "write_bedgraph": True,
        "spmr": False,
    }
    assert config["atac_consensus"] == {
        "enabled": True,
        "conditions": [
            {
                "id": "embryo",
                "label": "embryo",
                "samples": ["atac_rep1", "atac_rep2"],
            }
        ],
        "minimum_replicates": 2,
        "replicate_overlap_fraction": 0.5,
    }
    assert sample["parameters"]["trimming"]["adapter_preset"] == "nextera"
    preparation = config["reference"]["preparation"]
    assert preparation["mode"] == "download"
    assert preparation["fasta"]["checksum"].startswith("md5:")
    assert preparation["annotation"]["url"].endswith("dm6.ncbiRefSeq.gtf.gz")


def test_identical_config_generation_does_not_replace_file(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    sheet = (
        HEADER
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\t"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
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


def test_atac_callpeak_override_applies_to_every_context_replicate(tmp_path):
    accessions = ("SRR123456", "SRR123457")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    header = HEADER + "\tmacs3_qvalue\tmacs3_shift\tmacs3_extsize"
    sheet = (
        header
        + "\nSRR123456\tatac_rep1\tatac\tembryo\ttreatment\t\tcallpeak\t0.01\t-75\t150"
        + "\nSRR123457\tatac_rep2\tatac\tembryo\ttreatment\t\tcallpeak\t0.01\t-75\t150\n"
    )
    output = _generate(tmp_path, plans, sheet)[0]
    peaks = [sample["peak_caller"] for sample in yaml.safe_load(output.read_text())["samples"]]

    assert all(peak["format"] == "BED" for peak in peaks)
    assert all(peak["qvalue"] == 0.01 for peak in peaks)
    assert all(peak["shift"] == -75 and peak["extsize"] == 150 for peak in peaks)
    assert all(peak["spmr"] is False for peak in peaks)


def test_h3k27ac_alias_is_broad_and_resolves_matched_input(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR100001", "SRR100001"),
        _run_plan(tmp_path / "raw", "SRR100002", "SRR100002"),
    ]
    sheet = (
        HEADER
        + "\nSRR100001\th3_rep1\th3k27ac\teye\ttreatment\teye_input\t"
        + "\nSRR100002\teye_input\th3k27ac\teye\tcontrol\t\t\n"
    )
    output = _generate(tmp_path, plans, sheet, genome="hg38")[0]
    config = yaml.safe_load(output.read_text())
    treatment = config["samples"][0]

    assert config["assay"] == "chip_histone"
    assert config["reference"]["name"] == "hg38"
    assert treatment["control"] == "eye_input"
    assert treatment["peak_caller"]["command"] == "callpeak"
    assert treatment["peak_caller"]["broad"] is True
    assert treatment["peak_caller"]["broad_cutoff"] == 0.1
    assert "peak_caller" not in config["samples"][1]


def test_matched_ip_and_input_must_have_same_resolved_layout(tmp_path):
    plans = [
        _run_plan(tmp_path / "raw", "SRR100001", "SRR100001", layout="PAIRED"),
        _run_plan(tmp_path / "raw", "SRR100002", "SRR100002", layout="SINGLE"),
    ]
    sheet = (
        HEADER
        + "\nSRR100001\th3_rep1\th3k27ac\teye\ttreatment\teye_input\t"
        + "\nSRR100002\teye_input\th3k27ac\teye\tcontrol\t\t\n"
    )

    with pytest.raises(AcquisitionError, match="same read layout"):
        _generate(tmp_path, plans, sheet)


def test_mixed_atlas_table_generates_separate_atac_and_h3k27ac_configs(tmp_path):
    accessions = ("SRR100001", "SRR100002", "SRR100003", "SRR100004")
    plans = [_run_plan(tmp_path / "raw", accession, accession) for accession in accessions]
    sheet = (
        HEADER
        + "\nSRR100001\teye_atac_rep1\tatac\teye\ttreatment\t\t"
        + "\nSRR100002\teye_atac_rep2\tatac\teye\ttreatment\t\t"
        + "\nSRR100003\teye_h3_rep1\th3k27ac\teye\ttreatment\teye_input\t"
        + "\nSRR100004\teye_input\th3k27ac\teye\tcontrol\t\t\n"
    )

    outputs = _generate(tmp_path, plans, sheet)
    configs = {yaml.safe_load(path.read_text())["assay"]: path for path in outputs}

    assert set(configs) == {"atac", "chip_histone"}
    assert configs["atac"].name == "test-project.atac.dm6.yaml"
    assert configs["chip_histone"].name == "test-project.chip_histone.dm6.yaml"
    chip = yaml.safe_load(configs["chip_histone"].read_text())
    assert chip["samples"][0]["control"] == "eye_input"


def test_h3k27ac_can_run_ip_only_from_four_columns(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001")
    sheet = (
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\th3_rep1\th3k27ac\teye\n"
    )
    output = _generate(tmp_path, [plan], sheet)[0]
    treatment = yaml.safe_load(output.read_text())["samples"][0]

    assert "control" not in treatment
    assert treatment["peak_caller"]["format"] == "BAMPE"
    assert treatment["peak_caller"]["broad"] is True
    assert treatment["parameters"]["trimming"]["adapter_preset"] == "truseq"


def test_atac_context_requires_two_biological_libraries(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001")
    sheet = (
        "accession\tlibrary_id\tassay\tcontext\n"
        "SRR100001\tatac_rep1\tatac\teye\n"
    )
    with pytest.raises(AcquisitionError, match="below minimum_replicates"):
        _generate(tmp_path, [plan], sheet)


def test_explicit_hmmratac_rejects_single_end_atac(tmp_path):
    plan = _run_plan(tmp_path / "raw", "SRR100001", "SRR100001", layout="SINGLE")
    sheet = (
        HEADER
        + "\nSRR100001\tatac_rep1\tatac\teye\ttreatment\t\thmmratac\n"
    )
    with pytest.raises(AcquisitionError, match="requires paired-end"):
        _generate(tmp_path, [plan], sheet)
