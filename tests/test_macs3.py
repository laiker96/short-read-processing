from pathlib import Path

from short_read_processing.macs3 import callpeak_arguments, callpeak_bedgraph_paths


def test_callpeak_argv_and_bedgraph_outputs():
    config = {
        "command": "callpeak",
        "format": "BAM",
        "qvalue": 0.01,
        "broad": False,
        "nomodel": True,
        "shift": -75,
        "extsize": 150,
        "write_bedgraph": True,
        "spmr": True,
    }
    arguments = callpeak_arguments(
        config,
        treatment_bam=Path("treatment.bam"),
        control_bam=Path("input.bam"),
        name="atac_rep1",
        genome_size="dm",
        output_dir=Path("peaks"),
    )

    assert arguments[:2] == ["macs3", "callpeak"]
    assert arguments[arguments.index("-q") + 1] == "0.01"
    assert "-B" in arguments
    assert "--SPMR" in arguments
    assert arguments[arguments.index("-c") + 1] == "input.bam"
    assert arguments[arguments.index("--shift") + 1] == "-75"
    assert arguments[arguments.index("--extsize") + 1] == "150"
    assert callpeak_bedgraph_paths(Path("peaks"), "atac_rep1") == (
        Path("peaks/atac_rep1_treat_pileup.bdg"),
        Path("peaks/atac_rep1_control_lambda.bdg"),
    )


def test_broad_callpeak_omits_external_control_when_ip_only():
    config = {
        "command": "callpeak",
        "format": "BAMPE",
        "qvalue": 0.01,
        "broad": True,
        "broad_cutoff": 0.1,
        "nomodel": False,
        "shift": None,
        "extsize": None,
        "write_bedgraph": True,
        "spmr": True,
    }
    arguments = callpeak_arguments(
        config,
        treatment_bam=Path("h3k27ac.bam"),
        control_bam=None,
        name="h3k27ac_rep1",
        genome_size="dm",
        output_dir=Path("peaks"),
    )

    assert "-c" not in arguments
    assert "--broad" in arguments
    assert arguments[arguments.index("--broad-cutoff") + 1] == "0.1"
