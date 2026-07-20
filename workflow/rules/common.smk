import shlex

from short_read_processing.macs3 import callpeak_arguments


def aria2_checksum(source):
    algorithm, digest = source["checksum"].split(":", 1)
    if algorithm == "sha256":
        algorithm = "sha-256"
    return f"{algorithm}={digest}"


def raw_fastq(wildcards):
    return RAW_FASTQ_BY_UNIT[wildcards.unit]


def trimmed_fastq(wildcards):
    return TRIMMED_FASTQ_BY_UNIT[wildcards.unit]


def lane_input(wildcards, mate):
    return LANES[(wildcards.sample, wildcards.lane)][mate]


def trim_parameter(wildcards, section, name):
    return SAMPLES[wildcards.sample]["parameters"][section][name]


def adapter_value(wildcards, mate):
    trimming = SAMPLES[wildcards.sample]["parameters"]["trimming"]
    preset = trimming["adapter_preset"]
    if preset == "custom":
        return f"file:{trimming['adapter_fasta']}"
    if preset == "nextera":
        return NEXTERA_ADAPTER
    return TRUSEQ_R1_ADAPTER if mate == "r1" else TRUSEQ_R2_ADAPTER


def trimmed_lane_reads(wildcards):
    key = (wildcards.sample, wildcards.lane)
    reads = [TRIMMED_R1[key]]
    if SAMPLES[wildcards.sample]["layout"] == "paired":
        reads.append(TRIMMED_R2[key])
    return reads


def bowtie_lane_arguments(wildcards):
    key = (wildcards.sample, wildcards.lane)
    r1 = shlex.quote(TRIMMED_R1[key])
    if SAMPLES[wildcards.sample]["layout"] == "paired":
        r2 = shlex.quote(TRIMMED_R2[key])
        return f"-1 {r1} -2 {r2}"
    return f"-U {r1}"


def sample_lane_bams(wildcards):
    return [
        LANE_BAMS[(sample, lane)]
        for sample, lane in LANES
        if sample == wildcards.sample
    ]


def worker_threads(wildcards, threads):
    return max(1, threads - 1)


def bowtie_layout_arguments(wildcards):
    sample = SAMPLES[wildcards.sample]
    if sample["layout"] == "paired":
        maximum = sample["parameters"]["alignment"]["maximum_fragment_length"]
        return f"--no-mixed --no-discordant -X {maximum}"
    return ""


def bowtie_preset(wildcards):
    return "--" + SAMPLES[wildcards.sample]["parameters"]["alignment"]["preset"]


def required_flags(wildcards):
    return "-f 2" if SAMPLES[wildcards.sample]["layout"] == "paired" else ""


def excluded_flags(wildcards):
    excluded = 4 + 256 + 512 + 2048
    if SAMPLES[wildcards.sample]["layout"] == "paired":
        excluded += 8
    if SAMPLES[wildcards.sample]["parameters"]["filtering"]["remove_duplicates"]:
        excluded += 1024
    return excluded


def peak_control_input(wildcards):
    control = SAMPLES[wildcards.sample].get("control")
    return [FINAL_BAMS[control]] if control else []


def callpeak_argv(wildcards):
    sample = wildcards.sample
    control = SAMPLES[sample].get("control")
    arguments = callpeak_arguments(
        SAMPLES[sample]["peak_caller"],
        treatment_bam=Path(FINAL_BAMS[sample]),
        control_bam=Path(FINAL_BAMS[control]) if control else None,
        name=sample,
        genome_size=REFERENCE["macs3_genome_size"],
        output_dir=Path(f"{RESULT_ROOT}/peaks/{sample}"),
    )
    return shlex.join(arguments)


def atac_refinement_callpeak_argv(wildcards):
    refinement = ATAC_REFINEMENT
    peak_config = {
        "command": "callpeak",
        "format": "BAM",
        "qvalue": refinement["macs3_qvalue"],
        "broad": False,
        "nomodel": True,
        "shift": refinement["macs3_shift"],
        "extsize": refinement["macs3_extsize"],
        "write_bedgraph": True,
        "spmr": True,
    }
    arguments = callpeak_arguments(
        peak_config,
        treatment_bam=Path(ATAC_SHORT_BAMS[wildcards.sample]),
        control_bam=None,
        name=wildcards.sample,
        genome_size=REFERENCE["macs3_genome_size"],
        output_dir=Path(f"{ATAC_REFINEMENT_ROOT}/macs3/{wildcards.sample}"),
    )
    arguments.extend(["--keep-dup", "all"])
    return shlex.join(arguments)
