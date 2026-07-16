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


def bowtie_read_arguments(wildcards):
    sample = wildcards.sample
    lanes = [lane for candidate, lane in LANES if candidate == sample]
    r1 = [TRIMMED_R1[(sample, lane)] for lane in lanes]
    if SAMPLES[sample]["layout"] == "paired":
        r2 = [TRIMMED_R2[(sample, lane)] for lane in lanes]
        return "-1 " + shlex.quote(",".join(r1)) + " -2 " + shlex.quote(",".join(r2))
    return "-U " + shlex.quote(",".join(r1))


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
