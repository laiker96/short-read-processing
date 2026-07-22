"""Semantic validation for resolved workflow YAML files."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from .accessions import AcquisitionError


SAFE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
CHECKSUM_RE = re.compile(r"^(?:md5:[0-9a-f]{32}|sha256:[0-9a-f]{64})$")
ATAC_QPOIS_FIELDS = {
    "fragment_maximum",
    "minimum_exponent",
    "maximum_exponent",
    "merge_gap",
    "minimum_length",
    "maximum_length",
}
ATAC_CONSENSUS_FIELDS = {
    "enabled",
    "conditions",
    "minimum_replicates",
    "replicate_overlap_fraction",
}
REFERENCE_FIELDS = {
    "name",
    "fasta",
    "bowtie2_index",
    "chrom_sizes",
    "blacklist_bed",
    "tss_bed",
    "autosomes_file",
    "mitochondrial_contig",
    "effective_genome_size",
    "macs3_genome_size",
}


def wildcard_regex(values: list[str]) -> str:
    return "(?:" + "|".join(re.escape(value) for value in values) + ")" if values else r"(?!)"


def aria2_checksum(source: dict[str, Any]) -> str:
    """Translate a catalog checksum into aria2's ALGORITHM=DIGEST form."""

    algorithm, digest = str(source["checksum"]).split(":", 1)
    return f"{'sha-256' if algorithm == 'sha256' else algorithm}={digest}"


def _required(mapping: dict[str, Any], fields: set[str], label: str) -> None:
    missing = sorted(field for field in fields if field not in mapping or mapping[field] is None)
    if missing:
        raise AcquisitionError(f"{label} is missing: " + ", ".join(missing))


def validate_workflow_config(config: dict[str, Any]) -> None:
    """Reject inconsistent resolved configurations before Snakemake builds the DAG."""

    _required(config, {"project", "run_id", "output_dir", "assay", "reference", "samples"}, "Config")
    for field in ("project", "run_id"):
        if not SAFE_ID_RE.fullmatch(str(config[field])):
            raise AcquisitionError(f"Invalid {field}: {config[field]!r}")
    if config["assay"] not in {"atac", "chip_tf", "chip_histone"}:
        raise AcquisitionError(f"Unsupported assay: {config['assay']!r}")
    qpois = config.get("atac_qpois")
    if config["assay"] == "atac" and not isinstance(qpois, dict):
        raise AcquisitionError("ATAC configurations require atac_qpois parameters")
    if qpois is not None:
        if config["assay"] != "atac" or not isinstance(qpois, dict):
            raise AcquisitionError("atac_qpois is only valid for ATAC configurations")
        _required(qpois, ATAC_QPOIS_FIELDS, "ATAC qpois")
        if (
            int(qpois["fragment_maximum"]) < 2
            or int(qpois["minimum_exponent"]) < 0
            or int(qpois["maximum_exponent"]) < int(qpois["minimum_exponent"])
            or int(qpois["merge_gap"]) < 0
            or int(qpois["minimum_length"]) < 1
            or int(qpois["maximum_length"]) < int(qpois["minimum_length"])
        ):
            raise AcquisitionError("ATAC qpois parameters are invalid")
    consensus = config.get("atac_consensus")
    if consensus is not None:
        if config["assay"] != "atac" or not isinstance(consensus, dict):
            raise AcquisitionError("atac_consensus is only valid for ATAC configurations")
        _required(consensus, ATAC_CONSENSUS_FIELDS, "ATAC consensus")
        if not isinstance(consensus["enabled"], bool):
            raise AcquisitionError("ATAC consensus enabled must be true or false")
        if (
            int(consensus["minimum_replicates"]) < 2
            or not 0 < float(consensus["replicate_overlap_fraction"]) <= 1
            or not isinstance(consensus["conditions"], list)
            or not consensus["conditions"]
        ):
            raise AcquisitionError("ATAC consensus parameters are invalid")
    if not isinstance(config["reference"], dict):
        raise AcquisitionError("reference must be a mapping")
    _required(config["reference"], REFERENCE_FIELDS, "Reference")
    if int(config["reference"]["effective_genome_size"]) < 1:
        raise AcquisitionError("effective_genome_size must be positive")
    preparation = config["reference"].get("preparation")
    if preparation is not None:
        if not isinstance(preparation, dict) or preparation.get("mode") != "download":
            raise AcquisitionError("Reference preparation mode must be 'download'")
        _required(
            preparation,
            {"mode", "fasta", "annotation", "blacklist", "autosomes"},
            "Reference preparation",
        )
        for source_name in ("fasta", "annotation", "blacklist"):
            source = preparation[source_name]
            if not isinstance(source, dict):
                raise AcquisitionError(f"Reference {source_name} source must be a mapping")
            _required(source, {"url", "checksum"}, f"Reference {source_name} source")
            if not str(source["url"]).startswith("https://"):
                raise AcquisitionError(f"Reference {source_name} URL must use HTTPS")
            if not CHECKSUM_RE.fullmatch(str(source["checksum"])):
                raise AcquisitionError(
                    f"Reference {source_name} checksum must be md5:<hex> or sha256:<hex>"
                )
        autosomes = preparation["autosomes"]
        if not isinstance(autosomes, list) or not autosomes or any(not item for item in autosomes):
            raise AcquisitionError("Reference autosomes must be a non-empty list")

    samples = config["samples"]
    if not isinstance(samples, list) or not samples:
        raise AcquisitionError("samples must be a non-empty list")
    sample_ids = [str(sample.get("id", "")) for sample in samples]
    if any(not SAFE_ID_RE.fullmatch(sample_id) for sample_id in sample_ids):
        raise AcquisitionError("Every sample must have a safe non-empty id")
    if len(sample_ids) != len(set(sample_ids)):
        raise AcquisitionError("Sample IDs must be unique")
    role_by_id = {str(sample["id"]): sample.get("role") for sample in samples}
    context_by_id = {str(sample["id"]): str(sample.get("context", "")) for sample in samples}

    for sample in samples:
        sample_id = str(sample["id"])
        _required(
            sample,
            {"accessions", "context", "role", "layout", "r1", "parameters"},
            f"Sample {sample_id}",
        )
        if not SAFE_ID_RE.fullmatch(str(sample["context"])):
            raise AcquisitionError(f"Sample {sample_id}: invalid context")
        if sample["role"] not in {"treatment", "control"}:
            raise AcquisitionError(f"Sample {sample_id}: invalid role")
        if sample["layout"] not in {"single", "paired"}:
            raise AcquisitionError(f"Sample {sample_id}: invalid layout")
        if not isinstance(sample["r1"], list) or not sample["r1"]:
            raise AcquisitionError(f"Sample {sample_id}: r1 must contain at least one FASTQ")
        if sample["layout"] == "paired":
            if not isinstance(sample.get("r2"), list) or len(sample["r1"]) != len(sample["r2"]):
                raise AcquisitionError(f"Sample {sample_id}: paired r1/r2 lane counts differ")
        elif sample.get("r2"):
            raise AcquisitionError(f"Sample {sample_id}: single-end input must not contain r2")

        if sample["role"] == "treatment":
            peak = sample.get("peak_caller")
            if not isinstance(peak, dict) or peak.get("command") not in {"callpeak", "hmmratac"}:
                raise AcquisitionError(f"Sample {sample_id}: treatment requires a peak caller")
            if peak["command"] == "hmmratac":
                if config["assay"] != "atac" or sample["layout"] != "paired":
                    raise AcquisitionError(
                        f"Sample {sample_id}: HMMRATAC requires paired-end ATAC-seq"
                    )
            else:
                _required(
                    peak,
                    {"format", "qvalue", "broad", "nomodel", "write_bedgraph", "spmr"},
                    f"Sample {sample_id} callpeak",
                )
                if config["assay"] == "atac":
                    if (
                        peak.get("mode") != "tn5_qpois"
                        or peak["format"] != "BED"
                        or peak["broad"]
                        or not peak["nomodel"]
                        or not peak["write_bedgraph"]
                        or peak["spmr"]
                        or peak.get("shift") is None
                        or peak.get("extsize") is None
                    ):
                        raise AcquisitionError(
                            f"Sample {sample_id}: invalid two-ended Tn5 qpois configuration"
                        )
                elif not peak["write_bedgraph"] or not peak["spmr"]:
                    raise AcquisitionError(
                        f"Sample {sample_id}: ChIP callpeak must write -B --SPMR bedGraphs"
                    )

        if config["assay"].startswith("chip") and sample["role"] == "treatment":
            control = str(sample.get("control") or "")
            if control and role_by_id.get(control) != "control":
                raise AcquisitionError(f"Sample {sample_id}: invalid matched ChIP control")
            if control and context_by_id[control] != str(sample["context"]):
                raise AcquisitionError(
                    f"Sample {sample_id}: treatment and control contexts differ"
                )
        if config["assay"] == "atac" and sample["role"] != "treatment":
            raise AcquisitionError("ATAC configurations cannot contain control samples")

    if consensus and consensus["enabled"] and config["assay"] != "atac":
        raise AcquisitionError("ATAC consensus requires an ATAC configuration")
    if consensus and consensus["enabled"]:
        from .consensus import condition_specs

        specifications = condition_specs(
            consensus["conditions"],
            sample_ids=[
                sample_id
                for sample_id in sample_ids
                if role_by_id[sample_id] == "treatment"
            ],
            minimum_replicates=int(consensus["minimum_replicates"]),
        )
        for specification in specifications:
            contexts = {
                context_by_id[sample_id] for sample_id in specification.samples
            }
            if contexts != {specification.condition_id}:
                raise AcquisitionError(
                    f"ATAC condition {specification.condition_id!r} does not match "
                    "its sample contexts"
                )


def resolve_input_paths(config: dict[str, Any], base: Path) -> None:
    """Resolve relative FASTQ/reference paths against a config's launch directory in place."""

    reference = config["reference"]
    for key in ("fasta", "bowtie2_index", "chrom_sizes", "blacklist_bed", "tss_bed", "autosomes_file"):
        path = Path(reference[key])
        reference[key] = str(path if path.is_absolute() else (base / path).resolve())
    for sample in config["samples"]:
        for key in ("r1", "r2"):
            if key in sample:
                sample[key] = [
                    str(path if path.is_absolute() else (base / path).resolve())
                    for value in sample[key]
                    for path in [Path(value)]
                ]
        adapter = sample["parameters"]["trimming"].get("adapter_fasta")
        if adapter:
            path = Path(adapter)
            sample["parameters"]["trimming"]["adapter_fasta"] = str(
                path if path.is_absolute() else (base / path).resolve()
            )
