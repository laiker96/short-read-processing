# short-read-processing

Accession-first Snakemake pipeline for bulk ATAC-seq, TF ChIP-seq, and histone
ChIP-seq. A CSV or TSV sample sheet is the only user data input: the pipeline
resolves SRR/SRX/ERR/ERX accessions, downloads FASTQs, prepares `dm6` or `hg38`,
aligns and filters reads, calls peaks, creates signal tracks, and summarizes QC.

Paired-end ATAC defaults to two complementary peak products in one run:

- MACS3 HMMRATAC accessible regions;
- lenient MACS3 candidates from `<150 bp` fragments, followed by CPM-based
  refinement.

![short-read-processing workflow DAG](docs/workflow-dag.svg)

The figure is the Snakemake rule graph for a representative paired-end ATAC
sample; independent samples and technical lanes expand these same branches in
parallel.

The public entry point is `src/run_pipeline.py`. Acquisition and resolved-config
generation are restart-safe pre-DAG phases; all reference preparation and read
processing steps are in `workflow/Snakefile`.

## Requirements

- Linux or a compatible HPC compute node;
- Mamba or Micromamba;
- enough disk for compressed FASTQs, references, BAMs, and rule environments;
- outbound HTTPS for accession/reference downloads, unless inputs and caches
  have already been staged.

Run commands from the repository root. All environments and workflow caches are
kept inside the repository:

```text
.venv/                 orchestration environment
.micromamba/           optional repo-local Mamba root/cache
.snakemake/conda/      per-rule bioinformatics environments
```

These paths are ignored by Git.

## Install the environment

With Mamba:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
mamba env create --prefix "$PWD/.venv" --file environment.yml
mamba run --prefix "$PWD/.venv" snakemake --version
```

With Micromamba, replace `mamba` with `micromamba`. To update an existing
environment after `environment.yml` changes:

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"
mamba env update --prefix "$PWD/.venv" --file environment.yml --prune
```

The orchestration environment contains Snakemake, download tools, and test
dependencies. Bowtie2, MACS3, SAMtools, deepTools, FastQC, Cutadapt, MultiQC,
and the R-based ChIP QC tools are isolated by task under `workflow/envs/`.
Snakemake creates those environments automatically below `.snakemake/conda` on
the first run. R packages are not installed into `.venv`.

## Prepare the input sample sheet

The input must be a comma- or tab-separated file following
[`schemas/sample-sheet.schema.yaml`](schemas/sample-sheet.schema.yaml). Do not
list local FASTQ paths. One row represents one public run or experiment
accession; SRX and ERX accessions expand to all associated runs.

Required columns:

| Column | Meaning |
|---|---|
| `accession` | SRR, SRX, ERR, or ERX identifier |
| `sample_id` | Biological library ID; repeat for technical runs of the same library |
| `assay` | `atac`, `chip_tf`, or `chip_histone` |
| `genome` | `dm6` or `hg38` |
| `role` | `treatment` or `control` |
| `control_id` | Matched ChIP control `sample_id`; blank for ATAC, controls, or IP-only ChIP |
| `replicate` | Positive biological-replicate number |
| `peak_caller` | `hmmratac`, `callpeak`, or blank for the assay default |

Minimal paired-end ATAC input:

```text
accession  sample_id       assay  genome  role       control_id  replicate  peak_caller
ERR3975804 e5_atac_rep1    atac   dm6    treatment              1          hmmratac
ERR3975777 e5_atac_rep1    atac   dm6    treatment              1          hmmratac
ERR3975789 e5_atac_rep2    atac   dm6    treatment              2          hmmratac
```

The first two rows are technical runs because they share `sample_id`; they are
trimmed and aligned independently, then merged before duplicate marking. The
third row is a separate biological replicate and remains a separate sample.

ChIP with an explicit input control uses a control row and links the treatment
through `control_id`:

```text
accession  sample_id  assay    genome  role       control_id  replicate  peak_caller
SRR100001  tf_rep1    chip_tf  hg38    treatment  input_rep1  1          callpeak
SRR100002  input_rep1 chip_tf  hg38    control                1
```

IP-only ChIP is supported by leaving `control_id` blank. MACS3 then runs without
`-c`; control-dependent fingerprint QC is omitted.

Optional typed columns control MACS3/HMMRATAC, trimming, alignment, filtering,
and adapters. Important examples include:

```text
macs3_format=BAM
macs3_nomodel=true
macs3_shift=-75
macs3_extsize=150
adapter_preset=custom
adapter_fasta=resources/adapters/nextera.fa
```

Free-form shell arguments are deliberately unsupported. Use the named columns
defined by the schema. Preset defaults are Nextera for ATAC, TruSeq for ChIP,
Bowtie2 `very-sensitive`, MAPQ 30, duplicate removal, and mitochondrial-read
removal.

Validate a sheet without downloading data:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/validate_sample_sheet.py samples.tsv
```

## Run the pipeline locally

This command validates the sheet, downloads FASTQs concurrently, writes the
resolved YAML, prepares the reference, and runs Snakemake:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py samples.tsv \
  --project chromatin-study \
  --run-id baseline \
  --output-dir data/raw/chromatin-study \
  --reference-root references \
  --cores 24 \
  --file-jobs 8 \
  --connections 8
```

The main locations are then:

```text
data/raw/chromatin-study/               downloaded FASTQs and manifest
configs/chromatin-study.yaml            fully resolved workflow config
references/<genome>/                    downloaded/prepared reference
results/chromatin-study/baseline/       final outputs and logs
work/chromatin-study/baseline/          restartable intermediates
```

Independent accessions, lanes, and samples run concurrently. `--cores` is the
aggregate local CPU limit; each rule also declares its own threads and memory.
For SRA Toolkit fallbacks, `--sra-jobs` controls simultaneous conversions and
`--threads` is divided among those jobs. For direct ENA transfers,
`--file-jobs` controls simultaneous files and `--connections` controls HTTP
range connections per file.

Useful execution boundaries:

```bash
# Resolve and download only
python src/run_pipeline.py samples.tsv --download-only --output-dir data/raw/project

# Generate YAML from an existing manifest without running Snakemake
python src/run_pipeline.py samples.tsv --skip-download \
  --manifest data/raw/project/download_manifest.tsv --config-only

# Inspect the complete DAG without executing jobs
python src/run_pipeline.py samples.tsv --skip-download \
  --manifest data/raw/project/download_manifest.tsv --snakemake-dry-run
```

Run those commands through `mamba run --prefix "$PWD/.venv"` as above.

### Run the dm6 atlas inputs

The curated atlas sheets are ready to use:

```bash
# ATAC: 23 accessions / 22 biological libraries
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_atac_selected.sample_sheet.tsv \
  --project atlas-atac-dm6 --run-id baseline \
  --output-dir data/raw/atlas_atac --cores 24

# H3K27ac: 15 IP-only runs
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py resources/atlas_h3k27ac_ip_only.sample_sheet.tsv \
  --project atlas-h3k27ac-dm6 --run-id baseline \
  --output-dir data/raw/atlas_h3k27ac --cores 24
```

The source selection metadata and deterministic sample-sheet generators are
documented in [`resources/README.md`](resources/README.md).

### Restart or resume

Re-run the same command with the same `project` and `run-id`:

- aria2 resumes managed ENA partial downloads and verifies reported MD5 sums;
- SRA conversion writes to staging and promotes FASTQs only after completion;
- manifests and resolved YAMLs are atomic and unchanged files are not replaced;
- Snakemake skips complete outputs and reruns incomplete jobs;
- raw FASTQs and completed canonical BAMs are never overwritten by downstream
  rules.

To reuse completed downloads explicitly, add:

```bash
--skip-download --manifest data/raw/project/download_manifest.tsv
```

Use a new `run-id` when changing scientific parameters so previous results are
preserved. The ATAC short-fragment refinement defaults are recorded in the
resolved YAML under `atac_refinement`.

## Reference preparation

Generated `dm6` and `hg38` configurations include checksum-pinned reference
sources. The Snakemake DAG downloads and prepares:

```text
references/<genome>/<genome>.fa
references/<genome>/<genome>.fa.fai
references/<genome>/<genome>.chrom.sizes
references/<genome>/<genome>.blacklist.bed
references/<genome>/<genome>.tss.bed
references/<genome>/<genome>.autosomes.txt
references/<genome>/bowtie2/<genome>.*.bt2
references/<genome>/sources/*
```

A hand-written resolved YAML may omit `reference.preparation` and point to
existing local assets. FASTA indexing, chromosome sizes, and missing Bowtie2
indexes remain workflow targets.

## Outputs

Each run writes below `results/<project>/<run-id>/`:

```text
bam/
  <sample>.final.bam[.bai]             filtered canonical alignments

peaks/<sample>/
  *_accessible_regions.narrowPeak      HMMRATAC ATAC peaks
  *_peaks.narrowPeak                   MACS3 narrow peaks
  *_peaks.broadPeak                    MACS3 broad histone peaks
  *_treat_pileup.bdg                   callpeak treatment signal
  *_control_lambda.bdg                 callpeak local background

tracks/
  <sample>.CPM.bw                      unshifted CPM coverage
  <sample>.Tn5-shifted.CPM.bw          ATAC insertion-oriented coverage

atac_short_fragments/
  bam/*.fragments-lt150.bam[.bai]      proper pairs with 0 < |TLEN| < 150
  tracks/*.Tn5-shifted.CPM.bw          CPM within the retained subset
  macs3/<sample>/*_peaks.narrowPeak    lenient q=0.10 candidates
  macs3/<sample>/*_{treat_pileup,control_lambda}.bdg
  refined/*.CPM-refined.bed            50-400 bp signal-refined peaks
  refined/*.Excluded.bed               final unselected signal intervals
  refined/*.stats.json                 refinement counts and thresholds
  qc/*.fragment-filter.tsv             retained-fragment statistics

qc/
  fastqc/raw/ and fastqc/trimmed/       per-FASTQ FastQC reports
  cutadapt/                             trimming JSON
  alignment/                            flagstat, stats, and idxstats
  frip/                                 numerator, denominator, and FRiP
  tss/ and fragments/                   ATAC TSS and fragment-size QC
  chip/                                 ChIP fingerprint/cross-correlation
  metrics.tsv and metrics.json          stable machine-readable summary
  multiqc/multiqc_report.html           aggregate report

provenance/resolved_config.json         resolved run configuration
logs/                                   commands and tool logs by stage
```

MACS3 `callpeak` uses `-B --SPMR`, so both treatment and control-lambda
bedGraphs are declared outputs. Paired-end ATAC defaults to HMMRATAC for its
primary peaks. The additional short-fragment branch uses MACS3 `-f BAM
--nomodel --shift -75 --extsize 150 -q 0.10 --keep-dup all`, followed by a mean
CPM floor of 2 and 50-400 bp geometry. CPM-refined scores are signal-derived;
they are not q-values or an independent FDR estimate.

### Build an IGV session

Create a portable session containing each ATAC short-fragment CPM track, its
lenient narrowPeak candidates, its refined peaks, and optional matching
H3K27ac CPM tracks:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/build_igv_session.py \
  results/atlas-atac-dm6/baseline/atac_short_fragments \
  --h3k27ac-tracks results/atlas-h3k27ac-dm6/baseline/tracks \
  --output results/atlas-atac-dm6/baseline/atlas.igv.xml \
  --genome dm6
```

The XML uses paths relative to the session file, so move the session and its
track files together.

## Run on SLURM

Workflow rules are executor-independent and expose `threads` and `mem_mb`.
Keep site-specific SBATCH launchers, accounts, partitions, paths, and Snakemake
SLURM profiles under the ignored `slurm/` directory. Never run downloads,
alignment, or peak calling on a cluster login node.

A local shared-filesystem profile can be passed with:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py samples.tsv \
  --workflow-profile slurm/profile \
  --jobs 50 --cores 200 --max-threads 16
```

Here, `--jobs` caps submitted/running jobs, `--cores` caps their aggregate CPU
requests, and `--max-threads` caps any one rule. Cluster-specific values belong
in the ignored profile rather than the portable workflow.

## Standalone acquisition and configuration

The one-command entry point is preferred, but phases can be run separately:

```bash
# One accession or experiment
python src/download_accession.py SRX017289 --output-dir data/raw

# Every accession in a canonical sheet
python src/download_batch.py samples.tsv --output-dir data/raw

# Resolve YAML from a completed manifest
python src/write_pipeline_configs.py samples.tsv \
  --manifest data/raw/download_manifest.tsv \
  --project chromatin-study --run-id baseline
```

Use the repository environment for each command.

## Tests and workflow validation

```bash
export MAMBA_ROOT_PREFIX="$PWD/.micromamba"

mamba run --prefix "$PWD/.venv" pytest -q

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint

mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --cores 8 --dry-run
```

Regenerate the embedded representative ATAC rule graph after changing workflow
dependencies:

```bash
XDG_CACHE_HOME="$PWD/.cache" \
  .venv/bin/snakemake --snakefile workflow/Snakefile \
  --configfile docs/workflow-dag.config.yaml --rulegraph \
  | .venv/bin/dot -Tsvg -o docs/workflow-dag.svg
```

See [`PLAN.md`](PLAN.md) for design decisions and [`AGENTS.md`](AGENTS.md) for
repository-specific contribution rules.
