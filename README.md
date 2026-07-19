# short-read-processing

Snakemake workflow for reproducible bulk ATAC-seq and ChIP-seq short-read
processing. Its public input is a CSV/TSV sample sheet containing
SRR/SRX/ERR/ERX accessions and typed experiment parameters. The pipeline caches
FASTQs internally, supports `dm6` and `hg38`, calls peaks with MACS3 `callpeak`
or `hmmratac`, and produces signal tracks and assay-aware QC summaries.

The workflow runs raw and trimmed FastQC, Cutadapt, Bowtie2, duplicate marking,
MAPQ/blacklist/mitochondrial filtering, MACS3 peak calling, BigWig generation,
FRiP and assay-aware QC, and MultiQC. See [PLAN.md](PLAN.md) for design details.

## Local environment

Create the orchestration and data-acquisition environment from the repository
root. An absolute prefix is used deliberately so Mamba cannot interpret
`.venv` as a named environment under its global `envs` directory.

```bash
mamba env create --prefix "$PWD/.venv" --file environment.yml
mamba run --prefix "$PWD/.venv" snakemake --version
```

Bioinformatics tools are defined separately under `workflow/envs/`. Snakemake's
local profile sets `conda-prefix: .snakemake/conda`, so all rule environments
are also created below the repository root when the workflow is run from here:

```bash
mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile configs/project.yaml \
  --workflow-profile profiles/local
```

Both `.venv` and `.snakemake` are intentionally ignored by Git. The R-based
ChIP cross-correlation dependencies are isolated in
`workflow/envs/chip_qc.yaml`; they are not installed into `.venv`.

## Run accession to peaks and tracks

The normal entry point performs validation, concurrent download, resolved YAML
generation, and Snakemake execution in one command:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py samples.tsv \
  --project chromatin-study \
  --run-id baseline \
  --reference-root references \
  --output-dir data/raw
```

Use `--skip-download --manifest data/raw/download_manifest.tsv` to reuse an
existing checksum-verified download. `--download-only`, `--config-only`, and
`--snakemake-dry-run` stop at the corresponding boundary. Samples and lanes are
scheduled independently subject to the local core and memory profile.

Each technical lane is trimmed and aligned independently with a four-core job.
The lane BAMs are coordinate-sorted, merged by biological library, and only then
duplicate-marked. Thus, a 12-core run can align three lanes concurrently without
assigning the same cores simultaneously to Bowtie2 and SAMtools sort.

## Run on SLURM

The generic `profiles/slurm/config.yaml` uses Snakemake's SLURM executor on a
shared filesystem. Every rule's `threads` value becomes its per-job CPU request;
SLURM determines which node runs it, while Snakemake submits independent lanes
and samples concurrently. The profile limits queued/running workflow jobs to 100
but deliberately does not impose an aggregate cluster-wide core limit.

Run with the cluster profile and optionally tune the submission limits and
per-rule CPUs from the launcher:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py samples.tsv \
  --project chromatin-study \
  --run-id baseline \
  --workflow-profile profiles/slurm \
  --jobs 50 \
  --cores 200 \
  --max-threads 16 \
  --set-threads align_lane=8 \
  --set-threads build_bowtie2_index=16
```

Here, `--jobs` caps simultaneous SLURM jobs, `--cores` caps their aggregate CPU
requests, and `--max-threads` caps any one job. `--set-threads` overrides a
specific rule and is repeatable. Omit `--cores` to let the job cap alone bound
submission. Add cluster-specific defaults such as an account or partition by
copying the profile and adding them under `default-resources`:

```yaml
default-resources:
  mem_mb: 2000
  runtime: 720
  slurm_account: my_account
  slurm_partition: compute
```

The orchestration environment includes the SLURM executor plugin. Re-run the
repository-local Mamba environment command after pulling this change. The
cluster nodes must see the repository, inputs, results, and `.snakemake/conda`
through the same shared filesystem paths.

`slurm/install_environment.sbatch` installs the orchestration environment at
`$HOME/short-read-processing/.venv`; `slurm/run_atlas_atac.sbatch` processes all
23 selected accessions in `resources/atlas_atac_selected.sample_sheet.tsv`
inside one 24-core compute allocation. Both scripts refuse to run on a host
named `cranex*`. The ATAC job reuses the staged, checksum-verified manifest and
does not perform network downloads. On an offline cluster, stage the Micromamba
package cache under `.cluster-bootstrap/root/pkgs`, packed rule environments
under `.cluster-bootstrap/conda-envs`, FASTQs under `data/raw/atlas_atac`, and
the prepared dm6 assets under `references/dm6` before submission.

`slurm/download_and_run_atlas_atac_dual.sbatch` is the direct-download variant.
Inside one compute allocation it installs the repo-local environment when
needed, downloads and checksum-verifies all selected ENA FASTQs, runs the
primary HMMRATAC workflow, and then runs MACS3 `callpeak` with `-f BAM
--nomodel --shift -75 --extsize 150 -B --SPMR` from the already-produced final
BAMs. MACS3 outputs are isolated below
`results/atlas-atac-dm6/macs3-shift-neg75-extsize-150/`; alignment is not run a
second time. The job is safe to resubmit: aria2 resumes managed partial files,
untracked size-mismatched partial copies are restarted, the manifest is written
atomically, and Snakemake reruns incomplete jobs. It also refuses `cranex*`,
checks compute-node ENA DNS before downloading, and uses a lock to prevent two
jobs from writing the same dataset concurrently. It also raises the per-process
open-file soft limit to 65,536 (or the cluster hard limit) so deepTools can
concatenate temporary ATAC-shift BAM chunks for the 1,870-contig dm6 assembly.
`slurm/download_and_run_atlas_h3k27ac_ip_only.sbatch` applies the same
restart-safe direct-download pattern to the 15 selected H3K27ac IP-only runs,
using broad MACS3 defaults without external input controls.
From a networked workstation, `slurm/stage_and_submit_atlas.sh` waits for the
completed local manifest, verifies a full dry-run, stages those assets, selects
an idle CPU node with at least 24 CPUs and 24 GiB RAM, and submits the install
and dependent processing jobs. Override `REMOTE` or `REMOTE_ROOT` when needed.

For generated `dm6` and `hg38` configs, reference preparation is part of the
Snakemake DAG. The workflow downloads checksum-pinned UCSC FASTA and NCBI
RefSeq GTF archives plus the ENCODE v2 blacklist, then creates:

```text
references/<genome>/<genome>.fa
references/<genome>/<genome>.blacklist.bed
references/<genome>/<genome>.tss.bed
references/<genome>/<genome>.autosomes.txt
references/<genome>/<genome>.fa.fai
references/<genome>/<genome>.chrom.sizes
references/<genome>/bowtie2/<genome>.*.bt2
```

Downloaded source archives are retained under
`references/<genome>/sources/`. Their HTTPS URLs and MD5/SHA-256 checksums are
written into each resolved config and therefore into run provenance. A custom
resolved config can omit `reference.preparation`; in that local-only mode the
listed FASTA, blacklist, TSS BED, and autosome file remain explicit inputs,
while FASTA/chromosome-size and Bowtie2 indexes are still built when missing.

Two real dm6 examples curated from `../ABC-map-dm6/metadata/` are included:

- `resources/abc_map_dm6_stage5_atac.sample_sheet.tsv`: two biological ATAC
  replicates, preserving replicate 1's two technical runs;
- `resources/abc_map_dm6_stage5_h3k27ac.sample_sheet.tsv`: two H3K27ac IP
  replicates with matched replicate-specific input controls.

For example, this starts the real stage-5 ATAC accession-to-peaks run and also
prepares dm6 automatically:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/run_pipeline.py \
  resources/abc_map_dm6_stage5_atac.sample_sheet.tsv \
  --project abc-map-dm6-stage5-atac \
  --run-id baseline \
  --output-dir data/raw
```

Results are written below `results/<project>/<run_id>/`, including:

- `bam/*.final.bam` and indexes;
- `peaks/` with narrow/broad/HMMRATAC peaks and MACS3 treatment/control
  bedGraphs for `callpeak`;
- `tracks/*.CPM.bw`, plus Tn5-shifted ATAC BigWigs;
- raw/trimmed FastQC, Cutadapt, alignment, FRiP, TSS, fragment-size, ChIP
  fingerprint/cross-correlation, metrics TSV/JSON, and MultiQC outputs.

## Download public FASTQs

The download commands accept SRR/SRX and ERR/ERX accessions. Experiment
accessions are expanded to every associated run. By default, the commands use
ENA's already-compressed FASTQ files and submit all files to one `aria2c` queue:
files download concurrently, and each large file uses multiple HTTP range
connections. Downloads resume and ENA MD5 checksums are verified. Runs without
a complete ENA FASTQ set fall back to concurrent `prefetch` + `fasterq-dump` +
`pigz` processing.

Download one run or experiment:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/download_accession.py SRX017289 \
  --output-dir data/raw \
  --file-jobs 8 \
  --connections 8
```

Inspect the resolution and expected transfer size without downloading:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/download_accession.py SRR123456 --dry-run
```

For batch download, use the canonical sample sheet. The complete machine-readable
contract is in `schemas/sample-sheet.schema.yaml`; a runnable template is in
`resources/sample_sheet.example.tsv`. The required columns are:

| column | meaning |
|---|---|
| `accession` | SRR/SRX/ERR/ERX; experiment IDs expand to all runs |
| `sample_id` | biological library ID; repeat it for technical runs |
| `assay` | `atac`, `chip_tf`, or `chip_histone` |
| `genome` | `dm6` or `hg38` |
| `role` | `treatment` or `control` |
| `control_id` | matched ChIP control `sample_id`; blank otherwise |
| `replicate` | positive biological replicate number |
| `peak_caller` | `hmmratac`, `callpeak`, or blank for the assay default |

Validate without network access:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/validate_sample_sheet.py resources/sample_sheet.example.tsv
```

Then download all accessions concurrently:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/download_batch.py samples.tsv \
  --output-dir data/raw
```

The stable `data/raw/download_manifest.tsv` records requested accessions,
expanded runs, layouts, backends, local FASTQs, and checksums. Repeated commands
update matching rows rather than discarding unrelated manifest entries. Raw
files are stored as `data/raw/<run accession>/*.fastq.gz`. FASTQ paths below the
manifest directory are stored relatively and resolved against that directory
when read, so the directory can be moved to another host without path edits.

Force `--backend ena` to disallow conversion fallbacks or `--backend sra` to
force SRA Toolkit. `--file-jobs` controls simultaneous FASTQ files,
`--connections` controls range connections per file, and `--sra-jobs` plus
`--threads` control concurrent SRA conversions. The defaults are intentionally
aggressive but bounded; reduce them on shared filesystems or narrow links.

## Generate pipeline YAML files

Generate one resolved config per assay/genome group from the same sample sheet
and completed internal manifest:

```bash
mamba run --prefix "$PWD/.venv" \
  python src/write_pipeline_configs.py \
  samples.tsv \
  --manifest data/raw/download_manifest.tsv \
  --project chromatin-study \
  --output-dir configs
```

The sheet can override peak and preprocessing parameters with typed optional
columns. For example, this selects MACS3 insertion-site-style ATAC calling
instead of the paired-end HMMRATAC default:

```text
peak_caller=callpeak
macs3_format=BAM
macs3_nomodel=true
macs3_shift=-75
macs3_extsize=150
```

Paired-end ATAC defaults to HMMRATAC. TF and histone ChIP default to MACS3
`callpeak`, with narrow and broad peaks respectively. Single-end ATAC must set
`peak_caller=callpeak`. All `callpeak` configurations enable `-B --SPMR` and
declare both MACS3 bedGraph products: treatment pileup and control lambda.
Other defaults include Nextera adapters for ATAC, TruSeq for ChIP, Bowtie2
`very-sensitive`, MAPQ 30, and duplicate and mitochondrial filtering.

ChIP input libraries use `role=control`; a treatment names its matched
control's `sample_id` in `control_id`. IP-only ChIP is supported by leaving
`control_id` blank; MACS3 then runs without `-c` and control-dependent
fingerprint QC is omitted. Invalid or cross-assay control references fail
sample-sheet validation. Generated reference paths and their pinned preparation
sources are placed below `references/dm6` or `references/hg38`; they do not need
to exist beforehand.

## Tests

Run the unit tests, all assay-branch Snakemake dry-runs, and the loopback aria2
integration check inside the repo-local environment:

```bash
mamba run --prefix "$PWD/.venv" pytest -q
```

Lint and dry-run the workflow directly:

```bash
mamba run --prefix "$PWD/.venv" \
  snakemake --snakefile workflow/Snakefile \
  --configfile tests/fixtures/workflow_config.yaml --lint
```
