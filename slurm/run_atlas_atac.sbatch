#!/usr/bin/env bash
#SBATCH --job-name=atlas-atac
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/slurm/%x-%j.out

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]] || [[ "$(hostname -s)" == cranex* ]]; then
    echo "Refusing to process data outside a SLURM compute-node allocation." >&2
    exit 97
fi

REPO_ROOT="${REPO_ROOT:-$HOME/short-read-processing}"
MICROMAMBA="$HOME/.local/bin/micromamba"
SAMPLE_SHEET="$REPO_ROOT/resources/atlas_atac_selected.sample_sheet.tsv"
MANIFEST="$REPO_ROOT/data/raw/atlas_atac/download_manifest.tsv"
CORES="${SLURM_CPUS_PER_TASK:-1}"
MEM_MB="${SLURM_MEM_PER_NODE:-24000}"
if (( MEM_MB > 2048 )); then
    MEM_MB=$((MEM_MB - 1024))
fi

cd "$REPO_ROOT"
mkdir -p logs/slurm
exec > >(tee -a "logs/atlas-atac-${SLURM_JOB_ID}.log") 2>&1

if [[ ! -x "$MICROMAMBA" ]] || [[ ! -d "$REPO_ROOT/.venv/conda-meta" ]]; then
    echo "Missing repo environment; complete install_environment.sbatch first." >&2
    exit 2
fi

JOB_TMP="${TMPDIR:-/tmp}/short-read-processing-$SLURM_JOB_ID"
mkdir -p "$JOB_TMP"
trap 'rm -rf "$JOB_TMP"' EXIT
export TMPDIR="$JOB_TMP"
export XDG_CACHE_HOME="$REPO_ROOT/.cache"
export MAMBA_ROOT_PREFIX="$REPO_ROOT/.cluster-bootstrap/root"

echo "[$(date --iso-8601=seconds)] host=$(hostname -s) job=$SLURM_JOB_ID"
echo "repo=$REPO_ROOT cores=$CORES schedulable_mem_mb=$MEM_MB tmp=$TMPDIR"

"$MICROMAMBA" run --prefix "$REPO_ROOT/.venv" \
    python "$REPO_ROOT/src/run_pipeline.py" "$SAMPLE_SHEET" \
    --project atlas-atac-dm6 \
    --run-id hmmratac \
    --config-dir "$REPO_ROOT/configs/cluster" \
    --reference-root "$REPO_ROOT/references" \
    --output-dir "$REPO_ROOT/data/raw/atlas_atac" \
    --skip-download \
    --manifest "$MANIFEST" \
    --workflow-profile "$REPO_ROOT/profiles/local" \
    --cores "$CORES" \
    --max-threads "$CORES" \
    --snakemake-arg=--resources \
    --snakemake-arg="mem_mb=$MEM_MB"

echo "[$(date --iso-8601=seconds)] atlas ATAC pipeline completed"
