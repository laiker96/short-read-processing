#!/usr/bin/env bash
#SBATCH --job-name=atlas-h3k27ac
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --output=logs/slurm/%x-%j.out

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]] || [[ "$(hostname -s)" == cranex* ]]; then
    echo "Refusing to download or process outside a SLURM compute-node allocation." >&2
    exit 97
fi

REPO_ROOT="${REPO_ROOT:-$HOME/short-read-processing}"
MICROMAMBA="$HOME/.local/bin/micromamba"
SAMPLE_SHEET="$REPO_ROOT/resources/atlas_h3k27ac_ip_only.sample_sheet.tsv"
RAW_ROOT="$REPO_ROOT/data/raw/atlas_h3k27ac"
MANIFEST="$RAW_ROOT/download_manifest.tsv"
CORES="${SLURM_CPUS_PER_TASK:-1}"
MEM_MB="${SLURM_MEM_PER_NODE:-24000}"
DOWNLOAD_ATTEMPTS="${DOWNLOAD_ATTEMPTS:-4}"
DOWNLOAD_FILE_JOBS="${DOWNLOAD_FILE_JOBS:-8}"
DOWNLOAD_CONNECTIONS="${DOWNLOAD_CONNECTIONS:-1}"
NOFILE_LIMIT="${NOFILE_LIMIT:-65536}"

if (( MEM_MB > 2048 )); then
    MEM_MB=$((MEM_MB - 1024))
fi

cd "$REPO_ROOT"
mkdir -p logs/slurm "$RAW_ROOT" configs/cluster .cache
exec > >(tee -a "logs/atlas-h3k27ac-${SLURM_JOB_ID}.log") 2>&1

exec 9>"$REPO_ROOT/.atlas-h3k27ac.lock"
if ! flock -n 9; then
    echo "Another atlas H3K27ac job holds $REPO_ROOT/.atlas-h3k27ac.lock" >&2
    exit 75
fi

hard_nofile="$(ulimit -Hn)"
target_nofile="$NOFILE_LIMIT"
if [[ "$hard_nofile" != "unlimited" ]] && (( target_nofile > hard_nofile )); then
    target_nofile="$hard_nofile"
fi
ulimit -Sn "$target_nofile"
soft_nofile="$(ulimit -Sn)"
if [[ "$soft_nofile" != "unlimited" ]] && (( soft_nofile < 4096 )); then
    echo "Open-file limit $soft_nofile is too low for dm6 processing." >&2
    exit 74
fi
echo "open_file_limit_soft=$soft_nofile open_file_limit_hard=$hard_nofile"

JOB_TMP="${TMPDIR:-/tmp}/short-read-processing-$SLURM_JOB_ID"
mkdir -p "$JOB_TMP"
trap 'rm -rf "$JOB_TMP"' EXIT
export TMPDIR="$JOB_TMP"
export XDG_CACHE_HOME="$REPO_ROOT/.cache"
export MAMBA_ROOT_PREFIX="$REPO_ROOT/.cluster-bootstrap/root"

echo "[$(date --iso-8601=seconds)] host=$(hostname -s) job=$SLURM_JOB_ID"
echo "repo=$REPO_ROOT cores=$CORES schedulable_mem_mb=$MEM_MB tmp=$TMPDIR"
echo "download_file_jobs=$DOWNLOAD_FILE_JOBS connections_per_file=$DOWNLOAD_CONNECTIONS"

for host in www.ebi.ac.uk ftp.sra.ebi.ac.uk; do
    if ! getent ahostsv4 "$host" >/dev/null; then
        echo "Compute node cannot resolve $host; direct ENA download is unavailable." >&2
        exit 69
    fi
done

echo "Ensuring the repository and packed workflow environments are complete."
bash "$REPO_ROOT/slurm/install_environment.sbatch"
if [[ ! -x "$MICROMAMBA" ]] || [[ ! -d "$REPO_ROOT/.venv/conda-meta" ]]; then
    echo "Repository environment installation did not complete." >&2
    exit 2
fi

run_in_env=("$MICROMAMBA" run --prefix "$REPO_ROOT/.venv")
download_command=(
    "${run_in_env[@]}"
    python "$REPO_ROOT/src/run_pipeline.py" "$SAMPLE_SHEET"
    --project atlas-h3k27ac-dm6
    --run-id macs3-broad-ip-only
    --config-dir "$REPO_ROOT/configs/cluster"
    --reference-root "$REPO_ROOT/references"
    --output-dir "$RAW_ROOT"
    --manifest "$MANIFEST"
    --backend ena
    --resolve-jobs 16
    --file-jobs "$DOWNLOAD_FILE_JOBS"
    --connections "$DOWNLOAD_CONNECTIONS"
    --threads "$CORES"
    --download-only
)

download_complete=false
for ((attempt = 1; attempt <= DOWNLOAD_ATTEMPTS; attempt++)); do
    echo "[$(date --iso-8601=seconds)] ENA download attempt $attempt/$DOWNLOAD_ATTEMPTS"
    if "${download_command[@]}"; then
        download_complete=true
        break
    fi
    if (( attempt < DOWNLOAD_ATTEMPTS )); then
        delay=$((attempt * 60))
        echo "Download attempt failed; retaining aria2 partial state and retrying in ${delay}s."
        sleep "$delay"
    fi
done
if [[ "$download_complete" != true ]]; then
    echo "Direct download did not complete after $DOWNLOAD_ATTEMPTS attempts." >&2
    echo "Resubmit this same SBATCH script to resume safely." >&2
    exit 3
fi

echo "[$(date --iso-8601=seconds)] running IP-only H3K27ac broad MACS3 workflow"
"${run_in_env[@]}" \
    python "$REPO_ROOT/src/run_pipeline.py" "$SAMPLE_SHEET" \
    --project atlas-h3k27ac-dm6 \
    --run-id macs3-broad-ip-only \
    --config-dir "$REPO_ROOT/configs/cluster" \
    --reference-root "$REPO_ROOT/references" \
    --output-dir "$RAW_ROOT" \
    --skip-download \
    --manifest "$MANIFEST" \
    --workflow-profile "$REPO_ROOT/profiles/local" \
    --cores "$CORES" \
    --max-threads "$CORES" \
    --snakemake-arg=--resources \
    --snakemake-arg="mem_mb=$MEM_MB" \
    --snakemake-arg=--retries \
    --snakemake-arg=2

echo "[$(date --iso-8601=seconds)] IP-only H3K27ac workflow completed"
