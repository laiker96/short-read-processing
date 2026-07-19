#!/usr/bin/env bash
set -euo pipefail

if [[ "$(hostname -s)" == cranex* ]]; then
    echo "Refusing to stage or orchestrate from the cranex login node." >&2
    exit 97
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOTE="${REMOTE:-ilaiker@cluster.qb.fcen.uba.ar}"
REMOTE_ROOT="${REMOTE_ROOT:-short-read-processing}"
KNOWN_HOSTS="${KNOWN_HOSTS:-/tmp/qb_cluster_known_hosts}"
MANIFEST="$REPO_ROOT/data/raw/atlas_atac/download_manifest.tsv"
SSH_COMMAND=(
    ssh
    -o BatchMode=yes
    -o ConnectTimeout=20
    -o StrictHostKeyChecking=accept-new
    -o "UserKnownHostsFile=$KNOWN_HOSTS"
)
RSYNC_SSH="${SSH_COMMAND[*]}"

echo "Waiting for the checksum-verified portable manifest: $MANIFEST"
while [[ ! -s "$MANIFEST" ]]; do
    sleep 30
done
while tmux has-session -t cluster_stage_static 2>/dev/null; do
    sleep 30
done

cd "$REPO_ROOT"
export PATH="$REPO_ROOT/.venv/bin:$PATH"
export XDG_CACHE_HOME="$REPO_ROOT/.cache"

python src/validate_sample_sheet.py resources/atlas_atac_selected.sample_sheet.tsv
python src/run_pipeline.py resources/atlas_atac_selected.sample_sheet.tsv \
    --project atlas-atac-dm6 \
    --run-id hmmratac \
    --config-dir configs/cluster \
    --reference-root references \
    --output-dir data/raw/atlas_atac \
    --skip-download \
    --manifest "$MANIFEST" \
    --workflow-profile profiles/local \
    --cores 24 \
    --max-threads 24 \
    --snakemake-dry-run \
    --snakemake-arg=--resources \
    --snakemake-arg=mem_mb=23000

rsync -az --partial --info=progress2 \
    --exclude=.git/ \
    --exclude=.venv/ \
    --exclude=.snakemake/ \
    --exclude=.cluster-bootstrap/ \
    --exclude=.cache/ \
    --exclude=data/ \
    --exclude=references/ \
    --exclude=results/ \
    --exclude=work/ \
    --exclude=logs/ \
    --exclude=tests/workflow_results/ \
    -e "$RSYNC_SSH" \
    ./ "$REMOTE:$REMOTE_ROOT/"

rsync -aR --partial --info=progress2 \
    -e "$RSYNC_SSH" \
    .cluster-bootstrap/tool/bin/micromamba \
    .cluster-bootstrap/root/pkgs \
    .cluster-bootstrap/conda-envs \
    references/dm6 \
    "$REMOTE:$REMOTE_ROOT/"

rsync -aR --partial --info=progress2 \
    -e "$RSYNC_SSH" \
    data/raw/atlas_atac \
    "$REMOTE:$REMOTE_ROOT/"

node_table="$("${SSH_COMMAND[@]}" "$REMOTE" \
    "sinfo -N -h -p cpu -t idle -o '%N|%c|%m' | sort -u")"
echo "Idle CPU nodes:"
printf '%s\n' "$node_table"
node="$(awk -F '|' '$2 >= 24 && $3 >= 24576 {print $1; exit}' <<< "$node_table")"
if [[ -z "$node" ]]; then
    echo "No idle CPU node with at least 24 CPUs and 24 GiB RAM was found." >&2
    exit 3
fi

submission="$("${SSH_COMMAND[@]}" "$REMOTE" \
    "cd '$REMOTE_ROOT' && mkdir -p logs/slurm && env_job=\$(sbatch --parsable --nodelist='$node' slurm/install_environment.sbatch) && run_job=\$(sbatch --parsable --dependency=afterok:\$env_job --nodelist='$node' slurm/run_atlas_atac.sbatch) && printf '%s|%s\n' \"\$env_job\" \"\$run_job\"")"
env_job="${submission%%|*}"
run_job="${submission##*|}"
echo "Submitted node=$node environment_job=$env_job pipeline_job=$run_job"
"${SSH_COMMAND[@]}" "$REMOTE" \
    "squeue -j '$env_job,$run_job' -o '%.18i %.12j %.9T %.10M %.6D %R'"
