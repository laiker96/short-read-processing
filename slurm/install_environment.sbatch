#!/usr/bin/env bash
#SBATCH --job-name=srp-install
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/slurm/%x-%j.out

set -euo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]] || [[ "$(hostname -s)" == cranex* ]]; then
    echo "Refusing to install outside a SLURM compute-node allocation." >&2
    exit 97
fi

REPO_ROOT="${REPO_ROOT:-$HOME/short-read-processing}"
MICROMAMBA="$HOME/.local/bin/micromamba"
BUNDLED_MICROMAMBA="$REPO_ROOT/.cluster-bootstrap/tool/bin/micromamba"
PACKED_CONDA_ENVS="$REPO_ROOT/.cluster-bootstrap/conda-envs"
MAMBA_ROOT_PREFIX="$REPO_ROOT/.cluster-bootstrap/root"
export MAMBA_ROOT_PREFIX

cd "$REPO_ROOT"
mkdir -p logs/slurm "$HOME/.local/bin" "$MAMBA_ROOT_PREFIX"
exec > >(tee -a "logs/install-environment-${SLURM_JOB_ID}.log") 2>&1

echo "[$(date --iso-8601=seconds)] host=$(hostname -s) job=$SLURM_JOB_ID"

if [[ ! -x "$MICROMAMBA" ]] && [[ -x "$BUNDLED_MICROMAMBA" ]]; then
    install -m 0755 "$BUNDLED_MICROMAMBA" "$MICROMAMBA"
elif [[ ! -x "$MICROMAMBA" ]]; then
    archive_dir="$(mktemp -d "${TMPDIR:-/tmp}/micromamba.XXXXXX")"
    trap 'rm -rf "$archive_dir"' EXIT
    curl --fail --location --retry 5 \
        https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C "$archive_dir" bin/micromamba
    install -m 0755 "$archive_dir/bin/micromamba" "$MICROMAMBA"
fi

micromamba_network=()
if [[ -d "$MAMBA_ROOT_PREFIX/pkgs/cache" ]]; then
    micromamba_network=(--offline)
    echo "Using staged offline Micromamba package cache."
fi

if [[ -d "$REPO_ROOT/.venv/conda-meta" ]]; then
    "$MICROMAMBA" install "${micromamba_network[@]}" --yes --prefix "$REPO_ROOT/.venv" \
        --file "$REPO_ROOT/environment.yml" --channel-priority strict
else
    "$MICROMAMBA" create "${micromamba_network[@]}" --yes --prefix "$REPO_ROOT/.venv" \
        --file "$REPO_ROOT/environment.yml" --channel-priority strict
fi

"$MICROMAMBA" run --prefix "$REPO_ROOT/.venv" python --version
"$MICROMAMBA" run --prefix "$REPO_ROOT/.venv" conda --version
"$MICROMAMBA" run --prefix "$REPO_ROOT/.venv" snakemake --version

if compgen -G "$PACKED_CONDA_ENVS/*.tar.gz" > /dev/null; then
    mkdir -p "$REPO_ROOT/.snakemake/conda"
    for archive in "$PACKED_CONDA_ENVS"/*.tar.gz; do
        env_name="$(basename "$archive" .tar.gz)"
        target="$REPO_ROOT/.snakemake/conda/$env_name"
        if [[ -f "$target/.snakemake_timestamp" ]]; then
            echo "Workflow environment already ready: $env_name"
            continue
        fi
        temporary="$REPO_ROOT/.snakemake/conda/.${env_name}.tmp-${SLURM_JOB_ID}"
        rm -rf "$temporary" "$target"
        mkdir -p "$temporary"
        echo "Unpacking workflow environment: $env_name"
        tar -xzf "$archive" -C "$temporary"
        mv "$temporary" "$target"
        env PATH="$target/bin:$PATH" "$target/bin/conda-unpack"
        touch "$target/.snakemake_timestamp"
    done
fi

echo "[$(date --iso-8601=seconds)] environment ready: $REPO_ROOT/.venv"
