from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_local_profile_keeps_environments_local_and_restarts_incomplete_jobs():
    profile = yaml.safe_load((REPO_ROOT / "profiles/local/config.yaml").read_text())

    assert profile["conda-prefix"] == ".snakemake/conda"
    assert profile["rerun-incomplete"] is True
