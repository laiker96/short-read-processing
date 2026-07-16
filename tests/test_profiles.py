from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_slurm_profile_uses_rule_threads_as_cpu_requests():
    profile = yaml.safe_load((REPO_ROOT / "profiles/slurm/config.yaml").read_text())

    assert profile["executor"] == "slurm"
    assert profile["jobs"] > 0
    assert "cores" not in profile
    assert profile["default-resources"]["runtime"] > 0
