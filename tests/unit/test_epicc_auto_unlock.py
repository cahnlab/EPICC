"""Unit tests for `epicc run --auto-unlock`.

snakemake parses the Snakefile before it unlocks, and the Snakefile reads the
sample sheet, so the unlock call has to carry the same config as the run.
"""

import importlib.util
import subprocess
import sys
import types
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EPICC_SRC = _REPO_ROOT / "epicc"


def _load_epicc():
    loader = SourceFileLoader("epicc_cli_unlock", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli_unlock", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli_unlock", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


@pytest.fixture
def args(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.touch()
    options = tmp_path / "opts.yaml"
    options.touch()
    return types.SimpleNamespace(
        samples=str(samples), output_dir=None, genome_dir=None,
        keep_intermediates=None, use_node_tmpdir=False, options=str(options),
        profile=None, cores=1, verbose=False, quiet=False,
        workflow_profile=None, no_rerun_incomplete=True, auto_unlock=True,
    )


def _config_values(cmd):
    idx = cmd.index("--config")
    return cmd[idx + 1:]


def test_unlock_cmd_carries_sample_file_override(args):
    cmd = epicc.build_snakemake_cmd(args, extra_args=[], unlock=True)
    assert "--unlock" in cmd
    assert f"sample_file={args.samples}" in _config_values(cmd)
    # --config is greedy, so --unlock must come before it.
    assert cmd.index("--unlock") < cmd.index("--config")


def test_unlock_cmd_keeps_passthrough_config(args):
    cmd = epicc.build_snakemake_cmd(
        args, extra_args=["--config", "te_analysis=false"], unlock=True)
    values = _config_values(cmd)
    assert f"sample_file={args.samples}" in values
    assert "te_analysis=false" in values


def test_failed_unlock_stops_the_run(args, tmp_path, monkeypatch):
    monkeypatch.setattr(epicc, "_lock_present", lambda repo: True)
    monkeypatch.setattr(epicc, "_is_installed_mode", lambda: False)
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 1)

    monkeypatch.setattr(epicc.subprocess, "run", fake_run)
    with pytest.raises(SystemExit) as exc:
        epicc.maybe_clear_stale_lock(args, [], tmp_path)
    assert "--auto-unlock failed" in str(exc.value)
    assert "--unlock" in calls[0]


def test_no_unlock_without_flag(args, tmp_path, monkeypatch):
    args.auto_unlock = False
    monkeypatch.setattr(epicc, "_lock_present", lambda repo: True)
    monkeypatch.setattr(epicc.subprocess, "run",
                        lambda *a, **k: pytest.fail("unlock ran"))
    epicc.maybe_clear_stale_lock(args, [], tmp_path)
