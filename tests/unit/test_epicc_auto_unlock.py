"""Unit tests for `epicc run --auto-unlock` and `epicc validate --auto-unlock`.

snakemake parses the Snakefile before it unlocks, and the Snakefile reads the
sample sheet, so the unlock call has to carry the same config as the run.
"""

import argparse
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


# ---------------------------------------------------------------------------
# epicc validate: accepts --auto-unlock, and reports a lock it leaves in place
# ---------------------------------------------------------------------------

def test_validate_accepts_auto_unlock():
    parse = epicc.build_parser().parse_args
    assert parse(["validate", "--auto-unlock"]).auto_unlock is True
    assert parse(["validate"]).auto_unlock is False
    assert parse(["run", "--auto-unlock"]).auto_unlock is True


def _validate_env(tmp_path, monkeypatch, *flags, locked=True):
    """Run cmd_validate in installed mode from tmp_path with snakemake stubbed."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(epicc, "_is_installed_mode", lambda: True)
    monkeypatch.setattr(epicc, "find_repo_root", lambda: _REPO_ROOT)
    if locked:
        locks = tmp_path / ".snakemake" / "locks"
        locks.mkdir(parents=True)
        (locks / "0.input.lock").write_text("")
    samples = tmp_path / "samples.tsv"
    samples.write_text("")
    options = tmp_path / "opts.yaml"
    options.write_text("analysis_name: t\ngenomes: {}\n")
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        if "--unlock" in cmd:
            for f in (tmp_path / ".snakemake" / "locks").glob("*"):
                f.unlink()
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(epicc.subprocess, "run", fake_run)
    args = epicc.build_parser().parse_args(
        ["validate", "--options", str(options), "--samples", str(samples),
         "--cores", "1", *flags])
    rc = epicc.cmd_validate(args, [])
    return rc, calls


def test_validate_auto_unlock_never_unlocks(tmp_path, monkeypatch, capsys):
    # validate runs beside live runs; a lock file can't tell live from dead.
    rc, calls = _validate_env(tmp_path, monkeypatch, "--auto-unlock")
    assert rc == 0
    assert not any("--unlock" in c for c in calls)
    assert any((tmp_path / ".snakemake" / "locks").iterdir())
    err = capsys.readouterr().err
    assert "working directory is locked" in err
    assert "validate never unlocks" in err
    assert "epicc run --auto-unlock" in err


def test_validate_reports_lock_without_flag(tmp_path, monkeypatch, capsys):
    rc, calls = _validate_env(tmp_path, monkeypatch)
    assert rc == 0
    assert not any("--unlock" in c for c in calls)
    err = capsys.readouterr().err
    assert "Validation passed" in err
    assert "working directory is locked" in err
    assert "epicc unlock" in err


def test_validate_silent_when_unlocked(tmp_path, monkeypatch, capsys):
    rc, _ = _validate_env(tmp_path, monkeypatch, locked=False)
    assert rc == 0
    assert "working directory is locked" not in capsys.readouterr().err


def test_validate_verbose_reports_lock(tmp_path, monkeypatch, capsys):
    rc, calls = _validate_env(tmp_path, monkeypatch, "--verbose")
    assert rc == 0
    assert not any("--unlock" in c for c in calls)
    err = capsys.readouterr().err
    assert "Validation passed" in err
    assert "working directory is locked" in err


def test_auto_unlock_help_differs_for_validate():
    def help_for(sub):
        p = epicc.build_parser()
        sp = next(a for a in p._actions
                  if isinstance(a, argparse._SubParsersAction))
        act = next(a for a in sp.choices[sub]._actions if a.dest == "auto_unlock")
        return act.help
    assert "default: halt" in help_for("run")
    assert "never clears a lock" in help_for("validate")
    assert "default: halt" not in help_for("validate")
