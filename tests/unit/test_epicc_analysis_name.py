"""Unit tests for `--analysis-name` on `epicc run` / `epicc validate`.

analysis_name is a flat top-level key, so the flag reaches snakemake through
the single --config block, and the matching --unlock call carries it too.
"""

import ast
import importlib.util
import re
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_EPICC_SRC = _REPO_ROOT / "epicc"


def _load_epicc():
    loader = SourceFileLoader("epicc_cli_analysis_name", str(_EPICC_SRC))
    spec = importlib.util.spec_from_loader("epicc_cli_analysis_name", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli_analysis_name", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


@pytest.fixture
def files(tmp_path):
    samples = tmp_path / "samples.tsv"
    samples.touch()
    options = tmp_path / "opts.yaml"
    options.write_text("analysis_name: fromyaml\ngenomes: {}\n")
    return samples, options


def _parse(subcommand, files, *flags):
    samples, options = files
    return epicc.build_parser().parse_args(
        [subcommand, "--options", str(options), "--samples", str(samples),
         "--cores", "1", *flags])


def _config_values(cmd):
    idx = cmd.index("--config")
    return cmd[idx + 1:]


@pytest.mark.parametrize("subcommand", ["run", "validate"])
def test_override_reaches_config(subcommand, files):
    args = _parse(subcommand, files, "--analysis-name", "myrun")
    cmd = epicc.build_snakemake_cmd(args, [])
    assert "analysis_name=myrun" in _config_values(cmd)


def test_override_on_unlock_cmd(files):
    args = _parse("run", files, "--auto-unlock", "--analysis-name", "myrun")
    cmd = epicc.build_snakemake_cmd(args, [], unlock=True)
    assert "--unlock" in cmd
    assert "analysis_name=myrun" in _config_values(cmd)


def test_auto_unlock_passes_override(files, tmp_path, monkeypatch):
    args = _parse("run", files, "--auto-unlock", "--analysis-name", "myrun")
    monkeypatch.setattr(epicc, "_lock_present", lambda repo: True)
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(epicc.subprocess, "run", fake_run)
    epicc.maybe_clear_stale_lock(args, [], tmp_path)
    assert "--unlock" in calls[0]
    assert "analysis_name=myrun" in _config_values(calls[0])


@pytest.mark.parametrize("subcommand", ["run", "validate"])
def test_absent_when_not_given(subcommand, files):
    args = _parse(subcommand, files)
    assert args.analysis_name is None
    cmd = epicc.build_snakemake_cmd(args, [])
    assert not any(v.startswith("analysis_name=") for v in _config_values(cmd))


@pytest.mark.parametrize("bad", [
    "my__run", "a/b", "my run", "tab\there", "_lead", "trail_", "",
    "x;y", "x$y", "x*y", "x#y",
])
@pytest.mark.parametrize("subcommand", ["run", "validate", "output"])
def test_invalid_name_rejected(subcommand, bad, capsys):
    with pytest.raises(SystemExit) as exc:
        epicc.build_parser().parse_args([subcommand, "--analysis-name", bad])
    assert exc.value.code == 2
    assert "--analysis-name" in capsys.readouterr().err


@pytest.mark.parametrize("good", ["myrun", "Run-2", "a_b", "v1.2", "WT_vs_mut"])
def test_valid_names_accepted(good):
    parse = epicc.build_parser().parse_args
    assert parse(["run", "--analysis-name", good]).analysis_name == good


def test_unsafe_pattern_matches_sample_sheet_validation():
    # The CLI copies the sample-sheet unsafe-char set rather than importing it
    # (samplefile_validation pulls in pandas); keep the two identical.
    src = (_REPO_ROOT / "workflow" / "scripts"
           / "samplefile_validation.py").read_text()
    m = re.search(r"_UNSAFE_CHARS = re\.compile\((r'.*')\)\n", src)
    assert m, "_UNSAFE_CHARS definition not found"
    assert ast.literal_eval(m.group(1)) == epicc._ANALYSIS_NAME_UNSAFE


def test_validate_summary_shows_override(files, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(epicc, "_is_installed_mode", lambda: True)
    monkeypatch.setattr(epicc, "find_repo_root", lambda: _REPO_ROOT)
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(epicc.subprocess, "run", fake_run)
    args = _parse("validate", files, "--analysis-name", "myrun")
    assert epicc.cmd_validate(args, []) == 0
    err = capsys.readouterr().err
    assert "Analysis name: myrun" in err
    assert "fromyaml" not in err
    assert any("analysis_name=myrun" in c for c in calls[0])
