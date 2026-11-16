"""`epicc output` heatmaps and metaplots: one methylation context per plot.

The pipeline draws each context (mCG, mCHG, mCHH) in a plot of its own, so
`--env mC` for a heatmap or metaplot asks for one target per configured
context. The browser keeps every context, one track each, in one plot.
"""

import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest
import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_epicc():
    loader = SourceFileLoader("epicc_cli_output_mc", str(_REPO_ROOT / "epicc"))
    spec = importlib.util.spec_from_loader("epicc_cli_output_mc", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli_output_mc", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


class TestOutputEnvs:
    def test_mc_expands_per_context(self):
        assert epicc.output_envs("heatmap", "mC", ["CG", "CHG", "CHH"]) == ["mCG", "mCHG", "mCHH"]
        assert epicc.output_envs("metaplot", "mC", ["CG"]) == ["mCG"]

    def test_browser_keeps_mc(self):
        assert epicc.output_envs("browser", "mC", ["CG", "CHG", "CHH"]) == ["mC"]

    def test_single_context_and_other_envs_unchanged(self):
        assert epicc.output_envs("heatmap", "mCHH", ["CG", "CHG", "CHH"]) == ["mCHH"]
        assert epicc.output_envs("heatmap", "most", ["CG"]) == ["most"]

    def test_unconfigured_context_rejected(self):
        with pytest.raises(SystemExit, match="CHH.*not in methylation_contexts"):
            epicc.output_envs("metaplot", "mCHH", ["CG"])


@pytest.mark.parametrize("plot_type,prefix", [("heatmap", "Heatmap"), ("metaplot", "Profile")])
def test_cmd_output_targets_one_plot_per_context(tmp_path, monkeypatch, plot_type, prefix):
    monkeypatch.chdir(tmp_path)
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t10\t20\tg1\t0\t+\n")
    opts = tmp_path / "opts.yaml"
    opts.write_text(yaml.safe_dump({"analysis_name": "run", "output_dir": "res",
                                    "methylation_contexts": ["CG", "CHH"]}))
    args = epicc.build_parser().parse_args([
        "output", "--plot-type", plot_type, "--input-file", str(bed),
        "--ref-genome", "ColCEN", "--env", "mC", "--matrix", "tss",
        "--options", str(opts), "--cores", "1",
    ])
    calls = []
    monkeypatch.setattr(epicc, "find_repo_root", lambda: tmp_path)
    monkeypatch.setattr(epicc.subprocess, "run",
                        lambda cmd, *a, **k: calls.append(cmd) or type("R", (), {"returncode": 0})())
    assert epicc.cmd_output(args, []) == 0
    cmd = calls[-1]
    targets = cmd[cmd.index("--") + 1:]
    assert targets == [
        f"res/combined/plots/{prefix}__tss__mCG__run__ColCEN__genes.pdf",
        f"res/combined/plots/{prefix}__tss__mCHH__run__ColCEN__genes.pdf",
    ]
