"""Tests for issue #98 parts 2-3: browser extra tracks, --te-file, --ref-genome.

Covers the workflow helpers (normalize_extra_tracks, restrict_to_genome,
check_extra_output_files on browser_extra_tracks) and the epicc CLI plumbing
that turns the new flags into config.
"""

import sys
import types
from importlib.machinery import SourceFileLoader
import importlib.util
from pathlib import Path

import pandas as pd
import pytest
import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "workflow"))

from scripts.browser_extra_tracks import (  # noqa: E402
    EXTRA_TRACK_PALETTE, label_from_path, normalize_extra_tracks,
)
from scripts.sample_sheet import restrict_to_genome  # noqa: E402
from scripts.samplefile_validation import check_extra_output_files  # noqa: E402


def _load_epicc():
    loader = SourceFileLoader("epicc_cli_browser", str(_REPO_ROOT / "epicc"))
    spec = importlib.util.spec_from_loader("epicc_cli_browser", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("epicc_cli_browser", mod)
    loader.exec_module(mod)
    return mod


epicc = _load_epicc()


# ---------------------------------------------------------------------------
# normalize_extra_tracks
# ---------------------------------------------------------------------------

class TestNormalizeExtraTracks:
    def test_empty(self):
        assert normalize_extra_tracks(None) == []
        assert normalize_extra_tracks([]) == []

    def test_bare_paths_get_label_and_palette_colour(self):
        tracks = normalize_extra_tracks(["data/peaks.bed", "x/dmrs.bed.gz"])
        assert tracks == [
            {"file": "data/peaks.bed", "label": "peaks", "color": EXTRA_TRACK_PALETTE[0]},
            {"file": "x/dmrs.bed.gz", "label": "dmrs", "color": EXTRA_TRACK_PALETTE[1]},
        ]

    def test_palette_wraps(self):
        n = len(EXTRA_TRACK_PALETTE) + 1
        tracks = normalize_extra_tracks([f"t{i}.bed" for i in range(n)])
        assert tracks[-1]["color"] == EXTRA_TRACK_PALETTE[0]

    def test_mapping_keeps_given_label_and_colour(self):
        (t,) = normalize_extra_tracks([{"file": "a.bed", "label": "My peaks", "color": "#B2182B"}])
        assert t == {"file": "a.bed", "label": "My peaks", "color": "#B2182B"}

    def test_r_colour_name_accepted(self):
        (t,) = normalize_extra_tracks([{"file": "a.bed", "color": "darkorange"}])
        assert t["color"] == "darkorange"

    @pytest.mark.parametrize("label", ["a\tb", 'say "hi"', "$HOME", "a`b`", "a\\b", "a\nb"])
    def test_unsafe_label_rejected(self, label):
        with pytest.raises(ValueError, match="label"):
            normalize_extra_tracks([{"file": "a.bed", "label": label}])

    @pytest.mark.parametrize("color", ["#12345", "red;rm", "#GGGGGG", "1red"])
    def test_bad_colour_rejected(self, color):
        with pytest.raises(ValueError, match="color"):
            normalize_extra_tracks([{"file": "a.bed", "color": color}])

    def test_missing_file_rejected(self):
        with pytest.raises(ValueError, match="'file' is required"):
            normalize_extra_tracks([{"label": "x"}])

    def test_unknown_key_rejected(self):
        with pytest.raises(ValueError, match="unknown key"):
            normalize_extra_tracks([{"file": "a.bed", "colour": "red"}])

    def test_not_a_list_rejected(self):
        with pytest.raises(ValueError, match="must be a list"):
            normalize_extra_tracks("a.bed")

    def test_label_from_path(self):
        assert label_from_path("/x/y/peaks.narrowPeak.bed.gz") == "peaks"


# ---------------------------------------------------------------------------
# check_extra_output_files: browser_extra_tracks
# ---------------------------------------------------------------------------

class TestExtraTracksValidation:
    def test_valid_bed_passes(self, tmp_path):
        bed = tmp_path / "peaks.bed"
        bed.write_text("chr1\t10\t20\tp1\nchr1\t30\t40\tp2\n")
        check_extra_output_files({"browser_extra_tracks": [str(bed)]})

    def test_missing_file_reported(self, tmp_path):
        with pytest.raises(ValueError, match=r"browser_extra_tracks\[0\].*does not exist"):
            check_extra_output_files({"browser_extra_tracks": [str(tmp_path / "nope.bed")]})

    def test_bad_coordinates_reported(self, tmp_path):
        bed = tmp_path / "bad.bed"
        bed.write_text("chr1\t50\t20\n")
        with pytest.raises(ValueError, match="invalid coordinates"):
            check_extra_output_files({"browser_extra_tracks": [str(bed)]})

    def test_bad_label_reported_without_path_checks(self):
        with pytest.raises(ValueError, match="label"):
            check_extra_output_files(
                {"browser_extra_tracks": [{"file": "a.bed", "label": "$x"}]},
                check_paths=False,
            )


# ---------------------------------------------------------------------------
# restrict_to_genome (only_genome)
# ---------------------------------------------------------------------------

def _sheet():
    return pd.DataFrame({
        "Sample_ID": ["s1", "s1", "s2", "in1", "in1"],
        "Genome": ["B73", "W22", "B73", "B73", "W22"],
    })


class TestRestrictToGenome:
    def test_empty_keeps_all(self):
        assert len(restrict_to_genome(_sheet(), "")) == 5

    def test_keeps_one_genome_and_its_controls(self):
        out = restrict_to_genome(_sheet(), "W22")
        assert list(out["Sample_ID"]) == ["s1", "in1"]
        assert set(out["Genome"]) == {"W22"}
        assert list(out.index) == [0, 1]

    def test_unknown_genome_names_the_present_ones(self):
        with pytest.raises(ValueError, match="matches no sample.*B73, W22"):
            restrict_to_genome(_sheet(), "Mo17")


# ---------------------------------------------------------------------------
# epicc CLI: --ref-genome / --te-file for run and validate
# ---------------------------------------------------------------------------

def _write_options(tmp_path, genomes, **extra):
    cfg = {"analysis_name": "t", "genomes": {g: {"fasta_file": "x.fa", "gff_file": "x.gff"} for g in genomes}}
    cfg.update(extra)
    path = tmp_path / "opts.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def _write_sheet(tmp_path, genome_cells):
    lines = ["Sample_ID\tAssay\tGenome"]
    lines += [f"s{i}\tChIP\t{g}" for i, g in enumerate(genome_cells)]
    lines.append("#parked\tChIP\tCommentedOut")
    path = tmp_path / "samples.tsv"
    path.write_text("\n".join(lines) + "\n")
    return path


def _merged(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


class TestRunGenomeFlags:
    def test_parser_accepts_flags(self, tmp_path):
        te = tmp_path / "te.bed"
        args = epicc.build_parser().parse_args(
            ["run", "--ref-genome", "B73", "--te-file", str(te)]
        )
        assert args.ref_genome == "B73"
        assert args.te_file == str(te.resolve())

    def test_te_file_relative_path_made_absolute(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        args = epicc.build_parser().parse_args(["run", "--te-file", "te.bed"])
        assert args.te_file == str(tmp_path / "te.bed")

    def test_te_file_url_left_alone(self):
        url = "https://example.org/te.gff3.gz"
        args = epicc.build_parser().parse_args(["run", "--te-file", url])
        assert args.te_file == url

    def test_ref_genome_sets_only_genome(self):
        args = types.SimpleNamespace(ref_genome="B73")
        assert "only_genome=B73" in epicc._collect_config_overrides(args)

    def test_sheet_genomes_splits_and_skips_comments(self, tmp_path):
        sheet = _write_sheet(tmp_path, ["B73", "B73,W22"])
        assert epicc._sheet_genomes(sheet) == ["B73", "W22"]

    def test_te_file_uses_ref_genome(self, tmp_path):
        te = tmp_path / "te.bed"
        te.write_text("chr1\t1\t2\n")
        opts = _write_options(tmp_path, ["B73", "W22.v2"])
        args = types.SimpleNamespace(options=str(opts), te_file=str(te), ref_genome="W22.v2", samples=None)
        cfg = _merged(epicc._prepare_options_file(args))
        # Dotted genome names must not be split into nested keys
        # alt_te_file, so the options-file te_file and its outputs stay put
        assert cfg["genomes"]["W22.v2"]["alt_te_file"] == str(te)
        assert "te_file" not in cfg["genomes"]["W22.v2"]
        assert cfg["genomes"]["W22.v2"]["fasta_file"] == "x.fa"
        assert "alt_te_file" not in cfg["genomes"]["B73"]

    def test_te_file_infers_single_genome_from_sheet(self, tmp_path):
        te = tmp_path / "te.bed"
        te.write_text("chr1\t1\t2\n")
        opts = _write_options(tmp_path, ["B73"])
        sheet = _write_sheet(tmp_path, ["B73", "B73"])
        args = types.SimpleNamespace(options=str(opts), te_file=str(te), ref_genome=None, samples=str(sheet))
        assert _merged(epicc._prepare_options_file(args))["genomes"]["B73"]["alt_te_file"] == str(te)

    def test_te_file_multi_genome_sheet_needs_ref_genome(self, tmp_path):
        te = tmp_path / "te.bed"
        te.write_text("chr1\t1\t2\n")
        opts = _write_options(tmp_path, ["B73", "W22"])
        sheet = _write_sheet(tmp_path, ["B73,W22"])
        args = types.SimpleNamespace(options=str(opts), te_file=str(te), ref_genome=None, samples=str(sheet))
        with pytest.raises(SystemExit, match="Pass --ref-genome"):
            epicc._prepare_options_file(args)

    def test_te_file_unknown_genome(self, tmp_path):
        te = tmp_path / "te.bed"
        te.write_text("chr1\t1\t2\n")
        opts = _write_options(tmp_path, ["B73"])
        args = types.SimpleNamespace(options=str(opts), te_file=str(te), ref_genome="Mo17", samples=None)
        with pytest.raises(SystemExit, match="'Mo17' has no entry"):
            epicc._prepare_options_file(args)

    def test_te_file_missing(self, tmp_path):
        opts = _write_options(tmp_path, ["B73"])
        args = types.SimpleNamespace(options=str(opts), te_file=str(tmp_path / "nope.bed"), ref_genome="B73", samples=None)
        with pytest.raises(SystemExit, match="--te-file not found"):
            epicc._prepare_options_file(args)

    def test_no_te_file_returns_options_unchanged(self, tmp_path):
        opts = _write_options(tmp_path, ["B73"])
        args = types.SimpleNamespace(options=str(opts), te_file=None, ref_genome="B73", samples=None)
        assert epicc._prepare_options_file(args) == str(opts)


# ---------------------------------------------------------------------------
# epicc output: --extra-tracks[-labels|-colors]
# ---------------------------------------------------------------------------

def _output_args(tmp_path, *extra):
    return epicc.build_parser().parse_args([
        "output", "--plot-type", "browser", "--input-file", "loci.bed",
        "--ref-genome", "B73", "--env", "all", *extra,
    ])


class TestOutputExtraTracks:
    def test_parse_and_entries(self, tmp_path):
        a, b = tmp_path / "a.bed", tmp_path / "b.bed.gz"
        a.write_text("")
        b.write_text("")
        args = _output_args(
            tmp_path, "--extra-tracks", f"{a},{b}",
            "--extra-tracks-labels", "Peaks,", "--extra-tracks-colors", "red,#00FF00",
        )
        entries = epicc._extra_track_entries(args)
        assert entries == [
            {"file": str(a), "label": "Peaks", "color": "red"},
            {"file": str(b), "label": "b", "color": "#00FF00"},
        ]

    def test_entries_reach_merged_options(self, tmp_path):
        a = tmp_path / "a.bed"
        a.write_text("")
        opts = _write_options(tmp_path, ["B73"], browser_extra_tracks=["old.bed"])
        args = _output_args(tmp_path, "--options", str(opts), "--extra-tracks", str(a))
        cfg = _merged(epicc._prepare_options_file(args))
        assert cfg["browser_extra_tracks"] == [
            {"file": str(a), "label": "a", "color": EXTRA_TRACK_PALETTE[0]}
        ]

    def test_label_count_mismatch(self, tmp_path):
        a = tmp_path / "a.bed"
        a.write_text("")
        args = _output_args(tmp_path, "--extra-tracks", str(a), "--extra-tracks-labels", "x,y")
        with pytest.raises(SystemExit, match="2 value"):
            epicc._extra_track_entries(args)

    def test_missing_file(self, tmp_path):
        args = _output_args(tmp_path, "--extra-tracks", str(tmp_path / "nope.bed"))
        with pytest.raises(SystemExit, match="not found"):
            epicc._extra_track_entries(args)

    def test_bad_colour_rejected_before_snakemake(self, tmp_path):
        a = tmp_path / "a.bed"
        a.write_text("")
        args = _output_args(tmp_path, "--extra-tracks", str(a), "--extra-tracks-colors", "not a colour")
        with pytest.raises(SystemExit, match="color"):
            epicc._extra_track_entries(args)

    def test_extra_tracks_rejected_for_other_plot_types(self, tmp_path, monkeypatch):
        args = epicc.build_parser().parse_args([
            "output", "--plot-type", "heatmap", "--input-file", "loci.bed",
            "--extra-tracks", "a.bed",
        ])
        monkeypatch.setattr(epicc.os, "chdir", lambda p: None)
        with pytest.raises(SystemExit, match="only apply to --plot-type browser"):
            epicc.cmd_output(args, [])

    def test_labels_without_tracks_rejected(self, tmp_path, monkeypatch):
        args = _output_args(tmp_path, "--extra-tracks-labels", "x")
        monkeypatch.setattr(epicc.os, "chdir", lambda p: None)
        with pytest.raises(SystemExit, match="need --extra-tracks"):
            epicc.cmd_output(args, [])


# ---------------------------------------------------------------------------
# prep_browser_on_region: the extra-track shell block, taken from the .smk
# ---------------------------------------------------------------------------

import re  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import gzip  # noqa: E402

_SMK = _REPO_ROOT / "workflow" / "rules" / "combined_analysis.smk"


def _render_extras_block(trackfolder, extras_out, locus, extra_tracks):
    """Cut the extra-track block out of prep_browser_on_region, expanded like Snakemake."""
    text = _SMK.read_text()
    start = text.index("rule prep_browser_on_region:")
    block = text[start:text.index("rule make_single_loci_browser_plot:", start)]
    body = re.search(r'shell:\s*"""(.*?)"""', block, re.DOTALL).group(1)
    body = body[body.index("### Extra user BED tracks"):]
    body = body[: body.index('done <<< "{params.extra_tracks}"') + len('done <<< "{params.extra_tracks}"')]
    for placeholder, value in (
        ("{params.trackfolder}", str(trackfolder)),
        ("{output.extras}", str(extras_out)),
        ("{output.templocus}", str(locus)),
        ("{params.extra_tracks}", extra_tracks),
    ):
        body = body.replace(placeholder, value)
    return "set -euo pipefail\n" + body.replace("{{", "{").replace("}}", "}")


@pytest.mark.skipif(shutil.which("bedtools") is None, reason="bedtools missing")
class TestExtraTracksShell:
    def _run(self, tmp_path, tracks):
        folder = tmp_path / "tracks"
        folder.mkdir()
        locus = tmp_path / "locus.bed"
        locus.write_text("chr1\t1000\t2000\n")
        extras = tmp_path / "extras.txt"
        spec = "\n".join(f"{label}\t{color}\t{path}" for label, color, path in tracks)
        proc = subprocess.run(
            ["bash", "-c", _render_extras_block(folder, extras, locus, spec)],
            capture_output=True, text=True,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        rows = [l.split("\t") for l in extras.read_text().splitlines()]
        return rows, folder

    def test_no_tracks_writes_header_only(self, tmp_path):
        rows, folder = self._run(tmp_path, [])
        assert rows == [["Label", "Color", "Path"]]
        assert list(folder.iterdir()) == []

    def test_plain_and_gzipped_tracks_cut_to_locus(self, tmp_path):
        plain = tmp_path / "peaks.bed"
        plain.write_text(
            "track name=peaks\n"
            "#comment\n"
            "chrom\tstart\tend\tname\n"            # column header
            "chr1\t900\t1100\tedge\tNA\t+\n"       # overlaps; non-numeric score dropped
            "chr1\t1500\t1600\tinside\n"
            "chr1\t5000\t6000\toutside\n"
            "chr2\t1500\t1600\tother_chrom\n"
        )
        gz = tmp_path / "dmrs.bed.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write("chr1\t1200\t1300\n")          # 3-column BED
        rows, folder = self._run(tmp_path, [
            ("Peaks", "#E78AC3", str(plain)),
            ("My DMRs", "darkred", str(gz)),
        ])
        assert rows[0] == ["Label", "Color", "Path"]
        assert rows[1] == ["Peaks", "#E78AC3", str(folder / "extra_1.bed")]
        assert rows[2] == ["My DMRs", "darkred", str(folder / "extra_2.bed")]
        # Clipped to the locus, so a feature hanging over the edge can't widen
        # the plotted range
        assert (folder / "extra_1.bed").read_text() == (
            "chr1\t1000\t1100\tedge\nchr1\t1500\t1600\tinside\n"
        )
        assert (folder / "extra_2.bed").read_text() == "chr1\t1200\t1300\t.\n"

    def test_feature_spanning_locus_is_clipped(self, tmp_path):
        bed = tmp_path / "wide.bed"
        bed.write_text("chr1\t10\t900000\twide\n")
        _, folder = self._run(tmp_path, [("Wide", "red", str(bed))])
        assert (folder / "extra_1.bed").read_text() == "chr1\t1000\t2000\twide\n"

    def test_track_with_nothing_in_locus_is_kept_empty(self, tmp_path):
        bed = tmp_path / "far.bed"
        bed.write_text("chr1\t9000\t9100\n")
        rows, folder = self._run(tmp_path, [("Far", "red", str(bed))])
        assert len(rows) == 2
        assert (folder / "extra_1.bed").read_text() == ""
