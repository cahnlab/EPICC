"""Unit tests for workflow/scripts/samplefile_validation.py

Focuses on the file-existence rules added for Read_files paths and genome
config file paths. SRA accessions and HTTP(S) URLs are skipped — only
local-filesystem paths are probed.
"""

import os
import sys

import pandas as pd
import pytest

# samplefile_validation imports `from scripts.sample_sheet import ...`,
# which only resolves when the parent of the workflow/ tree is on sys.path
# (mirroring how Snakemake adds workflow.basedir's parent at runtime).
_REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "workflow"))

from scripts.samplefile_validation import (
    check_table, check_genome_config, check_extra_output_files,
)


def _row(sample_id, assay, read_files, layout="SE", **overrides):
    base = {
        "Sample_ID": sample_id,
        "Assay": assay,
        "Genome": "test_genome",
        "Levels": "genotype:WT",
        "Replicate_ID": "rep1",
        "Read_files": read_files,
        "Read_layout": layout,
        "IP_target": "",
        "Control": "",
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# Read_files: local path existence
# ---------------------------------------------------------------------------

class TestReadFilesPathExistence:
    def test_existing_local_path_passes(self, tmp_path):
        f = tmp_path / "reads.fq.gz"
        f.write_bytes(b"")
        df = pd.DataFrame([_row("s1", "RNAseq", str(f))])
        check_table(df)  # no raise

    def test_missing_local_path_raises(self, tmp_path):
        missing = tmp_path / "missing.fq.gz"
        df = pd.DataFrame([_row("s1", "RNAseq", str(missing))])
        with pytest.raises(ValueError) as excinfo:
            check_table(df)
        assert "does not exist" in str(excinfo.value)
        assert str(missing) in str(excinfo.value)

    def test_missing_pe_mate_raises(self, tmp_path):
        r1 = tmp_path / "r1.fq.gz"
        r1.write_bytes(b"")
        r2 = tmp_path / "r2.fq.gz"  # never created
        rf = f"{r1},{r2}"
        df = pd.DataFrame([_row("s1", "RNAseq", rf, layout="PE")])
        with pytest.raises(ValueError) as excinfo:
            check_table(df)
        msg = str(excinfo.value)
        assert "r2.fq.gz" in msg
        assert "does not exist" in msg

    def test_missing_merge_component_raises(self, tmp_path):
        a = tmp_path / "a.fq.gz"
        a.write_bytes(b"")
        b = tmp_path / "b.fq.gz"  # never created
        df = pd.DataFrame([_row("s1", "RNAseq", f"{a}+{b}")])
        with pytest.raises(ValueError) as excinfo:
            check_table(df)
        assert "b.fq.gz" in str(excinfo.value)

    def test_sra_accession_not_probed(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "SRR1234567")])
        check_table(df)  # no raise — SRA IDs are never probed on disk

    def test_ena_drr_not_probed(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "DRR400324", layout="PE")])
        check_table(df)

    def test_http_url_not_probed(self):
        df = pd.DataFrame([
            _row("s1", "dmC", "https://example.org/data.bam"),
        ])
        check_table(df)

    def test_check_paths_false_skips_existence(self, tmp_path):
        missing = tmp_path / "missing.fq.gz"
        df = pd.DataFrame([_row("s1", "RNAseq", str(missing))])
        # Bypass: caller signals that inputs may not be staged yet
        check_table(df, check_paths=False)


# ---------------------------------------------------------------------------
# Read_files: '+'-merge type compatibility (SRA + FASTQ mergeable; BAM,
# bedMethyl, and mixed-type merges rejected)
# ---------------------------------------------------------------------------

class TestMergeTypeCompatibility:
    def test_sra_merge_ok(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "SRR111+SRR222")])
        check_table(df, check_paths=False)  # no raise

    def test_fastq_se_merge_ok(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "/d/a.fq.gz+/d/b.fq.gz")])
        check_table(df, check_paths=False)  # no raise

    def test_fastq_pe_merge_ok(self):
        rf = "/d/a_R1.fq.gz,/d/a_R2.fq.gz+/d/b_R1.fq.gz,/d/b_R2.fq.gz"
        df = pd.DataFrame([_row("s1", "RNAseq", rf, layout="PE")])
        check_table(df, check_paths=False)  # no raise

    def test_fastq_url_merge_ok(self):
        rf = "https://h/a.fq.gz+https://h/b.fq.gz"
        df = pd.DataFrame([_row("s1", "RNAseq", rf)])
        check_table(df, check_paths=False)  # no raise

    def test_bam_merge_rejected(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "/d/a.bam+/d/b.bam")])
        with pytest.raises(ValueError) as e:
            check_table(df, check_paths=False)
        assert "not supported for bam" in str(e.value)

    def test_bedmethyl_merge_rejected(self):
        df = pd.DataFrame([_row("s1", "dmC", "/d/a.bed.gz+/d/b.bed.gz")])
        with pytest.raises(ValueError) as e:
            check_table(df, check_paths=False)
        assert "not supported for bedmethyl" in str(e.value)

    def test_mixed_type_merge_rejected(self):
        df = pd.DataFrame([_row("s1", "RNAseq", "SRR111+/d/b.fq.gz")])
        with pytest.raises(ValueError) as e:
            check_table(df, check_paths=False)
        assert "same type" in str(e.value)


# ---------------------------------------------------------------------------
# Genome config: file existence
# ---------------------------------------------------------------------------

def _minimal_sample_df(genome="test_genome", assay="RNAseq"):
    return pd.DataFrame([_row("s1", assay, "SRR1234567")]).assign(Genome=genome)


class TestGenomeConfigPathExistence:
    def _genome_cfg(self, **fields):
        cfg = {"genomes": {"test_genome": fields}}
        return cfg

    def test_missing_fasta_raises(self, tmp_path):
        gff = tmp_path / "annot.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(
            fasta_file=str(tmp_path / "absent.fa"),
            gff_file=str(gff),
        )
        with pytest.raises(ValueError) as excinfo:
            check_genome_config(_minimal_sample_df(), cfg)
        assert "fasta_file" in str(excinfo.value)
        assert "does not exist" in str(excinfo.value)

    def test_existing_fasta_and_gff_pass(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_bytes(b"")
        gff = tmp_path / "g.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(fasta_file=str(fa), gff_file=str(gff))
        check_genome_config(_minimal_sample_df(), cfg)

    def test_auto_sentinel_skipped(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_bytes(b"")
        gff = tmp_path / "g.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(
            fasta_file=str(fa), gff_file=str(gff),
            gtf_file="<auto>",
        )
        check_genome_config(_minimal_sample_df(), cfg)

    def test_url_not_probed_for_existence(self, tmp_path):
        cfg = self._genome_cfg(
            fasta_file="https://example.org/g.fa",
            gff_file="https://example.org/g.gff3",
        )
        check_genome_config(_minimal_sample_df(), cfg)

    def test_te_file_existence(self, tmp_path):
        fa = tmp_path / "g.fa"
        fa.write_bytes(b"")
        gff = tmp_path / "g.gff3"
        gff.write_bytes(b"")
        cfg = self._genome_cfg(
            fasta_file=str(fa), gff_file=str(gff),
            te_file=str(tmp_path / "absent.te.bed"),
        )
        with pytest.raises(ValueError) as excinfo:
            check_genome_config(_minimal_sample_df(), cfg)
        assert "te_file" in str(excinfo.value)

    def test_check_paths_false_skips_existence(self, tmp_path):
        cfg = self._genome_cfg(
            fasta_file=str(tmp_path / "absent.fa"),
            gff_file=str(tmp_path / "absent.gff3"),
        )
        check_genome_config(_minimal_sample_df(), cfg, check_paths=False)


# ---------------------------------------------------------------------------
# Extra output target file format validation
# ---------------------------------------------------------------------------

class TestBrowserTargetFile:
    """Browser target file: chrom/start/end/label/binsize[/htstart/htwidth]."""

    def _write(self, path, rows, header=None):
        with open(path, "w") as fh:
            if header:
                fh.write("\t".join(header) + "\n")
            for r in rows:
                fh.write("\t".join(str(c) for c in r) + "\n")

    def test_minimal_valid_5col(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "promoterX", 50)])
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        check_extra_output_files(cfg)

    def test_with_optional_ht(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "loc1", 50, "150,250", "20,30")])
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        check_extra_output_files(cfg)

    def test_with_header_row(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(
            f, [("chrI", 100, 500, "loc1", 50)],
            header=["chrom", "start", "end", "name", "bs"],
        )
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        check_extra_output_files(cfg)

    def test_label_starting_with_dash_rejected(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "-flag", 50)])
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "must not start with '-'" in str(exc.value)

    def test_binsize_zero_rejected(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "loc1", 0)])
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "binsize must be >= 1" in str(exc.value)

    def test_binsize_non_integer_rejected(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "loc1", "abc")])
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "binsize" in str(exc.value)

    def test_too_few_columns_rejected(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "loc1")])  # missing binsize
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "at least 5" in str(exc.value)

    def test_invalid_coords(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 500, 100, "loc1", 50)])  # end < start
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "invalid coordinates" in str(exc.value)

    def test_ht_pair_required(self, tmp_path):
        f = tmp_path / "browser.bed"
        self._write(f, [("chrI", 100, 500, "loc1", 50, "150", "")])
        cfg = {"full_analysis": True, "browser_target_file": str(f)}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "htstart" in str(exc.value)

    def test_skipped_when_full_analysis_off(self, tmp_path):
        f = tmp_path / "browser.bed"
        # Even malformed rows should slip through when the analysis is disabled.
        self._write(f, [("chrI", 100, 500, "-flag", 0)])
        cfg = {"full_analysis": False, "browser_target_file": str(f)}
        check_extra_output_files(cfg)

    def test_default_placeholder_path_skipped(self):
        # The shipped default value points to a path that doesn't exist;
        # users who don't customize it shouldn't see a startup error.
        cfg = {"full_analysis": True,
               "browser_target_file": "data/target_loci.bed"}
        check_extra_output_files(cfg)

    def test_missing_customized_path_errors(self, tmp_path):
        cfg = {"full_analysis": True,
               "browser_target_file": str(tmp_path / "absent.bed")}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "does not exist" in str(exc.value)


class TestOtherExtraTargetFiles:
    def test_heatmap_bed_validates(self, tmp_path):
        f = tmp_path / "regions.bed"
        f.write_text("chrI\t100\t500\tname\n")
        cfg = {"full_analysis": True, "heatmap_target_file": str(f)}
        check_extra_output_files(cfg)

    def test_heatmap_bed_too_few_cols(self, tmp_path):
        f = tmp_path / "regions.bed"
        f.write_text("chrI\t100\n")
        cfg = {"full_analysis": True, "heatmap_target_file": str(f)}
        with pytest.raises(ValueError):
            check_extra_output_files(cfg)

    def test_motif_bed_validated_only_if_motifs_on(self, tmp_path):
        f = tmp_path / "motifs.bed"
        f.write_text("chrI\t100\n")  # malformed
        cfg_off = {"motifs": False, "motif_target_file": str(f)}
        check_extra_output_files(cfg_off)  # gate is off → no validation
        cfg_on = {"motifs": True, "motif_target_file": str(f)}
        with pytest.raises(ValueError):
            check_extra_output_files(cfg_on)

    def test_rnaseq_target_geneids_validates(self, tmp_path):
        f = tmp_path / "genes.txt"
        f.write_text("AT1G01010\nAT1G01020\tlabel2\n")
        cfg = {"rnaseq_target_file": str(f)}
        check_extra_output_files(cfg)

    def test_rnaseq_target_empty_first_col_rejected(self, tmp_path):
        f = tmp_path / "genes.txt"
        f.write_text("\tlabel\n")
        cfg = {"rnaseq_target_file": str(f)}
        with pytest.raises(ValueError):
            check_extra_output_files(cfg)

    def test_go_background_default_sentinel_skipped(self, tmp_path):
        cfg = {"GO": True, "rnaseq_background_file": "default"}
        check_extra_output_files(cfg)

    def test_go_background_path_existence_required_when_set(self, tmp_path):
        cfg = {"GO": True,
               "rnaseq_background_file": str(tmp_path / "missing.txt")}
        with pytest.raises(ValueError) as exc:
            check_extra_output_files(cfg)
        assert "rnaseq_background_file" in str(exc.value)

    def test_check_paths_false_bypass(self, tmp_path):
        cfg = {"full_analysis": True,
               "browser_target_file": str(tmp_path / "absent.bed")}
        check_extra_output_files(cfg, check_paths=False)


# ---------------------------------------------------------------------------
# Control chain depth (dual-role samples)
# ---------------------------------------------------------------------------

class TestControlChainDepth:
    """A control may declare its own Control (depth 2) but no deeper.

    Depth 2 is the *dual-role* case: a sample that is another row's control while
    also being analysed itself (e.g. an H3 ChIP used as H3K9me2's control and
    normalized against Input). Peak calling only resolves one step, so deeper
    chains cannot be expressed — and bounding the depth also makes cycles
    impossible.
    """

    _counter = [0]

    @classmethod
    def _ip(cls, sample_id, ip_target, control=""):
        # Unique accession per row: Read_files must not repeat across samples.
        cls._counter[0] += 1
        return _row(sample_id, "ChIP", f"SRR100000{cls._counter[0]}",
                    IP_target=ip_target, Control=control, Peak_type="broad")

    def test_depth_one_classic_pair_passes(self):
        df = pd.DataFrame([
            self._ip("Input_rep1", "Input"),
            self._ip("K9_rep1", "H3K9me2", "Input_rep1"),
        ])
        check_table(df, check_paths=False)  # no raise

    def test_depth_two_dual_role_passes(self):
        # H3K9me2 -> H3 -> Input: H3 is both a control and an analysis target.
        df = pd.DataFrame([
            self._ip("Input_rep1", "Input"),
            self._ip("H3_rep1", "H3", "Input_rep1"),
            self._ip("K9_rep1", "H3K9me2", "H3_rep1"),
        ])
        check_table(df, check_paths=False)  # no raise

    def test_depth_three_rejected(self):
        df = pd.DataFrame([
            self._ip("A", "Input"),
            self._ip("B", "H3", "A"),
            self._ip("C", "H3K9me2", "B"),
            self._ip("D", "H3K27me3", "C"),
        ])
        with pytest.raises(ValueError) as excinfo:
            check_table(df, check_paths=False)
        assert "too deep" in str(excinfo.value)

    def test_two_cycle_rejected(self):
        # Bounding depth at 2 makes loops impossible: A is its own grandparent
        # and still declares a Control, so the depth check fires.
        df = pd.DataFrame([
            self._ip("A", "H3", "B"),
            self._ip("B", "H3K9me2", "A"),
        ])
        with pytest.raises(ValueError) as excinfo:
            check_table(df, check_paths=False)
        assert "too deep" in str(excinfo.value)

    def test_self_control_rejected(self):
        df = pd.DataFrame([self._ip("A", "H3", "A")])
        with pytest.raises(ValueError) as excinfo:
            check_table(df, check_paths=False)
        assert "refers to the sample itself" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Peak_type: required only where peaks are actually called (issue #66)
# ---------------------------------------------------------------------------

class TestPeakTypeRequiredOnlyForPeakCallers:
    """A control-only row calls no peaks, so it needs no Peak_type.

    Requiring one forced users to invent a meaningless value on every
    Input/WCE/IgG row. The value never reached an output path either: control
    artefacts are named from the bare Sample_ID.
    """

    def _sheet(self, tmp_path, ip_peak_type, ctrl_peak_type):
        # Distinct files per row: a shared Read_files entry trips the
        # cross-row duplicate check, which is a different rule.
        ip_fq = tmp_path / "ip.fq.gz"
        ip_fq.write_bytes(b"")
        ctrl_fq = tmp_path / "ctrl.fq.gz"
        ctrl_fq.write_bytes(b"")
        return pd.DataFrame([
            _row("ip1", "ChIP", str(ip_fq), IP_target="H3K9me2",
                 Control="ctrl1", Peak_type=ip_peak_type),
            _row("ctrl1", "ChIP", str(ctrl_fq), IP_target="Input",
                 Control="", Peak_type=ctrl_peak_type),
        ])

    def test_control_only_row_may_omit_peak_type(self, tmp_path):
        check_table(self._sheet(tmp_path, "broad", ""))  # no raise

    def test_control_only_row_may_still_carry_one(self, tmp_path):
        # Existing sheets put an arbitrary value here; they must keep working.
        check_table(self._sheet(tmp_path, "broad", "narrow"))  # no raise

    def test_peak_caller_still_requires_peak_type(self, tmp_path):
        with pytest.raises(ValueError) as e:
            check_table(self._sheet(tmp_path, "", "broad"))
        assert "Peak_type is required" in str(e.value)
        assert "ip1" in str(e.value)

    def test_invalid_value_still_rejected_on_peak_caller(self, tmp_path):
        with pytest.raises(ValueError) as e:
            check_table(self._sheet(tmp_path, "wide", "broad"))
        assert "not in ['broad', 'narrow']" in str(e.value)

    def test_invalid_value_still_rejected_on_control_only_row(self, tmp_path):
        # Blank is allowed there, but a typo must not slip through.
        with pytest.raises(ValueError) as e:
            check_table(self._sheet(tmp_path, "broad", "wide"))
        assert "not in ['broad', 'narrow']" in str(e.value)
        assert "ctrl1" in str(e.value)

    def test_ip_without_control_may_omit_peak_type(self, tmp_path):
        # Not peak-callable either (is_peak_call_target drops it), so the same
        # relaxation applies.
        f = tmp_path / "reads.fq.gz"
        f.write_bytes(b"")
        df = pd.DataFrame([
            _row("orphan", "ChIP", str(f), IP_target="H3K9me2",
                 Control="", Peak_type=""),
        ])
        check_table(df)  # no raise

    def test_non_pulldown_assay_still_must_be_blank(self, tmp_path):
        f = tmp_path / "reads.fq.gz"
        f.write_bytes(b"")
        df = pd.DataFrame([_row("r1", "RNAseq", str(f), Peak_type="broad")])
        with pytest.raises(ValueError) as e:
            check_table(df)
        assert "must be blank" in str(e.value)
