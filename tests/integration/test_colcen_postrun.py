"""
Post-run validation tests for the A. thaliana ColCEN integration run.

The ColCEN test config sets ``GO: true``, so a completed run builds and
installs the per-genome GO annotation package
``genomes_test_colcen/ColCEN/GO/org.Athaliana.eg.db`` via
``create_GO_database`` (makeOrgPackage -> install.packages) and then runs
``perform_GO_on_target_file``. This is the only integration coverage of that
code path — see issue #33, where the lazy-load step of the install
(``writeLines`` -> ``file(con, "w")``) failed on a user's cluster.

These tests require an existing completed ColCEN run and are skipped
automatically when none is found.

Run with: pytest tests/integration/test_colcen_postrun.py -v -m slow
"""

import shutil
import subprocess
import pytest
from pathlib import Path

from tests.integration.conftest import REPO_ROOT, load_output_dir, load_genome_dir

# ---------------------------------------------------------------------------
# Constants: results / genomes dirs and the expected GO package
# ---------------------------------------------------------------------------
RESULTS = REPO_ROOT / load_output_dir("test_options_colcen.yaml")
GENOMES = REPO_ROOT / load_genome_dir("test_options_colcen.yaml")

ANALYSIS_NAME = "test_colcen"
REF_GENOME = "ColCEN"
DBNAME = "org.Athaliana.eg.db"          # genus[0] + species => A + thaliana
GO_DB_DIR = GENOMES / REF_GENOME / "GO" / DBNAME


def _file_exists_nonempty(path):
    """Check that a file exists and is non-empty."""
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


# ---------------------------------------------------------------------------
# Guard fixture: skip all tests if no completed ColCEN run exists
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def results_exist():
    """Skip if the ColCEN run hasn't produced RNA output (GO runs off RNA)."""
    rna_dir = RESULTS / "RNA" / "mapped"
    if not rna_dir.exists() or not any(rna_dir.glob("final__*.bam")):
        pytest.skip("No completed ColCEN run found (results_test_colcen/RNA)")


# ===========================================================================
# TestGODatabase  — regression coverage for issue #33
# ===========================================================================

@pytest.mark.slow
@pytest.mark.integration
class TestGODatabase:
    """Validate the GO annotation package built by create_GO_database.

    A failed install (the #33 symptom) ends with R removing the package
    directory, so a completed run that produced these artifacts proves the
    makeOrgPackage -> install.packages -> lazy-load path succeeded end to end.
    """

    def test_go_package_installed(self, results_exist):
        """The installed org.db package directory and its key files exist."""
        assert GO_DB_DIR.is_dir(), f"GO package not installed: {GO_DB_DIR}"
        # Files written by the lazy-load step that fails in #33.
        for rel in (
            "DESCRIPTION",
            "Meta/package.rds",
            f"R/{DBNAME}.rdb",
            f"R/{DBNAME}.rdx",
        ):
            assert _file_exists_nonempty(GO_DB_DIR / rel), \
                f"Missing or empty installed-package file: {GO_DB_DIR / rel}"

    def test_go_source_tables_exist(self, results_exist):
        """The gaf/gene_info tables emitted alongside the package exist."""
        go_dir = GO_DB_DIR.parent
        for suffix in ("_gaf_file.tab", "_gene_info.tab"):
            tab = go_dir / f"{DBNAME}{suffix}"
            assert _file_exists_nonempty(tab), f"Missing or empty: {tab}"

    def test_go_package_lazyload(self, results_exist):
        """The package's lazy-load DB deserializes with base R.

        This is the direct regression guard for #33: the install step that
        failed there writes the R/<pkg>.rdb / .rdx lazy-load database, so
        loading it back exercises the exact artifact. base::lazyLoad needs
        only base R (no AnnotationDbi), so it runs in the epicc env. Skips
        if no Rscript is on PATH.
        """
        if shutil.which("Rscript") is None:
            pytest.skip("Rscript not on PATH; cannot exercise lazy-load DB")
        filebase = GO_DB_DIR / "R" / DBNAME
        rscript = (
            "e <- new.env(); "
            f'lazyLoad("{filebase}", envir = e); '
            "stopifnot(length(ls(e)) > 0)"
        )
        result = subprocess.run(
            ["Rscript", "-e", rscript],
            capture_output=True, text=True, timeout=120,
        )
        assert result.returncode == 0, (
            f"Loading GO package lazy-load DB failed:\n{result.stderr}"
        )

    def test_go_analysis_done(self, results_exist):
        """perform_GO_on_target_file produced its checkpoint for the DEGs."""
        done = RESULTS / "RNA" / "GO" / \
            f"TopGO__{ANALYSIS_NAME}__{REF_GENOME}__unique_DEGs.done"
        assert done.exists(), f"Missing GO analysis checkpoint: {done}"

    def test_go_outputs_carry_analysis_name(self, results_exist):
        """GO tables and treemaps are named per analysis and genome (issue #101).

        R_GO_analysis.R writes a table only when it has enriched terms, so an
        empty result is allowed; what must not happen is a run producing only
        the pre-#101 ``topGO_<label>_<ont>_*`` names, which a second analysis
        in the same output_dir would overwrite.
        """
        prefix = f"topGO__{ANALYSIS_NAME}__{REF_GENOME}__"
        go_dir = RESULTS / "RNA" / "GO"
        plots_dir = RESULTS / "RNA" / "plots"
        new = sorted(go_dir.glob(f"{prefix}*.txt")) + \
            sorted(plots_dir.glob(f"{prefix}*_treemap.pdf"))
        legacy = [
            p for d in (go_dir, plots_dir) for p in d.glob("topGO_*")
            if not p.name.startswith("topGO__")
        ]
        if not new and not legacy:
            pytest.skip("GO analysis wrote no enrichment tables (no enriched terms)")
        assert new, (
            f"Only legacy GO file names without the analysis name found: "
            f"{[p.name for p in legacy]}"
        )
