"""Tests for workflow/scripts/mc_contexts.py: which contexts a DMR pair runs in."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))

from scripts.mc_contexts import (  # noqa: E402
    dmc_contexts,
    dmr_contexts_for_pair,
    normalize_contexts,
)

ALL = ["CG", "CHG", "CHH"]


class TestNormalizeContexts:
    def test_orders_cg_chg_chh(self):
        assert normalize_contexts(["CHH", "CG"], "k") == ["CG", "CHH"]

    def test_drops_duplicates(self):
        assert normalize_contexts(["CG", "CG"], "k") == ["CG"]

    @pytest.mark.parametrize("bad", [[], "CG", None, {"CG": 1}])
    def test_rejects_non_list_or_empty(self, bad):
        with pytest.raises(ValueError, match="non-empty list"):
            normalize_contexts(bad, "dmc_methylation.contexts")

    def test_rejects_subcontext_and_names_key(self):
        with pytest.raises(ValueError, match="dmc_methylation.contexts entry 'CAG'"):
            normalize_contexts(["CG", "CAG"], "dmc_methylation.contexts")


class TestDmcContexts:
    def test_unset_follows_run(self):
        assert dmc_contexts({}, ALL) == ALL
        assert dmc_contexts({"dmc_methylation": {"pileup": {}}}, ALL) == ALL

    def test_null_follows_run(self):
        assert dmc_contexts({"dmc_methylation": {"contexts": None}}, ALL) == ALL
        assert dmc_contexts({"dmc_methylation": None}, ALL) == ALL

    def test_cpg_only(self):
        assert dmc_contexts({"dmc_methylation": {"contexts": ["CG"]}}, ALL) == ["CG"]

    def test_limited_by_run_contexts(self):
        cfg = {"dmc_methylation": {"contexts": ["CG", "CHG", "CHH"]}}
        assert dmc_contexts(cfg, ["CG"]) == ["CG"]

    def test_disjoint_from_run_is_empty(self):
        cfg = {"dmc_methylation": {"contexts": ["CHH"]}}
        assert dmc_contexts(cfg, ["CG"]) == []

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            dmc_contexts({"dmc_methylation": {"contexts": ["5mC"]}}, ALL)


class TestDmrContextsForPair:
    def test_bisulfite_pair_gets_every_run_context(self):
        assert dmr_contexts_for_pair(ALL, ["CG"], False, False) == ALL

    @pytest.mark.parametrize("s1,s2", [(True, False), (False, True), (True, True)])
    def test_pair_with_dmc_side_limited_to_dmc_contexts(self, s1, s2):
        assert dmr_contexts_for_pair(ALL, ["CG"], s1, s2) == ["CG"]

    def test_all_context_dmc_unchanged(self):
        assert dmr_contexts_for_pair(ALL, ALL, True, False) == ALL

    def test_empty_when_nothing_shared(self):
        assert dmr_contexts_for_pair(["CG"], [], True, False) == []

    def test_keeps_run_order(self):
        assert dmr_contexts_for_pair(ALL, ["CHH", "CG"], True, False) == ["CG", "CHH"]
