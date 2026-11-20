"""Which methylation contexts each DMR pair is called in.

``methylation_contexts`` sets the contexts for the whole run: bigwigs, DMRs
and PCA. Direct-methylation (dmC) data can carry fewer. Most ONT basecalling
models call 5mC only at CpG (Dorado's ``5mCG_5hmCG`` models), so their
modkit bedMethyl has no CHG or CHH position at all. Paired with a bisulfite
sample in CHG, such a group has nothing to test, and metilene cannot be given
an empty group: it segfaults. The DMR script therefore refuses one, which is
right for a broken sample and wrong for data that never had the context.

Which contexts a dmC sample carries cannot be read before the DAG is built.
A bedMethyl has no context column (``bedmethyl_to_cx_report.py`` derives it
from the reference), and a modBAM's MM tags don't say which model wrote them.
So it is declared: ``dmc_methylation.contexts`` lists the contexts dmC data
holds, and a pair with a dmC side is called only in the contexts both sides
have. Left unset, dmC follows ``methylation_contexts`` and nothing changes.
"""

ALL_CONTEXTS = ("CG", "CHG", "CHH")


def normalize_contexts(value, key):
    """Validate a context list and return it in CG, CHG, CHH order.

    ``key`` names the option in the error message.
    """
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError(
            f"{key} must be a non-empty list, e.g. [\"CG\", \"CHG\", \"CHH\"] "
            "or [\"CG\"]"
        )
    for c in value:
        if c not in ALL_CONTEXTS:
            raise ValueError(
                f"{key} entry {c!r} not in {list(ALL_CONTEXTS)} "
                "(CAG/CAA and other subcontexts are not yet supported)"
            )
    return [c for c in ALL_CONTEXTS if c in value]


def dmc_contexts(config, run_contexts):
    """The contexts dmC data is analysed in.

    ``dmc_methylation.contexts`` intersected with the run's contexts, so a
    run-wide ``[CG]`` still wins. Unset (or null) means the run's contexts.
    """
    declared = (config.get("dmc_methylation") or {}).get("contexts")
    if declared is None:
        return list(run_contexts)
    declared = normalize_contexts(declared, "dmc_methylation.contexts")
    return [c for c in run_contexts if c in declared]


def dmr_contexts_for_pair(run_contexts, dmc_ctx, side1_is_dmc, side2_is_dmc):
    """Contexts a pair of analysis groups is compared in.

    A side whose replicates are dmC contributes only ``dmc_ctx``; a
    bisulfite side contributes every run context.
    """
    allowed = set(run_contexts)
    if side1_is_dmc or side2_is_dmc:
        allowed &= set(dmc_ctx)
    return [c for c in run_contexts if c in allowed]
