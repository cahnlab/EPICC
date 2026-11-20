"""Startup checks on the options file.

Kept separate from samplefile_validation.py, which validates sample sheets.
"""


def _walk(node, path, in_flags, problems):
    if isinstance(node, dict):
        for key, val in node.items():
            here = f"{path}.{key}" if path else str(key)
            flags = in_flags or key == "params" or str(key).endswith(
                ("_params", "_flags"))
            _walk(val, here, flags, problems)
    elif in_flags and isinstance(node, str) and ('"' in node or "'" in node):
        problems.append((path, node))


def find_quoted_params(config):
    """Locate quote characters in the free-form tool-parameter strings.

    Returns a list of (config path, offending value).

    These strings are pasted into rule bodies and expanded unquoted, so bash
    word-splits them but never strips quotes: "--keep-dup 'all'" reaches the
    tool as the literal argument 'all', which macs2 rejects. Quoting is neither
    needed nor supported, so callers should refuse it at startup rather than
    let a rule fail an hour into a run.
    """
    problems = []
    _walk(config, "", False, problems)
    return problems


def format_quoted_params_error(problems):
    """Render find_quoted_params output as an actionable error message."""
    lines = ["ERROR: quote characters in tool-parameter settings:"]
    for path, value in problems:
        stripped = value.replace('"', "").replace("'", "")
        lines.append(f"  {path}:")
        lines.append(f"      found:   {value}")
        lines.append(f"      use:     {stripped}")
    lines.append(
        "These strings are split into words but not unquoted, so the quotes "
        "would reach the tool as part of the argument."
    )
    return "\n".join(lines)


# --- plot theme --------------------------------------------------------------

PROFILE_TYPES = ("lines", "fill", "se", "std")
HEATMAP_INTERPOLATIONS = ("bilinear", "nearest")

# Data types the theme's per-type maps name (plot_theme.R's heatmap_colors /
# heatmap_interpolation, kept in step by a test), plus the mC contexts the
# display table uses as types. A map entry may also name a single mark.
THEME_TYPES = ("default", "ChIP", "TF", "ATAC", "RNA", "RNA_RAMPAGE", "sRNA", "mC",
               "mCG", "mCHG", "mCHH")

# Options of the deeptools-drawn plots, removed once heatmaps and profiles
# were drawn in R, and what replaces each.
REMOVED_PLOT_OPTIONS = {
    "heatmaps_plot_params": "plot_theme_overrides: heatmap_colors, heatmap_interpolation "
                            "and heatmap_missing_color (see config/epicc-options.yaml)",
    "profiles_plot_params": "plot_theme_overrides.profile_type (lines | fill | se | std)",
    "heatmap_sort_mc_after_others": "nothing: Heatmap__...__all__ sorts its rows on the "
                                    "non-mC columns; for mC sorted on itself, "
                                    "epicc output --plot-type heatmap --env mC",
}

# Theme keys that are maps from data type (or mark) to a value.
_TYPE_MAPS = ("heatmap_colors", "heatmap_interpolation")


def _colormap_names():
    import json
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_colormaps.json")
    with open(path) as fh:
        return set(json.load(fh)["maps"])


def _known_colormap(name, names):
    return name in names or (name.endswith("_r") and name[:-2] in names)


def _check_type_maps(overrides):
    """Refuse a malformed heatmap_colors / heatmap_interpolation entry at
    startup rather than when the first heatmap is drawn. Colour names in a
    list are left to R, which knows them."""
    names = None
    for key in _TYPE_MAPS:
        value = overrides.get(key)
        if value is None:
            continue
        if not isinstance(value, dict):
            raise ValueError(f"plot_theme_overrides.{key} must map data types or marks to values")
        for entry, v in value.items():
            at = f"plot_theme_overrides.{key}.{entry}"
            if key == "heatmap_interpolation":
                if v not in HEATMAP_INTERPOLATIONS:
                    raise ValueError(f"{at} must be one of {', '.join(HEATMAP_INTERPOLATIONS)}, got {v!r}")
            elif isinstance(v, str):
                names = names or _colormap_names()
                if not _known_colormap(v, names):
                    raise ValueError(
                        f"{at}: unknown colour map {v!r} (add _r to reverse one). "
                        f"Known maps: {', '.join(sorted(names))}")
            elif not (isinstance(v, list) and len(v) >= 2 and all(isinstance(c, str) for c in v)):
                raise ValueError(f"{at} must be a colour-map name or a list of at least two colours")


def plot_theme_spec(config):
    """The plot-theme choice handed to the R plotting scripts, as JSON.

    `plot_theme` names a built-in theme of workflow/scripts/plot_theme.R and
    `plot_theme_overrides` replaces individual keys of it. The shape is checked
    here, and the per-type maps' entries; plot_theme.R owns the key list and
    refuses unknown keys and wrong value types when a plot is drawn.

    Raises ValueError on a malformed setting, or on an option of the old
    deeptools plots (REMOVED_PLOT_OPTIONS), so the Snakefile can refuse it at
    startup.
    """
    import json

    removed = [k for k in REMOVED_PLOT_OPTIONS if k in config]
    if removed:
        raise ValueError(
            "Option(s) of the old deeptools plots, removed now that heatmaps and profiles "
            "are drawn in R: " + "; ".join(
                f"{k} -> {REMOVED_PLOT_OPTIONS[k]}" for k in removed)
            + ". Remove them from the options file.")
    name = config.get("plot_theme") or "default"
    overrides = config.get("plot_theme_overrides") or {}
    if not isinstance(name, str):
        raise ValueError(f"plot_theme must be a theme name, got {name!r}")
    if not isinstance(overrides, dict):
        raise ValueError(
            "plot_theme_overrides must be a mapping of theme keys to values, "
            f"got {type(overrides).__name__}")
    _check_type_maps(overrides)
    return json.dumps({"theme": name, "overrides": overrides}, sort_keys=True)


def unmatched_theme_entries(config, marks):
    """heatmap_colors / heatmap_interpolation entries that name neither a data
    type of the theme nor one of `marks`, the marks this run plots. Such an
    entry is never looked up -- a typo, most likely -- and its columns fall
    back to their data type's entry."""
    overrides = config.get("plot_theme_overrides") or {}
    known = set(THEME_TYPES) | set(marks)
    return [f"plot_theme_overrides.{key}.{entry}"
            for key in _TYPE_MAPS
            for entry in (overrides.get(key) or {})
            if entry not in known]
