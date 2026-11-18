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


# --- deprecated deeptools plot flags -----------------------------------------
#
# profiles_plot_params holds deeptools flags that set how plotProfile drew.
# Profiles are drawn in R now, so its flags are translated onto plot-theme
# keys with a deprecation warning, and any flag with no theme equivalent stops
# the run. The Snakefile lists the families drawn in R in
# R_DRAWN_PLOT_FAMILIES. heatmaps_plot_params still goes to plotHeatmap
# untouched; its translation (--colorMap, --colorList, --missingDataColor,
# --interpolationMethod) returns with the R heatmaps in phase 3, from branch
# claude/r-plot-legacy-flags (commit fc1130c).

PROFILE_TYPES = ("lines", "fill", "se", "std")


def _flag_values(words):
    """Split ['--a', 'x', 'y', '--b=z'] into [('--a', ['x', 'y']), ('--b', ['z'])].

    A long flag may carry its value after '=' (argparse's --flag=value form,
    which deeptools accepted)."""
    out = []
    for w in words:
        if w.startswith("--"):
            flag, eq, value = w.partition("=")
            out.append((flag, [value] if eq else []))
        elif w.startswith("-") and len(w) == 2 and w[1].isalpha():
            out.append((w, []))
        elif out:
            out[-1][1].append(w)
        else:
            out.append((None, [w]))
    return out


def legacy_plot_params(config, r_drawn=()):
    """Translate profiles_plot_params onto theme keys.

    `r_drawn` names the plot families already drawn in R; only "profiles" is
    translated so far. Returns (overrides, warnings): theme overrides to apply
    under the user's own plot_theme_overrides, and one deprecation warning per
    key translated. Raises ValueError for a flag or value with no theme
    equivalent.

      --plotType T   -> profile_type (lines | fill | se | std)
    """
    overrides, warnings = {}, []
    if "profiles" in r_drawn and config.get("profiles_plot_params"):
        for flag, vals in _flag_values(str(config["profiles_plot_params"]).split()):
            if flag == "--plotType":
                if len(vals) != 1 or vals[0] not in PROFILE_TYPES:
                    raise ValueError(
                        f"profiles_plot_params: --plotType {' '.join(vals)} is not drawn in R; "
                        f"use one of {', '.join(PROFILE_TYPES)} (plot_theme_overrides.profile_type).")
                # The default is dropped, so the old shipped value adds
                # nothing to the theme JSON (and reruns nothing).
                if vals[0] != "lines":
                    overrides["profile_type"] = vals[0]
            else:
                raise ValueError(
                    f"profiles_plot_params: '{flag or vals[0]}' has no equivalent now that "
                    "profiles are drawn in R. Remove it; set the look with plot_theme / "
                    "plot_theme_overrides (keys: profile_type, line_width, panel_width, "
                    "panel_height).")
        warnings.append(
            "profiles_plot_params is deprecated: profiles are drawn in R and take their look "
            "from plot_theme / plot_theme_overrides. "
            + (f"Translated to plot_theme_overrides: {overrides_summary(overrides)}. "
               if overrides else "Its settings match the theme defaults. ")
            + "Remove profiles_plot_params from the options file.")
    return overrides, warnings


def overrides_summary(overrides):
    import json
    return ", ".join(f"{k}={json.dumps(v)}" for k, v in sorted(overrides.items()))


def plot_theme_spec(config, r_drawn=()):
    """The plot-theme choice handed to the R plotting scripts, as JSON.

    `plot_theme` names a built-in theme of workflow/scripts/plot_theme.R and
    `plot_theme_overrides` replaces individual keys of it. Only the shape is
    checked here; plot_theme.R owns the key list and refuses unknown keys and
    wrong value types when a plot is drawn.

    The deprecated deeptools flags of the families in `r_drawn` are translated
    (legacy_plot_params) and applied beneath plot_theme_overrides, which wins
    key by key.

    Returns (json, warnings). Raises ValueError on a malformed setting so the
    Snakefile can refuse it at startup.
    """
    import json

    name = config.get("plot_theme") or "default"
    overrides = config.get("plot_theme_overrides") or {}
    if not isinstance(name, str):
        raise ValueError(f"plot_theme must be a theme name, got {name!r}")
    if not isinstance(overrides, dict):
        raise ValueError(
            "plot_theme_overrides must be a mapping of theme keys to values, "
            f"got {type(overrides).__name__}")
    legacy, warnings = legacy_plot_params(config, r_drawn)
    merged = {**legacy, **overrides}
    return json.dumps({"theme": name, "overrides": merged}, sort_keys=True), warnings
