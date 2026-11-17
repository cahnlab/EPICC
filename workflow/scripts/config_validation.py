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
# heatmaps_plot_params (per env) and profiles_plot_params hold deeptools flags
# that set how plotHeatmap / plotProfile draw. Once a plot family is drawn in
# R, its flags are translated onto plot-theme keys with a deprecation warning,
# and any flag with no theme equivalent stops the run. Until then deeptools
# still reads them as before and nothing is translated. The Snakefile lists
# the families already drawn in R in R_DRAWN_PLOT_FAMILIES.

# Translations equal to the theme's defaults are dropped, so the shipped
# options file adds nothing to the theme JSON (and reruns nothing).
_DEFAULT_COLORMAP = {"mC": "Oranges"}          # every other env: seismic
_DEFAULT_MISSING = {"0.9", "grey90", "gray90", "#e6e6e6", "#e5e5e5"}
PROFILE_TYPES = ("lines", "fill", "se", "std")


def _colormap_names():
    import json
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_colormaps.json")
    with open(path) as fh:
        return set(json.load(fh)["maps"])


def _known_colormap(name, names):
    return name in names or (name.endswith("_r") and name[:-2] in names)


def _flag_values(words):
    """Split ['--a', 'x', 'y', '--b', 'z'] into [('--a', ['x', 'y']), ('--b', ['z'])]."""
    out = []
    for w in words:
        if w.startswith("--") or (w.startswith("-") and len(w) == 2 and w[1].isalpha()):
            out.append((w, []))
        elif out:
            out[-1][1].append(w)
        else:
            out.append((None, [w]))
    return out


def _missing_color(value):
    """deeptools takes a matplotlib colour; a bare number in [0, 1] is a grey."""
    try:
        f = float(value)
    except ValueError:
        return value
    if not 0 <= f <= 1:
        raise ValueError(f"--missingDataColor {value}: a grey level must be in [0, 1]")
    level = round(f * 255)
    return "#{0:02x}{0:02x}{0:02x}".format(level)


def legacy_plot_params(config, r_drawn=()):
    """Translate heatmaps_plot_params / profiles_plot_params onto theme keys.

    `r_drawn` names the plot families already drawn in R ("heatmaps",
    "profiles"); only their flags are translated. Returns (overrides, warnings):
    theme overrides to apply under the user's own plot_theme_overrides, and one
    deprecation warning per setting translated. Raises ValueError for a flag
    or value with no theme equivalent.

      --colorMap NAME         -> heatmap_colors.<env>: NAME (a matplotlib map)
      --colorList C1,C2,...   -> heatmap_colors.<env>: [C1, C2, ...]
      --missingDataColor X    -> heatmap_missing_color
      --interpolationMethod   -> dropped: R bins rows to the output resolution
      --plotType T            -> profile_type
    """
    overrides, warnings = {}, []
    if "heatmaps" in r_drawn and config.get("heatmaps_plot_params"):
        hp = config["heatmaps_plot_params"]
        if not isinstance(hp, dict):
            raise ValueError("heatmaps_plot_params must map each env to a flag string")
        names = _colormap_names()
        colors, missing = {}, {}
        for env, flags in hp.items():
            where = f"heatmaps_plot_params.{env}"
            for flag, vals in _flag_values(str(flags or "").split()):
                if flag == "--colorMap":
                    if len(vals) != 1:
                        raise ValueError(
                            f"{where}: --colorMap with {len(vals)} maps; the R heatmaps take one "
                            "colour map per data type. Set plot_theme_overrides.heatmap_colors.")
                    if not _known_colormap(vals[0], names):
                        raise ValueError(f"{where}: unknown colour map '{vals[0]}'")
                    if vals[0] != _DEFAULT_COLORMAP.get(env, "seismic"):
                        colors[env] = vals[0]
                elif flag == "--colorList":
                    if len(vals) != 1:
                        raise ValueError(
                            f"{where}: --colorList with {len(vals)} lists; the R heatmaps take one "
                            "list per data type. Set plot_theme_overrides.heatmap_colors.")
                    colors[env] = [c for c in vals[0].split(",") if c]
                    if len(colors[env]) < 2:
                        raise ValueError(f"{where}: --colorList needs at least two colours")
                elif flag == "--missingDataColor":
                    if len(vals) != 1:
                        raise ValueError(f"{where}: --missingDataColor takes one colour")
                    missing[env] = vals[0]
                elif flag == "--interpolationMethod":
                    pass
                else:
                    raise ValueError(
                        f"{where}: '{flag or vals[0]}' has no equivalent now that heatmaps are "
                        "drawn in R. Remove it; set the look with plot_theme / "
                        "plot_theme_overrides (keys: heatmap_colors, heatmap_missing_color, "
                        "heatmap_column_width, heatmap_height, profile_strip_height).")
        if colors:
            overrides["heatmap_colors"] = colors
        distinct = {v.lower() for v in missing.values()} - _DEFAULT_MISSING
        if len(distinct) > 1:
            raise ValueError(
                "heatmaps_plot_params: --missingDataColor differs between envs "
                f"({', '.join(f'{e}={v}' for e, v in sorted(missing.items()))}); the theme has one "
                "heatmap_missing_color. Set plot_theme_overrides.heatmap_missing_color.")
        if distinct:
            overrides["heatmap_missing_color"] = _missing_color(
                next(v for v in missing.values() if v.lower() in distinct))
        warnings.append(
            "heatmaps_plot_params is deprecated: heatmaps are drawn in R and take their look "
            "from plot_theme / plot_theme_overrides. "
            + (f"Translated to plot_theme_overrides: {overrides_summary(overrides)}. "
               if overrides else "Its settings match the theme defaults. ")
            + "Remove heatmaps_plot_params from the options file.")
    if "profiles" in r_drawn and config.get("profiles_plot_params"):
        before = dict(overrides)
        for flag, vals in _flag_values(str(config["profiles_plot_params"]).split()):
            if flag == "--plotType":
                if len(vals) != 1 or vals[0] not in PROFILE_TYPES:
                    raise ValueError(
                        f"profiles_plot_params: --plotType {' '.join(vals)} is not drawn in R; "
                        f"use one of {', '.join(PROFILE_TYPES)} (plot_theme_overrides.profile_type).")
                if vals[0] != "lines":
                    overrides["profile_type"] = vals[0]
            else:
                raise ValueError(
                    f"profiles_plot_params: '{flag or vals[0]}' has no equivalent now that "
                    "profiles are drawn in R. Remove it; set the look with plot_theme / "
                    "plot_theme_overrides (keys: profile_type, line_width, panel_width, "
                    "panel_height).")
        added = {k: v for k, v in overrides.items() if k not in before}
        warnings.append(
            "profiles_plot_params is deprecated: profiles are drawn in R and take their look "
            "from plot_theme / plot_theme_overrides. "
            + (f"Translated to plot_theme_overrides: {overrides_summary(added)}. "
               if added else "Its settings match the theme defaults. ")
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
    key by key (heatmap_colors entry by entry).

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
    merged = dict(legacy)
    for key, value in overrides.items():
        if key == "heatmap_colors" and isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = {**merged[key], **value}
        else:
            merged[key] = value
    return json.dumps({"theme": name, "overrides": merged}, sort_keys=True), warnings
