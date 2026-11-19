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
# that set how plotHeatmap / plotProfile drew. Both families are drawn in R
# now, so their flags are translated onto plot-theme keys with a deprecation
# warning, and any flag with no theme equivalent stops the run. The Snakefile
# lists the families drawn in R in R_DRAWN_PLOT_FAMILIES; only their flags are
# translated.

PROFILE_TYPES = ("lines", "fill", "se", "std")
HEATMAP_INTERPOLATIONS = ("bilinear", "nearest")

# Translations equal to the theme's defaults are dropped, so the shipped
# options file adds nothing to the theme JSON (and reruns nothing).
_DEFAULT_COLORMAP = {"mC": "Oranges"}           # every other type: seismic
_DEFAULT_INTERPOLATION = {"mC": "nearest"}      # every other type: bilinear
_DEFAULT_MISSING = {"0.9", "grey90", "gray90", "#e6e6e6", "#e5e5e5"}

# The other values matplotlib's imshow takes for --interpolationMethod, and
# the one of the two R modes each is closer to.
_MATPLOTLIB_INTERPOLATIONS = {
    "none": "nearest",
    **{name: "bilinear" for name in (
        "antialiased", "auto", "bicubic", "spline16", "spline36", "hanning", "hamming",
        "hermite", "kaiser", "quadric", "catrom", "gaussian", "bessel", "mitchell", "sinc",
        "lanczos", "blackman")},
}

# heatmaps_plot_params is keyed by env; the theme maps are keyed by data type.
# 'all' and 'most' hold every type, so they set the fallback entry only.
_ENV_TO_TYPE = {"all": "default", "most": "default"}


def _colormap_names():
    import json
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_colormaps.json")
    with open(path) as fh:
        return set(json.load(fh)["maps"])


def _known_colormap(name, names):
    return name in names or (name.endswith("_r") and name[:-2] in names)


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


def _set_type_entry(store, dtype, value, where, key):
    """One theme-map entry per data type; 'all' and 'most' both land on
    'default', so they must agree."""
    if dtype in store and store[dtype][0] != value:
        raise ValueError(
            f"{where}: {key} {value!r} differs from {store[dtype][1]}'s {store[dtype][0]!r}; "
            f"both set heatmap_{key}.{dtype}. Set plot_theme_overrides.heatmap_{key}.")
    store[dtype] = (value, where)


def _translate_heatmap_params(hp):
    """heatmaps_plot_params -> (overrides, notes)."""
    if not isinstance(hp, dict):
        raise ValueError("heatmaps_plot_params must map each env to a flag string")
    names = _colormap_names()
    colors, interp, missing, notes = {}, {}, {}, []
    for env, flags in hp.items():
        where = f"heatmaps_plot_params.{env}"
        dtype = _ENV_TO_TYPE.get(env, env)
        for flag, vals in _flag_values(str(flags or "").split()):
            if flag == "--colorMap":
                if len(vals) != 1:
                    raise ValueError(
                        f"{where}: --colorMap with {len(vals)} maps; the R heatmaps take one "
                        "colour map per data type. Set plot_theme_overrides.heatmap_colors.")
                if not _known_colormap(vals[0], names):
                    raise ValueError(f"{where}: unknown colour map '{vals[0]}'")
                _set_type_entry(colors, dtype, vals[0], where, "colors")
            elif flag == "--colorList":
                if len(vals) != 1:
                    raise ValueError(
                        f"{where}: --colorList with {len(vals)} lists; the R heatmaps take one "
                        "list per data type. Set plot_theme_overrides.heatmap_colors.")
                stops = [c for c in vals[0].split(",") if c]
                if len(stops) < 2:
                    raise ValueError(f"{where}: --colorList needs at least two colours")
                _set_type_entry(colors, dtype, stops, where, "colors")
            elif flag == "--missingDataColor":
                if len(vals) != 1:
                    raise ValueError(f"{where}: --missingDataColor takes one colour")
                missing[env] = vals[0]
            elif flag == "--interpolationMethod":
                if len(vals) != 1:
                    raise ValueError(f"{where}: --interpolationMethod takes one method")
                method = vals[0]
                if method not in HEATMAP_INTERPOLATIONS:
                    if method not in _MATPLOTLIB_INTERPOLATIONS:
                        raise ValueError(
                            f"{where}: --interpolationMethod {method} is not a matplotlib "
                            "interpolation method")
                    sub = _MATPLOTLIB_INTERPOLATIONS[method]
                    notes.append(f"{where}: --interpolationMethod {method} is drawn as "
                                 f"'{sub}' (the R heatmaps do bilinear or nearest)")
                    method = sub
                _set_type_entry(interp, dtype, method, where, "interpolation")
            else:
                raise ValueError(
                    f"{where}: '{flag or vals[0]}' has no equivalent now that heatmaps are "
                    "drawn in R. Remove it; set the look with plot_theme / "
                    "plot_theme_overrides (keys: heatmap_colors, heatmap_interpolation, "
                    "heatmap_missing_color, heatmap_column_width, heatmap_height, "
                    "profile_strip_height).")
    overrides = {}
    colors = {t: v for t, (v, _) in colors.items()
              if v != _DEFAULT_COLORMAP.get(t, "seismic")}
    if colors:
        overrides["heatmap_colors"] = colors
    interp = {t: v for t, (v, _) in interp.items()
              if v != _DEFAULT_INTERPOLATION.get(t, "bilinear")}
    if interp:
        overrides["heatmap_interpolation"] = interp
    distinct = {v.lower() for v in missing.values()} - _DEFAULT_MISSING
    if len(distinct) > 1:
        raise ValueError(
            "heatmaps_plot_params: --missingDataColor differs between envs "
            f"({', '.join(f'{e}={v}' for e, v in sorted(missing.items()))}); the theme has one "
            "heatmap_missing_color. Set plot_theme_overrides.heatmap_missing_color.")
    if distinct:
        overrides["heatmap_missing_color"] = _missing_color(
            next(v for v in missing.values() if v.lower() in distinct))
    return overrides, notes


def legacy_plot_params(config, r_drawn=()):
    """Translate heatmaps_plot_params / profiles_plot_params onto theme keys.

    `r_drawn` names the plot families drawn in R ("heatmaps", "profiles");
    only their flags are translated. Returns (overrides, warnings): theme
    overrides to apply under the user's own plot_theme_overrides, and one
    deprecation warning per option translated. Raises ValueError for a flag or
    value with no theme equivalent.

      --colorMap NAME           -> heatmap_colors.<type>: NAME (a matplotlib map)
      --colorList C1,C2,...     -> heatmap_colors.<type>: [C1, C2, ...]
      --missingDataColor X      -> heatmap_missing_color
      --interpolationMethod M   -> heatmap_interpolation.<type>: bilinear | nearest
                                   (other matplotlib methods: the closer of the
                                   two, with a note in the warning)
      --plotType T              -> profile_type (lines | fill | se | std)

    <type> is the env key of heatmaps_plot_params, except 'all' and 'most',
    which set the 'default' entry (the fallback for every type).
    """
    overrides, warnings = {}, []
    if "heatmaps" in r_drawn and config.get("heatmaps_plot_params"):
        added, notes = _translate_heatmap_params(config["heatmaps_plot_params"])
        overrides.update(added)
        warnings.append(
            "heatmaps_plot_params is deprecated: heatmaps are drawn in R and take their look "
            "from plot_theme / plot_theme_overrides. "
            + (f"Translated to plot_theme_overrides: {overrides_summary(added)}. "
               if added else "Its settings match the theme defaults. ")
            + "".join(f"{n}. " for n in notes)
            + "Remove heatmaps_plot_params from the options file.")
    if "heatmaps" in r_drawn and "heatmap_sort_mc_after_others" in config:
        warnings.append(
            "heatmap_sort_mc_after_others is deprecated and ignored: the default heatmap "
            "(Heatmap__<matrix>__all__...) draws every data type in one figure, its rows "
            "sorted on the non-mC columns. For mC sorted on itself, make the mC-only heatmap: "
            "epicc output --plot-type heatmap --env mC. Remove the option from the options file.")
    if "profiles" in r_drawn and config.get("profiles_plot_params"):
        added = {}
        for flag, vals in _flag_values(str(config["profiles_plot_params"]).split()):
            if flag == "--plotType":
                if len(vals) != 1 or vals[0] not in PROFILE_TYPES:
                    raise ValueError(
                        f"profiles_plot_params: --plotType {' '.join(vals)} is not drawn in R; "
                        f"use one of {', '.join(PROFILE_TYPES)} (plot_theme_overrides.profile_type).")
                # The default is dropped, so the old shipped value adds
                # nothing to the theme JSON (and reruns nothing).
                if vals[0] != "lines":
                    added["profile_type"] = vals[0]
            else:
                raise ValueError(
                    f"profiles_plot_params: '{flag or vals[0]}' has no equivalent now that "
                    "profiles are drawn in R. Remove it; set the look with plot_theme / "
                    "plot_theme_overrides (keys: profile_type, line_width, panel_width, "
                    "panel_height).")
        overrides.update(added)
        warnings.append(
            "profiles_plot_params is deprecated: profiles are drawn in R and take their look "
            "from plot_theme / plot_theme_overrides. "
            + (f"Translated to plot_theme_overrides: {overrides_summary(added)}. "
               if added else "Its settings match the theme defaults. ")
            + "Remove profiles_plot_params from the options file.")
    return overrides, warnings


# Theme keys that are maps from data type to a value (plot_theme.R's open maps).
_TYPE_MAPS = ("heatmap_colors", "heatmap_interpolation")


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
    key by key (the per-type maps heatmap_colors and heatmap_interpolation
    entry by entry).

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
        if key in _TYPE_MAPS and isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = {**merged[key], **value}
        else:
            merged[key] = value
    return json.dumps({"theme": name, "overrides": merged}, sort_keys=True), warnings
