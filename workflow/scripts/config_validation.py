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


def plot_theme_spec(config):
    """The plot-theme choice handed to the R plotting scripts, as JSON.

    `plot_theme` names a built-in theme of workflow/scripts/plot_theme.R and
    `plot_theme_overrides` replaces individual keys of it. Only the shape is
    checked here; plot_theme.R owns the key list and refuses unknown keys and
    wrong value types when a plot is drawn.

    Raises ValueError on a malformed setting so the Snakefile can refuse it at
    startup.
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
    return json.dumps({"theme": name, "overrides": overrides}, sort_keys=True)
