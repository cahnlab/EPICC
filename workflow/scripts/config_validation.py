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
