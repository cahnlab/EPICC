"""Normalize the ``browser_extra_tracks`` option into plottable track specs.

Each entry is either a bare BED path or a mapping with ``file`` and optional
``label`` and ``color``:

    browser_extra_tracks:
      - data/peaks.bed
      - {file: data/dmrs.bed.gz, label: DMRs, color: "#B2182B"}

Missing labels come from the file name, missing colours from a fixed palette
in entry order, so the same options always give the same plot. Labels and
colours end up in a TSV that the browser rule reads in a shell loop, which is
why tabs, newlines and shell-active characters are rejected rather than
escaped.
"""

import os
import re

# Set2-derived; skips the light green the TE track already uses.
EXTRA_TRACK_PALETTE = ("#E78AC3", "#8DA0CB", "#FC8D62", "#A6D854", "#FFD92F", "#E5C494")

_UNSAFE_LABEL_CHARS = set('\t\n\r"$`\\')
_COLOR_RE = re.compile(r"^(#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?|[A-Za-z][A-Za-z0-9]*)$")


def label_from_path(path):
    """'data/peaks.bed.gz' -> 'peaks'."""
    name = os.path.basename(str(path))
    while True:
        stem, ext = os.path.splitext(name)
        if not ext:
            break
        name = stem
    return name or "track"


def normalize_extra_tracks(entries):
    """Return ``[{"file", "label", "color"}, ...]``, raising ValueError on bad input.

    ``entries`` may be None or empty (no extra tracks).
    """
    if not entries:
        return []
    if not isinstance(entries, list):
        raise ValueError(
            "browser_extra_tracks must be a list of BED paths or "
            "{file, label, color} entries"
        )
    tracks = []
    for i, entry in enumerate(entries):
        where = f"browser_extra_tracks[{i}]"
        if isinstance(entry, str):
            entry = {"file": entry}
        if not isinstance(entry, dict):
            raise ValueError(f"{where}: expected a path or a mapping, got {entry!r}")
        unknown = set(entry) - {"file", "label", "color"}
        if unknown:
            raise ValueError(f"{where}: unknown key(s) {', '.join(sorted(unknown))}")
        path = str(entry.get("file") or "").strip()
        if not path:
            raise ValueError(f"{where}: 'file' is required")
        label = str(entry.get("label") or label_from_path(path)).strip()
        if not label or _UNSAFE_LABEL_CHARS & set(label):
            raise ValueError(
                f"{where}: label {label!r} is empty or contains a tab, newline, "
                f"quote, '$', backtick or backslash"
            )
        color = str(entry.get("color") or EXTRA_TRACK_PALETTE[i % len(EXTRA_TRACK_PALETTE)]).strip()
        if not _COLOR_RE.match(color):
            raise ValueError(
                f"{where}: color {color!r} is not a hex colour (#RRGGBB) or an R colour name"
            )
        tracks.append({"file": path, "label": label, "color": color})
    return tracks
