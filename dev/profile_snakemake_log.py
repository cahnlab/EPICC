#!/usr/bin/env python3
"""Parse a Snakemake log file and produce an execution profile report.

Usage:
    python dev/profile_snakemake_log.py .snakemake/log/<logfile>.snakemake.log
    python dev/profile_snakemake_log.py --html report.html <logfile>
    python dev/profile_snakemake_log.py                   # aggregate latest run
    python dev/profile_snakemake_log.py --single          # newest single log only

By default, all logs from the same resumed run are aggregated. Run identity is
derived from log content (output_dir + analysis_name wildcards), not timestamps.
"""

import argparse
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timedelta
from pathlib import Path

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

TIMESTAMP_RE = re.compile(r"^\[(.+?)\]$")
RULE_START_RE = re.compile(r"^(?:localrule|rule) (\S+):$")
JOBID_RE = re.compile(r"^\s+jobid:\s*(\d+)")
WILDCARDS_RE = re.compile(r"^\s+wildcards:\s*(.+)")
THREADS_RE = re.compile(r"^\s+threads:\s*(\d+)")
FINISHED_RE = re.compile(r"^Finished jobid:\s*(\d+)\s*\(Rule:\s*(\S+)\)")
STEPS_RE = re.compile(r"^(\d+) of (\d+) steps")
# SLURM-mode submission line: encodes rule and wildcard path in the slurm log path.
SLURM_SUBMIT_RE = re.compile(
    r"^Job (\d+) has been submitted with SLURM jobid \d+ \(log: .*?/rule_(\S+?)/(\S+?)/[^/]+\.log\)"
)

# Rule-name prefix → phase mapping (used for analysis and Gantt coloring)
PHASE_MAP = {
    "get_fastq": "Download", "process_fastq": "Download",
    "get_available_bam": "Download", "get_modbam": "Download",
    "get_bedmethyl": "Download",
    "filter_bam": "Alignment", "filter_rna": "Alignment",
    "STAR_map": "Alignment", "map_chromap": "Alignment",
    "bowtie2_map": "Alignment", "bismark_map": "Alignment",
    "align_modbam": "Alignment", "shortstack_map": "Alignment",
    "filter_structural_rna": "Alignment", "dispatch_srna_fastq": "Alignment",
    "merge_replicates": "Merge", "merging_": "Merge",
    "making_pseudo_replicates": "Merge",
    "calling_peaks": "Peak calling", "best_peaks": "Peak calling",
    "idr_analysis": "Peak calling", "select_peaks": "Peak calling",
    "perform_pairwise_diff": "Differential",
    "call_all_DEGs": "Differential", "call_all_differential": "Differential",
    "call_DMR": "Differential",
    "make_bigwig": "Tracks", "make_rna_stranded": "Tracks",
    "make_rna_unstranded": "Tracks", "make_srna_stranded": "Tracks",
    "make_mc_bigwig": "Tracks", "make_coverage": "Tracks",
    "computing_matrix": "Visualization", "making_stranded_matrix": "Visualization",
    "merging_matrix": "Visualization", "computing_matrix_scales": "Visualization",
    "sort_heatmap": "Visualization", "plot": "Visualization",
    "prep_browser": "Visualization", "prep_chromosomes": "Visualization",
    "make_single_loci": "Visualization", "merge_region_browser": "Visualization",
    "summarize_tracks": "Visualization",
    "make_fingerprint": "QC", "make_mapping_stats": "QC",
    "make_rna_stats": "QC", "make_mc_stats": "QC",
    "make_peak_stats": "QC", "make_srna_size": "QC",
    "has_header": "QC", "is_stranded": "QC", "run_fastqc": "QC",
    "modkit_summary": "QC",
    "filter_size_srna": "sRNA", "analyze_all_srna": "sRNA",
    "make_cluster": "sRNA", "combine_cluster": "sRNA",
    "deduplicate_srna": "sRNA",
    "check_fasta": "Setup", "check_gff": "Setup", "check_gtf": "Setup",
    "check_chrom_sizes": "Setup", "compute_genome_stats": "Setup",
    "resolve_taxid": "Setup", "prep_region_file": "Setup",
    "download_rfam": "Setup", "build_structural_rna_db": "Setup",
    "download_rfam_family_table": "Setup", "select_structural_rna_classes": "Setup",
    "check_te_file": "Setup", "make_bt2_indices": "Setup",
    "make_chromap_index": "Setup", "make_bismark_indices": "Setup",
    "make_STAR_indices": "Setup", "make_bowtie1_indices": "Setup",
    "create_GO_database": "Setup",
}

# Colorblind-friendly palette (Okabe-Ito derived)
PHASE_PALETTE = {
    "Download":      "#0072B2",  # blue
    "Alignment":     "#D55E00",  # vermillion
    "Merge":         "#009E73",  # bluish green
    "Peak calling":  "#CC79A7",  # reddish purple
    "Differential":  "#E69F00",  # orange
    "Tracks":        "#56B4E9",  # sky blue
    "Visualization": "#F0E442",  # yellow
    "QC":            "#555555",  # dark gray
    "sRNA":          "#882255",  # wine (Paul Tol extension)
    "Setup":         "#999999",  # gray
    "Other":         "#AA4499",  # purple (Paul Tol extension)
}


def _rule_to_phase(rule_name):
    """Map a rule name to its phase via prefix matching."""
    for prefix, phase in PHASE_MAP.items():
        if rule_name.startswith(prefix):
            return phase
    return "Other"


def parse_log(path):
    """Parse Snakemake log, return list of completed job records."""
    jobs_pending = {}  # jobid -> {rule, start, wildcards, threads}
    jobs_done = []
    current_ts = None
    current_rule = None

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")

            m = TIMESTAMP_RE.match(line)
            if m:
                try:
                    current_ts = datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y")
                except ValueError:
                    pass
                continue

            m = RULE_START_RE.match(line)
            if m:
                current_rule = m.group(1)
                continue

            m = JOBID_RE.match(line)
            if m and current_rule and current_ts:
                jobid = int(m.group(1))
                jobs_pending[jobid] = {
                    "rule": current_rule,
                    "start": current_ts,
                    "wildcards": "",
                    "threads": 1,
                }
                current_rule = None
                continue

            m = WILDCARDS_RE.match(line)
            if m:
                # attach to most recently seen jobid
                for jid in reversed(list(jobs_pending)):
                    if "wildcards" in jobs_pending[jid]:
                        jobs_pending[jid]["wildcards"] = m.group(1)
                        break
                continue

            m = THREADS_RE.match(line)
            if m:
                for jid in reversed(list(jobs_pending)):
                    jobs_pending[jid]["threads"] = int(m.group(1))
                    break
                continue

            m = SLURM_SUBMIT_RE.match(line)
            if m and current_ts:
                jobid = int(m.group(1))
                # Don't clobber a pending entry — if the verbose rule block already
                # registered this jobid (with proper wildcards), keep that.
                if jobid not in jobs_pending:
                    jobs_pending[jobid] = {
                        "rule": m.group(2),
                        "start": current_ts,
                        "wildcards": m.group(3),
                        "threads": 1,
                    }
                continue

            m = FINISHED_RE.match(line)
            if m and current_ts:
                jobid = int(m.group(1))
                if jobid in jobs_pending:
                    rec = jobs_pending.pop(jobid)
                    rec["end"] = current_ts
                    rec["duration"] = (rec["end"] - rec["start"]).total_seconds()
                    rec["jobid"] = jobid
                    jobs_done.append(rec)
                continue

    return jobs_done


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def fmt_duration(seconds):
    """Format seconds as human-readable duration."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        m, s = divmod(seconds, 60)
        return f"{int(m)}m {int(s)}s"
    else:
        h, rem = divmod(seconds, 3600)
        m, s = divmod(rem, 60)
        return f"{int(h)}h {int(m)}m {int(s)}s"


def analyze(jobs):
    """Compute summary statistics from job records."""
    if not jobs:
        return {}

    total_start = min(j["start"] for j in jobs)
    total_end = max(j["end"] for j in jobs)
    wall_time = (total_end - total_start).total_seconds()

    # Per-rule aggregates
    by_rule = defaultdict(lambda: {"count": 0, "total_sec": 0.0, "max_sec": 0.0, "jobs": []})
    for j in jobs:
        r = by_rule[j["rule"]]
        r["count"] += 1
        r["total_sec"] += j["duration"]
        r["max_sec"] = max(r["max_sec"], j["duration"])
        r["jobs"].append(j)

    # Sort by total time descending
    rules_ranked = sorted(by_rule.items(), key=lambda x: -x[1]["total_sec"])

    # Top 10 longest individual jobs
    top_jobs = sorted(jobs, key=lambda j: -j["duration"])[:15]

    # Phase analysis: group rules into logical phases
    phases = defaultdict(lambda: {"total_sec": 0.0, "count": 0})
    for j in jobs:
        phase = _rule_to_phase(j["rule"])
        phases[phase]["total_sec"] += j["duration"]
        phases[phase]["count"] += j["count"] if "count" in j else 1

    return {
        "total_jobs": len(jobs),
        "wall_time": wall_time,
        "total_start": total_start,
        "total_end": total_end,
        "cpu_time": sum(j["duration"] * j["threads"] for j in jobs),
        "rules_ranked": rules_ranked,
        "top_jobs": top_jobs,
        "phases": dict(sorted(phases.items(), key=lambda x: -x[1]["total_sec"])),
    }


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def report_markdown(stats):
    lines = []
    lines.append("# Snakemake Execution Profile\n")
    if stats["total_start"].date() == stats["total_end"].date():
        period = (f"{stats['total_start'].strftime('%Y-%m-%d %H:%M:%S')} — "
                  f"{stats['total_end'].strftime('%H:%M:%S')}")
    else:
        period = (f"{stats['total_start'].strftime('%Y-%m-%d %H:%M:%S')} — "
                  f"{stats['total_end'].strftime('%Y-%m-%d %H:%M:%S')}")
    parallelism = (f"{stats['cpu_time'] / stats['wall_time']:.1f}x"
                   if stats["wall_time"] > 0 else "n/a")
    if stats.get("num_logs", 1) > 1:
        lines.append(f"_Aggregated across {stats['num_logs']} resumed-run logs; "
                     f"wall time excludes idle gaps between resumptions._\n")
    lines.append("| Run period | Wall time | Total jobs | Total CPU time | Avg parallelism |")
    lines.append("|------------|-----------|------------|----------------|-----------------|")
    lines.append(f"| {period} | {fmt_duration(stats['wall_time'])} | "
                 f"{stats['total_jobs']} | {fmt_duration(stats['cpu_time'])} | "
                 f"{parallelism} |")
    lines.append("")

    # Phase summary — % is share of total job time, so it sums to 100%
    # regardless of parallelism. (Job time = sum of durations across all jobs.)
    total_job_sec = sum(d["total_sec"] for d in stats["phases"].values()) or 1.0
    lines.append("## Phase Summary\n")
    lines.append("| Phase | Jobs | Total time | % of total |")
    lines.append("|-------|------|-----------|------------|")
    for phase, data in stats["phases"].items():
        pct = data["total_sec"] / total_job_sec * 100
        lines.append(f"| {phase} | {data['count']} | "
                     f"{fmt_duration(data['total_sec'])} | {pct:.1f}% |")

    # Top individual jobs
    lines.append("\n## Slowest Individual Jobs\n")
    lines.append("| # | Rule | Duration | Wildcards |")
    lines.append("|---|------|----------|-----------|")
    for i, j in enumerate(stats["top_jobs"], 1):
        wc = j["wildcards"][:60] if j["wildcards"] else ""
        lines.append(f"| {i} | {j['rule']} | {fmt_duration(j['duration'])} | {wc} |")

    # Per-rule table
    lines.append("\n## Per-Rule Breakdown\n")
    lines.append("| Rule | Count | Total | Mean | Max |")
    lines.append("|------|-------|-------|------|-----|")
    for rule, data in stats["rules_ranked"]:
        mean = data["total_sec"] / data["count"]
        lines.append(f"| {rule} | {data['count']} | "
                     f"{fmt_duration(data['total_sec'])} | "
                     f"{fmt_duration(mean)} | "
                     f"{fmt_duration(data['max_sec'])} |")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def _md_to_html(md_text):
    """Convert markdown to HTML, returning (summary, details) tuple.

    Summary = everything before the first ## heading (title + stats table).
    Details = everything from the first ## heading onward.
    """
    rows = md_text.split("\n")
    summary_parts = []
    detail_parts = []
    in_table = False
    seen_h2 = False

    for row in rows:
        target = detail_parts if seen_h2 else summary_parts

        if row.startswith("# "):
            if in_table:
                target.append("</table>")
                in_table = False
            target.append(f"<h1>{row[2:]}</h1>")
        elif row.startswith("## "):
            if in_table:
                target.append("</table>")
                in_table = False
            if not seen_h2:
                seen_h2 = True
                target = detail_parts
            target.append(f"<h2>{row[3:]}</h2>")
        elif row.startswith("**"):
            target.append(f"<p>{row}</p>")
        elif row.startswith("|") and not row.startswith("|--"):
            cells = [c.strip() for c in row.split("|")[1:-1]]
            if not in_table:
                target.append('<table>')
                tag = "th"
                in_table = True
            else:
                tag = "td"
            target.append(
                "<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in cells) + "</tr>"
            )
        elif row.startswith("|--"):
            continue
        elif row.strip() == "":
            if in_table:
                target.append("</table>")
                in_table = False
        else:
            target.append(f"<p>{row}</p>")

    if in_table:
        (detail_parts if seen_h2 else summary_parts).append("</table>")

    return "".join(summary_parts), "".join(detail_parts)


def report_html(stats):
    md = report_markdown(stats)
    summary_html, details_html = _md_to_html(md)
    gantt = _build_gantt_svg(stats)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Snakemake Execution Profile</title>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 1200px;
         margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }}
  h1 {{ border-bottom: 2px solid #2563eb; padding-bottom: 0.3rem; }}
  h2 {{ color: #2563eb; margin-top: 2rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: 0.5rem 0 1.5rem; }}
  th {{ background: #2563eb; color: white; text-align: left; padding: 0.5rem 0.75rem; }}
  td {{ padding: 0.4rem 0.75rem; border-bottom: 1px solid #e5e7eb; }}
  tr:hover td {{ background: #f0f4ff; }}
  p {{ line-height: 1.6; }}
  .gantt {{ margin: 1rem 0; overflow-x: auto; }}
  svg text {{ font-size: 11px; font-family: monospace; }}
</style>
</head>
<body>
{summary_html}
<h2>Execution Timeline</h2>
<div class="gantt">{gantt}</div>
{details_html}
</body>
</html>"""


def _abbreviate_rule(rule):
    """Shorten a rule name to fit inside a Gantt bar."""
    # Drop common prefixes/suffixes to keep it compact
    abbrev = rule
    for prefix in ("make_", "calling_", "perform_", "compute_", "create_",
                   "prep_", "check_", "build_", "download_", "resolve_"):
        if abbrev.startswith(prefix):
            abbrev = abbrev[len(prefix):]
            break
    # Truncate long names
    if len(abbrev) > 20:
        abbrev = abbrev[:18] + ".."
    return abbrev


def _text_color_for_bg(hex_color):
    """Return white or black text depending on background luminance."""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    # Relative luminance (sRGB approximation)
    lum = 0.299 * r + 0.587 * g + 0.114 * b
    return "#ffffff" if lum < 140 else "#000000"


def _build_gantt_svg(stats):
    """Build a simple SVG Gantt chart of job execution, colored by phase."""
    jobs = []
    for rule, data in stats["rules_ranked"]:
        for j in data["jobs"]:
            jobs.append(j)

    if not jobs:
        return "<p>No job data available.</p>"

    t0 = stats["total_start"]
    wall = stats["wall_time"]
    if wall == 0:
        return "<p>Zero wall time.</p>"

    # Collect phases and rules-per-phase (in order of first appearance)
    phase_rules = {}  # phase -> list of rule names (unique, ordered)
    for j in sorted(jobs, key=lambda x: x["start"]):
        phase = _rule_to_phase(j["rule"])
        if phase not in phase_rules:
            phase_rules[phase] = []
        if j["rule"] not in phase_rules[phase]:
            phase_rules[phase].append(j["rule"])

    # Layout — SVG uses viewBox for responsive scaling to container width
    row_height = 18
    left_margin = 10
    chart_width = 1160
    svg_width = chart_width + left_margin + 10

    # Greedy lane assignment
    lanes = []
    job_rows = []
    sorted_jobs = sorted(jobs, key=lambda j: j["start"])
    for j in sorted_jobs:
        s = (j["start"] - t0).total_seconds()
        placed = False
        for li, lane_end in enumerate(lanes):
            if s >= lane_end:
                lanes[li] = (j["end"] - t0).total_seconds()
                job_rows.append((j, li))
                placed = True
                break
        if not placed:
            lanes.append((j["end"] - t0).total_seconds())
            job_rows.append((j, len(lanes) - 1))

    num_lanes = len(lanes)
    svg_height = max(num_lanes * row_height + 40, 80)

    parts = [f'<svg width="100%" viewBox="0 0 {svg_width} {svg_height}" '
             f'xmlns="http://www.w3.org/2000/svg">']

    # Time axis ticks
    for pct in range(0, 101, 25):
        x = left_margin + (pct / 100) * chart_width
        t_sec = pct / 100 * wall
        parts.append(f'<line x1="{x}" y1="0" x2="{x}" y2="{svg_height - 20}" '
                     f'stroke="#e5e7eb" stroke-width="1"/>')
        parts.append(f'<text x="{x}" y="{svg_height - 5}" text-anchor="middle" '
                     f'fill="#6b7280">{fmt_duration(t_sec)}</text>')

    # Job bars — colored by phase, with rule abbreviation inside wider bars
    char_width = 6.6  # approximate monospace character width at 11px
    for j, lane in job_rows:
        s = (j["start"] - t0).total_seconds()
        e = (j["end"] - t0).total_seconds()
        x = left_margin + (s / wall) * chart_width
        w = max(((e - s) / wall) * chart_width, 2)
        y = lane * row_height + 2
        h = row_height - 4
        phase = _rule_to_phase(j["rule"])
        color = PHASE_PALETTE.get(phase, PHASE_PALETTE["Other"])
        wc = j["wildcards"][:40] if j["wildcards"] else ""
        title = f'{j["rule"]} [{phase}] ({fmt_duration(j["duration"])})\n{wc}'
        parts.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{h}" '
                     f'fill="{color}" rx="2" opacity="0.85">'
                     f'<title>{title}</title></rect>')

        # Add rule abbreviation inside bars wide enough to fit text
        abbrev = _abbreviate_rule(j["rule"])
        text_width = len(abbrev) * char_width
        if w > text_width + 4:
            txt_color = _text_color_for_bg(color)
            tx = x + 3
            ty = y + h - 3
            parts.append(f'<text x="{tx:.1f}" y="{ty}" fill="{txt_color}" '
                         f'font-size="10px" font-family="monospace">'
                         f'<title>{title}</title>{abbrev}</text>')

    # Legend: phase colors with constituent rules listed
    legend_y_start = svg_height
    legend_lines = []
    line_height = 16
    for phase, rules in phase_rules.items():
        color = PHASE_PALETTE.get(phase, PHASE_PALETTE["Other"])
        # Phase header
        legend_lines.append(("phase", phase, color, None))
        # Rules under this phase (compact, multiple per row)
        rule_str = ", ".join(rules)
        # Wrap at ~100 chars
        while rule_str:
            chunk = rule_str[:100]
            if len(rule_str) > 100:
                # Break at last comma before 100 chars
                idx = chunk.rfind(",")
                if idx > 0:
                    chunk = chunk[:idx + 1]
            legend_lines.append(("rules", chunk.strip(), None, None))
            rule_str = rule_str[len(chunk):].lstrip(", ")

    legend_height = len(legend_lines) * line_height + 20
    total_height = svg_height + legend_height
    parts[0] = parts[0].replace(f'0 0 {svg_width} {svg_height}',
                                f'0 0 {svg_width} {total_height}')

    ly = legend_y_start + 16
    for kind, text, color, _ in legend_lines:
        if kind == "phase":
            parts.append(f'<rect x="10" y="{ly - 10}" width="12" height="12" '
                         f'fill="{color}" rx="2"/>')
            parts.append(f'<text x="28" y="{ly}" fill="#1a1a1a" '
                         f'font-weight="bold" font-size="12px">{text}</text>')
        else:
            parts.append(f'<text x="34" y="{ly}" fill="#555" '
                         f'font-size="10px" font-family="monospace">{text}</text>')
        ly += line_height

    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Multi-section HTML report
# ---------------------------------------------------------------------------

def report_multi_html(sections):
    """Generate a tabbed HTML report from multiple (label, stats) pairs."""
    tab_buttons = []
    tab_panels = []
    for i, (label, stats) in enumerate(sections):
        tab_id = f"tab{i}"
        active_cls = " active" if i == 0 else ""
        tab_buttons.append(f'<button class="tab-btn{active_cls}" '
                           f'onclick="showTab(\'{tab_id}\')">{label}</button>')

        md = report_markdown(stats)
        summary_html, details_html = _md_to_html(md)
        gantt = _build_gantt_svg(stats)
        display = "block" if i == 0 else "none"
        tab_panels.append(
            f'<div id="{tab_id}" class="tab-panel" style="display:{display}">'
            f'{summary_html}'
            f'<h2>Execution Timeline</h2><div class="gantt">{gantt}</div>'
            f'{details_html}'
            f'</div>'
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>EPICC Pipeline Execution Profiles</title>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 1200px;
         margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }}
  h1 {{ border-bottom: 2px solid #2563eb; padding-bottom: 0.3rem; }}
  h2 {{ color: #2563eb; margin-top: 2rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: 0.5rem 0 1.5rem; }}
  th {{ background: #2563eb; color: white; text-align: left; padding: 0.5rem 0.75rem; }}
  td {{ padding: 0.4rem 0.75rem; border-bottom: 1px solid #e5e7eb; }}
  tr:hover td {{ background: #f0f4ff; }}
  p {{ line-height: 1.6; }}
  .gantt {{ margin: 1rem 0; overflow-x: auto; }}
  svg text {{ font-size: 11px; font-family: monospace; }}
  .tab-bar {{ display: flex; gap: 0; border-bottom: 2px solid #2563eb; margin-bottom: 1.5rem; }}
  .tab-btn {{ background: #f0f4ff; border: 1px solid #d1d5db; border-bottom: none;
              padding: 0.6rem 1.2rem; cursor: pointer; font-size: 0.95rem;
              border-radius: 6px 6px 0 0; margin-bottom: -2px; }}
  .tab-btn.active {{ background: white; border-color: #2563eb; border-bottom: 2px solid white;
                     font-weight: 600; color: #2563eb; }}
  .tab-btn:hover:not(.active) {{ background: #e0e7ff; }}
</style>
<script>
function showTab(id) {{
  document.querySelectorAll('.tab-panel').forEach(p => p.style.display = 'none');
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(id).style.display = 'block';
  event.target.classList.add('active');
}}
</script>
</head>
<body>
<div class="tab-bar">{"".join(tab_buttons)}</div>
{"".join(tab_panels)}
</body>
</html>"""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def find_latest_log():
    log_dir = Path(".snakemake/log")
    if not log_dir.exists():
        sys.exit("No .snakemake/log directory found. Run from the Snakemake project root.")
    logs = sorted(log_dir.glob("*.snakemake.log"))
    if not logs:
        sys.exit("No Snakemake log files found.")
    return logs[-1]


# Output-dir capture: the path component appearing before /<env>/, where <env>
# is one of EPICC's per-datatype results subdirs.
_RESULT_PREFIX_RE = re.compile(
    r"(?:^|[\s'\"])(\S+?)/(?:ChIP|ATAC|RNA|sRNA|mC|combined)/"
)
_ANALYSIS_NAME_RE = re.compile(r"analysis_name=([^,\s]+)")


def extract_run_signature(path, max_hits=200):
    """Identify the run a log belongs to, from log content.

    Returns (output_dir, frozenset_of_analysis_names) or None for an empty log.
    output_dir is the most common path prefix appearing before an EPICC results
    subdirectory; analysis_names is the set of 'analysis_name=...' wildcards.
    """
    prefix_counts = Counter()
    analysis_names = set()
    try:
        with open(path) as fh:
            for line in fh:
                for m in _RESULT_PREFIX_RE.finditer(line):
                    prefix_counts[m.group(1)] += 1
                m = _ANALYSIS_NAME_RE.search(line)
                if m:
                    analysis_names.add(m.group(1))
                if sum(prefix_counts.values()) >= max_hits:
                    break
    except OSError:
        return None
    if not prefix_counts:
        return None
    output_dir = prefix_counts.most_common(1)[0][0]
    return (output_dir, frozenset(analysis_names))


def find_run_logs(output_dir=None):
    """Return (logs_in_run, signature) for the most recently active resumed run.

    Walks logs newest-first to find the most recent log with an extractable
    signature, then collects all logs sharing the same output_dir and either
    overlapping or empty analysis_name sets. Falls back to the newest log alone
    if no signature can be extracted.

    If ``output_dir`` is given, anchor on that output_dir instead of the most
    recent run: collect every log whose captured output_dir matches it. This is
    how a run is selected by its results directory (``epicc perf RESULTS_DIR``).
    """
    log_dir = Path(".snakemake/log")
    if not log_dir.exists():
        sys.exit("No .snakemake/log directory found. Run from the Snakemake project root.")
    logs = sorted(log_dir.glob("*.snakemake.log"))
    if not logs:
        sys.exit("No Snakemake log files found.")

    if output_dir is not None:
        want = str(Path(output_dir))
        matched = []
        names_seen = set()
        available = set()
        for p in logs:
            sig = extract_run_signature(p)
            if sig is None:
                continue
            out_dir, names = sig
            available.add(out_dir)
            if str(Path(out_dir)) != want:
                continue
            matched.append(p)
            names_seen |= names
        if not matched:
            hint = ("\nOutput dirs seen in .snakemake/log: "
                    + ", ".join(sorted(available))) if available else ""
            sys.exit(f"No logs found for output_dir '{output_dir}'.{hint}")
        return matched, (want, frozenset(names_seen))

    target_sig = None
    for p in reversed(logs):
        sig = extract_run_signature(p)
        if sig is not None:
            target_sig = sig
            break
    if target_sig is None:
        return [logs[-1]], None

    target_dir, target_names = target_sig
    matched = []
    for p in logs:
        sig = extract_run_signature(p)
        if sig is None:
            continue
        out_dir, names = sig
        if out_dir != target_dir:
            continue
        # If both logs have analysis_names, require overlap; otherwise accept.
        if target_names and names and target_names.isdisjoint(names):
            continue
        matched.append(p)
    return matched, target_sig


def main():
    parser = argparse.ArgumentParser(description="Profile a Snakemake run from its log file.")
    parser.add_argument("logfile", nargs="?",
                        help="Path to a specific .snakemake.log file, or a run's "
                             "results directory (output_dir) to aggregate all of "
                             "that run's logs")
    parser.add_argument("--latest", action="store_true",
                        help="(Default) Aggregate logs from the most recent resumed run")
    parser.add_argument("--single", action="store_true",
                        help="Profile only the newest single log, without aggregating")
    parser.add_argument("--html", metavar="FILE", help="Write HTML report to FILE")
    parser.add_argument("--multi", nargs="+", metavar="LABEL=LOG",
                        help="Generate multi-section HTML: 'Label=path/to/log' ...")
    args = parser.parse_args()

    if args.multi:
        if not args.html:
            sys.exit("--multi requires --html <output>")
        sections = []
        for spec in args.multi:
            if "=" not in spec:
                sys.exit(f"--multi entries must be LABEL=LOGFILE, got: {spec}")
            label, logpath = spec.split("=", 1)
            logpath = Path(logpath)
            if not logpath.exists():
                sys.exit(f"Log file not found: {logpath}")
            print(f"Parsing {logpath} ({label}) ...", file=sys.stderr)
            jobs = parse_log(logpath)
            if not jobs:
                print(f"  Warning: no completed jobs in {logpath}", file=sys.stderr)
                continue
            sections.append((label, analyze(jobs)))
        html = report_multi_html(sections)
        Path(args.html).write_text(html)
        print(f"Multi-section HTML report written to {args.html}", file=sys.stderr)
        return

    signature = None
    if args.logfile:
        logpath = Path(args.logfile)
        if not logpath.exists():
            sys.exit(f"Log file not found: {logpath}")
        if logpath.is_dir():
            # A run's results directory (output_dir) was given — select every
            # log belonging to that run rather than treating it as a logfile.
            logs_to_parse, signature = find_run_logs(output_dir=logpath)
        else:
            logs_to_parse = [logpath]
    elif args.single:
        logs_to_parse = [find_latest_log()]
    else:
        logs_to_parse, signature = find_run_logs()

    if len(logs_to_parse) == 1:
        print(f"Parsing {logs_to_parse[0]} ...", file=sys.stderr)
    else:
        out_dir = signature[0] if signature else "?"
        names = sorted(signature[1]) if signature and signature[1] else []
        names_label = f", analysis_name={'|'.join(names)}" if names else ""
        print(f"Aggregating {len(logs_to_parse)} logs from resumed run "
              f"(output_dir={out_dir}{names_label}):", file=sys.stderr)
        for p in logs_to_parse:
            print(f"  {p.name}", file=sys.stderr)

    log_jobs = [parse_log(p) for p in logs_to_parse]
    all_jobs = [j for jobs in log_jobs for j in jobs]

    if not all_jobs:
        sys.exit("No completed jobs found in log(s).")

    stats = analyze(all_jobs)

    # When aggregating, replace wall_time with the sum of per-log spans so
    # idle gaps between manual resumptions don't dilute parallelism.
    if len(logs_to_parse) > 1:
        active = sum(
            (max(j["end"] for j in jobs) - min(j["start"] for j in jobs)).total_seconds()
            for jobs in log_jobs if jobs
        )
        if active > 0:
            stats["wall_time"] = active
        stats["num_logs"] = len(logs_to_parse)

    if args.html:
        html = report_html(stats)
        Path(args.html).write_text(html)
        print(f"HTML report written to {args.html}", file=sys.stderr)
    else:
        print(report_markdown(stats))


if __name__ == "__main__":
    main()
