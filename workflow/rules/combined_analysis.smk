
def return_log_combined(analysis_name, genome, types):
    return os.path.join(REPO_FOLDER, RESULTS_DIR,"combined","logs",f"tmp__{analysis_name}__{genome}__{types}.log")

def define_combined_target_file(wildcards):
    heatname = config['heatmap_target_file_label']
    browsername = config['browser_target_file_label']
    target_name = wildcards.target_name
    ref_genome = wildcards.ref_genome
    
    if target_name == heatname:
        return config['heatmap_target_file']
    elif target_name == browsername:
        return config['browser_target_file']
    elif target_name == "full_chromosomes":
        return f"{RESULTS_DIR}/combined/bedfiles/full_chromosomes__{ref_genome}.bed"
    elif target_name.startswith("combined_peaks"):
        file = f"{RESULTS_DIR}/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("combined_clusters"):
        file = f"{RESULTS_DIR}/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("combined_TSS"):
        file = f"{RESULTS_DIR}/combined/bedfiles/{target_name}__{ref_genome}.bed"
    elif target_name.startswith("all_genes") or target_name.startswith("protein_coding_genes"):
        file = f"{RESULTS_DIR}/combined/bedfiles/{ref_genome}__{target_name}.bed"
    elif te_suffix_from_target(target_name) is not None:
        # all_TEs -> <G>__TE_file.bed, all_TEs_<tag> -> <G>__TE_file_<tag>.bed
        file = te_bed_path(GENOMES_DIR, ref_genome, te_suffix_from_target(target_name))
    else:
        raise ValueError(   
            f"{target_name} does not match possible files. It can be 'combined_peaks', 'combined_clusters', 'all_genes', 'all_TEs'" 
            "'full_chromosomes' or the value of 'heatmap_target_file_label' or 'browser_target_file_label' in the options file"
        )
    
    return file

def get_heatmap_param(matrix, key):
    override = config.get(key)
    if override is not None:
        return override

    return config['heatmaps'][matrix][key]

def get_matrix_inputs(wildcards):
    stranded_heatmaps = config['stranded_heatmaps']
    bedfile = define_combined_target_file(wildcards)
    prefix = f"{RESULTS_DIR}/combined/matrix/matrix_{wildcards.matrix_param}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}__{wildcards.target_name}"
    with checkpoints.is_stranded.get(bedfile=bedfile).output[0].open() as f:
        verdict = f.read().strip()
    if not stranded_heatmaps or verdict == "unstranded":
        return [ f"{prefix}__unstranded.gz" ]
    if verdict == "stranded_mixed":
        # The '.'/'?' rows become their own region group -- see is_stranded.
        return [ f"{prefix}__plus.gz", f"{prefix}__minus.gz", f"{prefix}__nostrand.gz" ]
    return [ f"{prefix}__plus.gz", f"{prefix}__minus.gz" ]

def define_sort_options(wildcards):
    sort_options = config['heatmaps_sort_options']
    if sort_options == "no":
        return "--sortRegions keep"
    elif sort_options == "mean":
        return "--sortRegions descend --sortUsing mean"
    elif sort_options == "median":
        return "--sortRegions descend --sortUsing median"
    else:
        print("unclear option: no sorting done")
        return "--sortRegions keep"

def define_samples_for_upset(wildcards, string):
    # One entry per UpSet sample: its label components (build_plot_labels),
    # file and type. type is None for sRNA, whose plotted columns are the label
    # per size class.
    entries = []
    types = set()
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    allreps = config['upset_allreps']
    if globenv == "all_chip":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'].isin(["ChIP","ATAC"])) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    elif globenv == "RAMPAGE":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['sample_type'] == "RAMPAGE") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        if row.env in ["ChIP", "ATAC"]:
            if allreps:
                for sid in get_replicate_sample_ids(row['sample_name'], samples):
                    peaktype = get_peaktype_for_env(row['sample_name'], row.env)
                    paired = get_sample_info_from_name(row['sample_name'], analysis_samples, 'paired')
                    if row.env == "ATAC":
                        prefix = "peaks_atac"
                    elif paired == "PE":
                        prefix = "peaks_pe"
                    else:
                        prefix = "peaks_se"
                    file = f"{RESULTS_DIR}/{row.env}/peaks/{prefix}__final__{sid}_peaks.{peaktype}Peak"
                    rep = parse_sample_name(sid)['replicate']
                    entries.append(dict(levels=row.Levels, mark=row.sample_type, assay=row.Assay, replicate=rep,
                                        sid=sid, file=file, type=row.sample_type))
                    types.add(row.sample_type)
            else:
                file = f"{RESULTS_DIR}/{row.env}/peaks/selected_peaks__{row['mapped_name']}.bedPeak"
                entries.append(dict(levels=row.Levels, mark=row.sample_type, assay=row.Assay,
                                    sid=row['sample_name'], file=file, type=row.sample_type))
                types.add(row.sample_type)
        elif globenv == "RAMPAGE":
            if allreps:
                for sid in get_replicate_sample_ids(row['sample_name'], samples):
                    file = f"{RESULTS_DIR}/RNA/TSS/TSS__final__{sid}_peaks.narrowPeak"
                    rep = parse_sample_name(sid)['replicate']
                    entries.append(dict(levels=row.Levels, mark=row.sample_type, assay=row.Assay, replicate=rep,
                                        sid=sid, file=file, type=row.levels_label))
                    types.add(f"{row.levels_label}")
            else:
                file = f"{RESULTS_DIR}/RNA/TSS/TSS__merged__{row['sample_name']}_peaks.narrowPeak"
                entries.append(dict(levels=row.Levels, mark=row.sample_type, assay=row.Assay,
                                    sid=row['sample_name'], file=file, type=row.levels_label))
                types.add(f"{row.levels_label}")
        elif row.env == "sRNA":
            for sid in get_replicate_sample_ids(row['sample_name'], samples):
                file = f"{RESULTS_DIR}/sRNA/mapped/{sid}/clusters.bed"
                rep = parse_sample_name(sid)['replicate']
                entries.append(dict(levels=row.Levels, mark=row.sample_type, assay=row.Assay, replicate=rep,
                                    sid=sid, file=file, type=None))

    labels = build_plot_labels(entries)
    check_unique_labels(labels, f"the {globenv} UpSet plot")
    files = [e["file"] for e in entries]
    # The merged peak/cluster/TSS file names each sample by its sample id, not
    # its label, so the expensive merge never depends on how labels are
    # spelled; plotting_upset_regions swaps the ids for labels ('relabel').
    keys = [e["sid"] for e in entries]
    names = [f"{key}:{file}" for key, file in zip(keys, files)]
    if globenv == "sRNA":
        srna_min = config['srna_min_size']
        srna_max = config['srna_max_size']
        types = [f"{s}nt" for s in range(srna_min, srna_max + 1)]
        types += ["MIRNA", "Others"]
        ordered = ":".join(types)
        # combine_clusterfiles suffixes each label with the cluster's size class
        # (clusters.bed column 4), so the plotted column is '{label}_{type}'.
        label_types = {f"{label}_{t}": t for label in labels for t in types}
    else:
        ordered = ":".join(sorted(types))
        label_types = {label: e["type"] for label, e in zip(labels, entries)}

    if string == "pairs":
        return names
    elif string == "labels":
        return labels
    elif string == "files":
        return files
    elif string == "types":
        return ordered
    elif string == "label_types":
        # Explicit label->type map for the UpSet scripts: a type can be a
        # substring of another (H3K9me / H3K9me2), so the scripts can't
        # recover it by matching on the label text.
        return ",".join(f"{label}={t}" for label, t in label_types.items())
    elif string == "relabel":
        # "id=label,id=label" for the Samples column of the merged file
        # (recode_sample_ids in upset_common.R; sRNA's '{id}_{class}' columns
        # are resolved there).
        return ",".join(f"{key}={label}" for key, label in zip(keys, labels))

def define_upset_script(wildcards):
    globenv = wildcards.env
    if globenv in ["all_chip", "ChIP", "ATAC"]:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot_peaks.R")
    elif globenv == "sRNA":
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot_clusters.R")
    elif globenv == "RAMPAGE":
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_Upset_plot_TSS.R")

    return script

import colorsys
from types import SimpleNamespace

# Matplotlib tab20 and Set2 palettes as hex — avoids runtime matplotlib dependency
_PALETTES = {
    "tab20": [
        "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
        "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
        "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
        "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
    ],
    "Set2": [
        "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
        "#ffd92f", "#e5c494", "#b3b3b3",
    ],
}

def assign_colors(keys, cmap_name="tab20"):
    pal = _PALETTES[cmap_name]
    colors = {}
    for i, key in enumerate(sorted(set(keys))):
        colors[key] = pal[i % len(pal)]
    return colors

def make_it_lighter(palette, factor):
    new_palette = {}
    for k, c in palette.items():
        r, g, b = int(c[1:3], 16) / 255, int(c[3:5], 16) / 255, int(c[5:7], 16) / 255
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        l = min(0.9, l * factor)
        rr, gg, bb = colorsys.hls_to_rgb(h, l, s)
        new_palette[k] = f"#{int(rr*255+.5):02x}{int(gg*255+.5):02x}{int(bb*255+.5):02x}"
    return new_palette

def define_key_for_plots(wildcards, string):
    if string == "display":
        # Display rows for the heatmap and profile plots, "key<TAB>label<TAB>
        # colour<TAB>track<TAB>series" per matrix column, looked up by key at plot time
        # (matrix_sample_labels.py). The plot rules don't know whether the
        # matrix came from the unstranded pass or the plus/minus/nostrand
        # passes, which carry different columns, so both are listed; keys they
        # share map to the same row. One pass per strand layout builds every
        # column at once.
        rows = {}
        for s in ("unstranded", "plus"):
            c = _plot_columns(SimpleNamespace(ref_genome=wildcards.ref_genome, env=wildcards.env, strand=s))
            for key, *rest in zip(c["keys"], c["labels"], c["colors"], c["tracks"], c["series"]):
                rows.setdefault(key, rest)
        return "\n".join("\t".join([key] + rest) for key, rest in rows.items())
    c = _plot_columns(wildcards)
    if string == "table":
        # What prep_browser_on_region extracts: no labels or colours, so the
        # extraction never reruns over a label change.
        table_name = f"{RESULTS_DIR}/combined/matrix/sample_table__{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}.tab"
        os.makedirs(os.path.dirname(table_name), exist_ok=True)
        tab = pd.DataFrame({
            "bigwigs": c["bigwigs"],
            "keys": c["keys"],
            "marks": c["tracks"],
        })
        tab.to_csv(table_name, sep="\t", index=False, header=False)
        return table_name
    if string == "browser_display":
        # "key<TAB>label<TAB>back<TAB>track<TAB>fill+<TAB>fill-" per browser
        # track, joined onto the extracted tracks by make_single_loci_browser_plot.
        return "\n".join("\t".join(r) for r in zip(c["keys"], c["labels"], c["colors"], c["trackcolors"],
                                                      c["fillcolorsplus"], c["fillcolorsminus"]))
    return c[string]

def _plot_columns(wildcards):
    """Every per-column list of one plot, in column order, built in one pass:
    bigwigs, keys, labels, series, tracks (each column's mark: mCG, H3K9me2,
    RNAseq_plus, ...), colors (tab20 keyed on levels_label, the browser
    backgrounds too, so a genotype keeps one colour everywhere -- deeptools 4's
    default line colours sample a sequential colormap and leave consecutive
    samples near-identical), and the browser's track/fill colours."""
    unique_chip = set()
    unique_atac = set()
    unique_rna = set()
    unique_srna = set()
    unique_mc = set()
    # One entry per plotted track, grouped by mark. Labels are made unique only
    # once every entry is in, so everything keyed by label is built afterwards
    # from these entries, in the same order as the bigwigs.
    grouped = defaultdict(list)
    srna_sizes = config['srna_heatmap_sizes']
    plot_allreps = config['plot_allreps']
    mc_contexts = get_methylation_contexts()
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    strand = getattr(wildcards, "strand", "unstranded")
    # The 'nostrand' pass covers regions whose strand is '.' or '?', so sense is
    # undefined for them. rbind requires every matrix it merges to carry the same
    # samples in the same order, so this pass cannot emit both strands the way
    # 'unstranded' does -- it takes the plus track, the coordinate-forward one,
    # which is also the orientation computeMatrix gives an unstranded region.
    if strand == "nostrand":
        strand = "plus"

    def add(group, row, sid, bw, mark, rep="", track=None, suffix=""):
        track = track or mark
        # The key names a column in every data file (matrix, browser scratch
        # tracks) so they never depend on the label. Its last '__' field is the
        # track, which computing_matrix_scales reads back, and a minus-strand
        # key ends in '_minus', which prep_browser_on_region negates on.
        grouped[group].append({"bw": bw, "levels": row.Levels, "mark": mark,
                               "assay": row.Assay, "replicate": rep, "sid": sid,
                               "suffix": suffix, "track": track, "key": f"{sid}__{track}",
                               "type": row.levels_label})

    def add_pair(group, row, sid, bw1, bw2, mark, rep=""):
        add(group, row, sid, bw1, mark, rep, f"{mark}_plus", "_plus")
        add(group, row, sid, bw2, mark, rep, f"{mark}_minus", "_minus")

    if globenv == "all":
        filtered_analysis_samples = analysis_samples[ analysis_samples['ref_genome'] == ref_genome ].copy()
    elif globenv == "most":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] != "mC") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        name = row['sample_name']
        rep_ids = get_replicate_sample_ids(name, samples)
        if row.env == "ChIP":
            unique_chip.add(row.sample_type)
            if not plot_allreps:
                bw = f"{RESULTS_DIR}/{row.env}/tracks/FC__merged__{row['mapped_name']}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/FC__final__{rep_ids[0]}.bw"
                add(f"chip_{row.sample_type}", row, name, bw, row.sample_type)
            else:
                for sid in rep_ids:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/FC__final__{sid}.bw"
                    rep = parse_sample_name(sid)['replicate']
                    add(f"chip_{row.sample_type}", row, sid, bw, row.sample_type, rep)

        elif row.env == "RNA":
            unique_rna.add(row.data_type)
            strandedness = config['rna_tracks'][row.data_type]['strandedness']
            if strandedness == "unstranded":
                if not plot_allreps:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__unstranded.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__unstranded.bw"
                    add(row.data_type, row, name, bw, row.data_type)
                else:
                    for sid in rep_ids:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__unstranded.bw"
                        rep = parse_sample_name(sid)['replicate']
                        add(row.data_type, row, sid, bw, row.data_type, rep)
            elif strand == "unstranded":
                if not plot_allreps:
                    bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__plus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__plus.bw"
                    bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__minus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__minus.bw"
                    add_pair(row.data_type, row, name, bw1, bw2, row.data_type)
                else:
                    for sid in rep_ids:
                        bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__plus.bw"
                        bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__minus.bw"
                        rep = parse_sample_name(sid)['replicate']
                        add_pair(row.data_type, row, sid, bw1, bw2, row.data_type, rep)
            else:
                if not plot_allreps:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{strand}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{strand}.bw"
                    add(row.data_type, row, name, bw, row.data_type)
                else:
                    for sid in rep_ids:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{strand}.bw"
                        rep = parse_sample_name(sid)['replicate']
                        add(row.data_type, row, sid, bw, row.data_type, rep)

        elif row.env == "sRNA":
            for size in srna_sizes:
                mark = f"sRNA_{size}nt"
                unique_srna.add(mark)
                if strand == "unstranded":
                    if not plot_allreps:
                        bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{size}nt__plus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{size}nt__plus.bw"
                        bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{size}nt__minus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{size}nt__minus.bw"
                        add_pair(mark, row, name, bw1, bw2, mark)
                    else:
                        for sid in rep_ids:
                            bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{size}nt__plus.bw"
                            bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{size}nt__minus.bw"
                            rep = parse_sample_name(sid)['replicate']
                            add_pair(mark, row, sid, bw1, bw2, mark, rep)
                else:
                    if not plot_allreps:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{size}nt__{strand}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{size}nt__{strand}.bw"
                        add(mark, row, name, bw, mark)
                    else:
                        for sid in rep_ids:
                            bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{size}nt__{strand}.bw"
                            rep = parse_sample_name(sid)['replicate']
                            add(mark, row, sid, bw, mark, rep)

        elif row.env == "mC":
            if not plot_allreps:
                for context in mc_contexts:
                    unique_mc.add(f"m{context}")
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{context}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{context}.bw"
                    add(f"m{context}", row, name, bw, f"m{context}")
            else:
                for sid in rep_ids:
                    for context in mc_contexts:
                        unique_mc.add(f"m{context}")
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{context}.bw"
                        rep = parse_sample_name(sid)['replicate']
                        add(f"m{context}", row, sid, bw, f"m{context}", rep)

        elif row.env == "ATAC":
            unique_atac.add("ATAC")
            if not plot_allreps:
                bw = f"{RESULTS_DIR}/ATAC/tracks/coverage__merged__{row['mapped_name']}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/ATAC/tracks/coverage__final__{rep_ids[0]}.bw"
                add("atac", row, name, bw, "ATAC")
            else:
                for sid in rep_ids:
                    bw = f"{RESULTS_DIR}/ATAC/tracks/coverage__final__{sid}.bw"
                    rep = parse_sample_name(sid)['replicate']
                    add("atac", row, sid, bw, "ATAC", rep)

    order = ( [f"chip_{chip}" for chip in sorted(unique_chip)] + ["atac"] + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )
    entries = sum([grouped.get(group, []) for group in order], [])
    bigwigs = [e["bw"] for e in entries]
    keys = [e["key"] for e in entries]
    check_unique_labels(keys, f"the {globenv} plot keys")
    # A label shows each component that varies across this plot's samples
    # (build_plot_labels). Colour stays keyed on the full levels_label, so a
    # genotype keeps one colour across plots whatever its label shows.
    labels = build_plot_labels(entries)
    # The line name across the methylation-context panels of one plot
    # (plot_profile_by_context in plot_profile.R): the label without its mark.
    series = build_plot_labels(entries, omit=("mark",))
    check_unique_labels(labels, f"the {globenv} plots")
    marksforbrowser = ( sorted(unique_chip) + sorted(unique_atac) + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )

    types = sorted(filtered_analysis_samples["levels_label"].tolist())

    back_palette = assign_colors(types, "tab20")
    track_palette = assign_colors(marksforbrowser, "Set2")
    plus_palette = make_it_lighter(track_palette, 1.1)
    minus_palette = make_it_lighter(track_palette, 1.4)
    for m in unique_rna | unique_srna:
        minus_palette[m] = plus_palette[m]

    backcolors = [back_palette[e["type"]] for e in entries]
    trackcolors = [track_palette[e["mark"]] for e in entries]
    fillcolorsplus = [plus_palette[e["mark"]] for e in entries]
    fillcolorsminus = [minus_palette[e["mark"]] for e in entries]
    alignedmarks = [e["track"] for e in entries]

    return {
        "bigwigs": bigwigs, "keys": keys, "labels": labels, "series": series,
        "tracks": alignedmarks, "colors": backcolors, "trackcolors": trackcolors,
        "fillcolorsplus": fillcolorsplus, "fillcolorsminus": fillcolorsminus,
    }

def define_individual_browser_plots(wildcards):
    ref_genome = wildcards.ref_genome
    analysis_name = wildcards.analysis_name    
    env = wildcards.env
    target_name = wildcards.target_name
    target_file = define_combined_target_file(wildcards)
    
    checkpoints.is_stranded.get(bedfile=target_file)
    
    with open(target_file) as f:
        first_line = f.readline().strip().split("\t")
    try:
        header = "no" if (int(first_line[1]) >=0 and int(first_line[2]) >=0) else "yes"
    except (ValueError, IndexError):
        header = "yes"

    regions = pd.read_csv(target_file, sep="\t", header=0 if header=="yes" else None)
    
    files = []
    starter = 2 if header == "yes" else 1
    for row_number, _ in enumerate(regions.itertuples(index=False),start=starter):
        files.append(f"{RESULTS_DIR}/combined/plots/single_browser__{target_name}__line{row_number}__{env}__{analysis_name}__{ref_genome}.pdf")    
    return files

# User BED files drawn as extra annotation rows on every browser plot
# (browser_extra_tracks; validated at startup by check_extra_output_files).
from scripts.browser_extra_tracks import normalize_extra_tracks
BROWSER_EXTRA_TRACKS = [
    dict(t, file=resolve_repo_data(t["file"]))
    for t in normalize_extra_tracks(config.get("browser_extra_tracks"))
]

def get_TE_file(wildcards):
    ref_genome = wildcards.ref_genome
    if config['browser_TE_file'] and te_source_for(ref_genome):
        return te_bed_for(ref_genome)
    return []

def define_input_for_pca(wildcards, string):
    # `indexes` returns the .bam.bai companions of `tracks` for the BAM-based
    # envs (empty for the mC contexts, which summarize bigwigs and need no
    # index). summarize_tracks_pca declares them as inputs so the index is a
    # tracked dependency: multiBamSummary requires it, but without the edge
    # Snakemake neither rebuilds a missing index nor counts the PCA job as a
    # consumer — so a temp() index can be reclaimed while the BAM it belongs to
    # is still alive. See the note on that rule.
    tracks = []
    indexes = []
    # keys name the samples in the summary matrix; labels are applied only by
    # plotPCA, so a label change reruns the plot, not the summary.
    keys = []
    entries = []
    groups = []
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env

    if globenv in ["mCG", "mCHG", "mCHH"]:
        filtered_samples = samples[ (samples['env'] == "mC") & (samples['ref_genome'] == ref_genome) ]
        context = globenv[1:]
        for _, row in filtered_samples.iterrows():
            tracks.append(f"{RESULTS_DIR}/mC/tracks/{row['mapped_name']}__{context}.bw")
            keys.append(row['mapped_name'])
            entries.append(dict(levels=row.Levels, mark=globenv, assay=row.Assay,
                                replicate=row.replicate, sid=row.Sample_ID))
            groups.append(f"{row.levels_label}")
    elif globenv in ["ChIP", "ATAC", "all_chip"]:
        envs = ["ChIP", "ATAC"] if globenv == "all_chip" else [globenv]
        filtered_samples = samples[ (samples['env'].isin(envs)) & (samples['ref_genome'] == ref_genome) ]
        for _, row in filtered_samples.iterrows():
            bam = f"{RESULTS_DIR}/{row.env}/mapped/final__{row['mapped_name']}.bam"
            tracks.append(bam)
            indexes.append(bam + ".bai")
            keys.append(row['mapped_name'])
            entries.append(dict(levels=row.Levels, mark=row.sample_type, assay=row.Assay,
                                replicate=row.replicate, sid=row.Sample_ID))
            groups.append(f"{row.sample_type}_{row.levels_label}")
    # PCA points are replicates, controls included, so two libraries can share
    # Levels, mark, assay and replicate (separate Inputs for two marks); the
    # builder falls back to the Sample_ID for those.
    labels = build_plot_labels(entries)
    check_unique_labels(labels, f"the {globenv} PCA")

    palette = assign_colors(groups, "tab20")
    colors = [palette[group] for group in groups]

    if string == "tracks":
        return tracks
    elif string == "indexes":
        return indexes
    elif string == "keys":
        return keys
    elif string == "labels":
        return labels
    elif string == "colors":
        return colors

def define_final_stats_output():
    aligned_bams = config['aligned_bams']
    stat_files = []    
    if not aligned_bams:
        stat_files += expand(f"{RESULTS_DIR}/combined/plots/mapping_stats_{{analysis_name}}_{{env}}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","ATAC"]])

    stat_files += expand(f"{RESULTS_DIR}/combined/plots/mapping_stats_{{analysis_name}}_{{env}}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["RNA","mC"]])
    stat_files += expand(f"{RESULTS_DIR}/combined/plots/srna_sizes_stats_{{analysis_name}}_{{env}}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["sRNA"]])
    stat_files += expand(f"{RESULTS_DIR}/combined/plots/peak_stats_{{analysis_name}}_{{env}}.pdf", analysis_name = analysis_name, env=[env for env in UNIQUE_ENVS if env in ["ChIP","ATAC"]])

    return stat_files

def define_final_combined_output(ref_genome):
    qc_option = config["QC_option"]
    full_analysis = config['full_analysis']
    te_analysis = config['te_analysis']
    analysis_name = config['analysis_name']
    mc_sort = config['heatmap_sort_mc_after_others']
    mc_contexts = get_methylation_contexts()
    plot_files = []
    te_plots = []
    te_target = te_target_for(ref_genome)
    
    all_analysis_samples = analysis_samples[ analysis_samples['ref_genome'] == ref_genome ].copy()
    chip_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    atac_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ATAC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    mc_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'mC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    rna_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'RNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    srna_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'sRNA') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    rampage_analysis_samples = analysis_samples[ (analysis_samples['data_type'] == 'RAMPAGE') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
        
    if len(all_analysis_samples) >=1:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Browser_full_chromosomes__all__{analysis_name}__{ref_genome}.pdf")
    
    if len(chip_analysis_samples) >=2:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Upset_combined_peaks__ChIP__{analysis_name}__{ref_genome}.pdf")

    peak_envs = sum(1 for s in [chip_analysis_samples, atac_analysis_samples] if len(s) >= 1)
    if peak_envs >= 2:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Upset_combined_peaks__all_chip__{analysis_name}__{ref_genome}.pdf")
    
    if len(srna_analysis_samples) >=1:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Upset_combined_clusters__sRNA__{analysis_name}__{ref_genome}.pdf")

    if len(atac_analysis_samples) >=2:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Upset_combined_peaks__ATAC__{analysis_name}__{ref_genome}.pdf")
    
    if len(rampage_analysis_samples) >=2:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Upset_combined_TSS__RAMPAGE__{analysis_name}__{ref_genome}.pdf")

    if len(mc_analysis_samples) >=1:
        if len(all_analysis_samples) > len(mc_analysis_samples) and mc_sort:
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")            
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__regions__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__tss__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__tes__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
            
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
        else:
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
        
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__{te_target}.pdf")
        
    else:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
    
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__{te_target}.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__{te_target}.pdf")

    mc_rep_samples = samples[ (samples['env'] == "mC") & (samples['ref_genome'] == ref_genome) ].copy()
    if len(mc_rep_samples) >=3:
        for ctxt in mc_contexts:
            plot_files.append(f"{RESULTS_DIR}/combined/plots/PCA__m{ctxt}__{analysis_name}__{ref_genome}.pdf")
    
    all_chip_samples = samples[ (samples['env'].isin(["ChIP","ATAC"])) & (samples['ref_genome'] == ref_genome) ].copy()
    if len(all_chip_samples) >=3:
            plot_files.append(f"{RESULTS_DIR}/combined/plots/PCA__all_chip__{analysis_name}__{ref_genome}.pdf")

    for env in [e for e in UNIQUE_ENVS if e in ["ChIP", "ATAC"]]:
        env_rep_samples = samples[ (samples['env'] == env) & (samples['ref_genome'] == ref_genome) ].copy()
        if len(env_rep_samples) >=3:
            plot_files.append(f"{RESULTS_DIR}/combined/plots/PCA__{env}__{analysis_name}__{ref_genome}.pdf")

    results = []
    
    if full_analysis:
        results += plot_files
        
    if te_analysis:
        results += te_plots   
    
    return results

###
# rules to look for header or strandedness of bedfile
rule has_header:
    input:
        bedfile = "{bedfile}"
    output:
        file = temp("{bedfile}.header")
    localrule: True
    run:
        with open(input.bedfile) as f:
            first_line = f.readline().strip().split('\t')
            try:
                res = "no" if (int(first_line[1]) >=0 and int(first_line[2]) >=0) else "yes"
            except (ValueError, IndexError):
                res = "yes"
        
        with open(output.file, "w") as out:
            out.write(res + "\n")
              
checkpoint is_stranded:
    input:
        bedfile = "{bedfile}",
        header = "{bedfile}.header"
    output:
        file = temp("{bedfile}.stranded")
    localrule: True
    run:
        with open(input.header) as h:
            header = h.read().strip()
        
        has_strand = False
        strand_values = set()
        with open(input.bedfile) as f:
            if header == "yes":
                next(f)
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 6:
                    has_strand = True
                    strand_values.add(cols[5])

        # A file with any '+'/'-' at all is stranded: annotations routinely carry
        # a few '.' or '?' rows (GFF3 uses '?' for "stranded but unknown"), and
        # demoting the whole file on account of them costs every stranded track
        # its sense-strand orientation. 631 of ColCEN's 42,927 TEs were enough to
        # do that. Those rows go to their own 'nostrand' pass and stay a separate
        # region group downstream, since sense is undefined for them.
        oriented = strand_values & {"+", "-"}
        with open(output.file, "w") as out:
            if not has_strand or not oriented:
                out.write("unstranded\n")
            elif strand_values - {"+", "-"}:
                out.write("stranded_mixed\n")
            else:
                out.write("stranded\n")

###
# Rules to prep and then plot the mapping stats:
rule prepping_mapping_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"{RESULTS_DIR}/{wildcards.env}/reports/summary_{wildcards.env}_{get_sample_info_from_name(sample_name, samples, 'paired')}_mapping_stats_{sample_name}.txt"
                                                for sample_name in get_mapped_names_by_env(wildcards.env, samples) ]
    output:
        temp_stat_file = temp(f"{RESULTS_DIR}/combined/reports/temp_summary_mapping_stats_{{analysis_name}}_{{env}}.txt"),
        stat_file = f"{RESULTS_DIR}/combined/reports/summary_mapping_stats_{{analysis_name}}_{{env}}.txt"
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "prep_mapping_stats"))
    shell:
        """
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk -F "\t" -v OFS="\t" 'NR>1 {{print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_mapping_stats:
    input:
        summary_stats = f"{RESULTS_DIR}/combined/reports/summary_mapping_stats_{{analysis_name}}_{{env}}.txt"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/mapping_stats_{{analysis_name}}_{{env}}.pdf"
    params:
        analysis_name = lambda wildcards: wildcards.analysis_name,
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_mapping_stats.R")
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_mapping_stats"))
    conda: CONDA_ENV
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}" "{output.plot}"
        """
        
###
# Rules to prep and then plot the peak stats:
rule prepping_peak_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"{RESULTS_DIR}/{wildcards.env}/reports/summary_{wildcards.env}_peak_stats_{sample_name}.txt" for sample_name in get_mapped_names_by_env(wildcards.env, analysis_samples) ]
    output:
        temp_stat_file = temp(f"{RESULTS_DIR}/combined/reports/temp_summary_peak_stats_{{analysis_name}}_{{env}}.txt"),
        stat_file = f"{RESULTS_DIR}/combined/reports/summary_peak_stats_{{analysis_name}}_{{env}}.txt"
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "prep_peak_stats"))
    shell:
        """
        printf "Group\tLevels\tMark\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tPeaks_in_merged\tPeaks_in_pseudo_reps\tPeaks_in_idr\tSelected_peaks\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk 'NR>1' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_peak_stats:
    input:
        summary_stats = f"{RESULTS_DIR}/combined/reports/summary_peak_stats_{{analysis_name}}_{{env}}.txt"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/peak_stats_{{analysis_name}}_{{env}}.pdf"
    params:
        analysis_name = lambda wildcards: wildcards.analysis_name,
        env = lambda wildcards: wildcards.env,
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_peak_stats.R")
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_peak_stats"))
    conda: CONDA_ENV
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}" "{output.plot}" "{params.env}"
        """

###
# Rules to prep and then plot the sizes stats for sRNA
rule prepping_srna_sizes_stats:
    input:
        sample_stat_files = lambda wildcards: [ f"{RESULTS_DIR}/sRNA/reports/sizes_stats__{sample_name}.txt"
                                                for sample_name in get_mapped_names_by_env(wildcards.env, samples) ]
    output:
        temp_stat_file = temp(f"{RESULTS_DIR}/combined/reports/temp_summary_sizes_stats_{{analysis_name}}_{{env}}.txt"),
        stat_file = f"{RESULTS_DIR}/combined/reports/summary_sizes_stats_{{analysis_name}}_{{env}}.txt"
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "prep_srna_sizes"))
    shell:
        """
        printf "Sample\tType\tSize\tCount\n" > "{output.stat_file}"
        for f in {input.sample_stat_files}
        do
            awk -F "\t" -v OFS="\t" 'NR>1' $f >> "{output.temp_stat_file}"
        done
        sort {output.temp_stat_file} -u >> "{output.stat_file}"
        """
    
rule plotting_srna_sizes_stats:
    input:
        summary_stats = f"{RESULTS_DIR}/combined/reports/summary_sizes_stats_{{analysis_name}}_{{env}}.txt"
    output:
        plot1 = f"{RESULTS_DIR}/combined/plots/srna_sizes_stats_{{analysis_name}}_{{env}}.pdf",
        plot2 = f"{RESULTS_DIR}/combined/plots/srna_sizes_stats_zoom_{{analysis_name}}_{{env}}.pdf"
    params:
        analysis_name = lambda wildcards: wildcards.analysis_name,
        script=os.path.join(REPO_FOLDER,"workflow","scripts","R_sizes_stats.R"),
        zoommin=config['srna_min_size'],
        zoommax=config['srna_max_size']
    log:
        temp(return_log_combined("{analysis_name}", "{env}", "plot_srna_sizes"))
    conda: CONDA_ENV
    shell:
        """
        Rscript "{params.script}" "{input.summary_stats}" "{params.analysis_name}" "{params.zoommin}" "{params.zoommax}" "{output.plot1}" "{output.plot2}"
        """

###
# Rules to prep and plot ChIP upset plots
rule combine_clusterfiles:
    input:
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/chrom.sizes",
        clusterfiles = lambda wildcards: define_samples_for_upset(wildcards, "files")
    output:
        temp1_file = temp(f"{RESULTS_DIR}/combined/bedfiles/temp1_combined_clusters__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        temp2_file = temp(f"{RESULTS_DIR}/combined/bedfiles/temp2_combined_clusters__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        merged_file = f"{RESULTS_DIR}/combined/bedfiles/combined_clusters__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samples_for_upset(wildcards, "pairs"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_clusters_{env}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Merging peakfiles for {params.env} {params.analysis_name} {params.ref_genome}\n"
        for pair in {params.names}; do
            label=$(echo ${{pair}} | cut -d":" -f1)
            file=$(echo ${{pair}} | cut -d":" -f2)
            awk -v OFS="\t" -v l=${{label}} '{{print $1,$2,$3,l"_"$4}}' ${{file}} >> {output.temp1_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        printf "Chr\tStart\tStop\tClusterID\tSamples\n" > {output.merged_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v e={params.env} -v a={params.analysis_name} '{{print $1,$2,$3,"combined_clusters_"e"_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule combine_peakfiles:
    input:
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/chrom.sizes",
        peakfiles = lambda wildcards: define_samples_for_upset(wildcards, "files")
    output:
        temp1_file = temp(f"{RESULTS_DIR}/combined/bedfiles/temp1_combined_peaks__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        temp2_file = temp(f"{RESULTS_DIR}/combined/bedfiles/temp2_combined_peaks__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        merged_file = f"{RESULTS_DIR}/combined/bedfiles/combined_peaks__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samples_for_upset(wildcards, "pairs"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_peaks_{env}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Merging peakfiles for {params.env} {params.analysis_name} {params.ref_genome}\n"
        for pair in {params.names}; do
            label=$(echo ${{pair}} | cut -d":" -f1)
            file=$(echo ${{pair}} | cut -d":" -f2)
            awk -v OFS="\t" -v l=${{label}} '{{print $1,$2,$3,l}}' ${{file}} >> {output.temp1_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        printf "Chr\tStart\tStop\tPeakID\tSamples\n" > {output.merged_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v e={params.env} -v a={params.analysis_name} '{{print $1,$2,$3,"combined_peaks_"e"_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """
        
rule combine_TSS:
    input:
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/chrom.sizes",
        peakfiles = lambda wildcards: define_samples_for_upset(wildcards, "files")
    output:
        temp1_file = temp(f"{RESULTS_DIR}/combined/bedfiles/temp1_combined_TSS__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        temp2_file = temp(f"{RESULTS_DIR}/combined/bedfiles/temp2_combined_TSS__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        merged_file = f"{RESULTS_DIR}/combined/bedfiles/combined_TSS__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        names = lambda wildcards: define_samples_for_upset(wildcards, "pairs"),
        analysis_name = config['analysis_name']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "combined_TSS_{env}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Merging TSS files for {params.env} {params.analysis_name} {params.ref_genome}\n"
        for pair in {params.names}; do
            label=$(echo ${{pair}} | cut -d":" -f1)
            file=$(echo ${{pair}} | cut -d":" -f2)
            awk -v OFS="\t" -v l=${{label}} '{{print $1,$2,$3,l}}' ${{file}} >> {output.temp1_file}
        done
        sort -k1,1 -k2,2n {output.temp1_file} > {output.temp2_file}
        printf "Chr\tStart\tStop\tTSSID\tSamples\n" > {output.merged_file}
        bedtools merge -i {output.temp2_file} -c 4 -o distinct | bedtools sort -g {input.chrom_sizes} | awk -v OFS="\t" -v a={params.analysis_name} '{{print $1,$2,$3,"combined_TSS_"a"_"NR,$4}}' >> {output.merged_file}
        }} 2>&1 | tee -a "{log}"
        """

rule get_annotations_for_bedfile:
    input:
        bedfile = lambda wildcards: define_combined_target_file(wildcards),
        region_file = lambda wildcards: f"{RESULTS_DIR}/combined/bedfiles/{wildcards.ref_genome}__all_genes.bed",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/chrom.sizes",
        header = lambda wildcards: f"{define_combined_target_file(wildcards)}.header"
    output:
        temp_bedfile = temp(f"{RESULTS_DIR}/combined/bedfiles/temp__{{target_name}}__{{ref_genome}}.bed"),
        annotated_file = f"{RESULTS_DIR}/combined/bedfiles/annotated__{{target_name}}__{{ref_genome}}.bed"
    params:
        target_name = lambda wildcards: wildcards.target_name,
    log:
        temp(return_log_combined("{target_name}", "{ref_genome}", "annotate_bedfile"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Annotating {params.target_name} to the closest genes\n"
        header=$(cat {input.header})
        if [[ "${{header}}" == "no" ]]; then
            awk -v OFS="\t" -v n={params.target_name} '{{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile}
        else
            awk -v OFS="\t" -v n={params.target_name} 'NR>1 {{if ($4=="") $4=n"_"NR; print $1,$2,$3,$4}}' {input.bedfile} > {output.temp_bedfile}
        fi
        printf "Chr\tStart\tStop\tRegionID\tDistance\tGene_strand\tGID\tCategory\n" > {output.annotated_file}
        bedtools closest -a {output.temp_bedfile} -b {input.region_file} -g {input.chrom_sizes} -D ref | awk -v OFS="\t" '{{if ($10=="+") print $1,$2,$3,$4,$11,$10,$8; else print $1,$2,$3,$4,-$11,$10,$8}}' | awk -F"[=;]" -v OFS="\t" '{{print $1,$2}}' | sed 's/gene://' | awk -v OFS="\t" '{{if ($5<-2000) {{d="Distal_downstream"}} else if ($5<0) {{d="Terminator"}} else if ($5==0) {{d="Gene_body"}} else if ($5>2000) {{d="Distal_upstream"}} else {{d="Promoter"}}; print $1,$2,$3,$4,$5,$6,$8,d}}' >> {output.annotated_file}
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_upset_regions:
    input:
        mergedfile = f"{RESULTS_DIR}/combined/bedfiles/{{target_name}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed",
        annotatedfile = f"{RESULTS_DIR}/combined/bedfiles/annotated__{{target_name}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/Upset_{{target_name}}__{{env}}__{{analysis_name}}__{{ref_genome}}.pdf"
    params:
        env = lambda wildcards: wildcards.env,
        types = lambda wildcards: define_samples_for_upset(wildcards, "types"),
        label_types = lambda wildcards: define_samples_for_upset(wildcards, "label_types"),
        relabel = lambda wildcards: define_samples_for_upset(wildcards, "relabel"),
        script = lambda wildcards: define_upset_script(wildcards)
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_upset_{target_name}_{env}"))
    conda: CONDA_ENV
    shell:
        """
        # The merged file names samples by id (see define_samples_for_upset);
        # the script swaps in the plot labels, so a label change reruns only
        # this plot.
        Rscript "{params.script}" "{input.mergedfile}" "{input.annotatedfile}" "{params.env}" "{params.types}" "{params.label_types}" "{params.relabel}" "{output.plot}"
        """

###
# rule to plot heatmaps
rule making_stranded_matrix_on_targetfile:
    input:
        bigwigs = lambda wildcards: define_key_for_plots(wildcards, "bigwigs"),
        target_file = lambda wildcards: define_combined_target_file(wildcards),
        header = lambda wildcards: f"{define_combined_target_file(wildcards)}.header"
    output:
        temp = temp(f"{RESULTS_DIR}/combined/matrix/temp_file_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}_{{strand}}.bed"),
        matrix = temp(f"{RESULTS_DIR}/combined/matrix/matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}__{{strand}}.gz")
    wildcard_constraints:
        strand = "plus|minus|nostrand|unstranded"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        # Columns are named by stable keys, not labels: labels and colours are
        # applied by the plot rules, so changing them never rebuilds a matrix.
        sample_keys = lambda wildcards: define_key_for_plots(wildcards, "keys"),
        matrix = lambda wildcards: wildcards.matrix_param,
        strand = lambda wildcards: wildcards.strand,
        base = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'base'),
        bs = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'bs'),
        base_mc = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'base_mc'),
        bs_mc = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'bs_mc'),
        before = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'before'),
        after = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'after'),
        middle = lambda wildcards: get_heatmap_param(wildcards.matrix_param, 'middle')
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "making_matrix_{matrix_param}_{target_name}_{strand}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        header="$(cat {input.header})"
        if [[ "{params.strand}" == "unstranded" ]]; then
            if [[ "${{header}}" == "no" ]]; then
                cat {input.target_file} > {output.temp}
            else
                awk 'NR>1' {input.target_file} > {output.temp}
            fi
        else
            case "{params.strand}" in
                plus)     awk '$6=="+"' {input.target_file} > {output.temp};;
                minus)    awk '$6=="-"' {input.target_file} > {output.temp};;
                # Unlike the two above, this predicate matches a header line
                # (its column 6 is neither '+' nor '-'), so drop it explicitly.
                nostrand) awk -v h="${{header}}" 'h=="yes" && NR==1 {{next}} $6!="+" && $6!="-"' {input.target_file} > {output.temp};;
            esac
        fi
        # computeMatrix accepts only '+', '-' and '.' in the strand column and
        # panics on anything else ("Strand should either be + or - or . \"?\" is
        # not supported"). GFF3-derived annotations legitimately carry '?' for
        # unknown strand, which is what '.' means here, so normalise rather
        # than drop the regions.
        awk -v OFS="\t" '
            {{ if (NF >= 6 && $6 != "+" && $6 != "-" && $6 != ".") {{ n++; $6 = "." }} print }}
            END {{ if (n > 0) print "WARNING: " n " region(s) had an unaccepted strand character in column 6; replaced with a dot" > "/dev/stderr" }}
        ' {output.temp} > "$TMPDIR/target_strand_fixed.bed"
        mv "$TMPDIR/target_strand_fixed.bed" {output.temp}
        printf "Making {params.strand} strand {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
        if [[ "{params.env}" == "mC" ]]; then
            computeMatrix {params.base_mc} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.sample_keys} -bs {params.bs_mc} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        else
            computeMatrix {params.base} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.sample_keys} -bs {params.bs} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        fi
        }} 2>&1 | tee -a "{log}"
        """
                
rule merging_matrix:
    input:
        get_matrix_inputs
    output:
        f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "merging_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        nfile=$(echo {input} | wc -w)
        if [[ ${{nfile}} -eq 1 ]]; then
            # Relabel rather than copy: computeMatrix names the group after its
            # input BED, which here is an internal temp file, and that name is
            # what ends up printed on the plot as the region label.
            computeMatrixOperations relabel -m {input} \
                --groupLabels "{params.target_name}" -o {output}
        else
            printf "\nMerging stranded matrices aligned by {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
            # rbind merges region groups that share a label and keeps apart any
            # that differ. computeMatrix names each group after its input BED, so
            # the plus and minus halves arrive labelled '..._plus' and '..._minus'
            # and survive the merge as two groups -- while every rule downstream
            # (plotHeatmap, plotProfile, sort_heatmap's relabel) passes one
            # --regionsLabel per group and dies on a count mismatch. Relabelling
            # decides the grouping outright: plus and minus take the target name
            # and fold into one group, and the 'nostrand' pass -- regions whose
            # strand is '.' or '?', for which sense is undefined -- takes its own
            # label and stays a group of its own.
            #
            # deeptools 3 labelled every single-BED group 'genes' regardless of
            # filename, so plus and minus merged by accident there. Don't rely on
            # that default again.
            relabelled=()
            i=0
            for matrix in {input}; do
                i=$((i+1))
                case "${{matrix}}" in
                    *__nostrand.gz) label="{params.target_name}_nostrand";;
                    *)              label="{params.target_name}";;
                esac
                computeMatrixOperations relabel -m "${{matrix}}" \
                    --groupLabels "${{label}}" -o "$TMPDIR/relabelled_${{i}}.gz"
                relabelled+=("$TMPDIR/relabelled_${{i}}.gz")
            done
            computeMatrixOperations rbind -m "${{relabelled[@]}}" -o {output}
        fi
        }} 2>&1 | tee -a "{log}"
        """

def get_all_matrix_inputs(wildcards):
    # The 'all' matrix is the join of the 'most' (non-mC) and mC matrices, each
    # built with its own computeMatrix options; either is absent when the
    # genome has no samples of that kind.
    prefix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{wildcards.matrix_param}"
    suffix = f"{wildcards.analysis_name}__{wildcards.ref_genome}__{wildcards.target_name}.gz"
    envs = set(analysis_samples.loc[analysis_samples['ref_genome'] == wildcards.ref_genome, 'env'])
    inputs = {}
    if envs - {"mC"}:
        inputs["most"] = f"{prefix}__most__{suffix}"
    if "mC" in envs:
        inputs["mc"] = f"{prefix}__mC__{suffix}"
    return inputs

# merging_matrix's output pattern also matches env 'all'; that one is joined.
ruleorder: joining_all_matrix > merging_matrix

rule joining_all_matrix:
    # Not computed from the bigwigs: computeMatrix takes one set of options per
    # run, and the non-mC ones (--missingDataAsZero --skipZeros, bin size 'bs')
    # read an uncovered cytosine as 0 % methylation. Joining the two matrices
    # keeps NaN in the mC columns and each half's own values; the join puts every
    # column on one bin size (the gcd of 'bs' and 'bs_mc', bins repeated).
    input:
        unpack(get_all_matrix_inputs)
    output:
        f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__all__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        args = lambda wildcards, input: " ".join(f"--{k} {v}" for k, v in input.items()),
        script = os.path.join(REPO_FOLDER, "workflow", "scripts", "matrix_join_all.py")
    log:
        temp(return_log_combined("{analysis_name}", "all_{ref_genome}", "joining_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Joining most and mC {params.matrix} matrices for all {params.target_name} on {params.ref_genome}\n"
        python3 "{params.script}" -o {output} {params.args}
        }} 2>&1 | tee -a "{log}"
        """

rule computing_matrix_scales:
    input:
        matrix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz"
    output:
        params_heatmap = temp(f"{RESULTS_DIR}/combined/matrix/params_heatmap_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"),
        params_profile = temp(f"{RESULTS_DIR}/combined/matrix/params_profile_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"),
        params_regions = temp(f"{RESULTS_DIR}/combined/matrix/params_regions_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        scales = config['heatmaps_scales'],
        profile = config['profiles_scale'],
        labeller = os.path.join(REPO_FOLDER, "workflow", "scripts", "matrix_region_labels.py"),
        scales_script = os.path.join(REPO_FOLDER, "workflow", "scripts", "R_matrix_scales.R"),
        cg_scale = config['heat_mcg'],
        chg_scale = config['heat_mchg'],
        chh_scale = config['heat_mchh']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "getting_scales_matrix_{matrix_param}_{target_name}"))
    conda: CONDA_ENV_PLOTS
    shell:
        """
        {{
        # One label per region group, counts taken from the matrix itself --
        # a stranded target can carry a separate 'nostrand' group, and
        # --skipZeros drops regions, so the target file cannot tell us either.
        python3 "{params.labeller}" {input.matrix} > {output.params_regions}

        # One --zMin/--zMax and --yMin/--yMax value per matrix column, in
        # column order (heatmaps_scales: type = one range per mark, the key's
        # last '__' field matched exactly; sample = one per column; default =
        # empty files). Computed in R from the matrix (plot_scales.R): the
        # 10th/90th percentiles dataRange gave, the heat_m* / 0.005 fallbacks,
        # and y from the average profile (profiles_scale) padded by a tenth of
        # its range -- the same averaging code the R profiles draw with.
        printf "Getting scales ({params.scales}) for {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
        Rscript "{params.scales_script}" {input.matrix} "{params.scales}" "{params.profile}" \
            "{params.cg_scale}" "{params.chg_scale}" "{params.chh_scale}" \
            {output.params_heatmap} {output.params_profile}
        }} 2>&1 | tee -a "{log}"
        """

rule plotting_heatmap_on_targetfile:
    input:
        matrix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz",
        params_regions = f"{RESULTS_DIR}/combined/matrix/params_regions_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt",
        params_heatmap = f"{RESULTS_DIR}/combined/matrix/params_heatmap_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt",
        params_profile = f"{RESULTS_DIR}/combined/matrix/params_profile_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/Heatmap__{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf",
        sorted_regions = f"{RESULTS_DIR}/combined/matrix/Heatmap__{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}_sorted_regions.bed"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        env = lambda wildcards: wildcards.env,
        plot_params = lambda wildcards: config['heatmaps_plot_params'][wildcards.env],
        sort = lambda wildcards: define_sort_options(wildcards),
        display = lambda wildcards: define_key_for_plots(wildcards, "display"),
        sample_labeller = os.path.join(REPO_FOLDER, "workflow", "scripts", "matrix_sample_labels.py"),
        plot_heatmap = os.path.join(REPO_FOLDER, "workflow", "scripts", "plotheatmap_lowmem.py")
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "plot_heatmap_{matrix_param}_{target_name}"))
    # Sized from the matrix (heatmap_memory.py); the profile sets only runtime.
    resources:
        mem_mb = lambda wildcards, input, attempt: heatmap_mem_mb(input.matrix, attempt)
    conda: CONDA_ENV
    shell:
        """
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        reg="$(cat {input.params_regions})"
        heat="$(cat {input.params_heatmap})"
        prof="$(cat {input.params_profile})"
        printf '%s\n' "{params.display}" > "$TMPDIR/display.tsv"
        lab="--samplesLabel $(python3 "{params.sample_labeller}" {input.matrix} "$TMPDIR/display.tsv" | cut -f2 | tr '\n' ' ')"
        printf "Plotting heatmap {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        # plotHeatmap through a wrapper that stops matplotlib keeping a
        # full-figure canvas per heatmap in the PDF (plotheatmap_lowmem.py):
        # memory grew with the square of the column count.
        python3 "{params.plot_heatmap}" -m {input.matrix} -out {output.plot} {params.plot_params} {params.sort} ${{reg}} ${{heat}} ${{prof}} ${{lab}} ${{add}} --outFileSortedRegions {output.sorted_regions}
        """

rule sort_heatmap:
    input: 
        matrix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz",
        sorted_regions = f"{RESULTS_DIR}/combined/matrix/Heatmap__{{matrix_param}}__most__{{analysis_name}}__{{ref_genome}}__{{target_name}}_sorted_regions.bed",
        params_regions = f"{RESULTS_DIR}/combined/matrix/params_regions_final_matrix_{{matrix_param}}__most__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"
    output:
        temp_matrix = temp(f"{RESULTS_DIR}/combined/matrix/temp_sorted_final_matrix_{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz"),
        matrix = f"{RESULTS_DIR}/combined/matrix/sorted_final_matrix_{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param
    log:
        temp(return_log_combined("{analysis_name}", "mC_{ref_genome}", "sort_heatmap_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    shell:
        """
        printf "Sorting heatmap {params.matrix} for mC {params.target_name} on {params.ref_genome}\n"
        # Everything after the flag: one label per region group, not just the first.
        labels="$(cat {input.params_regions} | cut -d" " -f 2-)"
        computeMatrixOperations relabel -m {input.matrix} --groupLabels ${{labels}} -o {output.temp_matrix}
        computeMatrixOperations sort -m {output.temp_matrix} -R {input.sorted_regions} -o {output.matrix}
        """

rule plotting_sorted_heatmap_on_targetfile:
    input:
        matrix = f"{RESULTS_DIR}/combined/matrix/sorted_final_matrix_{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz",
        params_heatmap = f"{RESULTS_DIR}/combined/matrix/params_heatmap_final_matrix_{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt",
        params_profile = f"{RESULTS_DIR}/combined/matrix/params_profile_final_matrix_{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__{{matrix_param}}__mC__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        plot_params = lambda wildcards: config['heatmaps_plot_params']['mC'],
        display = lambda wildcards: define_key_for_plots(SimpleNamespace(ref_genome=wildcards.ref_genome, env="mC"), "display"),
        sample_labeller = os.path.join(REPO_FOLDER, "workflow", "scripts", "matrix_sample_labels.py"),
        plot_heatmap = os.path.join(REPO_FOLDER, "workflow", "scripts", "plotheatmap_lowmem.py")
    log:
        temp(return_log_combined("{analysis_name}", "mC_{ref_genome}", "plot_sorted_heatmap_{matrix_param}_{target_name}"))
    resources:
        mem_mb = lambda wildcards, input, attempt: heatmap_mem_mb(input.matrix, attempt)
    conda: CONDA_ENV
    shell:
        """
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        heat="$(cat {input.params_heatmap})"
        prof="$(cat {input.params_profile})"
        printf '%s\n' "{params.display}" > "$TMPDIR/display.tsv"
        lab="--samplesLabel $(python3 "{params.sample_labeller}" {input.matrix} "$TMPDIR/display.tsv" | cut -f2 | tr '\n' ' ')"
        printf "Plotting heatmap {params.matrix} for mC {params.target_name} on {params.ref_genome}\n"
        python3 "{params.plot_heatmap}" -m {input.matrix} -out {output.plot} {params.plot_params} --sortRegions 'keep' ${{heat}} ${{prof}} ${{lab}} ${{add}}
        """

rule plotting_profile_on_targetfile:
    input:
        matrix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz"
    output:
        plot1 = f"{RESULTS_DIR}/combined/plots/Profile__{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf",
        plot2 = f"{RESULTS_DIR}/combined/plots/Profile_pergroup__{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        env = lambda wildcards: wildcards.env,
        scales = config['heatmaps_scales'],
        average = config['profiles_scale'],
        display = lambda wildcards: define_key_for_plots(wildcards, "display"),
        sample_labeller = os.path.join(REPO_FOLDER, "workflow", "scripts", "matrix_sample_labels.py"),
        profile_script = os.path.join(REPO_FOLDER, "workflow", "scripts", "R_profile.R"),
        theme = PLOT_THEME_JSON
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "plot_profile_{matrix_param}_{target_name}"))
    conda: CONDA_ENV_PLOTS
    shell:
        """
        {{
        printf "Plotting profiles {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        # Labels and colours, one row per matrix column in column order, looked
        # up by the column keys (key, label, colour, track, series).
        printf '%s\n' "{params.display}" > "$TMPDIR/display.tsv"
        python3 "{params.sample_labeller}" {input.matrix} "$TMPDIR/display.tsv" > "$TMPDIR/columns.tsv"
        printf '%s\n' {params.theme:q} > "$TMPDIR/plot_theme.json"
        # Profile__ (a panel per column) and Profile_pergroup__ (a panel per
        # region group; for mC a panel per group x context), drawn from the
        # matrix. The y-limits are computed from the curves drawn (plot_scales.R,
        # the y half of what computing_matrix_scales does). heat_mcg & co. only
        # set heatmap colour ranges, so this rule doesn't take them: changing
        # one reruns no profile.
        Rscript "{params.profile_script}" {input.matrix} "$TMPDIR/columns.tsv" "{params.matrix}" "{params.env}" "{params.scales}" "{params.average}" {output.plot1} {output.plot2} "$TMPDIR/plot_theme.json"
        }} 2>&1 | tee -a "{log}"
        """

###
# rules to plot browser shots
rule prep_chromosomes_for_browser:
    input: 
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/chrom.sizes"
    output:
        bedfile = f"{RESULTS_DIR}/combined/bedfiles/full_chromosomes__{{ref_genome}}.bed"
    params:
        chromosome_bs = config['chromosome_bs']
    log:
        temp(return_log_combined("bedfile", "{ref_genome}", "prep_chromosomes"))
    conda: CONDA_ENV
    shell:
        """
        {{
        awk -v OFS="\t" -v c={params.chromosome_bs} 'NR <= 50 {{if ($2/c > 1) b=c; else b=int($2/500); print $1,"1",$2,$1,b}}' {input.chrom_sizes} > {output.bedfile}
        }} 2>&1 | tee -a "{log}" 
        """
        
rule prep_browser_on_region:
    input:
        bigwigs = lambda wildcards: define_key_for_plots(wildcards, "bigwigs"),
        target_file = lambda wildcards: define_combined_target_file(wildcards),
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/chrom.sizes",
        gff = lambda wildcards: f"{GENOMES_DIR}/{wildcards.ref_genome}/{wildcards.ref_genome}.gff",
        all_genes = lambda wildcards: f"{RESULTS_DIR}/combined/bedfiles/{wildcards.ref_genome}__all_genes.bed",
        TE_file = lambda wildcards: get_TE_file(wildcards),
        extra_tracks = [t["file"] for t in BROWSER_EXTRA_TRACKS]
    output:
        filenames = f"{RESULTS_DIR}/combined/matrix/filenames__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        genes = f"{RESULTS_DIR}/combined/matrix/genes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.gff",
        tes = f"{RESULTS_DIR}/combined/matrix/tes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed",
        htstart = f"{RESULTS_DIR}/combined/matrix/highlight_start__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        htwidth = f"{RESULTS_DIR}/combined/matrix/highlight_width__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        name = f"{RESULTS_DIR}/combined/matrix/name__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        extras = f"{RESULTS_DIR}/combined/matrix/extra_tracks__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        tempgenes = temp(f"{RESULTS_DIR}/combined/matrix/genes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt"),
        templocus = temp(f"{RESULTS_DIR}/combined/matrix/locus_{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        temparray = temp(f"{RESULTS_DIR}/combined/matrix/array__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.npz"),
        tempvalues = temp(f"{RESULTS_DIR}/combined/matrix/values__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.tab")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        # Tracks are extracted under their stable keys; labels and colours are
        # joined on at plot time by make_single_loci_browser_plot.
        sample_keys = lambda wildcards: define_key_for_plots(wildcards, "keys"),
        sample_table = lambda wildcards: define_key_for_plots(wildcards, "table"),
        trackfolder = lambda wildcards: f"{RESULTS_DIR}/combined/matrix/tracks_{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}",
        regionID = lambda wildcards: wildcards.regionID,
        # One "label<TAB>color<TAB>path" line per extra track; labels can't
        # hold tabs, quotes or '$' (normalize_extra_tracks), so this is safe
        # to read back in a shell loop.
        extra_tracks = "\n".join(f"{t['label']}\t{t['color']}\t{t['file']}" for t in BROWSER_EXTRA_TRACKS),
        browser_scales = config['browser_scales'],
        mc_scales = config['fixed_mc_scales'],
        cg_scale = config['fixed_mcg'],
        chg_scale = config['fixed_mchg'],
        chh_scale = config['fixed_mchh']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "prep_files_{target_name}_{regionID}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Extracting values for {params.regionID}\n"
        tmpregion="{params.regionID}"
        line_nb=$(echo "${{tmpregion}}" | sed 's/line//')
        chr=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $1}}')
        start=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $2}}')
        end=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $3}}')
        printf "${{chr}}\t${{start}}\t${{end}}\n" > {output.templocus}
        region="${{chr}}:${{start}}:${{end}}"

        regionID=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $4}}')
        printf "${{regionID}}\n" > {output.name}
        
        binsize=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $5}}')
        
        htstart=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $6}}')
        htwidth=$(cat {input.target_file} | awk -v r=${{line_nb}} 'NR==r {{print $7}}')
        if [[ ${{htstart}} != "" ]]; then
            printf "${{htstart}}\n" | awk -F"," '{{for (i=1;i<=NF;i++) print $i}}' > "{output.htstart}"
            printf "${{htwidth}}\n" | awk -F"," '{{for (i=1;i<=NF;i++) print $i}}' > "{output.htwidth}"
        else
            touch {output.htstart}
            touch {output.htwidth}
        fi
        
        ### To get genes in the region
        bedtools intersect -a {input.all_genes} -b {output.templocus} | awk '{{print $4}}' > {output.tempgenes}
        if [[ "{params.target_name}" == "full_chromosomes" ]]; then
            printf "Do not include genes in whole chromosomes\n"
            touch {output.genes}
        elif [[ -s "{output.tempgenes}" ]]; then
            printf "Getting gene track\n"
            bedtools intersect -wb -a {input.gff} -b {output.templocus} | awk -v OFS="\t" '{{if ($7!="+" && $7!="-") $7="*"; print $1,$2,$3,$4,$5,$6,$7,$8,$9}}' > {output.genes}
        else
            printf "No genes in this region\n"
            touch {output.genes}
        fi
        ### To get the bed files of TEs. For now relying on a bed file of TEs (only one, needing to match the species).
        if [[ -n "{input.TE_file}" && -s "{input.TE_file}" && "{params.target_name}" != "full_chromosomes" ]]; then
            printf "Getting TE track\n"
            bedtools intersect -a {input.TE_file} -b {output.templocus} | awk -v OFS="\t" '{{if ($6!="+" && $6!="-") $6="*"; print $0}}' > {output.tes}
        else
            printf "No TE file to be included in browser\n"
            touch {output.tes}
        fi
        
        printf "Testing if all bw have data in this region\n"
        filelist2=()
        rm -rf "{params.trackfolder}" && mkdir "{params.trackfolder}"

        ### Extra user BED tracks, clipped to the locus. Clipped, not just
        # filtered (no -u): Gviz draws the union of every track's range, so one
        # feature reaching past the locus would zoom the whole plot out. The
        # TE and gene tracks are clipped the same way.
        # zcat -f reads plain and gzipped files alike; comment, track/browser
        # and column-header lines are dropped so bedtools sees only records,
        # and only chrom/start/end/name are kept since the track is drawn as
        # plain boxes (a non-numeric score column would break the BED import).
        printf "Label\tColor\tPath\n" > {output.extras}
        n=0
        while IFS=$'\t' read -r label color src; do
            [[ -z "${{src}}" ]] && continue
            n=$(( n + 1 ))
            path="{params.trackfolder}/extra_${{n}}.bed"
            zcat -f "${{src}}" \
                | awk -F"\t" -v OFS="\t" '/^(#|track|browser)/ {{next}} $2 !~ /^[0-9]+$/ {{next}} {{print $1,$2,$3,(NF>=4 ? $4 : ".")}}' \
                | bedtools intersect -a stdin -b {output.templocus} > "${{path}}"
            printf "Extra track %s: %s features in the region\n" "${{label}}" "$(wc -l < "${{path}}")"
            printf "%s\t%s\t%s\n" "${{label}}" "${{color}}" "${{path}}" >> {output.extras}
        done <<< "{params.extra_tracks}"
        while read bw key mark
        do
            path="{params.trackfolder}/${{key}}_empty"
            bigWigToBedGraph -chrom=${{chr}} -start=${{start}} -end=${{end}} ${{bw}} "${{path}}.bg"
            if [[ -s "${{path}}.bg" ]]; then
                printf "${{key}} has data on ${{chr}}:${{start}}-${{end}}\n"
                filelist2+=("${{bw}}")
            else
                printf "${{key}} is empty on ${{chr}}:${{start}}-${{end}}\n"
                awk -v OFS="\t" -v c="${{chr}}" '$1==c {{print $1,"1",$2,"0"}}' {input.chrom_sizes} | bedtools sort -i - > "${{path}}.bg"
                bedGraphToBigWig "${{path}}.bg" {input.chrom_sizes} "${{path}}.bw"
                filelist2+=("${{path}}.bw")
            fi
            rm -f "${{path}}.bg"            
        done < {params.sample_table}
        printf "Summarize bigwigs in binsize of ${{binsize}} bp on {params.regionID}\n"
        multiBigwigSummary bins -b ${{filelist2[@]}} -l {params.sample_keys} -r ${{region}} -p {threads} -bs ${{binsize}} -out {output.temparray} --outRawCounts {output.tempvalues}

        # --outRawCounts writes one column per bigwig, in -b order, after
        # chr/start/end; -b was built from the sample table's rows in order, so
        # row n is column n+3.
        ncols=$(head -1 {output.tempvalues} | awk '{{print NF}}')
        if [[ ${{ncols}} -ne $(( ${{#filelist2[@]}} + 3 )) ]]; then
            printf "Error: {output.tempvalues} has ${{ncols}} columns for ${{#filelist2[@]}} bigwigs\n"
            exit 1
        fi
        col=3
        while read bw key mark
        do
            path="{params.trackfolder}/${{key}}"
            printf "Making bw for ${{key}}\n"
            col=$(( col + 1 ))
            if [[ "${{key}}" == *_minus ]]; then
                awk -v OFS="\t" -v a=${{col}} 'NR>1 {{if ($a == "nan") b=0; else b=-$a; print $1,$2,$3,b}}' {output.tempvalues} | bedtools sort -g {input.chrom_sizes} > "${{path}}.bedGraph"
            else
                awk -v OFS="\t" -v a=${{col}} 'NR>1 {{if ($a == "nan") b=0; else b=$a; print $1,$2,$3,b}}' {output.tempvalues} | bedtools sort -g {input.chrom_sizes} > "${{path}}.bedGraph"
            fi
            bedGraphToBigWig "${{path}}.bedGraph" {input.chrom_sizes} "${{path}}.bw"
        done < {params.sample_table}
        
        # Name is the track key; make_single_loci_browser_plot joins the label
        # and colours onto it.
        printf "Name\tPath\tYmin\tYmax\n" > {output.filenames}
        while read bw key mark
        do
            path="{params.trackfolder}/${{key}}"
            if [[ {params.browser_scales} == "sample" ]]; then
                ymin=$(cat "${{path}}.bedGraph" | awk 'BEGIN {{a=9999}} {{if ($4<a) a=$4;}} END {{if (a<0) b=a*1.2; else b=a*0.8; print b}}')
                ymax=$(cat "${{path}}.bedGraph" | awk 'BEGIN {{a=-9999}} {{if ($4>a) a=$4;}} END {{if (a>0) b=a*1.2; else b=a*0.8; print b}}')
            elif [[ {params.browser_scales} == "type" ]]; then
                ymin=$(cat {params.trackfolder}/*__${{mark}}.bedGraph | awk 'BEGIN {{a=9999}} {{if ($4<a) a=$4;}} END {{if (a<0) b=a*1.2; else b=a*0.8; print b}}')
                ymax=$(cat {params.trackfolder}/*__${{mark}}.bedGraph | awk 'BEGIN {{a=-9999}} {{if ($4>a) a=$4;}} END {{if (a>0) b=a*1.2; else b=a*0.8; print b}}')
            fi
            if [[ {params.mc_scales} == "True" ]]; then
                if [[ ${{mark}} == "mCG" ]]; then
                    ymin=0
                    ymax={params.cg_scale}
                elif [[ ${{mark}} == "mCHG" ]]; then
                    ymin=0
                    ymax={params.chg_scale}
                elif [[ ${{mark}} == "mCHH" ]]; then
                    ymin=0
                    ymax={params.chh_scale}
                fi
            fi
            printf "${{key}}\t${{path}}.bw\t${{ymin}}\t${{ymax}}\n" >> {output.filenames}
        done < {params.sample_table}
        rm -f {params.trackfolder}/*.bedGraph
        
        }} 2>&1 | tee -a "{log}" 
        """

rule make_single_loci_browser_plot:
    input:
        filenames = f"{RESULTS_DIR}/combined/matrix/filenames__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        genes = f"{RESULTS_DIR}/combined/matrix/genes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.gff",
        tes = f"{RESULTS_DIR}/combined/matrix/tes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed",
        name = f"{RESULTS_DIR}/combined/matrix/name__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        htstart = f"{RESULTS_DIR}/combined/matrix/highlight_start__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        htwidth = f"{RESULTS_DIR}/combined/matrix/highlight_width__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        extras = f"{RESULTS_DIR}/combined/matrix/extra_tracks__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt"
    output:
        browser_plot = temp(f"{RESULTS_DIR}/combined/plots/single_browser__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.pdf")
    params:
        regionID = lambda wildcards: wildcards.regionID,
        title = lambda wildcards: f"{RESULTS_DIR}/combined/plots/single_browser__{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}.pdf",
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_browser_plot.R"),
        display = lambda wildcards: define_key_for_plots(wildcards, "browser_display"),
        trackfolder = lambda wildcards: f"{RESULTS_DIR}/combined/matrix/tracks_{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}"
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "single_browser_{target_name}_{regionID}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        name=$(cat {input.name})
        # Empty highlight files and a header-only extras table mean "none";
        # the script checks, so the call is the same either way.
        printf "\nPlotting browser on {params.regionID}\n\n"
        # The extracted tracks are named by key; join on each track's label and
        # colours, keeping the extraction's row order.
        printf '%s\n' "{params.display}" > "$TMPDIR/display.tsv"
        awk -F"\t" -v OFS="\t" 'NR==FNR {{d[$1]=$2 OFS $3 OFS $4 OFS $5 OFS $6; next}}
            FNR==1 {{print "Name", "Path", "Backcolor", "Trackcolor", "Fillcolorplus", "Fillcolorminus", "Ymin", "Ymax"; next}}
            !($1 in d) {{print "No label for browser track " $1 > "/dev/stderr"; exit 1}}
            {{split(d[$1], v, OFS); print v[1], $2, v[2], v[3], v[4], v[5], $3, $4}}' \
            "$TMPDIR/display.tsv" {input.filenames} > "$TMPDIR/filenames.tsv"
        Rscript "{params.script}" "$TMPDIR/filenames.tsv" "{input.genes}" "{input.tes}" "${{name}}" "{params.title}" "{input.htstart}" "{input.htwidth}" "{input.extras}"
        rm -rf {params.trackfolder}
        }} 2>&1 | tee -a "{log}"
        """

rule merge_region_browser_plots:
    input: 
        plots = lambda wildcards: define_individual_browser_plots(wildcards)
    output:
        merged_plots = f"{RESULTS_DIR}/combined/plots/Browser_{{target_name}}__{{env}}__{{analysis_name}}__{{ref_genome}}.pdf"
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "merging_browser_{target_name}"))
    conda: CONDA_ENV
    shell:
        """
        pdfunite {input.plots} {output.merged_plots}
        """

############## To plot PCA

rule summarize_tracks_pca:
    input:
        tracks = lambda wildcards: define_input_for_pca(wildcards, "tracks"),
        # multiBamSummary reads the .bai alongside each BAM, so the indexes must
        # be declared even though the shell command never names them. Declaring
        # them makes the dependency real in two ways Snakemake cannot infer:
        # a missing index gets rebuilt (dispatch_final_bam produces both), and
        # a temp() index is not reclaimed while this job still needs it. Empty
        # for the mC contexts, which summarize bigwigs instead.
        indexes = lambda wildcards: define_input_for_pca(wildcards, "indexes")
    output:
        array = f"{RESULTS_DIR}/combined/matrix/pca_matrix__{{env}}__{{analysis_name}}__{{ref_genome}}.npz"
    params:
        # Samples are named by mapped_name here and labelled by plotPCA, so a
        # label change reruns the plot only.
        sample_keys = lambda wildcards: define_input_for_pca(wildcards, "keys"),
        bs = config['pca_bs'],
        step = config['pca_step']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "summarize_tracks_pca_{env}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        if [[ {wildcards.env} == "mCG" || {wildcards.env} == "mCHG" || {wildcards.env} == "mCHH" ]]; then
            printf "Summarizing bigwigs for {wildcards.analysis_name} {wildcards.ref_genome} for mC samples in {wildcards.env} sequence context\n"
            multiBigwigSummary bins -b {input.tracks} -o {output.array} -l {params.sample_keys} -bs {params.bs} -n {params.step} -p {threads}
        else
            printf "Summarizing bams for {wildcards.analysis_name} {wildcards.ref_genome} for {wildcards.env} samples\n"
            multiBamSummary bins -b {input.tracks} -o {output.array} -l {params.sample_keys} -bs {params.bs} -n {params.step} -p {threads}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule plot_PCA_correlation:
    input:
        array = f"{RESULTS_DIR}/combined/matrix/pca_matrix__{{env}}__{{analysis_name}}__{{ref_genome}}.npz"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/PCA__{{env}}__{{analysis_name}}__{{ref_genome}}.pdf"
    params:
        labels = lambda wildcards: define_input_for_pca(wildcards, "labels"),
        colors = lambda wildcards: define_input_for_pca(wildcards, "colors"),
        bs = config['pca_bs']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_PCA_correlation_{env}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Plotting PCA for {wildcards.analysis_name} {wildcards.ref_genome} for {wildcards.env} samples\n"
        if ! plotPCA -in {input.array} -T "PCA for {wildcards.env} in {params.bs}bp bins" -o {output.plot} --labels {params.labels:q} --colors {params.colors:q} --transpose; then
            printf "WARNING: plotPCA failed (insufficient data for PCA) — creating placeholder\n"
            touch {output.plot}
        fi
        }} 2>&1 | tee -a "{log}"
        """

###
# final rule
rule all_combined:
    input:
        stats = define_final_stats_output(),
        final = lambda wildcards: define_final_combined_output(wildcards.ref_genome)
    output:
        touch = f"{RESULTS_DIR}/combined/chkpts/combined_analysis__{{analysis_name}}__{{ref_genome}}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """        
