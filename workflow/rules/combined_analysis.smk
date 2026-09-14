
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
    elif target_name.startswith("all_TEs"):
        file = f"{GENOMES_DIR}/{ref_genome}/{ref_genome}__TE_file.bed"
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
        if f.read().strip() == "stranded" and stranded_heatmaps:
            return [ f"{prefix}__plus.gz", f"{prefix}__minus.gz" ]
        else:
            return [ f"{prefix}__unstranded.gz" ]

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
    names = []
    files = []
    types = set()
    ref_genome = wildcards.ref_genome
    srna_sizes = config['srna_heatmap_sizes']
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
                    label = f"{row.levels_label}_{row.sample_type}_{rep}"
                    names.append(f"{label}:{file}")
                    files.append(file)
                    types.add(row.sample_type)
            else:
                file = f"{RESULTS_DIR}/{row.env}/peaks/selected_peaks__{row['mapped_name']}.bedPeak"
                label = f"{row.levels_label}_{row.sample_type}"
                names.append(f"{label}:{file}")
                files.append(file)
                types.add(row.sample_type)
        elif globenv == "RAMPAGE":
            if allreps:
                for sid in get_replicate_sample_ids(row['sample_name'], samples):
                    file = f"{RESULTS_DIR}/RNA/TSS/TSS__final__{sid}_peaks.narrowPeak"
                    rep = parse_sample_name(sid)['replicate']
                    label = f"{row.levels_label}_{rep}"
                    names.append(f"{label}:{file}")
                    files.append(file)
                    types.add(f"{row.levels_label}")
            else:
                file = f"{RESULTS_DIR}/RNA/TSS/TSS__merged__{row['sample_name']}_peaks.narrowPeak"
                label = f"{row.levels_label}"
                names.append(f"{label}:{file}")
                files.append(file)
                types.add(f"{row.levels_label}")
        elif row.env == "sRNA":
            for sid in get_replicate_sample_ids(row['sample_name'], samples):
                file = f"{RESULTS_DIR}/sRNA/mapped/{sid}/clusters.bed"
                rep = parse_sample_name(sid)['replicate']
                label = f"{row.levels_label}_{rep}"
                names.append(f"{label}:{file}")
                files.append(file)

    if globenv == "sRNA":
        srna_min = config['srna_min_size']
        srna_max = config['srna_max_size']
        types = [f"{s}nt" for s in range(srna_min, srna_max + 1)]
        types += ["MIRNA", "Others"]
        ordered = ":".join(types)
    else:
        ordered = ":".join(sorted(types))
    
    if string == "pairs":
        return names
    elif string == "files":
        return files
    elif string == "types":
        return ordered

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
    bigwigs = []
    labels = []
    marks = []
    unique_chip = set()
    unique_atac = set()
    unique_rna = set()
    unique_srna = set()
    unique_mc = set()
    grouped_bw = defaultdict(list)
    grouped_labs = defaultdict(list)
    label_to_mark = {}
    label_to_type = {}
    label_to_track = {}
    srna_sizes = config['srna_heatmap_sizes']
    plot_allreps = config['plot_allreps']
    mc_contexts = get_methylation_contexts()
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    strand = getattr(wildcards, "strand", "unstranded")
    
    if globenv == "all":
        filtered_analysis_samples = analysis_samples[ analysis_samples['ref_genome'] == ref_genome ].copy()
    elif globenv == "most":
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] != "mC") & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    else:
        filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == globenv) & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        rep_ids = get_replicate_sample_ids(row['sample_name'], samples)
        if row.env == "ChIP":
            if not plot_allreps:
                bw = f"{RESULTS_DIR}/{row.env}/tracks/FC__merged__{row['mapped_name']}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/FC__final__{rep_ids[0]}.bw"
                label = f"{row.levels_label}_{row.sample_type}"
                grouped_bw[f"chip_{row.sample_type}"].append(bw)
                grouped_labs[f"chip_{row.sample_type}"].append(label)
                unique_chip.add(row.sample_type)
                label_to_mark[label] = row.sample_type
                label_to_type[label] = f"{row.levels_label}"
                label_to_track[label] = row.sample_type
            else:
                for sid in rep_ids:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/FC__final__{sid}.bw"
                    rep = parse_sample_name(sid)['replicate']
                    label = f"{row.levels_label}_{row.sample_type}_{rep}"
                    grouped_bw[f"chip_{row.sample_type}"].append(bw)
                    grouped_labs[f"chip_{row.sample_type}"].append(label)
                    unique_chip.add(row.sample_type)
                    label_to_mark[label] = row.sample_type
                    label_to_type[label] = f"{row.levels_label}"
                    label_to_track[label] = row.sample_type

        elif row.env == "RNA":
            strandedness = config['rna_tracks'][row.data_type]['strandedness']
            if strandedness == "unstranded":
                if not plot_allreps:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__unstranded.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__unstranded.bw"
                    label = f"{row.levels_label}_{row.sample_type}"
                    grouped_bw[f"{row.data_type}"].append(bw)
                    grouped_labs[f"{row.data_type}"].append(f"{label}")
                    unique_rna.add(row.data_type)
                    label_to_mark[label] = row.data_type
                    label_to_type[label] = f"{row.levels_label}"
                    label_to_track[label] = row.data_type
                else:
                    for sid in rep_ids:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__unstranded.bw"
                        rep = parse_sample_name(sid)['replicate']
                        label = f"{row.levels_label}_{row.sample_type}_{rep}"
                        grouped_bw[f"{row.data_type}"].append(bw)
                        grouped_labs[f"{row.data_type}"].append(f"{label}")
                        unique_rna.add(row.data_type)
                        label_to_mark[label] = row.data_type
                        label_to_type[label] = f"{row.levels_label}"
                        label_to_track[label] = row.data_type
            elif strand == "unstranded":
                if not plot_allreps:
                    bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__plus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__plus.bw"
                    bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__minus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__minus.bw"
                    label = f"{row.levels_label}_{row.sample_type}"
                    grouped_bw[f"{row.data_type}"].extend([bw1, bw2])
                    grouped_labs[f"{row.data_type}"].extend([f"{label}_plus", f"{label}_minus"])
                    unique_rna.add(row.data_type)
                    label_to_mark[f"{label}_plus"] = row.data_type
                    label_to_mark[f"{label}_minus"] = row.data_type
                    label_to_type[f"{label}_plus"] = f"{row.levels_label}"
                    label_to_type[f"{label}_minus"] = f"{row.levels_label}"
                    label_to_track[f"{label}_plus"] = f"{row.data_type}_plus"
                    label_to_track[f"{label}_minus"] = f"{row.data_type}_minus"
                else:
                    for sid in rep_ids:
                        bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__plus.bw"
                        bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__minus.bw"
                        rep = parse_sample_name(sid)['replicate']
                        label = f"{row.levels_label}_{row.sample_type}_{rep}"
                        grouped_bw[f"{row.data_type}"].extend([bw1, bw2])
                        grouped_labs[f"{row.data_type}"].extend([f"{label}_plus", f"{label}_minus"])
                        unique_rna.add(row.data_type)
                        label_to_mark[f"{label}_plus"] = row.data_type
                        label_to_mark[f"{label}_minus"] = row.data_type
                        label_to_type[f"{label}_plus"] = f"{row.levels_label}"
                        label_to_type[f"{label}_minus"] = f"{row.levels_label}"
                        label_to_track[f"{label}_plus"] = f"{row.data_type}_plus"
                        label_to_track[f"{label}_minus"] = f"{row.data_type}_minus"
            else:
                if not plot_allreps:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{strand}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{strand}.bw"
                    label = f"{row.levels_label}_{row.sample_type}"
                    grouped_bw[f"{row.data_type}"].append(bw)
                    grouped_labs[f"{row.data_type}"].append(f"{label}")
                    unique_rna.add(row.data_type)
                    label_to_mark[label] = row.data_type
                    label_to_type[label] = f"{row.levels_label}"
                    label_to_track[label] = row.data_type
                else:
                    for sid in rep_ids:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{strand}.bw"
                        rep = parse_sample_name(sid)['replicate']
                        label = f"{row.levels_label}_{row.sample_type}_{rep}"
                        grouped_bw[f"{row.data_type}"].append(bw)
                        grouped_labs[f"{row.data_type}"].append(f"{label}")
                        unique_rna.add(row.data_type)
                        label_to_mark[label] = row.data_type
                        label_to_type[label] = f"{row.levels_label}"
                        label_to_track[label] = row.data_type
                        
        elif row.env == "sRNA":
            for size in srna_sizes:
                if strand == "unstranded":
                    if not plot_allreps:
                        bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{size}nt__plus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{size}nt__plus.bw"
                        bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{size}nt__minus.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{size}nt__minus.bw"
                        label = f"{row.levels_label}_sRNA_{size}nt"
                        grouped_bw[f"sRNA_{size}nt"].extend([bw1, bw2])
                        grouped_labs[f"sRNA_{size}nt"].extend([f"{label}_plus", f"{label}_minus"])
                        unique_srna.add(f"sRNA_{size}nt")
                        label_to_mark[f"{label}_plus"] = f"sRNA_{size}nt"
                        label_to_mark[f"{label}_minus"] = f"sRNA_{size}nt"
                        label_to_type[f"{label}_plus"] = f"{row.levels_label}"
                        label_to_type[f"{label}_minus"] = f"{row.levels_label}"
                        label_to_track[f"{label}_plus"] = f"sRNA_{size}nt_plus"
                        label_to_track[f"{label}_minus"] = f"sRNA_{size}nt_minus"
                    else:
                        for sid in rep_ids:
                            bw1 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{size}nt__plus.bw"
                            bw2 = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{size}nt__minus.bw"
                            rep = parse_sample_name(sid)['replicate']
                            label = f"{row.levels_label}_sRNA_{rep}_{size}nt"
                            grouped_bw[f"sRNA_{size}nt"].extend([bw1, bw2])
                            grouped_labs[f"sRNA_{size}nt"].extend([f"{label}_plus", f"{label}_minus"])
                            unique_srna.add(f"sRNA_{size}nt")
                            label_to_mark[f"{label}_plus"] = f"sRNA_{size}nt"
                            label_to_mark[f"{label}_minus"] = f"sRNA_{size}nt"
                            label_to_type[f"{label}_plus"] = f"{row.levels_label}"
                            label_to_type[f"{label}_minus"] = f"{row.levels_label}"
                            label_to_track[f"{label}_plus"] = f"sRNA_{size}nt_plus"
                            label_to_track[f"{label}_minus"] = f"sRNA_{size}nt_minus"
                else:
                    if not plot_allreps:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{size}nt__{strand}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{size}nt__{strand}.bw"
                        label = f"{row.levels_label}_sRNA_{size}nt"
                        grouped_bw[f"sRNA_{size}nt"].append(bw)
                        grouped_labs[f"sRNA_{size}nt"].append(f"{label}")
                        unique_srna.add(f"sRNA_{size}nt")
                        label_to_mark[label] = f"sRNA_{size}nt"
                        label_to_type[label] = f"{row.levels_label}"
                        label_to_track[label] = f"sRNA_{size}nt"
                    else:
                        for sid in rep_ids:
                            bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{size}nt__{strand}.bw"
                            rep = parse_sample_name(sid)['replicate']
                            label = f"{row.levels_label}_sRNA_{rep}_{size}nt"
                            grouped_bw[f"sRNA_{size}nt"].append(bw)
                            grouped_labs[f"sRNA_{size}nt"].append(f"{label}")
                            unique_srna.add(f"sRNA_{size}nt")
                            label_to_mark[label] = f"sRNA_{size}nt"
                            label_to_type[label] = f"{row.levels_label}"
                            label_to_track[label] = f"sRNA_{size}nt"
                        
        elif row.env == "mC":
            if not plot_allreps:
                for context in mc_contexts:
                    bw = f"{RESULTS_DIR}/{row.env}/tracks/{row['mapped_name']}__{context}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/{row.env}/tracks/{rep_ids[0]}__{context}.bw"
                    label = f"{row.levels_label}_m{context}"
                    grouped_bw[f"m{context}"].append(bw)
                    grouped_labs[f"m{context}"].append(f"{label}")
                    unique_mc.add(f"m{context}")
                    label_to_mark[label] = f"m{context}"
                    label_to_type[label] = f"{row.levels_label}"
                    label_to_track[label] = f"m{context}"
            else:
                for sid in rep_ids:
                    for context in mc_contexts:
                        bw = f"{RESULTS_DIR}/{row.env}/tracks/{sid}__{context}.bw"
                        rep = parse_sample_name(sid)['replicate']
                        label = f"{row.levels_label}_{rep}_m{context}"
                        grouped_bw[f"m{context}"].append(bw)
                        grouped_labs[f"m{context}"].append(f"{label}")
                        unique_mc.add(f"m{context}")
                        label_to_mark[label] = f"m{context}"
                        label_to_type[label] = f"{row.levels_label}"
                        label_to_track[label] = f"m{context}"

        elif row.env == "ATAC":
            if not plot_allreps:
                bw = f"{RESULTS_DIR}/ATAC/tracks/coverage__merged__{row['mapped_name']}.bw" if len(rep_ids) >=2 else f"{RESULTS_DIR}/ATAC/tracks/coverage__final__{rep_ids[0]}.bw"
                label = f"{row.levels_label}_{row.sample_type}"
                grouped_bw["atac"].append(bw)
                grouped_labs["atac"].append(label)
                unique_atac.add("ATAC")
                label_to_mark[label] = "ATAC"
                label_to_type[label] = f"{row.levels_label}"
                label_to_track[label] = "ATAC"
            else:
                for sid in rep_ids:
                    bw = f"{RESULTS_DIR}/ATAC/tracks/coverage__final__{sid}.bw"
                    rep = parse_sample_name(sid)['replicate']
                    label = f"{row.levels_label}_{row.sample_type}_{rep}"
                    grouped_bw["atac"].append(bw)
                    grouped_labs["atac"].append(label)
                    unique_atac.add("ATAC")
                    label_to_mark[label] = "ATAC"
                    label_to_type[label] = f"{row.levels_label}"
                    label_to_track[label] = "ATAC"

    bigwigs = (
        sum([grouped_bw.get(f"chip_{chip}", []) for chip in sorted(unique_chip)], []) +
        sum([grouped_bw.get("atac", [])], []) +
        sum([grouped_bw.get(f"{rna}", []) for rna in sorted(unique_rna)], []) +
        sum([grouped_bw.get(f"{srna}", []) for srna in sorted(unique_srna)], []) +
        sum([grouped_bw.get(f"{mc}", []) for mc in sorted(unique_mc)], [])
    )
    labels = (
        sum([grouped_labs.get(f"chip_{chip}", []) for chip in sorted(unique_chip)], []) +
        sum([grouped_labs.get("atac", [])], []) +
        sum([grouped_labs.get(f"{rna}", []) for rna in sorted(unique_rna)], []) +
        sum([grouped_labs.get(f"{srna}", []) for srna in sorted(unique_srna)], []) +
        sum([grouped_labs.get(f"{mc}", []) for mc in sorted(unique_mc)], [])
    )
    marks = ( sorted(unique_chip) + sorted(unique_atac) + [f"{rna}_{strand}" for rna in sorted(unique_rna) for strand in ["plus", "minus"]] + [f"{srna}_{strand}" for srna in sorted(unique_srna) for strand in ["plus", "minus"]] + sorted(unique_mc) ) if strand == "unstranded" else ( sorted(unique_chip) + sorted(unique_atac) + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )

    marksforbrowser = ( sorted(unique_chip) + sorted(unique_atac) + sorted(unique_rna) + sorted(unique_srna) + sorted(unique_mc) )
    
    types = sorted(filtered_analysis_samples["levels_label"].tolist())
    
    back_palette = assign_colors(types, "tab20")
    track_palette = assign_colors(marksforbrowser, "Set2")
    plus_palette = make_it_lighter(track_palette, 1.1)
    minus_palette = make_it_lighter(track_palette, 1.4)
    for m in unique_rna | unique_srna:
        minus_palette[m] = plus_palette[m]
    
    backcolors = [back_palette[label_to_type[lab]] for lab in labels]
    trackcolors = [track_palette[label_to_mark[lab]] for lab in labels]
    fillcolorsplus = [plus_palette[label_to_mark[lab]] for lab in labels]
    fillcolorsminus = [minus_palette[label_to_mark[lab]] for lab in labels]
    alignedmarks = [label_to_track[lab] for lab in labels]
    
    if string == "bigwigs":
        return bigwigs
    elif string == "labels":
        return labels
    elif string == "marks":
        return marks
    elif string == "table":
        table_name = f"{RESULTS_DIR}/combined/matrix/sample_table__{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}.tab"
        os.makedirs(os.path.dirname(table_name), exist_ok=True)
        tab = pd.DataFrame({
            "bigwigs": bigwigs,
            "labels": labels,
            "backcolors": backcolors,
            "trackcolors": trackcolors,
            "fillcolorsplus": fillcolorsplus,
            "fillcolorsminus": fillcolorsminus,
            "marks": alignedmarks,
        })
        tab.to_csv(table_name, sep="\t", index=False, header=False)
        return table_name

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

def get_TE_file(wildcards):
    ref_genome = wildcards.ref_genome
    te_file = config['genomes'][ref_genome].get('te_file')
    if config['browser_TE_file'] and te_file:
        return f"{GENOMES_DIR}/{ref_genome}/{ref_genome}__TE_file.bed"
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
    labels = []
    unique_group = set()
    label_to_group = {}
    ref_genome = wildcards.ref_genome
    globenv = wildcards.env
    
    if globenv in ["mCG", "mCHG", "mCHH"]:
        filtered_samples = samples[ (samples['env'] == "mC") & (samples['ref_genome'] == ref_genome) ].copy()
        context = globenv[1:]
        for _, row in filtered_samples.iterrows():
            bw = f"{RESULTS_DIR}/mC/tracks/{row['mapped_name']}__{context}.bw"
            label = f"{row.levels_label}_{row.replicate}"
            group = f"{row.levels_label}"
            tracks.append(bw)
            labels.append(label)
            unique_group.add(group)
            label_to_group[label] = group
    elif globenv in ["ChIP", "ATAC"]:
        filtered_samples = samples[ (samples['env'] == globenv) & (samples['ref_genome'] == ref_genome) ].copy()
        for _, row in filtered_samples.iterrows():
            bam = f"{RESULTS_DIR}/{globenv}/mapped/final__{row['mapped_name']}.bam"
            label = f"{row.sample_type}_{row.levels_label}_{row.data_type}_{row.replicate}"
            group = f"{row.sample_type}_{row.levels_label}"
            tracks.append(bam)
            indexes.append(bam + ".bai")
            labels.append(label)
            unique_group.add(group)
            label_to_group[label] = group
    elif globenv == "all_chip":
        filtered_samples = samples[ (samples['env'].isin(["ChIP","ATAC"])) & (samples['ref_genome'] == ref_genome) ].copy()
        for _, row in filtered_samples.iterrows():
            bam = f"{RESULTS_DIR}/{row.env}/mapped/final__{row['mapped_name']}.bam"
            label = f"{row.sample_type}_{row.levels_label}_{row.data_type}_{row.replicate}"
            group = f"{row.sample_type}_{row.levels_label}"
            tracks.append(bam)
            indexes.append(bam + ".bai")
            labels.append(label)
            unique_group.add(group)
            label_to_group[label] = group
    
    palette = assign_colors(unique_group, "tab20")
    colors = [palette[label_to_group[lab]] for lab in labels]
    
    if string == "tracks":
        return tracks
    elif string == "indexes":
        return indexes
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
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__regions__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__tss__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap_sorted__tes__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        else:
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
            te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__mC__{analysis_name}__{ref_genome}__all_TEs.pdf")
        
    else:
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__regions__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tss__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Heatmap__tes__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
    
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        plot_files.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_genes.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__regions__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tss__most__{analysis_name}__{ref_genome}__all_TEs.pdf")
        te_plots.append(f"{RESULTS_DIR}/combined/plots/Profile__tes__most__{analysis_name}__{ref_genome}__all_TEs.pdf")

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

        with open(output.file, "w") as out:
            if has_strand and strand_values.issubset({"+", "-"}):
                out.write("stranded\n")
            else:
                out.write("unstranded\n")

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
        script = lambda wildcards: define_upset_script(wildcards)
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_upset_{target_name}_{env}"))
    conda: CONDA_ENV
    shell:
        """
        Rscript "{params.script}" "{input.mergedfile}" "{input.annotatedfile}" "{params.env}" "{params.types}" "{output.plot}"
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
        strand = "plus|minus|unstranded"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        labels = lambda wildcards: define_key_for_plots(wildcards, "labels"),
        marks = lambda wildcards: define_key_for_plots(wildcards, "marks"),
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
                plus)   sign="+";;
                minus)  sign="-";;
            esac
            awk -v s=${{sign}} '$6==s' {input.target_file} > {output.temp}
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
        echo "{params.labels}" | xargs -n1 > "{config[output_dir]}/combined/matrix/labels_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt"
        echo "{params.marks}" | xargs -n1 > "{config[output_dir]}/combined/matrix/marks_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt"
        printf "Making {params.strand} strand {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
        if [[ "{params.env}" == "mC" ]]; then
            computeMatrix {params.base_mc} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.labels} -bs {params.bs_mc} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
        else
            computeMatrix {params.base} -R {output.temp} -S {input.bigwigs} --samplesLabel {params.labels} -bs {params.bs} -b {params.before} -a {params.after} {params.middle} -p {threads} -o {output.matrix}
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
        if [[ ${{nfile}} -eq 2 ]]; then
            printf "\nMerging stranded matrices aligned by {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations rbind -m {input} -o {output}
        else
            cp {input} {output}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule computing_matrix_scales:
    input:
        matrix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz",
        target_file = lambda wildcards: define_combined_target_file(wildcards),
        header = lambda wildcards: f"{define_combined_target_file(wildcards)}.header"
    output:
        params_heatmap = temp(f"{RESULTS_DIR}/combined/matrix/params_heatmap_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"),
        params_profile = temp(f"{RESULTS_DIR}/combined/matrix/params_profile_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"),
        params_regions = temp(f"{RESULTS_DIR}/combined/matrix/params_regions_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"),
        temp_values = temp(f"{RESULTS_DIR}/combined/matrix/temp_values_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"),
        temp_profile = temp(f"{RESULTS_DIR}/combined/matrix/temp_profile_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf"),
        temp_profile_values = temp(f"{RESULTS_DIR}/combined/matrix/temp_profile_values_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        env = lambda wildcards: wildcards.env,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        scales = config['heatmaps_scales'],
        profile = config['profiles_scale'],
        cg_scale = config['heat_mcg'],
        chg_scale = config['heat_mchg'],
        chh_scale = config['heat_mchh']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "getting_scales_matrix_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    shell:
        """
        {{        
        header="$(cat {input.header})"
        count=$(wc -l {input.target_file} | cut -d' ' -f 1)
        if [[ "${{header}}" == "yes" ]]; then
            count=$((count-1))
        fi
        awk -v ORS="" -v r=${{count}} -v n={params.target_name} 'BEGIN {{print "--regionsLabel "n"("r")"}}' > {output.params_regions}

        if [[ "{params.scales}" == "default" ]]; then
            touch {output.params_heatmap}
            touch {output.params_profile}
            touch {output.temp_values}
            touch {output.temp_profile}
            touch {output.temp_profile_values}
            
        elif [[ "{params.scales}" == "type" ]]; then
            printf "Getting scales per type for {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations dataRange -m {input.matrix} > {output.temp_values}
            plotProfile -m {input.matrix} -out {output.temp_profile} --averageType {params.profile} --outFileNameData {output.temp_profile_values}
            
            mins=()
            maxs=()
            ymins=()
            ymaxs=()
            while read mark
            do
                zmini=$(grep "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=999999}} {{a=$5; if (a != "nan" && a<m) m=a;}} END {{if (m != 999999) print m; else print 0}}')
                zmaxi=$(grep "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=-999999}} {{a=$6; if (a != "nan" && a>m) m=a;}} END {{if (m != -999999) print m; else print 0}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" && ${{mark}} == "mCG" ]]; then
                    zmini="0"
                    zmaxi="{params.cg_scale}"
                elif [[ "${{test}}" == "yes" && ${{mark}} == "mCHG" ]]; then
                    zmini="0"
                    zmaxi="{params.chg_scale}"
                elif [[ "${{test}}" == "yes" && ${{mark}} == "mCHH" ]]; then
                    zmini="0"
                    zmaxi="{params.chh_scale}"
                elif [[ "${{test}}" == "yes" ]]; then
                    zmini="0"
                    zmaxi="0.005"
                fi
                ymini=$(grep "${{mark}}" {output.temp_profile_values} | awk '{{m=$3; for (i=3;i<=NF;i++) if ($i<m) m=$i; print m}}' | awk 'BEGIN {{m=99999}} {{if ($1<m) m=$1}} END {{if (m<0) a=m*1.2; else a=m*0.8; print a}}')
                ymaxi=$(grep "${{mark}}" {output.temp_profile_values} | awk '{{m=$3; for (i=3;i<=NF;i++) if ($i>m) m=$i; print m}}' | awk 'BEGIN {{m=-99999}} {{if ($1>m) m=$1}} END {{if (m<0) a=m*0.8; else a=m*1.2; print a}}')
                test=$(awk -v a=${{ymini}} -v b=${{ymaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ ${{test}} == "yes" ]]; then
                    ymini=("0")
                    ymaxi=("0.01")
                fi
                num=$(grep "${{mark}}" {output.temp_values} | wc -l)
                for i in $(seq 1 ${{num}})
                do
                    zmins+=("$zmini")
                    zmaxs+=("$zmaxi")
                    ymins+=("$ymini")
                    ymaxs+=("$ymaxi")
                done
            done < {config[output_dir]}/combined/matrix/marks_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt
            
            awk -v ORS="" -v a="${{zmins[*]}}" -v b="${{zmaxs[*]}}" 'BEGIN {{print "--zMin "a" --zMax "b}}' > {output.params_heatmap}
            awk -v ORS="" -v c="${{ymins[*]}}" -v d="${{ymaxs[*]}}" 'BEGIN {{print "--yMin "c" --yMax "d}}' > {output.params_profile}
        
        elif [[ "{params.scales}" == "sample" ]]; then
            printf "Getting scales per sample for {params.matrix} matrix for {params.env} {params.target_name} on {params.ref_genome}\n"
            computeMatrixOperations dataRange -m {input.matrix} > {output.temp_values}
            plotProfile -m {input.matrix} -out {output.temp_profile} --averageType {params.profile} --outFileNameData {output.temp_profile_values}
            
            zmins=()
            zmaxs=()
            ymins=()
            ymaxs=()
            while read sample
            do
                zmini=$(grep "${{sample}}" {output.temp_values} | awk '{{if ($5 != "nan") print $5; else print 0}}')
                zmaxi=$(grep "${{sample}}" {output.temp_values} | awk '{{if ($6 != "nan") print $6; else print 0}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" && "${{sample}}" =~ mCG ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.cg_scale}")
                elif [[ "${{test}}" == "yes" && "${{sample}}" =~ mCHG ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.chg_scale}")
                elif [[ "${{test}}" == "yes" && "${{sample}}" =~ mCHH ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.chh_scale}")
                elif [[ "${{test}}" == "yes" ]]; then
                    zmins+=("0")
                    zmaxs+=("0.005")
                else
                    zmins+=("$zmini")
                    zmaxs+=("$zmaxi")
                fi
                
                ymini=$(grep "${{sample}}" {output.temp_profile_values} | awk '{{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}}' | awk 'BEGIN {{m=99999}} {{if ($1<m) m=$1}} END {{if (m<0) a=m*1.2; else a=m*0.8; print a}}')
                ymaxi=$(grep "${{sample}}" {output.temp_profile_values} | awk '{{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}}' | awk 'BEGIN {{m=-99999}} {{if ($1>m) m=$1}} END {{if (m<0) a=m*0.8; else a=m*1.2; print a}}')
                test=$(awk -v a=${{ymini}} -v b=${{ymaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" ]]; then
                    ymins+=("0")
                    ymaxs+=("0.01")
                else
                    ymins+=("$ymini")
                    ymaxs+=("$ymaxi")
                fi
            done < {config[output_dir]}/combined/matrix/labels_{params.matrix}__{params.env}__{params.analysis_name}__{params.ref_genome}__{params.target_name}.txt
            
            awk -v ORS="" -v a="${{zmins[*]}}" -v b="${{zmaxs[*]}}" 'BEGIN {{print "--zMin "a" --zMax "b}}' > {output.params_heatmap}
            awk -v ORS="" -v c="${{ymins[*]}}" -v d="${{ymaxs[*]}}" 'BEGIN {{print "--yMin "c" --yMax "d}}' > {output.params_profile}
        else
            printf "{params.scales} unknown. Returning default\n"
            touch {output.params_heatmap}
            touch {output.params_profile}
            touch {output.temp_values}
            touch {output.temp_profile}
            touch {output.temp_profile_values}
        fi
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
        sort = lambda wildcards: define_sort_options(wildcards)
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "plot_heatmap_{matrix_param}_{target_name}"))
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
        printf "Plotting heatmap {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        plotHeatmap -m {input.matrix} -out {output.plot} {params.plot_params} {params.sort} ${{reg}} ${{heat}} ${{prof}} ${{add}} --outFileSortedRegions {output.sorted_regions}
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
        label="$(cat {input.params_regions} | cut -d" " -f 2)"
        computeMatrixOperations relabel -m {input.matrix} --groupLabels ${{label}} -o {output.temp_matrix}
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
        plot_params = lambda wildcards: config['heatmaps_plot_params']['mC']
    log:
        temp(return_log_combined("{analysis_name}", "mC_{ref_genome}", "plot_sorted_heatmap_{matrix_param}_{target_name}"))
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
        printf "Plotting heatmap {params.matrix} for mC {params.target_name} on {params.ref_genome}\n"
        plotHeatmap -m {input.matrix} -out {output.plot} {params.plot_params} --sortRegions 'keep' ${{heat}} ${{prof}} ${{add}}
        """

rule plotting_profile_on_targetfile:
    input:
        matrix = f"{RESULTS_DIR}/combined/matrix/final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.gz",
        params_regions = f"{RESULTS_DIR}/combined/matrix/params_regions_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt",
        params_profile = f"{RESULTS_DIR}/combined/matrix/params_profile_final_matrix_{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.txt"
    output:
        plot1 = f"{RESULTS_DIR}/combined/plots/Profile__{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf",
        plot2 = f"{RESULTS_DIR}/combined/plots/Profile_pergroup__{{matrix_param}}__{{env}}__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf"
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        matrix = lambda wildcards: wildcards.matrix_param,
        env = lambda wildcards: wildcards.env,
        plot_params = config['profiles_plot_params']
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "plot_profile_{matrix_param}_{target_name}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        if [[ "{params.matrix}" == "tes" ]]; then
            add="--refPointLabel end"
        elif [[ "{params.matrix}" == "tss" ]]; then
            add="--refPointLabel start"
        else
            add="--startLabel start --endLabel end"
        fi
        printf "Plotting profile {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        reg="$(cat {input.params_regions})"
        prof="$(cat {input.params_profile})"
        plotProfile -m {input.matrix} -out {output.plot1} {params.plot_params} ${{reg}} ${{prof}} ${{add}}
        
        printf "Plotting per group profile {params.matrix} for {params.env} {params.target_name} on {params.ref_genome}\n"
        ymin=$(cat {input.params_profile} | awk 'BEGIN {{y=99999}} {{for (i=1; i<=NF; i++) {{if ($i == "--yMin") {{for (j=i+1; j<=NF && $j !~ /^--/; j++) {{if ($j<y) y=$j}} break}} }} }} END {{print y}}' )
        ymax=$(cat {input.params_profile} | awk 'BEGIN {{y=-99999}} {{for (i=1; i<=NF; i++) {{if ($i == "--yMax") {{for (j=i+1; j<=NF && $j !~ /^--/; j++) {{if ($j>y) y=$j}} break}} }} }} END {{print y}}' )
        plotProfile -m {input.matrix} -out {output.plot2} {params.plot_params} ${{reg}} --yMin ${{ymin}} --yMax ${{ymax}} ${{add}} --perGroup
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
        TE_file = lambda wildcards: get_TE_file(wildcards)
    output:
        filenames = f"{RESULTS_DIR}/combined/matrix/filenames__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        genes = f"{RESULTS_DIR}/combined/matrix/genes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.gff",
        tes = f"{RESULTS_DIR}/combined/matrix/tes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed",
        htstart = f"{RESULTS_DIR}/combined/matrix/highlight_start__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        htwidth = f"{RESULTS_DIR}/combined/matrix/highlight_width__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        name = f"{RESULTS_DIR}/combined/matrix/name__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt",
        tempgenes = temp(f"{RESULTS_DIR}/combined/matrix/genes_in_locus__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt"),
        templocus = temp(f"{RESULTS_DIR}/combined/matrix/locus_{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.bed"),
        temparray = temp(f"{RESULTS_DIR}/combined/matrix/array__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.npz"),
        tempvalues = temp(f"{RESULTS_DIR}/combined/matrix/values__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.tab")
    params:
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name,
        labels = lambda wildcards: define_key_for_plots(wildcards, "labels"),
        sample_table = lambda wildcards: define_key_for_plots(wildcards, "table"),
        trackfolder = lambda wildcards: f"{RESULTS_DIR}/combined/matrix/tracks_{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}",
        regionID = lambda wildcards: wildcards.regionID,
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
        while read bw lab back track plus minus mark
        do
            path="{params.trackfolder}/${{lab}}_empty"
            bigWigToBedGraph -chrom=${{chr}} -start=${{start}} -end=${{end}} ${{bw}} "${{path}}.bg"
            if [[ -s "${{path}}.bg" ]]; then
                printf "${{lab}} has data on ${{chr}}:${{start}}-${{end}}\n"
                filelist2+=("${{bw}}")
            else
                printf "${{lab}} is empty on ${{chr}}:${{start}}-${{end}}\n"
                awk -v OFS="\t" -v c="${{chr}}" '$1==c {{print $1,"1",$2,"0"}}' {input.chrom_sizes} | bedtools sort -i - > "${{path}}.bg"
                bedGraphToBigWig "${{path}}.bg" {input.chrom_sizes} "${{path}}.bw"
                filelist2+=("${{path}}.bw")
            fi
            rm -f "${{path}}.bg"            
        done < {params.sample_table}
        printf "Summarize bigwigs in binsize of ${{binsize}} bp on {params.regionID}\n"
        multiBigwigSummary bins -b ${{filelist2[@]}} -l {params.labels} -r ${{region}} -p {threads} -bs ${{binsize}} -out {output.temparray} --outRawCounts {output.tempvalues}
        
        while read bw lab back track plus minus mark
        do
            path="{params.trackfolder}/${{lab}}_${{mark}}"
            printf "Making bw for ${{lab}}\n"
            col=($(awk -v ORS=" " -v t=${{lab}} 'NR==1 {{for (i=1;i<=NF;i++) if ($i~t) print i}}' {output.tempvalues}))
            if [[ "${{lab}}" == *_minus ]]; then
                awk -v OFS="\t" -v a=${{col}} 'NR>1 {{if ($a == "nan") b=0; else b=-$a; print $1,$2,$3,b}}' {output.tempvalues} | bedtools sort -g {input.chrom_sizes} > "${{path}}.bedGraph"
            else
                awk -v OFS="\t" -v a=${{col}} 'NR>1 {{if ($a == "nan") b=0; else b=$a; print $1,$2,$3,b}}' {output.tempvalues} | bedtools sort -g {input.chrom_sizes} > "${{path}}.bedGraph"
            fi
            bedGraphToBigWig "${{path}}.bedGraph" {input.chrom_sizes} "${{path}}.bw"
        done < {params.sample_table}
        
        printf "Name\tPath\tBackcolor\tTrackcolor\tFillcolorplus\tFillcolorminus\tYmin\tYmax\n" > {output.filenames}
        while read bw lab back track plus minus mark
        do
            path="{params.trackfolder}/${{lab}}_${{mark}}"
            if [[ {params.browser_scales} == "sample" ]]; then
                ymin=$(cat "${{path}}.bedGraph" | awk 'BEGIN {{a=9999}} {{if ($4<a) a=$4;}} END {{if (a<0) b=a*1.2; else b=a*0.8; print b}}')
                ymax=$(cat "${{path}}.bedGraph" | awk 'BEGIN {{a=-9999}} {{if ($4>a) a=$4;}} END {{if (a>0) b=a*1.2; else b=a*0.8; print b}}')
            elif [[ {params.browser_scales} == "type" ]]; then
                ymin=$(cat {params.trackfolder}/*_${{mark}}.bedGraph | awk 'BEGIN {{a=9999}} {{if ($4<a) a=$4;}} END {{if (a<0) b=a*1.2; else b=a*0.8; print b}}')
                ymax=$(cat {params.trackfolder}/*_${{mark}}.bedGraph | awk 'BEGIN {{a=-9999}} {{if ($4>a) a=$4;}} END {{if (a>0) b=a*1.2; else b=a*0.8; print b}}')
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
            printf "${{lab}}\t${{path}}.bw\t${{back}}\t${{track}}\t${{plus}}\t${{minus}}\t${{ymin}}\t${{ymax}}\n" >> {output.filenames}
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
        htwidth = f"{RESULTS_DIR}/combined/matrix/highlight_width__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.txt"
    output:
        browser_plot = temp(f"{RESULTS_DIR}/combined/plots/single_browser__{{target_name}}__{{regionID}}__{{env}}__{{analysis_name}}__{{ref_genome}}.pdf")
    params:
        regionID = lambda wildcards: wildcards.regionID,
        title = lambda wildcards: f"{RESULTS_DIR}/combined/plots/single_browser__{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}.pdf",
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_browser_plot.R"),
        trackfolder = lambda wildcards: f"{RESULTS_DIR}/combined/matrix/tracks_{wildcards.target_name}__{wildcards.regionID}__{wildcards.env}__{wildcards.analysis_name}__{wildcards.ref_genome}"
    log:
        temp(return_log_combined("{analysis_name}", "{env}_{ref_genome}", "single_browser_{target_name}_{regionID}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        name=$(cat {input.name})
        if [[ -s {input.htstart} ]]; then
            printf "\nPlotting browser on {params.regionID} with higlights\n\n"
            Rscript "{params.script}" "{input.filenames}" "{input.genes}" "{input.tes}" "${{name}}" "{params.title}" "{input.htstart}" "{input.htwidth}"
        else
            printf "\nPlotting browser on {params.regionID} without higlights\n\n"
            Rscript "{params.script}" "{input.filenames}" "{input.genes}" "{input.tes}" "${{name}}" "{params.title}"
        fi
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
        labels = lambda wildcards: define_input_for_pca(wildcards, "labels"),
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
            multiBigwigSummary bins -b {input.tracks} -o {output.array} -l {params.labels} -bs {params.bs} -n {params.step} -p {threads}
        else
            printf "Summarizing bams for {wildcards.analysis_name} {wildcards.ref_genome} for {wildcards.env} samples\n"
            multiBamSummary bins -b {input.tracks} -o {output.array} -l {params.labels} -bs {params.bs} -n {params.step} -p {threads}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule plot_PCA_correlation:
    input:
        array = f"{RESULTS_DIR}/combined/matrix/pca_matrix__{{env}}__{{analysis_name}}__{{ref_genome}}.npz"
    output:
        plot = f"{RESULTS_DIR}/combined/plots/PCA__{{env}}__{{analysis_name}}__{{ref_genome}}.pdf"
    params:
        colors = lambda wildcards: define_input_for_pca(wildcards, "colors"),
        bs = config['pca_bs']
    log:
        temp(return_log_combined("{analysis_name}", "{ref_genome}", "plot_PCA_correlation_{env}"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Plotting PCA for {wildcards.analysis_name} {wildcards.ref_genome} for {wildcards.env} samples\n"
        if ! plotPCA -in {input.array} -T "PCA for {wildcards.env} in {params.bs}bp bins" -o {output.plot} --colors {params.colors:q} --transpose; then
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
