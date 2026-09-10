CONDA_ENV_CHIP=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_chip.yaml")
CONDA_ENV_IDR=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_idr.yaml")

def return_log_chip(env, sample_name, step, paired):
    return os.path.join(REPO_FOLDER, RESULTS_DIR,env,"logs",f"tmp__{sample_name}__{step}__{paired}.log")

def get_mapping_strategy(env):
    """Select mapping strategy based on environment."""
    if env == "ATAC":
        return config.get('atac_mapping_strategy', 'default')
    return config.get('chip_mapping_strategy', 'default')

def get_chip_aligner(env):
    """Get configured aligner for ChIP or ATAC."""
    if env == "ATAC":
        return config.get('atac_aligner', 'bowtie2')
    return config.get('chip_aligner', 'bowtie2')

def get_effective_aligner(env):
    """Get effective aligner, auto-falling back to bowtie2 for repeat/repeatall strategies."""
    aligner = get_chip_aligner(env)
    strategy = get_mapping_strategy(env)
    if aligner == "chromap" and strategy in ("repeat", "repeatall"):
        return "bowtie2"
    return aligner

def get_mapping_index(wildcards):
    """Return the appropriate index path based on effective aligner."""
    ref_genome = parse_sample_name(wildcards.sample_name)['ref_genome']
    aligner = get_effective_aligner(wildcards.env)
    if aligner == "chromap":
        return f"{GENOMES_DIR}/{ref_genome}/chromap_index/{ref_genome}.index"
    return f"{GENOMES_DIR}/{ref_genome}/bt2_index"

def get_mapping_fasta(wildcards):
    """Return FASTA path (needed by chromap; always available from environment_setup)."""
    ref_genome = parse_sample_name(wildcards.sample_name)['ref_genome']
    return f"{GENOMES_DIR}/{ref_genome}/{ref_genome}.fa"

def get_mapq_filter(env):
    """Extract MAPQ threshold from mapping strategy filter string."""
    import re
    strategy = get_mapping_strategy(env)
    filter_str = config['bt2_mapping_strategy'][strategy]['filter']
    match = re.search(r'-q\s+(\d+)', filter_str)
    return int(match.group(1)) if match else 0

def get_peaktype_for_env(sample_name_or_assay, env):
    """Return peaktype using the appropriate config for the environment.

    Accepts an internal assay token (e.g. "ChIP_broad") or a sample/analysis
    name, which _resolve_assay looks up. ATAC takes its type from
    atac_callpeaks.peaktype so the setting stays user-overridable.
    """
    if env == "ATAC":
        return config["atac_callpeaks"]["peaktype"]
    assay = _resolve_assay(sample_name_or_assay)
    return _peaktype_for_assay(assay)

def _resolve_peakcaller(sample_name_or_assay):
    """Return the peak-calling tool for a sample: ``macs2``, ``seacr``, or ``epic2``.

    For CUT&RUN/CUT&Tag samples, the family-level default in
    ``cut_callpeaks`` is consulted (``broad_caller`` for *_broad assays,
    ``narrow_caller`` for *_narrow). All other peak-calling assays
    (ChIP_*, ATAC) currently default to MACS2.
    """
    assay = _resolve_assay(sample_name_or_assay)
    if assay.startswith(("CUT_RUN", "CUT_TAG")):
        cut_cfg = config.get("cut_callpeaks", {})
        if assay.endswith("_broad"):
            return cut_cfg.get("broad_caller", "epic2")
        return cut_cfg.get("narrow_caller", "seacr")
    return "macs2"


def _resolve_assay(sample_name_or_assay):
    """Resolve an internal assay token or a sample/analysis name to its Assay.

    Per-replicate callers pass a genome-qualified ``mapped_name``, so the
    per-replicate lookup goes through match_sample_rows, which accepts the bare
    and qualified forms both. Analysis names are matched directly, since they
    already carry the genome.

    Raises rather than returning the input unchanged: an unresolved name used to
    fall through to a "broad"/macs2 default, which silently mis-called peaks
    instead of failing.
    """
    if sample_name_or_assay in ASSAY_TO_PEAKTYPE:
        return sample_name_or_assay
    match = match_sample_rows(samples, sample_name_or_assay)
    if not match.empty:
        return match.iloc[0]["Assay"]
    match2 = analysis_samples.loc[analysis_samples["sample_name"] == sample_name_or_assay]
    if not match2.empty:
        return match2.iloc[0]["Assay"]
    raise ValueError(
        f"Cannot resolve '{sample_name_or_assay}' to an assay: it is not an "
        f"internal assay token, a Sample_ID/mapped_name in the sample sheet, "
        f"or a known analysis name."
    )


def _peaktype_for_assay(assay):
    """broad/narrow for an internal combined assay token (or ATAC).

    Raises on an unknown token. The previous ``.get(assay, "broad")`` default
    meant a resolution failure produced broad-called peaks written to a
    narrow-named path, with nothing in the log to say so.
    """
    try:
        return ASSAY_TO_PEAKTYPE[assay]
    except KeyError:
        raise ValueError(
            f"No peak type for assay '{assay}'. Expected one of "
            f"{sorted(ASSAY_TO_PEAKTYPE)}; pulldown rows need Assay + Peak_type "
            f"(e.g. Assay=ChIP, Peak_type=broad)."
        ) from None


def _checked_peaktype(sample_name, env, requested):
    """Return the peak type from the output filename, checked against the sheet.

    The filename carries the peak type, so returning ``requested`` is what makes
    the name and the contents the same value — they cannot drift. The check
    catches the other direction: a path requested by hand that contradicts the
    sample's Peak_type fails instead of being quietly honoured.
    """
    expected = get_peaktype_for_env(sample_name, env)
    if requested != expected:
        raise ValueError(
            f"Peak type mismatch for '{sample_name}': the requested file is "
            f"'{requested}Peak' but its Peak_type resolves to '{expected}'. "
            f"Ask for the '{expected}Peak' path, or change Peak_type in the "
            f"sample sheet."
        )
    return requested

def assign_mapping_paired(wildcards, rulename, outputfile):
    sname = wildcards.sample_name
    env = get_sample_info_from_name(sname, samples, 'env')
    paired = get_sample_info_from_name(sname, samples, 'paired')
    aligned_bams = config['aligned_bams']
    if paired == "PE":
        rule_obj = getattr(rules, f"{rulename}_pe")
    elif paired == "SE":
        rule_obj = getattr(rules, f"{rulename}_se")

    return getattr(rule_obj.output, outputfile).format(sample_name=sname, env=env)

def assign_bam_file(wildcards):
    sname = wildcards.sample_name
    env = get_sample_info_from_name(sname, samples, 'env')
    aligned_bams = config['aligned_bams']
    new_bam = assign_mapping_paired(wildcards, "filter_bam", "bamfile")
    if aligned_bams:
        return f"{RESULTS_DIR}/{env}/mapped/copied__{sname}.bam"
    else:
        return new_bam

def assign_chip_input(wildcards):
    """Get the control sample name for a ChIP sample.

    For merged/pseudo file_types, finds any replicate's control.
    For per-replicate file_types, uses direct Control column lookup.
    """
    sname = wildcards.sample_name
    if hasattr(wildcards, 'file_type') and wildcards.file_type in ['merged', 'pseudo1', 'pseudo2']:
        # For merged/pseudo peak calling, find any replicate's control
        # The sample_name here is an analysis_name; look up replicate Sample_IDs
        rep_sids = get_replicate_sample_ids(sname, samples)
        if rep_sids:
            for sid in rep_sids:
                ctrl = get_control_sample_id(sid, samples)
                if ctrl:
                    return ctrl
        # Fallback: try direct lookup in case it's a per-replicate name
        ctrl = get_control_sample_id(sname, samples)
        if ctrl:
            return ctrl
        raise ValueError(f"No control found for {sname}")
    else:
        ctrl = get_control_sample_id(sname, samples)
        if ctrl is None:
            raise ValueError(f"No control found for sample '{sname}'")
        return ctrl

def _peak_prefix_for_sample(sid, env):
    """Return the peak file prefix for a per-replicate sample."""
    if env == "ATAC":
        return "peaks_atac"
    paired = get_sample_field(sid, samples, 'paired')
    return "peaks_pe" if paired == "PE" else "peaks_se"


def assign_peak_files_for_idr(wildcards):
    """Return list of per-replicate peak files for IDR analysis."""
    sname = wildcards.sample_name  # analysis_name
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)
    rep_sids = get_replicate_sample_ids(sname, samples)
    return [f"{RESULTS_DIR}/{env}/peaks/{_peak_prefix_for_sample(sid, env)}__final__{sid}_peaks.{peaktype}Peak"
            for sid in rep_sids]

def input_peak_files_for_best_peaks(wildcards):
    """Return peak files for selecting best peaks via pseudorep consistency."""
    sname = wildcards.sample_name  # analysis_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)

    # Analysis-level prefix for merged/pseudo peaks (uses analysis-level paired)
    if env == "ATAC":
        prefix = "peaks_atac"
    elif paired == "PE":
        prefix = "peaks_pe"
    else:
        prefix = "peaks_se"

    rep_sids = get_replicate_sample_ids(sname, samples)

    if len(rep_sids) >= 2:
        result = [f"{RESULTS_DIR}/{env}/peaks/{prefix}__merged__{sname}_peaks.{peaktype}Peak",
                  f"{RESULTS_DIR}/{env}/peaks/{prefix}__pseudo1__{sname}_peaks.{peaktype}Peak",
                  f"{RESULTS_DIR}/{env}/peaks/{prefix}__pseudo2__{sname}_peaks.{peaktype}Peak",
                  f"{RESULTS_DIR}/{env}/peaks/idr_peaks__{sname}.bed"]
    else:
        one_sid = rep_sids[0]
        rep_prefix = _peak_prefix_for_sample(one_sid, env)
        result = [f"{RESULTS_DIR}/{env}/peaks/{rep_prefix}__final__{one_sid}_peaks.{peaktype}Peak",
                  f"{RESULTS_DIR}/{env}/peaks/{rep_prefix}__pseudo1__{one_sid}_peaks.{peaktype}Peak",
                  f"{RESULTS_DIR}/{env}/peaks/{rep_prefix}__pseudo2__{one_sid}_peaks.{peaktype}Peak",
                  f"{RESULTS_DIR}/empty.txt"]

    return result

def get_replicate_name(wildcards, pos):
    """Return the peak file path for replicate at position `pos`."""
    sname = wildcards.sample_name  # analysis_name
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)
    rep_sids = get_replicate_sample_ids(sname, samples)

    if pos >= len(rep_sids):
        return "missingrep"
    else:
        sid = rep_sids[pos]
        return f"{RESULTS_DIR}/{env}/peaks/{_peak_prefix_for_sample(sid, env)}__final__{sid}_peaks.{peaktype}Peak"

def get_replicate_pairs(wildcards):
    """Return colon-separated pairs of replicate Sample_IDs for IDR."""
    sname = wildcards.sample_name  # analysis_name
    rep_sids = get_replicate_sample_ids(sname, samples)
    pairs = []
    for i in range(0, len(rep_sids)):
        for j in range(i + 1, len(rep_sids)):
            pairs.append(f"{rep_sids[i]}:{rep_sids[j]}")
    return pairs

def define_chipseq_target_file(wildcards):
    tarname = config['motif_target_file_label']
    env = wildcards.env
    peak_file = wildcards.peak_file
    parts = peak_file.split("__")
    file_type = parts[0]
    if file_type == "selected_peaks":
        # analysis_name is everything after "selected_peaks__"
        spname = "__".join(parts[1:])
        inputfile = f"{RESULTS_DIR}/{env}/peaks/{file_type}__{spname}.bedPeak"
        ref_genome = get_sample_info_from_name(spname, analysis_samples, 'ref_genome')
        fasta = f"{GENOMES_DIR}/{ref_genome}/{ref_genome}.fa"
        if any(analysis_samples['sample_name'] == spname):
            return [inputfile, fasta]
    elif file_type == "idr_peaks":
        spname = "__".join(parts[1:])
        inputfile = f"{RESULTS_DIR}/{env}/peaks/{file_type}__{spname}.bed"
        ref_genome = get_sample_info_from_name(spname, analysis_samples, 'ref_genome')
        fasta = f"{GENOMES_DIR}/{ref_genome}/{ref_genome}.fa"
        if any(analysis_samples['sample_name'] == spname):
            return [inputfile, fasta]
    elif file_type.startswith("peaks_"):
        # Format: peaks_{pe|se|atac}__{file_cat}__{sample_name}_peaks
        # parts[0] = "peaks_pe" (or se/atac), parts[1] = file_cat (final/merged/pseudo),
        # rest = sample_name + "_peaks" suffix
        file_cat = parts[1]
        # Rejoin remaining parts and strip "_peaks" suffix
        rest = "__".join(parts[2:])
        if rest.endswith("_peaks"):
            sname = rest[:-len("_peaks")]
        else:
            sname = rest
        assay = _resolve_assay(sname)
        peaktype = get_peaktype_for_env(assay, env)
        inputfile = f"{RESULTS_DIR}/{env}/peaks/{file_type}__{file_cat}__{sname}_peaks.{peaktype}Peak"
        parsed = parse_sample_name(sname)
        ref_genome = parsed['ref_genome']
        fasta = f"{GENOMES_DIR}/{ref_genome}/{ref_genome}.fa"
        if any(samples['sample_name'] == sname) or any(analysis_samples['sample_name'] == sname):
            return [inputfile, fasta]
    elif peak_file == tarname:
        ref_genome = config['motif_ref_genome']
        fasta = f"{GENOMES_DIR}/{ref_genome}/{ref_genome}.fa"
        inputfile = config['motif_target_file']
        return [inputfile, fasta]
    else:
        return ValueError(
            f"{wildcards.peak_file} is unknown."
            "Options are either peakfiles generated by the pipeline"
            "or the value of 'motifs_target_file_name' in the options file"
            )

def define_input_manorm(wildcards, string):
    """Resolve inputs for pairwise differential peak analysis (MAnorm).

    sample1 and sample2 are analysis_names (wildcards).
    """
    sample1 = wildcards.sample1
    env1 = get_sample_info_from_name(sample1, analysis_samples, 'env')
    paired1 = get_sample_info_from_name(sample1, analysis_samples, 'paired')

    sample2 = wildcards.sample2
    env2 = get_sample_info_from_name(sample2, analysis_samples, 'env')
    paired2 = get_sample_info_from_name(sample2, analysis_samples, 'paired')

    if env1 == "ChIP":
        params = config['diffpeaks_params']['chip_pe'] if paired1 == "PE" and paired2 == "PE" else config['diffpeaks_params']['chip_se']
    elif env1 == "ATAC":
        params = config['diffpeaks_params']['ATAC_pe'] if paired1 == "PE" and paired2 == "PE" else config['diffpeaks_params']['ATAC_se']
    else:
        # Fallback for any other env
        params = config['diffpeaks_params'].get('chip_se', '')

    assay1 = _resolve_assay(sample1)
    assay2 = _resolve_assay(sample2)
    peaktype1 = get_peaktype_for_env(assay1, env1)
    peaktype2 = get_peaktype_for_env(assay2, env2)
    if peaktype1 != peaktype2:
        raise ValueError(f"{sample1} and {sample2} have different peaktypes.")
    else:
        peaktype = peaktype1

    peakfile1 = f"{RESULTS_DIR}/{env1}/peaks/selected_peaks__{sample1}.bedPeak"
    peakfile2 = f"{RESULTS_DIR}/{env2}/peaks/selected_peaks__{sample2}.bedPeak"

    rep_sids1 = get_replicate_sample_ids(sample1, samples)
    add1 = "shifted_" if env1 == "ATAC" else ""
    if len(rep_sids1) >= 2:
        bamfile1 = f"{RESULTS_DIR}/{env1}/mapped/{add1}merged__{sample1}.bam"
    else:
        bamfile1 = f"{RESULTS_DIR}/{env1}/mapped/{add1}final__{rep_sids1[0]}.bam"

    rep_sids2 = get_replicate_sample_ids(sample2, samples)
    add2 = "shifted_" if env2 == "ATAC" else ""
    if len(rep_sids2) >= 2:
        bamfile2 = f"{RESULTS_DIR}/{env2}/mapped/{add2}merged__{sample2}.bam"
    else:
        bamfile2 = f"{RESULTS_DIR}/{env2}/mapped/{add2}final__{rep_sids2[0]}.bam"

    if string == "peaks1":
        return peakfile1
    elif string == "peaks2":
        return peakfile2
    elif string == "reads1":
        return bamfile1
    elif string == "reads2":
        return bamfile2
    elif string == "params":
        return params
    elif string == "format":
        return peaktype

def define_logs_final_input(wildcards):
    """Collect all log files needed for peak stats summary."""
    log_files = []
    sname = wildcards.sample_name  # analysis_name
    paired = get_sample_info_from_name(sname, analysis_samples, 'paired')
    env = get_sample_info_from_name(sname, analysis_samples, 'env')
    assay = _resolve_assay(sname)
    peaktype = get_peaktype_for_env(assay, env)
    rep_sids = get_replicate_sample_ids(sname, samples)

    # ATAC's single peak-calling rule uses a BED-format input and does not
    # encode PE/SE in its log filename; ChIP's PE- and SE-specific rules do.
    peakcall_paired = "" if env == "ATAC" else None

    for sid in rep_sids:
        rep_paired = peakcall_paired if peakcall_paired is not None else get_sample_field(sid, samples, 'paired')
        log_files.append(return_log_chip(env, sid, f"final__{peaktype}peak_calling", rep_paired))
        log_files.append(return_log_chip(env, sid, "making_bigwig_final", ""))
        if env != "ATAC":
            log_files.append(return_log_chip(env, sid, "making_fingerprint_final", ""))

    merged_paired = peakcall_paired if peakcall_paired is not None else paired
    if len(rep_sids) >= 2:
        log_files.append(return_log_chip(env, sname, "IDR", ""))
        log_files.append(return_log_chip(env, sname, "merging_reps", ""))
        log_files.append(return_log_chip(env, sname, f"merged__{peaktype}peak_calling", merged_paired))
        log_files.append(return_log_chip(env, sname, f"pseudo1__{peaktype}peak_calling", merged_paired))
        log_files.append(return_log_chip(env, sname, f"pseudo2__{peaktype}peak_calling", merged_paired))
        log_files.append(return_log_chip(env, sname, "splitting_pseudreps", ""))
    else:
        one_sid = rep_sids[0]
        log_files.append(return_log_chip(env, one_sid, "splitting_pseudreps", ""))
        log_files.append(return_log_chip(env, one_sid, f"pseudo1__{peaktype}peak_calling", merged_paired))
        log_files.append(return_log_chip(env, one_sid, f"pseudo2__{peaktype}peak_calling", merged_paired))

    log_files.append(return_log_chip(env, sname, "selecting_best_peaks", ""))

    return log_files

def define_final_chip_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    motifs = config['motifs']
    motifs_allreps = config['motifs_allreps']
    trimmed_fastqs = config['trimmed_fastqs']
    aligned_bams = config['aligned_bams']
    map_files = []
    stat_files = []
    qc_files = []
    peak_files = []
    bigwig_files = []
    motif_files = []
    allrep_files = []

    # All ChIP samples (including controls) for mapping/QC
    filtered_rep_samples = samples[(samples['env'] == 'ChIP') & (samples['ref_genome'] == ref_genome)].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = row['sample_name']
        mname = row['mapped_name']
        paired = row['paired']
        env = row['env']
        if paired == "PE" and not aligned_bams:
            qc_files.append(f"{RESULTS_DIR}/{env}/reports/trim__{sname}__R1_fastqc.html")
            qc_files.append(f"{RESULTS_DIR}/{env}/reports/trim__{sname}__R2_fastqc.html")
            map_files.append(f"{RESULTS_DIR}/{env}/logs/process_chip_pe_sample__{mname}.log")
            if not trimmed_fastqs:
                qc_files.append(f"{RESULTS_DIR}/{env}/reports/raw__{sname}__R1_fastqc.html")
                qc_files.append(f"{RESULTS_DIR}/{env}/reports/raw__{sname}__R2_fastqc.html")
        elif paired == "SE" and not aligned_bams:
            qc_files.append(f"{RESULTS_DIR}/{env}/reports/trim__{sname}__R0_fastqc.html")
            map_files.append(f"{RESULTS_DIR}/{env}/logs/process_chip_se_sample__{mname}.log")
            if not trimmed_fastqs:
                qc_files.append(f"{RESULTS_DIR}/{env}/reports/raw__{sname}__R0_fastqc.html")

    # Samples that get peak calls and FC bigwigs. The test is structural — a
    # pulldown row must declare a Control — rather than "is anyone using this row
    # as their control?". That drops what cannot be peak-called (an Input/WCE/IgG
    # row no IP points at, i.e. an "orphan" control, and an IP with no Control;
    # both would make assign_chip_input raise "No control found" during DAG
    # build) while still admitting *dual-role* samples: a row may serve as
    # another sample's control AND be analysed in its own right, as long as it
    # has a Control of its own (e.g. an H3 ChIP used as its own target and as the
    # control for H3K9me2). See is_peak_call_target.
    filtered_rep_samples_no_input = peak_callable_rows(filtered_rep_samples)
    for _, row in filtered_rep_samples_no_input.iterrows():
        assay = row['Assay']
        peaktype = _peaktype_for_assay(assay)
        sname = row['sample_name']
        mname = row['mapped_name']
        paired = row['paired']
        env = row['env']
        bigwig_files.append(f"{RESULTS_DIR}/{env}/tracks/FC__final__{mname}.bw")
        if config.get("chip_fingerprint_plots", True):
            stat_files.append(f"{RESULTS_DIR}/{env}/plots/Fingerprint__final__{mname}.png")
        if paired == "PE":
            peak_files.append(f"{RESULTS_DIR}/{env}/peaks/peaks_pe__final__{mname}_peaks.{peaktype}Peak")
        else:
            peak_files.append(f"{RESULTS_DIR}/{env}/peaks/peaks_se__final__{mname}_peaks.{peaktype}Peak")

    # Analysis-level outputs (merged, IDR, best peaks)
    filtered_analysis_samples = analysis_samples[(analysis_samples['env'] == 'ChIP') & (analysis_samples['ref_genome'] == ref_genome)].copy()
    for _, row in filtered_analysis_samples.iterrows():
        aname = row['sample_name']  # analysis_name
        env = row['env']
        peak_files.append(f"{RESULTS_DIR}/{env}/peaks/selected_peaks__{aname}.bedPeak")
        rep_sids = get_replicate_sample_ids(aname, samples)
        if len(rep_sids) >= 2:
            bigwig_files.append(f"{RESULTS_DIR}/{env}/tracks/FC__merged__{aname}.bw")
            stat_files.append(f"{RESULTS_DIR}/{env}/chkpts/idr__{aname}.done")

    # Motif discovery for narrow-peak (TF-like) IP assays. The pre-refactor
    # pipeline keyed this on a dedicated "TF" env; that env is gone, so we now
    # select narrow assays (ChIP_narrow / CUT_RUN_narrow / CUT_TAG_narrow) via
    # ASSAY_TO_PEAKTYPE. Motif-around-summit is meaningful for narrow
    # (point-source) peaks, not broad histone domains.
    def _is_narrow(assay):
        return _peaktype_for_assay(assay) == "narrow"

    narrow_analysis = filtered_analysis_samples[
        filtered_analysis_samples['Assay'].map(_is_narrow)
    ]
    for _, row in narrow_analysis.iterrows():
        aname = row['sample_name']
        env = row['env']
        motif_files.append(f"{RESULTS_DIR}/{env}/chkpts/motifs__selected_peaks__{aname}.done")
        if len(get_replicate_sample_ids(aname, samples)) >= 2:
            motif_files.append(f"{RESULTS_DIR}/{env}/chkpts/motifs__idr_peaks__{aname}.done")

    # Per-replicate motif discovery (allreps path) for narrow rep samples.
    narrow_rep_samples = filtered_rep_samples_no_input[
        filtered_rep_samples_no_input['Assay'].map(_is_narrow)
    ]
    for _, row in narrow_rep_samples.iterrows():
        sname = row['sample_name']
        mname = row['mapped_name']
        env = row['env']
        pe = "pe" if row['paired'] == "PE" else "se"
        allrep_files.append(f"{RESULTS_DIR}/{env}/chkpts/motifs__peaks_{pe}__final__{mname}_peaks.done")

    # Pairwise differential peaks (MAnorm) for ChIP
    chip_analysis_samples = analysis_samples[(analysis_samples["env"] == "ChIP") & (analysis_samples["ref_genome"] == ref_genome)].copy()
    for ip_target, group in chip_analysis_samples.groupby("sample_type"):
        if len(group) >= 2:
            for a, b in combinations(group.itertuples(index=False), 2):
                sample1 = a.sample_name
                sample2 = b.sample_name
                assay1 = a.Assay
                assay2 = b.Assay
                peaktype1 = _peaktype_for_assay(assay1)
                peaktype2 = _peaktype_for_assay(assay2)
                if peaktype1 == peaktype2:
                    peak_files.append(f"{RESULTS_DIR}/ChIP/peaks/{sample1}_vs_{sample2}/{sample1}_vs_{sample2}_all_MAvalues.xls")

    results = map_files + bigwig_files

    if qc_option == "all":
        results += qc_files

    if analysis:
        results += peak_files + stat_files

    if motifs:
        results += motif_files

    if motifs_allreps:
        results += allrep_files

    return results

rule make_bt2_indices:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa"
    output:
        indices = directory(f"{GENOMES_DIR}/{{ref_genome}}/bt2_index")
    log:
        temp(os.path.join(REPO_FOLDER, RESULTS_DIR,"logs","bowtie_index_{ref_genome}.log"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nBuilding Bowtie2 index for {wildcards.ref_genome}\n"
        mkdir {config[genome_dir]}/{wildcards.ref_genome}/bt2_index
        bowtie2-build --threads {threads} "{input.fasta}" "{output.indices}/{wildcards.ref_genome}"
        }} 2>&1 | tee -a "{log}"
        """

rule make_chromap_index:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa"
    output:
        index = f"{GENOMES_DIR}/{{ref_genome}}/chromap_index/{{ref_genome}}.index"
    log:
        temp(os.path.join(REPO_FOLDER, RESULTS_DIR,"logs","chromap_index_{ref_genome}.log"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nBuilding Chromap index for {wildcards.ref_genome}\n"
        mkdir -p {config[genome_dir]}/{wildcards.ref_genome}/chromap_index
        chromap -i -t {threads} -r "{input.fasta}" -o "{output.index}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_bam_pe:
    input:
        fastq1 = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/fastq/trim__{base_sample(wildcards.sample_name)}__R1.fastq.gz",
        fastq2 = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/fastq/trim__{base_sample(wildcards.sample_name)}__R2.fastq.gz",
        index = get_mapping_index,
        fasta = get_mapping_fasta
    output:
        bamfile = temp(f"{RESULTS_DIR}/{{env}}/mapped/mapped_pe__{{sample_name}}.bam"),
        metrics_map = f"{RESULTS_DIR}/{{env}}/reports/bt2_pe__{{sample_name}}.txt",
        metrics_dup = f"{RESULTS_DIR}/{{env}}/reports/markdup_pe__{{sample_name}}.txt",
        metrics_flag = f"{RESULTS_DIR}/{{env}}/reports/flagstat_pe__{{sample_name}}.txt"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        aligner = lambda wildcards: get_effective_aligner(wildcards.env),
        map_option = lambda wildcards: get_mapping_strategy(wildcards.env),
        mapping_params = lambda wildcards: config['bt2_mapping_strategy'][get_mapping_strategy(wildcards.env)]['map_pe'],
        mapq_filter = lambda wildcards: get_mapq_filter(wildcards.env),
        sort_mem = lambda wildcards, resources, threads: f"{max(768, int(resources.mem_mb * 0.75 / threads))}M",
        chromap_extra = config.get('chromap_extra_flags', '')
    log:
        temp(return_log_chip("{env}","{sample_name}", "map_filter", "PE"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        set -o pipefail
        aligner="{params.aligner}"
        printf "\nMapping {params.sample_name} to {params.ref_genome} with $aligner ({params.map_option}) and filtering with samtools\n"
        printf "samtools sort memory: {params.sort_mem} per thread\n"
        samtools --version | head -1
        # Sort spills and the pre-markdup BAM stay in the per-job $TMPDIR.
        # samtools sort will not overwrite an existing temp chunk, yet exits 0
        # and silently drops those reads, so a retry must never find one.
        sorted_bam="$TMPDIR/sorted.bam"

        # NOTE: samtools view -q (MAPQ filter) must run *after* fixmate, not
        # before. fixmate requires name-collated input where mates are
        # adjacent. Filtering on MAPQ first can drop one mate but keep the
        # other, so fixmate then sees a stream where adjacent records are
        # from different pairs and incorrectly pairs them — visible as
        # large unexplained read loss on highly repetitive references
        # (e.g. >90% loss on ColCEN PE).
        if [[ "$aligner" == "chromap" ]]; then
            chromap --version
            chromap -t {threads} \
                -l 2000 --SAM -q 0 \
                {params.chromap_extra} \
                -r "{input.fasta}" -x "{input.index}" \
                -1 "{input.fastq1}" -2 "{input.fastq2}" \
                -o /dev/stdout 2> "{output.metrics_map}" \
            | samtools view -@ 2 -bh -F 256 \
            | samtools fixmate -@ 2 -m - - \
            | samtools view -@ 2 -bh -q {params.mapq_filter} \
            | samtools sort -@ {threads} -m {params.sort_mem} -T "$TMPDIR/sorted.sort" -o "$sorted_bam"
        else
            bowtie2 --version
            bowtie2 -p {threads} {params.mapping_params} \
                -x "{input.index}/{params.ref_genome}" \
                -1 "{input.fastq1}" -2 "{input.fastq2}" \
                2> "{output.metrics_map}" \
            | samtools view -@ 2 -bh -F 256 \
            | samtools fixmate -@ 2 -m - - \
            | samtools view -@ 2 -bh -q {params.mapq_filter} \
            | samtools sort -@ 2 -m {params.sort_mem} -T "$TMPDIR/sorted.sort" -o "$sorted_bam"
        fi

        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} \
            "$sorted_bam" "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f "$sorted_bam"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_bam_se:
    input:
        fastq = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/fastq/trim__{base_sample(wildcards.sample_name)}__R0.fastq.gz",
        index = get_mapping_index,
        fasta = get_mapping_fasta
    output:
        bamfile = temp(f"{RESULTS_DIR}/{{env}}/mapped/mapped_se__{{sample_name}}.bam"),
        metrics_map = f"{RESULTS_DIR}/{{env}}/reports/bt2_se__{{sample_name}}.txt",
        metrics_dup = f"{RESULTS_DIR}/{{env}}/reports/markdup_se__{{sample_name}}.txt",
        metrics_flag = f"{RESULTS_DIR}/{{env}}/reports/flagstat_se__{{sample_name}}.txt"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        aligner = lambda wildcards: get_effective_aligner(wildcards.env),
        map_option = lambda wildcards: get_mapping_strategy(wildcards.env),
        mapping_params = lambda wildcards: config['bt2_mapping_strategy'][get_mapping_strategy(wildcards.env)]['map_se'],
        mapq_filter = lambda wildcards: get_mapq_filter(wildcards.env),
        sort_mem = lambda wildcards, resources, threads: f"{max(768, int(resources.mem_mb * 0.75 / threads))}M",
        chromap_extra = config.get('chromap_extra_flags', '')
    log:
        temp(return_log_chip("{env}","{sample_name}", "map_filter", "SE"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        set -o pipefail
        aligner="{params.aligner}"
        printf "\nMapping {params.sample_name} to {params.ref_genome} with $aligner ({params.map_option}) and filtering with samtools\n"
        printf "samtools sort memory: {params.sort_mem} per thread\n"
        samtools --version | head -1
        # Sort spills and the pre-markdup BAM stay in the per-job $TMPDIR.
        # samtools sort will not overwrite an existing temp chunk, yet exits 0
        # and silently drops those reads, so a retry must never find one.
        sorted_bam="$TMPDIR/sorted.bam"

        if [[ "$aligner" == "chromap" ]]; then
            chromap --version
            # No samtools fixmate here: this is the SE branch (single-end has
            # no mates). markdup works on coord-sorted SE input directly.
            chromap -t {threads} \
                --SAM -q 0 \
                {params.chromap_extra} \
                -r "{input.fasta}" -x "{input.index}" \
                -1 "{input.fastq}" \
                -o /dev/stdout 2> "{output.metrics_map}" \
            | samtools view -@ 2 -bh -q {params.mapq_filter} -F 256 \
            | samtools sort -@ {threads} -m {params.sort_mem} -T "$TMPDIR/sorted.sort" -o "$sorted_bam"
        else
            bowtie2 --version
            bowtie2 -p {threads} {params.mapping_params} \
                -x "{input.index}/{params.ref_genome}" \
                -U "{input.fastq}" \
                2> "{output.metrics_map}" \
            | samtools view -@ 2 -bh -q {params.mapq_filter} -F 256 \
            | samtools sort -@ 2 -m {params.sort_mem} -T "$TMPDIR/sorted.sort" -o "$sorted_bam"
        fi

        samtools markdup -r -s -f "{output.metrics_dup}" -@ {threads} \
            "$sorted_bam" "{output.bamfile}"
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.bamfile}" > "{output.metrics_flag}"
        rm -f "$sorted_bam"
        }} 2>&1 | tee -a "{log}"
        """

rule make_mapping_stats_pe:
    input:
        metrics_trim = lambda wildcards: f"{RESULTS_DIR}/{{env}}/reports/trim_pe__{base_sample(wildcards.sample_name)}.json",
        metrics_map = f"{RESULTS_DIR}/{{env}}/reports/bt2_pe__{{sample_name}}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.env, (base_sample(wildcards.sample_name) if step in ("downloading", "trimming") else wildcards.sample_name), step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "map_filter"] ]
    output:
        stat_file = f"{RESULTS_DIR}/{{env}}/reports/summary_{{env}}_PE_mapping_stats_{{sample_name}}.txt",
        log = f"{RESULTS_DIR}/{{env}}/logs/process_chip_pe_sample__{{sample_name}}.log"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        group = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'line'),
        levels = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'levels_label'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'sample_type'),
        replicate = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'replicate'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'ref_genome'),
        trimmed_fastq = config['trimmed_fastqs']
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(python3 -c "import json; print(json.load(open('{input.metrics_trim}'))['summary']['before_filtering']['total_reads'] // 2)")
        else
            tot=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        fi
        # Detect aligner from metrics format
        if grep -q "overall alignment rate" "{input.metrics_map}"; then
            # bowtie2 format (numbers are read pairs for PE)
            filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
            multi=$(grep "aligned concordantly >1 times" "{input.metrics_map}" | awk '{{print $1}}')
            single=$(grep "aligned concordantly exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        else
            # chromap format — numbers have trailing dots, and count individual reads (not pairs)
            filt=$(grep "^Number of reads:" "{input.metrics_map}" | awk '{{sub(/\\.$/, ""); print $NF}}')
            mapped=$(grep "^Number of mapped reads:" "{input.metrics_map}" | awk '{{sub(/\\.$/, ""); print $NF}}')
            multi=$(grep "^Number of reads have multi-mappings:" "{input.metrics_map}" | awk '{{sub(/\\.$/, ""); print $NF}}')
            multi=${{multi:-0}}
            # Convert individual reads to pairs for PE
            filt=$((filt / 2))
            mapped=$((mapped / 2))
            multi=$((multi / 2))
            single=$((mapped - multi))
        fi
        allmap=$((multi+single))
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule make_mapping_stats_se:
    input:
        metrics_trim = lambda wildcards: f"{RESULTS_DIR}/{{env}}/reports/trim_se__{base_sample(wildcards.sample_name)}.json",
        metrics_map = f"{RESULTS_DIR}/{{env}}/reports/bt2_se__{{sample_name}}.txt",
        logs = lambda wildcards: [ return_log_chip(wildcards.env, (base_sample(wildcards.sample_name) if step in ("downloading", "trimming") else wildcards.sample_name), step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "map_filter"] ]
    output:
        stat_file = f"{RESULTS_DIR}/{{env}}/reports/summary_{{env}}_SE_mapping_stats_{{sample_name}}.txt",
        log = f"{RESULTS_DIR}/{{env}}/logs/process_chip_se_sample__{{sample_name}}.log"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        group = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'line'),
        levels = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'levels_label'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'sample_type'),
        replicate = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'replicate'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'ref_genome'),
        trimmed_fastq = config['trimmed_fastqs']
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(python3 -c "import json; print(json.load(open('{input.metrics_trim}'))['summary']['before_filtering']['total_reads'])")
        else
            tot=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
        fi
        # Detect aligner from metrics format
        if grep -q "overall alignment rate" "{input.metrics_map}"; then
            # bowtie2 format
            filt=$(grep "reads" "{input.metrics_map}" | awk '{{print $1}}')
            multi=$(grep "aligned >1 times" "{input.metrics_map}" | awk '{{print $1}}')
            single=$(grep "aligned exactly 1 time" "{input.metrics_map}" | awk '{{print $1}}')
        else
            # chromap format — numbers have trailing dots
            filt=$(grep "^Number of reads:" "{input.metrics_map}" | awk '{{sub(/\\.$/, ""); print $NF}}')
            mapped=$(grep "^Number of mapped reads:" "{input.metrics_map}" | awk '{{sub(/\\.$/, ""); print $NF}}')
            multi=$(grep "^Number of reads have multi-mappings:" "{input.metrics_map}" | awk '{{sub(/\\.$/, ""); print $NF}}')
            multi=${{multi:-0}}
            single=$((mapped - multi))
        fi
        allmap=$((multi+single))
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        """

rule dispatch_final_bam:
    input:
        assign_bam_file
    output:
        bam = maybe_temp(f"{RESULTS_DIR}/{{env}}/mapped/final__{{sample_name}}.bam", config.get('keep_final_bams', True)),
        bai = maybe_temp(f"{RESULTS_DIR}/{{env}}/mapped/final__{{sample_name}}.bam.bai", config.get('keep_final_bams', True)),
        touch = f"{RESULTS_DIR}/{{env}}/chkpts/map_{{env}}__{{sample_name}}.done"
    wildcard_constraints:
        env = "ChIP|ATAC"
    conda: CONDA_ENV_CHIP
    shell:
        """
        cp {input} {output.bam}
        samtools index -@ {threads} "{output.bam}"
        touch {output.touch}
        """

rule make_coverage_chip:
    input:
        bamfile = f"{RESULTS_DIR}/{{env}}/mapped/final__{{sample_name}}.bam",
        bambai = f"{RESULTS_DIR}/{{env}}/mapped/final__{{sample_name}}.bam.bai"
    output:
        bigwigcov = f"{RESULTS_DIR}/{{env}}/tracks/coverage__{{sample_name}}.bw"
    wildcard_constraints:
        env = "ChIP"
    params:
        binsize = config['chip_tracks']['binsize']
    conda: CONDA_ENV_CHIP
    shell:
        """
        bamCoverage -b {input.bamfile} -o {output.bigwigcov} -bs {params.binsize} -p {threads}
        """

rule make_bigwig_chip:
    input:
        ipfile = f"{RESULTS_DIR}/{{env}}/mapped/{{file_type}}__{{sample_name}}.bam",
        ipbai = f"{RESULTS_DIR}/{{env}}/mapped/{{file_type}}__{{sample_name}}.bam.bai",
        inputfile = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam",
        inputbai = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam.bai"
    output:
        bigwigfile = f"{RESULTS_DIR}/{{env}}/tracks/FC__{{file_type}}__{{sample_name}}.bw"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged"
    params:
        ipname = lambda wildcards: f"{wildcards.file_type}__{wildcards.sample_name}",
        inputname = lambda wildcards: f"{wildcards.file_type}__{assign_chip_input(wildcards)}",
        binsize = config['chip_tracks']['binsize'],
        params = config['chip_tracks']['params']
    log:
        temp(return_log_chip("{env}","{sample_name}", "making_bigwig_{file_type}", ""))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nMaking bigwig files for {params.ipname} (vs {params.inputname}) with deeptools version:\n"
        deeptools --version
        bamCompare -b1 {input.ipfile} -b2 {input.inputfile} -o {output.bigwigfile} -p {threads} --binSize {params.binsize} {params.params}
        }} 2>&1 | tee -a "{log}"
        """

rule make_fingerprint_plot:
    input:
        ipfile = f"{RESULTS_DIR}/{{env}}/mapped/final__{{sample_name}}.bam",
        ipbai = f"{RESULTS_DIR}/{{env}}/mapped/final__{{sample_name}}.bam.bai",
        inputfile = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/final__{assign_chip_input(wildcards)}.bam",
        inputbai = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/final__{assign_chip_input(wildcards)}.bam.bai",
        genome_stats = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/genome_stats.json"
    output:
        pngplot = f"{RESULTS_DIR}/{{env}}/plots/Fingerprint__final__{{sample_name}}.png"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        ipname = lambda wildcards: f"final__{wildcards.sample_name}",
        inputname = lambda wildcards: f"final__{assign_chip_input(wildcards)}",
        n_samples = config.get("fingerprint_samples", "auto")
    log:
        temp(return_log_chip("{env}","{sample_name}", "making_fingerprint_final", ""))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nPlotting fingerprint for {params.ipname} (vs {params.inputname}) with deeptools version:\n"
        deeptools --version

        # Auto-scale numberOfSamples to 50% of genomic bins
        n_samples="{params.n_samples}"
        if [[ "$n_samples" == "auto" ]]; then
            genome_size=$(python3 -c "import json; print(json.load(open('{input.genome_stats}'))['effective_size'])")
            bin_size=500
            n_samples=$(python3 -c "
total_bins = int($genome_size / $bin_size)
print(int(total_bins * 0.5))
")
            printf "Auto-scaled --numberOfSamples to %s (genome: %s bp, bin: %s bp)\n" "$n_samples" "$genome_size" "$bin_size"
        fi

        plotFingerprint -b {input.ipfile} {input.inputfile} \
            -o {output.pngplot} -p {threads} \
            -n "$n_samples" \
            -l {params.ipname} {params.inputname}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_pe:
    input:
        ipfile = f"{RESULTS_DIR}/{{env}}/mapped/{{file_type}}__{{sample_name}}.bam",
        inputfile = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam",
        genome_stats = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/genome_stats.json",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes",
        converter = os.path.join(REPO_FOLDER, "workflow", "scripts", "convert_peaks.py")
    output:
        peakfile = f"{RESULTS_DIR}/{{env}}/peaks/peaks_pe__{{file_type}}__{{sample_name}}_peaks.{{peaktype}}Peak"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged|pseudo1|pseudo2",
        peaktype = "broad|narrow"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        inputname = lambda wildcards: assign_chip_input(wildcards),
        peaktype = lambda wildcards: _checked_peaktype(
            wildcards.sample_name, wildcards.env, wildcards.peaktype),
        caller = lambda wildcards: _resolve_peakcaller(wildcards.sample_name),
        is_cut = lambda wildcards: "1" if _resolve_assay(wildcards.sample_name).startswith(("CUT_RUN", "CUT_TAG")) else "0",
        file_type = lambda wildcards: wildcards.file_type,
        env = lambda wildcards: wildcards.env,
        chip_macs2_params = config["chip_callpeaks"]["params"],
        cut_macs2_params_pe = config.get("cut_callpeaks", {}).get("macs2_params_pe", "--keep-dup all --nomodel"),
        epic2_params = config.get("cut_callpeaks", {}).get("epic2_params", "--bin-size 150 --gaps-allowed 2"),
        seacr_norm = config.get("cut_callpeaks", {}).get("seacr_norm", "non"),
        seacr_threshold = config.get("cut_callpeaks", {}).get("seacr_threshold", "stringent"),
        genomesize_override = lambda wildcards: config["genomes"][parse_sample_name(wildcards.sample_name)['ref_genome']].get('genomesize', '')
    log:
        temp(return_log_chip("{env}","{sample_name}", "{file_type}__{peaktype}peak_calling", "PE"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        set -euo pipefail
        if [[ -n "{params.genomesize_override}" ]]; then
            gsize="{params.genomesize_override}"
        else
            gsize=$(python3 -c "import json; print(json.load(open('{input.genome_stats}'))['effective_size'])")
        fi
        outdir="{config[output_dir]}/{params.env}/peaks"
        prefix="peaks_pe__{params.file_type}__{params.sample_name}"
        printf "\\nCalling {params.peaktype} peaks (caller={params.caller}) for paired-end {params.sample_name} (vs {params.inputname})\\n"
        case "{params.caller}" in
        macs2)
            if [[ "{params.is_cut}" == "1" ]]; then
                macs_params="{params.cut_macs2_params_pe}"
            else
                macs_params="{params.chip_macs2_params}"
            fi
            extra=""
            if [[ "{params.peaktype}" == "broad" ]]; then extra="--broad"; fi
            macs2 --version
            macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAMPE -g $gsize $macs_params -n $prefix --outdir $outdir $extra
            ;;
        epic2)
            if [[ "{params.peaktype}" != "broad" ]]; then
                echo "[X] epic2 only supports broad peak calling; pick seacr or macs2 for narrow." >&2
                exit 2
            fi
            efrac=$(python3 -c "import json; s=json.load(open('{input.genome_stats}')); print(s['effective_size']/s['total_bases'])")
            raw="$outdir/${{prefix}}.epic2.tsv"
            epic2 --version
            epic2 -t {input.ipfile} -c {input.inputfile} \
                  --chromsizes {input.chrom_sizes} \
                  --effective-genome-fraction "$efrac" \
                  {params.epic2_params} \
                  --output "$raw"
            python3 {input.converter} --caller epic2 --peaktype broad "$raw" "{output.peakfile}"
            rm -f "$raw"
            ;;
        seacr)
            tmpd=$(mktemp -d)
            trap 'rm -rf "$tmpd"' EXIT
            for tag in ip ctrl; do
                if [[ "$tag" == "ip" ]]; then bam="{input.ipfile}"; else bam="{input.inputfile}"; fi
                samtools sort -n -@ {threads} -O bam -T "$tmpd/sort_$tag" "$bam" \
                  | bedtools bamtobed -bedpe -i - 2>/dev/null \
                  | awk 'BEGIN{{OFS="\\t"}} $1==$4 && $6-$2 < 1000 {{print $1, $2, $6}}' \
                  | sort -k1,1 -k2,2n -k3,3n \
                  | bedtools genomecov -bg -i - -g {input.chrom_sizes} > "$tmpd/$tag.bg"
            done
            raw="$outdir/${{prefix}}.seacr.bed"
            SEACR_1.3.sh "$tmpd/ip.bg" "$tmpd/ctrl.bg" {params.seacr_norm} {params.seacr_threshold} "$outdir/$prefix"
            mv "$outdir/${{prefix}}.{params.seacr_threshold}.bed" "$raw"
            python3 {input.converter} --caller seacr --peaktype {params.peaktype} "$raw" "{output.peakfile}"
            rm -f "$raw"
            ;;
        *)
            echo "[X] Unknown peak caller: {params.caller}" >&2
            exit 1
            ;;
        esac
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_se:
    input:
        ipfile = f"{RESULTS_DIR}/{{env}}/mapped/{{file_type}}__{{sample_name}}.bam",
        inputfile = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/{wildcards.file_type}__{assign_chip_input(wildcards)}.bam",
        genome_stats = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/genome_stats.json",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes",
        converter = os.path.join(REPO_FOLDER, "workflow", "scripts", "convert_peaks.py")
    output:
        peakfile = f"{RESULTS_DIR}/{{env}}/peaks/peaks_se__{{file_type}}__{{sample_name}}_peaks.{{peaktype}}Peak"
    wildcard_constraints:
        env = "ChIP",
        file_type = "final|merged|pseudo1|pseudo2",
        peaktype = "broad|narrow"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        inputname = lambda wildcards: assign_chip_input(wildcards),
        peaktype = lambda wildcards: _checked_peaktype(
            wildcards.sample_name, wildcards.env, wildcards.peaktype),
        caller = lambda wildcards: _resolve_peakcaller(wildcards.sample_name),
        is_cut = lambda wildcards: "1" if _resolve_assay(wildcards.sample_name).startswith(("CUT_RUN", "CUT_TAG")) else "0",
        file_type = lambda wildcards: wildcards.file_type,
        env = lambda wildcards: wildcards.env,
        chip_macs2_params = config["chip_callpeaks"]["params"],
        cut_macs2_params_se = config.get("cut_callpeaks", {}).get("macs2_params_se", "--keep-dup all --nomodel --shift -75 --extsize 150"),
        epic2_params = config.get("cut_callpeaks", {}).get("epic2_params", "--bin-size 150 --gaps-allowed 2"),
        seacr_norm = config.get("cut_callpeaks", {}).get("seacr_norm", "non"),
        seacr_threshold = config.get("cut_callpeaks", {}).get("seacr_threshold", "stringent"),
        genomesize_override = lambda wildcards: config["genomes"][parse_sample_name(wildcards.sample_name)['ref_genome']].get('genomesize', '')
    log:
        temp(return_log_chip("{env}","{sample_name}", "{file_type}__{peaktype}peak_calling", "SE"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        set -euo pipefail
        if [[ -n "{params.genomesize_override}" ]]; then
            gsize="{params.genomesize_override}"
        else
            gsize=$(python3 -c "import json; print(json.load(open('{input.genome_stats}'))['effective_size'])")
        fi
        outdir="{config[output_dir]}/{params.env}/peaks"
        prefix="peaks_se__{params.file_type}__{params.sample_name}"
        printf "\\nCalling {params.peaktype} peaks (caller={params.caller}) for single-end {params.sample_name} (vs {params.inputname})\\n"
        case "{params.caller}" in
        macs2)
            if [[ "{params.is_cut}" == "1" ]]; then
                macs_params="{params.cut_macs2_params_se}"
            else
                macs_params="{params.chip_macs2_params}"
            fi
            extra=""
            if [[ "{params.peaktype}" == "broad" ]]; then extra="--broad"; fi
            macs2 --version
            macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAM -g $gsize $macs_params -n $prefix --outdir $outdir $extra
            ;;
        epic2)
            if [[ "{params.peaktype}" != "broad" ]]; then
                echo "[X] epic2 only supports broad peak calling; pick seacr or macs2 for narrow." >&2
                exit 2
            fi
            efrac=$(python3 -c "import json; s=json.load(open('{input.genome_stats}')); print(s['effective_size']/s['total_bases'])")
            raw="$outdir/${{prefix}}.epic2.tsv"
            epic2 --version
            epic2 -t {input.ipfile} -c {input.inputfile} \
                  --chromsizes {input.chrom_sizes} \
                  --effective-genome-fraction "$efrac" \
                  {params.epic2_params} \
                  --output "$raw"
            python3 {input.converter} --caller epic2 --peaktype broad "$raw" "{output.peakfile}"
            rm -f "$raw"
            ;;
        seacr)
            tmpd=$(mktemp -d)
            trap 'rm -rf "$tmpd"' EXIT
            for tag in ip ctrl; do
                if [[ "$tag" == "ip" ]]; then bam="{input.ipfile}"; else bam="{input.inputfile}"; fi
                bedtools genomecov -bg -ibam "$bam" > "$tmpd/$tag.bg"
            done
            raw="$outdir/${{prefix}}.seacr.bed"
            SEACR_1.3.sh "$tmpd/ip.bg" "$tmpd/ctrl.bg" {params.seacr_norm} {params.seacr_threshold} "$outdir/$prefix"
            mv "$outdir/${{prefix}}.{params.seacr_threshold}.bed" "$raw"
            python3 {input.converter} --caller seacr --peaktype {params.peaktype} "$raw" "{output.peakfile}"
            rm -f "$raw"
            ;;
        *)
            echo "[X] Unknown peak caller: {params.caller}" >&2
            exit 1
            ;;
        esac
        }} 2>&1 | tee -a "{log}"
        """

rule idr_analysis_replicates:
    input:
        peak_file = lambda wildcards: assign_peak_files_for_idr(wildcards)
    output:
        touch = f"{RESULTS_DIR}/{{env}}/chkpts/idr__{{sample_name}}.done",
        merged = f"{RESULTS_DIR}/{{env}}/peaks/idr_peaks__{{sample_name}}.bed"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        peaktype = lambda wildcards: get_peaktype_for_env(_resolve_assay(wildcards.sample_name), wildcards.env),
        paired = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'paired'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome'),
        env = lambda wildcards: wildcards.env,
        replicate_pairs = lambda wildcards: get_replicate_pairs(wildcards),
        rep_sids = lambda wildcards: get_replicate_sample_ids(wildcards.sample_name, samples),
        sid_prefixes = lambda wildcards: " ".join(
            f"{sid}={_peak_prefix_for_sample(sid, wildcards.env)}"
            for sid in get_replicate_sample_ids(wildcards.sample_name, samples)
        )
    log:
        temp(return_log_chip("{env}","{sample_name}", "IDR", ""))
    conda: CONDA_ENV_IDR
    shell:
        """
        {{
        printf "\nLooping over each unique pair of biological replicates for {params.sname} to perform IDR with:\n"
        idr --version
        # Build associative array mapping sample_id → peak prefix (pe/se/atac)
        declare -A sid_pre
        for entry in {params.sid_prefixes}; do
            sid_pre[${{entry%%=*}}]=${{entry#*=}}
        done
        # Analysis-level prefix for merged/pseudo files
        if [[ "{params.env}" == "ATAC" ]]; then
            analysis_pre="peaks_atac"
        elif [[ "{params.paired}" == "PE" ]]; then
            analysis_pre="peaks_pe"
        else
            analysis_pre="peaks_se"
        fi
        temp="{config[output_dir]}/{params.env}/peaks/temp_idr_peaks__{params.sname}.bed"
        while read chr max; do
            printf "${{chr}}\t1\t${{max}}\n" >> "${{temp}}"
        done < "{config[genome_dir]}/{params.ref_genome}/chrom.sizes"
        mkdir -p {config[output_dir]}/{params.env}/plots/
        idr_failed=false
        for pair in {params.replicate_pairs}; do
            rep1_sid=$(echo ${{pair}} | cut -d":" -f1)
            rep2_sid=$(echo ${{pair}} | cut -d":" -f2)
            pre1=${{sid_pre[$rep1_sid]}}
            pre2=${{sid_pre[$rep2_sid]}}
            file1="{config[output_dir]}/{params.env}/peaks/${{pre1}}__final__${{rep1_sid}}_peaks.{params.peaktype}Peak"
            file2="{config[output_dir]}/{params.env}/peaks/${{pre2}}__final__${{rep2_sid}}_peaks.{params.peaktype}Peak"
            outfile="{config[output_dir]}/{params.env}/peaks/idr_${{analysis_pre}}__{params.sname}__${{rep1_sid}}_vs_${{rep2_sid}}_peaks.{params.peaktype}Peak"

            n1=$(wc -l < "${{file1}}")
            n2=$(wc -l < "${{file2}}")
            printf "\nPerforming IDR for ${{rep1_sid}} (%d peaks) vs ${{rep2_sid}} (%d peaks)\n" "${{n1}}" "${{n2}}"

            if idr --input-file-type {params.peaktype}Peak --output-file-type {params.peaktype}Peak \
                    --samples ${{file1}} ${{file2}} -o ${{outfile}} \
                    -l {config[output_dir]}/{params.env}/reports/idr_{params.sname}.log --plot 2>&1; then
                mv "${{outfile}}.png" {config[output_dir]}/{params.env}/plots/ 2>/dev/null || true
                filtered="${{outfile}}.filtered"
                awk -v OFS="\t" '$5>=540 {{print $1,$2,$3}}' ${{outfile}} | sort -k1,1 -k2,2n > "${{filtered}}"
                new="${{temp}}.new"
                bedtools intersect -a ${{temp}} -b ${{filtered}} > "${{new}}"
                mv "${{new}}" "${{temp}}"
            else
                idr_failed=true
                printf "\n"
                printf "========================================================================\n"
                printf "WARNING: IDR failed for %s vs %s\n" "${{rep1_sid}}" "${{rep2_sid}}"
                printf "This is an analytical outcome, not a pipeline error.\n"
                printf "Common causes: too few peaks (<20 after merge), low replicate concordance.\n"
                printf "Rep1 peaks: %d, Rep2 peaks: %d\n" "${{n1}}" "${{n2}}"
                printf "IDR-based peak filtering will be skipped for this comparison.\n"
                printf "Peak selection will rely on merged + pseudoreplicate consistency.\n"
                printf "========================================================================\n"
            fi
        done

        if [[ "${{idr_failed}}" == "true" ]]; then
            printf "\nIDR analysis incomplete for {params.sname} — creating empty IDR output\n"
            touch {output.merged}
            touch "{config[output_dir]}/{params.env}/peaks/idr_peaks__{params.sname}.{params.peaktype}Peak"
        else
            cat ${{temp}} > {output.merged}
            bedtools intersect -a ${{file2}} -b ${{temp}} -u > "{config[output_dir]}/{params.env}/peaks/idr_peaks__{params.sname}.{params.peaktype}Peak"
        fi
        rm -f ${{temp}} "${{temp}}.new" "${{outfile}}.filtered"
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule merging_bam_replicates:
    input:
        bamfiles = lambda wildcards: [
            f"{RESULTS_DIR}/{wildcards.env}/mapped/final__{sid}.bam"
            for sid in get_replicate_sample_ids(wildcards.sample_name, samples)
        ]
    output:
        mergefile = maybe_temp(f"{RESULTS_DIR}/{{env}}/mapped/merged__{{sample_name}}.bam", config.get('keep_merged_bams', False)),
        mergebai = maybe_temp(f"{RESULTS_DIR}/{{env}}/mapped/merged__{{sample_name}}.bam.bai", config.get('keep_merged_bams', False))
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env
    log:
        temp(return_log_chip("{env}","{sample_name}", "merging_reps", ""))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
        samtools merge -u -@ {threads} - {input.bamfiles} | samtools sort -@ {threads} -T "$TMPDIR/merge.sort" -o {output.mergefile}
        samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """

rule making_pseudo_replicates:
    input:
        bamfile = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'final'}__{wildcards.sample_name}.bam"
    output:
        temp_pseudo1 = temp(f"{RESULTS_DIR}/{{env}}/mapped/temp_pseudo1__{{sample_name}}.bam"),
        temp_pseudo2 = temp(f"{RESULTS_DIR}/{{env}}/mapped/temp_pseudo2__{{sample_name}}.bam"),
        pseudo1 = temp(f"{RESULTS_DIR}/{{env}}/mapped/pseudo1__{{sample_name}}.bam"),
        pseudo1_bai = temp(f"{RESULTS_DIR}/{{env}}/mapped/pseudo1__{{sample_name}}.bam.bai"),
        pseudo2 = temp(f"{RESULTS_DIR}/{{env}}/mapped/pseudo2__{{sample_name}}.bam"),
        pseudo2_bai = temp(f"{RESULTS_DIR}/{{env}}/mapped/pseudo2__{{sample_name}}.bam.bai")
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env
    log:
        temp(return_log_chip("{env}","{sample_name}", "splitting_pseudreps", ""))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nSplitting {params.sname} in two pseudo-replicates\n"
        samtools view -b -h -s 1.5 -@ {threads} -U {output.temp_pseudo2} -o {output.temp_pseudo1} {input.bamfile}
		samtools sort -@ {threads} -T "$TMPDIR/pseudo1.sort" -o {output.pseudo1} {output.temp_pseudo1}
		samtools sort -@ {threads} -T "$TMPDIR/pseudo2.sort" -o {output.pseudo2} {output.temp_pseudo2}
        # Index alongside the BAM so downstream rules (e.g. atac_shift_bam,
        # which calls alignmentSieve) can require a fresh .bai as input.
        samtools index -@ {threads} {output.pseudo1}
        samtools index -@ {threads} {output.pseudo2}
        }} 2>&1 | tee -a "{log}"
        """

rule create_empty_file:
    output:
        f"{RESULTS_DIR}/empty.txt"
    localrule: True
    shell:
        "touch {output}"

rule best_peaks_pseudoreps:
    input:
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome')}/chrom.sizes",
        peakfiles = input_peak_files_for_best_peaks
    output:
        bestpeaks = f"{RESULTS_DIR}/{{env}}/peaks/selected_peaks__{{sample_name}}.bedPeak",
        stats_pseudoreps = temp(f"{RESULTS_DIR}/{{env}}/reports/stats_pseudoreps__{{sample_name}}.txt")
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        env = lambda wildcards: wildcards.env,
        peaktype = lambda wildcards: get_peaktype_for_env(_resolve_assay(wildcards.sample_name), wildcards.env)
    log:
        temp(return_log_chip("{env}","{sample_name}", "selecting_best_peaks", ""))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        printf "\nIntersecting merged peaks ({input.peakfiles[0]}), and both pseudo replicates to select the best peaks for {params.sname}\n"
        awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[0]} | sort -k1,1 -k2,2n -u > "{config[output_dir]}/{params.env}/peaks/temp_{params.sname}_merged.bed"
		awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[1]} | sort -k1,1 -k2,2n -u > "{config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudo1.bed"
		awk -v OFS="\t" '{{print $1,$2,$3}}' {input.peakfiles[2]} | sort -k1,1 -k2,2n -u > "{config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudo2.bed"
		bedtools intersect -a {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudo1.bed -b {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudo2.bed > "{config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudos.bed"
		bedtools intersect -a {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_merged.bed -b {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudo1.bed -u > "{config[output_dir]}/{params.env}/peaks/temp_{params.sname}_selected.bed"
		bedtools intersect -a {input.peakfiles[0]} -b {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_selected.bed -u > "{config[output_dir]}/{params.env}/peaks/selected_peaks__{params.sname}.{params.peaktype}Peak"
        printf "\nGetting best quality peaks peaks\n"
        ## Note: If broadpeak, an additional "summit" column will be added for potential downstream processes, which only represent the middle of the peak, not its actual summit.
        sort -k1,1 -k2,2n -k5nr {config[output_dir]}/{params.env}/peaks/selected_peaks__{params.sname}.{params.peaktype}Peak | awk -v OFS="\t" -v t={params.peaktype} '{{a=$1":"$2":"$3; if (a!=n) {{if (t=="broad") $10=int(($3-$2)/2); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}} n=$1":"$2":"$3}}' | bedtools sort -g {input.chrom_sizes} > "{output.bestpeaks}"
        printf "\nExtracting peak stats for {params.sname}\n"
        merged=$(wc -l {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_merged.bed | cut -d" " -f1)
		pseudos=$(awk '{{print $1,$2,$3}}' {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_pseudos.bed | sort -k1,1 -k2,2n -u | wc -l)
		selected=$(cat {config[output_dir]}/{params.env}/peaks/temp_{params.sname}_selected.bed | sort -k1,1 -k2,2n -u | wc -l)
        if [[ "{input.peakfiles[3]}" == "{config[output_dir]}/empty.txt" ]]; then
            idr="0"
        else
            idr=$(wc -l {input.peakfiles[3]} | cut -d" " -f1)
        fi
		printf "Merged=${{merged}}\nPseudos=${{pseudos}}\nIDR=${{idr}}\nSelected=${{selected}}\n" > "{output.stats_pseudoreps}"
        rm -f "{config[output_dir]}/{params.env}/peaks/temp_{params.sname}"*
        }} 2>&1 | tee -a "{log}"
        """

rule make_peak_stats:
    input:
        logs = lambda wildcards: define_logs_final_input(wildcards),
        stats_pseudoreps = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/reports/stats_pseudoreps__{wildcards.sample_name}.txt"
        ## maybe a better solution is to append a stat file with wc -l as they are generated, or to create a new stat file for each file, then accessible in bash by regex on the samplename
    output:
        stat_file = f"{RESULTS_DIR}/{{env}}/reports/summary_{{env}}_peak_stats_{{sample_name}}.txt",
        log = f"{RESULTS_DIR}/{{env}}/logs/called_peaks__{{sample_name}}.log"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        sname = lambda wildcards: wildcards.sample_name,
        paired = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'paired'),
        peaktype = lambda wildcards: get_peaktype_for_env(_resolve_assay(wildcards.sample_name), wildcards.env),
        group = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'line'),
        levels = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'levels_label'),
        sample_type = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'sample_type'),
        ref_genome = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'ref_genome'),
        env = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'env'),
        tf_name = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, analysis_samples, 'extra_info'),
        rep1 = lambda wildcards: get_replicate_name(wildcards, 0),
        rep2 = lambda wildcards: get_replicate_name(wildcards, 1)
    shell:
        """
        nrep1=$(awk '{{print $1,$2,$3}}' {params.rep1} | sort -k1,1 -k2,2n -u | wc -l)
        if [[ "{params.rep2}" == "missingrep" ]]; then
            nrep2=0
        else
            nrep2=$(awk '{{print $1,$2,$3}}' {params.rep2} | sort -k1,1 -k2,2n -u | wc -l)
        fi
        merged=$(grep "Merged" {input.stats_pseudoreps} | cut -d"=" -f2)
        pseudos=$(grep "Pseudos" {input.stats_pseudoreps} | cut -d"=" -f2)
        idr=$(grep "IDR" {input.stats_pseudoreps} | cut -d"=" -f2)
        selected=$(grep "Selected" {input.stats_pseudoreps} | cut -d"=" -f2)
        printf "Group\tLevels\tMark\tReference_genome\tPeaks_in_Rep1\tPeaks_in_Rep2\tPeaks_in_merged\tPeaks_in_pseudo_reps\tPeaks_in_idr\tSelected_peaks\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v m={params.sample_type} -v r={params.ref_genome} -v a=${{nrep1}} -v b=${{nrep2}} -v c=${{merged}} -v d=${{pseudos}} -v e=${{idr}} -v f=${{selected}} 'BEGIN {{if (c==0) {{x=a}} else {{x=c}}; if (x==0) print l,t,m,r,a,b,c,d,e,f" (no%)"; else print l,t,m,r,a,b,c,d,e,f" ("f/x*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        """

rule find_motifs_in_file:
    input:
        define_chipseq_target_file
    output:
        temp_bed = temp(f"{RESULTS_DIR}/{{env}}/motifs/temp_regions_{{peak_file}}.bed"),
        temp_fa = temp(f"{RESULTS_DIR}/{{env}}/motifs/temp_regions_{{peak_file}}.fa"),
        touch = f"{RESULTS_DIR}/{{env}}/chkpts/motifs__{{peak_file}}.done"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        env = lambda wildcards: wildcards.env,
        peak_file = lambda wildcards: wildcards.peak_file,
        jaspar_db = resolve_repo_data(config['jaspar_db']),
        min_region_size = config.get('motif_params', {}).get('min_region_size', 8),
        n_motifs = config.get('motif_params', {}).get('n_motifs', 10),
        min_width = config.get('motif_params', {}).get('min_width', 6),
        max_width = config.get('motif_params', {}).get('max_width', 20)
    log:
        temp(return_log_chip("{env}","{peak_file}", "motifs", ""))
    conda: CONDA_ENV_IDR
    shell:
        """
        {{
        peakfile="{input[0]}"
        ext=${{peakfile##*.}}
        if [[ "${{ext}}" == "narrowPeak" ]]; then
            printf "\nGetting peak fasta sequences around the summit for narrowPeak file: ${{peakfile}}\n"
            sort -k1,1 -k2,2n -k5,5nr ${{peakfile}} | awk -v OFS="\t" '{{a=$1":"$2":"$3; if (a!=n) {{s=$2+$10; print $1,s-200,s+200,$4;}} n=$1":"$2":"$3}}' > {output.temp_bed}
        elif [[ "${{ext}}" == "broadPeak" ]]; then
            printf "\nGetting peak fasta sequences around the middle for broadPeak file: ${{peakfile}}\n"
            sort -k1,1 -k2,2n -k5,5nr ${{peakfile}} | awk -v OFS="\t" '{{s=int(($2+$3)/2); t=($3-$2); if (t<500) print $1,$2,$3,$4; else print $1,s-200,s+200,$4}}' > {output.temp_bed}
        elif [[ "${{ext}}" == "bedPeak" ]]; then
            printf "\nGetting peak fasta sequences for bed file: ${{peakfile}}\n"
            cat ${{peakfile}} | awk -v OFS="\t" '{{s=$2+$10; print $1,s-200,s+200,$4}}' > {output.temp_bed}
        else
            printf "\nGetting peak fasta sequences for unknown file format: ${{peakfile}}\n"
            cat ${{peakfile}} | awk -v OFS="\t" '{{s=int(($2+$3)/2); t=($3-$2); if (t<500) print $1,$2,$3,$4; else print $1,s-200,s+200,$4}}' > {output.temp_bed}
        fi
        head {output.temp_bed}
        # Clamp negative starts (peaks near a contig edge) and drop regions
        # below the configured minimum size. meme-chip errors when sequences
        # are shorter than the motif width, so filter degenerate regions up front.
        awk -v OFS="\t" -v min={params.min_region_size} '{{ if ($2<0) $2=0; if (($3-$2) >= min) print }}' {output.temp_bed} > {output.temp_bed}.filt
        mv {output.temp_bed}.filt {output.temp_bed}
        nregions=$(wc -l < {output.temp_bed})
        printf "\n%s input regions >= {params.min_region_size}bp for {params.peak_file}\n" "$nregions"
        bedtools getfasta -name -fi {input[1]} -bed {output.temp_bed} > {output.temp_fa}
        head {output.temp_fa}
        if [[ -s {output.temp_fa} ]]; then
            printf "\nGetting motifs for {params.peak_file} with meme version:\n"
            meme -version
            meme-chip -oc {config[output_dir]}/{params.env}/motifs/{params.peak_file}/meme -meme-p {threads} \
                -meme-nmotifs {params.n_motifs} -streme-nmotifs {params.n_motifs} \
                -minw {params.min_width} -maxw {params.max_width} {output.temp_fa}
            if [[ -s {config[output_dir]}/{params.env}/motifs/{params.peak_file}/meme/combined.meme ]]; then
                printf "\nLooking for similar motifs in JASPAR database with tomtom\n"
                tomtom -oc {config[output_dir]}/{params.env}/motifs/{params.peak_file}/tomtom/ {config[output_dir]}/{params.env}/motifs/{params.peak_file}/meme/combined.meme {params.jaspar_db}
            fi
        else
            printf "\nWARNING: no input regions >= {params.min_region_size}bp for {params.peak_file}; skipping motif discovery.\n"
        fi
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule perform_pairwise_diff_peaks:
    input:
        peak_file1 = lambda wildcards: define_input_manorm(wildcards, "peaks1"),
        peak_file2 = lambda wildcards: define_input_manorm(wildcards, "peaks2"),
        read_file1 = lambda wildcards: define_input_manorm(wildcards, "reads1"),
        read_file2 = lambda wildcards: define_input_manorm(wildcards, "reads2")
    output:
        result = f"{RESULTS_DIR}/{{env}}/peaks/{{sample1}}_vs_{{sample2}}/{{sample1}}_vs_{{sample2}}_all_MAvalues.xls"
    wildcard_constraints:
        env = "ChIP|ATAC"
    params:
        diffpeaks = lambda wildcards: define_input_manorm(wildcards, "params"),
        peakformat = lambda wildcards: define_input_manorm(wildcards, "format"),
        output_folder = lambda wildcards: f"{RESULTS_DIR}/{wildcards.env}/peaks/{wildcards.sample1}_vs_{wildcards.sample2}"
    log:
        temp(return_log_chip("{env}","{sample1}_vs_{sample2}", "diff_peaks", ""))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        n1=$(wc -l < {input.peak_file1})
        n2=$(wc -l < {input.peak_file2})
        min_peaks=2

        if [[ "$n1" -lt "$min_peaks" || "$n2" -lt "$min_peaks" ]]; then
            printf "\n"
            printf "========================================================================\n"
            printf "WARNING: Skipping differential peak analysis for {wildcards.sample1} vs {wildcards.sample2}\n"
            printf "Insufficient peaks for meaningful comparison (minimum %d required per sample).\n" "$min_peaks"
            printf "{wildcards.sample1}: %d peaks, {wildcards.sample2}: %d peaks\n" "$n1" "$n2"
            printf "This is an analytical outcome, not a pipeline error.\n"
            printf "========================================================================\n"
            mkdir -p {params.output_folder}
            printf "# MAnorm skipped: insufficient peaks ({wildcards.sample1}=%d, {wildcards.sample2}=%d)\n" "$n1" "$n2" > {output.result}
        else
            if [[ "{params.peakformat}" == "narrow" ]]; then
                pf="narrowpeak"
            else
                pf="broadpeak"
            fi
            printf "\nComparing {wildcards.sample1} with {wildcards.sample2} with .%sPeak files with MAnorm version:\n" "{params.peakformat}"
            manorm --version
            # MAnorm's Huber regressor can fail when the two peak sets have
            # too few overlapping regions (e.g. sparse low-coverage ChIP).
            # That's an analytical outcome, not a pipeline bug — catch it and
            # emit a placeholder so downstream rules still succeed.
            if ! manorm --p1 {input.peak_file1} --p2 {input.peak_file2} --r1 {input.read_file1} --r2 {input.read_file2} --n1 {wildcards.sample1} --n2 {wildcards.sample2} -o {params.output_folder} --rf "bam" --pf "$pf" {params.diffpeaks}; then
                printf "\n========================================================================\n"
                printf "WARNING: MAnorm failed for {wildcards.sample1} vs {wildcards.sample2}\n"
                printf "(likely insufficient overlapping peaks for Huber regression).\n"
                printf "This is an analytical outcome, not a pipeline error.\n"
                printf "========================================================================\n"
                mkdir -p {params.output_folder}
                printf "# MAnorm failed: insufficient overlapping peaks ({wildcards.sample1}=%d, {wildcards.sample2}=%d)\n" "$n1" "$n2" > {output.result}
            fi
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule all_chip:
    input:
        final = lambda wildcards: define_final_chip_output(wildcards.ref_genome)
    output:
        touch = f"{RESULTS_DIR}/{{env}}/chkpts/{{env}}_analysis__{{analysis_name}}__{{ref_genome}}.done"
    wildcard_constraints:
        env = "ChIP"
    localrule: True
    shell:
        """
        touch {output.touch}
        """
