CONDA_ENV_ATAC=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_chip.yaml")

def define_final_atac_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    trimmed_fastqs = config['trimmed_fastqs']
    aligned_bams = config['aligned_bams']
    motifs = config['motifs']
    motifs_allreps = config['motifs_allreps']
    map_files = []
    stat_files = []
    qc_files = []
    peak_files = []
    bigwig_files = []
    motif_files = []
    allrep_files = []

    filtered_rep_samples = samples[ (samples['env'] == 'ATAC') & (samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_rep_samples.iterrows():
        sname = row['sample_name']
        mname = row['mapped_name']
        paired = get_sample_info_from_name(sname, samples, 'paired')
        env = "ATAC"
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

    # ATAC has no Input samples to filter out
    peaktype = config["atac_callpeaks"]["peaktype"]
    for _, row in filtered_rep_samples.iterrows():
        sname = row['sample_name']
        mname = row['mapped_name']
        paired = get_sample_info_from_name(sname, samples, 'paired')
        bigwig_files.append(f"{RESULTS_DIR}/ATAC/tracks/coverage__final__{mname}.bw")
        peak_files.append(f"{RESULTS_DIR}/ATAC/peaks/peaks_atac__final__{mname}_peaks.{peaktype}Peak")

    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'ATAC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        spname = row['sample_name']
        peak_files.append(f"{RESULTS_DIR}/ATAC/peaks/selected_peaks__{spname}.bedPeak")
        reps = get_replicate_sample_ids(row['sample_name'], samples)
        if len(reps) >= 2:
            bigwig_files.append(f"{RESULTS_DIR}/ATAC/tracks/coverage__merged__{row['mapped_name']}.bw")
            stat_files.append(f"{RESULTS_DIR}/ATAC/chkpts/idr__{spname}.done")
            
    for a, b in combinations(filtered_analysis_samples.itertuples(index=False), 2):
        sample1 = a._asdict()['sample_name']
        sample2 = b._asdict()['sample_name']
        peak_files.append(f"{RESULTS_DIR}/ATAC/peaks/{sample1}_vs_{sample2}/{sample1}_vs_{sample2}_all_MAvalues.xls")

    # Motif discovery. ATAC peaks are point-source (open-chromatin summits), so
    # motif-around-summit is meaningful when peaks are called narrow. The
    # shared find_motifs_in_file rule (ChIPseq.smk) handles env=ATAC.
    if peaktype == "narrow":
        for _, row in filtered_analysis_samples.iterrows():
            spname = row['sample_name']
            motif_files.append(f"{RESULTS_DIR}/ATAC/chkpts/motifs__selected_peaks__{spname}.done")
            if len(get_replicate_sample_ids(spname, samples)) >= 2:
                motif_files.append(f"{RESULTS_DIR}/ATAC/chkpts/motifs__idr_peaks__{spname}.done")
        for _, row in filtered_rep_samples.iterrows():
            sname = row['sample_name']
            mname = row['mapped_name']
            allrep_files.append(f"{RESULTS_DIR}/ATAC/chkpts/motifs__peaks_atac__final__{mname}_peaks.done")

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


rule atac_shift_bam:
    input:
        # alignmentSieve requires an indexed input BAM; depend on the .bai
        # explicitly so Snakemake (re)builds it via the producing rule
        # (dispatch_final_bam / merging_bam_replicates / making_pseudo_replicates)
        # whenever the BAM itself changes.
        bamfile = f"{RESULTS_DIR}/ATAC/mapped/{{file_type}}__{{sample_name}}.bam",
        bambai = f"{RESULTS_DIR}/ATAC/mapped/{{file_type}}__{{sample_name}}.bam.bai"
    output:
        shifted_bam = maybe_temp(f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bam", config.get('keep_shifted_bams', False)),
        shifted_bai = maybe_temp(f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bam.bai", config.get('keep_shifted_bams', False))
    wildcard_constraints:
        file_type = "final|merged|pseudo1|pseudo2"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        file_type = lambda wildcards: wildcards.file_type
    log:
        temp(return_log_chip("ATAC","{sample_name}", "atac_shift_{file_type}", ""))
    conda: CONDA_ENV_ATAC
    shell:
        """
        {{
        printf "\nApplying Tn5 shift for {params.file_type}__{params.sample_name}\n"
        # Single-threaded alignmentSieve (stdout-safe) piped into parallel
        # samtools sort: Tn5 shift is IO-bound and fast, and alignmentSieve's
        # multithreaded mode emits records out of position order anyway. If
        # this becomes a wall-time bottleneck on very large BAMs, revisit.
        alignmentSieve --ATACshift -b {input.bamfile} -p 1 -o /dev/stdout \
            | samtools sort -@ {threads} -T "$TMPDIR/shift.sort" -o {output.shifted_bam} -
        samtools index -@ {threads} {output.shifted_bam}
        }} 2>&1 | tee -a "{log}"
        """

rule atac_bam_to_bed:
    input:
        bamfile = f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bam"
    output:
        bedfile = temp(f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bed.gz")
    wildcard_constraints:
        file_type = "final|merged|pseudo1|pseudo2"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        file_type = lambda wildcards: wildcards.file_type
    log:
        temp(return_log_chip("ATAC","{sample_name}", "atac_bed_{file_type}", ""))
    conda: CONDA_ENV_ATAC
    shell:
        """
        {{
        printf "\nConverting shifted BAM to BED for {params.file_type}__{params.sample_name}\n"
        bedtools bamtobed -i {input.bamfile} | pigz -p {threads} > {output.bedfile}
        }} 2>&1 | tee -a "{log}"
        """

rule calling_peaks_atac:
    input:
        bedfile = f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bed.gz",
        genome_stats = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/genome_stats.json"
    output:
        peakfile = f"{RESULTS_DIR}/ATAC/peaks/peaks_atac__{{file_type}}__{{sample_name}}_peaks.{{peaktype}}Peak"
    wildcard_constraints:
        file_type = "final|merged|pseudo1|pseudo2",
        peaktype = "broad|narrow"
    params:
        ipname = lambda wildcards: wildcards.sample_name,
        filetype = lambda wildcards: wildcards.file_type,
        params = config["atac_callpeaks"]['params'],
        peaktype = lambda wildcards: _checked_peaktype(
            wildcards.sample_name, "ATAC", wildcards.peaktype),
        broad_flag = lambda wildcards: "--broad" if wildcards.peaktype == "broad" else "",
        genomesize_override = lambda wildcards: config["genomes"][parse_sample_name(wildcards.sample_name)['ref_genome']].get('genomesize', '')
    log:
        temp(return_log_chip("ATAC","{sample_name}", "{file_type}__{peaktype}peak_calling", ""))
    conda: CONDA_ENV_ATAC
    shell:
        """
        {{
        if [[ -n "{params.genomesize_override}" ]]; then
            gsize="{params.genomesize_override}"
        else
            gsize=$(python3 -c "import json; print(json.load(open('{input.genome_stats}'))['effective_size'])")
        fi
        printf "\nCalling {params.peaktype} peaks for ATAC-seq {params.ipname} using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.bedfile} -f BED \
            -g $gsize {params.params} {params.broad_flag} \
            -n peaks_atac__{params.filetype}__{params.ipname} \
            --outdir {config[output_dir]}/ATAC/peaks/
        }} 2>&1 | tee -a "{log}"
        """

rule make_coverage_atac:
    input:
        bamfile = f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bam",
        bai = f"{RESULTS_DIR}/ATAC/mapped/shifted_{{file_type}}__{{sample_name}}.bam.bai"
    output:
        bigwig = f"{RESULTS_DIR}/ATAC/tracks/coverage__{{file_type}}__{{sample_name}}.bw"
    wildcard_constraints:
        file_type = "final|merged"
    params:
        binsize = config['atac_tracks']['binsize'],
        params = config['atac_tracks']['params']
    log:
        temp(return_log_chip("ATAC","{sample_name}", "making_bigwig_{file_type}", ""))
    conda: CONDA_ENV_ATAC
    shell:
        """
        {{
        printf "\nMaking coverage bigwig for ATAC-seq\n"
        bamCoverage -b {input.bamfile} -o {output.bigwig} -bs {params.binsize} -p {threads} {params.params}
        }} 2>&1 | tee -a "{log}"
        """

rule all_atac:
    input:
        final = lambda wildcards: define_final_atac_output(wildcards.ref_genome)
    output:
        touch = f"{RESULTS_DIR}/ATAC/chkpts/ATAC_analysis__{{analysis_name}}__{{ref_genome}}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """
