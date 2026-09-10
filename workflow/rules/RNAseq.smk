CONDA_ENV_RNA=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_rnaseq.yaml")

def return_log_rna(sample_name, step, paired):
    return os.path.join(REPO_FOLDER, RESULTS_DIR,"RNA","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def define_RNA_input_for_degs(ref_genome):
    filtered_samples = samples[(samples['Assay'] == 'RNAseq') & (samples['Genome'] == ref_genome)]
    # Per-replicate counts come out of alignment, so they carry the genome.
    return [f"{RESULTS_DIR}/RNA/DEG/counts__{mname}.tab" for mname in filtered_samples['mapped_name']]

def define_rnaseq_target_file(wildcards):
    tname = config['rnaseq_target_file_label']
    if wildcards.target_name == tname:
        return config['rnaseq_target_file']
    elif wildcards.target_name == "unique_DEGs":
        return f"{RESULTS_DIR}/RNA/DEG/unique_DEGs__{wildcards.analysis_name}__{wildcards.ref_genome}.txt"
    else:
        raise ValueError(
            f"{wildcards.target_name} does not match possible files."
            "It can be 'unique_DEGs', or the value of "
            "'rnaseq_target_file_name' in the options file"
        )

def define_rnaseq_background_file(wildcards):
    tname = config['rnaseq_target_file_label']
    bgfile = config['rnaseq_background_file']
    if wildcards.target_name == "unique_DEGs":
        return f"{RESULTS_DIR}/RNA/DEG/counts__{wildcards.analysis_name}__{wildcards.ref_genome}.txt"
    elif wildcards.target_name == tname and os.path.exists(bgfile):
        return config['rnaseq_background_file']
    else:
        return f"{RESULTS_DIR}/combined/bedfiles/{wildcards.ref_genome}__all_genes.bed"

def get_go_database(ref_genome):
    species=config["genomes"][ref_genome]['species']
    genus=config["genomes"][ref_genome]['genus']
    return f"{GENOMES_DIR}/{ref_genome}/GO/org.{genus[0]}{species}.eg.db"

def assign_rna_input(wildcards):
    """Find the RNA control bam for a RAMPAGE sample.

    Uses the Control column when available, otherwise falls back to
    matching by the full factor-level combination among RNAseq samples.
    """
    sname = wildcards.sample_name
    # Check if there is an explicit control in the Control column
    ctrl = get_control_sample_id(sname, samples)
    if ctrl:
        ctrl_parsed = parse_sample_name(ctrl)
        if ctrl_parsed['replicate'] == 'merged':
            return f"merged__{ctrl}"
        else:
            return f"final__{ctrl}"
    # Fallback: find matching RNAseq sample by the full factor-level combination
    parsed = parse_sample_name(sname)
    same_name = samples[
        (samples['Assay'] == 'RNAseq') &
        (samples['Genome'] == parsed['ref_genome']) &
        (samples['levels_label'] == parsed['levels_label'])
    ]
    if len(same_name) == 1:
        return f"final__{same_name.iloc[0]['sample_name']}"
    elif len(same_name) >= 2:
        # Multiple replicates: find the merged analysis_name for this group
        for _, r in analysis_samples.iterrows():
            if (r['Assay'] == 'RNAseq' and r['levels_label'] == parsed['levels_label'] and
                r['Genome'] == parsed['ref_genome']):
                return f"merged__{r['sample_name']}"
    raise ValueError(f"\nSample '{sname}' does not have corresponding RNA control for calling TSS")

def define_final_rna_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    analysis_name = config['analysis_name']
    go_analysis = config['GO']
    trimmed_fastqs = config['trimmed_fastqs']
    map_files = []
    bigwig_files = []
    qc_files = []
    deg_files = []
    tss_files = []
    filtered_rep_samples = samples[(samples['env'] == 'RNA') & (samples['Genome'] == ref_genome)].copy()

    for _, row in filtered_rep_samples.iterrows():
        sname = row['sample_name']
        mname = row['mapped_name']
        paired = row['paired']
        if paired == "PE":
            qc_files.append(f"{RESULTS_DIR}/RNA/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
            qc_files.append(f"{RESULTS_DIR}/RNA/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
            map_files.append(f"{RESULTS_DIR}/RNA/logs/process_rna_pe_sample__{mname}.log")
            if not trimmed_fastqs:
                qc_files.append(f"{RESULTS_DIR}/RNA/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
                qc_files.append(f"{RESULTS_DIR}/RNA/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
        elif paired == "SE":
            qc_files.append(f"{RESULTS_DIR}/RNA/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
            map_files.append(f"{RESULTS_DIR}/RNA/logs/process_rna_se_sample__{mname}.log")
            if not trimmed_fastqs:
                qc_files.append(f"{RESULTS_DIR}/RNA/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file

        strand = config['rna_tracks'][row['Assay']]['strandedness']
        if strand == "unstranded":
            bigwig_files.append(f"{RESULTS_DIR}/RNA/tracks/{mname}__unstranded.bw")
        else:
            bigwig_files.append(f"{RESULTS_DIR}/RNA/tracks/{mname}__plus.bw")
            bigwig_files.append(f"{RESULTS_DIR}/RNA/tracks/{mname}__minus.bw")

    filtered_analysis_samples = analysis_samples[(analysis_samples['env'] == 'RNA') & (analysis_samples['Genome'] == ref_genome)].copy()
    for _, row in filtered_analysis_samples.iterrows():
        strand = config['rna_tracks'][row['Assay']]['strandedness']
        akey = build_analysis_key(row)
        if len(analysis_to_replicates.get(akey, [])) >= 2:
            aname = row['sample_name']
            if strand == "unstranded":
                bigwig_files.append(f"{RESULTS_DIR}/RNA/tracks/{aname}__unstranded.bw")
            else:
                bigwig_files.append(f"{RESULTS_DIR}/RNA/tracks/{aname}__plus.bw")
                bigwig_files.append(f"{RESULTS_DIR}/RNA/tracks/{aname}__minus.bw")

    filtered_samples2 = samples[(samples['Assay'] == 'RNAseq') & (samples['Genome'] == ref_genome)].copy()
    # Group by the full factor-level combination (all factors), not just the
    # first two — a 3+-factor design must not collapse distinct groups.
    filtered_samples2['Sample'] = filtered_samples2['levels_label']
    if len(filtered_samples2['Sample'].drop_duplicates()) >= 2:
        deg_files.append(f"{RESULTS_DIR}/RNA/chkpts/calling_DEGs__{analysis_name}__{ref_genome}.done")
        deg_files.append(f"{RESULTS_DIR}/RNA/DEG/genes_rpkm__{analysis_name}__{ref_genome}.txt")
        deg_files.append(f"{RESULTS_DIR}/RNA/plots/plot_expression__{analysis_name}__{ref_genome}__unique_DEGs.pdf")

        if go_analysis:
            deg_files.append(f"{RESULTS_DIR}/RNA/GO/TopGO__{analysis_name}__{ref_genome}__unique_DEGs.done")

    elif len(filtered_samples2['Sample'].drop_duplicates()) == 1:
        deg_files.append(f"{RESULTS_DIR}/RNA/DEG/genes_rpkm__{analysis_name}__{ref_genome}.txt")

    filtered_samples3 = samples[(samples['Assay'] == 'RAMPAGE') & (samples['Genome'] == ref_genome)].copy()
    filtered_samples3['Sample'] = filtered_samples3['levels_label']
    valid_samples = set(filtered_samples2['Sample'])
    for _, row in filtered_samples3.iterrows():
        if row['Sample'] in valid_samples:
            tss_files.append(f"{RESULTS_DIR}/RNA/TSS/TSS__final__{row['mapped_name']}_peaks.narrowPeak")

    filtered_analysis_samples2 = analysis_samples[(analysis_samples['Assay'] == 'RAMPAGE') & (analysis_samples['Genome'] == ref_genome)].copy()
    filtered_analysis_samples2['Sample'] = filtered_analysis_samples2['levels_label']
    for _, row in filtered_analysis_samples2.iterrows():
        if row['Sample'] in valid_samples:
            tss_files.append(f"{RESULTS_DIR}/RNA/TSS/TSS__merged__{row['sample_name']}_peaks.narrowPeak")

    results = map_files + bigwig_files

    if qc_option == "all":
        results += qc_files

    if analysis:
        results += deg_files + tss_files

    return results

rule make_STAR_indices:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa",
        gtf = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gtf",
        genome_stats = f"{GENOMES_DIR}/{{ref_genome}}/genome_stats.json"
    output:
        indices = directory(f"{GENOMES_DIR}/{{ref_genome}}/STAR_index")
    params:
        star_index_override = lambda wildcards: config["genomes"][wildcards.ref_genome].get('star_index', ''),
        # Explicit STAR scratch — STAR insists the directory not exist
        # at startup, so we rm -rf before the call. Routing into the
        # genomes/ tree keeps temp work off any tmpfs/$TMPDIR.
        outtmpdir = f"{GENOMES_DIR}/{{ref_genome}}/STAR_index_tmp"
    log:
        temp(os.path.join(REPO_FOLDER, RESULTS_DIR,"logs","STAR_index_{ref_genome}.log"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        if [[ -n "{params.star_index_override}" ]]; then
            star_idx="{params.star_index_override}"
        else
            nbases=$(python3 -c "import json; print(json.load(open('{input.genome_stats}'))['star_sa_index_nbases'])")
            star_idx="--genomeSAindexNbases $nbases"
        fi
        printf "\nBuilding STAR index directory for {wildcards.ref_genome} ($star_idx)\n"
        mkdir "{output.indices}"
        rm -rf "{params.outtmpdir}"
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir "{output.indices}" --genomeFastaFiles "{input.fasta}" --sjdbGTFfile "{input.gtf}" --outTmpDir "{params.outtmpdir}" $star_idx
        rm -rf "{params.outtmpdir}"
        }} 2>&1 | tee -a "{log}"
        """

rule STAR_map_pe:
    input:
        fastq1 = lambda wildcards: f"{RESULTS_DIR}/RNA/fastq/trim__{base_sample(wildcards.sample_name)}__R1.fastq.gz",
        fastq2 = lambda wildcards: f"{RESULTS_DIR}/RNA/fastq/trim__{base_sample(wildcards.sample_name)}__R2.fastq.gz",
        indices = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        bamfile = temp(f"{RESULTS_DIR}/RNA/mapped/star_pe__{{sample_name}}_Aligned.out.bam"),
        count_file = temp(f"{RESULTS_DIR}/RNA/mapped/star_pe__{{sample_name}}_ReadsPerGene.out.tab"),
        metrics_map = f"{RESULTS_DIR}/RNA/reports/star_pe__{{sample_name}}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        file_order = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['file_order'],
        prefix = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/star_pe__{wildcards.sample_name}_",
        # STAR's intermediate scratch. STAR insists the directory not
        # exist at startup, so we rm -rf before the call. Routing to
        # results/ keeps temp work off any tmpfs/$TMPDIR.
        outtmpdir = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/star_pe_tmp__{wildcards.sample_name}"
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "PE"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        if [[ "{params.file_order}" == "rampage" ]]; then
            printf "Input file order for RAMPAGE (R2 R1)\n"
            input='"{input.fastq2}" "{input.fastq1}"'
        else
            printf "Input file order for RNAseq (R1 R2)\n"
            input='"{input.fastq1}" "{input.fastq2}"'
        fi
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
        STAR --version
        rm -rf "{params.outtmpdir}"
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn ${{input}} --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outTmpDir "{params.outtmpdir}" --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 20 --quantMode GeneCounts
        rm -rf "{params.outtmpdir}"
        mv "{config[output_dir]}/RNA/mapped/star_pe__{params.sample_name}_Log.final.out" "{output.metrics_map}"
        rm -f {params.prefix}Log*
        }} 2>&1 | tee -a "{log}"
        """

rule STAR_map_se:
    input:
        fastq0 = lambda wildcards: f"{RESULTS_DIR}/RNA/fastq/trim__{base_sample(wildcards.sample_name)}__R0.fastq.gz",
        indices = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/STAR_index"
    output:
        bamfile = temp(f"{RESULTS_DIR}/RNA/mapped/star_se__{{sample_name}}_Aligned.out.bam"),
        count_file = temp(f"{RESULTS_DIR}/RNA/mapped/star_se__{{sample_name}}_ReadsPerGene.out.tab"),
        metrics_map = f"{RESULTS_DIR}/RNA/reports/star_se__{{sample_name}}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/star_se__{wildcards.sample_name}_",
        outtmpdir = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/star_se_tmp__{wildcards.sample_name}"
    log:
        temp(return_log_rna("{sample_name}", "mappingSTAR", "SE"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        printf "\nMapping {params.sample_name} to {params.ref_genome} with STAR version:\n"
        STAR --version
        rm -rf "{params.outtmpdir}"
        STAR --runMode alignReads --genomeDir "{input.indices}" --readFilesIn "{input.fastq0}" --readFilesCommand zcat --runThreadN {threads} --genomeLoad NoSharedMemory --outMultimapperOrder Random --outFileNamePrefix "{params.prefix}" --outTmpDir "{params.outtmpdir}" --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --quantMode GeneCounts
        rm -rf "{params.outtmpdir}"
        mv "{config[output_dir]}/RNA/mapped/star_se__{params.sample_name}_Log.final.out" "{output.metrics_map}"
        rm -f {params.prefix}Log*
        }} 2>&1 | tee -a "{log}"
        """

rule filter_rna_pe:
    input:
        bamfile = f"{RESULTS_DIR}/RNA/mapped/star_pe__{{sample_name}}_Aligned.out.bam"
    output:
        sorted_file=temp(f"{RESULTS_DIR}/RNA/mapped/star_pe__{{sample_name}}_Processed.sorted.out.bam"),
        metrics_flag = f"{RESULTS_DIR}/RNA/reports/flagstat_pe__{{sample_name}}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        # Deduplicate RAMPAGE always (5'-tag TSS libraries), mRNA-seq only when
        # the user opts in (no UMI handling yet — see rna_deduplicate).
        deduplicate = lambda wildcards: (parse_sample_name(wildcards.sample_name)['data_type'] == 'RAMPAGE') or bool(config.get('rna_deduplicate', False)),
        # ENCODE RAMPAGE folds 15 bases of mate 2 (the 5' tag region) into the
        # duplicate signature; plain mRNA-seq dedup does not.
        dedup_extra = lambda wildcards: "--bamRemoveDuplicatesMate2basesN 15" if parse_sample_name(wildcards.sample_name)['data_type'] == 'RAMPAGE' else "",
        mrkdup = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/star_pe__{wildcards.sample_name}_Processed.out.bam"
    log:
        temp(return_log_rna("{sample_name}", "filteringRNA", "PE"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        ## Read deduplication is OFF by default for mRNA-seq: without UMIs, PCR/optical
        ## duplicates are indistinguishable from genuine high-expression reads, so
        ## deduplicating biases quantification. RAMPAGE is deduplicated following ENCODE
        ## (UniqueIdentical + Mate2basesN 15). Errors can happen because of limitBAMsortRAM,
        ## which seem to happen when bam files are sorted by coordinates (now removed from
        ## mapping step).
        if [[ "{params.deduplicate}" == "True" ]]; then
            printf "\nMarking duplicates\n"
            STAR --runMode inputAlignmentsFromBAM --inputBAMfile "{input.bamfile}" --bamRemoveDuplicatesType UniqueIdentical {params.dedup_extra} --outFileNamePrefix "{config[output_dir]}/RNA/mapped/star_pe__{params.sample_name}_"
            to_sort="{params.mrkdup}"
        else
            printf "\nDeduplication disabled; sorting aligned reads as-is\n"
            to_sort="{input.bamfile}"
        fi
        printf "\nSorting bam file\n"
        samtools sort -@ {threads} -T "{output.sorted_file}.sort" "${{to_sort}}" -o "{output.sorted_file}"
        [[ "{params.deduplicate}" == "True" ]] && rm -f "{params.mrkdup}"
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "{output.sorted_file}"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.sorted_file}" > "{output.metrics_flag}"
        }} 2>&1 | tee -a "{log}"
        """

rule filter_rna_se:
    input:
        bamfile = f"{RESULTS_DIR}/RNA/mapped/star_se__{{sample_name}}_Aligned.out.bam"
    output:
        sorted_file=temp(f"{RESULTS_DIR}/RNA/mapped/star_se__{{sample_name}}_Aligned.sorted.out.bam"),
        metrics_flag = f"{RESULTS_DIR}/RNA/reports/flagstat_se__{{sample_name}}.txt"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        # SE deduplication only on explicit opt-in (no UMI handling yet). ENCODE
        # RAMPAGE dedup is paired-end (needs mate 2), so RAMPAGE SE — an unusual
        # case — keeps its prior no-dedup behavior.
        deduplicate = lambda wildcards: bool(config.get('rna_deduplicate', False)),
        mrkdup = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/star_se__{wildcards.sample_name}_Processed.out.bam"
    log:
        temp(return_log_rna("{sample_name}", "filteringRNA", "SE"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        if [[ "{params.deduplicate}" == "True" ]]; then
            printf "\nMarking duplicates\n"
            STAR --runMode inputAlignmentsFromBAM --inputBAMfile "{input.bamfile}" --bamRemoveDuplicatesType UniqueIdentical --outFileNamePrefix "{config[output_dir]}/RNA/mapped/star_se__{params.sample_name}_"
            to_sort="{params.mrkdup}"
        else
            printf "\nDeduplication disabled; sorting aligned reads as-is\n"
            to_sort="{input.bamfile}"
        fi
        #### Sorting bam file
        printf "\nSorting bam file\n"
        samtools sort -@ {threads} -T "{output.sorted_file}.sort" "${{to_sort}}" -o "{output.sorted_file}"
        [[ "{params.deduplicate}" == "True" ]] && rm -f "{params.mrkdup}"
        #### Indexing bam file
        printf "\nIndexing bam file\n"
        samtools index -@ {threads} "{output.sorted_file}"
        #### Getting stats from bam file
        printf "\nGetting some stats\n"
        samtools flagstat -@ {threads} "{output.sorted_file}" > "{output.metrics_flag}"
        }} 2>&1 | tee -a "{log}"
        """

rule make_rna_stats_pe:
    input:
        metrics_trim = lambda wildcards: f"{RESULTS_DIR}/RNA/reports/trim_pe__{base_sample(wildcards.sample_name)}.json",
        metrics_map = f"{RESULTS_DIR}/RNA/reports/star_pe__{{sample_name}}.txt",
        logs = lambda wildcards: [ return_log_rna((base_sample(wildcards.sample_name) if step in ("downloading", "trimming") else wildcards.sample_name), step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingSTAR", "filteringRNA"] ]
    output:
        stat_file = f"{RESULTS_DIR}/RNA/reports/summary_RNA_PE_mapping_stats_{{sample_name}}.txt",
        log = f"{RESULTS_DIR}/RNA/logs/process_rna_pe_sample__{{sample_name}}.log"
    params:
        group = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        levels = lambda wildcards: parse_sample_name(wildcards.sample_name)['levels_label'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        trimmed_fastq = config['trimmed_fastqs']
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(python3 -c "import json; print(json.load(open('{input.metrics_trim}'))['summary']['before_filtering']['total_reads'] // 2)")
        else
            tot=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        fi
        filt=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Number of reads mapped to multiple loci" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Uniquely mapped reads number" "{input.metrics_map}" | awk '{{print $NF}}')
        allmap=$((multi+single))
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule make_rna_stats_se:
    input:
        metrics_trim = lambda wildcards: f"{RESULTS_DIR}/RNA/reports/trim_se__{base_sample(wildcards.sample_name)}.json",
        metrics_map = f"{RESULTS_DIR}/RNA/reports/star_se__{{sample_name}}.txt",
        logs = lambda wildcards: [ return_log_rna((base_sample(wildcards.sample_name) if step in ("downloading", "trimming") else wildcards.sample_name), step, get_sample_info_from_name(wildcards.sample_name, samples, 'paired')) for step in ["downloading", "trimming", "mappingSTAR", "filteringRNA"] ]
    output:
        stat_file = f"{RESULTS_DIR}/RNA/reports/summary_RNA_SE_mapping_stats_{{sample_name}}.txt",
        log = f"{RESULTS_DIR}/RNA/logs/process_rna_se_sample__{{sample_name}}.log"
    params:
        group = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        levels = lambda wildcards: parse_sample_name(wildcards.sample_name)['levels_label'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        trimmed_fastq = config['trimmed_fastqs']
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(python3 -c "import json; print(json.load(open('{input.metrics_trim}'))['summary']['before_filtering']['total_reads'])")
        else
            tot=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        fi
        filt=$(grep "Number of input reads" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Number of reads mapped to multiple loci" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Uniquely mapped reads number" "{input.metrics_map}" | awk '{{print $NF}}')
        allmap=$((multi+single))
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\n" > {output.stat_file}
        awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v m={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v a=${{tot}} -v b=${{filt}} -v c=${{allmap}} -v d=${{single}} 'BEGIN {{print l,t,m,r,g,a,b" ("b/a*100"%)",c" ("c/a*100"%)",d" ("d/a*100"%)"}}' >> "{output.stat_file}"
        cat {input.logs} > "{output.log}"
        rm -f {input.logs}
        """

rule pe_or_se_rna_dispatch:
    input:
        bamfile = lambda wildcards: assign_mapping_paired(wildcards, "filter_rna", "sorted_file"),
        countfile = lambda wildcards: assign_mapping_paired(wildcards, "STAR_map", "count_file")
    output:
        bam_file = maybe_temp(f"{RESULTS_DIR}/RNA/mapped/final__{{sample_name}}.bam", config.get('keep_final_bams', True)),
        count_file = f"{RESULTS_DIR}/RNA/DEG/counts__{{sample_name}}.tab",
        touch = f"{RESULTS_DIR}/RNA/chkpts/map_RNA__{{sample_name}}.done"
    localrule: True
    shell:
        """
        mv {input.bamfile} {output.bam_file}
        mv {input.bamfile}.bai {output.bam_file}.bai
        mv {input.countfile} {output.count_file}
        touch {output.touch}
        """

rule merging_rna_replicates:
    input:
        bamfiles = lambda wildcards: [
            f"{RESULTS_DIR}/RNA/mapped/final__{sid}.bam"
            for sid in get_replicate_sample_ids(wildcards.sample_name, samples)
        ]
    output:
        mergefile = maybe_temp(f"{RESULTS_DIR}/RNA/mapped/merged__{{sample_name}}.bam", config.get('keep_merged_bams', False))
    params:
        sname = lambda wildcards: wildcards.sample_name
    log:
        temp(return_log_rna("{sample_name}", "merging_rna_reps", ""))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
        samtools merge -u -@ {threads} - {input.bamfiles} | samtools sort -@ {threads} -T {output.mergefile}.sort -o {output.mergefile}
        samtools index -@ {threads} {output.mergefile}
        }} 2>&1 | tee -a "{log}"
        """

rule make_rna_stranded_bigwigs:
    input:
        bamfile = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'final'}__{wildcards.sample_name}.bam",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        bw_plus = f"{RESULTS_DIR}/RNA/tracks/{{sample_name}}__plus.bw",
        bw_minus = f"{RESULTS_DIR}/RNA/tracks/{{sample_name}}__minus.bw"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        param_bg = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['param_bg'],
        strandedness = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['strandedness'],
        multimap = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['multimap']
    log:
        temp(return_log_rna("{sample_name}", "making_bigiwig", ""))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --runThreadN {threads} --inputBAMfile "{input.bamfile}" --outWigStrand Stranded {params.param_bg} --outFileNamePrefix "{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_"
        ### Converting to bigwig files
        printf "\nConverting bedGraphs to bigWigs\n"
        if [[ "{params.multimap}" == "multiple" ]]; then
            bed1="{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str1.out.bg"
            bed2="{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str2.out.bg"
        elif [[ "{params.multimap}" == "unique" ]]; then
            bed1="{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_Signal.Unique.str1.out.bg"
            bed2="{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_Signal.Unique.str2.out.bg"
        fi
        bedSort ${{bed1}} "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg"
        bedSort ${{bed2}} "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str2.out.bg"
        if [[ "{params.strandedness}" == "forward" ]]; then
            bedGraphToBigWig "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
            bedGraphToBigWig "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str2.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
        elif [[ "{params.strandedness}" == "reverse" ]]; then
            bedGraphToBigWig "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg" "{input.chrom_sizes}" "{output.bw_minus}"
            bedGraphToBigWig "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str2.out.bg" "{input.chrom_sizes}" "{output.bw_plus}"
        fi
        rm -f {config[output_dir]}/RNA/tracks/bg_{params.sample_name}_* {config[output_dir]}/RNA/tracks/{params.sample_name}_Signal*
        }} 2>&1 | tee -a "{log}"
        """

rule make_rna_unstranded_bigwigs:
    input:
        bamfile = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/{'merged' if parse_sample_name(wildcards.sample_name)['replicate'] == 'merged' else 'final'}__{wildcards.sample_name}.bam",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        bw_unstranded = f"{RESULTS_DIR}/RNA/tracks/{{sample_name}}__unstranded.bw"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        param_bg = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['param_bg'],
        strandedness = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['strandedness'],
        multimap = lambda wildcards: config['rna_tracks'][parse_sample_name(wildcards.sample_name)['sample_type']]['multimap']
    log:
        temp(return_log_rna("{sample_name}", "making_bigiwig", ""))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        ### Making BedGraph files
        printf "\nMaking bedGraph files\n"
        STAR --runMode inputAlignmentsFromBAM --runThreadN {threads} --inputBAMfile "{input.bamfile}" --outWigStrand Unstranded {params.param_bg} --outFileNamePrefix "{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_"
        printf "\nConverting bedGraphs to bigWigs\n"
        if [[ "{params.multimap}" == "multiple" ]]; then
            bed1="{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_Signal.UniqueMultiple.str1.out.bg"
        elif [[ "{params.multimap}" == "unique" ]]; then
            bed1="{config[output_dir]}/RNA/tracks/bg_{params.sample_name}_Signal.Unique.str1.out.bg"
        fi
        bedSort ${{bed1}} "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg"
        bedGraphToBigWig "{config[output_dir]}/RNA/tracks/{params.sample_name}_Signal.sorted.str1.out.bg" "{input.chrom_sizes}" "{output.bw_unstranded}"
        rm -f {config[output_dir]}/RNA/tracks/bg_{params.sample_name}_* {config[output_dir]}/RNA/tracks/{params.sample_name}_Signal*
        }} 2>&1 | tee -a "{log}"
        """

rule prep_files_for_DEGs:
    input:
        lambda wildcards: define_RNA_input_for_degs(wildcards.ref_genome)
    output:
        rna_samples = temp(f"{RESULTS_DIR}/RNA/DEG/samples__{{analysis_name}}__{{ref_genome}}.txt"),
        rna_counts = temp(f"{RESULTS_DIR}/RNA/DEG/counts__{{analysis_name}}__{{ref_genome}}.txt")
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        strand = config['rna_tracks']['RNAseq']['strandedness']
    log:
        temp(return_log_rna("{ref_genome}", "prep_for_DEGs", "{analysis_name}"))
    run:
        filtered_samples = samples[(samples['Assay'] == 'RNAseq') & (samples['Genome'] == params.ref_genome)].copy()
        # Group by the full factor-level combination (see define_rnaseq_outputs).
        filtered_samples['Sample'] = filtered_samples['levels_label']
        filtered_samples['Replicate'] = filtered_samples['Sample'] + "__" + filtered_samples['replicate'].astype(str)

        RNA_samples = filtered_samples[['Replicate','Sample']].drop_duplicates()
        RNA_samples = RNA_samples.sort_values(by=['Sample', 'Replicate'],ascending=[True, True]).reset_index(drop=True)
        RNA_samples['Color'] = pd.factorize(RNA_samples['Sample'])[0] + 1

        RNA_samples.to_csv(output.rna_samples, sep="\t", index=False)

        RNA_counts = None
        # mapped_name, not sample_name: the per-replicate counts are a
        # post-alignment product, so they carry the genome.
        replicates = filtered_samples[['mapped_name', 'Replicate']].drop_duplicates()
        for mname, rep in replicates.values:
            file_path = f"{RESULTS_DIR}/RNA/DEG/counts__{mname}.tab"
            if params.strand == "reverse":
                temp = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 3])
            elif params.strand == "forward":
                temp = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 2])
            elif params.strand == "unstranded":
                temp = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1])
            else:
                print("Unknown strandedness option, defaulting to unstranded")
                temp = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1])

            temp.columns = ['GID', rep]

            if RNA_counts is None:
                RNA_counts = temp
            else:
                RNA_counts = pd.merge(RNA_counts, temp, on='GID', how='outer')

        replicate_order = RNA_samples['Replicate'].tolist()
        column_order = ['GID'] + replicate_order
        RNA_counts = RNA_counts[column_order]
        RNA_counts.to_csv(output.rna_counts, sep="\t", index=False)

rule call_all_DEGs:
    input:
        samples = f"{RESULTS_DIR}/RNA/DEG/samples__{{analysis_name}}__{{ref_genome}}.txt",
        counts = f"{RESULTS_DIR}/RNA/DEG/counts__{{analysis_name}}__{{ref_genome}}.txt",
        region_file = f"{RESULTS_DIR}/combined/bedfiles/{{ref_genome}}__all_genes.bed"
    output:
        rdata = temp(f"{RESULTS_DIR}/RNA/DEG/ReadyToPlot__{{analysis_name}}__{{ref_genome}}.RData"),
        unique_degs = f"{RESULTS_DIR}/RNA/DEG/unique_DEGs__{{analysis_name}}__{{ref_genome}}.txt",
        mds_plot = f"{RESULTS_DIR}/combined/plots/MDS_RNAseq_{{analysis_name}}_{{ref_genome}}_d12.pdf",
        touch = f"{RESULTS_DIR}/RNA/chkpts/calling_DEGs__{{analysis_name}}__{{ref_genome}}.done"
    params:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_call_DEGs.R"),
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_rna("{ref_genome}", "call_DEGs", "{analysis_name}"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        printf "running edgeR for all samples in {params.ref_genome}\n"
        Rscript "{params.script}" "{input.counts}" "{input.samples}" "{params.analysis_name}" "{params.ref_genome}" "{input.region_file}" "{config[output_dir]}"
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule gather_gene_expression_rpkm:
    input:
        samples = f"{RESULTS_DIR}/RNA/DEG/samples__{{analysis_name}}__{{ref_genome}}.txt",
        counts = f"{RESULTS_DIR}/RNA/DEG/counts__{{analysis_name}}__{{ref_genome}}.txt",
        region_file = f"{RESULTS_DIR}/combined/bedfiles/{{ref_genome}}__all_genes.bed"
    output:
        rpkm = f"{RESULTS_DIR}/RNA/DEG/genes_rpkm__{{analysis_name}}__{{ref_genome}}.txt"
    params:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_gene_expression_rpkm.R"),
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_rna("{ref_genome}", "gene_expression", "{analysis_name}"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        printf "Gathering gene expression levels for samples from {params.analysis_name} mapping to {params.ref_genome}\n"
        Rscript "{params.script}" "{input.counts}" "{input.samples}" "{params.analysis_name}" "{params.ref_genome}" "{input.region_file}" "{config[output_dir]}"
        }} 2>&1 | tee -a "{log}"
        """

rule plot_expression_levels:
    input:
        rdata = f"{RESULTS_DIR}/RNA/DEG/ReadyToPlot__{{analysis_name}}__{{ref_genome}}.RData",
        target_file = lambda wildcards: define_rnaseq_target_file(wildcards)
    output:
        plot = f"{RESULTS_DIR}/RNA/plots/plot_expression__{{analysis_name}}__{{ref_genome}}__{{target_name}}.pdf"
    params:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_plot_expression_level.R"),
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name
    log:
        temp(return_log_rna("{ref_genome}", "plot_expression_{target_name}", "{analysis_name}"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        printf "running plot expression levels for {input.target_file} (from {params.analysis_name} and {params.ref_genome})\n"
        Rscript "{params.script}" "{params.analysis_name}" "{params.ref_genome}" "{input.target_file}" "{params.target_name}" "{config[output_dir]}"
        }} 2>&1 | tee -a "{log}"
        """

rule create_GO_database:
    input:
        taxid_file = f"{GENOMES_DIR}/{{ref_genome}}/taxid.json"
    output:
        godb = directory(f"{GENOMES_DIR}/{{ref_genome}}/GO/{{dbname}}"),
        gaf = f"{GENOMES_DIR}/{{ref_genome}}/GO/{{dbname}}_gaf_file.tab",
        geneinfo = f"{GENOMES_DIR}/{{ref_genome}}/GO/{{dbname}}_gene_info.tab"
    params:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_build_GO_database.R"),
        ref_genome = lambda wildcards: wildcards.ref_genome,
        species = lambda wildcards: config["genomes"][wildcards.ref_genome]['species'],
        genus = lambda wildcards: config["genomes"][wildcards.ref_genome]['genus'],
        gaffile = lambda wildcards: resolve_repo_data(config["genomes"][wildcards.ref_genome]['gaf_file']),
        geneinfofile = lambda wildcards: resolve_repo_data(config["genomes"][wildcards.ref_genome]['gene_info_file'])
    log:
        temp(return_log_rna("{ref_genome}", "build_GO", "{dbname}"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        rm -rf {output.godb}
        # gaf_file / gene_info_file may be HTTP(S) URLs (validation permits them,
        # like the other genome-config fields). Fetch to a temp file first if so,
        # then the content-based gzip-or-plain handling below proceeds on the
        # local copy. curl flags mirror the genome-download rules in
        # environment_setup.smk.
        gaf_src="{params.gaffile}"
        geneinfo_src="{params.geneinfofile}"
        gaf_tmp="" ; geneinfo_tmp=""
        trap 'rm -f "$gaf_tmp" "$geneinfo_tmp" 2>/dev/null || true' EXIT
        if [[ "$gaf_src" == http://* || "$gaf_src" == https://* ]]; then
            printf "\nDownloading GO gaf_file from URL: $gaf_src\n"
            gaf_tmp=$(mktemp --suffix=.gaf.dl)
            curl --fail --silent --show-error --location --max-redirs 5 \
                 --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 1800 \
                 --proto '=https,http' -o "$gaf_tmp" "$gaf_src"
            gaf_src="$gaf_tmp"
        fi
        if [[ "$geneinfo_src" == http://* || "$geneinfo_src" == https://* ]]; then
            printf "\nDownloading GO gene_info_file from URL: $geneinfo_src\n"
            geneinfo_tmp=$(mktemp --suffix=.geneinfo.dl)
            curl --fail --silent --show-error --location --max-redirs 5 \
                 --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 1800 \
                 --proto '=https,http' -o "$geneinfo_tmp" "$geneinfo_src"
            geneinfo_src="$geneinfo_tmp"
        fi
        if file "$gaf_src" | grep -q 'gzip compressed'; then
            gunzip -c "$gaf_src" > {output.gaf}
        else
            cp "$gaf_src" {output.gaf}
        fi
        if file "$geneinfo_src" | grep -q 'gzip compressed'; then
            gunzip -c "$geneinfo_src" > {output.geneinfo}
        else
            cp "$geneinfo_src" {output.geneinfo}
        fi
        ncbi_taxid=$(python3 -c "import json; print(json.load(open('{input.taxid_file}'))['ncbi_taxid'])")
        printf "Creating GO database for {params.ref_genome} (TaxId: $ncbi_taxid)\n"
        Rscript "{params.script}" "{output.gaf}" "{output.geneinfo}" "{params.ref_genome}" "{params.genus}" "{params.species}" "$ncbi_taxid"
        }} 2>&1 | tee -a "{log}"
        """

rule perform_GO_on_target_file:
    input:
        godb = lambda wildcards: directory(get_go_database(wildcards.ref_genome)),
        gaf = lambda wildcards: f"{get_go_database(wildcards.ref_genome)}_gaf_file.tab",
        geneinfo = lambda wildcards: f"{get_go_database(wildcards.ref_genome)}_gene_info.tab",
        target_file = lambda wildcards: define_rnaseq_target_file(wildcards),
        background_file = lambda wildcards: define_rnaseq_background_file(wildcards)
    output:
        touch = f"{RESULTS_DIR}/RNA/GO/TopGO__{{analysis_name}}__{{ref_genome}}__{{target_name}}.done"
    params:
        script = os.path.join(REPO_FOLDER,"workflow","scripts","R_GO_analysis.R"),
        dbname = lambda wildcards: os.path.basename(get_go_database(wildcards.ref_genome)),
        analysis_name = config['analysis_name'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        target_name = lambda wildcards: wildcards.target_name
    log:
        temp(return_log_rna("{ref_genome}", "GO_{target_name}", "{analysis_name}"))
    conda: CONDA_ENV_RNA
    shell:
        """
        {{
        printf "running GO analysis for {input.target_file} (from {params.analysis_name} and {params.ref_genome})\n"
        Rscript "{params.script}" "{params.dbname}" "{params.analysis_name}" "{params.ref_genome}" "{input.target_file}" "{input.background_file}" "{params.target_name}" "{config[output_dir]}" "{config[genome_dir]}"
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

rule call_rampage_TSS:
    input:
        ipfile = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/{wildcards.file_type}__{wildcards.sample_name}.bam",
        inputfile = lambda wildcards: f"{RESULTS_DIR}/RNA/mapped/{assign_rna_input(wildcards)}.bam",
        genome_stats = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/genome_stats.json"
    output:
        peakfile = f"{RESULTS_DIR}/RNA/TSS/TSS__{{file_type}}__{{sample_name}}_peaks.narrowPeak"
    wildcard_constraints:
        env = "RNA",
        file_type = "final|merged"
    params:
        ipname = lambda wildcards: wildcards.sample_name,
        inputname = lambda wildcards: assign_rna_input(wildcards),
        filetype = lambda wildcards: wildcards.file_type,
        params = config["rampage_calltss"]['params'],
        genomesize_override = lambda wildcards: config["genomes"][parse_sample_name(wildcards.sample_name)['ref_genome']].get('genomesize', '')
    log:
        temp(return_log_rna("{sample_name}", "{file_type}__TSS_calling", "SE"))
    conda: CONDA_ENV_CHIP
    shell:
        """
        {{
        if [[ -n "{params.genomesize_override}" ]]; then
            gsize="{params.genomesize_override}"
        else
            gsize=$(python3 -c "import json; print(json.load(open('{input.genome_stats}'))['effective_size'])")
        fi
        printf "\nCalling TSS (narrow peaks) for {params.ipname} (vs {params.inputname}) using macs2 version:\n"
        macs2 --version
        macs2 callpeak -t {input.ipfile} -c {input.inputfile} -f BAM -g $gsize {params.params} -n TSS__{params.filetype}__{params.ipname} --outdir {config[output_dir]}/RNA/TSS/
        }} 2>&1 | tee -a "{log}"
        """

rule all_rna:
    input:
        final = lambda wildcards: define_final_rna_output(wildcards.ref_genome)
    output:
        touch = f"{RESULTS_DIR}/RNA/chkpts/RNA_analysis__{{analysis_name}}__{{ref_genome}}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """
