import json as _json

CONDA_ENV_MC=os.path.join(REPO_FOLDER,"workflow","envs","epibutton_mc.yaml")

# Build wildcard constraint patterns for dmC vs bisulfite rule routing.
# dmC samples (Assay == "dmC") use modkit/ONT pipeline; all others use Bismark.
# Both the bare Sample_ID (pre-alignment paths) and the genome-qualified
# mapped_name (post-alignment paths) must match, or a qualified token silently
# fails the constraint and its producing rule never fires.
_dmc_ids = sorted(set(samples.loc[samples["Assay"] == "dmC", "Sample_ID"]).union(
    samples.loc[samples["Assay"] == "dmC", "mapped_name"]
    if "mapped_name" in samples.columns else []
))
if _dmc_ids:
    _DMC_WC = "(?:" + "|".join(re.escape(s) for s in _dmc_ids) + ")"
    _NON_DMC_WC = "(?!(?:" + "|".join(re.escape(s) for s in _dmc_ids) + ")$).*"
else:
    _DMC_WC = "(?!x)x"  # matches nothing
    _NON_DMC_WC = ".*"

# dmC-specific rules take priority over bisulfite rules when both match
ruleorder: make_mc_stats_dmc > make_mc_stats_se

def return_log_mc(sample_name, step, paired):
    return os.path.join(REPO_FOLDER, RESULTS_DIR,"mC","logs",f"tmp__{sample_name}__{step}__{paired}.log")

def parameters_for_mc(sample_name):
    assay = get_sample_info_from_name(sample_name, samples, 'Assay')
    options = {"WGBS", "WGBS_nd", "PBAT", "EMseq", "dmC"}
    return assay if assay in options else "default"

def is_dmc_sample(sample_name):
    """Check if a sample uses direct methylation (dmC) workflow (vs bisulfite)."""
    return get_sample_info_from_name(sample_name, samples, 'Assay') == "dmC"

def define_cx_report_input(wildcards):
    """Get CX_report path for a sample (used by bigwig generation and replicate merging).

    Routes to appropriate path based on sample_type:
    - Bismark samples: {config[output_dir]}/mC/methylcall/...deduplicated.CX_report.txt.gz
    - dmC (direct methylation): {config[output_dir]}/mC/dmc/cx_report__...CX_report.txt.gz
    - Merged replicates: {config[output_dir]}/mC/methylcall/...merged.CX_report.txt.gz
    """
    sname = wildcards.sample_name
    parsed = parse_sample_name(sname)
    if parsed['replicate'] == "merged":
        return f"{RESULTS_DIR}/mC/methylcall/{sname}.merged.CX_report.txt.gz"
    elif is_dmc_sample(sname):
        # dmC samples: use converted CX_report from bedMethyl
        return f"{RESULTS_DIR}/mC/dmc/cx_report__{sname}.CX_report.txt.gz"
    else:
        # Bismark samples: use deduplicated CX_report
        return f"{RESULTS_DIR}/mC/methylcall/{sname}.deduplicated.CX_report.txt.gz"

def cx_report_for_rep_sid(rep_sid):
    """Per-replicate CX_report path (Bismark deduplicated vs dmC converted)."""
    if is_dmc_sample(rep_sid):
        return f"{RESULTS_DIR}/mC/dmc/cx_report__{rep_sid}.CX_report.txt.gz"
    return f"{RESULTS_DIR}/mC/methylcall/{rep_sid}.deduplicated.CX_report.txt.gz"

def define_DMR_samples(sample_name):
    """Get CX_report files for DMR analysis.

    For Bismark samples: returns deduplicated CX_report files
    For dmC samples: returns converted CX_report files (from bedMethyl)
    Uses get_replicate_sample_ids() to find per-replicate Sample_IDs.
    """
    rep_sids = get_replicate_sample_ids(sample_name, samples)
    if not rep_sids:
        # sample_name is itself a per-replicate name, not an analysis name
        rep_sids = [sample_name]
    return [cx_report_for_rep_sid(sid) for sid in rep_sids]

def rep_rds_for_group(group_name, mc_context):
    """Per-(replicate, context) RDS cache paths for one analysis group."""
    rep_sids = get_replicate_sample_ids(group_name, samples)
    if not rep_sids:
        rep_sids = [group_name]
    return [f"{RESULTS_DIR}/mC/pools/per_rep/{sid}__{mc_context}.rds" for sid in rep_sids]


def call_dmrs_pair_script():
    """Resolve the per-(pair, context) DMR-call script based on dmr_caller config.

    Both scripts share the same CLI signature, so the rule's shell line
    is identical regardless of caller — only the script path differs.
    """
    caller = config.get('dmr_caller', 'metilene')
    if caller == 'metilene':
        return os.path.join(REPO_FOLDER, "workflow", "scripts", "R_call_DMRs_pair_metilene.R")
    if caller == 'dmrcaller':
        return os.path.join(REPO_FOLDER, "workflow", "scripts", "R_call_DMRs_pair.R")
    raise ValueError(
        f"Unknown dmr_caller: {caller!r}. Expected 'metilene' or 'dmrcaller'."
    )


_DMR_THRESHOLD_DEFAULTS = {
    'min_diff':      {'CG': 0.3, 'CHG': 0.2, 'CHH': 0.1},
    'min_cytosines': 5,
    'bin_size':      200,
    'q_value':       0.01,
    'min_gap':       200,
    'min_size':      50,
    'min_reads':     3,
    'maxdist':       300,
    'valley':        {'CG': 0.7, 'CHG': 0.7, 'CHH': 0.3},
    'maxseg':        {'CG': -1,  'CHG': -1,  'CHH': 10000},
}

def get_dmr_threshold(key, context=None):
    """Return a DMR threshold from config, falling back to plant-tuned defaults."""
    thresholds = config.get('dmr_thresholds', {})
    val = thresholds.get(key, _DMR_THRESHOLD_DEFAULTS[key])
    if isinstance(val, dict):
        default_per_ctx = _DMR_THRESHOLD_DEFAULTS[key]
        return val.get(context, default_per_ctx.get(context)) if context else val
    return val


# True when the user has requested data-driven min_diff calibration.
_dmr_min_diff_auto = (
    str(config.get('dmr_thresholds', {}).get('min_diff', '')).lower() == 'auto'
)


def get_methylation_contexts():
    """Return the list of methylation contexts to analyze (subset of CG, CHG, CHH).

    Configured via ``methylation_contexts: ["CG", "CHG", "CHH"]`` — a list
    of contexts to generate bigwigs for, call DMRs in, and run PCA on.
    Animal genomes can pass ``["CG"]`` to skip empty mCHG/mCHH plots.
    Defaults to the full plant set when the key is not set. Output order
    is normalized to CG → CHG → CHH.

    Subcontexts beyond CG/CHG/CHH (CAG, CAA, etc.) are not yet supported
    and raise an error.
    """
    valid = ["CG", "CHG", "CHH"]
    contexts = config.get("methylation_contexts", valid)
    if not isinstance(contexts, list) or not contexts:
        raise ValueError(
            "methylation_contexts must be a non-empty list, e.g. "
            "[\"CG\", \"CHG\", \"CHH\"] or [\"CG\"] for animal genomes"
        )
    valid_set = set(valid)
    for c in contexts:
        if c not in valid_set:
            raise ValueError(
                f"methylation_contexts entry '{c}' not in {valid} "
                "(CAG/CAA and other subcontexts are not yet supported; "
                "support deferred to a future PR)"
            )
    return [c for c in valid if c in contexts]

def define_final_mC_output(ref_genome):
    qc_option = config["QC_option"]
    analysis = config['full_analysis']
    trimmed_fastqs = config['trimmed_fastqs']
    map_files = []
    dmr_files = []
    bigwig_files = []
    qc_files = []
    ont_files = []
    filtered_rep_samples = samples[ (samples['env'] == 'mC') & (samples['ref_genome'] == ref_genome) ].copy()

    for _, row in filtered_rep_samples.iterrows():
        sname = row['sample_name']
        mname = row['mapped_name']
        paired = row['paired']
        assay = row['Assay']

        # dmC samples use direct methylation workflow
        if assay == "dmC":
            bigwig_files.append(f"{RESULTS_DIR}/mC/chkpts/bigwig__{mname}.done")
            ont_files.append(f"{RESULTS_DIR}/mC/dmc/summary__{mname}.txt")  # modkit summary
        else:
            # Bismark workflow
            bigwig_files.append(f"{RESULTS_DIR}/mC/chkpts/bigwig__{mname}.done")
            if paired == "PE":
                map_files.append(f"{RESULTS_DIR}/mC/reports/final_report_pe__{mname}.html")
                qc_files.append(f"{RESULTS_DIR}/mC/reports/trim__{sname}__R1_fastqc.html") # fastqc of trimmed Read1 fastq files
                qc_files.append(f"{RESULTS_DIR}/mC/reports/trim__{sname}__R2_fastqc.html") # fastqc of trimmed Read2 fastq files
                if not trimmed_fastqs:
                    qc_files.append(f"{RESULTS_DIR}/mC/reports/raw__{sname}__R1_fastqc.html") # fastqc of raw Read1 fastq file
                    qc_files.append(f"{RESULTS_DIR}/mC/reports/raw__{sname}__R2_fastqc.html") # fastqc of raw Read2 fastq file
            else:
                map_files.append(f"{RESULTS_DIR}/mC/reports/final_report_se__{mname}.html")
                qc_files.append(f"{RESULTS_DIR}/mC/reports/trim__{sname}__R0_fastqc.html") # fastqc of trimmed (Read0) fastq files
                if not trimmed_fastqs:
                    qc_files.append(f"{RESULTS_DIR}/mC/reports/raw__{sname}__R0_fastqc.html") # fastqc of raw (Read0) fastq file

    filtered_analysis_samples = analysis_samples[ (analysis_samples['env'] == 'mC') & (analysis_samples['ref_genome'] == ref_genome) ].copy()
    for _, row in filtered_analysis_samples.iterrows():
        aname = row['sample_name']
        rep_sids = get_replicate_sample_ids(aname, samples)
        if len(rep_sids) >= 2:
            bigwig_files.append(f"{RESULTS_DIR}/mC/chkpts/bigwig__{aname}.done") # merged bigwig files

    # DMR analysis: all sample types use DMRcaller via unified CX_report format
    for a, b in combinations(filtered_analysis_samples.itertuples(index=False), 2):
        a_dict = a._asdict()
        b_dict = b._asdict()
        sample1 = a_dict['sample_name']
        sample2 = b_dict['sample_name']
        dmr_files.append(f"{RESULTS_DIR}/mC/DMRs/summary__{sample1}__vs__{sample2}__DMRs.txt")

    results = map_files + bigwig_files + ont_files

    if qc_option == "all":
        results += qc_files

    if analysis:
        results += dmr_files

    return results

rule make_bismark_indices:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa"
    output:
        indices = directory(f"{GENOMES_DIR}/{{ref_genome}}/Bisulfite_Genome")
    params:
        limthreads = lambda wildcards, threads: max(1, threads // 2)
    log:
        temp(os.path.join(REPO_FOLDER, RESULTS_DIR,"logs","bismark_index_{ref_genome}.log"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        printf "\nBuilding bismark index directory for {wildcards.ref_genome}\n"
        if [[ {params.limthreads} -gt 1 ]]; then
            bismark_genome_preparation --parallel {params.limthreads} --bowtie2 --genomic_composition {config[genome_dir]}/{wildcards.ref_genome}
        else
            bismark_genome_preparation --bowtie2 --genomic_composition {config[genome_dir]}/{wildcards.ref_genome}
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule bismark_map_pe:
    input:
        fastq1 = lambda wildcards: f"{RESULTS_DIR}/mC/fastq/trim__{base_sample(wildcards.sample_name)}__R1.fastq.gz",
        fastq2 = lambda wildcards: f"{RESULTS_DIR}/mC/fastq/trim__{base_sample(wildcards.sample_name)}__R2.fastq.gz",
        indices = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/Bisulfite_Genome"
    output:
        temp_bamfile = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/trim__{{sample_name}}__R1_bismark_bt2_pe.bam"),
        bamfile = maybe_temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/PE__{{sample_name}}.deduplicated.bam", config.get('keep_final_bams', True)),
        cx_report = maybe_temp(f"{RESULTS_DIR}/mC/mapped/PE__{{sample_name}}.deduplicated.CX_report.txt.gz", config.get('keep_cx_reports', False)),
        metrics_alignement = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/trim__{{sample_name}}__R1_bismark_bt2_PE_report.txt"),
        metrics_dedup = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/PE__{{sample_name}}.deduplication_report.txt"),
        # Companion reports consumed by bismark2report in make_mc_stats_pe.
        # Declared here so Snakemake guarantees/tracks them (a missing one now
        # fails this rule with a clear MissingOutputException instead of a
        # cryptic bismark2report crash downstream). splitting/M-bias land at the
        # methylation-extractor's -o dir (mapped/, no subdir); nucleotide_stats
        # is produced by the explicit bam2nuc step below (Bismark 3.x Rust
        # aligner recognises --nucleotide_coverage but does not yet implement
        # it, so we compute it separately from the deduplicated BAM).
        splitting_report = temp(f"{RESULTS_DIR}/mC/mapped/PE__{{sample_name}}.deduplicated_splitting_report.txt"),
        mbias_report = temp(f"{RESULTS_DIR}/mC/mapped/PE__{{sample_name}}.deduplicated.M-bias.txt"),
        nucleotide_stats = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/PE__{{sample_name}}.deduplicated.nucleotide_stats.txt")
    wildcard_constraints:
        sample_name = _NON_DMC_WC  # Exclude dmC (direct methylation) samples
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome_path = lambda wildcards: os.path.join(REPO_FOLDER, GENOMES_DIR,parse_sample_name(wildcards.sample_name)['ref_genome']),
        base = lambda wildcards: base_sample(wildcards.sample_name),
        mapping = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['map_pe'],
        process = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['process_pe'],
        prefix = lambda wildcards: f"{RESULTS_DIR}/mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: max(1, threads // 2)
    log:
        temp(return_log_mc("{sample_name}", "mapping", "PE"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        # Unset OMP_NUM_THREADS: Snakemake sets this to the thread count, but
        # bowtie2 uses it for internal parallelism which conflicts with bismark's
        # own --multicore process management, causing near-zero mapping rates.
        unset OMP_NUM_THREADS
        printf "\nAligning {params.sample_name} with bismark/bowtie2\n"
        # gzip output omitted: bismark rejects gzipped output in PBAT mode,
        # and the flag only affects intermediate temp files (not the final
        # BAM), so the disk cost is small and consistent across all mC types.
        bismark --genome {params.ref_genome_path} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --temp_dir {params.prefix} -1 {input.fastq1} -2 {input.fastq2}
        # bismark names its outputs after the input FASTQ, which is genome-free
        # (one trim serves every reference), while our paths carry the genome.
        mv "{params.prefix}/trim__{params.base}__R1_bismark_bt2_pe.bam" "{output.temp_bamfile}"
        mv "{params.prefix}/trim__{params.base}__R1_bismark_bt2_PE_report.txt" "{output.metrics_alignement}"
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -p --output_dir {params.prefix}/ -o "PE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nComputing nucleotide-coverage stats (bam2nuc)\n"
        bam2nuc --genome_folder {params.ref_genome_path} --dir {params.prefix} {output.bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -p --comprehensive -o {config[output_dir]}/mC/mapped/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder {params.ref_genome_path} {output.bamfile}
        rm -f {config[output_dir]}/mC/mapped/C*context_PE__{params.sample_name}*
        rm -f {config[output_dir]}/mC/mapped/PE__{params.sample_name}*bismark.cov*
        rm -f {config[output_dir]}/mC/mapped/PE__{params.sample_name}*bedGraph*
        }} 2>&1 | tee -a "{log}"
        """

rule bismark_map_se:
    input:
        fastq0 = lambda wildcards: f"{RESULTS_DIR}/mC/fastq/trim__{base_sample(wildcards.sample_name)}__R0.fastq.gz",
        indices = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/Bisulfite_Genome"
    output:
        temp_bamfile = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/trim__{{sample_name}}__R0_bismark_bt2.bam"),
        bamfile = maybe_temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/SE__{{sample_name}}.deduplicated.bam", config.get('keep_final_bams', True)),
        cx_report = maybe_temp(f"{RESULTS_DIR}/mC/mapped/SE__{{sample_name}}.deduplicated.CX_report.txt.gz", config.get('keep_cx_reports', False)),
        metrics_map = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/trim__{{sample_name}}__R0_bismark_bt2_SE_report.txt"),
        metrics_dedup = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/SE__{{sample_name}}.deduplication_report.txt"),
        # Companion reports consumed by bismark2report in make_mc_stats_se (see
        # the matching note in bismark_map_pe). nucleotide_stats comes from the
        # explicit bam2nuc step (Bismark 3.x aligner does not yet implement
        # --nucleotide_coverage).
        splitting_report = temp(f"{RESULTS_DIR}/mC/mapped/SE__{{sample_name}}.deduplicated_splitting_report.txt"),
        mbias_report = temp(f"{RESULTS_DIR}/mC/mapped/SE__{{sample_name}}.deduplicated.M-bias.txt"),
        nucleotide_stats = temp(f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/SE__{{sample_name}}.deduplicated.nucleotide_stats.txt")
    wildcard_constraints:
        sample_name = _NON_DMC_WC  # Exclude dmC (direct methylation) samples
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome_path = lambda wildcards: os.path.join(REPO_FOLDER, GENOMES_DIR,parse_sample_name(wildcards.sample_name)['ref_genome']),
        base = lambda wildcards: base_sample(wildcards.sample_name),
        mapping = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['map_se'],
        process = lambda wildcards: config["mC_mapping"][parameters_for_mc(wildcards.sample_name)]['process_se'],
        prefix = lambda wildcards: f"{RESULTS_DIR}/mC/mapped/{wildcards.sample_name}",
        limthreads = lambda wildcards, threads: max(1, threads // 2)
    log:
        temp(return_log_mc("{sample_name}", "mapping", "SE"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        # Unset OMP_NUM_THREADS: Snakemake sets this to the thread count, but
        # bowtie2 uses it for internal parallelism which conflicts with bismark's
        # own --multicore process management, causing near-zero mapping rates.
        unset OMP_NUM_THREADS
        printf "\nAligning {params.sample_name} with bismark/bowtie2\n"
        # --gzip omitted for symmetry with bismark_map_pe (see that rule).
        bismark --genome {params.ref_genome_path} {params.mapping} --local --multicore {params.limthreads} -o {params.prefix} --temp_dir {params.prefix} {input.fastq0}
        # bismark names its outputs after the input FASTQ, which is genome-free
        # (one trim serves every reference), while our paths carry the genome.
        mv "{params.prefix}/trim__{params.base}__R0_bismark_bt2.bam" "{output.temp_bamfile}"
        mv "{params.prefix}/trim__{params.base}__R0_bismark_bt2_SE_report.txt" "{output.metrics_map}"
        printf "\nDeduplicating with bismark\n"
        deduplicate_bismark -s --output_dir {params.prefix} -o "SE__{params.sample_name}" --bam {output.temp_bamfile}
        printf "\nComputing nucleotide-coverage stats (bam2nuc)\n"
        bam2nuc --genome_folder {params.ref_genome_path} --dir {params.prefix} {output.bamfile}
        printf "\nCalling mC for {params.sample_name}"
        bismark_methylation_extractor -s --comprehensive -o {config[output_dir]}/mC/mapped/ {params.process} --gzip --multicore {params.limthreads} --cytosine_report --CX --genome_folder {params.ref_genome_path} {output.bamfile}
        rm -f {config[output_dir]}/mC/mapped/C*context_SE__{params.sample_name}*
        rm -f {config[output_dir]}/mC/mapped/SE__{params.sample_name}*bismark.cov*
        rm -f {config[output_dir]}/mC/mapped/SE__{params.sample_name}*bedGraph*
        }} 2>&1 | tee -a "{log}"
        """

rule pe_or_se_mc_dispatch:
    input:
        lambda wildcards: assign_mapping_paired(wildcards, "bismark_map", "cx_report")
    output:
        cx_report = f"{RESULTS_DIR}/mC/methylcall/{{sample_name}}.deduplicated.CX_report.txt.gz",
        touch = f"{RESULTS_DIR}/mC/chkpts/map_mC__{{sample_name}}.done"
    wildcard_constraints:
        sample_name = _NON_DMC_WC  # Exclude dmC (direct methylation) samples
    localrule: True
    shell:
        """
        mv {input} {output.cx_report}
        touch {output.touch}
        """

rule make_mc_stats_pe:
    input:
        metrics_trim = lambda wildcards: f"{RESULTS_DIR}/mC/reports/trim_pe__{base_sample(wildcards.sample_name)}.json",
        metrics_map = f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/trim__{{sample_name}}__R1_bismark_bt2_PE_report.txt",
        metrics_dedup = f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/PE__{{sample_name}}.deduplication_report.txt",
        cx_report = f"{RESULTS_DIR}/mC/methylcall/{{sample_name}}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes",
        splitting_report = f"{RESULTS_DIR}/mC/mapped/PE__{{sample_name}}.deduplicated_splitting_report.txt",
        mbias_report = f"{RESULTS_DIR}/mC/mapped/PE__{{sample_name}}.deduplicated.M-bias.txt",
        nucleotide_stats = f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/PE__{{sample_name}}.deduplicated.nucleotide_stats.txt"
    wildcard_constraints:
        sample_name = _NON_DMC_WC  # Exclude dmC (direct methylation) samples
    output:
        stat_file = f"{RESULTS_DIR}/mC/reports/summary_mC_PE_mapping_stats_{{sample_name}}.txt",
        reportfile = f"{RESULTS_DIR}/mC/reports/final_report_pe__{{sample_name}}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        group = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        levels = lambda wildcards: parse_sample_name(wildcards.sample_name)['levels_label'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"{RESULTS_DIR}/mC/mapped/{wildcards.sample_name}",
        trimmed_fastq = config['trimmed_fastqs']
    log:
        temp(return_log_mc("{sample_name}", "making_stats", "PE"))
    conda: CONDA_ENV_MC
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(python3 -c "import json; print(json.load(open('{input.metrics_trim}'))['summary']['before_filtering']['total_reads'] // 2)")
        else
            tot=$(grep "Sequence pairs analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        fi
        filt=$(grep "Sequence pairs analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Sequence pairs did not map uniquely" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Number of paired-end alignments with a unique best hit" "{input.metrics_map}" | awk '{{print $NF}}')
        uniq=$(grep "Total count of deduplicated leftover sequences" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > {output.stat_file}
        ## Can change the name of the plastid chromosome to calculate non-conversion rate
        zcat {input.cx_report} | awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v s={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v x=${{tot}} -v y=${{filt}} -v z=${{allmap}} -v u=${{uniq}} '{{a+=1; b=$4+$5; i+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {{m+=$4; n+=b;}}; if (b>0) {{c+=1; d+=b;}}; if (b>2) e+=1}} END {{if (n>0) {{o=m/n*100;}} else o="NA"; print l,t,s,r,g,x,y" ("y/x*100"%)",z" ("z/x*100"%)",u" ("u/x*100"%)",c/a*100,e/a*100,i/a,d/c,o}}' >> "{output.stat_file}"

        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_pe__{params.sample_name}.html" --dir {config[output_dir]}/mC/reports/ --alignment_report {input.metrics_map} --dedup_report {input.metrics_dedup} --splitting_report {input.splitting_report} --mbias_report {input.mbias_report} --nucleotide_report {input.nucleotide_stats}
        cp {config[output_dir]}/mC/mapped/PE__"{params.sample_name}"*.txt {config[output_dir]}/mC/reports/
        cp {params.prefix}/trim__"{params.sample_name}"*.txt {config[output_dir]}/mC/reports/
        """

rule make_mc_stats_se:
    input:
        metrics_trim = lambda wildcards: f"{RESULTS_DIR}/mC/reports/trim_se__{base_sample(wildcards.sample_name)}.json",
        metrics_map = f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/trim__{{sample_name}}__R0_bismark_bt2_SE_report.txt",
        metrics_dedup = f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/SE__{{sample_name}}.deduplication_report.txt",
        cx_report = f"{RESULTS_DIR}/mC/methylcall/{{sample_name}}.deduplicated.CX_report.txt.gz",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes",
        splitting_report = f"{RESULTS_DIR}/mC/mapped/SE__{{sample_name}}.deduplicated_splitting_report.txt",
        mbias_report = f"{RESULTS_DIR}/mC/mapped/SE__{{sample_name}}.deduplicated.M-bias.txt",
        nucleotide_stats = f"{RESULTS_DIR}/mC/mapped/{{sample_name}}/SE__{{sample_name}}.deduplicated.nucleotide_stats.txt"
    wildcard_constraints:
        sample_name = _NON_DMC_WC  # Exclude dmC (direct methylation) samples
    output:
        stat_file = f"{RESULTS_DIR}/mC/reports/summary_mC_SE_mapping_stats_{{sample_name}}.txt",
        reportfile = f"{RESULTS_DIR}/mC/reports/final_report_se__{{sample_name}}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        group = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        levels = lambda wildcards: parse_sample_name(wildcards.sample_name)['levels_label'],
        sample_type = lambda wildcards: parse_sample_name(wildcards.sample_name)['sample_type'],
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        prefix = lambda wildcards: f"{RESULTS_DIR}/mC/mapped/{wildcards.sample_name}",
        trimmed_fastq = config['trimmed_fastqs']
    log:
        temp(return_log_mc("{sample_name}", "making_stats", "SE"))
    conda: CONDA_ENV_MC
    shell:
        """
        printf "\nMaking mapping statistics summary\n"
        if [[ "{params.trimmed_fastq}" == "False" ]]; then
            tot=$(python3 -c "import json; print(json.load(open('{input.metrics_trim}'))['summary']['before_filtering']['total_reads'])")
        else
            tot=$(grep "Sequences analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        fi
        filt=$(grep "Sequences analysed in total" "{input.metrics_map}" | awk '{{print $NF}}')
        multi=$(grep "Sequences did not map uniquely" "{input.metrics_map}" | awk '{{print $NF}}')
        single=$(grep "Number of alignments with a unique best hit" "{input.metrics_map}" | awk '{{print $NF}}')
        uniq=$(grep "Total count of deduplicated leftover sequences" {input.metrics_dedup} | awk -v FS=":" 'END {{print $2}}' | awk '{{print $1}}')
        allmap=$((single+multi))
        printf "Group\tLevels\tSample\tRep\tReference_genome\tTotal_reads\tPassing_filtering\tAll_mapped_reads\tUniquely_mapped_reads\tPercentage_covered\tPercentage_covered_min3reads\tAverage_coverage_all\tAverage_coverage_covered\tNon_conversion_rate(Pt/Lambda)\n" > {output.stat_file}
        ## Can change the name of the plastid chromosome to calculate non-conversion rate
        zcat {input.cx_report} | awk -v OFS="\t" -v l={params.group} -v t={params.levels} -v s={params.sample_type} -v r={params.replicate} -v g={params.ref_genome} -v x=${{tot}} -v y=${{filt}} -v z=${{allmap}} -v u=${{uniq}} '{{a+=1; b=$4+$5; i+=b; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {{m+=$4; n+=b;}}; if (b>0) {{c+=1; d+=b;}}; if (b>2) e+=1}} END {{if (n>0) {{o=m/n*100;}} else o="NA"; print l,t,s,r,g,x,y" ("y/x*100"%)",z" ("z/x*100"%)",u" ("u/x*100"%)",c/a*100,e/a*100,i/a,d/c,o}}' >> "{output.stat_file}"

        printf "\nMaking final html report for {params.sample_name}\n"
        bismark2report -o "final_report_se__{params.sample_name}.html" --dir {config[output_dir]}/mC/reports/ --alignment_report {input.metrics_map} --dedup_report {input.metrics_dedup} --splitting_report {input.splitting_report} --mbias_report {input.mbias_report} --nucleotide_report {input.nucleotide_stats}
        cp {config[output_dir]}/mC/mapped/SE__"{params.sample_name}"*.txt {config[output_dir]}/mC/reports/
        cp {params.prefix}/trim__"{params.sample_name}"*.txt {config[output_dir]}/mC/reports/
        """

def get_cx_reports_for_merging(wildcards):
    """Get CX_report paths for all replicates of a sample, for merging.

    Routes to appropriate paths based on sample_type:
    - Bismark samples: {config[output_dir]}/mC/methylcall/...deduplicated.CX_report.txt.gz
    - dmC (direct methylation): {config[output_dir]}/mC/dmc/cx_report__...CX_report.txt.gz
    """
    aname = wildcards.sample_name
    rep_sids = get_replicate_sample_ids(aname, samples)

    result = []
    for sid in rep_sids:
        if is_dmc_sample(sid):
            result.append(f"{RESULTS_DIR}/mC/dmc/cx_report__{sid}.CX_report.txt.gz")
        else:
            result.append(f"{RESULTS_DIR}/mC/methylcall/{sid}.deduplicated.CX_report.txt.gz")
    return result

rule merging_mc_replicates:
    input:
        report_files = get_cx_reports_for_merging
    output:
        bedfile = temp(f"{RESULTS_DIR}/mC/methylcall/{{sample_name}}.bed"),
        tempmergefile = temp(f"{RESULTS_DIR}/mC/methylcall/{{sample_name}}.merged.CX_report.txt"),
        mergefile = temp(f"{RESULTS_DIR}/mC/methylcall/{{sample_name}}.merged.CX_report.txt.gz")
    params:
        sname = lambda wildcards: wildcards.sample_name
    log:
        temp(return_log_mc("{sample_name}", "merging_reps", ""))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        printf "\nMerging replicates of {params.sname}\n"
        zcat {input.report_files} | sort -k1,1 -k2,2n | awk -v OFS="\t" '{{print $1,$2-1,$2,$3,$4,$5,$6,$7}}' > {output.bedfile}
		bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i {output.bedfile} | awk -v OFS="\t" '{{print $1,$3,$4,$5,$6,$7,$8}}' > {output.tempmergefile}
        pigz -p {threads} "{output.tempmergefile}" -c > "{output.mergefile}"
        }} 2>&1 | tee -a "{log}"
        """

rule make_mc_bigwig_files:
    """Generate bigwig files from CX_report data.

    Handles all sample types (Bismark and dmC) via unified CX_report format.
    """
    input:
        cx_report = define_cx_report_input,
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        bigwigcg = f"{RESULTS_DIR}/mC/tracks/{{sample_name}}__CG.bw",
        bigwigchg = f"{RESULTS_DIR}/mC/tracks/{{sample_name}}__CHG.bw",
        bigwigchh = f"{RESULTS_DIR}/mC/tracks/{{sample_name}}__CHH.bw",
        touch = f"{RESULTS_DIR}/mC/chkpts/bigwig__{{sample_name}}.done"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome'],
        # Space-separated list of contexts to actually generate bigwig
        # data for. Contexts NOT in the list still get an empty placeholder
        # bigwig (downstream rules require all three .bw outputs to exist).
        contexts = " ".join(get_methylation_contexts())
    log:
        temp(return_log_mc("{sample_name}", "bigwig", ""))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        outdir="{config[output_dir]}/mC/tracks"
        sname="{params.sample_name}"
        active_contexts="{params.contexts}"

        # Pre-create empty bedGraph files so awk redirection doesn't leave
        # missing files when a context/strand has no data (e.g. combined-
        # strand dmC input has no minus-strand calls).
        for context in CG CHG CHH; do
            > "$outdir/${{sname}}__${{context}}.bedGraph"
            for strand in plus minus; do
                > "$outdir/${{sname}}__${{context}}__${{strand}}.bedGraph"
            done
        done

        # Demux methylation calls into per-context bedGraphs (and per-strand
        # variants). The CX_report has all three contexts interleaved; the
        # awk output redirection writes each line to the matching context
        # file. Cheap to do unconditionally — empty bedGraphs are dropped
        # below for inactive contexts.
        zcat {input.cx_report} | awk -v OFS="\t" -v s=$sname '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 >> "'"$outdir"'/"s"__CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 >> "'"$outdir"'/"s"__CHG.bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 >> "'"$outdir"'/"s"__CG.bedGraph"}}'
        for strand in plus minus; do
            case "$strand" in
                plus)  sign="+";;
                minus) sign="-";;
            esac
            zcat {input.cx_report} | awk -v n="$sign" '$3==n' | awk -v OFS="\t" -v s=$sname -v d=$strand '($4+$5)>0 {{a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 >> "'"$outdir"'/"s"__CHH__"d".bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 >> "'"$outdir"'/"s"__CHG__"d".bedGraph"; else if ($6=="CG") print $1,$2-1,$2,$4/a*100 >> "'"$outdir"'/"s"__CG__"d".bedGraph"}}'
        done

        for context in CG CHG CHH; do
            # Skip contexts the user didn't ask for: emit a 1-bp placeholder
            # bigwig so snakemake's declared outputs all exist (downstream
            # heatmaps/PCA rules gate on the presence of these files).
            if ! [[ " $active_contexts " == *" $context "* ]]; then
                printf "\nContext $context not in methylation_contexts; emitting empty placeholder bigwig\n"
                chrom=$(head -1 {input.chrom_sizes} | cut -f1)
                printf "%s\t0\t1\t0\n" "$chrom" > "$outdir/empty__${{sname}}__${{context}}.bg"
                bedGraphToBigWig "$outdir/empty__${{sname}}__${{context}}.bg" {input.chrom_sizes} "$outdir/${{sname}}__${{context}}.bw"
                rm -f "$outdir/empty__${{sname}}__${{context}}.bg"
                for strand in plus minus; do
                    printf "%s\t0\t1\t0\n" "$chrom" > "$outdir/empty__${{sname}}__${{context}}__${{strand}}.bg"
                    bedGraphToBigWig "$outdir/empty__${{sname}}__${{context}}__${{strand}}.bg" {input.chrom_sizes} "$outdir/${{sname}}__${{context}}__${{strand}}.bw"
                    rm -f "$outdir/empty__${{sname}}__${{context}}__${{strand}}.bg"
                done
                continue
            fi

            printf "\nMaking bigwig files of $context context for $sname\n"
            if [[ -s "$outdir/${{sname}}__${{context}}.bedGraph" ]]; then
                LC_COLLATE=C sort -k1,1 -k2,2n "$outdir/${{sname}}__${{context}}.bedGraph" > "$outdir/sorted__${{sname}}__${{context}}.bedGraph"
                bedGraphToBigWig "$outdir/sorted__${{sname}}__${{context}}.bedGraph" {input.chrom_sizes} "$outdir/${{sname}}__${{context}}.bw"
            else
                printf "No data for $context context — creating empty bigwig\n"
                chrom=$(head -1 {input.chrom_sizes} | cut -f1)
                printf "%s\t0\t1\t0\n" "$chrom" > "$outdir/empty__${{sname}}__${{context}}.bg"
                bedGraphToBigWig "$outdir/empty__${{sname}}__${{context}}.bg" {input.chrom_sizes} "$outdir/${{sname}}__${{context}}.bw"
                rm -f "$outdir/empty__${{sname}}__${{context}}.bg"
            fi
            for strand in plus minus; do
                printf "\nMaking $strand strand bigwig files of $context context for $sname\n"
                if [[ -s "$outdir/${{sname}}__${{context}}__${{strand}}.bedGraph" ]]; then
                    LC_COLLATE=C sort -k1,1 -k2,2n "$outdir/${{sname}}__${{context}}__${{strand}}.bedGraph" > "$outdir/sorted__${{sname}}__${{context}}__${{strand}}.bedGraph"
                    bedGraphToBigWig "$outdir/sorted__${{sname}}__${{context}}__${{strand}}.bedGraph" {input.chrom_sizes} "$outdir/${{sname}}__${{context}}__${{strand}}.bw"
                else
                    printf "No data for $context $strand strand — creating empty bigwig\n"
                    chrom=$(head -1 {input.chrom_sizes} | cut -f1)
                    printf "%s\t0\t1\t0\n" "$chrom" > "$outdir/empty__${{sname}}__${{context}}__${{strand}}.bg"
                    bedGraphToBigWig "$outdir/empty__${{sname}}__${{context}}__${{strand}}.bg" {input.chrom_sizes} "$outdir/${{sname}}__${{context}}__${{strand}}.bw"
                    rm -f "$outdir/empty__${{sname}}__${{context}}__${{strand}}.bg"
                fi
            done
        done

        rm -f "$outdir/"*"$sname"*bedGraph*
        touch {output.touch}
        }} 2>&1 | tee -a "{log}"
        """

if config.get('custom_script_dmrs', False):
    # Custom parameter-sweep mode — keep the single-job rule that
    # iterates parameter combinations inside R_call_DMRs_custom.R.
    rule call_DMRs_pairwise:
        """Call DMRs between two samples using DMRcaller's custom-sweep script."""
        input:
            sample1 = lambda wildcards: define_DMR_samples(wildcards.sample1),
            sample2 = lambda wildcards: define_DMR_samples(wildcards.sample2),
            chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}/chrom.sizes"
        output:
            dmr_summary = f"{RESULTS_DIR}/mC/DMRs/summary__{{sample1}}__vs__{{sample2}}__DMRs.txt"
        params:
            script = os.path.join(REPO_FOLDER, "workflow", "scripts", "R_call_DMRs_custom.R"),
            contexts = ",".join(get_methylation_contexts()),
            sample1 = lambda wildcards: wildcards.sample1,
            sample2 = lambda wildcards: wildcards.sample2,
            nb_sample1 = lambda wildcards: len(define_DMR_samples(wildcards.sample1)),
            nb_sample2 = lambda wildcards: len(define_DMR_samples(wildcards.sample2)),
            min_cytosines = get_dmr_threshold('min_cytosines'),
            q_value = get_dmr_threshold('q_value'),
            min_gap = get_dmr_threshold('min_gap'),
            min_size = get_dmr_threshold('min_size'),
            min_reads = get_dmr_threshold('min_reads'),
            # Per-context minProportionDifference. The custom sweep does not
            # support "auto" calibration (that lives in the default path), so
            # fall back to the plant-tuned defaults when min_diff is "auto".
            min_diff_spec = ",".join(
                f"{ctx}={_DMR_THRESHOLD_DEFAULTS['min_diff'][ctx] if _dmr_min_diff_auto else get_dmr_threshold('min_diff', ctx)}"
                for ctx in get_methylation_contexts()
            )
        log:
            temp(return_log_mc("{sample1}__vs__{sample2}", "DMRs", ""))
        conda: CONDA_ENV_MC
        shell:
            """
            {{
            printf "running DMRcaller for {params.sample1} vs {params.sample2}\n"
            Rscript "{params.script}" "{threads}" "{input.chrom_sizes}" "{params.contexts}" "{params.sample1}" "{params.sample2}" "{params.nb_sample1}" "{params.nb_sample2}" "{config[output_dir]}" "{params.min_cytosines}" "{params.q_value}" "{params.min_gap}" "{params.min_size}" "{params.min_reads}" "{params.min_diff_spec}" {input.sample1} {input.sample2}
            }} 2>&1 | tee -a "{log}"
            """
else:
    # Two-stage default: per-(replicate, context) caches keep replicates
    # separate so the per-pair rule can call computeDMRsReplicates with a
    # proper condition vector (beta-regression test with biological
    # variance). Pairs where either group has N<2 fall back to pooled
    # computeDMRs inside the R script. Adding a replicate to an existing
    # group creates 3 new caches and invalidates the ~9 pair rows that
    # involve that group (correct: the variance estimate changes);
    # adding a new group leaves the existing 45 pair outputs valid.
    rule cache_mc_replicate_for_context:
        """Cache one replicate's CX report filtered to one methylation
        context as an RDS, for the per-(pair, context) DMR rule."""
        input:
            cx_report = lambda wildcards: cx_report_for_rep_sid(wildcards.rep_sid)
        output:
            rds = f"{RESULTS_DIR}/mC/pools/per_rep/{{rep_sid}}__{{mc_context}}.rds"
        wildcard_constraints:
            mc_context = "CG|CHG|CHH"
        params:
            script = os.path.join(REPO_FOLDER, "workflow", "scripts", "R_cache_mc_replicate_for_context.R")
        log:
            temp(return_log_mc("{rep_sid}", "cache_mc", "{mc_context}"))
        conda: CONDA_ENV_MC
        shell:
            """
            {{
            printf "Caching %s for context %s\n" "{wildcards.rep_sid}" "{wildcards.mc_context}"
            Rscript "{params.script}" "{wildcards.mc_context}" "{output.rds}" "{input.cx_report}"
            }} 2>&1 | tee -a "{log}"
            """

    rule calibrate_dmr_min_diff:
        """Estimate the noise-floor sigma for the auto min_diff calibration.

        Computes the SD of per-position between-group methylation differences.
        At most positions no real DMR exists, so the SD is measurement-noise
        dominated and naturally scales with the context's baseline methylation
        (CG > CHG > CHH). Threshold = sigma_n * sigma, floored at 0.05.
        Only triggered when dmr_thresholds.min_diff: "auto".
        """
        input:
            reps1 = lambda wildcards: rep_rds_for_group(wildcards.sample1, wildcards.mc_context),
            reps2 = lambda wildcards: rep_rds_for_group(wildcards.sample2, wildcards.mc_context),
        output:
            json = temp(f"{RESULTS_DIR}/mC/DMRs/.calib__{{sample1}}__vs__{{sample2}}__{{mc_context}}.json")
        wildcard_constraints:
            mc_context = "CG|CHG|CHH"
        params:
            script  = os.path.join(REPO_FOLDER, "workflow", "scripts", "R_calibrate_dmr_min_diff.R"),
            sigma_n = config.get('dmr_thresholds', {}).get('sigma_n', 3.0),
            n1      = lambda wildcards: len(rep_rds_for_group(wildcards.sample1, wildcards.mc_context)),
        log:
            temp(return_log_mc("{sample1}__vs__{sample2}", "DMR_calib", "{mc_context}"))
        conda: CONDA_ENV_MC
        shell:
            """
            {{
            printf "Calibrating min_diff for %s vs %s (%s)\n" "{wildcards.sample1}" "{wildcards.sample2}" "{wildcards.mc_context}"
            Rscript "{params.script}" "{wildcards.mc_context}" "{output.json}" "{params.sigma_n}" "{params.n1}" {input.reps1} {input.reps2}
            }} 2>&1 | tee -a "{log}"
            """

    rule call_DMRs_for_pair_context:
        """Call DMRs between two groups for one methylation context.

        Caller is selected by the `dmr_caller` config key (default
        "metilene"); the per-replicate RDS inputs from
        cache_mc_replicate_for_context are the same regardless. The
        metilene driver writes a TSV, calls the metilene binary, and
        parses its output; the dmrcaller driver uses
        computeDMRsReplicates (beta-regression with replicate variance,
        with pooled computeDMRs fallback when either group has N<2).
        Both scripts emit the same per-pair-context DMR table + counts
        tsv format, so the aggregator below is caller-agnostic.
        When dmr_thresholds.min_diff: "auto", the calibrate_dmr_min_diff
        rule runs first and its JSON output is read to set min_diff.
        """
        input:
            reps1 = lambda wildcards: rep_rds_for_group(wildcards.sample1, wildcards.mc_context),
            reps2 = lambda wildcards: rep_rds_for_group(wildcards.sample2, wildcards.mc_context),
            chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{get_sample_info_from_name(wildcards.sample1, analysis_samples, 'ref_genome')}/chrom.sizes",
            calib_json = lambda wildcards: (
                f"{RESULTS_DIR}/mC/DMRs/.calib__{wildcards.sample1}__vs__{wildcards.sample2}__{wildcards.mc_context}.json"
                if _dmr_min_diff_auto else []
            ),
        output:
            dmrs   = f"{RESULTS_DIR}/mC/DMRs/{{sample1}}__vs__{{sample2}}__{{mc_context}}_DMRs.txt",
            counts = temp(f"{RESULTS_DIR}/mC/DMRs/.counts__{{sample1}}__vs__{{sample2}}__{{mc_context}}.tsv")
        wildcard_constraints:
            mc_context = "CG|CHG|CHH"
        params:
            script         = call_dmrs_pair_script(),
            caller         = config.get('dmr_caller', 'metilene'),
            n1             = lambda wildcards: len(rep_rds_for_group(wildcards.sample1, wildcards.mc_context)),
            min_diff       = lambda wildcards, input: (
                _json.load(open(str(input.calib_json)))['min_diff']
                if _dmr_min_diff_auto
                else get_dmr_threshold('min_diff', wildcards.mc_context)
            ),
            min_cytosines  = get_dmr_threshold('min_cytosines'),
            bin_size       = get_dmr_threshold('bin_size'),
            q_value        = get_dmr_threshold('q_value'),
            min_gap        = get_dmr_threshold('min_gap'),
            min_size       = get_dmr_threshold('min_size'),
            min_reads      = get_dmr_threshold('min_reads'),
            maxdist        = get_dmr_threshold('maxdist'),
            valley         = lambda wildcards: get_dmr_threshold('valley',        wildcards.mc_context),
            maxseg         = lambda wildcards: get_dmr_threshold('maxseg',        wildcards.mc_context),
        log:
            temp(return_log_mc("{sample1}__vs__{sample2}", "DMRs", "{mc_context}"))
        conda: CONDA_ENV_MC
        shell:
            """
            {{
            printf "running %s for %s vs %s (%s)\n" "{params.caller}" "{wildcards.sample1}" "{wildcards.sample2}" "{wildcards.mc_context}"
            Rscript "{params.script}" "{threads}" "{wildcards.mc_context}" "{input.chrom_sizes}" "{output.dmrs}" "{output.counts}" \
              "{params.min_diff}" "{params.min_cytosines}" "{params.bin_size}" "{params.q_value}" "{params.min_gap}" "{params.min_size}" "{params.min_reads}" \
              "{params.maxdist}" "{params.valley}" "{params.maxseg}" \
              "{params.n1}" {input.reps1} {input.reps2}
            }} 2>&1 | tee -a "{log}"
            """

    rule aggregate_DMR_summaries:
        """Merge the per-context hyper/hypo counts into the per-pair
        summary file the downstream pipeline expects."""
        input:
            counts = lambda wildcards: [
                f"{RESULTS_DIR}/mC/DMRs/.counts__{wildcards.sample1}__vs__{wildcards.sample2}__{ctx}.tsv"
                for ctx in get_methylation_contexts()
            ]
        output:
            summary = f"{RESULTS_DIR}/mC/DMRs/summary__{{sample1}}__vs__{{sample2}}__DMRs.txt"
        localrule: True
        run:
            import pandas as pd
            sname = f"{wildcards.sample1}_vs_{wildcards.sample2}"
            merged = None
            for path in input.counts:
                df = pd.read_csv(path, sep='\t')
                merged = df if merged is None else merged.merge(df, on='Type')
            if merged is None:
                merged = pd.DataFrame({'Type': ['hyper', 'hypo']})
            merged.insert(0, 'Sample', sname)
            merged.to_csv(output.summary, sep='\t', index=False)

rule all_mc:
    input:
        final = lambda wildcards: define_final_mC_output(wildcards.ref_genome)
    output:
        touch = f"{RESULTS_DIR}/mC/chkpts/mC_analysis__{{analysis_name}}__{{ref_genome}}.done"
    localrule: True
    shell:
        """
        touch {output.touch}
        """

################################################################################
# Direct Methylation (dmC) Rules
# Handles both modBAM (direct methylation basecalls) and pre-computed bedMethyl inputs
################################################################################

MODKIT_VERSION = "0.6.1"
MODKIT_BIN = os.path.join(REPO_FOLDER, "workflow", "bin", "modkit")

rule download_modkit:
    """Download modkit binary from GitHub releases (not available via conda)."""
    output:
        binary = os.path.join(REPO_FOLDER, "workflow", "bin", "modkit")
    params:
        version = MODKIT_VERSION,
        bin_dir = os.path.join(REPO_FOLDER, "workflow", "bin")
    log:
        temp(os.path.join(REPO_FOLDER, RESULTS_DIR, "logs", "download_modkit.log"))
    shell:
        """
        {{
        mkdir -p {params.bin_dir}

        # Download modkit release (u16 = Ubuntu 16 build, compatible with CentOS 7+)
        MODKIT_URL="https://github.com/nanoporetech/modkit/releases/download/v{params.version}/modkit_v{params.version}_u16_x86_64.tar.gz"
        printf "Downloading modkit v{params.version} from $MODKIT_URL\n"

        curl -fSL "$MODKIT_URL" -o /tmp/modkit.tar.gz
        tar -xzf /tmp/modkit.tar.gz -C {params.bin_dir}
        rm -f /tmp/modkit.tar.gz

        # Move modkit binary from extracted subdirectory to bin/
        # The tarball extracts to dist_modkit_v<version>_<hash>/modkit
        extracted_modkit=$(find {params.bin_dir} -name "modkit" -type f | head -1)
        if [[ -n "$extracted_modkit" && "$extracted_modkit" != "{output.binary}" ]]; then
            mv "$extracted_modkit" {output.binary}
            # Cleanup extracted directory
            rm -rf {params.bin_dir}/dist_modkit_*
        fi

        # Make executable
        chmod +x {output.binary}

        # Verify installation
        {output.binary} --version

        printf "modkit installed successfully at {output.binary}\n"
        }} > {log} 2>&1
        """

rule get_dmc_input:
    """Acquire and validate a direct methylation input file (modBAM or bedMethyl).

    Automatically detects the input type based on file extension and content:
    - modBAM: BAM files with MM/ML methylation tags
    - bedMethyl: pre-computed methylation calls in BED format

    Supports both:
    - Direct file paths: /path/to/sample.bam
    - Directory paths: /path/to/dir (finds file matching seq_id)

    When searching directories, bedMethyl files are preferred over modBAM if both
    exist with the same seq_id prefix (bedMethyl is pre-computed and ready to use).
    A warning is logged when a modBAM is skipped in favor of bedMethyl.

    Creates a marker file indicating the detected type for downstream rules.
    """
    input:
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt",
        # validated holds the actual data for URL inputs (downloaded directly)
        # and a symlink to the user's data for local-file/directory inputs.
        # Either way it's marked as a candidate for cleanup once downstream
        # consumers (prepare_modbam_for_pileup or copy_bedmethyl_input) are
        # done, mirroring how raw__*.fastq.gz works in sample_download.smk.
        # Set keep_dmc_inputs=True in the options to retain across runs.
        validated = maybe_temp(f"{RESULTS_DIR}/mC/dmc/validated__{{sample_name}}.input", config.get('keep_dmc_inputs', False))
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        dmc_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'fastq_path'),
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, 'seq_id'),
        validate_script = os.path.join(REPO_FOLDER,"workflow","scripts","validate_dmc_input.py")
    log:
        temp(return_log_mc("{sample_name}", "get_dmc_input", "dmC"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        printf "\nDetecting and validating dmC input for {params.sample_name}\n"

        # Resolve input path: can be a file or a directory
        dmc_path="{params.dmc_path}"
        seq_id="{params.seq_id}"

        if [[ "$dmc_path" == http://* || "$dmc_path" == https://* ]]; then
            # URL provided — download directly to {output.validated}. The
            # `.input` extension hides modBAM/bedMethyl, but downstream type
            # detection in validate_dmc_input.py uses BAM magic bytes and
            # gzip+text content checks, not extension matching, so the loss
            # of the original suffix is fine.
            printf "Downloading dmC input from URL: $dmc_path\n"
            mkdir -p "$(dirname {output.validated})"
            curl --fail --show-error --location --max-redirs 5 \
                 --retry 3 --connect-timeout 30 --max-time 7200 \
                 --proto '=https,http' -o "{output.validated}" "$dmc_path"
            input_file="{output.validated}"
            printf "Downloaded to: $input_file\n"
        elif [[ -f "$dmc_path" ]]; then
            # Direct file path provided
            input_file="$dmc_path"
            printf "Using direct file path: $input_file\n"
        elif [[ -d "$dmc_path" ]]; then
            # Directory path provided - find file matching seq_id
            # Uses *seq_id* pattern consistent with sample_download.smk
            # Prefers bedMethyl over modBAM if both exist (bedMethyl is pre-computed)
            printf "Searching directory for seq_id '$seq_id'...\n"

            # Check for bedMethyl files first (preferred)
            # Note: Use {{ || true; }} pattern to handle pipefail when ls finds no matches
            bedmethyl_file=""
            for ext in .bed.gz .bedmethyl .bed; do
                match_count=$( {{ ls -1 "$dmc_path"/*"$seq_id"*"$ext" 2>/dev/null || true; }} | wc -l)
                if [[ "$match_count" -eq 1 ]]; then
                    bedmethyl_file=$(ls "$dmc_path"/*"$seq_id"*"$ext")
                    break
                elif [[ "$match_count" -gt 1 ]]; then
                    printf "Error: Multiple bedMethyl files found matching seq_id '$seq_id' with extension '$ext':\n"
                    ls -1 "$dmc_path"/*"$seq_id"*"$ext"
                    printf "Please use a more specific seq_id or provide a direct file path.\n"
                    exit 1
                fi
            done

            # Check for modBAM files
            # Note: Use {{ || true; }} pattern to handle pipefail when ls finds no matches
            modbam_file=""
            modbam_count=$( {{ ls -1 "$dmc_path"/*"$seq_id"*.bam 2>/dev/null || true; }} | wc -l)
            if [[ "$modbam_count" -eq 1 ]]; then
                modbam_file=$(ls "$dmc_path"/*"$seq_id"*.bam)
            elif [[ "$modbam_count" -gt 1 ]]; then
                printf "Error: Multiple modBAM files found matching seq_id '$seq_id':\n"
                ls -1 "$dmc_path"/*"$seq_id"*.bam
                printf "Please use a more specific seq_id or provide a direct file path.\n"
                exit 1
            fi

            # Select input file: prefer bedMethyl over modBAM
            if [[ -n "$bedmethyl_file" ]]; then
                input_file="$bedmethyl_file"
                printf "Found bedMethyl: $input_file\n"
                if [[ -n "$modbam_file" ]]; then
                    printf "WARNING: modBAM also found ($modbam_file) but using bedMethyl (pre-computed calls preferred)\n"
                fi
            elif [[ -n "$modbam_file" ]]; then
                input_file="$modbam_file"
                printf "Found modBAM: $input_file\n"
            else
                printf "Error: No dmC file found matching seq_id '$seq_id' in '$dmc_path'\n"
                printf "Looked for patterns: *$seq_id*.bed.gz, *$seq_id*.bedmethyl, *$seq_id*.bed, *$seq_id*.bam\n"
                printf "Available files:\n"
                ls -la "$dmc_path"
                exit 1
            fi
        else
            printf "Error: Path does not exist: $dmc_path\n"
            exit 1
        fi

        # Auto-detect input type
        input_type=$(python {params.validate_script} detect "$input_file")
        printf "Detected input type: $input_type\n"

        # Write type marker for downstream rules
        echo "$input_type" > {output.type_marker}

        # Validate based on detected type
        python {params.validate_script} "$input_type" "$input_file" {input.chrom_sizes}

        # For non-URL paths, link the user's data into the validated location
        # (zero-copy). For URL paths the download already wrote there.
        if [[ "$input_file" != "{output.validated}" ]]; then
            ln -sf "$(realpath "$input_file")" "{output.validated}"
        fi

        printf "\nInput validated successfully\n"
        }} 2>&1 | tee -a "{log}"
        """

def get_dmc_input_type(wildcards):
    """Get the input type (modBAM or bedMethyl) for a dmC sample by reading the marker file."""
    sname = wildcards.sample_name
    marker_file = f"{RESULTS_DIR}/mC/dmc/input_type__{sname}.txt"
    # This function is called during DAG building, marker file may not exist yet
    # Return a checkpoint-compatible path
    return marker_file

checkpoint dmc_input_checkpoint:
    """Checkpoint to determine dmC input type for branching workflow."""
    input:
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt"
    output:
        touch = touch(f"{RESULTS_DIR}/mC/dmc/checkpoint__{{sample_name}}.done")
    wildcard_constraints:
        sample_name = _DMC_WC
    localrule: True

def get_pileup_input_for_dmc(wildcards):
    """Determine pileup input based on detected input type."""
    checkpoint_output = checkpoints.dmc_input_checkpoint.get(**wildcards).output[0]
    sname = wildcards.sample_name
    type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{sname}.txt"
    with open(type_marker) as f:
        input_type = f.read().strip()
    if input_type == "modBAM":
        return f"{RESULTS_DIR}/mC/dmc/pileup_modbam__{sname}.bedmethyl.gz"
    else:
        return f"{RESULTS_DIR}/mC/dmc/pileup_bedmethyl__{sname}.bedmethyl.gz"

rule prepare_modbam_for_pileup:
    """Prepare modBAM input: index and optionally realign."""
    input:
        validated = f"{RESULTS_DIR}/mC/dmc/validated__{{sample_name}}.input",
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt",
        fasta = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        aligned_bam = maybe_temp(f"{RESULTS_DIR}/mC/dmc/aligned__{{sample_name}}.bam", config.get('keep_dmc_intermediates', False)),
        aligned_bai = maybe_temp(f"{RESULTS_DIR}/mC/dmc/aligned__{{sample_name}}.bam.bai", config.get('keep_dmc_intermediates', False))
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        preset = config.get('dmc_methylation', {}).get('alignment', {}).get('preset', 'lr:hqae')
    log:
        temp(return_log_mc("{sample_name}", "prepare_modbam", "dmC"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        # Check if input is modBAM
        input_type=$(cat {input.type_marker})
        if [[ "$input_type" != "modBAM" ]]; then
            printf "Skipping - input is not modBAM (is $input_type)\n"
            touch {output.aligned_bam} {output.aligned_bai}
            exit 0
        fi

        printf "\nChecking alignment status for {params.sample_name}\n"

        # Check if BAM is aligned by looking for @SQ headers
        has_sq=$(samtools view -H {input.validated} | grep -c "^@SQ" || true)

        # Check if aligned to correct reference
        needs_realign=false

        if [[ "$has_sq" -eq 0 ]]; then
            printf "BAM is unaligned, will align to reference\n"
            needs_realign=true
        else
            # Check chromosome overlap with reference
            bam_chroms=$(samtools view -H {input.validated} | grep "^@SQ" | cut -f2 | sed 's/SN://' | sort | head -20)
            ref_chroms=$(cut -f1 {input.chrom_sizes} | sort | head -20)
            n_ref_chroms=$(wc -l < {input.chrom_sizes})
            overlap=$(comm -12 <(echo "$bam_chroms") <(echo "$ref_chroms") | wc -l)
            min_overlap=5
            if [[ "$n_ref_chroms" -lt "$min_overlap" ]]; then
                min_overlap=$n_ref_chroms
            fi
            if [[ "$overlap" -lt "$min_overlap" ]]; then
                printf "Low chromosome overlap ($overlap/$n_ref_chroms) with reference, will realign\n"
                needs_realign=true
            fi
        fi

        if [[ "$needs_realign" == "true" ]]; then
            ref_genome=$(echo {input.fasta} | sed 's|.*/||; s|\\.fa$||')
            printf "Aligning modBAM to $ref_genome with mm2plus\n"
            samtools fastq -T MM,ML {input.validated} | \
                mm2plus -ax {params.preset} -t {threads} -y {input.fasta} - | \
                samtools sort -@ {threads} -T "$TMPDIR/align.sort" -o {output.aligned_bam} -
            samtools index -@ {threads} {output.aligned_bam}
        else
            printf "BAM is already aligned to compatible reference, claiming\n"
            # Move (don't symlink) the input into this rule's output location
            # so aligned_bam is a single-owner file. For URL-downloaded data,
            # this is a rename of the real BAM. For local-file/directory
            # inputs where validated is a symlink, mv renames the symlink
            # itself (POSIX default) — the user's source file is untouched.
            mv {input.validated} {output.aligned_bam}
            samtools index -@ {threads} {output.aligned_bam}
        fi

        printf "\nAlignment complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule modkit_pileup_dmc:
    """Generate bedMethyl file from aligned modBAM using modkit pileup."""
    input:
        bam = f"{RESULTS_DIR}/mC/dmc/aligned__{{sample_name}}.bam",
        bai = f"{RESULTS_DIR}/mC/dmc/aligned__{{sample_name}}.bam.bai",
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt",
        fasta = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa",
        modkit = MODKIT_BIN
    output:
        bedmethyl = maybe_temp(f"{RESULTS_DIR}/mC/dmc/pileup_modbam__{{sample_name}}.bedmethyl.gz", config.get('keep_dmc_intermediates', False))
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        combine_mods = "--combine-mods" if config.get('dmc_methylation', {}).get('pileup', {}).get('combine_mods', True) else ""
    log:
        temp(return_log_mc("{sample_name}", "modkit_pileup", "dmC"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        printf "\nRunning modkit pileup for {params.sample_name}\n"

        {input.modkit} pileup \
            --threads {threads} \
            --ref {input.fasta} \
            {params.combine_mods} \
            --modified-bases C \
            --filter-threshold 0.75 \
            {input.bam} \
            /dev/stdout | pigz -p {threads} > {output.bedmethyl}

        printf "\nPileup complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule copy_bedmethyl_input:
    """Copy pre-computed bedMethyl to pileup location for consistent downstream processing."""
    input:
        validated = f"{RESULTS_DIR}/mC/dmc/validated__{{sample_name}}.input",
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt"
    output:
        bedmethyl = maybe_temp(f"{RESULTS_DIR}/mC/dmc/pileup_bedmethyl__{{sample_name}}.bedmethyl.gz", config.get('keep_dmc_intermediates', False))
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        temp(return_log_mc("{sample_name}", "copy_bedmethyl", "dmC"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        printf "\nCopying bedMethyl input for {params.sample_name}\n"

        # Copy/compress the validated file
        # Check actual file content, not path extension (validated input path
        # is always .input but may be a symlink to a .gz file)
        if file -bL {input.validated} | grep -q gzip; then
            cp {input.validated} {output.bedmethyl}
        else
            pigz -p {threads} -c {input.validated} > {output.bedmethyl}
        fi

        printf "\nbedMethyl copied\n"
        }} 2>&1 | tee -a "{log}"
        """

rule merge_pileup_sources:
    """Create unified pileup path regardless of input type."""
    input:
        pileup = get_pileup_input_for_dmc
    output:
        bedmethyl = f"{RESULTS_DIR}/mC/dmc/pileup__{{sample_name}}.bedmethyl.gz"
    wildcard_constraints:
        sample_name = _DMC_WC
    localrule: True
    shell:
        """
        mv {input.pileup} {output.bedmethyl}
        """

rule modkit_summary_dmc:
    """Generate QC statistics from modBAM using modkit summary."""
    input:
        bam = f"{RESULTS_DIR}/mC/dmc/aligned__{{sample_name}}.bam",
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt",
        modkit = MODKIT_BIN
    output:
        summary = f"{RESULTS_DIR}/mC/dmc/summary__{{sample_name}}.txt"
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        sample_name = lambda wildcards: wildcards.sample_name
    log:
        temp(return_log_mc("{sample_name}", "modkit_summary", "dmC"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        input_type=$(cat {input.type_marker})

        if [[ "$input_type" == "modBAM" ]]; then
            printf "\nGenerating modkit summary for {params.sample_name}\n"
            {input.modkit} summary --threads {threads} {input.bam} > {output.summary}
        else
            printf "\nGenerating summary for pre-computed bedMethyl {params.sample_name}\n"
            printf "Sample: {params.sample_name}\n" > {output.summary}
            printf "Input type: pre-computed bedMethyl\n" >> {output.summary}
            printf "Note: Limited statistics available for bedMethyl input\n" >> {output.summary}
        fi

        printf "\nSummary complete\n"
        }} 2>&1 | tee -a "{log}"
        """

rule make_mc_stats_dmc:
    """Generate mapping stats for dmC samples in Bismark-compatible format.

    Produces stats file compatible with prepping_mapping_stats rule for combined reports.
    Uses CX_report file (unified format) for coverage statistics.
    """
    input:
        cx_report = f"{RESULTS_DIR}/mC/dmc/cx_report__{{sample_name}}.CX_report.txt.gz",
        type_marker = f"{RESULTS_DIR}/mC/dmc/input_type__{{sample_name}}.txt",
        chrom_sizes = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/chrom.sizes"
    output:
        stat_file = f"{RESULTS_DIR}/mC/reports/summary_mC_SE_mapping_stats_{{sample_name}}.txt"
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        group = lambda wildcards: parse_sample_name(wildcards.sample_name)['line'],
        levels = lambda wildcards: parse_sample_name(wildcards.sample_name)['levels_label'],
        sample_type = "dmC",
        replicate = lambda wildcards: parse_sample_name(wildcards.sample_name)['replicate'],
        ref_genome = lambda wildcards: parse_sample_name(wildcards.sample_name)['ref_genome']
    log:
        temp(return_log_mc("{sample_name}", "making_stats", "dmC"))
    conda: CONDA_ENV_MC
    shell:
        """
        {{
        printf "\nMaking mapping statistics summary for dmC sample {params.sample_name}\n"

        input_type=$(cat {input.type_marker})

        # Get genome size and coverage stats from CX_report (CG sites only for consistency)
        genome_size=$(awk '{{sum+=$2}} END {{print sum}}' {input.chrom_sizes})
        cov_stats=$(zcat {input.cx_report} | awk -v gs="$genome_size" '
            BEGIN {{total_cov=0; n_sites=0; n_sites_3x=0}}
            $6 == "CG" {{
                cov = $4 + $5;
                total_cov += cov;
                n_sites += 1;
                if (cov >= 3) n_sites_3x += 1;
            }}
            END {{
                if (n_sites > 0) {{
                    pct_cov = n_sites / gs * 100;
                    pct_cov_3x = n_sites_3x / gs * 100;
                    avg_cov_all = total_cov / gs;
                    avg_cov_covered = total_cov / n_sites;
                }} else {{
                    pct_cov = 0; pct_cov_3x = 0; avg_cov_all = 0; avg_cov_covered = 0;
                }}
                printf "%.4f\\t%.4f\\t%.4f\\t%.4f", pct_cov, pct_cov_3x, avg_cov_all, avg_cov_covered;
            }}
        ')

        pct_cov=$(echo "$cov_stats" | cut -f1)
        pct_cov_3x=$(echo "$cov_stats" | cut -f2)
        avg_cov_all=$(echo "$cov_stats" | cut -f3)
        avg_cov_covered=$(echo "$cov_stats" | cut -f4)

        # Write header
        printf "Group\\tLevels\\tSample\\tRep\\tReference_genome\\tTotal_reads\\tPassing_filtering\\tAll_mapped_reads\\tUniquely_mapped_reads\\tPercentage_covered\\tPercentage_covered_min3reads\\tAverage_coverage_all\\tAverage_coverage_covered\\tNon_conversion_rate(Pt/Lambda)\\n" > {output.stat_file}

        # For modBAM input, we can get read counts from the aligned BAM
        if [[ "$input_type" == "modBAM" ]]; then
            bam_file="{config[output_dir]}/mC/dmc/aligned__{params.sample_name}.bam"
            if [[ -f "$bam_file" ]]; then
                flagstat=$(samtools flagstat "$bam_file")
                tot=$(echo "$flagstat" | grep "in total" | awk '{{print $1}}')
                mapped=$(echo "$flagstat" | grep "mapped (" | head -1 | awk '{{print $1}}')
                uniq=$(samtools view -c -q 20 "$bam_file")
                if [ "$tot" -gt 0 ]; then
                    pct_mapped=$(awk "BEGIN {{printf \\"%.2f\\", $mapped/$tot*100}}")
                    pct_uniq=$(awk "BEGIN {{printf \\"%.2f\\", $uniq/$tot*100}}")
                else
                    pct_mapped="0.00"
                    pct_uniq="0.00"
                fi
                printf "{params.group}\\t{params.levels}\\t{params.sample_type}\\t{params.replicate}\\t{params.ref_genome}\\t$tot\\t$tot (${{pct_mapped}}%%)\\t$mapped (${{pct_mapped}}%%)\\t$uniq (${{pct_uniq}}%%)\\t$pct_cov\\t$pct_cov_3x\\t$avg_cov_all\\t$avg_cov_covered\\tNA\\n" >> {output.stat_file}
            else
                printf "{params.group}\\t{params.levels}\\t{params.sample_type}\\t{params.replicate}\\t{params.ref_genome}\\tNA\\tNA\\tNA\\tNA\\t$pct_cov\\t$pct_cov_3x\\t$avg_cov_all\\t$avg_cov_covered\\tNA\\n" >> {output.stat_file}
            fi
        else
            # For bedMethyl input, no BAM stats available
            printf "{params.group}\\t{params.levels}\\t{params.sample_type}\\t{params.replicate}\\t{params.ref_genome}\\tNA\\tNA\\tNA\\tNA\\t$pct_cov\\t$pct_cov_3x\\t$avg_cov_all\\t$avg_cov_covered\\tNA\\n" >> {output.stat_file}
        fi

        printf "\ndmC stats complete\\n"
        }} 2>&1 | tee -a "{log}"
        """

rule convert_bedmethyl_to_cx_report:
    """Convert bedMethyl format to Bismark CX_report format.

    Determines context (CG/CHG/CHH) from the reference genome and outputs
    a single CX_report file matching the Bismark output format. This enables
    unified downstream processing with merging_mc_replicates, make_mc_bigwig_files,
    and call_DMRs_pairwise.

    When methylation_contexts is a strict subset of {CG, CHG, CHH}, the
    output is filtered to only those contexts.
    """
    input:
        bedmethyl = f"{RESULTS_DIR}/mC/dmc/pileup__{{sample_name}}.bedmethyl.gz",
        fasta = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa",
        fai = lambda wildcards: f"{GENOMES_DIR}/{parse_sample_name(wildcards.sample_name)['ref_genome']}/{parse_sample_name(wildcards.sample_name)['ref_genome']}.fa.fai"
    output:
        cx_report = maybe_temp(f"{RESULTS_DIR}/mC/dmc/cx_report__{{sample_name}}.CX_report.txt.gz", config.get('keep_cx_reports', False))
    wildcard_constraints:
        sample_name = _DMC_WC
    params:
        script = os.path.join(REPO_FOLDER, "workflow", "scripts", "bedmethyl_to_cx_report.py"),
        sample_name = lambda wildcards: wildcards.sample_name,
        # Pipe-separated regex alternatives ("CG|CHG|CHH" or "CG", etc.) for
        # the post-processing context filter. When all three are active no
        # filtering is needed.
        context_filter = "|".join(get_methylation_contexts())
    log:
        temp(return_log_mc("{sample_name}", "bedmethyl_to_cx", "dmC"))
    conda: CONDA_ENV_MC
    threads: 1
    shell:
        """
        {{
        printf "Converting bedMethyl to CX_report format for {params.sample_name}...\n"
        printf "Active methylation contexts: {params.context_filter}\n"

        # Convert bedMethyl to CX_report (context determined from reference)
        python {params.script} {input.bedmethyl} {input.fasta} /dev/stdout > {config[output_dir]}/mC/dmc/tmp__{params.sample_name}.cx

        # Filter to active contexts if the user requested a strict subset.
        # The full set is "CG|CHG|CHH" — when that's what we got, skip the
        # filter pass.
        if [[ "{params.context_filter}" != "CG|CHG|CHH" ]]; then
            printf "Filtering CX_report to contexts: {params.context_filter}\n"
            awk -F'\\t' -v p='^({params.context_filter})$' '$6 ~ p' {config[output_dir]}/mC/dmc/tmp__{params.sample_name}.cx > {config[output_dir]}/mC/dmc/tmp__{params.sample_name}_filtered.cx
            mv {config[output_dir]}/mC/dmc/tmp__{params.sample_name}_filtered.cx {config[output_dir]}/mC/dmc/tmp__{params.sample_name}.cx
        fi

        # Sort by chromosome and position, then compress
        printf "Sorting and compressing...\n"
        sort -k1,1 -k2,2n {config[output_dir]}/mC/dmc/tmp__{params.sample_name}.cx | pigz -p {threads} > {output.cx_report}
        rm -f {config[output_dir]}/mC/dmc/tmp__{params.sample_name}.cx

        printf "Conversion complete\n"
        }} 2>&1 | tee -a "{log}"
        """
