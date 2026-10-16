#!/usr/bin/env bash
# subset_test_data.sh — SLURM-based chromosome/region subsetting for test data prep
#
# Self-resubmitting controller that orchestrates three phases:
#   Phase 1 (index):    Build genome indexes (bowtie2, STAR, bismark)
#   Phase 2 (sample):   Spawn per-sample SLURM jobs (download, align, subset)
#   Phase 3 (gather):   Controller resubmits itself to wait and summarize
#
# Usage:
#   bash scripts/subset_test_data.sh <manifest.tsv> <target_region> [options]
#
# Manifest TSV columns (tab-separated, header required):
#   sample_id    assay    srr_accessions    read_layout    reference_fasta
#
#   sample_id:        Unique name for the output files
#   assay:            One of: ChIP, RNA, RAMPAGE, sRNA, WGBS, dmC
#   srr_accessions:   SRA run(s), +-separated for merging (e.g. SRR123+SRR456)
#                     OR a URL/path to a pre-aligned BAM/CRAM (s3://, https://, local)
#                     OR a local path / s3:// URI to BAM/modBAM/bedmethyl (for dmC)
#                     When a .bam/.cram URL is given for any assay, the script
#                     streams the file and extracts reads in the target region
#                     directly — no local alignment needed.
#   read_layout:      SE or PE
#   reference_fasta:  Path to full reference genome FASTA (indexed)
#
# Target region: a chromosome name (e.g. chr21) or samtools region string
#   (e.g. chr21:1-10000000). Must match sequence names in the reference.
#
# Options:
#   --outdir DIR      Output directory (default: test-data-prep/subset)
#   --keep-aligned    Keep intermediate full-genome BAMs (default: delete)
#   --dry-run         Print sbatch commands without submitting
#   --local           Run everything locally (no SLURM), serial execution
#   --clean           Remove intermediate files (downloads/, merged/, jobs/,
#                     status/) after verifying all samples completed. Keeps
#                     final output FASTQs/BAMs and logs.
#
# Output per sample:
#   {outdir}/{sample_id}_R0.fastq.gz              (SE)
#   {outdir}/{sample_id}_R1.fastq.gz + _R2.fastq.gz  (PE)
#   {outdir}/{sample_id}.bam                      (dmC modBAM)
#   {outdir}/{sample_id}.bedmethyl.gz             (dmC bedmethyl)
#
# Requirements: samtools, bowtie2, STAR, bismark, wget, fasterq-dump (sra-tools),
#               pigz, awscli (for s3:// dmC sources)
# Only the tools needed for the assay types in the manifest must be installed.

set -euo pipefail

# ---------------------------------------------------------------------------
# SLURM parameters
# ---------------------------------------------------------------------------
# Partition and QoS are site-specific; set SUBSET_PARTITION / SUBSET_QOS to
# pass them, or leave them unset to use the cluster's defaults.
PARTITION="${SUBSET_PARTITION:-}"
QOS="${SUBSET_QOS:-}"
SITE_FLAGS=()
[[ -n "$PARTITION" ]] && SITE_FLAGS+=(--partition="$PARTITION")
[[ -n "$QOS" ]] && SITE_FLAGS+=(--qos="$QOS")
ALIGN_THREADS=24           # alignment + samtools sort threads
ALIGN_MEM="80G"            # per alignment job (~12GB per bismark --parallel instance)
ALIGN_TMP="500G"
ALIGN_TIME="24:00:00"
INDEX_THREADS=16
INDEX_MEM="48G"            # STAR index needs ~32G for human
INDEX_TMP="96G"
INDEX_TIME="4:00:00"
DOWNLOAD_THREADS=8
DOWNLOAD_MEM="16G"
DOWNLOAD_TMP="48G"
DOWNLOAD_TIME="12:00:00"
CTRL_MEM="4G"              # controller resubmission
CTRL_TIME="2:00:00"
CONDA_ENV="epicc"          # conda env with aligners + samtools

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
OUTDIR="test-data-prep/subset"
KEEP_ALIGNED=false
DRY_RUN=false
LOCAL_MODE=false
CLEAN_MODE=false

# Internal: phase tracking (set by self-resubmission)
_PHASE="${_SUBSET_PHASE:-auto}"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
usage() {
    sed -n '2,/^[^#]/p' "$0" | head -n -1 | sed 's/^# \?//'
    exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --keep-aligned) KEEP_ALIGNED=true; shift ;;
        --dry-run)      DRY_RUN=true; shift ;;
        --local)        LOCAL_MODE=true; shift ;;
        --clean)        CLEAN_MODE=true; shift ;;
        -h|--help)      usage ;;
        *)              POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]}"

if [[ ${#POSITIONAL[@]} -lt 2 ]]; then
    echo "Error: requires <manifest.tsv> and <target_region>"
    usage
fi

MANIFEST="$(realpath "$1")"
TARGET_REGION="$2"
SCRIPT_PATH="$(realpath "$0")"
OUTDIR="$(realpath -m "$OUTDIR")"

if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: manifest file not found: $MANIFEST"
    exit 1
fi

# Directories
DL_DIR="${OUTDIR}/downloads"
MERGE_DIR="${OUTDIR}/merged"
LOGS_DIR="${OUTDIR}/logs"
JOBS_DIR="${OUTDIR}/jobs"    # per-sample job scripts
STATUS_DIR="${OUTDIR}/status" # completion markers
mkdir -p "$OUTDIR" "$DL_DIR" "$MERGE_DIR" "$LOGS_DIR" "$JOBS_DIR" "$STATUS_DIR"

# ---------------------------------------------------------------------------
# Logging — always to stderr so stdout can carry job IDs in subshells
# ---------------------------------------------------------------------------
log() { printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$1" >&2; }

# ---------------------------------------------------------------------------
# Helper: detect pre-aligned BAM/CRAM input
# ---------------------------------------------------------------------------
is_bam_or_cram() { [[ "$1" == *.bam || "$1" == *.cram ]]; }

# ---------------------------------------------------------------------------
# Conda activation preamble for generated job scripts
# ---------------------------------------------------------------------------
SLURM_TMPDIR_PREAMBLE='# Set up per-job TMPDIR (profile.d scripts do not run in non-login shells)
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    export SLURM_TMPDIR="/tmp/slurm_tmp/$SLURM_JOB_ID"
    export TMPDIR="$SLURM_TMPDIR"
    mkdir -p "$TMPDIR"
fi'

CONDA_PREAMBLE="$(cat << 'CONDAEOF'
# Activate conda environment (set +u to work around conda activate scripts
# that reference unset variables, e.g. activate-binutils_linux-64.sh ADDR2LINE)
set +u
eval "$(conda shell.bash hook)"
CONDAEOF
)
conda activate ${CONDA_ENV}
set -u"

# Combined preamble for all generated job scripts
JOB_PREAMBLE="${SLURM_TMPDIR_PREAMBLE}
${CONDA_PREAMBLE}"

# ---------------------------------------------------------------------------
# Utility: submit or run
# ---------------------------------------------------------------------------
submit_job() {
    # Usage: submit_job <job_name> <threads> <mem> <tmp> <time> <script> [dependency]
    local name="$1" threads="$2" mem="$3" tmp="$4" time="$5" script="$6"
    local dep="${7:-}"

    if $LOCAL_MODE; then
        log "  [local] Running $name ..."
        bash "$script" 2>&1 | tee -a "${LOGS_DIR}/${name}.log"
        return
    fi

    local dep_flag=""
    if [[ -n "$dep" ]]; then
        dep_flag="--dependency=afterok:${dep}"
    fi

    local cmd=(
        sbatch
        ${SITE_FLAGS[@]+"${SITE_FLAGS[@]}"}
        --job-name="subset-${name}"
        --cpus-per-task="$threads"
        --mem="$mem"
        --tmp="$tmp"
        --time="$time"
        --output="${LOGS_DIR}/${name}_%j.log"
        --error="${LOGS_DIR}/${name}_%j.log"
        --export=ALL
    )
    [[ -n "$dep_flag" ]] && cmd+=("$dep_flag")
    cmd+=("$script")

    if $DRY_RUN; then
        echo "  [dry-run] ${cmd[*]}" >&2
        echo "DRY_RUN_JOBID"
    else
        local output
        output=$("${cmd[@]}")
        local jobid
        jobid=$(echo "$output" | grep -oP '\d+$')
        log "  Submitted $name: job $jobid"
        echo "$jobid"
    fi
}

# ---------------------------------------------------------------------------
# Parse manifest into arrays
# ---------------------------------------------------------------------------
declare -a SAMPLE_IDS ASSAYS SRR_COLS LAYOUTS REF_FASTAS
while IFS=$'\t' read -r sid assay srrs layout ref; do
    [[ -z "$sid" || "$sid" == \#* || "$sid" == "sample_id" ]] && continue
    SAMPLE_IDS+=("$sid")
    ASSAYS+=("$assay")
    SRR_COLS+=("$srrs")
    LAYOUTS+=("$layout")
    # Resolve relative ref paths to absolute (jobs may run from different cwd)
    if [[ "$ref" != /* ]]; then
        ref="$(realpath -m "$ref")"
    fi
    REF_FASTAS+=("$ref")
done < "$MANIFEST"

N_SAMPLES=${#SAMPLE_IDS[@]}
log "Manifest: $MANIFEST ($N_SAMPLES samples)"
log "Target region: $TARGET_REGION"
log "Output: $OUTDIR"
log "Mode: $(if $LOCAL_MODE; then echo local; else echo SLURM; fi)"

# ---------------------------------------------------------------------------
# Determine which indexes are needed
# ---------------------------------------------------------------------------
needs_bt2=false
needs_star=false
needs_bismark=false
BT2_REF=""
STAR_REF=""
BISMARK_REF=""

for i in $(seq 0 $((N_SAMPLES - 1))); do
    ref="${REF_FASTAS[$i]}"
    case "${ASSAYS[$i]}" in
        ChIP|sRNA)   needs_bt2=true;     BT2_REF="$ref" ;;
        RNA|RAMPAGE)  needs_star=true;    STAR_REF="$ref" ;;
        WGBS|EMseq)
            if ! is_bam_or_cram "${SRR_COLS[$i]}"; then
                needs_bismark=true; BISMARK_REF="$ref"
            fi
            ;;
    esac
done

# ===========================================================================
# PHASE 1: Build indexes
# ===========================================================================
build_indexes() {
    log "=== Phase 1: Building indexes ==="
    local index_jobs=()

    if $needs_bt2; then
        local ref_base="${BT2_REF%.fa}"
        ref_base="${ref_base%.fasta}"
        if [[ ! -f "${ref_base}.1.bt2" ]]; then
            local script="${JOBS_DIR}/index_bt2.sh"
            cat > "$script" << IDXEOF
#!/usr/bin/env bash
set -euo pipefail
${JOB_PREAMBLE}
echo "Building bowtie2 index: ${BT2_REF}"
bowtie2-build --threads ${INDEX_THREADS} "${BT2_REF}" "${ref_base}"
touch "${STATUS_DIR}/index_bt2.done"
IDXEOF
            chmod +x "$script"
            local jid
            jid=$(submit_job "index_bt2" "$INDEX_THREADS" "$INDEX_MEM" "$INDEX_TMP" "$INDEX_TIME" "$script")
            index_jobs+=("$jid")
        else
            log "  bowtie2 index exists, skipping"
            touch "${STATUS_DIR}/index_bt2.done"
        fi
    fi

    if $needs_star; then
        local genome_dir
        genome_dir="$(dirname "$STAR_REF")/STAR_index"
        if [[ ! -f "${genome_dir}/SA" ]]; then
            local script="${JOBS_DIR}/index_star.sh"
            cat > "$script" << IDXEOF
#!/usr/bin/env bash
set -euo pipefail
${JOB_PREAMBLE}
echo "Building STAR index: ${STAR_REF}"
mkdir -p "${genome_dir}"
STAR --runMode genomeGenerate --runThreadN ${INDEX_THREADS} \
    --genomeDir "${genome_dir}" --genomeFastaFiles "${STAR_REF}" \
    --genomeSAindexNbases 14
touch "${STATUS_DIR}/index_star.done"
IDXEOF
            chmod +x "$script"
            local jid
            jid=$(submit_job "index_star" "$INDEX_THREADS" "$INDEX_MEM" "$INDEX_TMP" "$INDEX_TIME" "$script")
            index_jobs+=("$jid")
        else
            log "  STAR index exists, skipping"
            touch "${STATUS_DIR}/index_star.done"
        fi
    fi

    if $needs_bismark; then
        local genome_dir
        genome_dir="$(dirname "$BISMARK_REF")"
        if [[ ! -d "${genome_dir}/Bisulfite_Genome" ]]; then
            local script="${JOBS_DIR}/index_bismark.sh"
            cat > "$script" << IDXEOF
#!/usr/bin/env bash
set -euo pipefail
${JOB_PREAMBLE}
echo "Building bismark index: ${genome_dir}"
bismark_genome_preparation --parallel ${INDEX_THREADS} "${genome_dir}"
touch "${STATUS_DIR}/index_bismark.done"
IDXEOF
            chmod +x "$script"
            local jid
            jid=$(submit_job "index_bismark" "$INDEX_THREADS" "$INDEX_MEM" "$INDEX_TMP" "$INDEX_TIME" "$script")
            index_jobs+=("$jid")
        else
            log "  bismark index exists, skipping"
            touch "${STATUS_DIR}/index_bismark.done"
        fi
    fi

    # Return colon-separated job IDs for dependency chaining
    local IFS=':'
    echo "${index_jobs[*]:-}"
}

# ===========================================================================
# PHASE 2: Spawn per-sample jobs
# ===========================================================================
spawn_sample_jobs() {
    local index_dep="$1"
    log "=== Phase 2: Spawning per-sample jobs ==="

    local sample_jobs=()

    for i in $(seq 0 $((N_SAMPLES - 1))); do
        local sid="${SAMPLE_IDS[$i]}"
        local assay="${ASSAYS[$i]}"
        local srrs="${SRR_COLS[$i]}"
        local layout="${LAYOUTS[$i]}"
        local ref="${REF_FASTAS[$i]}"

        # Skip if output already exists
        if sample_complete "$sid" "$assay" "$layout"; then
            log "  $sid: output exists, skipping"
            touch "${STATUS_DIR}/${sid}.done"
            continue
        fi

        local script="${JOBS_DIR}/sample_${sid}.sh"
        write_sample_script "$sid" "$assay" "$srrs" "$layout" "$ref" "$script"

        # Determine dependency: index jobs if this assay needs one
        local dep=""
        if [[ -n "$index_dep" ]]; then
            dep="$index_dep"
        fi

        # Choose resources based on assay and input type
        local threads mem tmp time
        if is_bam_or_cram "$srrs" && [[ "$assay" != "dmC" ]]; then
            # Pre-aligned BAM/CRAM: stream + subset — no local alignment
            threads="$DOWNLOAD_THREADS"; mem="$DOWNLOAD_MEM"; tmp="$DOWNLOAD_TMP"; time="$DOWNLOAD_TIME"
            dep=""  # no index dependency
        else
            case "$assay" in
                dmC)
                    threads=4; mem="$DOWNLOAD_MEM"; tmp="$DOWNLOAD_TMP"; time="$DOWNLOAD_TIME"
                    dep=""  # no index dependency for dmC
                    ;;
                *)
                    threads="$ALIGN_THREADS"; mem="$ALIGN_MEM"; tmp="$ALIGN_TMP"; time="$ALIGN_TIME"
                    ;;
            esac
        fi

        local jid
        jid=$(submit_job "$sid" "$threads" "$mem" "$tmp" "$time" "$script" "$dep")
        sample_jobs+=("$jid")
    done

    local IFS=':'
    echo "${sample_jobs[*]:-}"
}

# ---------------------------------------------------------------------------
# Check if a sample's output already exists
# ---------------------------------------------------------------------------
sample_complete() {
    local sid="$1" assay="$2" layout="$3"
    if [[ "$assay" == "dmC" ]]; then
        [[ -f "${OUTDIR}/${sid}.bam" || -f "${OUTDIR}/${sid}.bedmethyl.gz" ]]
    elif [[ "$layout" == "PE" ]]; then
        [[ -f "${OUTDIR}/${sid}_R1.fastq.gz" && -f "${OUTDIR}/${sid}_R2.fastq.gz" ]]
    else
        [[ -f "${OUTDIR}/${sid}_R0.fastq.gz" ]]
    fi
}

# ---------------------------------------------------------------------------
# Write a self-contained per-sample job script
# ---------------------------------------------------------------------------
write_sample_script() {
    local sid="$1" assay="$2" srrs="$3" layout="$4" ref="$5" script="$6"

    cat > "$script" << SAMPLEEOF
#!/usr/bin/env bash
set -euo pipefail
${JOB_PREAMBLE}
SAMPLEEOF

    # Inject variables
    cat >> "$script" << VARSEOF
SAMPLE_ID="${sid}"
ASSAY="${assay}"
SRR_ACCESSIONS="${srrs}"
LAYOUT="${layout}"
REF_FASTA="${ref}"
TARGET_REGION="${TARGET_REGION}"
OUTDIR="${OUTDIR}"
DL_DIR="${DL_DIR}"
MERGE_DIR="${MERGE_DIR}"
STATUS_DIR="${STATUS_DIR}"
THREADS=${ALIGN_THREADS}
KEEP_ALIGNED=${KEEP_ALIGNED}
VARSEOF

    cat >> "$script" << 'BODYEOF'

log() { printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$1"; }

log "=== Processing ${SAMPLE_ID} (${ASSAY}, ${LAYOUT}) ==="

# -------------------------------------------------------------------
# dmC: download and subset modBAM or bedmethyl
# -------------------------------------------------------------------
if [[ "$ASSAY" == "dmC" ]]; then
    input_src="$SRR_ACCESSIONS"

    # Detect input type by extension
    if [[ "$input_src" == *.bedmethyl.gz || "$input_src" == *.bedmethyl ]]; then
        # --- bedmethyl path ---
        out_file="${OUTDIR}/${SAMPLE_ID}.bedmethyl.gz"
        if [[ -f "$out_file" ]]; then
            log "Output already exists, skipping."
            touch "${STATUS_DIR}/${SAMPLE_ID}.done"
            exit 0
        fi

        local_bed="$input_src"
        if [[ "$input_src" == s3://* ]]; then
            local_bed="${DL_DIR}/${SAMPLE_ID}_source.bedmethyl.gz"
            if [[ ! -f "$local_bed" ]]; then
                log "Downloading bedmethyl from S3 ..."
                aws s3 --no-sign-request cp "$input_src" "$local_bed"
            fi
        fi

        # Extract target chromosome from region (e.g. chr21 from chr21:1-1000000)
        target_chrom="${TARGET_REGION%%:*}"
        log "Subsetting bedmethyl to ${target_chrom} ..."
        zcat "$local_bed" | awk -v chr="$target_chrom" '$1 == chr' \
            | gzip > "$out_file"

        n_lines=$(zcat "$out_file" | wc -l)
        log "Subset lines: $n_lines"
        log "Done: $out_file"
    else
        # --- modBAM path ---
        out_file="${OUTDIR}/${SAMPLE_ID}.bam"
        if [[ -f "$out_file" ]]; then
            log "Output already exists, skipping."
            touch "${STATUS_DIR}/${SAMPLE_ID}.done"
            exit 0
        fi

        input_bam="$input_src"
        if [[ "$input_bam" == s3://* ]]; then
            local_bam="${DL_DIR}/${SAMPLE_ID}_source.bam"
            if [[ ! -f "$local_bam" ]]; then
                log "Downloading modBAM from S3 ..."
                aws s3 --no-sign-request cp "$input_bam" "$local_bam"
                aws s3 --no-sign-request cp "${input_bam}.bai" "${local_bam}.bai" 2>/dev/null || true
            fi
            input_bam="$local_bam"
        fi

        if [[ ! -f "${input_bam}.bai" && ! -f "${input_bam%.bam}.bai" ]]; then
            log "Indexing source BAM ..."
            samtools index -@ 4 "$input_bam"
        fi

        log "Subsetting modBAM to ${TARGET_REGION} ..."
        samtools view -b -h -@ 4 "$input_bam" "$TARGET_REGION" \
            | samtools sort -@ 4 -o "$out_file"
        samtools index "$out_file"
        log "Done: $out_file"
    fi

    touch "${STATUS_DIR}/${SAMPLE_ID}.done"
    exit 0
fi

# -------------------------------------------------------------------
# Pre-aligned BAM/CRAM: stream, filter target region, extract FASTQs
# -------------------------------------------------------------------
if [[ "$SRR_ACCESSIONS" == *.bam || "$SRR_ACCESSIONS" == *.cram ]]; then
    log "Pre-aligned BAM/CRAM detected: $SRR_ACCESSIONS"

    # Stream a BAM/CRAM from s3://, https://, or local path
    stream_bam() {
        local src="$1"
        if [[ "$src" == s3://* ]]; then
            aws s3 cp "$src" - --no-sign-request 2>/dev/null
        elif [[ "$src" == http://* || "$src" == https://* ]]; then
            curl -fsSL "$src"
        else
            cat "$src"
        fi
    }

    target_chrom="${TARGET_REGION%%:*}"
    view_threads=$(( THREADS / 2 > 0 ? THREADS / 2 - 1 : 0 ))
    sort_threads=$(( THREADS / 2 > 0 ? THREADS / 2 - 1 : 0 ))

    if [[ "$LAYOUT" == "PE" ]]; then
        log "Streaming and subsetting to ${target_chrom} → PE FASTQs ..."
        stream_bam "$SRR_ACCESSIONS" \
            | samtools view -@ "$view_threads" -h -b \
                --expr "rname == \"${target_chrom}\"" - \
            | samtools sort -@ "$sort_threads" -n -m 2G - \
            | samtools fastq -@ 1 \
                -1 "${OUTDIR}/${SAMPLE_ID}_R1.fastq.gz" \
                -2 "${OUTDIR}/${SAMPLE_ID}_R2.fastq.gz" \
                -0 /dev/null -s /dev/null -
    else
        log "Streaming and subsetting to ${target_chrom} → SE FASTQ ..."
        stream_bam "$SRR_ACCESSIONS" \
            | samtools view -@ "$view_threads" -h -b \
                --expr "rname == \"${target_chrom}\"" - \
            | samtools sort -@ "$sort_threads" -n -m 2G - \
            | samtools fastq -@ 1 \
                -0 "${OUTDIR}/${SAMPLE_ID}_R0.fastq.gz" -
    fi

    log "Done: ${SAMPLE_ID}"
    touch "${STATUS_DIR}/${SAMPLE_ID}.done"
    exit 0
fi

# -------------------------------------------------------------------
# Download SRA runs (ENA first, fasterq-dump fallback)
# -------------------------------------------------------------------
IFS='+' read -ra SRR_ARRAY <<< "$SRR_ACCESSIONS"

# Build ENA HTTPS base URL for a given accession
# Matches logic in workflow/scripts/ena_download.sh
ena_base_url() {
    local acc="$1"
    local letters="${acc%%[0-9]*}"
    local digits="${acc#"$letters"}"
    local ndigits="${#digits}"
    local first6="${acc:0:6}"
    local intermediate=""

    case "$ndigits" in
        6) intermediate="" ;;
        7) intermediate="00${digits: -1}" ;;
        8) intermediate="0${digits: -2}" ;;
        9) intermediate="${digits: -3}" ;;
        *) log "WARNING: unexpected accession format '$acc' ($ndigits digits)"; return 1 ;;
    esac

    if [[ -n "$intermediate" ]]; then
        echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${intermediate}/${acc}"
    else
        echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${acc}"
    fi
}

# Download and verify a single gzipped file; returns 0 on success
ena_fetch() {
    local url="$1" dest="$2"
    if ! wget -q -O "$dest" "$url" 2>/dev/null; then
        rm -f "$dest"; return 1
    fi
    if ! gzip -t "$dest" 2>/dev/null; then
        log "  gzip integrity check failed: $(basename "$dest")"
        rm -f "$dest"; return 1
    fi
}

# Try downloading from ENA; returns 0 on success
download_ena() {
    local srr="$1" layout="$2" dl_dir="$3"
    local base
    base=$(ena_base_url "$srr") || return 1

    if [[ "$layout" == "PE" ]]; then
        log "  Trying ENA (PE): $srr ..."
        if ! ena_fetch "${base}/${srr}_1.fastq.gz" "${dl_dir}/${srr}_1.fastq.gz" \
            || ! ena_fetch "${base}/${srr}_2.fastq.gz" "${dl_dir}/${srr}_2.fastq.gz"; then
            rm -f "${dl_dir}/${srr}_1.fastq.gz" "${dl_dir}/${srr}_2.fastq.gz"
            return 1
        fi
    else
        log "  Trying ENA (SE): $srr ..."
        if ! ena_fetch "${base}/${srr}.fastq.gz" "${dl_dir}/${srr}.fastq.gz"; then
            # Some SE datasets use _1.fastq.gz naming on ENA
            if ! ena_fetch "${base}/${srr}_1.fastq.gz" "${dl_dir}/${srr}.fastq.gz"; then
                return 1
            fi
        fi
    fi
}

# Fallback: fasterq-dump from SRA
download_sra() {
    local srr="$1" layout="$2" dl_dir="$3" threads="$4" tmpdir="$5"
    log "  Trying SRA (fasterq-dump): $srr ..."
    fasterq-dump --threads "$threads" --temp "$tmpdir" \
        --outdir "$dl_dir" --split-files --prefer-sra-lite false "$srr"
    # Compress
    if [[ "$layout" == "PE" ]]; then
        pigz -p "$threads" "${dl_dir}/${srr}_1.fastq" &
        pigz -p "$threads" "${dl_dir}/${srr}_2.fastq" &
        wait
    else
        if [[ -f "${dl_dir}/${srr}_1.fastq" && ! -f "${dl_dir}/${srr}_2.fastq" ]]; then
            mv "${dl_dir}/${srr}_1.fastq" "${dl_dir}/${srr}.fastq"
        fi
        pigz -p "$threads" "${dl_dir}/${srr}.fastq"
    fi
}

for srr in "${SRR_ARRAY[@]}"; do
    if [[ -f "${DL_DIR}/${srr}_1.fastq.gz" ]] || [[ -f "${DL_DIR}/${srr}.fastq.gz" ]]; then
        log "Already downloaded: $srr"
        continue
    fi
    log "Downloading $srr ..."
    if ! download_ena "$srr" "$LAYOUT" "$DL_DIR"; then
        log "  ENA failed, falling back to SRA"
        download_sra "$srr" "$LAYOUT" "$DL_DIR" "$THREADS" "$TMPDIR"
    fi
done

# -------------------------------------------------------------------
# Merge multiple SRR FASTQs
# -------------------------------------------------------------------
if [[ ${#SRR_ARRAY[@]} -eq 1 ]]; then
    srr="${SRR_ARRAY[0]}"
    if [[ "$LAYOUT" == "PE" ]]; then
        ln -sf "$(realpath "${DL_DIR}/${srr}_1.fastq.gz")" "${MERGE_DIR}/${SAMPLE_ID}_R1.fastq.gz"
        ln -sf "$(realpath "${DL_DIR}/${srr}_2.fastq.gz")" "${MERGE_DIR}/${SAMPLE_ID}_R2.fastq.gz"
    else
        ln -sf "$(realpath "${DL_DIR}/${srr}.fastq.gz")" "${MERGE_DIR}/${SAMPLE_ID}_R0.fastq.gz"
    fi
else
    log "Merging ${#SRR_ARRAY[@]} runs ..."
    if [[ "$LAYOUT" == "PE" ]]; then
        fqs_r1=()
        fqs_r2=()
        for srr in "${SRR_ARRAY[@]}"; do
            fqs_r1+=("${DL_DIR}/${srr}_1.fastq.gz")
            fqs_r2+=("${DL_DIR}/${srr}_2.fastq.gz")
        done
        cat "${fqs_r1[@]}" > "${MERGE_DIR}/${SAMPLE_ID}_R1.fastq.gz"
        cat "${fqs_r2[@]}" > "${MERGE_DIR}/${SAMPLE_ID}_R2.fastq.gz"
    else
        fqs=()
        for srr in "${SRR_ARRAY[@]}"; do fqs+=("${DL_DIR}/${srr}.fastq.gz"); done
        cat "${fqs[@]}" > "${MERGE_DIR}/${SAMPLE_ID}_R0.fastq.gz"
    fi
fi

# -------------------------------------------------------------------
# Align to full genome
# -------------------------------------------------------------------
bam_full="${OUTDIR}/tmp_${SAMPLE_ID}_full.bam"
sort_mem="$(( ${THREADS} > 0 ? 1536 : 768 ))M"

case "$ASSAY" in
    ChIP|sRNA)
        ref_base="${REF_FASTA%.fa}"
        ref_base="${ref_base%.fasta}"
        if [[ "$LAYOUT" == "PE" ]]; then
            log "Aligning (PE, bowtie2) ..."
            bowtie2 -p "$THREADS" -x "$ref_base" \
                -1 "${MERGE_DIR}/${SAMPLE_ID}_R1.fastq.gz" \
                -2 "${MERGE_DIR}/${SAMPLE_ID}_R2.fastq.gz" \
                2>"${OUTDIR}/tmp_${SAMPLE_ID}_bt2.log" \
                | samtools view -bh -F 4 -@ 2 \
                | samtools sort -@ "$THREADS" -m "$sort_mem" -o "$bam_full"
        else
            log "Aligning (SE, bowtie2) ..."
            bowtie2 -p "$THREADS" -x "$ref_base" \
                -U "${MERGE_DIR}/${SAMPLE_ID}_R0.fastq.gz" \
                2>"${OUTDIR}/tmp_${SAMPLE_ID}_bt2.log" \
                | samtools view -bh -F 4 -@ 2 \
                | samtools sort -@ "$THREADS" -m "$sort_mem" -o "$bam_full"
        fi
        cat "${OUTDIR}/tmp_${SAMPLE_ID}_bt2.log"
        rm -f "${OUTDIR}/tmp_${SAMPLE_ID}_bt2.log"
        ;;

    RNA|RAMPAGE)
        genome_dir="$(dirname "$REF_FASTA")/STAR_index"
        star_tmp="${OUTDIR}/_STARtmp_${SAMPLE_ID}"
        if [[ "$LAYOUT" == "PE" ]]; then
            log "Aligning (PE, STAR) ..."
            STAR --runThreadN "$THREADS" --genomeDir "$genome_dir" \
                --readFilesIn "${MERGE_DIR}/${SAMPLE_ID}_R1.fastq.gz" "${MERGE_DIR}/${SAMPLE_ID}_R2.fastq.gz" \
                --readFilesCommand zcat --outSAMtype BAM Unsorted \
                --outTmpDir "$star_tmp" \
                --outFileNamePrefix "${OUTDIR}/tmp_${SAMPLE_ID}_" \
                --outStd BAM_Unsorted 2>/dev/null \
                | samtools sort -@ "$THREADS" -m "$sort_mem" -o "$bam_full"
        else
            log "Aligning (SE, STAR) ..."
            STAR --runThreadN "$THREADS" --genomeDir "$genome_dir" \
                --readFilesIn "${MERGE_DIR}/${SAMPLE_ID}_R0.fastq.gz" \
                --readFilesCommand zcat --outSAMtype BAM Unsorted \
                --outTmpDir "$star_tmp" \
                --outFileNamePrefix "${OUTDIR}/tmp_${SAMPLE_ID}_" \
                --outStd BAM_Unsorted 2>/dev/null \
                | samtools sort -@ "$THREADS" -m "$sort_mem" -o "$bam_full"
        fi
        rm -rf "${star_tmp}" "${OUTDIR}/tmp_${SAMPLE_ID}_"Log* "${OUTDIR}/tmp_${SAMPLE_ID}_"SJ*
        ;;

    WGBS|EMseq)
        genome_dir="$(dirname "$REF_FASTA")"
        # bismark --parallel (--multicore) is incompatible with --basename,
        # so we use default output naming: {input}_bismark_bt2[_pe].bam
        # Each --parallel instance uses ~4 cores (bowtie2 -p 1 default + bismark
        # overhead) and ~10-12GB RAM for human genome. See:
        # https://github.com/FelixKrueger/Bismark/issues/96
        bm_parallel=$(( THREADS / 4 ))
        [[ $bm_parallel -lt 1 ]] && bm_parallel=1
        [[ $bm_parallel -gt 8 ]] && bm_parallel=8
        if [[ "$LAYOUT" == "PE" ]]; then
            log "Aligning (PE, bismark, --parallel $bm_parallel) ..."
            bismark --parallel "$bm_parallel" --genome "$genome_dir" --temp_dir "$TMPDIR" \
                -1 "${MERGE_DIR}/${SAMPLE_ID}_R1.fastq.gz" \
                -2 "${MERGE_DIR}/${SAMPLE_ID}_R2.fastq.gz" \
                --output_dir "${OUTDIR}"
            # Default PE output: {SAMPLE_ID}_R1_bismark_bt2_pe.bam
            bm_bam="${OUTDIR}/${SAMPLE_ID}_R1_bismark_bt2_pe.bam"
            samtools sort -@ "$THREADS" -m "$sort_mem" -o "$bam_full" "$bm_bam"
            rm -f "$bm_bam" "${OUTDIR}/${SAMPLE_ID}_R1_bismark_bt2"*report* "${OUTDIR}/${SAMPLE_ID}_R1_bismark_bt2"*nucleotide*
        else
            log "Aligning (SE, bismark, --parallel $bm_parallel) ..."
            bismark --parallel "$bm_parallel" --genome "$genome_dir" --temp_dir "$TMPDIR" \
                "${MERGE_DIR}/${SAMPLE_ID}_R0.fastq.gz" \
                --output_dir "${OUTDIR}"
            # Default SE output: {SAMPLE_ID}_R0_bismark_bt2.bam
            bm_bam="${OUTDIR}/${SAMPLE_ID}_R0_bismark_bt2.bam"
            samtools sort -@ "$THREADS" -m "$sort_mem" -o "$bam_full" "$bm_bam"
            rm -f "$bm_bam" "${OUTDIR}/${SAMPLE_ID}_R0_bismark_bt2"*report* "${OUTDIR}/${SAMPLE_ID}_R0_bismark_bt2"*nucleotide*
        fi
        ;;

    *)
        echo "ERROR: Unknown assay type '$ASSAY'"
        exit 1
        ;;
esac

# -------------------------------------------------------------------
# Subset BAM to target region → FASTQ
# -------------------------------------------------------------------
bam_subset="${OUTDIR}/tmp_${SAMPLE_ID}_subset.bam"

log "Indexing full BAM ..."
samtools index -@ "$THREADS" "$bam_full"

log "Subsetting to ${TARGET_REGION} ..."
samtools view -b -h -@ 4 "$bam_full" "$TARGET_REGION" -o "$bam_subset"

n_reads=$(samtools view -c "$bam_subset" 2>/dev/null || echo "?")
log "Subset reads: $n_reads"

if [[ "$LAYOUT" == "PE" ]]; then
    samtools sort -n -@ "$THREADS" -m "$sort_mem" -o "${bam_subset%.bam}.nsort.bam" "$bam_subset"
    samtools fastq -@ 4 \
        -1 "${OUTDIR}/${SAMPLE_ID}_R1.fastq.gz" \
        -2 "${OUTDIR}/${SAMPLE_ID}_R2.fastq.gz" \
        -0 /dev/null -s /dev/null \
        "${bam_subset%.bam}.nsort.bam"
    rm -f "${bam_subset%.bam}.nsort.bam"
else
    samtools sort -n -@ "$THREADS" -m "$sort_mem" -o "${bam_subset%.bam}.nsort.bam" "$bam_subset"
    samtools fastq -@ 4 \
        -0 "${OUTDIR}/${SAMPLE_ID}_R0.fastq.gz" \
        "${bam_subset%.bam}.nsort.bam"
    rm -f "${bam_subset%.bam}.nsort.bam"
fi

# Cleanup
rm -f "$bam_subset" "${bam_subset}.bai"
if [[ "$KEEP_ALIGNED" == "false" ]]; then
    rm -f "$bam_full" "${bam_full}.bai"
fi

log "Done: ${SAMPLE_ID}"
touch "${STATUS_DIR}/${SAMPLE_ID}.done"
BODYEOF

    chmod +x "$script"
}

# ===========================================================================
# PHASE 3: Gather — wait for all jobs and report
# ===========================================================================
gather() {
    log "=== Phase 3: Gather ==="

    local total=0 done=0 failed=0
    for i in $(seq 0 $((N_SAMPLES - 1))); do
        local sid="${SAMPLE_IDS[$i]}"
        local assay="${ASSAYS[$i]}"
        local layout="${LAYOUTS[$i]}"
        total=$((total + 1))

        if [[ -f "${STATUS_DIR}/${sid}.done" ]]; then
            if sample_complete "$sid" "$assay" "$layout"; then
                done=$((done + 1))
            else
                failed=$((failed + 1))
                log "  FAILED: $sid (status marker exists but output missing)"
            fi
        else
            failed=$((failed + 1))
            log "  FAILED: $sid (no completion marker)"
        fi
    done

    log ""
    log "========================================="
    log "  Total:     $total"
    log "  Completed: $done"
    log "  Failed:    $failed"
    log "========================================="

    if [[ $failed -gt 0 ]]; then
        log "Check logs in: ${LOGS_DIR}/"
        log "Re-run with the same command to retry failed samples (idempotent)."
        exit 1
    fi

    # Summary of output sizes
    log ""
    log "Output files:"
    du -sh "${OUTDIR}"/*.fastq.gz "${OUTDIR}"/*.bam "${OUTDIR}"/*.bedmethyl.gz 2>/dev/null | while read -r sz fn; do
        log "  $sz  $(basename "$fn")"
    done

    log ""
    log "All $total samples complete. Outputs in: $OUTDIR"
}

# ===========================================================================
# Clean intermediate files
# ===========================================================================
clean_intermediates() {
    log "=== Cleaning intermediate files ==="

    # Verify all samples completed before deleting anything
    local total=0 failed=0
    for i in $(seq 0 $((N_SAMPLES - 1))); do
        local sid="${SAMPLE_IDS[$i]}"
        local assay="${ASSAYS[$i]}"
        local layout="${LAYOUTS[$i]}"
        total=$((total + 1))
        if ! sample_complete "$sid" "$assay" "$layout"; then
            failed=$((failed + 1))
            log "  INCOMPLETE: $sid"
        fi
    done

    if [[ $failed -gt 0 ]]; then
        log "ERROR: $failed/$total samples incomplete — refusing to clean."
        log "Re-run without --clean to finish, then retry."
        exit 1
    fi

    local freed=0
    for dir in "$DL_DIR" "$MERGE_DIR" "$JOBS_DIR" "$STATUS_DIR"; do
        if [[ -d "$dir" ]]; then
            local sz
            sz=$(du -sh "$dir" 2>/dev/null | cut -f1)
            log "  Removing $dir ($sz)"
            rm -rf "$dir"
        fi
    done

    log "Done. Kept: output files + logs ($LOGS_DIR)"
}

# ===========================================================================
# Main controller logic
# ===========================================================================
main() {
    # Standalone modes
    if $CLEAN_MODE; then
        clean_intermediates
        exit $?
    fi

    if [[ "$_PHASE" == "gather" ]]; then
        gather
        exit $?
    fi

    # Phase 1: indexes
    log "Phase 1/3: Building genome indexes"
    local index_dep
    index_dep=$(build_indexes)
    log "  Index jobs: ${index_dep:-none needed}"

    # Phase 2: per-sample jobs
    log ""
    log "Phase 2/3: Spawning per-sample jobs"
    local sample_dep
    sample_dep=$(spawn_sample_jobs "$index_dep")
    log "  Sample jobs: ${sample_dep:-none needed}"

    if $LOCAL_MODE; then
        # In local mode everything ran inline — go straight to gather
        log ""
        gather
        exit $?
    fi

    if $DRY_RUN; then
        log ""
        log "[dry-run] Would resubmit controller for Phase 3 (gather)"
        exit 0
    fi

    # Phase 3: resubmit self as gather job
    if [[ -n "$sample_dep" ]]; then
        log ""
        log "Phase 3/3: Resubmitting controller for gather"

        local keep_flag=""
        $KEEP_ALIGNED && keep_flag="--keep-aligned"

        local gather_script="${JOBS_DIR}/gather.sh"
        cat > "$gather_script" << GATHEREOF
#!/usr/bin/env bash
set -euo pipefail
export _SUBSET_PHASE=gather
exec bash "${SCRIPT_PATH}" "${MANIFEST}" "${TARGET_REGION}" --outdir "${OUTDIR}" ${keep_flag}
GATHEREOF
        chmod +x "$gather_script"

        local gather_jid
        gather_jid=$(
            sbatch \
                ${SITE_FLAGS[@]+"${SITE_FLAGS[@]}"} \
                --job-name="subset-gather" \
                --cpus-per-task=1 \
                --mem="$CTRL_MEM" \
                --time="$CTRL_TIME" \
                --output="${LOGS_DIR}/gather_%j.log" \
                --error="${LOGS_DIR}/gather_%j.log" \
                --dependency="afterany:${sample_dep}" \
                "$gather_script" | grep -oP '\d+$'
        )
        log "  Gather job: $gather_jid (runs after all sample jobs)"
    else
        log "  No sample jobs submitted — running gather inline"
        gather
    fi

    log ""
    log "All jobs submitted. Monitor with: squeue -u \$USER -n 'subset-*'"
    log "Logs: ${LOGS_DIR}/"
}

main
