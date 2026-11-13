def return_log_env(ref_genome, step):
    return os.path.join(REPO_FOLDER, RESULTS_DIR,"logs",f"tmp_{step}_{ref_genome}.log")

rule prepare_reference:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa",
        gff = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gff",
        gtf = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gtf",
        chrom_sizes = f"{GENOMES_DIR}/{{ref_genome}}/chrom.sizes",
        genome_stats = f"{GENOMES_DIR}/{{ref_genome}}/genome_stats.json",
        taxid = f"{GENOMES_DIR}/{{ref_genome}}/taxid.json",
        region_files = [f"{RESULTS_DIR}/combined/bedfiles/{{ref_genome}}__protein_coding_genes.bed", f"{RESULTS_DIR}/combined/bedfiles/{{ref_genome}}__all_genes.bed"],
        logs = lambda wildcards: [ return_log_env(wildcards.ref_genome, step) for step in ["fasta", "gff", "gtf", "chrom_sizes", "genome_stats", "taxid", "region_file"] ]
    output:
        chkpt = f"{RESULTS_DIR}/combined/chkpts/ref__{{ref_genome}}.done",
        log = os.path.join(REPO_FOLDER, RESULTS_DIR,"logs","ref_prep__{ref_genome}.log")
    localrule: True
    shell:
        """
        cat {input.logs} > {output.log}
        rm {input.logs}
        touch {output.chkpt}
        """

rule check_fasta:
    output:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa"
    params:
        fasta = lambda wildcards: config["genomes"][wildcards.ref_genome]['fasta_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "fasta"))
    conda: CONDA_ENV
    shell:
        """
        {{
        fasta_src="{params.fasta}"
        if [[ "$fasta_src" == http://* || "$fasta_src" == https://* ]]; then
            printf "\nDownloading fasta from URL: $fasta_src\n"
            tmpfile=$(mktemp --suffix=.fasta.dl)
            trap 'rm -f "$tmpfile"' EXIT
            curl --fail --silent --show-error --location --max-redirs 5 \
                 --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 1800 \
                 --proto '=https,http' -o "$tmpfile" "$fasta_src"
            url_path="${{fasta_src%%\\?*}}"
            if [[ "$url_path" == *.fa.gz || "$url_path" == *.fna.gz || "$url_path" == *.fasta.gz ]]; then
                pigz -p {threads} -dc "$tmpfile" > {output.fasta}
            elif [[ "$url_path" == *.fa || "$url_path" == *.fna || "$url_path" == *.fasta ]]; then
                mv "$tmpfile" {output.fasta}
            else
                printf "\nExtension of downloaded fasta unknown: $fasta_src\n"
                exit 1
            fi
            rm -f "$tmpfile"
        elif [[ ! -s "$fasta_src" ]]; then
            printf "\nFasta file for {params.ref_genome} does not exist:\n$fasta_src\n"
            exit 1
        elif [[ "$fasta_src" == *.fa.gz || "$fasta_src" == *.fna.gz || "$fasta_src" == *.fasta.gz ]]; then
            printf "\nGzipped fasta file found: $fasta_src\n"
            pigz -p {threads} -dc "$fasta_src" > {output.fasta}
        elif [[ "$fasta_src" == *.fa || "$fasta_src" == *.fna || "$fasta_src" == *.fasta ]]; then
            printf "\nUnzipped fasta file found: $fasta_src\n"
            cp "$fasta_src" {output.fasta}
        else
            printf "\nExtension of fasta file unknown, should be .fasta(.gz) or .fa(.gz):\n $fasta_src\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """
        
rule check_gff:
    output:
        gff = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gff"
    params:
        gff = lambda wildcards: config["genomes"][wildcards.ref_genome]['gff_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "gff"))
    conda: CONDA_ENV
    shell:
        """
        {{
        gff_src="{params.gff}"
        if [[ "$gff_src" == http://* || "$gff_src" == https://* ]]; then
            printf "\nDownloading GFF from URL: $gff_src\n"
            tmpfile=$(mktemp --suffix=.gff.dl)
            trap 'rm -f "$tmpfile"' EXIT
            curl --fail --silent --show-error --location --max-redirs 5 \
                 --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 1800 \
                 --proto '=https,http' -o "$tmpfile" "$gff_src"
            url_path="${{gff_src%%\\?*}}"
            if [[ "$url_path" == *.gff*.gz ]]; then
                pigz -p {threads} -dc "$tmpfile" > {output.gff}
            elif [[ "$url_path" == *.gff* ]]; then
                mv "$tmpfile" {output.gff}
            else
                printf "\nExtension of downloaded GFF unknown: $gff_src\n"
                exit 1
            fi
            rm -f "$tmpfile"
        elif [[ ! -s "$gff_src" ]]; then
            printf "\nGFF file for {params.ref_genome} does not exist:\n$gff_src\n"
            exit 1
        elif [[ "$gff_src" == *.gff*.gz ]]; then
            printf "\nGzipped gff file found: $gff_src\n"
            pigz -p {threads} -dc "$gff_src" > {output.gff}
        elif [[ "$gff_src" == *.gff* ]]; then
            printf "\nUnzipped gff file found: $gff_src\n"
            cp "$gff_src" {output.gff}
        else
            printf "\nExtension of gff file unknown, should be .gff*(.gz):\n $gff_src\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule check_gtf:
    input:
        gff = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gff"
    output:
        gtf = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gtf"
    params:
        gtf = lambda wildcards: config["genomes"][wildcards.ref_genome].get('gtf_file', '<auto>'),
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "gtf"))
    conda: CONDA_ENV
    shell:
        """
        {{
        override="{params.gtf}"
        if [[ "$override" != "<auto>" && -n "$override" ]]; then
            if [[ ! -s "$override" ]]; then
                printf "\nGTF file for {params.ref_genome} does not exist:\n$override\n"
                exit 1
            elif [[ "$override" == *.gtf.gz ]]; then
                printf "\nGzipped gtf file found: $override\n"
                pigz -p {threads} -dc "$override" > {output.gtf}
            elif [[ "$override" == *.gtf ]]; then
                printf "\nUnzipped gtf file found: $override\n"
                cp "$override" {output.gtf}
            else
                printf "\nExtension of gtf file unknown, should be .gtf(.gz):\n $override\n"
                exit 1
            fi
        else
            printf "\nDeriving GTF from GFF for {params.ref_genome} using gffread\n"
            gffread -T {input.gff} -o {output.gtf}
            if [[ ! -s {output.gtf} ]]; then
                printf "\nERROR: gffread produced empty GTF from {input.gff}. Supply a GTF explicitly via gtf_file.\n"
                exit 1
            fi
        fi
        }} 2>&1 | tee -a "{log}"
        """
        
rule check_chrom_sizes:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa"
    output:
        fasta_index = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa.fai",
        chrom_sizes = f"{GENOMES_DIR}/{{ref_genome}}/chrom.sizes"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome
    log:
        temp(return_log_env("{ref_genome}", "chrom_sizes"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "\nMaking chrom.sizes file for {params.ref_genome}\n"
        samtools faidx {input.fasta}
        cut -f1,2 {output.fasta_index} > {output.chrom_sizes}
        }} 2>&1 | tee -a "{log}"
        """

rule compute_genome_stats:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa",
        fai = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa.fai"
    output:
        stats = f"{GENOMES_DIR}/{{ref_genome}}/genome_stats.json"
    log:
        temp(return_log_env("{ref_genome}", "genome_stats"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "\nComputing genome stats for {wildcards.ref_genome}\n"
        total=$(awk '{{s+=$2}} END {{print s}}' {input.fai})
        n_bases=$(awk '!/^>/ {{gsub(/[^Nn]/,""); n+=length}} END {{print n+0}}' {input.fasta})
        effective=$((total - n_bases))
        star_nbases=$(python3 -c "import math; print(min(14, int(math.log2($total)/2 - 1)))")
        printf '{{"total_bases":%d,"n_bases":%d,"effective_size":%d,"star_sa_index_nbases":%d}}\n' \
            "$total" "$n_bases" "$effective" "$star_nbases" > {output.stats}
        printf "Genome stats: total=%d, N=%d, effective=%d, STAR_SAindexNbases=%d\n" \
            "$total" "$n_bases" "$effective" "$star_nbases"
        }} 2>&1 | tee -a "{log}"
        """

rule resolve_taxid:
    output:
        taxid_file = f"{GENOMES_DIR}/{{ref_genome}}/taxid.json"
    params:
        genus = lambda wildcards: config["genomes"][wildcards.ref_genome].get("genus", ""),
        species = lambda wildcards: config["genomes"][wildcards.ref_genome].get("species", ""),
        override = lambda wildcards: str(config["genomes"][wildcards.ref_genome].get("ncbi_taxid", ""))
    log:
        temp(return_log_env("{ref_genome}", "taxid"))
    conda: CONDA_ENV
    shell:
        """
        {{
        override="{params.override}"
        if [[ -n "$override" && "$override" != "<auto>" ]]; then
            printf "Using user-provided NCBI TaxId: %s\n" "$override"
            printf '{{"ncbi_taxid":"%s"}}\n' "$override" > {output.taxid_file}
        else
            printf "Looking up NCBI TaxId for {params.genus} {params.species}\n"
            taxid=$(datasets summary taxonomy taxon "{params.genus} {params.species}" \
                | python3 -c "import sys,json; d=json.load(sys.stdin); print(d['reports'][0]['taxonomy']['tax_id'])" 2>/dev/null) || true
            if [[ -n "$taxid" ]]; then
                printf "Resolved NCBI TaxId: %s\n" "$taxid"
                printf '{{"ncbi_taxid":"%s"}}\n' "$taxid" > {output.taxid_file}
            else
                printf "WARNING: Could not resolve NCBI TaxId for {params.genus} {params.species}. GO analysis requiring TaxId may fail.\n" >&2
                printf '{{"ncbi_taxid":null}}\n' > {output.taxid_file}
            fi
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule prep_region_file:
    input:
        chrom_sizes = f"{GENOMES_DIR}/{{ref_genome}}/chrom.sizes",
        gff = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.gff"
    output:
        region_file1 = f"{RESULTS_DIR}/combined/bedfiles/{{ref_genome}}__protein_coding_genes.bed",
        region_file2 = f"{RESULTS_DIR}/combined/bedfiles/{{ref_genome}}__all_genes.bed"
    params:
        ref_genome = lambda wildcards: wildcards.ref_genome,
        to_ascii   = os.path.join(REPO_FOLDER, "workflow", "scripts", "to_ascii.py")
    log:
        temp(return_log_env("{ref_genome}", "region_file"))
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "\nMaking a bed file with gene coordinates from {params.ref_genome}\n" >> {log} 2>&1
        awk -v OFS="\t" '$3=="gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | python {params.to_ascii} | bedtools sort -g {input.chrom_sizes} > {output.region_file1}
        awk -v OFS="\t" '$3~"gene" {{print $1,$4-1,$5,$9,".",$7}}' {input.gff} | python {params.to_ascii} | bedtools sort -g {input.chrom_sizes} > {output.region_file2}
        }} 2>&1 | tee -a "{log}"
        """
        
rule download_rfam:
    output:
        cm = f"{GENOMES_DIR}/rfam/Rfam.cm",
        clanin = f"{GENOMES_DIR}/rfam/Rfam.clanin",
        pressed = touch(f"{GENOMES_DIR}/rfam/Rfam.cm.pressed")
    log:
        os.path.join(REPO_FOLDER, RESULTS_DIR,"logs","download_rfam.log")
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Downloading Rfam covariance models and clan file\\n"
        mkdir -p {config[genome_dir]}/rfam
        curl -fSL https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz \
            | pigz -dc > {output.cm}
        curl -fSL https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin \
            > {output.clanin}

        printf "Pressing covariance model database\\n"
        cmpress -F {output.cm}
        }} 2>&1 | tee -a "{log}"
        """

rule build_structural_rna_db:
    input:
        fasta = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa",
        fai = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}.fa.fai",
        rfam_cm = f"{GENOMES_DIR}/rfam/Rfam.cm",
        rfam_clanin = f"{GENOMES_DIR}/rfam/Rfam.clanin",
        rfam_pressed = f"{GENOMES_DIR}/rfam/Rfam.cm.pressed"
    output:
        structural_fa = f"{GENOMES_DIR}/structural_RNAs/{{ref_genome}}/structural_rnas.fa",
        tblout = f"{GENOMES_DIR}/structural_RNAs/{{ref_genome}}/infernal.tblout"
    params:
        threshold = config.get("infernal_threshold", "ga"),
        bin_size_bp = 1000000
    log:
        return_log_env("{ref_genome}", "structural_rna")
    conda: CONDA_ENV
    shell:
        """
        {{
        printf "Building structural RNA database for {wildcards.ref_genome} via Infernal\\n"

        workdir=$(mktemp -d)
        trap 'rm -rf "$workdir"' EXIT
        chroms_dir="$workdir/chroms"
        mkdir -p "$chroms_dir"

        # Compute search space: total bases * 2 (both strands), in Mb
        search_space=$(awk '{{s+=$2}} END {{printf "%d", s*2/1000000}}' {input.fai})
        printf "Search space (-Z): %s Mb\\n" "$search_space"

        # Threshold flag
        threshold="{params.threshold}"
        if [[ "$threshold" == "ga" ]]; then
            threshold_flag="--cut_ga"
        else
            threshold_flag="-E $threshold --incE $threshold"
        fi
        printf "Infernal threshold: %s (flag: %s)\\n" "$threshold" "$threshold_flag"

        # Split genome: large contigs (>= bin_size) get their own file,
        # small contigs (< bin_size) are binned together
        bin_size={params.bin_size_bp}
        bin_idx=0
        bin_total=0
        bin_file="$chroms_dir/bin_${{bin_idx}}.fa"

        while IFS=$'\\t' read -r chrom size rest; do
            if [[ "$size" -ge "$bin_size" ]]; then
                samtools faidx {input.fasta} "$chrom" > "$chroms_dir/${{chrom}}.fa"
            else
                if [[ "$bin_total" -gt 0 && $(( bin_total + size )) -gt "$bin_size" ]]; then
                    bin_idx=$(( bin_idx + 1 ))
                    bin_total=0
                    bin_file="$chroms_dir/bin_${{bin_idx}}.fa"
                fi
                samtools faidx {input.fasta} "$chrom" >> "$bin_file"
                bin_total=$(( bin_total + size ))
            fi
        done < {input.fai}

        n_chunks=$(find "$chroms_dir" -name '*.fa' | wc -l)
        printf "Split genome into %d chunks (bin_size=%d bp)\\n" "$n_chunks" "$bin_size"

        # Determine parallelism: allocate threads across chunks
        if [[ "$n_chunks" -le 0 ]]; then
            printf "ERROR: No genome chunks produced\\n" >&2
            exit 1
        fi
        threads_per=$(( {threads} / n_chunks ))
        if [[ "$threads_per" -lt 1 ]]; then threads_per=1; fi
        n_parallel=$(( {threads} / threads_per ))
        if [[ "$n_parallel" -lt 1 ]]; then n_parallel=1; fi
        printf "Running %d parallel cmscan jobs with %d threads each\\n" "$n_parallel" "$threads_per"

        # Run cmscan per chunk
        find "$chroms_dir" -name '*.fa' | \
            xargs -P "$n_parallel" -I {{}} bash -c '
                cmscan --cpu '"$threads_per"' \
                    -Z '"$search_space"' \
                    '"$threshold_flag"' \
                    --rfam --nohmmonly --fmt 2 \
                    --tblout "{{}}.tblout" \
                    --clanin {input.rfam_clanin} \
                    {input.rfam_cm} \
                    {{}} > /dev/null
            '

        # Merge tblout results: keep header from first file, concatenate data lines,
        # remove lower-scoring clan overlaps (marked with " = " in olp column)
        printf "Merging and filtering Infernal results\\n"
        head -2 "$(find "$chroms_dir" -name '*.tblout' | head -1)" > {output.tblout}
        cat "$chroms_dir"/*.tblout | grep -v '^#' | grep -v ' = ' \
            | sort -k16,16g >> {output.tblout} || true

        # Count significant hits
        n_hits=$(grep -vc '^#' {output.tblout} || echo 0)
        printf "Found %d structural RNA loci\\n" "$n_hits"

        if [[ "$n_hits" -eq 0 ]]; then
            printf "WARNING: No structural RNAs found for {wildcards.ref_genome}. Creating empty FASTA.\\n" >&2
            touch {output.structural_fa}
        else
            # Convert tblout to BED (0-based coords)
            awk '$0 !~ /^#/ {{
                chrom = $4; start = $10; end = $11; name = $2; strand = $12
                if (strand == "+") {{ s = start - 1; e = end }}
                else {{ s = end - 1; e = start }}
                printf "%s\\t%d\\t%d\\t%s::%s:%d-%d(%s)\\t0\\t%s\\n",
                    chrom, s, e, name, chrom, s, e, strand, strand
            }}' {output.tblout} \
                | sort -k1,1 -k2,2n > "$workdir/structural.bed"

            # Extract sequences
            bedtools getfasta -fi {input.fasta} -bed "$workdir/structural.bed" \
                -s -name -fo {output.structural_fa}
        fi

        printf "Structural RNA database built: %s\\n" "{output.structural_fa}"
        }} 2>&1 | tee -a "{log}"
        """

def alt_te_source(wildcards):
    """Source of <G>__TE_file_<tag>.bed: the run's alt_te_file, if that's its tag."""
    g = wildcards.ref_genome
    if te_suffix(g) != "_" + wildcards.te_tag:
        raise ValueError(
            f"{g}__TE_file_{wildcards.te_tag}.bed was requested, but this run's "
            f"alternative TE file for {g} is "
            f"{config['genomes'].get(g, {}).get(ALT_TE_KEY) or 'not set'}"
            f"{' (tag ' + te_suffix(g)[1:] + ')' if te_suffix(g) else ''}."
        )
    return config["genomes"][g][ALT_TE_KEY]

# The TE-preparation shell is shared by the options-file annotation and an
# alternative one given for a run (alt_te_file, `epicc --te-file`). It is kept
# as a string, byte-identical to what check_te_file always ran, so moving it
# here does not trip Snakemake's code rerun-trigger on existing TE files.
TE_PREP_SHELL = """
        {{
        te_src="{params.te_file}"

        # If URL, download first to a temp file, then treat as local
        if [[ "$te_src" == http://* || "$te_src" == https://* ]]; then
            printf "\nDownloading TE file from URL: $te_src\n"
            url_path="${{te_src%%\\?*}}"
            dl_ext="${{url_path##*/}}"
            tmpfile=$(mktemp --suffix=".$dl_ext")
            trap 'rm -f "$tmpfile"' EXIT
            curl --fail --silent --show-error --location --max-redirs 5 \
                 --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 1800 \
                 --proto '=https,http' -o "$tmpfile" "$te_src"
            te_src="$tmpfile"
        fi

        # Determine format and convert to BED6
        if [[ "$te_src" == *.gff3.gz || "$te_src" == *.gff.gz ]]; then
            printf "\nGFF3 TE annotation found, converting to BED6: $te_src\n"
            pigz -p {threads} -dc "$te_src" | awk -F'\t' 'BEGIN{{OFS="\t"}} /^[^#]/ && NF>=9 {{
                id = ""
                n = split($9, attrs, ";")
                for (i=1; i<=n; i++) {{
                    if (attrs[i] ~ /^ID=/) {{ sub(/^ID=/, "", attrs[i]); id = attrs[i] }}
                }}
                if (id == "") id = "TE_" NR
                score = ($6 == "." ? "0" : $6)
                print $1, $4-1, $5, id, score, $7
            }}' > {output.te_file}
        elif [[ "$te_src" == *.gff3 || "$te_src" == *.gff ]]; then
            printf "\nUncompressed GFF3 TE annotation found, converting to BED6: $te_src\n"
            awk -F'\t' 'BEGIN{{OFS="\t"}} /^[^#]/ && NF>=9 {{
                id = ""
                n = split($9, attrs, ";")
                for (i=1; i<=n; i++) {{
                    if (attrs[i] ~ /^ID=/) {{ sub(/^ID=/, "", attrs[i]); id = attrs[i] }}
                }}
                if (id == "") id = "TE_" NR
                score = ($6 == "." ? "0" : $6)
                print $1, $4-1, $5, id, score, $7
            }}' "$te_src" > {output.te_file}
        elif [[ ! -s "$te_src" ]]; then
            printf "\nThe TE file for {wildcards.ref_genome} does not exist:\n $te_src\n"
            exit 1
        elif [[ "$te_src" == *.bed.gz ]]; then
            printf "\nGzipped TE file found: $te_src\n"
            pigz -p {threads} -dc "$te_src" > {output.te_file}
        elif [[ "$te_src" == *.bed ]]; then
            printf "\nUnzipped TE file found: $te_src\n"
            cp "$te_src" {output.te_file}
        else
            printf "\nExtension of TE file unknown, should be .bed(.gz) or .gff3(.gz):\n $te_src\n"
            exit 1
        fi

        # Transliterate non-ASCII characters (e.g. Greek letters in TE names)
        # to avoid deeptools computeMatrix assertion errors on region files.
        python {params.to_ascii} < {output.te_file} > {output.te_file}.tmp \
            && mv {output.te_file}.tmp {output.te_file}

        # Validate uniqueness of TE names (column 4)
        tot=$(cat {output.te_file} | wc -l)
        unique=$(cat {output.te_file} | cut -f4 | sort -u | wc -l)
        if [[ ${{unique}} -ne ${{tot}} ]]; then
            printf "\nNot all the names of TEs are unique. This is required for follow-up. Remove redundant rows or add unique identifiers.\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule check_te_file:
    output:
        te_file = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}__TE_file.bed"
    params:
        te_file    = lambda wildcards: config["genomes"][wildcards.ref_genome]['te_file'],
        ref_genome = lambda wildcards: wildcards.ref_genome,
        to_ascii   = os.path.join(REPO_FOLDER, "workflow", "scripts", "to_ascii.py")
    log:
        temp(return_log_env("{ref_genome}", "TEs"))
    conda: CONDA_ENV
    shell:
        TE_PREP_SHELL

rule check_alt_te_file:
    # <G>__TE_file_<tag>.bed for genomes.<G>.alt_te_file; the .source record
    # lets a later run with a different file of the same name pick a new tag
    # (te_files.resolve_te_tag) instead of overwriting this one.
    output:
        te_file = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}__TE_file_{{te_tag}}.bed",
        source  = f"{GENOMES_DIR}/{{ref_genome}}/{{ref_genome}}__TE_file_{{te_tag}}.source"
    wildcard_constraints:
        te_tag = r"[A-Za-z0-9_-]+"
    params:
        te_file    = lambda wildcards: alt_te_source(wildcards),
        ref_genome = lambda wildcards: wildcards.ref_genome,
        to_ascii   = os.path.join(REPO_FOLDER, "workflow", "scripts", "to_ascii.py")
    log:
        temp(return_log_env("{ref_genome}", "TEs_{te_tag}"))
    conda: CONDA_ENV
    shell:
        TE_PREP_SHELL + """
        printf '%s\\n' "{params.te_file}" > {output.source}
        """

