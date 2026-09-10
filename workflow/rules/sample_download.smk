def return_log_sample(data_type, sample_name, step, paired):
    return os.path.join(REPO_FOLDER, RESULTS_DIR,data_type,"logs",f"tmp__{sample_name}__{step}__{paired}.log")
    
rule get_fastq_pe:
    output:
        fastq1 = temp(f"{RESULTS_DIR}/{{data_type}}/fastq/raw__{{sample_name}}__R1.fastq.gz"),
        fastq2 = temp(f"{RESULTS_DIR}/{{data_type}}/fastq/raw__{{sample_name}}__R2.fastq.gz")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        trimmed_fastqs = config['trimmed_fastqs'],
        exist_fastq1 = lambda wildcards: f"{RESULTS_DIR}/{wildcards.data_type}/fastq/trim__{wildcards.sample_name}__R1.fastq.gz",
        exist_fastq2 = lambda wildcards: f"{RESULTS_DIR}/{wildcards.data_type}/fastq/trim__{wildcards.sample_name}__R2.fastq.gz",
        ena_script = os.path.join(REPO_FOLDER, "workflow", "scripts", "ena_download.sh")
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "downloading", "PE"))
    conda: CONDA_ENV
    retries: 3
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" && -e "{params.exist_fastq1}" && -e "{params.exist_fastq2}" ]]; then
            printf "Fastqs already exist for PE {params.sample_name}\n"
            cp {params.exist_fastq1} {output.fastq1}
            cp {params.exist_fastq2} {output.fastq2}
        elif [[ "{params.fastq_path}" == "SRA" ]]; then
            numbers=$(echo "{params.seq_id}" | sed 's/,/ /g')
            outdir="{config[output_dir]}/{params.data_type}/fastq"
            ena_ok=true

            # Try ENA first (pre-compressed .fastq.gz)
            printf "Attempting ENA download for PE {params.sample_name} ({params.seq_id})\n"
            for nb in ${{numbers}}; do
                if ! bash "{params.ena_script}" "${{nb}}" "${{outdir}}" "PE"; then
                    ena_ok=false
                    break
                fi
            done

            if [[ "${{ena_ok}}" == true ]]; then
                # ENA succeeded — concatenate per-accession files into final outputs
                ena_files_r1=()
                ena_files_r2=()
                for nb in ${{numbers}}; do
                    ena_files_r1+=("${{outdir}}/${{nb}}_1.fastq.gz")
                    ena_files_r2+=("${{outdir}}/${{nb}}_2.fastq.gz")
                done
                cat "${{ena_files_r1[@]}}" > {output.fastq1}
                cat "${{ena_files_r2[@]}}" > {output.fastq2}
                rm -f "${{ena_files_r1[@]}}" "${{ena_files_r2[@]}}"
                printf "ENA download complete for {params.sample_name}\n"
            else
                # ENA failed — clean up partial downloads and fall back to fasterq-dump
                printf "ENA download failed, falling back to fasterq-dump for PE {params.sample_name}\n"
                for nb in ${{numbers}}; do
                    rm -f "${{outdir}}/${{nb}}_1.fastq.gz" "${{outdir}}/${{nb}}_2.fastq.gz"
                done
                fastq_files_r1=()
                fastq_files_r2=()
                for nb in ${{numbers}}; do
                    fasterq-dump -e {threads} --temp "${{TMPDIR:-/tmp}}" --outdir "${{outdir}}" "${{nb}}"
                    fastq_files_r1+=("${{outdir}}/${{nb}}_1.fastq")
                    fastq_files_r2+=("${{outdir}}/${{nb}}_2.fastq")
                done
                printf "\n{params.sample_name} ({params.seq_id}) downloaded via fasterq-dump\nCompressing R1 and R2 in parallel\n"
                half_threads=$(( {threads} / 2 ))
                [[ "${{half_threads}}" -lt 1 ]] && half_threads=1
                cat "${{fastq_files_r1[@]}}" | pigz -p "${{half_threads}}" > {output.fastq1} &
                pid_r1=$!
                cat "${{fastq_files_r2[@]}}" | pigz -p "${{half_threads}}" > {output.fastq2} &
                pid_r2=$!
                wait "${{pid_r1}}" "${{pid_r2}}"
                rm -f "${{fastq_files_r1[@]}}" "${{fastq_files_r2[@]}}"
            fi
        elif [[ "{params.seq_id}" == "URL" ]]; then
            # One or more '+'-separated R1,R2 URL pairs; merge by concatenation.
            printf "Downloading PE fastqs from URL(s) for {params.sample_name}\n  {params.fastq_path}\n"
            IFS='+' read -ra pe_parts <<< "{params.fastq_path}"
            r1_files=()
            r2_files=()
            n=0
            for pair in "${{pe_parts[@]}}"; do
                r1="${{pair%%,*}}"
                r2="${{pair#*,}}"
                t1="{output.fastq1}.part${{n}}"
                t2="{output.fastq2}.part${{n}}"
                curl --fail --show-error --location --max-redirs 5 \
                     --retry 3 --connect-timeout 30 --max-time 7200 \
                     --proto '=https,http' -o "${{t1}}" "${{r1}}" &
                pid_r1=$!
                curl --fail --show-error --location --max-redirs 5 \
                     --retry 3 --connect-timeout 30 --max-time 7200 \
                     --proto '=https,http' -o "${{t2}}" "${{r2}}" &
                pid_r2=$!
                wait "${{pid_r1}}" "${{pid_r2}}"
                r1_files+=("${{t1}}")
                r2_files+=("${{t2}}")
                n=$((n+1))
            done
            cat "${{r1_files[@]}}" > {output.fastq1}
            cat "${{r2_files[@]}}" > {output.fastq2}
            rm -f "${{r1_files[@]}}" "${{r2_files[@]}}"
            printf "URL download complete for PE {params.sample_name}\n"
        elif [[ "{params.seq_id}" == "EXPLICIT" ]]; then
            # One or more '+'-separated local R1,R2 pairs; merge by concatenation.
            printf "Copying explicit PE fastqs for {params.sample_name}\n  {params.fastq_path}\n"
            IFS='+' read -ra pe_parts <<< "{params.fastq_path}"
            : > "{output.fastq1}"
            : > "{output.fastq2}"
            for pair in "${{pe_parts[@]}}"; do
                r1="${{pair%%,*}}"
                r2="${{pair#*,}}"
                if [[ "${{r1}}" == *.gz ]]; then cat "${{r1}}" >> "{output.fastq1}"; else pigz -p {threads} -c "${{r1}}" >> "{output.fastq1}"; fi
                if [[ "${{r2}}" == *.gz ]]; then cat "${{r2}}" >> "{output.fastq2}"; else pigz -p {threads} -c "${{r2}}" >> "{output.fastq2}"; fi
            done
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q.gz 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q.gz 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying PE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q.gz "{output.fastq1}"
            cp "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q.gz "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying and gzipping PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*R1*f*q -c > "{output.fastq1}"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*R2*f*q -c > "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q.gz 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q.gz 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying PE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q.gz "{output.fastq1}"
            cp "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q.gz "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q 2>/dev/null | wc -l) -eq 1 ]] && [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q 2>/dev/null | wc -l) -eq 1 ]]; then
            printf "Copying and gzipping PE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*_1.f*q -c > "{output.fastq1}"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*_2.f*q -c > "{output.fastq2}"
        elif [[ $(ls -1 "{params.fastq_path}"/*"{params.seq_id}"*1*f*q 2>/dev/null | wc -l) -gt 1 ]]; then
            printf "Error: Too many fastqs found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\nThe seq_id used {params.seq_id} is likely not unique.\n"
            exit 1
        else
            printf "Error: No PE fastqs found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule get_fastq_se:
    output:
        fastq0 = temp(f"{RESULTS_DIR}/{{data_type}}/fastq/raw__{{sample_name}}__R0.fastq.gz")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        fastq_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        trimmed_fastqs = config['trimmed_fastqs'],
        exist_fastq0 = lambda wildcards: f"{RESULTS_DIR}/{wildcards.data_type}/fastq/trim__{wildcards.sample_name}__R0.fastq.gz",
        ena_script = os.path.join(REPO_FOLDER, "workflow", "scripts", "ena_download.sh")
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "downloading", "SE"))
    conda: CONDA_ENV
    retries: 3
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" && -e "{params.exist_fastq0}" ]]; then
            printf "Fastq already existing for SE {params.sample_name}\n"
            cp {params.exist_fastq0} {output.fastq0}
        elif [[ "{params.fastq_path}" == "SRA" ]]; then
            numbers=$(echo "{params.seq_id}" | sed 's/,/ /g')
            outdir="{config[output_dir]}/{params.data_type}/fastq"
            ena_ok=true

            # Try ENA first (pre-compressed .fastq.gz)
            printf "Attempting ENA download for SE {params.sample_name} ({params.seq_id})\n"
            for nb in ${{numbers}}; do
                if ! bash "{params.ena_script}" "${{nb}}" "${{outdir}}" "SE"; then
                    ena_ok=false
                    break
                fi
            done

            if [[ "${{ena_ok}}" == true ]]; then
                # ENA succeeded — concatenate per-accession files into final output
                ena_files=()
                for nb in ${{numbers}}; do
                    ena_files+=("${{outdir}}/${{nb}}.fastq.gz")
                done
                cat "${{ena_files[@]}}" > {output.fastq0}
                rm -f "${{ena_files[@]}}"
                printf "ENA download complete for {params.sample_name}\n"
            else
                # ENA failed — clean up partial downloads and fall back to fasterq-dump
                printf "ENA download failed, falling back to fasterq-dump for SE {params.sample_name}\n"
                for nb in ${{numbers}}; do
                    rm -f "${{outdir}}/${{nb}}.fastq.gz"
                done
                fastq_files=()
                for nb in ${{numbers}}; do
                    fasterq-dump -e {threads} --temp "${{TMPDIR:-/tmp}}" --outdir "${{outdir}}" "${{nb}}"
                    fastq_files+=("${{outdir}}/${{nb}}.fastq")
                done
                printf "\n{params.sample_name} ({params.seq_id}) downloaded via fasterq-dump\nCompressing files\n"
                cat "${{fastq_files[@]}}" | pigz -p {threads} > {output.fastq0}
                rm -f "${{fastq_files[@]}}"
            fi
        elif [[ "{params.seq_id}" == "URL" ]]; then
            # One or more '+'-separated URLs; merge by concatenation.
            printf "Downloading SE fastq(s) from URL for {params.sample_name}\n  {params.fastq_path}\n"
            IFS='+' read -ra url_parts <<< "{params.fastq_path}"
            tmp_files=()
            n=0
            for url in "${{url_parts[@]}}"; do
                tmp="{output.fastq0}.part${{n}}"
                curl --fail --show-error --location --max-redirs 5 \
                     --retry 3 --connect-timeout 30 --max-time 7200 \
                     --proto '=https,http' -o "${{tmp}}" "${{url}}"
                tmp_files+=("${{tmp}}")
                n=$((n+1))
            done
            cat "${{tmp_files[@]}}" > {output.fastq0}
            rm -f "${{tmp_files[@]}}"
            printf "URL download complete for SE {params.sample_name}\n"
        elif [[ "{params.seq_id}" == "EXPLICIT" ]]; then
            # One or more '+'-separated local FASTQs; merge by concatenation
            # (gzip members concatenate cleanly; non-gz inputs are compressed).
            printf "Copying explicit SE fastq(s) for {params.sample_name}\n  {params.fastq_path}\n"
            IFS='+' read -ra fq_parts <<< "{params.fastq_path}"
            : > "{output.fastq0}"
            for fq in "${{fq_parts[@]}}"; do
                if [[ "${{fq}}" == *.gz ]]; then
                    cat "${{fq}}" >> "{output.fastq0}"
                else
                    pigz -p {threads} -c "${{fq}}" >> "{output.fastq0}"
                fi
            done
        elif ls "{params.fastq_path}"/*"{params.seq_id}"*q.gz 1> /dev/null 2>&1; then
            printf "\nCopying SE gzipped fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            cp "{params.fastq_path}"/*"{params.seq_id}"*q.gz "{output.fastq0}"
        elif ls "{params.fastq_path}"/*"{params.seq_id}"*q 1> /dev/null 2>&1; then
            printf "\nCopying and gzipping SE fastq for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            pigz -p {threads} "{params.fastq_path}"/*"{params.seq_id}"*q -c > "{output.fastq0}"
        else
            printf "Error: No SE fastq found for {params.sample_name} ({params.seq_id} in {params.fastq_path})\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule run_fastqc:
    input:
        fastq = f"{RESULTS_DIR}/{{data_type}}/fastq/{{step}}__{{sample_name}}__{{read}}.fastq.gz"
    output:
        fastqc = f"{RESULTS_DIR}/{{data_type}}/reports/{{step}}__{{sample_name}}__{{read}}_fastqc.html"
    params:
        data_type = lambda wildcards: wildcards.data_type,
        step = lambda wildcards: wildcards.step,
        sample_name = lambda wildcards: wildcards.sample_name,
        read = lambda wildcards: wildcards.read
    conda: CONDA_ENV
    threads: 1
    shell:
        """
        fastqc -o "{config[output_dir]}/{params.data_type}/reports/" "{input.fastq}"
        """

rule process_fastq_pe:
    input:
        raw_fastq1 = f"{RESULTS_DIR}/{{data_type}}/fastq/raw__{{sample_name}}__R1.fastq.gz",
        raw_fastq2 = f"{RESULTS_DIR}/{{data_type}}/fastq/raw__{{sample_name}}__R2.fastq.gz"
    output:
        fastq1 = maybe_temp(f"{RESULTS_DIR}/{{data_type}}/fastq/trim__{{sample_name}}__R1.fastq.gz", config.get('keep_trimmed_fastqs', False)),
        fastq2 = maybe_temp(f"{RESULTS_DIR}/{{data_type}}/fastq/trim__{{sample_name}}__R2.fastq.gz", config.get('keep_trimmed_fastqs', False)),
        metrics = f"{RESULTS_DIR}/{{data_type}}/reports/trim_pe__{{sample_name}}.json",
        html_report = f"{RESULTS_DIR}/{{data_type}}/reports/trim_pe__{{sample_name}}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        adapter1 = lambda wildcards: config['adapter1'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        adapter2 = lambda wildcards: config['adapter2'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        quality_threshold = lambda wildcards: config['quality_threshold'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        min_read_length = lambda wildcards: config['min_read_length'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trim_front = lambda wildcards: config['trim_front'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimmed_fastqs = config['trimmed_fastqs']
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "trimming", "PE"))
    conda: CONDA_ENV
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" ]]; then
            printf "\nFastq for {params.sample_name} is already trimmed\n"
            cp {input.raw_fastq1} {output.fastq1}
            cp {input.raw_fastq2} {output.fastq2}
            printf '{{}}' > {output.metrics}
            touch {output.html_report}
        else
            printf "\nTrimming adapters for {params.sample_name} with fastp version:\n"
            fastp --version 2>&1

            fastp_args=""
            [[ "{params.adapter1}" != "auto" ]] && fastp_args+=" --adapter_sequence {params.adapter1}"
            [[ "{params.adapter2}" != "auto" ]] && fastp_args+=" --adapter_sequence_r2 {params.adapter2}"
            [[ {params.trim_front} -gt 0 ]] && fastp_args+=" --trim_front1 {params.trim_front} --trim_front2 {params.trim_front}"

            fastp --thread {threads} \
                --cut_tail --cut_tail_mean_quality {params.quality_threshold} \
                --length_required {params.min_read_length} \
                --detect_adapter_for_pe \
                $fastp_args \
                --in1 "{input.raw_fastq1}" --in2 "{input.raw_fastq2}" \
                --out1 "{output.fastq1}" --out2 "{output.fastq2}" \
                --json "{output.metrics}" --html "{output.html_report}"
        fi
        }} 2>&1 | tee -a "{log}"
        """
        
rule process_fastq_se:
    input:
        raw_fastq = f"{RESULTS_DIR}/{{data_type}}/fastq/raw__{{sample_name}}__R0.fastq.gz"
    output:
        fastq = maybe_temp(f"{RESULTS_DIR}/{{data_type}}/fastq/trim__{{sample_name}}__R0.fastq.gz", config.get('keep_trimmed_fastqs', False)),
        metrics = f"{RESULTS_DIR}/{{data_type}}/reports/trim_se__{{sample_name}}.json",
        html_report = f"{RESULTS_DIR}/{{data_type}}/reports/trim_se__{{sample_name}}.html"
    params:
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        adapter1 = lambda wildcards: config['adapter1'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        quality_threshold = lambda wildcards: config['quality_threshold'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        min_read_length = lambda wildcards: config['min_read_length'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trim_front = lambda wildcards: config['trim_front'][get_sample_info_from_name(wildcards.sample_name, samples, 'env')],
        trimmed_fastqs = config['trimmed_fastqs']
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "trimming", "SE"))
    conda: CONDA_ENV
    shell:
        """
        {{
        if [[ "{params.trimmed_fastqs}" == "True" ]]; then
            printf "\nFastq for {params.sample_name} is already trimmed\n"
            cp {input.raw_fastq} {output.fastq}
            printf '{{}}' > {output.metrics}
            touch {output.html_report}
        else
            printf "\nTrimming adapters for {params.sample_name} with fastp version:\n"
            fastp --version 2>&1

            fastp_args=""
            [[ "{params.adapter1}" != "auto" ]] && fastp_args+=" --adapter_sequence {params.adapter1}"
            [[ {params.trim_front} -gt 0 ]] && fastp_args+=" --trim_front1 {params.trim_front}"

            fastp --thread {threads} \
                --cut_tail --cut_tail_mean_quality {params.quality_threshold} \
                --length_required {params.min_read_length} \
                $fastp_args \
                --in1 "{input.raw_fastq}" \
                --out1 "{output.fastq}" \
                --json "{output.metrics}" --html "{output.html_report}"
        fi
        }} 2>&1 | tee -a "{log}"
        """

rule get_available_bam:
    output: 
        bam = temp(f"{RESULTS_DIR}/{{data_type}}/mapped/copied__{{sample_name}}.bam")
    params:
        seq_id = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "seq_id"),
        bam_path = lambda wildcards: get_sample_info_from_name(wildcards.sample_name, samples, "fastq_path"),
        sample_name = lambda wildcards: wildcards.sample_name,
        data_type = lambda wildcards: wildcards.data_type,
        aligned_bams = config['aligned_bams'],
        exist_bam = lambda wildcards: f"{RESULTS_DIR}/{wildcards.data_type}/mapped/final__{wildcards.sample_name}.bam"
    log:
        temp(return_log_sample("{data_type}","{sample_name}", "copy_bam", "either"))
    conda: CONDA_ENV
    shell:
        """
        {{
        if [[ "{params.aligned_bams}" == "True" && -e "{params.exist_bam}" ]]; then
            printf "\nFinal bam file already exists for {params.sample_name}\n"
            cp {params.exist_bam} {output.bam}
        elif ls "{params.bam_path}"/*"{params.seq_id}"*.bam 1> /dev/null 2>&1; then
            printf "\nCopying bam file for {params.sample_name} ({params.seq_id} in {params.bam_path})\n"
            samtools sort -@ {threads} -T "$TMPDIR/download.sort" -o "{output.bam}" "{params.bam_path}"/*"{params.seq_id}"*.bam
        elif ls "{params.bam_path}"/*"{params.seq_id}"*.sam 1> /dev/null 2>&1; then
            printf "\nSorting and converting sam file for {params.sample_name} ({params.seq_id} in {params.bam_path})\n"
            samtools sort -@ {threads} -T "$TMPDIR/download.sort" -o "{output.bam}" "{params.bam_path}"/*"{params.seq_id}"*.sam
        else
            printf "Error: No bam or sam file found for {params.sample_name} ({params.seq_id} in {params.bam_path})\n"
            exit 1
        fi
        }} 2>&1 | tee -a "{log}"
        """