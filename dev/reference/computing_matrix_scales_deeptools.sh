# The shell block of combined_analysis.smk's computing_matrix_scales as it was
# before the scales moved to R (R_matrix_scales.R), kept verbatim -- Snakemake
# placeholders and doubled braces included -- as the reference the R port is
# compared against (dev/scales_rule.py, tests/unit/test_plot_scales.py,
# dev/plot_compare.py --family scales). Needs deeptools (dataRange,
# plotProfile) on PATH. Do not edit.
        {{
        # The track (mark) of a column is the last '__' field of its key
        # (define_key_for_plots), so rows of one mark are matched exactly
        # rather than by substring: H3K9me must not pick up H3K9me2.
        rows_of_mark() {{ awk -F"\t" -v m="$1" '{{n=split($1, f, "__"); if (n>1 && f[n]==m) print}}' "$2"; }}
        # One label per region group, counts taken from the matrix itself --
        # a stranded target can carry a separate 'nostrand' group, and
        # --skipZeros drops regions, so the target file cannot tell us either.
        python3 "{params.labeller}" {input.matrix} > {output.params_regions}

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
            
            # One range per mark, then one value per matrix column in column
            # order. Columns are not grouped by mark: an unstranded matrix of a
            # stranded track interleaves each sample's columns (s1__X_plus,
            # s1__X_minus, s2__X_plus, ...), so emitting a block per mark would
            # hand plus columns the minus scale and vice versa.
            declare -A zlo zhi ylo yhi
            while read mark
            do
                zmini=$(rows_of_mark "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=999999}} {{a=$5; if (a != "nan" && a<m) m=a;}} END {{if (m != 999999) print m; else print 0}}')
                zmaxi=$(rows_of_mark "${{mark}}" {output.temp_values} | awk 'BEGIN {{m=-999999}} {{a=$6; if (a != "nan" && a>m) m=a;}} END {{if (m != -999999) print m; else print 0}}')
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
                # Pad the axis by a tenth of the data's RANGE, not of its value.
                # Multiplying the endpoints (min*0.8, max*1.2) scales the padding
                # with the baseline, so data sitting high above zero in a narrow
                # band gets a huge axis: ColCEN TE mCG spans 41.6-50.1 and was
                # drawn on 33.3-60.1, filling a third of the plot. Genes were
                # unaffected only because their mCG dips near zero.
                ybounds=$(rows_of_mark "${{mark}}" {output.temp_profile_values} | awk '{{ for (i=3;i<=NF;i++) if ($i+0==$i) {{ if (n==0 || $i+0<lo) lo=$i+0; if (n==0 || $i+0>hi) hi=$i+0; n++ }} }} END {{ if (n==0) {{ print 0, 0; exit }} pad=(hi-lo)*0.1; if (pad<=0) pad=(hi<0?-hi:hi)*0.1; if (pad<=0) pad=0.01; print lo-pad, hi+pad }}')
                ymini=$(echo "${{ybounds}}" | cut -d" " -f1)
                ymaxi=$(echo "${{ybounds}}" | cut -d" " -f2)
                test=$(awk -v a=${{ymini}} -v b=${{ymaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ ${{test}} == "yes" ]]; then
                    ymini=("0")
                    ymaxi=("0.01")
                fi
                zlo[${{mark}}]=${{zmini}}; zhi[${{mark}}]=${{zmaxi}}
                ylo[${{mark}}]=${{ymini}}; yhi[${{mark}}]=${{ymaxi}}
            done < <(awk -F"\t" 'NR>1 {{n=split($1, f, "__"); if (!(f[n] in seen)) {{seen[f[n]]=1; print f[n]}}}}' {output.temp_values})

            zmins=()
            zmaxs=()
            ymins=()
            ymaxs=()
            # dataRange lists the columns in matrix order; each one's mark is
            # its key's last '__' field.
            while read mark
            do
                zmins+=("${{zlo[${{mark}}]}}")
                zmaxs+=("${{zhi[${{mark}}]}}")
                ymins+=("${{ylo[${{mark}}]}}")
                ymaxs+=("${{yhi[${{mark}}]}}")
            done < <(awk -F"\t" 'NR>1 {{n=split($1, f, "__"); print f[n]}}' {output.temp_values})
            
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
                zmini=$(awk -F"\t" -v s="${{sample}}" '$1==s' {output.temp_values} | awk '{{if ($5 != "nan") print $5; else print 0}}')
                zmaxi=$(awk -F"\t" -v s="${{sample}}" '$1==s' {output.temp_values} | awk '{{if ($6 != "nan") print $6; else print 0}}')
                test=$(awk -v a=${{zmini}} -v b=${{zmaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" && "${{sample}}" == *__mCG ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.cg_scale}")
                elif [[ "${{test}}" == "yes" && "${{sample}}" == *__mCHG ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.chg_scale}")
                elif [[ "${{test}}" == "yes" && "${{sample}}" == *__mCHH ]]; then
                    zmins+=("0")
                    zmaxs+=("{params.chh_scale}")
                elif [[ "${{test}}" == "yes" ]]; then
                    zmins+=("0")
                    zmaxs+=("0.005")
                else
                    zmins+=("$zmini")
                    zmaxs+=("$zmaxi")
                fi
                
                # Pad the axis by a tenth of the data's RANGE, not of its value.
                # Multiplying the endpoints (min*0.8, max*1.2) scales the padding
                # with the baseline, so data sitting high above zero in a narrow
                # band gets a huge axis: ColCEN TE mCG spans 41.6-50.1 and was
                # drawn on 33.3-60.1, filling a third of the plot. Genes were
                # unaffected only because their mCG dips near zero.
                ybounds=$(awk -F"\t" -v s="${{sample}}" '$1==s' {output.temp_profile_values} | awk '{{ for (i=3;i<=NF;i++) if ($i+0==$i) {{ if (n==0 || $i+0<lo) lo=$i+0; if (n==0 || $i+0>hi) hi=$i+0; n++ }} }} END {{ if (n==0) {{ print 0, 0; exit }} pad=(hi-lo)*0.1; if (pad<=0) pad=(hi<0?-hi:hi)*0.1; if (pad<=0) pad=0.01; print lo-pad, hi+pad }}')
                ymini=$(echo "${{ybounds}}" | cut -d" " -f1)
                ymaxi=$(echo "${{ybounds}}" | cut -d" " -f2)
                test=$(awk -v a=${{ymini}} -v b=${{ymaxi}} 'BEGIN {{if (a==0 && b==0) c="yes"; else c="no"; print c}}')
                if [[ "${{test}}" == "yes" ]]; then
                    ymins+=("0")
                    ymaxs+=("0.01")
                else
                    ymins+=("$ymini")
                    ymaxs+=("$ymaxi")
                fi
            done < <(awk -F"\t" 'NR>1 {{print $1}}' {output.temp_values})
            
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
        