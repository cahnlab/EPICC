# Big Refactor TODO

## Pending

### Documentation

* [x] **high priority** Update README (and CLAUDE.md) with full functionality, especially about how to run epicc CLI, making sure to capture all changes on the big-refactor branch.

* [x] **high priority** Update Read the docs page to also reflect all changes on the big-refactor.

* [x] **high priority** Integrate README + Read the docs + epicc-builder (local HTML) for concerted changes

* [x] We should add a list of recommendations for which of ChIP_broad/narrow to use based on the histone mark, and/or refer users to a discussion of this distinction elsewhere (perhaps in the MACS documentation).

### Analysis

* [x] Expose DMR-call thresholds as `config/epicc-options.yaml` options. **Phase 1 done**: `dmr_thresholds:` block in options YAML exposes all per-context and per-caller parameters for both metilene and DMRcaller. Shared: `min_diff` (per-context), `min_cytosines`. DMRcaller-only: `bin_size`, `p_value`, `min_gap`, `min_size`, `min_reads`. Metilene-only: `max_cpgs`, `valley` (per-context), `maxseg` (per-context). Defaults preserve Arabidopsis-leaf-tuned behavior. Both R scripts read from CLI args; Snakemake rule passes all 10 args to both. Documented in README + CLAUDE.md.
   - **Phase 2 done**: `min_diff: "auto"` triggers a `calibrate_dmr_min_diff` rule that estimates `threshold = max(0.05, sigma_n × σ)` per-(pair, context), where σ = SD of per-position between-group methylation differences (noise-floor dominated; naturally scales CG > CHG > CHH with organism methylation levels). Writes a JSON consumed by `call_DMRs_for_pair_context`. Not the default — opt-in. `sigma_n` (default 3.0) is tunable.

### UI/UX

#### Configuration fileset

* [x] Parameters in the options file that are sub-settings in the yaml cannot be fed directly into snakemake command line. **Done**: `epicc run`/`validate` now accept named flags for the most commonly-tuned settings. Top-level keys go through `--config` directly; nested keys (`dmr_thresholds.*`, `chip_callpeaks.params`, `cut_callpeaks.*`) are applied via a temporary merged config (dot-path merge, sibling keys preserved). New flags: `--chip/atac-aligner`, `--chip/atac-mapping-strategy`, `--[no-]full-analysis`, `--[no-]te-analysis`, `--[no-]go`, `--methylation-contexts`, `--dmr-caller`, `--dmr-min-diff` (supports `auto`/single-float/`CG=f,CHG=f,CHH=f`), `--dmr-min-cytosines`, `--dmr-p-value`, `--dmr-sigma-n`, `--chip-callpeaks-params`, `--cut-broad-caller`, `--cut-narrow-caller`.

### Distribution

* [~] **High priority** Prepare for bioconda distribution - sticking point: snakemake execution profiles, investigate snakemake default search paths for these profiles, update README with information on how users can adapt to their cluster. **Done**: Added `epicc init-profile TYPE [NAME]` subcommand (copies bundled profile template to `~/.config/snakemake/<name>/config.yaml`, the XDG path Snakemake searches when `--profile <name>` is used). `find_repo_root()` now detects conda-installed layout (`$PREFIX/share/epicc/`); all command handlers skip `os.chdir(repo)` in installed mode and inject `--snakefile` so Snakemake finds the workflow without needing CWD = repo root. `check_profile_placeholders()` resolves XDG-named profiles. README cluster-setup section updated with both setup paths (init-profile and in-place edit). Packaging infrastructure: `pyproject.toml` (metadata + note on pip-install gap), `conda-recipe/meta.yaml` (Bioconda-ready recipe skeleton). Remaining gap: pip wheel install requires renaming `epicc` → `epicc.py` for `console_scripts`; deferred until Bioconda submission is confirmed working.

### Performance/Resource Usage

#### Full code review of snakemake rules

* [~] Do a full code review of the complete workflow with parallel subagents, and identify any potential inefficiencies for all sample type analysis rules and the combined analysis. **Bug-fix pass done** (2026-06-09): 6 parallel review agents over the rule modules + scripts; findings validated against the ColCEN run and the pombe dry-run. Fixed confirmed bugs: SAM-input `samtools sort -b` flag; SE `trimmed_fastqs` shortcut pointing at its own output; fail-loud on missing download inputs; ChIP peak-stats header/column misalignment (confirmed in ColCEN output); STAR cleanup glob deleting concurrent samples' files; bowtie2 aligner default fallback; `#`-truncation of sample-sheet fields; PE small-RNA guard; methylation calls on reference-absent chromosomes miscounted as CG; browser placeholder chrom substring match; empty-safe `seq_len()` loops; factor-consistency error row numbers. Design decisions implemented: N-factor RNA DEG grouping (was line+tissue only); custom-DMR sweep now honors `dmr_thresholds` + all contexts; `+`-merge now supported for FASTQ (SE+PE, local+URL) as well as SRA, with BAM/bedMethyl/mixed merges rejected (smoke-tested); prominent auto-unlock warning. Switched ColCEN test config to the bowtie2 default.
  - **Auto-unlock reverted** (2026-06-12): the "prominent auto-unlock warning" above was a stopgap on top of clearing locks by default. Reverted to snakemake's behavior — `epicc run` now halts on a stale lock and reports it; the `--no-auto-unlock` opt-out became an `--auto-unlock` opt-in. Also fixed a long-standing wart: `snakemake --unlock` / `epicc unlock` no longer require a sample sheet (unlock parses the whole Snakefile, which used to trip the no-`sample_file` guard). Unlocking now builds an empty `samples` frame and skips validation; `add_compat_columns`/`get_analysis_samples` made empty-DataFrame-safe via `_row_apply`. Verified: bare `snakemake --unlock` clears the lock with no sample sheet; 71 unit + 28 pombe-dryrun tests pass.
  - **Deferred follow-ups (need real-run verification / are perf not bugs):**
    - BAM input handling: `get_available_bam` globs the file path as a directory (`ls "$bam_path"/*"$seq_id"*.bam` with `bam_path` = the full file path) and has no URL-download branch — so local/URL BAM inputs through this rule look broken (untested; no test exercises it). Fix that before adding BAM `+`-merge (the reason BAM merge was left out of the FASTQ-merge feature).
    - Substring-matching cluster: UpSet `grep(type, colnames)` in `R_Upset_plot_{peaks,TSS,clusters}.R` and the `multiBigwigSummary` column-pick `$i~t` in `combined_analysis.smk` — column names embed the type as a substring, so the match is load-bearing; needs a real run to verify the exact (quoted) header format before anchoring. Bites only when one mark/label name is a substring of another.
    - Perf (out of scope for this bug pass): CX report re-read 3x per replicate (`R_cache_mc_replicate_for_context.R` — see the existing mC-bigwig TODO); bowtie2 sort hardcoded `-@ 2` while chromap uses `{threads}` (`ChIPseq.smk`); structural-RNA cmscan thread math underutilizes cores (`environment_setup.smk`); `dispatch_final_bam` cp vs hardlink.
    - RPKM: `R_gene_expression_rpkm.R` reports reads/kb without library-depth normalization while labeling the column RPKM — decide normalize vs rename. Decicion: Normalize.
    - Conversion-control chromosomes hardcoded to plastid names in `mC.smk` (no lambda/pUC spike-in support) — make a config key.
    - Snakemake `--dag`/`--rulegraph` (and thus `epicc validate --dag`) crash with this snakemake version on checkpoint DAGs (`AttributeError: 'NoneType'.edit_notebook` in `printdag`); pre-existing, unrelated to pipeline content. The 3 failing pombe-dryrun `TestDAGStructure` tests are this.

* [ ] Remove all functions/scripts used for backward compatibility with previous config. 

### Testing

#### ColCEN

* [x] colCEN integration test has become too large. Pare down, maybe remove some ChIP libraries that exercise the same code paths. **Done** (2026-05-29): removed rdr126ddm1hp5 CenH3/H3K9me2 (duplicate genotype), 3 extra dmC groups, rdr6 EMseq, tarp1tarp2 PBAT, and nrpd1 ATAC. Active sheet is 41 samples; removed rows commented out at bottom of TSV for reference.

#### H. sapiens test case

* [ ] Re-examine hg38 chr21 test data quality. The current dataset has limitations:
  * **sRNA**: EN-TEx sRNA libraries are rRNA-depleted total RNA, not true miRNA-seq. No miRNAs detected in the chr21 subset (zero reads in 21-24 nt range). Consider sourcing dedicated miRNA-seq libraries (e.g. TCGA miRNA-seq) or adding a second organism's sRNA test that exercises the miRNA discovery path.
  * **dmC**: Only a pre-computed bedMethyl sample (HG002). Adding a modBAM sample would exercise the `modkit pileup` path (alignment, pileup, `--combine-mods` handling) which is currently untested in integration.
  * **mC coverage**: ~~WGBS samples have very low chr21 mapping rates (~0.001%), resulting in sparse methylation data that causes PCA failures and uninformative DMR calls. May be inherent to the chr21 subset approach with whole-genome libraries.~~ **Fixed**: Root cause was Snakemake's `OMP_NUM_THREADS` export breaking bowtie2 inside bismark (see commit caa0809). WGBS samples now map at 99.8-99.9%. mCHG/mCHH PCA plots are still empty/corrupt because animals lack asymmetric methylation — those contexts contain only background false-positive calls.

### Known Issues/Bugs

* [x] Chromap PE rarely sets the `0x2` proper-pair flag on highly repetitive references (observed 0.02% on ColCEN CenH3, 33% on ColCEN ATAC, vs. typical ~80%+ for bowtie2). **Resolved**: A/B test on S. pombe PE ChIP confirmed root cause — chromap maps ends independently; when one mate multi-maps (MAPQ 0) and the other maps uniquely, `samtools view -q 10` keeps only the unique mate, producing orphan reads with no proper-pair flag. bowtie2 99.4–99.8% vs chromap 11–28% properly paired on non-repetitive pombe data. Fix: switched default aligner to bowtie2 (`chip_aligner: "bowtie2"`, `atac_aligner: "bowtie2"`). chromap remains available as a speed-optimized option for SE-only datasets.

## Epic-builder

* [x] Can deletion of rows require confirmation, or is there a way to have a "undo" arrow? (deleting a row or the whole table can be very quick - 1 misclick between duplicate and delete) **Done**: Added visible Undo/Redo toolbar buttons (wired to Tabulator's existing `history:true` — previously keyboard-only via Ctrl+Z/Y, undiscoverable). Added `confirm()` guards to both destructive paths: the context-menu "Remove row" (names the Sample_ID) and the bulk "Delete Selected" (shows the count); both note the action is undoable.

* [x] Remove tabs and white spaces from all columns before exporting sample sheet (can trigger column misalignement and errors in sample sheet loading if invisble tab present from copying SRR number for example) **Done**: New `sanitizeCell()` strips embedded tabs/CR/LF and trims each value; applied to every cell via the shared `buildExportMatrix()` used by both exporters. Verified end-to-end (a value pasted as `"  WT_rep1\t"` exports as `WT_rep1`).

* [x] When duplicating multiple rows, can the rows be as one block at the bottom of the table rather than each copied row under it's original? Easier to modifiy by sample type this way. **Done**: "Duplicate Selected" now orders copies by table position and appends them as one block at the bottom (was inserting each copy directly under its original). Single-row context-menu "Duplicate row" still inserts below the original.

* [x] Could an export option be a code block (i.e. printf "..." > sample_sheet.tsv) that can be pasted directly on cluster? with the sample file label chosen on top which automatically populates the epicc-options field? **Done**: New "Copy as command" button opens a modal with a paste-safe `printf '%s\t...' ... > <file>` command (real tabs come from the format string, so the pasted text carries literal `\t`; every cell is single-quoted with `'\''` escaping). A toolbar filename field names the target file, is reused by Export TSV, and populates the Options "Sample file path" (`sample_file`). Validated: the generated command reproduces the exact TSV in bash, incl. empty cells, single quotes, spaces, and `%`.

* [x] Motifs analysis can be merged with the other output options **Done**: The Motifs toggle, JASPAR database, and motifs-on-all-replicates controls moved into the "Output Options" section (after Gene Ontology); the standalone "Motifs Analysis" section was removed. Export/import keys unchanged.

* [x] Some details/explanations would be useful when selecting some parameters (WGBS directionality, RNAseq strandedness, ...) **Done**: Added `help` support to `buildSelectRow` and help text to RNAseq strandedness, mapped-reads/multimap, ChIP/ATAC aligner + mapping strategy, and CUT&RUN/CUT&Tag broad/narrow callers. For mC "directionality": discovered the builder's "Library prep method" (`mC_method`) select was **dead** — never exported, never read by the workflow (directionality is set per-sample by the Assay → `mC_mapping[Assay]`); replaced it with an accurate note. **Bonus bug fixed while here**: the config export wrote `rna_tracks` with only the `RNAseq` entry, omitting `RAMPAGE` — a builder-generated options file would `KeyError` at `config['rna_tracks']['RAMPAGE']` for any RAMPAGE dataset. Added round-trip-safe RAMPAGE strandedness/multimap controls and export/import.

## Deferred

### Documentation

* [ ] **Defer** Update README continuously after each additional modifications.

### Analysis

#### UMI handling (quantitative sequencing)

* [ ] **Defer** (important functionality gap; can close big-refactor without it) The pipeline does not extract or use UMIs (unique molecular identifiers). This matters for any quantitative protocol that carries them — UMI-based RNA-seq/RAMPAGE, and dedup-dependent assays — where PCR/optical duplicates should be collapsed *by UMI*, not by coordinate. Until this lands, mRNA-seq read deduplication is OFF by default (`rna_deduplicate: false`), because coordinate-based dedup without UMIs biases quantification (drops genuine high-expression reads). Scope: an optional sample-sheet UMI spec (location/length, or a separate index read), UMI extraction at the trimming step (e.g. `umi_tools extract` / fastp `--umi`), and UMI-aware dedup post-alignment (e.g. `umi_tools dedup` or `picard UmiAwareMarkDuplicatesWithMateCigar`). RAMPAGE currently follows ENCODE coordinate+Mate2basesN dedup (no true UMI), which would also be revisited.

#### (Differential) splicing analysis for RNAseq

* [ ] **deferred for future PR** heavy pipeline. Here is what was used in my MBD paper:
To look for novel splicing changes that occurred within the mC reader mutants, the reads were mapped with STAR (Dobin et al. 2013) were processed by StringTie and merged together (Pertea et al. 2015) into a master novel transcriptome comprising splicing events from TAIR10 and ones uniquely identified withinthis study. Then reads underwent lightweight alignment using Salmon version 1.4.0 (Patro et al. 2017) against the novel transcriptome from StringTie. Novel and known transcripts belonging to the same genewere analysed for splicing events by SUPPA2 (Trincado et al. 2018). Differential alternative splicing (DAS)was calculated for each event based on abundance of transcripts with and without inclusion of those eventsby SUPPA2 (Trincado et al. 2018).

#### Generic pre-computed bedMethyl input support

* **Defer** [ ] Support bedMethyl as a generic pre-computed methylation input format for any mC assay (WGBS, EMseq, dmC), not just dmC. Currently bedMethyl handling is exclusively gated behind dmC wildcard constraints (`_DMC_WC`). The underlying scripts (`validate_dmc_input.py`, `bedmethyl_to_cx_report.py`) are already assay-agnostic — the limitation is purely in rule routing.
  * **Current state**: A sample with `Assay: WGBS` and `Read_files: /path/to/precomputed.bed.gz` fails because (1) excluded from dmC rules by `_DMC_WC`, (2) `define_cx_report_input()` routes it to Bismark expecting FASTQs.
  * **Approach**: Add `bedMethyl` as a valid Assay in `VALID_ASSAYS` (maps to `mC` env). Modify `define_cx_report_input()` and wildcard constraints to route bedMethyl samples to the conversion pipeline. Reuse existing `get_dmc_input` → `copy_bedmethyl_input` → `convert_bedmethyl_to_cx_report` chain.
  * **Alternative**: Add an optional `Input_format` column (FASTQ/BAM/bedMethyl) to separate assay type from input format. More flexible but larger change.

#### ATAC-seq input sample support

* [ ] **Defer** Support calling ATAC peaks with Input

#### CUT&RUN/CUT&Tag follow-ups

* [ ] **Defer for later PR** Spike-in normalization for CUT&RUN (E. coli) and CUT&Tag (Drosophila). Requires a second alignment to the spike-in genome and a scale-factor injection into bigwig generation; useful for cross-sample normalization but not strictly necessary when using IgG-aware peak callers (SEACR with non-norm + stringent threshold, epic2 with control input).

* [ ] **Defer for later PR** CUT&RUN fragment-size splitting. The original Skene & Henikoff CUT&RUN protocol distinguishes TF/short-fragment peaks (<120bp) from histone/long-fragment peaks (>150bp) by filtering BEDPE intervals before bedgraph generation. Currently the entire fragment-size distribution is used. A flag like `cut_callpeaks.fragment_filter: short|long|all` could expose this for narrow-mark TF analyses.

* [ ] **Defer** CUT&Tag and CUT&RUN narrow (TF) ColCEN integration test data — currently only CUT_RUN_broad H3K27me3 from GSE123602 is in the test sheet. Identify a comparable public Arabidopsis CUT&Tag dataset (candidates: GSE201962 Ouyang 2022, PRJNA940156 Zhao 2023) and a CUT&RUN narrow/TF dataset.

#### mC contexts and analyses

* [ ] **defer** Add options to analyze subcontexts (e.g. CAG, CAA). This would require expanding the mC contexts options, sorting the cx-report based on the new patterns (potentially checking early that the pattern exists), allowing specific base nomenclature (W, N), checking whether DMRcaller works with custom contexts, and adapting the correct parameters for the corresponding type (e.g. CAA based on CHH parameters, CAG based on CHG).
 
* [ ] **defer** Rework the make_mc_bigwig_files to only have named output files matching the required contexts, without placeholders files for not necessary contexts (i.e. CHG and CHH if only ["CG"] present in the context list). This is better snakemake practice, and would allow sub-contexts bigwig files to be generated as well. This should still be done in a single rule however, to limit the number of time cx-report needs to be read, and filter each required contexts with the same awk line.

* [ ] **defer** Related to the configuration fileset above, the parameters for DMR analysis should be more easily customizable, to allow for different . While the custom script allows this, being able to pass the DMR parameters more explicitly could be useful, including as a list of different combinations of values to try (e.g. 50, 100, 150bp windows; min_c cutoff: 3,5,10; diff_cutoff: 0.1, 0.2). To do this, these parameters should either be top-levels in epicc-options (and allow for lists of parameters), or potentially managed through a specific CLI argument --DMR_parameters.


#### Differential peaks

* [ ] **defer** Summaries differential peaks from MA norm.

* [ ] **defer** Replace MAnorm with more recent approaches (DiffBind) which can incorporate replicates and adapt the code to use them appropriately. Provide a toggle to choose beteen using selected peaks, idr peaks or replicated peaks for differential analysis.

### UI/UX

#### custom adapter handling

* [ ] **defer for later PR** Sequencing adapters could vary on a per-sample basis. Maybe there should be an optional sample file column for custom adapters and we remove the global params from epicc-options.yaml.

#### Explicitly handle repeats vs coding gene annotations?

* [ ] **Defer for later consideration** Do we currently have a way to specify whether one or the other or both should be used in the analysis? Answer: Currently, analysis over gene is always on and there is a `te_analysis` toggle in the config to also do TE. Analysis over genes can be turned off with `full_analysis` toggle, but no TE analysis performed either in this case. Do we want to change this behavior into `full_analysis`=gene+TE (if TE file given), `gene_analysis_only`, `te_analysis_only` (if TE file given) and `no_analysis` options?

### Plotting

N.B. we'll work on plotting improvements in a separate branch after the Big Refactor is complete

* [ ] **Defer** See if we can improve browser plot sample label readability

* [ ] **Defer** In Plot Expression script, DOWN/UP color coding is based on "unique DEGS", i.e. genes differentially expressed in one sample vs all others, but does not show pairwise comparisons. Could be good to implement pairwise t-test with p-value/significance showing between all pairs of samples instead (or as a toggle option).

* [ ] **Defer** In plotting peak stats, for now only the first 2 reps are used (empty if not, and idr only between these 2). Would be best to allow for a flexible output where all reps are shown, and all pairwise idr too. Need refactoring the way stats are compiled.

* [ ] **Defer** Consider adding correlation matrix of coverage or pairwise plots between selected samples for more QC output

* [ ] **Defer** Add helper script, epicc command or a specific argument in epicc ouput to help generate input files for the different plot-type of epicc output from different output files of epicc run. For example, to make a files for --plot-type browser (bed-like with optional highlights) of the 24nt clusters DEG in one sample. Might be tricky to fully generalize, but should be relatively easy to do from main output (DMRs, peaks, RNAseq DEGs, sRNA clusters).

* [ ]  **Defer** Improve cutomization of plots (colors, number of metaplots per line, tick labels, etc..) through CLI parameters; maybe exclusively for epicc output; mimicking all/most/some deeptools options (colors may be done through exisitng define_key_for_plots), adding options to pass deeptools parameters directly to epicc command (e.g. epicc output --plot-type metaplot --deeptools startLabel=Start endLabel="")

### Codebase Hygiene

* [ ] **Defer** Improve logging system (naming, concatenating, and cleaning if chosen) - currently lots of redundancy.

* [ ] **Defer** Merging rules to call peaks for ChIP and TSS for RAMPAGE (both macs2), for merging regions (clusters, peaks, TSS)? Check existing conditions.

### Performance/Resource Usage

#### Default behavior adapted to number of samples

* [ ] **Defer** For all analysis or plotting which can include many samples, impose a threshold over which not to perform this analysis by default (which could still be forced with a specific flag). For example, do not perform heatmaps/metaplots on all genes if there are more than 20 samples in total. Suggest the list of analyses/plots to gate this way and a number of samples based on intensity of computation requirements for validation before installing these limitations.   

#### Speeding up bigwig conversion

* [ ] **Defer** Investigate faster libraries than UCSC - [bigtools](https://github.com/jackh726/bigtools), others? Should be actively maintained.

## Complete

### Documentation

* [x] Update README and CLAUDE.md to suggest creating a conda environment named epicc instead of smk9. Rename the config file epicc-env.txt. **Done**: Renamed `config/smk9.txt` → `config/epicc-env.txt`, updated all references in README.md, CLAUDE.md, validate_pombe.sh, tests/unit/README.md, test_rule_commands.py.

### Analysis

#### CUT&RUN / CUT&Tag support

* [x] Support CUT&RUN / CUT&Tag nomenclature + analysis. **Done**: Added four new assays (`CUT_RUN_broad`, `CUT_RUN_narrow`, `CUT_TAG_broad`, `CUT_TAG_narrow`) routed to the ChIP env via the new `IP_PEAK_ASSAYS` set (sample_sheet.py + samplefile_validation.py + epicc-builder). Peak calling defaults are peak-shape-aware: epic2 for `*_broad`, SEACR for `*_narrow`, MACS2 as explicit fallback (toggleable via `cut_callpeaks.{broad,narrow}_caller`). Both new tools' outputs are converted to UCSC broadPeak/narrowPeak (`workflow/scripts/convert_peaks.py`) so all downstream rules (idr, best_peaks_pseudoreps, peak_stats, motifs) work unchanged. ChIP rules `calling_peaks_macs2_pe/se` renamed to `calling_peaks_pe/se` and refactored as caller-aware dispatchers. Conda env gains `seacr` + `epic2`. Test coverage: 73 unit tests (converter + sample-sheet vocabulary), 17 dry-run integration tests covering all four assays + caller overrides. ColCEN test sheet adds three CUT&RUN samples from Zheng & Gehring 2019 (GSE123602) — H3K27me3 reps 1/2 sharing one IgG Control. Builder updated with all four assays, peaktype-aware known-marks recommendation, and per-family Control-source dropdowns. See sample-sheet-spec.md for the full spec.

### UI/UX

#### Handling IDR failures

* [x] How can we better handle IDR failures like in the case of the S. pombe test case? It shouldn't be a failure in the pipeline sense - it's an analytical outcome, and here reflects the biology of the S. pombe samples. If we can't work around the minimum requirement of 20 peaks by parameterizing IDR differently, we should also choose a different mark for the S. pombe test case that will be more likely to give a larger number of peaks. In any case, an IDR failure of this kind should not be a pipeline failure, and we should gracefully recover/continue and report. **Done**: `idr_analysis_replicates` now checks IDR exit code instead of blanket `|| true`. On failure: emits a clear warning banner (cause, peak counts, explanation that this is analytical not a bug), creates empty IDR output files, and continues. Downstream `best_peaks_pseudoreps` already handles IDR=0 correctly — peak selection uses merged ∩ pseudoreplicate consistency, IDR is informational only. `mv` of PNG plot is now conditional on IDR success. Postrun test updated to allow empty IDR peak files.

#### Configuration fileset

* [x] Change all references to the "config file" to refer to it as the "options file". This is currently config.yaml. It should become epicc-options.yaml. This will be an extensive rename - Make sure this change is uniformly applied throughout the workflow, documentation, test cases, and the builder. Keep the config/ directory name as-is. A full run configuration is actually the composite of a sample sheet and the information contained in the options file. **Done**: Renamed `config/config.yaml` → `config/epicc-options.yaml`, all test config files → `test_options_*.yaml`. Updated all references across workflow, rules, scripts, tests, builder, and documentation.

#### refactor sample sheet

##### New format

* [x] New sample sheet format streamlines and clarifies input specifications. Controls can now be arbitrary samples (WCE for yeast ChIP, RNA-seq for RAMPAGE, etc.) referenced by Sample_ID rather than a boolean flag. **Done**: New format with 9 columns (Sample_ID, Assay, Genome, Levels, Replicate_ID, Read_files, Read_layout, IP_target, Control). Central logic in `workflow/scripts/sample_sheet.py`. Migration script at `scripts/migrate_sample_sheet.py`.
  * [x] Update Snakefile sample-sheet parsing to read the new columns and build sample metadata accordingly. **Done**: `read_sample_sheet()` + `add_compat_columns()` in Snakefile.
  * [x] Update rule files: automatically generated filenames will change (e.g. "Input" -> the control's Sample_ID). **Done**: All 8 rule files migrated (ChIPseq, RNAseq, ATACseq, smallRNA, mC, combined_analysis, sample_download, environment_setup). TF env eliminated.
  * [x] Update documentation (README, Read the Docs, example sample sheets) to reflect the new format. **Done**: README.md and CLAUDE.md updated.
  * [x] Update test sample sheets and test code. **Done**: All test sample sheets converted (pombe, colcen, hg38_chr21, dmc). Unit tests for sample_sheet.py (57 tests).
  * [x] Validate with a full run of the S. pombe integration test. **Done**: 257 pipeline steps, all completed successfully.

  | Sample_ID | Assay | Genome | Levels | Replicate_ID | Read_files | Read_layout | IP_target | Control |
  |-----------|-------|--------|--------|--------------|------------|-------------|-----------|---------|
  | [freetext] | [ChIP_broad, ChIP_narrow, ATAC, RNAseq, RAMPAGE, sRNA, WGBS, WGBS_nd, PBAT, EMseq, dmC] | [freetext] | [freetext] | [freetext] | FASTQ SE:[/path/to/file/name.r1.fq], FASTQ PE:[/path/to/file/name.r1.fq,/path/to/file/name.r2.fq], BAM SE or PE: [/path/to/file/name.bam], SRA: [SRRxxxxx], SRA merge multiple:   [SRRxxxxx+SRRxxxxx+SRRxxxxx] | [SE or PE] | [freetext name of IP target or control, e.g. H3K9me2, or Input, WCE, etc.] required for ChIP_broad/ChIP_narrow | [valid sample ID] |

  Sample_ID: a name that uniquely identifies this sample. Will be used to track the sample internally, and can be used to assign controls to ChIP_broad, ChIP_narrow, and RAMPAGE Assays. We will not enforce any format, other than uniqueness, but the epicc-builder app should suggest a concise ID (see epicc-builder specification).
  Assay: controlled vocabulary, replaces data_type/sample_type and provides the menu of accepted assay   types for analysis.
  Genome: Reference genome name
  Levels: Comma-separated list of factor:level pairs describing the experimental conditions for this sample. Uses the syntax `factor_name:level_value` (e.g. `genotype:WT,tissue:root`). Factors can be any experimental variable — genotype (B73, Mo17), tissue (root, leaf), temperature (37deg, 24deg), time point (T0, T1), etc. If multiple factors are specified, multifactorial comparisons will be performed. Factor names are not currently used in pipeline logic (only levels are used for constructing comparisons and plot labels), but they improve readability and enable epicc-builder to parse and validate entries progressively. **All samples must have the same number of factor:level pairs, including controls.** Controls should still specify meaningful levels where possible (e.g. `genotype:WT,tissue:root`).
  Replicate_ID: identifies biological or technical replicates (e.g. rep1, rep2, repA, repB, 1, 2). Replicates are treated independently in analysis and merged only for specific downstream steps like peak calling. Note: `+` in Read_files is for concatenating multi-file inputs from the same library (e.g. multiple lanes or runs) before any processing — this is distinct from replicate handling.
  Read_files: Path to FASTA, FASTQ, BAM files, or bare SRA IDs. In this last case, read files will be downloaded from SRA. For paired-end FASTQs with separate mate files, use a comma to separate the R1 and R2 (/path/to/file.R1.fastq.gz,/path/to/file.R2.fastq.gz). If multiple read files or SRA IDs are separated by a "+", they will be merged before any processing.
  Read_layout: controlled vocabulary, single-end or paired-end sequencing (SE or PE)
  IP_target: Required for all ChIP_broad and ChIP_narrow samples, including controls. Describes what was pulled down — e.g. `H3K9me2` for an IP sample, or `Input`, `WCE`, `IgG` for a control sample.
  Control: Must be a valid Sample_ID. Used only for normalizing ChIP_broad, ChIP_narrow, and RAMPAGE samples. Validation rules: (1) the referenced Sample_ID must exist in the sheet, (2) a control sample must not itself reference another control (no chaining), (3) it is an error to specify a Control for non-ChIP/RAMPAGE assays, (4) multiple IP samples may share the same control.

###### Input validation

* [x] Per-field validation rules are defined in [`dev/docs/sample-sheet-spec.md`](sample-sheet-spec.md) (the canonical specification) and implemented in `workflow/scripts/samplefile_validation.py`. The epicc-builder app should implement the same rules. See the spec file for the full list of per-field constraints, cross-field checks, and derived name definitions.

* [x] Add a check to ensure multiple samples don't share the same inputs files or SRA accessions. This is almost certainly a user data entry error, can't think of a reasonable use case. Should be validated both in the builder and in the workflow. **Done**: Cross-row duplicate check in `samplefile_validation.py` and `epicc-builder.html`. Documented in `sample-sheet-spec.md`.

* [x] Make fatal the non-fatal warnings for PE with a single (non-BAM, non-SRA) path (i.e. just one FASTQ) and SE with more than one path passed with comma separation. **Done**: Changed from warnings to errors in both `samplefile_validation.py` and `epicc-builder.html`. Updated `sample-sheet-spec.md` severity.

* [x] epicc-builder: don't allow the user to generate an invalid sample sheet. The export button should be inactive if the sheet is in a validation failed state, and we should instruct the user to fix errors in the sheet to enable export (could be a hover tool tip, let's try to be screenreader friendly). **Done**: Export button disabled when errors exist, with tooltip and `aria-disabled` attribute.

* [x] epicc-builder: IP_target should not be editable if Assay is not chip_broad or chip_narrow. **Done**: Conditional `editable` function on IP_target column checks `CHIP_ASSAYS.has(assay)`.

* [x] epicc-builder: Unless it has been edited by the user, the Sample_ID suggestion should be updated automatically and continuously as the user changes other fields for the sample. **Done**: `_sid_user_edited` flag tracks manual edits; `applySuggestions()` auto-updates cell value when not user-edited; clearing Sample_ID resumes auto-suggestion.

* [x] If full path(s) are given for samples, verify that the files exist before starting the run and return detail error. **Done**: `check_table()` and `check_genome_config()` in `samplefile_validation.py` now probe local Read_files paths (PE mates, `+`-merged components) and per-genome reference paths (`fasta_file`, `gff_file`, `gtf_file`, `te_file`, `structural_rna_fafile`, `gaf_file`, `gene_info_file`). SRA accessions, HTTP(S) URLs, and the `<auto>` sentinel are skipped. Both validators gain a `check_paths` flag for callers that need a bypass. Integration tests get a session-autouse conftest fixture that stages the missing mock placeholders so the dry-run tests stay self-contained without committing binary stubs. Documented in `dev/docs/sample-sheet-spec.md`.

* [x] Preemptive input checks for extra-output target files. **Done**: `check_extra_output_files()` in `samplefile_validation.py` validates the optional target files referenced in the options file (`browser_target_file`, `heatmap_target_file`, `motif_target_file`, `rnaseq_target_file`, `srna_target_file`, `rnaseq_background_file`). Each check is gated on the corresponding analysis flag (`full_analysis`, `motifs`, `GO`) so users who haven't enabled an analysis aren't required to supply a real file. Default placeholder paths (e.g. `data/target_loci.bed`) are recognized and skipped to avoid false positives on out-of-the-box configs. The `browser_target_file` validator enforces the extended-BED schema described in this TODO entry: chrom + start + end + label (non-empty, must not start with `-`) + binsize (integer ≥ 1) + optional comma-separated `htstart`/`htwidth` (paired or absent). Wired into the Snakefile alongside `check_table` / `check_genome_config`. Tests in `tests/unit/test_samplefile_validation.py`.

##### Harmonize analysis_samplefiles

* [x] Analysis samplefiles need to be updated to follow the conventions of the new input samplefile format. **Done**: `analysis_samplefile_*.tsv` now uses new-format columns (Assay, Genome, Levels, IP_target, env, sample_name).
  - **Obsolete** (2026-09-09): the file is no longer written. Nothing ever read it back — it was a MaizeCode remnant from when an analysis could be started *from* an analysis_samplefile, and re-adding that entry point is not planned. The write also targeted `REPO_FOLDER/config/`, which is wrong for conda installs where the repo sits under `$PREFIX/share/epicc/`.

##### Remove old format detection in samplefile_validation.py

* [x] Old format TSVs didn't always have headers, so we can't rely on column names for detection. Let's just get rid of the detection routine and make sure the README.md informs users that they can run migrate_sample_sheet.py if they're using the old format. **Done**: Removed old-format detection from `samplefile_validation.py` and `sample_sheet.py`. Removed `OLD_COLNAMES` constant (migration script defines its own). Added usage example to README.md and sample-sheet-spec.md.

* [x] Ensure migrate_sample_sheet.py is up to date and remember to update it if any changes are made to the new sample sheet format. **Done**: Verified — handles all 9 new columns, all assay types, Control linking, and Sample_ID generation. **Superseded**: the script was removed ahead of 1.0.0 along with legacy assay-token support; there is no longer a migration path to maintain.

##### Derived Names

* [x] For samples without {IP_target}, we should avoid "____" in the filenames. Since many samples could lack an IP_target, best if we can do this by just eliminating the term and spacer, avoiding a placeholder like __NA__. This may already be the case, in which case we should just document the fact in sample-sheet-spec.md. **Done**: `build_analysis_name()` now filters empty parts before joining. E.g. `RNAseq__WT_cell__Spombe` (not `RNAseq__WT_cell____Spombe`). Documented in sample-sheet-spec.md.

##### epicc-builder

* [x] Currently, epicc-builder (as referenced in the README) is a standalone web app hosted remotely, and currently broken. We should develop a new version of epicc-builder as a self-contained HTML5/javascript app that helps users create a valid sample sheet with a tabular GUI. It will be deployed as a single HTML file that can be opened offline in any modern browser. **Done**: `tools/epicc-builder.html` (Tabulator 6.x, ~1165 lines), symlinked at repo root.

###### [x] Implementation Stage 1: Sample sheet preparation helper

**Prerequisite research**: Identify a suitable JS library for table drawing andwidgets (e.g. Handsontable, AG Grid, or a lightweight alternative that can bebundled into a single HTML file).  **Decision reached: use [tabulator] <https:/tabulator.info/>**.

**I/O**:

* [x] Output format: TSV with the column header as the first row, matching thepipeline's expected input format. **Done**.
* [x] Import an existing sample sheet (TSV) for editing. **Done**.

**UX**:

* [x] Accessibility (to e.g. screen readers) is an important concern. **Done**: aria labels, roles, keyboard-contained Tab navigation.
* [x] When focus (keyboard or mouse) is on a column, the description of that columnshould be shown to the user via a dynamic text display. **Done**: description panel updates on header/cell click.
* [x] Hamburger menu on each sample row allows for common actions like "Add areplicate", "Insert duplicate below", "Remove sample", etc. **Done**: context menu with Add replicate, Duplicate, Move up/down, Remove.
* [x] Sample rows are reorderable. **Done**: Tabulator movableRows + Move up/down in context menu.
* [x] Example text is shown for the freetext fields until a user enters information.For example, epicc-builder has the opportunity to progressively suggest uniquesample IDs as a user fills out the other fields for that sample. The user can editthe suggestion. **Done**: per-row example data (dimmed, individually cleared on edit), progressive Sample_ID suggestion with auto-accept.
* [x] Controlled-vocabulary fields are represented as a drop-down menu. **Done**: Assay, Read_layout, Control as list editors.
* [x] The drop-down for Control should be filtered based on the Assay of the currentrow to show only compatible types. For ChIP_broad and ChIP_narrow, both are validAssay sources for Control. For RAMPAGE, only RNAseq Assays are valid Controls. Ifcontrols are available, instead of an empty dropdown list, show an unselectablemessage item "None available" or similar. **Done**: CONTROL_SOURCE_ASSAYS constant drives filtering.
* [x] Replace underscores in the header row with spaces for display purposes only. **Done**.
* [x] Add a selection widget to each row to support multi-row selection and actions like delete and duplicate. **Done**: Checkbox column with `rowSelection` formatter, `selectableRows: true`, "Delete Selected" and "Duplicate Selected" toolbar buttons.
* [x] Change the EPICC branding home button in the upper left corner to EPICCbuilder (with builder italicized and in a complementary color). **Done**: Brand text changed to EPICC*builder* with "builder" in `--accent2` teal, italicized, lighter weight.

**Validation** (see [`dev/docs/sample-sheet-spec.md`](sample-sheet-spec.md) for thecanonical rules):

* [x] Perform the same input validation as the pipeline code, and give usersfeedback through diagnostic messages. **Done**: full validation engine ported from samplefile_validation.py, diagnostics panel.
* [x] Continuously evaluate user input as the table is filled out, opportunisticallyassigning defaults to sample column entries when there is sufficient input to do so. **Done**: editing-aware validation scheduling, auto IP_target clearing on Assay change.
* [x] Cells with unresolved warnings should show an indicator of this in the UI.Possibilities are yellow cell border or small yellow circle or exclamation mark inthe cell. Same for errors, but with red. **Done**: red left border for errors, yellow for warnings.
* [x] When a user enters an IP_target, we should compare with the list of knownmarks, and if the user has entered a ChIP_broad or ChIP_narrow value that conflictswith the recommendation in the known marks, we should warn. **Done**: `KNOWN_MARKS` constant with regex patterns from `config/epicc-options.yaml` peaktype mappings. Warning emitted during validation when assay conflicts with recommendation.

###### [x] Implementation Stage 2: Config file helper and concordance with original hosted app

* [x] Augment the sample sheet builder with another interactive section just above allowing for the creation, editing, removal of the factors that will be used in Levels in the sample file. This can be something like an New Factor text field and button. Committed factors will then appear as tiles showing their name with an "x" inside to enable removal. Then, in the sample sheet builder, instead of showing the Levels as a table column, we should automatically populate/update the table GUI with columns having the names of the factors, so the users will just enter the level value there. The TSV will keep the same factor:level format for the Levels column, so epicc-builder should handle this transformation on import/export. **Done**: Factor panel with add/remove tiles, dynamic `_factor_*` columns replacing Levels, import/export round-trips correctly.

* [x] Implement all currently unimplemented functionality from [the original epicc-builder](https://epicc-builder.streamlit.app/), but keep our updated sample sheet builder. **Done**: Config file form with all sections from original Streamlit app (required params, per-genome/species, output/input options, motifs, advanced ChIP/ATAC/RNA/mC/sRNA/plotting), YAML import/export via js-yaml@4.
  * [x] We will need to alter the app presentation so we can show different sections (sample sheet and epicc configuration file form settings). This could be implemented with upper nav bar with "EPICC" as the home item in the menu, followed by Sample Sheet and Config File. Clicking EPICC resolves to the first menu item (Sample Sheet). Links to documentation and the EPICC/epigeneticbutton [Github repo](https://github.com/cahnlab/EPICC) should also appear in this main nav menu. **Done**: Top nav bar with EPICC brand, Sample Sheet / Config File tabs, Docs / GitHub links.
  * [x] The dedicated sample file check button probably isn't necessary since we're doing real time validation. **Done**: Real-time validation replaces the button.
  * [x] Keep all of the jokey prompts in the options file form. **Done**: All jokey prompts from the original Streamlit app preserved (full_analysis, te_analysis, GO, trimmed_fastqs, aligned_bams, motifs, nextflex, structural RNA, stranded heatmaps, browser TE, etc.).

###### [x] Deployment

* [x] Host the new epicc-builder on Github Pages so users who are primarily working with a remote installation of epicc don't necessarily have to download the HTML file and run it locally. **Done**: GitHub Actions workflow `.github/workflows/deploy-builder.yml` deploys `tools/epicc-builder.html` as `index.html` on push to main. Requires enabling GitHub Pages (Source: GitHub Actions) in repo Settings > Pages.

###### [x] Tweaks

* [x] Show help text for column headers on mouse over of headers specifically, instead of just on click. Retain current show-on-click behavior for header row and value cells. **Done**: Added `headerMouseEnter` callbacks alongside existing `headerClick` and `cellClick` in `descHooks()` and factor column inline hooks.

* [x] Add an explanation of what Factors are to the right of the +Add button. Something like: "Experimental variables used for grouping samples for comparative analysis. All rows must have the same factors, and may have different levels (values)." **Done**: Added `.factor-help` span with muted styling after the +Add button.

* [x] Promote the Reference Genomes section out of the current Config file form to a separate top-level menu tab appearing just to the right of the EPICC brand (and becoming the new link destination for EPICC). The new layout becomes Reference Genomes | Samples | Options. The data from the Reference Genomes page will constrain which Reference Genomes are available to choose from in the sample sheet (Genomes field becomes a dropdown list), and will still be used as it is now to populate the relevant section of the options file. **Done**: New Reference Genomes tab with genome name management UI (add/remove tiles), per-genome config forms, Genome column in Samples tab constrained to dropdown (freetext still allowed). Genome data flows into YAML export via `configState.genomes`. EPICC brand links to genomes tab.

###### [x] Bugs

* [x] Sample sheet example should validate, be exportable. **Done**: Fixed PE example rows to include R1,R2 comma-separated paths. Removed `_example` skip guards from validation and export. Example rows now validate normally and are included in TSV export.

* [x] Help text for the factor levels columns should be better. Maybe should say: Levels   of factor {factor}. Stored as comma-separated factor:level pairs (e.g. genotype:WT,tissue:root) in the sample sheet TSV 'Levels' column. All rows must have the same factors. Level values form the 'levels_label' in analysis names. **Done**: Factor columns now show factor-specific description via `showFactorDescription()` instead of generic "Levels" text.

* [x] Fix rows selected on click in any column - should only be selected when clicking on the leftmost column selection box. **Done**: Added `selectableRowsCheck` returning false to block automatic row-click selection; checkbox column's explicit `toggleSelect()` call bypasses the check.

#### species-specific parameters

* [x] Reference Genome Species (as in Species-dependent parameters in epicc-options.yaml) should probably be defined as their binomial like Zea_mays to avoid collisions. Could we just get rid of that section altogether and stick ncbiID and go_database along with the params for each reference genome? We can just compute the genome size of the reference and not bother with it in the config, same with —genomeSAindexNbases, no? **Done**: Species params inlined into genome entries (earlier refactor). `genomesize` and `star_index` are now auto-computed from the reference FASTA by `compute_genome_stats` rule; user values in options file are optional overrides. Removed from test configs; kept in main config as override examples.

* [x] Let's move the Genus prompt in the Reference Genomes builder form up to the top above Species. For the Species example, make it thaliana (no mays, confusing to have more than one example). Remove all "e.g." prefixes from the example text for the Reference Genome fields. They're obviously examples from the styling. However, will screen-reader users have a way to understand these are examples? **Done**: Genus moved above Species in builder. Species label changed to "Species (specific epithet)" with placeholder "thaliana". All "e.g." prefixes removed. `aria-label` attributes added to inputs with example placeholders for screen reader accessibility.

* [x] Add a lookup using the genus and species name to get ncbi TaxId - The ncbi-datasets-cli conda package should do the trick with something like "datasets summary taxonomy taxon "Homo sapiens" | jq '.reports[0].taxonomy.tax_id'". NCBI species ID field (should rename this TaxId throughout the code base) should have a default value of "\<auto\>", but can be overridden. **Done**: Renamed `ncbiID` → `ncbi_taxid` throughout codebase. New `resolve_taxid` rule in `environment_setup.smk` uses `ncbi-datasets-cli` to auto-resolve TaxId from genus+species, with user override support. Builder shows `<auto>` placeholder. `ncbi-datasets-cli` added to `epibutton.yaml` conda env.

* [x] Derive the custom GO database name from user-inputted fields. Probably need to add the ID as a suffix to avoid collisions in the case of multiple ref genomes with the same binomial. **Done**: GO database name auto-derived as `org.<G><species>_<GenomeName>.eg.db` (e.g. `org.Athaliana_ColCEN.eg.db`). `go_database` config field removed entirely. `get_go_database()` in `RNAseq.smk` computes the name. `R_build_GO_database.R` takes `ref_genome` as 7th argument. Builder no longer shows GO database field. `samplefile_validation.py` no longer requires `go_database` or `ncbiID` for GO analysis.

  If we cannot find a match (and the option is set to auto), gracefully warn and skip the dependent analysis. If the user overrides, and their value doesn't match an actual OrgDb, then Error instead of warn.

#### Eliminate redundant requirement for GTF

* [x] Currently, users must supply both a GFF annotation file and GTF transcript annotation file. We provide instructions for deriving the latter from the former in the README.md, but we should instead simply try to create the GTF without requiring it from the user. If GTF creation fails we can raise a clear error message and ask the user to supply one explicitly and re-run. **Done**: `gtf_file` defaults to `<auto>`, auto-derived from GFF via `gffread` in `check_gtf` rule. User-provided paths still work as overrides. `gffread` added to `epibutton.yaml`. Validation no longer requires `gtf_file`. Builder shows `<auto>` placeholder and skips it on YAML export.

#### Use Infernal workflow for building structural_rna_depletion FASTA database

* [x] Current suggested approach is cumbersome and results in a FASTA database that isn't ref genome specific. If we instead just run Infernal (with lots of threads and paralellized by chromosome if on the cluster) and filter overlapping hits, determine an e-value threshold, we could incorporate this into the pipeline in the event that the user doesn't specify an override file. Research whether this approach will have significantly different outcomes than the currently prescribed approach. **Done**: `structural_rna_fafile` defaults to `<auto>`. New `download_rfam` rule fetches and presses Rfam CMs (shared across genomes). New `build_structural_rna_db` rule runs `cmscan` per chromosome/contig (parallelized via `xargs -P`; contigs < 1 Mbp binned together), filters clan overlaps, extracts FASTA via `bedtools getfasta`. Uses `--cut_ga` (Rfam gathering thresholds) by default; configurable via `infernal_threshold` option (`"ga"` or numeric E-value). `infernal` added to `epibutton.yaml`. Validation no longer requires `structural_rna_fafile`. epicc-builder updated with `<auto>` placeholder and `infernal_threshold` in sRNA options.

#### Configurable output directory

* [x] Allow specifying an output directory other than the hardcoded `results/`. **Done**: Added `output_dir` and `genome_dir` config keys with defaults `"results"` and `"genomes"`. Python context uses `RESULTS_DIR`/`GENOMES_DIR` variables; shell blocks use `{config[output_dir]}`/`{config[genome_dir]}` Snakemake substitution. 827 path references updated across 9 files via `dev/refactor_paths.py`. CLI wrapper wired with `--output-dir` and `--genome-dir` flags. Validated with pombe dry-run using both default and custom directories.

#### Provide a front-end CLI executable script

* [x] Currently, users must call snakemake to run the pipeline. Instead, we should build a project-specific executable wrapper script (called epicc) that will expose the necessary command line parameters for runtime configuration of the pipeline and calling snakemake. **Done**: `epicc` Python script at repo root with subcommands: `run`, `dry-run`, `validate`, `perf`, `unlock`, `clean`. Auto-detects SLURM vs local execution. Passthrough to snakemake via `--` separator. `--output-dir` and `--genome-dir` flags for configurable output directories. See `dev/docs/design-decisions.md` for rationale.

* [x] Add `DAG` and `rulegraph` generation options to subcommand `dry-run` and `validate` (output to png with dot). **Done**: `epicc dry-run` and `epicc validate` accept `--dag FILE` / `--rulegraph FILE`; format inferred from extension (png/svg/pdf via graphviz `dot`, or `.dot` for raw text).

* [x] Add specific command (`epicc output`) to ease the targeting of expression plots, browsers, GOs, etc. **Done**: `epicc output --plot-type TYPE --input-file PATH` with seven types: `rnaseq-histogram`, `go`, `motifs`, `srna-clusters`, `heatmap`, `metaplot`, `browser`. Per-type registry in `epicc` (OUTPUT_TYPES) maps each type to its target template, snakemake `--config` overrides, required args (`--ref-genome`, `--env`), and an input-format validator (TSV gene IDs / BED / loci). Optional `--plot-label` overrides the auto-derived label, `--analysis-name` overrides the options-file value, `--matrix {regions,tss,tes}` for heatmap/metaplot. `--dry-run` resolves the plan without executing. Unit tests in `tests/unit/test_epicc_output.py` cover label derivation, registry shape, and the format validators.

* [x] Catch snakemake/slurm warnings: "No wall time information given..." and "No SLURM account given, trying to guess. No account was given, not able to get a SLURM account via sacct: sacct: invalid option -- '1'".
  * [x] Wall time: every rule pulls `runtime=config["resources"][<rule>]["runtime"]` into its `resources:` block (109 rules across 8 `.smk` files), with a `default-resources: runtime=60` safety net in `profiles/slurm/config.yaml`.
  * [x] SLURM account: `profiles/slurm/config.yaml` now sets `slurm_account="martienssenlab"` via `default-resources`; annotated as site-specific with an `EDIT:` marker for users copying the profile to a different cluster.

* [x] maybe rename the `profile` sub-command to `perf-profile`, since profile is a specific option for snakemake. Also, does not seem to work ("Parsing .snakemake/log/2026-04-24T105453.098500.snakemake.log ... No completed jobs found in log." when some jobs were completed) **Done**: Renamed to `epicc perf` (chose `perf` over `perf-profile` for brevity and consistency with other single-word subcommands). Default mode now aggregates all logs from the same resumed run, identified by `output_dir` + `analysis_name` wildcards extracted from log content. Fixed the "No completed jobs found" bug — the parser was missing SLURM-submitted jobs (which use `Job N has been submitted with SLURM jobid M (log: .../rule_RULE/WC/...)` instead of the verbose `rule X:` block).

* [x] for 'epicc validate', capture the snakemake output silently (all the rules that would run, the files that would be removed, etc..) and only print out the diagnostic/results of the command. **Done**: `cmd_validate` now captures snakemake's dry-run output and only emits a concise summary (Job stats table + warnings + files-to-be-removed). Full output is dumped on dry-run failure or when `--verbose` is passed. New `_summarize_dryrun` helper extracts the relevant sections without parsing the per-job rule blocks.

* [x] remove 'epicc dry-run' option, which is redundant with validate (not as informative), and can be easily run with 'epicc run -- --dry-run' anyway. **Done**: Removed the `dry-run` subcommand (parser, dispatcher, and cmd_dry_run function). Header docstring updated to point users at `epicc validate` for the friendly path or `epicc run -- --dry-run` for the raw snakemake output.

* [x] Clean up all entries with personal information (cluster paths and email addresses). **Done**: Removed the lab-specific `/grid/martienssen/home/jcahn/...` paths from the B73_v5/la8011/CabSauv genome blocks in `config/epicc-options.yaml`; replaced with a generic commented-out skeleton (`MyGenome`) so users see the schema without picking up Martienssen-lab cluster paths. `profiles/slurm/config.yaml` slurm_account changed from the `martienssenlab` literal to a `<your-slurm-account>` placeholder. `profiles/geno/config.yaml` is correctly scoped as a site-specific example and untouched.

* [x] change default paths to genome files in epicc-options, to remove hardcoded paths from Elzar. Should we remove genomes that are not fetchable online, and/or add some that are? **Done**: Default options now ship one fully-fetchable working genome (ColCEN with public schatzlab URLs, used by the ColCEN integration test) plus three commented-out templates: Spombe (PomBase URLs we already validate against in the pombe integration test), TAIR10 (Ensembl Plants URLs — release-XX placeholder noted), and a generic `MyGenome` skeleton. Lab-specific Elzar paths (B73_v5, la8011, CabSauv) removed in the same configuration cleanup pass.

#### Skip mCHG/mCHH analysis for animal genomes

* [x] Animals lack asymmetric (non-CpG) DNA methylation, so mCHG and mCHH contexts contain only background false-positive calls. PCA, DMR, and browser plots for these contexts are empty or corrupt for animal genomes. **Done**: Replaced the legacy boolean `mC_context: "all"|"CG-only"` with a `methylation_contexts: ["CG", "CHG", "CHH"]` list (default for plants; animals can pass `["CG"]` to skip empty asymmetric-context outputs). New `get_methylation_contexts()` helper at the top of `mC.smk` resolves the list with backward-compat support for the legacy key (deprecation warning emitted once). Wired through `make_mc_bigwig_files` (skipped contexts get 1-bp placeholder bigwigs to keep the snakemake DAG happy), `convert_bedmethyl_to_cx_report` (CX_report filtered to active contexts via awk regex), `call_DMRs_pairwise` (R script now takes a comma-separated context list and loops), and `combined_analysis.smk` PCA gating (`PCA__mCG/mCHG/mCHH` plots emitted only when those contexts are active). `R_call_DMRs.R` rewritten with a `call_dmrs_for_context()` + `summarize_dmrs()` pair and a per-context loop, replacing the hardcoded CG-then-conditional-CHG-CHH branches. Subcontexts beyond CG/CHG/CHH (CAG, CAA, etc.) are validated against and rejected — deferred to a future PR. All four shipped test configs (chr5, colcen, dmc, mC) updated to the new key.

### Codebase Hygiene

* [x] Rename shared routines to generic names, e.g. merging_chip_replicates → merging_bam_replicates. **Done**: 8 chip-specific rule names renamed to generic names across all rule files, resource config, tests, and profiling scripts (filter_chip_pe → filter_bam_pe, make_chip_stats_pe → make_mapping_stats_pe, pe_or_se_chip_dispatch → dispatch_final_bam, merging_chip_replicates → merging_bam_replicates, prepping_chip_peak_stats → prepping_peak_stats, plotting_peaks_stats_chip_tf → plotting_peak_stats, plus SE variants).

### Performance/Resource Usage

#### Data acquisition and preparation

* [x] Switch to direct fastq.gz downloads from ENA for download speed, better transitory disk space usage? Maybe add alternative fastq_path=ENA, or try ENA first and fall back to SRA. Look into storing SRA downloads as compressed FASTQs to avoid writing huge uncompressed data to disk, and the post-hoc wait for compression. **Done**: ENA-first downloads via `workflow/scripts/ena_download.sh` with automatic fallback to fasterq-dump. ENA provides pre-compressed `.fastq.gz`, eliminating uncompressed intermediates. `fasterq-dump --temp "$TMPDIR"` uses SLURM scratch instead of `/tmp`.

* [x] it looks like mate file compression for PE SRA accessions after fasterq-dump happens serially - the R1 file must complete before R2. I don't think there's any reason for this constraint provided sufficient resources are available. I noticed this on a local run, not sure if it's true on the cluster. **Done**: PE mate compression now runs in parallel (background `&` + `wait`), each using half the available threads.

* [x] pigz appears to be limited to 4 threads for at least local runs, which can bottleneck the pipeline when there are few samples, but with a high read volume. **Done**: New `download` resource preset with 8 threads (up from 4 in `heavy`). Download rules (`get_fastq_pe`, `get_fastq_se`) now use this preset.

#### Disk usage

* [x] See if we can refactor the conda envs to decrease required disk space - current config uses ~32GiB. Are there any low-hanging dependency fruits that can be removed? Would consolidation of the majority of the analysis packages into a single environment while maintaining packages with problematic dependencies or dependency version conflicts be an optimal solution to save disk space and initial installation time by eliminating core package redundancy? **Done**: Consolidated 8 → 5 envs. Merged `epibutton_upset` into `epibutton` (removed install_github hack), `epibutton_srna` into `epibutton_rnaseq` (shared R/Bioconductor stack), `epibutton_dmc` into `epibutton_mc` (independent aligners, no conflicts). Estimated ~3-4.5 GB savings.

* [x] Should snakemake --conda-cleanup-envs be added to the pipeline to get rid of old envs? Can this be run inside of the pipeline snakemake run? **Done**: Added `conda-cleanup-pkgs: tarballs` to SLURM profile. Documented `snakemake --sdm conda --conda-cleanup-envs` as a maintenance tip in README.md (standalone command, not inside pipeline runs).

* [x] add option to keep all intermediate files, default to using pipelining and cleanup to avoid storing large intermediates like processed FASTQs, BAMs, etc. Also includes intermediate files for plotting (tracks, heatmap parameters, ...) **Done**: Added `keep_intermediates` tiered option (`none`/`standard`/`custom`/`all`) with per-category booleans (`keep_trimmed_fastqs`, `keep_final_bams`, `keep_merged_bams`, `keep_shifted_bams`, `keep_cx_reports`). `maybe_temp()` helper in Snakefile conditionally wraps outputs with `temp()`. epicc-builder exposes tier selector and custom toggles.

* [x] Review to ensure all non-retainable (via the granular retention options) intermediate files are only ever written to temp storage. **Done**: Audited all rule files; applied unconditional `temp()` to ATAC shifted BED, sRNA clean FASTQs, RNA DEG intermediates (samples/counts/RData), and GO database temp files. Added `keep_dmc_intermediates` config option for dmC aligned modBAMs and bedMethyl pileups.

* [x] Refactor merge rules to pipe `samtools merge -u` into `samtools sort`, eliminating temp intermediate BAMs. **Done**: ChIPseq `merging_chip_replicates` and RNAseq `merging_rna_replicates` now pipe directly, removing `temp_merge`/`temp` output declarations.

* [x] ChIPseq.smk - we should not be writing SAM files to disk. Wasteful of both network storage I/O (slow) and disk space. **Done**: Merged bowtie2_map + filter_chip rules into single piped rules (filter_chip_pe, filter_chip_se). Bowtie2 stdout pipes directly into samtools view → fixmate → sort, eliminating the SAM-to-disk step and reducing intermediate BAMs from 3 to 1.

#### Read trimming

* [x] Explore faster options than cutadapt, like fastp or skewer, which supports automatic Illumina adapter detection (while still allowing explicit overrides in the options file). There are a number of different use cases for read trimming - let's investigate all to make sure a substitution would cover all of them. **Done**: Replaced cutadapt with fastp across the pipeline. Config restructured: `trimming_quality` (CLI strings) replaced with tool-agnostic `quality_threshold`, `min_read_length`, `trim_front` keys. Standard Illumina adapters set to `"auto"` (fastp auto-detection); non-standard adapters (NextFlex, Nextera) kept explicit. Trimming metrics changed from `.txt` to `.json`; HTML QC reports added. All 6 downstream stats rules updated to parse fastp JSON. epicc-builder config form updated. TF keys removed.

#### Read Mapping

* [x] (ChIP/ATAC):  look at adding option to use [Chromap](https://github.com/haowenz/chromap) for ~10X speedup (and set as default), consider supporting different sensitivity levels if possible as with bt2. **Done**: Chromap added as default aligner for ChIP/ATAC (`chip_aligner`/`atac_aligner` config keys). Auto-falls back to bowtie2 for `repeat`/`repeatall` strategies (chromap lacks `-k 100` multi-mapping). MAPQ filter extracted from config instead of hardcoded. Dual-format metrics parsing in stats rules. epicc-builder updated with aligner selection UI.

#### Local

* [x] Check snakemake log times for the S. pombe integration test to profile each stage of the DAG. **Done**: Use `dev/profile_snakemake_log.py` to profile runs. Top bottleneck was `plotFingerprint` (45% CPU); addressed via auto-scaled `--numberOfSamples`, optional `chip_fingerprint_plots` toggle, and per-rep-only restriction. Also added `samtools sort -m` to `filter_bam_pe`/`se`. Full pipeline completes in ~25 min on 56 cores (down from ~42 min pre-optimization).

#### Cluster

* [x] Reconsider the resource request delineations, and maybe add more fine-grained/task-specific options like proc-intensive/himem-lowproc/mapping, etc. **Partially done**: Bismark rules now compute `tmp_mb` dynamically from input FASTQ sizes instead of using a static tier. Bismark resource tier bumped to `*max` (16 threads). `--multicore` factor changed from `threads//3` to `threads//2` for better CPU utilization. Current tier system (`low`/`standard`/`download`/`heavy`/`heavier`/`max`/`single`) covers the main use cases; further per-rule tuning can be done as profiling data comes in.

* [x] Examine wall clock times on Elzar and estimate reasonable requests with slop - most steps should be O(n+k) for sequence inputs, e.g. trimming and mapping, but downstream analysis might differ. **Done**: See runtime item below.

* [x] Make better use of temporary storage on the cluster nodes to reduce NFS I/O bottlenecks and minimize temporary disk usage bloat. **Done**: Added `precommand` to SLURM profile that sets up per-job `$TMPDIR` → `/tmp/slurm_tmp/$SLURM_JOB_ID` when `$SLURM_TMPDIR` is not already set by the cluster prolog. `fasterq-dump` already uses `${TMPDIR:-/tmp}`. Bismark temp files stay on NFS intentionally (sequential I/O pattern, not worth consuming limited local scratch).
  - **Superseded** (2026-08-19): the `shell.prefix` TMPDIR routing in `workflow/Snakefile` (`{output_dir}/.tmp/<job>.<pid>`, `use_node_tmpdir` to opt out) landed later and runs *after* the precommand, so the profile's `export TMPDIR=` had no effect on any rule. It also pointed at `/tmp`, which is the tmpfs-sized-off-`mem_mb` trap on `job_container/tmpfs` sites — the thing the shell.prefix change exists to avoid. The precommand is now commented out in `profiles/slurm/config.yaml` and documented as a site-specific escape hatch only.

* [x] ~~Small one: mem_mb and tmp_mb should be changed to mem_mib and tmp_mib to meet expectations of binary byte counting.~~ **Won't do**: `mem_mb` is a Snakemake built-in resource name used for scheduling decisions and SLURM `--mem` mapping. Renaming would break Snakemake integration. `tmp_mb` follows the same convention for consistency. SLURM interprets `--mem` and `--tmp` values as megabytes regardless of the resource name.

* [x] Add `runtime` to all resource tiers and the SLURM profile to eliminate the "No wall time information given" warning. **Done**: Added `runtime` (minutes) to all 7 resource tiers in `epicc-options.yaml` and `time: "{resources.runtime}"` to SLURM profile sbatch section. Defaults: low=60m, standard=2h, download=12h, heavy=8h, heavier=8h, max=48h, single=8h. Based on hg38 chr21 profiling data with generous margins for full-genome production runs.

* [x] Resolve slurm issues with QOSMaxSubmitJobPerUserLimit reached sometimes (when it should be limited to 16 in the profile (specific to CSHL cluster, but could be helpful for other environments in case it' a shared bug). **Done**: switched back to qos=slow_nice for all jobs.

* [x] Split workflow-specific resource tuning from cluster-specific execution settings, following Snakemake 8+ conventions. **Done**: workflow profile at `profiles/default/config.yaml` carries the per-rule `set-resources` / `set-threads` maps (167 rules covered, with named tiers via YAML anchors). `profiles/slurm/config.yaml` is now a cluster-specific *example* (executor, account, partition, precommand) with `EDIT:` markers and an empty `set-resources:` stub for site overrides. Rule bodies no longer carry `qos=` or `tmp_mb=`; the only `disk_mb=lambda` remaining is the dynamic Bismark sizing, kept in the rule body because `set-resources:` does not accept free-form lambdas. The `resources:` block has been removed from `epicc-options.yaml` (no rule references `config["resources"]`). `epicc` CLI passes `--workflow-profile profiles/default` automatically. README/RTD/builder doc updates folded into the deferred docs pass.

### Testing

#### Methylation assay types and testing

* [x] Add named assay types for non-directional WGBS (`WGBS_nd`) and PBAT (`PBAT`) with correct bismark flags (`--non_directional`, `--pbat`). **Done**: WGBS_nd and PBAT added throughout codebase (sample_sheet.py, mC.smk, config, builder, tests). See design-decisions.md for rationale.

* [x] Create comprehensive mC dry-run test covering all methylation code paths. **Done**: `tests/integration/test_mC_dryrun.py` (52 tests) validates WGBS, WGBS_nd, PBAT, EMseq, dmC modBAM, and dmC bedMethyl workflows including parameter routing, replicate merging, DMR calling, and CX report conversion.

#### Schizosaccharomyces pombe test case

* [x] Add S. pombe integration test for faster development, user installation validation, and local single-host execution as well as cluster execution. **Done**: 17 samples (10 ChIP, 4 RNA-seq, 3 sRNA), 259 pipeline steps, ~1h 11m on gemmule with 56 threads. See `tests/integration/data/test_options_pombe.yaml` and `tests/integration/data/pombe_design.md`.

* [x] Gather all necessary genome reference resources (fasta, gff, gtf) from [Pombase.org](https://www.pombase.org/monthly_releases/2026/pombase-2026-02-01/). Derive an appropriate test config and test samplefile. **Done**: PomBase Feb 2026 FASTA/GFF3, gffread-derived GTF, Infernal/Rfam-15.0 structural RNA FASTA (261 loci). Files in `tests/integration/data/Spombe/`.

* [x] Let’s use Hyun Soo Kim’s ChIP-seq (H3K9me2, H3K9me3, sRNA) + Ekwall lab H3K4me3 + Chang/Rct1 WCE + Martienssen 2025 RNA-seq:
  * Kim et al. 2024 GSE156069: H3K9me2/me3 WT+dcr1 PE ChIP, sRNA SE
  * Ekwall lab GSE280066: H3K4me3 SE ChIP + Input
  * Chang et al. 2017 GSE97746: WCE control PE
  * Martienssen 2025 GSE278839: RNA-seq WT+dcr1 SE

* [x] Search only R.A. Martienssen publication datasets for additional RNA-seq and ChIP libraries. **Done**: used 3 Martienssen lab datasets + 1 Ekwall lab dataset.

* [x] Any necessary data caching for the development of this test case should be done in the untracked test-data-prep directory. **Done**: structural RNA build intermediates in `test-data-prep/pombe-structural-rna/`.

* [x] Add targeted rule tests executing real samtools commands on synthetic SAM data (filter_chip_pe, filter_chip_se, filter_rna_se, bamCoverage). **Done**: `tests/unit/test_rule_commands.py` (14 tests, 12 pass + 2 skip when bamCoverage unavailable).

* [x] Add post-run validation tests checking all pipeline outputs for existence, format integrity, and correct structure. **Done**: `tests/integration/test_pombe_postrun.py` (29 tests, marked `@slow`, auto-skips if no completed run).

* [x] Add validate-pombe orchestration script and Claude Code skill. **Done**: `scripts/validate_pombe.sh` (`--dry`/`--full`/`--check`/`--all`), `.claude/commands/validate-pombe.md` (`/validate-pombe` skill).

* [x] Remove tissue:cell from test_samples_pombe sample sheet - no reason to have this factor as it is not differentiated among the samples. **Done**: Reduced Levels to single `genotype:X` factor, updated all Sample_IDs, analysis names, and test assertions across 5 files.

* [x] Fix pombe test sample metadata: rename Ekwall `veg` → `WT` (HU3112 is WT), fix Kim/Chang WCE controls from `IP_target: Input` → `WCE` (SRR5445712 is whole-cell extract, not Input). **Done**: Updated test_samples_pombe.tsv, test_pombe_dryrun.py, test_pombe_postrun.py, test_sample_sheet.py.

#### Complete A. thaliana ColCEN test case

* [x] Remove files associated with the Chr5 subset. Let's just target the full ColCEN genome. **Done**: `git rm -r tests/integration/data/ColCEN_Chr5/`, `test_options_chr5.yaml`, `test_samples_chr5.tsv`.

* [x] Complete the A. thaliana ColCEN test case: **Done**: 38-sample test sheet with SRA/DDBJ accessions for ChIP (CenH3, H3K9me2), ATAC, EMseq, and PBAT; lemna.org URLs for dmC modBAMs. Self-contained `test_options_colcen.yaml` with GitHub URLs for genome reference (Col-CEN v1.2 FASTA, GFF3, EDTA TE GFF3). Required pipeline enhancements: URL support for Read_files and genome config fields, GFF3→BED auto-conversion for TE file, DRR/ERR accession support. See `tests/integration/data/colcen_design.md` for full dataset provenance.

  * [x] Find the corresponding SRA entries for the [Shimada Nat Plants 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11410651/) ChIP samples (CENH3, Input in genotypes WT, rdr126ddm1, rdr126ddm1hp5), GEO series: [GSE132005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132005) and replace local file paths with these. The ONT methylation samples from the same GEO series are also in SRA (mislabeled as bisulfite). They should've been deposited as modBAMs, need to verify this. **Done**: ChIP CenH3 SRR28453410-21 (12 samples), H3K9me2 SRR28453402-09 (8 samples), dmC modBAMs via lemna.org URLs.

  * [x] Find A. thaliana PBAT datasets for the ColCEN test case. **Done**: 4 PBAT PE samples (Col-0 WT × 2 reps, tarp1tarp2 × 2 reps) from [Takei et al. 2024 Plant Physiology](https://academic.oup.com/plphys/article/195/2/1333/7623186) PRJDB14218 (DRR400324-27). 10-day seedlings. Required expanding SRA accession regex to accept DRR/ERR prefixes. WGBS_nd (Pico Methyl-Seq) still TBD.

  * [x] Add A. thaliana EMseq data from [Trasser et al. 2024 EMBO reports](https://pmc.ncbi.nlm.nih.gov/articles/PMC11624286/) to the ColCEN test case. **Done**: 4 EMseq PE samples (Col-0 WT × 2 reps, rdr6 × 2 reps) from PRJNA1111825 (SRR29036835-40). Inflorescence tissue.

  * [x] Document the provenance of all datasets as previously done with S. pombe and H. sapiens test cases. **Done**: `tests/integration/data/colcen_design.md` and `tests/integration/data/pombe_design.md`.

#### H. sapiens test case

* [x] Add H. sapiens Chr21 test case.

#### Expand rule-level dry-run tests

* [x] Expand dmC dry-run tests to cover the full `mC.smk` rule (all bismark workflows + dmC). **Done**: `tests/integration/test_mC_dryrun.py` (52 tests) covers WGBS, WGBS_nd, PBAT, EMseq, dmC modBAM, and dmC bedMethyl including parameter routing, replicate merging, DMR calling, and CX report conversion. The original `test_dmc_dryrun.py` was merged into this and removed.
* [x] Add similar lightweight dry-run test modules for other rule files (ChIP, RNA, sRNA, ATAC) using mock inputs and a fake genome, following the mC test pattern. **Done**: `test_ChIP_dryrun.py` (22 tests), `test_ATAC_dryrun.py` (24 tests), `test_RNA_dryrun.py` (27 tests), `test_sRNA_dryrun.py` (24 tests). PBAT/EMseq already covered in `test_mC_dryrun.py`. Combined_analysis and folding dmc into mC still TBD.

#### Add test dataset documentation

* [x] Add a test design doc similar to hg38_chr21_design.md for S. pombe and A. thaliana test cases. **Done**: `tests/integration/data/pombe_design.md` and `tests/integration/data/colcen_design.md`.

### Known Unknowns

* [x] ~~For now, ChIPseq replicates are only properly merged if same paired information (all PE or all SE). Not sure what happens if both PE and SE reps are available with the same line+tissue name.~~ Mixed PE/SE analysis groups force SE treatment for merged peak calling with a startup warning. Per-replicate peaks use each sample's own layout. See `dev/docs/design-decisions.md`.

### Known Issues/Bugs

* [x] PlotPCA can fail if no dimensions found. check npz results before starting PCA? **Done**: `plot_PCA_correlation` rule now catches `plotPCA` failures (insufficient data for PCA) and creates a placeholder output instead of failing the pipeline.

* [x] ~~Chromap dropping >90% reads from ColCEN PE test data.~~ **Done** (commit `c78cd8f`): root cause was `samtools view -q` running before `samtools fixmate -m` in the chromap PE pipe. Filtering on MAPQ first dropped one mate of a pair while keeping the other, desyncing the name-collated stream so fixmate then mis-paired adjacent records. Reordered to `view -F → fixmate -m → view -q → sort → markdup`. Confirmed against the colcen run: CenH3 PE 47.4M → 20.4M (57% loss, mostly multi-mappers + duplicates), comparable to bowtie2 on the same data.

* [x] ~~Chromap SE output not getting properly piped/handled downstream by samtools (see ColCEN test).~~ **Done** (commit `c78cd8f`): a stray `samtools fixmate -m` had been added to the SE chromap branch in `5016558` to handle "chromap sets paired flags on some SE reads". `c78cd8f` removed it (PE-only operation; markdup works on coord-sorted SE input directly). Confirmed clean SE flag distribution on the colcen run (50/50 flag 0/16, no paired-flag contamination).
