# EPICC Design Decisions

A reference for architectural choices in the pipeline. Intended for contributors who need to understand why things are structured the way they are. Update this file when making major design decisions.

---

## Sample Sheet Design

### Sample_ID as the primary filesystem identifier

**Decision**: Each row in the sample sheet carries a user-assigned `Sample_ID` that serves as the filename stem for all per-replicate outputs (e.g. `final__WT_H3K9me2_rep1.bam`).

**Rationale**: The previous scheme constructed a compound name from six fixed columns (`data_type`, `line`, `tissue`, `sample_type`, `replicate`, `ref_genome`) and then deconstructed it with regex throughout the pipeline. This created fragility: any field containing a double-underscore or unexpected character could corrupt the name, and adding a new field required changing every rule that parsed the name. Using `Sample_ID` as an opaque token eliminates construction/parsing entirely. Wildcards in rule patterns bind directly to Sample_ID.

**Alternatives considered**: Keeping the compound name but centralizing its construction and parsing in one module. Rejected because it still requires all six fields to be present for every assay type (many fields are meaningless for RNA-seq or WGBS), and regex-based parsing is inherently fragile.

### Analysis-level names are derived, not stored

**Decision**: Merged-replicate (analysis-level) filenames are computed from group fields as `{Assay}__{levels_label}__{IP_target}__{Genome}`, with empty parts omitted so non-ChIP assays do not produce `____` in their paths.

**Rationale**: Analysis-level names need to be deterministic and human-readable, but they are fully derivable from sample metadata and should never be specified by the user. Deriving them keeps the sample sheet minimal while ensuring that replicates with identical metadata always resolve to the same analysis group.

**Alternatives considered**: Storing an explicit analysis group column. Rejected because it would require users to keep this column consistent with the other metadata fields, creating a validation burden and a new source of user error.

### `Levels` column for arbitrary experimental factors

**Decision**: Experimental conditions are encoded as comma-separated `factor:level` pairs in a single `Levels` column (e.g. `genotype:WT,tissue:root`). All samples in a run must have the same factor names and count.

**Rationale**: The old format had fixed `line` and `tissue` columns, which could not represent designs with a single factor, three or more factors, or factor names other than genotype and tissue. The `Levels` approach is schema-free while remaining machine-parseable. The factor names improve readability and allow the epicc-builder app to render dynamic per-factor columns for data entry.

**Alternatives considered**: Multiple fixed columns with a naming convention (e.g. `Factor1`, `Factor2`). Rejected because it shifts the schema negotiation problem from the pipeline to the user while providing no additional flexibility.

### Explicit control linking via a `Control` column

**Decision**: The `Control` column in each IP sample row must contain the `Sample_ID` of its corresponding input or whole-cell extract sample. Control chaining (a control pointing to another control) is forbidden.

**Rationale**: The previous pipeline assigned controls by matching `sample_type == "Input"` within the same group. This made it impossible to use different controls for different marks, or to use a whole-cell extract (WCE) with a name other than "Input". Explicit linking also means the relationship is self-documenting and survives any renaming of samples.

**Alternatives considered**: Keeping implicit matching by name convention. Rejected because it cannot handle multiple control types in the same experiment (Input, WCE, IgG) or different controls per mark.

### Peak type from an explicit column, not from the IP_target name

**Decision**: Peak type comes from the `Peak_type` column (`broad` for histone-domain marks, `narrow` for transcription factors and sharp marks). The IP_target name is not inspected to determine it.

**Rationale**: The original approach pattern-matched `IP_target` strings against a regex list in the options file. This required users to maintain the regex list and could silently mis-classify marks not in it. An explicit column makes the decision explicit at sample-sheet entry time and eliminates the regex config entirely.

**Alternatives considered**: (1) Retaining the regex approach. Rejected because it is error-prone for novel marks and obscures the setting from the sample sheet. (2) Encoding it in the Assay token (`ChIP_broad`), which is what shipped first. Rejected because it conflates an analytical parameter with an experimental invariant: the assay records what happened at the bench and cannot change, while peak type is a choice that can be revisited on the same libraries.

### No backward compatibility with old-format sheets

**Decision**: The sample-sheet format is a clean break from the old one; old-format detection was removed from the validation code. A migration script was provided for a time and then dropped ahead of 1.0.0 — see "No legacy assay tokens" below.

**Rationale**: Maintaining dual-format support would require the validation and parsing code to branch on format detection, and the old format could not be reliably detected by column names (old sheets did not always include headers).

---

## Conda Environment Strategy

### Consolidate to the minimum number of environments consistent with dependency conflicts

**Decision**: Eight conda environments were reduced to five by merging environments that shared an R/Bioconductor stack and had no package version conflicts.

| Env file | Serves |
|---|---|
| `epibutton.yaml` | base utilities, combined analysis, upset plots |
| `epibutton_chip.yaml` | ChIP-seq and ATAC-seq |
| `epibutton_rnaseq.yaml` | RNA-seq and small RNA-seq |
| `epibutton_mc.yaml` | bisulfite mC and direct methylation (dmC) |
| `epibutton_idr.yaml` | IDR only |

**Rationale**: The original eight environments contained significant redundancy: `r-base` installed three times, `bedtools` five times, `samtools` four times, and the common R visualization stack (dplyr, ggplot2, RColorBrewer) four or five times. The estimated total footprint was ~32 GB. Consolidating environments with compatible dependency trees reduces the installation footprint by roughly 3-4.5 GB and shortens initial setup time.

**Alternatives considered**: Merging all analysis tools into a single environment. Rejected because IDR requires `numpy <= 1.19` (due to deprecated `np.int`), which would conflict with numpy versions required by other tools. Merging ChIP tools into the base environment was also considered but rejected because bowtie2 + MACS2 + MaNorm represent a clean functional grouping with minimal overlap gain.

### IDR environment remains isolated

**Decision**: `epibutton_idr.yaml` is kept as a standalone environment pinned to `numpy <= 1.19`.

**Rationale**: The `idr` package depends on a numpy API (`np.int`) removed in numpy 1.20. Patching IDR upstream or using a compatibility shim was considered less portable than simply isolating the environment. Since IDR is run in a single rule, the isolation cost is low.

### `r-remotes::install_github` hack eliminated

**Decision**: The R upset plot scripts previously called `remotes::install_github("krassowski/complex-upset#212")` at runtime to work around a ComplexUpset incompatibility with ggplot2 4.0. This was removed; the conda-distributed `r-complexupset` package already incorporates the fix.

**Rationale**: Runtime installation of R packages inside a Snakemake rule is fragile (requires internet access from compute nodes, fails on air-gapped clusters, and makes the environment non-reproducible). The workaround was only necessary because the conda package lagged behind the fix; once the conda package was updated, the hack had no remaining purpose.

### Conda package tarball cleanup in SLURM profile

**Decision**: `conda-cleanup-pkgs: tarballs` is set in `profiles/slurm/config.yaml`. Orphaned environment cleanup (`--conda-cleanup-envs`) is documented as a manual maintenance step, not automated inside pipeline runs.

**Rationale**: Cleanup of orphaned environments changes which environments Snakemake considers valid and cannot be safely run concurrently with an active pipeline execution. Tarball cleanup is safe to automate because it only removes cached source archives after installation is complete.

---

## Mixed PE/SE Replicates in ChIP-seq

### PE reads treated as SE when an analysis group has mixed layouts

**Decision**: When an analysis group (same Assay + Levels + IP_target + Genome) contains both PE and SE replicates, the merged/pseudoreplicate peak calling forces SE treatment (`-f BAM` instead of `-f BAMPE`). Per-replicate peak calls still use each sample's own layout. A warning is printed at startup.

**Rationale**: Merging PE and SE BAMs produces a mixed BAM. MACS2 with `-f BAMPE` on such a BAM crashes (`fragment size = 0.0`, `ZeroDivisionError`) because SE reads have no mate pair information. The SE format (`-f BAM`) treats each read as an independent tag extended from its 5' end by the estimated fragment size `d` -- for PE reads this means both mates are treated as separate SE reads rather than using the actual insert size as fragment boundaries. This is slightly suboptimal for PE data but perfectly functional for peak detection. Per-replicate calls retain full PE information where available, so IDR comparison still benefits from PE fragment resolution on PE replicates.

**Alternatives considered**: (1) Rejecting mixed PE/SE at validation time. Rejected because it's a legitimate, if uncommon, scenario (e.g. ENCODE experiments from different eras). (2) Stripping pairing flags before merging. Rejected as unnecessary -- MACS2 `-f BAM` already handles both read types correctly. (3) Separate merged peak calling per layout. Rejected as over-complex for a rare edge case.

---

## Performance and Disk Usage

### Configurable intermediate file retention with tiered presets

**Decision**: The `keep_intermediates` options key accepts four tiers: `none`, `standard`, `all`, and `custom`. The default is `standard`, which keeps per-replicate final BAMs and deletes trimmed FASTQs, merged BAMs, ATAC shifted BAMs, and bisulfite CX reports. A `maybe_temp()` helper in the Snakefile conditionally wraps output paths with Snakemake's `temp()` based on the resolved tier.

**Rationale**: Large intermediate files (trimmed FASTQs can reach 50 GB each; merged BAMs are also substantial) accumulate quickly in multi-sample experiments. Snakemake's `temp()` mechanism handles automatic deletion after downstream rules consume the file, but it was previously all-or-nothing. Users running on NFS-backed HPC storage have strong incentives to minimize disk footprint, while users doing exploratory analysis may want to retain intermediates for inspection. Tiered presets provide sensible defaults without requiring per-file configuration.

**Per-file defaults rationale**:
- Trimmed FASTQs deleted by default: re-trimming from raw data is fast; raw FASTQs are the authoritative inputs and are typically archived separately.
- Per-replicate final BAMs kept by default: these are frequently loaded in IGV or used for custom downstream analysis outside the pipeline.
- Merged BAMs and shifted ATAC BAMs deleted by default: consumed entirely by downstream bigwig/peak rules; large and easily regenerated.

**Alternatives considered**: A single boolean `keep_all_intermediates`. Rejected because it provides no middle ground — the most common need (keep BAMs, discard FASTQs) cannot be expressed.

### SAM-to-disk eliminated in ChIP and ATAC mapping

**Decision**: The `bowtie2_map` and `filter_chip` rules were merged into single piped shell commands. Bowtie2 writes to stdout, which is piped through `samtools view`, `samtools fixmate`, and `samtools sort` before producing the first BAM written to disk.

**Rationale**: SAM format is uncompressed and roughly 3-5x larger than an equivalent BAM. Writing a SAM to NFS storage and immediately reading it back for conversion is wasteful in both I/O and disk space. On high-latency NFS mounts, eliminating this round-trip provides a measurable speedup. The piped approach also reduces the number of intermediate BAMs from three to one.

**Alternatives considered**: Keeping separate rules but marking the SAM output as `temp()`. Rejected because the SAM file would still be written and read through NFS storage; the file simply gets deleted afterward. Piping avoids the I/O entirely.

### Merge rules pipe `samtools merge -u` into `samtools sort`

**Decision**: BAM merging rules for ChIP and RNA-seq pipe uncompressed merge output directly into sort, eliminating a temporary merged BAM on disk.

**Rationale**: Intermediate unsorted merged BAMs are large (sum of all replicate BAM sizes) and immediately consumed by sort. Writing them to disk provides no benefit and doubles the transient disk footprint of the merge operation.

### ENA-first SRA downloads with fasterq-dump fallback

**Decision**: SRA download rules (`get_fastq_pe`, `get_fastq_se`) first attempt to download pre-compressed `.fastq.gz` files from ENA via HTTPS, falling back to `fasterq-dump` only on failure. The ENA URL is constructed deterministically from the accession using a helper script (`workflow/scripts/ena_download.sh`). The fasterq-dump fallback uses `--temp "${TMPDIR:-/tmp}"` for scratch storage, and PE mate compression runs in parallel.

**Rationale**: `fasterq-dump` writes uncompressed FASTQs (often 10-50 GB per mate) to disk, then `pigz` compresses them sequentially. ENA provides the same data pre-compressed, eliminating both the uncompressed intermediate and the compression step. This reduces peak disk usage by 2-3x and removes the I/O bottleneck of writing and reading large uncompressed files on NFS storage. The fallback ensures robustness: ENA may be unavailable, or accessions may not yet be mirrored (ENA mirrors lag SRA by hours to days). Using `$TMPDIR` for fasterq-dump scratch avoids filling shared `/tmp` and makes use of SLURM-allocated node-local storage. A dedicated `download` resource preset (8 threads) replaces the `heavy` preset (4 threads) for download rules, since both ENA downloads and pigz compression are I/O-bound and benefit from additional threads.

**Alternatives considered**: (1) Aspera (ascp) protocol for faster ENA transfers. Rejected because the IBM Aspera client is a non-trivial dependency that many HPC sites do not install. (2) ENA file-report API for URL/checksum validation. Rejected because URL construction from accession numbers is deterministic and well-documented; a failed curl cleanly signals fallback. (3) Always using fasterq-dump with `--include-technical` and piping directly to gzip. Rejected because fasterq-dump does not support stdout output, so the uncompressed intermediate cannot be avoided.

---

## Testing Strategy

### S. pombe as the primary integration test organism

**Decision**: A full multi-assay integration test uses publicly available S. pombe data (18 samples: ChIP-seq H3K9me2/me3/H3K4me3, RNA-seq, sRNA-seq) rather than a synthetic or downsampled A. thaliana dataset.

**Rationale**: S. pombe has a small genome (~12.5 Mb, 3 chromosomes) that makes full end-to-end pipeline runs practical on a laptop or single workstation in approximately 70 minutes. Using real published data from established labs (Kim et al. 2024 GSE156069, Ekwall lab GSE280066, Chang et al. 2017 GSE97746, Martienssen 2025 GSE278839) ensures the test exercises realistic data characteristics — mapping rates, peak shapes, differential expression patterns — rather than idealized synthetic data.

**Alternatives considered**: Subsetting an Arabidopsis dataset to a single chromosome. Rejected because chromosome-subset data introduces edge cases in reference genome indexing and annotation processing that do not reflect normal usage, and the Arabidopsis genome is larger enough that even a single chromosome produces a slower test than the complete pombe genome.

### Three-tier test structure: unit / dry-run / post-run

**Decision**: Tests are organized into unit tests (`tests/unit/`), integration dry-run tests (`tests/integration/test_pombe_dryrun.py`), and post-run validation tests (`tests/integration/test_pombe_postrun.py`). Post-run tests are marked `@slow` and auto-skip when no completed pipeline run is present.

**Rationale**: Dry-run tests (Snakemake DAG resolution without executing rules) catch the majority of configuration and wiring errors quickly without requiring compute resources or data downloads. Unit tests for `sample_sheet.py` functions catch logic errors in the parsing and validation layer without invoking Snakemake. Post-run tests verify output file existence, format integrity, and expected content after a real execution. Separating these tiers allows rapid iteration during development while preserving thorough validation for release.

### Rule command tests use synthetic SAM data

**Decision**: `tests/unit/test_rule_commands.py` constructs minimal synthetic SAM records in memory, runs the actual samtools pipeline commands extracted from rule shell blocks, and validates the BAM output.

**Rationale**: The most failure-prone rules after the SAM-to-disk elimination are the piped alignment rules (filter_chip_pe, filter_chip_se, filter_rna_se). These are difficult to test via dry-run because Snakemake cannot verify that a shell pipeline produces correctly sorted, filtered output. Testing the shell commands directly on synthetic data catches samtools flag errors, sort order issues, and pipe failures without requiring a real alignment run.

---

## Pipeline Architecture

### Sample-sheet logic centralized in `workflow/scripts/sample_sheet.py`

**Decision**: All functions for reading, validating, and querying sample metadata live in a single module (`sample_sheet.py`). The Snakefile imports from this module and calls a small set of well-defined functions rather than implementing parsing logic inline.

**Rationale**: Before the refactor, approximately 130 lines of sample-sheet parsing logic were embedded in the Snakefile, interleaved with rule definitions. This made the logic untestable in isolation and difficult to reason about. Centralizing it in an importable module enables unit testing of every parsing and validation function against in-memory DataFrames, independent of Snakemake.

### Modular rule files organized by assay type

**Decision**: Rules are split across eight `.smk` files, one per assay type, plus `environment_setup.smk` and `sample_download.smk`. The main Snakefile orchestrates inclusion and does no analysis work itself.

**Rationale**: A monolithic Snakefile becomes unnavigable as the number of assays grows. Per-assay rule files make it straightforward to locate rules, understand dependencies within an assay, and modify one assay's logic without risk of affecting others. The main Snakefile serves as an index and configuration point.

### Results directory organized by environment, not by assay

**Decision**: Output files go to `results/{env}/` where `env` is a coarser grouping than assay: `ChIP` serves both `ChIP_broad` and `ChIP_narrow`; `RNA` serves both `RNAseq` and `RAMPAGE`; `mC` serves `WGBS`, `WGBS_nd`, `PBAT`, `EMseq`, and `dmC`.

**Rationale**: Assays within the same environment share tools, conda environments, and downstream analysis steps (peak calling for both broad and narrow ChIP, differential expression for both RNAseq and RAMPAGE). Grouping their outputs together simplifies combined analysis rules that aggregate across replicates or marks, and reduces the number of top-level directories a user needs to navigate.

### Checkpoint files for analysis-level rerun control

**Decision**: Analysis-level steps (peak calling, differential expression, DMR calling, combined plots) write empty sentinel files to `results/{env}/chkpts/`. Deleting a checkpoint file forces that analysis to rerun on the next pipeline invocation.

**Rationale**: Snakemake re-runs rules when inputs are newer than outputs. For analysis-level rules that aggregate across many samples, the inputs include all mapped BAMs, all peak files, and configuration parameters. Any upstream change — adding a new sample, adjusting a parameter — would otherwise trigger rerunning all downstream combined analysis. Checkpoint files decouple the re-run decision from file timestamps, giving users explicit control over which analyses to regenerate.

### epicc-builder as a self-contained, offline HTML5 application

**Decision**: The sample sheet and options file builder is a single HTML file (`tools/epicc-builder.html`) that can be opened directly in any modern browser without a server, internet connection, or Python environment. It is also deployed to GitHub Pages for users who prefer not to download the file.

**Rationale**: The original epicc-builder was a Streamlit web application that required a running Python server and became broken when the hosting environment changed. A self-contained HTML file has no runtime dependencies, works on air-gapped HPC login nodes, and can be committed to the repository and versioned alongside the pipeline code it describes. Tabulator (a JavaScript library that can be bundled into a single file) provides the interactive table functionality.

**Alternatives considered**: Rebuilding as a Streamlit or Flask app. Rejected because server-side apps require deployment infrastructure and introduce a runtime dependency that can fail independently of the pipeline. A CLI tool was also considered but rejected because tabular data entry is substantially more usable with a graphical interface.

---

## Configuration Naming

### "Options file" replaces "config file" for the pipeline YAML

**Decision**: The pipeline's main YAML configuration file was renamed from `config/config.yaml` to `config/epicc-options.yaml`, and all references to "config file" in documentation, code comments, and error messages were changed to "options file". Test config files were similarly renamed (e.g. `test_config_pombe.yaml` → `test_options_pombe.yaml`). The `config/` directory name, Snakemake's `--configfile` CLI flag, and profile config files (`profiles/slurm/config.yaml`) were left unchanged.

**Rationale**: "Config" is overloaded in the Snakemake ecosystem — the pipeline's main YAML, the SLURM profile YAML, and Snakemake's own configuration all use `config.yaml`. Calling the pipeline's YAML the "options file" clarifies that a complete run configuration is the composite of a sample sheet and an options file, and eliminates confusion with profile configs. The `epicc-` prefix in the filename makes the file self-identifying when it appears outside the repository context.

**Alternatives considered**: Renaming to `pipeline.yaml` or `settings.yaml`. Rejected because "options" better conveys user-tunable parameters, and the `epicc-` prefix ties it to the pipeline identity. Renaming the `config/` directory was considered but rejected because it would break Snakemake's default config path resolution and add unnecessary churn.

---

## Read Trimming

### fastp replaces cutadapt for adapter removal and quality trimming

**Decision**: Replace cutadapt with fastp across all assay types (ChIP, ATAC, RNA, sRNA, mC). The config structure was changed from tool-specific CLI strings (`trimming_quality: "-q 10 -m 20"`) to tool-agnostic keys (`quality_threshold`, `min_read_length`, `trim_front`). Standard Illumina adapters are set to `"auto"` for fastp auto-detection; non-standard adapters (NextFlex for sRNA, Nextera for ATAC) remain explicit. Trimming metrics output changed from cutadapt `.txt` to fastp `.json`, and HTML QC reports are now generated alongside.

**Rationale**: fastp is 2-5x faster than cutadapt, actively maintained (v1.0 released 2025), and provides automatic adapter detection for both SE and PE reads, built-in HTML/JSON QC reports, and polyG tail trimming. The tool-agnostic config keys decouple the configuration from any specific trimming tool, making future tool changes easier. fastp also applies default read-level quality filtering (discard reads with >40% bases below Q15) that cutadapt does not have, which is generally beneficial for catching truly low-quality reads. The `trim_front` key replaces the previously unimplemented manual Pico/EMseq override comment with a functional config option.

**Alternatives considered**: (1) skewer: unmaintained since ~2016, lacks QC reporting, would require a separate QC tool. (2) Keeping cutadapt with parallelized JSON output: cutadapt lacks built-in JSON metrics and auto-detection, so the benefit would be marginal. (3) trim_galore (a cutadapt wrapper): adds a dependency layer without significant speed improvement over cutadapt itself.

### Chromap as default aligner for ChIP-seq and ATAC-seq

**Decision**: Add Chromap as the default aligner for ChIP-seq and ATAC-seq, with automatic fallback to bowtie2 for multi-mapping strategies. Options keys `chip_aligner`/`atac_aligner` select the aligner (default `"chromap"`); `chip_mapping_strategy`/`atac_mapping_strategy` replace the old `*_mapping_option` keys. When chromap is selected but the mapping strategy is `repeat` or `repeatall`, the pipeline automatically falls back to bowtie2 because chromap lacks bowtie2's `-k 100` multi-mapping mode. The MAPQ filter is extracted from the existing `bt2_mapping_strategy` config rather than hardcoded, making it aligner-agnostic. Stats rules auto-detect the metrics format (bowtie2 vs chromap stderr) for parsing.

**Rationale**: Chromap is a minimizer-based aligner specifically designed for ChIP-seq and ATAC-seq data. Benchmarks show ~99.8% peak concordance with bowtie2 at 10-40x the speed. For large-genome organisms (e.g. maize at 2.2 Gb), this translates to substantial wall-clock savings on both local and cluster runs. Chromap is available on bioconda and actively maintained. The automatic fallback mechanism means users get the speed benefit by default while retaining bowtie2's multi-mapping capability when needed.

**Alternatives considered**: (1) Replacing bowtie2 entirely: rejected because chromap cannot replicate `-k 100` multi-mapping mode needed for repeat analysis. (2) Using chromap's BED output for ATAC (with `--preset atac` and Tn5 shift): rejected because downstream bigwig generation requires BAM input, and the existing `atac_shift_bam` rule (deeptools `alignmentSieve --ATACshift`) handles the shift correctly on BAM. (3) Adding chromap sensitivity tuning via `-e` (error threshold): deferred as a future enhancement since chromap has no preset sensitivity modes analogous to bowtie2's `--very-sensitive`.

### Named assay types for bisulfite library prep methods

**Decision**: Non-directional WGBS (Zymo Pico Methyl-Seq, Swift Accel-NGS, etc.) and PBAT libraries are represented as distinct Assay types (`WGBS_nd`, `PBAT`) rather than parameterized via a generic `directionality` field on the existing `WGBS` type. Each assay type maps to a named entry in the `mC_mapping` config block, which specifies the bismark alignment and extraction flags.

| Assay type | Bismark alignment flag | `--ignore_r2 2` | Directionality |
|---|---|---|---|
| `WGBS` | (none) | Yes (end-repair artifact) | Directional (strands 1+2) |
| `WGBS_nd` | `--non_directional` | No | All 4 strands |
| `PBAT` | `--pbat` | No | Strands 3+4 only |
| `EMseq` | (none) | No | Directional (strands 1+2) |

**Rationale**: Users think in terms of library prep method, not bismark flags. The set of bisulfite library preps is bounded (directional WGBS, non-directional WGBS, PBAT, EMseq) and unlikely to grow. Named types make the sample sheet self-documenting, eliminate the need for users to understand bismark's `--non_directional` vs `--pbat` distinction, and allow the config to provide per-type defaults for all mapping and extraction parameters. The `mC_mapping` config block already supports per-type overrides, so adding new types requires no structural changes.

**Alternatives considered**: (1) A `directionality` field in the sample sheet (values: `directional`, `non_directional`, `pbat`). Rejected because it would be empty for all non-mC samples, adding column noise; it conflates a pipeline parameter with sample metadata; and it requires users to know which directionality their library prep produces. (2) A generic `WGBS` type with a per-sample `library_prep` annotation. Rejected for similar reasons — it pushes pipeline routing decisions into a secondary field rather than the primary Assay type that already determines all other routing.

### Genome config namespaced under `genomes:` with inlined species parameters

**Decision**: All reference genome blocks in `epicc-options.yaml` are namespaced under a `genomes:` key. Species-level parameters (`star_index`, `genomesize`, `ncbiID`, `genus`, `go_database`) that were previously in separate top-level species blocks are inlined into each genome entry. Rule files access genome config via `config["genomes"][ref_genome][field]`. A backward-compatibility shim in the Snakefile auto-migrates old bare-key format configs at startup with a deprecation warning. Startup validation (`check_genome_config()`) verifies that every genome referenced in the sample sheet has the required fields for its assay types.

**Rationale**: The previous layout placed genome blocks as bare top-level keys alongside unrelated config (e.g. `ColCEN:` at the same level as `quality_threshold:`). Species parameters required a fragile two-hop lookup pattern: `config[config[ref_genome]['species']][field]`. This was error-prone (a typo in the `species` field silently broke all downstream lookups), hard to validate at startup, and confusing for new users. Namespacing under `genomes:` makes the config self-documenting, inlining eliminates the indirection, and startup validation catches missing fields before rule execution.

**Alternatives considered**: (1) Keeping the two-level genome + species separation but namespacing both. Rejected because the indirection provided no real benefit — each genome already needed a unique species mapping, and most species blocks were used by exactly one genome. (2) Removing the `species` field entirely. Rejected because it still serves a documentation purpose and is used in the GO database naming convention.

### Auto-computed genomesize and STAR index from reference FASTA

**Decision**: `genomesize` (MACS2 `-g` effective genome size) and `star_index` (`--genomeSAindexNbases`) are computed automatically from the reference FASTA by a `compute_genome_stats` rule in `environment_setup.smk`. The rule produces `genomes/{ref_genome}/genome_stats.json` containing total bases, N bases, effective size, and STAR SA index nbases. Consumer rules (MACS2 peak calling, STAR index building) read the computed values from the JSON file at execution time. User-provided values in the options file take precedence (override check in shell blocks via params).

**Formulas**:
- `genomesize` = total bases − hard-masked (N/n) bases. This is appropriate for MACS2 `-g`, which needs the effective (non-N) genome size. Softmasked (lowercase) bases are NOT subtracted: they represent mappable repeats and should be included in the effective size.
- `star_index` = `min(14, floor(log2(totalBases) / 2 - 1))` per the STAR manual recommendation.

**Rationale**: Both values are fully derivable from the reference FASTA. Requiring users to specify them manually was error-prone (the S. pombe test case had an incorrect `star_index` of 11 instead of the correct 10) and created an unnecessary barrier for adding new genomes. The override mechanism preserves backward compatibility and allows expert users to tune values.

**Design**: The stats are computed in a dedicated rule (not a params lambda) because Snakemake params lambdas run at DAG time before execution, so they cannot read files produced by other rules on a fresh run. The `genome_stats.json` file is added as an input to consumer rules, ensuring proper dependency ordering. The shell block checks for a config override first (via a params lambda using `.get()`), falling back to the JSON file read only when no override is provided.

**Alternatives considered**: (1) Computing uniquely-mappable-regions for genomesize. Rejected because it is computationally heavyweight, read-length-dependent, and provides only marginal benefit for MACS2 peak calling. (2) Computing in a Python script imported at DAG time. Rejected because it would require the FASTA to exist before DAG resolution, breaking clean builds. (3) Using a checkpoint rule. Rejected as unnecessarily complex; a standard rule with input dependencies achieves the same effect.

### Auto-resolved NCBI TaxId and auto-derived GO database name

**Decision**: The NCBI TaxId (`ncbi_taxid`, renamed from `ncbiID`) is auto-resolved at pipeline runtime from genus+species using `ncbi-datasets-cli`. The GO database name is auto-derived as `org.<G><species>_<GenomeName>.eg.db` (e.g. `org.Athaliana_ColCEN.eg.db`). Both the `ncbi_taxid` and `go_database` config fields are no longer required — `ncbi_taxid` can be set as an optional override, and `go_database` is removed entirely.

**Rationale**: Both values are fully derivable from genus, species, and genome name. Requiring users to look up TaxIds manually and construct GO database names was error-prone and created unnecessary data entry. The genome-name suffix in the GO database name prevents collisions when multiple reference genomes share the same binomial (e.g. two maize assemblies). The `resolve_taxid` rule follows the same pattern as `compute_genome_stats`: produce a JSON file consumed as input by downstream rules, with user-provided values taking precedence.

**Design**: A new `resolve_taxid` rule in `environment_setup.smk` produces `genomes/{ref_genome}/taxid.json`. When `ncbi_taxid` is set in config (and is not `<auto>`), the override is written directly. Otherwise, `datasets summary taxonomy taxon "{genus} {species}"` is called. If the lookup fails, the JSON contains `null` and a warning is emitted — dependent analysis (GO) will fail at the R script stage with a clear error. The `create_GO_database` rule reads the TaxId from the JSON file at execution time.

**Alternatives considered**: (1) Looking up TaxId at DAG time via a params lambda calling `subprocess`. Rejected for the same reason as genome stats — it would require network access during DAG resolution and fail on air-gapped nodes. (2) Keeping `go_database` as an optional override. Rejected because no use case was identified for a user-specified name that deviates from the convention, and removing it simplifies the config surface.

### Auto-derived GTF from GFF via gffread

**Decision**: When `gtf_file` is omitted or set to `<auto>` in the genome config, the `check_gtf` rule automatically derives the GTF from the GFF using `gffread -T`. User-provided GTF paths continue to work as overrides. `gtf_file` is no longer a required field in startup validation.

**Rationale**: The GTF is only consumed by STAR indexing (`make_STAR_indices`), and is fully derivable from the GFF annotation that users must already provide. Requiring both files was a common source of setup friction — users had to manually run `gffread` before starting the pipeline. Since the conversion is deterministic and fast, auto-deriving it follows the same pattern as `genomesize`, `star_index`, and `ncbi_taxid`: compute by default, allow override.

**Design**: The `check_gtf` rule takes the GFF as an input (creating a dependency on `check_gff`), reads the `gtf_file` config value via `.get('gtf_file', '<auto>')`, and branches: user-provided paths are validated and copied/decompressed as before; `<auto>` triggers `gffread -T`. An empty-output guard catches GFF files that gffread cannot convert (e.g. non-standard GFF formats) and directs users to supply a GTF explicitly. `gffread` (bioconda) was added to the `epibutton.yaml` conda environment.

**Alternatives considered**: (1) Using a checkpoint rule. Rejected as unnecessary — a standard rule with an input dependency on the GFF achieves the same ordering. (2) Running gffread unconditionally and ignoring user GTF files. Rejected because some users have curated GTF files with modifications not present in their GFF (e.g. additional transcript annotations).

### Auto-derive structural RNA FASTA via Infernal

**Decision**: When `structural_rna_fafile` is omitted or set to `<auto>`, the pipeline automatically builds a genome-specific structural RNA FASTA by running Infernal `cmscan` against the Rfam covariance model database. The Rfam database is downloaded once and shared across genomes. User-provided FASTA paths continue to work as overrides.

**Rationale**: The previous approach required users to manually query the Rfam MySQL database, extract sequences with easel tools, and configure the path — a 9-step process documented in `Help/Structural_RNAs_Rfam.md`. This was the single highest barrier to enabling sRNA structural RNA depletion. Worse, users often reused a non-species-specific FASTA (e.g. Zea mays structural RNAs for Arabidopsis), reducing depletion specificity. Infernal `cmscan` with Rfam CMs produces genome-specific annotations of all structural RNA loci (rRNA, tRNA, snoRNA, snRNA, etc.) directly from the reference FASTA.

**Design**: Two new rules in `environment_setup.smk`. `download_rfam` fetches `Rfam.cm.gz` and `Rfam.clanin` from EBI FTP and runs `cmpress`. `build_structural_rna_db` splits the genome into per-chromosome FASTAs (contigs < 1 Mbp are binned together to avoid excessive parallelism with fragmented assemblies), runs `cmscan` in parallel via `xargs -P`, merges results, filters clan overlaps (lines marked ` = ` in cmscan `--fmt 2` output), and extracts sequences with `bedtools getfasta`. The `-Z` flag is set to the total genome size (both strands, in Mb) for proper E-value calibration across per-chromosome runs. Default threshold is `--cut_ga` (Rfam per-family gathering thresholds); configurable via `infernal_threshold` option (`"ga"` or a numeric E-value). `get_structural_rna_fasta()` in `smallRNA.smk` routes to the auto-derived or user-provided FASTA.

**Alternatives considered**: (1) Snakemake checkpoint scatter/gather for per-chromosome parallelization. Rejected as over-complex — checkpoint rules require DAG reconstruction and the number of chromosomes isn't known until runtime; internal `xargs -P` parallelization within a single rule is simpler and equally effective. (2) Per-sequence E-value threshold instead of `--cut_ga`. Rejected as default because gathering thresholds are per-family calibrated by Rfam curators and avoid both over- and under-filtering; E-value is available as an advanced override for divergent organisms.

---

## CLI Wrapper (`epicc`)

### Python script with subcommands, not a bash wrapper

**Decision**: The `epicc` CLI wrapper is a Python script using `argparse` subcommands (`run`, `dry-run`, `validate`, `profile`, `clean`). It constructs and invokes snakemake as a subprocess, mapping user-friendly flags to snakemake's `--config`, `--configfile`, `--profile`, and other arguments. Arbitrary snakemake arguments can be passed through via `--` separator.

**Rationale**: Python is consistent with the rest of the codebase (Snakefile, sample_sheet.py, profiler). `argparse` provides structured subcommands with auto-generated help text and type validation — all cumbersome in bash. The wrapper needs to read YAML config files (PyYAML already in the conda env) and may integrate more tightly with pipeline internals over time. Using `subprocess.run()` (not `os.exec`) allows the wrapper to capture exit codes and print post-run summaries.

**Alternatives considered**: (1) Bash wrapper script. Rejected due to poor YAML support, verbose argument parsing, and divergence from the Python codebase. (2) Click instead of argparse. Rejected as an unnecessary dependency — argparse is in the standard library and sufficient for this use case. (3) `--smk` flag for snakemake passthrough. Rejected in favor of the standard POSIX `--` separator, which is more intuitive and doesn't require a custom flag.

### Configurable output directories via dual-context path replacement

**Decision**: Hardcoded `results/` and `genomes/` path prefixes were replaced across all 9 pipeline files (827 lines) using two different patterns depending on context. In Python context (input/output directives, expand calls, params lambdas), a `RESULTS_DIR`/`GENOMES_DIR` variable derived from `config["output_dir"]`/`config["genome_dir"]` is used with f-string wrapping. In shell blocks (triple-quoted `shell: """..."""`), Snakemake's native `{config[output_dir]}`/`{config[genome_dir]}` substitution is used. Defaults are `"results"` and `"genomes"` for full backward compatibility.

**Rationale**: Shell blocks in Snakemake are not Python code — they are template strings where `{...}` patterns are substituted by Snakemake at runtime. Converting shell blocks to Python f-strings (prefixing with `f` and doubling all existing `{wildcards.x}`, `{params.x}`, `{input}`, `{output}`, etc.) would require hundreds of additional changes and make shell blocks harder to read. The `{config[key]}` syntax is natively supported in Snakemake shell blocks alongside other substitutions, providing a clean solution with no extra infrastructure (no `params:` injection needed).

**Alternatives considered**: (1) Symlinks (`results/` → target directory). Rejected: fragile, breaks on move, doesn't support concurrent runs with different output directories. (2) `workdir:` directive. Rejected: changes CWD for all rules, breaks relative paths in config. (3) Passing directory via `params:` in each rule. Rejected: would require adding `results_dir=RESULTS_DIR` to params in every rule with a shell block — more invasive and error-prone than `{config[key]}`. (4) Converting shell blocks to f-strings. Rejected: requires doubling all existing `{...}` patterns in shell blocks, hurting readability for minimal benefit.

---

## Data Acquisition

### URL support for Read_files and genome config fields

**Decision**: Read_files in the sample sheet and genome config fields (`fasta_file`, `gff_file`, `te_file`) accept HTTP(S) URLs in addition to local paths and SRA accessions. URLs are downloaded via curl at rule execution time with retry logic, redirect following, and scheme restriction (`--proto '=https,http'`).

**Rationale**: Making integration test cases self-contained and reproducible requires that all inputs can be fetched from public sources without pre-staging data on local filesystems. The ColCEN test case uses GitHub-hosted genome files and lemna.org-hosted modBAMs, which cannot be addressed as SRA accessions. URL support also benefits users who store data on institutional web servers or cloud storage rather than SRA.

**Design**: URL detection occurs at two levels:
1. **Sample sheet parsing** (`sample_sheet.py`): `get_seq_id_and_path()` detects `http://`/`https://` prefixes. BAM/bedMethyl URLs behave like local paths (URL string IS the path). FASTQ URLs use `seq_id="URL"` as a dispatch sentinel, analogous to `"EXPLICIT"` for local paths.
2. **Rule shell blocks**: Reference file rules (`check_fasta`, `check_gff`, `check_te_file`) and sample input rules (`get_fastq_pe/se`, `get_dmc_input`) check for URL prefixes before local file checks. Downloads use `curl --fail --location --max-redirs 5 --retry 3 --proto '=https,http'` for security and reliability.

**Alternatives considered**: (1) A pre-download step that resolves all URLs before DAG execution. Rejected because it would require network access during DAG resolution and would not integrate with Snakemake's retry mechanism. (2) Using `urllib`/`requests` in Python. Rejected because shell-level curl is consistent with the existing download patterns (ENA downloads) and available in all conda environments.

### GFF3-to-BED conversion for TE annotations

**Decision**: The `check_te_file` rule accepts `.gff3(.gz)` and `.gff(.gz)` files in addition to `.bed(.gz)`. GFF3 input is automatically converted to BED6 format using awk, extracting the `ID=` attribute from column 9 as the name (column 4). The existing BED uniqueness validation runs on the output regardless of input format.

**Rationale**: The Col-CEN reference genome distributes TE annotations in EDTA GFF3 format (`t2t-col.20201227.fasta.mod.EDTA.TEanno.gff3.gz`). Requiring users to pre-convert GFF3 to BED is an unnecessary manual step. Since EDTA GFF3 uses `ID=TE_homo_NNN` (already unique), the conversion is deterministic and lossless for the fields the pipeline uses (coordinates, name, strand).

**Alternatives considered**: (1) Requiring BED format only and documenting the conversion command. Rejected because it creates a manual step that blocks self-contained test cases. (2) Using `bedtools` or `gffread` for conversion. Rejected because awk is simpler, has no additional dependencies, and handles the specific GFF3→BED6 mapping needed here.

### Execution mode auto-detection

**Decision**: When neither `--profile` nor `--cores` is specified, the wrapper auto-detects SLURM by checking for `sbatch` in PATH. If found, it uses `profiles/slurm`; otherwise it runs locally with half of `nproc` cores.

**Rationale**: Users shouldn't need to remember different invocations for different environments. The explicit `--profile` and `--cores` flags serve as overrides for non-standard setups. This pattern matches the existing `validate_pombe.sh` script.

### No legacy assay tokens

**Decision**: Ahead of 1.0.0, the combined assay tokens (`ChIP_broad`, `CUT_RUN_narrow`, …) are no longer accepted in a sheet's `Assay` column, the builder no longer splits them on import, and `scripts/migrate_sample_sheet.py` is removed. `Assay` takes the method only; `Peak_type` carries broad/narrow.

The combined tokens remain as an **internal** representation: `add_compat_columns` folds `(Assay, Peak_type)` into one via `combine_assay_peaktype` and writes it back over the `Assay` column, and `ASSAY_TO_ENV`, `ASSAY_TO_PEAKTYPE` and the analysis-name builder all key on that form. They are listed in `_COMBINED_PEAK_ASSAYS`, not in `VALID_ASSAYS`.

**Rationale**: Dual-form input meant two ways to express one thing, and the validation code carried a branch whose only job was to tell users which form to switch to. A 1.0.0 release is the point to require the canonical form. Keeping the combined token internally is unrelated to back-compat — it is the key that output paths and env lookups are built from, and removing it would rename every ChIP output.

**Alternatives considered**: Removing the combined token throughout and keying lookups on `(Assay, Peak_type)` pairs. Rejected for this release: it renames every post-alignment ChIP path on top of the rename already coming from the Genome split, for no user-visible gain.

### One genome per run, and per-run TE files (`--ref-genome`, `--te-file`)

**Decision**: `epicc run`/`validate --ref-genome GENOME` sets a new `only_genome` option that drops every exploded sample-sheet row for any other genome (`restrict_to_genome`, applied before `check_table` and `check_genome_config`). `epicc output` sets it from its own `--ref-genome`. `--te-file PATH` writes `genomes.<GENOME>.alt_te_file` into the merged options file; the genome is `--ref-genome`, or the sheet's only genome, and anything else is an error. An `alt_te_file` is prepared as `<GENOME>__TE_file_<tag>.bed` and analysed under the target name `all_TEs_<tag>`, the tag taken from its file name (`workflow/scripts/te_files.py`); the options-file `te_file` keeps `<GENOME>__TE_file.bed` and `all_TEs`. The alternative rule records its source in `<GENOME>__TE_file_<tag>.source`, and a later file whose name gives the same tag but whose path differs gets the tag plus 8 hex digits of the path's SHA-1. An `alt_te_file` equal to `te_file` gets no tag. `check_te_file`'s shell moved into a shared `TE_PREP_SHELL` string, byte-identical to before, so existing TE files are not re-prepared by the code rerun-trigger. `--te-file` does not turn on `te_analysis`, but in `epicc output --plot-type browser` it does turn on `browser_TE_file`.

**Rationale**: A TE annotation is species-specific, so a CLI TE file can only mean one genome, and it is only unambiguous when the run is about one genome. Filtering before validation means a sheet listing `B73,W22` can be run for B73 without a W22 entry in the options file, which is the point of asking for one genome (issue #98). The genome key is written as a tuple path, not a dotted one, because genome names may contain dots. The CLI TE file is tagged rather than written over `te_file` because each genome has one prepared TE file and every TE output (about 20 heatmaps and profiles, the sRNA counts and checkpoint) names the target `all_TEs` with no trace of which annotation it used: swapping the source makes Snakemake re-prepare the file, so a run with `--te-file` followed by one without it would redo all TE analyses twice, which takes hours on large genomes, and leave whichever results ran last under the shared names.

**Alternatives considered**: (1) A `--te-file GENOME=PATH` form taking several genomes. Rejected as more syntax for a case the options file already covers. (2) Filtering targets, not samples, by genome. Rejected: every genome's config would still be validated and parsed, so an unconfigured genome would still stop the run. (3) Writing the CLI TE file over `te_file`, or tagging only the prepared file with a fixed suffix such as `_cli`. Rejected: a fixed suffix stops the options-file TE file being rebuilt but two different CLI files still share one name, and the downstream `all_TEs` outputs would still rerun whenever their input changed. (4) Always tagging with a hash of the path. Rejected for unreadable names (`all_TEs_3fa1c2d9`); the hash is only the fallback for a name clash. (5) Offering `--te-file` on `epicc output --plot-type browser` only. Rejected: issue #98 asked for it on `run`.

### Extra browser tracks are a list option, set from the CLI through the merged options file

**Decision**: `browser_extra_tracks` is a list of BED paths or `{file, label, color}` mappings, normalized by `workflow/scripts/browser_extra_tracks.py` (labels from file names, colours from a fixed palette). `epicc output --extra-tracks a.bed,b.bed [--extra-tracks-labels …] [--extra-tracks-colors …]` replaces the list for that run through the temp merged options file, since `--config` cannot express a list of mappings. `prep_browser_on_region` clips each file to the locus and writes an `extra_tracks__*.txt` table that `R_browser_plot.R` reads as its 8th argument. The plot script now always gets the highlight files, empty when there are none.

**Rationale**: The browser rules are also reached from a full run, so the setting belongs in the options file, with the CLI flags as a per-run override; one normalizer keeps the CLI's early errors and the workflow's startup check identical. Features are clipped, not just filtered, because Gviz plots the union of every track's range, so one feature reaching past the locus zoomed the plot out (seen with a multi-megabase sRNA cluster on hg38 chr21). Labels go through a shell loop, so tabs, quotes, `$`, backticks and backslashes are rejected instead of escaped.

**Alternatives considered**: (1) CLI flags only, with no options key. Rejected: full runs could not use the tracks. (2) Passing tracks to the R script as more positional arguments. Rejected: a manifest file keeps the argument list fixed whatever the number of tracks.

### Commented samples in the builder

**Decision**: `tools/epicc-builder.html` classifies a `#` line the way `read_sample_sheet` skips it (first non-whitespace character `#`). After the header, a `#` line with at least as many tab-separated fields as the required columns, and not a commented-out header, becomes a commented sample row: the prefix (`#`, `# `, `  ## `, …) is stored and the row is greyed out but editable. Every other `#` line, and every blank line, is a note kept verbatim: before the header in `sheetPreamble`, after the last row in `sheetTrailer`, otherwise on the row that follows it. Commented rows skip all validation and are left out of the cross-row checks, with two exceptions: an active row whose `Control` names a commented sample is an error, and a commented row sharing an active row's `Sample_ID` is a warning. A commented row's `Levels` is written back verbatim unless its factor cells were edited, and the sheet's factors come from active rows first (issue #65).

**Rationale**: The pipeline never sees a commented row, so errors on it would block export for nothing. The Control case is an error, not a warning: `samplefile_validation.py` rejects a Control that matches no Sample_ID, and a commented control is exactly that. Keeping the prefix and the notes is what makes import→export byte-identical. A parked row may come from an earlier design with different factors, so rebuilding its Levels from the current factor columns would silently drop or blank some of them.

**Alternatives considered**: (1) Treating every `#` line after the header as a sample. Rejected: section headings and other prose would become broken rows. (2) Dropping non-sample comment lines on import. Rejected: the round trip would lose them. (3) A warning for an active row pointing at a commented control. Rejected because the pipeline refuses that sheet. (4) Writing notes in "Copy as command". Not done: its printf format needs a fixed number of columns, so only Export TSV writes notes; commented rows are written by both.
