# Sample Sheet Specification

This document defines the canonical format and validation rules for EPICC sample
sheets. It is the single source of truth for both the pipeline
(`workflow/scripts/samplefile_validation.py`) and the epicc-builder app.

## Format

Tab-separated values (TSV), UTF-8 encoding. First row is the column header.
Trailing whitespace on any field is stripped during parsing.

## Columns

| # | Column | Required | Type |
|---|--------|----------|------|
| 1 | Sample_ID | Yes | Controlled freetext |
| 2 | Assay | Yes | Controlled vocabulary |
| 3 | Genome | Yes | Freetext |
| 4 | Levels | Yes | Structured freetext |
| 5 | Replicate_ID | Yes | Freetext |
| 6 | Read_files | Yes | Path / SRA ID |
| 7 | Read_layout | Yes | Controlled vocabulary |
| 8 | IP_target | Conditional | Freetext |
| 9 | Control | Conditional | Reference |
| 10 | Peak_type | Conditional | Controlled vocabulary |
| 11 | Comments | Optional | Freetext |

`Peak_type` is last because it applies to only three assays; column lookup is by
header name on both sides, so the order is ergonomic rather than semantic.

## Per-Field Validation Rules

### Sample_ID

Unique identifier for each sample row. Used as a filesystem path component and
as a reference target for the Control column.

| Rule | Detail |
|------|--------|
| Required | Must be non-empty |
| Unique | No two rows may share the same Sample_ID |
| No double underscores | Must not contain `__` (reserved as a delimiter in analysis-level names) |
| Filesystem-safe | Must not contain any of: `/\` whitespace `'";&\|<>$`` !{}()[]?*~#` |

There are no other format requirements. Sample_IDs can be fully arbitrary
strings (the pipeline never parses them for metadata). However, concise
descriptive IDs like `WT_leaf_H3K9me2_rep1` are recommended for readability.
The epicc-builder app should suggest IDs progressively as the user fills in
other fields, but always allow manual editing.

### Assay

Controlled vocabulary defining the experimental **method** only. The
peak-calling type (broad/narrow) is a separate parameter — see **Peak_type**
below. Determines pipeline routing and environment mapping.

**Valid values:**

| Assay | Env folder | Description |
|-------|------------|-------------|
| `ChIP` | `ChIP` | ChIP-seq (histone marks or TFs; set Peak_type) |
| `CUT_RUN` | `ChIP` | CUT&RUN (set Peak_type; default caller broad→epic2, narrow→SEACR) |
| `CUT_TAG` | `ChIP` | CUT&Tag (set Peak_type; default caller broad→epic2, narrow→SEACR) |
| `ATAC` | `ATAC` | ATAC-seq (narrow peaks, fixed; no Peak_type / IP_target) |
| `RNAseq` | `RNA` | RNA-seq (mRNA/total RNA) |
| `RAMPAGE` | `RNA` | RAMPAGE TSS profiling |
| `sRNA` | `sRNA` | Small RNA-seq |
| `WGBS` | `mC` | Whole-genome bisulfite sequencing |
| `WGBS_nd` | `mC` | Non-directional WGBS (e.g. Zymo Pico) |
| `PBAT` | `mC` | Post-bisulfite adapter tagging |
| `EMseq` | `mC` | Enzymatic methyl-seq |
| `dmC` | `mC` | Direct methylation (nanopore modBAM or pre-computed bedMethyl) |

These values are defined in `workflow/scripts/sample_sheet.py:VALID_ASSAYS`.

**Pulldown family.** `ChIP`/`CUT_RUN`/`CUT_TAG` are collected in
`IP_PEAK_ASSAYS` and share the same semantics: `IP_target` is required (for both
IPs and their controls), `Peak_type` is required on rows that call peaks (those
with a `Control`) and optional on control-only rows, and `Control` may reference
another sample. Default peak caller by peak shape: ChIP → MACS2, CUT&* broad →
epic2, CUT&* narrow → SEACR. Override via `cut_callpeaks.broad_caller` /
`narrow_caller` in the options file (`epic2`, `seacr`, or `macs2`).

### Peak_type

The peak-calling type — a separate analytical parameter from Assay. Controlled
vocabulary: `broad` or `narrow`.

| Rule | Detail |
|------|--------|
| Required | Non-empty (`broad`/`narrow`) for `ChIP`/`CUT_RUN`/`CUT_TAG` |
| Blank otherwise | Must be empty for all other assays (ATAC's type is fixed; the rest have no peaks) |

`broad` suits histone-domain marks (H3K9me2, H3K27me3, H3K36me3, …); `narrow`
suits transcription factors and sharp marks (H3K4me3, H3K27ac, …).

### Genome

Reference genome name. Must be non-empty. This value is used to locate genome
files under `genomes/{Genome}/` and to group samples for comparative analyses.
Examples: `ColCEN`, `Spombe`, `hg38`.

**Multiple references per sample.** A comma-separated list maps the same reads
to several genomes — e.g. `B73,W22`. The row is exploded internally into one
`(Sample_ID, Genome)` entry per reference, so `Sample_ID` remains globally
unique and does not encode the genome.

| Rule | Detail |
|------|--------|
| Required | Must be non-empty |
| Each entry must exist | Every listed genome needs an entry under `genomes:` in the options file |
| No `__` | A genome name must not contain `__`; it delimits the genome in post-alignment filenames |
| No duplicates | `A,A` is rejected (`explode_genomes`) |
| Empties tolerated | `A,,B` parses as two references; the builder flags the stray comma as a warning, not an error |

Read download, trimming and raw/trim QC happen **once per sample** and are
genome-independent; alignment and everything after it run **once per genome**.
See "Per-replicate naming" below.

In the builder, the Genome cell is a tick-box picker (`genomeMultiEditor`)
listing every genome known from the Reference Genomes tab and from the other
rows, plus an inline field for adding a new one. Ticks are written back as the
comma-separated list; the cell also stays typeable so paste still works.

### Levels

Comma-separated list of `factor:level` pairs describing experimental conditions.

**Format**: `factor1:level1,factor2:level2,...`

**Examples**:
- `genotype:WT,tissue:root`
- `genotype:dcr1,tissue:cell`
- `temperature:37deg`

| Rule | Detail |
|------|--------|
| Required | Must be non-empty |
| Pair format | Each entry must contain exactly one `:` separating a non-empty factor name from a non-empty level value |
| Consistent factor count | All rows must have the same number of `factor:level` pairs |
| Consistent factor names | All rows must use the same factor names in the same order (e.g. if row 1 has `genotype,tissue`, all rows must have `genotype,tissue`) |

**Derived values**: The level values (not factor names) are joined with `_` to
form the `levels_label` used in analysis-level filenames. For example,
`genotype:WT,tissue:root` produces `WT_root`.

**Controls**: Control samples (Input, WCE, IgG, etc.) must also specify Levels
with meaningful values where possible (e.g. `genotype:WT,tissue:root`), since
Levels determine how samples are grouped for analysis.

### Replicate_ID

Identifier for biological or technical replicates. Must be non-empty.

Common values: `rep1`, `rep2`, `repA`, `repB`, `1`, `2`.

Replicates are processed independently through mapping and initial analysis,
then merged for downstream steps like peak calling and IDR.

**Note**: `+` in Read_files concatenates multi-file inputs from the same library
(e.g. multiple sequencing lanes). This is distinct from replicate handling.

### Read_files

Path(s) to input data files or SRA accession(s). Must be non-empty.

**Supported formats**:

| Format | Syntax | Example |
|--------|--------|---------|
| SRA accession | `SRRnnnnnnn` | `SRR12345678` |
| SRA merge | `SRR...+SRR...` | `SRR111+SRR222+SRR333` |
| Local FASTQ (SE) | `/path/to/file.fq.gz` | `/data/reads.fastq.gz` |
| Local FASTQ (PE) | `/path/r1.fq.gz,/path/r2.fq.gz` | `/data/S1_R1.fq.gz,/data/S1_R2.fq.gz` |
| Local BAM | `/path/to/file.bam` | `/data/aligned.bam` |
| bedMethyl (dmC) | `/path/to/file.bed.gz` | `/data/methylation.bed.gz` |
| FASTQ merge (SE) | `f1.fq.gz+f2.fq.gz` | `/data/L1.fq.gz+/data/L2.fq.gz` |
| FASTQ merge (PE) | `a_R1,a_R2+b_R1,b_R2` | `/d/A_R1.fq.gz,/d/A_R2.fq.gz+/d/B_R1.fq.gz,/d/B_R2.fq.gz` |

`+`-merge concatenates components after download. Supported for **SRA
accessions** and **FASTQ files** (SE or PE; local or URL) — the components
are fetched/copied and concatenated (gzip members concatenate cleanly).
**BAM** (would need `samtools merge`) and **bedMethyl** (per-position
counts must be summed, not concatenated) cannot be `+`-merged; merge those
upstream. All components of a merge must be the same type.

| Rule | Detail |
|------|--------|
| Required | Must be non-empty |
| SRA regex | SRA accessions must match `^[SDE]RR\d+$` |
| No mixing | Within a single `+`-separated component, all entries must be the same type (all SRA or all paths) |
| Mergeable types | A multi-component `+`-merge must be all the same type and either all SRA accessions or all FASTQ files; BAM/bedMethyl/mixed `+`-merges are rejected (merge those upstream) |
| PE comma pair | For PE local FASTQs, each component should have a comma-separated R1,R2 pair |

| Rule | Detail |
|------|--------|
| PE comma pair required | PE layout with a single non-BAM path per component is an error |
| SE single path required | SE layout with multiple comma-separated paths is an error |
| No cross-row duplicates | The same file path or SRA accession must not appear in more than one sample's Read_files |
| Local paths must exist | Each local-path entry must resolve to an existing file on disk before the run starts. SRA accessions, HTTP(S) URLs and `s3://` URIs are not probed. |
| S3 URIs well-formed | An `s3://` entry must name both a bucket and a key (`s3://bucket/key`). Public (authentication-free) objects only — credentials and request signing are not supported. |

### Read_layout

Sequencing layout. Must be exactly `SE` (single-end) or `PE` (paired-end).

### IP_target

Describes what was immunoprecipitated (or equivalent). The value appears in
analysis-level names and plot labels.

| Rule | Detail |
|------|--------|
| Required for IP-peak assays | Must be non-empty for any assay in `IP_PEAK_ASSAYS` (`ChIP_broad/narrow`, `CUT_RUN_broad/narrow`, `CUT_TAG_broad/narrow`), **including control samples** (e.g. `Input`, `WCE`, `IgG`) |
| Blank for others | Must be empty/absent for all non-IP assays (`ATAC`, `RNAseq`, `RAMPAGE`, `sRNA`, `WGBS`, `EMseq`, `dmC`) |

**Examples**: `H3K9me2`, `H3K4me3`, `H3K27me3`, `CTCF`, `Input`, `WCE`, `IgG`

### Control

References the Sample_ID of the control sample used for normalization.

| Rule | Detail |
|------|--------|
| Valid reference | Must match an existing Sample_ID in the sheet |
| Allowed assays only | May only be specified for IP-peak assays (`ChIP_*`, `CUT_RUN_*`, `CUT_TAG_*`) and `RAMPAGE` samples |
| Not itself | A sample may not list itself as its own Control |
| Chain depth ≤ 2 | The referenced control **may** declare its own Control (the *dual-role* case below), but that one must not. Deeper chains are rejected; the bound also makes control cycles impossible |
| Sharing allowed | Multiple IP samples may reference the same control (typical CUT&RUN convention: one IgG per batch) |
| Optional | Blank/absent is valid (sample has no associated control) |

**Examples**:
- ChIP: `WT_leaf_Input_rep1` (`Assay: ChIP_broad`, `IP_target: Input`)
  referenced as `Control` by H3K9me2 IP samples.
- CUT&RUN: `WT_endosperm_IgG_rep1` (`Assay: CUT_RUN_broad`, `IP_target: IgG`)
  shared as `Control` by both H3K27me3 replicates in the same family. The
  builder constrains the Control dropdown to same-family assays
  (`ChIP_*`, `CUT_RUN_*`, or `CUT_TAG_*`), but freetext entry is allowed
  for cross-family controls if a study uses one.

**Dual-role samples.** A sample can serve as another row's control *and* be
analysed in its own right, as long as it declares a `Control` of its own. The
motivating case is an H3 ChIP used both to examine raw H3 distribution (its own
peaks, tracks, merged-replicate analysis, IDR and plots — normalized against
Input) and as the peak-calling control for H3K9me2:

```
Input_rep1     IP_target: Input      Control: (blank)
H3_rep1        IP_target: H3         Control: Input_rep1   <- dual role
H3K9me2_rep1   IP_target: H3K9me2    Control: H3_rep1
```

Dual-role samples receive the full analysis treatment. Qualification is
structural — *does this row declare a Control?* — not *is anyone using this row
as their control?*; see `is_peak_call_target` / `peak_callable_rows` in
`workflow/scripts/sample_sheet.py`. A pulldown row with **no** `Control` is not
analysable (peak calling requires one) and is dropped from peak/FC and
analysis-level targets, which is also what excludes an unreferenced
Input/WCE/IgG row.

### Comments

Free-text annotation for the user's own reference — batch, library prep,
sequencing run, caveats, anything. **The pipeline never reads it.**

| Rule | Detail |
|------|--------|
| Optional | May be blank, and the column may be absent entirely — sheets written before it existed still load (filled with `""`) |
| Unconstrained | Any text; no controlled vocabulary, no length limit, not assay-gated (unlike `Peak_type`/`IP_target`/`Control`, it is editable on every row) |
| Inert | Never contributes to analysis keys, analysis names, merge keys, or output paths |
| Tab/newline safe | The builder strips `\t`, `\r` and `\n` on export (`sanitizeCell`) so a pasted multi-line note cannot break the TSV |

A `#` inside the field is preserved: `read_sample_sheet` drops only *full-line*
comments, never truncating at a mid-field `#` (the same guarantee that keeps
URLs with fragments intact in `Read_files`).

## Cross-Field Validation

| Check | Severity | Detail |
|-------|----------|--------|
| Read_layout vs Read_files (PE with single path) | Error | PE layout but Read_files has only one non-BAM path per component |
| Read_layout vs Read_files (SE with multiple paths) | Error | SE layout but Read_files has multiple comma-separated paths |
| Duplicate Read_files entries | Error | Same file path or SRA accession used by more than one Sample_ID |

## Derived Names

These are not user-specified columns but are computed from the sample sheet:

- **levels_label**: `_`-joined level values from Levels (e.g. `WT_root`)
- **analysis_name**: Components `(Assay, levels_label, IP_target, Genome)` joined
  with `__`, omitting empty parts. For ChIP samples (with IP_target):
  `ChIP_broad__WT_root__H3K9me2__ColCEN`. For non-ChIP (blank IP_target):
  `RNAseq__WT_root__ColCEN`.
- **env**: Pipeline environment folder, mapped from Assay via `ASSAY_TO_ENV`
- **peak_type**: `broad` or `narrow`, mapped from Assay via `ASSAY_TO_PEAKTYPE`

## Implementation References

- Validation: `workflow/scripts/samplefile_validation.py`
- Constants and helpers: `workflow/scripts/sample_sheet.py`
