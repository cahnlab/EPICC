# Epigenetic Pipeline for Integrative Chromatin Characterization (EPICC) 
# aka Epigenetic Button

A Snakemake-based pipeline for analyzing and integrating various types of (epi)genomics datasets, including histone and transcription factor ChIP-seq, RNA-seq, RAMPAGE, small RNA-seq, and methylC-seq.

## Complete documentation: [Read the docs](https://epicc-documentation.readthedocs.io/en/latest/)

## Overview

EpigeneticButton is a comprehensive pipeline that processes and analyzes multiple types of genomics data. It provides an automated workflow for:
- Data preprocessing and quality control
- Read mapping and alignment
- Peak calling and differential expression analysis
- Data integration and visualization

## Features

- **Multiple Data Types Support**:
  - ChIP-seq — histone marks (broad peaks) and transcription factors (narrow peaks)
  - CUT&RUN / CUT&Tag — broad and narrow peak calling
  - ATAC-seq
  - RNA-seq
  - RAMPAGE
  - small RNA-seq
  - MethylC-seq — Bisulfite sequencing via Bismark (WGBS, WGBS_nd, PBAT, EMseq)
  - Direct Methylation (dmC) — Long-read native methylation (ONT, PacBio)

- **Automated Analysis**:
  - Reference genome preparation
  - Sample-specific processing
  - Data type-specific analysis
  - Combined analysis across samples
  - Quality control and reporting
  - Additional output options such as heatmaps, metaplots and browsers

- **Flexible Configuration**:
  - Self-contained local HTML builder (`tools/epicc-builder.html`) to build and validate sample sheets and options files
  - Customizable mapping parameters
  - Configurable analysis options
  - Resource management
  - Parallel processing

## Installation

1. Clone the repository:
```bash
git clone https://github.com/cahnlab/EPICC.git
```
or for ssh connection
```bash
git clone git@github.com:cahnlab/EPICC.git
```
```bash
cd EPICC
```

2. Install snakemake and other required packages from the depency file:
```bash
conda create -n epicc -y --file config/epicc-env.txt
conda activate epicc
```
If you want to run the pipeline on a different platform than locally or slurm, you will need to also install the corresponding snakemake-executor-plugin

## Usage

### Configuration

For new users, it is recommended to use the local HTML builder at `tools/epicc-builder.html` (open directly in a browser — self-contained, no internet connection required) to build and validate your sample metadata file and choose analysis options.

1. Prepare your sample metadata file (start from the documented template at `config/example_samples.tsv` and pass yours via `epicc run --samples ...`) with the required columns below (see Input requirements for more details specific to each data-type):
   - `Sample_ID`: Unique identifier for the sample (e.g. `WT_leaf_H3K9me2_rep1`). Must be filesystem-safe.
   - `Assay`: Type of assay — the experimental method only [`ChIP` | `CUT_RUN` | `CUT_TAG` | `ATAC` | `RNAseq` | `RAMPAGE` | `sRNA` | `WGBS` | `WGBS_nd` | `PBAT` | `EMseq` | `dmC`].
   - `Genome`: Reference genome name (e.g. `ColCEN`, `Spombe`), or a comma-separated list (e.g. `B73,W22`) to map the same reads to several references
   - `Levels`: Experimental conditions as comma-separated `factor:level` pairs (e.g. `genotype:WT,tissue:leaf`)
   - `Replicate_ID`: Replicate identifier (e.g. `rep1`, `rep2`)
   - `Read_files`: Path to FASTQ/BAM files, SRA accession (e.g. `SRR12345678`), HTTP(S) URL, or public `s3://bucket/key` URI. For PE FASTQs, comma-separate R1 and R2 paths. Use `+` to merge multiple SRA accessions or FASTQ files (e.g. `SRR111+SRR222`, or `a.fq.gz+b.fq.gz`); BAM/bedMethyl inputs must be merged upstream.
   - `Read_layout`: [`PE` | `SE`]
   - `IP_target`: Required for pulldown assays (`ChIP`/`CUT_RUN`/`CUT_TAG`) — the IP target (e.g. `H3K9me2`) or control type (e.g. `Input`, `WCE`, `IgG`). Leave blank for other assays.
   - `Control`: Sample_ID of the control sample for this IP (e.g. the Input sample's Sample_ID). Required for pulldown and RAMPAGE IP samples. Leave blank for plain controls and non-pulldown assays.
   - `Peak_type`: [`broad` | `narrow`] — the peak-calling type, kept separate from the assay. Required on `ChIP`/`CUT_RUN`/`CUT_TAG` rows that call peaks (those with a `Control`); leave blank on a control-only row and on every other assay.
   - `Comments`: Optional free text for your own notes. Never read by the pipeline. The column may be omitted entirely.

2. Update `config/epicc-options.yaml` with your paths and parameters:
   - Sample file: pass via `epicc run --samples` on the command line, or set `sample_file:` in the options file
   - Reference genome files: define each genome under the `genomes:` namespace (see [below](#common-to-all-types-of-samples) for format details)
   - Analysis parameters / options
   - Species-specific parameters
   - Resources allocation
   
3. If you are running the pipeline on a cluster, set up a cluster execution profile. EPICC ships templates for SLURM (`profiles/slurm/config.yaml`) and UGE/SGE (`profiles/uge/config.yaml`). Two setup paths:

   **Option A — `epicc init-profile` (recommended; works for both git-clone and conda-install):**
   ```bash
   epicc init-profile slurm my-cluster    # SLURM
   epicc init-profile uge my-cluster      # UGE/SGE
   ```
   This copies the template to `~/.config/snakemake/my-cluster/config.yaml`. Open that file and edit every line marked `EDIT:` for your site (SLURM account, TMPDIR routing, per-rule resource overrides). Then run the pipeline with `--profile my-cluster`.

   **Option B — edit the bundled template in place (git-clone only):**
   Edit `profiles/slurm/config.yaml` directly and pass `--profile profiles/slurm` to `epicc run`. Convenient for single-user setups; not portable across conda-installed instances.

   The SLURM profile defaults to 16 parallel jobs (`jobs: 16`). Per-rule resource overrides go in `set-resources:` within the profile. Memory units are in MB.

4. Default options: 
   - Full analysis: By default, a full analysis is performed form raw data to analysis plots. Change `full_analysis` in the options file ([see below](#main-output-options)).
   - Limited QC output: By default, some QC options are not performed to limit the time and amount of output files. Change `QC_option` in the options file ([see below](#main-output-options)).
   - No Gene Ontology analysis: Due to the difficulty in automating building a GO database, this option is OFF by default. Change `GO` option in the options file. Please refer to Additional output options #2 below and [Help GO](Help/Gene_Ontology.md) before setting it to `true` as it requires 2 other files. These files are available for Arabidopsis thaliana (Tair10 / ColCEN assembly) and Maize B73 (v5 or NAM assembly) in the `data` folder.
   - No TE analysis: By default, no analysis on transposable elements is performed. Change `te_analysis` in the options file ([see below](#main-output-options)).
   - For ChIP-seq and ATAC-seq: the default aligner is bowtie2. Set `chip_aligner`/`atac_aligner` to `"chromap"` for ~10x faster mapping — recommended for SE-only datasets; PE data shows only 11–28% properly-paired reads with chromap vs 99%+ with bowtie2 (chromap maps ends independently, causing one mate to be MAPQ-filtered). Mapping strategy options are available via `chip_mapping_strategy` ([see below](#chip-mapping-parameters)). When using `repeat` or `repeatall` strategies, bowtie2 is used automatically regardless of the aligner setting.
   - For sRNA-seq: the default is not based on Netflex v3 library preparation. If your data was made with this kit, an additional deduplication and read trimming is required. To turn it ON, change the `Netflex_v3_deduplication` in the options file. See [Known issues #3](#known-potential-issues) below if you have mixed libraries.
   - For sRNA-seq: the default is not to filter structural RNAs prior to shortstack analysis. Change `structural_rna_depletion` in the options file.  While this step is recommended for small interfering RNA analysis, it requires a pre-build database of fasta files. Please refer to the [Help structural RNAs](Help/Structural_RNAs_Rfam.md) before setting it to `true`. This file is available for Maize in the `data` folder.
   - For sRNA-seq: the default is to only perform *de novo* micro RNA identification (`--dn_mirna` argument in ShortStack). If you also want the known microRNAs, download the fasta file from [miRbase](https://www.mirbase.org), filter it for your species of interest, and add to the `srna_mapping_params` entry in the options file `--known_miRNAs <path/to/known_miRNA_file.fa>`.

### Running the Pipeline

The main entry point is the `epicc` CLI wrapper at the repository root. It handles configuration validation, TMPDIR routing, and cluster profile detection automatically.

1. To run the pipeline locally:
```bash
epicc run --samples config/my_samples.tsv --cores 12
```

2. To run the pipeline on a SLURM cluster:
```bash
epicc run --samples config/my_samples.tsv --profile profiles/slurm
```

To reduce terminal output, redirect to a log and run in the background:
```bash
epicc run --samples config/my_samples.tsv --profile profiles/slurm > epicc.log 2>&1 &
```

3. To run the pipeline on a UGE cluster (using qsub):
```bash
epicc run --samples config/my_samples.tsv --profile profiles/uge
```

*If using a profile, review and adapt the profile config to your cluster before running.*

4. Additional `epicc run` flags:
   - `--options FILE` — options YAML (default: `config/epicc-options.yaml`)
   - `--output-dir DIR` — results directory prefix (default: `results`)
   - `--genome-dir DIR` — genome directory prefix (default: `genomes`)
   - `--keep-intermediates TIER` — intermediate file retention: `none`, `standard` (default), `custom`, `all`
   - `--use-node-tmpdir` — skip the workflow's TMPDIR override and use the cluster's default (see TMPDIR routing below)
   - `--no-rerun-incomplete` — skip re-running Snakemake-flagged incomplete jobs from prior interrupted runs
   - `--forcerun RULE [RULE ...]` — force re-execution of specific rules
   - `-- SNAKEMAKE_ARGS` — anything after `--` is forwarded directly to Snakemake

5. Optional: validate configuration and generate a workflow graph before running:
```bash
epicc validate --samples config/my_samples.tsv --dag dag.png
epicc validate --samples config/my_samples.tsv --rulegraph rules.png
```

*For a full list of epicc subcommands run `epicc --help` or `epicc <subcommand> --help`. For Snakemake's own options: https://snakemake.readthedocs.io/en/stable/*

### Conda environment maintenance

Rule-specific conda environments are defined under `workflow/envs/` and are created automatically on first use by Snakemake. When YAML files are updated (e.g. after a pipeline upgrade), orphaned environment directories may accumulate under `.snakemake/conda/`. To clean them up:
```bash
snakemake --sdm conda --conda-cleanup-envs
```

The shipped SLURM profile (`profiles/slurm/config.yaml`) sets `conda-cleanup-pkgs: tarballs`, which frees downloaded package archives after each environment is built without removing the installed environment itself.

#### Pre-building envs before sbatch-wrapped runs

On some clusters, `conda env create --prefix` fails when invoked from inside a SLURM job allocation, while the same command works from a login node. If you launch `epicc run` via `sbatch`, build all rule envs once from a login or dev node first:

```bash
epicc validate --build-envs --options config/epicc-options.yaml --samples your_samples.tsv
```

This runs the standard configuration checks and dry-run, then calls `snakemake --conda-create-envs-only` to populate `.snakemake/conda/`. Subsequent sbatch-wrapped runs will reuse the pre-built envs and skip the failing creation step.

### TMPDIR routing

By default, every pipeline job sets `TMPDIR` to a per-job subdirectory under `{output_dir}/.tmp/` (e.g. `results/.tmp/<SLURM_JOB_ID>.<PID>`). Tools that spill large temporary data through `TMPDIR` — such as samtools sort, STAR, fasterq-dump, and deeptools — therefore write to the project filesystem rather than the cluster's `/tmp`. This avoids `ENOSPC` errors on sites where `/tmp` is a tmpfs sized to the job's RAM allocation.

To disable this override, pass `--use-node-tmpdir` on the command line or set `use_node_tmpdir: true` in `config/epicc-options.yaml`. This does **not** mean "use `/tmp`" — it means EPICC sets no `TMPDIR` at all, and jobs inherit one from the environment. Precedence is then:

1. `TMPDIR` exported in the shell you launch from — it reaches the jobs, because the SLURM plugin submits with `sbatch --export=ALL`
2. `TMPDIR` set by the site: a SLURM prolog, `job_container/tmpfs`, or an environment module
3. `TMPDIR` set by the SLURM profile's `precommand` — this is the case where it finally takes effect, since nothing overwrites it afterwards
4. Each tool's own default, which is usually `/tmp`, only if nothing above set one

So directing temp work at node-local NVMe scratch, or at a filesystem that isn't under metadata pressure, is just `export TMPDIR=...` before the run. `fasterq-dump` is explicit about this (`--temp "${TMPDIR:-/tmp}"`), and sort, samtools, STAR and deeptools follow `TMPDIR` the same way.

Two rules opt out in either mode: `STAR_map_pe`/`STAR_map_se` deliberately keep temp work in the `genomes/` and `results/` trees to stay off any tmpfs, and Bismark uses its own `--temp_dir` under `results/mC/mapped/<sample>/`. So `use_node_tmpdir` reduces project-filesystem load without eliminating it.

Each scratch directory is removed by an `EXIT` trap when its job finishes, including on most failures — but `SIGKILL` escapes the trap, so jobs killed by a SLURM timeout, an OOM kill, or a node failure leak one each. To reclaim them:

```bash
epicc clean --tmp
```

This only removes a directory when both safety conditions hold: it does not belong to a SLURM job that is still queued or running, and it has not been modified for `--tmp-min-age` hours (default 6). Directories it leaves alone are listed with the reason. If `squeue` is unavailable, every SLURM-named directory is skipped rather than assumed dead.

Note that this routing trades cluster `/tmp` pressure for metadata load on the project filesystem — one directory created and removed per job, across as many jobs as run concurrently. On a shared parallel filesystem that is worth knowing about: if the filesystem administrators report metadata contention, `use_node_tmpdir: true` and a lower `jobs:` setting in the profile are the two levers that reduce EPICC's contribution.

### Launching from a login node, not from inside a job

Start `epicc run` from a login or dev node rather than wrapping it in `sbatch`. When Snakemake detects it is already inside a SLURM allocation, `snakemake-executor-plugin-slurm` deletes every `SLURM_*` environment variable — `SLURM_CONF` included. On clusters where `SLURM_CONF` is the only route to `slurm.conf`, every submitted job then fails with `srun: fatal: Could not establish a configuration source`. If submitting from within a job is unavoidable, restore the path via the profile's `precommand`; see the comments in `profiles/slurm/config.yaml` and `dev/docs/upstream_blockers.md`.

The same failure has a second cause that a login-node launch does *not* avoid: `snakemake-executor-plugin-slurm-jobstep` 0.6.1 strips `SLURM_CONF` from the environment it hands the inner `srun`, so every job fails no matter where the run was started, and `precommand` cannot restore it. `config/epicc-env.txt` pins the plugin below that release. If you built your environment before this pin was added, fix an existing one with `conda install -n epicc 'snakemake-executor-plugin-slurm-jobstep<0.6.1'`.

Note that this interacts with the pre-building step above: if `conda env create` fails inside allocations on your cluster, run `epicc validate --build-envs` from the login node and then launch the run from the login node too.

## Sample file configuration

### Overview

The sample metadata file is a tab-separated file (TSV) with 11 columns. Each row defines one biological sample.

| Sample_ID | Assay | Genome | Levels | Replicate_ID | Read_files | Read_layout | IP_target | Control | Peak_type | Comments |
|-----------|-------|--------|--------|--------------|------------|-------------|-----------|---------|-----------|----------|

The four trailing columns (`IP_target`, `Control`, `Peak_type`, `Comments`) are optional and may be left blank or omitted from the right-hand end of a row — a WGBS row can stop after `Read_layout`. See the documented template at `config/example_samples.tsv`.

### Common to all types of samples

- **Sample_ID**: A unique identifier for this sample. Used in output filenames. Must be filesystem-safe (no `__`, `/`, whitespace, or shell metacharacters).
- **Levels**: Comma-separated `factor:level` pairs describing experimental conditions (e.g. `genotype:WT,tissue:leaf` or `genotype:Col0,treatment:control`). The combination of levels is used for comparisons. All samples must have the same number of factors with the same factor names.
- **Replicate_ID**: Any value to identify replicates (e.g. `rep1`, `repA`, `1`). Replicates with the same Assay, Levels, IP_target, and Genome are merged for downstream analysis.
- **Read_files**: Path to input data. Supports:
  - SRA accession: `SRR27821931` (will be downloaded automatically)
  - Multiple SRA accessions to merge: `SRR27821931+SRR27821932` (separated by `+`)
  - Local FASTQ path: `/archive/fastq/sample_R1.fq.gz` (SE) or `/archive/fastq/sample_R1.fq.gz,/archive/fastq/sample_R2.fq.gz` (PE, comma-separated)
  - Local BAM path for dmC: `/archive/bams/sample.bam`
  - HTTP(S) URL: `https://example.org/data/sample_R1.fq.gz`
  - Public S3 URI: `s3://bucket/key/sample.bam`. **Public (authentication-free) objects only** — there is no request signing, no credentials, and no dependency on the `aws` CLI. Requester-Pays buckets (including SRA's Original-Format buckets such as `s3://sra-pub-src-*`) therefore will not work; use the plain SRA accession instead.
- **Read_layout**: `PE` for paired-end data or `SE` for single-end data.
- **Peak_type**: `broad` or `narrow`, for the three pulldown assays only (`ChIP`, `CUT_RUN`, `CUT_TAG`). This is an *analytical* choice — how peaks are called — kept separate from the assay, which records what was done at the bench. The same libraries can therefore be re-analyzed under the other peak type by editing one column. Leave blank for every other assay: ATAC's peak type is fixed, and the rest call no peaks. Peak type still appears in output filenames.
- **Comments**: Free text for your own notes (e.g. `resequenced 2026-03; higher duplication than rep2`). Completely inert — never read by the pipeline and never part of analysis names or output paths. Tabs and newlines are stripped on builder export so the TSV stays intact.
- **Genome**: Name of the reference genome (e.g. `ColCEN`, `Spombe`). To map one sample to several references, comma-separate them (e.g. `B73,W22`) — or tick them in the builder's Genome cell. The reads are downloaded and trimmed once, then aligned to each reference; every output from alignment onward is tagged with the genome (`final__WT_H3K9me2_rep1__B73.bam`), so the two analyses never collide. A genome name must not contain `__`. Each genome is defined under the `genomes:` namespace in `config/epicc-options.yaml`:
```yaml
genomes:
  ColCEN:
    genus: "Arabidopsis"
    species: "thaliana"
    fasta_file: "path/to/ColCEN.fa.gz"       # local path or HTTP(S) URL; .fa/.fasta(.gz)
    gff_file: "path/to/ColCEN_genes.gff3.gz" # local path or HTTP(S) URL; .gff*(.gz)
    #gtf_file: "path/to/ColCEN.gtf"          # optional; auto-derived from GFF via gffread if omitted
    te_file: "path/to/ColCEN_TE.gff3.gz"     # optional; .bed(.gz) or .gff3(.gz)
    gaf_file: "data/ColCEN_infoGO.tab.gz"    # only required when GO: true
    gene_info_file: "data/ColCEN_genes_info.tab.gz" # only required when GO: true
    ncbi_taxid: "3702"
    #genomesize: 1.3e8        # optional; auto-computed from FASTA if omitted
    #star_index: "--genomeSAindexNbases 12"  # optional; auto-computed if omitted
    #structural_rna_fafile: "<auto>"         # optional; auto-derived via Infernal/Rfam if omitted
```
Several fields are auto-derived at runtime and do not need to be specified in most cases: `gtf_file` (generated from GFF via gffread), `genomesize` (total non-N bases in FASTA), `star_index` parameters, `structural_rna_fafile` (fetched from Rfam via Infernal), and `ncbi_taxid` (looked up from NCBI Datasets CLI). Any value provided explicitly in the options file overrides the computed value.

All `*_file` fields accept local paths (absolute or relative to the repo root) or HTTP(S) URLs; gzipped inputs are handled transparently. The `te_file` field accepts `.bed(.gz)` (used as-is) or `.gff3(.gz)` (auto-converted to BED6 using the `ID=` attribute as the name column).

Reference genomes must be defined under the `genomes:` namespace shown above.

### ChIP-seq (Histones and Transcription Factors)

- **Assay**: `ChIP`.
- **Peak_type**: `broad` for histone marks with broad peaks (e.g. H3K9me2, H3K27me3) or `narrow` for marks with narrow peaks (e.g. H3K4me3) and transcription factors (e.g. TB1, FLC). The [ENCODE histone ChIP-seq target categorization](https://www.encodeproject.org/chip-seq/histone/) is a good starting point: broad-domain marks (H3K27me3, H3K36me3, H3K9me1/2, H3K79me2/3, H4K20me1, H3K4me1) use `broad`; punctate marks and TFs (H3K4me2/3, H3K27ac, H3K9ac, H2AFZ) use `narrow`. **Note that H3K27ac is narrow** despite being on the same residue as the broad mark H3K27me3. H3K9me3 is nominally broad but heavily enriched in repeats and benefits from a focused analysis — see ENCODE's note. Some targets profit from running both passes and comparing: RNA Pol II has sharp TSS peaks plus broader gene-body coverage, and H3K27me3 Polycomb domains can contain internal islands. Because peak type is its own column, running both passes means copying the rows and changing `Peak_type` — the assay label stays `ChIP`.
- **IP_target**: Required for all ChIP samples including controls. The name of what was pulled down — e.g. `H3K9me2` or `TB1` for IP, or `Input`/`WCE`/`IgG` for controls.
- **Control**: For IP samples, the Sample_ID of the control sample. Multiple IP samples can share the same control. A control may itself declare a `Control` — the *dual-role* case, where a sample is both another row's control and an analysis target in its own right (e.g. H3 serving as H3K9me2's control while itself being normalized against Input). Such a sample gets the full analysis treatment: peaks, fold-change tracks, and a merged analysis group. Chaining stops there: the control's control must not declare one. A pulldown row with no `Control` at all is not peak-callable and is dropped from peak and analysis targets.

Histone example:
```
WT_leaf_H3K9me2_rep1	ChIP	ColCEN	genotype:WT,tissue:leaf	rep1	SRR12345	PE	H3K9me2	WT_leaf_Input_rep1	broad
WT_leaf_Input_rep1	ChIP	ColCEN	genotype:WT,tissue:leaf	rep1	SRR12346	PE	Input		broad
```

Transcription factor example:
```
WT_leaf_TB1_rep1	ChIP	ColCEN	genotype:WT,tissue:leaf	rep1	SRR12347	PE	TB1	WT_leaf_Input_rep1	narrow
WT_leaf_Input_rep1	ChIP	ColCEN	genotype:WT,tissue:leaf	rep1	SRR12348	PE	Input		narrow
```

To use different controls for different marks (e.g. H3 for one mark, H4 for another), simply assign the appropriate control Sample_ID in the `Control` column.

- Option: Differential nucleosome sensitivity (DNS-seq) can be analyzed with `Assay=ChIP` and `Peak_type=broad`, using `MNase` as IP_target for the light digest and `Input` for the heavy digest.

### CUT&RUN / CUT&Tag

- **Assay**: `CUT_RUN` or `CUT_TAG`. Both route through the ChIP env.
- **Peak_type**: `broad` for diffuse marks (H3K27me3, H3K9me2, etc.); `narrow` for sharp marks and TFs (H3K4me3, CTCF, etc.).
- **IP_target**: Required for all CUT&x samples including controls. Use the antibody target (e.g. `H3K27me3`, `CTCF`) for IPs and `IgG` for controls.
- **Control**: For IP samples, the Sample_ID of the IgG control sample. Multiple IPs commonly share a single IgG (the typical CUT&RUN convention is one IgG per batch, not per replicate). Because control-row replicate merging keys on `(Levels, IP_target, Genome)` and ignores the assay label, one biological IgG can serve broad and narrow IPs — and CUT&RUN and CUT&Tag IPs — at once.
- **Peak callers**: defaults are peak-shape-aware — `broad` → epic2, `narrow` → SEACR. MACS2 is available as a fallback. Override via `cut_callpeaks.broad_caller` / `narrow_caller` in `config/epicc-options.yaml`.

CUT&RUN example (sharing a single IgG across reps):
```
WT_endo_H3K27me3_rep1	CUT_RUN	ColCEN	genotype:WT,tissue:endosperm	rep1	SRR8310960	PE	H3K27me3	WT_endo_IgG_rep1	broad
WT_endo_H3K27me3_rep2	CUT_RUN	ColCEN	genotype:WT,tissue:endosperm	rep2	SRR8310958	PE	H3K27me3	WT_endo_IgG_rep1	broad
WT_endo_IgG_rep1	CUT_RUN	ColCEN	genotype:WT,tissue:endosperm	rep1	SRR8310961	PE	IgG		broad
```

CUT&Tag for a TF (single-end, defaulting to SEACR for narrow peak calling):
```
WT_leaf_CTCF_rep1	CUT_TAG	hg38	genotype:WT,tissue:leaf	rep1	SRR12345	SE	CTCF	WT_leaf_IgG_rep1	narrow
WT_leaf_IgG_rep1	CUT_TAG	hg38	genotype:WT,tissue:leaf	rep1	SRR12346	SE	IgG		narrow
```

### RNA-seq

- **Assay**: `RNAseq`
- **IP_target**: Leave blank.
- **Control**: Leave blank.

### RAMPAGE

- **Assay**: `RAMPAGE`
- **IP_target**: Leave blank.
- **Control**: Sample_ID of the corresponding RNAseq sample (used for normalization).

### small RNA-seq

- **Assay**: `sRNA`
- **IP_target**: Leave blank.
- **Control**: Leave blank.

### Whole Genome Bisulfite Sequencing

- **Assay**: `WGBS`, `WGBS_nd`, `PBAT`, or `EMseq`. These labels identify the conversion chemistry and inform Bismark alignment parameters:
  - `WGBS` — standard directional WGBS (e.g. TruSeq, Swift HTP)
  - `WGBS_nd` — non-directional WGBS (e.g. Zymo Pico Methyl-Seq, Swift Accel-NGS)
  - `PBAT` — post-bisulfite adapter tagging
  - `EMseq` — enzymatic methyl-seq (NEB EM-seq kit)
- **IP_target**: Leave blank.
- **Control**: Leave blank.

### Direct Methylation (Long-Read Sequencing)

- **Assay**: `dmC`. For samples with native base modifications (Oxford Nanopore, PacBio) that have not undergone bisulfite conversion.
- **Read_files**: Path to input file or directory. Supports:
  - **modBAM**: BAM files with MM/ML methylation tags from basecalling (e.g., Dorado, Guppy)
  - **bedMethyl**: Pre-computed methylation calls in bedMethyl format (e.g., from modkit pileup)
- **Read_layout**: Only `SE` is currently supported (long-read sequencing).
- **IP_target**: Leave blank.
- **Control**: Leave blank.
- Note: The pipeline automatically detects modBAM vs bedMethyl format. modBAM files are aligned and processed through modkit pileup. Both formats are converted to Bismark-compatible CX_report format for downstream analysis.

## Configuration Options

More details can be found in the `config/epicc-options.yaml` file (all options are commented inline) or via the local HTML builder at `tools/epicc-builder.html`.

### Main output options
- `full_analysis`: When `false`, only the mapping and the bigwigs will occur. When `true` (default), will also be performed: single-data analyses (e.g. peak calling for ChIP, differential expression for RNAseq, DMRs for mC) and combined analyses (e.g. Upset plots for ChIP, heatmaps and metaplots on all genes).
- `te_analysis`: When `true`, generates additional combined-analysis heatmaps and metaplots over the configured `te_file` annotation in addition to the standard gene-centered plots. Requires `te_file` to be set under the corresponding genome block in the options file (BED6 or GFF3, gzipped or not — see Sample file configuration). The name of each feature (4th column of the BED, or the `ID=` attribute for GFF3) must be unique. Default is `false`.
- `QC_option`: Controls FastQC reporting on raw and trimmed FASTQ files. `"none"` (default) skips FastQC entirely; `"all"` runs FastQC on all samples.

### Intermediate input formats
- `trimmed_fastqs`: When `false` (default), the analysis runs from raw, untrimmed fastq files and performs adapter trimming. If you already have trimmed fastqs, you can switch this config entry to `true` and no additional trimming will be performed (still compatible with nextflex_v3 deduplication and structural RNAs filtering for small RNAs).
- `aligned_bams`: When `true`, the pipeline expects pre-aligned ChIP-seq BAM/SAM files in the `Read_files` column of the sample sheet rather than FASTQs (still applies pipeline-wide, so all ChIP samples must follow the same convention). No mapping-stats plot is available when starting from BAMs this way. Default is `false`. Currently only supported for ChIP-seq assays.
- Note: These settings are applied to *all* samples in the analysis. If you have some samples to analyze from scratch and other already in an intermediate file:
	- 1) run the pipeline once with the samples to run from scratch - potentially switching `full_analysis` to `false` for less output.
	- 2) add the samples you already have intermediate files for in the samplefile and change the corresponding parameters in the options file.
	- 3) run the pipeline normally again.
	These steps can be repeated if you have raw data, trimmed fastqs and bam files, first creating all the fastq files and then the bam files.

###  Intermediate Target Rules
- `map_only`: Only performs the alignement of all samples. It returns bam files, QC files and mapping metrics.
- `coverage_chip`: Creates bigwig files of coverage for all ChIP samples. The binsize is by default 1bp (can be updated in the options file `chip_tracks: binsize: 1`).

###  Plotting parameters
- `plot_allreps`: When `true`, all individual replicates are shown on heatmaps, metaplots and browsers (can be heavy). When `false` (default), one sample with all merged replicates is used for each sample.

### ChIP Mapping Parameters
- `default`: Standard mapping parameters
- `repeat`: Centromere-specific mapping (more sensitive)
- `repeatall`: Centromere mapping with relaxed MAPQ
- `all`: Relaxed mapping parameters

### DMRs parameters
- By default, DNA methylation data is analyzed in all three sequence contexts (CG, CHG, and CHH, where H = A, T, or C). This is controlled by the `methylation_contexts` list in `config/epicc-options.yaml` (default: `["CG", "CHG", "CHH"]`). For animal genomes where non-CpG methylation is negligible, set `methylation_contexts: ["CG"]` to skip the empty CHG/CHH bigwigs, DMR calls, and PCA plots. Subcontexts (CAG, CAA, etc.) are not currently supported.
- DMR calling thresholds are configured via the `dmr_thresholds:` block in `config/epicc-options.yaml`. Defaults are tuned for Arabidopsis leaves ([Catoni et al. 2018](https://academic.oup.com/nar/article/46/19/e114/5050634)); cross-species and cross-tissue methylation varies widely (e.g. maize CHH is very low; pollen differs from leaf), so adjust as needed:
  ```yaml
  dmr_thresholds:
    min_diff:       {CG: 0.3, CHG: 0.2, CHH: 0.1}  # min mean methylation difference
    min_cytosines:  5      # min CpG/cytosine count per DMR window (both callers)
    q_value:        0.01   # FDR q-value significance cutoff (both callers)
    min_size:       50     # min DMR length in bp (both callers)
    bin_size:       200    # DMRcaller: window size (bp)
    min_gap:        200    # DMRcaller: min gap between DMRs (bp)
    min_reads:      3      # DMRcaller: min reads per cytosine
    maxdist:        300    # metilene: max distance between CpGs in a segment, nt (-M)
    valley:         {CG: 0.7, CHG: 0.7, CHH: 0.3}  # metilene: valley filter (-v)
    maxseg:         {CG: -1, CHG: -1, CHH: 10000}   # metilene: max segment cap (-G)
  ```
  `min_diff` and `min_cytosines` apply to both callers; the remaining keys are caller-specific and ignored by the other.
- To let the pipeline estimate `min_diff` automatically from the data, set `min_diff: "auto"` (string). For each comparison pair and context a short calibration step runs first and computes `sigma_n × σ`, where σ is the SD of per-position between-group methylation differences (noise-floor dominated at most genome positions). This scales naturally with baseline methylation: CG positions are typically highly methylated (high σ → higher threshold) while CHH positions are near zero (low σ → lower threshold), so per-context differentiation emerges without hard-coded values. A floor of 0.05 prevents pathological values at extreme coverage. Tune `sigma_n` (default 3.0) to tighten or relax sensitivity. `auto` is recommended when working outside Arabidopsis leaf tissue — maize or duckweed somatic CHH methylation is far lower than Arabidopsis, and pollen/endosperm methylation landscapes differ substantially from leaf.
- A deeper analysis is available to try different parameters and methods to call the DMRs. Toggle the `use custom_script_dmrs` on the options file to use it. Feel free to edit it as well for different parameters.

##  Additional output options
Below is a list of *cool* outputs that can be generated once the whole pipeline has run at least once. Use `epicc output --plot-type TYPE --input-file PATH [options]` to generate them, or pass the raw Snakemake target file directly via `snakemake --cores 1 <target>` for advanced use.

### **1. Plotting RNAseq expression levels on target genes**
Given a list of genes (and optional labels), it will plot the expression levels in all the different samples in the samplefile and analysis name defined. Genes uniquely differentially regulated in one sample versus one or more samples are color coded. It is based on a Rdata file created during the Differential Expression analysis.\
To run it, provide the target gene list file (one column of gene IDs matching the GTF of the reference genome used; an optional second column provides gene labels) and run:
```bash
epicc output --plot-type rnaseq-histogram \
  --input-file data/target_genes.txt \
  --plot-label my_genes_of_interest \
  --ref-genome TAIR10
```
Output is a single pdf file named `results/RNA/plots/plot_expression__<analysis_name>__<ref_genome>__<label>.pdf` where each gene of the list is on an individual page. The separator between variables is two underscores `__`.

### **2. Performing GO analysis on target genes**
Given a file containing a list of genes to do GO analysis on, and optionally a background file (default to all genes in the reference genome), it will perform Gene Ontology analysis.\
By default, GO is not performed since it requires manual input to build a database. To activate it, `GO` needs to be switched to `true` in the options file, and `gaf_file` and `gene_info_file` must be defined under the corresponding genome block in the options file. See [Gene_Ontology.md](Help/Gene_Ontology.md) for more details on how to create the GO database.\
To run it, provide the gene list file and run:
```bash
epicc output --plot-type go \
  --input-file data/target_genes.txt \
  --plot-label my_genes_of_interest \
  --ref-genome ColCEN
```
By default the enrichment background (the reference gene set tested against) is
all genes in the genome. To restrict it — e.g. to only the genes expressed in
your experiment — pass a gene-ID file with `--background`:
```bash
epicc output --plot-type go \
  --input-file data/target_genes.txt \
  --plot-label my_genes_of_interest \
  --ref-genome ColCEN \
  --background data/expressed_genes.txt
```
(`--background` sets `rnaseq_background_file` for this run; it is a tab-delimited
file whose first column is gene IDs matching the reference GFF.)
Output includes two pdf treemaps (`topGO_<label>_BP_treemap.pdf` for biological process and `topGO_<label>_MF_treemap.pdf` for molecular function) and corresponding enrichment tables under `results/RNA/GO/`.

### **3. Finding motifs on target regions**
Given a bed file containing different regions, it will perform a motifs analysis with meme.\
By default motifs analysis is only performed on the final selected TF peak files (`motifs: true` in the options file). Edit to `motifs_allreps: true` in the options file for motifs analysis to be performed on all replicates and pairwise IDR peaks if available. A plant motifs database is used by default for tomtom. Download the appropriate file from JASPAR and replace its name in the options file `jaspar_db` and change the `motif_ref_genome` to match the samples.\
To run the analysis:
```bash
epicc output --plot-type motifs \
  --input-file data/target_peaks.bed \
  --plot-label my_regions_of_interest \
  --ref-genome ColCEN
```
Output is the folder `results/ChIP/<label>` containing a subdirectory called `meme` and potentially one called `tomtom` with all the results, as described in https://meme-suite.org/meme/index.html.\
Use a reference genome that has already been used in a prior run for the motif analysis. For the target file, if the regions are over 500bp, only the middle 400bp will be used.

### **4. Performing sRNA differential analysis on regions**
Given a bed or gff file, it will perform the small RNA analysis with shortstack followed by differential analysis with edgeR, using all the samples from the sample file but limiting the mapping and counts to the loci in the target file. To run the analysis:
```bash
epicc output --plot-type srna-clusters \
  --input-file data/miRNA.gff \
  --plot-label miRNAs \
  --ref-genome ColCEN
```
Output is the results folder from Shortstack limited to this loci file, followed by the differential cluster analysis with edgeR. Output path: `results/sRNA/clusters/<analysis_name>__<ref_genome>__on_<label>/Counts.txt`.

The bed or gff file of regions **MUST HAVE** a header with a column called "Name" (the 4th column of a bed file or the 9th column of a gff3).

### **5. Plotting heatmap on regions**
Given a bed file, it will plot a heatmap using deeptools.
Edit `heatmap_target_file` and `heatmap_target_file_label` in the options file, or pass them on the command line. To run the analysis:
```bash
epicc output --plot-type heatmap \
  --input-file data/target_genes.bed \
  --plot-label interesting_genes \
  --ref-genome ColCEN \
  --matrix regions \
  --env most
```
- `--matrix` can be `regions` (scaled features), `tss` (reference point on TSS), or `tes` (reference point on TES).
- `--env` selects the data types to include: `most` includes all data types except mC; use `mC`, `mCG`, `mCHG`, or `mCHH` for methylation-specific heatmaps; other choices are `ChIP`, `ATAC`, `RNA`, `sRNA`.

Since mC requires different deeptools parameters it is handled independently. If you want the mC heatmap sorted by the same region order as the other samples, use the `Heatmap_sorted__` raw Snakemake target (see Snakemake docs for advanced usage).

Output is a pdf file, or two if sorted heatmap for mC samples was generated.\
By default, heatmaps are scaled by data type (`heatmaps_scales: "type"` in the options file; each ChIP mark, TF, RNA, sRNA size, and mC context on its own scale). Change to `"default"` for a single scale or `"sample"` for per-sample scaling.\
By default, regions are sorted by mean signal (`heatmaps_sort_options: "mean"`). Change to `"median"` or `"no"` to preserve the bedfile order.\
If the bedfile is stranded, heatmaps will split plus/minus strand for properly stranded data types (RNAseq, sRNA). Disable with `stranded_heatmaps: false` in the options file.\
Color scheme defaults: `seismic` for all non-mC samples, `Oranges` for mC. Change via `heatmaps_plot_params` in the options file.
Window sizes (`before`, `after`, scaled-region `middle`) and bin size (`binsize`) are configurable in `heatmaps` per matrix type in the options file.

### **6. Plotting metaplot profiles on regions**
Given a bed file, it will plot a metaplot profile using deeptools.
Edit `heatmap_target_file` and `heatmap_target_file_label` in the options file (the same target file fields as heatmap). To run the analysis:
```bash
epicc output --plot-type metaplot \
  --input-file data/target_genes.bed \
  --plot-label interesting_genes \
  --ref-genome ColCEN \
  --matrix regions \
  --env all
```
`--matrix` and `--env` options are identical to the heatmap above. Use `--env all` to include all samples including mC.\
Output is two pdf files, where the samples are grouped by regions or not.\
By default, profiles are scaled by data type (`heatmaps_scales: "type"`; change to `"default"` or `"sample"`).
By default, profiles represent the mean across all regions (`profiles_scale: "mean"`; change to `"median"`).
Plot type defaults to `"lines"`. See deeptools documentation for other options; edit `profiles_plot_params` in the options file.
Window sizes and bin size are configured in `heatmaps` per matrix type in the options file (same as heatmap).

### **7. Plotting browser screenshots on regions**
Given a region file, it will plot a browser screenshot using R packages.
Edit `browser_target_file` and `browser_target_file_label` in the options file. To run the analysis:
```bash
epicc output --plot-type browser \
  --input-file data/target_loci.bed \
  --plot-label target_loci \
  --ref-genome ColCEN \
  --env all
```
The target file is a bed-like file with columns: Chr Start End ID Binsize Highlight_starts Highlight_widths.\
Each region will be printed individually and merged into a final PDF.\
Highlight columns are optional; they mark regions of the browser view with colored boxes. As many highlights as needed can be provided as comma-separated lists — the first highlight is blue and all others are red. For example, for a region chr1:1000-5000, col6=3000,4000 and col7=50,200 will produce a blue box at chr1:3000-3050 and a red box at chr1:4000-4200.\
Use `--env all` to include all samples, `most` for all data types except mC, or any single environment `[ChIP, ATAC, RNA, sRNA, mC]`.\
By default, no TE file is used. To add TE annotations, supply a bed file in the options file under `browser_TE_file`.

### **8. Rerunning a specific analysis**
To force Snakemake to recreate a specific output, use `epicc run -- <target_file> --forcerun <rule>` or pass the target via `epicc run -- --force <target_file>`. For example:
```bash
epicc run -- results/combined/plots/srna_sizes_stats_test_snakemake_sRNA.pdf --force
```
If only the combined analysis needs to be redone (not the data-type-specific steps), delete the checkpoint files in `results/combined/chkpts/` and in the relevant environment directory `results/<env>/chkpts/<env>_analysis__<analysis_name>__<ref_genome>.done`. Alternatively, use `epicc clean --intermediates` to remove all checkpoints and trigger full reanalysis on the next run.

## Output Structure

```
EPICC/
├── config/			# Main options file and recommended location for sample files and target files
├── data/			# Test material and examples (e.g. zm_structural_RNAs.fa.gz)
├── Help/			# Help files (e.g. Structural_RNAs_Rfam.md, Gene_Ontology.md)
├── profiles/
│	├── default/		# Workflow-level per-rule resource/thread defaults
│	├── geno/		# Example profile for additional scheduler types
│	├── slurm/		# Config file to run the pipeline on a SLURM cluster
│	└── uge/		# Config file to run the pipeline on a UGE cluster (qsub)
├── tools/
│	└── epicc-builder.html	# Self-contained HTML5 sample sheet and options builder
├── workflow/
│	├── envs/		# Conda environment YAML files per analysis type
│	├── rules/		# Snakemake rule files by data type
│	├── scripts/		# R and Python scripts for analysis and plotting
│	└── Snakefile		# Main Snakefile
├── genomes/			# Genome directories created upon run
│	└── {ref_genome}/	# Reference genome directories with sequence, annotation and indexes
└── results/			# Results directories created upon run
	├── .tmp/		# Per-job TMPDIR scratch space (auto-cleaned after each job)
	├── combined/		# Combined analysis results
	│	├── bedfiles/	# Peak calling results
	│	├── chkpts/	# Checkpoint files; delete to trigger rerunning the corresponding analysis
	│	├── logs/	# Log files
	│	├── matrix/	# Data matrices
	│	├── plots/	# Visualization plots
	│	└── reports/	# Analysis reports 
	└── <env>/		# Data type specific directories (ChIP, ATAC, RNA, sRNA, mC)
		├── chkpts/	# Checkpoint files; delete to trigger rerunning the corresponding analysis
		├── fastq/	# Processed FASTQ files
		├── logs/	# Log files
		├── mapped/	# Mapped reads (BAM)
		├── plots/	# Data type specific plots
		├── reports/	# QC reports
		├── tracks/	# Track files (bigwigs)
		└── */		# Data-specific directories (e.g. peaks/ for ChIP, DEG/ for RNA, DMRs/ and methylcall/ for mC, clusters/ for sRNA)
```

## Known potential issues

1. Relationship between IP and Input for ChIP-seq samples\
Control samples are linked via the `Control` column in the sample sheet (must be a valid Sample_ID). IP and Control samples must use the same Read_layout (both PE or both SE).

2. small RNA-seq libraries\
Different small RNAseq libraries have different chemistry and might need to be trimmed differently. For now, the code only works if all your samples were done using the same library preparation, either netflex v3 or not. If you have a mix of libraries, you should run the pipeline with each kind separately, and then rerun the analysis with all the samples you want to anlayze together.

3. idr/numpy version\
IDR relies on an older version of numpy to work (due to deprecated np.int) and needs to be loaded as a seperate environment. Not best practice, but more portable than patching idr (np.int=int).

4. SLURM wall-time limits\
The shipped `profiles/slurm/config.yaml` sets a default runtime of 60 minutes per job; individual rules override this with their own `runtime` resource. If your cluster has tighter limits, increase the `default-resources: runtime` value in the profile, or add per-rule overrides under `set-resources`. If some jobs require a specific Quality-of-Service (QoS) setting on your cluster, add `slurm_extra: "'--qos=<your_qos>'"` under the relevant rule in `set-resources`.

5. Help for local fastq files naming convention\
If using local fastq files for paired-end data, provide comma-separated R1 and R2 paths in the `Read_files` column (e.g. `/path/sample_R1.fq.gz,/path/sample_R2.fq.gz`). Files can use extensions `.fq` or `.fastq` and may be gzipped (`.gz`).

## FAQ

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License — a permissive license that
allows use, modification, and distribution (including commercially), provided
the copyright notice and license text are retained.

See the [LICENSE](LICENSE) file for full details.

## Citation

If you use EpigeneticButton in your research, please cite:

```
Cahn, J., Regulski, M., Lynn, J. et al. MaizeCODE reveals bi-directionally expressed enhancers that harbor molecular signatures of maize domestication. Nat Commun 15, 10854 (2024). https://doi.org/10.1038/s41467-024-55195-w
```

## References 

This pipeline is only a combination of great tools developped by others. A non-exhaustive list of packages used are listed below. Please refer to them for more details.
- [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [ComplexUpset](https://krassowski.github.io/complex-upset/index.html)
- [Conda](https://anaconda.org/anaconda/conda)
- [fastp](https://github.com/OpenGene/fastp)
- [deepTools](https://deeptools.readthedocs.io/en/develop/index.html)
- [DMRcaller](https://bioconductor.org/packages/release/bioc/html/DMRcaller.html)
- [edgeR](https://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [IDR](https://github.com/nboley/idr)
- [MACS2](https://pypi.org/project/MACS2/#description)
- [Python](https://www.python.org/)
- [R](https://www.r-project.org/)
- [Samtools](https://www.htslib.org/)
- [ShortStack](https://github.com/MikeAxtell/ShortStack)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [SRA-Toolkit](https://github.com/ncbi/sra-tools)
- [STAR](https://github.com/alexdobin/STAR)
- [The MEME suite](https://meme-suite.org/meme/doc/meme-chip.html?man_type=web)
- [topGO](https://bioconductor.org/packages/release/bioc/html/topGO.html)
- [UCSC-GenomeBrowser-kent](https://github.com/ucscGenomeBrowser/kent/)

## Contact

For questions or support, please open an issue in the GitHub repository.
