# S. pombe Test Case Design

## Objective

A multi-omics fission yeast test case covering 4 assay types: ChIP_broad (PE),
ChIP_narrow (SE), RNAseq (SE), and sRNA (SE). Uses the full S. pombe genome
(~12.6 Mbp, 3 chromosomes). Reference files are bundled in the repository
under `tests/integration/data/Spombe/`. This is the primary dry-run test case
used for CI regression testing.

## Reference Genome

**Assembly**: ASM294v2 (PomBase)
**Source**: Bundled in `tests/integration/data/Spombe/`

| File | Description |
|------|-------------|
| `Schizosaccharomyces_pombe_all_chromosomes.fa.gz` | Genome FASTA |
| `Schizosaccharomyces_pombe_all_chromosomes.gff3.gz` | Gene annotations |
| `Schizosaccharomyces_pombe_all_chromosomes.gtf.gz` | GTF (pre-computed) |
| `sp_structural_RNAs.fa.gz` | Structural RNA sequences for sRNA filtering |

## Dataset 1: ChIP_broad (H3K9me2, H3K9me3) — Kim et al. 2024

**Source**: Kim HS, Roche B, Martienssen RA (2024) "Clr4SUV39H1 ubiquitination
and non-coding RNA mediate transcriptional silencing via heterochromatic phase
transitions." PubMed: 39477922.

**GEO**: [GSE156069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156069)
(76 samples total; subset used here)

### Experimental design

- **Factor**: `genotype` — WT vs dcr1
- **Marks**: H3K9me2, H3K9me3
- **Control**: WCE (whole-cell extract)
- **Layout**: PE
- **Replicates**: 2 for WT, 1 for dcr1

| Sample_ID | Genotype | Mark | SRR |
|-----------|----------|------|-----|
| WT_ChIP_H3K9me2_rep1 | WT | H3K9me2 | SRR20678305 |
| WT_ChIP_H3K9me2_rep2 | WT | H3K9me2 | SRR20678333 |
| WT_ChIP_H3K9me3_rep1 | WT | H3K9me3 | SRR20678311 |
| WT_ChIP_H3K9me3_rep2 | WT | H3K9me3 | SRR20678336 |
| dcr1_ChIP_H3K9me2_rep1 | dcr1 | H3K9me2 | SRR20678308 |
| dcr1_ChIP_H3K9me3_rep1 | dcr1 | H3K9me3 | SRR20678314 |
| WT_ChIP_WCE_rep1 | WT | WCE (control) | SRR5445712 |

### Design note

The WCE control (SRR5445712) is shared across all ChIP_broad samples. This is
standard practice for fission yeast where a single input control serves all
ChIP experiments from the same strain background. The dcr1 samples use the WT
WCE as control since dcr1 and WT share the same genome.

## Dataset 2: ChIP_narrow (H3K4me3) — Zeng & Ekwall 2024

**Source**: Zeng S, Ekwall K (2024) "Epigenome mapping in quiescent cells
reveals a key role for H3K4me3 in regulation of RNA polymerase II activity."
PubMed: 39449363.

**GEO**: [GSE280066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280066)

### Experimental design

- **Factor**: `genotype` — WT only
- **Mark**: H3K4me3
- **Control**: Input
- **Layout**: SE
- **Replicates**: 2

| Sample_ID | Mark | SRR |
|-----------|------|-----|
| WT_ChIP_H3K4me3_rep1 | H3K4me3 | SRR31073789 |
| WT_ChIP_H3K4me3_rep2 | H3K4me3 | SRR31073812 |
| WT_ChIP_Input_rep1 | Input (control) | SRR31073785 |

### Design rationale

H3K4me3 is a narrow-peak mark (active promoters), exercising the ChIP_narrow
code path with MACS2 narrow peak calling. Using SE data from a different lab
than the broad-peak ChIP data tests both read layout paths.

## Dataset 3: RNA-seq — Cheng & Martienssen 2025

**Source**: Cheng T, Martienssen RA. "Transcription-Replication Conflict
Resolution by Nuclear RNA Interference."

**GEO**: [GSE278839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278839)

### Experimental design

- **Factor**: `genotype` — WT vs dcr1
- **Layout**: SE
- **Replicates**: 2 per genotype

| Sample_ID | Genotype | SRR |
|-----------|----------|-----|
| WT_RNA_rep1 | WT | SRR30889044 |
| WT_RNA_rep2 | WT | SRR30889043 |
| dcr1_RNA_rep1 | dcr1 | SRR30889039 |
| dcr1_RNA_rep2 | dcr1 | SRR30889038 |

## Dataset 4: sRNA-seq — Kim et al. 2024

**Source**: Same as Dataset 1 (GSE156069).

### Experimental design

- **Factor**: `genotype` — WT vs dcr1
- **Layout**: SE
- **Replicates**: 2 for WT, 1 for dcr1

| Sample_ID | Genotype | SRR |
|-----------|----------|-----|
| WT_sRNA_rep1 | WT | SRR20678362 |
| WT_sRNA_rep2 | WT | SRR20678373 |
| dcr1_sRNA_rep1 | dcr1 | SRR20678363 |

### Design note

dcr1 is a Dicer mutant — sRNA production is severely impaired. This is
biologically expected and exercises the pipeline's handling of low-complexity
sRNA libraries.

## Assay Coverage Matrix

| Assay | Genotypes | Layout | Reps (WT / dcr1) | Samples | Source |
|-------|-----------|--------|-------------------|---------|--------|
| ChIP_broad (H3K9me2) | WT, dcr1 | PE | 2 / 1 | 3 | GSE156069 |
| ChIP_broad (H3K9me3) | WT, dcr1 | PE | 2 / 1 | 3 | GSE156069 |
| ChIP_broad (WCE ctrl) | WT | PE | 1 | 1 | GSE156069 |
| ChIP_narrow (H3K4me3) | WT | SE | 2 | 2 | GSE280066 |
| ChIP_narrow (Input ctrl) | WT | SE | 1 | 1 | GSE280066 |
| RNAseq | WT, dcr1 | SE | 2 / 2 | 4 | GSE278839 |
| sRNA | WT, dcr1 | SE | 2 / 1 | 3 | GSE156069 |

**Total**: 17 samples across 4 assay types (+ controls).

## Known Design Considerations

### Shared WCE control

A single WT WCE sample controls all 6 ChIP_broad IP samples (both WT and dcr1,
both H3K9me2 and H3K9me3). The pipeline handles shared controls correctly —
the control BAM is processed once and referenced by all IP samples.

### Asymmetric replicates

dcr1 has only 1 replicate for ChIP_broad and sRNA (WT has 2). This exercises
the pipeline's handling of analysis groups with different replicate counts.
IDR analysis requires 2 replicates and will only run for WT.

### Small genome advantage

At ~12.6 Mbp, S. pombe is ideal for fast dry-run and integration testing.
Full pipeline runs complete in minutes rather than hours. The genome is
bundled in the repository to avoid download dependencies during CI.
