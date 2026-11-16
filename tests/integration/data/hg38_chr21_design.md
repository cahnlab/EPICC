# Human Chr21 Test Case Design

## Objective

A multi-omics human test case covering all pipeline assay types: ChIP_broad,
ChIP_narrow, RNAseq, RAMPAGE, sRNA, WGBS, and dmC. All data is subsetted to
chromosome 21 (~46.7 Mbp, ~1.5% of the genome) to keep file sizes manageable.

## Core Dataset: EN-TEx (ENCODE)

**Source**: Rozowsky et al. (2023) "The EN-TEx resource of multi-tissue personal
epigenomes and variant-impact models." *Cell* 186(6), 1493-1511.
doi:10.1016/j.cell.2023.02.018

The EN-TEx project profiled ~30 tissues across 4 human donors with 13+ assay
types per tissue. All data is publicly available via ENCODE and SRA.

### Experimental design

- **Factor**: `tissue:stomach` vs `tissue:colon` (transverse colon)
- **Replicates**: Different donors serve as biological replicates
  - D1 = ENCDO271OUW (51F)
  - D2 = ENCDO793LXB (53F)
  - D3 = ENCDO451RUA (54M) — RNA-seq only (3rd replicate)
- **Genome**: hg38 (GRCh38), subsetted to chr21

### Design rationale

Using donors as replicates captures biological variation (genetics, age, sex)
within the tissue comparison. This is appropriate for a pipeline test case and
mirrors real-world multi-donor study designs. D3 provides a 3rd RNA-seq
replicate for robust DEG analysis; other assays have 2 replicates per tissue.

## Nanopore Methylation: ONT Open Data (GIAB HG002)

**Source**: Kolmogorov et al. (2023) "Scalable Nanopore sequencing of human
genomes provides a comprehensive view of haplotype-resolved variation and
methylation." *Nature Methods* 20, 1483-1492.

- **Sample**: HG002 (GM24385, Ashkenazi son, GIAB reference)
- **S3 path**: `s3://ont-open-data/giab_2025.01/flowcells/hg002/`
- **Chemistry**: R10.4.1, Kit V14 (SQK-LSK114), Dorado v5.0.0
- **Format**: Modbase-tagged BAM with MM/ML tags (5mC + 5hmC)
- **Coverage**: ~37X whole-genome

This dataset does not share the tissue factor with the EN-TEx samples
(tissue:cell_line). It exercises the dmC pipeline path as a standalone sample.

The DeepMod2 dataset (PRJNA951714) was evaluated and rejected: RRMS (not WGS),
raw signal only (no modBAM), and insufficient Chr21 coverage.

## Assay Coverage Matrix

| Assay | Stomach | Colon | Layout | Reps per tissue |
|-------|---------|-------|--------|-----------------|
| ChIP_broad (H3K27ac) | D1, D2 | D1, D2 | SE 76bp | 2 |
| ChIP_narrow (CTCF) | D1 (PE), D2 (SE) | D1 (PE), D2 (SE) | Mixed | 2 |
| Control ChIP (H3K27ac) | D1, D2 | D1, D2 | SE 76bp | 2 |
| Control ChIP (CTCF) | D1 (PE) | D1 (PE) | PE 100bp | 1 |
| RNAseq | D1, D2, D3 | D1, D2, D3 | PE 101bp | 3 |
| RAMPAGE | D1, D2 | D1, D2 | PE 101bp | 2 |
| sRNA | D1, D2 | D1, D2 | SE 76/101bp | 2 |
| WGBS | D1, D3 | D2, D3 | PE 151bp | 2 |
| dmC (ONT) | — | — | SE (long reads) | 1 (cell_line) |

**Total**: 33 samples across 7 assay types (+ controls).

## ENCODE Experiment Accessions

### Donor 1 (D1): ENCDO271OUW (51F)

| Assay | Tissue | ENCODE Exp | Control Exp | SRR(s) |
|-------|--------|------------|-------------|--------|
| H3K27ac ChIP | stomach | ENCSR751BHO | ENCSR948DNK | SRR10838362+363+364+365 |
| H3K27ac ChIP | colon | ENCSR792VLP | ENCSR384KIB | SRR10838452+453+454+455 (rep1 only) |
| CTCF ChIP | stomach | ENCSR173AIR | ENCSR955OCL | SRR6212897+SRR6212898 (PE) |
| CTCF ChIP | colon | ENCSR558HTE | ENCSR457NQR | SRR6213835 (PE) |
| Control (H3K27ac) | stomach | ENCSR948DNK | — | SRR10838592+593+594+595 (rep1) |
| Control (H3K27ac) | colon | ENCSR384KIB | — | SRR10837838+839+840+841 (rep1) |
| Control (CTCF) | stomach | ENCSR955OCL | — | SRR6214083 (PE) |
| Control (CTCF) | colon | ENCSR457NQR | — | SRR6213596 (PE) |
| RNA-seq | stomach | ENCSR853WOM | — | SRR4422373 (PE) |
| RNA-seq | colon | ENCSR403SZN | — | SRR4421756 (PE) |
| RAMPAGE | stomach | ENCSR497BYB | ENCSR853WOM | SRR4421872 (PE) |
| RAMPAGE | colon | ENCSR297QGC | ENCSR403SZN | SRR5111381 (PE) |
| sRNA | stomach | ENCSR573RGL | — | SRR4421975 (SE 76bp) |
| sRNA | colon | ENCSR183YQY | — | SRR4421452 (SE 101bp) |
| WGBS | stomach | ENCSR598OBF | — | SRR8659907 (PE 151bp) |

### Donor 2 (D2): ENCDO793LXB (53F)

| Assay | Tissue | ENCODE Exp | Control Exp | SRR(s) |
|-------|--------|------------|-------------|--------|
| H3K27ac ChIP | stomach | ENCSR133NBJ | ENCSR539HNK | SRR5832694+SRR5832695 |
| H3K27ac ChIP | colon | ENCSR208QRN | ENCSR828XQV | SRR5832686+SRR5832687 (rep1 only) |
| CTCF ChIP | stomach | ENCSR185CCV | ENCSR539HNK | SRR10837564+565+566+567 (SE) |
| CTCF ChIP | colon | ENCSR236YGF | ENCSR828XQV | SRR10837644+645+646+647 (SE, rep1) |
| Control (shared) | stomach | ENCSR539HNK | — | SRR5832845+846+847+848 (SE) |
| Control (shared) | colon | ENCSR828XQV | — | SRR5833123+124+125+126 (SE, rep1) |
| RNA-seq | stomach | ENCSR752UNJ | — | SRR4422210 (PE) |
| RNA-seq | colon | ENCSR800WIY | — | SRR4422293 (PE) |
| RAMPAGE | stomach | ENCSR686VRS | ENCSR752UNJ | SRR5111956 (PE) |
| RAMPAGE | colon | ENCSR855GRG | ENCSR800WIY | SRR4422374 (PE) |
| sRNA | stomach | ENCSR003SZJ | — | SRR4422000 (SE 76bp) |
| sRNA | colon | ENCSR316SXP | — | SRR4421636 (SE 101bp) |
| WGBS | colon | ENCSR668FTZ | — | SRR8659912 (PE 151bp) |

### Donor 3 (D3): ENCDO451RUA (54M) — RNA-seq + WGBS only

| Assay | Tissue | ENCODE Exp | SRR(s) |
|-------|--------|------------|--------|
| RNA-seq | stomach | ENCSR296PMS | SRR4421597 (PE) |
| RNA-seq | colon | ENCSR630VJN | SRR4422057 (PE) |
| WGBS | stomach | ENCSR267SNS | SRR8659894 (PE 151bp) |
| WGBS | colon | ENCSR156JXJ | SRR8659891 (PE 151bp) |

### HG002 (GIAB) — dmC only

| Assay | Source | Format |
|-------|--------|--------|
| dmC | ONT Open Data S3 | Modbase-tagged BAM (MM/ML tags) |

S3 download: `aws s3 --no-sign-request cp s3://ont-open-data/giab_2025.01/...`
Chr21 extraction: `samtools view -b aligned.bam chr21 > HG002_chr21.bam`

## Known Design Considerations

### Mixed PE/SE within CTCF

D1's CTCF experiments are PE 100bp (ENCODE 2016-era Illumina HiSeq), while D2's
are SE 76bp (ENCODE 2019-era). This exercises both pipeline code paths
(filter_bam_pe and filter_bam_se) for the same mark within the same analysis.

### WGBS donor overlap

D3 (ENCDO451RUA) appears in both tissue groups for WGBS. This is a paired
design (same donor, different tissues) rather than independent replicates. For a
test case exercising the pipeline this is acceptable; the DMR results should not
be interpreted as a properly powered independent comparison.

### D2 controls shared across H3K27ac and CTCF

D2's H3K27ac and CTCF experiments in each tissue reference the same control
experiment (ENCSR539HNK for stomach, ENCSR828XQV for colon). This is standard
ENCODE practice — a single input control serves all ChIP experiments from the
same biosample. D1 has separate controls for the SE (H3K27ac) and PE (CTCF)
experiments due to different library preparation protocols.

### EN-TEx sRNA libraries are not true miRNA-seq

The EN-TEx "small RNA-seq" libraries (ENCSR573RGL, ENCSR003SZJ, ENCSR183YQY,
ENCSR316SXP) are rRNA-depleted total RNA size-selected to <200 nt — NOT true
miRNA/sRNA-seq with adapter ligation. After trimming, reads remain 76-101 bp
with zero signal in the 18-30 nt small RNA range. ShortStack finds no de novo
miRNA candidates, triggering a crash in ShortStack 4.1.2 (`zip(*denovo_mloci1)`
on an empty list).

**Workaround**: De novo miRNA discovery is disabled for this test case via
`shortstack_dn_mirna: false` in `test_options_hg38_chr21.yaml`. ShortStack
still performs cluster analysis and size profiling. The sRNA results should not
be interpreted as a proper miRNA analysis.

### sRNA read lengths differ between donors

D1 sRNA is 76bp, D2 sRNA is 101bp. Both are SE. The pipeline handles variable
read lengths; fastp will trim accordingly.

### dmC does not participate in tissue comparison

The HG002 dmC sample has `tissue:cell_line` and will not be grouped with the
stomach/colon samples for differential analysis. It exercises the dmC pipeline
path (modBAM input, modkit pileup, bedMethyl conversion) independently.

## Chr21 Subsetting Strategy

All SRA data will be downloaded as full-genome FASTQs, aligned to hg38, then
subsetted to chr21 reads. This is a one-time preparation step:

1. Download FASTQs from SRA
2. Align to full hg38 (bowtie2 for ChIP/sRNA, STAR for RNA, bismark for WGBS)
3. `samtools view -b aligned.bam chr21 > chr21_subset.bam`
4. `samtools sort -n chr21_subset.bam | samtools fastq` to regenerate FASTQs
5. For dmC: `samtools view -b modbase.bam chr21 > HG002_chr21.bam`

The chr21-subsetted FASTQs and modBAM will be the test inputs stored in
`tests/integration/data/hg38_chr21/` or cached in `test-data-prep/`.

## Sample Sheet Preview

The Sample_ID column is shown for reference only: sheets no longer carry it, and the
pipeline derives each name as `{levels_label}_{IP_target|assay_short}_{Replicate_ID}` (#91).
The PE CTCF controls take IP_target `CTCFInput` so they derive a name distinct from the
SE Input of the same donor, and form their own control merge group.

```
Sample_ID	Assay	Genome	Levels	Replicate_ID	Read_files	Read_layout	IP_target	Control
stomach_H3K27ac_D1	ChIP_broad	hg38_chr21	tissue:stomach	D1	...	SE	H3K27ac	stomach_Input_D1
stomach_H3K27ac_D2	ChIP_broad	hg38_chr21	tissue:stomach	D2	...	SE	H3K27ac	stomach_Input_D2
colon_H3K27ac_D1	ChIP_broad	hg38_chr21	tissue:colon	D1	...	SE	H3K27ac	colon_Input_D1
colon_H3K27ac_D2	ChIP_broad	hg38_chr21	tissue:colon	D2	...	SE	H3K27ac	colon_Input_D2
stomach_CTCF_D1	ChIP_narrow	hg38_chr21	tissue:stomach	D1	...	PE	CTCF	stomach_CTCFInput_D1
stomach_CTCF_D2	ChIP_narrow	hg38_chr21	tissue:stomach	D2	...	SE	CTCF	stomach_Input_D2
colon_CTCF_D1	ChIP_narrow	hg38_chr21	tissue:colon	D1	...	PE	CTCF	colon_CTCFInput_D1
colon_CTCF_D2	ChIP_narrow	hg38_chr21	tissue:colon	D2	...	SE	CTCF	colon_Input_D2
stomach_Input_D1	ChIP_broad	hg38_chr21	tissue:stomach	D1	...	SE	Input
stomach_Input_D2	ChIP_broad	hg38_chr21	tissue:stomach	D2	...	SE	Input
stomach_CTCFInput_D1	ChIP_narrow	hg38_chr21	tissue:stomach	D1	...	PE	CTCFInput
colon_Input_D1	ChIP_broad	hg38_chr21	tissue:colon	D1	...	SE	Input
colon_Input_D2	ChIP_broad	hg38_chr21	tissue:colon	D2	...	SE	Input
colon_CTCFInput_D1	ChIP_narrow	hg38_chr21	tissue:colon	D1	...	PE	CTCFInput
stomach_RNA_D1	RNAseq	hg38_chr21	tissue:stomach	D1	...	PE
stomach_RNA_D2	RNAseq	hg38_chr21	tissue:stomach	D2	...	PE
stomach_RNA_D3	RNAseq	hg38_chr21	tissue:stomach	D3	...	PE
colon_RNA_D1	RNAseq	hg38_chr21	tissue:colon	D1	...	PE
colon_RNA_D2	RNAseq	hg38_chr21	tissue:colon	D2	...	PE
colon_RNA_D3	RNAseq	hg38_chr21	tissue:colon	D3	...	PE
stomach_RAMPAGE_D1	RAMPAGE	hg38_chr21	tissue:stomach	D1	...	PE		stomach_RNA_D1
stomach_RAMPAGE_D2	RAMPAGE	hg38_chr21	tissue:stomach	D2	...	PE		stomach_RNA_D2
colon_RAMPAGE_D1	RAMPAGE	hg38_chr21	tissue:colon	D1	...	PE		colon_RNA_D1
colon_RAMPAGE_D2	RAMPAGE	hg38_chr21	tissue:colon	D2	...	PE		colon_RNA_D2
stomach_sRNA_D1	sRNA	hg38_chr21	tissue:stomach	D1	...	SE
stomach_sRNA_D2	sRNA	hg38_chr21	tissue:stomach	D2	...	SE
colon_sRNA_D1	sRNA	hg38_chr21	tissue:colon	D1	...	SE
colon_sRNA_D2	sRNA	hg38_chr21	tissue:colon	D2	...	SE
stomach_WGBS_D1	WGBS	hg38_chr21	tissue:stomach	D1	...	PE
stomach_WGBS_D3	WGBS	hg38_chr21	tissue:stomach	D3	...	PE
colon_WGBS_D2	WGBS	hg38_chr21	tissue:colon	D2	...	PE
colon_WGBS_D3	WGBS	hg38_chr21	tissue:colon	D3	...	PE
cell_line_dmC_HG002	dmC	hg38_chr21	tissue:cell_line	HG002	...	SE
```

33 samples total.
