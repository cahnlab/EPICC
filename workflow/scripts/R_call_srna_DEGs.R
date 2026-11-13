#!/usr/bin/env Rscript

library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)
library(gplots)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

genecount<-read.delim(args[1], header = TRUE, row.names = "Name", check.names = FALSE)
keep.exprs<-rowSums(cpm(genecount)>1)>=2
filtered<-genecount[keep.exprs,]

targets<-read.delim(args[2], header = TRUE, check.names = FALSE)
samples<-as.factor(targets$Sample)
reps<-as.factor(targets$Replicate)
genotypes<-unique(samples)

targets$Color<-as.numeric(targets$Color)

qual_col_pals<-brewer.pal.info[brewer.pal.info$category == 'qual',]
colors<-unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #limited to 74 different samples. Could change this to include more if needed.
color_samples<-colors[targets$Color]

analysisname<-args[3]
refgenome<-args[4]
targetname<-args[5]
filename<-args[6]
output_dir<-args[7]
genome_dir<-args[8]

# The prepared TE file is <G>__TE_file.bed, or <G>__TE_file_<tag>.bed for an
# alternative annotation (epicc --te-file); both are headerless BED6.
is_te_file <- dirname(filename) == paste0(genome_dir,"/",refgenome) &&
	startsWith(basename(filename), paste0(refgenome,"__TE_file")) &&
	grepl("^__TE_file(_[A-Za-z0-9_-]+)?\\.bed$", substring(basename(filename), nchar(refgenome) + 1))
if ( is_te_file || filename %in% c(paste0(output_dir,"/combined/bedfiles/",refgenome,"__all_genes.bed"), paste0(output_dir,"/combined/bedfiles/",refgenome,"__protein_coding_genes.bed")) ) {
	region_file<-read.delim(filename, header = FALSE, col.names = c("Chr","Start","Stop","Name","Value","Strand"))
} else {
	region_file<-read.delim(filename, header = TRUE)
}

# EdgeR analysis

y<-DGEList(counts=filtered, group = samples)
y<-calcNormFactors(y)

pdf(paste0(output_dir,"/combined/plots/MDS_sRNA_",analysisname,"_",refgenome,"__on_",targetname,"_d12.pdf"),10,8)
plotMDS(y, col=color_samples, pch=16)
dev.off()

pdf(paste0(output_dir,"/combined/plots/MDS_sRNA_",analysisname,"_",refgenome,"__on_",targetname,"_d12_labs.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps)
dev.off()

pdf(paste0(output_dir,"/combined/plots/MDS_sRNA_",analysisname,"_",refgenome,"__on_",targetname,"_d23.pdf"),10,8)
plotMDS(y, col=color_samples, pch=16, dim.plot=c(2,3))
dev.off()

pdf(paste0(output_dir,"/combined/plots/MDS_sRNA_",analysisname,"_",refgenome,"__on_",targetname,"_d23_labs.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps, dim.plot=c(2,3))
dev.off()

y<-estimateCommonDisp(y, verbose = TRUE)
y<-estimateTagwiseDisp(y)

pdf(paste0(output_dir,"/combined/plots/BCV_sRNA_",analysisname,"_",refgenome,"__on_",targetname,".pdf"),10,8)
plotBCV(y)
dev.off()

#### Function to create a fold change table for all genes between a pair of samples

create.FC.table<-function(sample1, sample2, table) {
  et<-exactTest(table, pair = c(sample2,sample1))
  out<-topTags(et, n=Inf, adjust.method="BH")
  table<-as.data.frame(out)
  table$Name<-rownames(table)
  table<-mutate(table, Sample=paste0(sample1,"_vs_",sample2))
  table
}

#### Function to create a table of DEGs between a pair of samples

create.DEG.table<-function(sample1, sample2, y) {
  et<-exactTest(y, pair = c(sample2,sample1))
  out<-topTags(et, n=Inf, adjust.method="BH")
  table<-filter(out$table, FDR<=0.05, logFC<=-2 | logFC>=2)
  table$Name<-rownames(table)
  table<-mutate(table, DEG=as.factor(ifelse(logFC<0, "DOWN", "UP")),Sample=paste0(sample1,"_vs_",sample2))
  table
}

#### Applying above function to all possible pairs of samples, and creating temp files for unique DEGs

allDEG<-data.frame()
for (i in seq_len(length(genotypes)-1)) {
  sample1<-genotypes[i]
  for (j in (i+1):length(genotypes)) {
	sample2<-genotypes[j]
	FCtable<-create.FC.table(sample1,sample2,y)
	if ( nrow(FCtable) > 0 ) {
		FCtable<-merge(region_file,FCtable,by=c("Name"))
		write.table(FCtable,paste0(output_dir,"/sRNA/clusters/",analysisname,"__",refgenome,"__on_",targetname,"/FC_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	DEGtable<-create.DEG.table(sample1,sample2,y)
	if ( nrow(DEGtable) > 0 ) {
		DEGtable<-merge(region_file,DEGtable,by=c("Name"))
		write.table(DEGtable,paste0(output_dir,"/sRNA/clusters/",analysisname,"__",refgenome,"__on_",targetname,"/DEG_",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
		temptable<-mutate(DEGtable, firstsample = sample1, secondsample = sample2) %>%
					select(Name, DEG, firstsample, secondsample)
		allDEG<-rbind(allDEG,temptable)
	}
  }
}

#### To create a table of DEGs unique to each sample

uniqueUP<-data.frame()
uniqueDOWN<-data.frame()
if ( nrow(allDEG) > 0 ) {
	for (sample1 in genotypes) {
		tempUP<-filter(allDEG, (DEG=="UP" & firstsample==sample1) | (DEG=="DOWN" & secondsample==sample1)) %>%
				mutate(Sample=sample1) %>%
				select(Name, Sample)
		if ( nrow(tempUP) > 0 ) {
			uniqueUP<-rbind(uniqueUP, tempUP)
		}
		tempDOWN<-filter(allDEG, (DEG=="UP" & secondsample==sample1) | (DEG=="DOWN" & firstsample==sample1)) %>%
				mutate(Sample=sample1) %>%
				select(Name, Sample)
		if ( nrow(tempDOWN) > 0 ) {
			uniqueDOWN<-rbind(uniqueDOWN, tempDOWN)
		}
	}
}

if ( nrow(uniqueUP) > 0 ) {
	uniqueUP<-unique(uniqueUP) %>%
			group_by(Name) %>%
			filter(n() == 1) %>%
			ungroup() %>%
			mutate(DEG="UP")
}

if ( nrow(uniqueDOWN) > 0 ) {
	uniqueDOWN<-unique(uniqueDOWN) %>%
			group_by(Name) %>%
			filter(n() == 1) %>%
			ungroup() %>%
			mutate(DEG="DOWN")
}

uniqueDEGs<-data.frame()
if ( nrow(uniqueUP) > 0 && nrow(uniqueDOWN) > 0 ) {
	uniqueDEGs<-rbind(uniqueUP, uniqueDOWN) %>%
				arrange(Sample)
	write.table(uniqueDEGs,paste0(output_dir,"/sRNA/clusters/",analysisname,"__",refgenome,"__on_",targetname,"/unique_DEGs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
} else if ( nrow(uniqueUP) > 0 ) {
	uniqueDEGs<-uniqueUP
	write.table(uniqueDEGs,paste0(output_dir,"/sRNA/clusters/",analysisname,"__",refgenome,"__on_",targetname,"/unique_DEGs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
} else if ( nrow(uniqueDOWN) > 0 ) {
	uniqueDEGs<-uniqueDOWN
	write.table(uniqueDEGs,paste0(output_dir,"/sRNA/clusters/",analysisname,"__",refgenome,"__on_",targetname,"/unique_DEGs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

#### To create a summary table of number of DEGs

if ( nrow(allDEG) > 0 ) {
	stat_table<-group_by(allDEG, firstsample, secondsample, DEG) %>%
				summarize(Nb=n())
	if ( nrow(uniqueDEGs) > 0 ) {
		for (sample1 in genotypes) {
			nunique<-filter(uniqueDEGs, Sample == sample1) %>%
				mutate(secondsample="Unique") %>%
				select(firstsample=Sample, secondsample, DEG) %>%
				group_by(firstsample, secondsample, DEG) %>%
				summarize(Nb=n())
			if ( nrow(nunique) > 0 ) {
				stat_table<-rbind(stat_table, nunique)
			}
		}
	}
	write.table(stat_table,paste0(output_dir,"/sRNA/reports/summary_DEG_stats__",analysisname,"__",refgenome,"__on_",targetname,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
}

#### To create heatmaps over all DEGs (by count per million and z-score)

if ( nrow(allDEG) > 0 ) {
	keepDEG<-unique(allDEG$Name)

	if (length(keepDEG) >= 2) {
		logcounts<-cpm(y, log=TRUE)
		lcpm<-logcounts[keepDEG,]

		pdf(paste0(output_dir,"/combined/plots/Heatmap_sRNA_cpm__",analysisname,"__",refgenome,"__on_",targetname,".pdf"),10,15)
		heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
				main=paste0("Differential sRNA in ",refgenome," from ",analysisname," mapping to ",targetname),
				margins=c(12,2),cexCol=1, labRow = "", col="bluered", srtCol=45,
				lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "", key.xlab = "log(cpm)")
		dev.off()

		pdf(paste0(output_dir,"/combined/plots/Heatmap_sRNA_zscore__",analysisname,"__",refgenome,"__on_",targetname,".pdf"),10,15)
		heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
				main=paste0("Differential sRNA in ",refgenome," from ",analysisname," mapping to ",targetname),
				margins=c(12,2),cexCol=1, labRow = "", col="bluered", srtCol=45, scale="row",
				lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "")
		dev.off()
	}
}
