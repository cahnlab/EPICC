#!/usr/bin/env Rscript

library(limma)
library(edgeR)
library(dplyr)
library(tidyr)
library(stringr)
library(gplots)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)

genecount<-read.delim(args[1], header = TRUE, row.names = "GID", check.names = FALSE)
genecount<-genecount[!grepl("^N_", rownames(genecount)), ]
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

ref_genes<-read.delim(args[5], header = FALSE, 
                      col.names = c("Chr","Start","Stop","Name","Value","Strand"))
ref_genes<-mutate(ref_genes, GID=str_replace(ref_genes$Name, pattern = ".*ID=(gene:)?([^;]+).*", replacement = "\\2")) %>%
  select(-Name, -Value)
ref_genes$GID<-str_remove_all(ref_genes$GID, pattern = "_.$")

output_dir<-args[6]

# placeholder (default values) for user-defined parameters to filter DEG
xint<-0.5 #min absolute log2FC
yint<-0.05 #min FDR

# EdgeR analysis

y<-DGEList(counts=filtered, group = samples)
y<-calcNormFactors(y)

pdf(paste0(output_dir,"/combined/plots/MDS_RNAseq_",analysisname,"_",refgenome,"_d12.pdf"),10,8)
plotMDS(y, col=color_samples, pch=16)
dev.off()

pdf(paste0(output_dir,"/combined/plots/MDS_RNAseq_",analysisname,"_",refgenome,"_d12_labs.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps)
dev.off()

pdf(paste0(output_dir,"/combined/plots/MDS_RNAseq_",analysisname,"_",refgenome,"_d23.pdf"),10,8)
plotMDS(y, col=color_samples, pch=16, dim.plot=c(2,3))
dev.off()

pdf(paste0(output_dir,"/combined/plots/MDS_RNAseq_",analysisname,"_",refgenome,"_d23_labs.pdf"),10,8)
plotMDS(y, col=color_samples, labels=reps, dim.plot=c(2,3))
dev.off()

y<-estimateCommonDisp(y, verbose = TRUE)
y<-estimateTagwiseDisp(y)

pdf(paste0(output_dir,"/combined/plots/BCV_RNAseq_",analysisname,"_",refgenome,".pdf"),10,8)
plotBCV(y)
dev.off()

#### Function to create a fold change table for all genes between a pair of samples

create.FC.table<-function(sample1, sample2, table) {
  et<-exactTest(table, pair = c(sample2,sample1))
  out<-topTags(et, n=Inf, adjust.method="BH")
  table<-as.data.frame(out)
  table$GID<-rownames(table)
  table<-mutate(table, Sample=paste0(sample1,"_vs_",sample2))
  table
}

#### Function to create a table of DEGs between a pair of samples

create.DEG.table<-function(sample1, sample2, y) {
  et<-exactTest(y, pair = c(sample2,sample1))
  out<-topTags(et, n=Inf, adjust.method="BH")
  table<-filter(out$table, FDR<=yint, logFC<=-xint | logFC>=xint)
  table$GID<-rownames(table)
  table<-mutate(table, DEG=as.factor(ifelse(logFC<0, "DOWN", "UP")),Sample=paste0(sample1,"_vs_",sample2))
  table
}

#### Function to create a volcano plot to visualize number and distribution of DEGs between a pair of samples

plotvolcano<-function(sample1,sample2,table,xint,yint) {
  volcano<-ggplot(table, aes(logFC,-log2(FDR))) + 
    geom_point(aes(color=DEG, alpha=abs(logFC), size=-log2(FDR))) + 
    scale_color_manual(values = c(No="grey", UP="red", DOWN="blue")) +
    scale_alpha_continuous(range = c(0.1, 1)) +
    scale_size_continuous(range = c(0.3, 3)) +
    geom_vline(xintercept=-xint, linewidth=0.5, color="grey", linetype = 2) +
    geom_vline(xintercept=xint, linewidth=0.5, color="grey", linetype = 2) +
    geom_hline(yintercept=-log2(yint), linewidth=0.5, color="grey", linetype = 2) +
    labs(x=paste0("log2FC(",sample1," vs ",sample2,")")) +
    theme_bw() +
    theme(legend.position = "none")
  print(volcano)
}

#### Applying above function to all possible pairs of samples, and creating temp files for unique DEGs

allDEG<-data.frame()
for (i in seq_len(length(genotypes)-1)) {
  sample1<-genotypes[i]
  for (j in (i+1):length(genotypes)) {
	sample2<-genotypes[j]
	FCtable<-create.FC.table(sample1,sample2,y)
	FCtable<-merge(ref_genes,FCtable,by=c("GID")) %>%
		select(Chr,Start,Stop,GID,logFC,Strand,logCPM,PValue,FDR,Sample) %>%
		arrange(Chr,Start)
	write.table(FCtable,paste0(output_dir,"/RNA/DEG/FC_",analysisname,"__",refgenome,"__",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	tmptable<-mutate(FCtable, DEG=ifelse(FDR<=yint & logFC<=-xint, "DOWN", ifelse(FDR<=yint & logFC>=xint, "UP", "No")))
    tmptable$DEG<-as.factor(tmptable$DEG)
    pdf(paste0(output_dir,"/RNA/plots/Volcano__",analysisname,"__",refgenome,"__",sample1,"_vs_",sample2,".pdf"),6,6)
	plotvolcano(sample1,sample2,table2,xint,yint)
    dev.off()
	DEGtable<-create.DEG.table(sample1,sample2,y)
	DEGtable<-merge(ref_genes,DEGtable,by=c("GID")) %>%
		select(Chr,Start,Stop,GID,logFC,Strand,logCPM,PValue,FDR,Sample,DEG) %>%
		arrange(DEG,Chr,Start)
	write.table(DEGtable,paste0(output_dir,"/RNA/DEG/DEG_",analysisname,"__",refgenome,"__",sample1,"_vs_",sample2,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	temptable<-mutate(DEGtable, firstsample = sample1, secondsample = sample2) %>%
				select(GID, DEG, firstsample, secondsample)
	allDEG<-rbind(allDEG,temptable)
  }
}

#### To create a table of DEGs unique to each sample

uniqueUP<-data.frame()
uniqueDOWN<-data.frame()
for (sample1 in genotypes) {
	tempUP<-filter(allDEG, (DEG=="UP" & firstsample==sample1) | (DEG=="DOWN" & secondsample==sample1)) %>%
			mutate(Sample=sample1) %>%
			select(GID, Sample)
	tempUP
	uniqueUP<-rbind(uniqueUP, tempUP)
	tempDOWN<-filter(allDEG, (DEG=="UP" & secondsample==sample1) | (DEG=="DOWN" & firstsample==sample1)) %>%
			mutate(Sample=sample1) %>%
			select(GID, Sample)
	uniqueDOWN<-rbind(uniqueDOWN, tempDOWN)
	tempDOWN
}

uniqueUP<-unique(uniqueUP) %>%
		group_by(GID) %>%
		filter(n() == 1) %>%
		ungroup() %>%
		mutate(DEG="UP")

uniqueDOWN<-unique(uniqueDOWN) %>%
		group_by(GID) %>%
		filter(n() == 1) %>%
		ungroup() %>%
		mutate(DEG="DOWN")

uniqueDEGs<-rbind(uniqueUP, uniqueDOWN) %>%
			arrange(Sample)

write.table(uniqueDEGs,paste0(output_dir,"/RNA/DEG/unique_DEGs__",analysisname,"__",refgenome,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

#### To create a summary table of number of DEGs

stat_table<-group_by(allDEG, firstsample, secondsample, DEG) %>%
			summarize(Nb=n())
	
for (sample1 in genotypes) {
	nunique<-filter(uniqueDEGs, Sample == sample1) %>%
		mutate(secondsample="Unique") %>%
		select(firstsample=Sample, secondsample, DEG) %>%
		group_by(firstsample, secondsample, DEG) %>%
		summarize(Nb=n())
	
	stat_table<-rbind(stat_table, nunique)
}

write.table(stat_table,paste0(output_dir,"/RNA/DEG/summary_DEG_stats__",analysisname,"__",refgenome,".txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

#### To create heatmaps over all DEGs (by count per million and z-score)

keepDEG<-unique(allDEG$GID)

if (length(keepDEG) >= 2) {
	logcounts<-cpm(y, log=TRUE)
	lcpm<-logcounts[keepDEG,]

	pdf(paste0(output_dir,"/combined/plots/Heatmap_RNAseq_cpm__",analysisname,"__",refgenome,".pdf"),10,15)
	heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
			main=paste0("Differentially expressed genes in ",refgenome," from ",analysisname),
			margins=c(12,2),cexCol=1, labRow = "", col="bluered", srtCol=45,
			lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "", key.xlab = "log(cpm)")
	dev.off()

	pdf(paste0(output_dir,"/combined/plots/Heatmap_RNAseq_zscore__",analysisname,"__",refgenome,".pdf"),10,15)
	heatmap.2(lcpm,trace="none",ColSideColors = color_samples,
			main=paste0("Differentially expressed genes in ",refgenome," from ",analysisname),
			margins=c(12,2),cexCol=1, labRow = "", col="bluered", srtCol=45, scale="row",
			lwid=c(1,5),lhei=c(0.5,5,0.1), key.title = "")
	dev.off()
}

### To make R object for later plotting gene expression

norm<-cpm(y, normalized.lib.size=T)
genextable<-data.frame(norm, stringsAsFactors = FALSE, check.names = FALSE)
genextable<-mutate(genextable, GID=row.names(genextable))

plot.Expression <- function(gene, label) {
  
	dataline<-filter(genextable, GID==gene) %>% 
		pivot_longer(cols = -GID, names_to="Replicate", values_to="CountPerMillion")
  
	dataline<-merge(dataline, targets, by=c("Replicate")) %>%
			merge(uniqueDEGs, by=c("GID","Sample"), all.x = TRUE)
  
	dataline$Sample<-as.factor(dataline$Sample)
	dataline$CountPerMillion<-as.numeric(dataline$CountPerMillion)
	dataline<-group_by(dataline, Sample) %>% 
			mutate(Average = mean(CountPerMillion))
	dataline$Average[is.na(dataline$Average)]<-as.numeric(0)
	dataline$DEG<-as.factor(dataline$DEG)
	if ( label == "NoLabel" ) {
		plottitle<-paste0(gene)
	} else { 
		plottitle<-paste0(label," (",gene,")")
	}
  
	plot<-ggplot(dataline, aes(Sample,DEG)) + 
			geom_col(position="dodge", aes(y=Average, fill=DEG)) + 
			scale_fill_manual(values = c("0"="grey", "UP"="pink", "DOWN"="lightblue"),
							labels=c("0"="No", "UP"="Up", "DOWN"="Down")) +
			geom_point(aes(y=CountPerMillion), size=2, shape=3) + 
			labs(title = plottitle, y="cpm") + 
			theme(axis.title.y=element_text(size=10), axis.title.x=element_blank(),
				plot.title=element_text(size=15), 
				axis.text.x=element_text(size=10, angle = 90),
				panel.grid.major.y = element_blank(), 
				panel.grid.minor.y = element_blank(),
				panel.grid.major.x = element_blank(),
				panel.background = element_rect(fill = "white", colour = "black"))
	plot  
}

save(plot.Expression,genextable,targets,uniqueDEGs, file = paste0(output_dir,"/RNA/DEG/ReadyToPlot__",analysisname,"__",refgenome,".RData"))

