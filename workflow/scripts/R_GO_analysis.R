#!/usr/bin/env Rscript

library(edgeR)
library(AnnotationForge)
library(rrvgo)
library(dplyr)
library(topGO)
library(purrr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

dbname<-args[1]
analysisname<-args[2]
refgenome<-args[3]
targetfile<-args[4]
backgroundfile<-args[5]
targetname<-args[6]
output_dir<-args[7]
genome_dir<-args[8]

godir <- normalizePath(paste0("./",genome_dir,"/",refgenome,"/GO"), mustWork=FALSE)
# Prepend the per-genome GO directory to .libPaths() so requireNamespace,
# library, and rrvgo's orgdb-driven calls below all resolve dbname against
# this genome's installed package -- not whichever same-named package
# happens to be in the conda env's R library from an earlier run with
# a different reference genome of the same species.
.libPaths(c(godir, .libPaths()))
if (!requireNamespace(dbname, quietly = TRUE)) {
	stop(paste0("GO database package '", dbname, "' not found in ", godir,
	            ". Run epicc with GO: true so the create_GO_database rule builds it first."))
}
library(dbname, character.only = TRUE)
info<-read.delim(file.path(godir, paste0(dbname,"_gaf_file.tab")), header=FALSE)

fGO<-info[,c(1,6,10)]
colnames(fGO)<-c("GID","GO","EVIDENCE")
geneid2GO<-unique(fGO[,c(1,2)])
rn1<-paste(geneid2GO[,1], sep="")
gene2GO<-geneid2GO[,-1]
names(gene2GO)<-rn1

# Every GO output carries the analysis name and reference genome, like the
# other analysis-level RNA outputs (plot_expression__<analysis>__<genome>__...).
# Without them, two analyses sharing an output_dir and a target label (e.g. the
# same level name "WT" giving UP_in_WT) overwrite each other's tables (issue #101).
go_prefix<-function(name, ont) {
	paste0("topGO__",analysisname,"__",refgenome,"__",name,"__",ont)
}

getGO<-function(genelist, target, ont, name) {
	GOdata<-new("topGOdata", 
				ontology = ont, 
				allGenes = genelist,
				annot = annFUN.gene2GO, 
				gene2GO = gene2GO)
	resultFisher<-runTest(GOdata, algorithm = "weight01", statistic = "fisher")
	resultFisherSummary <- summary(attributes(resultFisher)$score <= 0.01)
	nSigTerms<-0
	if (length(resultFisherSummary) == 3) {
		nSigTerms<-as.integer(resultFisherSummary[[3]])
	}
	summary<-GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nSigTerms, numChar=1000)
	tab<-summary %>%
		rename_with(.cols = starts_with("apply"), .fn = ~ { if (length(.) > 0) { paste0("classicFisher") } else { . } }) %>%
		mutate(classicFisher = classicFisher %>% str_replace(pattern= "< *1e-30", replacement = "1e-30") %>% as.numeric())
	sigTerms<-tab$GO.ID
	genesInTerms<-genesInTerm(GOdata, sigTerms)
	genesInTerms2<-map(genesInTerms, ~ intersect(.x, myInterestedGenes) %>% paste(collapse = ","))
	tab2<-tab %>% 
		left_join(tibble(GO.ID = names(genesInTerms2), GID = genesInTerms2) %>% 
		tidyr::unnest(GID), by = "GO.ID")
	tab3<-tab %>%
		rename(GO=GO.ID) %>%
		merge(geneid2GO, by="GO") %>%
		merge(target, by="GID") %>%
		arrange(GO) %>%
		unique()
	if (nrow(tab2) > 1) {
		write.table(tab2,paste0(output_dir,"/RNA/GO/",go_prefix(name, ont),"_GOs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}
	if (nrow(tab3) > 0) {
		write.table(tab3,paste0(output_dir,"/RNA/GO/",go_prefix(name, ont),"_GIDs.txt"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
	}	
  
	scores<-setNames(-log10(as.numeric(tab$classicFisher)), tab$GO.ID)
	reducedTerms <- tab2 %>%
			rename("go" = "GO.ID", "term" = "Term") %>%
			mutate(parentTerm = term, score = scores)
  
	if (nrow(tab) > 1 ) {
		simMatrix<-calculateSimMatrix(tab$GO.ID,
									orgdb=dbname,
									ont=ont,
									method="Rel")
		if (! is.null(nrow(simMatrix))) {
			reducedTerms<-reduceSimMatrix(simMatrix,
										scores,
										threshold = 0.7,
										orgdb=dbname)
										
			pdf(paste0(output_dir,"/RNA/plots/",go_prefix(name, ont),"_treemap.pdf"), width=8, height=8)
			treemapPlot(reducedTerms, size = "score")
			dev.off()
		}
	}
}

if (startsWith(backgroundfile, paste0(output_dir,"/RNA/DEG/counts__"))) {
	genecount<-read.delim(backgroundfile, header = TRUE, row.names = "GID")
	target<-read.delim(targetfile, header = TRUE)

	genecount<-genecount[!grepl("^N_", rownames(genecount)), ]
	keep.exprs<-rowSums(cpm(genecount)>1)>=2
	filtered<-genecount[keep.exprs,]
	filtered$GID<-row.names(filtered)

	allGenes<-unique(unlist(filtered$GID))
	
	for ( samp in unique(target$Sample) ) {
		for ( deg in c("UP","DOWN") ) {
			samplename<-paste0(deg,"_in_",samp)
			sampletable<-filter(target, Sample==samp, DEG==deg)
			myInterestedGenes<-unique(unlist(sampletable$GID))
			geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
			names(geneList)<-allGenes
			
			if (length(levels(geneList)) == 2) {			
				for ( ont in c("BP","MF") ) {
					print(paste0("Getting ",ont," for ",samplename))
					getGO(geneList, sampletable, ont, samplename)
				}
			} else {
				print("Target genes empty or not in the background list")
			}
		}
	}

} else if (startsWith(backgroundfile, paste0(output_dir,"/combined/bedfiles/"))) {

	ref_genes<-read.delim(backgroundfile, header = FALSE, col.names = c("Chr","Start","Stop","Name","Value","Strand"))
	ref_genes<-mutate(ref_genes, GID=str_replace(ref_genes$Name, pattern = ".*ID=(gene:)?([^;]+).*", replacement = "\\2")) %>%
				select(-Name, -Value)
	ref_genes$GID<-str_remove_all(ref_genes$GID, pattern = "_.$")
	allGenes<-unique(unlist(ref_genes$GID))
	
	target<-read.delim(targetfile, header = TRUE)
	myInterestedGenes<-unique(unlist(target$GID))
	geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
	names(geneList)<-allGenes
	
	if (length(levels(geneList)) == 2) {
		for ( ont in c("BP","MF") ) {
			print(paste0("Getting ",ont," for ",targetname))
			getGO(geneList, target, ont, targetname)
		} 
	} else {
		print("Target genes empty or not in the background list")
	}

	
} else {

	background<-read.delim(backgroundfile, header = TRUE)
	allGenes<-unique(unlist(background$GID))
	
	target<-read.delim(targetfile, header = TRUE)
	myInterestedGenes<-unique(unlist(target$GID))
	geneList<-factor(as.integer(allGenes %in% myInterestedGenes))
	names(geneList)<-allGenes
	
	if (length(levels(geneList)) == 2) {
		for ( ont in c("BP","MF") ) {
			print(paste0("Getting ",ont," for ",targetname))
			getGO(geneList, target, ont, targetname)
		}
	} else {
		print("Target genes empty or not in the background list")
	}
}
