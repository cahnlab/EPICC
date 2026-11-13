#!/usr/bin/env Rscript

library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(txdbmaker)

args = commandArgs(trailingOnly=TRUE)
# 1 track table, 2 genes GFF, 3 TEs BED, 4 title, 5 output PDF,
# 6 highlight starts, 7 highlight widths (both may be empty files),
# 8 extra-track table (optional; Label/Color/Path, one row per BED track)

nonempty_path<-function(p) { file.exists(p) && file.info(p)$size > 0 }
nonempty<-function(i) { length(args) >= i && nonempty_path(args[i]) }

filenames<-read.delim(args[1], header=TRUE)
if ( nonempty(2) ) {
	genes<-txdbmaker::makeTxDbFromGFF(args[2], format="gff")
} else {
	genes<-c()
}

tes<-import(args[3], format="bed")

plotname<-args[4]
pdfname<-args[5]

tot<-nrow(filenames)

# The first highlighted region is blue and every further one red, however
# many the target file lists (issue #98: the old fixed six-colour vector left
# every region past the sixth without a colour).
htcol<-c()
htcol2<-c()
if ( nonempty(6) ) {
	htstart<-read.delim(args[6], header=FALSE)[[1]]
	htwidth<-read.delim(args[7], header=FALSE)[[1]]
	if ( length(htstart) != length(htwidth) ) {
		stop("Highlight starts (", length(htstart), ") and widths (", length(htwidth), ") differ in number")
	}
	first<-seq_along(htstart) == 1
	htcol<-ifelse(first, "#B7E2FD", "#fac0c7")
	htcol2<-ifelse(first, "#F6FBFE", "#fffafa")
}

# User BED files (browser_extra_tracks), each drawn as one annotation row
# below the TEs.
extratracks<-list()
if ( nonempty(8) ) {
	extras<-read.delim(args[8], header=TRUE, colClasses="character", quote="", comment.char="")
	for ( i in seq_len(nrow(extras)) ) {
		features<-if ( nonempty_path(extras$Path[i]) ) import(extras$Path[i], format="bed") else GRanges()
		extratracks[[i]]<-AnnotationTrack(features, name=extras$Label[i], stacking="dense", fill=extras$Color[i], col=extras$Color[i], shape="box", rotation.title=0, cex.title=0.5, lwd=0.1, fontsize=12, fontcolor.title="black")
	}
}

options(ucscChromosomeNames=FALSE)

tracksize<-c(1,1,0.5,rep(0.5,length(extratracks)))
tracklist<-list()
for ( i in c(1:tot) ) {
	label<-filenames$Name[i]
	path<-filenames$Path[i]
	backcolor<-filenames$Backcolor[i]
	trackcolor<-filenames$Trackcolor[i]
	fillcolorplus<-filenames$Fillcolorplus[i]
	fillcolorminus<-filenames$Fillcolorminus[i]
	ymin<-filenames$Ymin[i]
	ymax<-filenames$Ymax[i]
	ymintick<-sign(ymin)*((floor(abs(ymin)*100)/100))
	ymaxtick<-sign(ymax)*((floor(abs(ymax)*100)/100))
	tracksize<-c(tracksize,1)
	print(paste0("Importing bw for ",label))
	bw<-import(path)
	print(paste0("Creating track for ",label))
	lab<-gsub("_", "\n", label)
	track<-DataTrack(bw, type="polygon", baseline=0, name=lab, background.title = backcolor, col=trackcolor, fill.mountain=c(fillcolorminus,fillcolorplus), col.baseline="grey50", ylim=c(ymin,ymax), yTicksAt=c(ymintick,ymaxtick), cex.title=0.5, lwd=0.01, fontsize=12, fontcolor.title="black", col.axis="black", cex.axis=0.5)
	tracklist<-append(tracklist, track)
}

axistrack<-GenomeAxisTrack(scale=0.1, labelPos="above")
tetrack<-AnnotationTrack(tes, name="TEs", stacking = "dense", fill = "lightgreen", shape="box", rotation.title=0, cex.title=0.5, lwd=0.1, fontsize=12, fontcolor.title="black")

# Gene track layouts, most detailed first. Gviz refuses to draw a squished
# track whose rows don't fit the space it gets ("Too many stacks to draw"),
# which happens on gene-dense or wide regions (issue #98). Each fallback trades
# detail for rows: "meta" collapses every gene's transcripts onto one model so
# overlapping isoforms stop adding rows, and is labelled with its gene ID
# rather than one arbitrary transcript ID; the last layout
# draws all genes on a single row without names, which always fits.
gene_layouts<-list(
	full=list(shape="smallArrow", collapseTranscripts=FALSE, showId=TRUE, stacking="squish"),
	collapsed=list(shape="smallArrow", collapseTranscripts="meta", showId=TRUE, transcriptAnnotation="gene", stacking="squish"),
	dense=list(shape="box", collapseTranscripts=TRUE, showId=FALSE, stacking="dense")
)

make_gene_track<-function(layout) {
	do.call(GeneRegionTrack, c(list(genes, name="Genes", col="black", fill="grey60", rotation.title=0, cex.title=0.5, lwd=0.1, fontsize=12, fontcolor.title="black"), layout))
}

draw_browser<-function(genetrack) {
	if ( length(htcol) > 0 ) {
		httrack <- HighlightTrack(trackList=tracklist, start=htstart, width=htwidth, col=htcol, fill=htcol2)
		tracks<-c(list(axistrack, genetrack, tetrack), extratracks, list(httrack))
	} else {
		tracks<-c(list(axistrack, genetrack, tetrack), extratracks, tracklist)
	}
	# Extra annotation rows get half an inch each so they don't squeeze the
	# data tracks.
	pdf(pdfname, width = 12, height = tot + 0.5*length(extratracks))
	# Close the device whether or not plotTracks succeeds, so a failed attempt
	# leaves no half-written page behind for the next layout to append to.
	on.exit(dev.off())
	plotTracks(tracks, sizes=tracksize, main=plotname, cex.main = 1)
}

for ( layout_name in names(gene_layouts) ) {
	drawn<-tryCatch({
		draw_browser(make_gene_track(gene_layouts[[layout_name]]))
		TRUE
	}, error=function(e) {
		if ( !grepl("Too many stacks", conditionMessage(e)) ) stop(e)
		print(paste0("Gene track does not fit with the '", layout_name, "' layout, trying a more compact one"))
		FALSE
	})
	if ( drawn ) {
		print(paste0("Gene track drawn with the '", layout_name, "' layout"))
		break
	}
}
if ( !drawn ) stop("Gene track does not fit even with the dense layout")





