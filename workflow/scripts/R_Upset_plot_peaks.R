#!/usr/bin/env Rscript


library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ComplexUpset)
library(RColorBrewer)

# Resolved from this script's own path so it works from the repo and an install alike
script_dir<-dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])))
source(file.path(script_dir, "upset_common.R"))

args = commandArgs(trailingOnly=TRUE)

merged<-read.delim(args[1], header = TRUE)
annotated<-read.delim(args[2], header = TRUE) %>%
	mutate(distance=abs(Distance)+1) %>%
	select(PeakID=RegionID,Category,Gap=distance) %>%
	rename(Distance=Gap)
env<-args[3]
types<-unlist(strsplit(args[4], ":"))
label_types<-parse_label_types(args[5])
output<-args[6]

sampleslist<-unique(unlist(strsplit(as.character(merged$Samples), ",")))

# ComplexUpset needs at least two sets (samples with peaks) to build an
# intersection matrix; with 0 or 1 non-empty peak file it aborts with
# "Needs at least two indicator variables". Reachable as a corner case when
# all but one of the merged peak files are empty (no peaks called, or an
# upstream failure). Emit an informative placeholder PDF instead of failing
# the whole combined-analysis DAG.
if (length(sampleslist) < 2) {
	msg <- if (length(sampleslist) < 1) {
		"No sample has any peaks - UpSet plot skipped"
	} else {
		paste0("Only one sample has peaks (", sampleslist,
		       ") - UpSet plot needs at least two")
	}
	pdf(output, height = 10, width = 10)
	plot.new()
	text(0.5, 0.5, msg, cex = 1.3, col = "grey40")
	dev.off()
	quit(save = "no", status = 0)
}

figsize<-length(sampleslist)

mat<-separate_rows(merged, Samples, sep = ",") %>%
	mutate(value=1) %>%
	pivot_wider(names_from = Samples, values_from = value, values_fill = 0) %>%
	merge(annotated, by="PeakID")
	
mat$Category<-factor(mat$Category, levels=c("Distal_downstream","Terminator","Gene_body","Promoter","Distal_upstream"))

type_cols<-upset_type_cols(types, label_types, sampleslist)

## To create queries to color when the same mark is shared in the intersection matrix
qual_col_pals<-brewer.pal.info[brewer.pal.info$category == 'qual',]
colorlist<-unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
i<-1
queries<-c()
listcolor<-c()
for (sampletype in types) {
	setcols<-type_cols[[sampletype]]
	combos<-map(seq_along(setcols), ~ combn(setcols, ., FUN = c, simplify = FALSE)) %>% 
			unlist(recursive = FALSE)
	tmpqueries<-map(combos, ~ upset_query(intersect = .x, color = colorlist[i], 
                               fill = colorlist[i], only = TRUE, only_components = c('intersections_matrix')))
	queries<-append(queries, tmpqueries)
	listcolor<-append(listcolor, colorlist[i])
	i<-i+1
}
colmarks<-setNames(listcolor, types)

## To add a exclusive_mark column when all the sample of a the set contains the same mark in the violin plot

mat$exclusive_mark<-"Mix"
for (type in types) {
	type_cols_subset<-type_cols[[type]]
	exclu<-rowSums(mat[, type_cols_subset, drop=FALSE]) == rowSums(mat[, sampleslist, drop=FALSE])
	mat$exclusive_mark[exclu]<-type
}
mat <- mat %>% relocate(exclusive_mark, .after = Category)
colmarks["Mix"] <- "black"

## Make the Plot

plot<-upset(mat, sampleslist, name="Peaks", 
      mode='exclusive_intersection',
      n_intersections=30, 
      sort_sets=FALSE,
      height_ratio = 0.75,
      base_annotations = list(
        'Shared peaks'=intersection_size(
          counts=FALSE, mapping=aes(fill=Category)) +
          scale_fill_manual(values=c("Distal_downstream"="#B8B5B3","Terminator"="#B233FF",
                                     "Gene_body"="#3358FF","Promoter"="#FF33E0","Distal_upstream"="#2e2e2e"),
                            name="Distance category")
      ),
      annotations = list(
        'Distance to closest gene' = (
          ggplot(mapping = aes(x=intersection, y=Distance, fill = exclusive_mark)) +
            geom_violin(scale="width", na.rm=TRUE, color = "black") +
            scale_y_continuous(trans = "log10",
                               labels=scales::label_number(accuracy = 1, scale_cut = scales::cut_si("bp"))) +
            scale_fill_manual(values=colmarks, name="Exclusive marks"))
      ),
      queries = queries,
	  set_sizes = (upset_set_size() + ylab("Total peaks") +
        theme(axis.text.x = element_text(angle = 45))),
      matrix = (intersection_matrix(geom = geom_point(shape = "circle", size = 3),
          segment = geom_segment(linewidth = 1.5),
          outline_color = list(active = alpha("white", 0),inactive = alpha("white", 0))) +
          scale_color_manual(values = c("TRUE" = "black", "FALSE" = alpha("white", 0)),
            labels = c("TRUE" = "yes", "FALSE" = "no"),
            breaks = c("TRUE", "FALSE"),
            guide = "none") +
          theme(axis.ticks = element_blank(),
            panel.grid = element_blank())
      ),
      themes = upset_modify_themes(
        list(
          "default" = theme(
            panel.grid.major.x = element_blank(),
            axis.ticks.y = element_line(linewidth = unit(0.25, "pt"), color = "#2e2e2e")
          ),
          "intersections_matrix" = theme(
            panel.grid = element_blank(),
            panel.grid.major.y = element_line(color = c("#CFCCCF", "white"), linewidth = 5)
          ),
          "Intersection size" = theme(
            panel.grid = element_blank(),
          ),
          "overall_sizes" = theme(
            panel.grid = element_blank(),
            axis.ticks.x = element_line(linewidth = 0.25, color = "#2e2e2e")
          )
        )
      ),
      stripes = alpha("white", 0)
)

pdf(output,height=max(figsize/2,10),width=10)
print(plot)
dev.off()
