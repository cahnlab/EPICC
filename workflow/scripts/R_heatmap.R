#!/usr/bin/env Rscript
# Heatmap__ for one computeMatrix matrix (plotting_heatmap_on_targetfile),
# drawn in R with ComplexHeatmap from the matrix itself (plot_heatmap.R).
#
# Usage: R_heatmap.R <matrix.gz> <columns.tsv> <regions|tss|tes> <env>
#                    <type|sample|default|shared> <mean|median> <no|mean|median>
#                    <params_heatmap.txt> <Heatmap__.pdf> <sorted_regions.bed>
#                    [theme.json]
#
# columns.tsv: key, label, colour, track, series, type per matrix column, in
# column order (matrix_sample_labels.py). params_heatmap.txt: the
# "--zMin ... --zMax ..." computing_matrix_scales writes, one value per column
# in column order; empty under heatmaps_scales: default and shared, when the
# scales are worked out here (default_heatmap_zlim: one per data type, or one
# over every column for shared). The top profiles'
# y-limits come from the very curves drawn, as in Profile__.
#
# Rows sort within each region group (heatmaps_sort_options: mean/median =
# descend by that summary, no = keep), on every column, except in the all env:
# there on the non-mC columns, so mC is read against the other marks' order (a
# run with nothing but mC sorts on mC).

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 10:11) {
  stop("usage: R_heatmap.R matrix columns.tsv kind env mode average sort ",
       "params_heatmap heatmap.pdf sorted_regions.bed [theme.json]")
}
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("matrix_io.R", "plot_scales.R", "plot_theme.R", "plot_profile.R", "plot_heatmap.R")) {
  source(file.path(script_dir, f))
}

kind <- args[3]
env <- args[4]
mode <- scale_mode(args[5])
average <- args[6]
if (!average %in% c("mean", "median")) stop("profiles_scale must be mean or median, got ", average)
sort_opt <- args[7]
sort <- switch(sort_opt, no = c("keep", "mean"), mean = c("descend", "mean"),
               median = c("descend", "median"), {
                 message("unclear heatmaps_sort_options '", sort_opt, "': no sorting done")
                 c("keep", "mean")
               })
th <- epicc_theme(if (length(args) == 11) args[11] else NULL)

m <- read_matrix_table(args[1])
display <- read_display_table(args[2])

keys <- m$samples$label
zlim <- read_heatmap_zlim(args[8], keys)
if (is.null(zlim)) zlim <- default_heatmap_zlim(m, display, shared = mode == "shared")

prep <- prepare_heatmap(m, display, th, sort[1], sort[2], average, env)
# The table is no longer needed: everything the figure draws is in `prep`.
m$table <- NULL
invisible(gc())

# Top profiles on Profile__'s limits (each column its own range under default,
# one range over every column under shared).
ylim <- profile_limits(m, prep$d, mode, per_column = mode == "default")
write_sorted_regions(m, prep$ord, args[10])
fig <- build_heatmap(m, display, th, prep, zlim, ylim, kind)
draw_heatmap_pdf(fig, args[9])
