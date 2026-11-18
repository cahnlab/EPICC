#!/usr/bin/env Rscript
# Profile__ and Profile_pergroup__ for one computeMatrix matrix
# (plotting_profile_on_targetfile), drawn in R from the matrix itself.
#
# Usage: R_profile.R <matrix.gz> <columns.tsv> <regions|tss|tes> <env>
#                    <type|sample|default> <mean|median>
#                    <heat_mcg> <heat_mchg> <heat_mchh>
#                    <Profile__.pdf> <Profile_pergroup__.pdf> [theme.json]
#
# columns.tsv: key, label, colour, track, series per matrix column, in column
# order (matrix_sample_labels.py). The y-limits are worked out here with
# plot_scales.R from the same averages the curves are drawn with
# (heatmaps_scales, profiles_scale), so they are the ones
# computing_matrix_scales writes for the heatmaps.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 11:12) {
  stop("usage: R_profile.R matrix columns.tsv kind env mode average heat_mcg heat_mchg heat_mchh ",
       "profile.pdf profile_pergroup.pdf [theme.json]")
}
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("matrix_io.R", "plot_scales.R", "plot_theme.R", "plot_profile.R")) {
  source(file.path(script_dir, f))
}

kind <- args[3]
env <- args[4]
mode <- args[5]
if (!mode %in% c("type", "sample", "default")) {
  message(mode, " unknown. Returning default")
  mode <- "default"
}
average <- args[6]
if (!average %in% c("mean", "median")) stop("profiles_scale must be mean or median, got ", average)
mc <- c(mCG = as.numeric(args[7]), mCHG = as.numeric(args[8]), mCHH = as.numeric(args[9]))
th <- epicc_theme(if (length(args) == 12) args[12] else NULL)

m <- read_matrix(args[1])
display <- read.delim(args[2], header = FALSE, sep = "\t", quote = "", comment.char = "",
                      col.names = c("key", "label", "colour", "track", "series"),
                      colClasses = "character")
limits <- profile_limits(m, mode, average, mc)
ylab <- paste(average, if (env == "mC") "methylation" else "signal")

p1 <- plot_profile_per_column(m, display, th, limits, average, kind, ylab)
epicc_save(p1$plot, args[10], p1$width, p1$height)

p2 <- if (env == "mC") {
  # heatmaps_scales: default left each context panel to autoscale before; keep that.
  plot_profile_by_context(m, display, th, limits, average, kind, ylab,
                          use_limits = mode != "default")
} else {
  plot_profile_per_group(m, display, th, limits, average, kind, ylab)
}
epicc_save(p2$plot, args[11], p2$width, p2$height)
