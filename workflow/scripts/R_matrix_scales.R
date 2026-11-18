#!/usr/bin/env Rscript
# Per-column heatmap and profile limits for one computeMatrix matrix
# (computing_matrix_scales). Writes the deeptools flag files the plot rules
# read: "--zMin ... --zMax ..." and "--yMin ... --yMax ...", one value per
# matrix column in column order; both empty for heatmaps_scales: default.
#
# Usage: R_matrix_scales.R <matrix.gz> <type|sample|default> <mean|median>
#                          <heat_mcg> <heat_mchg> <heat_mchh>
#                          <params_heatmap.txt> <params_profile.txt>
#
# The logic is plot_scales.R's (see there); the y bounds use matrix_io.R's
# group_means, the averaging the R profiles draw with.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("usage: R_matrix_scales.R matrix mode average heat_mcg heat_mchg heat_mchh ",
       "params_heatmap params_profile")
}
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
source(file.path(script_dir, "matrix_io.R"))
source(file.path(script_dir, "plot_scales.R"))

mode <- scale_mode(args[2])
average <- args[3]
if (!average %in% c("mean", "median")) stop("profiles_scale must be mean or median, got ", average)
mc <- c(mCG = as.numeric(args[4]), mCHG = as.numeric(args[5]), mCHH = as.numeric(args[6]))

flags <- c(heatmap = "", profile = "")
if (mode != "default") {
  m <- read_matrix(args[1])
  flags <- scale_flags(matrix_scales(m, mode, average, mc))
}
# No trailing newline: the rules read these with $(cat ...).
cat(flags[["heatmap"]], file = args[7])
cat(flags[["profile"]], file = args[8])
