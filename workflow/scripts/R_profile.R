#!/usr/bin/env Rscript
# Profile__ and Profile_pergroup__ for one computeMatrix matrix
# (plotting_profile_on_targetfile), drawn in R from the matrix itself.
#
# Usage: R_profile.R <matrix.gz> <columns.tsv> <regions|tss|tes> <env>
#                    <type|sample|default> <mean|median>
#                    <Profile__.pdf> <Profile_pergroup__.pdf> [theme.json]
#
# columns.tsv: key, label, colour, track, series per matrix column, in column
# order (matrix_sample_labels.py). The y-limits are worked out here, with
# plot_scales.R, from the very curves drawn (heatmaps_scales,
# profiles_scale) -- the y values computing_matrix_scales writes for the
# heatmaps.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 8:9) {
  stop("usage: R_profile.R matrix columns.tsv kind env mode average ",
       "profile.pdf profile_pergroup.pdf [theme.json]")
}
script_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
for (f in c("matrix_io.R", "plot_scales.R", "plot_theme.R", "plot_profile.R")) {
  source(file.path(script_dir, f))
}

kind <- args[3]
env <- args[4]
mode <- scale_mode(args[5])
average <- args[6]
if (!average %in% c("mean", "median")) stop("profiles_scale must be mean or median, got ", average)
th <- epicc_theme(if (length(args) == 9) args[9] else NULL)

m <- read_matrix(args[1])
display <- read.delim(args[2], header = FALSE, sep = "\t", quote = "", comment.char = "",
                      col.names = c("key", "label", "colour", "track", "series"),
                      colClasses = "character")
# One pass over the matrix for the curves; the limits come from them.
d <- profile_data(m, average, .spread(th))
limits <- profile_limits(m, d, mode)
ylab <- paste(average, switch(env, mC = "methylation", all = "signal / methylation", "signal"))

# Under heatmaps_scales: default each Profile__ panel autoscales on its own
# curves, as the mC per-context panels below do.
limits1 <- if (mode == "default") profile_limits(m, d, mode, per_column = TRUE) else limits
p1 <- plot_profile_per_column(m, display, th, limits1, average, kind, ylab, data = d)
epicc_save(p1$plot, args[7], p1$width, p1$height)

p2 <- if (env %in% c("mC", "all")) {
  # A panel per mark (per context for mC): `all` mixes methylation (tens of
  # percent) with ChIP/RNA signal, which one shared range would flatten.
  # heatmaps_scales: default left each context panel to autoscale before; keep that.
  plot_profile_by_context(m, display, th, limits, average, kind, ylab,
                          use_limits = mode != "default", data = d)
} else {
  plot_profile_per_group(m, display, th, limits, average, kind, ylab, data = d)
}
epicc_save(p2$plot, args[8], p2$width, p2$height)
