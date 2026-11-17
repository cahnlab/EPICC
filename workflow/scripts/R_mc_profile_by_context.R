#!/usr/bin/env Rscript
# Per-group metaplot for methylation, faceted on sequence context.
#
# deeptools' `plotProfile --perGroup` draws one panel per region group with
# every matrix column as a line. For mC the columns are sample x context, and
# line colours are keyed on the sample's levels_label (so a genotype keeps one
# colour across plots), which made a sample's mCG, mCHG and mCHH lines the
# same colour in one panel. Here each context gets its own panel, side by side
# in one row (CG | CHG | CHH), with one line per sample in it. A target with a
# separate 'nostrand' region group gets a second row.
#
# Usage:
#   Rscript R_mc_profile_by_context.R <profile_data.tab> <samples.tsv> <out.pdf> <y_label>
#
# profile_data.tab : `plotProfile --outFileNameData` output for the mC matrix
# samples.tsv      : label <TAB> context (mCG/mCHG/mCHH) <TAB> colour <TAB>
#                    series, one row per matrix column, in matrix order (no
#                    header). The series names a sample's line in every context
#                    panel: its label without the context, built by
#                    build_plot_labels, which drops the context from the label
#                    itself when only one context is plotted.

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: R_mc_profile_by_context.R <profile_data.tab> <samples.tsv> <out.pdf> <y_label>")
}
data_file <- args[1]
sample_file <- args[2]
out_file <- args[3]
y_label <- args[4]

fields <- strsplit(readLines(data_file), "\t", fixed = TRUE)
keys <- vapply(fields, `[`, "", 1)
bin_labels <- fields[[which(keys == "bin labels")[1]]]
rows <- fields[!keys %in% c("bin labels", "bins")]
if (length(rows) == 0) stop("No profile rows in ", data_file)

nbins <- length(rows[[1]]) - 2
# The bin-label row repeats once per sample; the first block names the ticks.
# strsplit drops trailing empty fields, so pad before indexing.
first_block <- c(bin_labels[-(1:2)], rep("", nbins))[seq_len(nbins)]
tick_at <- which(first_block != "")
tick_label <- first_block[tick_at]

profile <- do.call(rbind, lapply(rows, function(r) {
  values <- suppressWarnings(as.numeric(r[-(1:2)]))
  data.frame(label = r[1], group = r[2], bin = seq_along(values), value = values,
             stringsAsFactors = FALSE)
}))

samples <- read.delim(sample_file, header = FALSE, sep = "\t",
                      col.names = c("label", "context", "color", "series"),
                      stringsAsFactors = FALSE, comment.char = "")
if (any(is.na(samples$series) | samples$series == "")) {
  stop("samples.tsv needs a series name (4th column) for every column")
}
missing <- setdiff(unique(profile$label), samples$label)
if (length(missing) > 0) {
  stop("Profile samples without a context: ", paste(missing, collapse = ", "))
}
profile <- merge(profile, samples, by = "label", sort = FALSE)
contexts <- unique(samples$context)
series <- unique(samples$series)
groups <- unique(vapply(rows, `[`, "", 2))
profile$context <- factor(profile$context, levels = contexts)
profile$series <- factor(profile$series, levels = series)
profile$group <- factor(profile$group, levels = groups)
palette <- setNames(samples$color[match(series, samples$series)], series)
# Colours are keyed on levels_label, so samples that differ only in assay or
# replicate share one. Those get distinct dash patterns, in series order, so
# the lines stay apart; a sample with a colour to itself stays solid.
dash_cycle <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
dash_rank <- ave(seq_along(series), palette, FUN = seq_along)
dashes <- setNames(dash_cycle[(dash_rank - 1) %% length(dash_cycle) + 1], series)

p <- ggplot(profile, aes(x = bin, y = value, colour = series, linetype = series)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  scale_colour_manual(values = palette, name = NULL) +
  scale_linetype_manual(values = dashes, name = NULL) +
  scale_x_continuous(breaks = tick_at, labels = tick_label, expand = c(0.02, 0)) +
  labs(x = NULL, y = y_label) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.spacing.x = unit(1.4, "lines"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.key.width = unit(2.2, "lines")) +
  guides(colour = guide_legend(ncol = ceiling(length(series) / 25)),
         linetype = guide_legend(ncol = ceiling(length(series) / 25)))

# Contexts sit at very different levels (CG tens of percent, CHH a few), so
# each panel gets its own y-axis.
if (length(groups) > 1) {
  p <- p + facet_wrap(vars(group, context), ncol = length(contexts), scales = "free_y")
} else {
  p <- p + facet_wrap(vars(context), nrow = 1, scales = "free_y")
}

legend_width <- 1.5 + 0.09 * max(nchar(series)) * ceiling(length(series) / 25)
width <- 3.2 * length(contexts) + legend_width
height <- max(3.2 * length(groups) + 0.6, 0.22 * min(length(series), 25) + 1)
ggsave(out_file, p, width = width, height = height, limitsize = FALSE)
