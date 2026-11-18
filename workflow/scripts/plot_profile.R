# Average-profile drawing shared by the R metaplots.
#
# Sourced after matrix_io.R, plot_scales.R and plot_theme.R (R_profile.R does
# all four). The pieces are kept apart so other plots can reuse them -- the
# phase 3 heatmap draws its top profile from profile_data() and
# profile_limits(), so its curves and axes are the ones Profile__ shows:
#
#   d   <- profile_data(m, "mean", spread = "none")   # one row per sample x group x bin
#   lim <- profile_limits(m, "type", "mean", mc_zmax) # per-column y range
#   p   <- profile_ggplot(d, ...)                     # the ggplot itself
#
# and the three EPICC profile plots on top of them:
#
#   plot_profile_per_column(m, display, ...)  Profile__: a panel per matrix
#                                             column, a line per region group
#   plot_profile_per_group(m, display, ...)   Profile_pergroup__: a panel per
#                                             region group, a line per column
#   plot_profile_by_context(m, display, ...)  Profile_pergroup__ for mC: a
#                                             panel per region group x context
#
# The curves and the y-limits come from the same group_means() call with the
# same `average`, so no curve leaves its own axis. (deeptools drew Profile__
# with its default mean even under profiles_scale: median, while the limits
# came from the median.)

suppressPackageStartupMessages({
  library(ggplot2)
})

# Average profile per matrix column and region group, in long form: sample
# (column key), group, x (0-based bin), value, and lo/hi for a shaded band
# when `spread` is "se" or "std" (NA otherwise). The spread is numpy's nanstd
# (population, missing values ignored) of each bin over the group's regions;
# se divides it by the square root of the regions with a value, as
# plotProfile --plotType se does.
profile_data <- function(m, average = "mean", spread = c("none", "se", "std")) {
  spread <- match.arg(spread)
  d <- group_means(m, average)
  d$lo <- NA_real_
  d$hi <- NA_real_
  if (spread == "none") return(d)
  k <- 0
  for (i in seq_len(nrow(m$samples))) {
    cols <- sample_columns(m, i)
    for (g in seq_len(nrow(m$groups))) {
      rows <- if (m$groups$n[g] == 0) integer(0) else m$groups$first[g]:m$groups$last[g]
      v <- m$values[rows, cols, drop = FALSE]
      n <- colSums(!is.na(v))
      mu <- colSums(v, na.rm = TRUE) / n
      dev <- sweep(v, 2, mu)
      sd <- sqrt(colSums(dev^2, na.rm = TRUE) / n)
      if (spread == "se") sd <- sd / sqrt(n)
      sd[n == 0] <- NA_real_
      idx <- k + seq_along(cols)
      d$lo[idx] <- d$value[idx] - sd
      d$hi[idx] <- d$value[idx] + sd
      k <- k + length(cols)
    }
  }
  d
}

# y range per matrix column, as data.frame(key, ymin, ymax) in column order.
# heatmaps_scales type/sample: matrix_scales() (one range per mark or per
# column, the tenth-of-range padding). default: the single range deeptools
# shared across panels (default_scales()) on every column.
profile_limits <- function(m, mode = "type", average = "mean",
                           mc_zmax = c(mCG = 100, mCHG = 50, mCHH = 20)) {
  keys <- m$samples$label
  if (mode %in% c("type", "sample")) {
    s <- matrix_scales(m, mode, average, mc_zmax)
    return(data.frame(key = keys, ymin = s$ymin, ymax = s$ymax, stringsAsFactors = FALSE))
  }
  d <- default_scales(m, average)
  data.frame(key = keys, ymin = d$ymin, ymax = d$ymax, stringsAsFactors = FALSE)
}

# Region-group strip labels as plotHeatmap/plotProfile's --regionsLabel gave
# them (matrix_region_labels.py): "label(regions)".
group_strip_labels <- function(m) {
  setNames(paste0(m$groups$label, "(", m$groups$n, ")"), m$groups$label)
}

# x-axis ticks for the matrix: what the rule passed deeptools per matrix kind
# (tss: --refPointLabel start; tes: --refPointLabel end; regions: --startLabel
# start --endLabel end).
profile_ticks <- function(m, kind = "regions") {
  ref <- if (kind == "tes") "end" else "start"
  matrix_ticks(m, 1, ref_label = ref, start_label = "start", end_label = "end")
}

# The ggplot for long profile data `d` (profile_data() joined with display
# columns). `colour`/`linetype` name columns of d; `colours`/`linetypes` are
# their manual scales (NULL linetypes = all solid). `facet` is a facet_*()
# spec; `limits` is a data.frame with the facet variables plus ymin/ymax,
# drawn as invisible points so every panel spans at least its own limits
# (pass NULL to let ggplot autoscale). `type` is the theme's profile_type.
profile_ggplot <- function(d, th, ticks, colour, colours, linetype = NULL, linetypes = NULL,
                           facet = NULL, limits = NULL, ylab = NULL, type = th$profile_type,
                           show_colour_legend = TRUE, legend_ncol = 1L) {
  aes_line <- if (is.null(linetype)) {
    aes(x = .data$x, y = .data$value, colour = .data[[colour]], group = .data$.line)
  } else {
    aes(x = .data$x, y = .data$value, colour = .data[[colour]], group = .data$.line,
        linetype = .data[[linetype]])
  }
  d$.line <- interaction(d$sample, d$group, drop = TRUE)
  p <- ggplot(d)
  if (type == "fill") {
    p <- p + geom_area(aes(x = .data$x, y = .data$value, fill = .data[[colour]], group = .data$.line),
                       position = "identity", alpha = 0.35, na.rm = TRUE, show.legend = FALSE)
  } else if (type %in% c("se", "std")) {
    p <- p + geom_ribbon(aes(x = .data$x, ymin = .data$lo, ymax = .data$hi,
                             fill = .data[[colour]], group = .data$.line),
                         alpha = 0.25, na.rm = TRUE, show.legend = FALSE)
  }
  p <- p + geom_line(aes_line, linewidth = th$line_width, na.rm = TRUE)
  if (!is.null(limits)) {
    lim <- rbind(transform(limits, y = limits$ymin), transform(limits, y = limits$ymax))
    lim <- lim[!is.na(lim$y), , drop = FALSE]
    p <- p + geom_blank(data = lim, aes(y = .data$y))
  }
  p <- p +
    scale_colour_manual(values = colours, name = NULL,
                        guide = if (show_colour_legend) guide_legend(ncol = legend_ncol) else "none") +
    scale_fill_manual(values = colours, guide = "none") +
    scale_x_continuous(breaks = ticks$x, labels = ticks$label, expand = c(0, 0)) +
    scale_y_continuous(expand = if (is.null(limits)) waiver() else c(0, 0)) +
    labs(x = NULL, y = ylab) +
    epicc_ggtheme(th)
  if (!is.null(linetype)) {
    p <- p + scale_linetype_manual(values = linetypes, name = NULL,
                                   guide = guide_legend(ncol = legend_ncol))
  }
  if (!is.null(facet)) p <- p + facet
  p
}

# Join display rows (key, label, colour, track, series; one per matrix column
# in column order, from matrix_sample_labels.py) onto profile data.
.with_display <- function(d, display) {
  i <- match(as.character(d$sample), display$key)
  if (anyNA(i)) stop("No display row for column(s) ",
                     paste(unique(as.character(d$sample)[is.na(i)]), collapse = ", "))
  d$label <- factor(display$label[i], levels = unique(display$label))
  d$colour <- display$colour[i]
  d$track <- display$track[i]
  d$series <- factor(display$series[i], levels = unique(display$series))
  d
}

.check_display <- function(m, display) {
  if (!identical(as.character(display$key), as.character(m$samples$label))) {
    stop("Display rows do not match the matrix columns (want ",
         paste(m$samples$label, collapse = ", "), ")")
  }
}

.panel_grid <- function(n, max_cols = 4L) {
  ncol <- min(n, max_cols)
  c(ncol = ncol, nrow = ceiling(n / ncol))
}

# Profile__: a panel per matrix column (titled with its label, drawn in its
# colour), a line per region group (told apart by dash pattern), each panel on
# its own y range from `limits` (profile_limits()).
plot_profile_per_column <- function(m, display, th, limits, average = "mean", kind = "regions",
                                    ylab = NULL) {
  .check_display(m, display)
  d <- .with_display(profile_data(m, average, .spread(th)), display)
  strips <- group_strip_labels(m)
  d$region <- factor(strips[as.character(d$group)], levels = strips)
  # A label can repeat when two keys share one (it shouldn't: check_unique_labels);
  # facet on the key and show the label.
  d$panel <- factor(as.character(d$sample), levels = m$samples$label)
  panel_labels <- setNames(display$label, display$key)
  lim <- data.frame(panel = factor(limits$key, levels = m$samples$label),
                    ymin = limits$ymin, ymax = limits$ymax)
  cols <- setNames(display$colour, display$key)
  d$colkey <- as.character(d$sample)
  nreg <- length(strips)
  dashes <- setNames(th$dash_cycle[(seq_len(nreg) - 1) %% length(th$dash_cycle) + 1], strips)
  grid <- .panel_grid(nrow(m$samples))
  p <- profile_ggplot(
    d, th, profile_ticks(m, kind), colour = "colkey", colours = cols,
    linetype = if (nreg > 1) "region" else NULL, linetypes = dashes,
    facet = facet_wrap(vars(.data$panel), ncol = grid[["ncol"]], scales = "free_y",
                       labeller = as_labeller(panel_labels)),
    limits = lim, ylab = ylab, show_colour_legend = FALSE)
  if (nreg == 1) p <- p + theme(legend.position = "none")
  legend_w <- if (nreg > 1) epicc_legend_width(th, strips) else 0
  list(plot = p,
       width = th$panel_width * grid[["ncol"]] + legend_w,
       height = th$panel_height * grid[["nrow"]] + 0.6)
}

# Profile_pergroup__ (non-mC): a panel per region group, a line per matrix
# column in its display colour, every panel on one y range: the lowest ymin
# and highest ymax over the columns (shared_ylim()).
plot_profile_per_group <- function(m, display, th, limits, average = "mean", kind = "regions",
                                   ylab = NULL) {
  .check_display(m, display)
  d <- .with_display(profile_data(m, average, .spread(th)), display)
  strips <- group_strip_labels(m)
  d$region <- factor(strips[as.character(d$group)], levels = strips)
  # One legend entry per column, named by its label.
  cols <- setNames(display$colour, display$label)
  dashes <- epicc_dashes(th, cols)
  ncol_leg <- epicc_legend_ncol(th, length(cols))
  y <- shared_ylim(limits)
  lim <- data.frame(region = factor(strips, levels = strips), ymin = y[1], ymax = y[2])
  p <- profile_ggplot(
    d, th, profile_ticks(m, kind), colour = "label", colours = cols,
    linetype = "label", linetypes = dashes,
    facet = facet_wrap(vars(.data$region), nrow = 1),
    limits = if (anyNA(y)) NULL else lim, ylab = ylab, legend_ncol = ncol_leg)
  list(plot = p,
       width = th$panel_width * length(strips) + epicc_legend_width(th, names(cols), ncol_leg),
       height = max(th$panel_height + 0.6, epicc_legend_height(th, length(cols))))
}

# Profile_pergroup__ for mC: matrix columns are sample x context, so a
# panel per context (CG | CHG | CHH in one row, a row per region group), a
# line per series (the label without its context) in its display colour.
# Each context panel spans the shared range of that context's columns.
plot_profile_by_context <- function(m, display, th, limits, average = "mean", kind = "regions",
                                    ylab = NULL, use_limits = TRUE) {
  .check_display(m, display)
  if (any(is.na(display$series) | display$series == "")) {
    stop("Display rows need a series name for every column")
  }
  d <- .with_display(profile_data(m, average, .spread(th)), display)
  contexts <- unique(display$track)
  d$context <- factor(d$track, levels = contexts)
  strips <- group_strip_labels(m)
  d$region <- factor(strips[as.character(d$group)], levels = strips)
  series <- levels(d$series)
  cols <- setNames(display$colour[match(series, display$series)], series)
  dashes <- epicc_dashes(th, cols)
  ncol_leg <- epicc_legend_ncol(th, length(series))
  lim <- NULL
  if (use_limits) {
    lim <- do.call(rbind, lapply(contexts, function(ctx) {
      y <- shared_ylim(limits[display$track == ctx, , drop = FALSE])
      data.frame(context = factor(ctx, levels = contexts), ymin = y[1], ymax = y[2])
    }))
    lim <- merge(lim, data.frame(region = factor(strips, levels = strips)))
    if (anyNA(lim$ymin)) lim <- NULL
  }
  # Contexts sit at very different levels (CG tens of percent, CHH a few), so
  # each panel gets its own y-axis.
  facet <- if (length(strips) > 1) {
    facet_wrap(vars(.data$region, .data$context), ncol = length(contexts), scales = "free_y")
  } else {
    facet_wrap(vars(.data$context), nrow = 1, scales = "free_y")
  }
  p <- profile_ggplot(
    d, th, profile_ticks(m, kind), colour = "series", colours = cols,
    linetype = "series", linetypes = dashes, facet = facet, limits = lim,
    ylab = ylab, legend_ncol = ncol_leg)
  list(plot = p,
       width = th$panel_width * length(contexts) + epicc_legend_width(th, series, ncol_leg),
       height = max(th$panel_height * length(strips) + 0.6, epicc_legend_height(th, length(series))))
}

.spread <- function(th) {
  if (th$profile_type %in% c("se", "std")) th$profile_type else "none"
}
