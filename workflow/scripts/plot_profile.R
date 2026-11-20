# Average-profile drawing shared by the R metaplots.
#
# Sourced after matrix_io.R, plot_scales.R and plot_theme.R (R_profile.R does
# all four). The pieces are kept apart so other plots can reuse them -- the
# phase 3 heatmap draws its top profile from profile_data() and
# profile_limits(), so its curves and axes are the ones Profile__ shows:
#
#   d   <- profile_data(m, "mean", spread = "none")   # one row per sample x group x bin
#   lim <- profile_limits(m, d, "type")              # per-column y range of those curves
#   p   <- profile_ggplot(d, ...)                     # the ggplot itself
#
# and the three EPICC profile plots on top of them:
#
#   plot_profile_per_column(m, display, ...)  Profile__: a panel per matrix
#                                             column, a line per region group
#   plot_profile_per_group(m, display, ...)   Profile_pergroup__: a panel per
#                                             region group, a line per column
#   plot_profile_by_context(m, display, ...)  Profile_pergroup__ for mC and
#                                             all: a panel per region group x
#                                             mark (the context for mC)
#
# The curves and the y-limits come from the same group_means() call with the
# same `average`, so no curve leaves its own axis. (deeptools drew Profile__
# with its default mean even under profiles_scale: median, while the limits
# came from the median.)

suppressPackageStartupMessages({
  library(ggplot2)
})

if (!exists("%||%", baseenv())) `%||%` <- function(a, b) if (is.null(a)) b else a

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

# y range per matrix column, as data.frame(key, ymin, ymax) in column order,
# from the profiles `d` (profile_data()) the plot draws. heatmaps_scales
# type/sample: profile_ylims() (one range per mark or per column, the
# tenth-of-range padding: the y half of matrix_scales(), so the same numbers
# computing_matrix_scales writes). default: the single range deeptools shared
# across panels (default_ylim()) on every column, or with `per_column` each
# column's own default_ylim() -- for Profile__, whose panels each autoscale
# under default, so a CHH or ChIP panel isn't drawn on mCG's range. shared:
# one default_ylim() over every curve, on every column.
profile_limits <- function(m, d, mode = "type", per_column = FALSE) {
  keys <- m$samples$label
  if (mode %in% c("type", "sample")) return(profile_ylims(keys, d, mode))
  if (per_column) {
    y <- vapply(keys, function(k) default_ylim(d[as.character(d$sample) == k, , drop = FALSE]),
                numeric(2))
    return(data.frame(key = keys, ymin = y[1, ], ymax = y[2, ], stringsAsFactors = FALSE,
                      row.names = NULL))
  }
  y <- default_ylim(d)
  data.frame(key = keys, ymin = y[1], ymax = y[2], stringsAsFactors = FALSE)
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
# spec; `limits` is a data.frame with the facet variables plus ymin/ymax:
# every panel spans exactly its own limits. They are drawn as invisible
# points with no axis expansion, and the se/std band and the fill are clipped
# to them (fill runs down to the panel's ymin, not to 0), as deeptools applied
# --yMin/--yMax and clipped -- otherwise a band or a 0 baseline would widen
# the axis and two panels of one mark would end up on different scales. The
# curves need no clipping: the limits are their padded range. With NULL limits
# ggplot autoscales and fill runs to the lowest curve value of its own panel
# (the facet variables). `type` is the theme's profile_type.
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
  if (!is.null(limits)) {
    by <- setdiff(intersect(names(limits), names(d)), c("ymin", "ymax"))
    if (length(by) == 0) stop("profile limits share no facet variable with the data")
    row_key <- function(x) do.call(paste, c(lapply(x[by], as.character), sep = "\r"))
    i <- match(row_key(d), row_key(limits))
    if (anyNA(i)) stop("profile data rows without limits")
    lo <- limits$ymin[i]; hi <- limits$ymax[i]
    d$base <- lo
    d$lo <- pmin(pmax(d$lo, lo), hi)
    d$hi <- pmin(pmax(d$hi, lo), hi)
  } else {
    # Each panel autoscales (free_y), so its fill runs to its own lowest value.
    by <- if (is.null(facet)) character(0) else intersect(names(facet$params$facets), names(d))
    panel <- if (length(by)) do.call(paste, c(lapply(d[by], as.character), sep = "\r"))
             else rep("", nrow(d))
    base <- tapply(d$value, panel, function(v) if (all(is.na(v))) 0 else min(v, na.rm = TRUE))
    d$base <- unname(base[panel])
  }
  p <- ggplot(d)
  if (type == "fill") {
    p <- p + geom_ribbon(aes(x = .data$x, ymin = .data$base, ymax = .data$value,
                             fill = .data[[colour]], group = .data$.line),
                         alpha = 0.35, na.rm = TRUE, show.legend = FALSE)
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
    # A little room at each end so the end ticks of neighbouring panels
    # ("1.0Kb" | "-1.0") do not run into each other.
    scale_x_continuous(breaks = ticks$x, labels = ticks$label, expand = expansion(mult = 0.04)) +
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

# The display table a plot rule writes through matrix_sample_labels.py: one
# row per matrix column in column order, key, label, colour, track, series and
# the data type (heatmap colours and interpolation; blank in older tables).
read_display_table <- function(path) {
  cols <- c("key", "label", "colour", "track", "series", "type")
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  f <- strsplit(lines, "\t", fixed = FALSE)
  f <- lapply(f, function(x) c(x, rep("", length(cols)))[seq_along(cols)])
  d <- as.data.frame(do.call(rbind, f), stringsAsFactors = FALSE)
  if (nrow(d) == 0) d <- as.data.frame(matrix(character(0), 0, length(cols)), stringsAsFactors = FALSE)
  names(d) <- cols
  d
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
                                    ylab = NULL, data = NULL) {
  .check_display(m, display)
  d <- .with_display(data %||% profile_data(m, average, .spread(th)), display)
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
                                   ylab = NULL, data = NULL) {
  .check_display(m, display)
  d <- .with_display(data %||% profile_data(m, average, .spread(th)), display)
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

# Profile_pergroup__ for mC and all: matrix columns are sample x mark (the
# key's last field: the context for mC), so a panel per mark (CG | CHG | CHH
# in one row, a row per region group), a line per series (the label without
# its mark) in its display colour. Each panel spans the shared range of that
# mark's columns, so mCHH or a ChIP mark is never drawn on mCG's range.
plot_profile_by_context <- function(m, display, th, limits, average = "mean", kind = "regions",
                                    ylab = NULL, use_limits = TRUE, data = NULL) {
  .check_display(m, display)
  if (any(is.na(display$series) | display$series == "")) {
    stop("Display rows need a series name for every column")
  }
  d <- .with_display(data %||% profile_data(m, average, .spread(th)), display)
  marks <- column_marks(as.character(display$key))
  contexts <- unique(marks)
  d$context <- factor(column_marks(as.character(d$sample)), levels = contexts)
  strips <- group_strip_labels(m)
  d$region <- factor(strips[as.character(d$group)], levels = strips)
  series <- levels(d$series)
  cols <- setNames(display$colour[match(series, display$series)], series)
  dashes <- epicc_dashes(th, cols)
  ncol_leg <- epicc_legend_ncol(th, length(series))
  if (use_limits) {
    lim <- do.call(rbind, lapply(contexts, function(ctx) {
      y <- shared_ylim(limits[marks == ctx, , drop = FALSE])
      data.frame(context = factor(ctx, levels = contexts), ymin = y[1], ymax = y[2])
    }))
    lim <- merge(lim, data.frame(region = factor(strips, levels = strips)))
  } else {
    # heatmaps_scales: default: each panel on its own curves' range with
    # ggplot's 5% margins, as the autoscaled plot had it -- spelled out so a
    # band or fill is clipped to it like everywhere else.
    lim <- do.call(rbind, lapply(split(d, list(d$region, d$context), drop = TRUE), function(p) {
      y <- default_ylim(p)
      data.frame(region = p$region[1], context = p$context[1], ymin = y[1], ymax = y[2])
    }))
  }
  if (anyNA(lim$ymin)) lim <- NULL
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
