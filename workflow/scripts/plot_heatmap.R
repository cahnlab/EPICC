# Heatmaps of computeMatrix matrices, drawn with ComplexHeatmap.
#
# Sourced after matrix_io.R, plot_scales.R, plot_theme.R and plot_profile.R
# (R_heatmap.R does all five). One ComplexHeatmap `Heatmap` per matrix column,
# joined into a HeatmapList that shares one row order and one split into
# region groups; each column has
#   - its own colour ramp, missing colour and interpolation for its data type
#     (theme: heatmap_colors, heatmap_missing_color, heatmap_interpolation),
#     over its own zMin/zMax (computing_matrix_scales, looked up by key),
#   - an average profile on top, from the same profile_data()/profile_limits()
#     Profile__ draws (plot_profile.R), so the two show the same curves on the
#     same axes,
#   - x ticks and a colour bar below.
#
# Pixels, not regions: a column's image is computed here at the size it is
# drawn (heatmap_raster_dpi) and placed with no further smoothing, so neither
# the PDF viewer nor ComplexHeatmap's raster filter decides how it looks, and
# memory depends on the figure, not on the number of regions:
#   rows   more regions than pixel rows are averaged in runs (bin_rows); fewer
#          are spread over the pixels like the bins below;
#   bins   bilinear: linear between bin centres (matplotlib's convention:
#          pixel centres mapped onto bin centres, edges clamped); nearest:
#          each bin a block of pixels; more bins than pixels are averaged;
#   NaN    a pixel whose interpolation would take in a missing value takes the
#          value of its nearest bin instead, so a missing mC bin stays a sharp
#          block of the missing colour and never bleeds into its neighbours.

suppressPackageStartupMessages({
  library(grid)
  library(ComplexHeatmap)
})
ht_opt$message <- FALSE

# --- resampling ---------------------------------------------------------------

# Resample the rows of `x` (a numeric matrix) to `n` rows. Fewer rows than
# there are: averaged in consecutive runs, NaN-aware (bin_rows). More: each
# output row's centre is mapped back onto the input rows' centres and takes
# the nearest row ("nearest") or the linear mix of the two around it
# ("bilinear"); when that mix would take in an NA, the nearest row's value.
resample_rows <- function(x, n, mode = c("bilinear", "nearest")) {
  mode <- match.arg(mode)
  nin <- nrow(x)
  if (nin == n) return(unname(x))
  if (nin == 0) return(matrix(NA_real_, n, ncol(x)))
  if (nin > n) {
    out <- bin_rows(x, n)
    attr(out, "rows") <- NULL
    return(out)
  }
  pos <- (seq_len(n) - 0.5) * nin / n - 0.5          # 0-based input coordinate
  pos <- pmin(pmax(pos, 0), nin - 1)
  i0 <- floor(pos)
  t <- pos - i0
  i1 <- pmin(i0 + 1, nin - 1)
  near <- ifelse(t < 0.5, i0, i1) + 1
  if (mode == "nearest") return(unname(x[near, , drop = FALSE]))
  out <- x[i0 + 1, , drop = FALSE] * (1 - t) + x[i1 + 1, , drop = FALSE] * t
  bad <- is.na(out)
  if (any(bad)) out[bad] <- x[near, , drop = FALSE][bad]
  unname(out)
}

# Resample a regions x bins matrix to an image of nr x nc pixels: rows first,
# then bins, each as resample_rows() describes.
resample_image <- function(x, nr, nc, mode = "bilinear") {
  x <- resample_rows(x, nr, mode)
  t(resample_rows(t(x), nc, mode))
}

# The per-column z limits computing_matrix_scales wrote ("--zMin a b ...
# --zMax c d ...", one value per matrix column in column order), as a scales
# table keyed by the matrix's column keys (scale_for_column looks them up);
# NULL when the file is empty (heatmaps_scales: default).
read_heatmap_zlim <- function(path, keys) {
  words <- scan(path, what = "", quiet = TRUE)
  if (length(words) == 0) return(NULL)
  grab <- function(flag) {
    i <- match(flag, words)
    if (is.na(i)) stop(path, ": no ", flag)
    as.numeric(words[i + seq_along(keys)])
  }
  z <- data.frame(key = keys, zmin = grab("--zMin"), zmax = grab("--zMax"),
                  ymin = NA_real_, ymax = NA_real_, stringsAsFactors = FALSE)
  if (anyNA(z$zmin) || anyNA(z$zmax)) {
    stop(path, ": expected ", length(keys), " --zMin/--zMax values (one per matrix column)")
  }
  z
}
# --- colours ------------------------------------------------------------------

# A function mapping values to colours for one column: the ramp's stops spread
# evenly from zmin to zmax (as matplotlib spreads a colour map over vmin..vmax),
# interpolated in RGB as matplotlib does, values beyond the ends clamped, NA
# to the missing colour.
heatmap_colour_fun <- function(ramp, zmin, zmax, missing) {
  if (!is.finite(zmin) || !is.finite(zmax)) { zmin <- 0; zmax <- 1 }
  if (zmax <= zmin) zmax <- zmin + max(abs(zmin) * 1e-6, 1e-12)
  breaks <- seq(zmin, zmax, length.out = length(ramp))
  rgb <- grDevices::col2rgb(ramp) / 255
  missing <- grDevices::rgb(t(grDevices::col2rgb(missing) / 255))
  f <- function(v) {
    out <- rep(missing, length(v))
    ok <- !is.na(v)
    if (any(ok)) {
      w <- pmin(pmax(v[ok], zmin), zmax)
      ch <- vapply(1:3, function(k) stats::approx(breaks, rgb[k, ], xout = w)$y,
                   numeric(length(w)))
      if (is.null(dim(ch))) ch <- matrix(ch, nrow = 1)
      out[ok] <- grDevices::rgb(ch[, 1], ch[, 2], ch[, 3])
    }
    out
  }
  attr(f, "breaks") <- breaks
  attr(f, "colors") <- ramp
  f
}

# --- the figure ---------------------------------------------------------------

# The data type of each column, for its colours and interpolation: the
# display table's type column (combined_analysis.smk), or, for a table
# without it, mC for the methylation contexts and 'default' for anything else.
column_types <- function(display) {
  type <- display$type
  if (is.null(type)) type <- rep("", nrow(display))
  fallback <- ifelse(grepl("^mC(G|HG|HH)$", column_marks(display$key)),
                     column_marks(display$key), "default")
  ifelse(is.na(type) | type == "", fallback, type)
}

# The columns the rows are sorted on: in the all env the non-mC ones (a run
# with nothing but mC sorts on mC), otherwise every column.
sort_columns <- function(display, env) {
  if (env != "all") return(NULL)
  mc <- grepl("^mC", column_types(display))
  if (all(mc)) return(NULL)
  display$key[!mc]
}

# deeptools' --outFileSortedRegions: a BED12-like table with the region group
# in a 13th column, in the order the heatmap draws the rows.
write_sorted_regions <- function(m, ord, file) {
  r <- m$regions[ord, , drop = FALSE]
  size <- r$end - r$start
  out <- data.frame(r$chrom, r$start, r$end, r$name, r$score, r$strand,
                    r$start, r$end, 0L, 1L, size, 0L, as.character(r$group))
  con <- file(file, "w")
  on.exit(close(con))
  writeLines(paste("#chrom", "start", "end", "name", "score", "strand", "thickStart",
                   "thickEnd", "itemRGB", "blockCount", "blockSizes", "blockStart",
                   "deepTools_group", sep = "\t"), con)
  utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Everything the figure needs from the matrix, one sample block at a time
# (matrix_block), so only one sample's values are copied out at once:
#   images  per column, per region group, the group's rows in heatmap order
#           averaged down to at most the pixel rows the group will get
#   d       profile_data() of every column (the curves Profile__ draws)
#   ord     the row order (for the sorted-regions BED)
prepare_heatmap <- function(m, display, th, sort_method = "descend", sort_using = "mean",
                            average = "mean", env = "most") {
  .check_display(m, display)
  ord <- row_order(m, sort_method, sort_using, sort_columns(display, env))
  total_px <- max(1L, as.integer(round(th$heatmap_height * th$heatmap_raster_dpi)))
  nreg <- nrow(m$regions)
  spread <- .spread(th)
  images <- vector("list", nrow(m$samples))
  d <- vector("list", nrow(m$samples))
  for (i in seq_len(nrow(m$samples))) {
    ms <- matrix_sample(m, i)
    d[[i]] <- profile_data(ms, average, spread)
    v <- ms$values
    images[[i]] <- lapply(seq_len(nrow(m$groups)), function(g) {
      if (m$groups$n[g] == 0) return(matrix(NA_real_, 0, ncol(v)))
      rows <- ord[m$groups$first[g]:m$groups$last[g]]
      target <- max(1L, as.integer(ceiling(total_px * m$groups$n[g] / max(nreg, 1))))
      out <- v[rows, , drop = FALSE]
      if (nrow(out) > target) {
        out <- bin_rows(out, target)
        attr(out, "rows") <- NULL
      }
      out
    })
    rm(ms, v)
  }
  d <- do.call(rbind, d)
  d$sample <- factor(as.character(d$sample), levels = m$samples$label)
  list(images = images, d = d, ord = ord)
}

# The custom annotation drawing one column's average profiles: a line per
# region group (dash pattern per group, as Profile__), in the column's
# colour, on [ymin, ymax]; x on the heatmap's bin axis, so a profile point
# sits above its bin's pixels.
.profile_annotation <- function(dk, nbins, ymin, ymax, colour, dashes, th, ticks, height) {
  fun <- function(index, k, n) {
    yr <- c(ymin, ymax)
    if (!all(is.finite(yr)) || yr[2] <= yr[1]) {
      v <- c(dk$value, dk$lo, dk$hi)
      v <- v[is.finite(v)]
      yr <- if (length(v)) range(v) else c(0, 1)
      if (yr[2] <= yr[1]) yr <- yr + c(-0.5, 0.5)
    }
    pushViewport(viewport(xscale = c(-0.5, nbins - 0.5), yscale = yr, clip = "off"))
    grid.rect(gp = gpar(col = "grey30", fill = NA, lwd = 0.6))
    for (x in ticks$x) {
      grid.lines(unit(c(x, x), "native"), unit(c(0, 1), "npc"),
                 gp = gpar(col = "grey88", lwd = 0.4))
    }
    pushViewport(viewport(xscale = c(-0.5, nbins - 0.5), yscale = yr, clip = "on"))
    groups <- levels(droplevels(dk$group))
    for (gi in seq_along(groups)) {
      p <- dk[as.character(dk$group) == groups[gi], , drop = FALSE]
      p <- p[order(p$x), , drop = FALSE]
      ok <- !is.na(p$value)
      if (!any(ok)) next
      if (th$profile_type == "fill") {
        base <- yr[1]                      # as Profile__: fill down to the panel's ymin
        q <- p[ok, ]
        grid.polygon(unit(c(q$x, rev(q$x)), "native"),
                     unit(c(pmin(pmax(q$value, yr[1]), yr[2]), rep(base, nrow(q))), "native"),
                     gp = gpar(fill = grDevices::adjustcolor(colour, 0.35), col = NA))
      } else if (th$profile_type %in% c("se", "std")) {
        q <- p[ok & !is.na(p$lo) & !is.na(p$hi), ]
        if (nrow(q) > 1) {
          grid.polygon(unit(c(q$x, rev(q$x)), "native"),
                       unit(c(pmin(pmax(q$hi, yr[1]), yr[2]), rev(pmin(pmax(q$lo, yr[1]), yr[2]))),
                            "native"),
                       gp = gpar(fill = grDevices::adjustcolor(colour, 0.25), col = NA))
        }
      }
      grid.lines(unit(p$x, "native"), unit(p$value, "native"),
                 gp = gpar(col = colour, lty = dashes[[groups[gi]]],
                           lwd = th$line_width * ggplot2::.pt))
    }
    popViewport()
    grid.yaxis(at = pretty(yr, 3)[pretty(yr, 3) >= yr[1] & pretty(yr, 3) <= yr[2]],
               gp = epicc_gpar(th, 0.6))
    popViewport()
  }
  AnnotationFunction(fun = fun, which = "column", height = height,
                     n = 1, subsettable = FALSE, show_name = FALSE)
}

# Below each column: the x ticks (as matrix_ticks places them on the bin axis)
# and a horizontal colour bar over the column's zmin..zmax.
.bottom_annotation <- function(nbins, ticks, cf, zmin, zmax, th, height) {
  fun <- function(index, k, n) {
    h <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
    tick_h <- 0.38
    pushViewport(viewport(y = unit(1, "npc"), height = unit(tick_h, "in"), just = "top",
                          xscale = c(-0.5, nbins - 0.5)))
    grid.segments(unit(ticks$x, "native"), unit(1, "npc"), unit(ticks$x, "native"),
                  unit(1, "npc") - unit(0.05, "in"), gp = gpar(col = "grey20", lwd = 0.6))
    keep <- nzchar(ticks$label)
    if (any(keep)) {
      grid.text(ticks$label[keep], unit(ticks$x[keep], "native"),
                unit(1, "npc") - unit(0.08, "in"), just = c("centre", "top"),
                gp = epicc_gpar(th, 0.6))
    }
    popViewport()
    bar_h <- 0.12
    zr <- c(zmin, zmax)
    if (!all(is.finite(zr)) || zr[2] <= zr[1]) zr <- range(attr(cf, "breaks"))
    pushViewport(viewport(y = unit(1, "npc") - unit(tick_h + 0.1, "in"), height = unit(bar_h, "in"),
                          width = unit(0.8, "npc"), just = "top", xscale = zr))
    grid.raster(matrix(cf(seq(zr[1], zr[2], length.out = 256)), nrow = 1),
                width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
    grid.rect(gp = gpar(col = "grey30", fill = NA, lwd = 0.5))
    at <- pretty(zr, 3)
    at <- at[at >= zr[1] - 1e-12 & at <= zr[2] + 1e-12]
    grid.xaxis(at = at, gp = epicc_gpar(th, 0.55), main = TRUE)
    popViewport()
  }
  AnnotationFunction(fun = fun, which = "column", height = height,
                     n = 1, subsettable = FALSE, show_name = FALSE)
}

# Build the HeatmapList for prepared data `prep` (prepare_heatmap) and the
# per-column z limits `zlim` (data.frame key, zmin, zmax) and profile limits
# `ylim` (data.frame key, ymin, ymax). Returns list(list, width, height,
# legends) with the page size in inches.
build_heatmap <- function(m, display, th, prep, zlim, ylim, kind = "regions") {
  .check_display(m, display)
  types <- column_types(display)
  strips <- group_strip_labels(m)
  ngroups <- length(strips)
  dashes <- setNames(th$dash_cycle[(seq_len(ngroups) - 1) %% length(th$dash_cycle) + 1],
                     m$groups$label)
  ticks <- profile_ticks(m, kind)
  nbins <- m$samples$nbins
  dpi <- th$heatmap_raster_dpi
  # Rows ComplexHeatmap lays out per group: the binned image rows, so the
  # groups' heights stay in proportion to their region counts.
  group_rows <- vapply(prep$images[[1]], nrow, 0L)
  split <- factor(rep(m$groups$label, group_rows), levels = m$groups$label)
  row_titles <- unname(strips[m$groups$label[group_rows > 0]])
  gap_in <- 0.08
  title_rel <- 0.85
  # Region-group titles run up the side of their slice, unless one of them is
  # longer than its slice is tall (a small nostrand group); then all are
  # written across, beside their slices.
  slice_in <- (th$heatmap_height - gap_in * (sum(group_rows > 0) - 1)) *
    group_rows[group_rows > 0] / max(sum(group_rows), 1)
  title_in <- nchar(row_titles) * th$base_size * 0.8 * 0.55 / 72
  row_rot <- if (all(title_in <= slice_in)) 90 else 0
  row_title_w <- if (row_rot == 90) 0.3 else max(title_in) + 0.15
  ht_list <- NULL
  for (i in seq_len(nrow(m$samples))) {
    key <- m$samples$label[i]
    z <- scale_for_column(zlim, key)
    y <- ylim[match(key, ylim$key), c("ymin", "ymax")]
    cf <- heatmap_colour_fun(epicc_heatmap_colors(th, types[i]), z$zmin, z$zmax,
                             th$heatmap_missing_color)
    mode <- epicc_heatmap_interpolation(th, types[i])
    imgs <- prep$images[[i]]
    nb <- nbins[i]
    layer <- local({
      imgs <- imgs; cf <- cf; mode <- mode; split <- split
      function(j, i, x, y, w, h, fill) {
        g <- as.integer(split[i[1]])
        img <- imgs[[g]]
        if (nrow(img) == 0) return(invisible())
        wpx <- max(1L, as.integer(round(convertWidth(unit(1, "npc"), "in", valueOnly = TRUE) * dpi)))
        hpx <- max(1L, as.integer(round(convertHeight(unit(1, "npc"), "in", valueOnly = TRUE) * dpi)))
        px <- resample_image(img, hpx, wpx, mode)
        grid.raster(matrix(cf(as.vector(px)), nrow = hpx), width = unit(1, "npc"),
                    height = unit(1, "npc"), interpolate = FALSE)
      }
    })
    dk <- prep$d[as.character(prep$d$sample) == key, , drop = FALSE]
    top <- HeatmapAnnotation(
      profile = .profile_annotation(dk, nb, y$ymin, y$ymax, display$colour[i], dashes, th,
                                    matrix_ticks(m, i, ref_label = if (kind == "tes") "end" else "start",
                                                 start_label = "start", end_label = "end"),
                                    unit(th$profile_strip_height, "in")),
      show_annotation_name = FALSE)
    bottom <- HeatmapAnnotation(
      axis = .bottom_annotation(nb, matrix_ticks(m, i, ref_label = if (kind == "tes") "end" else "start",
                                                 start_label = "start", end_label = "end"),
                                cf, z$zmin, z$zmax, th, unit(0.85, "in")),
      show_annotation_name = FALSE)
    ht <- Heatmap(
      matrix(0, length(split), 1), name = key, col = c("white"), rect_gp = gpar(type = "none"),
      layer_fun = layer, cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE,
      row_split = split, row_gap = unit(gap_in, "in"), row_title = row_titles,
      row_title_gp = epicc_gpar(th, 0.8), row_title_rot = row_rot, cluster_row_slices = FALSE,
      column_title = display$label[i], column_title_gp = epicc_gpar(th, title_rel, fontface = "bold"),
      top_annotation = top, bottom_annotation = bottom,
      width = unit(th$heatmap_column_width, "in"), height = unit(th$heatmap_height, "in"),
      border = "grey30", use_raster = FALSE)
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }
  legends <- list()
  if (ngroups > 1) {
    legends <- list(Legend(labels = unname(strips), type = "lines", title = NULL,
                           legend_gp = gpar(lty = unname(dashes), col = "grey20", lwd = 1.2),
                           labels_gp = epicc_gpar(th, 0.75)))
  }
  ncol <- nrow(m$samples)
  col_gap <- 0.55
  legend_w <- if (ngroups > 1) epicc_legend_width(th, strips) else 0
  width <- 0.6 + row_title_w + ncol * th$heatmap_column_width + (ncol - 1) * col_gap + legend_w + 0.3
  height <- 0.6 + th$profile_strip_height + 0.15 + th$heatmap_height + 0.85 + 0.4
  list(list = ht_list, legends = legends, width = width, height = height,
       gap = unit(col_gap, "in"))
}

# Draw the figure to a PDF.
draw_heatmap_pdf <- function(fig, file) {
  grDevices::pdf(file, width = fig$width, height = fig$height, useDingbats = FALSE)
  on.exit(grDevices::dev.off())
  draw(fig$list, ht_gap = fig$gap, annotation_legend_list = fig$legends,
       annotation_legend_side = "right", merge_legends = TRUE,
       padding = unit(c(0.2, 0.3, 0.2, 0.3), "in"))
}
