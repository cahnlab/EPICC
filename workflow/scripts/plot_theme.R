# Shared plot theme for EPICC's R plots.
#
# Every plotting script sources this file and takes its look from one theme
# object: fonts and sizes, line widths, facet strips, legend layout, the
# categorical palette, heatmap colour scales per data type, and PDF panel
# dimensions. Built-in themes are defined below; the options YAML picks one
# (`plot_theme`) and may replace individual keys (`plot_theme_overrides`).
# Snakemake hands the choice to a script as a small JSON file:
#
#   {"theme": "default", "overrides": {"base_size": 12,
#                                       "heatmap_colors": {"RNA": ["white", "darkgreen"]}}}
#
# Loading, from a script under workflow/scripts/:
#
#   script_dir <- dirname(normalizePath(sub("^--file=", "",
#     grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
#   source(file.path(script_dir, "plot_theme.R"))
#   th <- epicc_theme(theme_file)        # theme_file may be NULL or ""
#   p + epicc_ggtheme(th)
#
# Sample colours are not part of the theme: they arrive with each plot's
# display table (tab20 keyed on levels_label, assigned in combined_analysis.smk)
# so a sample keeps its colour across every plot. `palette` is the fallback for
# plots that have no display table.
#
# Override rules: a key the theme doesn't define is an error (a typo must not
# be silently ignored), and a value must have the type of the default it
# replaces. `heatmap_colors` is the one open map -- it may gain entries for
# data types the defaults don't name, and for single marks (an IP_target such
# as H3K9me2, RAMPAGE, sRNA_24nt), which win over their data type. A
# `heatmap_colors` entry is either a list of at least two colours or one
# matplotlib colour-map name ("viridis", "RdBu_r"), looked up in
# plot_colormaps.json beside this file.

suppressPackageStartupMessages({
  library(ggplot2)
})

# This file's directory, for plot_colormaps.json. source() records the path
# it was called with; Rscript scripts source this file by absolute path.
.EPICC_THEME_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                             error = function(e) NULL)

# matplotlib colour maps deeptools used, as hex, so the R heatmaps keep the
# look of the deeptools ones until a theme changes it.
.EPICC_SEISMIC <- c("#00004C", "#0000FF", "#FFFFFF", "#FF0000", "#7F0000")
.EPICC_ORANGES <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C",
                    "#F16913", "#D94801", "#A63603", "#7F2704")
# matplotlib tab20, the same list as _PALETTES["tab20"] in combined_analysis.smk.
.EPICC_TAB20 <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c",
                  "#98df8a", "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                  "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", "#7f7f7f",
                  "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5")

.epicc_default_theme <- list(
  # Text
  base_size = 11,            # pt; every other text size derives from it
  base_family = "",          # "" = the PDF device default (Helvetica)
  # Lines
  line_width = 0.6,          # mm (ggplot2 linewidth) for profile lines
  profile_type = "lines",    # lines | fill | se | std (deeptools --plotType)
  dash_cycle = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
  # Panels and facets
  grid_minor = FALSE,
  panel_spacing = 1.4,       # lines, between facet columns
  strip_fill = "grey92",
  strip_face = "bold",
  # Legend
  legend_position = "right",
  legend_key_width = 2.2,    # lines; wide enough to tell dash patterns apart
  legend_max_rows = 25,      # entries per legend column before wrapping
  legend_row_height = 0.22,  # in; height one legend entry needs
  legend_char_width = 0.09,  # in; width one legend character needs
  # Colours
  palette = .EPICC_TAB20,    # fallback categorical palette
  heatmap_colors = list(
    default = .EPICC_SEISMIC,
    ChIP = .EPICC_SEISMIC,
    TF = .EPICC_SEISMIC,
    ATAC = .EPICC_SEISMIC,
    RNA = .EPICC_SEISMIC,
    # TSS-capture counts, on a range of their own: a ramp of their own too,
    # so they don't read as RNA-seq next to it.
    RNA_RAMPAGE = "PRGn",
    sRNA = .EPICC_SEISMIC,
    mC = .EPICC_ORANGES
  ),
  heatmap_missing_color = "grey90",
  # How a heatmap column is resampled to its pixels, per data type (same
  # lookup as heatmap_colors): bilinear smooths between bins (continuous
  # tracks), nearest draws each bin as a block (mC's coarse, often-missing
  # bins). Done in R before drawing (plot_heatmap.R), so every viewer shows
  # the same image.
  heatmap_interpolation = list(
    default = "bilinear",
    ChIP = "bilinear",
    TF = "bilinear",
    ATAC = "bilinear",
    RNA = "bilinear",
    sRNA = "bilinear",
    mC = "nearest"
  ),
  heatmap_raster_dpi = 150,   # pixels per inch of the drawn heatmap images
  # PDF dimensions, inches
  panel_width = 3.2,         # one profile/facet panel
  panel_height = 3.2,
  heatmap_column_width = 1.6, # one sample column of a heatmap
  heatmap_height = 9,         # heatmap body, all region groups together
  profile_strip_height = 1.2  # average profile drawn above a heatmap column
)

EPICC_THEMES <- list(
  default = .epicc_default_theme,
  # Larger text and heavier lines for slides and posters.
  presentation = modifyList(.epicc_default_theme, list(
    base_size = 15,
    line_width = 1.0,
    legend_row_height = 0.3,
    legend_char_width = 0.12,
    panel_width = 4.2,
    panel_height = 3.8
  ))
)

# Keys whose value is an open map from data type to a value rather than a
# fixed record: they may gain entries for types the defaults don't name.
.EPICC_OPEN_MAPS <- c("heatmap_colors", "heatmap_interpolation")
.EPICC_INTERPOLATIONS <- c("bilinear", "nearest")

# Keys that take one of a fixed set of values.
.EPICC_CHOICES <- list(
  profile_type = c("lines", "fill", "se", "std"),
  legend_position = c("right", "bottom", "top", "left", "none")
)

.epicc_colormap_table <- NULL

# matplotlib colour maps by name, as hex stops (plot_colormaps.json).
epicc_colormaps <- function() {
  if (is.null(.epicc_colormap_table)) {
    path <- if (!is.null(.EPICC_THEME_DIR)) file.path(.EPICC_THEME_DIR, "plot_colormaps.json")
    if (is.null(path) || !file.exists(path)) {
      stop("plot_colormaps.json not found next to plot_theme.R")
    }
    .epicc_colormap_table <<- jsonlite::fromJSON(path)$maps
  }
  .epicc_colormap_table
}

# Hex stops of a named colour map; a `_r` suffix reverses it, as in matplotlib.
epicc_colormap <- function(name) {
  maps <- epicc_colormaps()
  base <- sub("_r$", "", name)
  if (!is.null(maps[[name]])) return(maps[[name]])
  if (base != name && !is.null(maps[[base]])) return(rev(maps[[base]]))
  stop("Unknown colour map '", name, "'. Known maps: ",
       paste(names(maps), collapse = ", "), " (add _r to reverse)")
}

.epicc_same_type <- function(default, value) {
  if (is.numeric(default)) return(is.numeric(value))
  if (is.logical(default)) return(is.logical(value))
  if (is.character(default)) return(is.character(value))
  if (is.list(default)) return(is.list(value))
  TRUE
}

# Replace keys of `base` with those of `overrides`, refusing unknown keys and
# type changes. `where` names the position for error messages.
epicc_merge_theme <- function(base, overrides, where = "plot_theme_overrides") {
  if (length(overrides) == 0) return(base)
  if (!is.list(overrides) || is.null(names(overrides)) || any(names(overrides) == "")) {
    stop(where, " must be a mapping of theme keys to values")
  }
  unknown <- setdiff(names(overrides), names(base))
  if (length(unknown) > 0) {
    stop("Unknown theme key(s) in ", where, ": ", paste(unknown, collapse = ", "),
         ". Known keys: ", paste(names(base), collapse = ", "))
  }
  for (key in names(overrides)) {
    value <- overrides[[key]]
    at <- paste0(where, ".", key)
    if (key %in% .EPICC_OPEN_MAPS) {
      if (!is.list(value) || is.null(names(value))) stop(at, " must be a mapping")
      for (entry in names(value)) {
        v <- value[[entry]]
        if (key == "heatmap_interpolation") {
          if (!is.character(v) || length(v) != 1 || !v %in% .EPICC_INTERPOLATIONS) {
            stop(at, ".", entry, " must be one of ", paste(.EPICC_INTERPOLATIONS, collapse = ", "))
          }
        } else {
          if (!is.character(v) || length(v) < 1) {
            stop(at, ".", entry, " must be a colour-map name or a list of at least two colours")
          }
          if (length(v) == 1) epicc_colormap(v)   # errors on an unknown name
        }
        base[[key]][[entry]] <- v
      }
    } else if (is.list(base[[key]])) {
      base[[key]] <- epicc_merge_theme(base[[key]], value, at)
    } else {
      if (!.epicc_same_type(base[[key]], value)) {
        stop(at, " must be ", class(base[[key]])[1], ", got ", class(value)[1])
      }
      allowed <- .EPICC_CHOICES[[key]]
      if (!is.null(allowed) && !(length(value) == 1 && value %in% allowed)) {
        stop(at, " must be one of ", paste(allowed, collapse = ", "), ", got ",
             paste(value, collapse = ", "))
      }
      base[[key]] <- value
    }
  }
  base
}

# The theme a script should use. `theme_file` is the JSON Snakemake wrote
# ({"theme": name, "overrides": {...}}); NULL or "" gives the default theme.
epicc_theme <- function(theme_file = NULL) {
  name <- "default"
  overrides <- list()
  if (!is.null(theme_file) && nzchar(theme_file)) {
    if (!file.exists(theme_file)) stop("Theme file not found: ", theme_file)
    spec <- jsonlite::fromJSON(theme_file, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                               simplifyMatrix = FALSE)
    if (!is.null(spec$theme) && nzchar(spec$theme)) name <- spec$theme
    if (!is.null(spec$overrides)) overrides <- spec$overrides
  }
  if (!name %in% names(EPICC_THEMES)) {
    stop("Unknown plot_theme '", name, "'. Built-in themes: ",
         paste(names(EPICC_THEMES), collapse = ", "))
  }
  epicc_merge_theme(EPICC_THEMES[[name]], overrides)
}

# The ggplot2 theme for a theme object.
epicc_ggtheme <- function(th) {
  g <- theme_bw(base_size = th$base_size, base_family = th$base_family) +
    theme(panel.spacing.x = unit(th$panel_spacing, "lines"),
          strip.background = element_rect(fill = th$strip_fill, colour = NA),
          strip.text = element_text(face = th$strip_face),
          legend.position = th$legend_position,
          legend.key.width = unit(th$legend_key_width, "lines"))
  if (!isTRUE(th$grid_minor)) g <- g + theme(panel.grid.minor = element_blank())
  g
}

# grid graphical parameters at the theme's font, for ComplexHeatmap and other
# grid-drawn plots that can't take a ggplot theme. `rel` scales the size.
epicc_gpar <- function(th, rel = 1, ...) {
  grid::gpar(fontsize = th$base_size * rel, fontfamily = th$base_family, ...)
}

# Number of legend columns for n entries.
epicc_legend_ncol <- function(th, n) {
  max(1L, as.integer(ceiling(n / th$legend_max_rows)))
}

# Width (in) a legend of these labels needs, at `ncol` columns.
epicc_legend_width <- function(th, labels, ncol = epicc_legend_ncol(th, length(labels))) {
  1.5 + th$legend_char_width * max(nchar(labels), 0) * ncol
}

# Height (in) a legend of n entries needs, before wrapping into more columns.
epicc_legend_height <- function(th, n) {
  th$legend_row_height * min(n, th$legend_max_rows) + 1
}

# Dash pattern per series. Colours are keyed on levels_label, so series that
# differ only in assay or replicate share a colour; those get distinct dash
# patterns, in series order, so the lines stay apart. A series with a colour
# to itself stays solid. `colours` is named by series.
epicc_dashes <- function(th, colours) {
  rank <- ave(seq_along(colours), unname(colours), FUN = seq_along)
  cycle <- th$dash_cycle
  setNames(cycle[(rank - 1) %% length(cycle) + 1], names(colours))
}

# Colour ramp (hex stops) for a heatmap column of one data type and mark: the
# mark's own entry (an IP_target such as H3K9me2, RAMPAGE, sRNA_24nt), else
# the type's, else its family's (mCG/mCHG/mCHH -> mC, RNA_RAMPAGE -> RNA),
# else 'default'. A colour-map name is expanded to its stops.
epicc_heatmap_colors <- function(th, type, mark = NULL) {
  v <- .epicc_type_entry(th$heatmap_colors, type, mark)
  if (length(v) == 1) epicc_colormap(v) else v
}

# bilinear | nearest for a heatmap column of one data type and mark, with the
# same fallbacks as epicc_heatmap_colors().
epicc_heatmap_interpolation <- function(th, type, mark = NULL) {
  .epicc_type_entry(th$heatmap_interpolation, type, mark)
}

# A per-type map's entry: the mark's own, else the type's, else its family
# (mCG/mCHG/mCHH -> mC, RNA_RAMPAGE -> RNA), else 'default'.
.epicc_type_entry <- function(maps, type, mark = NULL) {
  family <- sub("_.*$", "", type)
  if (length(mark) == 1 && !is.na(mark) && nzchar(mark) && !is.null(maps[[mark]])) maps[[mark]]
  else if (!is.null(maps[[type]])) maps[[type]]
  else if (grepl("^mC", type) && !is.null(maps[["mC"]])) maps[["mC"]]
  else if (!is.null(maps[[family]])) maps[[family]]
  else maps[["default"]]
}

# Write a ggplot to PDF at the given size in inches.
epicc_save <- function(plot, file, width, height) {
  ggsave(file, plot, width = width, height = height, limitsize = FALSE)
}
