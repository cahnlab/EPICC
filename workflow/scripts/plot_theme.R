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
# data types the defaults don't name.

suppressPackageStartupMessages({
  library(ggplot2)
})

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
    sRNA = .EPICC_SEISMIC,
    mC = .EPICC_ORANGES
  ),
  heatmap_missing_color = "grey90",
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

# Keys whose value is an open name -> value map rather than a fixed record.
.EPICC_OPEN_MAPS <- c("heatmap_colors")

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
        if (!is.character(value[[entry]]) || length(value[[entry]]) < 2) {
          stop(at, ".", entry, " must be a list of at least two colours")
        }
        base[[key]][[entry]] <- value[[entry]]
      }
    } else if (is.list(base[[key]])) {
      base[[key]] <- epicc_merge_theme(base[[key]], value, at)
    } else {
      if (!.epicc_same_type(base[[key]], value)) {
        stop(at, " must be ", class(base[[key]])[1], ", got ", class(value)[1])
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

# Colour ramp for a heatmap of one data type: the type's own entry, else its
# family (mCG/mCHG/mCHH -> mC, ChIP_broad -> ChIP), else 'default' (e.g. a
# mark name such as H3K9me2 that has no entry of its own).
epicc_heatmap_colors <- function(th, type) {
  maps <- th$heatmap_colors
  if (!is.null(maps[[type]])) return(maps[[type]])
  if (grepl("^mC", type) && !is.null(maps[["mC"]])) return(maps[["mC"]])
  family <- sub("_.*$", "", type)
  if (!is.null(maps[[family]])) return(maps[[family]])
  maps[["default"]]
}

# Write a ggplot to PDF at the given size in inches.
epicc_save <- function(plot, file, width, height) {
  ggsave(file, plot, width = width, height = height, limitsize = FALSE)
}
