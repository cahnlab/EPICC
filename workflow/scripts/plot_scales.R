# Heatmap colour scales and profile y-limits per matrix column, in R.
#
# This is computing_matrix_scales (combined_analysis.smk) ported from bash,
# computeMatrixOperations dataRange and plotProfile --outFileNameData, so the
# R heatmap and profile scripts can work out their limits from the matrix they
# draw, with the same averaging code (matrix_io.R's group_means) they draw it
# with. Sourced after matrix_io.R.
#
#   s <- matrix_scales(m, mode = "type", average = "mean",
#                      mc_zmax = c(mCG = 100, mCHG = 50, mCHH = 20))
#   s$zmin[i], s$zmax[i], s$ymin[i], s$ymax[i]   # limits of matrix column i
#
# One row per matrix column, in column order -- never one per mark. An
# unstranded matrix of a stranded track interleaves each sample's columns
# (s1__X_plus, s1__X_minus, s2__X_plus, ...), so scales must be looked up by
# column (scale_for_column), not handed out in mark blocks.
#
# What the bash rule does, and this reproduces:
#   heatmaps_scales: type    one range per mark; a column's mark is the last
#                            '__' field of its key, matched exactly (H3K9me
#                            does not take H3K9me2's columns)
#   heatmaps_scales: sample  one range per column
#   heatmaps_scales: default no per-column limits: deeptools then uses one
#                            scale for the whole plot (default_scales below)
# zMin/zMax are dataRange's 10th/90th percentiles of a column's values; a
# column holding any NaN gets NaN there (numpy's percentile propagates it), so
# it counts as having no range -- which is why mC columns, which always carry
# NaNs, end up on the heat_mcg/heat_mchg/heat_mchh fallback. A mark (or
# column) whose range is 0..0 or missing gets 0..heat_m<ctx> for mCG/mCHG/mCHH
# and 0..0.005 for anything else.
# y limits span the column's (or mark's) average profile over every region
# group, padded by a tenth of the range (or of |max| when the range is 0, or
# 0.01 when both are 0). A profile with no value at all gets 0..0.01.

# The mark (track) of each column key: its last '__' field.
column_marks <- function(keys) {
  vapply(strsplit(keys, "__", fixed = TRUE), function(f) f[length(f)], "")
}

# numpy's percentile (method "linear", R's quantile type 7), computed with
# numpy's own interpolation formula so the result is the same double, not
# merely the same to the last bit or two: a + (b - a) * t, or b - (b - a) *
# (1 - t) when t >= 0.5.
np_percentile <- function(v, q) {
  v <- sort(v)
  n <- length(v)
  vapply(q, function(qq) {
    h <- (n - 1) * (qq / 100)
    lo <- floor(h)
    t <- h - lo
    a <- v[lo + 1]
    b <- v[min(lo + 2, n)]
    d <- b - a
    if (t >= 0.5) b - d * (1 - t) else a + d * t
  }, 0)
}

# dataRange's 10th and 90th percentiles of each column's values (NA when the
# column holds a NaN, as numpy's percentile propagates it).
column_ranges <- function(m) {
  t(vapply(seq_len(nrow(m$samples)), function(i) {
    v <- m$values[, sample_columns(m, i)]
    if (length(v) == 0 || anyNA(v)) return(c(NA_real_, NA_real_))
    np_percentile(as.vector(v), c(10, 90))
  }, c(0, 0)))
}

# The bash rule's padding of an average profile's range.
.pad_range <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(c(0, 0.01))
  lo <- min(v); hi <- max(v)
  pad <- (hi - lo) * 0.1
  if (pad <= 0) pad <- abs(hi) * 0.1
  if (pad <= 0) pad <- 0.01
  c(lo - pad, hi + pad)
}

.z_fallback <- function(mark, mc_zmax) {
  if (mark %in% names(mc_zmax)) c(0, unname(mc_zmax[[mark]])) else c(0, 0.005)
}

# Per-column limits. `mode`: type | sample | default (anything else is
# treated as default, as the bash rule does). `average`: profiles_scale
# (mean | median). `mc_zmax`: heat_mcg/heat_mchg/heat_mchh by context name.
# Returns data.frame(key, mark, zmin, zmax, ymin, ymax); all NA for default.
matrix_scales <- function(m, mode = "type", average = "mean",
                          mc_zmax = c(mCG = 100, mCHG = 50, mCHH = 20)) {
  keys <- m$samples$label
  marks <- column_marks(keys)
  out <- data.frame(key = keys, mark = marks, zmin = NA_real_, zmax = NA_real_,
                    ymin = NA_real_, ymax = NA_real_, stringsAsFactors = FALSE)
  if (!mode %in% c("type", "sample")) return(out)
  rng <- column_ranges(m)
  prof <- group_means(m, average)
  prof_of <- function(cols) prof$value[as.character(prof$sample) %in% keys[cols]]
  # The bash matched a mark only in keys that have a '__' (every key the
  # pipeline builds does); in type mode a key without one found no data.
  if (mode == "type") {
    bare <- !grepl("__", keys, fixed = TRUE)
    rng[bare, ] <- NA_real_
    prof$value[as.character(prof$sample) %in% keys[bare]] <- NA_real_
  }
  units <- if (mode == "type") split(seq_along(keys), factor(marks, levels = unique(marks)))
           else as.list(seq_along(keys))
  for (cols in units) {
    lo <- rng[cols, 1]; hi <- rng[cols, 2]
    z <- c(if (all(is.na(lo))) 0 else min(lo, na.rm = TRUE),
           if (all(is.na(hi))) 0 else max(hi, na.rm = TRUE))
    if (z[1] == 0 && z[2] == 0) z <- .z_fallback(marks[cols[1]], mc_zmax)
    y <- .pad_range(prof_of(cols))
    out$zmin[cols] <- z[1]; out$zmax[cols] <- z[2]
    out$ymin[cols] <- y[1]; out$ymax[cols] <- y[2]
  }
  out
}

# The single scale deeptools falls back to when given no limits
# (heatmaps_scales: default): plotHeatmap colours everything between the 1st
# and 98th percentile of the whole matrix (NaN-propagating; when that gives
# NaN, or max <= min, matplotlib autoscales each image to its own data, here
# approximated by the matrix's NaN-aware min/max), and plotProfile shares one
# y axis across panels, autoscaled with matplotlib's 5% margins.
default_scales <- function(m, average = "mean") {
  v <- m$values
  z <- if (length(v) == 0 || anyNA(v)) c(NA_real_, NA_real_)
       else unname(stats::quantile(v, c(0.01, 0.98), type = 7, names = FALSE))
  if (anyNA(z) || z[2] <= z[1]) {
    z <- if (all(is.na(v))) c(0, 1) else range(v, na.rm = TRUE)
  }
  p <- group_means(m, average)$value
  p <- p[!is.na(p)]
  y <- if (length(p) == 0) c(0, 1) else range(p)
  d <- if (y[2] > y[1]) (y[2] - y[1]) * 0.05 else max(abs(y[2]) * 0.05, 0.01)
  list(zmin = z[1], zmax = z[2], ymin = y[1] - d, ymax = y[2] + d)
}

# The limits of one column, by key or index. Every R plot looks scales up
# through this, by its own column key, never by position within a mark.
scale_for_column <- function(scales, column) {
  i <- if (is.character(column)) match(column, scales$key) else column
  if (is.na(i) || i < 1 || i > nrow(scales)) stop("No scale for column ", column)
  scales[i, c("zmin", "zmax", "ymin", "ymax")]
}

# One y range for a plot that draws every column on one axis
# (Profile_pergroup__): the lowest ymin and highest ymax, as the profile rule
# collapses the per-column values today.
shared_ylim <- function(scales) {
  if (all(is.na(scales$ymin))) return(c(NA_real_, NA_real_))
  c(min(scales$ymin, na.rm = TRUE), max(scales$ymax, na.rm = TRUE))
}

# The deeptools flag strings computing_matrix_scales writes (params_heatmap /
# params_profile); empty for default. Formatted as the rule formats them: z at
# full precision (it passes dataRange's numbers through), y as awk prints
# numbers (%.6g).
scale_flags <- function(scales) {
  if (all(is.na(scales$zmin))) return(c(heatmap = "", profile = ""))
  z <- function(x) paste(vapply(x, .shortest, ""), collapse = " ")
  y <- function(x) paste(trimws(formatC(x, format = "g", digits = 6)), collapse = " ")
  c(heatmap = paste0("--zMin ", z(scales$zmin), " --zMax ", z(scales$zmax)),
    profile = paste0("--yMin ", y(scales$ymin), " --yMax ", y(scales$ymax)))
}

# A number with as few digits as identify its double (at most 17), in plain
# notation for the magnitudes scales take, as dataRange printed them. The text
# can still differ from Python's (R's own string-to-number parser is not
# always correctly rounded at 16-17 digits, so R may keep a digit Python
# drops), but deeptools and R read both back as the same double.
.shortest <- function(x) {
  if (x == 0) return("0")
  for (d in 1:17) {
    s <- trimws(formatC(x, format = "g", digits = d))
    if (as.numeric(s) == x) break
  }
  e <- floor(log10(abs(x)))
  if (e >= -4 && e < 16) {
    s <- formatC(x, format = "f", digits = max(d - 1 - e, 0))
    if (grepl(".", s, fixed = TRUE)) s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  trimws(s)
}
