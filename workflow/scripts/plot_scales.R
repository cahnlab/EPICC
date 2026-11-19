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
# One difference from the bash: its type-mode search started from the sentinels
# +-999999, so a mark whose 10th percentile was >= 999999 got zMin 0. Here the
# mark keeps its range.
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

# heatmaps_scales as the scripts take it: type | sample | default, anything
# else treated as default (as the bash rule did), with a message.
scale_mode <- function(mode) {
  if (mode %in% c("type", "sample", "default")) return(mode)
  message(mode, " unknown. Returning default")
  "default"
}

# The units a mode scales together: column indices per mark (type) or one per
# column (sample). Type mode needs keys of the form {sid}__{track}.
.scale_units <- function(keys, mode) {
  if (mode == "type" && any(!grepl("__", keys, fixed = TRUE))) {
    stop("heatmaps_scales: type needs column keys of the form <sample>__<mark>; got ",
         paste(keys[!grepl("__", keys, fixed = TRUE)], collapse = ", "))
  }
  marks <- column_marks(keys)
  if (mode == "type") split(seq_along(keys), factor(marks, levels = unique(marks)))
  else as.list(seq_along(keys))
}

# y limits per column from average profiles already computed (group_means()
# output: sample, group, x, value), for type | sample: the unit's profile over
# every region group, padded (.pad_range). The y half of matrix_scales(),
# without its percentile pass over the whole matrix, so a plot that already
# has the curves it draws gets the limits of exactly those curves.
profile_ylims <- function(keys, prof, mode = "type") {
  out <- data.frame(key = keys, ymin = NA_real_, ymax = NA_real_, stringsAsFactors = FALSE)
  if (!mode %in% c("type", "sample")) return(out)
  sample <- as.character(prof$sample)
  for (cols in .scale_units(keys, mode)) {
    y <- .pad_range(prof$value[sample %in% keys[cols]])
    out$ymin[cols] <- y[1]; out$ymax[cols] <- y[2]
  }
  out
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
  # Every key the pipeline builds is {sid}__{track} (define_key_for_plots);
  # without the '__' there is no mark to scale by (.scale_units stops).
  units <- .scale_units(keys, mode)
  rng <- column_ranges(m)
  y <- profile_ylims(keys, group_means(m, average), mode)
  out$ymin <- y$ymin; out$ymax <- y$ymax
  for (cols in units) {
    lo <- rng[cols, 1]; hi <- rng[cols, 2]
    # The bash searched a mark's range from the sentinels 999999 / -999999,
    # so in type mode a 10th percentile >= 999999 (or a 90th <= -999999) was
    # never taken and the mark fell back to 0 on that side. That was a bug,
    # not a choice (raw counts reach it); it is not reproduced.
    z <- c(if (all(is.na(lo))) 0 else min(lo, na.rm = TRUE),
           if (all(is.na(hi))) 0 else max(hi, na.rm = TRUE))
    if (z[1] == 0 && z[2] == 0) z <- .z_fallback(marks[cols[1]], mc_zmax)
    out$zmin[cols] <- z[1]; out$zmax[cols] <- z[2]
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
  z <- default_zlim(m$values)
  y <- default_ylim(group_means(m, average))
  list(zmin = z[1], zmax = z[2], ymin = y[1], ymax = y[2])
}

# The z half of default_scales() for the values `v` of the whole matrix.
default_zlim <- function(v) {
  z <- if (length(v) == 0 || anyNA(v)) c(NA_real_, NA_real_)
       else unname(stats::quantile(v, c(0.01, 0.98), type = 7, names = FALSE))
  if (anyNA(z) || z[2] <= z[1]) {
    z <- if (all(is.na(v))) c(0, 1) else range(v, na.rm = TRUE)
  }
  z
}

# default_zlim() over values held column by column (`get(j)` returns column
# j, j = 1..ncols; read_matrix_table's table), without gathering them into one
# vector: the same numbers, in memory proportional to one column. The 1st and
# 98th percentiles (R's type 7, as quantile() and numpy compute them) are
# found by counting the values into 2^16 buckets of their range, then sorting
# only the buckets that hold the wanted ranks.
default_zlim_columns <- function(get, ncols) {
  n <- 0; lo <- Inf; hi <- -Inf
  for (j in seq_len(ncols)) {
    v <- get(j)
    if (anyNA(v)) return(.range_columns(get, ncols))   # numpy's percentile gives NaN
    n <- n + length(v)
    if (length(v)) { lo <- min(lo, v); hi <- max(hi, v) }
  }
  if (n == 0) return(c(0, 1))
  index <- 1 + (n - 1) * c(0.01, 0.98)
  ranks <- sort(unique(c(floor(index), ceiling(index))))
  value_at <- if (hi == lo) setNames(rep(lo, length(ranks)), ranks) else {
    nb <- 65536
    breaks <- seq(lo, hi, length.out = nb + 1)
    counts <- numeric(nb)
    for (j in seq_len(ncols)) {
      counts <- counts + tabulate(findInterval(get(j), breaks, rightmost.closed = TRUE,
                                               all.inside = TRUE), nb)
    }
    cum <- cumsum(counts)
    bucket <- findInterval(ranks - 1, cum) + 1           # bucket holding each rank
    want <- unique(bucket)
    held <- unlist(lapply(seq_len(ncols), function(j) {
      v <- get(j)
      v[findInterval(v, breaks, rightmost.closed = TRUE, all.inside = TRUE) %in% want]
    }), use.names = FALSE)
    held <- sort(held)
    before <- c(0, cum)[want]                            # values below each wanted bucket
    first <- cumsum(c(0, counts[want]))[seq_along(want)] # where each bucket starts in `held`
    setNames(vapply(seq_along(ranks), function(k) {
      w <- match(bucket[k], want)
      held[first[w] + ranks[k] - before[w]]
    }, 0), ranks)
  }
  q <- vapply(index, function(i) {
    a <- value_at[[as.character(floor(i))]]
    b <- value_at[[as.character(ceiling(i))]]
    h <- i - floor(i)
    if (i > floor(i) && b != a) (1 - h) * a + h * b else a
  }, 0)
  if (q[2] <= q[1]) q <- c(lo, hi)
  q
}

.range_columns <- function(get, ncols) {
  lo <- Inf; hi <- -Inf
  for (j in seq_len(ncols)) {
    v <- get(j)
    v <- v[!is.na(v)]
    if (length(v)) { lo <- min(lo, v); hi <- max(hi, v) }
  }
  if (!is.finite(lo)) c(0, 1) else c(lo, hi)
}

# The shared y range of default_scales() from average profiles already
# computed: their range with matplotlib's 5% margins.
default_ylim <- function(prof) {
  p <- prof$value[!is.na(prof$value)]
  y <- if (length(p) == 0) c(0, 1) else range(p)
  d <- if (y[2] > y[1]) (y[2] - y[1]) * 0.05 else max(abs(y[2]) * 0.05, 0.01)
  c(y[1] - d, y[2] + d)
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
  y <- function(x) paste(vapply(x, .awk_number, ""), collapse = " ")
  c(heatmap = paste0("--zMin ", z(scales$zmin), " --zMax ", z(scales$zmax)),
    profile = paste0("--yMin ", y(scales$ymin), " --yMax ", y(scales$ymax)))
}

# A number as awk prints it (the old rule wrote the y bounds with awk's
# print): an integral value in full, anything else with OFMT %.6g.
.awk_number <- function(x) {
  if (x == round(x) && abs(x) < 2^53) return(formatC(x, format = "d", big.mark = ""))
  trimws(formatC(x, format = "g", digits = 6))
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
