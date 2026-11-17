# Read deeptools computeMatrix files in R.
#
# The R plotting scripts read the matrix that computeMatrix (and
# computeMatrixOperations) wrote, rather than text dumps from deeptools'
# plotters. A matrix file is gzipped: line 1 is '@' followed by a JSON header,
# then one row per region: six BED columns and one value per bin per sample,
# 'nan' where there was no signal.
#
#   source(file.path(script_dir, "matrix_io.R"))
#   m <- read_matrix("final_matrix.gz")
#   m <- sort_rows(m, "descend", "mean")      # deeptools' --sortRegions/--sortUsing
#   prof <- group_means(m, "mean")            # what plotProfile draws
#
# Coordinates follow deeptools: a sample's bins are numbered x = 0 .. nbins-1,
# and matrix_ticks() returns tick positions on that axis, so ticks, lines and
# heatmap cells line up the way plotProfile/plotHeatmap draw them.
#
# Summaries are NaN-aware (deeptools uses np.nanmean, np.nanmedian, ... and
# masked arrays): a bin or row with no value anywhere is NA, everything else
# ignores the missing entries.

suppressPackageStartupMessages({
  library(data.table)
})

# The matrix's JSON header as a list. Per-sample entries (upstream,
# downstream, body, bin size, ref point, unscaled 5/3 prime) are vectors with
# one element per sample.
read_matrix_header <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  first <- readLines(con, n = 1)
  if (length(first) == 0 || !startsWith(first, "@")) {
    stop(path, ": not a deeptools matrix (line 1 does not start with '@')")
  }
  h <- jsonlite::fromJSON(sub("^@", "", first), simplifyVector = TRUE)
  for (key in c("group_labels", "group_boundaries", "sample_labels", "sample_boundaries")) {
    if (is.null(h[[key]])) stop(path, ": matrix header has no '", key, "'")
  }
  if (length(h$group_boundaries) != length(h$group_labels) + 1) {
    stop(path, ": ", length(h$group_labels), " group labels but ",
         length(h$group_boundaries), " group boundaries")
  }
  if (length(h$sample_boundaries) != length(h$sample_labels) + 1) {
    stop(path, ": ", length(h$sample_labels), " sample labels but ",
         length(h$sample_boundaries), " sample boundaries")
  }
  h
}

# Per-sample value of a header field: deeptools 3+ stores a list with one
# entry per sample, older files a scalar.
.matrix_param <- function(h, key, i, default = 0) {
  v <- h[[key]]
  if (is.null(v) || length(v) == 0) return(default)
  v <- if (length(v) >= i) v[[i]] else v[[1]]
  if (is.null(v) || (length(v) == 1 && is.na(v))) default else v
}

# Read a matrix file. Returns a list:
#   header   the JSON header (read_matrix_header)
#   regions  data.frame: chrom, start, end, name, score, strand, group
#            (factor, levels in header order)
#   values   numeric matrix, one row per region, one column per bin of every
#            sample (samples side by side, in header order); NA = no value
#   samples  data.frame: label, first, last (column range in `values`),
#            nbins, and the bin layout (upstream, downstream, body, bin_size,
#            unscaled5, unscaled3, ref_point)
#   groups   data.frame: label, first, last (row range), n
read_matrix <- function(path) {
  h <- read_matrix_header(path)
  nrows <- tail(h$group_boundaries, 1)
  ncols <- tail(h$sample_boundaries, 1)
  if (nrows == 0) {
    body <- data.table(matrix(numeric(0), 0, 6 + ncols))
  } else {
    body <- fread(cmd = paste("gzip -dc", shQuote(path)), skip = 1, header = FALSE,
                  sep = "\t", na.strings = c("nan", "NaN", "-nan", "inf", "-inf"),
                  colClasses = c(rep("character", 6), rep("numeric", ncols)),
                  showProgress = FALSE)
  }
  if (nrow(body) != nrows) {
    stop(path, ": header announces ", nrows, " regions, file has ", nrow(body))
  }
  if (ncol(body) != 6 + ncols) {
    stop(path, ": header announces ", ncols, " value columns, rows have ", ncol(body) - 6)
  }
  gb <- h$group_boundaries
  sb <- h$sample_boundaries
  groups <- data.frame(label = h$group_labels, first = head(gb, -1) + 1,
                       last = tail(gb, -1), stringsAsFactors = FALSE)
  groups$n <- groups$last - groups$first + 1
  regions <- data.frame(
    chrom = body[[1]], start = as.integer(body[[2]]), end = as.integer(body[[3]]),
    name = body[[4]], score = body[[5]], strand = body[[6]],
    group = factor(rep(h$group_labels, groups$n), levels = unique(h$group_labels)),
    stringsAsFactors = FALSE)
  values <- as.matrix(body[, (7):(6 + ncols)])
  dimnames(values) <- NULL
  values[!is.finite(values)] <- NA_real_
  ns <- length(h$sample_labels)
  samples <- data.frame(
    label = h$sample_labels, first = head(sb, -1) + 1, last = tail(sb, -1),
    stringsAsFactors = FALSE)
  samples$nbins <- samples$last - samples$first + 1
  num <- function(key) vapply(seq_len(ns), function(i) as.numeric(.matrix_param(h, key, i)), 0)
  samples$upstream <- num("upstream")
  samples$downstream <- num("downstream")
  samples$body <- num("body")
  samples$bin_size <- num("bin size")
  samples$unscaled5 <- num("unscaled 5 prime")
  samples$unscaled3 <- num("unscaled 3 prime")
  samples$ref_point <- vapply(seq_len(ns), function(i) {
    v <- .matrix_param(h, "ref point", i, default = NA_character_)
    if (is.null(v)) NA_character_ else as.character(v)
  }, "")
  list(header = h, regions = regions, values = values, samples = samples, groups = groups)
}

# Column indices in m$values for sample i (index or label).
sample_columns <- function(m, i) {
  if (is.character(i)) {
    j <- match(i, m$samples$label)
    if (is.na(j)) stop("No sample '", i, "' in matrix")
    i <- j
  }
  m$samples$first[i]:m$samples$last[i]
}

.na_summary <- function(x, f) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else f(x)
}

.summary_fun <- function(using) {
  switch(using,
         mean = mean, median = stats::median, max = max, min = min, sum = sum,
         stop("Unsupported summary '", using, "' (mean, median, max, min, sum)"))
}

# NaN-aware summary of each row across the given columns.
row_summary <- function(values, using = "mean") {
  if (nrow(values) == 0) return(numeric(0))
  if (using == "mean") {
    n <- rowSums(!is.na(values))
    s <- rowSums(values, na.rm = TRUE)
    return(ifelse(n == 0, NA_real_, s / n))
  }
  if (using == "sum") {
    # np.nansum gives 0 for an all-NaN row
    return(rowSums(values, na.rm = TRUE))
  }
  f <- .summary_fun(using)
  apply(values, 1, .na_summary, f = f)
}

# NaN-aware summary of each column across the given rows.
col_summary <- function(values, using = "mean") {
  if (nrow(values) == 0) return(rep(NA_real_, ncol(values)))
  if (using == "mean") {
    n <- colSums(!is.na(values))
    s <- colSums(values, na.rm = TRUE)
    return(ifelse(n == 0, NA_real_, s / n))
  }
  f <- .summary_fun(using)
  apply(values, 2, .na_summary, f = f)
}

# Reorder the regions within each group, as plotHeatmap --sortRegions
# {keep,ascend,descend} --sortUsing {mean,median,max,min,sum} does: rank each
# row by the summary of all its values (every sample, or only `samples`, as
# --sortUsingSamples), NaN-aware, and order rows within their group.
#
# deeptools takes numpy's argsort (ascending, NaN last) and reverses it for
# 'descend', so a row with no value at all goes to the *top* of its group
# under descend and to the bottom under ascend, and tied rows come out in
# reverse order under descend. That is reproduced here with a stable sort; for
# groups larger than 16 rows numpy's quicksort is not stable, so the order of
# exactly tied rows can differ from deeptools there.
sort_rows <- function(m, method = "descend", using = "mean", samples = NULL) {
  method <- match.arg(method, c("keep", "no", "descend", "ascend"))
  if (method %in% c("keep", "no")) return(m)
  cols <- if (is.null(samples)) seq_len(ncol(m$values))
          else unlist(lapply(samples, sample_columns, m = m))
  score <- row_summary(m$values[, cols, drop = FALSE], using)
  ord <- unlist(lapply(seq_len(nrow(m$groups)), function(g) {
    if (m$groups$n[g] == 0) return(integer(0))
    idx <- m$groups$first[g]:m$groups$last[g]
    o <- idx[order(score[idx], na.last = TRUE, method = "radix")]
    if (method == "descend") rev(o) else o
  }))
  m$regions <- m$regions[ord, , drop = FALSE]
  rownames(m$regions) <- NULL
  m$values <- m$values[ord, , drop = FALSE]
  m
}

# Average profile per sample and region group: what plotProfile draws and
# writes with --outFileNameData, for --averageType `average`. Returns a long
# data.frame: sample, group, x (bin, 0-based), value.
group_means <- function(m, average = "mean") {
  out <- list()
  for (i in seq_len(nrow(m$samples))) {
    cols <- sample_columns(m, i)
    for (g in seq_len(nrow(m$groups))) {
      rows <- if (m$groups$n[g] == 0) integer(0) else m$groups$first[g]:m$groups$last[g]
      v <- col_summary(m$values[rows, cols, drop = FALSE], average)
      out[[length(out) + 1]] <- data.frame(
        sample = m$samples$label[i], group = m$groups$label[g],
        x = seq_along(cols) - 1, value = v, stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  res$sample <- factor(res$sample, levels = unique(m$samples$label))
  res$group <- factor(res$group, levels = unique(m$groups$label))
  res
}

# Shrink a block of rows to at most `n` rows by averaging consecutive rows
# (NaN-aware), so a heatmap of tens of thousands of regions is drawn from as
# many rows as the output can show. Row k of the result covers input rows
# breaks[k] .. breaks[k+1]-1; returns the matrix with attr "rows" listing the
# input row count behind each output row.
bin_rows <- function(values, n) {
  nr <- nrow(values)
  if (nr <= n) {
    attr(values, "rows") <- rep(1L, nr)
    return(values)
  }
  bin <- as.integer(floor((seq_len(nr) - 1) * n / nr)) + 1L
  s <- rowsum(ifelse(is.na(values), 0, values), bin, reorder = TRUE)
  k <- rowsum((!is.na(values)) * 1, bin, reorder = TRUE)
  out <- s / k
  out[k == 0] <- NA_real_
  dimnames(out) <- NULL
  attr(out, "rows") <- as.integer(tabulate(bin, n))
  out
}

# x-axis ticks for sample i, as plotProfile/plotHeatmap place them
# (deeptools heatmapper_utilities.getProfileTicks): positions on the 0-based
# bin axis and their labels. reference-point matrices: upstream distance,
# `ref_label`, downstream distance; scale-regions: upstream, `start_label`,
# `end_label`, downstream. deeptools writes the upstream distance without a
# unit ("-2.0") and the downstream one with it ("2.0Kb").
matrix_ticks <- function(m, i = 1, ref_label = "TSS", start_label = "TSS", end_label = "TES") {
  s <- m$samples[i, ]
  w <- s$bin_size; b <- s$upstream; a <- s$downstream; body <- s$body
  c5 <- s$unscaled5; c3 <- s$unscaled3
  if (b < 1e5) { q <- 1000; sym <- "Kb" } else { q <- 1e6; sym <- "Mb" }
  fmt <- function(v) formatC(v, format = "f", digits = 1)
  if (body == 0) {
    pos <- c(0, b - 0.5 * w, b + a - w) / w
    lab <- c(fmt(-(b / q)), ref_label, paste0(fmt(a / q), sym))
  } else {
    pos <- 0
    lab <- character(0)
    if (b > 0) { pos <- c(pos, b); lab <- c(lab, fmt(-(b / q))) }
    lab <- c(lab, start_label)
    if (c5 > 0) { pos <- c(pos, b + c5); lab <- c(lab, "") }
    if (c3 > 0) { pos <- c(pos, b + c5 + body); lab <- c(lab, "") }
    pos <- c(pos, b + c5 + body + c3 - w); lab <- c(lab, end_label)
    if (a > 0) {
      pos <- c(pos, b + c5 + body + c3 + a - w)
      lab <- c(lab, paste0(fmt(a / q), sym))
    }
    pos <- pmax(pos / w, 0)
  }
  data.frame(x = pos, label = lab, stringsAsFactors = FALSE)
}
