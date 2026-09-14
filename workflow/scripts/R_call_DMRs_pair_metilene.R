#!/usr/bin/env Rscript
# Metilene-backed driver for the call_DMRs_for_pair_context rule.
#
# Reads the per-(replicate, context) RDS caches produced by
# cache_mc_replicate_for_context, writes a metilene-format TSV (chr,
# pos, per-replicate methylation rate columns), calls the metilene
# binary, parses its BED-like output, and writes our standard DMR
# table + per-context hyper/hypo counts file.
#
# Per-context thresholds are passed as command-line args from the Snakemake
# rule and are configurable via dmr_thresholds: in epicc-options.yaml.
# Plant-tuned defaults: -d CG=0.3, CHG=0.2, CHH=0.1; -v CG/CHG=0.7, CHH=0.3
# (valley lowered for CHH because per-site CHH signal is intrinsically noisy).

suppressPackageStartupMessages({
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
threads         <- as.integer(args[1])
context         <- args[2]
chromsizes_path <- args[3]  # unused (kept for arg-order parity with DMRcaller path)
dmrs_out        <- args[4]
counts_out      <- args[5]
# Threshold args[6:15] — shared (min_diff, min_cytosines, q_value) + caller-specific
min_diff        <- as.numeric(args[6])   # -d: mean methylation difference threshold
min_cytosines   <- as.integer(args[7])   # -m: minimum CpGs per DMR
# args[8] (bin_size), args[10] (min_gap), args[12] (min_reads) are
# dmrcaller-specific with no metilene equivalent — ignored here. (metilene
# has no between-DMR merge step and consumes per-position rates, not coverage.)
q_value         <- as.numeric(args[9])   # significance cutoff applied to metilene's
                                         # FDR q-value (col 4); see filtering below
min_size        <- as.integer(args[11])  # minimum DMR length (nt); post-hoc filter,
                                         # equivalent to DMRcaller minSize and to
                                         # metilene_output.pl's -l (end-start >= length)
maxdist         <- as.integer(args[13])  # -M/--maxdist: maximum distance (nt) between
                                         # adjacent CpGs joinable into one segment
valley          <- as.numeric(args[14])  # -v: per-position valley filter
maxseg          <- as.integer(args[15])  # -G: max segment length (-1 = no cap)
n1              <- as.integer(args[16])
reps1_paths     <- args[17:(16 + n1)]
reps2_paths     <- args[(17 + n1):length(args)]
n2 <- length(reps2_paths)

# ---------------------------------------------------------------------------
# Per-context metilene parameters (read from command-line args, configured
# via dmr_thresholds: in epicc-options.yaml)
# ---------------------------------------------------------------------------
# maxseg (-G) caps the maximum continuous segment length to bound per-segment
# memory. Default -1 (no cap) works for CG and CHG. CHH on dense-coverage
# long chromosomes (e.g. ColCEN Chr1) produces single segments of ~8 million
# CpGs that drive metilene RSS to 400-700 GB and OOM-kill the job, so the
# default config caps CHH at 10000. See epicc-options.yaml for rationale.
p <- list(d = min_diff, m = min_cytosines, M = maxdist, v = valley, G = maxseg)
if (is.null(p$d)) stop(sprintf("Unknown context '%s'", context))

# ---------------------------------------------------------------------------
# Load per-replicate RDS caches and compute methylation rate per position
# ---------------------------------------------------------------------------
reps1 <- lapply(reps1_paths, readRDS)
reps2 <- lapply(reps2_paths, readRDS)
all_reps <- c(reps1, reps2)

cat(sprintf("R_call_DMRs_pair_metilene.R: n1=%d n2=%d context=%s threads=%d\n",
            n1, n2, context, threads))

# Prominent warning when either group has fewer than 2 replicates.
# Metilene's MWU and 2D-KS tests assume biological replication; with
# n=1 in either group there is no variance estimate and the reported
# q-values are NOT meaningful significance values (effectively a sign
# test on the single observation).
if (n1 < 2 || n2 < 2) {
  bar <- strrep("=", 78)
  msg <- paste(
    "", bar,
    sprintf("WARNING: metilene called with insufficient replicates (n1=%d, n2=%d)", n1, n2),
    bar,
    "  Metilene's MWU and 2D-KS tests assume biological replication.",
    "  With n=1 in either group there is no variance estimate, so the",
    "  reported q-values are NOT meaningful significance values --",
    "  treat the output as candidate DMR coordinates only, not as",
    "  statistically supported calls.",
    "",
    "  For meaningful q-values, rerun with >=2 biological replicates",
    "  in each group.",
    bar, "",
    sep = "\n"
  )
  cat(msg, file = stderr())
}

# Compute the union of (chr, pos) across all replicates.
#
# Bismark-aligned samples have CX reports covering every genomic position
# of the context, so all per-replicate caches have identical positions.
# dmC samples come from `modkit pileup → bedMethyl → CX`, where bedMethyl
# only contains positions with coverage — so different replicates have
# different position sets. Take the union and pad missing positions with
# NA, which metilene reads as "." (missing).
#
# Keep chr and pos as parallel vectors and dedupe via a packed numeric
# key (chr_idx * 1e10 + pos). The earlier paste/strsplit/do.call(rbind)
# pattern aborted R with SIGBUS on ~35M-key CHH inputs because rbind on
# a list of millions of length-2 char vectors blows up internally. 1e10
# is well above any plant chromosome length (<2e8) and well under 2^53,
# so the packed keys remain exact in double precision.
chr_list <- lapply(all_reps, function(gr) as.character(seqnames(gr)))
pos_list <- lapply(all_reps, function(gr) start(gr))
all_chr <- unlist(chr_list, use.names = FALSE)
all_pos <- unlist(pos_list, use.names = FALSE)
chr_levels <- unique(all_chr)
all_keys_num <- match(all_chr, chr_levels) * 1e10 + all_pos
keep <- !duplicated(all_keys_num)
union_chr <- all_chr[keep]
union_pos <- all_pos[keep]
union_keys_num <- all_keys_num[keep]
ord <- order(union_chr, union_pos)
union_chr <- union_chr[ord]
union_pos <- union_pos[ord]
union_keys_num <- union_keys_num[ord]
n_pos <- length(union_chr)

# Per-replicate methylation rates aligned to the union of positions.
rates <- lapply(all_reps, function(gr) {
  rep_keys_num <- match(as.character(seqnames(gr)), chr_levels) * 1e10 + start(gr)
  m <- mcols(gr)
  rep_rate <- m$readsM / m$readsN
  rep_rate[m$readsN == 0] <- NA
  rep_rate[match(union_keys_num, rep_keys_num)]
})

# Metilene segfaults (SIGSEGV, exit 139) when a group has no usable value at
# any position -- it still allocates and walks the per-group value arrays it
# never filled. That happens whenever a replicate's CX report is all-zero
# coverage, e.g. a sample whose FASTQ came back empty. Refuse the run here so
# the failure names the sample instead of surfacing as a core dump.
all_names <- c(paste0("g1_v", seq_len(n1)), paste0("g2_v", seq_len(n2)))
all_paths <- c(reps1_paths, reps2_paths)
usable <- vapply(rates, function(r) sum(!is.na(r)), numeric(1))
if (any(usable == 0)) {
  dead <- which(usable == 0)
  stop(sprintf(
    paste0("no covered %s position in %d of %d replicate(s): %s\n",
           "  Cached from: %s\n",
           "  These samples have zero coverage in this context, so metilene has ",
           "nothing to test. Check the mapping reports for them."),
    context, length(dead), length(rates),
    paste(all_names[dead], collapse = ", "),
    paste(all_paths[dead], collapse = ", ")))
}

# ---------------------------------------------------------------------------
# Write the metilene input TSV
# ---------------------------------------------------------------------------
sample_names <- all_names

tmpdir <- Sys.getenv("TMPDIR", "/tmp")
pid <- Sys.getpid()
input_tsv  <- file.path(tmpdir, sprintf("metilene_input_%s_%d.tsv",  context, pid))
output_tsv <- file.path(tmpdir, sprintf("metilene_output_%s_%d.tsv", context, pid))

# Build a data.frame keyed (chr, pos) with one column per replicate.
input_df <- data.frame(
  chr = union_chr,
  pos = union_pos,
  stringsAsFactors = FALSE
)
for (i in seq_along(rates)) {
  rate_col <- rates[[i]]
  rate_col_chr <- sprintf("%.4f", rate_col)
  rate_col_chr[is.na(rate_col)] <- "."
  input_df[[sample_names[i]]] <- rate_col_chr
}

# Header line (metilene -h: tab-separated; "chr\tpos\t<sample1>\t<sample2>\t...")
con <- file(input_tsv, "w")
writeLines(paste(c("chr", "pos", sample_names), collapse = "\t"), con)
write.table(input_df, con, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
close(con)
cat(sprintf("Wrote metilene input: %s (%d positions, %d samples)\n",
            input_tsv, n_pos, length(sample_names)))

# ---------------------------------------------------------------------------
# Run metilene
# ---------------------------------------------------------------------------
metilene_args <- c(
  "-a", "g1",
  "-b", "g2",
  "-d", as.character(p$d),
  "-m", as.character(p$m),
  "-M", as.character(p$M),
  # -G only when capping segment length. Passing `-G -1` as a "no cap"
  # sentinel breaks metilene's argparser (it reads `-1` as another flag);
  # metilene's internal default for --maxseg is already -1 (no cap), so
  # omitting -G is the correct way to express "no cap".
  if (p$G > 0L) c("-G", as.character(p$G)) else character(0),
  "-v", as.character(p$v),
  "-t", as.character(max(1L, threads)),
  "-X", "1",
  "-Y", "1",
  input_tsv
)
cat("Running: metilene", paste(metilene_args, collapse = " "), "\n")
rc <- system2("metilene", args = metilene_args, stdout = output_tsv)
cat(sprintf("metilene exit code: %d\n", rc))
if (rc != 0) stop(sprintf("metilene failed (exit %d)", rc))

# ---------------------------------------------------------------------------
# Parse metilene output and translate to our standard DMR table
# ---------------------------------------------------------------------------
# metilene de-novo (-f 1) output columns:
#   chr  start  end  q-value  meth-diff  #cpgs  pMWU  p2D-KS  mean-g1  mean-g2
# No header. start is 0-based, end is BED-style exclusive.
out_cols <- c("Chr", "Start", "End", "qValue", "meth_diff",
              "n_cpgs", "pMWU", "p2D_KS", "mean_g1", "mean_g2")
res <- tryCatch(
  read.table(output_tsv, sep = "\t", header = FALSE,
             col.names = out_cols, fill = TRUE,
             stringsAsFactors = FALSE),
  error = function(e) {
    cat(sprintf("Note: metilene output empty or unparseable (%s)\n",
                conditionMessage(e)))
    data.frame()
  }
)

# Post-hoc filters. metilene emits every candidate segment that clears the
# segmentation-stage effect-size/CpG filters (-d/-m), but does NOT apply a
# significance or DMR-length cutoff — that is the job of the second-stage
# metilene_output tool (`-q`/`-l`), which we replicate here. DMRcaller applies
# both internally (pValueThreshold on its BH-adjusted pValue, and minSize) and
# returns only passing regions, so filtering here brings the two callers into
# concordance: both _DMRs.txt tables are filtered at `q_value` and `min_size`.
# The strict `<` on qValue and `>=` on length match metilene_output.pl exactly.
n_called <- nrow(res)
if (n_called > 0) {
  res <- res[!is.na(res$qValue) &
               res$qValue < q_value &
               (res$End - res$Start) >= min_size, , drop = FALSE]
  cat(sprintf("metilene: %d candidate segments, %d pass qValue < %g & length >= %d\n",
              n_called, nrow(res), q_value, min_size))
}

if (nrow(res) > 0) {
  df_out <- data.frame(
    Chr          = res$Chr,
    Start        = res$Start,
    End          = res$End,
    firstsample  = res$mean_g1,
    secondsample = res$mean_g2,
    qValue       = res$qValue,
    Delta        = res$mean_g1 - res$mean_g2
  )
  write.table(df_out, dmrs_out, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  n_hyper <- sum(df_out$Delta > 0)
  n_hypo  <- sum(df_out$Delta <= 0)
  cat(sprintf("Wrote %d DMRs (hyper=%d, hypo=%d) to %s\n",
              nrow(df_out), n_hyper, n_hypo, dmrs_out))
} else {
  write.table(data.frame(Chr = character(), Start = integer(), End = integer(),
                         firstsample = numeric(), secondsample = numeric(),
                         qValue = numeric(), Delta = numeric()),
              dmrs_out, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  n_hyper <- 0L; n_hypo <- 0L
  cat(sprintf("No DMRs called; wrote empty header to %s\n", dmrs_out))
}

counts <- data.frame(Type = c("hyper", "hypo"))
counts[[paste0(context, "_DMRs")]] <- c(n_hyper, n_hypo)
write.table(counts, counts_out, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Cleanup TMPDIR scratch
suppressWarnings(file.remove(input_tsv, output_tsv))
