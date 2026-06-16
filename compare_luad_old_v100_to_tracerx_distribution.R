#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(x) {
  out <- list()
  for (a in x) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else ""
    out[[key]] <- val
  }
  out
}

cli <- parse_args(args)
old_root <- cli$old_root
tracerx_rda <- cli$tracerx_rda
out_dir <- if (!is.null(cli$out_dir) && nzchar(cli$out_dir)) cli$out_dir else old_root

if (is.null(old_root) || !dir.exists(old_root)) stop("Missing/invalid --old_root")
if (is.null(tracerx_rda) || !file.exists(tracerx_rda)) stop("Missing/invalid --tracerx_rda")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ccf_files <- list.files(old_root, pattern = "_old_v100_ccf\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(ccf_files) == 0) stop("No old_v100 CCF files found under old_root")

read_one <- function(p) {
  d <- read.delim(p, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if ("ccf_old_v100_raw" %in% names(d)) {
    d$ccf_for_compare <- as.numeric(d$ccf_old_v100_raw)
  } else if ("ccf_old_v100" %in% names(d)) {
    d$ccf_for_compare <- as.numeric(d$ccf_old_v100)
  } else if ("ccf_old_v100_valid" %in% names(d)) {
    d$ccf_for_compare <- as.numeric(d$ccf_old_v100_valid)
  } else {
    return(NULL)
  }

  d <- d[is.finite(d$ccf_for_compare), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d[, c("sample", "chr", "pos", "ccf_for_compare"), drop = FALSE]
}

old_list <- lapply(ccf_files, read_one)
old_list <- old_list[!vapply(old_list, is.null, logical(1))]
old_df <- do.call(rbind, old_list)

envx <- new.env()
load(tracerx_rda, envir = envx)
if (!exists("TRACERx_NEJM_2017", envir = envx)) stop("TRACERx_NEJM_2017 object not found in RDA")
tx <- envx$TRACERx_NEJM_2017
if (!"CCF" %in% names(tx)) stop("TRACERx object missing CCF column")

extract_ccf_vals <- function(s) {
  parts <- unlist(strsplit(as.character(s), ";", fixed = TRUE))
  vals <- suppressWarnings(as.numeric(sub("^.*:", "", parts)))
  vals[is.finite(vals)]
}

tx_vals <- unlist(lapply(tx$CCF, extract_ccf_vals))
old_vals <- as.numeric(old_df$ccf_for_compare)
old_vals <- old_vals[is.finite(old_vals)]

summ <- function(x, label) {
  data.frame(
    cohort = label,
    n = length(x),
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    p10 = as.numeric(quantile(x, 0.10, names = FALSE)),
    p90 = as.numeric(quantile(x, 0.90, names = FALSE)),
    stringsAsFactors = FALSE
  )
}

summary_df <- rbind(summ(old_vals, "TCGA_LUAD_old_dpclust3p_v100"), summ(tx_vals, "TRACERx_NEJM_2017"))

ks <- suppressWarnings(ks.test(old_vals, tx_vals))
ks_df <- data.frame(
  statistic_D = as.numeric(ks$statistic),
  p_value = as.numeric(ks$p.value),
  stringsAsFactors = FALSE
)

overlap_df <- data.frame(
  tcga_samples = length(unique(old_df$sample)),
  tracerx_patients = length(unique(tx$patientID)),
  direct_id_overlap = length(intersect(unique(old_df$sample), unique(tx$patientID))),
  stringsAsFactors = FALSE
)

write.table(summary_df, file = file.path(out_dir, "tracerx_vs_tcga_luad_ccf_distribution_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ks_df, file = file.path(out_dir, "tracerx_vs_tcga_luad_ccf_distribution_ks.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(overlap_df, file = file.path(out_dir, "tracerx_vs_tcga_luad_id_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("Done.\n")
cat(sprintf("old CCF values: %d\n", length(old_vals)))
cat(sprintf("TRACERx CCF values: %d\n", length(tx_vals)))
cat(sprintf("ID overlap: %d\n", overlap_df$direct_id_overlap[1]))
