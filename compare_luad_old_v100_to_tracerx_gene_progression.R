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

parse_tracerx_ccf_max <- function(x) {
  parts <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  vals <- suppressWarnings(as.numeric(sub("^.*:", "", parts)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA_real_)
  max(vals)
}

load_driver_genes <- function(driver_tsv) {
  if (is.null(driver_tsv) || !nzchar(driver_tsv)) return(NULL)
  if (!file.exists(driver_tsv)) stop("Missing/invalid --driver_tsv")

  d <- read.delim(driver_tsv, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  candidate_cols <- c("GENE_SYMBOL", "Gene Symbol", "gene_symbol", "gene")
  gene_col <- intersect(candidate_cols, names(d))[1]
  if (is.na(gene_col)) {
    stop("Cannot find gene symbol column in driver TSV. Expected one of: GENE_SYMBOL, Gene Symbol, gene_symbol, gene")
  }

  genes <- toupper(trimws(as.character(d[[gene_col]])))
  genes <- unique(genes[nzchar(genes) & genes != "NA"])
  genes
}

safe_quantile <- function(x, p) {
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE, na.rm = TRUE))
}

summarize_gene_progression <- function(df, gene_col, sample_col, ccf_col, cohort_label) {
  total_samples <- length(unique(df[[sample_col]]))

  split_idx <- split(seq_len(nrow(df)), df[[gene_col]])
  rows <- lapply(names(split_idx), function(g) {
    idx <- split_idx[[g]]
    ccf_vals <- df[[ccf_col]][idx]
    sample_vals <- unique(df[[sample_col]][idx])

    n_samples_gene <- length(sample_vals)
    prevalence <- if (total_samples > 0) n_samples_gene / total_samples else NA_real_
    median_ccf <- median(ccf_vals, na.rm = TRUE)
    mean_ccf <- mean(ccf_vals, na.rm = TRUE)
    high_ccf_rate <- mean(ccf_vals >= 0.9, na.rm = TRUE)
    # Weighted earliness score: frequent + high-CCF genes rank earlier.
    progression_score <- 0.6 * prevalence + 0.4 * median_ccf

    data.frame(
      cohort = cohort_label,
      gene = g,
      n_mutations = length(idx),
      n_samples_with_gene = n_samples_gene,
      prevalence = prevalence,
      median_ccf = median_ccf,
      mean_ccf = mean_ccf,
      p90_ccf = safe_quantile(ccf_vals, 0.9),
      high_ccf_rate = high_ccf_rate,
      progression_score = progression_score,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(-out$progression_score, -out$prevalence, -out$median_ccf, out$gene), , drop = FALSE]
  out$rank_progression <- seq_len(nrow(out))
  out
}

cli <- parse_args(args)
old_root <- cli$old_root
luad_fast_root <- cli$luad_fast_root
tracerx_rda <- cli$tracerx_rda
arg_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
if (is.na(arg_file) || !nzchar(arg_file)) arg_file <- "./compare_luad_old_v100_to_tracerx_gene_progression.R"
script_dir <- dirname(normalizePath(arg_file, winslash = "/", mustWork = FALSE))
default_candidates <- c(
  file.path(script_dir, "cosmic_lung_nsclc_driver_genes.tsv"),
  file.path(script_dir, "..", "cosmic_lung_nsclc_driver_genes.tsv"),
  file.path(script_dir, "..", "..", "cosmic_lung_nsclc_driver_genes.tsv")
)
default_hit <- default_candidates[file.exists(default_candidates)][1]
driver_tsv <- if (!is.null(cli$driver_tsv) && nzchar(cli$driver_tsv)) cli$driver_tsv else if (!is.na(default_hit)) normalizePath(default_hit, winslash = "/", mustWork = TRUE) else NULL
out_dir <- if (!is.null(cli$out_dir) && nzchar(cli$out_dir)) cli$out_dir else file.path(old_root, "tracerx_gene_progression")

if (is.null(old_root) || !dir.exists(old_root)) stop("Missing/invalid --old_root")
if (is.null(luad_fast_root) || !dir.exists(luad_fast_root)) stop("Missing/invalid --luad_fast_root")
if (is.null(tracerx_rda) || !file.exists(tracerx_rda)) stop("Missing/invalid --tracerx_rda")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

driver_genes <- load_driver_genes(driver_tsv)
is_driver_mode <- !is.null(driver_genes)
file_suffix <- if (is_driver_mode) "_driver_only" else ""

old_files <- list.files(old_root, pattern = "_old_v100_ccf\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(old_files) == 0) stop("No old-v100 CCF files found")

fast_files <- list.files(luad_fast_root, pattern = "_ccf\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(fast_files) == 0) stop("No LUAD fast CCF files found")

read_old_one <- function(p) {
  d <- read.delim(p, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  req <- c("sample", "chr", "pos")
  if (!all(req %in% names(d))) return(NULL)

  if ("ccf_old_v100_raw" %in% names(d)) {
    ccf_vals <- as.numeric(d$ccf_old_v100_raw)
  } else if ("ccf_old_v100" %in% names(d)) {
    ccf_vals <- as.numeric(d$ccf_old_v100)
  } else if ("ccf_old_v100_valid" %in% names(d)) {
    ccf_vals <- as.numeric(d$ccf_old_v100_valid)
  } else {
    return(NULL)
  }

  keep <- is.finite(ccf_vals)
  if (!any(keep)) return(NULL)

  out <- data.frame(
    sample = as.character(d$sample[keep]),
    chr = gsub("^chr", "", as.character(d$chr[keep]), ignore.case = TRUE),
    pos = as.integer(d$pos[keep]),
    ccf = ccf_vals[keep],
    stringsAsFactors = FALSE
  )
  out
}

read_fast_one <- function(p) {
  d <- read.delim(p, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  req <- c("chr", "pos", "gene")
  if (!all(req %in% names(d))) return(NULL)

  sample_id <- substr(basename(p), 1, 15)
  gene_vals <- toupper(trimws(as.character(d$gene)))
  keep <- nzchar(gene_vals) & gene_vals != "NA"
  if (!any(keep)) return(NULL)

  out <- data.frame(
    sample = sample_id,
    chr = gsub("^chr", "", as.character(d$chr[keep]), ignore.case = TRUE),
    pos = as.integer(d$pos[keep]),
    gene = gene_vals[keep],
    stringsAsFactors = FALSE
  )
  out
}

old_list <- lapply(old_files, read_old_one)
old_list <- old_list[!vapply(old_list, is.null, logical(1))]
if (length(old_list) == 0) stop("No valid old CCF rows found")
old_df <- do.call(rbind, old_list)

fast_list <- lapply(fast_files, read_fast_one)
fast_list <- fast_list[!vapply(fast_list, is.null, logical(1))]
if (length(fast_list) == 0) stop("No valid fast gene rows found")
fast_df <- do.call(rbind, fast_list)

# De-duplicate potential one-to-many mappings by keeping the first mapping per key.
fast_df <- fast_df[!duplicated(fast_df[, c("sample", "chr", "pos")]), , drop = FALSE]

luad_mut <- merge(old_df, fast_df, by = c("sample", "chr", "pos"), all.x = FALSE, all.y = FALSE)
if (nrow(luad_mut) == 0) stop("No overlap between old CCF rows and fast gene annotations")

if (is_driver_mode) {
  luad_mut <- luad_mut[luad_mut$gene %in% driver_genes, , drop = FALSE]
  if (nrow(luad_mut) == 0) stop("No LUAD mutations remain after driver-gene filtering")
}

luad_gene <- summarize_gene_progression(
  df = luad_mut,
  gene_col = "gene",
  sample_col = "sample",
  ccf_col = "ccf",
  cohort_label = "TCGA_LUAD_old_dpclust3p_v100_raw"
)

envx <- new.env()
load(tracerx_rda, envir = envx)
if (!exists("TRACERx_NEJM_2017", envir = envx)) stop("TRACERx_NEJM_2017 object not found in RDA")
tx <- envx$TRACERx_NEJM_2017

if (!all(c("patientID", "variantID", "CCF") %in% names(tx))) {
  stop("TRACERx_NEJM_2017 missing required columns: patientID, variantID, CCF")
}

tx_df <- data.frame(
  sample = as.character(tx$patientID),
  gene = toupper(trimws(as.character(tx$variantID))),
  ccf = vapply(tx$CCF, parse_tracerx_ccf_max, numeric(1)),
  stringsAsFactors = FALSE
)

tx_df <- tx_df[is.finite(tx_df$ccf) & nzchar(tx_df$gene) & tx_df$gene != "NA", , drop = FALSE]
if (nrow(tx_df) == 0) stop("No valid TRACERx rows after CCF parsing")

if (is_driver_mode) {
  tx_df <- tx_df[tx_df$gene %in% driver_genes, , drop = FALSE]
  if (nrow(tx_df) == 0) stop("No TRACERx mutations remain after driver-gene filtering")
}

tracerx_gene <- summarize_gene_progression(
  df = tx_df,
  gene_col = "gene",
  sample_col = "sample",
  ccf_col = "ccf",
  cohort_label = "TRACERx_NEJM_2017"
)

common_genes <- intersect(luad_gene$gene, tracerx_gene$gene)
luad_only <- setdiff(luad_gene$gene, tracerx_gene$gene)
tracerx_only <- setdiff(tracerx_gene$gene, luad_gene$gene)

shared <- merge(
  tracerx_gene[, c("gene", "rank_progression", "progression_score", "prevalence", "median_ccf", "n_samples_with_gene")],
  luad_gene[, c("gene", "rank_progression", "progression_score", "prevalence", "median_ccf", "n_samples_with_gene")],
  by = "gene",
  suffixes = c("_tracerx", "_luad")
)
shared <- shared[order(shared$rank_progression_tracerx), , drop = FALSE]

# Re-rank on common genes only so progression order comparison uses a shared universe.
tracerx_common <- tracerx_gene[tracerx_gene$gene %in% common_genes, , drop = FALSE]
luad_common <- luad_gene[luad_gene$gene %in% common_genes, , drop = FALSE]
tracerx_common <- tracerx_common[order(-tracerx_common$progression_score, -tracerx_common$prevalence, -tracerx_common$median_ccf, tracerx_common$gene), , drop = FALSE]
luad_common <- luad_common[order(-luad_common$progression_score, -luad_common$prevalence, -luad_common$median_ccf, luad_common$gene), , drop = FALSE]
if (nrow(tracerx_common) > 0) tracerx_common$rank_progression_common <- seq_len(nrow(tracerx_common))
if (nrow(luad_common) > 0) luad_common$rank_progression_common <- seq_len(nrow(luad_common))

shared_common_only <- merge(
  tracerx_common[, c("gene", "rank_progression_common", "progression_score", "prevalence", "median_ccf", "n_samples_with_gene")],
  luad_common[, c("gene", "rank_progression_common", "progression_score", "prevalence", "median_ccf", "n_samples_with_gene")],
  by = "gene",
  suffixes = c("_tracerx", "_luad")
)
shared_common_only <- shared_common_only[order(shared_common_only$rank_progression_common_tracerx), , drop = FALSE]

rank_metrics <- data.frame(
  shared_gene_count = length(common_genes),
  tracerx_gene_count = nrow(tracerx_gene),
  luad_gene_count = nrow(luad_gene),
  common_gene_fraction_tracerx = ifelse(nrow(tracerx_gene) > 0, length(common_genes) / nrow(tracerx_gene), NA_real_),
  common_gene_fraction_luad = ifelse(nrow(luad_gene) > 0, length(common_genes) / nrow(luad_gene), NA_real_),
  spearman_rho = ifelse(nrow(shared_common_only) > 1, suppressWarnings(cor(shared_common_only$rank_progression_common_tracerx, shared_common_only$rank_progression_common_luad, method = "spearman")), NA_real_),
  kendall_tau = ifelse(nrow(shared_common_only) > 1, suppressWarnings(cor(shared_common_only$rank_progression_common_tracerx, shared_common_only$rank_progression_common_luad, method = "kendall")), NA_real_),
  top20_overlap_count = length(intersect(head(tracerx_common$gene, 20), head(luad_common$gene, 20))),
  top20_jaccard = {
    a <- head(tracerx_common$gene, 20)
    b <- head(luad_common$gene, 20)
    u <- union(a, b)
    if (length(u) == 0) NA_real_ else length(intersect(a, b)) / length(u)
  },
  driver_filter_applied = is_driver_mode,
  driver_gene_reference_count = ifelse(is_driver_mode, length(driver_genes), NA_integer_),
  stringsAsFactors = FALSE
)

top20_shared <- data.frame(
  gene = intersect(head(tracerx_common$gene, 20), head(luad_common$gene, 20)),
  stringsAsFactors = FALSE
)

write.table(tracerx_gene, file = file.path(out_dir, paste0("tracerx_gene_progression", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(luad_gene, file = file.path(out_dir, paste0("luad_gene_progression", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(shared, file = file.path(out_dir, paste0("shared_genes_progression_comparison", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(shared_common_only, file = file.path(out_dir, paste0("shared_genes_progression_comparison_common_only", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = tracerx_only, stringsAsFactors = FALSE), file = file.path(out_dir, paste0("tracerx_only_genes", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = luad_only, stringsAsFactors = FALSE), file = file.path(out_dir, paste0("luad_only_genes", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(rank_metrics, file = file.path(out_dir, paste0("progression_rank_agreement_summary", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(top20_shared, file = file.path(out_dir, paste0("top20_overlap_genes", file_suffix, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

cat("Done.\n")
cat(sprintf("LUAD mapped mutations: %d\n", nrow(luad_mut)))
cat(sprintf("TRACERx parsed mutations: %d\n", nrow(tx_df)))
cat(sprintf("Shared genes: %d\n", length(common_genes)))
if (is_driver_mode) {
  cat(sprintf("Driver reference genes: %d\n", length(driver_genes)))
  cat("Driver-gene filter: applied\n")
}
cat(sprintf("Output dir: %s\n", out_dir))
