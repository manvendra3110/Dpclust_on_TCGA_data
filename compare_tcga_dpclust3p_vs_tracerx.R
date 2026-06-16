#!/usr/bin/env Rscript
# compare_tcga_dpclust3p_vs_tracerx.R
# End-to-end cohort-level comparison: TCGA LUAD dpclust3p CCF vs TRACERx LUAD CCF.

suppressPackageStartupMessages({
  library(ggplot2)
})

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else ""
    out[[key]] <- val
  }
  out
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE, na.rm = TRUE))
}

safe_sample <- function(x, max_n = 10000) {
  if (length(x) <= max_n) return(x)
  sample(x, max_n)
}

clean_ccf <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x < 0] <- NA_real_
  x
}

clean_chr <- function(x) {
  out <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
  numeric_like <- grepl("^[0-9]+\\.0*$", out)
  out[numeric_like] <- as.character(as.integer(as.numeric(out[numeric_like])))
  out[out == "23"] <- "X"
  out[out == "24"] <- "Y"
  out
}

parse_tracerx_ccf_max <- function(ccf_str) {
  if (is.na(ccf_str) || !nzchar(as.character(ccf_str))) return(NA_real_)
  parts <- strsplit(as.character(ccf_str), ";", fixed = TRUE)[[1]]
  vals <- suppressWarnings(as.numeric(sub("^R[0-9]+:", "", parts)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(NA_real_)
  max(vals)
}

parse_tracerx_n_regions <- function(ccf_str) {
  if (is.na(ccf_str) || !nzchar(as.character(ccf_str))) return(1L)
  length(strsplit(as.character(ccf_str), ";", fixed = TRUE)[[1]])
}

parse_tracerx <- function(rda_path) {
  cat("Loading TRACERx RDA...\n")
  env <- new.env(parent = emptyenv())
  load(rda_path, envir = env)

  obj_names <- ls(env)
  if (length(obj_names) == 0) {
    stop("No objects found in TRACERx RDA")
  }

  tx <- NULL
  if ("TRACERx_NEJM_2017" %in% obj_names) {
    tx <- env[["TRACERx_NEJM_2017"]]
    obj_name <- "TRACERx_NEJM_2017"
  } else {
    obj_name <- obj_names[1]
    tx <- env[[obj_name]]
  }

  if (!is.data.frame(tx)) {
    stop(sprintf("TRACERx object '%s' is not a data.frame", obj_name))
  }
  if (!"variantID" %in% names(tx)) stop("TRACERx data missing 'variantID'")
  if (!"patientID" %in% names(tx)) stop("TRACERx data missing 'patientID'")

  cat(sprintf("  Object: %s\n", obj_name))
  cat(sprintf("  Rows: %d\n", nrow(tx)))

  tx$chr <- sub("^CRUK[0-9]+:([^:]+):.*$", "\\1", tx$variantID)
  tx$chr <- clean_chr(tx$chr)
  tx$pos <- suppressWarnings(as.integer(sub("^CRUK[0-9]+:[^:]+:([0-9]+):.*$", "\\1", tx$variantID)))
  tx$alt <- sub("^CRUK[0-9]+:[^:]+:[0-9]+:([A-Za-z-]+)$", "\\1", tx$variantID)

  ccf_candidates <- c("CCF", "ccf", "PyCloneCCF", "pyclone_ccf")
  ccf_col <- intersect(ccf_candidates, names(tx))[1]
  if (is.na(ccf_col)) {
    for (col in names(tx)) {
      sample_val <- tx[[col]][which(!is.na(tx[[col]]) & nzchar(as.character(tx[[col]])))[1]]
      if (length(sample_val) == 1 && grepl("^R[0-9]+:[0-9.]+", as.character(sample_val))) {
        ccf_col <- col
        break
      }
    }
  }
  if (is.na(ccf_col)) stop("Cannot find TRACERx CCF column")

  cat(sprintf("  CCF column: %s\n", ccf_col))
  tx$ccf_max <- vapply(tx[[ccf_col]], parse_tracerx_ccf_max, numeric(1))
  tx$n_regions <- vapply(tx[[ccf_col]], parse_tracerx_n_regions, integer(1))
  tx$is_clonal_ccf <- tx$ccf_max >= 0.8

  if ("is.clonal" %in% names(tx)) {
    tx$is_clonal_pyclone <- as.logical(tx[["is.clonal"]])
  }

  cat(sprintf("  Parsed mutations: %d\n", nrow(tx)))
  cat(sprintf("  Patients: %d\n", length(unique(tx$patientID))))
  tx
}

load_mc3_annotation <- function(maf_path, luad_samples) {
  cat("Loading mc3 MAF for gene annotation...\n")
  if (!file.exists(maf_path)) stop("Missing/invalid maf_path")

  header <- read.delim(
    maf_path,
    nrows = 0,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "Tumor_Sample_Barcode")
  optional_cols <- c("Variant_Type")
  missing_cols <- setdiff(required_cols, names(header))
  if (length(missing_cols) > 0) {
    stop(sprintf("MAF missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  keep_cols <- c(required_cols, intersect(optional_cols, names(header)))
  col_classes <- rep("NULL", length(names(header)))
  names(col_classes) <- names(header)
  col_classes[keep_cols] <- NA

  maf <- read.delim(
    maf_path,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = col_classes
  )

  if ("Variant_Type" %in% names(maf)) {
    maf <- maf[maf$Variant_Type == "SNP", , drop = FALSE]
  }

  luad_prefixes <- unique(substr(luad_samples, 1, 15))
  maf$sample_15 <- substr(maf$Tumor_Sample_Barcode, 1, 15)
  maf <- maf[maf$sample_15 %in% luad_prefixes, , drop = FALSE]

  maf$chr_clean <- clean_chr(maf$Chromosome)
  maf$pos_clean <- suppressWarnings(as.integer(maf$Start_Position))
  maf$join_key <- paste(maf$chr_clean, maf$pos_clean, sep = ":")

  annot <- maf[!duplicated(maf$join_key), c("join_key", "Hugo_Symbol", "chr_clean", "pos_clean"), drop = FALSE]
  cat(sprintf("  Annotation rows: %d unique chr:pos positions\n", nrow(annot)))
  annot
}

read_tcga_old_ccf <- function(fp) {
  d <- read.delim(fp, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("sample", "chr", "pos")
  if (!all(req %in% names(d))) return(NULL)

  if ("ccf_old_v100_raw" %in% names(d)) {
    ccf <- clean_ccf(d$ccf_old_v100_raw)
  } else if ("ccf_old_v100" %in% names(d)) {
    ccf <- clean_ccf(d$ccf_old_v100)
  } else if ("ccf_old_v100_valid" %in% names(d)) {
    ccf <- as.numeric(d$ccf_old_v100_valid)
  } else {
    return(NULL)
  }

  keep <- is.finite(ccf)
  if (!any(keep)) return(NULL)

  data.frame(
    sample = as.character(d$sample[keep]),
    chr = clean_chr(d$chr[keep]),
    pos = suppressWarnings(as.integer(d$pos[keep])),
    ccf_raw = if ("ccf_old_v100_raw" %in% names(d)) as.numeric(d$ccf_old_v100_raw[keep]) else ccf[keep],
    ccf = ccf[keep],
    stringsAsFactors = FALSE
  )
}

read_tcga_dpclust_input <- function(fp) {
  d <- read.delim(fp, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("chr", "end", "subclonal.fraction")
  if (!all(req %in% names(d))) return(NULL)

  sample_id <- basename(fp)
  sample_id <- sub("_(dpclust_input.*|dpclust_input_old_v100)\\.txt$", "", sample_id)
  ccf_raw <- suppressWarnings(as.numeric(d$subclonal.fraction))
  ccf <- clean_ccf(ccf_raw)
  keep <- is.finite(ccf)
  if (!any(keep)) return(NULL)

  data.frame(
    sample = sample_id,
    chr = clean_chr(d$chr[keep]),
    pos = suppressWarnings(as.integer(d$end[keep])),
    ccf_raw = ccf_raw[keep],
    ccf = ccf[keep],
    stringsAsFactors = FALSE
  )
}

load_tcga_ccf <- function(ccf_dir, maf_annotation) {
  cat("Loading TCGA dpclust3p CCF files...\n")
  if (!dir.exists(ccf_dir)) stop("Missing/invalid ccf_dir")

  old_ccf_files <- list.files(
    ccf_dir,
    pattern = "_old_v100_ccf\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  dpclust_input_files <- list.files(
    ccf_dir,
    pattern = "_dpclust_input.*\\.txt$",
    recursive = TRUE,
    full.names = TRUE
  )

  tcga_rows <- list()
  if (length(old_ccf_files) > 0) {
    cat(sprintf("  Found %d old-v100 CCF TSV files\n", length(old_ccf_files)))
    tcga_rows <- lapply(old_ccf_files, read_tcga_old_ccf)
  } else if (length(dpclust_input_files) > 0) {
    cat(sprintf("  Found %d dpclust input files\n", length(dpclust_input_files)))
    tcga_rows <- lapply(dpclust_input_files, read_tcga_dpclust_input)
  } else {
    stop("No supported TCGA CCF files found in ccf_dir")
  }

  tcga_rows <- Filter(Negate(is.null), tcga_rows)
  if (length(tcga_rows) == 0) stop("No readable TCGA CCF rows found")

  tcga <- do.call(rbind, tcga_rows)
  tcga$join_key <- paste(tcga$chr, tcga$pos, sep = ":")
  tcga <- merge(tcga, maf_annotation, by = "join_key", all.x = TRUE)

  annotated_n <- sum(!is.na(tcga$Hugo_Symbol))
  cat(sprintf("  Loaded %d rows from %d samples\n", nrow(tcga), length(unique(tcga$sample))))
  cat(sprintf("  Gene-annotated: %d / %d (%.1f%%)\n", annotated_n, nrow(tcga), 100 * annotated_n / nrow(tcga)))
  tcga
}

compare_distributions <- function(tcga_ccf, tracerx_ccf, out_dir) {
  cat("Running distribution comparison...\n")
  tcga_vals <- tcga_ccf$ccf[is.finite(tcga_ccf$ccf)]
  tx_vals <- tracerx_ccf$ccf_max[is.finite(tracerx_ccf$ccf_max)]

  summarize_ccf <- function(x, label) {
    data.frame(
      cohort = label,
      n = length(x),
      mean = round(mean(x), 4),
      median = round(median(x), 4),
      sd = round(sd(x), 4),
      p10 = round(safe_quantile(x, 0.10), 4),
      p25 = round(safe_quantile(x, 0.25), 4),
      p75 = round(safe_quantile(x, 0.75), 4),
      p90 = round(safe_quantile(x, 0.90), 4),
      pct_clonal = round(100 * mean(x >= 0.8), 2),
      pct_subclonal = round(100 * mean(x < 0.5), 2),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- rbind(
    summarize_ccf(tcga_vals, "TCGA_LUAD_dpclust3p"),
    summarize_ccf(tx_vals, "TRACERx_LUAD_2017")
  )

  ks <- suppressWarnings(ks.test(tcga_vals, tx_vals))
  ks_df <- data.frame(
    D = round(as.numeric(ks$statistic), 4),
    p_value = as.numeric(ks$p.value),
    stringsAsFactors = FALSE
  )

  write.table(
    summary_df,
    file.path(out_dir, "ccf_distribution_summary.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    ks_df,
    file.path(out_dir, "ccf_distribution_ks_test.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  plot_df <- rbind(
    data.frame(cohort = "TCGA LUAD (dpclust3p)", ccf = safe_sample(tcga_vals)),
    data.frame(cohort = "TRACERx LUAD 2017", ccf = safe_sample(tx_vals))
  )
  p <- ggplot(plot_df, aes(x = ccf, fill = cohort)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("#378ADD", "#D85A30")) +
    labs(
      title = "CCF distribution: TCGA LUAD vs TRACERx LUAD",
      subtitle = sprintf("KS D=%.3f, p=%.2e", ks_df$D, ks_df$p_value),
      x = "Cancer Cell Fraction (CCF)",
      y = "Density",
      fill = "Cohort"
    ) +
    theme_minimal()
  ggsave(file.path(out_dir, "ccf_distribution_density.pdf"), p, width = 8, height = 5)

  cat(sprintf("  KS test: D=%.3f, p=%.2e\n", ks_df$D, ks_df$p_value))
  list(summary = summary_df, ks = ks_df)
}

summarize_gene_ccf <- function(df, ccf_col, sample_col, cohort_label, min_samples) {
  genes <- sort(unique(df$Hugo_Symbol[!is.na(df$Hugo_Symbol) & nzchar(df$Hugo_Symbol)]))
  rows <- lapply(genes, function(gene_name) {
    sub_df <- df[df$Hugo_Symbol == gene_name, , drop = FALSE]
    ccf_vals <- sub_df[[ccf_col]]
    ccf_vals <- ccf_vals[is.finite(ccf_vals)]
    n_samples <- length(unique(sub_df[[sample_col]]))
    if (n_samples < min_samples || length(ccf_vals) == 0) return(NULL)

    data.frame(
      cohort = cohort_label,
      gene = gene_name,
      n_mutations = length(ccf_vals),
      n_samples = n_samples,
      mean_ccf = round(mean(ccf_vals), 4),
      median_ccf = round(median(ccf_vals), 4),
      pct_clonal = round(100 * mean(ccf_vals >= 0.8), 2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

compare_driver_genes <- function(tcga_ccf, tracerx_ccf, out_dir, min_samples = 5) {
  cat("Running per-gene CCF comparison...\n")

  tx_drivers <- tracerx_ccf
  if ("is.driver" %in% names(tx_drivers)) {
    driver_flag <- tx_drivers[["is.driver"]]
    tx_drivers <- tx_drivers[driver_flag %in% c(TRUE, "TRUE", 1, "1"), , drop = FALSE]
    cat(sprintf("  TRACERx driver mutations: %d\n", nrow(tx_drivers)))
  }

  tcga_gene <- summarize_gene_ccf(
    tcga_ccf[!is.na(tcga_ccf$Hugo_Symbol), , drop = FALSE],
    "ccf",
    "sample",
    "TCGA_LUAD",
    min_samples
  )
  tx_gene <- summarize_gene_ccf(
    tx_drivers[!is.na(tx_drivers$Hugo_Symbol), , drop = FALSE],
    "ccf_max",
    "patientID",
    "TRACERx_LUAD",
    min_samples
  )

  if (is.null(tcga_gene) || is.null(tx_gene) || nrow(tcga_gene) == 0 || nrow(tx_gene) == 0) {
    cat("  No shared gene-level comparison available after filtering\n")
    return(list(tcga = tcga_gene, tracerx = tx_gene, shared = NULL))
  }

  shared_genes <- intersect(tcga_gene$gene, tx_gene$gene)
  cat(sprintf("  Genes in TCGA: %d, TRACERx: %d, shared: %d\n", nrow(tcga_gene), nrow(tx_gene), length(shared_genes)))

  shared <- merge(
    tcga_gene[, c("gene", "n_samples", "median_ccf", "pct_clonal")],
    tx_gene[, c("gene", "n_samples", "median_ccf", "pct_clonal")],
    by = "gene",
    suffixes = c("_tcga", "_tracerx")
  )
  shared <- shared[order(-shared$median_ccf_tracerx, shared$gene), , drop = FALSE]

  write.table(
    shared,
    file.path(out_dir, "per_gene_ccf_comparison.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  if (nrow(shared) > 1) {
    rho <- suppressWarnings(cor(shared$median_ccf_tcga, shared$median_ccf_tracerx, method = "spearman"))
    cat(sprintf("  Spearman rho (median CCF per gene): %.3f\n", rho))

    p <- ggplot(shared, aes(x = median_ccf_tcga, y = median_ccf_tracerx, label = gene)) +
      geom_point(aes(size = n_samples_tracerx), colour = "#378ADD", alpha = 0.7) +
      geom_text(size = 3, vjust = -0.5, check_overlap = TRUE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      labs(
        title = "Median CCF per driver gene: TCGA vs TRACERx",
        subtitle = sprintf("Spearman rho = %.3f | n = %d shared genes", rho, nrow(shared)),
        x = "Median CCF (TCGA LUAD)",
        y = "Median CCF (TRACERx LUAD)",
        size = "TRACERx samples"
      ) +
      theme_minimal()
    ggsave(file.path(out_dir, "per_gene_ccf_scatter.pdf"), p, width = 7, height = 7)
  }

  list(tcga = tcga_gene, tracerx = tx_gene, shared = shared)
}

compare_clonality <- function(tcga_ccf, tracerx_ccf, out_dir, clonal_threshold = 0.8, min_samples = 5) {
  cat("Running clonality comparison...\n")

  driver_genes <- c(
    "KRAS", "EGFR", "TP53", "STK11", "KEAP1", "SMARCA4", "RB1",
    "PIK3CA", "CDKN2A", "MET", "BRAF", "NF1"
  )

  get_clonality <- function(df, ccf_col, sample_col, cohort_name) {
    rows <- lapply(driver_genes, function(gene_name) {
      sub_df <- df[!is.na(df$Hugo_Symbol) & df$Hugo_Symbol == gene_name, , drop = FALSE]
      n_samples <- length(unique(sub_df[[sample_col]]))
      if (n_samples < min_samples) return(NULL)

      ccf_vals <- sub_df[[ccf_col]]
      ccf_vals <- ccf_vals[is.finite(ccf_vals)]
      if (length(ccf_vals) == 0) return(NULL)

      data.frame(
        cohort = cohort_name,
        gene = gene_name,
        n_samples = n_samples,
        pct_clonal = round(100 * mean(ccf_vals >= clonal_threshold), 1),
        median_ccf = round(median(ccf_vals), 3),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }

  tcga_cl <- get_clonality(tcga_ccf, "ccf", "sample", "TCGA LUAD")
  tx_cl <- get_clonality(tracerx_ccf, "ccf_max", "patientID", "TRACERx LUAD")

  combined <- rbind(tcga_cl, tx_cl)
  write.table(
    combined,
    file.path(out_dir, "clonality_per_driver_gene.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  shared_genes <- intersect(tcga_cl$gene, tx_cl$gene)
  if (length(shared_genes) > 0) {
    plot_df <- combined[combined$gene %in% shared_genes, , drop = FALSE]
    plot_df$gene <- factor(plot_df$gene, levels = tcga_cl$gene[order(-tcga_cl$pct_clonal)])

    p <- ggplot(plot_df, aes(x = gene, y = pct_clonal, fill = cohort)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("TCGA LUAD" = "#378ADD", "TRACERx LUAD" = "#D85A30")) +
      labs(
        title = "% clonal mutations per driver gene",
        subtitle = paste("Clonal threshold: CCF >=", clonal_threshold),
        x = "Driver gene",
        y = "% clonal mutations",
        fill = "Cohort"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "clonality_barplot.pdf"), p, width = 9, height = 5)
  }

  list(tcga = tcga_cl, tracerx = tx_cl)
}

main <- function(
  tracerx_rda = "TRACERx_NEJM_2017.rda",
  maf_path = "mc3.v0.2.8.PUBLIC.maf",
  ccf_dir = "TCGA/outputs_old/batch_output",
  sample_list = "TCGA/tcga_luad_samples.txt",
  out_dir = "TCGA/outputs_old/validation_final"
) {
  if (!file.exists(tracerx_rda)) stop("Missing/invalid tracerx_rda")
  if (!file.exists(sample_list)) stop("Missing/invalid sample_list")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tracerx <- parse_tracerx(tracerx_rda)

  luad_samples <- readLines(sample_list, warn = FALSE)
  luad_samples <- luad_samples[nzchar(luad_samples)]
  maf_annot <- load_mc3_annotation(maf_path, luad_samples)

  tracerx$join_key <- paste(tracerx$chr, tracerx$pos, sep = ":")
  tracerx <- merge(tracerx, maf_annot[, c("join_key", "Hugo_Symbol")], by = "join_key", all.x = TRUE)
  cat(sprintf("TRACERx mutations with gene annotation: %d / %d\n", sum(!is.na(tracerx$Hugo_Symbol)), nrow(tracerx)))

  tcga <- load_tcga_ccf(ccf_dir, maf_annot)

  distribution_res <- compare_distributions(tcga, tracerx, out_dir)
  gene_res <- compare_driver_genes(tcga, tracerx, out_dir)
  clonality_res <- compare_clonality(tcga, tracerx, out_dir)

  cat(sprintf("\nAll outputs written to: %s\n", out_dir))
  invisible(list(distribution = distribution_res, per_gene = gene_res, clonality = clonality_res))
}

if (sys.nframe() == 0) {
  cli <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  do.call(main, cli[intersect(names(cli), names(formals(main)))])
}

cat("Loaded compare_tcga_dpclust3p_vs_tracerx.R\n")
cat("Call: main(tracerx_rda, maf_path, ccf_dir, sample_list, out_dir)\n")
