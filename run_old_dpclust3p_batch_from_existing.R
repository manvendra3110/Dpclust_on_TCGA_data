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

input_root <- cli$input_root
output_root <- cli$output_root
sample_list_file <- cli$sample_list
clean_lib <- if (!is.null(cli$lib) && nzchar(cli$lib)) cli$lib else "C:/CloneHD_benchmarking/DPclust/.r_libs_dpclust3p_v100"
ref_tag <- if (!is.null(cli$dpclust3p_ref) && nzchar(cli$dpclust3p_ref)) cli$dpclust3p_ref else "dpclust3p-v1.0.0"

if (is.null(input_root) || !nzchar(input_root)) stop("Missing --input_root")
if (is.null(output_root) || !nzchar(output_root)) stop("Missing --output_root")
if (!dir.exists(input_root)) stop(sprintf("Input root does not exist: %s", input_root))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(clean_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(clean_lib, .libPaths()))

options(repos = c(CRAN = "https://cran.rstudio.com"))

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing)
}

install_if_missing(c("remotes"))

required_bioc <- c("VariantAnnotation", "GenomicRanges", "Rsamtools", "IRanges", "S4Vectors", "GenomeInfoDb")
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  stop(sprintf("Missing required Bioconductor packages in current R installation: %s", paste(missing_bioc, collapse = ", ")))
}

if (!requireNamespace("dpclust3p", quietly = TRUE)) {
  remotes::install_github("Wedge-lab/dpclust3p", ref = ref_tag, upgrade = "never", dependencies = TRUE)
}

suppressPackageStartupMessages(library(dpclust3p))

patch_dpclust3p_reader_bug <- function() {
  fixed <- function(loci_file, allele_frequencies_file) {
    subs.data <- tryCatch(
      read.table(loci_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE),
      error = function(e) NA
    )

    if (length(subs.data) == 1 && is.na(subs.data)) {
      return(NULL)
    }

    subs.data <- subs.data[order(subs.data[, 1], subs.data[, 2]), ]
    subs.data[, 3] <- apply(as.data.frame(subs.data[, 3]), 1, function(x) substring(x, 1, 1))
    subs.data[, 4] <- apply(as.data.frame(subs.data[, 4]), 1, function(x) substring(x, 1, 1))

    subs.data.gr <- GenomicRanges::GRanges(
      subs.data[, 1],
      IRanges::IRanges(subs.data[, 2], subs.data[, 2]),
      rep("*", nrow(subs.data))
    )
    S4Vectors::elementMetadata(subs.data.gr) <- subs.data[, c(3, 4)]

    alleleFrequencies <- read.delim(
      allele_frequencies_file,
      sep = "\t",
      header = TRUE,
      quote = NULL,
      stringsAsFactors = FALSE
    )
    alleleFrequencies <- alleleFrequencies[order(alleleFrequencies[, 1], alleleFrequencies[, 2]), ]

    alleleFrequencies.gr <- GenomicRanges::GRanges(
      alleleFrequencies[, 1],
      IRanges::IRanges(alleleFrequencies[, 2], alleleFrequencies[, 2]),
      rep("*", nrow(alleleFrequencies))
    )
    S4Vectors::elementMetadata(alleleFrequencies.gr) <- alleleFrequencies[, 3:7]

    overlap <- GenomicRanges::findOverlaps(subs.data.gr, alleleFrequencies.gr)
    alleleFrequencies <- alleleFrequencies[S4Vectors::subjectHits(overlap), ]

    nucleotides <- c("A", "C", "G", "T")
    ref.indices <- match(subs.data[, 3], nucleotides)
    alt.indices <- match(subs.data[, 4], nucleotides)

    WT.count <- as.numeric(unlist(sapply(1:nrow(alleleFrequencies), function(v, a, i) v[i, a[i] + 2], v = alleleFrequencies, a = ref.indices)))
    mut.count <- as.numeric(unlist(sapply(1:nrow(alleleFrequencies), function(v, a, i) v[i, a[i] + 2], v = alleleFrequencies, a = alt.indices)))

    combined.gr <- GenomicRanges::GRanges(GenomicRanges::seqnames(subs.data.gr), GenomicRanges::ranges(subs.data.gr), rep("*", nrow(subs.data)))
    S4Vectors::elementMetadata(combined.gr) <- data.frame(WT.count = WT.count, mut.count = mut.count)
    combined.gr
  }

  assignInNamespace("GetWTandMutCount", fixed, ns = "dpclust3p")

  fixed_cellularity <- function(rho_and_psi_file) {
    d <- read.table(rho_and_psi_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    if ("cellularity" %in% names(d)) {
      return(as.numeric(d$cellularity[1]))
    }
    if ("rho" %in% names(d) && "FRAC_GENOME" %in% rownames(d)) {
      return(as.numeric(d["FRAC_GENOME", "rho"]))
    }

    numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
    if (length(numeric_cols) > 0) {
      return(as.numeric(d[[numeric_cols[1]]][1]))
    }

    stop("Unable to parse cellularity from file.")
  }
  assignInNamespace("GetCellularity", fixed_cellularity, ns = "dpclust3p")

  if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    gdp <- get("GetDirichletProcessInfo", envir = asNamespace("dpclust3p"))
    gdp_body <- paste(deparse(body(gdp)), collapse = "\n")
    gdp_body <- gsub("sortSeqlevels\\(", "GenomeInfoDb::sortSeqlevels(", gdp_body)
    gdp_body <- gsub("subjectHits\\(", "S4Vectors::subjectHits(", gdp_body)
    gdp_body <- gsub("queryHits\\(", "S4Vectors::queryHits(", gdp_body)
    body(gdp) <- parse(text = gdp_body)[[1]]
    assignInNamespace("GetDirichletProcessInfo", gdp, ns = "dpclust3p")
  }
  invisible(TRUE)
}

patch_dpclust3p_reader_bug()

if (!is.null(sample_list_file) && nzchar(sample_list_file)) {
  samples <- unique(substr(readLines(sample_list_file, warn = FALSE), 1, 15))
  samples <- samples[nzchar(samples)]
} else {
  kids <- list.files(input_root, full.names = FALSE)
  sample_dirs <- kids[file.info(file.path(input_root, kids))$isdir]
  samples <- sample_dirs[grepl("^TCGA-", sample_dirs)]
}

samples <- sort(unique(samples))
cat(sprintf("Processing %d samples from %s\n", length(samples), input_root))

rows <- list()
fails <- list()

for (i in seq_along(samples)) {
  s <- samples[[i]]
  cat(sprintf("[%d/%d] %s\n", i, length(samples), s))

  in_dir <- file.path(input_root, s)
  if (!dir.exists(in_dir)) {
    fails[[length(fails) + 1]] <- data.frame(sample = s, stage = "input_dir", error = "missing input sample directory", stringsAsFactors = FALSE)
    next
  }

  out_dir <- file.path(output_root, s)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fast_file <- file.path(in_dir, paste0(s, "_ccf.tsv"))
  loci_file <- file.path(out_dir, paste0(s, "_loci_for_dpclust3p.txt"))
  allele_file <- file.path(in_dir, paste0(s, "_alleleFrequencies.txt"))
  cellularity_file <- file.path(in_dir, paste0(s, ".cellularity_ploidy"))
  subclone_file <- file.path(in_dir, paste0(s, ".battenberg"))
  out_file <- file.path(out_dir, paste0(s, "_dpclust_input_old_v100.txt"))
  out_ccf <- file.path(out_dir, paste0(s, "_old_v100_ccf.tsv"))

  req <- c(fast_file, allele_file, cellularity_file, subclone_file)
  miss <- req[!file.exists(req)]
  if (length(miss) > 0) {
    fails[[length(fails) + 1]] <- data.frame(sample = s, stage = "input_files", error = paste("missing:", paste(miss, collapse = "; ")), stringsAsFactors = FALSE)
    next
  }

  ok <- TRUE
  err <- ""

  tryCatch({
    fast_df <- read.delim(fast_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("chr", "pos", "ref", "alt") %in% names(fast_df))) {
      stop("fast ccf file missing chr/pos/ref/alt columns")
    }

    loci_df <- fast_df[, c("chr", "pos", "ref", "alt")]
    write.table(loci_df, file = loci_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    dpclust3p::runGetDirichletProcessInfo(
      loci_file = loci_file,
      allele_frequencies_file = allele_file,
      cellularity_file = cellularity_file,
      subclone_file = subclone_file,
      gender = "female",
      SNP.phase.file = "NA",
      mut.phase.file = "NA",
      output_file = out_file
    )

    d <- read.delim(out_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    if (all(c("chr", "end", "subclonal.fraction") %in% names(d))) {
      ccf_raw <- as.numeric(d$subclonal.fraction)
      ccf_valid <- ifelse(is.finite(ccf_raw) & ccf_raw >= 0 & ccf_raw <= 1, ccf_raw, NA_real_)
      ccf_df <- data.frame(
        sample = s,
        chr = d$chr,
        pos = as.integer(d$end),
        ccf_old_v100_raw = ccf_raw,
        ccf_old_v100_valid = ccf_valid,
        ccf_over_one_flag = ifelse(is.finite(ccf_raw) & ccf_raw > 1, TRUE, FALSE),
        stringsAsFactors = FALSE
      )
      write.table(ccf_df, file = out_ccf, sep = "\t", row.names = FALSE, quote = FALSE)
    }

    rows[[length(rows) + 1]] <- data.frame(
      sample = s,
      n_rows = nrow(d),
      n_non_na_ccf = sum(is.finite(as.numeric(d$subclonal.fraction))),
      n_ccf_gt1 = sum(is.finite(as.numeric(d$subclonal.fraction)) & as.numeric(d$subclonal.fraction) > 1),
      mean_ccf = mean(as.numeric(d$subclonal.fraction), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    ok <<- FALSE
    err <<- conditionMessage(e)
  })

  if (!ok) {
    fails[[length(fails) + 1]] <- data.frame(sample = s, stage = "run", error = err, stringsAsFactors = FALSE)
  }
}

sum_df <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(sample = character(), n_rows = integer(), n_non_na_ccf = integer(), n_ccf_gt1 = integer(), mean_ccf = numeric(), stringsAsFactors = FALSE)
fail_df <- if (length(fails) > 0) do.call(rbind, fails) else data.frame(sample = character(), stage = character(), error = character(), stringsAsFactors = FALSE)

write.table(sum_df, file = file.path(output_root, "old_v100_processing_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fail_df, file = file.path(output_root, "old_v100_failures.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Done. Success: %d, Failed: %d\n", nrow(sum_df), nrow(fail_df)))
cat(sprintf("Output root: %s\n", output_root))
