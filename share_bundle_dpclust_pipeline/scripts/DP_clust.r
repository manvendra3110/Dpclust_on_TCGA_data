options(repos = c(CRAN = "https://cran.rstudio.com"))

install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing)
  }
}

install_dpclust_dependencies <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    c("VariantAnnotation", "GenomicRanges", "Rsamtools", "IRanges", "S4Vectors"),
    ask = FALSE,
    update = FALSE
  )

  install_if_missing(c("optparse", "ggplot2", "reshape2", "remotes"))

  remotes::install_github("Wedge-lab/dpclust3p", upgrade = "never", dependencies = TRUE)
  remotes::install_github("Wedge-lab/dpclust", upgrade = "never", dependencies = TRUE)
}

pick_first_column <- function(df_names, candidates, label) {
  hit <- intersect(df_names, candidates)
  if (length(hit) == 0) {
    stop(sprintf("Could not find required %s column. Tried: %s", label, paste(candidates, collapse = ", ")))
  }
  hit[1]
}

prepare_dpclust3p_input <- function(maf_file, donor_id, output_prefix) {
  if (!file.exists(maf_file)) {
    stop(sprintf("MAF file not found: %s", maf_file))
  }

  out_dir <- dirname(output_prefix)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  maf <- read.delim(maf_file, stringsAsFactors = FALSE, comment.char = "#")

  sample_col <- pick_first_column(
    names(maf),
    c("Sample", "sample", "Tumor_Sample_Barcode"),
    "sample"
  )

  donor_maf <- maf[maf[[sample_col]] == donor_id, , drop = FALSE]
  cat(sprintf("Found %d mutations for %s\n", nrow(donor_maf), donor_id))
  if (nrow(donor_maf) == 0) {
    stop("No rows found for donor_id in MAF.")
  }

  chr_col <- pick_first_column(names(donor_maf), c("chr", "Chromosome", "chromosome"), "chromosome")
  pos_col <- pick_first_column(names(donor_maf), c("start", "Start_Position", "Start_position"), "position")
  ref_col <- pick_first_column(names(donor_maf), c("reference", "Reference_Allele", "ref"), "reference allele")
  alt_col <- pick_first_column(names(donor_maf), c("alt", "Tumor_Seq_Allele2", "ALT"), "alternate allele")
  t_alt_col <- pick_first_column(names(donor_maf), c("t_alt_count", "T_alt_count"), "alternate read count")
  t_ref_col <- pick_first_column(names(donor_maf), c("t_ref_count", "T_ref_count"), "reference read count")

  donor_maf$chr_clean <- gsub("^chr", "", donor_maf[[chr_col]], ignore.case = TRUE)

  is_snv <- nchar(donor_maf[[ref_col]]) == 1 &
    nchar(donor_maf[[alt_col]]) == 1 &
    donor_maf[[ref_col]] != "-" &
    donor_maf[[alt_col]] != "-"

  donor_maf <- donor_maf[is_snv, , drop = FALSE]
  cat(sprintf("After SNV filter: %d mutations\n", nrow(donor_maf)))
  if (nrow(donor_maf) == 0) {
    stop("No SNVs remained after filtering.")
  }

  loci <- data.frame(
    chr = donor_maf$chr_clean,
    pos = donor_maf[[pos_col]],
    ref = donor_maf[[ref_col]],
    alt = donor_maf[[alt_col]],
    stringsAsFactors = FALSE
  )

  loci_file <- paste0(output_prefix, "_loci.txt")
  write.table(
    loci,
    loci_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  allele_freq <- data.frame(
    `#CHR` = donor_maf$chr_clean,
    POS = donor_maf[[pos_col]],
    Count_A = integer(nrow(donor_maf)),
    Count_C = integer(nrow(donor_maf)),
    Count_G = integer(nrow(donor_maf)),
    Count_T = integer(nrow(donor_maf)),
    Good_depth = as.integer(donor_maf[[t_alt_col]]) + as.integer(donor_maf[[t_ref_col]]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(donor_maf))) {
    ref <- toupper(donor_maf[[ref_col]][i])
    alt <- toupper(donor_maf[[alt_col]][i])
    ref_count <- as.integer(donor_maf[[t_ref_col]][i])
    alt_count <- as.integer(donor_maf[[t_alt_col]][i])

    ref_col_name <- paste0("Count_", ref)
    alt_col_name <- paste0("Count_", alt)

    if (ref_col_name %in% names(allele_freq)) {
      allele_freq[[ref_col_name]][i] <- ref_count
    }
    if (alt_col_name %in% names(allele_freq)) {
      allele_freq[[alt_col_name]][i] <- alt_count
    }
  }

  allele_file <- paste0(output_prefix, "_alleleFreq.txt")
  write.table(
    allele_freq,
    allele_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  cat(sprintf("Written loci file: %s\n", loci_file))
  cat(sprintf("Written allele frequencies file: %s\n", allele_file))

  list(loci = loci_file, allele_freq = allele_file)
}

prepare_dpclust3p_input_from_mutect <- function(mutect_vcf_file, output_prefix) {
  if (!file.exists(mutect_vcf_file)) {
    stop(sprintf("Mutect file not found: %s", mutect_vcf_file))
  }

  out_dir <- dirname(output_prefix)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  vcf_lines <- readLines(mutect_vcf_file, warn = FALSE)
  header_idx <- grep("^#CHROM\\t", vcf_lines)
  if (length(header_idx) == 0) {
    stop("Could not find #CHROM header line in Mutect VCF file.")
  }

  vcf <- read.delim(
    text = paste(vcf_lines[header_idx:length(vcf_lines)], collapse = "\n"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  req_cols <- c("#CHROM", "POS", "REF", "ALT", "FORMAT", "tumor")
  if (!all(req_cols %in% names(vcf))) {
    stop(sprintf("Mutect file missing required columns. Required: %s", paste(req_cols, collapse = ", ")))
  }

  is_snv <- nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1 & vcf$REF != "-" & vcf$ALT != "-"
  snv <- vcf[is_snv, , drop = FALSE]
  cat(sprintf("Found %d SNVs in Mutect file\n", nrow(snv)))
  if (nrow(snv) == 0) {
    stop("No SNVs found in Mutect input.")
  }

  parse_ad <- function(format_str, tumor_str) {
    fmt <- strsplit(format_str, ":", fixed = TRUE)[[1]]
    tvals <- strsplit(tumor_str, ":", fixed = TRUE)[[1]]
    ad_idx <- match("AD", fmt)
    if (is.na(ad_idx) || ad_idx > length(tvals)) {
      return(c(NA_integer_, NA_integer_))
    }
    ad_vals <- strsplit(tvals[ad_idx], ",", fixed = TRUE)[[1]]
    if (length(ad_vals) < 2) {
      return(c(NA_integer_, NA_integer_))
    }
    c(as.integer(ad_vals[1]), as.integer(ad_vals[2]))
  }

  ad_counts <- t(mapply(parse_ad, snv$FORMAT, snv$tumor))
  snv$ref_count <- ad_counts[, 1]
  snv$alt_count <- ad_counts[, 2]

  bad_counts <- is.na(snv$ref_count) | is.na(snv$alt_count)
  if (any(bad_counts)) {
    cat(sprintf("Dropping %d SNVs with missing AD counts\n", sum(bad_counts)))
    snv <- snv[!bad_counts, , drop = FALSE]
  }
  if (nrow(snv) == 0) {
    stop("No SNVs with usable AD counts remained after filtering.")
  }

  chr_clean <- gsub("^chr", "", snv$`#CHROM`, ignore.case = TRUE)

  loci <- data.frame(
    chr = chr_clean,
    pos = snv$POS,
    ref = snv$REF,
    alt = snv$ALT,
    stringsAsFactors = FALSE
  )

  loci_file <- paste0(output_prefix, "_loci.txt")
  write.table(
    loci,
    loci_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  allele_freq <- data.frame(
    `#CHR` = chr_clean,
    POS = snv$POS,
    Count_A = integer(nrow(snv)),
    Count_C = integer(nrow(snv)),
    Count_G = integer(nrow(snv)),
    Count_T = integer(nrow(snv)),
    Good_depth = as.integer(snv$ref_count) + as.integer(snv$alt_count),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(snv))) {
    ref_col <- paste0("Count_", toupper(snv$REF[i]))
    alt_col <- paste0("Count_", toupper(snv$ALT[i]))

    if (ref_col %in% names(allele_freq)) {
      allele_freq[[ref_col]][i] <- as.integer(snv$ref_count[i])
    }
    if (alt_col %in% names(allele_freq)) {
      allele_freq[[alt_col]][i] <- as.integer(snv$alt_count[i])
    }
  }

  allele_file <- paste0(output_prefix, "_alleleFreq.txt")
  write.table(
    allele_freq,
    allele_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  cat(sprintf("Written loci file: %s\n", loci_file))
  cat(sprintf("Written allele frequencies file: %s\n", allele_file))

  list(loci = loci_file, allele_freq = allele_file)
}

validate_dpclust3p_inputs <- function(loci_file, allele_file) {
  req_af <- c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")

  loci <- read.delim(loci_file, stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
  af <- read.delim(allele_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(loci) < 4) {
    stop("Loci file must contain at least 4 columns: chr, pos, ref, alt (without header).")
  }
  if (!all(req_af %in% names(af))) {
    stop(sprintf("Allele frequency file missing columns. Required: %s", paste(req_af, collapse = ", ")))
  }
  if (nrow(loci) != nrow(af)) {
    stop("Loci and allele files have different row counts.")
  }

  same_order <- all(loci[[1]] == af$`#CHR`) && all(as.character(loci[[2]]) == as.character(af$POS))
  if (!same_order) {
    stop("Loci and allele files are not aligned row-by-row by chromosome and position.")
  }

  invisible(TRUE)
}

patch_dpclust3p_reader_bug <- function() {
  if (!requireNamespace("dpclust3p", quietly = TRUE)) {
    return(invisible(FALSE))
  }

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

    stop("Unable to parse cellularity from file. Expected a 'cellularity' or 'rho' value.")
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

run_dpclust3p_on_pcawg <- function(
  variant_file,
  cn_file,
  purity_file,
  donor_id,
  gender = "female",
  output_dir = "results",
  input_format = c("maf", "mutect_vcf")
) {
  input_format <- match.arg(input_format)

  if (!file.exists(cn_file)) {
    stop(sprintf("Subclone/Battenberg file not found: %s", cn_file))
  }
  if (!file.exists(purity_file)) {
    stop(sprintf("Cellularity file not found: %s", purity_file))
  }

  if (!requireNamespace("dpclust3p", quietly = TRUE)) {
    stop("Package dpclust3p is not installed. Run install_dpclust_dependencies() first.")
  }

  patch_dpclust3p_reader_bug()

  output_prefix <- file.path(output_dir, donor_id)

  cat("Step 1: Preparing loci and allele input files...\n")
  if (input_format == "maf") {
    files <- prepare_dpclust3p_input(variant_file, donor_id, output_prefix)
  } else {
    files <- prepare_dpclust3p_input_from_mutect(variant_file, output_prefix)
  }
  validate_dpclust3p_inputs(files$loci, files$allele_freq)

  cat("Step 2: Running dpclust3p preprocessing...\n")
  dpclust_input <- paste0(output_prefix, "_dpclust_input.txt")

  dpclust3p::runGetDirichletProcessInfo(
    loci_file = files$loci,
    allele_frequencies_file = files$allele_freq,
    cellularity_file = purity_file,
    subclone_file = cn_file,
    gender = gender,
    SNP.phase.file = "NA",
    mut.phase.file = "NA",
    output_file = dpclust_input
  )

  cat(sprintf("Done. dpclust3p output: %s\n", dpclust_input))
  dpclust_input
}

# CCF-only mode: DPClust clustering is intentionally commented out.
# run_dpclust_clustering <- function(dpclust_input_file, sample_id = "DO51110", cellularity = 0.885, output_dir = "results") {
#   if (!requireNamespace("DPClust", quietly = TRUE)) {
#     stop("Package DPClust is not installed. Run install_dpclust_dependencies() first.")
#   }
#   if (!file.exists(dpclust_input_file)) {
#     stop(sprintf("DPClust input file not found: %s", dpclust_input_file))
#   }
#
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }
#
#   # DPClust uses Sys.time() directly in filenames; sanitize for Windows paths.
#   rundp_fn <- get("RunDP", envir = asNamespace("DPClust"))
#   rundp_body <- paste(deparse(body(rundp_fn)), collapse = "\n")
#   rundp_body <- gsub(
#     "chartr\\(\" \", \"_\", Sys.time\\(\\)\\)",
#     "format(Sys.time(), \"%Y-%m-%d_%H-%M-%S\")",
#     rundp_body
#   )
#   body(rundp_fn) <- parse(text = rundp_body)[[1]]
#   assignInNamespace("RunDP", rundp_fn, ns = "DPClust")
#
#   run_params <- DPClust::make_run_params(
#     no.iters = 2000,
#     no.iters.burn.in = 1000,
#     mut.assignment.type = 1,
#     num_muts_sample = NA,
#     is.male = TRUE,
#     min_muts_cluster = 10
#   )
#   run_params$datpath <- normalizePath(dirname(dpclust_input_file), winslash = "/", mustWork = TRUE)
#
#   sample_params <- DPClust::make_sample_params(
#     datafiles = basename(dpclust_input_file),
#     cellularity = cellularity,
#     is.male = TRUE,
#     samplename = sample_id,
#     subsamples = sample_id
#   )
#
#   advanced_params <- DPClust::make_advanced_params(seed = 123)
#
#   cna_params <- list(
#     co_cluster_cna = FALSE,
#     cndatafiles = NA,
#     add.conflicts = FALSE,
#     cna.conflicting.events.only = FALSE,
#     num.clonal.events.to.add = 0,
#     min.cna.size = 0
#   )
#
#   DPClust::RunDP(
#     analysis_type = "nd_dp",
#     run_params = run_params,
#     sample_params = sample_params,
#     advanced_params = advanced_params,
#     outdir = file.path(output_dir, "dpclust_output"),
#     cna_params = cna_params
#   )
#
#   cat(sprintf("Done. DPClust output folder: %s\n", file.path(output_dir, "dpclust_output")))
# }

cat("Loaded DP_clust.r pipeline functions.\n")
cat("1) install_dpclust_dependencies()\n")
cat("2) run_dpclust3p_on_pcawg(variant_file, cn_file, purity_file, donor_id, gender, output_dir, input_format)\n")
cat("3) DPClust clustering section is commented out (CCF-only mode).\n")

print_cli_usage <- function() {
  cat("\nUsage:\n")
  cat("  Rscript DP_clust.r --donor_id=DO1000 --data_dir=C:/CloneHD_benchmarking/PCWGA/data --gender=male --input_format=mutect_vcf\n")
  cat("\n")
  cat("  OR provide explicit files:\n")
  cat("  Rscript DP_clust.r --donor_id=DO1000 --variant_file=... --cn_file=... --purity_file=... --output_dir=... --input_format=mutect_vcf\n")
  cat("\nOptions:\n")
  cat("  --donor_id=VALUE             Required donor/sample id (example: DO1000)\n")
  cat("  --data_dir=PATH              Base data folder containing donor subfolders\n")
  cat("  --variant_file=PATH          Variant input file (MAF or Mutect VCF-like)\n")
  cat("  --cn_file=PATH               Battenberg file\n")
  cat("  --purity_file=PATH           Cellularity/ploidy file\n")
  cat("  --output_dir=PATH            Output folder (default: <sample_dir>/results)\n")
  cat("  --input_format=maf|mutect_vcf  Input variant format (default: mutect_vcf)\n")
  cat("  --gender=male|female         Gender for dpclust3p (default: female)\n")
  cat("  --install_deps=true|false    Install packages before running (default: false)\n")
}

parse_cli_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) {
      next
    }
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else ""
    out[[key]] <- value
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || identical(x, "")) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

run_from_cli <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

  if (length(args) == 0 || as_bool(args$help, FALSE) || as_bool(args$h, FALSE)) {
    print_cli_usage()
    quit(save = "no", status = 0)
  }

  donor_id <- args$donor_id
  if (is.null(donor_id) || identical(donor_id, "")) {
    stop("Missing required argument: --donor_id")
  }

  input_format <- if (!is.null(args$input_format) && nzchar(args$input_format)) args$input_format else "mutect_vcf"
  gender <- if (!is.null(args$gender) && nzchar(args$gender)) args$gender else "female"

  variant_file <- args$variant_file
  cn_file <- args$cn_file
  purity_file <- args$purity_file
  output_dir <- args$output_dir

  if (!is.null(args$data_dir) && nzchar(args$data_dir)) {
    sample_dir <- file.path(args$data_dir, donor_id)

    if (is.null(variant_file) || !nzchar(variant_file)) {
      variant_file <- if (input_format == "maf") {
        file.path(sample_dir, paste0(donor_id, ".maf"))
      } else {
        file.path(sample_dir, paste0(donor_id, ".mutect"))
      }
    }
    if (is.null(cn_file) || !nzchar(cn_file)) {
      cn_file <- file.path(sample_dir, paste0(donor_id, ".battenberg"))
    }
    if (is.null(purity_file) || !nzchar(purity_file)) {
      purity_file <- file.path(sample_dir, paste0(donor_id, ".cellularity_ploidy"))
    }
    if (is.null(output_dir) || !nzchar(output_dir)) {
      output_dir <- file.path(sample_dir, "results")
    }
  }

  if (is.null(variant_file) || !nzchar(variant_file) ||
      is.null(cn_file) || !nzchar(cn_file) ||
      is.null(purity_file) || !nzchar(purity_file)) {
    stop("Provide either --data_dir with --donor_id, or explicit --variant_file, --cn_file, and --purity_file.")
  }

  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- "results"
  }

  if (as_bool(args$install_deps, FALSE)) {
    cat("Installing dependencies...\n")
    install_dpclust_dependencies()
  }

  cat("Running CCF preprocessing...\n")
  out <- run_dpclust3p_on_pcawg(
    variant_file = variant_file,
    cn_file = cn_file,
    purity_file = purity_file,
    donor_id = donor_id,
    gender = gender,
    output_dir = output_dir,
    input_format = input_format
  )

  cat(sprintf("\nCCF output file: %s\n", out))
}

if (sys.nframe() == 0) {
  run_from_cli()
}

