# convert_absolute_to_battenberg_v2.R
# Complete fix for ALL identified bugs in ABSOLUTE → dpclust3p CCF pipeline
# Covers: D1-D5 (data conversion), C1-C3 (code patches), N1-N2 (new bugs in v1)

CLONAL_CCF_THRESHOLD <- 0.95

# ─────────────────────────────────────────────────────────────────
# FIX N1: Robust chromosome string conversion
# Handles: float "1.0", integer "1", "X", "Y", "23"(=X), "24"(=Y)
# ─────────────────────────────────────────────────────────────────
clean_chromosome <- function(chr_vec) {
  chr_str <- as.character(chr_vec)
  chr_str <- gsub("^chr", "", chr_str, ignore.case = TRUE)
  # Convert floats like "1.0" → "1"
  is_numeric <- grepl("^[0-9]+\\.0*$", chr_str)
  chr_str[is_numeric] <- as.character(as.integer(as.numeric(chr_str[is_numeric])))
  # Some ABSOLUTE files encode X=23, Y=24
  chr_str[chr_str == "23"] <- "X"
  chr_str[chr_str == "24"] <- "Y"
  chr_str
}

# ─────────────────────────────────────────────────────────────────
# FIX N2: Correct diploid baseline based on genome doubling status
# For WGD samples (Genome_doublings >= 1), baseline is 2+2 not 1+1
# ─────────────────────────────────────────────────────────────────
get_baseline_cn <- function(genome_doublings) {
  # Returns list(major=, minor=) for the pre-event diploid baseline
  if (is.na(genome_doublings) || genome_doublings < 1) {
    list(major = 1L, minor = 1L)
  } else {
    list(major = 2L, minor = 2L)
  }
}

convert_segtab_to_battenberg <- function(segtab_file, tables_file, output_dir) {
  if (!file.exists(segtab_file)) {
    stop(sprintf("Segtab file not found: %s", segtab_file))
  }
  if (!file.exists(tables_file)) {
    stop(sprintf("Tables file not found: %s", tables_file))
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("Reading segtab:", segtab_file, "\n")
  seg <- read.delim(segtab_file, stringsAsFactors = FALSE, check.names = FALSE)

  cat("Reading tables:", tables_file, "\n")
  tab <- read.delim(tables_file, stringsAsFactors = FALSE, check.names = FALSE)

  # FIX D1 + N1: chromosome and coordinate types
  seg$Chromosome <- clean_chromosome(seg$Chromosome)
  seg$Start <- as.integer(seg$Start)
  seg$End   <- as.integer(seg$End)

  # Build per-sample genome_doublings lookup
  gd_col <- grep("Genome.doublings", names(tab), ignore.case = TRUE)[1]
  if (is.na(gd_col)) {
    cat("Warning: 'Genome doublings' column not found in tables; defaulting to 0\n")
    gd_lookup <- setNames(rep(0L, nrow(tab)), tab$array)
  } else {
    gd_lookup <- setNames(as.integer(tab[[gd_col]]), tab$array)
  }

  samples <- unique(seg$Sample)
  cat(sprintf("Converting %d samples...\n", length(samples)))

  success_count <- 0
  for (s in samples) {
    d <- seg[seg$Sample == s, , drop = FALSE]

    # FIX N2: get correct baseline for this sample
    gd <- if (s %in% names(gd_lookup)) gd_lookup[[s]] else 0L
    baseline <- get_baseline_cn(gd)
    base_maj <- baseline$major
    base_min <- baseline$minor

    # Determine which segments are clonal vs subclonal
    clonal_mask    <- d$Cancer_cell_frac_a2 >= CLONAL_CCF_THRESHOLD
    subclonal_mask <- !clonal_mask
    loh_clonal_mask <- (d$LOH == 1) & clonal_mask
    homdel_mask     <- d$Homozygous_deletion == 1

    n <- nrow(d)
    nMaj1 <- integer(n)
    nMin1 <- integer(n)
    frac1  <- numeric(n)
    nMaj2  <- rep(NA_integer_, n)
    nMin2  <- rep(NA_integer_, n)
    frac2  <- rep(NA_real_, n)

    # --- Clonal normal/gain/loss (not LOH, not homdel) ---
    cl <- clonal_mask & !loh_clonal_mask & !homdel_mask
    nMaj1[cl] <- as.integer(d$Modal_HSCN_2[cl])    # FIX D3: HSCN_2 = major
    nMin1[cl] <- as.integer(d$Modal_HSCN_1[cl])
    frac1[cl] <- 1.0
    # nMaj2/nMin2/frac2 stay NA

    # --- Clonal LOH ---
    nMaj1[loh_clonal_mask] <- as.integer(d$Modal_HSCN_2[loh_clonal_mask])
    nMin1[loh_clonal_mask] <- as.integer(d$Modal_HSCN_1[loh_clonal_mask])
    frac1[loh_clonal_mask] <- 1.0

    # --- Subclonal CN events ---
    # FIX D2 + N2: clonal state = baseline (1+1 or 2+2), subclonal state = modal HSCN
    sc <- subclonal_mask & !homdel_mask
    nMaj1[sc] <- base_maj
    nMin1[sc] <- base_min
    frac1[sc] <- 1.0 - d$Cancer_cell_frac_a2[sc]
    nMaj2[sc] <- as.integer(d$Modal_HSCN_2[sc])
    nMin2[sc] <- as.integer(d$Modal_HSCN_1[sc])
    frac2[sc] <- d$Cancer_cell_frac_a2[sc]

    # --- Homozygous deletion ---
    nMaj1[homdel_mask] <- 0L
    nMin1[homdel_mask] <- 0L
    frac1[homdel_mask] <- 1.0

    bb <- data.frame(
      chr      = d$Chromosome,
      startpos = d$Start,
      endpos   = d$End,
      BAF      = NA_real_,   # not in ABSOLUTE; placeholder
      pval     = NA_real_,
      LogR     = NA_real_,
      ntot     = as.integer(d$Modal_Total_CN),
      nMaj1    = nMaj1,
      nMin1    = nMin1,
      frac1    = frac1,
      nMaj2    = nMaj2,
      nMin2    = nMin2,
      frac2    = frac2,
      stringsAsFactors = FALSE
    )

    chr_order <- c(as.character(1:22), "X", "Y")
    bb <- bb[order(match(bb$chr, chr_order), bb$startpos), ]

    out_file <- file.path(output_dir, paste0(s, ".battenberg"))
    write.table(bb, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  %s → %d segments (WGD=%d, baseline %d+%d)\n",
                s, nrow(bb), gd, base_maj, base_min))
    success_count <- success_count + 1
  }

  cat(sprintf("Successfully converted %d samples.\n", success_count))
}

# FIX D4: Write cellularity files with explicit 'cellularity' column name
convert_tables_to_cellularity <- function(tables_file, output_dir) {
  if (!file.exists(tables_file)) {
    stop(sprintf("Tables file not found: %s", tables_file))
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  tab <- read.delim(tables_file, stringsAsFactors = FALSE, check.names = FALSE)

  if (!"array" %in% names(tab)) {
    stop("Tables file missing 'array' column")
  }
  if (!"purity" %in% names(tab)) {
    stop("Tables file missing 'purity' column")
  }

  success_count <- 0
  for (i in seq_len(nrow(tab))) {
    s   <- tab$array[i]
    pur <- as.numeric(tab$purity[i])
    plo <- if ("ploidy" %in% names(tab)) as.numeric(tab$ploidy[i]) else 2.0

    if (is.na(pur)) {
      cat(sprintf("  Skipping %s: purity is NA\n", s))
      next
    }

    cell_df <- data.frame(
      cellularity = pur,
      ploidy      = plo,
      stringsAsFactors = FALSE
    )

    out_file <- file.path(output_dir, paste0(s, ".cellularity_ploidy"))
    write.table(cell_df,
                file = out_file,
                sep = "\t", row.names = FALSE, quote = FALSE)
    success_count <- success_count + 1
  }

  cat(sprintf("Written %d cellularity files to %s\n", success_count, output_dir))
}

# ─────────────────────────────────────────────────────────────────
# FIX C1: Corrected GetWTandMutCount patch
# subs.data subsetted AFTER overlap to keep indices aligned
# ─────────────────────────────────────────────────────────────────
make_fixed_GetWTandMutCount <- function() {
  function(loci_file, allele_frequencies_file) {
    subs.data <- tryCatch(
      read.table(loci_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE),
      error = function(e) NA
    )
    if (length(subs.data) == 1 && is.na(subs.data)) return(NULL)

    subs.data <- subs.data[order(subs.data[, 1], subs.data[, 2]), ]
    subs.data[, 3] <- substring(subs.data[, 3], 1, 1)
    subs.data[, 4] <- substring(subs.data[, 4], 1, 1)

    subs.data.gr <- GenomicRanges::GRanges(
      subs.data[, 1],
      IRanges::IRanges(subs.data[, 2], subs.data[, 2]),
      rep("*", nrow(subs.data))
    )
    S4Vectors::elementMetadata(subs.data.gr) <- subs.data[, c(3, 4)]

    alleleFrequencies <- read.delim(
      allele_frequencies_file, sep = "\t",
      header = TRUE, quote = NULL, stringsAsFactors = FALSE
    )
    alleleFrequencies <- alleleFrequencies[order(alleleFrequencies[, 1], alleleFrequencies[, 2]), ]

    alleleFrequencies.gr <- GenomicRanges::GRanges(
      alleleFrequencies[, 1],
      IRanges::IRanges(alleleFrequencies[, 2], alleleFrequencies[, 2]),
      rep("*", nrow(alleleFrequencies))
    )
    S4Vectors::elementMetadata(alleleFrequencies.gr) <- alleleFrequencies[, 3:7]

    overlap <- GenomicRanges::findOverlaps(subs.data.gr, alleleFrequencies.gr)

    # FIX C1: subset BOTH tables to matched rows to keep indices aligned
    alleleFrequencies <- alleleFrequencies[S4Vectors::subjectHits(overlap), ]
    subs.data         <- subs.data[S4Vectors::queryHits(overlap), ]

    nucleotides <- c("A", "C", "G", "T")
    ref.indices <- match(subs.data[, 3], nucleotides)
    alt.indices <- match(subs.data[, 4], nucleotides)

    WT.count  <- as.numeric(vapply(
      seq_len(nrow(alleleFrequencies)),
      function(i) alleleFrequencies[i, ref.indices[i] + 2],
      numeric(1)
    ))
    mut.count <- as.numeric(vapply(
      seq_len(nrow(alleleFrequencies)),
      function(i) alleleFrequencies[i, alt.indices[i] + 2],
      numeric(1)
    ))

    # Rebuild GRanges from the now-subsetted subs.data
    combined.gr <- GenomicRanges::GRanges(
      subs.data[, 1],
      IRanges::IRanges(subs.data[, 2], subs.data[, 2]),
      rep("*", nrow(subs.data))
    )
    S4Vectors::elementMetadata(combined.gr) <- data.frame(
      WT.count = WT.count, mut.count = mut.count
    )
    combined.gr
  }
}

# ─────────────────────────────────────────────────────────────────
# FIX C2: gender-aware batch runner helper
# Reads sex from a clinical TSV; falls back to "female" with a warning
# ─────────────────────────────────────────────────────────────────
load_sample_sex <- function(clinical_file) {
  if (!file.exists(clinical_file)) {
    warning("Clinical file not found; defaulting all samples to gender='female'")
    return(NULL)
  }

  cl <- read.delim(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)

  gender_col   <- intersect(names(cl), c("gender", "Gender", "sex", "Sex"))[1]
  barcode_col  <- intersect(names(cl), c("bcr_patient_barcode", "submitter_id", "case_id"))[1]

  if (is.na(gender_col) || is.na(barcode_col)) {
    warning("Cannot find gender/barcode columns in clinical file")
    return(NULL)
  }

  sex_map <- setNames(
    tolower(cl[[gender_col]]),
    substr(cl[[barcode_col]], 1, 12)
  )
  sex_map
}

get_sample_gender <- function(sample_id, sex_map) {
  if (is.null(sex_map)) return("female")
  key <- substr(sample_id, 1, 12)
  if (key %in% names(sex_map)) {
    val <- sex_map[[key]]
    if (val %in% c("male", "m")) return("male") else return("female")
  }
  "female"
}

# ─────────────────────────────────────────────────────────────────
# FIX C3: Preserve original CCF values; only drop invalid negatives
# ─────────────────────────────────────────────────────────────────
clean_ccf <- function(ccf_raw) {
  ccf <- as.numeric(ccf_raw)
  # Negative values are invalid
  ccf[ccf < 0] <- NA_real_
  ccf
}

# ─────────────────────────────────────────────────────────────────
# FIX D5: MAF barcode matching (15-char prefix vs 28-char full barcode)
# ─────────────────────────────────────────────────────────────────
extract_sample_maf <- function(maf_df, sample_id,
                                barcode_col = "Tumor_Sample_Barcode") {
  maf_df[startsWith(maf_df[[barcode_col]], sample_id), , drop = FALSE]
}

# ─────────────────────────────────────────────────────────────────
# apply_all_patches: call once before running dpclust3p
# ─────────────────────────────────────────────────────────────────
apply_all_patches <- function() {
  if (!requireNamespace("dpclust3p", quietly = TRUE)) {
    stop("dpclust3p is not installed")
  }

  # C1: fixed GetWTandMutCount
  assignInNamespace("GetWTandMutCount",
                    make_fixed_GetWTandMutCount(),
                    ns = "dpclust3p")

  # Fixed GetCellularity: explicit 'cellularity' column lookup (D4)
  assignInNamespace("GetCellularity", function(rho_and_psi_file) {
    d <- read.table(rho_and_psi_file, header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
    if ("cellularity" %in% names(d)) return(as.numeric(d$cellularity[1]))
    if ("purity"      %in% names(d)) return(as.numeric(d$purity[1]))
    if ("rho"         %in% names(d)) return(as.numeric(d$rho[1]))
    stop("Cannot parse cellularity from file: ", rho_and_psi_file)
  }, ns = "dpclust3p")

  # Namespace qualifier patches for sortSeqlevels / queryHits / subjectHits
  if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    gdp <- get("GetDirichletProcessInfo", envir = asNamespace("dpclust3p"))
    gdp_body <- paste(deparse(body(gdp)), collapse = "\n")
    gdp_body <- gsub("sortSeqlevels\\(",  "GenomeInfoDb::sortSeqlevels(", gdp_body)
    gdp_body <- gsub("subjectHits\\(",    "S4Vectors::subjectHits(",      gdp_body)
    gdp_body <- gsub("queryHits\\(",      "S4Vectors::queryHits(",        gdp_body)
    body(gdp) <- parse(text = gdp_body)[[1]]
    assignInNamespace("GetDirichletProcessInfo", gdp, ns = "dpclust3p")
  }

  invisible(TRUE)
}

cat("Loaded convert_absolute_to_battenberg_v2.R\n")
cat("Available functions:\n")
cat("  convert_segtab_to_battenberg(segtab_file, tables_file, output_dir)\n")
cat("  convert_tables_to_cellularity(tables_file, output_dir)\n")
cat("  apply_all_patches()   # call before dpclust3p::runGetDirichletProcessInfo\n")
cat("  extract_sample_maf(maf_df, sample_id)\n")
cat("  load_sample_sex(clinical_file)\n")
cat("  get_sample_gender(sample_id, sex_map)\n")
cat("  clean_ccf(ccf_raw)\n")
cat("  clean_chromosome(chr_vec)\n")
