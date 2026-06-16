# run_absolute_to_dpclust_batch.R
# Batch converter: ABSOLUTE mastercalls → dpclust3p-compatible Battenberg/cellularity files
# Integrates all bug fixes from convert_absolute_to_battenberg_v2.R

library(optparse)

option_list <- list(
  make_option(c("--segtab"), type = "character", default = NULL,
              help = "TCGA ABSOLUTE segtabs file (e.g., TCGA_mastercalls.abs_segtabs.fixed.txt)"),
  make_option(c("--tables"), type = "character", default = NULL,
              help = "TCGA ABSOLUTE tables file (e.g., TCGA_mastercalls.abs_tables_JSedit.fixed.txt)"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory for Battenberg and cellularity files"),
  make_option(c("--sample_list"), type = "character", default = NULL,
              help = "Optional: sample list file (one sample ID per line) to subset conversion"),
  make_option(c("--clinical"), type = "character", default = NULL,
              help = "Optional: clinical TSV with gender info"),
  make_option(c("--help"), action = "store_true", default = FALSE,
              help = "Show this help message")
)

parser <- OptionParser(option_list = option_list, description = "
Convert TCGA ABSOLUTE mastercalls to dpclust3p-compatible format
with all known bug fixes (v2 converter).

Example:
  Rscript run_absolute_to_dpclust_batch.R \\
    --segtab='TCGA_mastercalls.abs_segtabs.fixed.txt' \\
    --tables='TCGA_mastercalls.abs_tables_JSedit.fixed.txt' \\
    --output_dir='TCGA/dpclust_batch_luad' \\
    --sample_list='TCGA/tcga_luad_samples.txt'
")

args <- parse_args(parser, positional_arguments = 0)

if (args$help || is.null(args$segtab) || is.null(args$tables) || is.null(args$output_dir)) {
  print_help(parser)
  quit(save = "no", status = 0)
}

# Source the converter (must be in same directory or in R search path)
converter_script <- "convert_absolute_to_battenberg_v2.R"
if (!file.exists(converter_script)) {
  stop("Cannot find ", converter_script, " in current directory")
}
source(converter_script)

cat("════════════════════════════════════════════════════════════\n")
cat("ABSOLUTE → dpclust3p Batch Converter\n")
cat("════════════════════════════════════════════════════════════\n\n")

# Create output directory structure
cat("Creating output directory structure...\n")
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}

# Load sample list if provided (to subset ABSOLUTE data)
sample_list <- NULL
if (!is.null(args$sample_list)) {
  if (!file.exists(args$sample_list)) {
    cat("Warning: sample list not found:", args$sample_list, "\n")
  } else {
    sample_list <- read.table(args$sample_list, stringsAsFactors = FALSE)[[1]]
    cat(sprintf("Loaded %d samples from sample list\n", length(sample_list)))
  }
}

# Load clinical data if provided (for gender info)
sex_map <- NULL
if (!is.null(args$clinical)) {
  cat("Loading clinical data...\n")
  sex_map <- load_sample_sex(args$clinical)
  if (!is.null(sex_map)) {
    cat(sprintf("Loaded gender info for %d patients\n", length(sex_map)))
  }
}

cat("\n")

# Step 1: Convert ABSOLUTE segtabs to Battenberg format
cat("Step 1/2: Converting ABSOLUTE segtabs to Battenberg format\n")
cat("────────────────────────────────────────────────────────────\n")

seg_df <- read.delim(args$segtab, stringsAsFactors = FALSE, check.names = FALSE)
if (!is.null(sample_list)) {
  seg_df <- seg_df[seg_df$Sample %in% sample_list, ]
  cat(sprintf("Filtered to %d samples from sample list\n", nrow(seg_df)))
}

convert_segtab_to_battenberg(args$segtab, args$tables, args$output_dir)

cat("\n")

# Step 2: Convert ABSOLUTE tables to cellularity files
cat("Step 2/2: Converting ABSOLUTE tables to cellularity files\n")
cat("────────────────────────────────────────────────────────────\n")

convert_tables_to_cellularity(args$tables, args$output_dir)

cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat("Conversion complete!\n")
cat("Output directory:", args$output_dir, "\n")
cat("\n")

# Report gender distribution if available
if (!is.null(sex_map)) {
  samples_with_gender <- unique(seg_df$Sample)[unique(seg_df$Sample) %in% names(sex_map)]
  male_count <- sum(sex_map[samples_with_gender] %in% c("male", "m"))
  female_count <- length(samples_with_gender) - male_count
  cat("Gender distribution in converted samples:\n")
  cat(sprintf("  Male:   %d\n", male_count))
  cat(sprintf("  Female: %d\n", female_count))
  cat("\n")
}

cat("Next steps:\n")
cat("1. Run dpclust3p on the output files:\n")
cat("   - Loci/allele files are in: ", args$output_dir, "/\n", sep = "")
cat("   - Use DP_clust.r with output_dir='", args$output_dir, "'\n", sep = "")
cat("2. Call apply_all_patches() before running dpclust3p::runGetDirichletProcessInfo\n")
cat("════════════════════════════════════════════════════════════\n")
