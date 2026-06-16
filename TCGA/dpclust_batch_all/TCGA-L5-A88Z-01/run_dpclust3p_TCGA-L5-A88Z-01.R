#!/usr/bin/env Rscript
# Auto-generated R script to run dpclust3p on TCGA-L5-A88Z-01
library(dpclust3p)

sample <- "TCGA-L5-A88Z-01"
output_dir <- "C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_all\TCGA-L5-A88Z-01"

loci_file <- file.path(output_dir, paste0(sample, "_loci.txt"))
allele_file <- file.path(output_dir, paste0(sample, "_alleleFrequencies.txt"))
cellularity_file <- file.path(output_dir, paste0(sample, ".cellularity_ploidy"))
subclone_file <- file.path(output_dir, paste0(sample, ".battenberg"))
output_file <- file.path(output_dir, paste0(sample, "_dpclust_input.txt"))

cat(sprintf("Running dpclust3p for %s\n", sample))

runGetDirichletProcessInfo(
    loci_file              = loci_file,
    allele_frequencies_file = allele_file,
    cellularity_file       = cellularity_file,
    subclone_file          = subclone_file,
    gender                 = "female",
    SNP.phase.file         = "NA",
    mut.phase.file         = "NA",
    output_file            = output_file
)

cat(sprintf("Done! Output: %s\n", output_file))

# Read and summarize
results <- read.delim(output_file)
cat(sprintf("Mutations: %d\n", nrow(results)))
cat(sprintf("Clonal (CCF>=0.85): %d (%.1f%%)\n",
    sum(results$subclonal.fraction >= 0.85, na.rm=TRUE),
    100 * mean(results$subclonal.fraction >= 0.85, na.rm=TRUE)))
cat(sprintf("Mean CCF: %.3f\n", mean(results$subclonal.fraction, na.rm=TRUE)))
