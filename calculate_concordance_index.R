#!/usr/bin/env Rscript

# Calculate Concordance Index for Driver-Only Gene Progression Analysis
# Compares gene ranking concordance between TracerX and TCGA-LUAD

# Set working directory
setwd('c:/CloneHD_benchmarking')

# Read driver-only gene progression data
driver_data <- read.csv("_dpclust_publish/TCGA/outputs_old/validation_gene_progression/shared_genes_progression_comparison_driver_only.tsv", sep="\t")

cat("=== Driver-Only Gene Progression Concordance Analysis ===\n\n")

# Extract data for each cohort
n_genes <- nrow(driver_data)

# Calculate Spearman correlation
spearman_result <- cor.test(driver_data$rank_progression_tracerx, 
                             driver_data$rank_progression_luad, 
                             method = "spearman")

# Calculate Kendall tau correlation
kendall_result <- cor.test(driver_data$rank_progression_tracerx, 
                            driver_data$rank_progression_luad, 
                            method = "kendall")

# Calculate Pearson correlation on progression scores
pearson_result <- cor.test(driver_data$progression_score_tracerx, 
                            driver_data$progression_score_luad, 
                            method = "pearson")

# Calculate concordance index (Harrell's C-index concept)
# Proportion of concordant pairs
concordant_pairs <- 0
discordant_pairs <- 0
tied_pairs <- 0

for (i in 1:(n_genes-1)) {
  for (j in (i+1):n_genes) {
    # Compare ranking differences in both cohorts
    tx_diff <- sign(driver_data$rank_progression_tracerx[i] - driver_data$rank_progression_tracerx[j])
    luad_diff <- sign(driver_data$rank_progression_luad[i] - driver_data$rank_progression_luad[j])
    
    if (tx_diff == luad_diff && tx_diff != 0) {
      concordant_pairs <- concordant_pairs + 1
    } else if (tx_diff != luad_diff && tx_diff != 0 && luad_diff != 0) {
      discordant_pairs <- discordant_pairs + 1
    } else if (tx_diff == 0 || luad_diff == 0) {
      tied_pairs <- tied_pairs + 1
    }
  }
}

total_pairs <- concordant_pairs + discordant_pairs
concordance_index <- concordant_pairs / total_pairs

cat("RANKING CONCORDANCE METRICS:\n")
cat("-----------------------------\n")
cat(sprintf("Genes analyzed: %d\n", n_genes))
cat(sprintf("Concordant pairs: %d\n", concordant_pairs))
cat(sprintf("Discordant pairs: %d\n", discordant_pairs))
cat(sprintf("Total comparable pairs: %d\n\n", total_pairs))

cat("CORRELATION COEFFICIENTS:\n")
cat("------------------------\n")
cat(sprintf("Spearman ρ: %.4f (p-value: %.2e)\n", spearman_result$estimate, spearman_result$p.value))
cat(sprintf("Kendall τ: %.4f (p-value: %.2e)\n", kendall_result$estimate, kendall_result$p.value))
cat(sprintf("Pearson r: %.4f (p-value: %.2e)\n\n", pearson_result$estimate, pearson_result$p.value))

cat("CONCORDANCE INDEX:\n")
cat("------------------\n")
cat(sprintf("Concordance Index: %.4f (%.2f%%)\n", concordance_index, concordance_index * 100))
cat(sprintf("Discordance Index: %.4f (%.2f%%)\n\n", 1 - concordance_index, (1 - concordance_index) * 100))

# Top genes by progression score agreement
cat("TOP 15 GENES BY PROGRESSION SCORE AGREEMENT:\n")
cat("---------------------------------------------\n")

score_agreement <- data.frame(
  gene = driver_data$gene,
  progression_score_tracerx = driver_data$progression_score_tracerx,
  progression_score_luad = driver_data$progression_score_luad,
  score_abs_diff = abs(driver_data$progression_score_tracerx - driver_data$progression_score_luad),
  rank_progression_tracerx = driver_data$rank_progression_tracerx,
  rank_progression_luad = driver_data$rank_progression_luad,
  rank_abs_diff = abs(driver_data$rank_progression_tracerx - driver_data$rank_progression_luad)
)

# Sort by score_abs_diff
score_agreement <- score_agreement[order(score_agreement$score_abs_diff), ]

print(head(score_agreement, 15))

cat("\n\nBOTTOM 5 GENES (MOST DISCORDANT):\n")
cat("----------------------------------\n")
print(tail(score_agreement, 5))

# Summary statistics
cat("\n\nDISCORDANCE DISTRIBUTION:\n")
cat("-------------------------\n")
cat(sprintf("Mean rank difference: %.2f\n", mean(score_agreement$rank_abs_diff)))
cat(sprintf("Median rank difference: %.2f\n", median(score_agreement$rank_abs_diff)))
cat(sprintf("Max rank difference: %d\n", max(score_agreement$rank_abs_diff)))

cat("\n=== Analysis Complete ===\n")
