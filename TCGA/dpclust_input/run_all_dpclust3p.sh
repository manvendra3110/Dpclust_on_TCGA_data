#!/usr/bin/env bash
# Run dpclust3p on all samples
set -euo pipefail

echo "Processing TCGA-05-4244-01..."
Rscript C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_input\TCGA-05-4244-01/run_dpclust3p_TCGA-05-4244-01.R

