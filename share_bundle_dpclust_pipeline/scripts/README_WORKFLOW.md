# DPclust Workflow: ABSOLUTE to TRACERx Comparison

This document captures the full workflow completed in this project:
1. **ABSOLUTE conversion**: TCGA ABSOLUTE mastercalls to Battenberg and cellularity files.
2. Batch execution of old dpclust3p v1.0.0 on TCGA LUAD.
3. Distribution-level comparison against TRACERx.
4. Gene-level progression order comparison against TRACERx.
5. Consolidated cohort-level comparison with plots and per-gene summaries.
6. Clean rerun into `TCGA/outputs_old` only.

## 1) Prerequisites

- Windows + PowerShell
- R 4.5.3 available as `Rscript.exe`
- TRACERx RDA file:
  `C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda`
- ABSOLUTE mastercalls files (segtabs + tables)
- Workspace root:
  `C:/CloneHD_benchmarking/DPclust`

## 2) Key Scripts

### ABSOLUTE Conversion (NEW)
- `convert_absolute_to_battenberg_v2.R`
  - Functions for converting ABSOLUTE segtabs/tables to dpclust3p-compatible format.
  - Includes all bug fixes (D1-D5, N1-N2, C1-C3).
- `run_absolute_to_dpclust_batch.R`
  - Batch runner: converts all ABSOLUTE samples at once.

### Core Workflow
- `DP_clust.r`
  - Core preprocessing and dpclust3p input generation (single sample).
- `run_old_dpclust3p_batch_from_existing.R`
  - Runs old dpclust3p v1.0.0 in batch mode on LUAD cohort.
- `compare_luad_old_v100_to_tracerx_distribution.R`
  - Cohort CCF distribution comparison (KS test and summaries).
- `compare_luad_old_v100_to_tracerx_gene_progression.R`
  - Gene mapping and progression rank comparison.
- `compare_tcga_dpclust3p_vs_tracerx.R`
  - Consolidated final comparison: distribution, per-gene CCF, and clonality plots/tables.
- `diag_dpclust3p.R`
  - Diagnostics helper for dpclust3p outputs/issues.

## 3) Output Structure

Organized under `TCGA/`:

- `outputs_old/batch_output/`
  - Old v1.0.0 per-sample results.
  - Also contains `old_v100_processing_summary.tsv` and `old_v100_failures.tsv`.
- `outputs_old/validation_ccf_distribution/`
  - Distribution comparison outputs.
- `outputs_old/validation_gene_progression/`
  - Gene progression comparison outputs.
- `outputs_old/validation_final/`
  - Final consolidated comparison tables and plots.

## 4) End-to-End Commands

For a clean rerun, use the PowerShell pipeline below. It deletes `TCGA/outputs_old` first and rebuilds all outputs from scratch.

### Cohort Selection (LUAD vs LUAD+LUSC)

The workflow script is now cohort-configurable via variables at the top of `WORKFLOW_PIPELINE.ps1`.

- Default (LUAD-only):
  - `$CohortName = 'LUAD'`
  - `$SampleList = "$Tcga/tcga_luad_samples.txt"`
  - `$DpclustInputRoot = "$Tcga/dpclust_batch_luad"`
- LUAD+LUSC:
  - `$CohortName = 'NSCLC'`
  - `$SampleList = "$Tcga/tcga_nsclc_samples.txt"`
  - `$DpclustInputRoot = "$Tcga/dpclust_batch_nsclc"`

Note: the pipeline requires a sample list file. To run LUAD+LUSC, provide `tcga_nsclc_samples.txt` (15-char sample IDs), otherwise the run remains LUAD-only.

```powershell
powershell -ExecutionPolicy Bypass -File .\WORKFLOW_PIPELINE.ps1
```

Use R executable:

```powershell
$R = 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe'
Set-Location 'C:/CloneHD_benchmarking/DPclust'
```

### Step 0: ABSOLUTE conversion (optional, if starting from mastercalls)

If you have ABSOLUTE segtabs and tables files:

```powershell
& $R 'run_absolute_to_dpclust_batch.R' `
  --segtab='TCGA_mastercalls.abs_segtabs.fixed.txt' `
  --tables='TCGA_mastercalls.abs_tables_JSedit.fixed.txt' `
  --output_dir='C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_luad' `
  --sample_list='C:/CloneHD_benchmarking/DPclust/TCGA/tcga_luad_samples.txt'
```

This converts ABSOLUTE data with all bug fixes:
- D1-D5: Data format corrections (chr types, coordinates, CCF semantics)
- N1-N2: Chromosome robustness, genome doubling baseline adjustment
- C1-C3: dpclust3p function patches (GetWTandMutCount, cellularity parsing, preserve original CCF)

### Step A: Old dpclust batch generation

```powershell
& $R 'run_old_dpclust3p_batch_from_existing.R' `
  --input_root='C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_luad' `
  --output_root='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --sample_list='C:/CloneHD_benchmarking/DPclust/TCGA/tcga_luad_samples.txt' `
  --lib='C:/CloneHD_benchmarking/DPclust/.r_libs_dpclust3p_v100' `
  --dpclust3p_ref='dpclust3p-v1.0.0'
```

### Step B: Distribution comparison (old vs TRACERx)

```powershell
& $R 'compare_luad_old_v100_to_tracerx_distribution.R' `
  --old_root='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --tracerx_rda='C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda' `
  --out_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/validation_ccf_distribution'
```

### Step C: Gene progression comparison

```powershell
& $R 'compare_luad_old_v100_to_tracerx_gene_progression.R' `
  --old_root='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --luad_fast_root='C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_luad' `
  --tracerx_rda='C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda' `
  --driver_tsv='C:/CloneHD_benchmarking/DPclust/cosmic_lung_nsclc_driver_genes.tsv' `
  --out_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/validation_gene_progression'
```

### Step D: Final consolidated comparison

```powershell
& $R 'compare_tcga_dpclust3p_vs_tracerx.R' `
  --tracerx_rda='C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda' `
  --maf_path='C:/CloneHD_benchmarking/DPclust/TCGA/mc3.v0.2.8.PUBLIC.maf/mc3.v0.2.8.PUBLIC.maf' `
  --ccf_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --sample_list='C:/CloneHD_benchmarking/DPclust/TCGA/tcga_luad_samples.txt' `
  --out_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/validation_final'
```

## 5) Result Files Produced

### Distribution validation

- `tracerx_vs_tcga_luad_ccf_distribution_ks.tsv`
- `tracerx_vs_tcga_luad_ccf_distribution_summary.tsv`
- `tracerx_vs_tcga_luad_id_overlap.tsv`

### Gene progression validation

- `tracerx_gene_progression.tsv`
- `luad_gene_progression.tsv`
- `shared_genes_progression_comparison.tsv`
- `progression_rank_agreement_summary.tsv`
- `top20_overlap_genes.tsv`
- `tracerx_only_genes.tsv`
- `luad_only_genes.tsv`

### Final consolidated comparison

- `ccf_distribution_summary.tsv`
- `ccf_distribution_ks_test.tsv`
- `ccf_distribution_density.pdf`
- `per_gene_ccf_comparison.tsv`
- `per_gene_ccf_scatter.pdf`
- `clonality_per_driver_gene.tsv`
- `clonality_barplot.pdf`

## 6) Notes

- Old dpclust3p v1.0.0 can produce values > 1 in raw subclonal fraction context.
- Workflow uses validated CCF handling for fair TRACERx comparisons.
- The clean rerun pipeline now resets `TCGA/outputs_old` before rebuilding outputs.
