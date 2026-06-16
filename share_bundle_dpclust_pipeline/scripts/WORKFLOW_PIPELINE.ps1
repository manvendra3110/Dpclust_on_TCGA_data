# End-to-end workflow runner for DPclust ABSOLUTE -> old-v100 to TRACERx comparisons
# Usage:
#   powershell -ExecutionPolicy Bypass -File .\WORKFLOW_PIPELINE.ps1

# Optional: Set $runAbsoluteConversion to $true to include Step 0
$runAbsoluteConversion = $false  # Set to $true if starting from ABSOLUTE mastercalls

$ErrorActionPreference = 'Stop'

$R = 'C:/Program Files/R/R-4.5.3/bin/Rscript.exe'
$Root = 'C:/CloneHD_benchmarking/DPclust'
$Tcga = 'C:/CloneHD_benchmarking/DPclust/TCGA'
$Tracerx = 'C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda'
$Mc3Maf = 'C:/CloneHD_benchmarking/DPclust/TCGA/mc3.v0.2.8.PUBLIC.maf/mc3.v0.2.8.PUBLIC.maf'

$OutputsRoot = "$Tcga/outputs_old"
$OldBatchOut = "$Tcga/outputs_old/batch_output"
$DistOut = "$Tcga/outputs_old/validation_ccf_distribution"
$GeneOut = "$Tcga/outputs_old/validation_gene_progression"
$FinalOut = "$Tcga/outputs_old/validation_final"

Write-Host 'Resetting outputs_old for a clean rerun...' -ForegroundColor Cyan
if (Test-Path $OutputsRoot) {
  Remove-Item $OutputsRoot -Recurse -Force
}

Write-Host 'Creating output folders...' -ForegroundColor Cyan
New-Item -ItemType Directory -Force -Path $OldBatchOut | Out-Null
New-Item -ItemType Directory -Force -Path $DistOut | Out-Null
New-Item -ItemType Directory -Force -Path $GeneOut | Out-Null
New-Item -ItemType Directory -Force -Path $FinalOut | Out-Null

Set-Location $Root

# Optional Step 0: ABSOLUTE conversion (if starting from mastercalls)
if ($runAbsoluteConversion) {
  Write-Host 'Step 0/4: Converting ABSOLUTE mastercalls to Battenberg format...' -ForegroundColor Yellow
  & $R 'run_absolute_to_dpclust_batch.R' `
    --segtab='TCGA_mastercalls.abs_segtabs.fixed.txt' `
    --tables='TCGA_mastercalls.abs_tables_JSedit.fixed.txt' `
    --output_dir='C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_luad' `
    --sample_list='C:/CloneHD_benchmarking/DPclust/TCGA/tcga_luad_samples.txt'
  if ($LASTEXITCODE -ne 0) { throw "Step 0 failed with exit code $LASTEXITCODE" }
  Write-Host ''
}

Write-Host 'Step 1/3: Running old dpclust3p batch export...' -ForegroundColor Yellow
& $R 'run_old_dpclust3p_batch_from_existing.R' `
  --input_root='C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_luad' `
  --output_root='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --sample_list='C:/CloneHD_benchmarking/DPclust/TCGA/tcga_luad_samples.txt' `
  --lib='C:/CloneHD_benchmarking/DPclust/.r_libs_dpclust3p_v100' `
  --dpclust3p_ref='dpclust3p-v1.0.0'
if ($LASTEXITCODE -ne 0) { throw "Step 1 failed with exit code $LASTEXITCODE" }

Write-Host 'Step 2/3: Running CCF distribution comparison...' -ForegroundColor Yellow
& $R 'compare_luad_old_v100_to_tracerx_distribution.R' `
  --old_root='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --tracerx_rda='C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda' `
  --out_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/validation_ccf_distribution'
if ($LASTEXITCODE -ne 0) { throw "Step 2 failed with exit code $LASTEXITCODE" }

Write-Host 'Step 3/3: Running gene progression comparison...' -ForegroundColor Yellow
& $R 'compare_luad_old_v100_to_tracerx_gene_progression.R' `
  --old_root='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --luad_fast_root='C:/CloneHD_benchmarking/DPclust/TCGA/dpclust_batch_luad' `
  --tracerx_rda='C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda' `
  --out_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/validation_gene_progression'
if ($LASTEXITCODE -ne 0) { throw "Step 3 failed with exit code $LASTEXITCODE" }

Write-Host 'Step 4/4: Running consolidated TCGA vs TRACERx comparison...' -ForegroundColor Yellow
& $R 'compare_tcga_dpclust3p_vs_tracerx.R' `
  --tracerx_rda='C:/Users/ManvendraSingh/Downloads/TRACERx_NEJM_2017.rda' `
  --maf_path='C:/CloneHD_benchmarking/DPclust/TCGA/mc3.v0.2.8.PUBLIC.maf/mc3.v0.2.8.PUBLIC.maf' `
  --ccf_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/batch_output' `
  --sample_list='C:/CloneHD_benchmarking/DPclust/TCGA/tcga_luad_samples.txt' `
  --out_dir='C:/CloneHD_benchmarking/DPclust/TCGA/outputs_old/validation_final'
if ($LASTEXITCODE -ne 0) { throw "Step 4 failed with exit code $LASTEXITCODE" }

Write-Host ''
Write-Host 'Workflow completed successfully.' -ForegroundColor Green
Write-Host "Outputs are in: $Tcga/outputs_old"
