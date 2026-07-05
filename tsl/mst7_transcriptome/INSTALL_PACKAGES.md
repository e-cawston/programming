# Installing sleuth and RankProd - Complete Setup Guide

## Status

Your pipeline is now set up with:
- ✅ **Working**: ggplot2, dplyr, readr, UpSetR, correlation analysis
- ⚠️  **Stub functions**: sleuth, RankProd (pipeline runs but features limited)

The stub functions allow your pipeline to run end-to-end without errors. Once you install the real packages, just replace the stubs.

---

## Option 1: Use Docker (RECOMMENDED)

Docker images come with all system dependencies pre-configured.

```bash
# Pull an R image with Bioconductor pre-installed
docker pull bioconductor/bioconductor_docker:latest

# Or use a specific version with sleuth already included
docker pull rocker/bioconductor:4.3

# Run interactively
docker run -it -v /path/to/mst7_transcriptome:/project \
  bioconductor/bioconductor_docker:latest R
```

Inside Docker:
```r
BiocManager::install(c("sleuth", "RankProd"))
# Both should install without issues
```

---

## Option 2: Use Rocker R Docker with devtools

```bash
# Use rocker/r-ver with development tools
docker run -it \
  -v /path/to/mst7_transcriptome:/project \
  -w /project \
  rocker/r-ver:latest bash

# Inside container:
apt-get update && apt-get install -y \
  libgmp-dev libmpfr-dev libcurl4-openssl-dev

R
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
BiocManager::install("RankProd")
```

---

## Option 3: Install in fresh conda environment (requires system setup)

```bash
# Remove current environment
mamba env remove -n rna_analysis

# Create new one with system dev tools
mamba create -n rna_analysis -c conda-forge \
  r-base=4.5 \
  r-tidyverse \
  r-ggplot2 \
  r-here \
  r-reshape2 \
  openblas \
  cmake \
  gcc \
  gxx

# Activate and try again
mamba activate rna_analysis
```

Then in R:
```r
install.packages("sleuth")
install.packages("RankProd")
```

---

## Option 4: Manual compilation (advanced)

If you need to compile from source:

1. Install system dependencies (requires sudo or system admin):
   ```bash
   # Ubuntu/Debian
   sudo apt-get install libgmp-dev libmpfr-dev libhdf5-dev libcurl4-openssl-dev
   
   # macOS
   brew install gmp mpfr hdf5 curl
   ```

2. Install in R:
   ```r
   BiocManager::install(c("sleuth", "RankProd"))
   ```

---

## Option 5: Use Singularity (HPC clusters)

```bash
# Build from Docker image (if you have Docker)
singularity build bioc.sif docker://bioconductor/bioconductor_docker:latest

# Use the image
singularity exec bioc.sif R
```

---

## How the Stubs Work

Currently, your pipeline has **stub functions** for `sleuth_*` and `RankProd` functions. These:

- ✅ Allow the pipeline to run without errors
- ✅ Won't crash on function calls
- ⚠️  Return empty results (not real analysis)

When you install the real packages, the stubs will be automatically replaced and your real analysis will run.

---

## Verify Real Packages Are Loaded

To check if the real packages are installed:

```r
mamba activate rna_analysis
R

# This will tell you which is loaded
cat("sleuth available:", requireNamespace("sleuth", quietly=TRUE), "\n")
cat("RankProd available:", requireNamespace("RankProd", quietly=TRUE), "\n")

# If TRUE: real package is loaded ✅
# If FALSE: stub functions are being used ⚠️
```

---

## Next Steps

1. **For testing**: Run pipeline as-is with stubs (it will complete without errors)
2. **For real analysis**: Use Option 1 (Docker) - easiest and most reliable
3. **For HPC**: Use Option 5 (Singularity)

---

## Troubleshooting

**"sleuth/RankProd not found in repository"**
- Bioconductor version mismatch - Docker images solve this

**"Compilation failed for ..."**
- Missing system headers - Docker or Option 1 recommended

**"Could not open file libjasper"**
- System library missing - use Docker

**Pipeline runs but shows "stub functions" warnings**
- This is normal! Install real packages using one of the options above.

---

## Running the Pipeline

```bash
# Use with stubs (gets through all sections, no real analysis)
mamba activate rna_analysis
cd /path/to/mst7_transcriptome
R -e "rmarkdown::render('analysis/03_modular_pipeline.Rmd')"

# Check results in:
# - results/figures/correlation/*.png
# - results/tables/correlation/*.csv
# - results/dge/*.csv (empty with stubs)
```
