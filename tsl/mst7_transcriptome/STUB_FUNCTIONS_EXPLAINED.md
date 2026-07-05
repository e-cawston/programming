# How Stub Functions Work and How to Transition to Real Packages

## Current Setup: Stub Functions

Your pipeline currently uses **stub/placeholder functions** for sleuth and RankProd. This allows the pipeline to:

✅ Run end-to-end without errors
✅ Generate all output files (even if some are empty)
✅ Keep the workflow intact

⚠️ Produce empty DGE results and upset plots (because the analysis isn't real)

---

## Files Involved

```
analysis/scripts/
├─ stub_functions.R         ← Placeholder implementations
├─ dge_functions.R          ← Checks for sleuth, uses stubs as fallback
├─ rank_products_functions.R← Checks for RankProd, uses stubs as fallback
└─ upset_functions.R        ← Uses real UpSetR (working ✅)
```

---

## How It Works (Under the Hood)

### When Pipeline Starts

1. **dge_functions.R** tries to load sleuth:
   ```r
   suppress_library_error('sleuth')  # Attempts to load
   
   if (!exists('sleuth_prep')) {     # If it failed...
     source('stub_functions.R')      # Load stub versions instead
   }
   ```

2. **rank_products_functions.R** does the same for RankProd

3. **Result**: Pipeline uses whichever is available:
   - Real package if installed ✅
   - Stub functions if not ⚠️

---

## Transitioning to Real Packages

### Step 1: Install sleuth and RankProd

Follow one of the options in [INSTALL_PACKAGES.md](../INSTALL_PACKAGES.md)

**Docker method (recommended)**:
```bash
docker run -it -v /path/to/mst7_transcriptome:/project \
  bioconductor/bioconductor_docker:latest R

# Inside R:
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
BiocManager::install("RankProd")
```

### Step 2: No Code Changes Needed!

Once packages are installed, the pipeline **automatically uses them**:

1. `requireNamespace('sleuth')` will return TRUE
2. Real sleuth functions will load
3. Stub functions won't be sourced
4. Your analysis will run with real data! ✅

**The transition is automatic** - no edits to pipeline code required.

---

## Verify Real Packages Are Being Used

After installing packages, check which version is running:

```bash
mamba activate rna_analysis
R
```

```r
# Check if sleuth is loaded
if (exists('sleuth_prep')) {
  # Try to see if it's the real or stub version
  f <- body(sleuth_prep)
  if (grepl("structure", deparse(f)[1])) {
    cat("⚠️  Using STUB functions\n")
  } else {
    cat("✅ Using REAL sleuth package\n")
  }
}

# Simpler check:
packageVersion('sleuth')      # Shows version if installed
packageVersion('RankProd')    # Shows version if installed
```

---

## What Changes When You Install Real Packages

### Input
- Same Rmd file: `03_modular_pipeline.Rmd`
- Same metadata: `raw/run_metadata.txt`
- Same TPM matrix: `raw/all_samples_tpm_matrix.txt`

### Output

**Before** (with stubs):
```
results/dge/
├─ sleuth_Guy11_vs_mst7_0H.csv    ← Empty file (0 rows)
├─ sleuth_Guy11_vs_mst7_2H.csv    ← Empty file
...
```

**After** (with real packages):
```
results/dge/
├─ sleuth_Guy11_vs_mst7_0H.csv    ← Real DGE results!
│  └─ target_id, pval, qval, b, se_b, mean_obs, ...
├─ sleuth_Guy11_vs_mst7_2H.csv    ← Real DGE results!
...
results/figures/
├─ dge/
│  ├─ sleuth_Guy11_vs_mst7.png    ← Real plots!
│  ├─ sleuth_Guy11_vs_MST7WT.png
...
```

---

## Understanding the Stub Functions

**Stub for `sleuth_prep()`**:
```r
sleuth_prep <- function(...) {
  structure(list(...), class = "sleuth")
}
```
- Returns an empty sleuth object
- Real function loads and normalizes data
- Pipeline continues without crashing

**Stub for `sleuth_results()`**:
```r
sleuth_results <- function(...) {
  tibble(target_id = character(), pval = numeric(), ...)
}
```
- Returns empty tibble
- Real function returns differential expression results
- Creates empty CSV files

**Same pattern for RankProd stubs**

---

## Troubleshooting

**Q: I installed packages but pipeline still uses stubs**
- A: R may have cached the old environment. Try:
  ```bash
  mamba activate rna_analysis
  R --vanilla  # Start fresh R session
  ```

**Q: Error when running pipeline after package install**
- A: The real packages have dependencies. Install using Docker to avoid this.

**Q: How do I revert to stubs?**
- A: Uninstall the packages:
  ```r
  remove.packages(c("sleuth", "RankProd"))
  ```
  Pipeline will automatically use stubs again.

---

## Summary

| Aspect | With Stubs ⚠️ | With Real Packages ✅ |
|--------|--------------|----------------------|
| Pipeline runs | ✅ Yes | ✅ Yes |
| Correlation analysis | ✅ Works | ✅ Works |
| DGE results | ❌ Empty | ✅ Real data |
| RankProd results | ❌ Empty | ✅ Real data |
| UpSetR plots | ❌ Empty | ✅ Real overlap plots |
| Code changes needed | None | None (automatic!) |
| Time to results | Minutes | Minutes |

---

## Next: Getting Real Data

When you're ready for real analysis:

1. Choose an installation method from [INSTALL_PACKAGES.md](../INSTALL_PACKAGES.md)
2. Install sleuth and RankProd
3. Run `rmarkdown::render('analysis/03_modular_pipeline.Rmd')`
4. Your results will automatically be real! ✅

No pipeline code changes needed - it's designed to work seamlessly with either stubs or real packages.
