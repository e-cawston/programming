# Infection Phosphoproteome Project — Copilot Instructions

## Purpose
This file tells GitHub Copilot how to behave inside this project.  
Copilot must use these instructions whenever editing or generating:

- `analysis/scripts/phosphoproteome_cleanup.r`
- `analysis/scripts/phosphoproteome_analysis.r`
- `analysis/infection_phosphoproteome.rmd`

Copilot may also use `chat.md` for additional context.

---

## Project structure

- `raw/pmk1/combined_modified_peptides.tsv`
- `raw/mst7/combined_modified_peptides.tsv`
- `analysis/infection_phosphoproteome.rmd`
- `analysis/scripts/phosphoproteome_cleanup.r`
- `analysis/scripts/phosphoproteome_analysis.r`
- `chat.md`
- `instructions.md`

---

## Dataset structure

Copilot must assume the TSV has at least:

- **Peptide Sequence**
- **Modified Sequence**
- **Assigned Modifications**
- **Protein**, **Protein ID**, **Gene**
- Spectral count columns like:
  - `guy11_0h_1 Spectral Count`, …, `mst7_6h_3 Spectral Count`
- Intensity columns
- MaxLFQ columns like:
  - `guy11_0h_1 MaxLFQ Intensity`, …, `mst7_6h_3 MaxLFQ Intensity`

Use:

- `Assigned Modifications` for phosphorylation detection (mass shift ≈ `79.9663`).
- `Peptide Sequence` as peptide ID.
- `Protein ID` (or `Protein` if needed) as protein ID.

---

## Cleanup script (`analysis/scripts/phosphoproteome_cleanup.r`)

Copilot must implement:

### `load_phospho_data(path)`
- Read TSV into a tibble.
- Return raw data.

### `filter_phospho_peptides_all(df, mod_col = "Assigned Modifications", mass_shift = 79, tolerance = 2)`
- Keep peptides with at least one phosphorylation (mass shift within `mass_shift ± tolerance`).
- Allow other modifications.

### `filter_phospho_peptides_only(df, mod_col = "Assigned Modifications", mass_shift = 79, tolerance = 2)`
- Keep peptides with phosphorylation and **no other modifications**.

### `filter_nonzero_spectral_counts(df, spec_cols)`
- Remove rows where **all** spectral count columns are zero or `NA`.

### `cleanup_phospho_dataset(path, mod_col, spec_cols, mass_shift = 79, tolerance = 2)`
- Load data.
- Create:
  - `phospho_all`
  - `phospho_only`
- Apply spectral count filtering to both.
- Return:

```r
list(
  raw = df,
  phospho_all = phospho_all,
  phospho_only = phospho_only
)
## renv environment

This project uses **renv**.  
Copilot must assume:

- All required R packages are already installed inside the renv environment.
- Copilot must **not** suggest `install.packages()`, `renv::install()`, or any installation commands.
- Copilot must write code that simply loads packages with `library()`.

If a package appears missing, Copilot must:
- Add a comment indicating the package should be added to renv by the user,
- But **must not** generate installation commands.
