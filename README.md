# smallarea-validation-design

Code for:

> Chen Y, Gao H, Ma TF, Olatosi B, Baker P, Zhang J. *Evaluation design determines model choice in small-area surveillance: comparing spatiotemporal models for county-level substance use disorder-related hospital encounters among people living with HIV in South Carolina, 2005–2023.*

The pipeline does four things:

1. Builds a county-year panel from person-level surveillance and hospital discharge records.
2. Fits crude estimates and four Bayesian small-area models in R-INLA:
   - a binomial GLMM;
   - a beta-binomial GLMM;
   - BYM2 + RW2;
   - BYM2 + RW2 with a county–year interaction.
3. Compares the five approaches under three holdout designs, one for each surveillance task:

   | Design | Held out | Surveillance task |
   |---|---|---|
   | A. Random county-year | 10 folds, stratified by denominator | fill in a suppressed county-year |
   | B. Forward in time | each year 2019–2023, trained on earlier years | release estimates for a new year |
   | C. Leave one county out | each county's full series | estimate a county with no data |

4. Applies a prespecified selection rule. Among the four models, it chooses the simplest one that meets all three criteria:
   - Design A RMSE within 0.25 percentage points of the best model;
   - a maximum calibration gap of at most 1.5 percentage points;
   - residual Moran's I significant (5% FDR) in no more than 2 of 19 years.

The South Carolina data are restricted and are **not** included. `run_toy.R` generates synthetic data with the same file layout and runs the complete pipeline on it.

## Requirements

- R ≥ 4.3
- These CRAN packages:

  ```r
  install.packages(c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
                     "tidyselect", "haven", "readxl", "sf", "spdep", "tigris", "ggplot2"))
  ```

- INLA, which is not on CRAN:

  ```r
  install.packages("INLA", repos = c(getOption("repos"),
                   INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
  ```

The code runs on macOS, Windows and Linux. County boundaries are downloaded with `tigris`, so the first run needs an internet connection. Offline, download the 2020 TIGER/Line county shapefile and set the environment variable `SAE_SHAPEFILE` to its path.

## Quick start

**RStudio:** open `smallarea-validation-design.Rproj`, then run

```r
source("run_toy.R")
```

**Terminal:** run this from the repository folder:

```bash
Rscript run_toy.R
```

The quick run takes about 4 minutes on a laptop (MacBook Air, Apple M3). It uses 5 folds and skips Design C. To run the full protocol used in the paper (10 folds, Designs A–C, 9,999 permutations), set `TOY_FULL=1` before running.

What the quick run does:

1. **Simulates** a person-level HIV case file and a hospital discharge file (`R/simulate_toy_data.R`). They have the same columns and code formats as the restricted sources. The simulation deliberately includes:
   - persons under 18 or with a missing county, who must be excluded;
   - tobacco, alcohol and unlisted dependence codes, which must not be counted.

   The true county-year proportion is saved in `toy_data/truth_county_year.csv`.
2. **Runs the pipeline** (`R/01_` to `R/04_`): cohort, ICD matching, panel, adjacency, Moran's I, models, holdout designs, selection rule, tables and figures.
3. **Scores each approach against the truth** and saves the result to `output_toy/toy_rmse_vs_truth.csv`. This check is not possible with real data.

All results are written to `output_toy/`. The hardware and software versions are recorded in `output_toy/logs/session_report_toy.txt`.

## Repository layout

```
run_toy.R                   complete pipeline on synthetic data (start here)
run_all.R                   analysis of the restricted data + 12 sensitivity variants
R/00_config.R               all options: paths, code list, priors, designs, thresholds
R/codelist.R                reads and checks the ICD code list
R/01_build_panel.R          cohort, ICD-9-CM/ICD-10-CM matching, county-year panel
R/02_descriptives_moran.R   adjacency, Table 1, Moran's I, Figure 1
R/03_models.R               models, holdout designs, metrics, selection rule
R/04_tables_figures.R       Tables 2–3, Figures 2–3, supplementary figures
R/simulate_toy_data.R       synthetic data generator with known truth
R/session_report.R          records R/package versions, CPU, memory, threads
data/ICD_Code_List_SUD.xlsx ICD code list (Supplementary Table S1)
tools/pi_betabinomial_M2.R  optional: beta-binomial predictive intervals for M2
```

## Using your own data

**Person-level files.** Point `cfg$file_cases` and `cfg$file_ub` in `R/00_config.R` to your files. Both SAS (`.sas7bdat`) and CSV files are read.

- The case file has one row per person, with these columns:
  - `RFA_ID`
  - `AGE_AT_HIV_DX`
  - `DATE_OF_HIV_DX`, `DATE_OF_AIDS_DX` and `DEATH_HARS` (format `YYYYMM`)
  - `COUNTY_AT_HIV_DX` and `CURRENT_COUNTY`
- The discharge file has one row per discharge, with these columns:
  - `RFA_ID`, `DISYEAR`, `DISMTH` and `TIME_DXDATE_DISD`
  - diagnosis fields `PDIAG`, `ADM_DIAG` and `SDIAG1`–`SDIAG14`
  - diagnosis fields `ADM_DIAG10` and `SDIAG10_1`–`SDIAG10_14`

  Codes are assigned to ICD-9-CM or ICD-10-CM by discharge date. Diagnosis fields that are absent are skipped.

**County-year panel.** If you already have a panel with columns `fips`, `county`, `year`, `deno` and `nume`, save it as `output/panel_core_surveillance.csv` and start from `R/02_descriptives_moran.R`.

**Another state.** Replace the county crosswalk and the boundary source in `R/00_config.R`.

## Computing environment of the published analysis

- **Machine:** MacBook Air, Apple M3 (8 cores: 4 performance, 4 efficiency), 16 GB RAM, macOS 27.2.
- **R and BLAS:** R 4.4.2 (aarch64-apple-darwin20) with Apple Accelerate BLAS.
- **Packages:**
  - INLA 24.12.11 (fmesher 0.7.0, Matrix 1.7-5)
  - sf 1.1-0, spdep 1.4-2, spData 2.3.4, tigris 2.2.1
  - haven 2.5.5, readxl 1.4.5
  - dplyr 1.2.0, tidyr 1.3.2, readr 2.2.0, purrr 1.2.1, stringr 1.6.0, tibble 3.3.1, ggplot2 4.0.2
- **Threading:** INLA used its default threading (`num.threads = "8:1"`). All other steps were single-threaded.

## Data availability

The linked HIV surveillance and all-payer hospital discharge records belong to the South Carolina Department of Public Health and the South Carolina Revenue and Fiscal Affairs Office. They are available only under data use agreements with those agencies.

## License

MIT. See `LICENSE`.
