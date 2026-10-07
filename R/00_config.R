## =====================================================================
## 00_config.R  -- paths, analysis options, county crosswalk, helpers
## Project: Validation design and model choice in small-area surveillance:
##          SUD-related hospital encounters among PLWH, South Carolina 2005-2023
## =====================================================================
## EDIT ONLY THIS FILE for a normal run. Everything downstream reads `cfg`.

## ---- 0. Packages -------------------------------------------------------
.needed <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "haven",
             "readxl", "sf", "spdep", "ggplot2", "tidyselect")
.missing <- .needed[!vapply(.needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("Missing R packages: ", paste(.missing, collapse = ", "), "\nInstall with:\n",
       "  install.packages(c(", paste0('"', .missing, '"', collapse = ", "), "))", call. = FALSE)
if (!requireNamespace("INLA", quietly = TRUE))
  message("Note: INLA is not installed; it is needed for 03_models.R. Install with:\n",
          '  install.packages("INLA", repos = c(getOption("repos"), ',
          'INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)')
suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(purrr); library(tibble); library(haven); library(readxl)
}))

cfg <- list()

## ---- 1. Data mode and paths --------------------------------------------
## Two ways to run the same pipeline:
##   "toy"        -- fully synthetic files produced by R/simulate_toy_data.R.
##                   No restricted data are needed; this is what run_toy.R uses.
##   "restricted" -- the linked South Carolina HIV surveillance (case) file and
##                   the all-payer hospital discharge (UB) file. These are held
##                   under a data use agreement and are NOT distributed. Point
##                   the two paths below at your local copies.
## The mode can also be set from the shell: SAE_DATA_MODE=restricted Rscript run_all.R
cfg$data_mode <- Sys.getenv("SAE_DATA_MODE", unset = "toy")
cfg$dir_code  <- "R"

if (cfg$data_mode == "restricted") {
  cfg$file_cases <- Sys.getenv("SAE_FILE_CASES", unset = "path/to/li_dhec_cases0725.sas7bdat")
  cfg$file_ub    <- Sys.getenv("SAE_FILE_UB",    unset = "path/to/li_ub_allpayer0825.sas7bdat")
} else {
  cfg$dir_toy    <- "toy_data"
  cfg$file_cases <- file.path(cfg$dir_toy, "toy_cases.csv.gz")
  cfg$file_ub    <- file.path(cfg$dir_toy, "toy_ub_discharges.csv.gz")
}

## SC county boundaries: a local shapefile, or tigris (needs internet).
## Offline? Download the TIGER/Line 2020 county shapefile once and point
## SAE_SHAPEFILE (or this line) at it.
cfg$shapefile   <- if (nzchar(Sys.getenv("SAE_SHAPEFILE"))) Sys.getenv("SAE_SHAPEFILE") else NULL
cfg$use_tigris  <- TRUE
cfg$tiger_year  <- 2020

cfg$dir_out <- if (cfg$data_mode == "restricted") "output" else "output_toy"
cfg$dir_log <- file.path(cfg$dir_out, "logs")

## ---- 2. THE CODE LIST -------------------------------------------------
## This block is the only thing you touch when the code list is revised.
## Point at the new workbook and rerun: the loader auto-detects the column
## layout, expands wildcards and ranges, hashes the file, and writes a diff
## against the previously used list (output/codelist_diff.csv) plus the
## resulting change in prevalence (output/panel_change_vs_previous.csv).
cfg$file_icd  <- "data/ICD_Code_List_SUD.xlsx"   # public; the codes are listed in Supplementary Table S1
cfg$icd_sheet <- "Substance use"   # NULL = first sheet

## Column detection. Leave NULL to auto-detect; override only if a future
## workbook uses names the detector cannot guess.
##   wide layout : an ICD-9 column and an ICD-10 column per row
##   long layout : a `code` column plus a `code_system` column
cfg$icd_cols <- list(icd9 = NULL, icd10 = NULL, code = NULL, system = NULL,
                     category = NULL, subcategory = NULL, label = NULL,
                     include = NULL)

## Tobacco (and, in the primary phenotype, alcohol) is excluded.
## Two independent guards, so a renamed or restructured code list cannot
## silently let nicotine codes back into the numerator:
##   (a) drop any category whose name matches this regex
cfg$exclude_category_regex <- "tobacc|nicotin|smok|cigar|vap|alcohol"
##   (b) drop these codes whatever any list says (normalized, prefix-matched)
cfg$tobacco_blocklist <- c("F17", "Z720", "Z716", "3051", "V1582")
cfg$alcohol_blocklist <- c("F10", "291", "303", "3050",          # use disorders
                           "K70", "K292", "K860", "G312", "G621", # organ sequelae
                           "G721", "I426", "5710", "5711", "5712",
                           "5713", "3575", "4255", "53530", "53531")

## Categories treated as alcohol-attributable ORGAN disease (sequelae rather
## than a coded use disorder). Excluded from `core`, kept in `broad`.
cfg$organ_category_regex <- "alcoholic_disease|organ|liver|hepat|cirrho"

## Phenotype. One of:
##   "core"          PRIMARY. listed codes, minus tobacco, minus organ disease
##   "listed_all"    every listed code minus tobacco (organ disease kept in)
##   "revised_core"  core + regex families the current list omits
##                   (ICD-9 304.xx / 292.x; ICD-10 F1x.2x dependence, F1x.9x)
##   "broad"         revised_core + alcohol-attributable organ disease
cfg$phenotype <- "core"

## Force codes in or out for one run without editing the workbook. Normalized,
## dot-free, prefix-matched. Handy for a quick "what if we add opioid
## dependence" check before the list is formally revised.
cfg$codes_force_in  <- character(0)   # e.g. c("F112", "3040")
cfg$codes_force_out <- character(0)

## ---- 3. Analysis options ---------------------------------------------
cfg$years <- 2005:2023

## Denominator: "surveillance" (PRIMARY) or "encounter"
cfg$denominator <- "surveillance"
## County assignment: "COUNTY_AT_HIV_DX" (primary) or "CURRENT_COUNTY"
cfg$county_var  <- "COUNTY_AT_HIV_DX"
## ICD matching: "exact" (a listed code matches only itself) or "prefix"
## (a diagnosis matches if it starts with a listed code). Wildcard entries in
## the workbook are always prefixes regardless of this setting.
cfg$match_mode  <- "exact"
## Require the code system to agree with the discharge date era (ICD-10 from
## Oct 2015). TRUE = primary; FALSE = sensitivity.
cfg$enforce_era <- TRUE
## Count encounters before the HIV diagnosis date but inside the diagnosis
## calendar year? FALSE = primary analysis.
cfg$require_post_dx <- FALSE

cfg$adjacency  <- "queen"      # or "rook"
cfg$moran_nsim <- 9999

## PC priors
cfg$pc_u <- 1; cfg$pc_alpha <- 0.01
cfg$phi_u <- 0.5; cfg$phi_alpha <- 2/3

## Holdout designs
cfg$cv_folds      <- 10
cfg$cv_seed       <- 20260820
cfg$forward_years <- 2019:2023
cfg$run_designB   <- TRUE
cfg$run_designC   <- TRUE      # leave-one-county-out: 46 refits per model, slow

## QUICK MODE -- first pass / preliminary result. Fewer folds, no Design C,
## fewer permutations; the whole pipeline finishes in a few minutes.
cfg$quick <- FALSE
if (isTRUE(cfg$quick)) {
  cfg$cv_folds <- 5; cfg$run_designC <- FALSE; cfg$moran_nsim <- 999
}

cfg$rmse_margin_pp <- 0.25     # practical-equivalence margin, selection rule

## ---- Computing ----------------------------------------------------------
## INLA threads. NULL keeps INLA's default, which is what the published run
## used (num.threads = "8:1" on an 8-core laptop). Set e.g. "4:1" to fix it.
## Fitting times depend on CPU, memory and this setting; R/session_report.R
## records all three next to the results.
cfg$inla_threads <- NULL

## ---- Selection-rule adequacy thresholds ---------------------------------
cfg$calib_gap_max_pp    <- 1.5  # max |observed - predicted| over predicted deciles, Design A
cfg$resid_sig_years_max <- 2    # residual Moran's I significant (FDR 5%) in at most 2 of 19 years

## ---- Additional analyses ----------------------------------------------
## Decomposition models separate the temporal from the spatial contribution
## and compare RW1 with RW2. They are reported in the supplement and are NOT
## part of the prespecified selection set (M0-M4).
cfg$run_decomposition <- TRUE
## Posterior predictive draws per held-out cell, for 95% prediction-interval
## coverage and width.
cfg$n_pred_draws <- 1000
## Predictive intervals for the beta-binomial model. FALSE (as published):
## Y ~ Binomial(N, p) for every model, so the beta-binomial overdispersion is
## not propagated. TRUE: Y ~ BetaBinomial(N, p, rho) with rho drawn from its
## posterior, for the beta-binomial model only. Point predictions, RMSE, MAE
## and the selection rule are unaffected by this switch.
cfg$pi_betabinomial <- identical(Sys.getenv("SAE_PI_BETABINOMIAL"), "1")
## Frequentist cross-check of M1/M2 with glmmTMB. Not reported in the paper;
## every model-based approach (M1-M4) is fit in INLA.
cfg$run_glmmtmb <- FALSE
## Model used for the county maps (Figure 2) and the shrinkage figure:
## "M3_bym2_rw2" (default) or "selected" (prespecified Design A choice)
cfg$map_model <- "M3_bym2_rw2"

## ---- 4. South Carolina county crosswalk (46 counties) ------------------
sc_counties <- tibble::tribble(
  ~county,        ~fips,
  "Abbeville",    "45001",  "Aiken",        "45003",  "Allendale",    "45005",
  "Anderson",     "45007",  "Bamberg",      "45009",  "Barnwell",     "45011",
  "Beaufort",     "45013",  "Berkeley",     "45015",  "Calhoun",      "45017",
  "Charleston",   "45019",  "Cherokee",     "45021",  "Chester",      "45023",
  "Chesterfield", "45025",  "Clarendon",    "45027",  "Colleton",     "45029",
  "Darlington",   "45031",  "Dillon",       "45033",  "Dorchester",   "45035",
  "Edgefield",    "45037",  "Fairfield",    "45039",  "Florence",     "45041",
  "Georgetown",   "45043",  "Greenville",   "45045",  "Greenwood",    "45047",
  "Hampton",      "45049",  "Horry",        "45051",  "Jasper",       "45053",
  "Kershaw",      "45055",  "Lancaster",    "45057",  "Laurens",      "45059",
  "Lee",          "45061",  "Lexington",    "45063",  "McCormick",    "45065",
  "Marion",       "45067",  "Marlboro",     "45069",  "Newberry",     "45071",
  "Oconee",       "45073",  "Orangeburg",   "45075",  "Pickens",      "45077",
  "Richland",     "45079",  "Saluda",       "45081",  "Spartanburg",  "45083",
  "Sumter",       "45085",  "Union",        "45087",  "Williamsburg", "45089",
  "York",         "45091"
) |> dplyr::mutate(county_key = toupper(county))
stopifnot(nrow(sc_counties) == 46L)

## ---- 5. Helpers -------------------------------------------------------
read_panel <- function(path) {
  readr::read_csv(path, col_types = readr::cols(fips = readr::col_character())) |>
    dplyr::mutate(fips = sprintf("%05s", as.character(fips)))
}
dir_init <- function() {
  dir.create(cfg$dir_out, showWarnings = FALSE, recursive = TRUE)
  dir.create(cfg$dir_log, showWarnings = FALSE, recursive = TRUE)
}
.log_file <- function() file.path(cfg$dir_log, "run_log.txt")
log_step <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n"); dir_init(); cat(msg, "\n", file = .log_file(), append = TRUE)
  invisible(msg)
}
check_that <- function(condition, msg) {
  if (!isTRUE(condition)) stop("CHECK FAILED: ", msg, call. = FALSE)
  log_step("check ok: ", msg); invisible(TRUE)
}

## Normalize an ICD code for matching: uppercase, drop dots and spaces.
## Wildcard characters are preserved so the loader can turn them into prefixes.
fix_excel_float <- function(x) {
  x <- as.character(x)
  bad <- !is.na(x) & grepl("\\.[0-9]{6,}$", x)
  if (any(bad)) x[bad] <- format(round(as.numeric(x[bad]), 3), trim = TRUE, scientific = FALSE)
  x
}

normalize_icd <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9%*-]", "", x)
  ifelse(x %in% c("", "NOCODE", "NA", "."), NA_character_, x)
}

clean_county <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+CO\\.?$", "", x)
  x <- gsub("\\s+COUNTY$", "", x)
  x <- gsub("[^A-Z ]", "", x)
  x <- gsub("\\s+", " ", trimws(x))
  ifelse(x == "", NA_character_, x)
}

parse_yyyymm <- function(x, default_month = 6L) {
  x  <- trimws(as.character(x))
  yr <- suppressWarnings(as.integer(substr(x, 1, 4)))
  mo <- suppressWarnings(as.integer(substr(x, 5, 6)))
  mo[is.na(mo) | mo < 1 | mo > 12] <- default_month
  mo[is.na(yr)] <- NA_integer_
  tibble(year = yr, month = mo)
}

lower_names <- function(d) { names(d) <- tolower(names(d)); d }

## Read a source file. SAS (.sas7bdat) for the restricted data; CSV (optionally
## gzipped) for the toy data. CSV columns are read as text so that ICD codes
## keep their leading zeros, then the few numeric fields are converted.
read_source <- function(path, n_max = Inf, col_select = NULL) {
  if (!file.exists(path)) stop("Data file not found: ", path,
    if (cfg$data_mode == "toy") "\nRun run_toy.R, which simulates the toy files first." else "")
  if (grepl("\\.sas7bdat$", path, ignore.case = TRUE)) {
    if (is.null(col_select)) return(haven::read_sas(path, n_max = n_max))
    return(haven::read_sas(path, n_max = n_max, col_select = tidyselect::all_of(col_select)))
  }
  d <- if (is.null(col_select))
    readr::read_csv(path, n_max = n_max, col_types = readr::cols(.default = readr::col_character()),
                    na = c("", "NA"), progress = FALSE)
  else
    readr::read_csv(path, n_max = n_max, col_types = readr::cols(.default = readr::col_character()),
                    col_select = tidyselect::all_of(col_select), na = c("", "NA"), progress = FALSE)
  num <- names(d)[toupper(names(d)) %in% c("AGE_AT_HIV_DX", "DISYEAR", "DISMTH", "TIME_DXDATE_DISD")]
  for (v in num) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  d
}

## South Carolina county boundaries (46 counties), as an sf object with
## columns fips and county, ordered by FIPS. Used by 02_ and by run_toy.R.
load_sc_boundaries <- function() {
  if (!is.null(cfg$shapefile)) {
    sc <- sf::st_read(cfg$shapefile, quiet = TRUE) |> dplyr::rename_with(tolower)
    fips_col <- intersect(c("geoid", "fips", "geoid20"), names(sc))[1]
    sc <- sc |> dplyr::mutate(fips = as.character(.data[[fips_col]])) |>
      dplyr::filter(substr(fips, 1, 2) == "45")
  } else if (isTRUE(cfg$use_tigris)) {
    if (!requireNamespace("tigris", quietly = TRUE))
      stop("tigris is not installed and cfg$shapefile is NULL. Either ",
           'install.packages("tigris") on a machine with internet, or download ',
           "the TIGER/Line county shapefile and set cfg$shapefile.")
    sc <- tigris::counties(state = "SC", year = cfg$tiger_year, cb = FALSE,
                           progress_bar = FALSE) |>
      dplyr::rename_with(tolower) |> dplyr::mutate(fips = as.character(geoid))
  } else stop("Set cfg$shapefile or cfg$use_tigris.")
  sc |> dplyr::select(fips) |>                       # geometry column is sticky
    dplyr::left_join(sc_counties |> dplyr::select(county, fips), by = "fips") |>
    dplyr::arrange(fips)
}

pc_prec <- function() list(prec = list(prior = "pc.prec",
                                       param = c(cfg$pc_u, cfg$pc_alpha)))
pc_bym2 <- function() list(prec = list(prior = "pc.prec",
                                       param = c(cfg$pc_u, cfg$pc_alpha)),
                           phi  = list(prior = "pc",
                                       param = c(cfg$phi_u, cfg$phi_alpha)))

dir_init()
log_step("config loaded | data=", cfg$data_mode, " | phenotype=", cfg$phenotype,
         " | denominator=", cfg$denominator,
         " | county_var=", cfg$county_var,
         " | match=", cfg$match_mode,
         " | quick=", cfg$quick)
