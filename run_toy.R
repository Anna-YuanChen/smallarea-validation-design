## =====================================================================
## run_toy.R -- the whole pipeline on synthetic data, no restricted files
##
##   Rscript run_toy.R            (from the repository root)
##
## 1. simulates a person-level HIV case file and a hospital discharge file
##    with the same layout as the restricted sources (R/simulate_toy_data.R)
## 2. runs 01_ (code list, cohort, ICD matching, county-year panel),
##    02_ (adjacency, Table 1, Moran's I), 03_ (models M1-M4 and crude,
##    holdout Designs A/B[/C], selection rule) and 04_ (tables, figures)
## 3. compares every model's full-data estimates with the known truth
##
## Quick mode (default) takes a few minutes on a laptop: 5 folds, no Design C,
## 999 permutations, no decomposition models. Set TOY_FULL=1 for the full
## protocol used in the paper (10 folds, Designs A-C, 9,999 permutations).
## =====================================================================
## ---- locate the repository folder (works from Rscript, RStudio, or source())
.find_root <- function() {
  if (file.exists(file.path("R", "00_config.R"))) return(invisible(TRUE))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(f)) f <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) character(0))
  if (!length(f) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    f <- rstudioapi::getSourceEditorContext()$path
  if (length(f) && nzchar(f[1])) setwd(dirname(normalizePath(f[1])))
  if (!file.exists(file.path("R", "00_config.R")))
    stop("Please run this script from the repository folder, e.g. open ",
         "smallarea-validation-design.Rproj in RStudio, or setwd() to the folder first.",
         call. = FALSE)
}
.find_root()
t_start <- Sys.time()
Sys.setenv(SAE_DATA_MODE = "toy")
source(file.path("R", "00_config.R"))

full <- identical(Sys.getenv("TOY_FULL"), "1")
cfg$quick <- !full
if (cfg$quick) {
  cfg$cv_folds <- 5; cfg$run_designC <- FALSE; cfg$moran_nsim <- 999
  cfg$run_decomposition <- FALSE
}
assign("cfg", cfg, envir = globalenv())

## ---- 1. synthetic source files -----------------------------------------
source(file.path("R", "codelist.R"))
source(file.path("R", "simulate_toy_data.R"))
cl <- load_code_list()                          # same code list as the paper
sc <- load_sc_boundaries()
stopifnot(identical(sc$fips, sc_counties$fips))
nb <- spdep::poly2nb(sc, queen = TRUE)
sim <- simulate_toy_data(out_dir    = cfg$dir_toy,
                         n_persons  = as.integer(Sys.getenv("TOY_N", "30000")),
                         sud_codes9  = cl$codes9_exact,
                         sud_codes10 = cl$codes10_exact,
                         nb          = nb)
log_step("toy data written to ", cfg$dir_toy, ": ", nrow(sim$cases), " persons, ",
         nrow(sim$ub), " discharge rows")

## ---- 2. the unchanged pipeline -------------------------------------------
source(file.path("R", "01_build_panel.R"),        local = new.env())
source(file.path("R", "02_descriptives_moran.R"), local = new.env())
source(file.path("R", "03_models.R"),             local = new.env())
source(file.path("R", "04_tables_figures.R"),     local = new.env())

## ---- 3. how close is each approach to the truth? -------------------------
res   <- readRDS(file.path(cfg$dir_out, sprintf("models_%s_%s.rds", cfg$phenotype, cfg$denominator)))
panel <- readr::read_csv(file.path(cfg$dir_out, sprintf("panel_%s_%s_indexed.csv",
                         cfg$phenotype, cfg$denominator)), show_col_types = FALSE,
                         col_types = readr::cols(fips = readr::col_character()))
truth <- readr::read_csv(file.path(cfg$dir_toy, "truth_county_year.csv"),
                         col_types = readr::cols(fips = readr::col_character()))
panel <- dplyr::left_join(panel, truth, by = c("fips", "county", "year"))
strat <- cut(panel$deno, c(-Inf, 50, 100, Inf), right = FALSE,
             labels = c("N<50", "50<=N<100", "N>=100"))
err <- function(est) 100 * (est - panel$p_true)
rows <- c(list(M0_crude = err(panel$encounter_rate)),
          lapply(res$full_fits, function(f) err(f$summary.fitted.values$mean)))
vs_truth <- purrr::imap_dfr(rows, function(e, nm) {
  tibble::tibble(model = nm, stratum = c("All", levels(strat)),
                 rmse_vs_truth_pp = c(sqrt(mean(e^2)),
                   vapply(levels(strat), function(s) sqrt(mean(e[strat == s]^2)), numeric(1))))
}) |> tidyr::pivot_wider(names_from = stratum, values_from = rmse_vs_truth_pp)
readr::write_csv(vs_truth, file.path(cfg$dir_out, "toy_rmse_vs_truth.csv"))
message("\nRMSE of full-data estimates against the known truth (percentage points):")
print(as.data.frame(vs_truth), digits = 3)
message("\nPrespecified (Design A) selection: ", res$selected)

runtime_min <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)
source(file.path("R", "session_report.R"))
write_session_report(file.path(cfg$dir_log, "session_report_toy.txt"),
                     extra = paste("Total wall-clock runtime (minutes):", runtime_min))
message("Total runtime: ", runtime_min, " minutes. Results in ", cfg$dir_out, "/")
