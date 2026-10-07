## =====================================================================
## run_all.R -- full analysis on the RESTRICTED data, then the sensitivity
##              variants. The restricted files are not distributed; to see
##              the pipeline work end to end without them, use run_toy.R.
## Usage (from the repository root):
##   SAE_DATA_MODE=restricted SAE_FILE_CASES=/path/cases.sas7bdat \
##   SAE_FILE_UB=/path/ub.sas7bdat Rscript run_all.R
## =====================================================================
## FIRST PASS: leave cfg$quick <- TRUE in 00_config.R and RUN_SENSITIVITY
## below at FALSE. That gives a complete preliminary result in a few minutes.
## Then set cfg$quick <- FALSE and RUN_SENSITIVITY <- TRUE for the real run.
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
RUN_SENSITIVITY <- TRUE

Sys.setenv(SAE_DATA_MODE = Sys.getenv("SAE_DATA_MODE", unset = "restricted"))
source(file.path("R", "00_config.R"))
assign("cfg", cfg, envir = globalenv())

source(file.path("R", "01_build_panel.R"),        local = new.env())
source(file.path("R", "02_descriptives_moran.R"), local = new.env())
source(file.path("R", "03_models.R"),             local = new.env())
source(file.path("R", "04_tables_figures.R"),     local = new.env())

## ---- sensitivity analyses ---------------------------------------------
## Each variant re-runs 01-03 with one option changed and files its own
## panel and comparison table under output/sens_<label>/.
## The first eight variants were prespecified; the last four were added after
## the primary analysis and are reported as exploratory.
sens_grid <- tibble::tribble(
  ~label,            ~opts,
  ## phenotype / code-list variants
  "pheno_with_alcohol", list(exclude_category_regex = "tobacc|nicotin|smok|cigar|vap",
                             alcohol_blocklist = character(0)),
  "pheno_revised",    list(phenotype = "revised_core"),
  "match_prefix",     list(match_mode = "prefix"),
  ## a quick "what if we add drug dependence codes" probe that does not
  ## require the workbook to be revised first
  "add_dependence",   list(codes_force_in = c("F112", "F122", "F132", "F142",
                                              "F152", "F162", "F182", "F192",
                                              "3040", "3041", "3042", "3043",
                                              "3044", "3045", "3046", "3047",
                                              "3048", "3049", "292")),
  ## design variants
  "denom_encounter",  list(denominator = "encounter"),
  "county_current",   list(county_var = "CURRENT_COUNTY"),
  "adj_rook",         list(adjacency = "rook"),
  "era_off",          list(enforce_era = FALSE),
  ## added after the primary analysis (exploratory)
  ## cohort timing: only encounters on or after the HIV diagnosis date count
  "post_dx_only",     list(require_post_dx = TRUE),
  ## penalized-complexity prior sensitivity (primary: Pr(sd > 1) = 0.01,
  ## Pr(phi < 0.5) = 2/3)
  "prior_sd_u05",     list(pc_u = 0.5),
  "prior_sd_u2",      list(pc_u = 2),
  "prior_phi_half",   list(phi_alpha = 0.5)
)

if (RUN_SENSITIVITY) {
  base_cfg <- cfg
  base_out <- cfg$dir_out
  for (i in seq_len(nrow(sens_grid))) {
    lab  <- sens_grid$label[i]
    opts <- sens_grid$opts[[i]]
    message("\n===== SENSITIVITY: ", lab, " =====")
    cfg <- modifyList(base_cfg, opts)
    cfg$dir_out <- file.path(base_out, paste0("sens_", lab))
    cfg$dir_log <- file.path(cfg$dir_out, "logs")
    cfg$run_designC <- FALSE
    cfg$run_decomposition <- FALSE     # decomposition models: primary analysis only
    assign("cfg", cfg, envir = globalenv())
    ok <- try({
      source(file.path("R", "01_build_panel.R"),        local = new.env())
      source(file.path("R", "02_descriptives_moran.R"), local = new.env())
      source(file.path("R", "03_models.R"),             local = new.env())
    }, silent = FALSE)
    if (inherits(ok, "try-error")) message("SENSITIVITY ", lab, " FAILED -- see its log")
  }
  assign("cfg", base_cfg, envir = globalenv())

  ## does the model ranking survive every variant?
  files <- list.files(base_out, pattern = "^model_comparison\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  rank_tab <- purrr::map_dfr(files, function(f) {
    variant <- basename(dirname(f))
    readr::read_csv(f, show_col_types = FALSE) |>
      mutate(variant = ifelse(variant == basename(base_out), "primary", variant)) |>
      arrange(rmse_pp_A) |> mutate(rank_designA = row_number()) |>
      select(variant, model, rank_designA, rmse_pp_A, waic)
  })
  readr::write_csv(rank_tab, file.path(base_out, "tableS3_ranking_across_variants.csv"))

  ## every reported metric, per variant and model (long format):
  ## Design A and B RMSE, WAIC, DIC, mean log CPO, Design B interval coverage
  metric_cols <- c("rmse_pp_A", "rmse_pp_B", "waic", "dic", "mean_log_cpo",
                   "cover95_B", "piwidth_pp_B")
  metrics_tab <- purrr::map_dfr(files, function(f) {
    variant <- basename(dirname(f))
    d <- readr::read_csv(f, show_col_types = FALSE)
    d |> mutate(variant = ifelse(variant == basename(base_out), "primary", variant)) |>
      select(variant, model, any_of(metric_cols))
  })
  readr::write_csv(metrics_tab, file.path(base_out, "tableS3_metrics_across_variants.csv"))

  ## how much does each variant move the headline prevalence?
  pfiles <- list.files(base_out, pattern = "^panel_.*_(surveillance|encounter)\\.csv$",
                       recursive = TRUE, full.names = TRUE)
  pfiles <- pfiles[!grepl("PREVIOUS", pfiles)]
  prev_tab <- purrr::map_dfr(pfiles, function(f) {
    variant <- basename(dirname(f))
    readr::read_csv(f, show_col_types = FALSE) |>
      group_by(year) |>
      summarise(wrate_pct = 100 * sum(nume) / sum(deno), .groups = "drop") |>
      mutate(variant = ifelse(variant == basename(base_out), "primary", variant))
  })
  readr::write_csv(prev_tab, file.path(base_out, "tableS4_encounter_rate_across_variants.csv"))
  message("\nRanking and prevalence stability tables written.")
}

runtime_min <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)
message("\nTotal runtime: ", runtime_min, " minutes")
source(file.path("R", "session_report.R"))
write_session_report(file.path(cfg$dir_log, "session_report_run_all.txt"),
                     extra = paste("Total wall-clock runtime (minutes):", runtime_min))
