## =====================================================================
## 05_revision_analyses.R -- analyses added after the primary comparison
##
## Run AFTER run_all.R (or run_toy.R) has finished. It reads the saved
## outputs of 01-03 and refits models only where new holdout settings are
## needed. Every file it writes holds aggregate results only (no person-level
## or county-year records), so the outputs can be shared with co-authors.
##
##   Checkpoint 2 (cohort and code-list checks; reads the source files)
##     person-years, diagnosis-year range, current vs diagnosis county,
##     on-list diagnosis entries, abuse codes absent from the list
##   (i)    persistence comparators under Design B (last value; last 3 years)
##   (ii)   Design B RMSE and mean error by target year
##   (iii)  paired county-cluster bootstrap of RMSE differences
##   (iv)   margin sweep for the selection rule, 0.05-0.30 pp
##   (v)    earlier forecast origins, targets 2010-2014
##   (vi)   two-year forecast horizon, targets 2019-2023
##   (vii)  log predictive density and noise-corrected MSE
##   (viii) level shift vs Design B gain across sensitivity variants
##
##   out: <dir_out>/revision/*.csv, figS6_level_shift.png,
##        revision_summary.txt   <- every number quoted in the paper
## Usage: open the project in RStudio and run  source("R/05_revision_analyses.R")
##        (or, in a flat folder,  source("05_revision_analyses.R"))
## =====================================================================

## ---- locate the code folder ----------------------------------------------
code_dir <- if (file.exists(file.path("R", "00_config.R"))) "R" else
  if (file.exists("00_config.R")) "." else
  stop("Run this from the project folder (the one containing 00_config.R or R/00_config.R).")
## Run on its own, the script always reads 00_config.R afresh (so a `cfg` left in
## the workspace by another run cannot be picked up by mistake); run_all.R and
## run_toy.R set .rev_use_existing_cfg so that their settings are used instead.
if (!isTRUE(get0(".rev_use_existing_cfg", envir = globalenv())) || !exists("cfg", envir = globalenv())) {
  if (code_dir == "R") Sys.setenv(SAE_DATA_MODE = Sys.getenv("SAE_DATA_MODE", unset = "restricted"))
  source(file.path(code_dir, "00_config.R"))
  assign("cfg", cfg, envir = globalenv())
}
cfg <- get("cfg", envir = globalenv())
suppressWarnings(suppressMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2) }))

## ---- options ---------------------------------------------------------------
REV <- list(
  seed        = 20261014,
  n_boot      = 2000,
  margins     = seq(0.05, 0.30, by = 0.05),
  early_years = 2010:2014,          # (v) forecast targets before the coding transition
  run_raw     = TRUE,               # Checkpoint 2 counts (reads the source files)
  run_refits  = TRUE,               # (v), (vi), (vii): INLA refits
  lpd_designC = TRUE,               # (vii) for Design C: 46 x 4 refits, the slowest step
  suppress_lt = 11                  # counts 1-10 are shown as "<11" in the summary
)
if (isTRUE(cfg$quick)) REV$lpd_designC <- FALSE

dir_rev <- file.path(cfg$dir_out, "revision")
dir.create(dir_rev, showWarnings = FALSE, recursive = TRUE)
summary_file <- file.path(dir_rev, "revision_summary.txt")
writeLines(c("Revision analyses -- aggregate results",
             paste("run:", format(Sys.time(), "%Y-%m-%d %H:%M"), "| data mode:",
                   if (is.null(cfg$data_mode)) "restricted" else cfg$data_mode), ""),
           summary_file)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n"); cat(txt, "\n", file = summary_file, append = TRUE)
}
supp <- function(n) ifelse(n > 0 & n < REV$suppress_lt, paste0("<", REV$suppress_lt), format(n, big.mark = ","))
f3 <- function(x) formatC(x, format = "f", digits = 3)
t0 <- Sys.time()

## ---- saved outputs of 01-03 -----------------------------------------------
pfile <- file.path(cfg$dir_out, sprintf("panel_%s_%s_indexed.csv", cfg$phenotype, cfg$denominator))
panel <- readr::read_csv(pfile, col_types = readr::cols(fips = readr::col_character())) |>
  mutate(county_id2 = county_id, year_id2 = year_id,
         st_id = as.integer(interaction(county_id, year_id, drop = TRUE)))
cv <- readr::read_csv(file.path(cfg$dir_out, "cv_predictions.csv"), show_col_types = FALSE)
strat <- function(n) as.character(cut(n, c(-Inf, 50, 100, Inf), right = FALSE,
                                      labels = c("N<50", "50<=N<100", "N>=100")))
APPROACH_ORDER <- c("M0_crude", "P1_last_value", "P3_last_3y", "M1_binomial_glmm",
                    "M2_betabinomial_glmm", "M3_bym2_rw2", "M4_bym2_rw2_typeI",
                    "D1_iid_rw2", "D2_bym2_iidyear", "D3_bym2_rw1")
SEL_SET <- c("M1_binomial_glmm", "M2_betabinomial_glmm", "M3_bym2_rw2", "M4_bym2_rw2_typeI")

## =====================================================================
## Checkpoint 2 -- cohort, denominator and code-list checks
## =====================================================================
say("== Checkpoint 2: cohort, denominator and code list ==")
say("Total person-years (sum of county-year denominators): ", format(sum(panel$deno), big.mark = ","))
pyf <- file.path(cfg$dir_out, sprintf("person_year_%s_%s.rds", cfg$phenotype, cfg$denominator))
py_ids <- NULL
if (file.exists(pyf)) {
  py <- readRDS(pyf)
  py_ids <- unique(py$rfa_id)
  dx <- py |> distinct(rfa_id, dx_year)
  say("Persons contributing person-years: ", format(nrow(dx), big.mark = ","),
      " | diagnosis years ", min(dx$dx_year), "-", max(dx$dx_year),
      " | diagnosed before 2005: ", format(sum(dx$dx_year < 2005), big.mark = ","),
      " (", sprintf("%.1f", 100 * mean(dx$dx_year < 2005)), "%)")
  say("Person-years from persons diagnosed before 2005: ",
      format(sum(py$dx_year < 2005), big.mark = ","))
  rm(py); invisible(gc())
}

if (isTRUE(REV$run_raw)) {
  rd <- if (exists("read_source")) read_source else function(path, ...) {
    if (grepl("\\.sas7bdat$", path, ignore.case = TRUE)) haven::read_sas(path, ...)
    else readr::read_csv(path, show_col_types = FALSE, ...)
  }
  ok <- try({
    ## ---- county at diagnosis vs current county, analytic cohort ----
    cases <- rd(cfg$file_cases) |> haven::zap_labels() |> lower_names()
    cases$rfa_id <- as.character(cases$rfa_id)
    if (!is.null(py_ids) && all(c("county_at_hiv_dx", "current_county") %in% names(cases))) {
      cc <- cases |> filter(rfa_id %in% py_ids) |> distinct(rfa_id, .keep_all = TRUE) |>
        mutate(k_dx = clean_county(county_at_hiv_dx), k_cur = clean_county(current_county),
               k_dx = ifelse(k_dx == "MC CORMICK", "MCCORMICK", k_dx),
               k_cur = ifelse(k_cur == "MC CORMICK", "MCCORMICK", k_cur),
               cur_ok = k_cur %in% sc_counties$county_key)
      n_diff <- sum(cc$cur_ok & cc$k_cur != cc$k_dx)
      say("Current county differs from county at diagnosis: ", format(n_diff, big.mark = ","),
          " of ", format(nrow(cc), big.mark = ","), " (", sprintf("%.1f", 100 * n_diff / nrow(cc)),
          "%); current county missing or outside SC: ", format(sum(!cc$cur_ok), big.mark = ","))
    }
    rm(cases); invisible(gc())

    ## ---- diagnosis entries: on the list, and abuse codes absent from it ----
    source(file.path(code_dir, "codelist.R"))
    cl <- suppressMessages(load_code_list())
    icd9_fields  <- c("pdiag", "adm_diag", paste0("sdiag", 1:14))
    icd10_fields <- c("adm_diag10", "pdiag10", paste0("sdiag10_", 1:14))
    head1 <- rd(cfg$file_ub, n_max = 1)
    actual <- names(head1); names(actual) <- tolower(actual)
    have9 <- intersect(icd9_fields, names(actual)); have10 <- intersect(icd10_fields, names(actual))
    meta <- intersect(c("rfa_id", "disyear", "dismth", "time_dxdate_disd"), names(actual))
    ub <- rd(cfg$file_ub, col_select = unname(actual[c(meta, have9, have10)])) |>
      haven::zap_labels() |> lower_names()
    ub$rfa_id <- as.character(ub$rfa_id); ub$disyear <- suppressWarnings(as.integer(ub$disyear))
    ub$dismth <- if ("dismth" %in% names(ub)) suppressWarnings(as.integer(ub$dismth)) else NA_integer_
    ub <- filter(ub, !is.na(disyear), disyear %in% cfg$years)
    if (!is.null(py_ids)) ub <- filter(ub, rfa_id %in% py_ids)
    era10 <- ub$disyear > 2015 | (ub$disyear == 2015 & !is.na(ub$dismth) & ub$dismth >= 10)
    ## abuse subcategories (F1x.1, x = 1-6, 8, 9; ICD-9-CM 305.2-305.9), alcohol and tobacco excluded
    RX10_ABUSE <- "^F1[1-689]1"
    RX9_ABUSE  <- "^305[2-9]"
    n_on <- 0L; cand <- list(); on_py <- list(); cand_py <- list()
    for (v in c(have9, have10)) {
      slot10 <- v %in% have10
      cn <- normalize_icd(ub[[v]]); okc <- !is.na(cn)
      if (!any(okc)) next
      h10 <- match_sud(cn, cl, "icd10"); h9 <- match_sud(cn, cl, "icd9")
      is_sud <- if (isTRUE(cfg$enforce_era)) (h10 & (slot10 | era10)) | (h9 & !slot10 & !era10) else h9 | h10
      is_sud[is.na(is_sud)] <- FALSE
      n_on <- n_on + sum(is_sud)
      if (any(is_sud)) on_py[[v]] <- distinct(tibble(rfa_id = ub$rfa_id[is_sud], year = ub$disyear[is_sud]))
      ab <- okc & !is_sud & ((((slot10 | era10)) & grepl(RX10_ABUSE, cn)) |
                               ((!slot10 & !era10) & grepl(RX9_ABUSE, cn)))
      ab[is.na(ab)] <- FALSE
      if (any(ab)) {
        cand[[v]] <- tibble(code = cn[ab], year = ub$disyear[ab])
        cand_py[[v]] <- distinct(tibble(rfa_id = ub$rfa_id[ab], year = ub$disyear[ab]))
      }
    }
    say("Diagnosis entries matched by the code list (primary rules): ", format(n_on, big.mark = ","))
    on_py <- distinct(bind_rows(on_py))
    say("Person-years with >=1 on-list entry (check: should equal the panel numerator total ",
        format(sum(panel$nume), big.mark = ","), "): ", format(nrow(on_py), big.mark = ","))
    cand <- bind_rows(cand)
    if (nrow(cand)) {
      ctab <- cand |> group_by(code) |>
        summarise(n_entries = n(), first_year = min(year), last_year = max(year), .groups = "drop") |>
        arrange(desc(n_entries))
      readr::write_csv(mutate(ctab, n_entries = supp(n_entries)),
                       file.path(dir_rev, "cp2_abuse_codes_not_on_list.csv"))
      add_py <- anti_join(distinct(bind_rows(cand_py)), on_py, by = c("rfa_id", "year"))
      say("Abuse-subcategory codes present but not on the list: ", nrow(ctab), " codes, ",
          format(sum(ctab$n_entries), big.mark = ","), " entries; person-years that adding them ",
          "would newly count: ", supp(nrow(add_py)), " (",
          sprintf("%.2f", 100 * nrow(add_py) / max(sum(panel$nume), 1)), "% of the current numerator)")
      for (k in seq_len(min(nrow(ctab), 25)))
        say("   ", ctab$code[k], ": ", supp(ctab$n_entries[k]), " entries (", ctab$first_year[k], "-",
            ctab$last_year[k], ")")
    } else say("Abuse-subcategory codes present but not on the list: none")
    rm(ub); invisible(gc())
  }, silent = TRUE)
  if (inherits(ok, "try-error"))
    say("Checkpoint 2 source-file checks FAILED: ", substr(as.character(ok), 1, 300))
}

## =====================================================================
## (i)-(iv) Design B from the saved predictions + persistence comparators
## =====================================================================
say("\n== Design B: persistence comparators, by year, bootstrap, margin sweep ==")
cvB <- cv |> filter(design == "B") |>
  select(model, row, year, county, deno, observed, pred, err_pp)
cvB$county_id <- panel$county_id[cvB$row]

## persistence comparators, built from training years only
persist <- function(target_rows, lag_from) {
  purrr::map_dfr(target_rows, function(r) {
    cid <- panel$county_id[r]; t <- panel$year[r]
    last <- panel[panel$county_id == cid & panel$year == t - lag_from, ]
    l3   <- panel[panel$county_id == cid & panel$year %in% (t - lag_from - 2):(t - lag_from), ]
    tibble(row = r, P1_last_value = last$nume / last$deno, P3_last_3y = sum(l3$nume) / sum(l3$deno))
  })
}
hold_B <- sort(unique(cvB$row))
pB <- persist(hold_B, 1) |>
  pivot_longer(-row, names_to = "model", values_to = "pred") |>
  mutate(year = panel$year[row], county = panel$county[row], deno = panel$deno[row],
         observed = panel$encounter_rate[row], err_pp = 100 * (pred - observed),
         county_id = panel$county_id[row])
cvB <- bind_rows(cvB, pB) |> mutate(stratum = strat(deno))

met <- function(d) d |> summarise(n = n(), rmse_pp = sqrt(mean(err_pp^2)), mae_pp = mean(abs(err_pp)),
                                  bias_pp = mean(err_pp), .groups = "drop")
S8a <- bind_rows(cvB |> group_by(model) |> met() |> mutate(stratum = "All"),
                 cvB |> group_by(model, stratum) |> met()) |>
  mutate(model = factor(model, APPROACH_ORDER)) |> arrange(stratum != "All", stratum, model)
readr::write_csv(S8a, file.path(dir_rev, "S8a_designB_with_persistence.csv"))
rB <- S8a |> filter(stratum == "All") |> mutate(model = as.character(model))
gR <- function(m) rB$rmse_pp[rB$model == m]
say("Design B RMSE (MAE): ", paste0(rB$model, " ", f3(rB$rmse_pp), " (", f3(rB$mae_pp), ")", collapse = "; "))
for (s in c("N<50", "50<=N<100", "N>=100")) {
  x <- S8a |> filter(stratum == s)
  say("   ", s, " (n=", x$n[1], "): ", paste0(x$model, " ", f3(x$rmse_pp), collapse = "; "))
}

## (ii) by target year
S8b <- cvB |> group_by(model, year) |> met() |>
  mutate(model = factor(model, APPROACH_ORDER)) |> arrange(model, year)
readr::write_csv(S8b, file.path(dir_rev, "S8b_designB_by_year.csv"))
for (m in c("M0_crude", "P1_last_value", "P3_last_3y", "M1_binomial_glmm", "M3_bym2_rw2", "M4_bym2_rw2_typeI")) {
  x <- S8b |> filter(model == m)
  say("By year ", m, ": RMSE ", paste0(x$year, "=", f3(x$rmse_pp), collapse = ", "),
      " | mean error ", paste0(x$year, "=", sprintf("%+.3f", x$bias_pp), collapse = ", "))
}

## (iii) paired county-cluster bootstrap
set.seed(REV$seed)
counties <- sort(unique(cvB$county_id))
E <- cvB |> select(model, county_id, year, err_pp) |>
  pivot_wider(names_from = model, values_from = err_pp)
idx_by_c <- split(seq_len(nrow(E)), E$county_id)
mods <- intersect(APPROACH_ORDER, names(E))
boot <- replicate(REV$n_boot, {
  ii <- unlist(idx_by_c[as.character(sample(counties, replace = TRUE))], use.names = FALSE)
  vapply(mods, function(m) sqrt(mean(E[[m]][ii]^2)), numeric(1))
})
rm_full <- vapply(mods, function(m) sqrt(mean(E[[m]]^2)), numeric(1))
pair_ci <- function(a, b) {
  d <- boot[a, ] - boot[b, ]; rr <- 100 * (1 - boot[a, ] / boot[b, ])
  tibble(approach = a, reference = b, rmse_diff_pp = rm_full[a] - rm_full[b],
         lo95 = unname(quantile(d, .025)), hi95 = unname(quantile(d, .975)),
         reduction_pct = 100 * (1 - rm_full[a] / rm_full[b]),
         red_lo95 = unname(quantile(rr, .025)), red_hi95 = unname(quantile(rr, .975)))
}
pairs <- c(lapply(setdiff(mods, "M1_binomial_glmm"), function(a) c(a, "M1_binomial_glmm")),
           unlist(lapply(c("M3_bym2_rw2", "M4_bym2_rw2_typeI", "D1_iid_rw2", "D3_bym2_rw1"), function(a)
             lapply(c("P1_last_value", "P3_last_3y", "M0_crude"), function(b) c(a, b))), recursive = FALSE))
pairs <- Filter(function(p) all(p %in% mods), pairs)
S8c <- purrr::map_dfr(pairs, function(p) pair_ci(p[1], p[2])) |> mutate(n_boot = REV$n_boot)
readr::write_csv(S8c, file.path(dir_rev, "S8c_bootstrap_rmse_differences.csv"))
say("Bootstrap (", REV$n_boot, " county-cluster replicates, seed ", REV$seed, "):")
for (k in seq_len(nrow(S8c)))
  say("   ", S8c$approach[k], " - ", S8c$reference[k], ": ", sprintf("%+.3f", S8c$rmse_diff_pp[k]),
      " pp (95% ", sprintf("%+.3f", S8c$lo95[k]), " to ", sprintf("%+.3f", S8c$hi95[k]),
      "); reduction ", sprintf("%.1f", S8c$reduction_pct[k]), "% (",
      sprintf("%.1f", S8c$red_lo95[k]), " to ", sprintf("%.1f", S8c$red_hi95[k]), ")")

## (iv) margin sweep (adequacy thresholds were met by every model, so the rule
## reduces to margin + complexity ordering; checked here, not assumed)
tab2 <- readr::read_csv(file.path(cfg$dir_out, "model_comparison.csv"), show_col_types = FALSE)
adequate <- tab2$model[tab2$model %in% SEL_SET &
                         tidyr::replace_na(tab2$calib_max_abs_gap_pp, 99) <= cfg$calib_gap_max_pp &
                         tidyr::replace_na(tab2$resid_moran_sig_years, 99) <= cfg$resid_sig_years_max]
say("Models meeting the calibration and residual-autocorrelation thresholds: ",
    paste(adequate, collapse = ", "), " | calibration gaps ",
    paste(f3(tab2$calib_max_abs_gap_pp[tab2$model %in% SEL_SET]), collapse = ", "))
sweep <- function(col, design) purrr::map_dfr(REV$margins, function(mg) {
  d <- tab2 |> filter(model %in% adequate)
  best <- min(d[[col]])
  ok  <- SEL_SET[SEL_SET %in% d$model[d[[col]] <= best + mg]]
  tibble(design = design, margin_pp = mg, best_rmse = best, selected = ok[1],
         eligible = paste(ok, collapse = ", "))
})
S8d <- bind_rows(sweep("rmse_pp_A", "A"), sweep("rmse_pp_B", "B"))
readr::write_csv(S8d, file.path(dir_rev, "S8d_margin_sweep.csv"))
for (k in seq_len(nrow(S8d)))
  say("   margin ", sprintf("%.2f", S8d$margin_pp[k]), " Design ", S8d$design[k], ": ", S8d$selected[k])
say("Design B: binomial GLMM misses the 0.25 margin by ",
    sprintf("%.4f", tab2$rmse_pp_B[tab2$model == "M1_binomial_glmm"] -
              min(tab2$rmse_pp_B[tab2$model %in% SEL_SET]) - 0.25), " pp")

## noise-corrected MSE for every saved design (vii, part 1)
## E[(pred - obs)^2] = (pred - p)^2 + p(1-p)/N ; p(1-p)/(N-1) estimated without bias by obs(1-obs)/(N-1)
nc <- function(d) d |> mutate(sv_pp2 = 1e4 * observed * (1 - observed) / pmax(deno - 1, 1)) |>
  summarise(n = n(), mse_pp2 = mean(err_pp^2), sampling_var_pp2 = mean(sv_pp2),
            nc_mse_pp2 = mean(err_pp^2 - sv_pp2), .groups = "drop") |>
  mutate(nc_rmse_pp = sqrt(pmax(nc_mse_pp2, 0)))
cv_all <- bind_rows(cv |> select(design, model, deno, observed, err_pp),
                    pB |> mutate(design = "B") |> select(design, model, deno, observed, err_pp)) |>
  mutate(stratum = strat(deno))
S8f_nc <- bind_rows(cv_all |> group_by(design, model) |> nc() |> mutate(stratum = "All"),
                    cv_all |> group_by(design, model, stratum) |> nc())
readr::write_csv(S8f_nc, file.path(dir_rev, "S8f_noise_corrected_mse.csv"))
for (dsg in sort(unique(S8f_nc$design))) for (s in c("All", "N<50", "50<=N<100", "N>=100")) {
  x <- S8f_nc |> filter(design == dsg, stratum == s, model %in% c("M0_crude", "P1_last_value", SEL_SET)) |>
    mutate(model = factor(model, APPROACH_ORDER)) |> arrange(model)
  if (nrow(x)) say("Noise-corrected MSE, Design ", dsg, ", ", s, ": ",
                   paste0(x$model, " ", sprintf("%.3f", x$nc_mse_pp2), collapse = "; "),
                   " | sampling variance ", sprintf("%.3f", x$sampling_var_pp2[1]), " pp^2")
}

## =====================================================================
## (v)-(vii) INLA refits: earlier origins, two-year horizon, log predictive density
## =====================================================================
if (isTRUE(REV$run_refits)) {
  if (!requireNamespace("INLA", quietly = TRUE)) stop("INLA is required for the refits.")
  library(INLA)
  if (!is.null(cfg$inla_threads)) INLA::inla.setOption(num.threads = cfg$inla_threads)
  spatial <- readRDS(file.path(cfg$dir_out, "sc_spatial.rds"))
  graph_file <- spatial$graph_file
  if (!file.exists(graph_file)) graph_file <- file.path(cfg$dir_out, basename(graph_file))

  ## same specifications as 03_models.R
  f_M1 <- nume ~ 1 + f(county_id, model = "iid", hyper = pc_prec()) +
    f(year_id, model = "iid", hyper = pc_prec())
  f_M3 <- nume ~ 1 + f(county_id, model = "bym2", graph = graph_file, scale.model = TRUE,
                       constr = TRUE, hyper = pc_bym2()) +
    f(year_id, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = pc_prec())
  f_M4 <- update(f_M3, . ~ . + f(st_id, model = "iid", hyper = pc_prec()))
  f_D1 <- nume ~ 1 + f(county_id, model = "iid", hyper = pc_prec()) +
    f(year_id, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = pc_prec())
  f_D3 <- nume ~ 1 + f(county_id, model = "bym2", graph = graph_file, scale.model = TRUE,
                       constr = TRUE, hyper = pc_bym2()) +
    f(year_id, model = "rw1", scale.model = TRUE, constr = TRUE, hyper = pc_prec())
  specs <- list(M1_binomial_glmm = list(f = f_M1, fam = "binomial"),
                M2_betabinomial_glmm = list(f = f_M1, fam = "betabinomial"),
                M3_bym2_rw2 = list(f = f_M3, fam = "binomial"),
                M4_bym2_rw2_typeI = list(f = f_M4, fam = "binomial"),
                D1_iid_rw2 = list(f = f_D1, fam = "binomial"),
                D3_bym2_rw1 = list(f = f_D3, fam = "binomial"))

  fit_pred <- function(sp, dat, pos, y, N) {
    r <- try(INLA::inla(sp$f, family = sp$fam, data = as.data.frame(dat), Ntrials = deno,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute = list(return.marginals.predictor = TRUE),
      control.fixed = list(mean.intercept = 0, prec.intercept = 0.001, mean = 0, prec = 0.001),
      control.inla = list(strategy = "adaptive", int.strategy = "auto"), verbose = FALSE), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    marg <- r$marginals.fitted.values[pos]
    rho <- NA_real_
    if (sp$fam == "betabinomial") {
      k <- grep("overdispersion", rownames(r$summary.hyperpar), ignore.case = TRUE)
      rho <- r$summary.hyperpar$mean[k[1]]
    }
    lpd <- vapply(seq_along(pos), function(j) {
      fn <- if (is.na(rho)) function(p) stats::dbinom(y[j], N[j], pmin(pmax(p, 1e-12), 1 - 1e-12)) else
        function(p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-12); a <- p * (1 - rho) / rho; b <- (1 - p) * (1 - rho) / rho
          exp(lchoose(N[j], y[j]) + lbeta(y[j] + a, N[j] - y[j] + b) - lbeta(a, b)) }
      log(max(INLA::inla.emarginal(fn, marg[[j]]), 1e-300))
    }, numeric(1))
    tibble(pred = r$summary.fitted.values$mean[pos], lpd = lpd)
  }

  ## generic holdout: hold_idx = held-out rows; train = rows kept in the data
  ## (held-out responses set to NA); lag_from = years between last training year and target
  run_holdout <- function(hold_idx, keep, label, models, lag_from = 1, comparators = TRUE) {
    dat <- panel; dat$nume[hold_idx] <- NA; dat <- dat[keep, , drop = FALSE]
    pos <- match(hold_idx, which(keep)); y <- panel$nume[hold_idx]; N <- panel$deno[hold_idx]
    out <- list()
    for (nm in models) {
      pr <- fit_pred(specs[[nm]], dat, pos, y, N)
      if (is.null(pr)) { message("  ", label, "/", nm, " FAILED"); next }
      out[[nm]] <- mutate(pr, model = nm, row = hold_idx)
    }
    if (comparators) {
      tr <- dat[!is.na(dat$nume), ]
      cty <- tr |> group_by(county_id) |> summarise(p = sum(nume) / sum(deno), .groups = "drop")
      pc <- tibble(county_id = panel$county_id[hold_idx]) |> left_join(cty, by = "county_id") |>
        mutate(p = tidyr::replace_na(p, sum(tr$nume) / sum(tr$deno))) |> pull(p)
      out[["M0_crude"]] <- tibble(model = "M0_crude", row = hold_idx, pred = pc, lpd = NA_real_)
      if (lag_from > 0) {
        pp <- persist(hold_idx, lag_from)
        out[["P1_last_value"]] <- tibble(model = "P1_last_value", row = hold_idx, pred = pp$P1_last_value, lpd = NA_real_)
        out[["P3_last_3y"]] <- tibble(model = "P3_last_3y", row = hold_idx, pred = pp$P3_last_3y, lpd = NA_real_)
      }
    }
    message("  ", label, " done (", format(Sys.time(), "%H:%M:%S"), ")")
    bind_rows(out) |> mutate(setting = label)
  }
  ## fit through year t - h and predict year t; responses of the years in between
  ## are unknown at release time, so they are set to NA as well
  forward <- function(targets, h, label, models) purrr::map_dfr(targets, function(t) {
    hold <- which(panel$year == t)
    gap  <- which(panel$year > t - h & panel$year < t)
    r <- run_holdout(c(hold, gap), panel$year <= t, paste0(label, " ", t), models, lag_from = h)
    filter(r, row %in% hold)
  })
  fwd_models <- names(specs)
  R_early <- forward(REV$early_years, 1, "early origin", fwd_models)
  R_h2    <- forward(cfg$forward_years, 2, "two-year horizon", fwd_models)
  R_B     <- forward(cfg$forward_years, 1, "Design B refit", SEL_SET)

  ## Design A folds rebuilt exactly as in 03_models.R
  set.seed(cfg$cv_seed)
  strata <- dplyr::ntile(panel$deno, 5); fold <- integer(nrow(panel))
  for (s in unique(strata)) { idx <- which(strata == s)
    fold[idx] <- sample(rep_len(seq_len(cfg$cv_folds), length(idx))) }
  R_A <- purrr::map_dfr(sort(unique(fold)), function(k)
    run_holdout(which(fold == k), rep(TRUE, nrow(panel)), paste0("Design A fold ", k), SEL_SET, 0, FALSE))
  R_C <- if (isTRUE(REV$lpd_designC)) purrr::map_dfr(sort(unique(panel$county_id)), function(c)
    run_holdout(which(panel$county_id == c), rep(TRUE, nrow(panel)), paste0("Design C county ", c),
                SEL_SET, 0, FALSE)) else NULL

  add_obs <- function(d) d |> mutate(observed = panel$encounter_rate[row], deno = panel$deno[row],
                                     year = panel$year[row], err_pp = 100 * (pred - observed),
                                     stratum = strat(deno))
  R_early <- add_obs(R_early); R_h2 <- add_obs(R_h2); R_B <- add_obs(R_B); R_A <- add_obs(R_A)
  if (!is.null(R_C)) R_C <- add_obs(R_C)

  ## check: the Design B refit reproduces the saved predictions
  chk <- R_B |> inner_join(cv |> filter(design == "B") |> select(model, row, pred0 = pred), by = c("model", "row"))
  say("\nCheck: Design B refit vs saved predictions, max |difference| = ",
      formatC(max(abs(chk$pred - chk$pred0)), format = "e", digits = 1), " (proportion scale)")

  ## (v) and (vi)
  S8e <- bind_rows(R_early |> group_by(setting = "early origins 2010-2014", model) |> met(),
                   R_h2 |> group_by(setting = "two-year horizon 2019-2023", model) |> met(),
                   R_early |> group_by(setting = "early origins 2010-2014", model, year) |> met(),
                   R_h2 |> group_by(setting = "two-year horizon 2019-2023", model, year) |> met()) |>
    mutate(model = factor(model, APPROACH_ORDER)) |> arrange(setting, is.na(year) == FALSE, year, model)
  readr::write_csv(S8e, file.path(dir_rev, "S8e_early_origins_two_year_horizon.csv"))
  for (st in unique(S8e$setting)) {
    x <- S8e |> filter(setting == st, is.na(year))
    ref <- x$rmse_pp[x$model == "M1_binomial_glmm"]
    say("\n", st, " RMSE: ", paste0(x$model, " ", f3(x$rmse_pp), collapse = "; "))
    say("   reduction vs binomial GLMM: ", paste0(x$model, " ", sprintf("%.1f", 100 * (1 - x$rmse_pp / ref)), "%",
                                              collapse = "; "))
  }
  ## bootstrap for the early-origin and horizon contrasts
  boot_contrast <- function(R, a, b) {
    W <- R |> filter(model %in% c(a, b)) |> select(model, row, err_pp) |>
      pivot_wider(names_from = model, values_from = err_pp) |> mutate(cid = panel$county_id[row])
    ix <- split(seq_len(nrow(W)), W$cid)
    set.seed(REV$seed)
    bb <- replicate(REV$n_boot, { ii <- unlist(ix[sample(names(ix), replace = TRUE)], use.names = FALSE)
      sqrt(mean(W[[a]][ii]^2)) - sqrt(mean(W[[b]][ii]^2)) })
    c(diff = sqrt(mean(W[[a]]^2)) - sqrt(mean(W[[b]]^2)), lo = unname(quantile(bb, .025)), hi = unname(quantile(bb, .975)))
  }
  for (pair in list(c("M3_bym2_rw2", "M1_binomial_glmm"), c("M3_bym2_rw2", "P1_last_value"),
                    c("M3_bym2_rw2", "P3_last_3y"), c("D3_bym2_rw1", "M3_bym2_rw2"))) {
    for (nmR in c("early", "h2")) {
      R <- if (nmR == "early") R_early else R_h2
      if (!all(pair %in% R$model)) next
      bc <- boot_contrast(R, pair[1], pair[2])
      say("   ", nmR, ": ", pair[1], " - ", pair[2], " = ", sprintf("%+.3f", bc["diff"]), " pp (95% ",
          sprintf("%+.3f", bc["lo"]), " to ", sprintf("%+.3f", bc["hi"]), ")")
    }
  }

  ## (vii) log predictive density by design and stratum
  L <- bind_rows(mutate(R_A, design = "A"), mutate(R_B, design = "B"),
                 if (!is.null(R_C)) mutate(R_C, design = "C")) |> filter(!is.na(lpd))
  S8f_lpd <- bind_rows(L |> group_by(design, model) |> summarise(n = n(), mean_lpd = mean(lpd), .groups = "drop") |>
                         mutate(stratum = "All"),
                       L |> group_by(design, model, stratum) |>
                         summarise(n = n(), mean_lpd = mean(lpd), .groups = "drop"))
  readr::write_csv(S8f_lpd, file.path(dir_rev, "S8f_log_predictive_density.csv"))
  say("")
  for (dsg in sort(unique(S8f_lpd$design))) for (s in c("All", "N<50", "50<=N<100", "N>=100")) {
    x <- S8f_lpd |> filter(design == dsg, stratum == s) |> mutate(model = factor(model, APPROACH_ORDER)) |> arrange(model)
    say("Mean log predictive density, Design ", dsg, ", ", s, ": ",
        paste0(x$model, " ", sprintf("%.4f", x$mean_lpd), collapse = "; "))
  }
}

## =====================================================================
## (viii) level shift vs Design B gain across sensitivity variants
## =====================================================================
s4 <- file.path(cfg$dir_out, "tableS4_encounter_rate_across_variants.csv")
s3 <- file.path(cfg$dir_out, "tableS3_metrics_across_variants.csv")
if (file.exists(s4) && file.exists(s3)) {
  lev <- readr::read_csv(s4, show_col_types = FALSE) |>
    mutate(per = ifelse(year <= 2018, "pre", ifelse(year >= 2019, "post", NA))) |>
    group_by(variant, per) |> summarise(ml = mean(qlogis(wrate_pct / 100)), .groups = "drop") |>
    pivot_wider(names_from = per, values_from = ml) |> mutate(level_shift = pre - post)
  gain <- readr::read_csv(s3, show_col_types = FALSE) |>
    filter(model %in% c("M1_binomial_glmm", "M4_bym2_rw2_typeI")) |>
    select(variant, model, rmse_pp_B) |> pivot_wider(names_from = model, values_from = rmse_pp_B) |>
    mutate(gain_pct = 100 * (1 - M4_bym2_rw2_typeI / M1_binomial_glmm))
  S6 <- inner_join(select(lev, variant, level_shift), select(gain, variant, gain_pct), by = "variant") |>
    mutate(variant = sub("^sens_", "", variant)) |> arrange(desc(level_shift))
  readr::write_csv(S6, file.path(dir_rev, "figS6_level_shift_vs_gain.csv"))
  say("\n(viii) level shift (mean logit 2005-2018 minus 2019-2023) and Design B gain of the interaction model:")
  for (k in seq_len(nrow(S6))) say("   ", S6$variant[k], ": shift ", f3(S6$level_shift[k]),
                                   " | gain ", sprintf("%.1f", S6$gain_pct[k]), "%")
  say("   Spearman correlation across variants: ",
      sprintf("%.2f", suppressWarnings(cor(S6$level_shift, S6$gain_pct, method = "spearman"))))
  ## variants that coincide with the primary analysis (priors, adjacency, matching) are not labelled
  lab <- c(primary = "Primary analysis*", pheno_with_alcohol = "Alcohol included",
           pheno_revised = "Dependence and unspecified in list", add_dependence = "Dependence forced in",
           denom_encounter = "Encounter denominator", county_current = "Current residence",
           era_off = "No coding-era enforcement", post_dx_only = "Encounters after diagnosis")
  S6$label <- ifelse(S6$variant %in% names(lab), lab[S6$variant], NA)
  g <- ggplot(S6, aes(level_shift, gain_pct)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(size = 2.2, colour = "#1f4e79") +
    geom_text(aes(label = label), size = 2.7, hjust = -0.1, vjust = 0.5, na.rm = TRUE) +
    scale_x_continuous(limits = c(0, 0.45)) +
    labs(x = "Shift in statewide series (mean logit, 2005-2018 minus 2019-2023)",
         y = "Design B RMSE reduction (%)\ninteraction model vs binomial GLMM",
         caption = "*Overlapping the primary analysis: rook contiguity, prefix matching, and the three prior variants.") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0, size = 7))
  ggsave(file.path(dir_rev, "figS6_level_shift.png"), g, width = 7, height = 4.6, dpi = 300, bg = "white")
} else say("\n(viii) skipped: run_all.R sensitivity tables not found")

say("\nRuntime: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
writeLines(capture.output(sessionInfo()), file.path(dir_rev, "sessionInfo_05.txt"))
message("\nDone. Send this file back: ", normalizePath(summary_file))
