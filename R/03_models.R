## =====================================================================
## 03_models.R -- fit M0-M4 (+ decomposition models), run holdout Designs
##                A/B/C, compute all comparison metrics
##
##   in : output/panel_<...>_indexed.csv, output/sc_adj_<adjacency>.graph
##   out: output/model_comparison.csv          (Table 2, primary models)
##        output/model_comparison_all.csv      (incl. decomposition models)
##        output/cv_predictions.csv
##        output/table2b_cv_by_stratum.csv
##        output/calibration_designA.csv
##        output/residual_moran.csv
##        output/selection_rule_trace.csv          (prespecified, Design A)
##        output/selection_rule_trace_designB.csv  (secondary, Design B)
##        output/designB_relative_reduction.csv
##        output/precision_summary.csv         (credible-interval width)
##        output/fixed_effects.csv             (intercept)
##        output/hyperparameters.csv           (precision scale)
##        output/variance_components_sd.csv    (SD scale)
##        output/phi_marginal.csv              (BYM2 mixing parameter)
##        output/pearson_residuals.csv
##        output/fitted_selected_model.csv
##        output/models_<phenotype>_<denominator>.rds
## =====================================================================
if (!exists("cfg", envir = globalenv())) source(file.path("R", "00_config.R"))
cfg <- get("cfg", envir = globalenv())
suppressWarnings(suppressMessages({ library(spdep); library(sf) }))

if (!requireNamespace("INLA", quietly = TRUE))
  stop("INLA is required.\n",
       'install.packages("INLA", repos = c(getOption("repos"), ',
       'INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)')
library(INLA)

## defaults, so the script also runs with an older 00_config.R
if (is.null(cfg$calib_gap_max_pp))    cfg$calib_gap_max_pp    <- 1.5
if (is.null(cfg$resid_sig_years_max)) cfg$resid_sig_years_max <- 2
if (is.null(cfg$run_decomposition))   cfg$run_decomposition   <- TRUE
if (is.null(cfg$n_pred_draws))        cfg$n_pred_draws        <- 1000
if (is.null(cfg$run_glmmtmb))         cfg$run_glmmtmb         <- FALSE
if (is.null(cfg$pi_betabinomial))     cfg$pi_betabinomial     <- FALSE
if (!is.null(cfg$inla_threads)) INLA::inla.setOption(num.threads = cfg$inla_threads)
log_step("INLA ", as.character(utils::packageVersion("INLA")),
         " | num.threads = ", INLA::inla.getOption("num.threads"),
         " | detected CPU cores = ", parallel::detectCores())

panel <- readr::read_csv(file.path(cfg$dir_out,
                                   sprintf("panel_%s_%s_indexed.csv", cfg$phenotype, cfg$denominator)),
                         col_types = readr::cols(fips = readr::col_character()))
spatial    <- readRDS(file.path(cfg$dir_out, "sc_spatial.rds"))
graph_file <- spatial$graph_file
nb <- spatial$nb; sc <- spatial$sc

check_that(all(panel$deno > 0), "every county-year has a positive denominator")

## index copies -- INLA needs a distinct index vector per f() term
panel <- panel |>
  mutate(county_id2 = county_id, year_id2 = year_id,
         st_id = as.integer(interaction(county_id, year_id, drop = TRUE)))

EXCEED_THRESHOLD <- 0.20

## =====================================================================
## Model specifications
## =====================================================================
## The intercept alpha is the only fixed effect in every model.

## ---- prespecified selection set (M1-M4; M0 is the crude comparator) ----
f_M1 <- nume ~ 1 +
  f(county_id, model = "iid", hyper = pc_prec()) +
  f(year_id,   model = "iid", hyper = pc_prec())

f_M2 <- f_M1                       # same linear predictor, beta-binomial likelihood

f_M3 <- nume ~ 1 +
  f(county_id, model = "bym2", graph = graph_file,
    scale.model = TRUE, constr = TRUE, hyper = pc_bym2()) +
  f(year_id, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = pc_prec())

f_M4 <- nume ~ 1 +
  f(county_id, model = "bym2", graph = graph_file,
    scale.model = TRUE, constr = TRUE, hyper = pc_bym2()) +
  f(year_id, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = pc_prec()) +
  f(st_id, model = "iid", hyper = pc_prec())          # Knorr-Held type I

model_specs <- list(
  M1_binomial_glmm      = list(formula = f_M1, family = "binomial"),
  M2_betabinomial_glmm  = list(formula = f_M2, family = "betabinomial"),
  M3_bym2_rw2           = list(formula = f_M3, family = "binomial"),
  M4_bym2_rw2_typeI     = list(formula = f_M4, family = "binomial")
)

## ---- decomposition models (supplementary; not in the selection set) ----
## D1: temporal structure only    -> is the Design B gain temporal?
## D2: spatial structure only     -> does BYM2 alone help?
## D3: BYM2 + RW1 instead of RW2  -> is second order needed?
f_D1 <- nume ~ 1 +
  f(county_id, model = "iid", hyper = pc_prec()) +
  f(year_id, model = "rw2", scale.model = TRUE, constr = TRUE, hyper = pc_prec())

f_D2 <- nume ~ 1 +
  f(county_id, model = "bym2", graph = graph_file,
    scale.model = TRUE, constr = TRUE, hyper = pc_bym2()) +
  f(year_id, model = "iid", hyper = pc_prec())

f_D3 <- nume ~ 1 +
  f(county_id, model = "bym2", graph = graph_file,
    scale.model = TRUE, constr = TRUE, hyper = pc_bym2()) +
  f(year_id, model = "rw1", scale.model = TRUE, constr = TRUE, hyper = pc_prec())

decomp_specs <- list(
  D1_iid_rw2      = list(formula = f_D1, family = "binomial"),
  D2_bym2_iidyear = list(formula = f_D2, family = "binomial"),
  D3_bym2_rw1     = list(formula = f_D3, family = "binomial")
)

all_specs <- if (isTRUE(cfg$run_decomposition)) c(model_specs, decomp_specs) else model_specs
model_set <- c(setNames(rep("primary", length(model_specs)), names(model_specs)),
               setNames(rep("decomposition", length(decomp_specs)), names(decomp_specs)),
               M0_crude = "primary")

## =====================================================================
## Fitting wrapper
## =====================================================================
fit_one <- function(spec, dat, compute_extras = TRUE, marginals = compute_extras) {
  dat <- as.data.frame(dat)
  t0 <- proc.time()[["elapsed"]]
  res <- try(INLA::inla(
    spec$formula, family = spec$family, data = dat, Ntrials = deno,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute   = list(dic = compute_extras, waic = compute_extras,
                             cpo = compute_extras, config = compute_extras,
                             return.marginals.predictor = marginals),
    control.fixed     = list(mean.intercept = 0, prec.intercept = 0.001,
                             mean = 0, prec = 0.001),
    control.inla      = list(strategy = "adaptive", int.strategy = "auto"),
    verbose = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) return(list(ok = FALSE, msg = as.character(res)))
  list(ok = TRUE, fit = res, seconds = proc.time()[["elapsed"]] - t0)
}

## 95% posterior predictive interval for an observed proportion Y/N:
## draw p from its posterior marginal, then Y ~ Binomial(N, p).
## For the beta-binomial model this ignores the extra-binomial term unless
## cfg$pi_betabinomial is TRUE. The omission is not uniformly negligible: the
## beta-binomial variance is the binomial variance times 1 + (N - 1) * rho, so
## with rho ~ 3e-4 the factor is <= 1.03 when N < 100 but ~2 for the largest
## county-years (N ~ 3,500). Only the beta-binomial model's intervals are
## affected; point predictions, RMSE/MAE and the selection rule are not.
pred_interval <- function(marg_list, N, ndraw, rho_marg = NULL) {
  out <- vapply(seq_along(marg_list), function(k) {
    p <- INLA::inla.rmarginal(ndraw, marg_list[[k]])
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    if (!is.null(rho_marg)) {
      rho <- pmin(pmax(INLA::inla.rmarginal(ndraw, rho_marg), 1e-8), 1 - 1e-6)
      p <- stats::rbeta(ndraw, p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho)
    }
    y <- stats::rbinom(ndraw, N[k], p)
    stats::quantile(y / N[k], c(.025, .975), names = FALSE)
  }, numeric(2))
  t(out)
}
rho_marginal <- function(fit) {
  hp <- fit$marginals.hyperpar
  k  <- grep("overdispersion", names(hp), ignore.case = TRUE)
  if (length(k)) hp[[k[1]]] else NULL
}

## =====================================================================
## PART 1 -- full-data fits
## =====================================================================
full_fits <- list(); ic_rows <- list()
for (nm in names(all_specs)) {
  log_step("fitting (full data): ", nm)
  r <- fit_one(all_specs[[nm]], panel)
  if (!r$ok) { log_step("  FAILED: ", substr(r$msg, 1, 200)); next }
  f <- r$fit
  ic_rows[[nm]] <- tibble(
    model = nm, set = model_set[[nm]],
    waic  = f$waic$waic, p_waic = f$waic$p.eff,
    dic   = f$dic$dic,   p_dic  = f$dic$p.eff,
    mean_log_cpo = mean(log(pmax(f$cpo$cpo, 1e-12)), na.rm = TRUE),
    cpo_failures = sum(f$cpo$failure > 0, na.rm = TRUE),
    n_params = length(f$summary.hyperpar$mean) + nrow(f$summary.fixed),
    seconds = r$seconds
  )
  full_fits[[nm]] <- f
}
if (length(full_fits) == 0)
  stop("Every model failed to fit. Rerun one by one: fit_one(all_specs[[1]], panel)$msg")
ic_tab <- bind_rows(ic_rows)

## ---- residual spatial autocorrelation, by year, for every fit ----------
## moran.mc permutes the residual values across the 46 counties within each
## year (spatial permutation within time).
lw <- spdep::nb2listw(nb, style = "W")
resid_moran <- purrr::map_dfr(names(full_fits), function(nm) {
  fitted_p <- full_fits[[nm]]$summary.fitted.values$mean
  d <- panel |> mutate(resid = encounter_rate - fitted_p)
  purrr::map_dfr(sort(unique(d$year)), function(y) {
    dy <- d |> filter(year == y) |> arrange(match(fips, sc$fips))
    mc <- spdep::moran.mc(dy$resid, lw, nsim = 999, zero.policy = TRUE)
    tibble(model = nm, year = y, morans_i = unname(mc$statistic), p = mc$p.value)
  })
}) |> group_by(model) |> mutate(p_fdr = p.adjust(p, "BH")) |> ungroup()
readr::write_csv(resid_moran, file.path(cfg$dir_out, "residual_moran.csv"))

resid_summary <- resid_moran |> group_by(model) |>
  summarise(resid_moran_sig_years = sum(p_fdr < 0.05),
            resid_moran_max = max(morans_i), .groups = "drop")

## ---- precision of fitted values: 95% credible-interval width -----------
strat <- function(n) cut(n, c(-Inf, 50, 100, Inf), right = FALSE,
                         labels = c("N<50", "50<=N<100", "N>=100"))
precision_summary <- purrr::imap_dfr(full_fits, function(f, nm) {
  sfv <- f$summary.fitted.values
  d <- tibble(model = nm, stratum = as.character(strat(panel$deno)),
              width_pp = 100 * (sfv$`0.975quant` - sfv$`0.025quant`))
  bind_rows(mutate(d, stratum = "All"), d) |>
    group_by(model, stratum) |>
    summarise(mean_cri_width_pp = mean(width_pp),
              median_cri_width_pp = stats::median(width_pp), .groups = "drop")
})
readr::write_csv(precision_summary, file.path(cfg$dir_out, "precision_summary.csv"))

## ---- fixed effect (intercept) and covariance parameters ----------------
purrr::imap_dfr(full_fits, ~ tibble::rownames_to_column(.x$summary.fixed, "param") |>
                  dplyr::mutate(model = .y)) |>
  readr::write_csv(file.path(cfg$dir_out, "fixed_effects.csv"))

purrr::imap_dfr(full_fits, ~ tibble::rownames_to_column(.x$summary.hyperpar, "param") |>
                  dplyr::mutate(model = .y)) |>
  readr::write_csv(file.path(cfg$dir_out, "hyperparameters.csv"))

## precisions re-expressed as standard deviations on the logit scale
sd_tab <- purrr::imap_dfr(full_fits, function(f, nm) {
  hp <- f$marginals.hyperpar
  hp <- hp[grepl("^Precision", names(hp))]
  purrr::imap_dfr(hp, function(m, pn) {
    z <- INLA::inla.zmarginal(INLA::inla.tmarginal(function(x) 1 / sqrt(x), m),
                              silent = TRUE)
    tibble(model = nm, param = sub("^Precision", "SD", pn),
           mean = z$mean, sd = z$sd, q025 = z$quant0.025,
           q50 = z$quant0.5, q975 = z$quant0.975)
  })
})
readr::write_csv(sd_tab, file.path(cfg$dir_out, "variance_components_sd.csv"))

## posterior density of the BYM2 mixing parameter (skewness check)
phi_tab <- purrr::imap_dfr(full_fits, function(f, nm) {
  m <- f$marginals.hyperpar[["Phi for county_id"]]
  if (is.null(m)) return(NULL)
  sm <- INLA::inla.smarginal(m)
  tibble(model = nm, x = sm$x, density = sm$y)
})
if (nrow(phi_tab)) readr::write_csv(phi_tab, file.path(cfg$dir_out, "phi_marginal.csv"))

## Pearson residuals: the binomial likelihood should absorb the funnel shape
## seen in raw shrinkage (sampling variance proportional to 1/N)
pearson <- purrr::imap_dfr(full_fits[intersect(names(full_fits), names(model_specs))],
  function(f, nm) {
    p <- f$summary.fitted.values$mean
    panel |> transmute(model = nm, county, year, deno, nume, fitted = p,
                       raw_pp  = 100 * (nume / deno - p),
                       pearson = (nume - deno * p) / sqrt(deno * p * (1 - p)))
  })
readr::write_csv(pearson, file.path(cfg$dir_out, "pearson_residuals.csv"))

## =====================================================================
## PART 2 -- holdout designs
## =====================================================================
make_folds <- function(design) {
  set.seed(cfg$cv_seed)
  if (design == "A") {                      # random county-year, denominator-stratified
    strata <- dplyr::ntile(panel$deno, 5)
    fold <- integer(nrow(panel))
    for (s in unique(strata)) {
      idx <- which(strata == s)
      fold[idx] <- sample(rep_len(seq_len(cfg$cv_folds), length(idx)))
    }
    split(seq_len(nrow(panel)), fold)
  } else if (design == "B") {               # forward in time
    setNames(lapply(cfg$forward_years, function(y) which(panel$year == y)),
             paste0("y", cfg$forward_years))
  } else if (design == "C") {               # leave one county out
    setNames(lapply(sort(unique(panel$county_id)),
                    function(c) which(panel$county_id == c)),
             paste0("c", sort(unique(panel$county_id))))
  }
}

## Design B trains only on years strictly before the target year.
train_mask <- function(design, hold_idx) {
  if (design != "B") return(rep(TRUE, nrow(panel)))
  target_year <- unique(panel$year[hold_idx])
  panel$year < target_year | seq_len(nrow(panel)) %in% hold_idx
}

## Crude comparator (M0). It uses TRAINING information only; the held-out
## numerator is never read. A held-out county-year is predicted by that
## county's pooled training proportion (sum of numerators / sum of
## denominators over its training county-years); when the county has no
## training data (Design C) the statewide training proportion is used.
##   Design A: county's other years (held-out cells masked)
##   Design B: county's years before the target year
##   Design C: statewide proportion from the other 45 counties
predict_crude <- function(train, hold_idx) {
  cty <- train |> filter(!is.na(nume)) |> group_by(county_id) |>
    summarise(p = sum(nume) / sum(deno), .groups = "drop")
  overall <- with(filter(train, !is.na(nume)), sum(nume) / sum(deno))
  tibble(row = hold_idx, county_id = panel$county_id[hold_idx]) |>
    left_join(cty, by = "county_id") |>
    mutate(pred = tidyr::replace_na(p, overall)) |> pull(pred)
}

run_design <- function(design) {
  folds <- make_folds(design)
  purrr::imap_dfr(folds, function(hold_idx, fold_name) {
    keep <- train_mask(design, hold_idx)
    out <- list()
    N_hold <- panel$deno[hold_idx]

    ## crude comparator; plug-in binomial interval around its prediction
    tr <- panel; tr$nume[hold_idx] <- NA
    tr <- tr[keep, , drop = FALSE]
    pc <- predict_crude(tr, hold_idx)
    out[["M0_crude"]] <- tibble(model = "M0_crude", row = hold_idx, pred = pc,
      pi_lo = stats::qbinom(.025, N_hold, pc) / N_hold,
      pi_hi = stats::qbinom(.975, N_hold, pc) / N_hold)

    ## model-based approaches. Held-out responses are set to NA but the rows
    ## stay in the data, so the adjacency graph and the time index are intact
    ## and INLA predicts the held-out cells from the posterior (no na.omit).
    for (nm in names(all_specs)) {
      dat <- panel; dat$nume[hold_idx] <- NA
      dat <- dat[keep, , drop = FALSE]
      r <- fit_one(all_specs[[nm]], dat, compute_extras = FALSE, marginals = TRUE)
      if (!r$ok) { log_step("  ", design, "/", fold_name, "/", nm, " FAILED"); next }
      pos <- match(hold_idx, which(keep))
      set.seed(cfg$cv_seed + match(fold_name, names(folds)) * 100 +
                 match(nm, names(all_specs)))
      rho_m <- if (isTRUE(cfg$pi_betabinomial) && all_specs[[nm]]$family == "betabinomial")
        rho_marginal(r$fit) else NULL
      pi <- pred_interval(r$fit$marginals.fitted.values[pos], N_hold, cfg$n_pred_draws, rho_m)
      out[[nm]] <- tibble(model = nm, row = hold_idx,
                          pred  = r$fit$summary.fitted.values$mean[pos],
                          pi_lo = pi[, 1], pi_hi = pi[, 2])
    }
    log_step("design ", design, " fold ", fold_name, " done")
    bind_rows(out) |> mutate(design = design, fold = fold_name)
  })
}

cv <- run_design("A")
if (isTRUE(cfg$run_designB)) cv <- bind_rows(cv, run_design("B"))
if (isTRUE(cfg$run_designC)) cv <- bind_rows(cv, run_design("C"))

cv <- cv |>
  mutate(observed = panel$encounter_rate[row], deno = panel$deno[row],
         year = panel$year[row], county = panel$county[row],
         err_pp = 100 * (pred - observed),
         covered = observed >= pi_lo & observed <= pi_hi,
         stratum = strat(deno),
         set = unname(model_set[model]))
readr::write_csv(cv, file.path(cfg$dir_out, "cv_predictions.csv"))

## =====================================================================
## PART 3 -- metrics and selection rules
## =====================================================================
metrics <- function(d) tibble(
  n          = nrow(d),
  rmse_pp    = sqrt(mean(d$err_pp^2)),
  mae_pp     = mean(abs(d$err_pp)),
  rmse_wt    = sqrt(stats::weighted.mean(d$err_pp^2, d$deno)),
  bias_pp    = mean(d$err_pp),
  cover95    = mean(d$covered),
  piwidth_pp = 100 * mean(d$pi_hi - d$pi_lo)
)

cv_overall <- cv |> group_by(design, model) |> group_modify(~ metrics(.x)) |> ungroup()
cv_stratum <- cv |> group_by(design, model, stratum) |> group_modify(~ metrics(.x)) |> ungroup()
readr::write_csv(cv_stratum, file.path(cfg$dir_out, "table2b_cv_by_stratum.csv"))

## calibration on Design A: observed vs predicted by predicted decile
calib <- cv |> filter(design == "A") |> group_by(model) |>
  mutate(dec = dplyr::ntile(pred, 10)) |> group_by(model, dec) |>
  summarise(pred_mean = 100 * mean(pred),
            obs_mean  = 100 * stats::weighted.mean(observed, deno),
            n = n(), .groups = "drop")
readr::write_csv(calib, file.path(cfg$dir_out, "calibration_designA.csv"))
calib_summary <- calib |> group_by(model) |>
  summarise(calib_max_abs_gap_pp = max(abs(obs_mean - pred_mean)), .groups = "drop")

wide <- cv_overall |>
  select(design, model, rmse_pp, mae_pp, cover95, piwidth_pp) |>
  pivot_wider(names_from = design, values_from = c(rmse_pp, mae_pp, cover95, piwidth_pp))

table_all <- wide |>
  mutate(set = unname(model_set[model])) |>
  left_join(select(ic_tab, -set), by = "model") |>
  left_join(resid_summary, by = "model") |>
  left_join(calib_summary, by = "model") |>
  left_join(precision_summary |> filter(stratum == "All") |>
              select(model, mean_cri_width_pp), by = "model") |>
  arrange(set != "primary", rmse_pp_A)
readr::write_csv(table_all, file.path(cfg$dir_out, "model_comparison_all.csv"))

table2 <- table_all |> filter(set == "primary") |> select(-set)
readr::write_csv(table2, file.path(cfg$dir_out, "model_comparison.csv"))
print(as.data.frame(table2))

## ---- selection rules ----------------------------------------------------
complexity <- c(M0_crude = 0, M1_binomial_glmm = 1, M2_betabinomial_glmm = 2,
                M3_bym2_rw2 = 3, M4_bym2_rw2_typeI = 4)

apply_rule <- function(rmse_col) {
  ## reference = best MODEL-BASED approach; the crude comparator is a benchmark
  ## only (in the published run it had the highest RMSE, so this is identical
  ## to taking the minimum over all five rows)
  best <- min(table2[[rmse_col]][table2$model != "M0_crude"], na.rm = TRUE)
  table2 |>
    mutate(best_rmse     = best,
           within_margin = .data[[rmse_col]] <= best + cfg$rmse_margin_pp,
           calib_ok      = tidyr::replace_na(calib_max_abs_gap_pp, 99) <= cfg$calib_gap_max_pp,
           resid_ok      = tidyr::replace_na(resid_moran_sig_years, 99) <= cfg$resid_sig_years_max,
           ## the crude comparator has no residual diagnostic, so it cannot pass
           eligible      = within_margin & calib_ok & resid_ok & model != "M0_crude",
           complexity    = complexity[model]) |>
    arrange(desc(eligible), complexity)
}

trace_A <- apply_rule("rmse_pp_A") |> mutate(analysis = "prespecified (Design A)")
selected <- trace_A$model[trace_A$eligible][1]
if (is.na(selected)) selected <- table2$model[which.min(table2$rmse_pp_A)]
readr::write_csv(trace_A, file.path(cfg$dir_out, "selection_rule_trace.csv"))

## Secondary operational analysis: the same margin applied to Design B.
## This was added after the primary comparison; it does not replace the rule
## above and is reported as secondary.
if ("rmse_pp_B" %in% names(table2)) {
  trace_B <- apply_rule("rmse_pp_B") |>
    mutate(analysis = "secondary, added during analysis (Design B)")
  readr::write_csv(trace_B, file.path(cfg$dir_out, "selection_rule_trace_designB.csv"))
  selected_B <- trace_B$model[trace_B$eligible][1]
  log_step("SECONDARY (Design B) -> preferred model: ", selected_B)

  ref <- table2$rmse_pp_B[table2$model == "M1_binomial_glmm"]
  red <- table_all |> filter(!is.na(rmse_pp_B)) |>
    transmute(model, set, rmse_pp_B,
              reduction_vs_binomial_glmm_pct = 100 * (1 - rmse_pp_B / ref))
  readr::write_csv(red, file.path(cfg$dir_out, "designB_relative_reduction.csv"))
}

log_step("SELECTION RULE (prespecified, Design A) -> ", selected,
         " | margin ", cfg$rmse_margin_pp, " pp | calibration <= ", cfg$calib_gap_max_pp,
         " pp | residual Moran significant in <= ", cfg$resid_sig_years_max, " years")

## =====================================================================
## PART 4 -- prediction file for the selected model
## =====================================================================
write_fitted <- function(mod, path) {
  f <- full_fits[[mod]]
  marg <- f$marginals.fitted.values
  exceed <- vapply(marg, function(m) 1 - INLA::inla.pmarginal(EXCEED_THRESHOLD, m), numeric(1))
  panel |>
    select(county, fips, year, deno, nume, crude_pct = encounter_pct) |>
    mutate(model = mod,
           smoothed_pct = 100 * f$summary.fitted.values$mean,
           lo95_pct     = 100 * f$summary.fitted.values$`0.025quant`,
           hi95_pct     = 100 * f$summary.fitted.values$`0.975quant`,
           exceed_prob  = exceed,
           exceed_threshold_pct = 100 * EXCEED_THRESHOLD,
           shrinkage_pp = crude_pct - smoothed_pct) |>
    readr::write_csv(path)
}
## the prespecified selection, plus BYM2+RW2 (preferred for prospective
## surveillance) so that maps can be drawn under either model
if (selected %in% names(full_fits))
  write_fitted(selected, file.path(cfg$dir_out, "fitted_selected_model.csv"))
for (mod in intersect(c("M1_binomial_glmm", "M3_bym2_rw2"), names(full_fits)))
  write_fitted(mod, file.path(cfg$dir_out, paste0("fitted_", mod, ".csv")))

suppressWarnings(saveRDS(list(full_fits = full_fits, table2 = table2, table_all = table_all,
             cv = cv, selected = selected, cfg = cfg),
        file.path(cfg$dir_out,
                  sprintf("models_%s_%s.rds", cfg$phenotype, cfg$denominator))))

## =====================================================================
## PART 5 -- optional frequentist cross-check (not reported)
## =====================================================================
if (isTRUE(cfg$run_glmmtmb) && requireNamespace("glmmTMB", quietly = TRUE)) {
  d <- panel |> mutate(county_f = factor(county_id), year_f = factor(year_id))
  m1 <- try(glmmTMB::glmmTMB(cbind(nume, deno - nume) ~ 1 + (1 | county_f) + (1 | year_f),
                             family = stats::binomial(), data = d), silent = TRUE)
  m2 <- try(glmmTMB::glmmTMB(cbind(nume, deno - nume) ~ 1 + (1 | county_f) + (1 | year_f),
                             family = glmmTMB::betabinomial(), data = d), silent = TRUE)
  sink(file.path(cfg$dir_log, "glmmTMB_crosscheck.txt"))
  if (!inherits(m1, "try-error")) print(summary(m1))
  if (!inherits(m2, "try-error")) print(summary(m2))
  sink()
}

## machine specification for the Methods (computing environment)
source(file.path(cfg$dir_code, "session_report.R"))
write_session_report(file.path(cfg$dir_log, "sessionInfo_03.txt"))
