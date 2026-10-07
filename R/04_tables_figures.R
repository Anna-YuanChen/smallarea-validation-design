## =====================================================================
## 04_tables_figures.R -- publication-ready tables and figures
##   in : output/model_comparison.csv, output/cv_predictions.csv,
##        output/fitted_selected_model.csv, output/sc_spatial.rds
##   out: output/table2_model_comparison_formatted.csv
##        output/table3_selected_estimates.csv
##        output/fig3_maps_crude_vs_smoothed.png
##        output/fig4_obs_vs_pred.png, output/fig4b_error_by_stratum.png
## =====================================================================
if (!exists("cfg", envir = globalenv())) source(file.path("R", "00_config.R"))
cfg <- get("cfg", envir = globalenv())
suppressWarnings(suppressMessages({ library(ggplot2); library(sf) }))

out <- cfg$dir_out
tab2 <- readr::read_csv(file.path(out, "model_comparison.csv"), show_col_types = FALSE)
cv   <- readr::read_csv(file.path(out, "cv_predictions.csv"),  show_col_types = FALSE)
spatial <- readRDS(file.path(out, "sc_spatial.rds")); sc <- spatial$sc

## Maps and the shrinkage figure are drawn under cfg$map_model. The default is
## BYM2+RW2, the model preferred for prospective surveillance; set it to
## "selected" to use the model chosen by the prespecified Design A rule.
map_model <- if (is.null(cfg$map_model)) "M3_bym2_rw2" else cfg$map_model
fit_path <- if (map_model == "selected") file.path(out, "fitted_selected_model.csv") else
  file.path(out, paste0("fitted_", map_model, ".csv"))
if (!file.exists(fit_path)) fit_path <- file.path(out, "fitted_selected_model.csv")
log_step("maps and shrinkage drawn under: ", basename(fit_path))
fitted <- if (file.exists(fit_path))
  readr::read_csv(fit_path, col_types = readr::cols(fips = readr::col_character())) else NULL

## ---- Table 2, formatted ----------------------------------------------
fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))

## reader-facing labels
model_labels <- c(M0_crude = "Crude", M1_binomial_glmm = "Binomial GLMM",
                  M2_betabinomial_glmm = "Beta-binomial GLMM", M3_bym2_rw2 = "BYM2 + RW2",
                  M4_bym2_rw2_typeI = "BYM2 + RW2 + interaction",
                  D1_iid_rw2 = "IID county + RW2", D2_bym2_iidyear = "BYM2 + IID year",
                  D3_bym2_rw1 = "BYM2 + RW1")
primary_models <- c("M0_crude", "M1_binomial_glmm", "M2_betabinomial_glmm",
                    "M3_bym2_rw2", "M4_bym2_rw2_typeI")
stratum_levels <- c("N<50", "50<=N<100", "N>=100")
stratum_labels <- c("N < 50", "50 \u2264 N < 100", "N \u2265 100")

## quick mode skips Design C, so only format the designs that actually ran
has <- function(col) col %in% names(tab2)
tab2_fmt <- tibble(Model = tab2$model)
if (has("waic"))         tab2_fmt$WAIC <- fmt(tab2$waic, 1)
if (has("dic"))          tab2_fmt$DIC  <- fmt(tab2$dic, 1)
if (has("mean_log_cpo")) tab2_fmt$`Mean log-CPO` <- fmt(tab2$mean_log_cpo, 3)
for (dsg in c("A", "B", "C")) {
  rc <- paste0("rmse_pp_", dsg); mc <- paste0("mae_pp_", dsg)
  if (has(rc)) tab2_fmt[[paste0("RMSE ", dsg, " (pp)")]] <- fmt(tab2[[rc]])
  if (has(mc)) tab2_fmt[[paste0("MAE ",  dsg, " (pp)")]] <- fmt(tab2[[mc]])
}
if (has("calib_max_abs_gap_pp"))
  tab2_fmt$`Max calibration gap (pp)` <- fmt(tab2$calib_max_abs_gap_pp, 2)
if (has("resid_moran_sig_years"))
  tab2_fmt$`Years with residual spatial autocorr.` <- tab2$resid_moran_sig_years
if (has("n_params")) tab2_fmt$Parameters <- tab2$n_params
if (has("n_hyper"))  tab2_fmt$Parameters <- tab2$n_hyper
if (has("seconds"))  tab2_fmt$`Fit time (s)`  <- fmt(tab2$seconds, 1)

readr::write_csv(tab2_fmt, file.path(out, "table2_model_comparison_formatted.csv"))

## ---- Table 3: selected-model estimates --------------------------------
if (!is.null(fitted)) {
  show_years <- c(min(cfg$years), 2015, max(cfg$years))
  tab3 <- fitted |>
    filter(year %in% show_years) |>
    transmute(County = county, Year = year, N = deno, Cases = nume,
              `Crude %` = round(crude_pct, 2),
              `Smoothed % (95% CrI)` = sprintf("%.2f (%.2f, %.2f)",
                                               smoothed_pct, lo95_pct, hi95_pct),
              `Shrinkage (pp)` = round(shrinkage_pp, 2),
              `P(prev > threshold)` = round(exceed_prob, 3)) |>
    arrange(Year, County)
  readr::write_csv(tab3, file.path(out, "table3_selected_estimates.csv"))

  ## ---- Figure 3: crude vs smoothed maps -------------------------------
  map_dat <- sc |>
    left_join(fitted |> filter(year %in% show_years), by = "fips") |>
    tidyr::pivot_longer(c(crude_pct, smoothed_pct),
                        names_to = "type", values_to = "pct") |>
    mutate(type = factor(type, c("crude_pct", "smoothed_pct"),
                         c("Crude", "Smoothed")))
  p3 <- ggplot(map_dat) +
    geom_sf(aes(fill = pct), colour = "white", linewidth = .15) +
    facet_grid(type ~ year) +
    scale_fill_viridis_c(option = "magma", direction = -1, name = "Encounter rate (%)") +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold")) +
    labs(title = NULL)
  ggsave(file.path(out, "fig3_maps_crude_vs_smoothed.png"), p3,
         width = 9, height = 6, dpi = 300)

  ## Shrinkage against denominator -- the visual argument for the whole paper
  p3b <- ggplot(fitted, aes(deno, shrinkage_pp)) +
    geom_hline(yintercept = 0, linewidth = .3) +
    geom_point(alpha = .35, size = 1) +
    scale_x_log10() +
    labs(x = "County-year denominator (log scale)",
         y = "Crude minus smoothed (percentage points)",
         title = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file.path(out, "fig3b_shrinkage_by_denominator.png"), p3b,
         width = 6.5, height = 4, dpi = 300)
}

## ---- Figure 4: observed vs predicted (Design A) ----------------------
designs_present <- unique(sub("^rmse_pp_", "", grep("^rmse_pp_", names(tab2), value = TRUE)))
log_step("holdout designs present in Table 2: ", paste(designs_present, collapse = ", "))

p4 <- cv |> filter(design == "A") |>
  ggplot(aes(100 * observed, 100 * pred)) +
  geom_abline(slope = 1, intercept = 0, linewidth = .3, colour = "grey40") +
  geom_point(alpha = .3, size = .8) +
  facet_wrap(~ model, nrow = 1) +
  coord_equal() +
  labs(x = "Observed encounter rate (%)", y = "Held-out predicted encounter rate (%)",
       title = "Design A: random county-year holdout") +
  theme_minimal(base_size = 10)
ggsave(file.path(out, "fig4_obs_vs_pred.png"), p4, width = 11, height = 3.4, dpi = 300)

p4b <- cv |>
  filter(model %in% primary_models) |>
  group_by(design, model, stratum) |>
  summarise(rmse = sqrt(mean(err_pp^2)), .groups = "drop") |>
  mutate(stratum = factor(stratum, stratum_levels, stratum_labels),
         model   = factor(model, primary_models, model_labels[primary_models])) |>
  ggplot(aes(stratum, rmse, fill = model)) +
  geom_col(position = position_dodge(.8), width = .75) +
  facet_wrap(~ design, labeller = labeller(design = c(
    A = "Design A: random county-year", B = "Design B: forward in time",
    C = "Design C: leave-one-county-out"))) +
  scale_fill_grey(start = 0.85, end = 0.15) +
  labs(x = "County-year denominator", y = "Holdout RMSE (percentage points)", fill = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "fig4b_error_by_stratum.png"), p4b, width = 9, height = 4.5, dpi = 300)


## ---- Figure S: posterior density of the BYM2 mixing parameter ---------
phi_path <- file.path(out, "phi_marginal.csv")
if (file.exists(phi_path)) {
  phi <- readr::read_csv(phi_path, show_col_types = FALSE) |>
    filter(model %in% c("M3_bym2_rw2", "M4_bym2_rw2_typeI")) |>
    mutate(model = factor(model, names(model_labels), model_labels))
  pphi <- ggplot(phi, aes(x, density, linetype = model)) +
    geom_line(linewidth = .7) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = expression(phi ~ "(spatially structured share of county variance)"),
         y = "Posterior density", linetype = NULL) +
    theme_classic(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(out, "figS_phi_posterior.png"), pphi, width = 6.2, height = 3.8, dpi = 600)
}

## ---- Figure S: raw versus Pearson residuals against the denominator -----
pr_path <- file.path(out, "pearson_residuals.csv")
if (file.exists(pr_path)) {
  pr <- readr::read_csv(pr_path, show_col_types = FALSE) |>
    filter(model == "M1_binomial_glmm") |>
    tidyr::pivot_longer(c(raw_pp, pearson), names_to = "type", values_to = "value") |>
    mutate(type = factor(type, c("raw_pp", "pearson"),
                         c("Raw residual (percentage points)", "Pearson residual")))
  ppr <- ggplot(pr, aes(deno, value)) +
    geom_hline(yintercept = 0, linewidth = .3) +
    geom_point(alpha = .35, size = .9) +
    scale_x_log10() +
    facet_wrap(~ type, scales = "free_y") +
    labs(x = "County-year denominator (log scale)", y = NULL) +
    theme_classic(base_size = 10)
  ggsave(file.path(out, "figS_pearson_residuals.png"), ppr, width = 8, height = 3.6, dpi = 600)
}

log_step("tables and figures written to ", out)
