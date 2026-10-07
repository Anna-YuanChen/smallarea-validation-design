## =====================================================================
## tools/pi_betabinomial_M2.R -- optional sensitivity check for Table S4
##
## In the published analysis, every model's 95% predictive interval drew
## Y ~ Binomial(N, p). For the beta-binomial GLMM (M2) this leaves out the
## overdispersion term. This script refits ONLY M2 under Designs A, B and C,
## using the same folds and seeds as 03_models.R. It reports coverage and width
## for both the binomial and the full beta-binomial predictive draws.
## Nothing else in the paper depends on it.
##
##   Rscript tools/pi_betabinomial_M2.R output
## (argument = the output folder of a completed run; takes a few minutes)
## =====================================================================
args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else "output"
suppressWarnings(suppressMessages({ library(dplyr); library(INLA) }))

cv_seed <- 20260820; cv_folds <- 10; forward_years <- 2019:2023; ndraw <- 1000
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

pf <- list.files(out, pattern = "^panel_.*_indexed\\.csv$", full.names = TRUE)[1]
panel <- readr::read_csv(pf, show_col_types = FALSE,
                         col_types = readr::cols(fips = readr::col_character()))
f_M2 <- nume ~ 1 + f(county_id, model = "iid", hyper = pc_prec) +
  f(year_id, model = "iid", hyper = pc_prec)

make_folds <- function(design) {
  set.seed(cv_seed)
  if (design == "A") {
    strata <- dplyr::ntile(panel$deno, 5); fold <- integer(nrow(panel))
    for (s in unique(strata)) {
      idx <- which(strata == s)
      fold[idx] <- sample(rep_len(seq_len(cv_folds), length(idx)))
    }
    split(seq_len(nrow(panel)), fold)
  } else if (design == "B") {
    setNames(lapply(forward_years, function(y) which(panel$year == y)), paste0("y", forward_years))
  } else {
    setNames(lapply(sort(unique(panel$county_id)), function(c) which(panel$county_id == c)),
             paste0("c", sort(unique(panel$county_id))))
  }
}
train_mask <- function(design, hold) {
  if (design != "B") return(rep(TRUE, nrow(panel)))
  panel$year < unique(panel$year[hold]) | seq_len(nrow(panel)) %in% hold
}

res <- list()
for (design in c("A", "B", "C")) {
  folds <- make_folds(design)
  for (fn in names(folds)) {
    hold <- folds[[fn]]; keep <- train_mask(design, hold)
    dat <- panel; dat$nume[hold] <- NA; dat <- as.data.frame(dat[keep, , drop = FALSE])
    fit <- INLA::inla(f_M2, family = "betabinomial", data = dat, Ntrials = deno,
                      control.predictor = list(compute = TRUE, link = 1),
                      control.compute = list(return.marginals.predictor = TRUE),
                      control.fixed = list(mean.intercept = 0, prec.intercept = 0.001),
                      control.inla = list(strategy = "adaptive", int.strategy = "auto"))
    pos <- match(hold, which(keep))
    rho_m <- fit$marginals.hyperpar[[grep("overdispersion", names(fit$marginals.hyperpar))[1]]]
    set.seed(cv_seed + match(fn, names(folds)) * 100 + 2)
    for (k in seq_along(pos)) {
      N <- panel$deno[hold[k]]; obs <- panel$encounter_rate[hold[k]]
      p   <- pmin(pmax(INLA::inla.rmarginal(ndraw, fit$marginals.fitted.values[[pos[k]]]), 1e-10), 1 - 1e-10)
      rho <- pmin(pmax(INLA::inla.rmarginal(ndraw, rho_m), 1e-8), 1 - 1e-6)
      pb  <- stats::rbeta(ndraw, p * (1 - rho) / rho, (1 - p) * (1 - rho) / rho)
      qb  <- stats::quantile(stats::rbinom(ndraw, N, p)  / N, c(.025, .975), names = FALSE)
      qbb <- stats::quantile(stats::rbinom(ndraw, N, pb) / N, c(.025, .975), names = FALSE)
      res[[length(res) + 1]] <- data.frame(design, fold = fn, deno = N, observed = obs,
        lo_binom = qb[1], hi_binom = qb[2], lo_bb = qbb[1], hi_bb = qbb[2])
    }
    message(design, " ", fn, " done")
  }
}
cell <- bind_rows(res)
summ <- cell |>
  mutate(stratum = cut(deno, c(-Inf, 50, 100, Inf), right = FALSE,
                       labels = c("N<50", "50<=N<100", "N>=100"))) |>
  (\(d) bind_rows(mutate(d, stratum = "All"), mutate(d, stratum = as.character(stratum))))() |>
  group_by(design, stratum) |>
  summarise(n = n(),
            cover_binomial = mean(observed >= lo_binom & observed <= hi_binom),
            width_binomial_pp = 100 * mean(hi_binom - lo_binom),
            cover_betabinomial = mean(observed >= lo_bb & observed <= hi_bb),
            width_betabinomial_pp = 100 * mean(hi_bb - lo_bb), .groups = "drop")
readr::write_csv(summ, file.path(out, "pi_betabinomial_M2.csv"))
print(as.data.frame(summ), digits = 3)
