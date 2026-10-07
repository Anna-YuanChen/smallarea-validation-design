## =====================================================================
## simulate_toy_data.R -- a fully synthetic stand-in for the restricted files
##
## The analysis used two linked person-level files that cannot be shared: an
## HIV surveillance case file (one row per person) and an all-payer hospital
## discharge file (one row per discharge, with separate ICD-9-CM and ICD-10-CM
## diagnosis slots). simulate_toy_data() writes two CSV files with the SAME
## column names and formats, so the unchanged pipeline (01_ to 04_) runs on
## them end to end. Nothing in these files is derived from the real data:
## county sizes, trends and effects are drawn from the parameters below.
##
## The data-generating process is a known "truth" for every county-year,
## written to truth_county_year.csv, so users can check how close each
## model's estimates come to it -- something the real data cannot offer.
##
##   logit p_it = alpha + b_i + gamma_t + delta_it
##     b_i      county effect: sd_county * (sqrt(1-phi) v_i + sqrt(phi) u_i),
##              u_i a scaled ICAR field on the county adjacency (if supplied)
##     gamma_t  smooth year trend + a step at the Oct-2015 ICD-10 transition
##     delta_it small county-year noise
##   person-year SUD-related encounter ~ Bernoulli(p_it)
## =====================================================================

.icar_field <- function(nb) {
  n <- length(nb)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) W[i, nb[[i]]] <- 1
  Q <- diag(rowSums(W)) - W
  e <- eigen(Q, symmetric = TRUE)
  keep <- e$values > 1e-8                       # drop the constant direction
  z <- e$vectors[, keep] %*% (stats::rnorm(sum(keep)) / sqrt(e$values[keep]))
  z <- as.numeric(z - mean(z))
  z / stats::sd(z)
}

simulate_toy_data <- function(out_dir        = "toy_data",
                              n_persons      = 30000,
                              years          = 2005:2023,
                              sud_codes9,                 # listed ICD-9-CM codes (dot-free)
                              sud_codes10,                # listed ICD-10-CM codes (dot-free)
                              nb             = NULL,      # spdep nb in sc_counties FIPS order
                              alpha          = -3.6,
                              sd_county      = 0.30,
                              phi            = 0.5,
                              trend_slope    = -0.03,
                              icd10_step     = -0.35,
                              sd_year        = 0.05,
                              sd_interaction = 0.04,
                              p_encounter    = 0.28,
                              seed           = 20261007) {
  stopifnot(length(sud_codes9) > 0, length(sud_codes10) > 0)
  set.seed(seed)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  K  <- nrow(sc_counties)
  ty <- years - min(years)

  ## ---- county sizes: synthetic, skewed like most US states --------------
  w <- stats::rlnorm(K, meanlog = 0, sdlog = 1.0)
  w <- pmax(w / sum(w), 0.002); w <- w / sum(w)   # every county gets some persons

  ## ---- true county-year proportions ------------------------------------
  v <- stats::rnorm(K)
  u <- if (!is.null(nb)) .icar_field(nb) else stats::rnorm(K)
  b <- sd_county * (sqrt(1 - phi) * as.numeric(scale(v)) + sqrt(phi) * u)
  gamma <- trend_slope * ty + icd10_step * ifelse(years >= 2016, 1, ifelse(years == 2015, 0.25, 0)) +
    stats::rnorm(length(years), 0, sd_year)
  truth <- tidyr::crossing(k = seq_len(K), t = seq_along(years)) |>
    dplyr::mutate(fips = sc_counties$fips[k], county = sc_counties$county[k],
                  year = years[t],
                  eta = alpha + b[k] + gamma[t] + stats::rnorm(dplyr::n(), 0, sd_interaction),
                  p_true = stats::plogis(eta))

  ## ---- persons ---------------------------------------------------------
  n <- n_persons
  dx_year <- sample(1985:max(years), n, replace = TRUE,
                    prob = stats::plogis((1985:max(years) - 1995) / 4))
  dx_mon  <- sample(1:12, n, replace = TRUE)
  k_dx    <- sample(seq_len(K), n, replace = TRUE, prob = w)
  moved   <- stats::runif(n) < 0.15
  k_cur   <- ifelse(moved, sample(seq_len(K), n, replace = TRUE, prob = w), k_dx)
  age     <- round(stats::rgamma(n, shape = 9, rate = 0.27), 0)    # mean ~33
  age[stats::runif(n) < 0.01] <- NA                               # missing age
  u18 <- stats::runif(n) < 0.02
  age[u18] <- sample(13:17, sum(u18), replace = TRUE)             # under 18, excluded
  ## annual mortality ~1.5%, so most diagnosed persons remain in the cohort
  yrs_to_death <- stats::rgeom(n, 0.015) + 1
  death_year   <- ifelse(dx_year + yrs_to_death <= max(years), dx_year + yrs_to_death, NA)
  death_mon    <- ifelse(is.na(death_year), NA, sample(1:12, n, replace = TRUE))
  has_aids     <- stats::runif(n) < 0.35
  aids_lag     <- stats::rgeom(n, 0.25)
  aids_year    <- ifelse(has_aids, pmin(dx_year + aids_lag, max(years)), NA)
  aids_mon     <- ifelse(has_aids, sample(1:12, n, replace = TRUE), NA)

  ymd <- function(y, m, unknown_share = 0) {
    mm <- sprintf("%02d", m)
    mm[stats::runif(length(y)) < unknown_share] <- ".."           # unknown month
    ifelse(is.na(y), NA, paste0(y, mm))
  }
  cty_string <- function(k) {
    s <- paste0(toupper(sc_counties$county[k]), " CO.")
    s[stats::runif(length(k)) < 0.004] <- NA                        # missing county
    s
  }
  cases <- tibble::tibble(
    RFA_ID           = sprintf("T%06d", seq_len(n)),
    AGE_AT_HIV_DX    = age,
    DATE_OF_HIV_DX   = ymd(dx_year, dx_mon, unknown_share = 0.04),
    DATE_OF_AIDS_DX  = ymd(aids_year, aids_mon, unknown_share = 0.04),
    DEATH_HARS       = ymd(death_year, death_mon),
    COUNTY_AT_HIV_DX = cty_string(k_dx),
    CURRENT_COUNTY   = cty_string(k_cur),
    SEX_HARS         = sample(c("M", "F"), n, TRUE, c(.72, .28)),
    RACE_HARS        = sample(c("Black", "White", "Hisp", "Other"), n, TRUE, c(.68, .24, .05, .03)),
    RISK_HARS        = sample(c("01) MSM", "05) Heterosexual", "02) Injecting Drug Use",
                                "03) MSM/Injecting Drug Use", "08) Other"), n, TRUE,
                              c(.55, .30, .07, .03, .05))
  )

  ## ---- person-years in the study window --------------------------------
  first <- pmax(dx_year, min(years)); last <- pmin(ifelse(is.na(death_year), max(years), death_year), max(years))
  ok <- first <= last
  idx <- rep(which(ok), times = (last - first + 1)[ok])
  yr  <- unlist(lapply(which(ok), function(i) first[i]:last[i]))
  py  <- tibble::tibble(i = idx, year = yr, k = k_dx[idx]) |>
    dplyr::left_join(dplyr::select(truth, k, year, p_true), by = c("k", "year"))
  py$sud <- stats::rbinom(nrow(py), 1, py$p_true)
  py$enc <- ifelse(py$sud == 1, 1L, stats::rbinom(nrow(py), 1, p_encounter))
  py <- dplyr::filter(py, enc == 1L)

  ## ---- discharges ------------------------------------------------------
  ndis <- 1L + stats::rpois(nrow(py), 0.5)
  d <- py[rep(seq_len(nrow(py)), ndis), ]
  d$first_of_py <- !duplicated(paste(d$i, d$year))
  d$month <- sample(1:12, nrow(d), replace = TRUE)
  d$era10 <- d$year > 2015 | (d$year == 2015 & d$month >= 10)
  d$disch_id <- seq_len(nrow(d))
  d$time <- round(((d$year - dx_year[d$i]) * 12 + (d$month - dx_mon[d$i])) * 30.4)

  bg9  <- c("042", "V08", "4019", "25000", "2724", "311", "53081", "78650", "49390", "27800")
  bg10 <- c("B20", "Z21", "I10", "E119", "E785", "F329", "K219", "R0789", "J45909", "Z7984")
  tob9 <- c("3051", "V1582");  tob10 <- c("F17210", "Z87891")     # tobacco: never counted
  alc9 <- c("30500", "30390"); alc10 <- c("F1010", "F1020")       # alcohol: excluded by definition
  dep9 <- c("30400", "30420", "30430"); dep10 <- c("F1120", "F1420", "F1290")  # not on the list

  codes_for <- function(era10, sud) {
    m  <- sample(2:8, 1)
    cs <- sample(if (era10) bg10 else bg9, m, replace = TRUE)
    if (stats::runif(1) < 0.30) cs <- c(cs, sample(if (era10) tob10 else tob9, 1))
    if (stats::runif(1) < 0.05) cs <- c(cs, sample(if (era10) alc10 else alc9, 1))
    if (stats::runif(1) < 0.02) cs <- c(cs, sample(if (era10) dep10 else dep9, 1))
    if (sud) cs <- append(cs, sample(if (era10) sud_codes10 else sud_codes9, 1),
                          after = sample(0:length(cs), 1))
    unique(cs)[seq_len(min(length(unique(cs)), 16))]
  }
  diag_list <- lapply(seq_len(nrow(d)), function(j)
    codes_for(d$era10[j], d$sud[j] == 1L && d$first_of_py[j]))

  ## Layout of the source file: 16 ICD-9-CM slots and 15 ICD-10-CM slots. There
  ## is no separate ICD-10 principal-diagnosis field; after October 2015 the
  ## principal diagnosis is carried in PDIAG, which is why 01_ matches codes by
  ## discharge date (era) rather than by slot.
  slots9  <- c("PDIAG", "ADM_DIAG", paste0("SDIAG", 1:14))
  slots10 <- c("ADM_DIAG10", paste0("SDIAG10_", 1:14))
  M9  <- matrix(NA_character_, nrow(d), length(slots9),  dimnames = list(NULL, slots9))
  M10 <- matrix(NA_character_, nrow(d), length(slots10), dimnames = list(NULL, slots10))
  for (j in seq_len(nrow(d))) {
    cs <- diag_list[[j]]
    if (d$era10[j]) {
      M9[j, "PDIAG"] <- cs[1]
      rest <- cs[-1][seq_len(min(length(cs) - 1, length(slots10)))]
      if (length(rest)) M10[j, seq_along(rest)] <- rest
    } else M9[j, seq_along(cs)] <- cs
  }
  ub <- dplyr::bind_cols(
    tibble::tibble(RFA_ID = cases$RFA_ID[d$i], DISYEAR = d$year, DISMTH = d$month,
                   TIME_DXDATE_DISD = d$time),
    tibble::as_tibble(M9), tibble::as_tibble(M10))
  ub <- ub[sample(nrow(ub)), ]                                     # no implied ordering

  ## ---- write ------------------------------------------------------------
  readr::write_csv(cases, file.path(out_dir, "toy_cases.csv.gz"), na = "")
  readr::write_csv(ub,    file.path(out_dir, "toy_ub_discharges.csv.gz"), na = "")
  readr::write_csv(dplyr::select(truth, fips, county, year, p_true),
                   file.path(out_dir, "truth_county_year.csv"))
  writeLines(c(
    "Synthetic data generated by simulate_toy_data(); not derived from any real records.",
    paste("seed:", seed), paste("persons:", n_persons),
    paste("discharge rows:", nrow(ub)),
    paste("alpha:", alpha, "| sd_county:", sd_county, "| phi:", phi,
          "| trend_slope:", trend_slope, "| icd10_step:", icd10_step,
          "| sd_year:", sd_year, "| sd_interaction:", sd_interaction,
          "| p_encounter:", p_encounter)),
    file.path(out_dir, "README_toy_data.txt"))
  invisible(list(cases = cases, ub = ub, truth = truth))
}
