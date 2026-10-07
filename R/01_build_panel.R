## =====================================================================
## 01_build_panel.R -- person-level source files -> county-year analytic panel
##   in : the HIV case file and the hospital discharge (UB) file named in cfg
##        (restricted SAS files, or the synthetic toy CSVs from run_toy.R),
##        the ICD code workbook named in cfg$file_icd
##   out: output/panel_<phenotype>_<denominator>.csv   (46 x 19 = 874 rows)
##        output/person_year_<...>.rds
##        output/tableS1_codelist_<phenotype>.csv, output/codelist_diff.csv
##        output/cohort_flow.csv, output/run_manifest.csv
##        output/panel_change_vs_previous.csv
##        logs/unlisted_substance_codes.csv   <- what to add when the list is revised
## =====================================================================
if (!exists("cfg", envir = globalenv())) source(file.path("R", "00_config.R"))
cfg <- get("cfg", envir = globalenv())
source(file.path(cfg$dir_code, "codelist.R"))

RUN_PREFIX_DIAGNOSTIC <- TRUE   # extra pass reporting what prefix matching would add

## =====================================================================
## PART 1 -- code list
## =====================================================================
cl <- load_code_list()
audit_code_list(cl)

## =====================================================================
## PART 2 -- cohort from the DHEC CASES file
## =====================================================================
flow <- tibble(step = character(), n = integer())
add_flow <- function(step, n) {
  flow <<- add_row(flow, step = step, n = as.integer(n)); log_step(step, ": n=", n)
}

cases <- read_source(cfg$file_cases) |> haven::zap_labels() |> lower_names()
log_step("cases file read: ", nrow(cases), " rows, ", ncol(cases), " cols")

county_col <- tolower(cfg$county_var)
need <- c("rfa_id", "age_at_hiv_dx", "date_of_hiv_dx", "date_of_aids_dx",
          "death_hars", county_col)
missing_cols <- setdiff(need, names(cases))
if (length(missing_cols))
  stop("CASES is missing expected columns: ", paste(missing_cols, collapse = ", "),
       "\nCheck the column names of the case file and update this block.")

cases <- distinct(cases, rfa_id, .keep_all = TRUE)
add_flow("unique RFA_ID in CASES", nrow(cases))

hiv  <- parse_yyyymm(cases$date_of_hiv_dx)
aids <- parse_yyyymm(cases$date_of_aids_dx)
dth  <- parse_yyyymm(cases$death_hars)

coh <- cases |>
  transmute(
    rfa_id = as.character(rfa_id),
    age_at_hiv_dx = as.numeric(age_at_hiv_dx),
    sex  = if ("sex_hars"  %in% names(cases)) as.character(sex_hars)  else NA_character_,
    race = if ("race_hars" %in% names(cases)) as.character(race_hars) else NA_character_,
    risk = if ("risk_hars" %in% names(cases)) as.character(risk_hars) else NA_character_,
    county_raw = .data[[county_col]],
    hiv_y = hiv$year,  hiv_m = hiv$month,
    aids_y = aids$year, aids_m = aids$month,
    death_y = dth$year, death_m = dth$month
  ) |>
  mutate(
    hiv_idx  = ifelse(is.na(hiv_y),  NA_real_, hiv_y  * 100 + hiv_m),
    aids_idx = ifelse(is.na(aids_y), NA_real_, aids_y * 100 + aids_m),
    dx_idx   = suppressWarnings(pmin(hiv_idx, aids_idx, na.rm = TRUE)),
    dx_idx   = ifelse(is.finite(dx_idx), dx_idx, NA_real_),
    dx_year  = as.integer(dx_idx %/% 100),
    death_idx = ifelse(is.na(death_y), NA_real_, death_y * 100 + death_m)
  )

coh <- filter(coh, !is.na(age_at_hiv_dx), age_at_hiv_dx >= 18)
add_flow("age at HIV dx >= 18 and non-missing", nrow(coh))

coh <- filter(coh, !is.na(dx_year))
add_flow("resolvable HIV/AIDS diagnosis date", nrow(coh))

coh <- coh |>
  mutate(county_key = clean_county(county_raw),
         county_key = ifelse(county_key == "MC CORMICK", "MCCORMICK", county_key)) |>
  left_join(sc_counties, by = "county_key")

unmatched <- coh |> filter(is.na(fips)) |> count(county_key, sort = TRUE)
if (nrow(unmatched)) {
  readr::write_csv(unmatched, file.path(cfg$dir_log, "unmatched_county_values.csv"))
  log_step("NOTE: ", nrow(unmatched), " distinct county strings did not match the SC ",
           "crosswalk (", sum(unmatched$n), " persons). See ",
           "logs/unmatched_county_values.csv.")
}
coh <- filter(coh, !is.na(fips))
add_flow("county maps to one of the 46 SC counties", nrow(coh))

coh <- filter(coh, is.na(death_idx) | death_idx > dx_idx)
add_flow("not deceased on or before index diagnosis", nrow(coh))

readr::write_csv(flow, file.path(cfg$dir_out, "cohort_flow.csv"))

## =====================================================================
## PART 3 -- diagnoses from the UB file
## =====================================================================
## TRAP: SDIAG10 is the *10th ICD-9* secondary diagnosis, while SDIAG10_1 ...
## SDIAG10_14 are the *ICD-10* secondaries. Never glob on "sdiag10".
icd9_fields  <- c("pdiag", "adm_diag", paste0("sdiag", 1:14))
icd10_fields <- c("adm_diag10", "pdiag10", paste0("sdiag10_", 1:14))
meta_fields  <- c("rfa_id", "disyear", "dismth", "time_dxdate_disd")

## Source column names may be upper case; col_select is case sensitive, so map
## the lower-case names we reason with back to the file's actual names.
ub_head  <- read_source(cfg$file_ub, n_max = 1)
actual   <- names(ub_head); names(actual) <- tolower(actual)
have9    <- intersect(icd9_fields,  names(actual))
have10   <- intersect(icd10_fields, names(actual))
have_meta<- intersect(meta_fields,  names(actual))
if (!all(c("rfa_id", "disyear") %in% have_meta))
  stop("UB file must contain RFA_ID and DISYEAR. Found: ",
       paste(names(actual), collapse = ", "))
if (!"dismth" %in% have_meta)
  log_step("WARNING: DISMTH absent -- the Oct-2015 ICD split falls back to whole-year 2015.")
log_step("UB diagnosis fields found: ", length(have9), " ICD-9 slots (",
         paste(have9, collapse = ","), "); ", length(have10), " ICD-10 slots (",
         paste(have10, collapse = ","), ")")
missing_expected <- setdiff(c(icd9_fields, icd10_fields), c(have9, have10))
if (length(missing_expected))
  log_step("NOTE: expected but absent diagnosis fields: ",
           paste(missing_expected, collapse = ", "),
           " (absent slots are skipped).")

sel <- unname(actual[c(have_meta, have9, have10)])
ub  <- read_source(cfg$file_ub, col_select = sel) |>
  haven::zap_labels() |> haven::zap_formats() |> lower_names()

ub$rfa_id  <- as.character(ub$rfa_id)
ub$disyear <- suppressWarnings(as.integer(ub$disyear))
ub$dismth  <- if ("dismth" %in% names(ub)) suppressWarnings(as.integer(ub$dismth)) else NA_integer_
log_step("UB rows read: ", nrow(ub), " | discharge years ",
         min(ub$disyear, na.rm = TRUE), "-", max(ub$disyear, na.rm = TRUE))

## *** CRITICAL DATA CHECK ***
yr_cov <- ub |> filter(!is.na(disyear)) |> count(disyear) |> arrange(disyear)
readr::write_csv(yr_cov, file.path(cfg$dir_log, "ub_discharges_per_year.csv"))
if (max(ub$disyear, na.rm = TRUE) < max(cfg$years))
  warning("UB discharges end in ", max(ub$disyear, na.rm = TRUE),
          " but cfg$years runs to ", max(cfg$years),
          ". Numerators for the final years will be structurally incomplete. ",
          "Truncate cfg$years or obtain a newer extract.")

ub <- filter(ub, !is.na(disyear), disyear %in% cfg$years, rfa_id %in% coh$rfa_id)
if (isTRUE(cfg$require_post_dx) && "time_dxdate_disd" %in% names(ub))
  ub <- filter(ub, !is.na(time_dxdate_disd), time_dxdate_disd >= 0)
log_step("UB rows retained for cohort members in study years: ", nrow(ub))

## person-year with any encounter = the encounter-based denominator
enc_py <- distinct(ub, rfa_id, year = disyear)

## ---- scan the diagnosis fields one column at a time -------------------
## A column loop rather than pivot_longer: 365k rows x 30 slots would be ~11M
## rows in memory, and nothing downstream needs the long table itself.
era10 <- ub$disyear > 2015 | (ub$disyear == 2015 & !is.na(ub$dismth) & ub$dismth >= 10)

## families we expect a substance-use code list to cover, used to report codes
## that appear in the data but are NOT on the current list
UNLISTED_RX10 <- "^F1[0-9]"
UNLISTED_RX9  <- "^(29[12]|30[345])"

hits <- vector("list", length(c(have9, have10)))
unlisted <- list(); era_rows <- list()
n_codes_seen <- 0L; n_prefix_extra <- 0L
prefix_rx10 <- if (RUN_PREFIX_DIAGNOSTIC && length(cl$codes10_exact))
  paste0("^(", paste(cl$codes10_exact, collapse = "|"), ")") else NA_character_

for (k in seq_along(c(have9, have10))) {
  v        <- c(have9, have10)[k]
  slot_sys <- if (v %in% have10) "icd10" else "icd9"
  cn       <- normalize_icd(ub[[v]])
  ok       <- !is.na(cn)
  if (!any(ok)) next
  n_codes_seen <- n_codes_seen + sum(ok)

  h10 <- match_sud(cn, cl, "icd10")
  h9  <- match_sud(cn, cl, "icd9")

  is_sud <- if (isTRUE(cfg$enforce_era)) {
    (h10 & (slot_sys == "icd10" | era10)) | (h9 & slot_sys == "icd9" & !era10)
  } else h9 | h10

  if (any(is_sud, na.rm = TRUE))
    hits[[k]] <- tibble(rfa_id = ub$rfa_id[which(is_sud)],
                        year   = ub$disyear[which(is_sud)])

  ## era-discordance audit: a code that matched the other era's code set
  disc <- ok & ((h9 & era10) | (h10 & !era10 & slot_sys == "icd9"))
  if (any(disc, na.rm = TRUE))
    era_rows[[length(era_rows) + 1L]] <-
      tibble(field = v, disyear = ub$disyear[which(disc)], code = cn[which(disc)]) |>
      count(field, disyear, code, name = "n")

  ## substance-family codes present in the data but not counted
  fam <- ok & (( (slot_sys == "icd10" | era10) & stringr::str_detect(cn, UNLISTED_RX10)) |
               ( (slot_sys == "icd9"  & !era10) & stringr::str_detect(cn, UNLISTED_RX9)))
  miss <- fam & !is_sud
  if (any(miss, na.rm = TRUE))
    unlisted[[length(unlisted) + 1L]] <-
      tibble(code = cn[which(miss)],
             era  = ifelse(era10[which(miss)], "icd10", "icd9")) |>
      count(code, era, name = "n")

  if (!is.na(prefix_rx10)) {
    pf <- ok & stringr::str_detect(cn, prefix_rx10)
    n_prefix_extra <- n_prefix_extra + sum(pf & !is_sud, na.rm = TRUE)
  }
}

sud_py <- bind_rows(hits) |> distinct(rfa_id, year) |> mutate(sud = 1L)
log_step("diagnosis entries scanned: ", n_codes_seen,
         " | person-years with >=1 qualifying SUD diagnosis: ", nrow(sud_py))

if (length(era_rows)) {
  readr::write_csv(bind_rows(era_rows) |>
                     count(field, disyear, code, wt = n, name = "n_entries"),
                   file.path(cfg$dir_log, "era_discordance_check.csv"))
  log_step("era-discordant matches written to logs/era_discordance_check.csv")
}

## *** the report to read when the code list is about to be revised ***
if (length(unlisted)) {
  unl <- bind_rows(unlisted) |> count(code, era, wt = n, name = "n_entries") |>
    arrange(desc(n_entries))
  readr::write_csv(unl, file.path(cfg$dir_log, "unlisted_substance_codes.csv"))
  log_step("substance-family codes present in the data but NOT counted by the ",
           "current list: ", nrow(unl), " distinct codes, ", sum(unl$n_entries),
           " diagnosis entries. Top: ",
           paste(head(unl$code, 12), collapse = ", "),
           " -- see logs/unlisted_substance_codes.csv")
}
if (RUN_PREFIX_DIAGNOSTIC && !is.na(prefix_rx10))
  log_step("prefix-vs-exact diagnostic: prefix matching would add ",
           n_prefix_extra, " ICD-10 diagnosis entries")

rm(ub); invisible(gc())

## =====================================================================
## PART 4 -- person-year file and county-year panel
## =====================================================================
py <- coh |>
  select(rfa_id, county, fips, dx_year, death_y) |>
  tidyr::crossing(year = cfg$years) |>
  filter(dx_year <= year, is.na(death_y) | death_y >= year) |>
  left_join(sud_py, by = c("rfa_id", "year")) |>
  mutate(sud = tidyr::replace_na(sud, 0L)) |>
  left_join(mutate(enc_py, any_encounter = 1L), by = c("rfa_id", "year")) |>
  mutate(any_encounter = tidyr::replace_na(any_encounter, 0L))

check_that(all(py$sud <= py$any_encounter),
           "every SUD person-year also has an encounter that year")

if (cfg$denominator == "encounter") py <- filter(py, any_encounter == 1L)

saveRDS(py, file.path(cfg$dir_out,
        sprintf("person_year_%s_%s.rds", cfg$phenotype, cfg$denominator)))

panel <- py |>
  group_by(fips, county, year) |>
  summarise(deno = n(), nume = sum(sud), .groups = "drop") |>
  right_join(tidyr::crossing(sc_counties |> select(county, fips), year = cfg$years),
             by = c("county", "fips", "year")) |>
  mutate(deno = tidyr::replace_na(deno, 0L), nume = tidyr::replace_na(nume, 0L),
         encounter_rate = ifelse(deno > 0, nume / deno, NA_real_),
         encounter_pct  = 100 * encounter_rate) |>
  arrange(county, year)

check_that(nrow(panel) == 46 * length(cfg$years), "panel is 46 counties x 19 years")
check_that(all(panel$nume <= panel$deno), "numerator never exceeds denominator")
check_that(n_distinct(panel$fips) == 46, "46 distinct FIPS codes")
if (any(panel$deno == 0))
  log_step("WARNING: ", sum(panel$deno == 0), " county-years have a zero denominator")

out_panel <- file.path(cfg$dir_out,
             sprintf("panel_%s_%s.csv", cfg$phenotype, cfg$denominator))

## ---- what changed since the last run of this same configuration -------
if (file.exists(out_panel)) {
  prev <- readr::read_csv(out_panel, show_col_types = FALSE)
  cmp <- prev |> group_by(year) |>
      summarise(prev_deno = sum(deno), prev_nume = sum(nume), .groups = "drop") |>
    full_join(panel |> group_by(year) |>
      summarise(new_deno = sum(deno), new_nume = sum(nume), .groups = "drop"),
      by = "year") |>
    mutate(rate_pct_old = 100 * prev_nume / prev_deno,
           rate_pct_new = 100 * new_nume / new_deno,
           change_pp    = rate_pct_new - rate_pct_old)
  readr::write_csv(cmp, file.path(cfg$dir_out, "panel_change_vs_previous.csv"))
  log_step("panel differs from the previous run by a mean of ",
           round(mean(abs(cmp$change_pp), na.rm = TRUE), 3),
           " encounter-rate pp per year -- see output/panel_change_vs_previous.csv")
  file.copy(out_panel, file.path(cfg$dir_out,
            sprintf("panel_%s_%s_PREVIOUS.csv", cfg$phenotype, cfg$denominator)),
            overwrite = TRUE)
}

readr::write_csv(panel, out_panel)
log_step("panel written: ", out_panel)

statewide <- panel |> group_by(year) |>
  summarise(deno = sum(deno), nume = sum(nume),
            wpate_pct = round(100 * nume / deno, 2), .groups = "drop")
print(as.data.frame(statewide))

writeLines(capture.output(sessionInfo()),
           file.path(cfg$dir_log, "sessionInfo_01.txt"))
