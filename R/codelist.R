## =====================================================================
## codelist.R -- load, normalize, expand and audit the ICD code list
## =====================================================================
## Everything phenotype-related lives here so that a revised code list is a
## one-line change in 00_config.R. The loader tolerates:
##   * renamed columns          (auto-detection, with manual override)
##   * wide or long layouts     (icd9/icd10 columns, or code + code_system)
##   * wildcard entries         "F11.*", "F11%", "304.x"
##   * range entries            "F11.10-F11.19"
##   * an include/exclude flag column
##   * codes written with or without dots, padded, or lower case
## and it hashes the file and diffs it against the list used last time.
## =====================================================================

## ---- column detection -------------------------------------------------
.detect <- function(nms, patterns, override = NULL) {
  if (!is.null(override)) {
    if (!override %in% nms) stop("cfg$icd_cols names a column not in the sheet: ", override)
    return(override)
  }
  for (p in patterns) {
    hit <- grep(p, nms, value = TRUE, perl = TRUE)
    if (length(hit)) return(hit[1])
  }
  NULL
}

## ---- wildcard / range expansion ---------------------------------------
## Returns a tibble(code, is_prefix). A range is expanded into its members
## when the two endpoints share a prefix and differ only in a numeric tail;
## otherwise both endpoints are kept and a warning is logged.
.expand_entry <- function(x) {
  if (is.na(x) || x == "") return(tibble(code = character(0), is_prefix = logical(0)))

  ## wildcard -> prefix
  if (grepl("[%*]", x)) {
    p <- sub("[%*]+$", "", gsub("[%*]", "", x))
    if (p == "") return(tibble(code = character(0), is_prefix = logical(0)))
    return(tibble(code = p, is_prefix = TRUE))
  }

  ## range
  if (grepl("-", x)) {
    parts <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(parts) == 2L && nchar(parts[1]) == nchar(parts[2])) {
      a <- parts[1]; b <- parts[2]
      same <- which(strsplit(a, "")[[1]] == strsplit(b, "")[[1]])
      k <- if (length(same) && all(seq_along(same) == same)) max(same) else 0L
      ta <- substring(a, k + 1); tb <- substring(b, k + 1)
      if (k > 0 && grepl("^[0-9]+$", ta) && grepl("^[0-9]+$", tb) &&
          as.integer(tb) >= as.integer(ta)) {
        tail_seq <- formatC(as.integer(ta):as.integer(tb), width = nchar(ta),
                            flag = "0", format = "d")
        return(tibble(code = paste0(substring(a, 1, k), tail_seq), is_prefix = FALSE))
      }
    }
    log_step("NOTE: could not expand range entry '", x, "' -- kept endpoints only")
    return(tibble(code = parts, is_prefix = FALSE))
  }

  ## trailing X used as a wildcard (e.g. "F112X", "304XX") but NOT a real code
  if (grepl("^[A-Z]?[0-9]+X+$", x)) {
    return(tibble(code = sub("X+$", "", x), is_prefix = TRUE))
  }

  tibble(code = x, is_prefix = FALSE)
}

## ---- regex families the current workbook omits ------------------------
## Applied for phenotype "revised_core" and "broad". Tested against the
## NORMALIZED (dot-free) code, so "F1120" and "30400" are the strings matched.
##   ICD-10: F10-F19 except F17, 4th char in {1 abuse, 2 dependence, 9 use}
##   ICD-9 : 291.x, 292.x, 303.xx, 304.xx, 305.0x, 305.2x-305.9x (not 305.1)
RX10_DISORDER <- "^F1[12345689][129]"        # F11-F16, F18, F19 (F10 alcohol and F17 tobacco excluded)
RX9_DISORDER  <- "^(292|304|305[2-9])"       # 292 drug-induced mental disorders, 304 dependence, 305.2-305.9 abuse
RX10_ORGAN    <- NULL                        # organ sequelae are all alcohol-attributable; not used
RX9_ORGAN     <- NULL

## =====================================================================
## load_code_list()
## =====================================================================
load_code_list <- function() {

  sheet <- cfg$icd_sheet
  if (is.null(sheet)) sheet <- readxl::excel_sheets(cfg$file_icd)[1]
  raw <- readxl::read_excel(cfg$file_icd, sheet = sheet) |> rename_with(tolower)
  raw <- raw |> mutate(across(where(is.character), trimws))
  nms <- names(raw)
  log_step("code list: '", basename(cfg$file_icd), "' sheet '", sheet, "' -- ",
           nrow(raw), " rows, columns: ", paste(nms, collapse = ", "))

  co <- cfg$icd_cols
  c_icd9  <- .detect(nms, c("icd.?9", "dx.?9", "code.?9"),  co$icd9)
  c_icd10 <- .detect(nms, c("icd.?10", "dx.?10", "code.?10"), co$icd10)
  c_code  <- .detect(nms, c("^code$", "^icd$", "^icd_code$", "^dx_code$", "^diag"), co$code)
  c_sys   <- .detect(nms, c("system", "version", "icd_?ver", "coding"), co$system)
  c_cat   <- .detect(nms, c("^name$", "^category$", "^group$", "^class$",
                            "condition", "disease_?group"), co$category)
  c_sub   <- .detect(nms, c("sub_?class", "sub_?name", "sub_?category",
                            "sub_?group"), co$subcategory)
  c_lab   <- .detect(nms, c("specific_disease_name", "description", "^label$",
                            "^desc"), co$label)
  c_inc   <- .detect(nms, c("^include$", "^exclude$", "^use$", "^keep$",
                            "^flag$", "^final$"), co$include)

  ## ---- reshape to one row per (code, system) --------------------------
  if (!is.null(c_icd9) || !is.null(c_icd10)) {
    layout <- "wide"
    parts <- list()
    if (!is.null(c_icd9))
      parts[[1]] <- tibble(code_raw = raw[[c_icd9]],  system = "icd9",  row = seq_len(nrow(raw)))
    if (!is.null(c_icd10))
      parts[[2]] <- tibble(code_raw = raw[[c_icd10]], system = "icd10", row = seq_len(nrow(raw)))
    long <- bind_rows(parts)
  } else if (!is.null(c_code)) {
    layout <- "long"
    sysv <- if (!is.null(c_sys)) tolower(as.character(raw[[c_sys]])) else NA_character_
    long <- tibble(code_raw = raw[[c_code]], row = seq_len(nrow(raw)),
                   system = dplyr::case_when(
                     grepl("10", sysv)          ~ "icd10",
                     grepl("9",  sysv)          ~ "icd9",
                     ## no system column: infer from the code itself
                     grepl("^[A-Z]", toupper(trimws(as.character(raw[[c_code]])))) &
                       !grepl("^[EV]", toupper(trimws(as.character(raw[[c_code]])))) ~ "icd10",
                     TRUE                        ~ "icd9"))
    if (is.null(c_sys))
      log_step("NOTE: no code-system column found; system inferred from code format. ",
               "Add an explicit column to the workbook if ICD-9 E/V codes are used.")
  } else {
    stop("Could not find any ICD code column in sheet '", sheet, "'.\n",
         "Columns present: ", paste(nms, collapse = ", "), "\n",
         "Set cfg$icd_cols$icd9 / $icd10 (or $code and $system) in 00_config.R.")
  }

  meta <- tibble(
    row         = seq_len(nrow(raw)),
    category    = if (!is.null(c_cat)) as.character(raw[[c_cat]]) else NA_character_,
    subcategory = if (!is.null(c_sub)) as.character(raw[[c_sub]]) else NA_character_,
    label       = if (!is.null(c_lab)) as.character(raw[[c_lab]]) else NA_character_,
    include_raw = if (!is.null(c_inc)) as.character(raw[[c_inc]]) else NA_character_
  )
  long <- left_join(long, meta, by = "row")

  ## ---- normalize and expand -------------------------------------------
  long <- long |>
    mutate(code_norm = normalize_icd(fix_excel_float(code_raw))) |>
    filter(!is.na(code_norm))

  expanded <- purrr::map2_dfr(long$code_norm, seq_len(nrow(long)), function(cd, i) {
    e <- .expand_entry(cd)
    if (!nrow(e)) return(NULL)
    bind_cols(e, long[rep(i, nrow(e)), c("system", "category", "subcategory",
                                         "label", "include_raw", "code_raw")])
  })
  if (is.null(expanded) || nrow(expanded) == 0)
    stop("No usable ICD codes were read from sheet '", sheet, "'. ",
         "Check cfg$icd_sheet and cfg$icd_cols.")
  n_expanded <- nrow(expanded) - nrow(long)
  if (n_expanded > 0)
    log_step("wildcard/range expansion added ", n_expanded, " code entries")

  ## ---- include/exclude flag column ------------------------------------
  if (!is.null(c_inc)) {
    neg <- grepl("^(n|no|0|f|false|exclude|drop)$", tolower(expanded$include_raw))
    if (tolower(c_inc) == "exclude")
      neg <- grepl("^(y|yes|1|t|true|x)$", tolower(expanded$include_raw))
    if (any(neg, na.rm = TRUE)) {
      log_step("include/exclude column '", c_inc, "' removed ",
               sum(neg, na.rm = TRUE), " entries")
      expanded <- expanded[!(neg %in% TRUE), ]
    }
  }

  ## ---- GUARD (a): drop tobacco/nicotine categories ---------------------
  cat_str <- paste(tolower(expanded$category), tolower(expanded$subcategory))
  is_tob_cat <- grepl(cfg$exclude_category_regex, cat_str)
  is_organ   <- grepl(cfg$organ_category_regex, cat_str)
  expanded <- expanded |> mutate(is_tobacco_cat = is_tob_cat, is_organ = is_organ)
  log_step("category flags: ", sum(is_tob_cat), " tobacco/alcohol rows, ",
           sum(is_organ), " alcohol-attributable organ-disease rows")
  if (sum(is_tob_cat) == 0 && is.null(c_cat))
    log_step("WARNING: no category column detected, so the tobacco exclusion ",
             "rests entirely on the code-level blocklist. Verify manually.")

  ## ---- GUARD (b): code-level tobacco blocklist ------------------------
  blk <- normalize_icd(c(cfg$tobacco_blocklist, cfg$alcohol_blocklist))
  blk <- blk[!is.na(blk)]
  blk_rx <- if (length(blk)) paste0("^(", paste(blk, collapse = "|"), ")") else NA_character_
  if (!is.na(blk_rx)) {
    hit <- grepl(blk_rx, expanded$code)
    if (any(hit))
      log_step("tobacco/alcohol block list removed ", sum(hit), " codes")
    expanded <- expanded[!hit, ]
  }

  ## ---- assemble the requested phenotype -------------------------------
  ph <- cfg$phenotype
  base <- switch(ph,
                 "core"         = filter(expanded, !is_tobacco_cat, !is_organ),
                 "listed_all"   = filter(expanded, !is_tobacco_cat),
                 "revised_core" = filter(expanded, !is_tobacco_cat, !is_organ),
                 "broad"        = filter(expanded, !is_tobacco_cat),
                 stop("unknown cfg$phenotype: ", ph))

  rx9 <- character(0); rx10 <- character(0)
  if (ph %in% c("revised_core", "broad")) { rx9 <- RX9_DISORDER; rx10 <- RX10_DISORDER }
  if (ph == "broad") { rx9 <- c(rx9, RX9_ORGAN); rx10 <- c(rx10, RX10_ORGAN) }

  split_codes <- function(sys) {
    d <- filter(base, system == sys)
    list(exact  = sort(unique(d$code[!d$is_prefix])),
         prefix = sort(unique(d$code[ d$is_prefix])))
  }
  s9 <- split_codes("icd9"); s10 <- split_codes("icd10")

  ## In prefix mode every listed code becomes a prefix.
  if (cfg$match_mode == "prefix") {
    s9  <- list(exact = character(0), prefix = sort(unique(c(s9$exact,  s9$prefix))))
    s10 <- list(exact = character(0), prefix = sort(unique(c(s10$exact, s10$prefix))))
  }

  ## Run-level manual overrides (prefix-matched, applied to both systems).
  fin  <- normalize_icd(cfg$codes_force_in);  fin  <- fin[!is.na(fin)]
  fout <- normalize_icd(cfg$codes_force_out); fout <- fout[!is.na(fout)]
  if (length(fin)) {
    s9$prefix  <- sort(unique(c(s9$prefix,  fin)))
    s10$prefix <- sort(unique(c(s10$prefix, fin)))
    log_step("cfg$codes_force_in added prefixes: ", paste(fin, collapse = ", "))
  }

  out <- list(
    codes9_exact = s9$exact,  codes9_prefix = s9$prefix,
    codes10_exact = s10$exact, codes10_prefix = s10$prefix,
    rx9 = rx9, rx10 = rx10,
    blocklist_rx = blk_rx,
    force_out_rx = if (length(fout)) paste0("^(", paste(fout, collapse = "|"), ")") else NA_character_,
    table = base, layout = layout, sheet = sheet,
    hash = tryCatch(tools::md5sum(cfg$file_icd)[[1]], error = function(e) NA_character_)
  )

  ## no code should be claimed by both systems
  overlap <- intersect(c(out$codes9_exact, out$codes9_prefix),
                       c(out$codes10_exact, out$codes10_prefix))
  if (length(overlap))
    log_step("WARNING: ", length(overlap), " codes appear under both ICD-9 and ",
             "ICD-10: ", paste(head(overlap, 10), collapse = ", "),
             ". Era enforcement will decide which era they count in.")

  log_step("phenotype '", ph, "': ",
           length(out$codes9_exact), " exact + ", length(out$codes9_prefix),
           " prefix ICD-9 codes; ",
           length(out$codes10_exact), " exact + ", length(out$codes10_prefix),
           " prefix ICD-10 codes; ",
           length(rx9) + length(rx10), " regex families")
  out
}

## =====================================================================
## audit_code_list() -- snapshot, hash, and diff against the previous run
## =====================================================================
audit_code_list <- function(cl) {
  snap <- cl$table |>
    transmute(code, system, is_prefix, category, subcategory, label, code_raw) |>
    distinct() |> arrange(system, code)

  snap_path <- file.path(cfg$dir_out, "codelist_snapshot.csv")
  prev <- if (file.exists(snap_path))
    readr::read_csv(snap_path, show_col_types = FALSE) else NULL

  if (!is.null(prev)) {
    key_new <- paste(snap$system, snap$code)
    key_old <- paste(prev$system, prev$code)
    added   <- snap[!key_new %in% key_old, ]
    removed <- prev[!key_old %in% key_new, ]
    diff <- bind_rows(mutate(added, change = "added"),
                      mutate(removed, change = "removed"))
    readr::write_csv(diff, file.path(cfg$dir_out, "codelist_diff.csv"))
    if (nrow(diff))
      log_step("CODE LIST CHANGED since the last run: ", nrow(added), " codes added, ",
               nrow(removed), " removed -- see output/codelist_diff.csv")
    else
      log_step("code list unchanged since the last run")
  }
  readr::write_csv(snap, snap_path)
  readr::write_csv(snap, file.path(cfg$dir_out,
                   paste0("tableS1_codelist_", cfg$phenotype, ".csv")))

  manifest <- tibble(
    run_time = format(Sys.time()), phenotype = cfg$phenotype,
    denominator = cfg$denominator, county_var = cfg$county_var,
    match_mode = cfg$match_mode, enforce_era = cfg$enforce_era,
    codelist_file = basename(cfg$file_icd), codelist_sheet = cl$sheet,
    codelist_layout = cl$layout, codelist_md5 = cl$hash,
    n_codes_icd9 = length(cl$codes9_exact) + length(cl$codes9_prefix),
    n_codes_icd10 = length(cl$codes10_exact) + length(cl$codes10_prefix))
  mpath <- file.path(cfg$dir_out, "run_manifest.csv")
  readr::write_csv(manifest, mpath,
                   append = file.exists(mpath), col_names = !file.exists(mpath))
  invisible(snap)
}

## =====================================================================
## match_sud() -- TRUE where a normalized diagnosis code is a qualifying code
## =====================================================================
match_sud <- function(code_norm, cl, system = c("icd9", "icd10")) {
  system <- match.arg(system)
  exact  <- if (system == "icd9") cl$codes9_exact  else cl$codes10_exact
  prefix <- if (system == "icd9") cl$codes9_prefix else cl$codes10_prefix
  rx     <- if (system == "icd9") cl$rx9           else cl$rx10

  hit <- code_norm %in% exact
  if (length(prefix))
    hit <- hit | stringr::str_detect(code_norm, paste0("^(", paste(prefix, collapse = "|"), ")"))
  if (length(rx))
    hit <- hit | stringr::str_detect(code_norm, paste0("(", paste(rx, collapse = "|"), ")"))

  ## Final guards, applied to the DIAGNOSIS code itself. Even if a future code
  ## list or regex family reaches a tobacco code, it can never count here.
  if (!is.na(cl$blocklist_rx))
    hit <- hit & !stringr::str_detect(code_norm, cl$blocklist_rx)
  if (!is.na(cl$force_out_rx))
    hit <- hit & !stringr::str_detect(code_norm, cl$force_out_rx)

  hit & !is.na(code_norm)
}
