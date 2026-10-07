## =====================================================================
## session_report.R -- record the computing environment next to results
##
## Fitting times depend on the CPU, memory, BLAS and on how many threads each
## package actually uses (INLA is multithreaded; sf, spdep and the tidyverse
## steps run on one thread). This writes everything a reader needs to judge
## reported clock times: R and package versions, OS, CPU model, core count,
## physical memory, BLAS/LAPACK and the INLA thread setting.
## =====================================================================

.safe <- function(expr) tryCatch(expr, error = function(e) NA_character_,
                                 warning = function(w) NA_character_)

machine_spec <- function() {
  sys <- Sys.info()[["sysname"]]
  cpu <- NA_character_; mem_gb <- NA_real_; model <- NA_character_
  if (sys == "Darwin") {
    cpu    <- .safe(system("sysctl -n machdep.cpu.brand_string", intern = TRUE))
    model  <- .safe(system("sysctl -n hw.model", intern = TRUE))
    mem_gb <- as.numeric(.safe(system("sysctl -n hw.memsize", intern = TRUE))) / 1024^3
    perf   <- .safe(system("sysctl -n hw.perflevel0.physicalcpu", intern = TRUE))
    eff    <- .safe(system("sysctl -n hw.perflevel1.physicalcpu", intern = TRUE))
    if (!is.na(perf)) cpu <- paste0(cpu, " (", perf, " performance + ",
                                    ifelse(is.na(eff), "?", eff), " efficiency cores)")
  } else if (sys == "Linux") {
    ci  <- .safe(readLines("/proc/cpuinfo"))
    cpu <- .safe(sub(".*:\\s*", "", grep("model name", ci, value = TRUE)[1]))
    mi  <- .safe(readLines("/proc/meminfo"))
    mem_gb <- as.numeric(.safe(gsub("\\D", "", grep("MemTotal", mi, value = TRUE)))) / 1024^2
  } else if (sys == "Windows") {
    cpu <- Sys.getenv("PROCESSOR_IDENTIFIER")
    m <- .safe(system("wmic ComputerSystem get TotalPhysicalMemory /value", intern = TRUE))
    m <- suppressWarnings(as.numeric(gsub("\\D", "", grep("=", m, value = TRUE))))
    if (length(m) && !is.na(m[1])) mem_gb <- m[1] / 1024^3
  }
  list(os = paste(sys, Sys.info()[["release"]]),
       hardware_model = model,
       cpu = cpu,
       logical_cores = parallel::detectCores(logical = TRUE),
       physical_cores = .safe(parallel::detectCores(logical = FALSE)),
       memory_gb = round(mem_gb, 1))
}

write_session_report <- function(path, extra = NULL) {
  ms <- machine_spec()
  si <- utils::sessionInfo()
  key_pkgs <- c("INLA", "fmesher", "Matrix", "sf", "spdep", "spData", "tigris",
                "haven", "readxl", "dplyr", "tidyr", "readr", "purrr", "stringr",
                "tibble", "ggplot2")
  vers <- vapply(key_pkgs, function(p) .safe(as.character(utils::packageVersion(p))),
                 character(1))
  inla_thr <- if (requireNamespace("INLA", quietly = TRUE))
    .safe(as.character(INLA::inla.getOption("num.threads"))) else NA_character_
  lines <- c(
    "## Computing environment",
    paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    paste("R:", R.version.string, "|", si$platform),
    paste("OS:", ms$os, if (!is.null(si$running)) paste0("(", si$running, ")") else ""),
    paste("Hardware model:", ms$hardware_model),
    paste("CPU:", ms$cpu),
    paste("Cores (logical / physical):", ms$logical_cores, "/", ms$physical_cores),
    paste("Memory (GB):", ms$memory_gb),
    paste("BLAS:", si$BLAS), paste("LAPACK:", si$LAPACK),
    paste("INLA num.threads:", inla_thr,
          "(INLA is the only multithreaded step; other packages ran on one thread)"),
    "", "## Package versions",
    paste0("  ", names(vers), " ", vers),
    if (!is.null(extra)) c("", extra),
    "", "## Full sessionInfo()", utils::capture.output(si))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, path)
  invisible(lines)
}
