## =====================================================================
## 02_descriptives_moran.R -- adjacency, Table 1, Moran's I by year, Fig 2
##   in : output/panel_<phenotype>_<denominator>.csv
##   out: output/sc_adj.graph (INLA graph file), output/sc_counties.rds,
##        output/table1_data_profile.csv, output/tableS2_moran_by_year.csv,
##        output/fig2_statewide_trend.png
## =====================================================================
if (!exists("cfg", envir = globalenv())) source(file.path("R", "00_config.R"))
cfg <- get("cfg", envir = globalenv())
suppressWarnings(suppressMessages({ library(sf); library(spdep); library(ggplot2) }))

panel <- readr::read_csv(file.path(cfg$dir_out,
                                   sprintf("panel_%s_%s.csv", cfg$phenotype, cfg$denominator)),
                         col_types = readr::cols(fips = readr::col_character()))

## ---- 1. Boundaries and adjacency --------------------------------------
sc <- load_sc_boundaries()

check_that(nrow(sc) == 46, "boundary file has 46 SC counties")
check_that(!any(is.na(sc$county)), "every boundary FIPS matches the crosswalk")
check_that(setequal(sc$fips, unique(panel$fips)), "boundary and panel FIPS agree")

nb <- spdep::poly2nb(sc, queen = (cfg$adjacency == "queen"), row.names = sc$fips)
n_islands <- sum(spdep::card(nb) == 0)
n_comp    <- spdep::n.comp.nb(nb)$nc
check_that(n_islands == 0, "no island counties in the adjacency graph")
check_that(n_comp == 1, "adjacency graph is fully connected")
log_step(cfg$adjacency, " adjacency: mean ", round(mean(spdep::card(nb)), 2),
         " neighbours (range ", min(spdep::card(nb)), "-", max(spdep::card(nb)), ")")

graph_file <- file.path(cfg$dir_out, paste0("sc_adj_", cfg$adjacency, ".graph"))
spdep::nb2INLA(graph_file, nb)
saveRDS(list(sc = sc, nb = nb, graph_file = graph_file),
        file.path(cfg$dir_out, "sc_spatial.rds"))

## county_id must follow the ORDER OF THE GRAPH FILE, not alphabetical order.
county_index <- tibble(fips = sc$fips, county_id = seq_len(nrow(sc)))
panel <- panel |> left_join(county_index, by = "fips") |>
  mutate(year_id = as.integer(factor(year, levels = sort(cfg$years))))
readr::write_csv(panel, file.path(cfg$dir_out,
  sprintf("panel_%s_%s_indexed.csv", cfg$phenotype, cfg$denominator)))

## ---- 2. Table 1: county-year data profile -----------------------------
q <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))
tab1 <- tibble(
  metric = c("County-years", "Counties", "Years",
             "Denominator: min / median / max",
             "Denominator: IQR",
             "County-years with N < 100 (n, %)",
             "County-years with N < 50 (n, %)",
             "Numerator: min / median / max",
             "County-years with numerator 0",
             "Crude encounter rate: min / median / max"),
  value = c(
    nrow(panel), n_distinct(panel$fips), n_distinct(panel$year),
    sprintf("%d / %d / %d", min(panel$deno), as.integer(stats::median(panel$deno)), max(panel$deno)),
    sprintf("%d-%d", as.integer(q(panel$deno, .25)), as.integer(q(panel$deno, .75))),
    sprintf("%d (%.1f%%)", sum(panel$deno < 100), 100 * mean(panel$deno < 100)),
    sprintf("%d (%.1f%%)", sum(panel$deno <  50), 100 * mean(panel$deno <  50)),
    sprintf("%d / %d / %d", min(panel$nume), as.integer(stats::median(panel$nume)), max(panel$nume)),
    sum(panel$nume == 0),
    sprintf("%.2f / %.2f / %.2f", min(panel$encounter_pct, na.rm = TRUE),
            stats::median(panel$encounter_pct, na.rm = TRUE), max(panel$encounter_pct, na.rm = TRUE))
  )
)
readr::write_csv(tab1, file.path(cfg$dir_out, "table1_data_profile.csv"))
print(tab1, n = 50)

statewide <- panel |> group_by(year) |>
  summarise(deno = sum(deno), nume = sum(nume),
            wrate_pct = 100 * nume / deno, .groups = "drop") |>
  mutate(lo = 100 * qbeta(.025, nume + .5, deno - nume + .5),
         hi = 100 * qbeta(.975, nume + .5, deno - nume + .5))
readr::write_csv(statewide, file.path(cfg$dir_out, "table1_statewide_trend.csv"))

## ---- 3. Moran's I by year (Supp Table S2) -----------------------------
lw <- spdep::nb2listw(nb, style = "W")
set.seed(cfg$cv_seed)
moran_tab <- purrr::map_dfr(sort(unique(panel$year)), function(y) {
  d <- panel |> filter(year == y) |> arrange(match(fips, sc$fips))
  mc <- spdep::moran.mc(d$encounter_pct, lw, nsim = cfg$moran_nsim, zero.policy = TRUE)
  tibble(year = y, morans_i = unname(mc$statistic),
         expected = -1 / (nrow(d) - 1), p_perm = mc$p.value)
}) |>
  mutate(p_fdr = p.adjust(p_perm, method = "BH"),
         significant = p_fdr < 0.05)
readr::write_csv(moran_tab, file.path(cfg$dir_out, "tableS2_moran_by_year.csv"))
log_step("Moran's I significant (FDR<0.05) in ", sum(moran_tab$significant),
         " of ", nrow(moran_tab), " years; I ranged ",
         round(min(moran_tab$morans_i), 3), " to ", round(max(moran_tab$morans_i), 3))

## ---- 4. Figure 2: statewide trend -------------------------------------
p2 <- ggplot(statewide, aes(year, wrate_pct)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15) +
  geom_line(linewidth = .8) + geom_point(size = 1.6) +
  scale_x_continuous(breaks = seq(min(cfg$years), max(cfg$years), 2)) +
  labs(x = NULL, y = "PLWH with a SUD-related hospital encounter (%)") +
  theme_minimal(base_size = 11)
ggsave(file.path(cfg$dir_out, "fig2_statewide_trend.png"), p2,
       width = 7, height = 4.2, dpi = 300)

writeLines(capture.output(sessionInfo()),
           file.path(cfg$dir_log, "sessionInfo_02.txt"))
