#!/usr/bin/env Rscript

## ============================================================================
## Figures for the regression simulation study (Section 3.4, Appendix E.2).
##
## Reads the tab-delimited output of sim_study_xy.R, averages each metric over
## replicates within a design cell, and writes one faceted bar chart per
## (k, clust_sep, sig_const) scenario, with metrics down the rows and sample
## sizes across the columns.  The layout matches plot_sim_y_paper.R so the two
## studies' figures are directly comparable.
##
## Usage:
##   RESULTS_DIR=/path/to/sim_results_xy Rscript plot_sim_xy_paper.R
##
## Run the quality check on the same files first --
##   Rscript code/check_results.R 50 <files>
## -- because a replicate lost to the tryCatch in ss_study_xy(), or a cell
## duplicated across runs, would silently distort the averages plotted here.
## ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

## ------------------------------------------------------------------ input ---
## Every result file in RESULTS_DIR is read and row-bound, so combining several
## runs is a matter of writing them to the same directory.
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = ".")
result_files <- Sys.glob(file.path(RESULTS_DIR, "study_xy_saveK_results*.txt"))
if (!length(result_files)) stop("no result files found under RESULTS_DIR=", RESULTS_DIR)
message("reading ", length(result_files), " result file(s):\n  ",
        paste(basename(result_files), collapse = "\n  "))

results <- do.call(rbind, lapply(result_files,
                                 function(f) data.frame(read.table(f, header = TRUE))))
results <- as.data.frame(results)

id_cols <- intersect(c("N","k","clust_sep","sig_const","rep","M"), names(results))
results <- results %>% mutate(N = as.integer(N))

## Metrics reported in the manuscript.  
metrics_keep <- c("mnARI","sdARI","essMean","essM2","KSmean","L2mean","speed")

## Column names are "<metric>_<method>", so a single split on "_" recovers both
## factors.  read.table() turns the hyphens in the method names into dots
## ("ds-ml" becomes "ds.ml"), which is why the recode table below carries both
## spellings.
long <- results %>%
  pivot_longer(
    cols = -all_of(id_cols),
    names_to = c("metric","method"),
    names_sep = "_",
    values_to = "value"
  )

## ------------------------------------------------- averaging and filtering --
summ <- long %>%
  filter(metric %in% metrics_keep) %>%
  group_by(N, k, clust_sep, sig_const, metric, method) %>%
  summarise(
    mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
    sd   = if (all(is.na(value))) NA_real_ else sd(value,   na.rm = TRUE),
    n    = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  ## A method that does not run in a cell -- the direct sampler above N = 13 --
  ## is all-NA there and is dropped rather than plotted as an empty bar.
  filter(!is.na(mean))

## --------------------------------------------------- relabel and order ------
## Translate the internal metric and method codes into the names used in the
## manuscript.  Both hyphen and dot spellings are accepted so the script reads
## result files whatever their column-name mangling.
metric_map <- c(
  essMean = "ESS Mean",
  essM2   = "ESS M2",
  mnARI   = "ARI Mean",
  sdARI   = "SD ARI",
  KSmean  = "KS",
  L2mean  = "L2",
  speed   = "Time"
)

method_map <- c(
  "mcmc"        = "MCMC",
  "ds"          = "DS",
  "ds-gibbs"    = "MC-MCMC",
  "ds.gibbs"    = "MC-MCMC",
  ## Earlier result files named MC-MCMC "ds-mcmc"; kept so they still read.
  "ds-mcmc"     = "MC-MCMC",
  "ds_mcmc"     = "MC-MCMC",
  "ds.mcmc"     = "MC-MCMC",
  
  ## Subset-based constrained sampler.
  "ds10-obs"    = "DS-Const 0.1",
  "ds25-obs"    = "DS-Const 0.25",
  "ds50-obs"    = "DS-Const 0.5",
  "ds10.obs"    = "DS-Const 0.1",
  "ds25.obs"    = "DS-Const 0.25",
  "ds50.obs"    = "DS-Const 0.5",
  
  ## DS-Const-MAP is run only in the intercept-only study; these entries let
  ## the same map serve both scripts.
  "ds10-obsmap" = "DS-Const-MAP 0.1",
  "ds25-obsmap" = "DS-Const-MAP 0.25",
  "ds50-obsmap" = "DS-Const-MAP 0.5",
  "ds10.obsmap" = "DS-Const-MAP 0.1",
  "ds25.obsmap" = "DS-Const-MAP 0.25",
  "ds50.obsmap" = "DS-Const-MAP 0.5",
  
  "ds-ml"       = "DS-ML",
  "ds.ml"       = "DS-ML"
)

## Colourblind-safe palette, shared with plot_sim_y_paper.R: DS-Const in blues,
## DS-Const-MAP in reds, the unsubsetted methods in green, purple, grey, black.
fill_vals <- c(
  "DS"               = "#009E73",
  "DS-Const 0.1"     = "#deebf7",
  "DS-Const 0.25"    = "#9ecae1",
  "DS-Const 0.5"     = "#3182bd",
  "DS-Const-MAP 0.1" = "#fee0d2",
  "DS-Const-MAP 0.25"= "#fc9272",
  "DS-Const-MAP 0.5" = "#cb181d",
  "DS-ML"            = "#6a51a3",
  "MC-MCMC"          = "#999999",
  "MCMC"             = "#000000"
)
## Opacity increases with the subset fraction; the alpha legend is suppressed
## because the fraction is already part of the method name.
alpha_vals <- c(`0.1` = 0.55, `0.25` = 0.75, `0.5` = 1.00)


## Plotting order.  Levels absent from a given run are dropped below.
method_levels <- c(
  "DS",
  "DS-Const 0.1","DS-Const 0.25","DS-Const 0.5",
  "DS-Const-MAP 0.1","DS-Const-MAP 0.25","DS-Const-MAP 0.5",
  "DS-ML","MC-MCMC","MCMC"
)

metric_levels <- c("ARI Mean","SD ARI","ESS Mean","ESS M2", "KS","L2","Time")

summ_clean <- summ %>%
  mutate(
    metric = dplyr::recode(metric, !!!metric_map, .default = metric),
    method = dplyr::recode(method, !!!method_map, .default = method),
    frac = dplyr::case_when(
      grepl("0\\.1$",  method) ~ 0.10,
      grepl("0\\.25$", method) ~ 0.25,
      grepl("0\\.5$",  method) ~ 0.50,
      TRUE ~ NA_real_
    ),
    fill_key = method
  ) %>%
  mutate(
    ## Restrict to the levels actually present so empty methods disappear.
    method = factor(method, levels = intersect(method_levels, unique(method))),
    metric = factor(metric, levels = intersect(metric_levels, unique(metric))),
    N      = factor(N)
  )

## One faceted bar chart for a single scenario: metrics down the rows, sample
## sizes across the columns, bars labelled with the mean over replicates.
##
## @param df_one  rows of summ_clean for one (k, clust_sep, sig_const) cell.
## @param k,clust_sep,sig_const  the design values, used for the title.
plot_scenario <- function(df_one, k, clust_sep, sig_const) {

  ## Separation codes as named in the manuscript.
  overlap_map <- c(`1` = "Moderate Overlap",
                   `2` = "Large Overlap",
                   `3` = "No Overlap")
  overlap_lab <- overlap_map[as.character(clust_sep)]
  if (is.na(overlap_lab)) overlap_lab <- as.character(clust_sep)
  
  ggplot(df_one, aes(x = method, y = mean)) +
    geom_col(
      aes(fill = fill_key, alpha = factor(frac, levels = c(0.1, 0.25, 0.5))),
      width = 0.7,
      color = "black",
      linewidth = 0.35
    ) +
    ## Value labels: ESS to the nearest integer, times to one decimal, and
    ## the remaining metrics to three.
    geom_text(
      aes(label = case_when(
        metric %in% c("ESS Mean","ESS M2") ~ as.character(round(mean)),
        metric == "Time"                  ~ sprintf("%.1f", mean),
        TRUE                              ~ sprintf("%.3f", mean)
      )),
      vjust = -0.45,
      size = 3.8,
      fontface = "bold"
    ) +
    facet_grid(
      rows = vars(metric),
      cols = vars(N),
      scales = "free_y",
      labeller = labeller(
        N = function(x) paste0("N = ", x)
      )
    ) +
    ## Headroom above the tallest bar, and no clipping, so the value labels fit.
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = fill_vals, breaks = names(fill_vals)) +
    scale_alpha_manual(values = alpha_vals, na.value = 1, guide = "none") +
    labs(
      title = sprintf("True K = %s, %s", k, overlap_lab),
      x = NULL, y = "Mean"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),

      strip.text.x = element_text(face = "bold", size = 15),
      strip.text.y = element_text(face = "bold", size = 12),

      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),

      panel.spacing = grid::unit(0.9, "lines"),
      plot.margin   = margin(8, 18, 8, 8),

      legend.position = "bottom",
      legend.box      = "vertical",
      legend.title    = element_text(size = 15, face = "bold"),
      legend.text     = element_text(size = 13),
      legend.key.size = grid::unit(0.9, "lines")
    ) +
    guides(
      fill = guide_legend(title = "Method", nrow = 2, byrow = TRUE)
    )
}

## ----------------------------------------------------------- output --------
## One figure per scenario, named k<k>_sep<clust_sep>_sig<sig_const>.png.
plots <- summ_clean %>%
  group_by(k, clust_sep, sig_const) %>%
  group_map(~ plot_scenario(.x, .y$k[[1]], .y$clust_sep[[1]], .y$sig_const[[1]]), .keep = TRUE)

names(plots) <- summ_clean %>%
  distinct(k, clust_sep, sig_const) %>%
  mutate(name = paste0("k", k, "_sep", clust_sep, "_sig", sig_const)) %>%
  pull(name)

dir.create("plots_metric_byN_xy_sims_with10000", showWarnings = FALSE)
purrr::iwalk(
  plots,
  ~ ggsave(file.path("plots_metric_byN_xy_sims_with10000", paste0(.y, ".png")),
           plot = .x, width = 16, height = 11.5, dpi = 300)
)

