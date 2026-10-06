#!/usr/bin/env Rscript

## ============================================================================
## Figures for the intercept-only simulation study (Section 3.3, Appendix E.1).
##
## Reads the tab-delimited output of sim_study_y.R, averages each metric over
## replicates within a design cell, and writes one faceted bar chart per
## (k, clust_sep, sig_const) scenario, with metrics down the rows and sample
## sizes across the columns.
##
## Usage:
##   RESULTS_DIR=/path/to/sim_results_y Rscript plot_sim_y_paper.R
##
## Run the quality check on the same files first --
##   Rscript code/check_results.R 50 <files>
## -- because a replicate lost to the tryCatch in ss_study(), or a cell
## duplicated across runs, would silently distort the averages plotted here.
## ============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(xtable)
library(stringr)
library(ggplot2)

## ------------------------------------------------------------------ input ---
## Every result file in RESULTS_DIR is read and row-bound, so combining several
## runs is a matter of writing them to the same directory.
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", unset = ".")
result_files <- Sys.glob(file.path(RESULTS_DIR, "study1_saveK_results*.txt"))
if (!length(result_files)) stop("no result files found under RESULTS_DIR=", RESULTS_DIR)
message("reading ", length(result_files), " result file(s):\n  ",
        paste(basename(result_files), collapse = "\n  "))

results <- do.call(rbind, lapply(result_files,
                                 function(f) data.frame(read.table(f, header = TRUE))))
results <- as.data.frame(results)

stopifnot(all(c("N","k","clust_sep","sig_const", "M") %in% names(results)))
results <- results %>% mutate(N = as.integer(N))

## Column names are "<metric>_<method>", so a single split on "_" recovers both
## factors.  read.table() has already turned the hyphens in the method names
## into dots ("ds-ml" becomes "ds.ml"), which is why the filters and the recode
## table below use the dotted forms.
long <- results %>%
  pivot_longer(
    cols = -c(N, k, clust_sep, sig_const, rep, M),
    names_to = c("metric", "method"),
    names_sep = "_",
    values_to = "value"
  )

## ------------------------------------------------- averaging and filtering --
## Drop the blocked-Gibbs reference sampler, which is not reported in the
## manuscript, along with the diagnostic-only columns.  ISE is dropped in
## favour of its square root, L2, which is what the tables report.
summ <- long %>%
  filter(!(method %in% c("ds.mcmc"))) %>%
  filter(!(metric %in% c("postMean", "postMed", "sdKS", "sdISE", "sdL2",
                         "ISEmed", "mnISE", "ISEmean"))) %>%
  group_by(N, k, clust_sep, sig_const, metric, method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    n    = sum(!is.na(value)),
    .groups = "drop"
  )

## --------------------------------------------------- relabel and order ------
## Translate the internal metric and method codes into the names used in the
## manuscript, and fix the plotting order.
summ_clean <- summ %>%
  mutate(
    metric_code = as.character(metric)   # retain the internal code for filtering
  ) %>%
  ## Of the density metrics, keep only those defined on the posterior MEAN
  ## curve: KSmean and L2mean are the Appendix E.1 definitions.  The per-draw
  ## averages (mnKS, mnL2) and the median-curve variants are diagnostics.
  filter(!metric_code %in% c("mnKS", "mnL2")) %>%
  filter(!metric_code %in% c("KSmed", "L2med")) %>%
  mutate(
    method = factor(method),
    N      = factor(N),
    
    metric = dplyr::case_when(
      metric_code == "essMean" ~ "ESS Mean",
      metric_code == "essM2"   ~ "ESS M2",
      metric_code == "mnARI"   ~ "ARI Mean",
      metric_code == "sdARI"   ~ "SD ARI",
      metric_code == "KSmean"  ~ "KS",
      metric_code == "L2mean"  ~ "L2",
      metric_code == "speed"   ~ "Time",
      TRUE                     ~ metric_code
    ),
    
    method = recode(method,
                    "mcmc"         = "MCMC",
                    "ds"           = "DS",
                    
                    "ds10.obs"     = "DS-Const 0.1",
                    "ds25.obs"     = "DS-Const 0.25",
                    "ds50.obs"     = "DS-Const 0.5",
                    
                    "ds10.obsmap"  = "DS-Const-MAP 0.1",
                    "ds25.obsmap"  = "DS-Const-MAP 0.25",
                    "ds50.obsmap"  = "DS-Const-MAP 0.5",
                    
                    "ds.gibbs"     = "MC-MCMC",
                    "ds.ml"        = "DS-ML"
    ),
    
    method = factor(method, levels = c(
      "DS",
      "DS-Const 0.1", "DS-Const 0.25", "DS-Const 0.5",
      "DS-Const-MAP 0.1", "DS-Const-MAP 0.25", "DS-Const-MAP 0.5",
      "DS-ML",
      "MC-MCMC",
      "MCMC"
    )),
    
    metric = factor(metric, levels = c(
      "ARI Mean", "SD ARI", "ESS Mean", "ESS M2",
      "KS", "L2", "Time"
    ))
  ) %>%
  ## Aesthetic keys: the subset fraction drives opacity, the method family
  ## drives the fill hue.
  mutate(
    frac = case_when(
      grepl("0.1$",  method) ~ 0.10,
      grepl("0.25$", method) ~ 0.25,
      grepl("0.5$",  method) ~ 0.50,
      TRUE ~ NA_real_
    ),
    group = case_when(
      grepl("^DS-Const-MAP", method) ~ "DS-Const-MAP",
      grepl("^DS-Const",     method) ~ "DS-Const",
      method == "MC-MCMC"            ~ "MC-MCMC",
      method == "MCMC"               ~ "MCMC",
      method == "DS-ML"              ~ "DS-ML",
      method == "DS"                 ~ "DS",
      TRUE                           ~ "Other"
    ),
    fill_key = method
  )

## ---------------------------------------------------------- aesthetics ------
## Colourblind-safe palette: one hue per method family, shaded by subset
## fraction -- DS-Const in blues, DS-Const-MAP in reds, and the three
## unsubsetted methods in green, purple, grey and black.
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
  
  ggplot2::ggplot(df_one, ggplot2::aes(x = method, y = mean)) +
    ggplot2::geom_col(
      ggplot2::aes(
        fill  = fill_key,
        alpha = factor(frac, levels = c(0.1, 0.25, 0.5))
      ),
      width = 0.7,
      color = "black",
      size  = 0.35
    ) +
    ## Value labels: ESS to the nearest integer, times to one decimal, and
    ## the remaining metrics to three.
    ggplot2::geom_text(
      ggplot2::aes(
        label = dplyr::case_when(
          metric %in% c("ESS Mean", "ESS M2") ~ as.character(round(mean)),
          metric == "Time"                    ~ sprintf("%.1f", mean),
          TRUE                                ~ sprintf("%.3f", mean)
        )
      ),
      vjust    = -0.55,
      size     = 4,
      fontface = "bold"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(metric),
      cols = ggplot2::vars(N),
      scales = "free_y",
      labeller = ggplot2::labeller(
        N = function(x) paste0("N = ", x)
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.18))
    ) +
    ggplot2::scale_fill_manual(values = fill_vals, breaks = names(fill_vals)) +
    ggplot2::scale_alpha_manual(values = alpha_vals, na.value = 1, guide = "none") +
    ggplot2::labs(
      title = sprintf("True K = %s, %s", k, overlap_lab),
      x = NULL, y = "Mean"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),

      strip.text.x = ggplot2::element_text(face = "bold", size = 15),
      strip.text.y = ggplot2::element_text(face = "bold", size = 12),

      legend.position = "bottom",
      legend.title    = ggplot2::element_text(size = 15, face = "bold"),
      legend.text     = ggplot2::element_text(size = 13),
      legend.key.size = grid::unit(0.9, "lines"),
      
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(title = "Method")
    )
}


## ----------------------------------------------------------- output --------
## One figure per scenario, named k<k>_sep<clust_sep>_sig<sig_const>.png.
plots_main <- summ_clean %>%
  group_by(k, clust_sep, sig_const) %>%
  group_map(
    ~ plot_scenario(.x,
                    k         = .y$k[[1]],
                    clust_sep = .y$clust_sep[[1]],
                    sig_const = .y$sig_const[[1]]),
    .keep = TRUE
  )

plot_names_main <- summ_clean %>%
  distinct(k, clust_sep, sig_const) %>%
  mutate(name = str_glue("k{k}_sep{clust_sep}_sig{sig_const}")) %>%
  pull(name)

names(plots_main) <- plot_names_main

dir.create("plots_metric_byN_M25_int_only", showWarnings = FALSE)
purrr::iwalk(
  plots_main,
  ~ ggplot2::ggsave(file.path("plots_metric_byN_M25_int_only", paste0(.y, ".png")),
                    plot = .x, width = 19, height = 11.5, dpi = 300)
)
