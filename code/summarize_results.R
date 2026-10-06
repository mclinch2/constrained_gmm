#!/usr/bin/env Rscript
## ============================================================================
## Read finished simulation output and print it in the manuscript's table
## layout, so a run can be compared against Table 1 / Table 2 directly.
##
## Base R only, so it runs wherever the study ran.
##
## Usage:
##   Rscript code/summarize_results.R <dir>                  # every result file in dir
##   Rscript code/summarize_results.R <file> [file ...]      # named files
##   Rscript code/summarize_results.R <dir> --reps 5         # expected replicates
##
## It prints, for each (k, clust_sep, N) cell and each method: the number of
## replicates actually contributing, mean ARI, ESS of the first and second
## moments, the KS and L2 distances, and CPU seconds -- each averaged over
## replicates, which is what the manuscript tables report.
##
## Run code/check_results.R first for the integrity checks (duplicated cells,
## short replicate counts); this script summarizes, it does not validate.
## ============================================================================

args <- commandArgs(TRUE)
if (!length(args)) stop("usage: summarize_results.R <dir|files...> [--reps N]")

reps_expected <- NA_integer_
if ("--reps" %in% args) {
  i <- match("--reps", args)
  reps_expected <- as.integer(args[i + 1L])
  args <- args[-c(i, i + 1L)]
}

files <- unlist(lapply(args, function(a) {
  if (dir.exists(a)) c(Sys.glob(file.path(a, "study1_saveK_results*.txt")),
                       Sys.glob(file.path(a, "study_xy_saveK_results*.txt")))
  else a
}))
files <- files[file.exists(files)]
if (!length(files)) stop("no result files found in: ", paste(args, collapse = ", "))

res <- do.call(rbind, lapply(files, function(f)
  read.table(f, header = TRUE, check.names = TRUE)))
cat("files read : ", length(files), "\n", sep = "")
cat("rows       : ", nrow(res), "\n", sep = "")

## Which study?  The method columns differ, and so do the manuscript tables.
is_xy <- any(grepl("^KSmean_", names(res))) && !any(grepl("^mnKS_", names(res)))
cat("study      : ", if (is_xy) "regression (Table 2)" else "intercept-only (Table 1)",
    "\n\n", sep = "")

## Internal method code -> manuscript name, in the tables' own order.
methods <- if (is_xy) {
  c("ds" = "DS", "ds.ml" = "DS-ML", "ds10.obs" = "DS-Const 0.1",
    "ds25.obs" = "DS-Const 0.25", "ds50.obs" = "DS-Const 0.5",
    "ds.gibbs" = "MC-MCMC", "mcmc" = "MCMC")
} else {
  c("ds" = "DS", "ds.ml" = "DS-ML", "ds10.obs" = "DS-Const 0.1",
    "ds25.obs" = "DS-Const 0.25", "ds50.obs" = "DS-Const 0.5",
    "ds10.obsmap" = "DS-Const-MAP 0.1", "ds25.obsmap" = "DS-Const-MAP 0.25",
    "ds50.obsmap" = "DS-Const-MAP 0.5", "ds.gibbs" = "MC-MCMC", "mcmc" = "MCMC")
}
## KSmean_/L2mean_ are the posterior-MEAN-curve distances -- the manuscript's
## definition.  mnKS_/mnL2_, written only by the intercept-only study, are the
## mean over draws of the per-draw distance, a different quantity, and are not
## reported here.
metrics <- c(ARI = "mnARI", sdARI = "sdARI", ESS1 = "essMean", ESS2 = "essM2",
             KS = "KSmean", L2 = "L2mean", CPUs = "speed")

sep_lab <- c("1" = "Moderate overlap", "2" = "Large overlap", "3" = "No overlap")

cells <- unique(res[, c("k", "clust_sep", "N")])
cells <- cells[order(cells$k, cells$clust_sep, cells$N), , drop = FALSE]

for (i in seq_len(nrow(cells))) {
  cl <- cells[i, ]
  d  <- res[res$k == cl$k & res$clust_sep == cl$clust_sep & res$N == cl$N, , drop = FALSE]
  lab <- sep_lab[as.character(cl$clust_sep)]
  if (is.na(lab)) lab <- paste("clust_sep", cl$clust_sep)
  cat(sprintf("=== True K = %d, %s, N = %d  (%d replicate rows) ===\n",
              cl$k, lab, cl$N, nrow(d)))

  cat(sprintf("  %-18s %5s %8s %8s %8s %8s %8s %8s %10s\n",
              "Method", "reps", "ARI", "sd ARI", "ESS(1)", "ESS(2)", "KS", "L2", "CPU (s)"))
  for (code in names(methods)) {
    cols <- paste0(metrics, "_", code)
    have <- cols %in% names(d)
    if (!any(have)) next
    v <- vapply(seq_along(cols), function(j)
      if (have[j]) mean(d[[cols[j]]], na.rm = TRUE) else NA_real_, numeric(1))
    names(v) <- names(metrics)
    n_ok <- if (paste0("mnARI_", code) %in% names(d))
      sum(!is.na(d[[paste0("mnARI_", code)]])) else 0L
    if (n_ok == 0L) next   # method not run in this cell (e.g. DS above N = 13)
    flag <- if (!is.na(reps_expected) && n_ok != reps_expected) "  <- incomplete" else ""
    cat(sprintf("  %-18s %5d %8.3f %8.3f %8.0f %8.0f %8.4f %8.4f %10.1f%s\n",
                methods[[code]], n_ok, v[["ARI"]], v[["sdARI"]],
                v[["ESS1"]], v[["ESS2"]], v[["KS"]], v[["L2"]], v[["CPUs"]], flag))
  }
  cat("\n")
}