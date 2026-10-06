#!/usr/bin/env Rscript

## ============================================================================
## Quality checks on simulation output before plotting or tabulating.
##
## This script exits with a nonzero status if it detects:
##
##   1. duplicated design-cell/replicate combinations;
##   2. an unexpected number of replicates in a design cell;
##   3. incomplete method results;
##   4. metrics outside their valid ranges.
##
## Because the analysis shell scripts use `set -e`, a failed check stops
## the workflow before summaries or figures are produced.
##
## Usage:
##
##   Rscript code/check_results.R \
##       <expected_reps> <file1.txt> [file2.txt ...]
##
## ============================================================================

args <- commandArgs(TRUE)

if (length(args) < 2L) {
    stop(
        "Usage: check_results.R ",
        "<expected_reps> <file1.txt> [file2.txt ...]",
        call. = FALSE
    )
}

expected <- suppressWarnings(as.integer(args[1]))
files <- args[-1]

if (is.na(expected) || expected < 1L) {
    stop(
        "expected_reps must be a positive integer. Received: ",
        args[1],
        call. = FALSE
    )
}

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0L) {
    stop(
        "The following result files do not exist:\n  ",
        paste(missing_files, collapse = "\n  "),
        call. = FALSE
    )
}


# ============================================================
# Read results
# ============================================================

res <- do.call(
    rbind,
    lapply(files, function(file_name) {
        result <- read.table(
            file_name,
            header = TRUE,
            check.names = TRUE
        )

        result$.src <- basename(file_name)
        result
    })
)

cat(
    "Files: ", length(files),
    "   Rows: ", nrow(res),
    "\n\n",
    sep = ""
)

id_columns <- c(
    "N",
    "k",
    "clust_sep",
    "sig_const",
    "M",
    "rep"
)

if (!all(id_columns %in% names(res))) {
    stop(
        "Missing ID columns: ",
        paste(
            setdiff(id_columns, names(res)),
            collapse = ", "
        ),
        call. = FALSE
    )
}

# Track whether any check fails.
check_failed <- FALSE

key <- do.call(
    paste,
    c(res[id_columns], sep = "|")
)

cell_columns <- setdiff(id_columns, "rep")

cell <- do.call(
    paste,
    c(res[cell_columns], sep = "|")
)


# ============================================================
# 1. Duplicate design-cell/replicate combinations
# ============================================================

duplicate_counts <- table(key)
duplicates <- duplicate_counts[duplicate_counts > 1L]

if (length(duplicates) > 0L) {
    check_failed <- TRUE

    cat(
        "1. FAILED: duplicated (design cell, replicate) ",
        "combinations were found.\n",
        sep = ""
    )

    print(utils::head(duplicates, 20L))

    duplicate_rows <- duplicated(key) | duplicated(key, fromLast = TRUE)

    cat("\nFiles containing duplicated rows:\n")

    for (file_name in unique(res$.src[duplicate_rows])) {
        cat("  ", file_name, "\n", sep = "")
    }

    cat(
        "\nThese rows would be counted more than once when ",
        "the result files are combined.\n\n",
        sep = ""
    )
} else {
    cat(
        "1. PASSED: no duplicated (design cell, replicate) ",
        "combinations.\n\n",
        sep = ""
    )
}


# ============================================================
# 2. Replicate count in each design cell
# ============================================================

replicate_counts <- tapply(
    res$rep,
    cell,
    function(x) length(unique(x))
)

bad_counts <- replicate_counts[
    replicate_counts != expected
]

if (length(bad_counts) > 0L) {
    check_failed <- TRUE

    cat(
        "2. FAILED: the following design cells do not have ",
        expected,
        " unique replicates:\n",
        sep = ""
    )

    print(bad_counts)
    cat("\n")
} else {
    cat(
        "2. PASSED: every design cell has exactly ",
        expected,
        " unique replicates.\n\n",
        sep = ""
    )
}


# ============================================================
# Identify metric and method columns
# ============================================================

value_columns <- setdiff(
    names(res),
    c(id_columns, ".src")
)

first_underscore <- regexpr(
    "_",
    value_columns,
    fixed = TRUE
)

metric <- ifelse(
    first_underscore > 0L,
    substring(value_columns, 1L, first_underscore - 1L),
    NA_character_
)

method <- ifelse(
    first_underscore > 0L,
    substring(value_columns, first_underscore + 1L),
    NA_character_
)

valid_value_column <- !is.na(metric) & !is.na(method)


# ============================================================
# 3. Usable replicate counts by method and metric
# ============================================================

incomplete <- character(0L)

for (cell_value in unique(cell)) {
    rows <- which(cell == cell_value)
    current_N <- res$N[rows[1L]]

    for (j in which(valid_value_column)) {
        usable <- sum(
            !is.na(res[rows, value_columns[j]])
        )

        # The enumerating direct sampler is intentionally not run
        # for N > 13.
        if (method[j] == "ds" && current_N > 13) {
            next
        }

        if (usable != expected) {
            incomplete <- c(
                incomplete,
                sprintf(
                    "%s  %s: %d/%d",
                    cell_value,
                    value_columns[j],
                    usable,
                    expected
                )
            )
        }
    }
}

if (length(incomplete) > 0L) {
    check_failed <- TRUE

    cat(
        "3. FAILED: method/metric entries with a usable ",
        "replicate count different from ",
        expected,
        ":\n",
        sep = ""
    )

    cat(
        paste0(
            "   ",
            utils::head(incomplete, 30L),
            collapse = "\n"
        ),
        "\n"
    )

    if (length(incomplete) > 30L) {
        cat(
            "   ... and ",
            length(incomplete) - 30L,
            " more.\n",
            sep = ""
        )
    }

    cat(
        "\nThese summaries would be based on a different number ",
        "of replicates than intended.\n\n",
        sep = ""
    )
} else {
    cat(
        "3. PASSED: all applicable methods have ",
        expected,
        " usable replicates in every design cell.\n\n",
        sep = ""
    )
}


# ============================================================
# 4. Metric ranges
# ============================================================

cat("4. Metric ranges:\n")

for (metric_name in sort(unique(metric[valid_value_column]))) {
    selected_columns <- value_columns[
        valid_value_column &
        metric == metric_name
    ]

    values <- unlist(
        res[selected_columns],
        use.names = FALSE
    )

    values <- values[!is.na(values)]

    if (length(values) > 0L) {
        cat(
            sprintf(
                "   %-10s [%12.4f, %12.4f]\n",
                metric_name,
                min(values),
                max(values)
            )
        )
    }
}

invalid_columns <- character(0L)

for (j in which(valid_value_column)) {
    values <- res[[value_columns[j]]]
    values <- values[!is.na(values)]

    if (length(values) == 0L) {
        next
    }

    if (any(!is.finite(values))) {
        invalid_columns <- c(
            invalid_columns,
            value_columns[j]
        )

        next
    }

    if (
        grepl("ARI", metric[j]) &&
        (min(values) < -1 || max(values) > 1)
    ) {
        invalid_columns <- c(
            invalid_columns,
            value_columns[j]
        )
    }

    if (
        metric[j] == "speed" &&
        min(values) < 0
    ) {
        invalid_columns <- c(
            invalid_columns,
            value_columns[j]
        )
    }

    if (
        grepl("^ess", metric[j]) &&
        min(values) < 0
    ) {
        invalid_columns <- c(
            invalid_columns,
            value_columns[j]
        )
    }

    if (
        grepl("^(KS|L2|ISE|mnKS|mnL2|mnISE|sdKS|sdL2|sdISE)",
              metric[j]) &&
        min(values) < 0
    ) {
        invalid_columns <- c(
            invalid_columns,
            value_columns[j]
        )
    }
}

invalid_columns <- unique(invalid_columns)

if (length(invalid_columns) > 0L) {
    check_failed <- TRUE

    cat(
        "\n4. FAILED: invalid or out-of-range values were found in:\n"
    )

    for (column_name in invalid_columns) {
        cat("   ", column_name, "\n", sep = "")
    }
} else {
    cat("\n4. PASSED: no invalid metric values were detected.\n")
}


# ============================================================
# Final status
# ============================================================

cat("\n============================================================\n")

if (check_failed) {
    cat("RESULT INTEGRITY CHECK FAILED\n")
    cat("Summaries and figures should not be produced from these files.\n")
    cat("============================================================\n")

    quit(
        save = "no",
        status = 1L
    )
}

cat("RESULT INTEGRITY CHECK PASSED\n")
cat("============================================================\n")

quit(
    save = "no",
    status = 0L
)