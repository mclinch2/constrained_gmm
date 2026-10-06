#!/bin/bash
#SBATCH --job-name=AnalyseY
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/analysey_%j.out
#SBATCH --error=logs/analysey_%j.err

set -euo pipefail

# ============================================================
# Y-only simulation analysis
# ============================================================
#
# Submit from the root of the cloned repository:
#
#   cd /path/to/fmm_manuscript
#   mkdir -p logs
#   sbatch code/run_analysis_y.sh
#
# Analyze final 50-replicate results:
#
#   sbatch \
#       --export=ALL,FMM_REPS=50,FMM_OUT=/path/to/sim_results_y \
#       code/run_analysis_y.sh
#
# To specify a particular R library:
#
#   sbatch \
#       --export=ALL,FMM_R_LIB=/path/to/R/library \
#       code/run_analysis_y.sh
#
# Analysis steps:
#
#   1. check_results.R
#      Checks for missing replicates, duplicate cells, and other
#      result-integrity problems.
#
#   2. summarize_results.R
#      Summarizes the results using the manuscript Table 1 layout.
#
#   3. plot_sim_y_paper.R
#      Creates figures if the required packages are installed.
#
# Expected repository structure:
#
#   fmm_manuscript/
#   ├── code/
#   │   ├── run_analysis_y.sh
#   │   ├── check_results.R
#   │   ├── summarize_results.R
#   │   └── plot_sim_y_paper.R
#   ├── logs/
#   └── test_results_y/
#
# ============================================================


# ============================================================
# Repository locations
# ============================================================

# FMM_REPO_ROOT can explicitly specify the repository.
#
# Under SLURM, the default is the directory from which sbatch was
# called. Therefore, submit this script from the repository root.
#
# When run directly with bash, the repository is determined from
# the location of this script, assuming it is inside code/.

if [[ -n "${FMM_REPO_ROOT:-}" ]]; then
    REPO_ROOT="${FMM_REPO_ROOT}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
fi

CODE_DIR="${REPO_ROOT}/code"

# Make relative output paths consistent.
cd "${REPO_ROOT}"


# ============================================================
# Results settings
# ============================================================
#
# Defaults correspond to the five-replicate test produced by
# run_sim_y.sh.
# ============================================================

export FMM_OUT="${FMM_OUT:-${REPO_ROOT}/test_results_y}"
FMM_REPS="${FMM_REPS:-5}"

if [[ ! "${FMM_REPS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: FMM_REPS must be a positive integer."
    echo "Received: ${FMM_REPS}"
    exit 1
fi

if [[ ! -d "${FMM_OUT}" ]]; then
    echo "ERROR: Results directory does not exist:"
    echo "       ${FMM_OUT}"
    echo
    echo "The simulation may still be running, or FMM_OUT may be wrong."
    exit 1
fi


# ============================================================
# R environment
# ============================================================
#
# R_MODULE defaults to R/4.4.0. It can be overridden:
#
#   sbatch --export=ALL,R_MODULE=R/4.5.0 \
#       code/run_analysis_y.sh
#
# A user-specific R library can be supplied using:
#
#   sbatch --export=ALL,FMM_R_LIB=/path/to/R/library \
#       code/run_analysis_y.sh
#
# When FMM_R_LIB is supplied, R_LIBS is cleared so an inherited
# incompatible library is not added to .libPaths().
# ============================================================

R_MODULE="${R_MODULE-R/4.4.0}"

if command -v module >/dev/null 2>&1 && [[ -n "${R_MODULE}" ]]; then
    module purge
    module load "${R_MODULE}"
fi

if [[ -n "${FMM_R_LIB:-}" ]]; then
    if [[ ! -d "${FMM_R_LIB}" ]]; then
        echo "ERROR: FMM_R_LIB does not exist:"
        echo "       ${FMM_R_LIB}"
        exit 1
    fi

    export R_LIBS_USER="${FMM_R_LIB}"
    unset R_LIBS || true
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found."
    echo
    echo "Load an appropriate R module or make Rscript available in PATH."
    exit 1
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1


# ============================================================
# Verify analysis scripts
# ============================================================

REQUIRED_SCRIPTS=(
    "check_results.R"
    "summarize_results.R"
    "plot_sim_y_paper.R"
)

for script_name in "${REQUIRED_SCRIPTS[@]}"; do
    if [[ ! -f "${CODE_DIR}/${script_name}" ]]; then
        echo "ERROR: Could not find the required analysis script:"
        echo "       ${CODE_DIR}/${script_name}"
        exit 1
    fi
done


# ============================================================
# Record the code version
# ============================================================

CODE_VERSION="<not a Git checkout>"

if command -v git >/dev/null 2>&1 &&
   git -C "${REPO_ROOT}" rev-parse --is-inside-work-tree \
       >/dev/null 2>&1
then
    CODE_VERSION="$(
        git -C "${REPO_ROOT}" rev-parse --short HEAD
    )"

    if [[ -n "$(git -C "${REPO_ROOT}" status --porcelain)" ]]; then
        CODE_VERSION="${CODE_VERSION} +uncommitted-edits"
    fi
fi


# ============================================================
# Job information
# ============================================================

echo "============================================================"
echo "Y-only simulation analysis"
echo "============================================================"
echo "Date:                $(date)"
echo "Host:                $(hostname)"
echo "SLURM job ID:        ${SLURM_JOB_ID:-not running under SLURM}"
echo
echo "Repository root:     ${REPO_ROOT}"
echo "Code directory:      ${CODE_DIR}"
echo "Code version:        ${CODE_VERSION}"
echo
echo "R executable:        $(command -v R)"
echo "Rscript executable:  $(command -v Rscript)"
echo "R module:            ${R_MODULE:-not loaded by this script}"
echo "R_LIBS_USER:         ${R_LIBS_USER:-R default}"
echo
echo "Results directory:   ${FMM_OUT}"
echo "Expected replicates: ${FMM_REPS}"
echo "============================================================"
echo


# ============================================================
# Locate result files
# ============================================================

shopt -s nullglob
FILES=("${FMM_OUT}"/study1_saveK_results*.txt)
shopt -u nullglob

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "ERROR: No result files were found matching:"
    echo "       ${FMM_OUT}/study1_saveK_results*.txt"
    echo
    echo "The simulation may still be running, or FMM_OUT may be wrong."
    echo
    echo "The default five-replicate test directory is:"
    echo "       ${REPO_ROOT}/test_results_y"
    exit 1
fi

echo "Found ${#FILES[@]} result file(s):"

for result_file in "${FILES[@]}"; do
    echo "  $(basename "${result_file}")"
done

echo


# ============================================================
# 1. Integrity check
# ============================================================

echo "============================================================"
echo "1. Integrity check"
echo "============================================================"
echo

Rscript --vanilla \
    "${CODE_DIR}/check_results.R" \
    "${FMM_REPS}" \
    "${FILES[@]}"


# ============================================================
# 2. Summary in the manuscript Table 1 layout
# ============================================================

echo
echo "============================================================"
echo "2. Summary"
echo "============================================================"
echo

Rscript --vanilla \
    "${CODE_DIR}/summarize_results.R" \
    "${FMM_OUT}" \
    --reps "${FMM_REPS}"


# ============================================================
# 3. Figures
# ============================================================
#
# Plotting requires:
#
#   dplyr
#   tidyr
#   purrr
#   ggplot2
#   stringr
#
# If these are unavailable, the integrity check and numerical
# summary remain complete, and only plotting is skipped.
# ============================================================

echo
echo "============================================================"
echo "3. Figures"
echo "============================================================"
echo

if Rscript --vanilla - <<'RSCRIPT'
plot_packages <- c(
    "dplyr",
    "tidyr",
    "purrr",
    "ggplot2",
    "stringr"
)

available <- vapply(
    plot_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
)

missing_packages <- plot_packages[!available]

if (length(missing_packages) > 0L) {
    cat(
        "Missing plotting packages: ",
        paste(missing_packages, collapse = ", "),
        "\n",
        sep = ""
    )

    quit(status = 1L)
}

cat("All required plotting packages are available.\n")
RSCRIPT
then
    RESULTS_DIR="${FMM_OUT}" \
        Rscript --vanilla "${CODE_DIR}/plot_sim_y_paper.R"

    echo
    echo "Figures written under:"
    echo "  ${REPO_ROOT}/plots_metric_byN_M25_int_only"
else
    echo
    echo "SKIPPED: The figures were not created because one or more"
    echo "required plotting packages are unavailable."

    if [[ -f "${CODE_DIR}/install_packages.R" ]]; then
        echo
        echo "The packages may be installed using:"
        echo "  Rscript ${CODE_DIR}/install_packages.R"
    fi
fi


# ============================================================
# Completion information
# ============================================================

echo
echo "============================================================"
echo "Y-ONLY ANALYSIS COMPLETE"
echo "Results read from:   ${FMM_OUT}"
echo "Expected replicates: ${FMM_REPS}"
echo "Finished:            $(date)"
echo "============================================================"