#!/bin/bash
#SBATCH --job-name=MixXY_test
#SBATCH --partition=statistics_q
#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=6-00:00:00
#SBATCH --output=logs/simxy_MixXY_test_%A_%a.out
#SBATCH --error=logs/simxy_MixXY_test_%A_%a.err

set -euo pipefail

# ============================================================
# Submission instructions
# ============================================================
#
# Submit from the root of the cloned repository:
#
#   cd /path/to/fmm_manuscript
#   mkdir -p logs
#   sbatch code/run_sim_xy.sh
#
# To specify a particular R library:
#
#   sbatch --export=ALL,FMM_R_LIB=/path/to/R/library \
#       code/run_sim_xy.sh
#
# Expected repository structure:
#
#   fmm_manuscript/
#   ├── code/
#   │   ├── run_sim_xy.sh
#   │   ├── sim_study_xy.R
#   │   └── sim_xy_functions.R
#   ├── logs/
#   └── test_results_xy/
#
# ============================================================


# ============================================================
# Repository locations
# ============================================================

# SLURM_SUBMIT_DIR is the directory from which sbatch was run.
# Users should therefore submit this job from the repository root.
#
# FMM_REPO_ROOT can be supplied to override this location.

REPO_ROOT="${FMM_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CODE_DIR="${REPO_ROOT}/code"

SIM_SCRIPT="${CODE_DIR}/sim_study_xy.R"
export FMM_SRC="${CODE_DIR}/sim_xy_functions.R"

# Run from the repository root so relative paths used by the R
# scripts behave consistently.
cd "${REPO_ROOT}"

if [[ ! -f "${SIM_SCRIPT}" ]]; then
    echo "ERROR: Could not find the simulation script:"
    echo "       ${SIM_SCRIPT}"
    echo
    echo "Submit from the repository root or set FMM_REPO_ROOT."
    exit 1
fi

if [[ ! -f "${FMM_SRC}" ]]; then
    echo "ERROR: Could not find the functions file:"
    echo "       ${FMM_SRC}"
    echo
    echo "Submit from the repository root or set FMM_REPO_ROOT."
    exit 1
fi


# ============================================================
# Simulation design
# ============================================================
#
# One SLURM array task is used for each combination of:
#
#   N = 10, 500, 1000, 10000
#   cluster separation = 1, 2, 3
#
# Task mapping:
#
#   task  1: N=10,    separation=1
#   task  2: N=10,    separation=2
#   task  3: N=10,    separation=3
#
#   task  4: N=500,   separation=1
#   task  5: N=500,   separation=2
#   task  6: N=500,   separation=3
#
#   task  7: N=1000,  separation=1
#   task  8: N=1000,  separation=2
#   task  9: N=1000,  separation=3
#
#   task 10: N=10000, separation=1
#   task 11: N=10000, separation=2
#   task 12: N=10000, separation=3
#
# Five replicate datasets are run within each task by default.
# ============================================================

N_VALUES=(10 500 1000 10000)

TASK_ID="${SLURM_ARRAY_TASK_ID:-${FMM_TASK_ID:-}}"

if [[ -z "${TASK_ID}" ]]; then
    echo "ERROR: No array task ID was provided."
    echo
    echo "Submit using sbatch or set FMM_TASK_ID for a direct test."
    exit 1
fi

if ! [[ "${TASK_ID}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Task ID must be an integer from 1 through 12."
    echo "Received: ${TASK_ID}"
    exit 1
fi

if (( TASK_ID < 1 || TASK_ID > 12 )); then
    echo "ERROR: Task ID must be between 1 and 12."
    echo "Received: ${TASK_ID}"
    exit 1
fi

# Convert task 1-12 into:
#
#   N_INDEX  = 0-3
#   CLUST_SEP = 1-3

N_INDEX=$(( (TASK_ID - 1) / 3 ))
CLUST_SEP=$(( (TASK_ID - 1) % 3 + 1 ))

FMM_N="${N_VALUES[$N_INDEX]}"


# ============================================================
# Simulation settings
# ============================================================
#
# These defaults can be overridden without editing this file.
#
# For example:
#
#   sbatch --export=ALL,FMM_REPS=50,PRIOR_M=25 \
#       code/run_sim_xy.sh
#
# ============================================================

FMM_REPS="${FMM_REPS:-5}"
K_TRUE="${K_TRUE:-3}"
SIG_CONST="${SIG_CONST:-0.25}"
PRIOR_M="${PRIOR_M:-2}"
CORES="${CORES:-${SLURM_CPUS_PER_TASK:-1}}"

# Keep the five-replicate test separate from final results.
export FMM_OUT="${FMM_OUT:-${REPO_ROOT}/test_results_xy}"

mkdir -p "${FMM_OUT}"


# ============================================================
# Unique output identifier
# ============================================================

ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID:-${TASK_ID}}"

export FMM_TAG="N${FMM_N}_sep${CLUST_SEP}_${ARRAY_JOB_ID}_${ARRAY_TASK_ID}"


# ============================================================
# R environment
# ============================================================
#
# R_MODULE defaults to R/4.4.0. Users on another cluster can
# override it:
#
#   sbatch --export=ALL,R_MODULE=R/4.5.0 \
#       code/run_sim_xy.sh
#
# To specify a clean personal R library:
#
#   sbatch --export=ALL,FMM_R_LIB=/path/to/R/library \
#       code/run_sim_xy.sh
#
# When FMM_R_LIB is supplied, R_LIBS is cleared to prevent an
# inherited incompatible library from entering .libPaths().
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

# Prevent oversubscription within each one-core SLURM task.
export OMP_NUM_THREADS="${CORES}"
export MKL_NUM_THREADS="${CORES}"
export OPENBLAS_NUM_THREADS="${CORES}"
export BLIS_NUM_THREADS="${CORES}"


# ============================================================
# Job information
# ============================================================

echo "============================================================"
echo "XY Gaussian mixture regression simulation test"
echo "============================================================"
echo "Date:                $(date)"
echo "Host:                $(hostname)"
echo "Repository root:     ${REPO_ROOT}"
echo "Code directory:      ${CODE_DIR}"
echo
echo "SLURM job ID:        ${SLURM_JOB_ID:-not running under SLURM}"
echo "SLURM array job ID:  ${SLURM_ARRAY_JOB_ID:-not available}"
echo "SLURM array task:    ${TASK_ID}"
echo
echo "R executable:        $(command -v R)"
echo "Rscript executable:  $(command -v Rscript)"
echo "R module:            ${R_MODULE:-not loaded by this script}"
echo "R_LIBS_USER:         ${R_LIBS_USER:-R default}"
echo
echo "Simulation script:   ${SIM_SCRIPT}"
echo "Functions file:      ${FMM_SRC}"
echo "Output directory:    ${FMM_OUT}"
echo "Output tag:          ${FMM_TAG}"
echo
echo "N:                   ${FMM_N}"
echo "K:                   ${K_TRUE}"
echo "Cluster separation:  ${CLUST_SEP}"
echo "Sigma constant:      ${SIG_CONST}"
echo "Prior M option:      ${PRIOR_M}"
echo "Replicates:          ${FMM_REPS}"
echo "Cores:               ${CORES}"
echo "============================================================"
echo


# ============================================================
# R package and environment check
# ============================================================

Rscript --vanilla - <<'RSCRIPT'

cat("\n================ R ENVIRONMENT CHECK ================\n")

cat("\nR version:\n")
cat(R.version.string, "\n")

cat("\nR_LIBS_USER:\n")
cat(Sys.getenv("R_LIBS_USER", unset = "<not set>"), "\n")

cat("\nR library paths:\n")
print(.libPaths())

required_packages <- c(
    "Matrix",
    "quantreg",
    "salso",
    "MCMCpack",
    "coda",
    "mvtnorm",
    "miscPack",
    "mclust",
    "kernlab",
    "cluster",
    "dbscan",
    "e1071",
    "matrixStats"
)

cat("\n================ PACKAGE CHECK ======================\n\n")

for (package_name in required_packages) {

    if (!requireNamespace(package_name, quietly = TRUE)) {
        stop(
            "Required R package is not installed: ",
            package_name,
            "\nInstall the required dependencies before running ",
            "the simulation.",
            call. = FALSE
        )
    }

    cat(
        sprintf(
            "%-12s version %-12s %s\n",
            package_name,
            as.character(packageVersion(package_name)),
            find.package(package_name)
        )
    )
}

cat("\n================ MISC PACK CHECK ====================\n\n")

if (!exists(
    "gaussian_mixture",
    where = asNamespace("miscPack"),
    inherits = FALSE
)) {
    stop(
        "miscPack::gaussian_mixture() was not found. ",
        "Confirm that the required version of miscPack is installed.",
        call. = FALSE
    )
}

cat("miscPack::gaussian_mixture() found successfully.\n")

cat("\n================ PACKAGE LOAD CHECK =================\n\n")

for (package_name in required_packages) {

    suppressPackageStartupMessages(
        library(package_name, character.only = TRUE)
    )

    cat("Loaded: ", package_name, "\n", sep = "")
}

cat("\n*****************************************************\n")
cat("PACKAGE CHECK PASSED\n")
cat("*****************************************************\n\n")

RSCRIPT


# ============================================================
# Run XY simulation
# ============================================================
#
# sim_study_xy.R arguments:
#
#   1. N
#   2. k
#   3. clust_sep
#   4. sig_const
#   5. prior_M
#   6. reps
#   7. cores
#
# ============================================================

echo
echo "============================================================"
echo "Starting XY simulation"
echo "============================================================"
echo

Rscript --vanilla "${SIM_SCRIPT}" \
    "${FMM_N}" \
    "${K_TRUE}" \
    "${CLUST_SEP}" \
    "${SIG_CONST}" \
    "${PRIOR_M}" \
    "${FMM_REPS}" \
    "${CORES}"

echo
echo "============================================================"
echo "XY SIMULATION COMPLETED SUCCESSFULLY"
echo "N:                   ${FMM_N}"
echo "Cluster separation:  ${CLUST_SEP}"
echo "Replicates:          ${FMM_REPS}"
echo "Results:             ${FMM_OUT}"
echo "Finished:            $(date)"
echo "============================================================"