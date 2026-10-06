#!/bin/bash
#SBATCH --job-name=MixY_stats_test
#SBATCH --array=7-12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p statistics_q
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/simy_stats_test_%A_%a.out
#SBATCH --error=logs/simy_stats_test_%A_%a.err

set -euo pipefail

# ============================================================
# One-replicate Y-only timing test on the Statistics partition
# ============================================================
#
# Submit this file from the root of the cloned repository:
#
#   cd /path/to/constrained_gmm
#   mkdir -p logs
#   sbatch code/run_sim_y_statistics_test.sh
#
# This batch script uses the same Statistics partition and
# resource-request style as the existing Statistics jobs:
#
#   partition: statistics_q
#   tasks:     1
#   CPUs:      1
#   memory:    128 GB
#   wall time: 4 days
#
# The script submits array tasks 7-12 from run_sim_y.sh:
#
#   task  7: N=1000,  separation=1
#   task  8: N=1000,  separation=2
#   task  9: N=1000,  separation=3
#   task 10: N=10000, separation=1
#   task 11: N=10000, separation=2
#   task 12: N=10000, separation=3
#
# FMM_REPS defaults to 1, so this runs one replicate for each
# of the six N/separation combinations: six data sets in total.
#
# This is a small launcher for run_sim_y.sh. The main script
# continues to perform the package checks, map array task IDs to
# N and separation, and call sim_study_y.R. Keeping that logic in
# one file prevents the two batch scripts from drifting apart.
# ============================================================


# ============================================================
# Repository locations
# ============================================================

# SLURM_SUBMIT_DIR is the directory from which sbatch was run.
# Submit from the repository root as shown above. Alternatively,
# set FMM_REPO_ROOT when submitting the job.
REPO_ROOT="${FMM_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
CODE_DIR="${REPO_ROOT}/code"
MAIN_JOB_SCRIPT="${CODE_DIR}/run_sim_y.sh"

if [[ ! -f "${MAIN_JOB_SCRIPT}" ]]; then
    echo "ERROR: Could not find the main Y-only batch script:"
    echo "       ${MAIN_JOB_SCRIPT}"
    echo
    echo "Submit from the repository root or set FMM_REPO_ROOT."
    exit 1
fi


# ============================================================
# Test settings
# ============================================================

# Use one replicate per selected array task unless explicitly
# overridden at submission.
export FMM_REPS="${FMM_REPS:-1}"

# Keep these timing-test results separate from the earlier
# genacc_q results and the final publication simulation results.
export FMM_OUT="${FMM_OUT:-${REPO_ROOT}/test_results_y_statistics}"

# Make the repository location explicit for run_sim_y.sh.
export FMM_REPO_ROOT="${REPO_ROOT}"

mkdir -p "${FMM_OUT}"


# ============================================================
# Job information
# ============================================================

echo "============================================================"
echo "Y-only timing test on the Statistics partition"
echo "============================================================"
echo "Date:                $(date)"
echo "Host:                $(hostname)"
echo "SLURM job ID:        ${SLURM_JOB_ID}"
echo "SLURM array job ID:  ${SLURM_ARRAY_JOB_ID}"
echo "SLURM array task:    ${SLURM_ARRAY_TASK_ID}"
echo "SLURM partition:     ${SLURM_JOB_PARTITION:-statistics_q}"
echo "Repository root:     ${REPO_ROOT}"
echo "Main batch script:   ${MAIN_JOB_SCRIPT}"
echo "Results directory:   ${FMM_OUT}"
echo "Replicates:          ${FMM_REPS}"
echo "============================================================"
echo


# ============================================================
# Run the existing Y-only simulation workflow
# ============================================================

# The #SBATCH directives inside run_sim_y.sh are comments when
# it is invoked with bash here. The current Statistics allocation
# and array task ID remain in effect.
exec bash "${MAIN_JOB_SCRIPT}"

