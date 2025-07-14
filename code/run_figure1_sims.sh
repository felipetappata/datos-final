#!/bin/bash

# ========================================================================
# PARALLEL EXECUTION SCRIPT FOR FIGURE 1 SIMULATIONS
# ========================================================================
# 
# This script runs simulations for Figure 1, varying N and rho values.
# For each (N, rho, model) combination, it runs endogenous2.do which 
# generates both endogenous selection estimates AND full sample estimates
# (corresponding to "AB all" and "system all" from Figure 1).
#
# Usage: bash run_figure1_sims.sh "N_values" "rho_values" "models"
#
# Example: bash run_figure1_sims.sh "200 400" "0.25 0.50" "A B"
#
# Output files per combination:
#   - endo_modelA_N200_rho0.25.csv (endogenous selection estimates)
#   - fullsample_modelA_N200_rho0.25.csv (full sample estimates - "AB all", "system all")
#
# Note: Each simulation runs as a background process. Monitor your system 
# resources and adjust the number of parallel processes accordingly.
# ========================================================================

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    echo "Error: Incorrect number of arguments"
    echo "Usage: bash run_figure1_sims.sh \"N_values\" \"rho_values\" \"models\""
    echo "Example: bash run_figure1_sims.sh \"200 400\" \"0.25 0.50\" \"A B\""
    exit 1
fi

# Parse command line arguments
N_VALUES=($1)
RHO_VALUES=($2)
MODELS=($3)

echo "========================================================================"
echo "STARTING FIGURE 1 SIMULATIONS"
echo "========================================================================"
echo "Starting simulations at: $(date)"
echo "N values: ${N_VALUES[@]}"
echo "Rho values: ${RHO_VALUES[@]}"
echo "Models: ${MODELS[@]}"
echo ""

# Change to the script directory
cd "$(dirname "$0")"

# Define Stata command (using the alias path)
STATA="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE"

# Create output directory if it doesn't exist
mkdir -p output/partial

# Counter for tracking number of jobs
job_count=0

# Start simulations for each combination
for N in "${N_VALUES[@]}"; do
    for rho in "${RHO_VALUES[@]}"; do
        for model in "${MODELS[@]}"; do
            echo "Starting simulation for N=$N, rho=$rho, model=$model"
            
            # Run endogenous2.do which generates both endogenous and full sample estimates
            $STATA -b endogenous2.do $N $model $rho &
            job_count=$((job_count + 1))
        done
    done
done

echo ""
echo "========================================================================"
echo "ALL $job_count PROCESSES STARTED"
echo "========================================================================"
echo "Expected output files will be in output/partial/ with names like:"
echo "  endo_modelA_N200_rho0.25.csv (endogenous selection estimates)"
echo "  fullsample_modelA_N200_rho0.25.csv (full sample estimates - AB all, system all)"
echo ""
echo "Monitor progress manually or use system monitoring tools."
echo "========================================================================"
