#!/bin/bash

# ==============================================================================
# RUN EXPERIMENT 2 (MORE SELECTION) - TABLE 2 REPLICATION 
# ==============================================================================
# This script runs Experiment 2 (More sample selection, 25%) for all combinations of:
# - Models: A (static), B (dynamic) 
# - Autoregressive parameters: 0.25, 0.50, 0.75
# For a specified sample size N.
#
# Usage: bash run_experiment2.sh N
# Example: bash run_experiment2.sh 500
# ==============================================================================

# Check if N parameter is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide sample size N as argument"
    echo "Usage: bash run_experiment2.sh N"
    echo "Example: bash run_experiment2.sh 500"
    exit 1
fi

N=$1

# Change to the script directory
cd "$(dirname "$0")"

# Define Stata command (using the full path)
STATA="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE"

echo "========================================================================"
echo "STARTING EXPERIMENT 2: MORE SAMPLE SELECTION (25%) - TABLE 2 REPLICATION"
echo "========================================================================"
echo "Running all combinations for N=$N endogenous selection sensitivity analysis"
echo "Time: $(date)"
echo ""

# Determine output directory based on N
if [ $N -eq 500 ]; then
    OUTPUT_DIR="output/partial_tab2"
elif [ $N -eq 5000 ]; then
    OUTPUT_DIR="output/partial_tab3"
else
    echo "Warning: N=$N not recognized. Using default output/partial_tab2/"
    OUTPUT_DIR="output/partial_tab2"
fi

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Start all simulations in parallel
echo "Starting all Experiment 2 specifications in parallel for N=$N..."
echo "Output directory: $OUTPUT_DIR"
echo ""

# Run all 6 combinations in parallel
$STATA -b experiment2.do $N A 0.25 $OUTPUT_DIR &
$STATA -b experiment2.do $N A 0.50 $OUTPUT_DIR &
$STATA -b experiment2.do $N A 0.75 $OUTPUT_DIR &
$STATA -b experiment2.do $N B 0.25 $OUTPUT_DIR &
$STATA -b experiment2.do $N B 0.50 $OUTPUT_DIR &
$STATA -b experiment2.do $N B 0.75 $OUTPUT_DIR &

echo "All 6 specifications started in parallel for N=$N."
echo "Jobs running in background - check $OUTPUT_DIR/ for results."

echo "========================================================================"
echo "EXPERIMENT 2 STARTED"
echo "========================================================================"
echo "All specifications running in parallel for N=$N"
echo "Results will be saved to: $OUTPUT_DIR/"
echo "File pattern: exp2_model[A|B]_N$N_rho[0.25|0.50|0.75].*"
echo "Time: $(date)"
echo "========================================================================"
