#!/bin/bash

# ========================================================================
# PARALLEL EXECUTION SCRIPT FOR N=500 SIMULATIONS
# ========================================================================
# 
# This script runs all N=500 combinations in parallel to speed up execution.
# Each simulation runs as a background process.
#
# Usage: ./run_parallel_N500.sh
#
# Note: This will start 12 parallel processes. Monitor your system resources.
# ========================================================================

echo "========================================================================"
echo "STARTING PARALLEL N=500 SIMULATIONS"
echo "========================================================================"
echo "Starting 12 parallel processes at: $(date)"
echo ""

# Change to the script directory
cd "$(dirname "$0")"

# Define Stata command (using the alias path)
STATA="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE"

# Create output directory if it doesn't exist
mkdir -p output/partial

# Start all N=500 simulations in parallel
echo "Starting simulations..."

# Non-endogenous simulations
$STATA -b nonendogenous.do 500 A 0.25 &
$STATA -b nonendogenous.do 500 A 0.50 &
$STATA -b nonendogenous.do 500 A 0.75 &
$STATA -b nonendogenous.do 500 B 0.25 &
$STATA -b nonendogenous.do 500 B 0.50 &
$STATA -b nonendogenous.do 500 B 0.75 &

# Endogenous simulations
$STATA -b endogenous.do 500 A 0.25 &
$STATA -b endogenous.do 500 A 0.50 &
$STATA -b endogenous.do 500 A 0.75 &
$STATA -b endogenous.do 500 B 0.25 &
$STATA -b endogenous.do 500 B 0.50 &
$STATA -b endogenous.do 500 B 0.75 &

echo "All 12 jobs started in parallel."

echo ""
echo "========================================================================"
echo "ALL 12 PROCESSES STARTED"
echo "========================================================================"
echo "Use 'wait' to wait for all jobs to complete, or monitor manually."
echo "========================================================================"
