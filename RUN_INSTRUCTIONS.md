# How to Run the Monte Carlo Experiment

## Prerequisites
- Stata software (version 14 or later recommended)
- The `xtabond2` package (will be auto-installed)

## Step-by-Step Instructions

### 1. Navigate to the Code Directory
```
cd "your_path/final/code"
```

### 2. Start with Test Mode (Recommended)
Open Stata and run:
```stata
do main.do test
```

This will:
- Run a quick test with small samples (N=100, 10 replications)
- Take approximately 5-10 minutes
- Verify that all code works correctly
- Generate test results in `output/`

### 3. Check Test Results
After the test completes, check:
- Console output for any errors
- `output/table1_basic_summary.csv` for results summary
- `output/partial/progress_log.txt` for execution log

### 4. Run Full Experiments
If the test works correctly, proceed with full experiments:

#### Option A: Run All Tables (12-16 hours)
```stata
do main.do full
```

#### Option B: Run Individual Tables
```stata
do main.do table1    // Table 1: Basic results (2-4 hours)
do main.do table2    // Table 2: Sensitivity N=500 (4-6 hours)  
do main.do table3    // Table 3: Sensitivity N=5000 (4-6 hours)
```

### 5. Resume Interrupted Runs
If a run is interrupted, you can resume from the last completed replication:
```stata
run_monte_carlo, table(1) resume_rep(150)
```

### 6. Monitor Progress
- Watch console output for progress updates
- Check `output/partial/progress_log.txt` for detailed logs
- Partial results are saved continuously in `output/partial/`

## Expected Output

### Table 1 Results
The implementation should produce results similar to Table 1 in the paper, showing:
- AB estimator bias decreases with sample size
- System GMM has small, stable bias (1-2.5%)
- Bias patterns match theoretical predictions

### Sensitivity Analysis
Tables 2 and 3 will show how results change under:
- Very short panels (T=4)
- Higher selection rates (25%)
- Different variance ratios
- Lower error correlations
- Non-stationary components

## Troubleshooting

### Common Issues

1. **Missing xtabond2 package**: The code will auto-install it
2. **Memory issues**: Reduce sample sizes in test mode first
3. **Long runtimes**: This is expected; use test mode for verification
4. **Convergence failures**: Some replications may fail; this is normal

### Getting Help

If you encounter issues:
1. Check the console output for specific error messages
2. Review `output/partial/progress_log.txt`
3. Try running in test mode first
4. Check that you have sufficient disk space and memory

## Validation

Compare your results with Table 1 from the paper:
- AB bias should be larger for small N, smaller for large N
- System GMM should show consistent small bias
- Standard errors should be reasonable
- Patterns should match across ρ values (0.25, 0.50, 0.75)

The exact numbers may differ due to randomness, but patterns should be consistent.
