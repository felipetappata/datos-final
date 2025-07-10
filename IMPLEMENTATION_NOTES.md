# Monte Carlo Replication of Sadoon et al. (2019)

## Implementation Summary

This code provides an exact replication of the Monte Carlo experiments from "Simple methods for consistent estimation of dynamic panel data sample selection models" by Sadoon, Jiménez-Martín, and Labeaga (2019).

### Key Features Implemented

1. **Exact Data Generating Process (DGP)**
   - AR(1) dynamic panel model with sample selection
   - Static (Model A) and dynamic (Model B) selection equations
   - Correlated error components (α_i, ε_it) with selection components (η_i, u_it)
   - Correlation coefficient = 0.447 as specified in equations (46-49)

2. **Estimators**
   - **Arellano-Bond (AB)**: First-differenced GMM using `xtabond2 y L.y, gmm(L.y, lag(2 .)) nolevel`
   - **System GMM**: Combined system using `xtabond2 y L.y, gmm(L.y, lag(2 .)) iv(D.L.y, equation(level))`

3. **Sample Selection**
   - 85% selection probability (15% attrition)
   - Requires 3 consecutive observations for GMM estimation
   - Substantial data loss as documented in the paper

4. **Experimental Design**
   - **Table 1**: Basic results for ρ ∈ {0.25, 0.50, 0.75}, N ∈ {500, 5000}
   - **Table 2**: Sensitivity analysis for N=500
   - **Table 3**: Sensitivity analysis for N=5000
   - 500 Monte Carlo replications each (10 in test mode)

5. **Sensitivity Experiments**
   - Short T (T=4 instead of T=7)
   - Higher selection (25% instead of 15%)
   - Different variance ratios (σ_η/σ_ε = 2)
   - Lower correlation (0.25 instead of 0.5)
   - Non-stationary time-varying components

### File Structure

```
code/
├── main.do                 # Main implementation file
├── input/                  # Input data directory (empty for simulation)
└── output/
    ├── tables/            # Final summary tables
    └── partial/           # Partial results for resuming
```

### Usage Instructions

#### Quick Test (Recommended First)
```stata
do main.do test
```
This runs Table 1 with small samples (N=100, 10 replications) to verify everything works.

#### Full Replication
```stata
do main.do full            # All tables (very time-consuming)
do main.do table1          # Just Table 1 (basic results)
do main.do table2          # Just Table 2 (sensitivity, N=500)
do main.do table3          # Just Table 3 (sensitivity, N=5000)
```

#### Programmatic Usage
```stata
run_monte_carlo, test_mode table(1)
run_monte_carlo, table(1 2 3)
```

### Expected Runtime

- **Test mode**: ~5-10 minutes
- **Table 1 full**: ~2-4 hours
- **Tables 2-3**: ~4-6 hours each
- **All tables**: ~12-16 hours

### Output Files

- `table1_basic_summary.csv`: Main results replicating Table 1
- `table2_sensitivity_summary.csv`: Sensitivity analysis for N=500
- `table3_sensitivity_summary.csv`: Sensitivity analysis for N=5000
- Progress logs and partial results for resuming interrupted runs

### Key Implementation Notes

1. **Exact Specification**: All parameters match the paper exactly
2. **State Management**: Can resume from interruptions
3. **Progress Tracking**: Console output and log files
4. **Error Handling**: Robust to estimation failures
5. **Test Mode**: Quick verification with small samples

### Verification

The implementation should reproduce the key findings from the paper:
- AB estimator is consistent regardless of selection endogeneity
- System GMM has small bias due to invalid moment conditions
- Bias is mainly from correlation between time-invariant components
- System GMM performs better in small samples despite bias

### Dependencies

- Stata with `xtabond2` package (auto-installed if missing)
- Sufficient disk space for large intermediate files
- Adequate RAM for large panel datasets (especially N=5000 case)
