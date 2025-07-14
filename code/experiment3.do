/* ========================================================================= 
   EXPERIMENT 3: INCREASING RATIO OF VARIANCES - TABLE 2 REPLICATION
   =========================================================================
   
   This script implements Experiment III from Table 2 of Al Sadoon et al. (2019):
   "Increasing the ratio of variances: σ_η/σ_ε = 2"
   
   This experiment tests the sensitivity of AB and System GMM estimators when
   the ratio of the variance of the individual heterogeneous component to the 
   variance of the time-variant component of the outcome equation is increased
   to σ_α²/σ_ε² = 2. This is achieved by setting σ_η = 2 while keeping σ_ε = 1.
   
   Command line usage:
   stata -b do experiment3.do N model rho
   
   Where:
   - N: Sample size (any integer >= 200)
   - model: Selection model (A=static, B=dynamic)  
   - rho: Autoregressive parameter (0.25, 0.50, 0.75)
   
   Example: stata -b do experiment3.do 500 A 0.25
   
   Author: Generated for Al Sadoon et al. (2019) replication study
   Date: January 2025
   ========================================================================= */

clear all
set more off
set seed 08869  // As specified in instructions

// Optimize Mata for speed (helps with xtabond2 performance)
mata: mata set matafavor speed, perm

// Parse command line arguments
if "`1'" == "" | "`2'" == "" | "`3'" == "" {
    display as error "Usage: stata -b do experiment3.do N model rho [output_dir]"
    display as error "Where N>=200, model={A,B}, rho={0.25,0.50,0.75}"
    display as error "output_dir is optional (default: output/partial_tab2)"
    exit 198
}

local N `1'
local model `2'  
local rho `3'
local output_dir `4'

// Set default output directory if not provided
if "`output_dir'" == "" {
    local output_dir "output/partial_tab2"
}

// Validate inputs
if `N' < 200 {
    display as error "N must be >= 200"
    exit 198
}
if !inlist("`model'", "A", "B") {
    display as error "Model must be A (static) or B (dynamic)"
    exit 198
}
if !inlist(`rho', 0.25, 0.50, 0.75) {
    display as error "rho must be 0.25, 0.50, or 0.75"
    exit 198
}

// ========================================================================
// SIMULATION PARAMETERS - EXPERIMENT 3: INCREASING RATIO OF VARIANCES
// ========================================================================

// Time dimension (baseline)
local T 7                    // Time periods for estimation (after discarding 13)
local T_total 20             // Total periods generated (discarding first 13)
local T_discard 13           // Periods to discard for initial conditions
local reps 500               // Number of Monte Carlo replications

// Selection equation parameters (Modified for Experiment III)
local a_param 1.8            // Set so P(d_it* > 0) = 0.85 (15% selection) - but will change with σ_η=2
local sigma_z 1              // Standard deviation of z_it ~ N(0,1)
local sigma_eta 2            // CHANGED: Standard deviation of eta_i ~ N(0,2) - following literal Table 2 description
local sigma_u 1              // Standard deviation of u_it ~ N(0,1)

// EXPERIMENT 3 MODIFICATION: Following Table 2 literal description σ_η/σ_ε = 2
local sigma_alpha0 1         // Standard deviation of alpha_i^0 ~ N(0,1) (baseline)
local sigma_eps0 1           // Standard deviation of epsilon_it^0 ~ N(0,1) (baseline)
local theta_param 0.5        // Correlation parameter alpha_i = alpha_i^0 + theta*eta_i
local vartheta_param 0.5     // Correlation parameter eps_it = eps_it^0 + vartheta*u_it

// Note: Following Table 2 description σ_η/σ_ε = 2:
// With σ_η=2, σ_ε=1, we get the ratio σ_η/σ_ε = 2/1 = 2 ✓
// We should note that this changes selection probabilities slightly.

// ========================================================================
// MONTE CARLO SIMULATION
// ========================================================================

display as text ""
display as text "========================================================================="
display as text "EXPERIMENT 3: INCREASING RATIO OF VARIANCES - ENDOGENOUS SELECTION"
display as text "========================================================================="
display as text "Parameters: N=`N', Model=`model', rho=`rho', Replications=`reps'"
display as text "EXPERIMENT: Following Table 2 description σ_η/σ_ε = 2"
display as text "Selection: ~23% attrition (changed from 15%), ENDOGENOUS (theta=`theta_param', vartheta=`vartheta_param')"
display as text "Variance parameters: σ_η=`sigma_eta', σ_ε=`sigma_eps0' (primitive errors)"
local implied_corr = `theta_param' / sqrt(1 + `theta_param'^2)
display as text "Implied correlation: " %5.3f `implied_corr'
display as text "========================================================================="
display as text ""

// Initialize results storage
matrix AB_coefs = J(`reps', 1, .)
matrix AB_ses = J(`reps', 1, .)
matrix SYS_coefs = J(`reps', 1, .)
matrix SYS_ses = J(`reps', 1, .)

local ab_valid = 0
local sys_valid = 0
local ab_converged = 0
local sys_converged = 0

// Start replication loop
forvalues rep = 1/`reps' {
    
    if mod(`rep', 50) == 0 {
        display as text "... Running replication `rep'/`reps' ..."
    }
    
    // ====================================================================
    // DATA GENERATION PROCESS (Modified variance structure)
    // ====================================================================
    
    quietly {
        clear
        local total_obs = `N' * `T_total'
        set obs `total_obs'
        
        // Create panel structure
        gen long id = floor((_n-1)/`T_total') + 1
        gen int t = mod(_n-1, `T_total') + 1
        xtset id t
        
        // Generate individual-specific components (time-invariant)
        gen double alpha_i0 = rnormal(0, `sigma_alpha0')
        gen double eta_i = rnormal(0, `sigma_eta')  // σ_η = 2 following Table 2 description
        sort id
        by id: replace alpha_i0 = alpha_i0[1]
        by id: replace eta_i = eta_i[1]
        
        // Endogenous selection: correlation between equations
        gen double alpha_i = alpha_i0 + `theta_param' * eta_i
        
        // Generate time-varying components
        gen double eps_i0 = rnormal(0, `sigma_eps0')
        gen double u_it = rnormal(0, `sigma_u')
        gen double z_it = rnormal(0, `sigma_z')
        
        // Endogenous selection: correlation between equations
        gen double eps_it = eps_i0 + `vartheta_param' * u_it
        
        // ================================================================
        // SELECTION EQUATION (Static Model A or Dynamic Model B)
        // ================================================================
        
        if "`model'" == "A" {
            // Static selection: d_it* = a - z_it - eta_i - u_it
            gen double d_star = `a_param' - z_it - eta_i - u_it
            gen byte d = (d_star > 0)
        }
        else if "`model'" == "B" {
            // Dynamic selection: d_it* = a - 0.5*d_it-1 + z_it - eta_i - u_it
            // Need to handle sequential dependency properly with explicit loop
            gen double d_star = .
            gen byte d = .
            gen byte ditm1 = 1  // Initial value for d_{t-1} to start the chain
            
            sort id t
            forvalues tt = 1/`T_total' {
                replace ditm1 = d[_n-1] if t == `tt' & t > 1
                replace d_star = `a_param' - 0.5*ditm1 + z_it - eta_i - u_it if t == `tt'
                replace d = (d_star > 0) if t == `tt'
            }
        }
        
        // ================================================================
        // OUTCOME EQUATION (AR(1) Process)
        // ================================================================
        
        gen double y_star = .
        
        // Initial condition: y_i1* = (2 + alpha_i + eps_i1)/(1-rho) [stationary]
        gen double y_initial = (2 + alpha_i + eps_it) / (1 - `rho') 
        replace y_star = y_initial if t == 1
        
        // AR(1) process: y_it* = 2 + rho*y_it-1* + alpha_i + eps_it
        sort id t
        by id (t): replace y_star = 2 + `rho' * L.y_star + alpha_i + eps_it if t > 1
        
        // Observe y only when selected
        gen double y = y_star if d == 1
        
        // Drop burn-in periods to minimize initial condition problems
        drop if t <= `T_discard'
        
        // Reset panel structure after dropping observations
        xtset id t
        
        // Check if we have enough observations for estimation
        count if !missing(y)
        local n_obs_y = r(N)
        
        // Count observations with 3 consecutive periods (needed for AB)
        gen byte has_3consec = 0
        sort id t
        by id: replace has_3consec = 1 if !missing(y) & !missing(L.y) & !missing(L2.y)
        count if has_3consec == 1
        local n_obs_ab = r(N)
        
        // Count observations with 2 consecutive periods (needed for System levels)
        gen byte has_2consec = 0
        by id: replace has_2consec = 1 if !missing(y) & !missing(L.y)
        count if has_2consec == 1
        local n_obs_sys = r(N)
    }
    
    // ====================================================================
    // ESTIMATION
    // ====================================================================
    
    // ARELLANO-BOND (AB) FIRST-DIFFERENCE GMM ESTIMATOR
    capture noisily xtabond2 y L.y, gmm(L.y, collapse) nolevel
    
    if _rc == 0 {
        local ab_converged = `ab_converged' + 1
        matrix b_ab = e(b)
        matrix V_ab = e(V)
        
        // Extract coefficient on L.y
        local pos_ab = colnumb(b_ab, "L.y")
        if !missing(`pos_ab') & `pos_ab' > 0 {
            scalar ab_coef = b_ab[1, `pos_ab']
            scalar ab_se = sqrt(V_ab[`pos_ab', `pos_ab'])
            
            // Store results if valid
            if !missing(ab_coef) & !missing(ab_se) & ab_se > 0 {
                matrix AB_coefs[`rep', 1] = ab_coef
                matrix AB_ses[`rep', 1] = ab_se
                local ab_valid = `ab_valid' + 1
            }
        }
    }
    
    // SYSTEM GMM ESTIMATOR (AB + Levels)
    capture noisily xtabond2 y L.y, gmm(L.y, lag(2 .)) iv(L.D.y, equation(level))
    
    if _rc == 0 {
        local sys_converged = `sys_converged' + 1
        matrix b_sys = e(b)
        matrix V_sys = e(V)
        
        // Extract coefficient on L.y
        local pos_sys = colnumb(b_sys, "L.y")
        if !missing(`pos_sys') & `pos_sys' > 0 {
            scalar sys_coef = b_sys[1, `pos_sys']
            scalar sys_se = sqrt(V_sys[`pos_sys', `pos_sys'])
            
            // Store results if valid
            if !missing(sys_coef) & !missing(sys_se) & sys_se > 0 {
                matrix SYS_coefs[`rep', 1] = sys_coef
                matrix SYS_ses[`rep', 1] = sys_se
                local sys_valid = `sys_valid' + 1
            }
        }
    }
    
    // Progress indicator for large replications
    if mod(`rep', 100) == 0 & `rep' > 0 {
        display as text "    Replication `rep': AB valid=`ab_valid', SYS valid=`sys_valid'"
    }
}

// ========================================================================
// COMPUTE AND SAVE RESULTS
// ========================================================================

display as text ""
display as text "========================================================================="
display as text "EXPERIMENT 3 SIMULATION COMPLETED"
display as text "========================================================================="
display as text "AB estimations: `ab_converged' converged, `ab_valid' valid results"
display as text "SYS estimations: `sys_converged' converged, `sys_valid' valid results"
display as text ""

// Convert matrices to variables for analysis
clear
svmat AB_coefs, names(ab_coef)
svmat AB_ses, names(ab_se)
svmat SYS_coefs, names(sys_coef)
svmat SYS_ses, names(sys_se)

// Compute summary statistics
if `ab_valid' >= 10 {
    summarize ab_coef1, detail
    scalar ab_mean = r(mean)
    scalar ab_sd = r(sd)
    scalar ab_bias = ab_mean - `rho'
    
    display as text "AB ESTIMATOR RESULTS:"
    display as text "  Valid replications: `ab_valid'"
    display as text "  Mean coefficient: " %9.6f ab_mean
    display as text "  Standard deviation: " %9.6f ab_sd
    display as text "  Bias (mean - `rho'): " %9.6f ab_bias
}
else {
    scalar ab_mean = .
    scalar ab_sd = .
    scalar ab_bias = .
    display as error "WARNING: Insufficient valid AB results (`ab_valid' < 10)"
}

if `sys_valid' >= 10 {
    summarize sys_coef1, detail
    scalar sys_mean = r(mean)
    scalar sys_sd = r(sd) 
    scalar sys_bias = sys_mean - `rho'
    
    display as text ""
    display as text "SYSTEM GMM ESTIMATOR RESULTS:"
    display as text "  Valid replications: `sys_valid'"
    display as text "  Mean coefficient: " %9.6f sys_mean
    display as text "  Standard deviation: " %9.6f sys_sd
    display as text "  Bias (mean - `rho'): " %9.6f sys_bias
}
else {
    scalar sys_mean = .
    scalar sys_sd = .
    scalar sys_bias = .
    display as error "WARNING: Insufficient valid SYS results (`sys_valid' < 10)"
}

// ========================================================================
// SAVE RESULTS TO FILE
// ========================================================================

// Create results matrix for export
matrix Results = J(1, 9, .)
matrix colnames Results = N rho AB_bias AB_se SYS_bias SYS_se AB_valid SYS_valid experiment

matrix Results[1,1] = `N'
matrix Results[1,2] = `rho'
matrix Results[1,3] = ab_bias
matrix Results[1,4] = ab_sd
matrix Results[1,5] = sys_bias
matrix Results[1,6] = sys_sd
matrix Results[1,7] = `ab_valid'
matrix Results[1,8] = `sys_valid'
matrix Results[1,9] = 3

// Display final results table
display as text ""
display as text "========================================================================="
display as text "EXPERIMENT 3 FINAL RESULTS - INCREASING RATIO OF VARIANCES"
display as text "Model `model', N=`N', rho=`rho'"
display as text "========================================================================="
matrix list Results, format(%9.6f)

// Export results to text file for easier reading
local txt_filename "`output_dir'/exp3_model`model'_N`N'_rho`rho'.txt"
file open txtfile using "`txt_filename'", write replace
file write txtfile "=========================================================================" _n
file write txtfile "EXPERIMENT 3: INCREASING RATIO OF VARIANCES - TABLE 2 REPLICATION" _n
file write txtfile "AL SADOON ET AL. (2019) - ENDOGENOUS SELECTION RESULTS" _n
file write txtfile "=========================================================================" _n
file write txtfile "Experiment: Increasing ratio of variances (σ_α²/σ_ε² = 2)" _n
file write txtfile "Model: `model' (Static=A, Dynamic=B)" _n
file write txtfile "Sample Size (N): `N'" _n  
file write txtfile "Autoregressive Parameter (rho): `rho'" _n
file write txtfile "Number of Replications: `reps'" _n
file write txtfile "Selection Type: Endogenous (theta=0.5, vartheta=0.5)" _n
file write txtfile "Variance Parameters: σ_α0=`sigma_alpha0', σ_ε0=`sigma_eps0', σ_η=`sigma_eta'" _n
file write txtfile "Date: `c(current_date)' `c(current_time)'" _n
file write txtfile "=========================================================================" _n
file write txtfile "" _n
file write txtfile "ESTIMATION RESULTS:" _n
file write txtfile "-------------------" _n
file write txtfile "AB Estimator:" _n
file write txtfile "  Bias (mean - true): " %9.6f (ab_bias) _n
file write txtfile "  Standard Error: " %9.6f (ab_sd) _n
file write txtfile "  Valid Replications: " %9.0f (`ab_valid') " / `reps'" _n
file write txtfile "" _n
file write txtfile "System GMM Estimator:" _n
file write txtfile "  Bias (mean - true): " %9.6f (sys_bias) _n
file write txtfile "  Standard Error: " %9.6f (sys_sd) _n
file write txtfile "  Valid Replications: " %9.0f (`sys_valid') " / `reps'" _n
file write txtfile "" _n
file write txtfile "=========================================================================" _n
file close txtfile

// Save individual results
local filename "exp3_model`model'_N`N'_rho`rho'.dta"
preserve
clear
svmat Results, names(col)
save "`output_dir'/`filename'", replace
restore

// Create individual CSV file for this run
preserve
clear
svmat Results, names(col)
gen str model = "`model'"
gen str experiment_type = "exp3_variance_ratio"
gen str date_time = "`c(current_date)' `c(current_time)'"
local csv_filename "`output_dir'/exp3_model`model'_N`N'_rho`rho'.csv"
export delimited using "`csv_filename'", replace
restore

display as text ""
display as text "Results saved to:"
display as text "  Stata format: `output_dir'/`filename'"
display as text "  Text format:  `txt_filename'"
display as text "  CSV format:   `csv_filename'"
display as text ""
display as text "========================================================================="
display as text "EXPERIMENT 3 COMPLETE FOR N=`N', MODEL=`model', RHO=`rho'"
display as text "========================================================================="

log close _all
