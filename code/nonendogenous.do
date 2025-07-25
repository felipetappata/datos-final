/* ========================================================================= 
   REPLICATION OF AL SADOON ET AL. (2019) TABLE 1 - NON-ENDOGENOUS SELECTION
   =========================================================================
   Command line usage:
   stata -b do nonendogenous.do N model rho
   
   Where:
   - N: Sample size (>= 200)
   - model: Selection model (A=static, B=dynamic)  
   - rho: Autoregressive parameter (>= 0)
   
   Example: stata -b do nonendogenous.do 500 A 0.25
   
   Author: Felipe I. Tappata
   Date: July 2025
   ========================================================================= */

clear all
set more off
set seed 08869  // As specified in instructions

// Optimize Mata for speed (helps with xtabond2 performance)
mata: mata set matafavor speed, perm

// Parse command line arguments
if "`1'" == "" | "`2'" == "" | "`3'" == "" {
    display as error "Usage: stata -b do nonendogenous.do N model rho"
    display as error "Where N>=200, model={A,B}, rho>=0"
    exit 198
}

local N `1'
local model `2'  
local rho `3'

// Validate inputs
if `N' < 200 {
    display as error "N must be >= 200"
    exit 198
}
if !inlist("`model'", "A", "B") {
    display as error "Model must be A (static) or B (dynamic)"
    exit 198
}
if `rho' < 0 {
    display as error "rho must be >= 0"
    exit 198
}

// ========================================================================
// SIMULATION PARAMETERS (Following Section 3 of Al Sadoon et al. 2019)
// ========================================================================

// Define sample sizes for experimentation (can be reduced for testing)
local N_values "200 500 1000 5000"  // Full range for comprehensive study
local N_experiment "50 100 200"     // Reduced values for quick experimentation

// Core simulation parameters exactly as described in the paper
local T 7                    // Time periods for estimation (after discarding 13)
local T_total 20             // Total periods generated (discarding first 13)
local T_discard 13           // Periods to discard for initial conditions
local reps 500               // Number of Monte Carlo replications

// Selection equation parameters  
local a_param 1.794            // Set so P(d_it* > 0) = 0.85 (15% selection), using using a=√3*Φ^−1(0.85)=1.794
local sigma_z 1              // Standard deviation of z_it ~ N(0,1)
local sigma_eta 1            // Standard deviation of eta_i ~ N(0,1)  
local sigma_u 1              // Standard deviation of u_it ~ N(0,1)

// Outcome equation parameters (NON-ENDOGENOUS SELECTION)
local sigma_alpha0 1         // Standard deviation of alpha_i^0 ~ N(0,1)
local sigma_eps0 1           // Standard deviation of epsilon_it^0 ~ N(0,1)
local theta_param 0          // Correlation parameter alpha_i = alpha_i^0 + theta*eta_i (=0 for non-endogenous)
local vartheta_param 0       // Correlation parameter eps_it = eps_it^0 + vartheta*u_it (=0 for non-endogenous)

// ========================================================================
// MONTE CARLO SIMULATION
// ========================================================================

display as text ""
display as text "========================================================================="
display as text "MONTE CARLO SIMULATION - NON-ENDOGENOUS SELECTION"
display as text "========================================================================="
display as text "Parameters: N=`N', Model=`model', rho=`rho', Replications=`reps'"
display as text "Selection: 15% attrition, NON-ENDOGENOUS (theta=`theta_param', vartheta=`vartheta_param')"
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
    // DATA GENERATION PROCESS (Exactly as in Section 3)
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
        gen double eta_i = rnormal(0, `sigma_eta')
        sort id
        by id: replace alpha_i0 = alpha_i0[1]
        by id: replace eta_i = eta_i[1]
        
        // Non-endogenous selection: no correlation between equations
        gen double alpha_i = alpha_i0 + `theta_param' * eta_i
        
        // Generate time-varying components
        gen double eps_i0 = rnormal(0, `sigma_eps0')
        gen double u_it = rnormal(0, `sigma_u')
        gen double z_it = rnormal(0, `sigma_z')
        
        // Non-endogenous selection: no correlation between equations
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
display as text "SIMULATION COMPLETED"
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
matrix Results = J(1, 8, .)
matrix colnames Results = N rho AB_bias AB_se SYS_bias SYS_se AB_valid SYS_valid

matrix Results[1,1] = `N'
matrix Results[1,2] = `rho'
matrix Results[1,3] = ab_bias
matrix Results[1,4] = ab_sd
matrix Results[1,5] = sys_bias
matrix Results[1,6] = sys_sd
matrix Results[1,7] = `ab_valid'
matrix Results[1,8] = `sys_valid'

// Display final results table
display as text ""
display as text "========================================================================="
display as text "FINAL RESULTS - NON-ENDOGENOUS SELECTION"
display as text "Model `model', N=`N', rho=`rho'"
display as text "========================================================================="
matrix list Results, format(%9.6f)

// Export results to text file for easier reading
local txt_filename "output/partial/nonendo_model`model'_N`N'_rho`rho'.txt"
file open txtfile using "`txt_filename'", write replace
file write txtfile "=========================================================================" _n
file write txtfile "AL SADOON ET AL. (2019) REPLICATION - NON-ENDOGENOUS SELECTION RESULTS" _n
file write txtfile "=========================================================================" _n
file write txtfile "Model: `model' (Static=A, Dynamic=B)" _n
file write txtfile "Sample Size (N): `N'" _n  
file write txtfile "Autoregressive Parameter (rho): `rho'" _n
file write txtfile "Number of Replications: `reps'" _n
file write txtfile "Selection Type: Non-endogenous (theta=0, vartheta=0)" _n
file write txtfile "Date: `c(current_date)' `c(current_time)'" _n
file write txtfile "=========================================================================" _n
file write txtfile "" _n
file write txtfile "ESTIMATION RESULTS:" _n
file write txtfile "---------—----------" _n
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
file write txtfile "Paper Predictions (Non-endogenous selection):" _n
file write txtfile "  AB: Should be approximately unbiased" _n
file write txtfile "  SYS: Should be approximately unbiased" _n
file write txtfile "" _n
file write txtfile "=========================================================================" _n
file close txtfile

// Save to file with descriptive name
local filename "nonendo_model`model'_N`N'_rho`rho'.dta"
preserve
clear
svmat Results, names(col)
save "output/partial/`filename'", replace
restore

// Also append to master results file
preserve
clear
svmat Results, names(col)
gen str model = "`model'"
gen str selection_type = "non_endogenous"
capture append using "output/partial/master_results_nonendogenous.dta"
save "output/partial/master_results_nonendogenous.dta", replace

// Export to CSV for easy reading
export delimited using "output/partial/master_results_nonendogenous.csv", replace
restore

// Create individual CSV file for this run
preserve
clear
svmat Results, names(col)
gen str model = "`model'"
gen str selection_type = "non_endogenous"
gen str date_time = "`c(current_date)' `c(current_time)'"
local csv_filename "output/partial/nonendo_model`model'_N`N'_rho`rho'.csv"
export delimited using "`csv_filename'", replace
restore

display as text ""
display as text "Results saved to:"
display as text "  Stata format: output/partial/`filename'"
display as text "  Text format:  `txt_filename'"
display as text "  CSV format:   `csv_filename'"
display as text "Results appended to: output/partial/master_results_nonendogenous.dta"
display as text "Master CSV updated: output/partial/master_results_nonendogenous.csv"
display as text ""
display as text "========================================================================="
display as text "SIMULATION COMPLETE FOR N=`N', MODEL=`model', RHO=`rho'"
display as text "========================================================================="

log close _all
