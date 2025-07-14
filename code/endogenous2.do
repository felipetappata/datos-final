/* ========================================================================= 
   REPLICATION OF AL SADOON ET AL. (2019) - ENDOGENOUS + FULL SAMPLE
   =========================================================================
   
   This script replicates both:
   1. "Endogenous selection" columns of Table 1 (selected sample estimates)
   2. "AB all" and "system all" from Figure 1 (full sample estimates)
   
   For each Monte Carlo replication, we estimate on both:
   - Selected sample: y = y_star if d == 1 (endogenous selection)
   - Full sample: y_full = y_star (no selection, complete NxT sample)
   
   display as text ""
display as text "Results saved to:"
display as text "SELECTED SAMPLE (endogenous selection):"
display as text "  Stata format: output/partial/`filename'"
display as text "  Text format:  `txt_filename'"
display as text "  CSV format:   `csv_filename'"
display as text "FULL SAMPLE (no selection):"
display as text "  Stata format: output/partial/`filename_full'"
display as text "  Text format:  `txt_filename_full'"
display as text "  CSV format:   `csv_filename_full'"
display as text ""
display as text "========================================================================="
display as text "SIMULATION COMPLETE FOR N=`N', MODEL=`model', RHO=`rho'"
display as text "Generated both endogenous selection and full sample estimates"
display as text "========================================================================"e usage:
   stata -b do endogenous2.do N model rho
   
   Where:
   - N: Sample size (>= 200)
   - model: Selection model (A=static, B=dynamic)  
   - rho: Autoregressive parameter (>= 0)
   
   Example: stata -b do endogenous2.do 500 A 0.25
   
   Output files:
   - endo_modelA_N500_rho0.25.csv (selected sample estimates)
   - fullsample_modelA_N500_rho0.25.csv (full sample estimates)
   
   Author: Generated for replication study
   Date: January 2025
   ========================================================================= */

clear all
set more off
set seed 08869  // As specified in instructions

// Optimize Mata for speed (helps with xtabond2 performance)
mata: mata set matafavor speed, perm

// Parse command line arguments
if "`1'" == "" | "`2'" == "" | "`3'" == "" {
    display as error "Usage: stata -b do endogenous.do N model rho"
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
local a_param 1.794            // Set so P(d_it* > 0) = 0.85 (15% selection)
local sigma_z 1              // Standard deviation of z_it ~ N(0,1)
local sigma_eta 1            // Standard deviation of eta_i ~ N(0,1)  
local sigma_u 1              // Standard deviation of u_it ~ N(0,1)

// Outcome equation parameters (ENDOGENOUS SELECTION)
local sigma_alpha0 1         // Standard deviation of alpha_i^0 ~ N(0,1)
local sigma_eps0 1           // Standard deviation of epsilon_it^0 ~ N(0,1)
local theta_param 0.5        // Correlation parameter alpha_i = alpha_i^0 + theta*eta_i
local vartheta_param 0.5     // Correlation parameter eps_it = eps_it^0 + vartheta*u_it

// Implied correlations: corr(eps_it, u_it) = corr(alpha_i, eta_i) = 0.5/sqrt(1+0.5^2) = 0.447

// ========================================================================
// MONTE CARLO SIMULATION
// ========================================================================

display as text ""
display as text "========================================================================="
display as text "MONTE CARLO SIMULATION - ENDOGENOUS SELECTION + FULL SAMPLE"
display as text "========================================================================="
display as text "Parameters: N=`N', Model=`model', rho=`rho', Replications=`reps'"
display as text "Selection: 15% attrition, ENDOGENOUS (theta=`theta_param', vartheta=`vartheta_param')"
local implied_corr = `theta_param' / sqrt(1 + `theta_param'^2)
display as text "Implied correlation: " %5.3f `implied_corr'
display as text "Estimating on: (1) Selected sample, (2) Full sample (no selection)"
display as text "========================================================================="
display as text ""

// Initialize results storage for SELECTED SAMPLE (endogenous selection)
matrix AB_coefs = J(`reps', 1, .)
matrix AB_ses = J(`reps', 1, .)
matrix SYS_coefs = J(`reps', 1, .)
matrix SYS_ses = J(`reps', 1, .)

// Initialize results storage for FULL SAMPLE (no selection)
matrix AB_full_coefs = J(`reps', 1, .)
matrix AB_full_ses = J(`reps', 1, .)
matrix SYS_full_coefs = J(`reps', 1, .)
matrix SYS_full_ses = J(`reps', 1, .)

local ab_valid = 0
local sys_valid = 0
local ab_converged = 0
local sys_converged = 0

local ab_full_valid = 0
local sys_full_valid = 0
local ab_full_converged = 0
local sys_full_converged = 0

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
        
        // Endogenous selection: correlation between equations via heterogeneity
        gen double alpha_i = alpha_i0 + `theta_param' * eta_i
        
        // Generate time-varying components
        gen double eps_i0 = rnormal(0, `sigma_eps0')
        gen double u_it = rnormal(0, `sigma_u')
        gen double z_it = rnormal(0, `sigma_z')
        
        // Endogenous selection: correlation between equations via shocks
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
        
        // Observe y only when selected (for endogenous selection estimates)
        gen double y = y_star if d == 1
        
        // Create full sample outcome (for "AB all" and "system all" estimates)
        gen double y_full = y_star
        
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
    // ESTIMATION (UNCORRECTED FOR SELECTION) - SELECTED SAMPLE
    // ====================================================================
    
    // ARELLANO-BOND (AB) FIRST-DIFFERENCE GMM ESTIMATOR - SELECTED SAMPLE
    // According to Al Sadoon et al., this should be consistent even under endogenous selection
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
    
    // SYSTEM GMM ESTIMATOR (AB + Levels) - SELECTED SAMPLE
    // According to Al Sadoon et al., this should have small bias under endogenous selection
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
    
    // ====================================================================
    // ESTIMATION - FULL SAMPLE (NO SELECTION) - "AB all" and "system all"
    // ====================================================================
    
    // ARELLANO-BOND (AB) FIRST-DIFFERENCE GMM ESTIMATOR - FULL SAMPLE
    capture noisily xtabond2 y_full L.y_full, gmm(L.y_full, collapse) nolevel
    
    if _rc == 0 {
        local ab_full_converged = `ab_full_converged' + 1
        matrix b_ab_full = e(b)
        matrix V_ab_full = e(V)
        
        // Extract coefficient on L.y_full
        local pos_ab_full = colnumb(b_ab_full, "L.y_full")
        if !missing(`pos_ab_full') & `pos_ab_full' > 0 {
            scalar ab_full_coef = b_ab_full[1, `pos_ab_full']
            scalar ab_full_se = sqrt(V_ab_full[`pos_ab_full', `pos_ab_full'])
            
            // Store results if valid
            if !missing(ab_full_coef) & !missing(ab_full_se) & ab_full_se > 0 {
                matrix AB_full_coefs[`rep', 1] = ab_full_coef
                matrix AB_full_ses[`rep', 1] = ab_full_se
                local ab_full_valid = `ab_full_valid' + 1
            }
        }
    }
    
    // SYSTEM GMM ESTIMATOR (AB + Levels) - FULL SAMPLE
    capture noisily xtabond2 y_full L.y_full, gmm(L.y_full, lag(2 .)) iv(L.D.y_full, equation(level))
    
    if _rc == 0 {
        local sys_full_converged = `sys_full_converged' + 1
        matrix b_sys_full = e(b)
        matrix V_sys_full = e(V)
        
        // Extract coefficient on L.y_full
        local pos_sys_full = colnumb(b_sys_full, "L.y_full")
        if !missing(`pos_sys_full') & `pos_sys_full' > 0 {
            scalar sys_full_coef = b_sys_full[1, `pos_sys_full']
            scalar sys_full_se = sqrt(V_sys_full[`pos_sys_full', `pos_sys_full'])
            
            // Store results if valid
            if !missing(sys_full_coef) & !missing(sys_full_se) & sys_full_se > 0 {
                matrix SYS_full_coefs[`rep', 1] = sys_full_coef
                matrix SYS_full_ses[`rep', 1] = sys_full_se
                local sys_full_valid = `sys_full_valid' + 1
            }
        }
    }
    
    // Progress indicator for large replications
    if mod(`rep', 100) == 0 & `rep' > 0 {
        display as text "    Replication `rep': Selected[AB=`ab_valid', SYS=`sys_valid'] Full[AB=`ab_full_valid', SYS=`sys_full_valid']"
    }
}

// ========================================================================
// COMPUTE AND SAVE RESULTS
// ========================================================================

display as text ""
display as text "========================================================================="
display as text "SIMULATION COMPLETED"
display as text "========================================================================="
display as text "SELECTED SAMPLE (endogenous selection):"
display as text "  AB estimations: `ab_converged' converged, `ab_valid' valid results"
display as text "  SYS estimations: `sys_converged' converged, `sys_valid' valid results"
display as text "FULL SAMPLE (no selection):"
display as text "  AB estimations: `ab_full_converged' converged, `ab_full_valid' valid results"
display as text "  SYS estimations: `sys_full_converged' converged, `sys_full_valid' valid results"
display as text ""

// Convert matrices to variables for analysis
clear
svmat AB_coefs, names(ab_coef)
svmat AB_ses, names(ab_se)
svmat SYS_coefs, names(sys_coef)
svmat SYS_ses, names(sys_se)
svmat AB_full_coefs, names(ab_full_coef)
svmat AB_full_ses, names(ab_full_se)
svmat SYS_full_coefs, names(sys_full_coef)
svmat SYS_full_ses, names(sys_full_se)

// Compute summary statistics for SELECTED SAMPLE
if `ab_valid' >= 10 {
    summarize ab_coef1, detail
    scalar ab_mean = r(mean)
    scalar ab_sd = r(sd)
    scalar ab_bias = ab_mean - `rho'
    
    display as text "SELECTED SAMPLE - AB ESTIMATOR RESULTS:"
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
    display as text "SELECTED SAMPLE - SYSTEM GMM ESTIMATOR RESULTS:"
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

// Compute summary statistics for FULL SAMPLE
if `ab_full_valid' >= 10 {
    summarize ab_full_coef1, detail
    scalar ab_full_mean = r(mean)
    scalar ab_full_sd = r(sd)
    scalar ab_full_bias = ab_full_mean - `rho'
    
    display as text ""
    display as text "FULL SAMPLE - AB ESTIMATOR RESULTS:"
    display as text "  Valid replications: `ab_full_valid'"
    display as text "  Mean coefficient: " %9.6f ab_full_mean
    display as text "  Standard deviation: " %9.6f ab_full_sd
    display as text "  Bias (mean - `rho'): " %9.6f ab_full_bias
}
else {
    scalar ab_full_mean = .
    scalar ab_full_sd = .
    scalar ab_full_bias = .
    display as error "WARNING: Insufficient valid AB full results (`ab_full_valid' < 10)"
}

if `sys_full_valid' >= 10 {
    summarize sys_full_coef1, detail
    scalar sys_full_mean = r(mean)
    scalar sys_full_sd = r(sd) 
    scalar sys_full_bias = sys_full_mean - `rho'
    
    display as text ""
    display as text "FULL SAMPLE - SYSTEM GMM ESTIMATOR RESULTS:"
    display as text "  Valid replications: `sys_full_valid'"
    display as text "  Mean coefficient: " %9.6f sys_full_mean
    display as text "  Standard deviation: " %9.6f sys_full_sd
    display as text "  Bias (mean - `rho'): " %9.6f sys_full_bias
}
else {
    scalar sys_full_mean = .
    scalar sys_full_sd = .
    scalar sys_full_bias = .
    display as error "WARNING: Insufficient valid SYS full results (`sys_full_valid' < 10)"
}

// ========================================================================
// SAVE RESULTS TO FILES - SELECTED SAMPLE (ENDOGENOUS SELECTION)
// ========================================================================

// Create results matrix for export - SELECTED SAMPLE
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

// Display final results table for SELECTED SAMPLE
display as text ""
display as text "========================================================================="
display as text "FINAL RESULTS - ENDOGENOUS SELECTION (SELECTED SAMPLE)"
display as text "Model `model', N=`N', rho=`rho'"
display as text "========================================================================="
matrix list Results, format(%9.6f)

// Export SELECTED SAMPLE results to text file
local txt_filename "output/partial/endo_model`model'_N`N'_rho`rho'.txt"
file open txtfile using "`txt_filename'", write replace
file write txtfile "=========================================================================" _n
file write txtfile "AL SADOON ET AL. (2019) REPLICATION - ENDOGENOUS SELECTION RESULTS" _n
file write txtfile "=========================================================================" _n
file write txtfile "Model: `model' (Static=A, Dynamic=B)" _n
file write txtfile "Sample Size (N): `N'" _n  
file write txtfile "Autoregressive Parameter (rho): `rho'" _n
file write txtfile "Number of Replications: `reps'" _n
file write txtfile "Selection Type: Endogenous (theta=0.5, vartheta=0.5)" _n
local implied_corr = 0.5 / sqrt(1 + 0.5^2)
file write txtfile "Implied Correlation: " %5.3f (`implied_corr') _n
file write txtfile "Date: `c(current_date)' `c(current_time)'" _n
file write txtfile "=========================================================================" _n
file write txtfile "" _n
file write txtfile "ESTIMATION RESULTS (SELECTED SAMPLE):" _n
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
file write txtfile "=========================================================================" _n
file close txtfile

// Save SELECTED SAMPLE to file with descriptive name
local filename "endo_model`model'_N`N'_rho`rho'.dta"
preserve
clear
svmat Results, names(col)
save "output/partial/`filename'", replace
restore

// Create individual CSV file for SELECTED SAMPLE
preserve
clear
svmat Results, names(col)
gen str model = "`model'"
gen str selection_type = "endogenous"
gen str date_time = "`c(current_date)' `c(current_time)'"
local csv_filename "output/partial/endo_model`model'_N`N'_rho`rho'.csv"
export delimited using "`csv_filename'", replace
restore

// ========================================================================
// SAVE RESULTS TO FILES - FULL SAMPLE (NO SELECTION)
// ========================================================================

// Create results matrix for export - FULL SAMPLE
matrix Results_full = J(1, 8, .)
matrix colnames Results_full = N rho AB_bias AB_se SYS_bias SYS_se AB_valid SYS_valid

matrix Results_full[1,1] = `N'
matrix Results_full[1,2] = `rho'
matrix Results_full[1,3] = ab_full_bias
matrix Results_full[1,4] = ab_full_sd
matrix Results_full[1,5] = sys_full_bias
matrix Results_full[1,6] = sys_full_sd
matrix Results_full[1,7] = `ab_full_valid'
matrix Results_full[1,8] = `sys_full_valid'

// Display final results table for FULL SAMPLE
display as text ""
display as text "========================================================================="
display as text "FINAL RESULTS - FULL SAMPLE (NO SELECTION)"
display as text "Model `model', N=`N', rho=`rho'"
display as text "========================================================================="
matrix list Results_full, format(%9.6f)

// Export FULL SAMPLE results to text file
local txt_filename_full "output/partial/fullsample_model`model'_N`N'_rho`rho'.txt"
file open txtfile_full using "`txt_filename_full'", write replace
file write txtfile_full "=========================================================================" _n
file write txtfile_full "AL SADOON ET AL. (2019) REPLICATION - FULL SAMPLE RESULTS" _n
file write txtfile_full "=========================================================================" _n
file write txtfile_full "Model: `model' (Static=A, Dynamic=B)" _n
file write txtfile_full "Sample Size (N): `N'" _n  
file write txtfile_full "Autoregressive Parameter (rho): `rho'" _n
file write txtfile_full "Number of Replications: `reps'" _n
file write txtfile_full "Selection Type: None (Full N x T sample)" _n
file write txtfile_full "Date: `c(current_date)' `c(current_time)'" _n
file write txtfile_full "=========================================================================" _n
file write txtfile_full "" _n
file write txtfile_full "ESTIMATION RESULTS (FULL SAMPLE):" _n
file write txtfile_full "---------—----------" _n
file write txtfile_full "AB Estimator:" _n
file write txtfile_full "  Bias (mean - true): " %9.6f (ab_full_bias) _n
file write txtfile_full "  Standard Error: " %9.6f (ab_full_sd) _n
file write txtfile_full "  Valid Replications: " %9.0f (`ab_full_valid') " / `reps'" _n
file write txtfile_full "" _n
file write txtfile_full "System GMM Estimator:" _n
file write txtfile_full "  Bias (mean - true): " %9.6f (sys_full_bias) _n
file write txtfile_full "  Standard Error: " %9.6f (sys_full_sd) _n
file write txtfile_full "  Valid Replications: " %9.0f (`sys_full_valid') " / `reps'" _n
file write txtfile_full "" _n
file write txtfile_full "=========================================================================" _n
file close txtfile_full

// Save FULL SAMPLE to file with descriptive name
local filename_full "fullsample_model`model'_N`N'_rho`rho'.dta"
preserve
clear
svmat Results_full, names(col)
save "output/partial/`filename_full'", replace
restore

// Create individual CSV file for FULL SAMPLE
preserve
clear
svmat Results_full, names(col)
gen str model = "`model'"
gen str selection_type = "full_sample"
gen str date_time = "`c(current_date)' `c(current_time)'"
local csv_filename_full "output/partial/fullsample_model`model'_N`N'_rho`rho'.csv"
export delimited using "`csv_filename_full'", replace
restore

display as text ""
display as text "Results saved to:"
display as text "  Stata format: output/partial/`filename'"
display as text "  Text format:  `txt_filename'"
display as text "  CSV format:   `csv_filename'"
display as text "Results appended to: output/partial/master_results_endogenous.dta"
display as text "Master CSV updated: output/partial/master_results_endogenous.csv"
display as text ""
display as text "========================================================================="
display as text "SIMULATION COMPLETE FOR N=`N', MODEL=`model', RHO=`rho'"
display as text "========================================================================="

log close _all
