clear all
set more off
set seed 46863

* ==== PARÁMETROS ====
local T 7
local reps 500
local a 1.8
local sigma_z 1
local sigma_eta 1
local sigma_u 1
local sigma_alpha0 1
local sigma_eps0 1
local N 500
local Tplus13 = `T' + 13
local obs = `N' * `Tplus13'

* ==== MATRIZ DE RESULTADOS FINAL ====
matrix Results = J(3, 4, .)
matrix rownames Results = rho025 rho050 rho075
matrix colnames Results = AB_bias AB_se SYS_bias SYS_se

local idx = 1

foreach rho in 0.25 0.5 0.75 {

    * Estims: columnas = [coef_AB, se_AB, coef_SYS, se_SYS]
    matrix Estims = J(`reps', 4, .)
    local ab_valid = 0
    local sys_valid = 0

    forvalues i = 1/`reps' {
        clear
        set obs `obs'

        gen id = floor((_n-1)/`Tplus13') + 1
        gen t = mod(_n-1, `Tplus13') + 1
        xtset id t

        gen alpha0 = rnormal(0,`sigma_alpha0')
        gen etai   = rnormal(0,`sigma_eta')
        sort id
        by id: replace alpha0 = alpha0[1]
        by id: replace etai = etai[1]
        gen alphai = alpha0
        gen eps0   = rnormal(0,`sigma_eps0')
        gen epsi   = eps0
        gen uit    = rnormal(0,`sigma_u')

        gen zit = rnormal(0,`sigma_z')

        * Selección no endógena (modelo A)
        gen d_star = `a' - zit - etai - uit
        gen d = (d_star > 0)

        * Variable latente
        gen ystar = .
        gen y0 = (2 + alphai + epsi) / (1 - `rho') 
        replace ystar = y0 if t == 1
        sort id t
        by id (t): replace ystar = 2 + `rho'*L.ystar + alphai + epsi if t > 1

        gen y = ystar if d == 1
        drop if t <= 13

        xtset id t

        * Estimación AB
        capture noisily xtabond2 y L.y, gmm(L.y, collapse) nolevel
        if _rc == 0 {
            matrix b = e(b)
            matrix V = e(V)
            local pos = colnumb(b, "L.y")
            if !missing(`pos') & `pos' > 0 {
                if b[1,`pos'] < . & V[`pos',`pos'] < . {
                    matrix Estims[`i',1] = b[1,`pos']
                    matrix Estims[`i',2] = sqrt(V[`pos',`pos'])
                    local ++ab_valid
                }
            }
        }

        * Estimación SYS
        capture noisily xtabond2 y L.y, gmm(L.y, lag(2 .)) iv(L.D.y, equation(level))
        if _rc == 0 {
            matrix b = e(b)
            matrix V = e(V)
            local pos = colnumb(b, "L.y")
            if !missing(`pos') & `pos' > 0 {
                if b[1,`pos'] < . & V[`pos',`pos'] < . {
                    matrix Estims[`i',3] = b[1,`pos']
                    matrix Estims[`i',4] = sqrt(V[`pos',`pos'])
                    local ++sys_valid
                }
            }
        }
    }

    * Procesamiento de resultados
    svmat double Estims, names(est)

    count if !missing(est1)
    local ab_count = r(N)
    count if !missing(est3)
    local sys_count = r(N)

    if `ab_count' >= 2 {
        summarize est1, meanonly
        local ab_bias = r(mean) - `rho'

        quietly summarize est1, detail
        local ab_sd = r(sd)
    }
    else {
        local ab_bias = .
        local ab_sd = .
    }

    if `sys_count' >= 2 {
        summarize est3, meanonly
        local sys_bias = r(mean) - `rho'

        quietly summarize est3, detail
        local sys_sd = r(sd)
    }
    else {
        local sys_bias = .
        local sys_sd = .
    }

    matrix Results[`idx',1] = `ab_bias'
    matrix Results[`idx',2] = `ab_sd'
    matrix Results[`idx',3] = `sys_bias'
    matrix Results[`idx',4] = `sys_sd'

    di ">>> rho = `rho': AB válidas = `ab_valid' (`ab_count' usadas), SYS válidas = `sys_valid' (`sys_count' usadas)"
    local ++idx
}

* ==== IMPRIMIR RESULTADOS ====
di as text "===> Resultados Monte Carlo - Modelo A (no endógena), N=500, T=7"
set cformat %9.6f
matrix list Results

