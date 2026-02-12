* =============================================================================
* Test suite for ivsensemakr
* Compares Stata output against R reference values from iv.sensemakr
* =============================================================================

clear all
set more off

* Load the package from repo root
adopath ++ "./"

* Load data
use card.dta, clear

* =========================================================================
* Test 1: Basic command with benchmarks
* =========================================================================
display _n "=========================================="
display "TEST 1: Basic command with benchmarks"
display "=========================================="

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa)

* =========================================================================
* Test 2: Verify ereturn scalars against R reference values
* =========================================================================
display _n "=========================================="
display "TEST 2: Verify ereturn scalars"
display "=========================================="

* R reference values
local r_iv_estimate = 0.131503836244589
local r_fs_coef     = 0.319898940091703
local r_fs_se       = 0.087863817795222
local r_rf_coef     = 0.042067937832637
local r_rf_se       = 0.018077600953643
local r_rho         = 0.361814301664605
local r_rv_iv       = 0.006666407438599
local r_xrv_iv      = 0.000523244341685
local r_rv_fs       = 0.030231294087421
local r_xrv_fs      = 0.003129076419604
local r_rv_rf       = 0.006666407438599
local r_xrv_rf      = 0.000523244341685
local r_ar_ci_lwr   = 0.024804835965035
local r_ar_ci_upr   = 0.284823593338720

* Check each value
local tol = 1e-6
local pass = 0
local fail = 0

foreach var in iv_estimate fs_coef fs_se rf_coef rf_se rho rv_iv xrv_iv rv_fs xrv_fs rv_rf xrv_rf {
    local stata_val = e(`var')
    local r_val = `r_`var''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol') {
        display as text "  PASS: `var' = " %12.10f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local pass = `pass' + 1
    }
    else {
        display as error "  FAIL: `var' = " %12.10f `stata_val' " vs R = " %12.10f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local fail = `fail' + 1
    }
}

* Check CI bounds (slightly looser tolerance)
local tol_ci = 1e-4
foreach var in ci_lwr ci_upr {
    if ("`var'" == "ci_lwr") local evar "ci_lwr"
    if ("`var'" == "ci_upr") local evar "ci_upr"
    if ("`var'" == "ci_lwr") local r_val = `r_ar_ci_lwr'
    if ("`var'" == "ci_upr") local r_val = `r_ar_ci_upr'
    local stata_val = e(`evar')
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_ci') {
        display as text "  PASS: `var' = " %12.10f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local pass = `pass' + 1
    }
    else {
        display as error "  FAIL: `var' = " %12.10f `stata_val' " vs R = " %12.10f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local fail = `fail' + 1
    }
}

display _n as text "Scalar tests: `pass' passed, `fail' failed"

* =========================================================================
* Test 3: Verify IV bounds matrix
* =========================================================================
display _n "=========================================="
display "TEST 3: Verify IV bounds"
display "=========================================="

matrix list e(bounds_iv)

* R reference: IV bounds for black (row 1) and smsa (row 2)
* Columns: kz, ky, r2zw_x, r2y0w_zx, lower_CI, upper_CI
local r_iv_black_r2zw   = 0.002214714828641
local r_iv_black_r2y0w  = 0.074999286087440
local r_iv_black_lwr    = -0.021215585945664
local r_iv_black_upr    = 0.401911999411675

local r_iv_smsa_r2zw    = 0.006394072343943
local r_iv_smsa_r2y0w   = 0.020182013342911
local r_iv_smsa_lwr     = -0.019230595280639
local r_iv_smsa_upr     = 0.395750574468739

local tol_bounds = 1e-4
local bpass = 0
local bfail = 0

matrix B = e(bounds_iv)

* Black row
foreach col_pair in "3 r2zw" "4 r2y0w" "5 lwr" "6 upr" {
    local col : word 1 of `col_pair'
    local name : word 2 of `col_pair'
    local stata_val = B[1, `col']
    local r_val = `r_iv_black_`name''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_bounds') {
        display as text "  PASS: IV black `name' = " %12.8f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local bpass = `bpass' + 1
    }
    else {
        display as error "  FAIL: IV black `name' = " %12.8f `stata_val' " vs R = " %12.8f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local bfail = `bfail' + 1
    }
}

* SMSA row
foreach col_pair in "3 r2zw" "4 r2y0w" "5 lwr" "6 upr" {
    local col : word 1 of `col_pair'
    local name : word 2 of `col_pair'
    local stata_val = B[2, `col']
    local r_val = `r_iv_smsa_`name''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_bounds') {
        display as text "  PASS: IV smsa `name' = " %12.8f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local bpass = `bpass' + 1
    }
    else {
        display as error "  FAIL: IV smsa `name' = " %12.8f `stata_val' " vs R = " %12.8f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local bfail = `bfail' + 1
    }
}

display _n as text "IV bounds tests: `bpass' passed, `bfail' failed"

* =========================================================================
* Test 4: Verify FS bounds matrix
* =========================================================================
display _n "=========================================="
display "TEST 4: Verify FS bounds"
display "=========================================="

matrix list e(bounds_fs)

local r_fs_black_r2zw  = 0.002214714828641
local r_fs_black_r2dw  = 0.033418338559507
local r_fs_black_lwr   = 0.108899742662464
local r_fs_black_upr   = 0.530898137520942

local r_fs_smsa_r2zw   = 0.006394072343943
local r_fs_smsa_r2dw   = 0.004981186551447
local r_fs_smsa_lwr    = 0.120248363377218
local r_fs_smsa_upr    = 0.519549516806188

local fspass = 0
local fsfail = 0

matrix BFS = e(bounds_fs)

* Black row
foreach col_pair in "3 r2zw" "4 r2dw" "5 lwr" "6 upr" {
    local col : word 1 of `col_pair'
    local name : word 2 of `col_pair'
    local stata_val = BFS[1, `col']
    local r_val = `r_fs_black_`name''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_bounds') {
        display as text "  PASS: FS black `name' = " %12.8f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local fspass = `fspass' + 1
    }
    else {
        display as error "  FAIL: FS black `name' = " %12.8f `stata_val' " vs R = " %12.8f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local fsfail = `fsfail' + 1
    }
}

* SMSA row
foreach col_pair in "3 r2zw" "4 r2dw" "5 lwr" "6 upr" {
    local col : word 1 of `col_pair'
    local name : word 2 of `col_pair'
    local stata_val = BFS[2, `col']
    local r_val = `r_fs_smsa_`name''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_bounds') {
        display as text "  PASS: FS smsa `name' = " %12.8f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local fspass = `fspass' + 1
    }
    else {
        display as error "  FAIL: FS smsa `name' = " %12.8f `stata_val' " vs R = " %12.8f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local fsfail = `fsfail' + 1
    }
}

display _n as text "FS bounds tests: `fspass' passed, `fsfail' failed"

* =========================================================================
* Test 5: Verify RF bounds matrix
* =========================================================================
display _n "=========================================="
display "TEST 5: Verify RF bounds"
display "=========================================="

matrix list e(bounds_rf)

local r_rf_black_r2zw  = 0.002214714828641
local r_rf_black_r2yw  = 0.065659478818227
local r_rf_black_lwr   = -0.004179559451282
local r_rf_black_upr   = 0.088315435116555

local r_rf_smsa_r2zw   = 0.006394072343943
local r_rf_smsa_r2yw   = 0.019733114638746
local r_rf_smsa_lwr    = -0.004291689581061
local r_rf_smsa_upr    = 0.088427565246335

local rfpass = 0
local rffail = 0

matrix BRF = e(bounds_rf)

* Black row
foreach col_pair in "3 r2zw" "4 r2yw" "5 lwr" "6 upr" {
    local col : word 1 of `col_pair'
    local name : word 2 of `col_pair'
    local stata_val = BRF[1, `col']
    local r_val = `r_rf_black_`name''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_bounds') {
        display as text "  PASS: RF black `name' = " %12.8f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local rfpass = `rfpass' + 1
    }
    else {
        display as error "  FAIL: RF black `name' = " %12.8f `stata_val' " vs R = " %12.8f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local rffail = `rffail' + 1
    }
}

* SMSA row
foreach col_pair in "3 r2zw" "4 r2yw" "5 lwr" "6 upr" {
    local col : word 1 of `col_pair'
    local name : word 2 of `col_pair'
    local stata_val = BRF[2, `col']
    local r_val = `r_rf_smsa_`name''
    local rdiff = reldif(`stata_val', `r_val')
    if (`rdiff' < `tol_bounds') {
        display as text "  PASS: RF smsa `name' = " %12.8f `stata_val' " (reldif = " %8.2e `rdiff' ")"
        local rfpass = `rfpass' + 1
    }
    else {
        display as error "  FAIL: RF smsa `name' = " %12.8f `stata_val' " vs R = " %12.8f `r_val' " (reldif = " %8.2e `rdiff' ")"
        local rffail = `rffail' + 1
    }
}

display _n as text "RF bounds tests: `rfpass' passed, `rffail' failed"

* =========================================================================
* Summary
* =========================================================================
local total_pass = `pass' + `bpass' + `fspass' + `rfpass'
local total_fail = `fail' + `bfail' + `fsfail' + `rffail'
local total = `total_pass' + `total_fail'

display _n "=========================================="
display "TOTAL: `total_pass'/`total' tests passed, `total_fail' failed"
display "=========================================="

if (`total_fail' > 0) {
    error 9
}
