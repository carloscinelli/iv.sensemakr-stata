*! version 1.0.0  07feb2026
*! Sensitivity Analysis for Instrumental Variable Estimates
*! Authors: Carlos Cinelli and Chad Hazlett
*! Reference: Cinelli and Hazlett (2025), Biometrika

program define ivsensemakr, eclass
version 13

syntax varlist(min=3 ts fv) [if] [in] [, ///
	Benchmark(varlist ts fv) ///
	kz(numlist min=1 > 0) ///
	ky(numlist min=1 > 0) ///
	kd(numlist min=1 > 0) ///
	r2zw(numlist min=1 >=0 <=1) ///
	r2y0w(numlist min=1 >=0 <=1) ///
	bound_label(string) ///
	q(real 1) ///
	alpha(real .05) ///
	h0(real 0) ///
	CONTOURplot ///
	TCONTOURplot ///
	Clim(numlist min=2 max=2) ///
	Clines(real 7) ///
	noMIN ///
	Suppress ///
]

marksample touse

// =========================================================================
// Parse variables: first 3 are depvar, treatvar, instrvar; rest are covariates
// =========================================================================
local depvar   : word 1 of `varlist'
local treatvar : word 2 of `varlist'
local instrvar : word 3 of `varlist'
local covars   : list varlist - depvar
local covars   : list covars - treatvar
local covars   : list covars - instrvar

// =========================================================================
// Options and error handling
// =========================================================================
if (`q' <= 0) {
	display as error "q must be greater than 0"
	exit 198
}

if (`alpha' >= 1 | `alpha' <= 0) {
	display as error "alpha must be between 0 and 1"
	exit 198
}

// min option
if ("`min'" == "nomin") {
	local use_min = 0
}
else {
	local use_min = 1
}

// kz defaults
if ("`kz'" == "") {
	local kz 1
}

// ky defaults to kz
if ("`ky'" == "") {
	local ky `kz'
}

// kd defaults to kz
if ("`kd'" == "") {
	local kd `kz'
}

// Validate kz/ky/kd lengths match
local count_kz : word count `kz'
local count_ky : word count `ky'
local count_kd : word count `kd'

if ("`count_kz'" != "`count_ky'") {
	display as error "ky must have the same number of elements as kz"
	exit 198
}
if ("`count_kz'" != "`count_kd'") {
	display as error "kd must have the same number of elements as kz"
	exit 198
}

// Manual bounds defaults
if ("`r2y0w'" == "" & "`r2zw'" != "") {
	local r2y0w `r2zw'
}

if ("`bound_label'" == "") {
	local bound_label "Manual Bound"
}

// Contour plot limits
if ("`clim'" == "") {
	local lim_lb = 0
	local lim_ub = .4
	local user_clim = 0
}
else {
	local lim_lb : word 1 of `clim'
	local lim_ub : word 2 of `clim'
	local user_clim = 1
}

if (`lim_ub' > 1) {
	display as error "Max upper limit for contour plot is 1"
	exit 198
}
if (`lim_lb' < 0) {
	display as error "Min lower limit for contour plot is 0"
	exit 198
}
if (`lim_lb' >= `lim_ub') {
	display as error "Lower limit must be less than upper limit"
	exit 198
}

// =========================================================================
// Step 1: Run three regressions
// =========================================================================

// First Stage: d ~ z + x
qui: reg `treatvar' `instrvar' `covars' if `touse'
local fs_coef = _b[`instrvar']
local fs_se   = _se[`instrvar']
local fs_t    = `fs_coef' / `fs_se'
local fs_dof  = e(df_r)
local fs_p    = 2*ttail(`fs_dof', abs(`fs_t'))
local fs_lwr  = `fs_coef' - invttail(`fs_dof', `alpha'/2) * `fs_se'
local fs_upr  = `fs_coef' + invttail(`fs_dof', `alpha'/2) * `fs_se'
estimates store __ivsm_fs

// Reduced Form: y ~ z + x
qui: reg `depvar' `instrvar' `covars' if `touse'
local rf_coef = _b[`instrvar']
local rf_se   = _se[`instrvar']
local rf_t    = `rf_coef' / `rf_se'
local rf_dof  = e(df_r)
local rf_p    = 2*ttail(`rf_dof', abs(`rf_t'))
local rf_lwr  = `rf_coef' - invttail(`rf_dof', `alpha'/2) * `rf_se'
local rf_upr  = `rf_coef' + invttail(`rf_dof', `alpha'/2) * `rf_se'
estimates store __ivsm_rf

// Anderson-Rubin: (y - h0*d) ~ z + x
tempvar ar_depvar
qui: gen double `ar_depvar' = `depvar' - `h0' * `treatvar' if `touse'
qui: reg `ar_depvar' `instrvar' `covars' if `touse'
local ar_coef = _b[`instrvar']
local ar_se   = _se[`instrvar']
local ar_t    = `ar_coef' / `ar_se'
local ar_dof  = e(df_r)
local ar_p    = 2*ttail(`ar_dof', abs(`ar_t'))

// =========================================================================
// Step 2: Compute rho (correlation of FS and RF residuals)
// =========================================================================
tempvar fs_resid rf_resid
qui: estimates restore __ivsm_fs
qui: predict double `fs_resid' if `touse', resid
qui: estimates restore __ivsm_rf
qui: predict double `rf_resid' if `touse', resid
qui: corr `fs_resid' `rf_resid' if `touse'
local rho = r(rho)

// =========================================================================
// Step 3: IV estimate and AR confidence interval
// =========================================================================
local iv_estimate = `rf_coef' / `fs_coef'

// Call Mata to compute AR confidence interval
mata: __ivsm_ar_confint(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `alpha')
// Results stored in locals: ar_ci_type, ar_ci_lwr, ar_ci_upr (and ar_ci_lwr2, ar_ci_upr2 for disjoint)

// =========================================================================
// Step 4: Compute sensitivity statistics
// =========================================================================

// Refit AR at h0 = (1-q)*iv_estimate for sensitivity
local h0_q = (1 - `q') * `iv_estimate'
tempvar ar_depvar_q
qui: gen double `ar_depvar_q' = `depvar' - `h0_q' * `treatvar' if `touse'
qui: reg `ar_depvar_q' `instrvar' `covars' if `touse'
local ar_t_q = _b[`instrvar'] / _se[`instrvar']

// RV and XRV for IV
mata: st_local("rv_iv",  strofreal(__ivsm_rv_iv(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `q', `alpha', `use_min')))
mata: st_local("xrv_iv", strofreal(__ivsm_xrv_iv(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `q', `alpha', `use_min')))

// RV and XRV for FS (standard OLS, q=1)
mata: st_local("rv_fs",  strofreal(__ivsm_rv_ols(`fs_t', `fs_dof', 1, `alpha')))
mata: st_local("xrv_fs", strofreal(__ivsm_xrv_ols(`fs_t', `fs_dof', 1, `alpha')))

// RV and XRV for RF (standard OLS, q=1)
mata: st_local("rv_rf",  strofreal(__ivsm_rv_ols(`rf_t', `rf_dof', 1, `alpha')))
mata: st_local("xrv_rf", strofreal(__ivsm_xrv_ols(`rf_t', `rf_dof', 1, `alpha')))

// =========================================================================
// Step 5: Print output tables
// =========================================================================

// Header
mata: printf("\n{txt}Sensitivity Analysis for Instrumental Variables\n")
mata: printf("{txt}(Anderson-Rubin Approach)\n")
mata: printf("{txt}{hline 65}\n")

// IV Estimates
mata: printf("\n{txt}IV Estimates:\n")
mata: printf("{txt}  Coef. Estimate: {res}%10.4f\n", `iv_estimate')
if ("`ar_ci_type'" == "bounded") {
	mata: printf("{txt}  Conf. Interval: {res}[%10.4f, %10.4f]\n", `ar_ci_lwr', `ar_ci_upr')
}
else if ("`ar_ci_type'" == "disjoint") {
	mata: printf("{txt}  Conf. Interval: {res}(-inf, %10.4f] U [%10.4f, +inf)\n", `ar_ci_lwr', `ar_ci_upr')
}
else if ("`ar_ci_type'" == "all_reals") {
	mata: printf("{txt}  Conf. Interval: {res}(-inf, +inf)\n")
}
else {
	mata: printf("{txt}  Conf. Interval: {res}Empty set\n")
}
mata: printf("{txt}  Note: H0 = %6.3f, alpha = %4.2f, df = %5.0f\n", `h0', `alpha', `ar_dof')

// Sensitivity Statistics
mata: printf("\n{txt}{hline 65}\n")
mata: printf("\n{txt}Sensitivity Statistics:\n")
mata: printf("{txt}  %16s {c |} %10s  %10s  %10s\n", "", "Estimate", "XRV_qa", "RV_qa")
mata: printf("{txt}  {hline 16}{c +}{hline 36}\n")
mata: printf("{txt}  %16s {c |} {res}%10.4f  %10.6f  %10.6f\n", "IV", `iv_estimate', `xrv_iv', `rv_iv')
mata: printf("{txt}  %16s {c |} {res}%10.4f  %10.6f  %10.6f\n", "First-Stage", `fs_coef', `xrv_fs', `rv_fs')
mata: printf("{txt}  %16s {c |} {res}%10.4f  %10.6f  %10.6f\n", "Reduced-Form", `rf_coef', `xrv_rf', `rv_rf')

if (`use_min' == 1) {
	mata: printf("\n{txt}  Note: q >= %3.2f, alpha = %4.2f\n", `q', `alpha')
}
else {
	mata: printf("\n{txt}  Note: q = %3.2f, alpha = %4.2f\n", `q', `alpha')
}

// First-Stage and Reduced-Form detail
if ("`suppress'" == "") {
	mata: printf("\n{txt}{hline 65}\n")
	mata: printf("\n{txt}First-Stage Estimates (D ~ Z | X):\n")
	mata: printf("{txt}  Coef. Estimate: {res}%10.4f\n", `fs_coef')
	mata: printf("{txt}  Standard Error: {res}%10.4f\n", `fs_se')
	mata: printf("{txt}  t-value:        {res}%10.4f\n", `fs_t')
	mata: printf("{txt}  p-value:        {res}%10.6f\n", `fs_p')
	mata: printf("{txt}  Conf. Interval: {res}[%10.4f, %10.4f]\n", `fs_lwr', `fs_upr')

	mata: printf("\n{txt}Reduced-Form Estimates (Y ~ Z | X):\n")
	mata: printf("{txt}  Coef. Estimate: {res}%10.4f\n", `rf_coef')
	mata: printf("{txt}  Standard Error: {res}%10.4f\n", `rf_se')
	mata: printf("{txt}  t-value:        {res}%10.4f\n", `rf_t')
	mata: printf("{txt}  p-value:        {res}%10.6f\n", `rf_p')
	mata: printf("{txt}  Conf. Interval: {res}[%10.4f, %10.4f]\n", `rf_lwr', `rf_upr')
}

mata: printf("{txt}{hline 65}\n")

// =========================================================================
// Step 6: Benchmark bounds
// =========================================================================

if ("`benchmark'" != "") {

	local max_strl = 15
	local bench_count = 0

	// Regressions for partial R2 of z with benchmarks
	// Regress z on covariates (excluding z itself)
	qui: reg `instrvar' `covars' if `touse'
	estimates store __ivsm_z_model

	foreach bench of local benchmark {
		local bench_count = `bench_count' + 1
	}

	// ---------------------------------------------------------------
	// IV Bounds
	// ---------------------------------------------------------------
	mata: printf("\n{txt}Bounds on Omitted Variable Bias (IV):\n")
	mata: printf("{txt}  %-15s {c |} %10s  %10s  %10s  %10s\n", "Bound Label", "R2zw.x", "R2y0w.zx", "Lower CI", "Upper CI")
	mata: printf("{txt}  {hline 15}{c +}{hline 48}\n")

	local bench_idx = 0
	local bench_row_idx = 0
	foreach bench of local benchmark {
		local bench_idx = `bench_idx' + 1

		// Get partial R2 of z with benchmark (for instrument)
		qui: estimates restore __ivsm_z_model
		local bench_t_z = _b[`bench'] / _se[`bench']
		local dof_z = e(df_r)
		local r2dxj_x = `bench_t_z' * `bench_t_z' / (`bench_t_z' * `bench_t_z' + `dof_z')

		// Residualize y, d, xj on z + remaining covariates (for potential outcome benchmark)
		local covars_minus_bench : list covars - bench

		tempvar y_zx d_zx xj_zx
		qui: reg `depvar' `instrvar' `covars_minus_bench' if `touse'
		qui: predict double `y_zx' if `touse', resid

		qui: reg `treatvar' `instrvar' `covars_minus_bench' if `touse'
		qui: predict double `d_zx' if `touse', resid

		qui: reg `bench' `instrvar' `covars_minus_bench' if `touse'
		qui: predict double `xj_zx' if `touse', resid

		// Compute maxR2 via Mata
		mata: st_local("r2yxj_x", strofreal(__ivsm_maxR2(st_data(., "`y_zx'", "`touse'"), st_data(., "`d_zx'", "`touse'"), st_data(., "`xj_zx'", "`touse'"))))

		drop `y_zx' `d_zx' `xj_zx'

		// Iterate over kz/ky values
		local ki = 0
		foreach kz_val of local kz {
			local ki = `ki' + 1
			local ky_val : word `ki' of `ky'

			// Compute bounds via Mata
			mata: __ivsm_compute_iv_bounds(`r2dxj_x', `r2yxj_x', `kz_val', `ky_val', ///
				`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `alpha')
			// Results in locals: bound_r2zw_x, bound_r2y0w_zx, bound_lwr, bound_upr

			// Label
			if (`kz_val' == `ky_val') {
				local blabel "`kz_val'x `bench'"
			}
			else {
				local blabel "`kz_val'/`ky_val'x `bench'"
			}
			local bench_row_idx = `bench_row_idx' + 1
			local __ivsm_blbl_`bench_row_idx' "`blabel'"

			// Print
			mata: printf("{txt}  %-15s {c |} {res}%10.5f  %10.5f  %10.4f  %10.4f\n", ///
				substr("`blabel'", 1, 15), `bound_r2zw_x', `bound_r2y0w_zx', `bound_lwr', `bound_upr')

			// Store in matrix
			if (`bench_idx' == 1 & `ki' == 1) {
				mat __ivsm_bounds_iv = (`kz_val', `ky_val', `bound_r2zw_x', `bound_r2y0w_zx', `bound_lwr', `bound_upr')
			}
			else {
				mat __ivsm_bounds_iv = __ivsm_bounds_iv \ (`kz_val', `ky_val', `bound_r2zw_x', `bound_r2y0w_zx', `bound_lwr', `bound_upr')
			}

			// Track limits for contour plot
			if (`user_clim' == 0) {
				if ((`bound_r2zw_x' + .1) > `lim_ub') {
					local lim_ub = `bound_r2zw_x' + .1
				}
				if ((`bound_r2y0w_zx' + .05) > `lim_ub') {
					local lim_ub = `bound_r2y0w_zx' + .05
				}
			}
		}
	}
	matrix colnames __ivsm_bounds_iv = "kz" "ky" "r2zw_x" "r2y0w_zx" "lower_CI" "upper_CI"

	// ---------------------------------------------------------------
	// FS Bounds (standard OLS sensemakr on FS model, treatment = z)
	// ---------------------------------------------------------------
	mata: printf("\n{txt}Bounds on Omitted Variable Bias (First-Stage):\n")
	mata: printf("{txt}  %-15s {c |} %10s  %10s  %10s  %10s\n", "Bound Label", "R2zw.x", "R2dw.zx", "Lower CI", "Upper CI")
	mata: printf("{txt}  {hline 15}{c +}{hline 48}\n")

	// For FS: "outcome" is d (treatment), "treatment" is z (instrument)
	// kd_fs = kz, ky_fs = kd (parameter mapping from R code)
	qui: estimates restore __ivsm_fs
	local bench_idx = 0
	foreach bench of local benchmark {
		local bench_idx = `bench_idx' + 1

		// Partial R2 of z with xj in FS model
		local bench_t_fs = _b[`bench'] / _se[`bench']
		local dof_fs_b = e(df_r)
		local r2yxj_fs = `bench_t_fs' * `bench_t_fs' / (`bench_t_fs' * `bench_t_fs' + `dof_fs_b')

		// Partial R2 of z with xj in z-model
		qui: estimates restore __ivsm_z_model
		local bench_t_z2 = _b[`bench'] / _se[`bench']
		local dof_z2 = e(df_r)
		local r2dxj_fs = `bench_t_z2' * `bench_t_z2' / (`bench_t_z2' * `bench_t_z2' + `dof_z2')

		qui: estimates restore __ivsm_fs

		local ki = 0
		foreach kz_val of local kz {
			local ki = `ki' + 1
			local kd_val : word `ki' of `kd'

			// OLS bounds: kd_fs = kz_val, ky_fs = kd_val
			mata: __ivsm_compute_ols_bounds(`r2dxj_fs', `r2yxj_fs', `kz_val', `kd_val', ///
				`fs_coef', `fs_se', `fs_dof', `alpha')

			if (`kz_val' == `kd_val') {
				local blabel = "`kz_val'x `bench'"
			}
			else {
				local blabel = "`kz_val'/`kd_val'x `bench'"
			}

			mata: printf("{txt}  %-15s {c |} {res}%10.5f  %10.5f  %10.4f  %10.4f\n", ///
				substr("`blabel'", 1, 15), `ols_r2dz_x', `ols_r2yz_dx', `ols_lwr', `ols_upr')

			if (`bench_idx' == 1 & `ki' == 1) {
				mat __ivsm_bounds_fs = (`kz_val', `kd_val', `ols_r2dz_x', `ols_r2yz_dx', `ols_lwr', `ols_upr')
			}
			else {
				mat __ivsm_bounds_fs = __ivsm_bounds_fs \ (`kz_val', `kd_val', `ols_r2dz_x', `ols_r2yz_dx', `ols_lwr', `ols_upr')
			}
		}
	}
	matrix colnames __ivsm_bounds_fs = "kz" "kd" "r2zw_x" "r2dw_zx" "lower_CI" "upper_CI"

	// ---------------------------------------------------------------
	// RF Bounds (standard OLS sensemakr on RF model, treatment = z)
	// ---------------------------------------------------------------
	mata: printf("\n{txt}Bounds on Omitted Variable Bias (Reduced-Form):\n")
	mata: printf("{txt}  %-15s {c |} %10s  %10s  %10s  %10s\n", "Bound Label", "R2zw.x", "R2yw.zx", "Lower CI", "Upper CI")
	mata: printf("{txt}  {hline 15}{c +}{hline 48}\n")

	// For RF: "outcome" is y, "treatment" is z
	// kd_rf = kz, ky_rf = ky (parameter mapping from R code)
	qui: estimates restore __ivsm_rf
	local bench_idx = 0
	foreach bench of local benchmark {
		local bench_idx = `bench_idx' + 1

		// Partial R2 of z with xj in RF model
		local bench_t_rf = _b[`bench'] / _se[`bench']
		local dof_rf_b = e(df_r)
		local r2yxj_rf = `bench_t_rf' * `bench_t_rf' / (`bench_t_rf' * `bench_t_rf' + `dof_rf_b')

		// Partial R2 of z with xj
		qui: estimates restore __ivsm_z_model
		local bench_t_z3 = _b[`bench'] / _se[`bench']
		local dof_z3 = e(df_r)
		local r2dxj_rf = `bench_t_z3' * `bench_t_z3' / (`bench_t_z3' * `bench_t_z3' + `dof_z3')

		qui: estimates restore __ivsm_rf

		local ki = 0
		foreach kz_val of local kz {
			local ki = `ki' + 1
			local ky_val : word `ki' of `ky'

			mata: __ivsm_compute_ols_bounds(`r2dxj_rf', `r2yxj_rf', `kz_val', `ky_val', ///
				`rf_coef', `rf_se', `rf_dof', `alpha')

			if (`kz_val' == `ky_val') {
				local blabel = "`kz_val'x `bench'"
			}
			else {
				local blabel = "`kz_val'/`ky_val'x `bench'"
			}

			mata: printf("{txt}  %-15s {c |} {res}%10.5f  %10.5f  %10.4f  %10.4f\n", ///
				substr("`blabel'", 1, 15), `ols_r2dz_x', `ols_r2yz_dx', `ols_lwr', `ols_upr')

			if (`bench_idx' == 1 & `ki' == 1) {
				mat __ivsm_bounds_rf = (`kz_val', `ky_val', `ols_r2dz_x', `ols_r2yz_dx', `ols_lwr', `ols_upr')
			}
			else {
				mat __ivsm_bounds_rf = __ivsm_bounds_rf \ (`kz_val', `ky_val', `ols_r2dz_x', `ols_r2yz_dx', `ols_lwr', `ols_upr')
			}
		}
	}
	matrix colnames __ivsm_bounds_rf = "kz" "ky" "r2zw_x" "r2yw_zx" "lower_CI" "upper_CI"

	mata: printf("\n{txt}{hline 65}\n")
}

// =========================================================================
// Step 7: Contour plots
// =========================================================================

if ("`contourplot'" != "" | "`tcontourplot'" != "") {

	if ("`contourplot'" != "") {
		preserve

		// Lower CI contour
		mata: __ivsm_contour_plot_iv(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `lim_lb', `lim_ub', `alpha', "lwr", `clines', `h0_q')
		matrix colnames __ivsm_contourgrid_lwr = "r2zw_x" "r2y0w_zx" "lwr_ci"

		capture: graph drop __ivsm_lwr_plot
		local round_h0 : di %6.3f `h0_q'
		local obs_val = scalar(__ivsm_obs_value)
		local round_obs : di %5.3f `obs_val'

		// Build single contourline with all levels + threshold
		local all_cuts ""
		local all_colors ""
		local all_patterns ""
		local all_widths ""
		local ncuts = rowsof(toprange_lwr)
		forvalues i = 1(1)`ncuts' {
			local cutvalue = toprange_lwr[`i', 1]
			local skip = 0
			if (abs(`cutvalue' - `h0_q') < 0.005) local skip = 1
			if (`skip' == 0) {
				local all_cuts "`all_cuts' `cutvalue'"
				local all_colors "`all_colors' black"
				local all_patterns "`all_patterns' solid"
				local all_widths "`all_widths' thin"
			}
		}
		// Add threshold as last level
		local all_cuts "`all_cuts' `h0_q'"
		local all_colors "`all_colors' red"
		local all_patterns "`all_patterns' dash"
		local all_widths "`all_widths' medthick"

		// Build scatter overlay for benchmarks
		local bench_scatter ""
		if ("`benchmark'" != "") {
			local dim = rowsof(__ivsm_bounds_iv)
			forvalues idx = 1(1)`dim' {
				local r2yz_val = __ivsm_bounds_iv[`idx', 4]
				local r2dz_val = __ivsm_bounds_iv[`idx', 3]
				local coef_val = __ivsm_bounds_iv[`idx', 5]
				local coef_val : di %5.3f `coef_val'
				local blbl "`__ivsm_blbl_`idx''"
				if (`r2yz_val' < `lim_ub' & `r2dz_val' < `lim_ub') {
					local bench_scatter "`bench_scatter' || scatteri `r2yz_val' `r2dz_val' `"`blbl' (`coef_val')"', mcolor(red) mfcolor(red) msize(small) msymbol(D) mlabcolor(black) mlabsize(vsmall) mlabposition(3)"
				}
			}
		}

		twoway contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
			name(__ivsm_lwr_plot, replace) ccuts(`all_cuts') ccolor(`all_colors') clpattern(`all_patterns') clwidths(`all_widths') nodraw ///
			xlab(, labsize(small)) ylab(, labsize(small)) ///
			ytitle("Partial R{superscript:2} of confounder(s) with pot. outcome", size(small)) ///
			xtitle("Partial R{superscript:2} of confounder(s) with instrument", size(small)) ///
			title("Lower CI Limit", size(medium)) ///
			legend(off) ///
			|| scatteri 0 0 "Observed (`round_obs')", mcolor(black) msize(small) msymbol(T) mlabcolor(black) mlabsize(vsmall) mlabposition(3) ///
			`bench_scatter'

		capture: drop sense_contour_*
		capture: mat drop toprange_lwr

		// Upper CI contour
		mata: __ivsm_contour_plot_iv(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `lim_lb', `lim_ub', `alpha', "upr", `clines', `h0_q')
		matrix colnames __ivsm_contourgrid_upr = "r2zw_x" "r2y0w_zx" "upr_ci"

		capture: graph drop __ivsm_upr_plot
		local obs_val = scalar(__ivsm_obs_value)
		local round_obs : di %5.3f `obs_val'

		local all_cuts ""
		local all_colors ""
		local all_patterns ""
		local all_widths ""
		local ncuts = rowsof(toprange_upr)
		forvalues i = 1(1)`ncuts' {
			local cutvalue = toprange_upr[`i', 1]
			local skip = 0
			if (abs(`cutvalue' - `h0_q') < 0.005) local skip = 1
			if (`skip' == 0) {
				local all_cuts "`all_cuts' `cutvalue'"
				local all_colors "`all_colors' black"
				local all_patterns "`all_patterns' solid"
				local all_widths "`all_widths' thin"
			}
		}
		local all_cuts "`all_cuts' `h0_q'"
		local all_colors "`all_colors' red"
		local all_patterns "`all_patterns' dash"
		local all_widths "`all_widths' medthick"

		local bench_scatter ""
		if ("`benchmark'" != "") {
			local dim = rowsof(__ivsm_bounds_iv)
			forvalues idx = 1(1)`dim' {
				local r2yz_val = __ivsm_bounds_iv[`idx', 4]
				local r2dz_val = __ivsm_bounds_iv[`idx', 3]
				local coef_val = __ivsm_bounds_iv[`idx', 6]
				local coef_val : di %5.3f `coef_val'
				local blbl "`__ivsm_blbl_`idx''"
				if (`r2yz_val' < `lim_ub' & `r2dz_val' < `lim_ub') {
					local bench_scatter "`bench_scatter' || scatteri `r2yz_val' `r2dz_val' `"`blbl' (`coef_val')"', mcolor(red) mfcolor(red) msize(small) msymbol(D) mlabcolor(black) mlabsize(vsmall) mlabposition(3)"
				}
			}
		}

		twoway contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
			name(__ivsm_upr_plot, replace) ccuts(`all_cuts') ccolor(`all_colors') clpattern(`all_patterns') clwidths(`all_widths') nodraw ///
			xlab(, labsize(small)) ylab(, labsize(small)) ///
			ytitle("Partial R{superscript:2} of confounder(s) with pot. outcome", size(small)) ///
			xtitle("Partial R{superscript:2} of confounder(s) with instrument", size(small)) ///
			title("Upper CI Limit", size(medium)) ///
			legend(off) ///
			|| scatteri 0 0 "Observed (`round_obs')", mcolor(black) msize(small) msymbol(T) mlabcolor(black) mlabsize(vsmall) mlabposition(3) ///
			`bench_scatter'

		capture: drop sense_contour_*
		capture: mat drop toprange_upr

		// Combine side by side
		graph combine __ivsm_lwr_plot __ivsm_upr_plot, name(__ivsm_contour, replace) nodraw
		graph display __ivsm_contour

		capture: graph drop __ivsm_lwr_plot
		capture: graph drop __ivsm_upr_plot

		restore
	}

	if ("`tcontourplot'" != "") {
		preserve

		mata: __ivsm_contour_plot_t(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', ///
			`lim_lb', `lim_ub', `alpha', `clines', `h0_q')

		capture: graph drop __ivsm_t_plot

		local t_crit = sign(`iv_estimate') * invttail(`ar_dof' - 1, `alpha'/2)
		local round_h0 : di %6.3f `h0_q'
		local obs_val = scalar(__ivsm_obs_value)
		local round_obs : di %5.3f `obs_val'

		local all_cuts ""
		local all_colors ""
		local all_patterns ""
		local all_widths ""
		local ncuts = rowsof(toprange_t)
		forvalues i = 1(1)`ncuts' {
			local cutvalue = toprange_t[`i', 1]
			local skip = 0
			if (abs(`cutvalue' - `t_crit') < 0.05) local skip = 1
			if (`skip' == 0) {
				local all_cuts "`all_cuts' `cutvalue'"
				local all_colors "`all_colors' black"
				local all_patterns "`all_patterns' solid"
				local all_widths "`all_widths' thin"
			}
		}
		local all_cuts "`all_cuts' `t_crit'"
		local all_colors "`all_colors' red"
		local all_patterns "`all_patterns' dash"
		local all_widths "`all_widths' medthick"

		// Build scatter overlay for benchmarks
		local bench_scatter ""
		if ("`benchmark'" != "") {
			local dim = rowsof(__ivsm_bounds_iv)
			forvalues idx = 1(1)`dim' {
				local r2yz_val = __ivsm_bounds_iv[`idx', 4]
				local r2dz_val = __ivsm_bounds_iv[`idx', 3]
				local blbl "`__ivsm_blbl_`idx''"
				local bench_lwr = __ivsm_bounds_iv[`idx', 5]
				local bench_upr = __ivsm_bounds_iv[`idx', 6]
				local bench_t : di %5.3f (`bench_lwr' + `bench_upr') / 2
				if (`r2yz_val' < `lim_ub' & `r2dz_val' < `lim_ub') {
					local bench_scatter "`bench_scatter' || scatteri `r2yz_val' `r2dz_val' `"`blbl' (`bench_t')"', mcolor(red) mfcolor(red) msize(small) msymbol(D) mlabcolor(black) mlabsize(vsmall) mlabposition(3)"
				}
			}
		}

		twoway contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
			name(__ivsm_t_plot, replace) ccuts(`all_cuts') ccolor(`all_colors') clpattern(`all_patterns') clwidths(`all_widths') nodraw ///
			xlab(, labsize(small)) ylab(, labsize(small)) ///
			ytitle("Partial R{superscript:2} of confounder(s) with pot. outcome", size(small)) ///
			xtitle("Partial R{superscript:2} of confounder(s) with instrument", size(small)) ///
			title("t-value for H0 = `round_h0'", size(medium)) ///
			legend(off) ///
			|| scatteri 0 0 "Observed (`round_obs')", mcolor(black) msize(small) msymbol(T) mlabcolor(black) mlabsize(vsmall) mlabposition(3) ///
			`bench_scatter'

		graph display __ivsm_t_plot

		capture: drop sense_contour_*
		capture: mat drop toprange_t

		restore
	}
}

// =========================================================================
// Step 8: ereturn storage
// =========================================================================

ereturn post, esample(`touse')

// Scalars
ereturn scalar iv_estimate = `iv_estimate'
ereturn scalar fs_coef = `fs_coef'
ereturn scalar fs_se = `fs_se'
ereturn scalar rf_coef = `rf_coef'
ereturn scalar rf_se = `rf_se'
ereturn scalar rho = `rho'
ereturn scalar ar_t = `ar_t'
ereturn scalar dof = `ar_dof'
ereturn scalar rv_iv = `rv_iv'
ereturn scalar xrv_iv = `xrv_iv'
ereturn scalar rv_fs = `rv_fs'
ereturn scalar xrv_fs = `xrv_fs'
ereturn scalar rv_rf = `rv_rf'
ereturn scalar xrv_rf = `xrv_rf'
ereturn scalar q = `q'
ereturn scalar alpha = `alpha'
ereturn scalar h0 = `h0'

// CI info
if ("`ar_ci_type'" == "bounded") {
	ereturn scalar ci_lwr = `ar_ci_lwr'
	ereturn scalar ci_upr = `ar_ci_upr'
}
else if ("`ar_ci_type'" == "disjoint") {
	ereturn scalar ci_lwr = `ar_ci_lwr'
	ereturn scalar ci_upr = `ar_ci_upr'
}
ereturn local ci_type = "`ar_ci_type'"

// Macros
ereturn local outcome    = "`depvar'"
ereturn local treatment  = "`treatvar'"
ereturn local instrument = "`instrvar'"
ereturn local covariates = "`covars'"

if ("`benchmark'" != "") {
	ereturn local bench = "`benchmark'"
	ereturn matrix bounds_iv = __ivsm_bounds_iv
	ereturn matrix bounds_fs = __ivsm_bounds_fs
	ereturn matrix bounds_rf = __ivsm_bounds_rf
}

if ("`contourplot'" != "") {
	capture: ereturn matrix contourgrid_lwr = __ivsm_contourgrid_lwr
	capture: ereturn matrix contourgrid_upr = __ivsm_contourgrid_upr
}

// Clean up
capture: estimates drop __ivsm_fs
capture: estimates drop __ivsm_rf
capture: estimates drop __ivsm_z_model

end


// =========================================================================
// Addplot helper (modified from Ben Jann)
// =========================================================================
program __ivsm_addplot
*! modified version of program written by Ben Jann
	version 9
	local caller : di _caller()
	_on_colon_parse `0'
	local plots `"`s(before)'"'
	local cmd   `"`s(after)'"'
	gettoken name : plots
	capt confirm name `name'
	if _rc == 0 {
		gettoken name plots : plots
	}
	else local name
	if `"`plots'"' != "" {
		numlist `"`plots'"', integer range(>=1)
		local plots `r(numlist)'
	}
	if `"`name'"' == "" {
		gr_current name :
	}
	else {
		capt classutil d `name'
	}
	if `"`plots'"' == "" {
		local grtype `"`.`name'.graphfamily'"'
		if `"`grtype'"' == "twoway" {
			version `caller': .`name'.parse `cmd'
			exit
		}
		local nplots `"`.`name'.n'"'
		capt confirm integer number `nplots'
		if _rc == 0 {
			forv plot = 1/`nplots' {
				if `"`.`name'.graphs[`plot'].graphfamily'"' == "twoway" {
					local plots `plots' `plot'
				}
			}
		}
	}
	else {
		foreach plot of local plots {
			capt classutil d `name'.graphs[`plot']
			local grtype `"`.`name'.graphs[`plot'].graphfamily'"'
		}
	}
	foreach plot of local plots {
		version `caller': .`name'.graphs[`plot'].parse `cmd'
	}
end


// =========================================================================
// MATA FUNCTIONS
// =========================================================================

version 13
mata:
mata set matastrict on
mata set matafavor speed


// -----------------------------------------------------------------
// ar_confint: Anderson-Rubin confidence interval
// Maps from: iv_functions.R:135-167
// -----------------------------------------------------------------
void __ivsm_ar_confint(real scalar fs_coef, real scalar fs_se,
					   real scalar rf_coef, real scalar rf_se,
					   real scalar rho, real scalar dof,
					   real scalar alpha) {
	real scalar ts, a, b, c, delta
	real scalar root1, root2, lo, hi

	// Critical value: sqrt(F_{1,dof}(alpha))
	ts = sqrt(invFtail(1, dof, alpha))

	// Quadratic coefficients
	a = fs_coef^2 - fs_se^2 * ts^2
	b = 2 * rho * rf_se * fs_se * ts^2 - 2 * rf_coef * fs_coef
	c = rf_coef^2 - rf_se^2 * ts^2
	delta = b^2 - 4 * a * c

	if (a < 0 & delta < 0) {
		// All reals: (-Inf, +Inf)
		st_local("ar_ci_type", "all_reals")
		st_local("ar_ci_lwr", ".")
		st_local("ar_ci_upr", ".")
		return
	}

	if (a < 0 & delta > 0) {
		// Disjoint: (-Inf, root1] U [root2, +Inf)
		root1 = (-b + sqrt(delta)) / (2 * a)
		root2 = (-b - sqrt(delta)) / (2 * a)
		lo = min((root1, root2))
		hi = max((root1, root2))
		st_local("ar_ci_type", "disjoint")
		st_local("ar_ci_lwr", strofreal(lo))
		st_local("ar_ci_upr", strofreal(hi))
		return
	}

	if (a > 0 & delta < 0) {
		// Empty set
		st_local("ar_ci_type", "empty")
		st_local("ar_ci_lwr", ".")
		st_local("ar_ci_upr", ".")
		return
	}

	if (a > 0 & delta >= 0) {
		// Bounded interval
		root1 = (-b + sqrt(delta)) / (2 * a)
		root2 = (-b - sqrt(delta)) / (2 * a)
		lo = min((root1, root2))
		hi = max((root1, root2))
		st_local("ar_ci_type", "bounded")
		st_local("ar_ci_lwr", strofreal(lo))
		st_local("ar_ci_upr", strofreal(hi))
		return
	}

	// Edge case: a == 0
	st_local("ar_ci_type", "all_reals")
	st_local("ar_ci_lwr", ".")
	st_local("ar_ci_upr", ".")
}


// -----------------------------------------------------------------
// xrv_ols: Extreme robustness value for OLS
// Maps from: sensemakr/sensitivity_stats.R:310-345
// -----------------------------------------------------------------
real scalar __ivsm_xrv_ols(real scalar t_stat, real scalar dof,
						   real scalar q, real scalar alpha) {
	real scalar fq2, f_crit2

	fq2     = (q * abs(t_stat) / sqrt(dof))^2
	f_crit2 = (abs(invttail(dof - 1, alpha / 2)) / sqrt(dof - 1))^2

	if (fq2 <= f_crit2) return(0)
	return((fq2 - f_crit2) / (1 + fq2))
}


// -----------------------------------------------------------------
// rv_ols: Robustness value for OLS
// Maps from: sensemakr/sensitivity_stats.R:162-215
// -----------------------------------------------------------------
real scalar __ivsm_rv_ols(real scalar t_stat, real scalar dof,
						  real scalar q, real scalar alpha) {
	real scalar fq, f_crit, fqa, rv, xrv_val

	fq     = q * abs(t_stat) / sqrt(dof)
	f_crit = abs(invttail(dof - 1, alpha / 2)) / sqrt(dof - 1)
	fqa    = fq - f_crit

	if (fqa < 0) return(0)

	// Numerically stable formula
	rv = 2 / (1 + sqrt(1 + 4 / fqa^2))

	// Check constraint not binding case
	xrv_val = __ivsm_xrv_ols(t_stat, dof, q, alpha)
	if (fqa > 0 & fq > 1 / f_crit) {
		return(xrv_val)
	}

	return(rv)
}


// -----------------------------------------------------------------
// rv_iv: Robustness value for IV estimate
// Maps from: iv.sensemakr/sensitivity_stats.R:99-134
// -----------------------------------------------------------------
real scalar __ivsm_rv_iv(real scalar fs_coef, real scalar fs_se,
						 real scalar rf_coef, real scalar rf_se,
						 real scalar rho, real scalar dof,
						 real scalar q, real scalar alpha,
						 real scalar use_min) {
	real scalar tau_r, tau0, phi_r, se_phi_r, t_phi_r, t_fs_r
	real scalar rv_phi, rv_fs

	tau_r    = rf_coef / fs_coef
	tau0     = (1 - q) * tau_r
	phi_r    = rf_coef - tau0 * fs_coef
	se_phi_r = sqrt(rf_se^2 + tau0^2 * fs_se^2 - 2 * tau0 * rho * fs_se * rf_se)
	t_phi_r  = abs(phi_r / se_phi_r)
	t_fs_r   = abs(fs_coef / fs_se)

	rv_phi = __ivsm_rv_ols(t_phi_r, dof, 1, alpha)
	rv_fs  = __ivsm_rv_ols(t_fs_r, dof, 1, alpha)

	if (use_min) {
		return(min((rv_phi, rv_fs)))
	}
	return(rv_phi)
}


// -----------------------------------------------------------------
// xrv_iv: Extreme robustness value for IV estimate
// -----------------------------------------------------------------
real scalar __ivsm_xrv_iv(real scalar fs_coef, real scalar fs_se,
						  real scalar rf_coef, real scalar rf_se,
						  real scalar rho, real scalar dof,
						  real scalar q, real scalar alpha,
						  real scalar use_min) {
	real scalar tau_r, tau0, phi_r, se_phi_r, t_phi_r, t_fs_r
	real scalar xrv_phi, xrv_fs

	tau_r    = rf_coef / fs_coef
	tau0     = (1 - q) * tau_r
	phi_r    = rf_coef - tau0 * fs_coef
	se_phi_r = sqrt(rf_se^2 + tau0^2 * fs_se^2 - 2 * tau0 * rho * fs_se * rf_se)
	t_phi_r  = abs(phi_r / se_phi_r)
	t_fs_r   = abs(fs_coef / fs_se)

	xrv_phi = __ivsm_xrv_ols(t_phi_r, dof, 1, alpha)
	xrv_fs  = __ivsm_xrv_ols(t_fs_r, dof, 1, alpha)

	if (use_min) {
		return(min((xrv_phi, xrv_fs)))
	}
	return(xrv_phi)
}


// -----------------------------------------------------------------
// adjusted_critical_value
// Maps from: sensemakr/bias_functions.R:32-56
// -----------------------------------------------------------------
real scalar __ivsm_adjusted_critical_value(real scalar r2dz_x,
										   real scalar r2yz_dx,
										   real scalar dof,
										   real scalar alpha) {
	real scalar t_crit, f_crit, sef_val, bf_val, ddof, t_dagger

	// Traditional critical value
	t_crit = sqrt(invFtail(1, dof - 1, alpha))
	f_crit = t_crit / sqrt(dof)

	// Max adjustment (worst-case)
	if (r2dz_x < f_crit^2 * (r2yz_dx / (1 - r2yz_dx))) {
		r2yz_dx = r2dz_x / (f_crit^2 + r2dz_x)
	}

	// SE factor and bias factor
	sef_val = sqrt(1 - r2yz_dx) / sqrt(1 - r2dz_x)
	bf_val  = sqrt(r2yz_dx * r2dz_x / (1 - r2dz_x))

	ddof = sqrt(dof / (dof - 1))
	t_dagger = t_crit * sef_val * ddof + bf_val * sqrt(dof)

	return(t_dagger)
}


// -----------------------------------------------------------------
// abc_builder: quadratic coefficients for IV CI
// Maps from: plots.R:841-858
// -----------------------------------------------------------------
real vector __ivsm_abc_builder(real scalar fs_coef, real scalar fs_se,
							   real scalar rf_coef, real scalar rf_se,
							   real scalar rho, real scalar crit_thr) {
	real scalar a, b, c

	a = fs_coef^2 - fs_se^2 * crit_thr^2
	b = 2 * rho * rf_se * fs_se * crit_thr^2 - 2 * rf_coef * fs_coef
	c = rf_coef^2 - rf_se^2 * crit_thr^2

	return((a, b, c))
}


// -----------------------------------------------------------------
// quad_lwr_limit / quad_upr_limit
// Maps from: plots.R:725-746
// -----------------------------------------------------------------
real scalar __ivsm_quad_lwr_limit(real scalar a, real scalar b, real scalar c) {
	real scalar delta
	delta = b^2 - 4 * a * c
	if (a <= 0) return(.)
	if (a > 0 & delta >= 0) return((-b - sqrt(delta)) / (2 * a))
	return(.)
}

real scalar __ivsm_quad_upr_limit(real scalar a, real scalar b, real scalar c) {
	real scalar delta
	delta = b^2 - 4 * a * c
	if (a < 0) return(.)
	if (a > 0 & delta > 0) return((-b + sqrt(delta)) / (2 * a))
	return(.)
}


// -----------------------------------------------------------------
// iv_adjusted_limit: adjusted CI limit at given R2 values
// Maps from: plots.R:566-601
// -----------------------------------------------------------------
real scalar __ivsm_iv_adjusted_limit(real scalar fs_coef, real scalar fs_se,
									 real scalar rf_coef, real scalar rf_se,
									 real scalar rho, real scalar dof,
									 real scalar r2zw_x, real scalar r2y0w_zx,
									 real scalar alpha,
									 string scalar ci_limit) {
	real scalar t_dagger
	real vector abc

	t_dagger = __ivsm_adjusted_critical_value(r2zw_x, r2y0w_zx, dof, alpha)
	abc = __ivsm_abc_builder(fs_coef, fs_se, rf_coef, rf_se, rho, t_dagger)

	if (ci_limit == "lwr") {
		return(__ivsm_quad_lwr_limit(abc[1], abc[2], abc[3]))
	}
	else {
		return(__ivsm_quad_upr_limit(abc[1], abc[2], abc[3]))
	}
}


// -----------------------------------------------------------------
// maxR2: find max cor(y_zx - tau0*d_zx, xj_zx)^2 over tau0
// Maps from: plots.R:710-717
// Uses golden-section search
// -----------------------------------------------------------------
real scalar __ivsm_maxR2(real colvector y_zx, real colvector d_zx,
						 real colvector xj_zx) {
	real scalar a, b, gr, c_pt, d_pt, fc, fd
	real scalar i, tau0, best_r2
	real colvector ytau0
	real matrix corr_mat

	a = -1e5
	b = 1e5
	gr = (sqrt(5) + 1) / 2

	for (i = 1; i <= 100; i++) {
		c_pt = b - (b - a) / gr
		d_pt = a + (b - a) / gr

		ytau0 = y_zx :- c_pt :* d_zx
		corr_mat = correlation((ytau0, xj_zx))
		fc = corr_mat[1, 2]^2

		ytau0 = y_zx :- d_pt :* d_zx
		corr_mat = correlation((ytau0, xj_zx))
		fd = corr_mat[1, 2]^2

		if (fc > fd) {
			b = d_pt
		}
		else {
			a = c_pt
		}
	}

	tau0 = (a + b) / 2
	ytau0 = y_zx :- tau0 :* d_zx
	corr_mat = correlation((ytau0, xj_zx))
	best_r2 = corr_mat[1, 2]^2

	return(best_r2)
}


// -----------------------------------------------------------------
// ovb_partial_r2_bound: compute R2 bounds from benchmark values
// Maps from: sensemakr ovb_bounds formula
// Returns (r2dz_x, r2yz_dx) as a row vector
// -----------------------------------------------------------------
real vector __ivsm_ovb_partial_r2_bound(real scalar r2dxj_x, real scalar r2yxj_x,
										real scalar kz, real scalar ky) {
	real scalar r2dz_x, r2yz_dx, r2zxj_xd

	r2dz_x = kz * (r2dxj_x / (1 - r2dxj_x))
	if (r2dz_x > 1) r2dz_x = 1

	r2zxj_xd = kz * (r2dxj_x^2) / ((1 - kz * r2dxj_x) * (1 - r2dxj_x))
	if (r2zxj_xd > 1) r2zxj_xd = 1

	r2yz_dx = ((sqrt(ky) + sqrt(r2zxj_xd)) / sqrt(1 - r2zxj_xd))^2 * (r2yxj_x / (1 - r2yxj_x))
	if (r2yz_dx > 1) r2yz_dx = 1

	return((r2dz_x, r2yz_dx))
}


// -----------------------------------------------------------------
// compute_iv_bounds: compute IV bounds for a single benchmark x k
// -----------------------------------------------------------------
void __ivsm_compute_iv_bounds(real scalar r2dxj_x, real scalar r2yxj_x,
							  real scalar kz, real scalar ky,
							  real scalar fs_coef, real scalar fs_se,
							  real scalar rf_coef, real scalar rf_se,
							  real scalar rho, real scalar dof,
							  real scalar alpha) {
	real scalar r2zw_x, r2y0w_zx, lwr, upr
	real vector bounds

	bounds = __ivsm_ovb_partial_r2_bound(r2dxj_x, r2yxj_x, kz, ky)
	r2zw_x   = bounds[1]
	r2y0w_zx = bounds[2]

	lwr = __ivsm_iv_adjusted_limit(fs_coef, fs_se, rf_coef, rf_se, rho, dof,
									r2zw_x, r2y0w_zx, alpha, "lwr")
	upr = __ivsm_iv_adjusted_limit(fs_coef, fs_se, rf_coef, rf_se, rho, dof,
									r2zw_x, r2y0w_zx, alpha, "upr")

	st_local("bound_r2zw_x", strofreal(r2zw_x))
	st_local("bound_r2y0w_zx", strofreal(r2y0w_zx))
	st_local("bound_lwr", strofreal(lwr))
	st_local("bound_upr", strofreal(upr))
}


// -----------------------------------------------------------------
// compute_ols_bounds: compute standard OLS bounds for FS/RF
// Maps from: sensemakr-stata iterate_bounds logic
// -----------------------------------------------------------------
void __ivsm_compute_ols_bounds(real scalar r2dxj_x, real scalar r2yxj_x,
							   real scalar kd_val, real scalar ky_val,
							   real scalar coef, real scalar se,
							   real scalar dof, real scalar alpha) {
	real scalar r2dz_x, r2zxj_xd, r2yz_dx
	real scalar t_dagger, lwr, upr

	r2dz_x = kd_val * (r2dxj_x / (1 - r2dxj_x))
	if (r2dz_x > 1) r2dz_x = 1

	r2zxj_xd = kd_val * (r2dxj_x^2) / ((1 - kd_val * r2dxj_x) * (1 - r2dxj_x))
	if (r2zxj_xd > 1) r2zxj_xd = 1

	r2yz_dx = ((sqrt(ky_val) + sqrt(r2zxj_xd)) / sqrt(1 - r2zxj_xd))^2 * (r2yxj_x / (1 - r2yxj_x))
	if (r2yz_dx > 1) r2yz_dx = 1

	// Adjusted critical value (t-dagger), then CI = estimate +/- t_dagger * se
	// This matches the R code: sensemakr.R lines 136-137
	t_dagger = __ivsm_adjusted_critical_value(r2dz_x, r2yz_dx, dof, alpha)
	lwr = coef - t_dagger * se
	upr = coef + t_dagger * se

	st_local("ols_r2dz_x", strofreal(r2dz_x))
	st_local("ols_r2yz_dx", strofreal(r2yz_dx))
	st_local("ols_lwr", strofreal(lwr))
	st_local("ols_upr", strofreal(upr))
}


// -----------------------------------------------------------------
// contour_plot_iv: generate contour grid for CI limits (vectorized)
// Maps from: plots.R ovb4iv_contour_plot.numeric
// -----------------------------------------------------------------
void __ivsm_contour_plot_iv(real scalar fs_coef, real scalar fs_se,
							real scalar rf_coef, real scalar rf_se,
							real scalar rho, real scalar dof,
							real scalar lim_lb, real scalar lim_ub,
							real scalar alpha,
							string scalar ci_limit,
							real scalar clines,
							real scalar threshold) {
	real colvector grid_x, grid_y, r2x_vec, r2y_vec
	real matrix contour_grid
	real scalar i, j, z, nx, ny, n_existing, n_total, val, tc, fc
	string scalar mat_name, toprange_name
	real colvector sef_vec, bf_vec, t_dagger
	real colvector a_vec, b_vec, c_vec, delta_vec, result
	real colvector r2y_adj, cap_mask

	nx = 50
	ny = 50
	n_total = nx * ny

	grid_x = rangen(lim_lb, lim_ub, nx)
	grid_y = rangen(lim_lb, lim_ub, ny)

	// Build flat grid vectors using Kronecker-style expansion (no loop)
	r2x_vec = J(ny, 1, 1) # grid_x  // each x repeated ny times
	r2y_vec = J(1, 1, 1) # J(nx, 1, 1) # grid_y[ny::1]
	// Actually use simple loop — Mata Kronecker may not work this way
	r2x_vec = J(n_total, 1, .)
	r2y_vec = J(n_total, 1, .)
	z = 1
	for (i = 1; i <= nx; i++) {
		for (j = ny; j >= 1; j--) {
			r2x_vec[z] = grid_x[i]
			r2y_vec[z] = grid_y[j]
			z++
		}
	}

	// Vectorized adjusted_critical_value (no loops)
	tc = sqrt(invFtail(1, dof - 1, alpha))
	fc = tc / sqrt(dof)

	// Cap r2y where r2x < fc^2 * (r2y/(1-r2y))
	cap_mask = r2x_vec :< (fc^2 :* (r2y_vec :/ (1 :- r2y_vec)))
	r2y_adj = r2y_vec
	r2y_adj = cap_mask :* (r2x_vec :/ (fc^2 :+ r2x_vec)) + (1 :- cap_mask) :* r2y_vec

	sef_vec = sqrt(1 :- r2y_adj) :/ sqrt(1 :- r2x_vec)
	bf_vec  = sqrt(r2y_adj :* r2x_vec :/ (1 :- r2x_vec))
	t_dagger = tc :* sef_vec :* sqrt(dof / (dof - 1)) :+ bf_vec :* sqrt(dof)

	// Vectorized abc_builder
	a_vec = fs_coef^2 :- fs_se^2 :* t_dagger:^2
	b_vec = 2 * rho * rf_se * fs_se :* t_dagger:^2 :- 2 * rf_coef * fs_coef
	c_vec = rf_coef^2 :- rf_se^2 :* t_dagger:^2

	// Vectorized quadratic solver (no loops)
	delta_vec = b_vec:^2 :- 4 :* a_vec :* c_vec
	result = J(n_total, 1, .)

	if (ci_limit == "lwr") {
		// lwr: a > 0 & delta >= 0
		for (i = 1; i <= n_total; i++) {
			if (a_vec[i] > 0 & delta_vec[i] >= 0) {
				result[i] = (-b_vec[i] - sqrt(delta_vec[i])) / (2 * a_vec[i])
			}
		}
	}
	else {
		// upr: a > 0 & delta > 0
		for (i = 1; i <= n_total; i++) {
			if (a_vec[i] > 0 & delta_vec[i] > 0) {
				result[i] = (-b_vec[i] + sqrt(delta_vec[i])) / (2 * a_vec[i])
			}
		}
	}

	contour_grid = (r2x_vec, r2y_vec, result)

	// Store matrix
	if (ci_limit == "lwr") {
		mat_name = "__ivsm_contourgrid_lwr"
		toprange_name = "toprange_lwr"
	}
	else {
		mat_name = "__ivsm_contourgrid_upr"
		toprange_name = "toprange_upr"
	}
	st_matrix(mat_name, contour_grid)

	// Create Stata variables for plotting
	n_existing = st_nobs()
	(void) st_addobs(n_total)
	(void) st_addvar("double", "sense_contour_x")
	(void) st_addvar("double", "sense_contour_y")
	(void) st_addvar("double", "sense_contour_z")

	st_store((n_existing + 1::n_existing + n_total), "sense_contour_x", contour_grid[., 1])
	st_store((n_existing + 1::n_existing + n_total), "sense_contour_y", contour_grid[., 2])
	st_store((n_existing + 1::n_existing + n_total), "sense_contour_z", contour_grid[., 3])

	// Compute contour levels using quantiles (matching R)
	real colvector finite_z, quantile_probs, raw_levels, cuts
	real scalar clines_actual, n_finite

	clines_actual = min((clines, nx - 2))

	finite_z = select(contour_grid[., 3], contour_grid[., 3] :< .)
	n_finite = length(finite_z)
	if (n_finite > 0) {
		finite_z = sort(finite_z, 1)
		quantile_probs = rangen(0, 1, max((2, clines_actual)))
		raw_levels = J(length(quantile_probs), 1, .)
		for (i = 1; i <= length(quantile_probs); i++) {
			raw_levels[i] = finite_z[max((1, ceil(quantile_probs[i] * n_finite)))]
		}
		raw_levels[1] = finite_z[max((1, ceil(0.01 * n_finite)))]
		raw_levels[length(raw_levels)] = finite_z[max((1, ceil(0.99 * n_finite)))]

		raw_levels = round(raw_levels, .001)
		raw_levels = raw_levels \ round(threshold, .001)
		raw_levels = sort(uniqrows(raw_levels), 1)
		cuts = raw_levels
	}
	else {
		cuts = J(1, 1, threshold)
	}
	st_matrix(toprange_name, cuts)

	// Store the observed value at (0,0) for labeling
	val = __ivsm_iv_adjusted_limit(fs_coef, fs_se, rf_coef, rf_se,
								   rho, dof, 0, 0, alpha, ci_limit)
	st_numscalar("__ivsm_obs_value", val)
}


// -----------------------------------------------------------------
// contour_plot_t: generate contour grid for t-values (AR approach)
// Uses OLS adjusted_t on the Anderson-Rubin regression
// -----------------------------------------------------------------
void __ivsm_contour_plot_t(real scalar fs_coef, real scalar fs_se,
						   real scalar rf_coef, real scalar rf_se,
						   real scalar rho, real scalar dof,
						   real scalar lim_lb, real scalar lim_ub,
						   real scalar alpha,
						   real scalar clines,
						   real scalar h0) {
	real colvector grid_x, grid_y, r2x_vec, r2y_vec
	real matrix contour_grid
	real scalar i, j, z, nx, ny, n_existing, n_total
	real scalar ar_coef, ar_se, clines_actual
	real colvector bias_vec, adj_e_vec, adj_se_vec, result

	ar_coef = rf_coef - h0 * fs_coef
	ar_se = sqrt(rf_se^2 + h0^2 * fs_se^2 - 2 * h0 * rho * fs_se * rf_se)

	grid_x = rangen(lim_lb, lim_ub, 100)
	grid_y = rangen(lim_lb, lim_ub, 100)
	nx = length(grid_x)
	ny = length(grid_y)
	n_total = nx * ny

	// Build flat grid vectors
	r2x_vec = J(n_total, 1, .)
	r2y_vec = J(n_total, 1, .)
	z = 1
	for (i = 1; i <= nx; i++) {
		for (j = ny; j >= 1; j--) {
			r2x_vec[z] = grid_x[i]
			r2y_vec[z] = grid_y[j]
			z++
		}
	}

	// Vectorized OLS-style adjusted t-value
	bias_vec = sqrt(r2y_vec :* r2x_vec :/ (1 :- r2x_vec)) :* (ar_se * sqrt(dof))
	adj_e_vec = sign(ar_coef) :* (abs(ar_coef) :- bias_vec)
	adj_se_vec = sqrt((1 :- r2y_vec) :/ (1 :- r2x_vec)) :* (ar_se * sqrt(dof / (dof - 1)))

	result = J(n_total, 1, .)
	for (i = 1; i <= n_total; i++) {
		if (adj_se_vec[i] > 0) {
			result[i] = adj_e_vec[i] / adj_se_vec[i]
		}
	}

	contour_grid = (r2x_vec, r2y_vec, result)
	st_matrix("__ivsm_contourgrid_t", contour_grid)

	n_existing = st_nobs()
	(void) st_addobs(n_total)
	(void) st_addvar("double", "sense_contour_x")
	(void) st_addvar("double", "sense_contour_y")
	(void) st_addvar("double", "sense_contour_z")

	st_store((n_existing + 1::n_existing + n_total), "sense_contour_x", contour_grid[., 1])
	st_store((n_existing + 1::n_existing + n_total), "sense_contour_y", contour_grid[., 2])
	st_store((n_existing + 1::n_existing + n_total), "sense_contour_z", contour_grid[., 3])

	// Compute contour levels using quantiles (matching R)
	real colvector finite_z, quantile_probs, raw_levels, cuts
	real scalar n_finite, t_crit

	clines_actual = min((clines, nx - 2))
	t_crit = sign(ar_coef) * invttail(dof - 1, alpha/2)

	finite_z = select(contour_grid[., 3], contour_grid[., 3] :< .)
	n_finite = length(finite_z)
	if (n_finite > 0) {
		finite_z = sort(finite_z, 1)
		quantile_probs = rangen(0, 1, max((2, clines_actual)))
		raw_levels = J(length(quantile_probs), 1, .)
		for (i = 1; i <= length(quantile_probs); i++) {
			raw_levels[i] = finite_z[max((1, ceil(quantile_probs[i] * n_finite)))]
		}
		raw_levels[1] = finite_z[max((1, ceil(0.01 * n_finite)))]
		raw_levels[length(raw_levels)] = finite_z[max((1, ceil(0.99 * n_finite)))]

		raw_levels = round(raw_levels, .001)
		raw_levels = raw_levels \ round(t_crit, .001)
		raw_levels = sort(uniqrows(raw_levels), 1)
		cuts = raw_levels
	}
	else {
		cuts = J(1, 1, t_crit)
	}
	st_matrix("toprange_t", cuts)

	// Store the observed t-value at (0,0)
	st_numscalar("__ivsm_obs_value", ar_coef / (ar_se * sqrt(dof / (dof - 1))))
}


end
