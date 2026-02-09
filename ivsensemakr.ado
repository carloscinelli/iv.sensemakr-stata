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

// AR CI at the adjusted h0
mata: __ivsm_ar_confint(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `alpha')

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
				local blabel = "`kz_val'x `bench'"
			}
			else {
				local blabel = "`kz_val'/`ky_val'x `bench'"
			}

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

		// Contour plot - lower limit
		capture: graph drop __ivsm_lwr_plot
		local round_est : di %6.3f `iv_estimate'
		local round_h0 : di %6.3f `h0_q'

		twoway contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
			name(__ivsm_lwr_plot, replace) ccuts(`h0_q') ccolor(red) nodraw ///
			xlab(, labsize(small)) ylab(, labsize(small)) ///
			ytitle("Partial R{superscript:2} of confounder(s) with pot. outcome", size(small)) ///
			xtitle("Partial R{superscript:2} of confounder(s) with instrument", size(small)) ///
			title("Lower CI Limit", size(medium)) ///
			text(0.001 0.001 "Unadjusted (`round_est')", size(vsmall) place(ne)) ///
			legend(off) ///
			|| scatteri 0 0, mcolor(black) msize(small) msymbol(T)

		// Add extra contour lines
		forvalues i = 1(1)`clines' {
			local cutvalue = toprange_lwr[`i', 3]
			local cutvalue_round : di %6.3f toprange_lwr[`i', 3]
			local cut_x = toprange_lwr[`i', 2] + .001
			capture: __ivsm_addplot __ivsm_lwr_plot: contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
				ccuts(`cutvalue') ccolor(gray) clwidths(vthin) ///
				text(`cut_x' `cut_x' "`cutvalue_round'", size(vsmall) bcolor(white) box place(ne))
		}

		// Add benchmark points if available
		if ("`benchmark'" != "") {
			local dim = rowsof(__ivsm_bounds_iv)
			local rownms : rown __ivsm_bounds_iv
			forvalues idx = 1(1)`dim' {
				local r2yz_val = __ivsm_bounds_iv[`idx', 4]
				local r2dz_val = __ivsm_bounds_iv[`idx', 3]
				local coef_val = __ivsm_bounds_iv[`idx', 5]
				local coef_val : di %5.3f `coef_val'
				local kz_v = __ivsm_bounds_iv[`idx', 1]
				local ky_v = __ivsm_bounds_iv[`idx', 2]
				if (`r2yz_val' < `lim_ub' & `r2dz_val' < `lim_ub') {
					capture: __ivsm_addplot __ivsm_lwr_plot: scatteri `r2yz_val' `r2dz_val' "`kz_v'x (`coef_val')", ///
						mcolor(black) msize(vsmall) mlabcolor(black) mlabsize(vsmall)
				}
			}
		}

		capture: drop sense_contour_*
		capture: mat drop toprange_lwr

		// Upper CI contour
		mata: __ivsm_contour_plot_iv(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', `lim_lb', `lim_ub', `alpha', "upr", `clines', `h0_q')
		matrix colnames __ivsm_contourgrid_upr = "r2zw_x" "r2y0w_zx" "upr_ci"

		capture: graph drop __ivsm_upr_plot

		twoway contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
			name(__ivsm_upr_plot, replace) ccuts(`h0_q') ccolor(red) nodraw ///
			xlab(, labsize(small)) ylab(, labsize(small)) ///
			ytitle("Partial R{superscript:2} of confounder(s) with pot. outcome", size(small)) ///
			xtitle("Partial R{superscript:2} of confounder(s) with instrument", size(small)) ///
			title("Upper CI Limit", size(medium)) ///
			text(0.001 0.001 "Unadjusted (`round_est')", size(vsmall) place(ne)) ///
			legend(off) ///
			|| scatteri 0 0, mcolor(black) msize(small) msymbol(T)

		forvalues i = 1(1)`clines' {
			local cutvalue = toprange_upr[`i', 3]
			local cutvalue_round : di %6.3f toprange_upr[`i', 3]
			local cut_x = toprange_upr[`i', 2] + .001
			capture: __ivsm_addplot __ivsm_upr_plot: contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
				ccuts(`cutvalue') ccolor(gray) clwidths(vthin) ///
				text(`cut_x' `cut_x' "`cutvalue_round'", size(vsmall) bcolor(white) box place(ne))
		}

		if ("`benchmark'" != "") {
			local dim = rowsof(__ivsm_bounds_iv)
			forvalues idx = 1(1)`dim' {
				local r2yz_val = __ivsm_bounds_iv[`idx', 4]
				local r2dz_val = __ivsm_bounds_iv[`idx', 3]
				local coef_val = __ivsm_bounds_iv[`idx', 6]
				local coef_val : di %5.3f `coef_val'
				local kz_v = __ivsm_bounds_iv[`idx', 1]
				if (`r2yz_val' < `lim_ub' & `r2dz_val' < `lim_ub') {
					capture: __ivsm_addplot __ivsm_upr_plot: scatteri `r2yz_val' `r2dz_val' "`kz_v'x (`coef_val')", ///
						mcolor(black) msize(vsmall) mlabcolor(black) mlabsize(vsmall)
				}
			}
		}

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

		// t-value contour uses the AR model (OLS-style adjusted t)
		mata: __ivsm_contour_plot_t(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', ///
			`lim_lb', `lim_ub', `alpha', `clines', `h0_q')

		capture: graph drop __ivsm_t_plot

		// Critical t-value threshold
		local t_crit = sign(`iv_estimate') * invttail(`ar_dof' - 1, `alpha'/2)
		local round_t : di %6.3f `ar_t_q'

		twoway contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
			name(__ivsm_t_plot, replace) ccuts(`t_crit') ccolor(red) nodraw ///
			xlab(, labsize(small)) ylab(, labsize(small)) ///
			ytitle("Partial R{superscript:2} of confounder(s) with pot. outcome", size(small)) ///
			xtitle("Partial R{superscript:2} of confounder(s) with instrument", size(small)) ///
			title("t-value for H0 = `round_h0'", size(medium)) ///
			text(0.001 0.001 "Unadjusted (`round_t')", size(vsmall) place(ne)) ///
			legend(off) ///
			|| scatteri 0 0, mcolor(black) msize(small) msymbol(T)

		forvalues i = 1(1)`clines' {
			local cutvalue = toprange_t[`i', 3]
			local cutvalue_round : di %6.3f toprange_t[`i', 3]
			local cut_x = toprange_t[`i', 2] + .001
			capture: __ivsm_addplot __ivsm_t_plot: contourline sense_contour_z sense_contour_y sense_contour_x if sense_contour_z != ., ///
				ccuts(`cutvalue') ccolor(gray) clwidths(vthin) ///
				text(`cut_x' `cut_x' "`cutvalue_round'", size(vsmall) bcolor(white) box place(ne))
		}

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
mata clear
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
// contour_plot_iv: generate contour grid for CI limits
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
	real colvector grid_x, grid_y
	real matrix contour_grid, cuts
	real scalar i, j, z, nx, ny, val
	string scalar mat_name, toprange_name

	grid_x = rangen(lim_lb, lim_ub, 51)
	grid_y = rangen(lim_lb, lim_ub, 51)
	nx = length(grid_x)
	ny = length(grid_y)

	contour_grid = J(nx * ny, 3, .)
	z = 1
	for (i = 1; i <= nx; i++) {
		for (j = ny; j >= 1; j--) {
			contour_grid[z, 1] = grid_x[i]
			contour_grid[z, 2] = grid_y[j]
			contour_grid[z, 3] = __ivsm_iv_adjusted_limit(fs_coef, fs_se, rf_coef, rf_se,
														   rho, dof, grid_x[i], grid_y[j],
														   alpha, ci_limit)
			z++
		}
	}

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
	(void) st_addobs(nx * ny)
	(void) st_addvar("double", "sense_contour_x")
	(void) st_addvar("double", "sense_contour_y")
	(void) st_addvar("double", "sense_contour_z")

	for (i = 1; i <= (nx * ny); i++) {
		st_store(i, "sense_contour_x", contour_grid[i, 1])
		st_store(i, "sense_contour_y", contour_grid[i, 2])
		st_store(i, "sense_contour_z", contour_grid[i, 3])
	}

	// Compute contour line cut values and label positions
	real colvector diag_vals
	real scalar clines_actual
	clines_actual = min((clines, nx - 2))
	diag_vals = rangen(lim_lb + (lim_ub - lim_lb) * 0.05, lim_ub * 0.95, clines_actual)
	cuts = J(clines_actual, 3, .)
	for (i = 1; i <= clines_actual; i++) {
		cuts[i, 1] = diag_vals[i]
		cuts[i, 2] = diag_vals[i]
		cuts[i, 3] = __ivsm_iv_adjusted_limit(fs_coef, fs_se, rf_coef, rf_se,
											   rho, dof, diag_vals[i], diag_vals[i],
											   alpha, ci_limit)
	}
	st_matrix(toprange_name, cuts)
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
	real colvector grid_x, grid_y
	real matrix contour_grid, cuts
	real scalar i, j, z, nx, ny
	real scalar ar_coef, ar_se, bias, adj_e, adj_se, adj_t
	real colvector diag_vals
	real scalar clines_actual

	// The AR regression tests H0: tau = h0
	// The coefficient is rf_coef - h0*fs_coef (approximately), se is composite
	// For t-value contour, we use the OLS adjusted_t approach on the AR model
	// The AR model: y - h0*d = z*phi + x*gamma + e
	// phi = rf_coef - h0*fs_coef
	// We compute the adjusted t-value treating the AR model as a standard OLS model
	ar_coef = rf_coef - h0 * fs_coef
	ar_se = sqrt(rf_se^2 + h0^2 * fs_se^2 - 2 * h0 * rho * fs_se * rf_se)

	grid_x = rangen(lim_lb, lim_ub, 51)
	grid_y = rangen(lim_lb, lim_ub, 51)
	nx = length(grid_x)
	ny = length(grid_y)

	contour_grid = J(nx * ny, 3, .)
	z = 1
	for (i = 1; i <= nx; i++) {
		for (j = ny; j >= 1; j--) {
			contour_grid[z, 1] = grid_x[i]
			contour_grid[z, 2] = grid_y[j]

			// OLS-style adjusted t-value
			bias = sqrt(grid_y[j] * grid_x[i] / (1 - grid_x[i])) * ar_se * sqrt(dof)
			adj_e = sign(ar_coef) * (abs(ar_coef) - bias)
			adj_se = sqrt((1 - grid_y[j]) / (1 - grid_x[i])) * ar_se * sqrt(dof / (dof - 1))
			if (adj_se > 0) {
				adj_t = adj_e / adj_se
			}
			else {
				adj_t = .
			}
			contour_grid[z, 3] = adj_t
			z++
		}
	}

	st_matrix("__ivsm_contourgrid_t", contour_grid)

	(void) st_addobs(nx * ny)
	(void) st_addvar("double", "sense_contour_x")
	(void) st_addvar("double", "sense_contour_y")
	(void) st_addvar("double", "sense_contour_z")

	for (i = 1; i <= (nx * ny); i++) {
		st_store(i, "sense_contour_x", contour_grid[i, 1])
		st_store(i, "sense_contour_y", contour_grid[i, 2])
		st_store(i, "sense_contour_z", contour_grid[i, 3])
	}

	clines_actual = min((clines, nx - 2))
	diag_vals = rangen(lim_lb + (lim_ub - lim_lb) * 0.05, lim_ub * 0.95, clines_actual)
	cuts = J(clines_actual, 3, .)
	for (i = 1; i <= clines_actual; i++) {
		cuts[i, 1] = diag_vals[i]
		cuts[i, 2] = diag_vals[i]
		bias = sqrt(diag_vals[i] * diag_vals[i] / (1 - diag_vals[i])) * ar_se * sqrt(dof)
		adj_e = sign(ar_coef) * (abs(ar_coef) - bias)
		adj_se = sqrt((1 - diag_vals[i]) / (1 - diag_vals[i])) * ar_se * sqrt(dof / (dof - 1))
		if (adj_se > 0) {
			cuts[i, 3] = adj_e / adj_se
		}
		else {
			cuts[i, 3] = .
		}
	}
	st_matrix("toprange_t", cuts)
}


end
