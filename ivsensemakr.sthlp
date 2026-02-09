{smcl}
{* *! version 1.0.0 07feb2026}

{cmd:help ivsensemakr}
{hline}

{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{hi:ivsensemakr} {hline 2} Sensitivity analysis for instrumental variable estimates}{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 18 2}
{cmdab:ivsensemakr} {it:depvar} {it:treatvar} {it:instrvar} [{it:covariates}]
{ifin}
[{cmd:,}
{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Main}
{synopt:{opt b:enchmark(varlist)}}specify benchmark covariates for bounding omitted variable bias{p_end}
{synopt:{opt kz(numlist)}}strength of confounder relative to benchmark in explaining the instrument. Default: 1{p_end}
{synopt:{opt ky(numlist)}}strength of confounder relative to benchmark in explaining the potential outcome. Default: kz{p_end}
{synopt:{opt kd(numlist)}}strength of confounder relative to benchmark in explaining the treatment. Default: kz{p_end}
{synopt:{opt r2zw(numlist)}}manual partial R2 of omitted variable(s) with the instrument{p_end}
{synopt:{opt r2y0w(numlist)}}manual partial R2 of omitted variable(s) with the potential outcome. Default: r2zw{p_end}
{synopt:{opt bound_label(string)}}label for manual bounds. Default: "Manual Bound"{p_end}
{synopt:{opt s:uppress}}suppress verbose output (first-stage and reduced-form details){p_end}

{syntab:Sensitivity Parameters}
{synopt:{opt q(real)}}fraction of the effect estimate to be explained away. Default: 1{p_end}
{synopt:{opt alpha(real)}}significance level. Default: 0.05{p_end}
{synopt:{opt h0(real)}}null hypothesis for the IV estimate. Default: 0{p_end}
{synopt:{opt nomin}}compute robustness value for bias of exactly q, rather than >= q{p_end}

{syntab:Graphing}
{synopt:{opt contourplot}}generate side-by-side contour plots for lower and upper CI limits{p_end}
{synopt:{opt tcontourplot}}generate a contour plot for the t-value{p_end}
{synopt:{opt clim(numlist)}}lower and upper limits for contour plot axes{p_end}
{synopt:{opt clines(real)}}number of contour lines. Default: 7{p_end}

{synoptline}
{p2colreset}{...}

{title:Description}

{pstd}
{opt ivsensemakr} implements sensitivity analysis tools for instrumental variable (IV) estimates,
as described in Cinelli and Hazlett (2025).

{pstd}
The command performs IV estimation using the Anderson-Rubin (AR) approach, which provides
confidence intervals with correct coverage regardless of instrument strength.
It then computes sensitivity statistics including robustness values (RV) and extreme
robustness values (XRV) that measure how strong omitted variable confounding must be
to change the research conclusions.

{pstd}
When benchmark covariates are specified, the command computes bounds on the strength of
omitted variables based on the observed explanatory power of the benchmarks.
Contour plots visualize how the IV confidence interval changes as a function of the
hypothetical partial R2 of omitted variables with the instrument and with the potential outcome.

{title:Required}

{phang}
{cmd:depvar} {it:{help varname}} specifies the outcome (dependent) variable.

{phang}
{cmd:treatvar} {it:{help varname}} specifies the treatment (endogenous) variable.

{phang}
{cmd:instrvar} {it:{help varname}} specifies the instrumental variable.

{phang}
{cmd:covariates} (optional) {it:{help varlist}} specifies the observed control variables.
These may include factor variables; see {help fvvarlist}.

{title:Options}

{dlgtab:Main}

{phang}
{opt benchmark(varlist)} specifies one or more covariates to be used as benchmarks
for bounding the plausible strength of omitted variables. The benchmark approach
uses the observed explanatory power of these covariates to provide reference points
for how strong the confounding would need to be.

{phang}
{opt kz(numlist)} parameterizes how many times stronger the omitted variables are
related to the instrument compared to the benchmark covariates. Default is 1.

{phang}
{opt ky(numlist)} parameterizes how many times stronger the omitted variables are
related to the potential outcome compared to the benchmark covariates. Default equals kz.

{phang}
{opt kd(numlist)} parameterizes how many times stronger the omitted variables are
related to the treatment compared to the benchmark covariates. Default equals kz.

{phang}
{opt suppress} suppresses the detailed first-stage and reduced-form output tables.

{dlgtab:Sensitivity Parameters}

{phang}
{opt q(real)} specifies the fraction of the effect estimate that would need to be
explained away to be problematic. Default is 1 (i.e., enough confounding to bring
the estimate to zero). Values between 0 and 1 represent partial reductions.

{phang}
{opt alpha(real)} specifies the significance level for hypothesis tests and
confidence intervals. Default is 0.05.

{phang}
{opt h0(real)} specifies the null hypothesis value for the IV estimate.
Default is 0.

{phang}
{opt nomin} computes robustness values for a bias of exactly q, rather than
biases of q or larger (the default).

{dlgtab:Graphing}

{phang}
{opt contourplot} generates two side-by-side contour plots showing the adjusted
lower and upper limits of the Anderson-Rubin confidence interval as a function
of the partial R2 of omitted variable(s) with the instrument (x-axis) and with
the potential outcome (y-axis). The red dashed line shows the critical threshold.

{phang}
{opt tcontourplot} generates a contour plot showing the adjusted t-value for
the Anderson-Rubin test as a function of the partial R2 of omitted variables.

{phang}
{opt clim(numlist)} sets the lower and upper limits for the contour plot axes.
Both values must be between 0 and 1. Default limits are determined automatically.

{phang}
{opt clines(real)} sets the number of contour lines to draw. Default is 7.

{title:Methodology}

{pstd}
The Anderson-Rubin (AR) approach constructs IV confidence intervals by inverting the
test of H0: tau = tau0. The AR point estimate equals the two-stage least squares (2SLS)
estimate (ratio of reduced-form to first-stage coefficients), but the confidence
interval is constructed via test inversion rather than the delta method. This provides
correct coverage even with weak instruments, at the cost of potentially unbounded intervals.

{pstd}
The sensitivity analysis extends the omitted variable bias (OVB) framework of
Cinelli and Hazlett (2020) to the IV setting. The key sensitivity parameters are:

{pstd}
{it:R2zw.x}: partial R2 of the omitted variable(s) W with the instrument Z, given X.

{pstd}
{it:R2y0w.zx}: partial R2 of the omitted variable(s) W with the potential outcome Y(0),
given Z and X.

{pstd}
The {it:Robustness Value (RV)} measures the minimum strength of confounding (in terms
of partial R2) that is sufficient to change the conclusions. The {it:Extreme Robustness
Value (XRV)} considers the worst-case scenario.

{title:Examples}

{pstd}Load example data{p_end}
{phang2}{cmd:. use card.dta, clear}{p_end}

{pstd}Basic IV sensitivity analysis{p_end}
{phang2}{cmd:. ivsensemakr lwage educ nearc4 exper expersq black south smsa reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66}{p_end}

{pstd}With benchmark covariates{p_end}
{phang2}{cmd:. ivsensemakr lwage educ nearc4 exper expersq black south smsa reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, benchmark(black smsa)}{p_end}

{pstd}With contour plots{p_end}
{phang2}{cmd:. ivsensemakr lwage educ nearc4 exper expersq black south smsa reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, benchmark(black smsa) contourplot clim(0 0.09)}{p_end}

{pstd}With custom sensitivity parameters{p_end}
{phang2}{cmd:. ivsensemakr lwage educ nearc4 exper expersq black south smsa reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, benchmark(black smsa) kz(1 2) ky(1 2) q(0.5) alpha(0.1)}{p_end}

{title:Saved results}

{p 4 8 2}
{cmd:ivsensemakr} ereturns the following results, which can be displayed by typing {cmd:ereturn list}.

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(iv_estimate)}}the IV point estimate{p_end}
{synopt:{cmd:e(fs_coef)}}first-stage coefficient of the instrument{p_end}
{synopt:{cmd:e(fs_se)}}first-stage standard error{p_end}
{synopt:{cmd:e(rf_coef)}}reduced-form coefficient of the instrument{p_end}
{synopt:{cmd:e(rf_se)}}reduced-form standard error{p_end}
{synopt:{cmd:e(rho)}}correlation of first-stage and reduced-form residuals{p_end}
{synopt:{cmd:e(ar_t)}}Anderson-Rubin t-statistic{p_end}
{synopt:{cmd:e(dof)}}degrees of freedom{p_end}
{synopt:{cmd:e(rv_iv)}}robustness value for IV estimate{p_end}
{synopt:{cmd:e(xrv_iv)}}extreme robustness value for IV estimate{p_end}
{synopt:{cmd:e(rv_fs)}}robustness value for first-stage{p_end}
{synopt:{cmd:e(xrv_fs)}}extreme robustness value for first-stage{p_end}
{synopt:{cmd:e(rv_rf)}}robustness value for reduced-form{p_end}
{synopt:{cmd:e(xrv_rf)}}extreme robustness value for reduced-form{p_end}
{synopt:{cmd:e(q)}}the value of q{p_end}
{synopt:{cmd:e(alpha)}}the significance level{p_end}
{synopt:{cmd:e(h0)}}the null hypothesis value{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(outcome)}}outcome variable name{p_end}
{synopt:{cmd:e(treatment)}}treatment variable name{p_end}
{synopt:{cmd:e(instrument)}}instrument variable name{p_end}
{synopt:{cmd:e(covariates)}}list of covariate names{p_end}
{synopt:{cmd:e(bench)}}list of benchmark covariates{p_end}
{synopt:{cmd:e(ci_type)}}type of AR confidence interval: "bounded", "disjoint", "all_reals", or "empty"{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(bounds_iv)}}bounds table for IV estimate (kz, ky, r2zw_x, r2y0w_zx, lower_CI, upper_CI){p_end}
{synopt:{cmd:e(bounds_fs)}}bounds table for first-stage{p_end}
{synopt:{cmd:e(bounds_rf)}}bounds table for reduced-form{p_end}

{title:References}

{p 4 8 2}
Cinelli, C. and Hazlett, C. (2025), "An Omitted Variable Bias Framework for Sensitivity
Analysis of Instrumental Variables." {it:Biometrika}.
{browse "https://doi.org/10.1093/biomet/asaf004":doi:10.1093/biomet/asaf004}

{p 4 8 2}
Cinelli, C. and Hazlett, C. (2020), "Making Sense of Sensitivity: Extending Omitted
Variable Bias." {it:Journal of the Royal Statistical Society, Series B} 82(1): 39-67.

{p 4 8 2}
Anderson, T.W. and Rubin, H. (1949), "Estimation of the Parameters of a Single Equation
in a Complete System of Stochastic Equations." {it:Annals of Mathematical Statistics} 20: 46-63.

{title:Authors}

      Carlos Cinelli
      University of Washington

      Chad Hazlett
      UCLA

{title:Also see}

{p 4 8 2}
R package: {browse "https://github.com/carloscinelli/iv.sensemakr":iv.sensemakr}

{p 4 8 2}
OLS version: {help sensemakr:sensemakr} (if installed)
