* =============================================================================
* Test contour plots and additional options
* =============================================================================

clear all
set more off

adopath ++ "./"
use card.dta, clear

* =========================================================================
* Test: Contour plot
* =========================================================================
display _n "TEST: Contour plot"

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa) contourplot clim(0 0.09)

graph export "misc/contourplot.png", replace width(1200)

* =========================================================================
* Test: t-Contour plot
* =========================================================================
display _n "TEST: t-Contour plot"

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa) tcontourplot clim(0 0.09)

graph export "misc/tcontourplot.png", replace width(1200)

* =========================================================================
* Test: No benchmarks (just sensitivity stats)
* =========================================================================
display _n "TEST: No benchmarks"

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66

* =========================================================================
* Test: suppress option
* =========================================================================
display _n "TEST: Suppress option"

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa) suppress

* =========================================================================
* Test: nomin option
* =========================================================================
display _n "TEST: nomin option"

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa) nomin

* =========================================================================
* Test: custom q and alpha
* =========================================================================
display _n "TEST: Custom q and alpha"

ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa) q(0.5) alpha(0.10)

display _n "ALL PLOT/OPTION TESTS COMPLETED"
