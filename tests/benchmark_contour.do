* Quick benchmark: Mata marching squares + twoway line at small grids
* Skip Approach 1 — already measured at 380 seconds for 100x100

clear all
set more off
adopath ++ "./"
use card.dta, clear

quietly ivsensemakr lwage educ nearc4 exper expersq black south smsa ///
    reg661 reg662 reg663 reg664 reg665 reg666 reg667 reg668 smsa66, ///
    benchmark(black smsa)

local fs_coef = e(fs_coef)
local fs_se = e(fs_se)
local rf_coef = e(rf_coef)
local rf_se = e(rf_se)
local rho = e(rho)
local ar_dof = e(dof)
local alpha = e(alpha)

display _n "=============================================="
display "BENCHMARK: Marching Squares Scaling"
display "Approach 1 (contourline+addplot 100x100): ~380 sec (measured previously)"
display "=============================================="

mata:

real matrix __bench_ms(real colvector gx, real colvector gy, real matrix Z, real scalar lev) {
    real scalar nx, ny, i, j, ns, ci
    real scalar zbl, zbr, ztl, ztr, x0, x1, y0, y1, f, zmid
    real scalar hb, ht, hl, hr, pb_x, pb_y, pt_x, pt_y, pl_x, pl_y, pr_x, pr_y
    real scalar pc, s1x, s1y, s2x, s2y
    real matrix segs

    nx = length(gx)
    ny = length(gy)
    segs = J((nx-1)*(ny-1)*2, 4, .)
    ns = 0

    for (i = 1; i <= nx-1; i++) {
        for (j = 1; j <= ny-1; j++) {
            zbl = Z[j,i]; zbr = Z[j,i+1]; ztl = Z[j+1,i]; ztr = Z[j+1,i+1]
            if (missing(zbl) | missing(zbr) | missing(ztl) | missing(ztr)) continue

            x0 = gx[i]; x1 = gx[i+1]; y0 = gy[j]; y1 = gy[j+1]
            ci = 8*(ztl>=lev) + 4*(ztr>=lev) + 2*(zbr>=lev) + (zbl>=lev)
            if (ci==0 | ci==15) continue

            hb = ((zbl>=lev) != (zbr>=lev))
            ht = ((ztl>=lev) != (ztr>=lev))
            hl = ((zbl>=lev) != (ztl>=lev))
            hr = ((zbr>=lev) != (ztr>=lev))

            pb_x=.; pb_y=.; pt_x=.; pt_y=.; pl_x=.; pl_y=.; pr_x=.; pr_y=.

            if (hb) { f=(lev-zbl)/(zbr-zbl); pb_x=x0+f*(x1-x0); pb_y=y0; }
            if (ht) { f=(lev-ztl)/(ztr-ztl); pt_x=x0+f*(x1-x0); pt_y=y1; }
            if (hl) { f=(lev-zbl)/(ztl-zbl); pl_x=x0; pl_y=y0+f*(y1-y0); }
            if (hr) { f=(lev-zbr)/(ztr-zbr); pr_x=x1; pr_y=y0+f*(y1-y0); }

            if (ci==5 | ci==10) {
                zmid = (zbl+zbr+ztl+ztr)/4
                if (ci==5) {
                    if (zmid>=lev) {
                        ns++; segs[ns,.] = (pb_x,pb_y,pr_x,pr_y)
                        ns++; segs[ns,.] = (pl_x,pl_y,pt_x,pt_y)
                    } else {
                        ns++; segs[ns,.] = (pb_x,pb_y,pl_x,pl_y)
                        ns++; segs[ns,.] = (pr_x,pr_y,pt_x,pt_y)
                    }
                } else {
                    if (zmid>=lev) {
                        ns++; segs[ns,.] = (pb_x,pb_y,pl_x,pl_y)
                        ns++; segs[ns,.] = (pr_x,pr_y,pt_x,pt_y)
                    } else {
                        ns++; segs[ns,.] = (pb_x,pb_y,pr_x,pr_y)
                        ns++; segs[ns,.] = (pl_x,pl_y,pt_x,pt_y)
                    }
                }
            } else {
                pc=0; s1x=.; s1y=.; s2x=.; s2y=.
                if (hb) { pc++; if (pc==1) { s1x=pb_x; s1y=pb_y; } else { s2x=pb_x; s2y=pb_y; } }
                if (hr) { pc++; if (pc==1) { s1x=pr_x; s1y=pr_y; } else { s2x=pr_x; s2y=pr_y; } }
                if (ht) { pc++; if (pc==1) { s1x=pt_x; s1y=pt_y; } else { s2x=pt_x; s2y=pt_y; } }
                if (hl) { pc++; if (pc==1) { s1x=pl_x; s1y=pl_y; } else { s2x=pl_x; s2y=pl_y; } }
                if (pc==2) { ns++; segs[ns,.] = (s1x,s1y,s2x,s2y); }
            }
        }
    }
    if (ns==0) return(J(0,4,.))
    return(segs[1::ns,.])
}


real matrix __bench_chain(real matrix segs) {
    real scalar n, i, j, found, nr, nc, clen
    real scalar sx, sy, ex, ey
    real colvector used
    real matrix res, ch
    real scalar tol

    n = rows(segs)
    if (n==0) return(J(0,2,.))

    tol = 1e-12
    used = J(n,1,0)
    res = J(3*n, 2, .)
    nr = 0; nc = 0

    for (i=1; i<=n; i++) {
        if (used[i]) continue
        used[i] = 1
        nc++
        if (nc > 1) {
            nr++
            res[nr,1] = .; res[nr,2] = .
        }

        ch = J(n+1, 2, .)
        ch[1,1] = segs[i,1]; ch[1,2] = segs[i,2]
        ch[2,1] = segs[i,3]; ch[2,2] = segs[i,4]
        clen = 2

        ex = ch[clen,1]; ey = ch[clen,2]
        found = 1
        while (found) {
            found = 0
            for (j=1; j<=n; j++) {
                if (used[j]) continue
                if (abs(segs[j,1]-ex)<tol & abs(segs[j,2]-ey)<tol) {
                    used[j]=1; clen++
                    ch[clen,1]=segs[j,3]; ch[clen,2]=segs[j,4]
                    ex=segs[j,3]; ey=segs[j,4]; found=1; break
                }
                if (abs(segs[j,3]-ex)<tol & abs(segs[j,4]-ey)<tol) {
                    used[j]=1; clen++
                    ch[clen,1]=segs[j,1]; ch[clen,2]=segs[j,2]
                    ex=segs[j,1]; ey=segs[j,2]; found=1; break
                }
            }
        }

        sx = ch[1,1]; sy = ch[1,2]
        found = 1
        while (found) {
            found = 0
            for (j=1; j<=n; j++) {
                if (used[j]) continue
                if (abs(segs[j,3]-sx)<tol & abs(segs[j,4]-sy)<tol) {
                    used[j]=1
                    ch = (segs[j,1],segs[j,2]) \ ch[1::clen,.]
                    clen++; sx=segs[j,1]; sy=segs[j,2]; found=1; break
                }
                if (abs(segs[j,1]-sx)<tol & abs(segs[j,2]-sy)<tol) {
                    used[j]=1
                    ch = (segs[j,3],segs[j,4]) \ ch[1::clen,.]
                    clen++; sx=segs[j,3]; sy=segs[j,4]; found=1; break
                }
            }
        }

        for (j=1; j<=clen; j++) { nr++; res[nr,1]=ch[j,1]; res[nr,2]=ch[j,2]; }
    }
    if (nr==0) return(J(0,2,.))
    return(res[1::nr,.])
}


void __bench_run(real scalar fs_coef, real scalar fs_se,
                 real scalar rf_coef, real scalar rf_se,
                 real scalar rho, real scalar dof,
                 real scalar lim_lb, real scalar lim_ub,
                 real scalar alpha, real scalar ngrid) {
    real scalar nx, ny, nt, i, j, z, tc, fc, k
    real colvector gx, gy, r2x, r2y, sef, bf, td, av, bv, cv, dv, result
    real colvector r2ya, cm
    real matrix Z, segs, chain
    real colvector cuts
    real scalar nl, tp, offset

    nx = ngrid; ny = ngrid; nt = nx*ny
    gx = rangen(lim_lb, lim_ub, nx)
    gy = rangen(lim_lb, lim_ub, ny)

    r2x = J(nt,1,.); r2y = J(nt,1,.)
    z = 1
    for (i=1; i<=nx; i++) {
        for (j=ny; j>=1; j--) {
            r2x[z]=gx[i]; r2y[z]=gy[j]; z++
        }
    }

    tc = sqrt(invFtail(1, dof-1, alpha))
    fc = tc/sqrt(dof)
    cm = r2x :< (fc^2 :* (r2y :/ (1:-r2y)))
    r2ya = cm :* (r2x :/ (fc^2 :+ r2x)) + (1:-cm) :* r2y
    sef = sqrt(1:-r2ya) :/ sqrt(1:-r2x)
    bf = sqrt(r2ya :* r2x :/ (1:-r2x))
    td = tc :* sef :* sqrt(dof/(dof-1)) :+ bf :* sqrt(dof)
    av = fs_coef^2 :- fs_se^2 :* td:^2
    bv = 2*rho*rf_se*fs_se :* td:^2 :- 2*rf_coef*fs_coef
    cv = rf_coef^2 :- rf_se^2 :* td:^2
    dv = bv:^2 :- 4 :* av :* cv
    result = J(nt,1,.)
    for (i=1; i<=nt; i++) {
        if (av[i]>0 & dv[i]>=0) result[i] = (-bv[i]-sqrt(dv[i]))/(2*av[i])
    }

    Z = J(ny, nx, .)
    z = 1
    for (i=1; i<=nx; i++) {
        for (j=ny; j>=1; j--) { Z[j,i]=result[z]; z++; }
    }

    // Contour levels
    real colvector fz, cz
    real scalar nf, clo, chi, nc2
    fz = select(result, result:<.)
    nf = length(fz)
    fz = sort(fz, 1)
    clo = fz[max((1,ceil(0.01*nf)))]
    chi = fz[max((1,ceil(0.99*nf)))]
    cz = select(fz, fz:>=clo :& fz:<=chi)
    nc2 = length(cz)

    real colvector qp, rl
    qp = rangen(0,1,7)
    rl = J(7,1,.)
    for (i=1; i<=7; i++) rl[i] = cz[max((1,ceil(qp[i]*nc2)))]
    cuts = round(sort(uniqrows(rl), 1), .001)
    nl = length(cuts)

    printf("  Grid: %g x %g, Levels: %g\n", nx, ny, nl)

    // Extract + chain contour lines
    tp = 0
    pointer(real matrix) colvector lc
    real colvector lnp
    lc = J(nl,1,NULL)
    lnp = J(nl,1,0)

    for (k=1; k<=nl; k++) {
        segs = __bench_ms(gx, gy, Z, cuts[k])
        if (rows(segs)>0) {
            chain = __bench_chain(segs)
            lc[k] = &chain
            lnp[k] = rows(chain)
            tp = tp + rows(chain)
        } else {
            chain = J(0,2,.)
            lc[k] = &chain
        }
    }
    tp = tp + nl - 1
    printf("  Total contour points: %g\n", tp)

    // Store as Stata variables
    real scalar ne
    ne = st_nobs()
    if (tp > ne) (void) st_addobs(tp - ne)
    (void) st_addvar("double", "__cl_x")
    (void) st_addvar("double", "__cl_y")
    (void) st_addvar("double", "__cl_lev")

    offset = 0
    for (k=1; k<=nl; k++) {
        chain = *lc[k]
        if (rows(chain)==0) continue
        if (offset > 0) {
            offset++
            st_store(offset, "__cl_x", .)
            st_store(offset, "__cl_y", .)
            st_store(offset, "__cl_lev", .)
        }
        for (i=1; i<=rows(chain); i++) {
            offset++
            st_store(offset, "__cl_x", chain[i,1])
            st_store(offset, "__cl_y", chain[i,2])
            if (missing(chain[i,1])) {
                st_store(offset, "__cl_lev", .)
            } else {
                st_store(offset, "__cl_lev", k)
            }
        }
    }

    st_numscalar("__bench_nlevels", nl)
    st_matrix("__bench_cuts", cuts)
}

end


* =========================================================================
* Test at 20, 50, 100, 200
* =========================================================================
foreach gridsize in 20 50 {
    display _n "--- marching squares + twoway line (`gridsize'x`gridsize') ---"

    preserve
    timer clear
    timer on 1

    mata: __bench_run(`fs_coef', `fs_se', `rf_coef', `rf_se', `rho', `ar_dof', 0, 0.09, `alpha', `gridsize')

    local nlevels = scalar(__bench_nlevels)
    local plot_cmd ""
    forvalues k = 1(1)`nlevels' {
        local plot_cmd `"`plot_cmd' (line __cl_y __cl_x if __cl_lev == `k', lcolor(black) lwidth(vthin))"'
    }

    twoway `plot_cmd', ///
        name(__bench, replace) ///
        legend(off) title("MS `gridsize'x`gridsize'") ///
        xtitle("R2zw.x") ytitle("R2y0w.zx")

    graph export "misc/bench_ms_`gridsize'.png", replace width(800)

    timer off 1
    timer list 1
    local time_ms = r(t1)

    capture graph drop __bench
    restore

    display "  => Time (`gridsize'x`gridsize'): `time_ms' seconds"
}

display _n "=============================================="
display "BENCHMARK COMPLETE"
display "=============================================="
display "Approach 1 (contourline+addplot 100x100): ~380 sec (measured)"
display "See above for Approach 2 scaling at 20, 50"
