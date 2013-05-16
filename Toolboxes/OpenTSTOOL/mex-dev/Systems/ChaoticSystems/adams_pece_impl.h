#ifndef ADAM_PECE_DECL_H
#define ADAM_PECE_DECL_H

#include <cmath>

#define TRUE_ (1)
#define FALSE_ (0)
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)

template<class ExplicitODE>
int Adams<ExplicitODE>::ode(double *t, double *tout, double *relerr, double *abserr, long *iflag)
{
    /* Initialized data */

    static long ialpha = 1;
    static long ih = 89;
    static long ihold = 90;
    static long istart = 91;
    static long itold = 92;
    static long idelsn = 93;
    static long ibeta = 13;
    static long isig = 25;
    static long iv = 38;
    static long iw = 50;
    static long ig = 62;
    static long iphase = 75;
    static long ipsi = 76;
    static long ix = 88;

    long iphi;
    long nornd, start, phase1;

	long ip, iypout, iyp, iwt, iyy;
	
	long* const iwork = iwork_c - 1;
	double* const y = f.v-1;
    double* const work = work_c - 1;
   
    /* Function Body */
    iyy = 100;
    iwt = iyy + neqn;
    ip = iwt + neqn;
    iyp = ip + neqn;
    iypout = iyp + neqn;
    iphi = iypout + neqn;
    if (abs(*iflag) == 1) {
		goto L1;
    }
    start = work[istart] > 0.;
    phase1 = work[iphase] > 0.;
    nornd = iwork[2] != -1;
L1:
    de(&y[1], t, tout, relerr, abserr, iflag, &work[iyy], &
	    work[iwt], &work[ip], &work[iyp], &work[iypout], &work[iphi], &
	    work[ialpha], &work[ibeta], &work[isig], &work[iv], &work[iw], &
	    work[ig], &phase1, &work[ipsi], &work[ix], &work[ih], &work[ihold]
	    , &start, &work[itold], &work[idelsn], &iwork[1], &nornd, &iwork[3], &iwork[4], &iwork[5]);
   
    work[istart] = -1.;
    if (start) {
		work[istart] = 1.;
    }
    work[iphase] = -1.;
    if (phase1) {
		work[iphase] = 1.;
    }
    iwork[2] = -1;
    if (nornd) {
		iwork[2] = 1;
    }
	
    return 0;
} /* ode_ */


template<class ExplicitODE>
int Adams<ExplicitODE>::de(double *y, double *t, 
	double *tout, double *relerr, double *abserr, long *
	iflag, double *yy, double *wt, double *p, double *yp, 
	double *ypout, double *phi, double *alpha, double *
	beta, double *sig, double *v, double *w, double *g, 
	long *phase1, double *psi, double *x, double *h__, 
	double *hold, long *start, double *told, double *
	delsgn, long *ns, long *nornd, long *k, long *kold, 
	long *isnold)
{
    /* Initialized data */

    const long maxnum = 1000;

    /* System generated locals */
    long phi_dim1, phi_offset, i__1;
    double d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    double tend;
    long l;
    long crash, stiff;

    double absdel, abseps, releps;
    long nostep;
    double del, eps;
    long isn, kle4;

    /* Parameter adjustments */
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1 * 1;
    phi -= phi_offset;
    --ypout;
    --yp;
    --p;
    --wt;
    --yy;
    --y;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;
    --psi;

    /* Function Body */

/*            ***            ***            *** */
/*   test for improper parameters */

    if (neqn < 1) {
		goto L10;
    }
    if (*t == *tout) {
		goto L10;
    }
    if (*relerr < 0. || *abserr < 0.) {
		goto L10;
    }
    eps = max(*relerr,*abserr);
    if (eps <= 0.) {
		goto L10;
    }
    if (*iflag == 0) {
		goto L10;
    }
    isn = i_sign(&c__1, iflag);
    *iflag = abs(*iflag);
    if (*iflag == 1) {
		goto L20;
    }
    if (*t != *told) {
		goto L10;
    }
    if (*iflag >= 2 && *iflag <= 5) {
		goto L20;
    }
L10:
    *iflag = 6;
    return 0;

/*   on each call set interval of integration and counter for number of */
/*   steps.  adjust input error tolerances to define weight vector for */
/*   subroutine  step */

L20:
    del = *tout - *t;
    absdel = abs(del);
    tend = *t + del * 10.;
    if (isn < 0) {
		tend = *tout;
    }
    nostep = 0;
    kle4 = 0;
    stiff = FALSE_;
    releps = *relerr / eps;
    abseps = *abserr / eps;
    if (*iflag == 1) {
		goto L30;
    }
    if (*isnold < 0) {
		goto L30;
    }
    if (*delsgn * del > 0.) {
		goto L50;
    }

/*   on start and restart also set work variables x and yy(*), store the */
/*   direction of integration and initialize the step size */

L30:
    *start = TRUE_;
    *x = *t;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
		yy[l] = y[l];
    }
    *delsgn = d_sign(&c_b11, &del);
	
/* Computing MAX */
    
	d__3 = (d__1 = *tout - *x, abs(d__1)), d__4 = fouru * abs(*x);
    d__2 = max(d__3,d__4);
    d__5 = *tout - *x;
    *h__ = d_sign(&d__2, &d__5);

/*   if already past output point, interpolate and return */

L50:
    if ((d__1 = *x - *t, abs(d__1)) < absdel) {
		goto L60;
    }
    intrp(*x, &yy[1], *tout, &y[1], &ypout[1], *kold, &phi[phi_offset], &psi[1]);
	
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   if cannot go past output point and sufficiently close, */
/*   extrapolate and return */

L60:
    if (isn > 0 || (d__1 = *tout - *x, abs(d__1)) >= fouru * abs(*x)) {
		goto L80;
    }
    *h__ = *tout - *x;
    f(x, &yy[1], &yp[1]);
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
		y[l] = yy[l] + *h__ * yp[l];
    }
    *iflag = 2;
    *t = *tout;
    *told = *t;
    *isnold = isn;
    return 0;

/*   test for too many steps */

L80:
    if (nostep < maxnum) {
		goto L100;
    }
    *iflag = isn << 2;
    if (stiff) {
		*iflag = isn * 5;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L90: */
		y[l] = yy[l];
    }
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   limit step size, set weight vector and take a step */

L100:
	/* Computing MIN */

    d__3 = abs(*h__), d__4 = (d__1 = tend - *x, abs(d__1));
    d__2 = min(d__3,d__4);
    *h__ = d_sign(&d__2, h__);
	
    i__1 = neqn;
	
    for (l = 1; l <= i__1; ++l) {
		wt[l] = releps * (d__1 = yy[l], abs(d__1)) + abseps;
    }
	
    step(x, &yy[1], h__, &eps, &wt[1], start, hold, k, kold, &
	    crash, &phi[phi_offset], &p[1], &yp[1], &psi[1], &alpha[1], &beta[
	    1], &sig[1], &v[1], &w[1], &g[1], phase1, ns, nornd);

/*   test for tolerances too small */

    if (! crash) {
		goto L130;
    }
    *iflag = isn * 3;
    *relerr = eps * releps;
    *abserr = eps * abseps;
    i__1 = neqn;
	
    for (l = 1; l <= i__1; ++l) {
		y[l] = yy[l];
    }
	
    *t = *x;
    *told = *t;
    *isnold = 1;
    return 0;

/*   augment counter on number of steps and test for stiffness */

L130:
    ++nostep;
    ++kle4;
    if (*kold > 4) {
	kle4 = 0;
    }
    if (kle4 >= 50) {
		stiff = TRUE_;
    }
    goto L50;
} /* de_ */


template<class ExplicitODE>
int Adams<ExplicitODE>::step(double *x, double *y, double *h__, double *eps, double *wt, long *start, double *hold, long *k, 
	long *kold, long *crash, double *phi, double *p, double *yp, double *psi, double *alpha, double *beta, double *sig, double *v, 
	double *w, double *g, long *phase1, long *ns, long *nornd)
{
    /* Initialized data */

    static double two[13] = { 2.,4.,8.,16.,32.,64.,128.,256.,512.,1024.,
	    2048.,4096.,8192. };
    static double gstr[13] = { .5,.0833,.0417,.0264,.0188,.0143,.0114,
	    .00936,.00789,.00679,.00592,.00524,.00468 };

    /* System generated locals */
    long phi_dim1, phi_offset, i__1, i__2;
    double d__1, d__2, d__3;

    /* Local variables */
    double absh, hnew;
    long knew;
    double xold, erkm1, erkm2, erkp1, temp1, temp2, temp3, temp4, temp5, temp6, p5eps;
    long i__, j, l;
    double r__;
    long ifail;
    double reali, round;

    long limit1, limit2, iq;
    double realns;
    long im1, km1, km2, ip1, kp1, kp2;
    double erk, err, tau, rho, sum;
    long nsm2, nsp1, nsp2;
    /* Parameter adjustments */
    --yp;
    --p;
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1 * 1;
    phi -= phi_offset;
    --wt;
    --y;
    --psi;
    --alpha;
    --beta;
    --sig;
    --v;
    --w;
    --g;

    /* Function Body */
/*     data g(1),g(2)/1.0,0.5/,sig(1)/1.0/ */

/*       ***     begin block 0     *** */
/*   check if step size or error tolerance is too small for machine */
/*   precision.  if first step, initialize phi array and estimate a */
/*   starting step size. */
/*                   *** */

/*   if step size is too small, determine an acceptable one */

    *crash = TRUE_;
    if (abs(*h__) >= fouru * abs(*x)) {
	goto L5;
    }
    d__1 = fouru * abs(*x);
    *h__ = d_sign(&d__1, h__);
    return 0;
L5:
    p5eps = *eps * .5;

/*   if error tolerance is too small, increase it to an acceptable value */

    round = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L10: */
/* Computing 2nd power */
	d__1 = y[l] / wt[l];
	round += d__1 * d__1;
    }
    round = twou * sqrt(round);
    if (p5eps >= round) {
	goto L15;
    }
    *eps = round * 2.f * (fouru + 1.);
    return 0;
L15:
    *crash = FALSE_;
    g[1] = 1.;
    g[2] = .5;
    sig[1] = 1.;
    if (! (*start)) {
	goto L99;
    }

/*   initialize.  compute appropriate step size for first step */

    f(x, &y[1], &yp[1]);
    sum = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + phi_dim1] = yp[l];
	phi[l + (phi_dim1 << 1)] = 0.;
/* L20: */
/* Computing 2nd power */
	d__1 = yp[l] / wt[l];
	sum += d__1 * d__1;
    }
    sum = sqrt(sum);
    absh = abs(*h__);
    if (*eps < sum * 16. * *h__ * *h__) {
	absh = sqrt(*eps / sum) * .25;
    }
/* Computing MAX */
    d__2 = absh, d__3 = fouru * abs(*x);
    d__1 = max(d__2,d__3);
    *h__ = d_sign(&d__1, h__);
    *hold = 0.;
    *k = 1;
    *kold = 0;
    *start = FALSE_;
    *phase1 = TRUE_;
    *nornd = TRUE_;
    if (p5eps > round * 100.) {
	goto L99;
    }
    *nornd = FALSE_;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L25: */
	phi[l + phi_dim1 * 15] = 0.;
    }
L99:
    ifail = 0;
/*       ***     end block 0     *** */

/*       ***     begin block 1     *** */
/*   compute coefficients of formulas for this step.  avoid computing */
/*   those quantities not changed when step size is not changed. */
/*                   *** */

L100:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    km2 = *k - 2;

/*   ns is the number of steps taken with size h, including the current */
/*   one.  when k.lt.ns, no coefficients change */

    if (*h__ != *hold) {
	*ns = 0;
    }
    if (*ns <= *kold) {
	++(*ns);
    }
    nsp1 = *ns + 1;
    if (*k < *ns) {
	goto L199;
    }

/*   compute those components of alpha(*),beta(*),psi(*),sig(*) which */
/*   are changed */

    beta[*ns] = 1.;
    realns = (double) (*ns);
    alpha[*ns] = 1. / realns;
    temp1 = *h__ * realns;
    sig[nsp1] = 1.;
    if (*k < nsp1) {
	goto L110;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	temp2 = psi[im1];
	psi[im1] = temp1;
	beta[i__] = beta[im1] * psi[im1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	reali = (double) i__;
/* L105: */
	sig[i__ + 1] = reali * alpha[i__] * sig[i__];
    }
L110:
    psi[*k] = temp1;

/*   compute coefficients g(*) */

/*   initialize v(*) and set w(*).  g(2) is set in data statement */

    if (*ns > 1) {
	goto L120;
    }
    i__1 = *k;
    for (iq = 1; iq <= i__1; ++iq) {
	temp3 = (double) (iq * (iq + 1));
	v[iq] = 1. / temp3;
/* L115: */
	w[iq] = v[iq];
    }
    goto L140;

/*   if order was raised, update diagonal part of v(*) */

L120:
    if (*k <= *kold) {
	goto L130;
    }
    temp4 = (double) (*k * kp1);
    v[*k] = 1. / temp4;
    nsm2 = *ns - 2;
    if (nsm2 < 1) {
	goto L130;
    }
    i__1 = nsm2;
    for (j = 1; j <= i__1; ++j) {
	i__ = *k - j;
/* L125: */
	v[i__] -= alpha[j + 1] * v[i__ + 1];
    }

/*   update v(*) and set w(*) */

L130:
    limit1 = kp1 - *ns;
    temp5 = alpha[*ns];
    i__1 = limit1;
    for (iq = 1; iq <= i__1; ++iq) {
	v[iq] -= temp5 * v[iq + 1];
/* L135: */
	w[iq] = v[iq];
    }
    g[nsp1] = w[1];

/*   compute the g(*) in the work vector w(*) */

L140:
    nsp2 = *ns + 2;
    if (kp1 < nsp2) {
	goto L199;
    }
    i__1 = kp1;
    for (i__ = nsp2; i__ <= i__1; ++i__) {
	limit2 = kp2 - i__;
	temp6 = alpha[i__ - 1];
	i__2 = limit2;
	for (iq = 1; iq <= i__2; ++iq) {
/* L145: */
	    w[iq] -= temp6 * w[iq + 1];
	}
/* L150: */
	g[i__] = w[1];
    }
L199:
/*       ***     end block 1     *** */

/*       ***     begin block 2     *** */
/*   predict a solution p(*), evaluate derivatives using predicted */
/*   solution, estimate local error at order k and errors at orders k, */
/*   k-1, k-2 as if constant step size were used. */
/*                   *** */

/*   change phi to phi star */

    if (*k < nsp1) {
	goto L215;
    }
    i__1 = *k;
    for (i__ = nsp1; i__ <= i__1; ++i__) {
	temp1 = beta[i__];
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
/* L205: */
	    phi[l + i__ * phi_dim1] = temp1 * phi[l + i__ * phi_dim1];
	}
/* L210: */
    }

/*   predict solution and differences */

L215:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1];
	phi[l + kp1 * phi_dim1] = 0.;
/* L220: */
	p[l] = 0.;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__ = kp1 - j;
	ip1 = i__ + 1;
	temp2 = g[i__];
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
	    p[l] += temp2 * phi[l + i__ * phi_dim1];
/* L225: */
	    phi[l + i__ * phi_dim1] += phi[l + ip1 * phi_dim1];
	}
/* L230: */
    }
    if (*nornd) {
	goto L240;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	tau = *h__ * p[l] - phi[l + phi_dim1 * 15];
	p[l] = y[l] + tau;
/* L235: */
	phi[l + (phi_dim1 << 4)] = p[l] - y[l] - tau;
    }
    goto L250;
L240:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L245: */
	p[l] = y[l] + *h__ * p[l];
    }
L250:
    xold = *x;
    *x += *h__;
    absh = abs(*h__);
    f(x, &p[1], &yp[1]);

/*   estimate errors at orders k,k-1,k-2 */

    erkm2 = 0.;
    erkm1 = 0.;
    erk = 0.;
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	temp3 = 1. / wt[l];
	temp4 = yp[l] - phi[l + phi_dim1];
	if (km2 < 0) {
	    goto L265;
	} else if (km2 == 0) {
	    goto L260;
	} else {
	    goto L255;
	}
L255:
/* Computing 2nd power */
	d__1 = (phi[l + km1 * phi_dim1] + temp4) * temp3;
	erkm2 += d__1 * d__1;
L260:
/* Computing 2nd power */
	d__1 = (phi[l + *k * phi_dim1] + temp4) * temp3;
	erkm1 += d__1 * d__1;
L265:
/* Computing 2nd power */
	d__1 = temp4 * temp3;
	erk += d__1 * d__1;
    }
    if (km2 < 0) {
	goto L280;
    } else if (km2 == 0) {
	goto L275;
    } else {
	goto L270;
    }
L270:
    erkm2 = absh * sig[km1] * gstr[km2 - 1] * sqrt(erkm2);
L275:
    erkm1 = absh * sig[*k] * gstr[km1 - 1] * sqrt(erkm1);
L280:
    temp5 = absh * sqrt(erk);
    err = temp5 * (g[*k] - g[kp1]);
    erk = temp5 * sig[kp1] * gstr[*k - 1];
    knew = *k;

/*   test if order should be lowered */

    if (km2 < 0) {
	goto L299;
    } else if (km2 == 0) {
	goto L290;
    } else {
	goto L285;
    }
L285:
    if (max(erkm1,erkm2) <= erk) {
	knew = km1;
    }
    goto L299;
L290:
    if (erkm1 <= erk * .5) {
	knew = km1;
    }

/*   test if step successful */

L299:
    if (err <= *eps) {
	goto L400;
    }
/*       ***     end block 2     *** */

/*       ***     begin block 3     *** */
/*   the step is unsuccessful.  restore  x, phi(*,*), psi(*) . */
/*   if third consecutive failure, set order to one.  if step fails more */
/*   than three times, consider an optimal step size.  double error */
/*   tolerance and return if estimated step size is too small for machine */
/*   precision. */
/*                   *** */

/*   restore x, phi(*,*) and psi(*) */

    *phase1 = FALSE_;
    *x = xold;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 1. / beta[i__];
	ip1 = i__ + 1;
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
/* L305: */
	    phi[l + i__ * phi_dim1] = temp1 * (phi[l + i__ * phi_dim1] - phi[
		    l + ip1 * phi_dim1]);
	}
/* L310: */
    }
    if (*k < 2) {
	goto L320;
    }
    i__1 = *k;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L315: */
	psi[i__ - 1] = psi[i__] - *h__;
    }

/*   on third failure, set order to one.  thereafter, use optimal step */
/*   size */

L320:
    ++ifail;
    temp2 = .5;
    if ((i__1 = ifail - 3) < 0) {
	goto L335;
    } else if (i__1 == 0) {
	goto L330;
    } else {
	goto L325;
    }
L325:
    if (p5eps < erk * .25) {
	temp2 = sqrt(p5eps / erk);
    }
L330:
    knew = 1;
L335:
    *h__ = temp2 * *h__;
    *k = knew;
    if (abs(*h__) >= fouru * abs(*x)) {
	goto L340;
    }
    *crash = TRUE_;
    d__1 = fouru * abs(*x);
    *h__ = d_sign(&d__1, h__);
    *eps += *eps;
    return 0;
L340:
    goto L100;
/*       ***     end block 3     *** */

/*       ***     begin block 4     *** */
/*   the step is successful.  correct the predicted solution, evaluate */
/*   the derivatives using the corrected solution and update the */
/*   differences.  determine best order and step size for next step. */
/*                   *** */
L400:
    *kold = *k;
    *hold = *h__;

/*   correct and evaluate */

    temp1 = *h__ * g[kp1];
    if (*nornd) {
	goto L410;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	rho = temp1 * (yp[l] - phi[l + phi_dim1]) - phi[l + (phi_dim1 << 4)];
	y[l] = p[l] + rho;
/* L405: */
	phi[l + phi_dim1 * 15] = y[l] - p[l] - rho;
    }
    goto L420;
L410:
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L415: */
	y[l] = p[l] + temp1 * (yp[l] - phi[l + phi_dim1]);
    }
L420:
    f(x, &y[1], &yp[1]);

/*   update differences for next step */

    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
	phi[l + kp1 * phi_dim1] = yp[l] - phi[l + phi_dim1];
/* L425: */
	phi[l + kp2 * phi_dim1] = phi[l + kp1 * phi_dim1] - phi[l + kp2 * 
		phi_dim1];
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = neqn;
	for (l = 1; l <= i__2; ++l) {
/* L430: */
	    phi[l + i__ * phi_dim1] += phi[l + kp1 * phi_dim1];
	}
/* L435: */
    }

/*   estimate error at order k+1 unless: */
/*     in first phase when always raise order, */
/*     already decided to lower order, */
/*     step size not constant so estimate unreliable */

    erkp1 = 0.;
    if (knew == km1 || *k == 12) {
	*phase1 = FALSE_;
    }
    if (*phase1) {
	goto L450;
    }
    if (knew == km1) {
	goto L455;
    }
    if (kp1 > *ns) {
	goto L460;
    }
    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
/* L440: */
/* Computing 2nd power */
	d__1 = phi[l + kp2 * phi_dim1] / wt[l];
	erkp1 += d__1 * d__1;
    }
    erkp1 = absh * gstr[kp1 - 1] * sqrt(erkp1);

/*   using estimated error at order k+1, determine appropriate order */
/*   for next step */

    if (*k > 1) {
	goto L445;
    }
    if (erkp1 >= erk * .5) {
	goto L460;
    }
    goto L450;
L445:
    if (erkm1 <= min(erk,erkp1)) {
	goto L455;
    }
    if (erkp1 >= erk || *k == 12) {
	goto L460;
    }

/*   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have */
/*   been lowered in block 2.  thus order is to be raised */

/*   raise order */

L450:
    *k = kp1;
    erk = erkp1;
    goto L460;

/*   lower order */

L455:
    *k = km1;
    erk = erkm1;

/*   with new order determine appropriate step size for next step */

L460:
    hnew = *h__ + *h__;
    if (*phase1) {
	goto L465;
    }
    if (p5eps >= erk * two[*k]) {
	goto L465;
    }
    hnew = *h__;
    if (p5eps >= erk) {
	goto L465;
    }
    temp2 = (double) (*k + 1);
    d__1 = p5eps / erk;
    d__2 = 1. / temp2;
    r__ = pow_dd(&d__1, &d__2);
/* Computing MAX */
    d__1 = .5, d__2 = min(.9,r__);
    hnew = absh * max(d__1,d__2);
/* Computing MAX */
    d__2 = hnew, d__3 = fouru * abs(*x);
    d__1 = max(d__2,d__3);
    hnew = d_sign(&d__1, h__);
L465:
    *h__ = hnew;
    return 0;
/*       ***     end block 4     *** */
} /* step_ */



template<class ExplicitODE>
int Adams<ExplicitODE>::intrp(const double x, const double *y, const double xout, double *yout, double *ypout, const long kold, 
	double *phi, double *psi)
{
    /* Initialized data */

    static double g[13] = { 1. };
    static double rho[13] = { 1. };

    /* System generated locals */
    long phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
    double term, temp1, temp2, temp3;
    long i__, j, l;
    double gamma, w[13];
    long limit1;
    double psijm1, hi;
    long ki, jm1;
    double eta;
    long kip1;


/*   the methods in subroutine  step  approximate the solution near  x */
/*   by a polynomial.  subroutine  intrp  approximates the solution at */
/*   xout  by evaluating the polynomial there.  information defining this */
/*   polynomial is passed from  step  so  intrp  cannot be used alone. */

/*   this code is completely explained and documented in the text, */
/*   computer solution of ordinary differential equations:  the initial */
/*   value problem  by l. f. shampine and m. k. gordon. */

/*   input to intrp -- */

/*   all floating point variables are double precision */
/*   the user provides storage in the calling program for the arrays in */
/*   the call list */
/*   and defines */
/*      xout -- point at which solution is desired. */
/*   the remaining parameters are defined in  step  and passed to  intrp */
/*   from that subroutine */

/*   output from  intrp -- */

/*      yout(*) -- solution at  xout */
/*      ypout(*) -- derivative of solution at  xout */
/*   the remaining parameters are returned unaltered from their input */
/*   values.  integration with  step  may be continued. */

    /* Parameter adjustments */
    phi_dim1 = neqn;
    phi_offset = 1 + phi_dim1 * 1;
    phi -= phi_offset;

    --ypout;
    --yout;
    --y;
    --psi;

    /* Function Body */

    hi = xout - x;
    ki = kold + 1;
    kip1 = ki + 1;

/*   initialize w(*) for computing g(*) */

    i__1 = ki;
    for (i__ = 1; i__ <= i__1; ++i__) {
		temp1 = (double) i__;
		w[i__ - 1] = 1. / temp1;
    }
    term = 0.;

/*   compute g(*) */

    i__1 = ki;
	
    for (j = 2; j <= i__1; ++j) {
		jm1 = j - 1;
		psijm1 = psi[jm1];
		gamma = (hi + term) / psijm1;
		eta = hi / psijm1;
		limit1 = kip1 - j;
		i__2 = limit1;
		for (i__ = 1; i__ <= i__2; ++i__) {
	    	w[i__ - 1] = gamma * w[i__ - 1] - eta * w[i__];
		}
		g[j - 1] = w[0];
		rho[j - 1] = gamma * rho[jm1 - 1];
		term = psijm1;
    }

/*   interpolate */

    i__1 = neqn;
    for (l = 1; l <= i__1; ++l) {
		ypout[l] = 0.;
		yout[l] = 0.;
    }
	
    i__1 = ki;
    
	for (j = 1; j <= i__1; ++j) {
		i__ = kip1 - j;
		temp2 = g[i__ - 1];
		temp3 = rho[i__ - 1];
		i__2 = neqn;
		for (l = 1; l <= i__2; ++l) {
	    	yout[l] += temp2 * phi[l + i__ * phi_dim1];
	    	ypout[l] += temp3 * phi[l + i__ * phi_dim1];
		}
    }

    i__1 = neqn;

    for (l = 1; l <= i__1; ++l) {
		yout[l] = y[l] + hi * yout[l];
    }
    return 0;
} /* intrp_ */

double NumericalAlgorithm::machine_precision()
{
	double halfu = 0.5;

	while(((double)(1.0 + halfu)) > 1.0) {
		halfu *= 0.5;
	}

	return 2 * halfu;
}	

#endif	
