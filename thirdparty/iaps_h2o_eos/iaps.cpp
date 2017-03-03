/*------------------------------------------------------------*- C -*-
# $Id: iaps.cpp,v 1.1.1.1 2006/08/09 13:48:38 apaluszn Exp $
#---------------------------------------------------------------------
# author(s):  Olaf Bauer
# maintainer: Ole Engel <engel@tu-harburg.de>
# copyright:  GNU LIBRARY GENERAL PUBLIC LICENSE (see COPYING.LIB-2.0)
#---------------------------------------------------------------------
# IAPS formulation 1984
# main part of PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/

/*********************************************************************
	      see file PROST4.TXT for help
**********************************************************************/

#include <cfloat>   /* DBL_EPSILON */
#include <cmath>    /* fabs(), sqrt(), pow(), exp(), log() */
#include <cstdlib>  /* NULL */
#include "iaps.h"

using namespace std;

/*-------------------- external functions -------------------------*/

void sat_t(double t, Prop *pliq, Prop *pvap)
/*
  Computes liquid and vapour properties for a given saturation temperature
*/
{
    double dl, dv, p;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if ((pliq == NULL) || (pvap == NULL)) { 
	if (pliq != NULL) { pliq->error = 1; }
	if (pvap != NULL) { pvap->error = 1; }
	return; 
    }

    if ((t > crit.t) || (t < tripl.t)) {
	pliq->error = 1;
	pvap->error = 1;
    } else {
	psat(t, &p, &dl, &dv, &MLiq, &MPro);
	format_pro(t, dv, &MPro, pvap);      /* MPro -> *pvap */
	load(&MLiq, &MPro);                  /* MLiq ->  MPro */
	format_pro(t, dl, &MPro, pliq);      /* MPro -> *pliq */
    }
}

void sat_p(double p, Prop *pliq, Prop *pvap)
/*
  Computes liquid and vapour properties for a given saturation pressure
*/
{
    double dl, dv, t;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if ((pliq == NULL) || (pvap == NULL)) { 
	if (pliq != NULL) { pliq->error = 1; }
	if (pvap != NULL) { pvap->error = 1; }
	return; 
    }

    p *= 1.0e-6;
    if ((p > crit.p) || (p < tripl.p)) {
	pliq->error = 1;
	pvap->error = 1;
    } else {
	tsat(p, &t, &dl, &dv, &MLiq, &MPro);
	format_pro(t, dv, &MPro, pvap);      /* MPro -> *pvap */
	load(&MLiq, &MPro);                  /* MLiq ->  MPro */
	format_pro(t, dl, &MPro, pliq);      /* MPro -> *pliq */
    }
}

void water_td(double t, double d, Prop *prop)
/*
  Computes properties for given temperature and density
*/
{
    double p, dl, dv, x;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    d *= 1.0e-3;
    if (! valid_td(t, d)) {
	prop->error = 1;
	return;
    }
    if ((t >= tripl.t) && (t <= crit.t)) {
	psat(t, &p, &dl, &dv, &MLiq, &MPro);
	if ((d > dv) && (d < dl)) {
	    x = (1.0 / d - 1.0 / dl) / (1.0 / dv - 1.0 / dl);
	    format_two(t, p, x, dl, dv, &MLiq, &MPro, prop);
	    return;
	}
    }
    td(t, d, &MPro, prop);
}

void water_tp(double t, double p, double d, double dp, Prop *prop)
/*
  Computes properties for given temperature and pressure in the one-
  phase-region within a given tolerance of pressure using a given
  approximation of density.
  No output in case of saturation, i.e. |1 - p / ps(T)| < 1e-6.
*/
{
    double reg, dl, dv;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    d *= 1.0e-3;
    prop->error = 1;
    if (! valid_tp(t, p)) {
	return;
    }
    if (t < crit.t) {
	reg = region_tp(t, p, &dl, &dv, &MLiq, &MPro);
	if (reg == SATURATED) {
	    return;                                 /* saturated */
	}
	if ((reg == LIQUID) && (d < crit.d)) {         /* liquid */
	    d = 1.01 * crit.d;
	} else if ((reg == VAPOUR) && (d > crit.d)) {  /* vapour */
	    d = 0.99 * crit.d;
	}
    }
    tp(t, p, &d, dp, &MPro, prop);
}

void water_ph(double p, double h, double t, double d,
	      double dp, double dh, Prop *prop)
/*
  Computes properties for given pressure and enthalpy within given
  tolerances of pressure and enthalpy using given approximations
  of temperature and density.
  For invalid (p,h) it extrapolates from bounds of validity.
*/
{
    double ts, dl, dv, x;
    int reg;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    h *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_ph(p, h)) {
	extra_ph(p, h, &t, &d, dp, dh, &MPro, prop);  /* extrapolation */
	prop->error = 1;
	return;
    }
    reg = region_ph(p, h, &ts, &dl, &dv, &MLiq, &MPro);
    if (reg == TWO) {
	x = (h - MLiq.spro.h) / (MPro.spro.h - MLiq.spro.h);
	format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
    } else {
	ph(p, h, &t, &d, dp, dh, &MPro, prop);
    }
}

void water_ps(double p, double s, double t, double d,
	      double dp, double ds, Prop *prop)
/*
  Computes properties for given pressure and entropy within given
  tolerances of pressure and entropy using given approximations
  of temperature and density.
*/
{
    double ts,dl,dv,sl,sv,x;
    S_mpro MPro;
    S_mliq MLiq;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    s *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_ps(p, s)) {
	prop->error = 1;
	return;
    }
    if ((p >= tripl.p) && (p <= crit.p)) {
	tsat(p, &ts, &dl, &dv, &MLiq, &MPro);
	sl = MLiq.spro.s;
	sv = MPro.spro.s;
	if ((s > sl) && (s < sv)) {
	    x = (s - sl) / (sv - sl);
	    format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
	    return;
	}
    }
    ps(p, s, &t, &d, dp, ds, &MPro, prop);
}

void water_hd(double h, double d, double t, double dh, Prop *prop)
/*
  Computes properties for given enthalpy h and density d within a given
  relative tolerance of enthalpy dh using a given approximation of
  temperature t.
*/
{
    int reg;
    double ts,p,dl,dv,x;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    h *= 1e-3;
    d *= 1e-3;
    if (! valid_hd(h, d)) {
	prop->error = 1;           /* invalid */
	return;
    }
    reg = region_hd(h, d, dh, &ts, &p, &dl, &dv, &x, &MLiq, &MPro);
    if (reg == TWO) {
	format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);   /* two-phase */
    } else {
	hd(h, d, &t, dh, &MPro, prop);               /* one-phase */
    }
}



/*-------------------- Isolated external functions --------------------*/

double speed(Prop *prop)
/*
  Speed of sound for a calculated prop structure,
  only in one-phase-region.
*/
{
    /* check pointer */
    if (prop == NULL) { return 0.0; }

    return sqrt(fabs(prop->cp * prop->dp->d_CT / prop->cv));
}

double wkappa(Prop *prop)
/*
  Isentropic expansion coefficient for a calculated prop structure,
  only in one-phase-region.
*/
{
    const double T   = prop->T;
    const double d   = prop->d;
    const double p   = prop->p;
    const double cv  = prop->cv;
    const double dpt = prop->dp->T_Cd;
    const double dpd = prop->dp->d_CT;
    return (d * d * cv * dpd + T * dpt * dpt) / (d * p * cv);
}

double wbetas(Prop *prop)
/*
  Isentropic pressure-temperature coefficient for a calculated
  prop-structure, only in one-phase-region.
*/
{
    /* check pointer */
    if (prop == NULL) { return 0.0; }

    {
	const double T   = prop->T;
	const double d   = prop->d;
	const double cv  = prop->cv;
	const double dpt = prop->dp->T_Cd;
	const double dpd = prop->dp->d_CT;
	return (T * dpt) / (dpd * cv * d * d + T * dpt * dpt);
    }
}

double viscos(Prop *prop)
/*
  Viscosity [Pas] for a calculated prop-structure
  Returns zero if region of validity is exceeded
  (specified by the initial if statement).
  Results do not agree entirely with the NBS/NRC values,
  because an improved formula is used.
  Code from "WaterC" by M.L. Demsey, University of Arizona.
*/
{
    int i, j;
    static const double a[4] = { 
	 0.0181583,  0.0177624,  0.0105287, -0.0036744 };
    static const double b[42]= { 
	 0.5132047,  0.3205656,  0.0,        0.0,
	-0.7782567,  0.1885447,  0.2151778,  0.7317883,
	 1.2410440,  1.4767830,  0.0,        0.0,
	-0.2818107, -1.0707860, -1.2631840,  0.0,
	 0.0,        0.0,        0.1778064,  0.4605040,
	 0.2340379, -0.4924179,  0.0,        0.0,
	-0.0417661,  0.0,        0.0,        0.1600435,
	 0.0,        0.0,        0.0,       -0.01578386,
	 0.0,        0.0,        0.0,        0.0,
	 0.0,        0.0,        0.0,       -0.003629481,
	 0.0,        0.0 };
    static const double tol = 0.01;
    double TdegC, T, D, sum, u0, u1, xt, u2;
    double Pbars, Dkgm3, betaPa;

    /* check pointer */
    if (prop == NULL) { return 0.0; }

    Pbars = prop->p * 1.0e-5;
    Dkgm3 = prop->d;
    betaPa = 1.0 / (prop->d * prop->dp->d_CT);
    TdegC = prop->T - 273.15;

    if ((Pbars > (5000.0 + tol)) 
	|| ((Pbars > (3500.0 + tol)) && (TdegC > (150.0 + tol))) 
	|| ((Pbars > (3000.0 + tol)) && (TdegC > (600.0 + tol))) 
	|| (TdegC > (900.0 + tol))) {
	return 0.0;
    }

    T = prop->T / 647.27;
    D = Dkgm3 / 317.763;
    sum = 0.0;
    for (i = 0; i < 4; i++) {
	sum += a[i] / pow(T, (double) i);
    }
    u0 = 1.0e-6 * sqrt(T) / sum;
    sum = 0.0;
    for (i = 0; i < 6; i++) {
	for (j = 0; j < 7; j++) {
	    sum += b[j * 6 + i] * 
		pow((1.0 / T - 1.0), (double) i) * pow((D - 1.0), (double) j);
	}
    }
    u1 = exp(D * sum);
    if ((0.997 <= T) && (T <= 1.0082) && (0.755 <= D) && (D <= 1.2900)) {
	xt = 22115000.0 / (100973.324169) * betaPa * Dkgm3 * Dkgm3;
	if (xt < 22.0) {
	    u2 = 1.0;
	} else {
	    u2 = 0.922 * pow(xt, 0.0263);
	}
    } else {
	u2 = 1.0;
    }
    return u0 * u1 * u2;
}

double thcond(Prop *prop)
/*
  Thermal conductivity [W/mK] for a calculated prop-structure
  Returns zero if region of validity is exceeded
  (specified by the initial if statement).
  Code from "WaterC" by M.L. Demsey, University of Arizona.
*/
{
    int i, j;
    double T, D, sum, L0, L1, L2, TdegC, u0, u1, xt, dPdT;
    static const double aL[4] = { 
	 2.02223,     14.11166,      
	 5.25597,     -2.0187 };
    static const double au[4] = { 
	 0.0181583,    0.0177624,    
	 0.0105287,   -0.0036744 };
    static const double bL[30] = { 
	 1.329304600, -0.404524370,  0.244094900,
	 0.018660751, -0.129610680,  0.044809953,
	 1.701836300, -2.215684500,  1.651105700,
	-0.767360020,  0.372833440, -0.112031600,
	 5.224615800, -10.12411100,  4.987468700,
	-0.272976940, -0.430833930,  0.133338490,
	 8.712767500, -9.500061100,  4.378660600,
	-0.917837820,  0.0,          0.0,
	-1.852599900,  0.934046900,  0.0,
	 0.0,          0.0,          0.0 };
    static const double bu[30] = { 
	 0.5019380,  0.2356220, -0.2746370,  0.1458310,
	-0.0270448,  0.1628880,  0.7893930, -0.7435390,
	 0.2631290, -0.0253093, -0.1303560,  0.6736650,
	-0.9594560,  0.3472470, -0.0267758,  0.9079190,
	 1.2075520, -0.6873430,  0.2134860, -0.0822904,
	-0.5511190,  0.0670665, -0.4970890,  0.1007540,
	 0.0602253,  0.1465430, -0.0843370,  0.1952860,
	-0.0329320, -0.0202595};
    static const double tol = .01;
    double Pbars, Dkgm3, alph, betaPa;

    /* check pointer */
    if (prop == NULL) { return 0.0; }

    Pbars = prop->p * 1.0e-5;
    Dkgm3 = prop->d;
    alph = prop->dp->T_Cd / prop->dp->d_CT / prop->d;
    betaPa = 1.0 / (prop->d * prop->dp->d_CT);
    TdegC = prop->T - 273.15;

    if ((Pbars > (4000.0 + tol))
        || ((Pbars > (2000.0 + tol)) && (TdegC > (125.0 + tol))) 
        || ((Pbars > (1500.0 + tol)) && (TdegC > (400.0 + tol)))
        || (TdegC > (800.0 + tol))) {
	return 0.0;
    }
    T = prop->T / 647.27;
    D = Dkgm3 / 317.763;
    sum = 0.0;
    for (i = 0; i < 4; i++) {
	sum += aL[i] / pow(T, (double) i);
    }
    L0 = sqrt(T) / sum;
    sum = 0.0;
    for (i = 0; i < 5; i++) {
	for (j = 0; j < 6; j++) {
	    sum += bL[i * 6 + j] * 
		pow((1.0 / T - 1.0), (double) i) * pow((D - 1.0), (double) j);
	}
    }
    L1 = exp(D * sum);
    sum = 0.0;
    for (i = 0; i < 4; i++) {
	sum += au[i] / pow(T, (double) i);
    }
    u0 = 1.0e-6 * sqrt(T) / sum;
    sum = 0.0;
    for (i = 0; i < 6; i++) {
	for (j = 0; j < 5; j++) {
	    sum += bu[i * 5 + j] 
		* pow((1.0 / T - 1.0), (double) i) * pow((D - 1.0), (double) j);
	}
    }
    u1 = exp(D * sum);
    xt = 22115000.0 / (100973.324169) * betaPa * Dkgm3 * Dkgm3;
    dPdT = (647.27 / 22115000.0) * alph / betaPa;
    L2 = 0.000000037711 / (u0 * u1) * T * T/ (D * D) 
	* dPdT * dPdT * pow(xt, 0.4678) * sqrt(D) 
	* exp(-18.66 * (T - 1.0) * (T - 1.0) - pow((D - 1.0), 4.0));
    return L0 * L1 + L2;
}



/*---------------- calculation of one-phase-properties -------------*/

void td(double t, double d, S_mpro *mPro, Prop *prop)
/*
  Computes thermodynamic functions for given temperature and density
  in the one-phase-region
*/
{
    calctd(t, d, mPro);
    format_pro(t, d, mPro, prop);
}

void tp(double t, double p, double *d, double delp, S_mpro *mPro, Prop *prop)
/*
  Computes thermodynamic functions for given temperature and pressure
  in the one-phase-region within a given tolerance of pressure using a
  given approximation of density.
  For t < 647.126 K region of state must be known and *d chosen as follows:
    p < ps(t) => vapour => (*d) < 0.32189 g/cm3
    p > ps(t) => liquid => (*d) > 0.32189 g/cm3
*/
{
    int i = 0;
    double dp, dmin, dmax;

    adjust_tp(t, *d, &dmin, &dmax);
    bb(t, mPro);
    ideal(t, mPro);
    do {
	if ((*d) < dmin) {
	    *d = dmin;
	} else if ((*d) > dmax) {
	    *d = dmax;
	}
	base(t, *d, mPro);
	resid(t, *d, mPro);
	props(t, *d, mPro);
	dp = mPro->spro.p - p;
	if (fabs(dp / p) <= delp) {
	    format_pro(t, *d, mPro, prop);
	    return;
	}
	*d -= dp / mPro->spro.dpd;
	i++;
    }  while (i < IMAX);
    prop->error = 1;
}

void ph(double p, double h, double *t, double *d,
	double delp, double delh, S_mpro *mPro, Prop *prop)
/*
  Computes thermodynamic functions for given enthalpy and pressure in the one-
  phase-region within given tolerances of enthalpy and pressure using given
  approximations of temperature and density (which are changed to yield h,p)
*/
{
    int i = 0;
    double dh, dp, det, deld, delt;
    static const double dhmin = 10.0e3;

    do {
	adjust_hsp(t, d);
	calctd(*t, *d, mPro);
	dh = mPro->spro.h - h;
	dp = mPro->spro.p - p;
	if ((fabs(dh) <= delh * (fabs(h) + dhmin)) && (fabs(dp / p) <= delp)) {
	    format_pro(*t, *d, mPro, prop);
	    return;
	}
	derive(*t, *d, mPro);
	det = mPro->spro.dht * mPro->spro.dpd - mPro->spro.dpt * mPro->spro.dhd;
	delt = (mPro->spro.dpd * dh - mPro->spro.dhd * dp) / det;
	deld = (mPro->spro.dht * dp - mPro->spro.dpt * dh) / det;
	*t -= delt;
	*d -= deld;
	i++;
    } while (i < IMAX);
    prop->error = 1;
}

void ps(double p, double s, double *t, double *d,
	double delp, double dels, S_mpro *mPro, Prop *prop)
/*
  Computes thermodynamic functions for given entropy and pressure in the one-
  phase-region within given tolerances of entropy and pressure using given
  approximations of temperature and density (which are changed to yield s, p)
*/
{
    int i = 0;
    double ds, dp, det, deld, delt;
    static const double dsmin = 0.1;

    do {
	adjust_hsp(t, d);
	calctd(*t, *d, mPro);
	ds = mPro->spro.s - s;
	dp = mPro->spro.p - p;
	if ((fabs(ds) <= dels * (fabs(s) + dsmin)) && (fabs(dp / p) <= delp)) {
	    format_pro(*t, *d, mPro, prop);
	    return;
	}
	derive(*t, *d, mPro);
	det = mPro->spro.dst * mPro->spro.dpd - mPro->spro.dpt * mPro->spro.dsd;
	delt = (mPro->spro.dpd * ds - mPro->spro.dsd * dp) / det;
	deld = (mPro->spro.dst * dp-mPro->spro.dpt * ds) / det;
	*t -= delt;
	*d -= deld;
	i++;
    } while(i < IMAX);
    prop->error = 1;
}

void hd(double h, double d, double *t, double delh, S_mpro *mPro, Prop *prop)
/*
  Computes thermodynamic functions for given enthalpy and density
  in the one-phase-region within a given tolerance of enthalpy using a
  given approximation of temperature.
*/
{
    int i = 0;
    double dh, dht, tmin, tmax;
    static const double dhmin = 100.0;

    adjust_hd(d, &tmin, &tmax);
    do {
	if ((*t) > tmax) {
	    *t = tmax;
	} else if ((*t) < tmin) {
	    *t = tmin;
	}
	calctd(*t, d, mPro);
	dh = mPro->spro.h - h;
	if (fabs(dh) <= delh * (fabs(h) + dhmin)) {
	    format_pro(*t, d, mPro, prop);
	    return;
	}
	dht = mPro->spro.cv + mPro->spro.dpt / d;
	*t -= dh / dht;
	i++;
    } while (i < IMAX);
    prop->error = 1;
}

/* ----------------- extrapolation functions -------------------- */

void extra_ph(double p, double h, double *t, double *d,
	      double dp, double dh, S_mpro *mPro, Prop *prop)
/*
   For invalid (p,h) extrapolates towards (p,h) from a point on
   validity border (specified in valid_ph) closest to (p,h).
   In some areas ph is used to calculate border-point. If interative
   process in ph is not successful, no extrapolation takes place.
*/
{
    double h0, p0, delh, delp;

    /* make sure that first derivatives by p and h will be computed */
    if ((prop->indep1 != 'p') || (prop->indep2 != 'h') || (prop->deriv < 1)) {
	prop->error = 1;
	return;
    }

    if (h > 9950.0) {
	if (p > 3000.0) {
	    /* h = 9950, p = 3000 */
	    td(2635.0940908011, 1.0070832064, mPro, prop);
    } else if (p < P_MIN) {
	    /* h = 9950, p = 1e-6 */
	    td(3122.1018913308, 6.94e-10, mPro, prop);
	} else {
	    ph(p, 9950.0, t, d, dp, dh, mPro, prop);
	}
    } else if (p < tripl.p) {
	if (h > 2460.0) {
        ph(P_MIN, h, t, d, dp, dh, mPro, prop);
    } else if (p < P_MIN) {
	    /* h = 2460, p = 1e-6 */
	    td(251.1302297950, 8.628e-9, mPro, prop);
	} else {
	    ph(p, 2460.0, t, d, dp, dh, mPro, prop);
	}
    } else if (p > 6000.0 - (h + 20.0) / 0.95) {
	if (h > 2830.0) {
	    ph(3000.0, h, t, d, dp, dh, mPro, prop);
	} else {
	    /* h = 2830, p = 3000 */
	    td(432.6356634498, 1.3818099902, mPro, prop);
	}
    } else if (p < 2.0 *tripl.p - (h + 20.0) / 0.95) {
	/* h = -19.9994, p = tripl.p */
	td(268.4580202963, 0.9991594944, mPro, prop);
    } else {
	p0 = 0.5 * (p + (h + 20.0) / 0.95);
	h0 = 0.95 * p0 - 20.0;
	ph(p0, h0, t, d, dp, dh, mPro, prop);
    }
    if (prop->error == 0) {
	p *= 1.0e6;
	h *= 1.0e3;
	delh = h - prop->h;
	delp = p - prop->p;
	prop->T += prop->dT->h_Cp * delh + prop->dT->p_Ch * delp;
	prop->d += prop->dd->h_Cp * delh + prop->dd->p_Ch * delp;
	prop->s += prop->ds->h_Cp * delh + prop->ds->p_Ch * delp;
	prop->u += prop->du->h_Cp * delh + prop->du->p_Ch * delp;
	prop->p = p;
	prop->h = h;
	prop->f = prop->u - prop->T * prop->s;
	prop->g = prop->h - prop->T * prop->s;
    }
}


/* --------------- Calculation of saturation properties -----------*/

void psat(double t, double *p, double *dl, double *dv, 
	  S_mliq *mLiq, S_mpro *mPro)
/*
  Computes liquid and vapour properties for a given saturation temperature
  Liquid properties are put into structure mLiq->spro
  Vapour properties will be found as structure mPro->spro
*/
{
    int i = 0;
    double psa, delg, delp;

    if (t <= creg.t) {
	bb(t, mPro);
	ideal(t, mPro);
	approx_dlv(t, dl, dv);
	do {
	    i++;
	    sat(t, *dl, *dv, &delg, mLiq, mPro);
	    delp = fabs(1.0 - mPro->spro.p / mLiq->spro.p);
	    if ((delp < 0.000001) && (delg < 0.00000001)) {
		i = IMAX;
	    } else {
		psa = (mLiq->spro.f - mPro->spro.f) 
		    / (1.0 / (*dv) - 1.0 / (*dl));
		*dl -= (mLiq->spro.p - psa) / mLiq->spro.dpd;
		*dv -= (mPro->spro.p - psa) / mPro->spro.dpd;
	    }
	} while(i < IMAX);
	*p = 0.5 * (mPro->spro.p + mLiq->spro.p);
    } else {
	psatc(t, p, dl, dv, mLiq, mPro);
    }
}

void psatc(double t, double *p, double *dl, double *dv,
	   S_mliq *mLiq, S_mpro *mPro)
/*  Saturation properties for 646.303775 K < t <= 647.126 K */
{
    double tt, dd, dc;

    bb(t, mPro);
    ideal(t, mPro);
    tt = 1.0 - (t / 647.126);
    dd = 0.657128 * pow(tt, 0.325);
    dc = 0.32189;
    *dl = dc + dd;
    *dv = dc - dd;
    base(t, *dl, mPro);
    resid(t, *dl, mPro);
    props(t, *dl, mPro);
    save(mLiq, mPro);    /* save liquid properties => mLiq->spro */
    base(t, *dv, mPro);
    resid(t, *dv, mPro);
    props(t, *dv, mPro);
    *p = 0.6 * mPro->spro.p + 0.4 * mLiq->spro.p; /* sic! */
}

void tsat(double p, double *t, double *dl, double *dv, 
	  S_mliq *mLiq, S_mpro *mPro)
/*
  Computes liquid and vapour properties for a given saturation pressure
  Liquid properties are put into structure mLiq->spro
  Vapour properties will be found as structure mPro->spro
*/
{
    int i = 0;
    double dt, delpl, delpv, delg;

    if (p <= creg.p) {
	*t = approx_ts(p);
	approx_dlv(*t, dl, dv);
	do {
	    i++;
	    bb(*t, mPro);
	    ideal(*t, mPro);
	    sat(*t, *dl, *dv, &delg, mLiq, mPro);
	    delpl = fabs(1.0 - mLiq->spro.p / p);
	    delpv = fabs(1.0 - mPro->spro.p / p);
	    if ((delpl < 0.000001) && (delpv < 0.000001) 
		&& (delg < 0.00000001)) {
		i = IMAX;
	    } else {
		dt = (mLiq->spro.f - mPro->spro.f 
		      + p * (1.0 / (*dl) - 1.0 / (*dv))) 
		    / (mLiq->spro.s - mPro->spro.s);
		*t += dt;
		*dl += (p - mLiq->spro.p - mLiq->spro.dpt * dt) 
		    / mLiq->spro.dpd;
		*dv += (p - mPro->spro.p - mPro->spro.dpt * dt) 
		    / mPro->spro.dpd;
	    }
	} while(i < IMAX);
    } else {
	tsatc(p, t, dl, dv, mLiq, mPro);
    }
}

void tsatc(double p, double *t, double *dl, double *dv,
	   S_mliq *mLiq, S_mpro *mPro)
/* Saturation properties for 21.839129 MPa <= p <= 22.054915 */
/* Don't touch anything, or we all go to hell ! */

{
    int i = 0;
    static const int imax1 = 10;
    static const int imax2 = 20;
    double dp, psa, dpsdt, t1, t2;

    if (p < 22.05485) {
	*t = 647.126;
	do {
	    i++;
	    psatc(*t, &psa, dl, dv, mLiq, mPro);
	    dp = psa - p;
	    if (fabs(dp) <= 0.000001) {
		i = imax1;
	    } else {
		approx_ps(*t, &dpsdt);
		*t -= dp / dpsdt;
	    }
	} while (i < imax1);
    } else {
	t1 = 647.1259;
	t2 = 647.126;
	do {
	    i++;
	    *t = 0.5 * (t1 + t2);
	    psatc(*t, &psa, dl, dv, mLiq, mPro);
	    if (psa > p) {
		t2 = (*t);
	    } else {
		t1 = (*t);
	    }
	} while((i < imax2) && (fabs(psa - p) > 0.000001));
    }
}

void sat(double t, double dl, double dv, double *delg,
	 S_mliq *mLiq, S_mpro *mPro)
/*
  Computes supposed liquid (mLiq->spro.x) and vapour (mPro->spro.x) properties
  and offers difference of gibbs-function delg = |(gl - gv) / (R T)|
  bb and ideal must have been called previously
*/
{
    base(t, dl, mPro);
    resid(t, dl, mPro);
    props(t, dl, mPro);
    save(mLiq, mPro);	/* save liquid properties => mLiq->spro */
    base(t, dv, mPro);
    resid(t, dv, mPro);
    props(t, dv, mPro);
    *delg = fabs((mLiq->spro.g - mPro->spro.g) / con.gas / t);
}


void hdsat(double h, double d, double delh, double *t,
	   double *p,double *dl, double *dv, double *x,
	   S_mliq *mLiq, S_mpro *mPro)
/*
  find saturation properties for (h, d) in two-phase-region;
  in case of success, enthalpy will be within delh * (h + 100);
  works for tripl.dv <= d <= 1.0 and tripl.t <= t <= creg.t only;

  0 <= x <= 1 : (h, d) is two phase point, use format_two for output,
  x < 0       : (h, d) is liquid point, use hd for calculation
  x > 1       : (h, d) is vapour point, use hd for calculation
*/
{
    int i = 0;
    static const int imax = 100;
    double vl, vv, psa, hll, hvv, dll, dvv, hh, tt, v;
    double cv, dht, cvl, cvv, dpdt, ddpl, ddpv, pss;
    double delp, delg, hx, dh, xx;
    static const double dhmin = 100.0;

    v = 1.0 / d;

    /* guess of t, dl, dv to fit (h, v) */
    *t = approx_thd(h, d);
    approx_dlv(*t, dl, dv);
    vl = 1.0 / (*dl);
    vv = 1.0 / (*dv);

    do {
	i++;

	/* results of the fundamental equation at (t, dl) and (t, dv) */
	calctd(*t, *dl, mPro);
	save(mLiq, mPro);  /* save liquid properties => mLiq->spro */
	base(*t, *dv, mPro);
	resid(*t, *dv, mPro);
	props(*t, *dv, mPro);

	/* examine conditions for two-phase-equilibrium */
	delp = fabs(1.0 - mLiq->spro.p / mPro->spro.p);
	delg = fabs(mLiq->spro.g - mPro->spro.g) / con.gas / (*t);
	*x = (v - vl) / (vv - vl);
	hx = mLiq->spro.h + (*x) * (mPro->spro.h - mLiq->spro.h);
	dh = fabs(hx - h);

	if ((delp < 1e-6) && (delg < 1e-8) && (dh < delh * (fabs(h) + dhmin))) {
	    i = imax;
	} else {
	    /* construct would-be-two-phase-derivation (dh / dt)d */
	    dpdt = (mLiq->spro.s - mPro->spro.s) / (vl - vv);
	    ddpl = dpdt - mLiq->spro.dpt;
	    ddpv = dpdt - mPro->spro.dpt;
	    cvl = mLiq->spro.cv + (*t) * ddpl * ddpl * vl * vl / mLiq->spro.dpd;
	    cvv = mPro->spro.cv + (*t) * ddpv * ddpv * vv * vv / mPro->spro.dpd;
	    cv = cvl + (*x) * (cvv - cvl);
	    dht = cv + v * dpdt;

	    /* better approximations of saturation properties for guess of t */
	    psa = (mLiq->spro.f -mPro->spro.f) / (vv - vl);
	    hll = mLiq->spro.h + (psa - mLiq->spro.p) 
		* vl * (1.0- vl * (*t) * mLiq->spro.dpt / mLiq->spro.dpd);
	    hvv = mPro->spro.h 
		+ (psa - mPro->spro.p) * vv 
		* (1.0 - vv * (*t) * mPro->spro.dpt / mPro->spro.dpd);
	    dll = *dl + (psa - mLiq->spro.p) / mLiq->spro.dpd;
	    dvv = *dv + (psa - mPro->spro.p) / mPro->spro.dpd;

	    /* use approximations and (dh / dt)d to construct better t to 
	       fit (h,v) 
	    */
	    vl = 1.0 / dll;
	    vv = 1.0 / dvv;
	    xx = (v - vl) / (vv - vl);
	    hh = hll + xx * (hvv - hll);
	    tt = (*t) + (h - hh) / dht;
	    if (tt > 646.304) {
		tt = 646.304;
	    } else if (tt < tripl.t) {
		tt = tripl.t;
	    }

	    /* better approximations of saturation properties for new t */
	    pss = psa + dpdt * (tt - (*t));
	    *dl = (*dl) + (pss - mLiq->spro.p 
			   - mLiq->spro.dpt * (tt - (*t))) / mLiq->spro.dpd;

	    *dv = (*dv) + (pss - mPro->spro.p 
			   - mPro->spro.dpt * (tt - (*t))) / mPro->spro.dpd;
	    if ((*dl) < 0.397) {
		*dl = 0.397;
	    } else if ((*dl) > 1.001) {
		*dl = 1.001;
	    }
	    if ((*dv) < 4.85e-6) {
		*dv = 4.85e-6;
	    } else if ((*dv) > 0.247) {
		*dv = 0.247;
	    }
	    vl = 1.0 / (*dl);
	    vv = 1.0 / (*dv);
	    *t = tt;
	}
    } while (i < imax);
    *p = 0.5 * (mPro->spro.p + mLiq->spro.p);
}


void hdsatc(double h, double d, double delh, double *t,
	    double *p, double *dl, double *dv, double *x,
	    S_mliq *mLiq, S_mpro *mPro)
/*
  find saturation properties for (h, d) in critical two-phase-region;
  in case of success, enthalpy will be within delh * (h + 100);
  works for 0.2439 <= d <= 0.4 and creg.t <= t < crit.t only;

  0 <= x <= 1 : (h, d) is two phase point, use format_two for output,
  x < 0       : (h, d) is liquid point, use hd for calculation
  x > 1       : (h, d) is vapour point, use hd for calculation
*/
{
    int i = 0;
    static const int imax = 100;
    double v, vl, vv, hh, dh;
    double dvlt, dvvt, dhlt, dhvt, dhtv;
    static const double dhmin = 100.0;

    v = 1.0 / d;
    *t = 646.303;   /* first guess of t */
    do {
	i++;
	/* calculate two-phase-properties for t */
	psatc(*t, p, dl, dv, mLiq, mPro);
	vl = 1.0 / (*dl);
	vv = 1.0 / (*dv);
	*x = (v - vl) / (vv - vl);
	hh = mLiq->spro.h + (*x) * (mPro->spro.h - mLiq->spro.h);
	/* compare enthalpies */
	dh = h - hh;
	if (fabs(dh) < delh * (fabs(h) + dhmin)) {
	    i = imax;
	} else {
	    /* construct two-phase derivation (dh / dt)v and use it to 
	       improve t */
	    dvlt = vl * 0.325 * (1.0 - 0.32189 * vl) / (crit.t - (*t));
	    dvvt = vv * 0.325 * (1.0 - 0.32189 * vv) / (crit.t - (*t));
	    dhlt = mLiq->spro.cv + vl * mLiq->spro.dpt 
		+ ((*t) * mLiq->spro.dpt - mLiq->spro.dpd / vl) * dvlt;
	    dhvt = mPro->spro.cv + vv * mPro->spro.dpt 
		+ ((*t) * mPro->spro.dpt - mPro->spro.dpd / vv) * dvvt;
	    dhtv = dhlt + (*x) * (dhvt - dhlt) 
		- (mPro->spro.h - mLiq->spro.h) / (vv - vl) 
		* (dvlt + (*x) * (dvvt - dvlt));
	    *t += dh / dhtv;
	    if ((*t) >= crit.t) {
		*t = 647.1259999999999; /* t = 647.126 causes div by zero */
	    } else if ((*t) < 646.303) {
		*t = 646.303;
	    }
	}
    } while (i < imax);
}



/*-------------- approximations for saturation properties ----------*/

double approx_ps(double t, double *dpsdt)
/*
  Approximation of saturation pressure and its derivation by temperature.
  Equation from program code in "NBS/NRC Steam-Tables".
  Below 646.3 K error of ps is within 0.015 % compared to solutions of psat
*/
{
    int i;
    static const double a[8] = { 
	-7.8889166,  2.5514255,
	-6.716169,   33.239495,
	-105.38479,  174.35319, 
	-148.39348,  48.631602};
    double v, w, w2, pa, pb, b, c, psa;

    if (t <= 314.0) {
	pa = 8858.843 / t;
	pb = 607.56335 * pow(t, -0.6);
	psa = 0.1 * exp(6.3573118 - pa + pb);
	*dpsdt = psa * (pa - 0.6 * pb) / t;
    } else {
	v = t / 647.25;
	w = fabs(1.0 - v);
	w2 = sqrt(w);
	b = c = 0.0;
	for (i = 7; i >= 0; i--) {
	    b = b * w2 + a[i];
	    c = c * w2 + a[i] * 0.5 * (i + 2);
	}
	b *= w / v;
	c = -(c + b) / t;
	psa = 22.093 * exp(b);
	*dpsdt = psa * c;
    }
    return psa;
}

double approx_ts(double p)
/*
  Approximation of saturation temperature for given pressure
  (by iterative inversion of approx_ps)
*/
{
    int i = 0;
    static const int imax = 9;
    double pl, t, psa, dpsdt;

    pl = 2.302585 + log(p);
    t = 372.83 
	+ pl * (27.7589 + pl * (2.3819 + pl * (0.24834 + pl * 0.0193855)));
    do {
	i++;
	if (t < 273.15) {
	    t = 273.15;
	} else if (t > 647.126) {
	    t = 647.126;
	}
	psa = approx_ps(t, &dpsdt);
	if (fabs(1.0 - psa / p) < 0.00001) {
	    i = imax;
	} else {
	    t -= (psa - p) / dpsdt;
	}
    } while(i < imax);
    return t;
}


void approx_dlv(double t ,double *dl, double *dv)
/*
  Approximations of liquid and vapour density for given temperature
*/
{
    int i;
    double ts, tt;
    double xl = 0.0;
    double xv = 0.0;
    static const double al1[11] = { 
	 0.3155901296,  -0.03592060636,
	 2.396625841,   -36.39240662,
	 413.6745246,	-2911.342409,
	 12844.66533,	-35543.55367,
	 59925.07856,	-56266.61248,
	 22585.58 };
    static const double av1[11] = { 
	 11.08333753,	-44.64350654,
	 121.0778198,	-278.5153762,
	 327.9464497,	1462.514145,
	-11593.55916,	38922.19673,
	-73229.67131,	74466.37149,
	-31962.35787 };
    static const double al2[10] = { 
	 1.0,	        0.643432247,
	-25.98811457,	192.8271795,
	-947.6312526,	3190.638964,
	-6805.842363,	8088.134131,
	-4034.574025,	0.0 };
    static const double av2[10] = { 
	 1.000000272,    4.026415669,
	-129.9023268,	1867.667009,
	-13845.24815,	61824.71587,
	-171560.6251,	290312.0606,
	-274292.5181,	111202.097 };

    ts = t / 647.3;
    if (t <= 623.15) {
	tt = ts - 273.15 / 647.3;
	for(i = 10; i >= 0; i--) {
	    xl = xl * tt + al1[i];
	    xv = xv * tt + av1[i];
	}
	xv = exp(xv);
    }
    else {
	tt = pow(1.0 - ts, 0.25);
	for (i = 9; i >= 0; i--) {
	    xl = xl * tt + al2[i];
	    xv = xv * tt + av2[i];
	}
    }
    *dl = 1.0 / (3.17 * xl);
    *dv = 1.0 / (3.17 * xv);
}

void approx_hlvp(double p,double *hl, double *hv)
/*
  Approximations for saturation enthalpies for tripl.p <= p <= creg.p.
  Error is within (+-)0.08 kJ/kg compared to solutions of psat.
  For p > 7 MPa error of hv is within (+-)0.03 kJ/kg.
*/
{
    int i;
    double hlr, hvr, pr;
    static const double hl1[10] = { 
	 0.8733214860,           0.3139052092,
	 0.09683295123,          0.02789725252,
	 0.006097758153, 	0.0009363482053,
	 0.00009698796863,	0.000006435851441,
	 0.0000002467695234,	0.000000004155400906 };
    static const double hv1[10] = { 
	 1.191055872,	        -0.237600508,
	-0.1495666133,	        -0.05331935637,
	-0.01282450167,	        -0.002106161251,
	-0.0002309312619,	-0.00001609217758,
	-0.000000642634847,	-0.00000001117696881 };
    static const double hl2[10] = { 

	 0.9506622997,	1.144019695,
	-8.981135832,	22.05395845,
	 10.71497685,	-183.1514846,
	 441.2001644,	-517.5797252,
	 310.7525469,	-76.7287569 };
    static const double hv2[10] = { 
	 0.8378327686,	2.500557338,
	-14.09731747,	46.9525944977,
	-85.16478445,	70.61754229,
	 19.13336011,	-94.08040182,
	 75.79109507,	-21.14194941 };

    hlr = hvr = 0.0;
    if (p < 7.0) {
	pr = log(p / 22.055);
	for(i = 9; i >= 0; i--) {
	    hlr = hlr * pr + hl1[i];
	    hvr = hvr * pr + hv1[i];
	}
    } else {
	pr = pow(1.0 - p / 22.055, 0.25);
	for(i = 9; i >= 0; i--) {
	    hlr = hlr * pr + hl2[i];
	    hvr = hvr * pr + hv2[i];
	}
    }
    *hl = hlr * 2086.0;
    *hv = hvr * 2086.0;
}


double approx_thd(double h, double d)
/*
  Offers approximation of t(h, d) in two-phase-region;
  (for tripl.t < t < creg.t only);
*/
{
    int i = 1;
    static const int imax = 20;
    double t, v, h0, h1, k, dt, h00;
    static const double dhv[20] = {
	0.0121393548,  0.0416114641,   0.1189528301,   0.2935390356,
	0.6423205107,  1.2729584726,   2.3237199593,   3.9607701193,
	6.3740647913,  9.7730565689,   14.3832258357,  20.4443912785,
	28.2120376294, 37.9636079295,  50.0130522183,  64.7397792345,
	82.6458631478, 104.4826727629, 131.6332904572, 169.3249912165 };

    v = 1.0 / d;
    k = 4.1478835;
    dt = 19.63914605;  /* (creg.t - tripl.t) / 19.0 */
    t = tripl.t;
    h00 = -0.0115017594;
    h0 = h00 + dhv[0] * v;     /* start from triple line */
    do {
	h1 = h00 + k * (t + dt - tripl.t) + dhv[i] * v;
	if (h1 > h) {
	    t += dt * (h - h0) / (h1 - h0);
	    i = imax;
	}
	else {
	    t += dt;
	    i++;
	    h0 = h1;
	}
    } while (i < imax);
    return t;
}


/*--------------------- input control ----------------------*/

int valid_td(double t, double d)
/*
   returns 1 if (t, p) falls within IAPS-bounds (fluid region of water)
   (upper bound of d is a bit exceeded - for d > 1.4 water is no longer fluid)
*/
{
    int val = 1;

    if ((t < 260.0) || (t > 2500.0) || (d < DBL_EPSILON) || (d > 1.8)) {
	val = 0;
    }
    return val;
}


int valid_tp(double t, double p)
/*
  returns 1 if (t, p) falls within IAPS-bounds (fluid region of water)
*/
{
    int val = 0;

    if ((t < 260.0) || (t > 2500.0) || (p < P_MIN) || (p > 3000.0)) {
	return val;
    }
    if (p < 209.9) {
	val = ((t > tripl.t) || (p <= psublm(t)) || (p >= pice1(t)));
    } else {
	val = ((t > 413.0) || (p <= pice(t)));
    }
    return val;
}

double psublm(double t)
/*

  sublimation pressure for 190 K <= T <= 273.16 K
*/
{
    double tn, pn, a1, a2, b1, b2, e1, e2, tt;
    
    tn = 273.16;
    pn = 0.000611657;
    a1 = -13.928169;
    a2 = 34.7078238;
    e1 = -1.5;
    e2 = -1.25;
    tt = t / tn;
    b1 = (1.0 - pow(tt, e1));
    b2 = (1.0 - pow(tt, e2));
    return pn * exp(a1 * b1 + a2 * b2);
}

double pice1(double t)
/*
  melting pressure for 251.165 K <= T <= 273.16 K and p < 209.9 MPa
*/
{
    double tt;

    tt = t / 273.16;
    return 0.000611657 * (1.0 - 626000.0 * (1.0 - pow(tt, -3.0)) + 
			  197135.0 * (1.0 - pow(tt, 21.2)));
}

double pice(double t)
/*
  melting pressure for 251.165 K <= T <= 715 K and p > 209.9 MPa
*/
{
    int i;
    double tt;
    static const double tn[4] = { 251.165, 256.164, 273.31, 355.0 };
    static const double pn[3] =  {209.9, 350.1, 632.4 };
    static const double a[3] = { 0.295252, 1.18721, 1.07476 };
    static const double e[3] = { 60.0, 8.0, 4.6 };

    if (t < tn[1]) {
	i = 0;
    } else if (t < tn[2]) {
	i = 1;
    } else if (t < tn[3]) {
	i = 2;
    } else {
	tt = t / tn[3];
	return 2216.0 * exp(1.73683 * (1.0 - 1.0 / tt) - 
			    0.0544606 * (1.0 - pow(tt, 5.0)) + 
			    0.0000000806106 * (1.0 - pow(tt, 22.0)));
    }
    tt = t / tn[i];
    return pn[i] * (1.0 - a[i] * (1.0 - pow(tt, e[i])));
}

int region_tp(double t, double p, double *dl, double *dv, 
	      S_mliq *mLiq, S_mpro *mPro)
/*
  For t < crit.t speed-optimized control of region of (t, p)

  (no control of validity).
  1: p > ps(t) * (1 + 1e-6) => LIQUID
  2: p < ps(t) * (1 - 1e-6) => VAPOUR
  3: p = ps * (1 +- 1e-6)   => SATURATED (d'=dl, d''=dv)
*/
{
    double psa, dpsdt;
    int reg = 0;

    if (t < tripl.t) {
	if (p > tripl.p) {
	    reg = LIQUID;
	} else {
	    reg = VAPOUR;
	}
    } else {
	psa = approx_ps(t, &dpsdt);        /* fast control of region */
	if (t <= creg.t) {
	    if (p > 1.00015 * psa) {
		reg = LIQUID;
	    } else if (p < 0.99985 * psa) {
		reg = VAPOUR;
	    }
	} else {
	    if (p > psa + 0.003) {
		reg = LIQUID;
	    } else if (p < psa - 0.007) {
		reg = VAPOUR;
	    }
	}
	if (reg == 0) {
	    psat(t, &psa, dl, dv, mLiq, mPro);   /* exact control of region */
	    if (t <= creg.t) {
		if (p > 1.000001 * psa) {
		    reg = LIQUID;
		} else if (p < 0.999999 * psa) {
		    reg = VAPOUR;
		} else {
		    reg = SATURATED;
		}
	    }else {
		if (p > mLiq->spro.p) {
		    reg = LIQUID;
		} else if (p < mPro->spro.p) {
		    reg = VAPOUR;

		} else {
		    reg = SATURATED;
		}
	    }
	}
    }
    return reg;
}

int valid_ph(double p, double h)
/*
   Returns 1 if (p, h) is valid.
   A clumsy control of validity of (p, h) with a guarantee, that
   outside two-phase-region only one (t, d) can be found to fit (p, h).
   This region of validity includes the one specified by valid_tp
   except for the ice-near-region above p = 100 MPa.
   It also includes boundaries of (t, d) that result from adjustement
   of (t, d) in adjust_hsp (when used), so that finding (t, d) is possible.
*/
{
    int val = 1;

    if ((p < P_MIN) || (p > 3000.0) || (h > 9950.0) || (h < 0.95 * p - 20.0)) {
	val = 0;
    } else if ((p < tripl.p) && (h < 2460.0)) {
	val = 0;
    }
    return val;
}

int region_ph(double p, double h, double *ts, double *dl, double *dv, 
	      S_mliq *mLiq, S_mpro *mPro)
/*
  checks region of (p, h)
  ONE = not two-phase , TWO = two-phase,
  for two-phase also calculates mLiq->spro properties
*/
{
    double hl, hv;
    int reg;

    if ((p < tripl.p) || (p > crit.p)) {
	reg = ONE;
    } else {
	if (p < creg.p) {
	    approx_hlvp(p, &hl, &hv);  /* approximations of hl, hv */
	} else {
	    hl = 1975.0;
	    hv = 2235.0;
	}
	if ((h < hl - 0.08) || (h > hv + 0.08)) {  /* fast control of region */
	    reg = ONE;
	} else {
	    tsat(p, ts, dl, dv, mLiq, mPro);
	    hl = mLiq->spro.h;
	    hv = mPro->spro.h;
	    if ((h < hl) || (h > hv)) {  /* exact control of region */
		reg = ONE;
	    } else {
		reg = TWO;
	    }
	}
    }
    return reg;
}


int valid_ps(double p, double s)
/*
   Returns 1 if (p, s) is valid.
   A clumsy control of validity of (p, s) with a guarantee, that
   outside two-phase-region only one (t, d) can be found to fit (p, s).
   This region of validity includes the one specified by valid_tp
   except for the ice-near-region above p = 100 MPa.
   It also includes boundaries of (t, d) that result from adjustement
   of (t, d) in adjust_hsp (when used), so that finding (t, d) is possible.
*/
{
    int val = 1;

    if ((p < P_MIN) || (p > 3000.0) || (s < 0.00045 * p - 0.28)) {
	val = 0;
    } else if (s > 11.0 - 0.47 * log(p)) {
	val = 0;
    } else if ((p < tripl.p) && (s < 5.62-0.461 * log(p))) {
	val = 0;
    }
    return val;
}


int valid_hd(double h, double d)
/*
   Returns 1 if (h, d) is valid.
   A clumsy control of validity of (h, d) with a guarantee, that
   outside two-phase-region only one t(h, d) can be found.
   This region of validity includes the one specified by valid_hd
   except for the ice-near-region above p = 100 MPa.
   It also includes boundaries of (t, d) that result from adjustement
   of t in adjust_hd, so that finding t is possible.
*/

{
    int val;
    double v;

    v = 1.0 / d;
    if ((h < 0.000640286) || (h > 9550.0) || (d < DBL_EPSILON) || (v < 0.723)) {
	val = 0;
    } else {
	if (v > 1.0 / tripl.dv) {
	    val = (h >= 2460.0);
	} else if (v > 1.0 / tripl.dl) {
	    val = (h >= -0.0115017594 + v * 0.0121393548); /* triple line */
	} else if (v>0.923) {
	    val = (h >= 1295.2500899117 - v * 1294.9621775858);
	} else if (v>0.882) {
	    val = (h >= 4602.4390243902 - v * 4878.0487804878);
	} else if (v>0.801) {
	    val = (h >= 9011.1111111111 - v * 9876.5432098766);
	} else {
	    val = (h >= 19173.8461538461 - v * 22564.1025641025);
	}
    }
    return val;
}

int region_hd(double h, double d, double dh, double *t,
	      double *p, double *dl, double *dv, double *x,
	      S_mliq *mLiq, S_mpro *mPro)
/*
  checks region of (h, d)
  ONE = not two-phase , TWO = two-phase,
  for two-phase also calculates mLiq->spro properties
  note that dh has an influence on phase classification
*/
{
    int reg = ONE;
    double v;

    v = 1.0 / d;

    if ((v > 1.0) && (v <= 1.0 / tripl.dv) && (h < 2803.3)
        && (h < 2990.0 - 40.0 * log(v))) {
	*x = -1.0;
	if (h < 1547.745404137 + 169.3249912165 * v) {  /* h(v, creg.t) */
	    hdsat(h, d, dh, t, p, dl, dv, x, mLiq, mPro);
	} else if ((v > 2.5) && (v < 4.1)
		   && (h < 1547.2357851199 + 173.4098851329 * v)) { 
	    /* tangent on critical point*/
	    hdsatc(h, d, dh, t, p, dl, dv, x, mLiq, mPro);
	}
	if ((*x >= 0.0) && (*x <= 1.0)) {
	    reg = TWO;
	}
    }
    return reg;
}


/*---------------- supervision of iterations -----------------*/

void adjust_tp(double t, double d, double *dmin, double *dmax)
/*
   Offers maximum and minimum density to adjust d
   during iteration to avoid multiple solutions of d(p, t).
   For t <= crit.t region of state must be known
   (d < 0.32189 means VAPOUR state, LIQUID otherwise)
*/
{
    int reg;
    double da, db;

    *dmin = DBL_EPSILON;
    *dmax = 1.8;
    if (t < crit.t) {
	if (d < crit.d) {
	    reg = VAPOUR;
	} else {
	    reg = LIQUID;
	}
    } else {
	return;
    }
    if (t <= 620.0) {
	if (reg == VAPOUR) {
	    *dmax = exp(-10.1 + 0.013 * t);
	} else {
	    *dmin = 1.2925 - t * 0.00117;
	}
    } else {
	da = 0.32189;
	db = 0.657128 * pow(1.0 - t / 647.126, 0.325) 
	    + 0.03 * (t - 647.126) / 27.126;
	if (reg == VAPOUR) {
	    *dmax = da - db;
	} else {
	    *dmin = da + db;
	}
    }
}

void adjust_hsp(double *t, double *d)
/*
   for systems (h, p) and (s, p) adjusts t and d during iteration
   to avoid multiple solutions
*/
{
    double tmin;

    if ((*d) < DBL_EPSILON) {
	*d = DBL_EPSILON;
    } else if ((*d) > 1.8) {
	*d = 1.8;
    }
    if ((*t) > 3150.0) {
	*t = 3150.0;
    } else if ((*t) < 615.0) {
	if ((*d) < 0.00106) {
	    tmin = 250.0;
	} else if ((*d) < 0.121846) {
	    tmin = (log(*d) + 10.1) / 0.013;
	} else if ((*d) < 0.57295) {
	    tmin = 615.0;
	} else if ((*d) < 1.0) {
	    tmin = (1.2925 - (*d)) / 0.00117;
	} else if ((*d) < 1.147) {
	    tmin = 250.0;
	} else {
	    tmin = ((*d) - 0.702) / 0.00178;
	}
	if ((*t) < tmin) {
	    (*t) = tmin;
	}
    }
}


void adjust_hd(double d, double *tmin, double *tmax)
/*
  offers minimum and maximum values of temperature
  for adjustment during iteration in hd to avoid multiple solutions
*/
{
    *tmax = 3000.0;
    if (d < 0.000170803) {
	*tmin = 250.0;
    } else if (d < 0.24549) {
	*tmin = (log(d) + 13.3) / 0.0185;
    } else if (d < 0.40406) {
	*tmin = 643.0;
    } else if (d < 1.025) {
	*tmin = (1.42 - d) / 0.00158;
    } else if (d < 1.1375) {
	*tmin = 250.0;
    } else {
	*tmin = (d - 0.755) / 0.00153;
    }
}


/*-------------- Calculation and Formation of Properties ---------------*/

void calctd(double t, double d, S_mpro *mPro)
/*
  Calculates properties and places them in structure pro.
*/

{
    bb(t, mPro);
    ideal(t, mPro);
    base(t, d, mPro);
    resid(t, d, mPro);
    props(t, d, mPro);
}

void format_pro (double t, double d, S_mpro *mPro, Prop *prop)
/*
  Formation of output structure prop for one-phase properties
*/
{
    prop->T = t;
    prop->d = d * 1.0e3;
    prop->p = mPro->spro.p * 1.0e6;
    prop->f = mPro->spro.f * 1.0e3;
    prop->g = mPro->spro.g * 1.0e3;
    prop->s = mPro->spro.s * 1.0e3;
    prop->u = mPro->spro.u * 1.0e3;
    prop->h = mPro->spro.h * 1.0e3;
    prop->cv = mPro->spro.cv * 1.0e3;
    prop->cp = mPro->spro.cp * 1.0e3;
    
    if (prop->deriv >= 1) {
      prop->dp->T_Cd = mPro->spro.dpt * 1.0e6;
      prop->dp->d_CT = mPro->spro.dpd * 1.0e3;
    }

    if (prop->deriv == 2) {
	third(t, d, mPro, &mPro->spro);
	prop->dp->dd_CT->d_CT = mPro->spro.dptt * 1.0e6;
	prop->dp->dd_CT->T_Cd = mPro->spro.dptd * 1.0e6;
	prop->dp->dT_Cd->d_CT = prop->dp->dd_CT->T_Cd;
	prop->dp->dT_Cd->T_Cd = mPro->spro.dpdd;
	prop->dcv->T_Cd = mPro->spro.dcvt * 1.0e3;
    }

    if ((prop->indep1 == 'p') && (prop->deriv >= 1)) {
	if (prop->indep2=='h') {
	    deriv_ph(t, d, mPro, prop);
	}
	if (prop->indep2=='s') {
	    deriv_ps(t, d, mPro, prop);
	}
    }
    prop->phase = ONE;
    prop->error = 0;
}

void deriv_ph(double t, double d, S_mpro *mPro, Prop *prop)
/*
  Calculates first and, if required, second derivatives of T, d, s, u
  by p(h = const) and h(p = const) from a formatted prop-struct
  in the one phase region.
*/
{
    double detPH, detPH_t, detPH_d;
    double dtph, dthp, ddph, ddhp;
    double dcv_t, dcv_d, dcp_t, dcp_d;
    double dthp_t, dthp_d, dtph_t, dtph_d;
    double ddhp_t, ddhp_d, ddph_t, ddph_d;
    double duhp_t, duhp_d, duph_t, duph_d;

    detPH = mPro->spro.cp * mPro->spro.dpd;
    mPro->spro.dht = mPro->spro.cv + mPro->spro.dpt / d;
    mPro->spro.dhd = (mPro->spro.dpd - t * mPro->spro.dpt / d) / d;


    dtph = -mPro->spro.dhd / detPH;
    dthp = mPro->spro.dpd / detPH;
    ddph = mPro->spro.dht / detPH;
    ddhp = -mPro->spro.dpt / detPH;

    prop->dT->p_Ch = dtph * 1.0e-6;
    prop->dT->h_Cp = dthp * 1.0e-3;
    prop->dd->p_Ch = ddph * 1.0e-3;
    prop->dd->h_Cp = ddhp;
    prop->ds->p_Ch = -1.0e-3 / t / d;
    prop->ds->h_Cp = 1.0 / t;
    prop->du->p_Ch = 1.0e-3 * (ddph * mPro->spro.p / d - 1.0) / d;
    prop->du->h_Cp = ddhp * mPro->spro.p / d / d + 1.0;

    if (prop->deriv == 2) {

	detPH_t = mPro->spro.dcvt * mPro->spro.dpd 
	    + mPro->spro.cv * mPro->spro.dptd 
	    + (mPro->spro.dpt + 2.0 * t * mPro->spro.dptt) 
	    * mPro->spro.dpt / d / d;
	detPH_d = mPro->spro.cv * mPro->spro.dpdd 
	    + (2.0 * mPro->spro.dpt * (mPro->spro.dptd - mPro->spro.dpt / d)
	       - mPro->spro.dptt * mPro->spro.dpd) * t / d / d;

	mPro->spro.dhtt = mPro->spro.dcvt + mPro->spro.dptt / d;
	mPro->spro.dhtd = (mPro->spro.dptd 
			   - (t * mPro->spro.dptt + mPro->spro.dpt) / d) / d;
	mPro->spro.dhdd = (((2.0 * t * mPro->spro.dpt / d
			     - (t * mPro->spro.dptd + mPro->spro.dpd)) / d
			    + mPro->spro.dpdd)) / d;

	dcv_t = mPro->spro.dcvt;
	dcv_d = -t * mPro->spro.dptt / d / d;

	dcp_t = (detPH_t - mPro->spro.cp * mPro->spro.dptd) / mPro->spro.dpd;
	dcp_d = (detPH_d - mPro->spro.cp * mPro->spro.dpdd) / mPro->spro.dpd;

	dthp_t = dthp * (mPro->spro.dptd / mPro->spro.dpd - detPH_t / detPH);
	dthp_d = dthp * (mPro->spro.dpdd / mPro->spro.dpd - detPH_d / detPH);
	dtph_t = dtph * (mPro->spro.dhtd / mPro->spro.dhd - detPH_t / detPH);
	dtph_d = dtph * (mPro->spro.dhdd / mPro->spro.dhd - detPH_d / detPH);

	ddhp_t = ddhp * (mPro->spro.dptt / mPro->spro.dpt - detPH_t / detPH);
	ddhp_d = ddhp * (mPro->spro.dptd / mPro->spro.dpt - detPH_d / detPH);
	ddph_t = ddph * (mPro->spro.dhtt / mPro->spro.dht - detPH_t / detPH);
	ddph_d = ddph * (mPro->spro.dhtd / mPro->spro.dht - detPH_d / detPH);

	duhp_t = (ddhp * mPro->spro.dpt + mPro->spro.p * ddhp_t) / d / d;
	duhp_d = (ddhp * (mPro->spro.dpd - 2.0 * mPro->spro.p / d)
		  + mPro->spro.p * ddhp_d) / d / d;
	duph_t = (ddph * mPro->spro.dpt + mPro->spro.p * ddph_t) / d / d;
	duph_d = (ddph * (mPro->spro.dpd - 2.0 * mPro->spro.p / d)
		  + mPro->spro.p * ddph_d + 1.0) / d / d;

	prop->dcv->h_Cp =  ddhp * dcv_d + dthp * dcv_t;
	prop->dcv->p_Ch = (ddph * dcv_d + dtph * dcv_t) * 1.0e-3;

	prop->dcp->h_Cp =  ddhp * dcp_d + dthp * dcp_t;
	prop->dcp->p_Ch = (ddph * dcp_d + dtph * dcp_t) * 1.0e-3;

	prop->dT->dh_Cp->h_Cp = (ddhp * dthp_d + dthp * dthp_t) * 1.0e-6;
	prop->dT->dh_Cp->p_Ch = (ddph * dthp_d + dtph * dthp_t) * 1.0e-9;
	prop->dT->dp_Ch->h_Cp = prop->dT->dh_Cp->p_Ch;
	prop->dT->dp_Ch->p_Ch = (ddph * dtph_d + dtph * dtph_t) * 1.0e-12;

	prop->dd->dh_Cp->h_Cp = (ddhp * ddhp_d + dthp * ddhp_t) * 1.0e-3;
	prop->dd->dh_Cp->p_Ch = (ddph * ddhp_d + dtph * ddhp_t) * 1.0e-6;
	prop->dd->dp_Ch->h_Cp = prop->dd->dh_Cp->p_Ch;
	prop->dd->dp_Ch->p_Ch = (ddph * ddph_d + dtph * ddph_t) * 1.0e-9;

	prop->du->dh_Cp->h_Cp = (ddhp * duhp_d + dthp * duhp_t) * 1.0e-3;
	prop->du->dh_Cp->p_Ch = (ddph * duhp_d + dtph * duhp_t) * 1.0e-6;
	prop->du->dp_Ch->h_Cp = prop->du->dh_Cp->p_Ch;
	prop->du->dp_Ch->p_Ch = (ddph * duph_d + dtph * duph_t) * 1.0e-9;

	prop->ds->dh_Cp->h_Cp = - dthp / t /t * 1.0e-3;
	prop->ds->dh_Cp->p_Ch = - dtph / t /t * 1.0e-6;
	prop->ds->dp_Ch->h_Cp = prop->ds->dh_Cp->p_Ch;
	prop->ds->dp_Ch->p_Ch = (ddph / d + dtph / t) / d / t * 1.0e-9;
    }
}

void deriv_ps(double t, double d, S_mpro *mPro, Prop *prop)
/*
  Calculates first derivatives of T, d, h, u by p(s = const) and 
  s(p = const) from a formatted prop-struct in the one phase region.
*/
{
    double detPH, detPH_t, detPH_d;
    double ddsp, ddps, dtsp, dtps;
    double dcv_t, dcv_d, dcp_t, dcp_d;

    detPH = mPro->spro.cp * mPro->spro.dpd;
    mPro->spro.dst = mPro->spro.cv / t;
    mPro->spro.dsd = -mPro->spro.dpt / d / d;

    dtsp = t * mPro->spro.dpd / detPH;
    dtps = -t * mPro->spro.dsd / detPH;
    ddsp = -t * mPro->spro.dpt / detPH;
    ddps = mPro->spro.cv / detPH;

    prop->dT->s_Cp = dtsp * 1.0e-3;
    prop->dT->p_Cs = dtps * 1.0e-6;
    prop->dd->s_Cp = ddsp;
    prop->dd->p_Cs = ddps * 1.0e-3;
    prop->dh->s_Cp = t;
    prop->dh->p_Cs = 1.0e-3 / d;
    prop->du->p_Cs = 1.0e-3 * ddps*mPro->spro.p / d / d;
    prop->du->s_Cp = ddsp * mPro->spro.p / d / d + t;

    if (prop->deriv == 2) {

	detPH_t = mPro->spro.dcvt * mPro->spro.dpd
	    + mPro->spro.cv * mPro->spro.dptd 
	    + (mPro->spro.dpt + 2.0 * t * mPro->spro.dptt) 
	    * mPro->spro.dpt / d / d;
	detPH_d = mPro->spro.cv * mPro->spro.dpdd 
	    + (2.0 * mPro->spro.dpt * (mPro->spro.dptd - mPro->spro.dpt / d)
	       - mPro->spro.dptt * mPro->spro.dpd) * t / d / d;

	dcv_t = mPro->spro.dcvt;
	dcv_d = -t * mPro->spro.dptt / d / d;

	dcp_t = (detPH_t - mPro->spro.cp * mPro->spro.dptd) / mPro->spro.dpd;
	dcp_d = (detPH_d - mPro->spro.cp * mPro->spro.dpdd) / mPro->spro.dpd;

	prop->dcv->s_Cp =  ddsp * dcv_d + dtsp * dcv_t;
	prop->dcv->p_Cs = (ddps * dcv_d + dtps * dcv_t) * 1.0e-3;

	prop->dcp->s_Cp =  ddsp * dcp_d + dtsp * dcp_t;
	prop->dcp->p_Cs = (ddps * dcp_d + dtps * dcp_t) * 1.0e-3;
    }
}


void format_two(double t, double p, double x, double dl, double dv,
		S_mliq *mLiq, S_mpro *mPro, Prop *prop)
/*
  Formation of output structure prop for two-phase properties
*/
{
    double d, dxv, dxd, dxt, dpT, dvTl, dvTv, duTl, duTv, cv;
    double dxdd, dxtd, dxtt, dpTT, dvTTl, dvTTv, duTTl, duTTv, dcvt;

    d = dl * dv / (dv + x * (dl - dv));
    dxv = dl * dv / (dl - dv);
    dxd = -dxv / d / d;
    dpT = (mPro->spro.s - mLiq->spro.s) * dxv;
    dvTl = (mLiq->spro.dpt - dpT) / mLiq->spro.dpd / dl / dl;
    dvTv = (mPro->spro.dpt - dpT) / mPro->spro.dpd / dv / dv;
    dxt = -dxv * (dvTl + x * (dvTv - dvTl));
    duTl = mLiq->spro.cv + (t * mLiq->spro.dpt - mLiq->spro.p) * dvTl;
    duTv = mPro->spro.cv + (t * mPro->spro.dpt - mPro->spro.p) * dvTv;
    cv = duTl + x * (duTv - duTl) + dxt * (mPro->spro.u - mLiq->spro.u);

    prop->x = x;
    prop->T = t;
    prop->d = d * 1.0e3;
    prop->p = p * 1.0e6;
    prop->s = (mLiq->spro.s + x * (mPro->spro.s - mLiq->spro.s)) * 1.0e3;
    prop->u = (mLiq->spro.u + x * (mPro->spro.u - mLiq->spro.u)) * 1.0e3;
    prop->h = (mLiq->spro.h + x * (mPro->spro.h - mLiq->spro.h)) * 1.0e3;
    prop->f = prop->u - t * prop->s;
    prop->g = prop->h - t * prop->s;
    prop->cv = cv * 1.0e3;
    prop->cp = 0.0; /* cp is infinite! */
    if (prop->deriv >= 1) {
	prop->dp->d_CT = 0.0;
	prop->dp->T_Cd = dpT * 1.0e6;
	prop->dx->T_Cd = dxt;
	prop->dx->d_CT = dxd * 1.0e-3;
    }

    if (prop->deriv == 2) {
	/* calculate vapour third derivatives */
	third(t, dv, mPro, &mPro->spro); 
	mPro->mp = mLiq->mp;  /* reload internal results for liquid props */
	mPro->rg = mLiq->rg;
	mPro->glb = mLiq->glb;
	mPro->loc = mLiq->loc;
	/* calculate liquid third derivatives */
	third(t, dl, mPro, &(mLiq->spro));

	dpTT = dxv * (mPro->spro.cv / t - mLiq->spro.cv / t 
		      + dvTv * (mPro->spro.dpt - dpT) 
		      - dvTl * (mLiq->spro.dpt-dpT));
	dxdd = 2.0 * dxv / d / d / d;
	dxtd = dxv * dxv * (dvTv - dvTl) / d / d;
	dvTTl = ((mLiq->spro.dptt - dpTT) / dl / dl 
		 + dvTl * (dl * dvTl 
			   * (2.0 * mLiq->spro.dpd + dl * mLiq->spro.dpdd)
			   - 2.0 * mLiq->spro.dptd)) / mLiq->spro.dpd;
	dvTTv = ((mPro->spro.dptt - dpTT) / dv / dv 
		 + dvTv * (dv * dvTv 
			   * (2.0 * mPro->spro.dpd + dv * mPro->spro.dpdd)
			   - 2.0 * mPro->spro.dptd)) / mPro->spro.dpd;
	dxtt = -dxv * (2.0 * dxt * (dvTv - dvTl) + dvTTl + x * (dvTTv - dvTTl));
	duTTl = mLiq->spro.dcvt 
	    + (mLiq->spro.dpt - dpT 
	       + t * (2.0 * mLiq->spro.dptt 
		      - dl * dl * mLiq->spro.dptd*dvTl)) * dvTl
	    + (t * mLiq->spro.dpt - mLiq->spro.p) * dvTTl;
	duTTv = mPro->spro.dcvt 
	    + (mPro->spro.dpt - dpT 
	       + t * (2.0 * mPro->spro.dptt 
		      - dv * dv * mPro->spro.dptd*dvTv)) * dvTv
	    + (t * mPro->spro.dpt - mPro->spro.p) * dvTTv;
	dcvt = duTTl + x * (duTTv - duTTl) 
	    + 2.0 * dxt * (duTv - duTl) + dxtt * (mPro->spro.u - mLiq->spro.u);

	prop->dp->dd_CT->d_CT = 0.0;
	prop->dp->dT_Cd->d_CT = 0.0;
	prop->dp->dd_CT->T_Cd = 0.0;
	prop->dp->dT_Cd->T_Cd = dpTT * 1.0e6;
	prop->dx->dd_CT->d_CT = dxdd * 1.0e-6;
	prop->dx->dT_Cd->d_CT = dxtd * 1.0e-3;
	prop->dx->dd_CT->T_Cd = prop->dx->dT_Cd->d_CT;
	prop->dx->dT_Cd->T_Cd = dxtt;

	prop->dcv->T_Cd = dcvt * 1.0e3;
    }
    if ((prop->indep1 == 'p') && (prop->deriv >= 1)) {
	if (prop->indep2=='h') {
	    deriv_ph2(prop);
	}
	if (prop->indep2=='s') {
	    deriv_ps2(prop);
	}
    }
    prop->phase = TWO;
    prop->error = 0;
}


void deriv_ph2(Prop *prop)
/*
  Calculate first and, if required, second derivatives of x, T, d, s, u
  by p(h = const) and h(p = const) from a formatted prop-struct
  in the two phase region.
*/
{
    double detPH, detPH_t_PH, detPH_d_PH;
    double dtph, ddph, ddhp;
    double dcv_t, dcv_d;
    double dxph_t, dxph_d;
    double dtph_t;
    double ddhp_t, ddhp_d, ddph_t, ddph_d;
    double duhp_t, duph_t, duph_d;
    double t, d, p, dpt, dptt, dht, dhd, dhtt, dhtd;

    t = prop->T;
    d = prop->d;
    p = prop->p;
    dpt = prop->dp->T_Cd;

    dht = prop->cv + dpt / d;
    dhd = -t * dpt / d / d;
    detPH = -dpt * dhd;

    dtph = 1.0 / dpt;
    ddph = dht / detPH;
    ddhp = -dpt / detPH;

    prop->dx->p_Ch = ddph * prop->dx->d_CT + dtph * prop->dx->T_Cd;
    prop->dx->h_Cp = ddhp * prop->dx->d_CT;
    prop->dT->p_Ch = dtph;
    prop->dT->h_Cp = 0.0;
    prop->dd->p_Ch = ddph;
    prop->dd->h_Cp = ddhp;
    prop->ds->p_Ch = -1.0 / t / d;
    prop->ds->h_Cp = 1.0 / t;
    prop->du->p_Ch = (ddph * p / d - 1.0) / d;
    prop->du->h_Cp = ddhp * p / d / d + 1.0;

    if (prop->deriv == 2) {
	dptt = prop->dp->dT_Cd->T_Cd;

	detPH_t_PH = 2.0 * dptt / dpt + 1.0 / t; /* = detPH_t / detPH */
	detPH_d_PH = -2.0 / d;                   /* = detPH_d / detPH */

	dhtt = prop->dcv->T_Cd + dptt / d;
	dhtd = -(t * dptt + dpt) / d / d;

	dcv_t = prop->dcv->T_Cd;
	dcv_d = -t * dptt / d / d;

	dtph_t = dtph * (dhtd / dhd - detPH_t_PH);

	ddhp_t = ddhp * (dptt / dpt - detPH_t_PH);
	ddhp_d = ddhp * (-detPH_d_PH);
	ddph_t = ddph * (dhtt / dht - detPH_t_PH);
	ddph_d = ddph * (dhtd / dht - detPH_d_PH);

	duhp_t = (ddhp * dpt + p * ddhp_t) / d / d;
	duph_t = (ddph * dpt + p * ddph_t) / d / d;
	duph_d = ((-2.0 * ddph / d + ddph_d) * p + 1.0) / d / d;

	prop->dcv->h_Cp = ddhp * dcv_d;
	prop->dcv->p_Ch = ddph * dcv_d + dtph * dcv_t;

	prop->dcp->h_Cp = 0.0;
	prop->dcp->p_Ch = 0.0;

	prop->dT->dh_Cp->h_Cp = 0.0;
	prop->dT->dh_Cp->p_Ch = 0.0;
	prop->dT->dp_Ch->h_Cp = 0.0;
	prop->dT->dp_Ch->p_Ch = dtph * dtph_t;

	prop->dd->dh_Cp->h_Cp = ddhp * ddhp_d;
	prop->dd->dp_Ch->h_Cp = ddhp * ddph_d;
	prop->dd->dh_Cp->p_Ch = prop->dd->dp_Ch->h_Cp;
	prop->dd->dp_Ch->p_Ch = ddph * ddph_d + dtph * ddph_t;

	prop->du->dh_Cp->h_Cp = 0.0;
	prop->du->dh_Cp->p_Ch = dtph * duhp_t;
	prop->du->dp_Ch->h_Cp = prop->du->dh_Cp->p_Ch;
	prop->du->dp_Ch->p_Ch = ddph * duph_d + dtph * duph_t;

	prop->ds->dh_Cp->h_Cp = 0.0;
	prop->ds->dp_Ch->h_Cp = ddhp / d / d / t;
	prop->ds->dh_Cp->p_Ch = prop->ds->dp_Ch->h_Cp;
	prop->ds->dp_Ch->p_Ch = (ddph / d + dtph / t) / d / t;

	dxph_t =  prop->dx->dT_Cd->d_CT * ddph + prop->dx->d_CT * ddph_t
	    +prop->dx->dT_Cd->T_Cd * dtph + prop->dx->T_Cd * dtph_t;
	dxph_d =  prop->dx->dd_CT->d_CT * ddph + prop->dx->d_CT * ddph_d
	    +prop->dx->dT_Cd->d_CT * dtph;

	prop->dx->dh_Cp->h_Cp = 0.0;
	prop->dx->dp_Ch->h_Cp = ddhp * dxph_d;
	prop->dx->dh_Cp->p_Ch = prop->dx->dp_Ch->h_Cp;
	prop->dx->dp_Ch->p_Ch = ddph * dxph_d + dtph * dxph_t;
    }
}


void deriv_ps2(Prop *prop)
/*
  Calculates first derivatives of x, T, d, h, u by p(s = const) and 
  s(p = const) from a formatted prop-struct in the two phase region.
*/
{
    double detPH, ddsp, ddps, dtps, dcv_t, dcv_d;
    double t, d, p, dpt, dptt;

    t = prop->T;
    d = prop->d;
    p = prop->p;
    dpt = prop->dp->T_Cd;

    detPH = t * dpt * dpt / d / d;
    dtps = 1.0 / dpt;
    ddsp = -t * dpt / detPH;
    ddps = prop->cv / detPH;

    prop->dx->p_Cs = ddps * prop->dx->d_CT + dtps * prop->dx->T_Cd;
    prop->dx->s_Cp = ddsp * prop->dx->d_CT;
    prop->dT->s_Cp = 0.0;
    prop->dT->p_Cs = dtps;
    prop->dd->s_Cp = ddsp;
    prop->dd->p_Cs = ddps;
    prop->dh->s_Cp = t;
    prop->dh->p_Cs = 1.0 / d;
    prop->du->p_Cs = ddps * p / d / d;
    prop->du->s_Cp = ddsp * p / d / d + t;

    if (prop->deriv == 2) {
	dptt = prop->dp->dT_Cd->T_Cd;

	dcv_t = prop->dcv->T_Cd;
	dcv_d = -t * dptt / d / d;

	prop->dcv->s_Cp = ddsp * dcv_d;
	prop->dcv->p_Cs = ddps * dcv_d + dtps * dcv_t;

	prop->dcp->s_Cp = 0.0;
	prop->dcp->p_Cs = 0.0;
    }
}


void save(S_mliq *mLiq, S_mpro *mPro)
/*
  Save liquid properties and some internal results
  for later calculation of liquid third derivatives
*/
{
    mLiq->spro = mPro->spro;
    mLiq->mp = mPro->mp;
    mLiq->rg = mPro->rg;
    mLiq->glb = mPro->glb;
    mLiq->loc = mPro->loc;
}

void load(S_mliq *mLiq, S_mpro *mPro)
/*
  Reload liquid properties and some internal results
  for output of structure pliq in sat_t and sat_p
*/
{
    mPro->spro = mLiq->spro;
    mPro->mp = mLiq->mp;
    mPro->rg = mLiq->rg;
    mPro->glb = mLiq->glb;
    mPro->loc = mLiq->loc;
}


/*---------------------------- IAPS-Formulation ---------------------*/


void props(double t, double d, S_mpro *mPro)
/*
   properties, put together from base, residual and ideal gas function
*/
{
    mPro->spro.f = mPro->bs.f + mPro->rs.f + mPro->id.f;
    mPro->spro.p = d * d * (mPro->bs.fd + mPro->rs.fd);
    mPro->spro.s = -(mPro->bs.ft + mPro->rs.ft + mPro->id.ft);
    mPro->spro.u = mPro->spro.f + t * mPro->spro.s;
    mPro->spro.h = mPro->spro.u + mPro->spro.p / d;
    mPro->spro.g = mPro->spro.f + mPro->spro.p / d;
    mPro->spro.dpd = 2.0 * mPro->spro.p / d 
	+ d * d * (mPro->bs.fdd + mPro->rs.fdd);
    mPro->spro.dpt = d * d * (mPro->bs.ftd + mPro->rs.ftd);
    mPro->spro.cv = -t * (mPro->bs.ftt + mPro->rs.ftt + mPro->id.ftt);
    mPro->spro.cp = mPro->spro.cv + 
	t * mPro->spro.dpt * mPro->spro.dpt / (d * d * mPro->spro.dpd);
}

void derive(double t, double d, S_mpro *mPro)
/*
  derivatives of properties, useful for numerical methods
  props must be called previously.
  Pressure derivations are computed in props already.
*/
{
    mPro->spro.dft = -mPro->spro.s;
    mPro->spro.dfd = mPro->spro.p / d / d;
    mPro->spro.dgt = mPro->spro.dpt / d - mPro->spro.s;
    mPro->spro.dgd = mPro->spro.dpd / d;
    mPro->spro.dst = mPro->spro.cv / t;
    mPro->spro.dsd = -mPro->spro.dpt / d / d;
    mPro->spro.dut = mPro->spro.cv;
    mPro->spro.dud = (mPro->spro.p - t * mPro->spro.dpt) / d / d;
    mPro->spro.dht = mPro->spro.cv + mPro->spro.dpt / d;
    mPro->spro.dhd = (mPro->spro.dpd - t * mPro->spro.dpt / d) / d;
}

void bb(double t, S_mpro *mPro)
/*
  molecular parameters
*/
{
    double tt, tt2, tt3;

    tt = con.tz / t;
    tt2 = tt * tt;
    tt3 = tt2 * tt;
    mPro->mp.b1 = Cm.b1[0] + Cm.b1[1] * log(1.0 / tt) 
	+ (Cm.b1[2] + Cm.b1[3] * tt2) * tt3;
    mPro->mp.b1t = (Cm.b1[1] - (3.0 * Cm.b1[2] + 5.0 * Cm.b1[3] * tt2) 
		    * tt3) / t;
    mPro->mp.b1tt = (-Cm.b1[1] 
		     + (12.0 * Cm.b1[2] + 30.0 * Cm.b1[3] * tt2) * tt3) / t / t;
    mPro->mp.b2 = Cm.b2[0] + (Cm.b2[1] + (Cm.b2[2] + Cm.b2[3] * tt2) * tt) * tt;
    mPro->mp.b2t = -(Cm.b2[1] 
		     + (2.0 * Cm.b2[2] + 4.0 * Cm.b2[3] * tt2) * tt) * tt / t;
    mPro->mp.b2tt = (2.0 * Cm.b2[1] 
		     + (6.0 * Cm.b2[2] + 20.0 * Cm.b2[3] * tt2) 
		     * tt) * tt / t / t;
}

void base(double t, double d, S_mpro *mPro)
/*
  Base function and its first and second derivatives
*/
{
    double x;

    x = 1.0 - 0.25 * mPro->mp.b1 * d;
    mPro->rg.z1 = -log(x) + (91.0 + (-260.0 + 169.0 / x) / x) / 6.0;
    mPro->rg.z2 = (1.0 + (-130.0 + 169.0 / x) / x / 3.0) / x;
    mPro->rg.z3 = (1.0 + (-260.0 / 3.0 + 169.0 / x) / x) / x / x;

    mPro->rg.k = mPro->mp.b2 - 3.5 * mPro->mp.b1;
    mPro->rg.kt = mPro->mp.b2t - 3.5 * mPro->mp.b1t;
    mPro->rg.ktt = mPro->mp.b2tt - 3.5 * mPro->mp.b1tt;

    mPro->bs.f  = con.gas * t 
	* (mPro->rg.z1 + d * mPro->rg.k + log(d * con.gas * t / con.pz));
    mPro->bs.fd = con.gas * t 
	* (mPro->rg.z2 * mPro->mp.b1 / 4.0 + mPro->rg.k + 1.0 / d);
    mPro->bs.fdd = con.gas * t 
	* (mPro->rg.z3 * mPro->mp.b1 * mPro->mp.b1 / 16.0 - 1.0 / d / d);
    mPro->bs.ft = mPro->bs.f / t + con.gas 
	* ((t * mPro->rg.z2 * mPro->mp.b1t / 4.0 + t * mPro->rg.kt) * d + 1.0);
    mPro->bs.ftt = con.gas 
	* ((((mPro->mp.b1t + t * mPro->mp.b1tt / 2.0) * mPro->rg.z2 
	     + t * mPro->rg.z3 * mPro->mp.b1t *mPro->mp.b1t
	     * d / 8.0) / 2.0 + (2.0 * mPro->rg.kt + t * mPro->rg.ktt)) 
	   * d + 1.0 / t);
    mPro->bs.ftd = con.gas 
	* (((mPro->mp.b1 + t * mPro->mp.b1t) * mPro->rg.z2 
	    + t * mPro->rg.z3 * mPro->mp.b1 * mPro->mp.b1t * d / 4.0) / 4.0
	   + (mPro->rg.k + t * mPro->rg.kt) + 1.0 / d);
}


void resid(double t, double d, S_mpro *mPro)
/*
  Residual function and its first and second derivatives
*/
{
    int j;
    double ed, dd, tt;
    double tau, del, dk, dl;

    tt = con.tz / t;
    ed = exp(-d);
    dd = (1.0 - ed);

    mPro->rs.f = mPro->rs.ft = mPro->rs.ftt = 0.0;
    mPro->rs.ftd = mPro->rs.fd = mPro->rs.fdd = 0.0;

    /* global parts */

    for (j = 8; j > -1; j--) {
	mPro->glb.k[j] = Cr.g[j][0] 
	    + (Cr.g[j][1] + (Cr.g[j][2] + (Cr.g[j][3]
	    + (Cr.g[j][4] + Cr.g[j][5] * tt * tt)
	    * tt) * tt) * tt) * tt;
	mPro->glb.kt[j] = (Cr.g[j][1] + (2.0 * Cr.g[j][2] + (3.0 * Cr.g[j][3]
	    + (4.0 * Cr.g[j][4] + 6.0 * Cr.g[j][5] * tt * tt) 
            * tt) * tt) * tt) * tt;
	mPro->glb.ktt[j] = (2.0 * Cr.g[j][1] + (6.0 * Cr.g[j][2] 
	    + (12.0 * Cr.g[j][3] + (20.0 * Cr.g[j][4] + 42.0 * Cr.g[j][5] * tt
            * tt) * tt) * tt) * tt) * tt;

	mPro->rs.f   = mPro->rs.f * dd + mPro->glb.k[j] / (j + 1.0);
	mPro->rs.fd  = mPro->rs.fd * dd + mPro->glb.k[j];
	mPro->rs.fdd = mPro->rs.fdd * dd + mPro->glb.k[j] * (j * ed / dd - 1.0);
	mPro->rs.ft  = mPro->rs.ft * dd + mPro->glb.kt[j] / (j + 1.0);
	mPro->rs.ftt = mPro->rs.ftt * dd + mPro->glb.ktt[j] / (j + 1.0);
	mPro->rs.ftd = mPro->rs.ftd * dd + mPro->glb.kt[j];
    }
    mPro->rs.f *= dd;
    mPro->rs.fd *= ed;
    mPro->rs.fdd *= ed;
    mPro->rs.ft  *= -dd / t;
    mPro->rs.ftt *= dd / t / t;
    mPro->rs.ftd *= -ed / t;

    /* local parts */

    for (j = 0; j < 4; j++) {
	mPro->loc.tau[j] = tau = (t - Cr.t[j]) / Cr.t[j];
	mPro->loc.del[j] = del = (d - Cr.d[j]) / Cr.d[j];
	if (fabs(del)<DBL_EPSILON) {
	    mPro->loc.del[j] = del = DBL_EPSILON;  /* avoid division by zero */
	}
	mPro->loc.dk[j] = dk = pow(del, Cr.k[j]);
	mPro->loc.dl[j] = dl = pow(del, Cr.l[j]);

	mPro->loc.k[j] = Cr.gg[j] * dl 
	    * exp(-Cr.a[j] * dk - Cr.b[j] * tau * tau);
	mPro->loc.kt[j] = -2.0 * Cr.b[j] * tau / Cr.t[j];
	mPro->loc.ktt[j] = (4.0 * tau * tau * Cr.b[j] - 2.0) 
	    * Cr.b[j] / Cr.t[j] / Cr.t[j];
	mPro->loc.kd[j] = (Cr.l[j] - Cr.a[j] * Cr.k[j] * dk) / Cr.d[j] / del;
	mPro->loc.kdd[j] = (Cr.k[j] * Cr.a[j] * dk * (1.0 - Cr.k[j]) 
			    - Cr.l[j]) / Cr.d[j] / Cr.d[j] / del / del;

	mPro->rs.f   += mPro->loc.k[j];
	mPro->rs.ft  += mPro->loc.k[j] * mPro->loc.kt[j];
	mPro->rs.ftt += mPro->loc.k[j] * mPro->loc.ktt[j];
	mPro->rs.ftd += mPro->loc.k[j] * mPro->loc.kt[j] * mPro->loc.kd[j];
	mPro->rs.fd  += mPro->loc.k[j] * mPro->loc.kd[j];
	mPro->rs.fdd += mPro->loc.k[j] * (mPro->loc.kdd[j] + mPro->loc.kd[j] 
					  * mPro->loc.kd[j]);
    }
}

void ideal(double t, S_mpro *mPro)
/*
  Ideal gas function and its first and second derivatives
*/
{
    double tr, tl;

    tr = t / 100.0;
    tl = log(tr);

    mPro->id.f = -con.gas *
	(t * ((Ci.c[0] / tr + Ci.c[1]) * tl
	    + ((Ci.c[2] / tr + Ci.c[3]) / tr + Ci.c[4]) / tr + Ci.c[5]
	    + (Ci.c[6] + (Ci.c[7] + (Ci.c[8] + (Ci.c[9] + (Ci.c[10] + (Ci.c[11]
	    + (Ci.c[12] + (Ci.c[13] + (Ci.c[14] + (Ci.c[15] + (Ci.c[16] 
            + Ci.c[17] * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr)
			  * tr) * tr) * tr + 1.0 - con.sref) + con.uref);

    mPro->id.ft = -con.gas * (Ci.c[0] / tr + Ci.c[1] * (tl + 1.0)
	 - (2.0 * Ci.c[2] / tr + Ci.c[3]) / tr / tr + Ci.c[5]
	 + (2.0 * Ci.c[6] + (3.0 * Ci.c[7] + (4.0 * Ci.c[8] + (5.0 * Ci.c[9]
	 + (6.0 * Ci.c[10] + (7.0 * Ci.c[11] + (8.0 * Ci.c[12] 
         + (9.0 * Ci.c[13] + (10.0 * Ci.c[14] + (11.0 * Ci.c[15] 
         + (12.0 * Ci.c[16] + 13.0 * Ci.c[17] * tr) * tr) * tr) * tr) * tr) 
         * tr) * tr) * tr) * tr) * tr) * tr) * tr + 1.0 - con.sref);

    mPro->id.ftt = -con.gas / t * (Ci.c[1] - Ci.c[0] / tr
	 + (6.0 * Ci.c[2] / tr + 2.0 * Ci.c[3]) / tr / tr
	 + (2.0 * Ci.c[6] + (6.0 * Ci.c[7] + (12.0 * Ci.c[8] 
         + (20.0 * Ci.c[9] + (30.0 * Ci.c[10] + (42.0 * Ci.c[11]
         + (56.0 * Ci.c[12] + (72.0 * Ci.c[13] + (90.0 * Ci.c[14]
         + (110.0 * Ci.c[15] + (132.0 * Ci.c[16] + 156.0 * Ci.c[17]
         * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) 
         * tr) * tr);
}


void third(double t, double d, S_mpro *mPro, S_pro *out)
/*

  Third derivatives of the entire fundamental equation.
  Uses internal results of bb(), base() and resid().
*/
{
    int j;
    double x, tt, ed, dd, kttt, kddd, tr;

    /* molecular parameters */

    tt = con.tz / t;
    tr = t / 100.0;

    mPro->mp.b1ttt = (2.0 * Cm.b1[1] - (60.0 * Cm.b1[2] 
                     + 210.0 * Cm.b1[3] * tt * tt) * tt * tt * tt) / t / t / t;
    mPro->mp.b2ttt = -(6.0 * Cm.b2[1] + (24.0 * Cm.b2[2] 
                     + 120.0 * Cm.b2[3] * tt * tt) * tt) * tt / t / t / t;

    /* base function */

    x = 1.0 - 0.25 * mPro->mp.b1 * d;
    mPro->rg.z4 = (2.0 + (-260.0 + 676.0 / x) / x) / x / x / x;
    mPro->rg.kttt = mPro-> mp.b2ttt - 3.5 * mPro->mp.b1ttt;

    mPro->bs.fddd = con.gas * t * (mPro->rg.z4 * mPro->mp.b1 
		      * mPro->mp.b1 * mPro->mp.b1 
		      / 64.0 + 2.0 / (d * d * d));
    mPro->bs.ftdd = con.gas * (((mPro->mp.b1 / 2.0 
		      + t * mPro->mp.b1t) * mPro->rg.z3
		      + t * mPro->rg.z4 * mPro->mp.b1 * mPro->mp.b1t * d / 8.0) 
		      * mPro->mp.b1 / 8.0 - 1.0 / (d * d));
    mPro->bs.fttd = con.gas * (((mPro->mp.b1t + t * mPro->mp.b1tt/2.0) 
		      * mPro->rg.z2 + (((mPro->mp.b1 + t * mPro->mp.b1t) 
		      * mPro->mp.b1t + t * mPro->mp.b1 * mPro->mp.b1tt / 2.0) 
		      * mPro->rg.z3 + t * mPro->rg.z4 * mPro->mp.b1 
		      * mPro->mp.b1t * mPro->mp.b1t * d / 8.0) * d / 4.0) / 2.0 
	              + (2.0 * mPro->rg.kt + t * mPro->rg.ktt));
    mPro->bs.fttt = con.gas *
	(((3.0 * mPro->mp.b1tt + t * mPro->mp.b1ttt) * mPro->rg.z2
	  + (3.0 * (mPro->mp.b1t + t * mPro->mp.b1tt)
	     * mPro->rg.z3 + t * mPro->rg.z4 * mPro->mp.b1t * 
	     mPro->mp.b1t * d / 4.0) * mPro->mp.b1t * d / 4.0) * d / 4.0
	 + (3.0 * mPro->rg.ktt + t * mPro->rg.kttt) * d - 1.0 / (t * t));

    /* residual function : global parts */

    tt = con.tz / t;
    ed = exp(-d);
    dd = (1.0 - ed);

    mPro->rs.fttt = mPro->rs.fttd = mPro->rs.ftdd = mPro->rs.fddd = 0.0;

    for (j = 8; j > -1; j--) {

	kttt = (6.0 * Cr.g[j][1] + (24.0 * Cr.g[j][2] + 
				   (60.0 * Cr.g[j][3] + 
				    (120.0 * Cr.g[j][4]
				     + 336.0 * Cr.g[j][5] * tt * tt) * tt)
				   * tt) * tt) *tt;

	mPro->rs.fttt = mPro->rs.fttt*dd + kttt / (j + 1.0);
	mPro->rs.fttd = mPro->rs.fttd*dd + mPro->glb.ktt[j];
	mPro->rs.ftdd = mPro->rs.ftdd*dd 
	    + mPro->glb.kt[j] * (j * ed / dd - 1.0);
	mPro->rs.fddd = mPro->rs.fddd*dd + 
	    mPro->glb.k[j] * (1.0 + j * ed / dd * ((j - 1.0) * ed / dd - 3.0));
    }
    mPro->rs.fttt *= -dd / t / t / t;
    mPro->rs.fttd *= ed / t / t;
    mPro->rs.ftdd *= -ed / t;
    mPro->rs.fddd *= ed;

    /* residual function : local parts */

    for(j = 0; j < 4; j++) {

	kddd = (Cr.a[j] * Cr.k[j] * mPro->loc.dk[j] 
		* (1.0 - Cr.k[j]) * (Cr.k[j] - 2.0) + 2.0 * Cr.l[j])
	    / Cr.d[j] / Cr.d[j] / Cr.d[j] / mPro->loc.del[j] / 
	    mPro->loc.del[j] / mPro->loc.del[j];
	mPro->rs.fddd += mPro->loc.k[j] * 
	    (kddd + mPro->loc.kd[j] * (3.0 * mPro->loc.kdd[j] 
				       + mPro->loc.kd[j] * mPro->loc.kd[j]));
	mPro->rs.fttt += mPro->loc.k[j] * (3.0 - 2.0 * Cr.b[j] 
		   * mPro->loc.tau[j] * mPro->loc.tau[j])
	    * 4.0 * Cr.b[j] * Cr.b[j] * mPro->loc.tau[j] 
	    / Cr.t[j] / Cr.t[j] / Cr.t[j];
	mPro->rs.fttd += mPro->loc.k[j] * mPro->loc.ktt[j] * mPro->loc.kd[j];
	mPro->rs.ftdd += mPro->loc.k[j] * mPro->loc.kt[j] * 
	    (mPro->loc.kdd[j] + mPro->loc.kd[j] * mPro->loc.kd[j]);
    }

    /* ideal function */

    mPro->id.fttt = -con.gas / t / t *
	(2.0 * Ci.c[0]/tr - Ci.c[1]
	 - (24.0 * Ci.c[2] / tr + 6.0 * Ci.c[3]) / tr / tr
	 + (6.0 * Ci.c[7] + (24.0 * Ci.c[8] +
			       (60.0 * Ci.c[9] + (120.0 * Ci.c[10]
		   + (210.0 * Ci.c[11]+ (336.0 * Ci.c[12] + 
		 (504.0 * Ci.c[13] + (720.0 * Ci.c[14]
	      + (990.0 * Ci.c[15] + (1320.0 * Ci.c[16] + 1716.0 * Ci.c[17]
	     * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) 
	      * tr) * tr * tr);

    /* properties' derivatives */

    out->dpdd = -6.0 * out->p / d / d + 4.0 * out->dpd / d 
	+ d * d * (mPro->bs.fddd + mPro->rs.fddd);
    out->dptd = 2.0 * out->dpt / d 
	+ d * d * (mPro->bs.ftdd + mPro->rs.ftdd);
    out->dptt = d * d * (mPro->bs.fttd + mPro->rs.fttd);
    out->dcvt = out->cv / t 
	- t * (mPro->bs.fttt + mPro->rs.fttt + mPro->id.fttt);
}
