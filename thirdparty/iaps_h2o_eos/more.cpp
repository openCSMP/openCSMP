// SPDX-FileCopyrightText: © 1998 Olaf Bauer <o.bauer@tu-harburg.de>
// SPDX-FileCopyrightText: © 1998 Ole Engel <engel@tu-harburg.de>
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: more.cpp,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Olaf Bauer
#             Ole Engel
# maintainer: Ole Engel <engel@tu-harburg.de>
#---------------------------------------------------------------------
# Part of PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/

/*
 * Modified by "The openCSMP project", 2026:
 * - Ported from more.c to more.cpp for C++ compilation
 * - Replaced C headers (stdio.h, stdlib.h, math.h, float.h) with
 *   C++ equivalents (cstdio, cstdlib, cmath, cfloat)
 * - Use of namespace std
 * - Fixed character encoding in comment (line 579)
 *   Gegenma\DFnahme -> Gegenmassnahme
 */

/* all possibilities:
       1. in iaps.c available: td, tp, hp, sp, hd
       2. in more.c available: st, ht, hs, pd, ud, sd, pu, us, ut
       3. not used:            ft, gt, fd, gd, pf, pg, fg, fh, fu, fs,
                   gh, gu, gs, hu,

   double names:
       pt, ph, ps
*/
#include <cstdio>   /* printf() */
#include <cstdlib>  /* rand(), RAND_MAX */
#include <cmath>    /* fabs() */
#include <cfloat>   /* DBL_EPSILON */
#include "iaps.h"
#include "more.h"

using namespace std;

/* external functions */

void water_th(double t, double h, double d, double dh, Prop *prop)
/*
  Computes properties for given temperature and enthalpy within a given
  tolerance of enthalpy using a given approximation of density.
*/
{
    double hl, hv, dl, dv, x, p;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    h *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_th(t, h)) {
    prop->error = 1;
    return;
    }
    if (t <= crit.t) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    hl = MLiq.spro.h;
    hv = MPro.spro.h;
    if ((h > hl) && (h < hv)) {
        x = (h - hl) / (hv - hl);
        format_two(t, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }
    }
    ht(h, t, &d, dh, &MPro, prop);
}

void water_tu(double t, double u, double d, double du, Prop *prop)
/*
  Computes properties for given temperature and internal energy within a
  given tolerance of internal energy using a given approximation of density.
*/
{
    double ul, uv, dl, dv, x, p;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    u *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_tu(t, u)) {
    prop->error = 1;
    return;
    }
    if (t <= crit.t) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    ul = MLiq.spro.u;
    uv = MPro.spro.u;
    if ((u > ul) && (u < uv)) {
        x = (u - ul) / (uv - ul);
        format_two(t, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }
    }
    ut(u, t, &d, du, &MPro, prop);
}

void water_ts(double t, double s, double d, double ds, Prop *prop)
/*
  Computes properties for given temperature and entropy within a given
  tolerance of entropy using a given approximation of density.
*/
{
    double sl, sv, dl, dv, x, p;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    s *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_ts(t, s)) {
    prop->error = 1;
    return;
    }
    if (t <= crit.t) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    sl = MLiq.spro.s;
    sv = MPro.spro.s;
    if ((s > sl) && (s < sv)) {
        x = (s - sl) / (sv - sl);
        format_two(t, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }
    }
    st(s, t, &d, ds, &MPro, prop);
}

void water_hs(double h, double s, double t, double d,
        double dh, double ds, Prop *prop)
/*
  Computes properties for given enthalpy and entropy within given
  tolerances of enthalpy and entropy using given approximations
  of temperature and density.
*/
{
    double p, dl, dv;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    h *= 1.0e-3;
    s *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_hs(h, s)) {
    prop->error = 1;
    return;
    }
    hs(h, s, &t, &d, dh, ds, &MPro, prop);
    if (t <= crit.t) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    if ((d < dl) && (d > dv)) {
        prop->error = 1;
    }
    }
}

void water_pu(double p, double u, double t, double d,
        double dp, double du, Prop *prop)
/*
  Computes properties for given pressure and internal energy within given
  tolerances of pressure and internal energy using given approximations
  of temperature and density.
*/
{
    double ts, dl, dv, ul, uv, x;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    u *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_pu(p, u)) {
    prop->error = 1;
    return;
    }
    if ((p >= tripl.p) && (p <= crit.p)) {
    tsat(p, &ts, &dl, &dv, &MLiq, &MPro);
    ul = MLiq.spro.u;
    uv = MPro.spro.u;
    if ((u > ul) && (u < uv)) {              /* exact control of region */
        x = (u - ul) / (uv - ul);
        format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }
    }
    pu(p, u, &t, &d, dp, du, &MPro, prop);
}

void water_us(double u, double s, double t, double d,
        double du, double ds, Prop *prop)
/*
  Computes properties for given internal energy and entropy within given
  tolerances of internal energy and entropy using given approximations
  of temperature and density.
*/
{
    double p, dl, dv;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    u *= 1.0e-3;
    s *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_us(u, s)) {
    prop->error = 1;
    return;
    }
    us(u, s, &t, &d, du, ds, &MPro, prop);
    if (t <= crit.t) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    if ((d < dl) && (d > dv)) {
        prop->error = 1;
    }
    }
}


void water_pd(double p, double d, double t, double dp, Prop *prop)
/*
  Computes properties for given pressure and density within given
  tolerance of pressure using given approximation
  of temperature.
*/
{
    double ts, dl, dv, x, v, vl, vv;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    d *= 1.0e-3;
    if (! valid_pd(p, d)) {
    prop->error = 1;
    return;
    }
    if ((p >= tripl.p) && (p <= crit.p)) {
    tsat(p, &ts, &dl, &dv, &MLiq, &MPro);
    v = 1.0 / d;
    vl = 1.0 / dl;
    vv = 1.0 / dv;
    if ((v > vl) && (v < vv)) {              /* exact control of region */
        x = (v - vl) / (vv - vl);
        format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }
    }
    pd(p, d, &t, dp, &MPro, prop);
}

void water_ud(double u, double d, double t, double du, Prop *prop)
/*
  Computes properties for given internal energy and density within given
  tolerance of internal energy using given approximation
  of temperature.
*/
{
    double p, dl, dv;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    u *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_ud(u, d)) {
    prop->error = 1;
    return;
    }
    ud(u, d, &t, du, &MPro, prop);
    if (t <= crit.t) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    if ((d < dl) && (d > dv)) {
        prop->error = 1;
    }
    }
}

void water_sd(double s, double d, double t, double ds, Prop *prop)
/*
  Computes properties for given entropy and density within given
  tolerance of entropy using given approximation
  of temperature.
*/
{
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    s *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_sd(s, d)) {
    prop->error = 1;
    return;
    }
    sd(s, d, &t, ds, prop, &MLiq, &MPro);
}

/* ------------ validate functions -------------- */

int valid_th(double t, double h)
{
    if (t) {} /* avoid 'unused parameter' warning */
    if (h) {} /* avoid 'unused parameter' warning */
    return 1;
}

int valid_tu(double t, double u)
{
    if (t) {} /* avoid 'unused parameter' warning */
    if (u) {} /* avoid 'unused parameter' warning */
    return 1;
}

int valid_ts(double t, double s)
{
    if (t) {} /* avoid 'unused parameter' warning */
    if (s) {} /* avoid 'unused parameter' warning */
    return 1;
}


int valid_pd(double p, double d)
{
    if (p) {} /* avoid 'unused parameter' warning */
    if (d) {} /* avoid 'unused parameter' warning */
    return 1;
}

int valid_ud(double u, double d)
{
    if (u) {} /* avoid 'unused parameter' warning */
    if (d) {} /* avoid 'unused parameter' warning */
    return 1;
}

int valid_us(double u, double s)
{
    if (u) {} /* avoid 'unused parameter' warning */
    if (s) {} /* avoid 'unused parameter' warning */
    return 1;
}

int valid_sd(double s, double d)
{
    if (s) {} /* avoid 'unused parameter' warning */
    if (d) {} /* avoid 'unused parameter' warning */
   return 1;
}

int valid_pu(double p, double u)
{
    if (p) {} /* avoid 'unused parameter' warning */
    if (u) {} /* avoid 'unused parameter' warning */
    return 1;
}

/*
  for adding more inversion functions

  a)  x = (p, f, g, h, u, s) and t are known, d must be found
  b)  x = (p, f, g, h, u, s) and d are known, t must be found
  c)  x = (p, f, g, h, u, s) and y = (p, f, g, h, u, s) are known,
      t and d be must be found

*/

/*--------------------------------------------------------------*/
/* a)  x = (p, f, g, h, u, s) and t are known, d must be found  */
/*--------------------------------------------------------------*/

/* replace x by known parameter: */

/*
void xt(double x, double t, double *d, double delx, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double dx;

    bb(t, mPro);
    ideal(t, mPro);
    do {
    base(t, *d, mPro);
    resid(t, *d, mPro);
    props(t, *d. mPro);
    dx = mPro->spro.x - x;
    if (fabs(dx) <= delx * fabs(x)) {
        format_pro(t, *d, mPro, prop);
        return;
    }
        *//* call of derive() is not necessary for x = p *//*
    derive(t, *d, mPro);
        *d -= dx / mPro->spro.dxd;
        i++;
    }  while (i < IMAX);
    prop->error = 1;
}
*/

void ht(double h, double t, double *d, double delh, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double dh;

    bb(t, mPro);
    ideal(t, mPro);
    do {
    base(t, *d, mPro);
    resid(t, *d, mPro);
    props(t, *d, mPro);
    dh = mPro->spro.h - h;
    if (fabs(dh) <= delh * fabs(h)) {
        format_pro(t, *d, mPro, prop);
        return;
    }
    derive(t, *d, mPro);
    *d -= dh / mPro->spro.dhd;
    if (*d < 0.0) {
        *d = 0.0000001;
    }
    i++;
    }  while (i < IMAX);
    prop->error = 1;
}

void ut(double u, double t, double *d, double delu, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double du;

    bb(t, mPro);
    ideal(t, mPro);
    do {
    base(t, *d, mPro);
    resid(t, *d, mPro);
    props(t, *d, mPro);
    du = mPro->spro.u - u;
    if (fabs(du) <= delu * fabs(u)) {
        format_pro(t, *d, mPro, prop);
        return;
    }
    derive(t, *d, mPro);
    *d -= du / mPro->spro.dud;
    if ((*d) < 0.0) {
        *d = 0.0000001;
    }
    i++;
    }  while (i < IMAX);
    prop->error = 1;
}

void st(double s, double t, double *d, double dels, S_mpro *mPro, Prop *prop)
{
    int i=0;
    double ds;
    const double dsmin = 1.0;

    if ((*d) <= 0.0) {
    *d = 1.0e-6;
    }
    bb(t, mPro);
    ideal(t, mPro);
    do {
    base(t, *d, mPro);
    resid(t, *d, mPro);
    props(t, *d, mPro);
    ds = mPro->spro.s - s;
    if (fabs(ds) <= dels * (fabs(s) + dsmin)) {
        format_pro(t, *d, mPro, prop);
        return;
    }
    derive(t, *d, mPro);
    ds = ds / mPro->spro.dsd;
    while ((*d - ds) < 0.0) {
        ds = ds * 0.5;
    }
    /* printf ("d = % .8e\t\tds = % .8e\n", *d, ds); */
    *d -= ds;
    i++;
    }  while (i < IMAX * 50);
    prop->error = 1;
}

/*
 */
/*---------------------------------------------------------------*/
/*  b)  x = (p, f, g, h, u, s) and d are known, t must be found  */
/*---------------------------------------------------------------*/
/* replace x by known parameter: */

/*
void xd(double x, double d, double *t, double delx, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double dx;

    do {
    calctd(*t, d, mPro);
    dx = mPro->spro.x - x;
    if (fabs(dx) <= delx * fabs(x)) {
        format_pro(*t, d, mPro, prop);
        return;
    }
        *//* call of deriv is not necessary for x=p *//*
    derive(*t, d, mPro);
        *t -= dx / mPro->spro.dxt;
        i++;
    }  while (i < IMAX);
    prop->error = 1;
}
*/

void pd(double p, double d, double *t, double delp, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double dp;

    do {
    calctd(*t, d, mPro);
    dp = mPro->spro.p - p;
    if (fabs(dp) <= delp * fabs(p)) {
        format_pro(*t, d, mPro, prop);
        return;
    }
    *t -= dp / mPro->spro.dpt;
    i++;
    }  while (i < IMAX);
    prop->error = 1;
}


void ud(double u, double d, double *t, double delu, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double du;
    const double dumin = 100.0;

    do {
    calctd(*t, d, mPro);
    du = mPro->spro.u - u;
    if (fabs(du) <= delu * (fabs(u) + dumin)) {
        format_pro(*t, d, mPro, prop);
        return;
    }
    derive(*t, d, mPro);
    *t -= du / mPro->spro.dut;
    i++;
    }  while (i < IMAX);
    prop->error = 1;
}

void sd(double s, double d, double *t, double dels, Prop *prop,
    S_mliq *mLiq, S_mpro *mPro)
/* Da calctd() nicht abfragt, ob ein Punkt im Nassdampfgebiet liegt,
   muss die Abfrage in dieser Funktion erfolgen. Die Berechnung von
   dht ist analog zu format_two(). Leider ist das damit gebildete
   Newton-Verfahren nicht stabil. Als Gegenmassnahme habe ich eine
   pseudozufaellige Daempfung eingefuehrt und gleichzeitig die Zahl der
   erlaubten Iterationen erhoeht. (Ole Engel)
*/
{
    int i = 0;
    double ds;
    int saturation;
    double p, dl, dv, x, sx;
    double cvl, cvv, cv;
    double dpt, ddpl, ddpv;
    double dst;
    double dt = 0.0;
    double rand1;
    const double dsmin = 1.0;
    /*  int show=0; */

    /*label: */
    do {
    /*    if (show) { */
    /*      fprintf (stderr, "*t = % .8e\n", *t); */
    /*    } */
    saturation = 0;
    if (((*t) >= tripl.t) && ((*t) <= crit.t)) {
        psat(*t, &p, &dl, &dv, mLiq, mPro);
        if ((d > dv) && (d < dl)) {
        saturation = 1;
        x = (1.0 / d - 1.0 / dl) / (1.0 / dv - 1.0 / dl);
        sx = mLiq->spro.s + x * (mPro->spro.s - mLiq->spro.s);
        ds = sx - s;
        if (fabs(ds) <= dels * (fabs(s) + dsmin)) {
            format_two(*t, p, x, dl, dv, mLiq, mPro, prop);
            return;
        }
        dpt = (mLiq->spro.s - mPro->spro.s) / (1.0 / dl - 1.0 / dv);
        ddpl = dpt - mLiq->spro.dpt;
        ddpv = dpt - mPro->spro.dpt;

        cvl = mLiq->spro.cv
            + (*t) * ddpl * ddpl / dl / dl / mLiq->spro.dpd;
        cvv = mPro->spro.cv
            + (*t) * ddpv * ddpv / dv / dv / mPro->spro.dpd;
        cv = cvl + x * (cvv - cvl);

        dst = cv / *t;
        dt  = ds / dst;
        }
    }
    if (saturation == 0) {
        calctd(*t, d, mPro);
        ds = mPro->spro.s - s;
        if (fabs(ds) <= dels * (fabs(s) + dsmin)) {
        format_pro(*t, d, mPro, prop);
        return;
        }
        derive(*t, d, mPro);
        dt = ds / mPro->spro.dst;
    }
    if (fabs(dt) < 1.0e-10 * (*t)) {
        /*      fprintf (stderr, "sd: return\n"); */
        format_pro(*t, d, mPro, prop);
        return;
    }
    /*    if (show) { */
    /*      fprintf (stderr, "dt = % .8e\n", dt); */
    /*    } */
    rand1 = (*t / 8.0) * rand() / (RAND_MAX + 1.0);
    if (dt > rand1) {
        *t -= rand1;
    } else if (dt < -rand1) {
        *t -= -rand1;
    } else {
        *t -= dt;
    }
    i++;
    }  while (i < 50 * IMAX);
    /*fprintf (stderr, "sd: maximum number of random search exceeded\n"); */
    /*show=1;  */
    /*i=0; */
    /*goto label; */
    prop->error = 1;
}

/*
 */
/*--------------------------------------------------------------------*/
/*  c)  x = (p, f, g, h, u, s) and y = (p, f, g, h, u, s) are known,  */
/*  t and d must be found                             */
/*--------------------------------------------------------------------*/

/* replace x and y by known parameters: */

/*
void xy(double x, double y, double *t, double *d,
    double delx, double dely, Prop *prop)
{
    int i = 0;
    double dx, dy, det, delt, deld;

    do {
    calctd(*t, *d);
    dx = mPro->spro.x - x;
    dy = mPro->spro.y - y;
    if ((fabs(dx) <= delx * fabs(x)) && (fabs(dy) <= dely * fabs(y))) {
        format_pro(*t, *d, S_mpro *mPro, prop);
        return;
    }
    derive(*t, *d);
    det = mPro->spro.dxt * mPro->spro.dyd - mPro->spro.dyt * mPro->spro.dxd;
    delt = (mPro->spro.dyd * dx - mPro->spro.dxd * dy) / det;
    deld = (mPro->spro.dxt * dy - mPro->spro.dyt * dx) / det;
    *t -= delt;
    *d -= deld;
    i++;
    } while(i < IMAX);
    prop->error = 1;
}
*/

void hs(double h, double s, double *t, double *d,
    double delh, double dels, S_mpro *mPro, Prop *prop)
/*
  computes thermodynamic functions for given enthalpy and entropy in the one-
  phase-region within given tolerances of enthalpy and entropy using given
  approximations of temperature and density (which are changed to yield h, s)
*/
{
    int i = 0;
    double dh, ds, det, delt, deld;
    const double dhmin = 100.0;
    const double dsmin = 1.0;

    do {
    adjust_hsp(t, d);
    calctd(*t, *d, mPro);
    dh = mPro->spro.h - h;
    ds = mPro->spro.s - s;
    if ((fabs(dh) <= delh * (fabs(h) + dhmin))
        && (fabs(ds) <= dels * (fabs(s) + dsmin))) {
        format_pro(*t, *d, mPro, prop);
        return;
    }
    derive(*t, *d, mPro);
    det = mPro->spro.dht * mPro->spro.dsd - mPro->spro.dst * mPro->spro.dhd;
    delt = (mPro->spro.dsd * dh - mPro->spro.dhd * ds) / det;
    deld = (mPro->spro.dht * ds - mPro->spro.dst * dh) / det;
    *t -= delt;
    *d -= deld;
    i++;
    } while (i < IMAX);
    prop->error = 1;
}


/* you may take the following as an input control for system (h, s)*/

int valid_hs(double h, double s)
/*
   Returns 1 if (h, s) is valid.
   A clumsy control of validity of (h, s) with a vague guarantee, that
   outside two-phase-region only one (t, d) can be found to fit (h, s).
   This region of validity includes the one specified by valid_tp
   except for the ice-near-region above p = 100 MPa.
*/
{
    int val = 1;
    double hv = 2500.5374565175789, hl = 0.000640282639798101265;
    double sv = 9.15410669412127653, sl = -0.00000377679288979491875;

    if ((h < -10.0) || (h > 9550.0) || (s > 19.5) || (s < 0.0005 * h - 0.275)) {
    val = 0;
    } else if (h < hv) {
    if (s > sv) {
        val = (h > 2460.0);
    } else {
        val = (h > hl + (s - sl) * (hv - hl) / (sv - sl));  /* p = tripl.p*/
    }
    }
    return val;
}

void pu(double p, double u, double *t, double *d,
    double delp, double delu, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double dp, du, det, delt, deld;
    const double dumin = 100.0;

    do {
    calctd(*t, *d, mPro);
    dp = mPro->spro.p - p;
    du = mPro->spro.u - u;
    if ((fabs(dp) <= delp * fabs(p))
        && (fabs(du) <= delu * (fabs(u) + dumin))) {
        format_pro(*t, *d, mPro, prop);
        return;
    }
    derive(*t, *d, mPro);
    det = mPro->spro.dpt * mPro->spro.dud - mPro->spro.dut * mPro->spro.dpd;
    delt = (mPro->spro.dud * dp - mPro->spro.dpd * du) / det;
    deld = (mPro->spro.dpt * du - mPro->spro.dut * dp) / det;
    *t -= delt;
    *d -= deld;
    i++;
    } while (i < IMAX);
    prop->error = 1;
}

void us(double u, double s, double *t, double *d,
    double delu, double dels, S_mpro *mPro, Prop *prop)
{
    int i = 0;
    double du, ds, det, delt, deld;
    const double dumin = 100.0;
    const double dsmin = 1.0;

    do {
    calctd(*t, *d, mPro);
    du = mPro->spro.u - u;
    ds = mPro->spro.s - s;
    if ((fabs(du) <= delu * (fabs(u) + dumin))
        && (fabs(ds) <= dels * (fabs(s) + dsmin))) {
        format_pro(*t, *d, mPro, prop);
        return;
    }
    derive(*t, *d, mPro);
    det = mPro->spro.dut * mPro->spro.dsd - mPro->spro.dst * mPro->spro.dud;
    delt = (mPro->spro.dsd * du - mPro->spro.dud * ds) / det;
    deld = (mPro->spro.dut * ds - mPro->spro.dst * du) / det;
    *t -= delt;
    *d -= deld;
    i++;
    } while (i < IMAX);
    prop->error = 1;
}

void water_dx0(double d, double deld, Prop *prop)
{
    int i=0;
    double t2 = tripl.t;
    double d2 = tripl.dl;
    double t1 = crit.t;
    double d1 = crit.d;
    double tm, dd;
    double dv, dl, p;
    const double x = 0.0;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    d *= 1.0e-3;

    prop->error=1;
    if ((d < d1) || (d > d2)) {
    return;
    }
    do {
    tm = (t1 + t2) * 0.5;
    psat(tm, &p, &dl, &dv, &MLiq, &MPro);

    dd = dl - d;
    if ((fabs(dd) <= deld * fabs(d))
        || (fabs((t2 - t1) / t2) <= DBL_EPSILON)) {
        format_two(tm, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }

    if (d < dl) {
        t2 = tm;
    } else {
        t1 = tm;
    }

    i++;
    } while (i < IMAX * 50);
}

void water_dx1(double d, double deld, Prop *prop)
{
    int i = 0;
    double t1 = tripl.t;
    double d1 = tripl.dv;
    double t2 = crit.t;
    double d2 = crit.d;
    double tm, dd;
    double dv, dl, p;
    const double x = 1.0;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    d *= 1.0e-3;

    prop->error = 1;
    if ((d < d1) || (d > d2)) {
    return;
    }
    do {
    tm = (t1 + t2) * 0.5;
    psat(tm, &p, &dl, &dv, &MLiq, &MPro);

    dd = dv - d;
    if ((fabs(dd) <= deld * fabs(d))
        || (fabs((t2 - t1) / t2) <= DBL_EPSILON)) {
        format_two(tm, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }

    if (d < dv) {
        t2 = tm;
    } else {
        t1 = tm;
    }

    i++;
    } while (i < IMAX * 50);
}

void water_dxm(double d, double x, double deld, Prop *prop)
{
    int i = 0;
    double t1 = tripl.t;
    double d1;
    double t2 = crit.t;
    double d2 = crit.d;
    double tm, dd;
    double dx, dv, dl, p;
    double vv, vl, vx;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    d *= 1.0e-3;

    prop->error = 1;
    if (x < 0.0 || x > 1.0) {
    return;
    }
    vv = 1.0 / tripl.dv;
    vl = 1.0 / tripl.dl;
    vx = vl + x * (vv - vl);
    d1 = 1.0 / vx;
    if ((d < d1) || (d > d2)) {
    return;
    }
    do {
    tm = (t1 + t2) * 0.5;
    psat(tm, &p, &dl, &dv, &MLiq, &MPro);
    vv = 1.0 / dv;
    vl = 1.0 / dl;
    vx = vl + x * (vv - vl);
    dx = 1.0 / vx;

    dd = dv - dx;
    if ((fabs(dd) <= deld * fabs(d))
        || (fabs((t2 - t1) / t2) <= DBL_EPSILON)) {
        format_two(tm, p, x, dl, dv, &MLiq, &MPro, prop);
        return;
    }

    if (d < dx) {
        t2 = tm;
    } else {
        t1 = tm;
    }

    i++;
    } while (i < IMAX * 50);
}

void water_dx(double d, double x, double t,
          double deld, Prop *prop)
{
    /* check pointer */
    if (prop == NULL) { return; }

    if (t) {} /* avoid 'unused parameter' warning */

    if (x == 0.0) {
    water_dx0(d, deld, prop);
    } else if (x == 1.0) {
    water_dx1(d, deld, prop);
    } else {
    water_dxm(d, x, deld, prop);
    }
}

void water_px(double p, double x, Prop *prop)
{
    double dl, dv, t;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    prop->error = 1;
    if ((p >= tripl.p) && (p <= crit.p)) {
    tsat(p, &t, &dl, &dv, &MLiq, &MPro);
    format_two(t, p, x, dl, dv, &MLiq, &MPro, prop);
    }
}

void water_tx(double t, double x, Prop *prop)
{
    double dl, dv, p;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    prop->error = 1;
    /* Damit kein Bereichsfehler auftritt, wenn t minimal kleiner als
       der Tripelpunkt ist, wird nach (tripl.t - 0.01) gefragt. */
    if ((t >= (tripl.t - 0.01)) && (t <= crit.t)) {
    psat(t, &p, &dl, &dv, &MLiq, &MPro);
    format_two(t, p, x, dl, dv, &MLiq, &MPro, prop);
    }
}
