// SPDX-FileCopyrightText: © 1998 Olaf Bauer <o.bauer@tu-harburg.de>
// SPDX-FileCopyrightText: © 1998 Ole Engel <engel@tu-harburg.de>
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: iaps.h,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Olaf Bauer
# maintainer: Ole Engel <engel@tu-harburg.de>
#---------------------------------------------------------------------
# IAPS formulation 1984
# main part of PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/


#ifndef _IAPS_H_
#define _IAPS_H_

#include "steam4.h"
#include "decl.h"

#define IMAX 20        /* maximum number of iterations */
#define P_MIN 0.000001  /* minimum pressure */


/**********************/
/*                    */
/* internal functions */
/*                    */
/**********************/

/*---------------- calculation of one-phase-properties -------------*/
/* internal functions doesn't check the pointers, they must no be NULL */

void td(double t, double d, S_mpro *mPro, Prop *prop);
/*
  Computes thermodynamic functions for given temperature and density
  in the one-phase-region
*/

void tp(double t, double p, double *d, double delp, S_mpro *mPro, Prop *prop);
/*
  Computes thermodynamic functions for given temperature and pressure
  in the one-phase-region within a given tolerance of pressure using a
  given approximation of density.
  For t < 647.126K region of state must be known and *d chosen as follows:
    p < ps(t) => vapour => (*d) < 0.32189 g/cm3
    p > ps(t) => liquid => (*d) > 0.32189 g/cm3
*/

void ph(double p, double h, double *t, double *d,
    double delp, double delh, S_mpro *mPro, Prop *prop);
/*
  Computes thermodynamic functions for given enthalpy and pressure in the one-
  phase-region within given tolerances of enthalpy and pressure using given
  approximations of temperature and density (which are changed to yield h,p)
*/

void ps(double p, double s, double *t, double *d,
    double delp, double dels, S_mpro *mPro, Prop *prop);
/*
  Computes thermodynamic functions for given entropy and pressure in the one-
  phase-region within given tolerances of entropy and pressure using given
  approximations of temperature and density (which are changed to yield s,p)
*/

void hd(double h, double d, double *t, double delh, S_mpro *mPro, Prop *prop);
/*
  Computes thermodynamic functions for given enthalpy and density
  in the one-phase-region within a given tolerance of enthalpy using a
  given approximation of temperature.
*/


/* ----------------- extrapolation function(s) -------------------- */

void extra_ph(double p, double h, double *t, double *d,
          double dp, double dh, S_mpro *mPro, Prop *prop);
/*
   For invalid (p, h) extrapolates towards (p, h) from a point on
   validity border (specified in valid_ph) closest to (p, h).
   In some areas ph is used to calculate border-point. If interative
   process in ph is not successful, no extrapolation takes place.
*/


/* --------------- Calculation of saturation properties -----------*/

void psat(double t, double *p, double *dl, double *dv,
      S_mliq *mLiq, S_mpro *mPro);
/*
  Computes liquid and vapour properties for a given saturation temperature
*/

void psatc(double t, double *p, double *dl, double *dv,
       S_mliq *mLiq, S_mpro *mPro);
/*
  Saturation properties for 646.303775 K < t <= 647.126 K
*/

void tsat(double p, double *t, double *dl, double *dv,
      S_mliq *mLiq, S_mpro *mPro);
/*
  Computes liquid and vapour properties for a given saturation pressure
*/

void tsatc(double p, double *t, double *dl, double *dv,
       S_mliq *mLiq, S_mpro *mPro);
/*
  Saturation properties for 21.839129 MPa <= p <= 22.054915
*/

void sat(double t, double dl, double dv, double *delg,
     S_mliq *mLiq, S_mpro *mPro);
/*
  Computes supposed liquid (sliq.x) and vapour (spro.x) properties
  and offers difference of gibbs-function delg = |(gl - gv) / (R T)|
*/

void hdsat(double h, double d, double delh, double *t,
       double *p,double *dl, double *dv, double *x,
       S_mliq *mLiq, S_mpro *mPro);
/*
  find saturation properties for (h, d) in two-phase-region;
  in case of success, enthalpy will be within delh * (h + 100);
  works for tripl.dv <= d <= 1.0 and tripl.t <= t <= creg.t only;
*/


void hdsatc(double h, double d, double delh, double *t,
        double *p, double *dl, double *dv, double *x,
        S_mliq *mLiq, S_mpro *mPro);
/*
  find saturation properties for (h, d) in critical two-phase-region;
  in case of success, enthalpy will be within delh * (h + 100);
  works for 0.2439 <= d <= 0.4 and creg.t <= t < crit.t only;
*/


void save(S_mliq *mLiq, S_mpro *mPro);
/*
  Save liquid properties and some internal results
  for later calculation of liquid third derivatives
*/

void load(S_mliq *mLiq, S_mpro *mPro);
/*
  Reload liquid properties and some internal results
*/


/*-------------- approximations for saturation properties ----------*/

double approx_ps(double t, double *dpsdt);
/*
  Approximation of saturation pressure and its derivation by temperature.
  Equation from program code in "NBS/NRC Steam-Tables".
  Below 646.3 K error of ps is within 0.01 5% compared to solutions of psat
*/

double approx_ts(double p);
/*
  Approximation of saturation temperature for given pressure
*/

void approx_dlv(double t ,double *dl, double *dv);
/*
  Approximations of liquid and vapour density for given temperature
*/

void approx_hlvp(double p,double *hl, double *hv);
/*
  Approximations for saturation enthalpies for tripl.p <= p <= creg.p.
  Error is within (+-)0.08 kJ/kg compared to solutions of psat.
  For p > 7 MPa error of hv is within (+-)0.03 kJ/kg.
*/

double approx_thd(double h, double d);
/*
  Offers approximation of t(h, d) in two-phase-region;
*/


/* ------------------------ input control ---------------------- */

int valid_td(double t, double d);
/*
   returns 1 if (t, p) falls within IAPS-bounds (fluid region of water)
   (upper bound of d is a bit exceeded - for d > 1.4 water is no longer fluid)
*/


int valid_tp(double t, double p);
/*
  returns 1 if (t, p) falls within IAPS-bounds (fluid region of water)
*/

double psublm(double t);
/*
  sublimation pressure for 190 K <= T <= 273.16 K
*/

double pice1(double t);
/*
  melting pressure for 251.165 K <= T <= 273.16 K and p < 209.9 MPa
*/

double pice(double t);
/*
  melting pressure for 251.165 K <= T <= 715 K and p > 209.9 MPa
*/

int region_tp(double t, double p, double *dl, double *dv,
          S_mliq *mLiq, S_mpro *mPro);
/*
  For t < crit.t speed-optimized control of region of (t, p)
*/

int valid_ph(double p, double h);
/*
   Returns 1 if (p, h) is valid.
*/

int region_ph(double p, double h, double *ts, double *dl, double *dv,
          S_mliq *mLiq, S_mpro *mPro);
/*
  checks region of (p, h)
*/

int valid_ps(double p, double s);
/*
   Returns 1 if (p, s) is valid.
*/

int valid_hd(double h, double d);
/*
   Returns 1 if (h, d) is valid.
*/

int region_hd(double h, double d, double dh, double *t,
          double *p, double *dl, double *dv, double *x,
          S_mliq *mLiq, S_mpro *mPro);
/*
  checks region of (h, d)
*/


/*---------------- supervision of iterations -----------------*/

void adjust_tp(double t, double d, double *dmin, double *dmax);
/*
   Offers maximum and minimum density to adjust d
   during iteration to avoid multiple solutions of d(p, t).
   For t <= crit.t region of state must be known
   (d < 0.32189 means vapour state, liquid otherwise)
*/

void adjust_hsp(double *t, double *d);
/*
   for systems (h, p) and (s, p) adjusts t and d during iteration
   to avoid multiple solutions
*/

void adjust_hd(double d, double *tmin, double *tmax);
/*
  offers minimum and maximum values of temperature
  for adjustment during iteration in hd to avoid multiple solutions
*/


/* -------------- Formation of output structure  --------------- */


void format_pro (double t, double d, S_mpro *mPro, Prop *prop);
/*
  Formation of output structure prop for one-phase properties
*/

void deriv_ph(double t, double d, S_mpro *mPro, Prop *prop);
/*
  Calculate first and, if required, second derivatives of T,d,s,u
  by p(h = const) and h(p = const) from a formatted pro-struct
  in the one phase region.
*/

void deriv_ps(double t, double d, S_mpro *mPro, Prop *prop);
/*
  Calculate first derivatives of T, d, h, u by p(s = const) and
  s(p = const) from a formatted pro-struct in the one phase region.
*/


void format_two(double t, double p, double x, double dl, double dv,
        S_mliq *mLiq, S_mpro *mPro, Prop *prop);
/*
  Formation of output structure prop for two-phase properties
*/

void deriv_ph2(Prop *prop);
/*
  Calculates first and, if required, second derivatives of x, T, d, s, u
  by p(h = const) and h(p = const) from a formatted prop-struct
  in the two phase region.
*/

void deriv_ps2(Prop *prop);
/*
  Calculates first derivatives of x, T, d, h, u by p(s = const) and
  s(p = const) from a formatted prop-struct in the two phase region.
*/


/* -------------------------- IAPS-Formulation ---------------------*/

void calctd(double t, double d, S_mpro *mPro);
/*
  Calculates properties and places them in internal structure pro
*/

void props(double t, double d, S_mpro *mPro);
/*
   properties, put together from base, residual and ideal gas function
*/

void derive(double t, double d, S_mpro *mPro);
/*
  derivatives of properties, useful for numerical methods
*/

void bb(double t, S_mpro *mPro);
/*
  molecular parameters
*/

void base(double t, double d, S_mpro *mPro);
/*
  Base function and its first and second derivatives
*/

void resid(double t, double d, S_mpro *mPro);
/*
  Residual function and its first and second derivatives
*/

void ideal(double t, S_mpro *mPro);
/*
  Ideal gas function and its first and second derivatives
*/

void third(double t, double d, S_mpro *mPro, S_pro *out);
/*
  Third derivatives of the entire fundamental equation.
*/


#endif /* _IAPS_H_ */
