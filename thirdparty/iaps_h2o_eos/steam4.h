// SPDX-FileCopyrightText: © 1998 Ole Engel <engel@tu-harburg.de>
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: steam4.h,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Ole Engel
# maintainer: Ole Engel <engel@tu-harburg.de>
#---------------------------------------------------------------------
# All functions needed for using PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/

/*
 * Modified by "The openCSMP project", 2026:
 * - Replaced headers (stdio.h) with (cstdio)
 * - Modified lines 109:
 *   extern void dumpProp(FILE *fp, Prop *prop); -> extern void dumpProp(std::FILE *fp, Prop *prop);
 */

/*********************************************************************
          see file PROST4.TXT for help
**********************************************************************/


#ifndef _STEAM4_H_
#define _STEAM4_H_

#include <cstdio> /* FILE */

#define ONE 1  /* one phase */
#define TWO 2  /* two phase */

#define LIQUID    3
#define SATURATED 4
#define VAPOUR    5

/* properties at the critical point */
#define TCRIT (crit.t)          /* 647.126 K */
#define PCRIT (crit.p * 1.0e6)  /* 22.054915e6 Pa */
#define DCRIT (crit.d * 1.0e3)  /* 0.32189e3 kg/m^3 */


typedef struct deriv_prop_struct
{
    double T_Cd;
    double d_CT;
    double h_Cp;
    double p_Ch;
    double p_Cs;
    double s_Cp;
    struct deriv_prop_struct *dT_Cd;
    struct deriv_prop_struct *dd_CT;
    struct deriv_prop_struct *dh_Cp;
    struct deriv_prop_struct *dp_Ch;
    struct deriv_prop_struct *dp_Cs;
    struct deriv_prop_struct *ds_Cp;
} d_Prop;

typedef struct prop_struct
{
    double x;
    double T;
    double d;
    double p;
    double f;
    double g;
    double s;
    double u;
    double h;
    double cv;
    double cp;
    d_Prop *dx;
    d_Prop *dT;
    d_Prop *dd;
    d_Prop *dp;
    d_Prop *df;
    d_Prop *dg;
    d_Prop *ds;
    d_Prop *du;
    d_Prop *dh;
    d_Prop *dcv;
    d_Prop *dcp;
    int phase;
    int error;
    char indep1;
    char indep2;
    int deriv;
} Prop;


/***********************************/
/*                                 */
/* external functions from steam.c */
/*                                 */
/***********************************/

extern Prop *newProp(int indep1, int indep2, int deriv);
/*
  Allocate memory for a Prop structure. Call newProp() with characters
  as argument for indep1 and indep2.
  For example:
      Prop *myProp = newProp('p', 'h', 2);
*/

extern Prop *freeProp(Prop *prop);
/*
  Free all the memory allocated for a Prop structure. If you just call
  free(), you will keep allocated memory.
*/

extern void dumpProp(std::FILE *fp, Prop *prop);
/*
  Print all values in the structure prop to file fp.
*/


/**********************************/
/*                                */
/* external functions from iaps.c */
/*                                */
/**********************************/
/****************************************
   normal functions:
   offer output value for prop->phase
*****************************************/

extern void sat_t(double t, Prop *pliq, Prop *pvap);
/*
  Computes liquid and vapour properties for a given saturation temperature
*/

extern void sat_p(double p, Prop *pliq, Prop *pvap);
/*
  Computes liquid and vapour properties for a given saturation pressure
*/

extern void water_td(double t, double d, Prop *prop);
/*
  Computes properties for given temperature and density
*/

extern void water_tp(double t, double p, double d, double dp, Prop *prop);
/*
  Computes properties for given temperature and pressure in the one-
  phase-region within a given tolerance of pressure using a given
  approximation of density.
  No output in case of saturation, i.e. |1 - p / ps(T)| < 1e-6.
*/

extern void water_ph(double p, double h, double t, double d,
             double dp, double dh, Prop *prop);
/*
  Computes properties for given pressure and enthalpy within given
  tolerances of pressure and enthalpy using given approximations
  of temperature and density.
  For invalid (p, h) it extrapolates from bounds of validity.
*/

extern void water_ps(double p, double s, double t, double d,
             double dp, double ds, Prop *prop);
/*
  Computes properties for given pressure and entropy within given
  tolerances of pressure and entropy using given approximations
  of temperature and density.
*/

extern void water_hd(double h, double d, double t, double dh, Prop *prop);
/*
  Computes properties for given enthalpy h and density d within a given
  relative tolerance of enthalpy dh using a given approximation of
  temperature t.
*/


/****************************************
   meta-functions:
   require input value for prop->phase
*****************************************/

extern double speed(Prop *prop);
/*
  Speed of sound for a calculated prop structure, only in one-phase-region.
*/

extern double wkappa(Prop *prop);
/*
  Isentropic expansion coefficient for a calculated prop structure,
  only in one-phase-region.
*/

extern double wbetas(Prop *prop);
/*
  Isentropic pressure-temperature coefficient for a calculated
  prop-structure, only in one-phase-region.
*/

extern double viscos(Prop *prop);
/*
  Viscosity [Pas] for a calculated prop-structure
  Returns zero if region of validity is exceeded
  (specified by the initial if statement).
  Results do not agree entirely with the NBS/NRC values,
  because an improved formula is used.
  Code from "WaterC" by M.L. Demsey, University of Arizona.
*/

extern double thcond(Prop *prop);
/*
  Thermal conductivity [W/mK] for a calculated prop-structure
  Returns zero if region of validity is exceeded
  (specified by the initial if statement).
  Code from "WaterC" by M.L. Demsey, University of Arizona.
*/


/**********************************/
/*                                */
/* external functions from meta.c */
/*                                */
/**********************************/

/*
  Modified versions of water_td, water_ph, water_ps and water_hd called
  meta_td, meta_ph, meta_ps and meta_hd. They produce a fictitious
  steadiness of the Helmholtz-function's derivatives on the saturation line.

  Functions require a logical argument prop->phase.
  Point will be treated as if it belongs to prop->phase.

  For phase != TWO functions directly use the fundamental equation,
  which in the two-phase region produces some kind of meta stable
  properties. If they fail (error = 1), no definit solution can be found.

  For phase == TWO functions simply ignore the fact that the
  steam quota can only be 0 <= x <= 1. If they fail (error = 1),
  no saturation-like properties can be offered.

  For invalid (p, h)-points meta_ph, just like water_ph,
  extrapolates from bounds of validity towards (p, h).
*/

void meta_td(double t, double d, Prop *prop);
/*
   Supposes that region of (p, h) is same as given prop->phase.
   prop->phase can be ONE or TWO.
   For phase == ONE the fundamental equation is used directly.
   In the two-phase-region this produces some kind of meta stable properties.
   For phase == TWO properties are calculated via saturation properties of t,
   using a steam quota x < 0 or x > 1, when (t, s) "really" is a one-phase
   point.

   This works only for tripl.t <= t <= crit.t (otherwise: error = 1).
*/

void meta_ph(double p, double h, double t, double d,
         double dp, double dh, Prop *prop);
/*
   Supposes that region of (p, h) is same as given prop->phase.
   prop->phase can be ONE or TWO.
   prop->phase does not have to agree with phase of (t, d),
   (t, d) are used only as starting values.

   For phase == ONE the fundamental equation is used directly.
   In the two-phase-region this produces some kind of meta stable properties.
   You cannot step too deeply inside it this way, because there are multiple
   solutions for t(p, h) and d(p, h) in it. The Newton-Iteration is modified in
   a way that these multiple solutions cannot be found (result: error = 1).
   However there should be no difficulties for
   h < h'(p) + 150 kJ/kg  or  h > h''(p) - 80 kJ/kg.
   For p > 16.5 MPa all enthalpies allowed.
   For h < 1500 kJ/kg all pressures allowed.

   For phase == TWO properties are calculated via saturation properties of p,
   using a steam quota x < 0 or x > 1, when (p, h) "really" is a one-phase
   point.

   This works only for tripl.p <= p <= crit.p (otherwise: error = 1).

   For invalid (p, h) extrapolates from bounds of validity and sets error = 1.
*/


void meta_ps(double p, double s, double t, double d,
         double dp, double ds, Prop *prop);
/*
   Supposes that region of (p, s) is same as given prop->phase.
   prop->phase can be ONE or TWO.
   prop->phase does not have to agree with phase of (t, d),
   (t, d) are used only as starting values.

   For phase == ONE the fundamental equation is used directly.
   In the two-phase-region this produces some kind of meta stable properties.
   You cannot step too deeply inside it this way, because there are multiple
   solutions for t(p, s) and d(p, s) in it. The Newton-Iteration is modified in
   a way that these multiple solutions cannot be found (result: error = 1).

   For phase == TWO properties are calculated via saturation properties of p,
   using a steam quota x < 0 or x > 1, when (p, s) "really" is a one-phase
   point.

   This works only for tripl.p <= p*1e-6 <= crit.p (otherwise: error = 1).
*/

void meta_hd(double h, double d, double t, double dh, Prop *prop);
/*
   Supposes that region of (h, d) is same as given prop->phase.
   prop->phase can be ONE or TWO.
   t is used only as starting value.

   For phase == ONE the fundamental equation is used directly.
   In the two-phase-region this produces some kind of meta stable properties.
   You cannot step too deeply inside it this way, because there are multiple
   solutions for t(h, d) in it. The Newton-Iteration is modified in
   a way that these multiple solutions cannot be found (result: error = 1)

   For phase == TWO properties are calculated via saturation properties of
   t fitting (h, d) using a steam quota x < 0 or x > 1, when (h, d) "really"
   is a one-phase point. This works only for:

   0.00485467589830342234 <= d <= 1000.0
   and
   -11.5017594+12139.3548/d <= h <= hmax
   with
   hmax = 1547235.7851199+173409885.1329 / d for 243.9 <= d <= 400.0
   and
   hmax = 1547745.4041370+169324991.2165 / d otherwise
*/


/**********************************/
/*                                */
/* external functions from more.c */
/*                                */
/**********************************/
/****************************************
   normal functions:
   offer output value for prop->phase
*****************************************/

extern void water_hs(double h, double s, double t, double d,
             double dh, double ds, Prop *prop);

extern void water_pd(double p, double d, double t,
             double dp, Prop *prop);

extern void water_pu(double p, double u, double t, double d,
             double dp, double du, Prop *prop);

extern void water_px(double p, double x, Prop *prop);

extern void water_sd(double s, double d, double t,
             double ds, Prop *prop);

extern void water_th(double t, double h, double d,
             double dh, Prop *prop);

extern void water_ts(double t, double s, double d,
             double ds, Prop *prop);

extern void water_tu(double t, double u, double d,
             double du, Prop *prop);

extern void water_tx(double t, double x, Prop *prop);

extern void water_ud(double u, double d, double t,
             double du, Prop *prop);

extern void water_us(double u, double s, double t, double d,
             double du, double ds, Prop *prop);

extern void water_dx(double d, double x, double t,
             double deld, Prop *prop);

extern void water_dx0(double d,
              double deld, Prop *prop);

extern void water_dx1(double d,
              double deld, Prop *prop);

extern void water_dxm(double d, double x,
              double deld, Prop *prop);


/****************************************
   normal functions:
   offer output value for prop->phase
*****************************************/

extern double wthcond(double t, double d);
/* thermal conductivity as in thcond, but checked for two-phase regions.
   if IsInTwoPhaseRegion returns true, the values from the phase regions
   are linearly interpolated.
   Not safe around the critical point!!!
   !!! TEST-VERSION !!!
*/

extern double wviscos(double t, double d);
/* thermal conductivity as in thcond, but checked for two-phase regions.
   if IsInTwoPhaseRegion returns true, the values from the phase regions
   are linearly interpolated.
   Not safe around the critical point!!!
   !!! TEST-VERSION !!!
*/

#endif /* _STEAM4_H_ */
