// SPDX-FileCopyrightText: © 1998 Ole Engel <engel@tu-harburg.de>
//
// SPDX-License-Identifier: LGPL-2.0-only
/*------------------------------------------------------------*- C -*-
# $Id: steam.cpp,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Ole Engel
#             Joerg Ungethuem
# maintainer: Ole Engel <engel@tu-harburg.de>
#---------------------------------------------------------------------
# All functions needed for using PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/

/*
 * Modified by "The openCSMP project", 2026:
 * - Ported from meta.c to meta.cpp for C++ compilation
 * - Replaced C headers (stdlib.h) with C++ equivalents (cstdlib)
 * - Use of namespace std
 * - Modified lines 73, 74:
 *   prop->indep1 = indep1; -> prop->indep1 = (char)indep1;
 *   prop->indep2 = indep2; -> prop->indep2 = (char)indep2;
 * - Modified line 198:
 *   " t = %g °C\t" -> " t = %g C\t"
 */


#include <cstdlib> /* NULL */
#include "steam4.h"

using namespace std;

/* local functions */
static d_Prop *new_dProp(void);
static d_Prop *free_dProp(d_Prop *d_prop);


/* ---- functions to create, delete and dump the Prop structure ---- */

static d_Prop *new_dProp(void)
{
    d_Prop *d_prop;

    /* allocate memory for the new data struct */
    if (! (d_prop = (d_Prop*) malloc(sizeof(d_Prop)))) { return NULL; }

    /* init the real variables to zero */
    d_prop->T_Cd = d_prop->d_CT = 0.0;
    d_prop->h_Cp = d_prop->p_Ch = 0.0;
    d_prop->p_Cs = d_prop->s_Cp = 0.0;

    /* init the derivation struct pointers to NULL */
    d_prop->dT_Cd = d_prop->dd_CT = NULL;
    d_prop->dh_Cp = d_prop->dp_Ch = NULL;
    d_prop->dp_Cs = d_prop->ds_Cp = NULL;

    return d_prop;
}

Prop *newProp(int indep1, int indep2, int deriv)
/*
  Allocate memory for a Prop structure. Call newProp() with characters
  as argument for indep1 and indep2.
  For example:
      Prop *myProp = newProp('p', 'h', 2);
*/
{
    int error = 0;
    Prop *prop = NULL;

    /* allocate memory for the new data struct */
    if (! (prop = (Prop*) malloc(sizeof(Prop)))) {
    return NULL;
    }

    prop->indep1 = (char)indep1;
    prop->indep2 = (char)indep2;
    prop->deriv = deriv;

    /* check depth of derivation to calculate */
    if ((deriv < 0) || (deriv > 2)) {
    fprintf(stderr, "newProp(): deriv out of range\n");
    exit(1);
    }

    /* init the real variables to zero */
    prop->x = prop->T = prop->d = prop->p = prop->f = prop->g = 0.0;
    prop->s = prop->u = prop->h = prop->cv = prop->cp = 0.0;

    /* init the derivation struct pointers to NULL */
    prop->dx = prop->dT = prop->dd = prop->dp = prop->df = prop->dg = NULL;
    prop->ds = prop->du = prop->dh = prop->dcv = prop->dcp = NULL;

    /* if we need derivations allocate the memory */
    if (deriv >= 1) {
    if (! (prop->dx  = new_dProp())) { error++; }
    if (! (prop->dp  = new_dProp())) { error++; }
    }
    if (deriv == 2) {
    if (! (prop->dcv = new_dProp())) { error++; }
    if (! (prop->dp->dT_Cd  = new_dProp())) { error++; }
    if (! (prop->dp->dd_CT  = new_dProp())) { error++; }
    if (! (prop->dx->dT_Cd  = new_dProp())) { error++; }
    if (! (prop->dx->dd_CT  = new_dProp())) { error++; }
    }

    /* allocate memory for some special derivations */
    if ((indep1 == 'p') && (indep2 == 'h')) {
    if (deriv >= 1) {
        if (! (prop->dT = new_dProp())) { error++; }
        if (! (prop->dd = new_dProp())) { error++; }
        if (! (prop->ds = new_dProp())) { error++; }
        if (! (prop->du = new_dProp())) { error++; }
    }
    if (deriv == 2) {
        if (! (prop->dcp = new_dProp())) { error++; }
        if (! (prop->dx->dp_Ch  = new_dProp())) { error++; }
        if (! (prop->dT->dp_Ch  = new_dProp())) { error++; }
        if (! (prop->dd->dp_Ch  = new_dProp())) { error++; }
        if (! (prop->ds->dp_Ch  = new_dProp())) { error++; }
        if (! (prop->du->dp_Ch  = new_dProp())) { error++; }
        if (! (prop->dx->dh_Cp  = new_dProp())) { error++; }
        if (! (prop->dT->dh_Cp  = new_dProp())) { error++; }
        if (! (prop->dd->dh_Cp  = new_dProp())) { error++; }
        if (! (prop->ds->dh_Cp  = new_dProp())) { error++; }
        if (! (prop->du->dh_Cp  = new_dProp())) { error++; }
    }
    } else if ((indep1 == 'p') && (indep2 == 's')) {
    if (deriv >= 1) {
        if (! (prop->dT = new_dProp())) { error++; }
        if (! (prop->dd = new_dProp())) { error++; }
        if (! (prop->dh = new_dProp())) { error++; }
        if (! (prop->du = new_dProp())) { error++; }
    }
    if (deriv == 2) {
        if (! (prop->dcp = new_dProp())) { error++; }
    }
    }

    /* if any error occured, free all the memory and return NULL  */
    if (error != 0) {
    prop = freeProp(prop);
    }
    return prop;
}

d_Prop *free_dProp(d_Prop *d_prop)
{
    /* if the memory is already freed just return NULL */
    if (d_prop == NULL) { return NULL; }

    /* now we free all the memory allocated for derivations */
    if (d_prop->dT_Cd != NULL) { d_prop->dT_Cd = free_dProp(d_prop->dT_Cd); }
    if (d_prop->dd_CT != NULL) { d_prop->dd_CT = free_dProp(d_prop->dd_CT); }
    if (d_prop->dh_Cp != NULL) { d_prop->dh_Cp = free_dProp(d_prop->dh_Cp); }
    if (d_prop->dp_Ch != NULL) { d_prop->dp_Ch = free_dProp(d_prop->dp_Ch); }
    if (d_prop->dp_Cs != NULL) { d_prop->dp_Cs = free_dProp(d_prop->dp_Cs); }
    if (d_prop->ds_Cp != NULL) { d_prop->ds_Cp = free_dProp(d_prop->ds_Cp); }

    free(d_prop);
    return NULL;
}

Prop *freeProp(Prop *prop)
/*
  Free all the memory allocated for a Prop structure. If you just call
  free(), you will keep allocated memory.
*/
{
    /* if the memory is already freed just return NULL */
    if (prop == NULL) { return NULL; }

    /* now we free all the memory allocated for derivations */
    if (prop->dx  != NULL) { prop->dx  = free_dProp(prop->dx); }
    if (prop->dT  != NULL) { prop->dT  = free_dProp(prop->dT); }
    if (prop->dd  != NULL) { prop->dd  = free_dProp(prop->dd); }
    if (prop->dp  != NULL) { prop->dp  = free_dProp(prop->dp); }
    if (prop->df  != NULL) { prop->df  = free_dProp(prop->df); }
    if (prop->dg  != NULL) { prop->dg  = free_dProp(prop->dg); }
    if (prop->ds  != NULL) { prop->ds  = free_dProp(prop->ds); }
    if (prop->du  != NULL) { prop->du  = free_dProp(prop->du); }
    if (prop->dh  != NULL) { prop->dh  = free_dProp(prop->dh); }
    if (prop->dcv != NULL) { prop->dcv = free_dProp(prop->dcv); }
    if (prop->dcp != NULL) { prop->dcp = free_dProp(prop->dcp); }

    free(prop);
    return NULL;
}

void dumpProp(FILE *fp, Prop *prop)
/*
  Print all values in the structure prop to file fp.
*/
{
    if ((fp == NULL) || (prop == NULL)) { return; }

    fprintf(fp,
        "independent vars: %c, %c, deriv grade: %d, phase: %d\n"
        "================================================\n"
        " T = %g K\t"
        " t = %g C\t"
        " d = %g kg/m3\n"
        " p = %g bar\t"
        " s = %g kJ/kgK\t"
        " u = %g kJ/kg\n"
        " h = %g kJ/kg\t"
        " f = %g kJ/kg\t"
        " g = %g kJ/kg\n"
        " cp = %g kJ/kgK\t"
        " cv = %g kJ/kgK",
        prop->indep1, prop->indep2, prop->deriv, prop->phase,
        prop->T,
        prop->error == 0 ? prop->T - 273.15 : 0.0,
        prop->d,
        prop->p*1.0e-5,
        prop->s*1.0e-3,
        prop->u*1.0e-3,
        prop->h*1.0e-3,
        prop->f*1.0e-3,
        prop->g*1.0e-3,
        prop->cp*1.0e-3,
        prop->cv*1.0e-3);
    if (prop->phase == TWO) {
    fprintf(fp,"\tx = %g\n",prop->x);
    } else {
    fprintf(fp,"\n");
    }

    if (prop->deriv >= 1) {
    fprintf(fp,
        "\nfirst derivatives (SI-units):\n"
        "dp/dT = %g\t"
        "dp/dd = %g\n",
        prop->dp->T_Cd,
        prop->dp->d_CT);
    if (prop->phase == TWO) {
        fprintf(fp,
            "dx/dT = %g\t"
            "dx/dd = %g\n",
            prop->dx->T_Cd,
            prop->dx->d_CT);
    }
    if ((prop->indep1 == 'p') && (prop->indep2 == 'h')) {
        fprintf(fp,
            "dT/dp = %g\t"
            "dT/dh = %g\n"
            "dd/dp = %g\t"
            "dd/dh = %g\n"
            "ds/dp = %g\t"
            "ds/dh = %g\n"
            "du/dp = %g\t"
            "du/dh = %g\n",
            prop->dT->p_Ch,
            prop->dT->h_Cp,
            prop->dd->p_Ch,
            prop->dd->h_Cp,
            prop->ds->p_Ch,
            prop->ds->h_Cp,
            prop->du->p_Ch,
            prop->du->h_Cp);
        if (prop->phase == TWO) {
        fprintf(fp,
            "dx/dp = %g\t"
            "dx/dh = %g\n",
            prop->dx->p_Ch,
            prop->dx->h_Cp);
        }
    } else if ((prop->indep1 == 'p') && (prop->indep2 == 's')) {
        fprintf(fp,
            "dT/dp = %g\t"
            "dT/ds = %g\n"
            "dd/dp = %g\t"
            "dd/ds = %g\n"
            "du/dp = %g\t"
            "du/ds = %g\n"
            "dh/dp = %g\t"
            "dh/ds = %g\n",
            prop->dT->p_Cs,
            prop->dT->s_Cp,
            prop->dd->p_Cs,
            prop->dd->s_Cp,
            prop->du->p_Cs,
            prop->du->s_Cp,
            prop->dh->p_Cs,
            prop->dh->s_Cp);
        if (prop->phase == TWO) {
        fprintf(fp,
            "dx/dp = %g\t"
            "dx/ds = %g\n",
            prop->dx->p_Cs,
            prop->dx->s_Cp);
        }
    }
    }

    if (prop->deriv == 2) {
    fprintf(fp, "\nsecond derivatives (SI-units):\n"
        "dcv/dT   = %g\n"
        "dp/dT/dd = %g    "
        "dp/dT/dT = %g\n"
        "dp/dd/dT = %g    "
        "dp/dd/dd = %g\n",
        prop->dcv->T_Cd,
        prop->dp->dT_Cd->d_CT,
        prop->dp->dT_Cd->T_Cd,
        prop->dp->dd_CT->T_Cd,
        prop->dp->dd_CT->d_CT);
    if (prop->phase == TWO) {
        fprintf(fp,
            "dx/dT/dd = %g    "
            "dx/dT/dT = %g\n"
            "dx/dd/dT = %g    "
            "dx/dd/dd = %g\n",
            prop->dx->dT_Cd->d_CT,
            prop->dx->dT_Cd->T_Cd,
            prop->dx->dd_CT->T_Cd,
            prop->dx->dd_CT->d_CT);
    }
    if ((prop->indep1 == 'p') && (prop->indep2 == 'h')) {
        fprintf(fp,
            "dT/dp/dh = %g    "
            "dT/dp/dp = %g\n"
            "dT/dh/dp = %g    "
            "dT/dh/dh = %g\n"
            "dd/dp/dh = %g    "
            "dd/dp/dp = %g\n"
            "dd/dh/dp = %g    "
            "dd/dh/dh = %g\n"
            "ds/dp/dh = %g    "
            "ds/dp/dp = %g\n"
            "ds/dh/dp = %g    "
            "ds/dh/dh = %g\n"
            "du/dp/dh = %g    "
            "du/dp/dp = %g\n"
            "du/dh/dp = %g    "
            "du/dh/dh = %g\n"
            "dcv/dp   = %g    "
            "dcv/dh   = %g\n"
            "dcp/dp   = %g    "
            "dcp/dh   = %g\n",
            prop->dT->dp_Ch->h_Cp,
            prop->dT->dp_Ch->p_Ch,
            prop->dT->dh_Cp->p_Ch,
            prop->dT->dh_Cp->h_Cp,
            prop->dd->dp_Ch->h_Cp,
            prop->dd->dp_Ch->p_Ch,
            prop->dd->dh_Cp->p_Ch,
            prop->dd->dh_Cp->h_Cp,
            prop->ds->dp_Ch->h_Cp,
            prop->ds->dp_Ch->p_Ch,
            prop->ds->dh_Cp->p_Ch,
            prop->ds->dh_Cp->h_Cp,
            prop->du->dp_Ch->h_Cp,
            prop->du->dp_Ch->p_Ch,
            prop->du->dh_Cp->p_Ch,
            prop->du->dh_Cp->h_Cp,
            prop->dcv->p_Ch,
            prop->dcp->h_Cp,
            prop->dcp->p_Ch,
            prop->dcp->h_Cp);
        if (prop->phase == TWO) {
        fprintf(fp,
            "dx/dp/dh = %g    "
            "dx/dp/dp = %g\n"
            "dx/dh/dp = %g    "
            "dx/dh/dh = %g\n",
            prop->dx->dp_Ch->h_Cp,
            prop->dx->dp_Ch->p_Ch,
            prop->dx->dh_Cp->p_Ch,
            prop->dx->dh_Cp->h_Cp);
        }
    } else if ((prop->indep1 == 'p') && (prop->indep2 == 's')) {
        fprintf(fp,
            "dcv/dp   = %g    "
            "dcv/ds   = %g\n"
            "dcp/dp   = %g    "
            "dcp/ds   = %g\n",
            prop->dcv->p_Cs,
            prop->dcp->s_Cp,
            prop->dcp->p_Cs,
            prop->dcp->s_Cp);
    }
    }

    if (! prop->error) {
    fprintf(fp, " ***   OK    ***\n");
    } else {
    fprintf(fp, " *** Failure ***\n");
    }
}


/*************************/
/* !!! TEST-VERSIONS !!! */
/*************************/

double wthcond(double t, double d)
/*
   thermal conductivity as in thcond, but checked for two-phase regions.
   if IsInTwoPhaseRegion returns true, the values from the phase regions
   are linearly interpolated.
   Not safe around the critical point!!!
   !!! TEST-VERSION !!!
*/
{
    Prop *prop = newProp('x', 'x', 0);
    double lambda;
    double lambda_l;
    double lambda_v;
    double x;

    water_td(t, d, prop);
    if (prop->phase == ONE) {
    lambda = thcond(prop);
    } else {
        Prop *pliq = newProp('x', 'x', 0);
        Prop *pvap = newProp('x', 'x', 0);
    sat_t(t, pliq, pvap);
    lambda_l = thcond(pliq);
    lambda_v = thcond(pvap);
    x = (1.0 / d - 1.0 / pliq->d) / (1.0 / pvap->d - 1.0 / pliq->d);
    lambda = lambda_l + x * (lambda_v - lambda_l);
        pliq = freeProp(pliq);
        pvap = freeProp(pvap);
    }
    prop = freeProp(prop);
    return lambda;
}

double wviscos(double t, double d)
/*
   thermal conductivity as in thcond, but checked for two-phase regions.
   if IsInTwoPhaseRegion returns true, the values from the phase regions
   are linearly interpolated.
   Not safe around the critical point!!!
   !!! TEST-VERSION !!!
*/
{
    Prop *prop = newProp('x', 'x', 0);
    double eta;
    double eta_l;
    double eta_v;
    double x;

    water_td(t, d, prop);
    if (prop->phase == ONE) {
    eta = viscos(prop);
    } else {
    Prop *pliq = newProp('x', 'x', 0);
    Prop *pvap = newProp('x', 'x', 0);
    sat_t(t, pliq, pvap);
    eta_l = viscos(pliq);
    eta_v = viscos(pvap);
    x = (1.0 / d - 1.0 / pliq->d) / (1.0 / pvap->d - 1.0 / pliq->d);
    eta = d * (x * eta_v / pvap->d + (1.0 - x) * eta_l / pliq->d);
    /* aus Baehr/Stefan Waerme- und Stoffuebertragung Seite 495 */
    pliq = freeProp(pliq);
    pvap = freeProp(pvap);
    }
    prop = freeProp(prop);
    return eta;
}
