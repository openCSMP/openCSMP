// SPDX-FileCopyrightText: (C) 1998 Olaf Bauer <o.bauer@tu-harburg.de>
// SPDX-FileCopyrightText: (C) 1998 Ole Engel <engel@tu-harburg.de>
// SPDX-FileContributor: Max Mustermann
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: meta.cpp,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Olaf Bauer
# maintainer: Ole Engel <engel@tu-harburg.de>
# copyright:  GNU LIBRARY GENERAL PUBLIC LICENSE (see COPYING.LIB-2.0)
#---------------------------------------------------------------------
# Part of PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/

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
  steam quota can only be 0 <= x <= 1. If they fail (error =1),
  no saturation-like properties can be offered.

  For invalid (p, h)-points meta_ph, just like water_ph,
  extrapolates from bounds of validity towards (p, h).
*/

#include <cstdlib> /* NULL */
#include "steam4.h"
#include "iaps.h"

using namespace std;

void meta_td(double t, double d, Prop *prop)
/*
   Uses prop->phase (TWO or !TWO) as an info of how the fundamental equation
   should be evaluated. Therefore, in case prop->phase does not agree with
   real state if offers meta stable properties.
*/
{
    double psa, dl, dv, x;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    d *= 1.0e-3;
    if (! valid_td(t, d)) {
    prop->error = 1;
    return;
    }
    if (prop->phase != TWO) {
    td(t, d, &MPro, prop);  /* calculate supposed 1-phase (t, d) */
    } else if ((t >= tripl.t) && (t <= crit.t)) {
    /* calculate supposed 2-phase (t, d) */
    psat(t, &psa, &dl, &dv, &MLiq, &MPro);
    x = (1.0 / d - 1.0 / dl) / (1.0 / dv - 1.0 / dl);
    format_two(t, psa, x, dl, dv, &MLiq, &MPro, prop);
    } else {
    prop->error = 1;    /* no two-phase props available for given t */
    }
}


void meta_ph(double p, double h, double t, double d,
         double dp, double dh, Prop *prop)
/*
   Uses prop->phase (TWO or !TWO) as an info of how the fundamental equation
   should be evaluated. Therefore, in case prop->phase does not agree with
   real state if offers meta stable properties.
   For invalid (p, h) it extrapolates from bounds of validity
*/
{
    double ts, dl, dv, x;
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
    if (prop->phase != TWO) {
    /* calculate supposed 1-phase (p, h) */
    ph(p, h, &t, &d, dp, dh, &MPro, prop);
    } else if ((p >= tripl.p) && (p <= crit.p)) {
    /* calculate supposed 2-phase (p, h) */
    tsat(p, &ts, &dl, &dv, &MLiq, &MPro);
    x = (h - MLiq.spro.h) / (MPro.spro.h - MLiq.spro.h);
    format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
    } else {
    prop->error = 1;    /* no two-phase props available for given p */
    }
}


void meta_ps(double p, double s, double t, double d,
         double dp, double ds, Prop *prop)
/*
   Uses prop->phase (TWO or !TWO) as an info of how the fundamental equation
   should be evaluated. Therefore, in case prop->phase does not agree with
   real state if offers meta stable properties.
*/
{
    double ts, dl, dv, x;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    p *= 1.0e-6;
    s *= 1.0e-3;
    d *= 1.0e-3;
    if (! valid_ps(p, s)) {
    prop->error = 1;
    return;
    }
    if (prop->phase != TWO) {
    /* calculate supposed 1-phase (p, s) */
    ps(p, s, &t, &d, dp, ds, &MPro, prop);
    } else if ((p >= tripl.p) && (p <= crit.p)) {
    /* calculate supposed 2-phase (p, s) */
    tsat(p, &ts, &dl, &dv, &MLiq, &MPro);
    x = (s - MLiq.spro.s) / (MPro.spro.s - MLiq.spro.s);
    format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
    } else {
    prop->error = 1;    /* no two-phase props available for given p */
    }
}


void meta_hd(double h, double d, double t,
         double dh, Prop *prop)
/*
   Uses prop->phase (TWO or !TWO) as an info of how the fundamental equation
   should be evaluated. Therefore, in case prop->phase does not agree with
   real state if offers meta stable properties.
*/
{
    double ts, p, dl, dv, x;
    S_mliq MLiq;
    S_mpro MPro;

    /* check pointer */
    if (prop == NULL) { return; }

    h *= 1e-3;
    d *= 1e-3;
    if (! valid_hd(h, d)) {
    prop->error = 1;
    return;
    }
    if (prop->phase != TWO) {
    hd(h, d, &t, dh, &MPro, prop); /* calculate supposed 1-phase (h, d) */
    } else if ((d >= tripl.dv) && (d <= 1.0)
           && (h < 1547.745404137 + 169.3249912165 / d)) {
    /* calculate supposed 2-phase (h, d) */
    hdsat(h, d, dh, &ts, &p, &dl, &dv, &x, &MLiq, &MPro);
    format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
    } else if ((d > 0.2439) && (d < 0.4)
           && (h < 1547.2357851199 + 173.4098851329 / d)) {
    /* calc supp critical 2-phase (h, d) */
    hdsatc(h, d, dh, &ts, &p, &dl, &dv, &x, &MLiq, &MPro);
    format_two(ts, p, x, dl, dv, &MLiq, &MPro, prop);
    } else {
    prop->error = 1;     /* no two-phase props available */
    }
}
