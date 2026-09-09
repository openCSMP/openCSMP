// SPDX-FileCopyrightText: (C) 1998 Olaf Bauer <o.bauer@tu-harburg.de>
// SPDX-FileCopyrightText: (C) 1998 Ole Engel <engel@tu-harburg.de>
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: decl.h,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Olaf Bauer
# maintainer: Ole Engel <engel@tu-harburg.de>
# copyright:  GNU LIBRARY GENERAL PUBLIC LICENSE (see COPYING.LIB-2.0)
#---------------------------------------------------------------------
# IAPS formulation 1984
# structs for PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/


#ifndef _DECL_H_
#define _DECL_H_

typedef struct {
    double p, dpt, dpd, f, dft, dfd, g, dgt, dgd, s, dst, dsd,
    u, dut, dud, h, dht, dhd, cv, cp, x, dxt, dxd, dxtt, dxtd, dxdd,
    dptt, dptd, dpdd, dcvt, dhtt, dhtd, dhdd;
} S_pro;             /* props-saving structure */

typedef struct {
    double f, ft, ftt, fd, fdd, ftd, fttt, fttd, ftdd, fddd;
} S_br;              /* base- and residualfunction results */

typedef struct {
    double f, ft, ftt, fttt;
} S_id;              /* idealgas function results */

typedef struct {
    double b1, b2, b1t, b2t, b1tt, b2tt, b1ttt, b2ttt;
} S_mp;              /* internal results (molecular parameters)*/

typedef struct {
    double z1, z2, z3, z4, k, kt, ktt, kttt;
} S_rg;              /* internal results (realgas properties) */

typedef struct {
    double k[9], kt[9], ktt[9];
} S_glb;             /* internal results (local part of resid) */

typedef struct {
    double tau[4], del[4], dk[4], dl[4], k[4], kt[4], ktt[4], kd[4], kdd[4];
} S_loc;             /* internal results (local part of resid) */

typedef struct {
    S_pro spro;
    S_br  bs;
    S_br  rs;
    S_id  id;
    S_mp  mp;
    S_rg  rg;
    S_glb glb;
    S_loc loc;
} S_mpro;

typedef struct {
    S_pro spro;
    S_mp mp;
    S_rg rg;
    S_glb glb;
    S_loc loc;
} S_mliq;


/****************************
*  constant global structs  *
****************************/

typedef struct {
    double t, p, dl, dv;
} S_tripl;           /* triple point values */
extern const S_tripl tripl;

typedef struct {
    double t, p, d;
} S_crit;            /* critical point values */
extern const S_crit crit;

typedef struct {
    double t, p;
} S_creg;            /* critical region */
extern const S_creg creg;

typedef struct {
    double gas, tz, pz, uref, sref;
} S_con;             /* gas constant, reference values */
extern const S_con con;

typedef struct {
    double b1[4], b2[4];
} S_m;               /* coefficients for molecular parameters */
extern const S_m Cm;

typedef struct {
    double g[9][6], k[4], l[4], gg[4], t[4], d[4], a[4], b[4];
} S_r;               /* coefficients for residual function */
extern const S_r Cr;

typedef struct {
    double c[18];
} S_ig;               /* coefficients for ideal gas function */
extern const S_ig Ci;

#endif /* _DECL_H_ */
