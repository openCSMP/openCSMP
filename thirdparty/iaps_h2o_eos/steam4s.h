// SPDX-FileCopyrightText: © 1998 Olaf Bauer <o.bauer@tu-harburg.de>
// SPDX-FileCopyrightText: © 1998 Ole Engel <engel@tu-harburg.de>
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: steam4s.h,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Ole Engel
# maintainer: Ole Engel <engel@tu-harburg.de>
#---------------------------------------------------------------------
# All functions needed for using PROST (PROperties of water and STeam)
# with SMILE. Long headers are in steam4.h.
#-------------------------------------------------------------------*/

/*********************************************************************
          see file PROST4.TXT for help
**********************************************************************/


#ifndef _STEAM4_H_
#define _STEAM4_H_

#define ONE 1  /* one phase */
#define TWO 2  /* two phase */

#define LIQUID    4
#define SATURATED 5
#define VAPOUR    6

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

extern Prop *newProp(int, int, int);
extern Prop *freeProp(Prop *);
/* FILE is not declared in Smile, we can use void instead. */
extern void dumpProp(void *, Prop *);

extern void sat_t(double, Prop *, Prop *);
extern void sat_p(double, Prop *, Prop *);

extern void water_td(double, double, Prop *);
extern void water_tp(double, double, double, double, Prop *);
extern void water_ph(double, double, double, double,
             double, double, Prop *);
extern void water_ps(double, double, double, double,
             double, double, Prop *);
extern void water_hd(double, double, double, double, Prop *);
extern double speed(Prop *);
extern double wkappa(Prop *);
extern double wbetas(Prop *);
extern double viscos(Prop *);
extern double thcond(Prop *);
void meta_td(double, double, Prop *);
void meta_ph(double, double, double, double,
         double, double, Prop *);
void meta_ps(double, double, double, double,
         double, double, Prop *);
void meta_hd(double, double, double, double, Prop *);
extern void water_hs(double, double, double, double,
             double, double, Prop *);
extern void water_pd(double, double, double,
             double, Prop *);
extern void water_pu(double, double, double, double,
             double, double, Prop *);
extern void water_px(double, double, Prop *);
extern void water_sd(double, double, double,
             double, Prop *);
extern void water_th(double, double, double,
             double, Prop *);
extern void water_ts(double, double, double,
             double, Prop *);
extern void water_tu(double, double, double,
             double, Prop *);
extern void water_tx(double, double, Prop *);

extern void water_ud(double, double, double,
             double, Prop *);
extern void water_us(double, double, double, double,
             double, double, Prop *);
extern void water_dx(double, double, double,
             double, Prop *);
extern void water_dx0(double,
              double, Prop *);
extern void water_dx1(double,
              double, Prop *);
extern void water_dxm(double, double,
              double, Prop *);

extern double wthcond(double, double);
/* !!! TEST-VERSION !!! */

extern double wviscos(double, double);
/* !!! TEST-VERSION !!! */

#endif /* _STEAM4_H_ */
