// SPDX-FileCopyrightText: © 1998 Olaf Bauer <o.bauer@tu-harburg.de>
// SPDX-FileContributor: Erika Mustermann
//
// SPDX-License-Identifier: LGPL-2.0-only

/*------------------------------------------------------------*- C -*-
# $Id: more.h,v 1.2 2006/10/11 19:31:13 mattha2k Exp $
#---------------------------------------------------------------------
# author(s):  Olaf Bauer
#             Ole Engel
# maintainer: Ole Engel <engel@tu-harburg.de>
# copyright:  GNU LIBRARY GENERAL PUBLIC LICENSE (see COPYING.LIB-2.0)
#---------------------------------------------------------------------
# Part of PROST (PROperties of water and STeam)
#-------------------------------------------------------------------*/


#ifndef _MORE_H_
#define _MORE_H_

#include "steam4.h"
#include "decl.h"

/**************************
 *                        *
 *   internal functions   *
 *                        *
 **************************/

/* validate functions aren't implemented yet */
int valid_hs(double h, double s);
int valid_pd(double p, double d);
int valid_pu(double p, double u);
int valid_sd(double s, double d);
int valid_th(double t, double h);
int valid_ts(double t, double s);
int valid_tu(double t, double u);
int valid_ud(double u, double d);
int valid_us(double u, double s);

/* internal functions doesn't check the pointers, they must no be NULL */
void hs(double h, double s, double *t, double *d,
    double delh, double dels, S_mpro *mPro, Prop *prop);
void ht(double h, double t, double *d, double delh, S_mpro *mPro, Prop *prop);
void pd(double p, double d, double *t, double delp, S_mpro *mPro, Prop *prop);
void pu(double p, double u, double *t, double *d,
    double delp, double delu, S_mpro *mPro, Prop *prop);
void sd(double s, double d, double *t, double dels, Prop *prop,
    S_mliq *mLiq, S_mpro *mPro);
void st(double s, double t, double *d, double dels, S_mpro *mPro, Prop *prop);
void ud(double u, double d, double *t, double delu, S_mpro *mPro, Prop *prop);
void us(double u, double s, double *t, double *d,
    double delu, double dels, S_mpro *mPro, Prop *prop);
void ut(double u, double t, double *d, double delu, S_mpro *mPro, Prop *prop);

#endif /* _MORE_H_ */
