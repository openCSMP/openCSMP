#ifndef __machine_h__
#define __machine_h__

/* Note special macros: ANSI_C (ANSI C syntax)
			SEGMENTED (segmented memory machine e.g. MS-DOS)
			MALLOCDECL (declared if malloc() etc have
					been declared) */

/* Platform dependent include and define statements added by SG */
#include "platform.h"

/*#define ANSI_C 1*/
/* #undef MALLOCDECL */
#define NOT_SEGMENTED 1
/* #undef HAVE_COMPLEX_H */
#define HAVE_MALLOC_H 1
#define STDC_HEADERS 1
#define HAVE_BCOPY 1
#define HAVE_BZERO 1
#define CHAR0ISDBL0 1
#define WORDS_BIGENDIAN 1
#define U_INT_DEF 1
#define PROTOTYPES_IN_STRUCT

/* number types on the Mac */
#ifdef PC
#define isascii(c)  (isprint(c) || iscntrl(c))
typedef unsigned int u_int;
#endif


/* for basic or larger versions */
#define COMPLEX 1
#define SPARSE 1

/* for loop unrolling */
/* #undef VUNROLL */
/* #undef MUNROLL */

/* for segmented memory */
#ifndef NOT_SEGMENTED
#define	SEGMENTED
#endif

/* if the system has malloc.h */
#ifdef HAVE_MALLOC_H
#define	MALLOCDECL	1
#endif

/* Check for ANSI C memmove and memset */
#ifdef STDC_HEADERS

/* standard copy & zero functions */
#define	MEM_COPY(from,to,size)	memmove((to),(from),(size))
#define	MEM_ZERO(where,size)	memset((where),'\0',(size))

#ifndef ANSI_C
#define ANSI_C 1
#endif

#endif

/* if have bcopy & bzero and no alternatives yet known, use them */
#ifdef HAVE_BCOPY
#ifndef MEM_COPY
/* nonstandard copy function */
#define	MEM_COPY(from,to,size)	bcopy((char *)(from),(char *)(to),(int)(size))
#endif
#endif

#ifdef HAVE_BZERO
#ifndef MEM_ZERO
/* nonstandard zero function */
#define	MEM_ZERO(where,size)	bzero((char *)(where),(int)(size))
#endif
#endif

/* if the system has complex.h */
#ifdef HAVE_COMPLEX_H
#include	<complex.h>
#endif

/* If prototypes are available & ANSI_C not yet defined, then define it,
	but don't include any header files as the proper ANSI C headers
        aren't here */
#define HAVE_PROTOTYPES 1
#ifdef HAVE_PROTOTYPES
#ifndef ANSI_C
#define ANSI_C  1
#endif
#endif

/* floating point precision */

/* you can choose single, double or long double (if available) precision */

#define FLOAT 		1
#define DOUBLE 		2
#define LONG_DOUBLE 	3

/* #undef REAL_FLT */
/* #undef REAL_DBL */

/* if nothing is defined, choose double precision */
#ifndef REAL_DBL
#ifndef REAL_FLT
#define REAL_DBL 1
#endif
#endif

/* single precision */
#ifdef REAL_FLT
#define  Real float
#define  LongReal float
#define REAL FLOAT
#define LONGREAL FLOAT
#endif

/* double precision */
#ifdef REAL_DBL
#define Real double
#define LongReal double
#define REAL DOUBLE
#define LONGREAL DOUBLE
#endif


/* machine epsilon or unit roundoff error */
/* This is correct on most IEEE Real precision systems */
#ifdef DBL_EPSILON
#if REAL == DOUBLE
#define	MACHEPS	DBL_EPSILON
#elif REAL == FLOAT
#define	MACHEPS	FLT_EPSILON
#elif REAL == LONGDOUBLE
#define MACHEPS LDBL_EPSILON
#endif
#endif

#define F_MACHEPS 1.19209e-07
#define D_MACHEPS 2.22045e-16

#ifndef MACHEPS
#if REAL == DOUBLE
#define	MACHEPS	D_MACHEPS
#elif REAL == FLOAT  
#define MACHEPS F_MACHEPS
#elif REAL == LONGDOUBLE
#define MACHEPS D_MACHEPS
#endif
#endif

/* #undef M_MACHEPS */

/********************
#ifdef DBL_EPSILON
#define	MACHEPS	DBL_EPSILON
#endif
#ifdef M_MACHEPS
#ifndef MACHEPS
#define MACHEPS	M_MACHEPS
#endif
#endif
********************/

#define	M_MAX_INT 2147483647
#ifdef	M_MAX_INT
#ifndef MAX_RAND
#define	MAX_RAND ((double)(M_MAX_INT))
#endif
#endif

/* for non-ANSI systems */
#ifdef HUGE_VAL
#ifndef HUGE
#define HUGE HUGE_VAL
#endif
#endif

/*#ifdef ANSI_C
extern	int	isatty(int);
#endif*/

double	__ip__(double	*dp1,
               double *dp2,
               int len );


void	__mltadd__(double *dp1,double *dp2,
                   double s, int len);

void	__smlt__(double *dp,double s, double	*out,
                 int len);

void	__add__(double	*dp1, double *dp2,double *out,
                int len);

void	__sub__(double	*dp1, double *dp2, double *out,
                int len);

void	__zero__(double *dp,int len);

void	__ip4__(double *v0,double *v1,double *v2,double *v3,double *w,double out[4], int len);

void	__lc4__(double *v0,double *v1,double *v2,double *v3,double *w,double	a[4], int len);

void	__ma4__(double *v0,double *v1,double *v2,double *v3,double *w,double	a[4], int len);

#endif

