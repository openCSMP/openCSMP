#ifndef MESCHACH_PLATFORM_H
#define MESCHACH_PLATFORM_H

/* Platform dependent include and define statements */
/* Added by SG 03/11/2005 */
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstddef>
#include <cstring>
#include <cfloat>
/* unix idenfication function for whether a terminal is used */
//extern	int isatty(int __fd) __THROW;
#include <unistd.h>
#endif

#ifdef __MWERKS__
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <cfloat>
#include <cstddef>
#include <cstring>
#include <unix.h>
#define PC 
#define CODE_WARRIOR 
extern	int	isatty(int);
#else
#ifndef __INTEL_COMPILER
#ifndef _WIN32
#include <unistd.h>
#endif
#ifdef _WIN32
#include <io.h>
#define isatty _isatty
#endif
//#include <ctype>
#include <ctype.h>
#define _fileno fileno
#endif
#endif

// _fileno (unix.h) is only defined on non-unix systems
#ifndef __MWERKS__
#ifndef _fileno
int _fileno(FILE* stream ) { return fileno(stream); }
#endif
#endif

#endif
