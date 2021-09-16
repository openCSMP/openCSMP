#ifndef CSMP_NUMBER_TYPES_H
#define CSMP_NUMBER_TYPES_H

/**
@file CSMP_number_types.h
@author S.K. Matthai

// ------------------------------------------------

// CSMP number types for current computing platform

SKM comments:

@todo Use double for all your calculations and temp variables.
Use float when you need to maintain an array of numbers - float[] (if precision is sufficient),
and you are dealing with over tens of thousands of float numbers.

Many/most C++ math functions or operators convert/return double,
and you don't want to cast the numbers back to float for any intermediate steps.

E.g. If you have an input of 100,000 numbers from a file or a stream and need to sort them, put the numbers in a float[].

Since you can fit twice as many floats in processor cache as you can with doubles,
and memory latency is likely to be the main bottleneck in many programs,
keeping a whole working set of floats warm in cache may be literally an order of magnitude faster than using doubles and having them spill to RAM.

// ------------------------------------------------

*/

namespace csmp { 

/**
@addtogroup CSMPglobalTypedefs
@{
*/

typedef char                 char8;
typedef unsigned char        uchar8;
typedef char                 int8;
typedef unsigned char        uint8;
typedef short                int16;
typedef unsigned short       uint16;
typedef int                  int32;
typedef unsigned int         uint32;
typedef long                 long64;
typedef unsigned long        ulong64;
typedef float                float32;
typedef double               double64;
typedef long double          double128;

/**
@}
*/

 } // end namespace csmp


#endif
