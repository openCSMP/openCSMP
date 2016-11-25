#ifndef CSMP_DEFINITIONS_H
#define CSMP_DEFINITIONS_H

/**
@file CSMP_definitions.h
@brief Global declarations and definitions
@author S.K. Matthai (1995)

Global declarations and definitions for the CSMP++ code for the simulation
of complex earth systems.
This file is included by all headers in the C++ namespace csmp which contains
the CSMP++ application programmer interface. To make
porting to different platforms easier, this file includes all necessary
C++ / STL header files.
*/


/*
=============================================================================
Intel math library if available
(because it affects the standard namespace this block must come first)
=============================================================================
*/

#include <cmath>


/*
==============================================================
Backward compatibility with C
==============================================================
*/

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstdarg>
#include <cctype>
#include <cstring>
#include <ctime>
#include <cfloat>
#include <climits>
#include <cassert>


/*
==============================================================
C++ ANSI Standard Compliance / Standard Template Library (STL)
==============================================================
*/
#include <limits>
#include <algorithm>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <istream>
#include <sstream>
#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <array>
#include <ciso646>

/*
===========================
g++ issues
===========================
*/

#ifndef to_string
template<typename T> std::string to_string( T& number ) {
     std::ostringstream  stream;
     stream << number;
     return stream.str();
  }
#endif




/*
============================================================================
Platform-independent integer and floating point types. NOTE:

!!! The 'CSMP_number_types.h' file must be generated for each compiler
    and platform by compiling and running the program:

     'CSMP_number_types_main.cpp'

    in the main-library source code directory. Only after this has been
    done, CSMP should be compiled.

The CSMP size_t variable will be equivalent to size_t, which is used
by the STL as index variable for random access arrays.

============================================================================
*/
#include "CSMP_number_types.h"

/*
=======================
JPEG Interface
=======================
*/
#define CSMP_WITH_IMAGE_OUTPUT


/**
@mainpage About

The Complex Systems Modeling Platform (CSMP++) is an application programmer interface (API) for the simulation of THMC processes and their interactions in
complex space and time domains trough applicaton of combinations of FEM and FVM methods.
It was created in 1995 by Stephan K. Matthai and Stephen G. Roberts at Stanford University.

CSMP has been designed for the simulation of complex physics in geometrically complex and scale-variant domains. 
Its modular structure facilitates the analysis of the emergent properties of the studied system via large-scale mathematical simulations with models 
incorporating fundamental laws of nature, integrated for scenarios based on observations, problems and possible solutions.
Governing PDEs are integrated using FEM and FVM methods and combinations thereof. 

CSMP can be licensed for academic, commercial and public sector  (NGO) use. The license can be obtained from the CSMP Originator group.
Currently (March 2016) the official contact point is the Technology Transfer of the ETHZ, Switzerland ("Meyns Silke (F&W)" <silke.meyns@sl.ethz.ch>).
*/


/**
@defgroup CSMPglobalTypedefs CSMP number types for current computing platform
*/

#include "Index.h"
#include "CSMP_global_enumerations.h"

#endif
