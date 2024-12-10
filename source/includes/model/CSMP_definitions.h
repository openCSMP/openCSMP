#ifndef CSMP_DEFINITIONS_H
#define CSMP_DEFINITIONS_H

/**
@file CSMP_definitions.h
@brief Global declarations and definitions
@author S.K. Matthai (1995)

Global declarations and definitions for the CSMP++ code for the simulation
of complex natural systems.
This file includes the most important C++ header files to make the use
of the library easier for the less experienced cross platform-, cross-compiler user and
It is included by most headers in the CSMP C++ library,
Last updated for use with C++17 in November 2021.
*/

/*
==============================================================
Backward compatibility with C
==============================================================
*/

#include <cmath>
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
#include <algorithm>
#include <array>
#include <chrono>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <type_traits>
#include <utility>
#include <vector>
#include <unordered_set>
#include <unordered_map>

/*
=====================================================
    Matthew Bentham's flexible storage container used
    to manage the storage of Elements, Faces, InterFaces and Nodes in MeshManager
    https://plflib.org/colony.htm
=====================================================
*/
//#include "plf_colony.h"


/*
=======================
JPEG Interface
to create images etc.
=======================
*/
#define CSMP_WITH_IMAGE_OUTPUT


#if defined(DEBUG) && defined(NDEBUG)
#error && "DEBUG and NDEBUG cannot be defined simultaneously"
#endif

/**
@mainpage About

The Complex Systems Modeling Platform (Open CSMP++) is an application programmer interface (API)
for the simulation of THMC processes and their interactions in complex space and time domains.
It primarily relies on combinations of finite element (FEM) and finite volue (FVM) methods for this purpose.
Open CSMP++  was created as CSP in 1995 by Stephan K. Matthai and Stephen G. Roberts at Stanford University.

Open CSMP++ has been designed for the simulation of complex physics in geometrically complex and scale-variant domains.
Its modular structure aims to support analysis of the emergent properties of the studied system via mathematical simulations.
Governing PDEs are integrated using FEM and FVM methods and combinations thereof. 

This core library of Open CSMP++ is  distributed under the L-GPL license for general and including commercial use.
Until Dezember 2021,  CSMP++ had to be licensed and the official contact point was the Technology Transfer of the ETHZ, Switzerland
("Meyns Silke (F&W)" <silke.meyns@sl.ethz.ch>).

*/


 /**
     Tests wheter a  size_t   or any other unsigned integer has been assigned a maximum value to signify that it has not been initialised.
     (example: raw size_t's in CSMP are initialised to UINT_MAX, just as floating point values are initialised to a quiet NaN).
 */
 template<typename T>
 bool constexpr isUninitialisedInteger( T i ) {
     static_assert( std::is_integral<T>::value, "hasIntegerUnitialisedDefaultValue: inappropriate argument type for this function." );
     static_assert( std::is_arithmetic<T>::value, "hasIntegerUnitialisedDefaultValue: inappropriate argument type for this function." );
     static_assert( std::is_unsigned<T>::value, "hasIntegerUnitialisedDefaultValue: inappropriate argument type for this function." );
     return ( i == std::numeric_limits<T>::max() );
  }


/**
@defgroup CSMPglobalTypedefs CSMP number types for current computing platform
*/

#include "Index.h"
#include "CSMP_global_enumerations.h"


#endif /* CSMP_DEFINITIONS_H */
