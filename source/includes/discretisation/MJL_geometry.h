#ifndef MJL_GEOMETRY_H
#define MJL_GEOMETRY_H

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <cmath>

namespace mjl {

/// @file CSMP_definitions.h

#ifndef DBL_MAX
#define DBL_MAX  1.7976931348623157E+308  /* max decimal value of a "double"*/
#endif

/**
  @addtogroup CSMPglobalEnums
  @{
*/

enum LOCATION { LEFT, RIGHT, BEYOND, BEHIND, BETWEEN, ORIGIN, DESTINATION };
enum INTERSECTION { COLLINEAR, PARALLEL, SKEW, SKEW_CROSS, SKEW_NO_CROSS };
enum POINT_CLASSIFICATION { OUTSIDE, INSIDE, BOUNDARY };
enum EDGE_CLASSIFICATION { TOUCHING, CROSSING, INESSENTIAL };
enum POLYGON_CLASSIFICATION { UNKNOWN, P_IS_INSIDE, Q_IS_INSIDE };

/**
  @}
*/

} // end namespace mjl

#endif
