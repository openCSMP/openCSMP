#ifndef MJL_CSMP_GLUE_H
#define MJL_CSMP_GLUE_H

#include "Point.h"
#include "MJL_Point.h"
#include "MJL_Point3D.h"
#include "CSMP_Definitions.h"

mjl::Point    assign( const csmp::Point<2U>& );
mjl::Point3D  assign( const csmp::Point<3U>& );

csmp::Point<2U>  assign( const mjl::Point& );
csmp::Point<3U>  assign( const mjl::Point3D& );

// definitions of inline functions
inline mjl::Point assign( const csmp::Point<2U>& p )
 {
    return mjl::Point( p[0], p[1] );
 }
 
inline mjl::Point3D  assign( const csmp::Point<3U>& p )
 {
    return mjl::Point3D( p[0], p[1], p[2] );
 }

inline csmp::Point<2U>  assign( const mjl::Point& p )
 {
	 return csmp::Point<2U>( p[0], p[1] );
 }
 
inline csmp::Point<3U>  assign( const mjl::Point3D& p )
 {
	 return csmp::Point<3U>( p[0], p[1], p[2] );
 }


#endif