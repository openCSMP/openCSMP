#include "CSMP_ElementSpecifications.h"
#include "Box.h"
#include "Exception.h"
#include "FiniteElement.h"

namespace csmp {

/**

Returns the CSMP finite element type name.
*/
std::string  CSMP_ElementSpecifications::CSMP_TypeName( int8_t etype )
{
    return parseFiniteElementType( etype );
} // end CSMP_TypeName( int, bool )


/** Returns the CSMP integer flag corresponding to the CSMP element name std::string.
*/
CSMP_FEM_TYPE  CSMP_ElementSpecifications::CSMP_Type( const std::string& csmp_FEtype )
{
    return parseFiniteElementType( csmp_FEtype );
} // end CSMP_Type


/**

Method deduces order of element interpolation functions from the
input finite-element type.

@section arguments Input Arguments

The CSMP type name of the finite element that shall be queried.

@return the interpolation order is returned as an unsigned integer value:
1 = linear, 2 = quadratic, 3 = cubic.

*/
size_t  CSMP_ElementSpecifications::InterpolationOrder( const std::string& etype )
 {
     if ( LinearElement( CSMP_Type(etype) ) )    return 1U;
     if ( QuadraticElement( CSMP_Type(etype) ) ) return 2U;
     if ( CubicElement( CSMP_Type(etype) ) )     return 3U;
     return 0U; // order not identified, probably constant as in a polygon
 }

size_t  CSMP_ElementSpecifications::InterpolationOrder( int8_t etype )
  {
    if ( LinearElement( etype ) )    return 1U;
    if ( QuadraticElement( etype ) ) return 2U;
    if ( CubicElement( etype ) )     return 3U;
    return 0U; // order not identified, probably constant as in a polygon
  }


/** Considers the element types.
 
    ISOPARAMETRIC_LINEAR_BAR,    						            // BAR_2            = 2,
    ISOPARAMETRIC_QUADRATIC_BAR, 						            // BAR_3            = 3,
    ISOPARAMETRIC_CUBIC_BAR,
    ISOPARAMETRIC_LINEAR_TRIANGLE,  					          // TRI_3            = 8,
    ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE, 		      // TRI_3_X          = 9,
    ISOPARAMETRIC_QUADRATIC_TRIANGLE, 					        // 2D & 3D TRI_6    = 10,
    ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE, 		  // TRI_6_X          = 11,
    ISOPARAMETRIC_CUBIC_TRIANGLE,
    ISOPARAMETRIC_LINEAR_TETRAHEDRON, 					        // TETRA_4          = 4,
    ISOPARAMETRIC_QUADRATIC_TETRAHEDRON, 			          // TETRA_10         = 5,
    ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON,    // TETRA_11
    ISOPARAMETRIC_CUBIC_TETRAHEDRON,
    ISOPARAMETRIC_LINEAR_PYRAMID,						            // PYRA_5           = 18,
    ISOPARAMETRIC_QUADRATIC_PYRAMID13, 					        // PYRA_13          = 24,
    ISOPARAMETRIC_QUADRATIC_PYRAMID14,     				      // PYRA_14          = 22,
    ISOPARAMETRIC_CUBIC_PYRAMID,
    ISOPARAMETRIC_LINEAR_PRISM,  						            // PENTA_6          = 12,
    ISOPARAMETRIC_QUADRATIC_PRISM15,         			      // PENTA_15         = 13,
    ISOPARAMETRIC_QUADRATIC_PRISM18,           			    // PENTA_18         = 21,
    ISOPARAMETRIC_CUBIC_PRISM,
    ISOPARAMETRIC_LINEAR_QUADRILATERAL,                 // QUAD_4           = 14,
    ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL,     // QUAD_4_X         = 15,
    ISOPARAMETRIC_QUADRATIC_QUADRILATERAL,              // QUAD_8           = 16,
    ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL,  // QUAD_8_X         = 17,
    ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9,  			      // QUAD_9           = 19,
    ISOPARAMETRIC_CUBIC_QUADRILATERAL,
    ISOPARAMETRIC_LINEAR_HEXAHEDRON,              		  // HEXA_8           = 6,
    ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20,       		    // HEXA_20          = 7,
    ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27, 				      // HEXA_27          = XX,
    ISOPARAMETRIC_CUBIC_HEXAHEDRON,
*/
bool CSMP_ElementSpecifications::UsesLocalCoordinates( int8_t etype )
 {
    // linear elements
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM ) return true;
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ) return true;
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return true;
    // quadratic elements
    if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return true;
    // cubic elements
    if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return true;
    return false;
 }






CSMP_FEM_TYPE CSMP_ElementSpecifications::CSMP_TypeUsingLocalCoordinates( int8_t etype )
 {
    // linear elements
    if ( etype == LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE3D ) return ISOPARAMETRIC_LINEAR_TRIANGLE;
    if ( etype == LINEAR_TETRAHEDRON ) return ISOPARAMETRIC_LINEAR_TETRAHEDRON;
    if ( etype == LINEAR_BAR ) return ISOPARAMETRIC_LINEAR_BAR;
    if ( etype == LINEAR_RECTANGLE ) return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
    if ( etype == LINEAR_CUBOID ) return ISOPARAMETRIC_LINEAR_HEXAHEDRON;
    if ( etype == BARYCENTRIC_LINEAR_TRIANGLE ) return ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE;
    return UNKNOWN;
    
    // TODO: handle the extra cases
    // quadratic elements
/*
    if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return true;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return true;
    // cubic elements
    if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return true;
    if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return true;
*/
 }



/**

Returns whether the CSMP element type specifier denotes a finite
element type that uses linear basis functions.

CSMP element types with such characteristics are:
@code
linear
==============
LINEAR_BAR
LINEAR_TRIANGLE
LINEAR_TRIANGLE3D
BARYCENTRIC_LINEAR_TRIANGLE
LINEAR_TETRAHEDRON
ISOPARAMETRIC_LINEAR_BAR
ISOPARAMETRIC_LINEAR_TRIANGLE
ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE
ISOPARAMETRIC_LINEAR_QUADRILATERAL
ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL
ISOPARAMETRIC_LINEAR_TETRAHEDRON
ISOPARAMETRIC_LINEAR_PYRAMID
ISOPARAMETRIC_LINEAR_PRISM
ISOPARAMETRIC_LINEAR_HEXAHEDRON

@endcode */
bool  CSMP_ElementSpecifications::LinearElement( int8_t etype )
 {
    if ( etype == LINEAR_BAR ||
         etype == LINEAR_TRIANGLE ||
         etype == LINEAR_TRIANGLE3D ||
		     etype == LINEAR_CUBOID ||
		     etype == LINEAR_RECTANGLE ||
         etype == BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == LINEAR_TETRAHEDRON ||
         etype == ISOPARAMETRIC_LINEAR_BAR ||
         etype == ISOPARAMETRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ||
         etype == ISOPARAMETRIC_LINEAR_PYRAMID ||
         etype == ISOPARAMETRIC_LINEAR_PRISM ||
         etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON)
        return true;
    return false;
 }




/**

Returns whether the CSMP element type specifier denotes a finite
element type that uses quadratic basis functions.
CSMP element types with such characteristics are:

@code
quadratic
==============
QUADRATIC_BAR
ISOPARAMETRIC_QUADRATIC_BAR
QUADRATIC_TRIANGLE
BARYCENTRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE
ISOPARAMETRIC_QUADRATIC_QUADRILATERAL
ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL
ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9
QUADRATIC_TETRAHEDRON
BARYCENTRIC_QUADRATIC_TETRAHEDRON
ISOPARAMETRIC_QUADRATIC_TETRAHEDRON
ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON
ISOPARAMETRIC_QUADRATIC_PYRAMID13
ISOPARAMETRIC_QUADRATIC_PYRAMID14
ISOPARAMETRIC_QUADRATIC_PRISM15
ISOPARAMETRIC_QUADRATIC_PRISM18
ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20
ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27

@endcode*/
bool  CSMP_ElementSpecifications::QuadraticElement( int8_t etype )
 {
    if ( etype == QUADRATIC_BAR ||
         etype == ISOPARAMETRIC_QUADRATIC_BAR ||
         etype == QUADRATIC_TRIANGLE ||
         etype == BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ||
         etype == QUADRATIC_TETRAHEDRON ||
         etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 )
        return true;
    return false;
 }



/**

Returns whether the CSMP element type specifier denotes a finite
element type that uses cubic basis functions.
CSMP element types with such characteristics are:

@code
cubic
==============
CUBIC_BAR
ISOPARAMETRIC_CUBIC_BAR
CUBIC_TRIANGLE
ISOPARAMETRIC_CUBIC_TRIANGLE
ISOPARAMETRIC_CUBIC_QUADRILATERAL
CUBIC_TETRAHEDRON
ISOPARAMETRIC_CUBIC_TETRAHEDRON
ISOPARAMETRIC_CUBIC_PYRAMID
ISOPARAMETRIC_CUBIC_PRISM
ISOPARAMETRIC_CUBIC_HEXAHEDRON

@endcode*/
bool  CSMP_ElementSpecifications::CubicElement( int8_t etype )
 {
    if ( etype == CUBIC_BAR ||
         etype == ISOPARAMETRIC_CUBIC_BAR ||
         etype == CUBIC_TRIANGLE ||
         etype == ISOPARAMETRIC_CUBIC_TRIANGLE ||
         etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ||
         etype == CUBIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_CUBIC_PYRAMID ||
         etype == ISOPARAMETRIC_CUBIC_PRISM ||
         etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON )
        return true;
    return false;
 }




/**

The element type used defines the spatial dimension of the model. It is
indicative only of the minimum dimension since lines and surfaces can
also exist in 3D space, which the name makes explicit.

@section arguments Input Arguments

Method expects a std::string argument that contains the ANSYS element type
to deduce the spatial dimension from.

@return size_t The minimum spatial dimension is returned as an unsigned integer ranging
from 1(=1D) to 3(=3D).

*/
uint32_t  CSMP_ElementSpecifications::MinimumSpatialDimension( const std::string& etype )
 {
    if ( LineElement(etype) )    return 1U;
    if ( SurfaceElement(etype) ) return 2U;
    return 3U;
 }
uint32_t  CSMP_ElementSpecifications::MinimumSpatialDimension( int8_t etype )
 {
    if ( LineElement(etype) )    return 1U;
    if ( SurfaceElement(etype) ) return 2U;
    return 3U;
 }


/**
Returns a list of line type elements (see FiniteElement.h).
"LINEAR_BAR"
"QUADRATIC_BAR"
"CUBIC_BAR"
"ISOPARAMETRIC_LINEAR_BAR"
"ISOPARAMETRIC_QUADRATIC_BAR"
"ISOPARAMETRIC_CUBIC_BAR"
*/
void  CSMP_ElementSpecifications::LineElements( std::list<std::string>& line_elements )
 {
    if ( !line_elements.empty() )
      line_elements.erase( line_elements.begin(), line_elements.end() );
    line_elements.push_back("LINEAR_BAR");
    line_elements.push_back("QUADRATIC_BAR");
    line_elements.push_back("CUBIC_BAR");
    // isoparametric elements
    line_elements.push_back("ISOPARAMETRIC_LINEAR_BAR");
    line_elements.push_back("ISOPARAMETRIC_QUADRATIC_BAR");
    line_elements.push_back("ISOPARAMETRIC_CUBIC_BAR");
 }

/** Returns true if the CSMP element name corresponds to a line element.
*/
bool  CSMP_ElementSpecifications::LineElement( const std::string& etype )
 {
    if ( etype == "LINEAR_BAR" ||
         etype == "QUADRATIC_BAR" ||
         etype == "CUBIC_BAR" ||
         etype == "ISOPARAMETRIC_LINEAR_BAR" ||
         etype == "ISOPARAMETRIC_QUADRATIC_BAR" ||
         etype == "ISOPARAMETRIC_CUBIC_BAR" )
        return true;
    return false;
 }

/** Returns true if the CSMP input element type is a line.
*/
bool  CSMP_ElementSpecifications::LineElement( int8_t etype )
 {
    if ( etype == LINEAR_BAR ||
         etype == QUADRATIC_BAR ||
         etype == CUBIC_BAR ||
         etype == ISOPARAMETRIC_LINEAR_BAR ||
         etype == ISOPARAMETRIC_QUADRATIC_BAR ||
         etype == ISOPARAMETRIC_CUBIC_BAR )
        return true;

    return false;
 }



/**
Returns a list of the CSMP surface element type names:
"LINEAR_TRIANGLE"
"LINEAR_TRIANGLE3D"
"BARYCENTRIC_LINEAR_TRIANGLE"
"QUADRATIC_TRIANGLE"
"BARYCENTRIC_QUADRATIC_TRIANGLE"
"CUBIC_TRIANGLE"
"ISOPARAMETRIC_LINEAR_TRIANGLE"
"ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE"
"ISOPARAMETRIC_QUADRATIC_TRIANGLE"
"ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE"
"ISOPARAMETRIC_CUBIC_TRIANGLE"
"ISOPARAMETRIC_LINEAR_QUADRILATERAL"
"ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL"
"ISOPARAMETRIC_QUADRATIC_QUADRILATERAL"
"ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL"
"ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9"
"ISOPARAMETRIC_CUBIC_QUADRILATERAL"
*/
void  CSMP_ElementSpecifications::SurfaceElements( std::list<std::string>& surf_elements )
 {
    if ( !surf_elements.empty() )
      surf_elements.erase( surf_elements.begin(), surf_elements.end() );

    surf_elements.push_back("LINEAR_TRIANGLE");
    surf_elements.push_back("LINEAR_TRIANGLE3D");
    surf_elements.push_back("BARYCENTRIC_LINEAR_TRIANGLE");
    surf_elements.push_back("QUADRATIC_TRIANGLE");
    surf_elements.push_back("BARYCENTRIC_QUADRATIC_TRIANGLE");
    surf_elements.push_back("CUBIC_TRIANGLE");
	  surf_elements.push_back("LINEAR_RECTANGLE");
    // isoparametric elements
    surf_elements.push_back("ISOPARAMETRIC_LINEAR_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_QUADRATIC_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_CUBIC_TRIANGLE");
    surf_elements.push_back("ISOPARAMETRIC_LINEAR_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_QUADRATIC_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL");
    surf_elements.push_back("ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9");
    surf_elements.push_back("ISOPARAMETRIC_CUBIC_QUADRILATERAL");
 }

/**
Returns true is the CSMP input element name is a surface element.
*/
bool  CSMP_ElementSpecifications::SurfaceElement( const std::string& etype )
 {
    if ( etype == "LINEAR_TRIANGLE" ||
		     etype == "LINEAR_RECTANGLE" ||
         etype == "LINEAR_TRIANGLE3D" ||
         etype == "BARYCENTRIC_LINEAR_TRIANGLE" ||
         etype == "QUADRATIC_TRIANGLE" ||
         etype == "BARYCENTRIC_QUADRATIC_TRIANGLE" ||
         etype == "CUBIC_TRIANGLE" ||
         etype == "ISOPARAMETRIC_LINEAR_TRIANGLE" ||
         etype == "ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE" ||
         etype == "ISOPARAMETRIC_QUADRATIC_TRIANGLE" ||
         etype == "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE" ||
         etype == "ISOPARAMETRIC_CUBIC_TRIANGLE" ||
         etype == "ISOPARAMETRIC_LINEAR_QUADRILATERAL" ||
         etype == "ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL" ||
         etype == "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL" ||
         etype == "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL" ||
         etype == "ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9" |
         etype == "ISOPARAMETRIC_CUBIC_QUADRILATERAL" )
      return true;
    return false;
 }

/** Returns true if the CSMP input element type is a surface element.
*/
bool  CSMP_ElementSpecifications::SurfaceElement( int8_t etype )
 {
    if ( etype == LINEAR_TRIANGLE ||
		     etype == LINEAR_RECTANGLE ||
         etype == LINEAR_TRIANGLE3D ||
         etype == BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == QUADRATIC_TRIANGLE ||
         etype == BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == CUBIC_TRIANGLE ||
         etype == ISOPARAMETRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_CUBIC_TRIANGLE ||
         etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ||
         etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL )
      return true;
    return false;
 }

/**
Returns a list of the CSMP element type names that denote volumetric elements:
"LINEAR_TETRAHEDRON"
"QUADRATIC_TETRAHEDRON"
"BARYCENTRIC_QUADRATIC_TETRAHEDRON"
"CUBIC_TETRAHEDRON"
"ISOPARAMETRIC_LINEAR_TETRAHEDRON"
"ISOPARAMETRIC_QUADRATIC_TETRAHEDRON"
"ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON"
"ISOPARAMETRIC_CUBIC_TETRAHEDRON"
"ISOPARAMETRIC_LINEAR_PYRAMID"
"ISOPARAMETRIC_QUADRATIC_PYRAMID13"
"ISOPARAMETRIC_QUADRATIC_PYRAMID14"
"ISOPARAMETRIC_CUBIC_PYRAMID"
"ISOPARAMETRIC_LINEAR_PRISM"
"ISOPARAMETRIC_QUADRATIC_PRISM15"
"ISOPARAMETRIC_QUADRATIC_PRISM18"
"ISOPARAMETRIC_CUBIC_PRISM"
"ISOPARAMETRIC_LINEAR_HEXAHEDRON"
"ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20"
"ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27"
"ISOPARAMETRIC_CUBIC_HEXAHEDRON"
*/
void  CSMP_ElementSpecifications::VolumeElements( std::list<std::string>& vol_elements )
 {
    if ( !vol_elements.empty() )
      vol_elements.erase( vol_elements.begin(), vol_elements.end() );
   
    vol_elements.push_back("LINEAR_TETRAHEDRON");
    vol_elements.push_back("QUADRATIC_TETRAHEDRON");
    vol_elements.push_back("BARYCENTRIC_QUADRATIC_TETRAHEDRON");
    vol_elements.push_back("CUBIC_TETRAHEDRON");
    vol_elements.push_back("LINEAR_CUBOID");
    // isoparametric elements
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_TETRAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_TETRAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_CUBIC_TETRAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_PYRAMID");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PYRAMID13");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PYRAMID14");
    vol_elements.push_back("ISOPARAMETRIC_CUBIC_PYRAMID");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_PRISM");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PRISM15");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_PRISM18");
    vol_elements.push_back("ISOPARAMETRIC_CUBIC_PRISM");
    vol_elements.push_back("ISOPARAMETRIC_LINEAR_HEXAHEDRON");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20");
    vol_elements.push_back("ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27");
    vol_elements.push_back("ISOPARAMETRIC_CUBIC_HEXAHEDRON");
 }


/**
Return true if the ANSYS input element name corresponds to a volumetric element.
*/
bool  CSMP_ElementSpecifications::VolumeElement( const std::string& etype )
 {
    if ( etype == "LINEAR_TETRAHEDRON" ||
		     etype == "LINEAR_CUBOID" ||
         etype == "QUADRATIC_TETRAHEDRON" ||
         etype == "BARYCENTRIC_QUADRATIC_TETRAHEDRON" ||
         etype == "CUBIC_TETRAHEDRON" ||
         etype == "ISOPARAMETRIC_LINEAR_TETRAHEDRON" ||
         etype == "ISOPARAMETRIC_QUADRATIC_TETRAHEDRON" ||
         etype == "ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON" ||
         etype == "ISOPARAMETRIC_CUBIC_TETRAHEDRON" ||
         etype == "ISOPARAMETRIC_LINEAR_PYRAMID" ||
         etype == "ISOPARAMETRIC_QUADRATIC_PYRAMID13" ||
         etype == "ISOPARAMETRIC_QUADRATIC_PYRAMID14" ||
         etype == "ISOPARAMETRIC_CUBIC_PYRAMID" ||
         etype == "ISOPARAMETRIC_LINEAR_PRISM" ||
         etype == "ISOPARAMETRIC_QUADRATIC_PRISM15" ||
         etype == "ISOPARAMETRIC_QUADRATIC_PRISM18" ||
         etype == "ISOPARAMETRIC_CUBIC_PRISM" ||
         etype == "ISOPARAMETRIC_LINEAR_HEXAHEDRON" ||
         etype == "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20" ||
         etype == "ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27" ||
         etype == "ISOPARAMETRIC_CUBIC_HEXAHEDRON" )
        return true;
    return false;
 }


/** Returns true if the ANSYS input element type is a volumetric element.
*/
bool  CSMP_ElementSpecifications::VolumeElement( int8_t etype )
 {
    if ( etype == LINEAR_TETRAHEDRON ||
		     etype == LINEAR_CUBOID ||
         etype == QUADRATIC_TETRAHEDRON ||
         etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == CUBIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_LINEAR_PYRAMID ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ||
         etype == ISOPARAMETRIC_CUBIC_PYRAMID ||
         etype == ISOPARAMETRIC_LINEAR_PRISM ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ||
         etype == ISOPARAMETRIC_CUBIC_PRISM ||
         etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ||
         etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON )
        return true;
    return false;
 }









/** 
    Returns the number of nodes of the CSMP element type as identified by
    the CSMP element integer code (enumeration CSMP_FEM_TYPE in FiniteElement.h).
*/
uint32_t CSMP_ElementSpecifications::NodesPerElementOfType( int8_t etype )
 {
    // bar
    if ( etype == ISOPARAMETRIC_LINEAR_BAR || etype == LINEAR_BAR )
        return 2U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_BAR || etype == QUADRATIC_BAR )
        return 3U;
    if ( etype == ISOPARAMETRIC_CUBIC_BAR || etype == CUBIC_BAR )
        return 4U;

    // tetrahedron
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON || etype == LINEAR_TETRAHEDRON )
        return 4U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON || etype == QUADRATIC_TETRAHEDRON )
        return 10U;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON || etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON )
        return 11U;
    if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON || etype == CUBIC_TETRAHEDRON )
        return 16U;

    // hexahedron
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON || etype == LINEAR_CUBOID )
        return 8U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 )
        return 20U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 )
        return 27U;
    if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON )
        return 32U;

    // quadrilateral
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL || etype == LINEAR_RECTANGLE )
        return 4U;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL )
        return 5U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL )
        return 8U;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL )
        return 9U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 )
        return 9U;
    if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL )
        return 12U;

    // triangle
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE3D )
        return 3U;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE || etype == BARYCENTRIC_LINEAR_TRIANGLE )
        return 4U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE || etype == QUADRATIC_TRIANGLE )
        return 6U;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE || etype == BARYCENTRIC_QUADRATIC_TRIANGLE )
        return 7U;
    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE || etype == CUBIC_TRIANGLE )
        return 9U;

    // prism
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM )
        return 6U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 )
        return 15U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 )
        return 18U;
    if ( etype == ISOPARAMETRIC_CUBIC_PRISM )
        return 24U;

    // pyramid
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID )
        return 5U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 )
        return 13U;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 )
        return 14U;
    if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID )
        return 21U;

    std::cerr <<"\nCSMP_ElementSpecifications::NodesPerElementOfType: ";
    std::cerr <<"unable to parse element type, returning 0"<< std::endl;
    return std::numeric_limits<uint32_t>::max(); // unknown element type

 } // end




/** Returns the  number of segments of the input CSMP finite element type
    as specified in CSMP_FEM_conventions.pdf   file in CSMP's documentation directory.
*/
uint32_t CSMP_ElementSpecifications::SegmentsPerElementOfType( int8_t etype )
 {
    // bar
    if ( etype == ISOPARAMETRIC_LINEAR_BAR || etype == LINEAR_BAR ||
         etype == ISOPARAMETRIC_QUADRATIC_BAR || etype == QUADRATIC_BAR ||
         etype == ISOPARAMETRIC_CUBIC_BAR || etype == CUBIC_BAR )
        return 2U;

    // tetrahedron
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON || etype == LINEAR_TETRAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON || etype == QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON || etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON || etype == CUBIC_TETRAHEDRON )
        return 6U;

    // hexahedron
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ||
         etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ||
		 etype == LINEAR_CUBOID)
        return 12U;

    // quadrilateral
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ||
         etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ||
		 etype == LINEAR_RECTANGLE)
        return 4U;

    // triangle
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE3D ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE || etype == BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE || etype == QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE || etype == BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_CUBIC_TRIANGLE || etype == CUBIC_TRIANGLE )
        return 3U;

    // prism
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ||
         etype == ISOPARAMETRIC_CUBIC_PRISM )
        return 9U;

    // pyramid
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ||
         etype == ISOPARAMETRIC_CUBIC_PYRAMID )
        return 8U;

    std::cerr <<"\nCSMP_ElementSpecifications::SegmentsPerElementOfType: ";
    std::cerr <<"unable to parse element type, returning 0"<< std::endl;
    return std::numeric_limits<uint32_t>::max(); // unknown element type
 } // end



uint32_t CSMP_ElementSpecifications::FacesPerElementOfType( int8_t CSMP_finite_element_type )
 {
    // since these are equivalent numbers
    return NeighborsPerElementOfType( CSMP_finite_element_type );
 } // end FacesPerElementOfType




/** Returns the maximum possible number of neighbors of the input CSMP finite element type
    as specified in CSMP_FEM_conventions.pdf in the documentation directory.
*/
uint32_t CSMP_ElementSpecifications::NeighborsPerElementOfType( int8_t etype )
 {
    // tetrahedron
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON || etype == LINEAR_TETRAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON || etype == QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON || etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON || etype == CUBIC_TETRAHEDRON )
        return 4U;

    // hexahedron
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ||
         etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ||
         etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ||
		 etype == LINEAR_CUBOID)
        return 6U;

    // quadrilateral
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ||
         etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ||
		 etype == LINEAR_RECTANGLE)
        return 4U;

    // triangle
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE3D ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE || etype == BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE || etype == QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE || etype == BARYCENTRIC_QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_CUBIC_TRIANGLE || etype == CUBIC_TRIANGLE )
        return 3U;

    // prism
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ||
         etype == ISOPARAMETRIC_CUBIC_PRISM )
        return 5U;

    // pyramid
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ||
         etype == ISOPARAMETRIC_CUBIC_PYRAMID )
        return 5U;

    // bar
    if ( etype == ISOPARAMETRIC_LINEAR_BAR || etype == LINEAR_BAR ||
         etype == ISOPARAMETRIC_QUADRATIC_BAR || etype == QUADRATIC_BAR ||
         etype == ISOPARAMETRIC_CUBIC_BAR || etype == CUBIC_BAR )
        return 2U;

    std::cerr <<"\nCSMP_ElementSpecifications::NeighborsPerElementOfType: ";
    std::cerr <<"unable to parse element type, returning 0"<< std::endl;
    return std::numeric_limits<uint32_t>::max(); // unknown element type
    
 } // end



/**

    Returns number of nodes that make up a particular face of the supplied 
    CSMP_FEM_TYPE enumeration as specified in FiniteElement.h
    
    Use this method to retrieve this information before any finite elements have been built.

    @attention 17/12/2015 SKM fixed method for quadratic elements
*/
uint32_t CSMP_ElementSpecifications::NodesPerFaceForElementOfType( int8_t etype,
                                                                   uint32_t face )
 {
    // bar
    // ---
    if ( etype == ISOPARAMETRIC_LINEAR_BAR || etype == LINEAR_BAR ||
         etype == ISOPARAMETRIC_QUADRATIC_BAR || etype == QUADRATIC_BAR ||
         etype == ISOPARAMETRIC_CUBIC_BAR || etype == CUBIC_BAR )
    {
        assert( face < 2U );
        return 1U;
    }

    // tetrahedron
    // -----------
    // linear elements
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON || etype == LINEAR_TETRAHEDRON  ) {
          assert( face < 4U );
          return 3U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON || etype == QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON || etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ) {
          assert( face < 4U );
          return 6U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON || etype == CUBIC_TETRAHEDRON ) {
          assert( face < 4U );
          return 8U;
      }

    // hexahedron
    // ----------
	if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON || etype == LINEAR_CUBOID ) {
          assert( face < 6U );
          return 4U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) {
          assert( face < 6U );
          return 8U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) {
          assert( face < 6U );
          return 9U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) {
          assert( face < 6U );
          return 13U;
      }

    // quadrilateral
    // -------------
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL || etype == LINEAR_RECTANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) {
         assert( face < 4U );
         return 2U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) {
         assert( face < 4U );
         return 3U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) {
         assert( face < 4U );
         return 12U;
      }

    // triangle
    // --------
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE3D ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE || etype == BARYCENTRIC_LINEAR_TRIANGLE ) {
          assert( face < 4U );
          return 2U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE || etype == QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE || etype == BARYCENTRIC_QUADRATIC_TRIANGLE ) {
          assert( face < 4U );
          return 3U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE || etype == CUBIC_TRIANGLE ) {
          assert( face < 4U );
          return 9U;
      }

    // prism
    // -----
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM ) {
          assert( face < 5U );
          if ( face==1U || face==2U || face==3U ) return 4U;
          return 3U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) {
          assert( face < 5U );
          if ( face==1U || face==2U || face==3U ) return 8U;
          return 6U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) {
          assert( face < 5U );
          if ( face==1U || face==2U || face==3U ) return 9U;
          return 6U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) {
          assert( face < 5U );
          if ( face==1U || face==2U || face==3U ) return 12U;
          return 9U;
      }

    // pyramid
    // -------
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ) {
          assert( face < 5U );
          if ( face == 4U ) return 4U;
          return 3U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) {
          assert( face < 5U );
          if ( face == 4U ) return 8U;
          return 6U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) {
          assert( face < 5U );
          if ( face == 4U ) return 9U;
          return 6U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) {
          assert( face < 5U );
          if ( face == 4U ) return 12U;
          return 8U;
      }

    std::cerr <<"\nCSMP_ElementSpecifications::NodesPerFaceForElementOfType: ";
    std::cerr <<"unable to parse element type, returning ULONG_MAX"<< std::endl;
    return std::numeric_limits<uint32_t>::max(); // unknown number of faces

 } // end NodesPerFaceForElementOfType







std::pair<uint32_t,uint32_t>  CSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType( int8_t CSMP_FE_type,
                                                                                                 uint32_t segm_id )
 {
    switch( CSMP_FE_type )
      {
         // tetrahedra
         case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
         case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON:
         case ISOPARAMETRIC_CUBIC_TETRAHEDRON:
         case LINEAR_TETRAHEDRON:
         case QUADRATIC_TETRAHEDRON:
         case BARYCENTRIC_QUADRATIC_TETRAHEDRON:
         case CUBIC_TETRAHEDRON :
           switch( segm_id ) {
               case 0: return std::make_pair( 0, 1 );
               case 1: return std::make_pair( 1, 2 );
               case 2: return std::make_pair( 2, 0 );
               case 3: return std::make_pair( 0, 3 );
               case 4: return std::make_pair( 1, 3 );
               case 5: return std::make_pair( 2, 3 );
               default:
                 std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
                 std::cerr <<"segment id="<< segm_id <<" out of range.\n";
                 return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
            }
         // hexahedra
         case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
         case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20:
         case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27:
         case ISOPARAMETRIC_CUBIC_HEXAHEDRON:
         case LINEAR_CUBOID :
           switch( segm_id ) {
               case 0: return std::make_pair( 0, 1 );
               case 1: return std::make_pair( 1, 2 );
               case 2: return std::make_pair( 2, 3 );
               case 3: return std::make_pair( 3, 0 );
               case 4: return std::make_pair( 0, 4 );
               case 5: return std::make_pair( 1, 5 );
               case 6: return std::make_pair( 2, 6 );
               case 7: return std::make_pair( 3, 7 );
               case 8: return std::make_pair( 4, 5 );
               case 9: return std::make_pair( 5, 6 );
               case 10: return std::make_pair( 6, 7 );
               case 11: return std::make_pair( 7, 4 );
               default:
                 std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
                 std::cerr <<"segment id="<< segm_id <<" out of range.\n";
                 return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
            }
         // prisms
         case ISOPARAMETRIC_LINEAR_PRISM:
         case ISOPARAMETRIC_QUADRATIC_PRISM15:
         case ISOPARAMETRIC_QUADRATIC_PRISM18:
         case ISOPARAMETRIC_CUBIC_PRISM :
           switch( segm_id ) {
               case 0: return std::make_pair( 0, 1 );
               case 1: return std::make_pair( 1, 2 );
               case 2: return std::make_pair( 2, 0 );
               case 3: return std::make_pair( 0, 3 );
               case 4: return std::make_pair( 1, 4 );
               case 5: return std::make_pair( 2, 5 );
               case 6: return std::make_pair( 3, 4 );
               case 7: return std::make_pair( 4, 5 );
               case 8: return std::make_pair( 5, 3 );
               default:
                 std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
                 std::cerr <<"segment id="<< segm_id <<" out of range.\n";
                 return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
            }
         // pyramids
         case ISOPARAMETRIC_LINEAR_PYRAMID:
         case ISOPARAMETRIC_QUADRATIC_PYRAMID13:
         case ISOPARAMETRIC_QUADRATIC_PYRAMID14:
         case ISOPARAMETRIC_CUBIC_PYRAMID :
           switch( segm_id ) {
               case 0: return std::make_pair( 0, 1 );
               case 1: return std::make_pair( 1, 2 );
               case 2: return std::make_pair( 2, 3 );
               case 3: return std::make_pair( 3, 0 );
               case 4: return std::make_pair( 0, 4 );
               case 5: return std::make_pair( 1, 4 );
               case 6: return std::make_pair( 2, 4 );
               case 7: return std::make_pair( 3, 4 );
               default:
                 std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
                 std::cerr <<"segment id="<< segm_id <<" out of range.\n";
                 return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
            }
         // triangles
         case ISOPARAMETRIC_LINEAR_TRIANGLE:
         case ISOPARAMETRIC_QUADRATIC_TRIANGLE:
         case ISOPARAMETRIC_CUBIC_TRIANGLE:
         case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE:
         case LINEAR_TRIANGLE:
         case LINEAR_TRIANGLE3D:
         case BARYCENTRIC_LINEAR_TRIANGLE :
           switch( segm_id ) {
               case 0: return std::make_pair( 1, 2 );
               case 1: return std::make_pair( 2, 0 );
               case 2: return std::make_pair( 0, 1 );
               default:
                 std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
                 std::cerr <<"segment id="<< segm_id <<" out of range.\n";
                 return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
            }
         // quadrilaterals
         case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
         case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL:
         case ISOPARAMETRIC_CUBIC_QUADRILATERAL:
         case LINEAR_RECTANGLE:
         case LINEAR_QUADRILATERAL :
           switch( segm_id ) {
               case 0: return std::make_pair( 0, 1 );
               case 1: return std::make_pair( 1, 2 );
               case 2: return std::make_pair( 2, 3 );
               case 3: return std::make_pair( 3, 4 );
               default:
                 std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
                 std::cerr <<"segment id="<< segm_id <<" out of range.\n";
                 return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
             }
         // line elements
         case ISOPARAMETRIC_LINEAR_BAR:
         case ISOPARAMETRIC_QUADRATIC_BAR:
         case ISOPARAMETRIC_CUBIC_BAR:
         case LINEAR_BAR:
         case QUADRATIC_BAR:
         case CUBIC_BAR :
           return std::make_pair( 0, 1 );
         default:
           std::cerr <<"\nCSMP_ElementSpecifications::CornerNodesPerSegmentForElementOfType: ";
           std::cerr <<"element type could not be identified.\n";
      }
      
    return std::pair<size_t,size_t>{ UINT_MAX, UINT_MAX };
    
 } // end CornerNodesPerSegmentForElementOfType






/**
    Returns the number of the node that is the n'th node of a particular face of the supplied
    element type as specified by CSMP_FEM_TYPE enumeration, see FiniteElement.h
    
    Use this method to retrieve this information before any finite elements have been built.
*/
uint32_t CSMP_ElementSpecifications::FaceNodeForElementOfType( int8_t etype,
                                                               uint32_t face,
                                                               uint32_t face_node )
 {
    // bar element
    // -----------
    if ( etype == ISOPARAMETRIC_LINEAR_BAR || etype == LINEAR_BAR )
      {
          assert( face < 2U );
          assert( face_node == 0U );
          // face 1 is adjacent to node 0 and 1 to node 1
          if ( face==0U ) return 0U;
          if ( face==1U ) return 1U;
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_BAR || etype == QUADRATIC_BAR )
      {
          assert( face < 2U );
          assert( face_node == 0U );
          // face 1 is adjacent to node 0 and 1 to node 1
          if ( face==0U ) return 0U;
          if ( face==1U ) return 1U;
      }
    if ( etype == ISOPARAMETRIC_CUBIC_BAR || etype == CUBIC_BAR )
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "elements with cubic interpolation functions are not handled yet");

    // tetrahedron
    // -----------
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON || etype == LINEAR_TETRAHEDRON )
      {
          assert( face < 4U );
          assert( face_node < 3U );
          if ( face==0U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 3U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 3U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 3U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 2U;
            }
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON || etype == QUADRATIC_TETRAHEDRON ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ||
         etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON )
      {
          assert( face < 4U );
          assert( face_node < 6U );
          if ( face==0U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 3U;
               // midside nodes
               if ( face_node==3U ) return 5U;
               if ( face_node==4U ) return 9U;
               if ( face_node==5U ) return 6U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 3U;
               // midside nodes
               if ( face_node==3U ) return 6U;
               if ( face_node==4U ) return 9U;
               if ( face_node==5U ) return 7U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 3U;
               // midside nodes
               if ( face_node==3U ) return 4U;
               if ( face_node==4U ) return 8U;
               if ( face_node==5U ) return 7U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 2U;
               // midside nodes
               if ( face_node==3U ) return 4U;
               if ( face_node==4U ) return 5U;
               if ( face_node==5U ) return 6U;
            }
      }
    if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON || etype == CUBIC_TETRAHEDRON )
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "elements with cubic interpolation functions are not handled yet");

    // hexahedron
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON || etype == LINEAR_CUBOID )
      {
          assert( face < 6U );
          assert( face_node < 4U );
          if ( face==0U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 2U;
               if ( face_node==3U ) return 1U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 5U;
               if ( face_node==3U ) return 4U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 6U;
               if ( face_node==3U ) return 5U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 7U;
               if ( face_node==3U ) return 6U;
            }
          if ( face==4U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 4U;
               if ( face_node==2U ) return 7U;
               if ( face_node==3U ) return 3U;
            }
          if ( face==5U ) {
               if ( face_node==0U ) return 4U;
               if ( face_node==1U ) return 5U;
               if ( face_node==2U ) return 6U;
               if ( face_node==3U ) return 7U;
            }
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 )
      {
          assert( face < 6U );
          assert( face_node < 8U );
          if ( face==0U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 2U;
               if ( face_node==3U ) return 1U;
               // mid-edge nodes
               if ( face_node==4U ) return 11U;
               if ( face_node==5U ) return 10U;
               if ( face_node==6U ) return 9U;
               if ( face_node==7U ) return 8U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 5U;
               if ( face_node==3U ) return 4U;
               // mid-edge nodes
               if ( face_node==4U ) return 8U;
               if ( face_node==5U ) return 10U;
               if ( face_node==6U ) return 16U;
               if ( face_node==7U ) return 12U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 6U;
               if ( face_node==3U ) return 5U;
               // mid-edge nodes
               if ( face_node==4U ) return 9U;
               if ( face_node==5U ) return 14U;
               if ( face_node==6U ) return 17U;
               if ( face_node==7U ) return 10U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 7U;
               if ( face_node==3U ) return 6U;
               // mid-edge nodes
               if ( face_node==4U ) return 13U;
               if ( face_node==5U ) return 15U;
               if ( face_node==6U ) return 18U;
               if ( face_node==7U ) return 14U;
            }
          if ( face==4U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 4U;
               if ( face_node==2U ) return 7U;
               if ( face_node==3U ) return 3U;
               // mid-edge nodes
               if ( face_node==4U ) return 11U;
               if ( face_node==5U ) return 12U;
               if ( face_node==6U ) return 19U;
               if ( face_node==7U ) return 15U;
            }
          if ( face==5U ) {
               if ( face_node==0U ) return 4U;
               if ( face_node==1U ) return 5U;
               if ( face_node==2U ) return 6U;
               if ( face_node==3U ) return 7U;
               // mid-edge nodes
               if ( face_node==4U ) return 19U;
               if ( face_node==5U ) return 16U;
               if ( face_node==6U ) return 17U;
               if ( face_node==7U ) return 18U;
            }
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 || etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON )
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "elements with cubic interpolation functions are not handled yet");

    // quadrilateral
    // -------------
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL || etype == LINEAR_RECTANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL )
      {
          assert( face < 4U );
          assert( face_node < 2U );
          if ( face==0U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 3U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 3U;
               if ( face_node==1U ) return 0U;
            }
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ||
         etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 )
      {
          assert( face < 4U );
          assert( face_node < 3U );
          if ( face==0U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 4U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 5U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 6U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 3U;
               if ( face_node==1U ) return 0U;
               if ( face_node==2U ) return 7U;
            }
      }
    if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL )
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "elements with cubic interpolation functions are not handled yet");

    // triangle
    // --------
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE || etype == LINEAR_TRIANGLE3D )
      {
          assert( face < 3U );
          assert( face_node < 2U );
          if ( face==0U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 0U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
            }
      }
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE || etype == BARYCENTRIC_LINEAR_TRIANGLE ||
         etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE || etype == QUADRATIC_TRIANGLE ||
         etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE || etype == BARYCENTRIC_QUADRATIC_TRIANGLE )
      {
          assert( face < 3U );
          assert( face_node < 3U );
          if ( face==0U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 4U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 0U;
               if ( face_node==2U ) return 5U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 3U;
            }
      }
    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE || etype == CUBIC_TRIANGLE )
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "elements with cubic interpolation functions are not handled yet");

    // prism
    // -----
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM )
      {
          assert( face < 5U );
          if ( face==0U ) {
               assert( face_node < 3U );
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 1U;
            }
          if ( face==1U ) {
               assert( face_node < 4U );
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 4U;
               if ( face_node==3U ) return 3U;
            }
          if ( face==2U ) {
               assert( face_node < 4U );
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 5U;
               if ( face_node==3U ) return 4U;
            }
          if ( face==3U ) {
               assert( face_node < 4U );
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 5U;
               if ( face_node==3U ) return 2U;
            }
          if ( face==4U ) {
               assert( face_node < 3U );
               if ( face_node==0U ) return 3U;
               if ( face_node==1U ) return 4U;
               if ( face_node==2U ) return 5U;
            }
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ||
         etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ||
         etype == ISOPARAMETRIC_CUBIC_PRISM )
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "prism elements with quadratic or cubic interpolation functions are not handled yet");

    // pyramid
    // -------
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID )
      {
          assert( face < 5U );
          if ( face==0U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 1U;
               if ( face_node==2U ) return 4U;
            }
          if ( face==1U ) {
               if ( face_node==0U ) return 1U;
               if ( face_node==1U ) return 2U;
               if ( face_node==2U ) return 4U;
            }
          if ( face==2U ) {
               if ( face_node==0U ) return 2U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 4U;
            }
          if ( face==3U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 4U;
               if ( face_node==2U ) return 3U;
            }
          if ( face==4U ) {
               if ( face_node==0U ) return 0U;
               if ( face_node==1U ) return 3U;
               if ( face_node==2U ) return 2U;
               if ( face_node==3U ) return 1U;
            }
      }
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ||
         etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ||
         etype == ISOPARAMETRIC_CUBIC_PYRAMID ) {
      std::cerr <<"\nCSMP_ElementSpecifications::FaceNodeForElementOfType: ";
      std::cerr <<"unable to parse element type, returning ULONG_MAX"<< std::endl;
      throw csmp::Exception( ERROR, "SMP_ElementSpecifications::FaceNodeForElementOfType:",
                            "pyramid elements with quadratic or cubic interpolation functions are not handled yet");
      }

    return std::numeric_limits<uint32_t>::max();

 } // end FaceNodeForElementOfType

} // end namespace csmp
