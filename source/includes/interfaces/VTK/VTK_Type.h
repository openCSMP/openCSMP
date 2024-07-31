#ifndef VTK_TYPE_H
#define VTK_TYPE_H

#include "FiniteElement.h"

namespace csmp {

/// polygonal data types supported by the Visual Tool Kit (VTK) from kitware.com
enum VTK_TYPE // keep default int because it is written to file
{ 
  VTK_VERTEX=1,
  VTK_POLY_VERTEX=2,
  VTK_LINE=3,
  VTK_POLYLINE=4,
  VTK_TRIANGLE=5,
  VTK_TRIANGLE_STRIP=6,
  VTK_POLYGON=7,
  VTK_QUAD=9,
  VTK_TETRA=10,
  VTK_HEXAHEDRON=12,
  VTK_WEDGE=13,
  VTK_PYRAMID=14,
  VTK_QUADRATIC_EDGE=21,
  VTK_QUADRATIC_TRIANGLE=22,
  VTK_QUADRATIC_QUAD=23,
  VTK_QUADRATIC_TETRA=24,
  VTK_QUADRATIC_HEXAHEDRON=25,
  VTK_QUADRATIC_WEDGE=26,
  VTK_QUADRATIC_PYRAMID=27,
  VTK_BIQUADRATIC_QUAD=28,
  VTK_BIQUADRATIC_TRIANGLE=34
};

inline std::string parseVTK_TYPE( int type_id ) {
    switch( type_id ) {
        case VTK_VERTEX:
      return "VTK_VERTEX";
        case VTK_POLY_VERTEX:
      return "VTK_POLY_VERTEX";
        case VTK_LINE:
      return "VTK_LINE";
        case VTK_POLYLINE:
      return "VTK_POLYLINE";
        case VTK_TRIANGLE:
      return "VTK_TRIANGLE";
        case VTK_TRIANGLE_STRIP:
      return "VTK_TRIANGLE_STRIP";
        case VTK_POLYGON:
      return "VTK_POLYGON";
        case VTK_QUAD:
      return "VTK_QUAD";
        case VTK_TETRA:
      return "VTK_TETRA";
        case VTK_HEXAHEDRON:
      return "VTK_HEXAHEDRON";
        case VTK_WEDGE:
      return "VTK_WEDGE";
        case VTK_PYRAMID:
      return "VTK_PYRAMID";
        case VTK_QUADRATIC_EDGE:
      return "VTK_QUADRATIC_EDGE";
        case VTK_QUADRATIC_TRIANGLE:
      return "VTK_QUADRATIC_TRIANGLE";
        case VTK_QUADRATIC_QUAD:
      return "VTK_QUADRATIC_QUAD";
        case VTK_QUADRATIC_TETRA:
      return "VTK_QUADRATIC_TETRA";
        case VTK_QUADRATIC_HEXAHEDRON:
      return "VTK_QUADRATIC_HEXAHEDRON";
        case VTK_QUADRATIC_WEDGE:
      return "VTK_QUADRATIC_WEDGE";
        case VTK_QUADRATIC_PYRAMID:
      return "VTK_QUADRATIC_PYRAMID";
        case VTK_BIQUADRATIC_QUAD:
      return "VTK_BIQUADRATIC_QUAD";
        case VTK_BIQUADRATIC_TRIANGLE:
      return "VTK_IQUADRATIC_TRIANGLE";
        default:
          std::string failure_message{"\nparseVTK_TYPE: type_id="};
          failure_message += std::to_string(type_id);
          failure_message += " out-of-range of VTK_TYPE identifiers.\n";
          return failure_message.c_str();
      }
    return "undefined";
 }
 
 
/** from FiniteElement.h
enum CSMP_FEM_TYPE : std::int8_t { UNKNOWN,
                    LINEAR_BAR,    										                  // BAR_2            = 2,
                    QUADRATIC_BAR, 										                  // BAR_3            = 3,
                    CUBIC_BAR, 										                      // BAR_4
                    LINEAR_TRIANGLE,
                    LINEAR_TRIANGLE3D, 									                // TRI_3            = 8,
					          LINEAR_QUADRILATERAL,
					          LINEAR_CUBOID,
					          LINEAR_RECTANGLE,
					          BARYCENTRIC_LINEAR_TRIANGLE, 			    	           	// TRI_3_X          = 9,
                    QUADRATIC_TRIANGLE, 								                // 2D & 3D TRI_6    = 10,
                    BARYCENTRIC_QUADRATIC_TRIANGLE, 					          // TRI_6_X          = 11,
                    CUBIC_TRIANGLE,
                    LINEAR_TETRAHEDRON, 								                // TETRA_4          = 4,
                    QUADRATIC_TETRAHEDRON, 								              // TETRA_10         = 5,
                    BARYCENTRIC_QUADRATIC_TETRAHEDRON,                  // PYRA_5           = 18,
                    CUBIC_TETRAHEDRON,
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
                    ZERO_DIMENSIONAL_FACE,                              // the face of a line element
                    POINT_ELEMENT,
                    POLYGONAL_ELEMENT,
                    POLYHEDRAL_ELEMENT,
                    EXPERIMENTAL_ELEMENT,
 */
 inline CSMP_FEM_TYPE csmpTypeFromVTK_TYPE( VTK_TYPE type_id )
  {
    switch( type_id ) {
        case VTK_VERTEX:
      return POINT_ELEMENT;
        case VTK_POLY_VERTEX:
      return POLYGONAL_ELEMENT;
        case VTK_LINE:
      return ISOPARAMETRIC_LINEAR_BAR;
        case VTK_POLYLINE:
//      return "VTK_POLYLINE";
        case VTK_TRIANGLE:
      return ISOPARAMETRIC_LINEAR_TRIANGLE;
        case VTK_TRIANGLE_STRIP:
//      return "VTK_TRIANGLE_STRIP";
        case VTK_POLYGON:
      return POLYGONAL_ELEMENT;
        case VTK_QUAD:
      return ISOPARAMETRIC_LINEAR_QUADRILATERAL;
        case VTK_TETRA:
      return ISOPARAMETRIC_LINEAR_TETRAHEDRON;
        case VTK_HEXAHEDRON:
      return ISOPARAMETRIC_LINEAR_HEXAHEDRON;
        case VTK_WEDGE:
      return ISOPARAMETRIC_LINEAR_PRISM;
        case VTK_PYRAMID:
      return ISOPARAMETRIC_LINEAR_PYRAMID;
        case VTK_QUADRATIC_EDGE:
      return ISOPARAMETRIC_QUADRATIC_BAR;
        case VTK_QUADRATIC_TRIANGLE:
      return ISOPARAMETRIC_QUADRATIC_TRIANGLE;
        case VTK_QUADRATIC_QUAD:
      return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;
        case VTK_QUADRATIC_TETRA:
      return ISOPARAMETRIC_QUADRATIC_TETRAHEDRON;
        case VTK_QUADRATIC_HEXAHEDRON:
      return ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20; // check whether this element matches correctly
        case VTK_QUADRATIC_WEDGE:
      return ISOPARAMETRIC_QUADRATIC_PRISM15;
        case VTK_QUADRATIC_PYRAMID:
      return ISOPARAMETRIC_QUADRATIC_PYRAMID13;
        case VTK_BIQUADRATIC_QUAD:
      return ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
        case VTK_BIQUADRATIC_TRIANGLE:
      return ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE;
        default:
          throw std::logic_error( "ERROR csmpTypeFromVTK_TYPE: no corresponding CSMP element found, returning UNKNOWN");
      }
    return UNKNOWN;
  } // end
 


inline int nodesPer_VTK_TYPE( VTK_TYPE type_id )
  {
    switch( type_id ) {
        case VTK_VERTEX:
      return 1;
        case VTK_POLY_VERTEX:
      return UNSPECIFIED;
        case VTK_LINE:
      return 2;
        case VTK_POLYLINE:
      return UNSPECIFIED;
        case VTK_TRIANGLE:
      return 3;
        case VTK_TRIANGLE_STRIP:
      return UNSPECIFIED;
        case VTK_POLYGON:
      return UNSPECIFIED;
        case VTK_QUAD:
      return 4;
        case VTK_TETRA:
      return 4;
        case VTK_HEXAHEDRON:
      return 8;
        case VTK_WEDGE:
      return 6;
        case VTK_PYRAMID:
      return 5;
        case VTK_QUADRATIC_EDGE:
      return 3;
        case VTK_QUADRATIC_TRIANGLE:
      return 6;
        case VTK_QUADRATIC_QUAD:
      return 8;
        case VTK_QUADRATIC_TETRA:
      return 8;
        case VTK_QUADRATIC_HEXAHEDRON:
      return 20;
        case VTK_QUADRATIC_WEDGE:
      return 15;
        case VTK_QUADRATIC_PYRAMID:
      return 13;
        case VTK_BIQUADRATIC_QUAD:
      return 9;
        case VTK_BIQUADRATIC_TRIANGLE:
      return 7;
        default:
          throw std::logic_error( "ERROR nodesPer_VTK_TYPE: VTK_Type not recognised");
      }
    return UNSPECIFIED;
    
  } // end nodesPer_VTK_TYPE



// TODO: put the conversions for the attributes integer and floating-point types here as well

} // csmp

#endif
