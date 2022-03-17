#include "CGNS_ElementSpecifications.h"

#include <cstdio>
/* cgnslib.h file must be located in directory specified by -I during compile: */
#include "cgnslib.h"

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

/**
   CGNS elemnt types
typedef enum {
  CGNS_ENUMV( ElementTypeNull  ) =CG_Null,
  CGNS_ENUMV( ElementTypeUserDefined ) =CG_UserDefined,
  CGNS_ENUMV( NODE ) =2,
  CGNS_ENUMV( BAR_2 ) =3,
  CGNS_ENUMV( BAR_3 ) =4,
  CGNS_ENUMV( TRI_3 ) =5,
  CGNS_ENUMV( TRI_6 ) =6,
  CGNS_ENUMV( QUAD_4 ) =7,
  CGNS_ENUMV( QUAD_8 ) =8,
  CGNS_ENUMV( QUAD_9 ) =9,
  CGNS_ENUMV( TETRA_4 ) =10,
  CGNS_ENUMV( TETRA_10 ) =11,
  CGNS_ENUMV( PYRA_5 ) =12,
  CGNS_ENUMV( PYRA_14 ) =13,
  CGNS_ENUMV( PENTA_6 ) =14,
  CGNS_ENUMV( PENTA_15 ) =15,
  CGNS_ENUMV( PENTA_18 ) =16,
  CGNS_ENUMV( HEXA_8 ) =17,
  CGNS_ENUMV( HEXA_20 ) =18,
  CGNS_ENUMV( HEXA_27 ) =19,
  CGNS_ENUMV( MIXED ) =20,
  CGNS_ENUMV( PYRA_13 ) =21,
  CGNS_ENUMV( NGON_n ) =22,
  CGNS_ENUMV( NFACE_n ) =23,
  CGNS_ENUMV( BAR_4 ) =24,
  CGNS_ENUMV( TRI_9 ) =25,
  CGNS_ENUMV( TRI_10 ) =26,
  CGNS_ENUMV( QUAD_12 ) =27,
  CGNS_ENUMV( QUAD_16 ) =28,
  CGNS_ENUMV( TETRA_16 ) =29,
  CGNS_ENUMV( TETRA_20 ) =30,
  CGNS_ENUMV( PYRA_21 ) =31,
  CGNS_ENUMV( PYRA_29 ) =32,
  CGNS_ENUMV( PYRA_30 ) =33,
  CGNS_ENUMV( PENTA_24 ) =34,
  CGNS_ENUMV( PENTA_38 ) =35,
  CGNS_ENUMV( PENTA_40 ) =36,
  CGNS_ENUMV( HEXA_32 ) =37,
  CGNS_ENUMV( HEXA_56 ) =38,
  CGNS_ENUMV( HEXA_64 ) =39
} CGNS_ENUMT( ElementType_t );
*/


namespace csmp {

CGNS_ElementSpecifications::CGNS_ElementSpecifications()
 {
 }
  
CGNS_ElementSpecifications::~CGNS_ElementSpecifications()
 {
 }
 

std::string  CGNS_ElementSpecifications::CGNS_TypeName( int cgns_FEtype ) const
{
    // bar elements
    if ( cgns_FEtype == CGNS_ENUMV( BAR_2  ) )    return "BAR_2";
    if ( cgns_FEtype == CGNS_ENUMV( BAR_3  ) )    return "BAR_3";
    if ( cgns_FEtype == CGNS_ENUMV( BAR_4  ) )    return "BAR_4";
    // trapezoidal elements (quads )
    if ( cgns_FEtype == CGNS_ENUMV( QUAD_4  ) )   return "QUAD_4";
    if ( cgns_FEtype == CGNS_ENUMV( QUAD_8  ) )   return "QUAD_8";
    if ( cgns_FEtype == CGNS_ENUMV( QUAD_9  ) )   return "QUAD_9";
    // triangles
    if ( cgns_FEtype == CGNS_ENUMV( TRI_3  ) )    return "TRI_3";
    if ( cgns_FEtype == CGNS_ENUMV( TRI_6  ) )    return "TRI_6";
    // tetrahedral elements
    if ( cgns_FEtype == CGNS_ENUMV( TETRA_4  ) )  return "TETRA_4";
    if ( cgns_FEtype == CGNS_ENUMV( TETRA_10  ) ) return "TETRA_10";
    // hexahedral elements
    if ( cgns_FEtype == CGNS_ENUMV( HEXA_8  ) )   return "HEXA_8";
    if ( cgns_FEtype == CGNS_ENUMV( HEXA_20  ) )  return "HEXA_20";
    if ( cgns_FEtype == CGNS_ENUMV( HEXA_27  ) )  return "HEXA_27";
    // pyramids
    if ( cgns_FEtype == CGNS_ENUMV( PYRA_5  ) )   return "PYRA_5";
    if ( cgns_FEtype == CGNS_ENUMV( PYRA_14  ) )  return "PYRA_14";
    // pentahedra
    if ( cgns_FEtype == CGNS_ENUMV( PENTA_6  ) )  return "PENTA_6";
    if ( cgns_FEtype == CGNS_ENUMV( PENTA_15  ) ) return "PENTA_15";
    if ( cgns_FEtype == CGNS_ENUMV( PENTA_18  ) ) return "PENTA_18";

    std::cout <<"\nCGNS_ElementSpecifications::CGNS_TypeName: ";
    std::cout <<" Element type appears different from CGNS types ?";
    std::cout <<"\nType: "<< cgns_FEtype << std::endl;

    return "UNKNOWN_TYPE";
}

/** Map between CGNS integer flag corresponding to the CGNS element name string.
*/
int  CGNS_ElementSpecifications::CGNS_Type( const std::string& cgns_FEtype ) const
 {
    // bar elements
    if ( cgns_FEtype == "BAR_2" )    return CGNS_ENUMV( BAR_2 );
    if ( cgns_FEtype == "BAR_3" )    return CGNS_ENUMV( BAR_3 );
    if ( cgns_FEtype == "BAR_4" )    return CGNS_ENUMV( BAR_4 );
    // trapezoidal elements (quads)
    if ( cgns_FEtype == "QUAD_4" )   return CGNS_ENUMV( QUAD_4 );
    if ( cgns_FEtype == "QUAD_8" )   return CGNS_ENUMV( QUAD_8 );
    if ( cgns_FEtype == "QUAD_9" )   return CGNS_ENUMV( QUAD_9 );
    // triangles
    if ( cgns_FEtype == "TRI_3" )    return CGNS_ENUMV( TRI_3 );
    if ( cgns_FEtype == "TRI_6" )    return CGNS_ENUMV( TRI_6 );
    // tetrahedral elements
    if ( cgns_FEtype == "TETRA_4" )  return CGNS_ENUMV( TETRA_4 );
    if ( cgns_FEtype == "TETRA_10" ) return CGNS_ENUMV( TETRA_10 );
    // hexahedral elements
    if ( cgns_FEtype == "HEXA_8" )   return CGNS_ENUMV( HEXA_8 );
    if ( cgns_FEtype == "HEXA_20" )  return CGNS_ENUMV( HEXA_20 );
    if ( cgns_FEtype == "HEXA_27" )  return CGNS_ENUMV( HEXA_27 );
    // pyramids
    if ( cgns_FEtype == "PYRA_5" )   return CGNS_ENUMV( PYRA_5 );
    if ( cgns_FEtype == "PYRA_14" )  return CGNS_ENUMV( PYRA_14 );
    // pentahedra
    if ( cgns_FEtype == "PENTA_6" )  return CGNS_ENUMV( PENTA_6 );
    if ( cgns_FEtype == "PENTA_15" ) return CGNS_ENUMV( PENTA_15 );
    if ( cgns_FEtype == "PENTA_18" ) return CGNS_ENUMV( PENTA_18 );

    std::cout <<"\nCGNS_ElementSpecifications::CGNS_Type: ";
    std::cout <<" Element type appears different from CGNS types ?";
    std::cout <<"\nType: "<< cgns_FEtype << std::endl;

    return 0;

} // end CGNS_Type



/**

Returns the CSP finite element type enum identifier that corresponds to the
CGNS type as identified by an integer identifier. The user needs to specify
whether an isoparametric or globally interpolated element shall be
returned because CGNS elements can be converted into both.

Note however that CGNS fits quadratic or cubic elements to the bounding
curves. Thus, global interpolation will give incorrect results as
it works only for straight-sided elements.
 */
csmp::CSMP_FEM_TYPE  CGNS_ElementSpecifications::CSMP_TypeFrom_CGNS_Type( int etype, bool isoparametric, uint32_t dim ) const
 {
    if ( isoparametric ) {
        // bar elements
        if ( etype ==  CGNS_ENUMV( BAR_2 ) ) return csmp::ISOPARAMETRIC_LINEAR_BAR;
        if ( etype ==  CGNS_ENUMV( BAR_4 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_BAR;
        // trapezoidal elements (quads)
        if ( etype ==  CGNS_ENUMV( QUAD_4 ) ) return csmp::ISOPARAMETRIC_LINEAR_QUADRILATERAL;
        if ( etype ==  CGNS_ENUMV( QUAD_8 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_QUADRILATERAL;
        if ( etype ==  CGNS_ENUMV( QUAD_9 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9;
        // hexahedral elements
        if ( etype ==  CGNS_ENUMV( HEXA_8 ) ) return csmp::ISOPARAMETRIC_LINEAR_HEXAHEDRON;
        if ( etype ==  CGNS_ENUMV( HEXA_20 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20;
        if ( etype ==  CGNS_ENUMV( HEXA_27 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27;
        // triangles
        if ( etype ==  CGNS_ENUMV( TRI_3 ) ) return csmp::ISOPARAMETRIC_LINEAR_TRIANGLE;
        if ( etype ==  CGNS_ENUMV( TRI_6 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_TRIANGLE;
        // tetrahedral elements
        if ( etype ==  CGNS_ENUMV( TETRA_4 ) )  return csmp::ISOPARAMETRIC_LINEAR_TETRAHEDRON;
        if ( etype ==  CGNS_ENUMV( TETRA_10 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_TETRAHEDRON;
        // pyramids
        if ( etype ==  CGNS_ENUMV( PYRA_5 ) ) return csmp::ISOPARAMETRIC_LINEAR_PYRAMID;
        if ( etype ==  CGNS_ENUMV( PYRA_14 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_PYRAMID14;
        if ( etype ==  CGNS_ENUMV( PYRA_13 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_PYRAMID13;
        // pentahedra
        if ( etype ==  CGNS_ENUMV( PENTA_6 ) ) return csmp::ISOPARAMETRIC_LINEAR_PRISM;
        if ( etype ==  CGNS_ENUMV( PENTA_15 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_PRISM15;
        if ( etype ==  CGNS_ENUMV( PENTA_18 ) ) return csmp::ISOPARAMETRIC_QUADRATIC_PRISM18;
      }

    // bar elements
    if ( etype ==  CGNS_ENUMV( BAR_2 ) ) return csmp::LINEAR_BAR;
    if ( etype ==  CGNS_ENUMV( BAR_4 ) ) return csmp::QUADRATIC_BAR;
    // trapezoidal elements (quads)
    if ( etype == CGNS_ENUMV( QUAD_4 ) ) return csmp::LINEAR_RECTANGLE;
//    if ( etype == CGNS_ENUMV( QUAD_8 ) ) return csmp::QUADRATIC_QUADRILATERAL;
//    if ( etype == CGNS_ENUMV( QUAD_9 ) ) return csmp::QUADRATIC_QUADRILATERAL9;
    // hexahedral elements
    if ( etype ==  CGNS_ENUMV( HEXA_8 ) )  return csmp::LINEAR_CUBOID;
//    if ( etype ==  CGNS_ENUMV( HEXA_20 ) ) return csmp::QUADRATIC_HEXAHEDRON20;
//    if ( etype ==  CGNS_ENUMV( HEXA_27 ) ) return csmp::QUADRATIC_HEXAHEDRON27;
    // triangles
    if ( etype ==  CGNS_ENUMV( TRI_3 ) && dim == 3U ) return csmp::LINEAR_TRIANGLE;
    if ( etype ==  CGNS_ENUMV( TRI_3 ) && dim == 2U ) return csmp::LINEAR_TRIANGLE3D;
    if ( etype ==  CGNS_ENUMV( TRI_6 ) ) return csmp::QUADRATIC_TRIANGLE;
    // tetrahedral elements
    if ( etype ==  CGNS_ENUMV( TETRA_4 ) )  return csmp::LINEAR_TETRAHEDRON;
    if ( etype ==  CGNS_ENUMV( TETRA_10 ) ) return csmp::QUADRATIC_TETRAHEDRON;
    // pyramids
//    if ( etype ==  CGNS_ENUMV( PYRA_5 ) )  return csmp::LINEAR_PYRAMID;
//    if ( etype ==  CGNS_ENUMV( PYRA_14 ) ) return csmp::QUADRATIC_PYRAMID14;
//    if ( etype ==  CGNS_ENUMV( PYRA_13 ) ) return csmp::QUADRATIC_PYRAMID13;
    // pentahedra
//    if ( etype ==  CGNS_ENUMV( PENTA_6 ) ) return csmp::LINEAR_PRISM;
//    if ( etype ==  CGNS_ENUMV( PENTA_15 ) ) return csmp::QUADRATIC_PRISM15;
//    if ( etype ==  CGNS_ENUMV( PENTA_18 ) ) return csmp::QUADRATIC_PRISM18;

    std::cout <<"\nCGNS_ElementSpecifications::CSMP_TypeFromCGNSType";
    if ( isoparametric ) std::cout <<" (isoparametric): ";
    else                 std::cout <<": ";
    std::cout <<" Element type specifier cannot be parsed ";
    std::cout <<"\nType specifier: "<< etype << std::endl;
    return csmp::UNKNOWN;

 }

csmp::CSMP_FEM_TYPE  CGNS_ElementSpecifications::CSMP_TypeFrom_CGNS_TypeName( const std::string& CGNS_element_type, bool isoparametric, uint32_t dim ) const
{
    int etype( CGNS_Type( CGNS_element_type ) );
    return CSMP_TypeFrom_CGNS_Type( etype, isoparametric, dim );
}

std::string  CGNS_ElementSpecifications::CSMP_TypeNameFrom_CGNS_Type( int etype, bool isoparametric, uint32_t dim ) const
{
    return std::string(csmp::parseFiniteElementType( CSMP_TypeFrom_CGNS_Type( etype, isoparametric, dim ) ) );
}

std::string  CGNS_ElementSpecifications::CSMP_TypeNameFrom_CGNS_TypeName( const std::string& CGNS_element_type, bool isoparametric, uint32_t dim ) const
{
    return std::string(csmp::parseFiniteElementType( CSMP_TypeFrom_CGNS_TypeName( CGNS_element_type, isoparametric, dim ) ) );
}

int  CGNS_ElementSpecifications::CGNS_TypeFrom_CSMP_Type( int etype  ) const
{
    // bar elements
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_BAR ) return CGNS_ENUMV( BAR_2 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_BAR ) return CGNS_ENUMV( BAR_4 );
    if ( etype == csmp::LINEAR_BAR) return CGNS_ENUMV( BAR_2 );
    if ( etype == csmp::QUADRATIC_BAR ) return CGNS_ENUMV( BAR_4 );

    // trapezoidal elements (quads)
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_QUADRILATERAL || etype == csmp::LINEAR_RECTANGLE ) return CGNS_ENUMV( QUAD_4 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return CGNS_ENUMV( QUAD_8 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return CGNS_ENUMV( QUAD_9 );
//    if ( etype == csmp::LINEAR_QUADRILATERAL ) return CGNS_ENUMV( QUAD_4 );
//    if ( etype == csmp::QUADRATIC_QUADRILATERAL ) return CGNS_ENUMV( QUAD_8 );
//    if ( etype == csmp::QUADRATIC_QUADRILATERAL9 ) return CGNS_ENUMV( QUAD_9 );

    // hexahedral elements
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_HEXAHEDRON || etype == csmp::LINEAR_CUBOID ) return CGNS_ENUMV( HEXA_8 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return CGNS_ENUMV( HEXA_20 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return CGNS_ENUMV( HEXA_27 );
//    if ( etype == csmp::LINEAR_HEXAHEDRON ) return CGNS_ENUMV( HEXA_8 );
//    if ( etype == csmp::QUADRATIC_HEXAHEDRON20 ) return CGNS_ENUMV( HEXA_20 );
//    if ( etype == csmp::QUADRATIC_HEXAHEDRON27 ) return CGNS_ENUMV( HEXA_27 );

    // triangles
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_TRIANGLE ) return CGNS_ENUMV( TRI_3 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return CGNS_ENUMV( TRI_6 );
    if ( etype == csmp::LINEAR_TRIANGLE ) return CGNS_ENUMV( TRI_3 );
    if ( etype == csmp::LINEAR_TRIANGLE3D ) return CGNS_ENUMV( TRI_3 );
    if ( etype == csmp::QUADRATIC_TRIANGLE ) return CGNS_ENUMV( TRI_6 );

    // tetrahedral elements
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return CGNS_ENUMV( TETRA_4 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return CGNS_ENUMV( TETRA_10 );
    if ( etype == csmp::LINEAR_TETRAHEDRON ) return CGNS_ENUMV( TETRA_4 );
    if ( etype == csmp::QUADRATIC_TETRAHEDRON ) return CGNS_ENUMV( TETRA_10 );

    // pyramids
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_PYRAMID ) return CGNS_ENUMV( PYRA_5 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return CGNS_ENUMV( PYRA_14 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return CGNS_ENUMV( PYRA_13 );
//    if ( etype == csmp::LINEAR_PYRAMID ) return CGNS_ENUMV( PYRA_5 );
//    if ( etype == csmp::QUADRATIC_PYRAMID14 ) return CGNS_ENUMV( PYRA_14 );
//    if ( etype == csmp::QUADRATIC_PYRAMID13) return CGNS_ENUMV( PYRA_13 );

    // pentahedra
    if ( etype == csmp::ISOPARAMETRIC_LINEAR_PRISM ) return CGNS_ENUMV( PENTA_6 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_PRISM15 ) return CGNS_ENUMV( PENTA_15 );
    if ( etype == csmp::ISOPARAMETRIC_QUADRATIC_PRISM18 ) return CGNS_ENUMV( PENTA_18 );
//    if ( etype == csmp::LINEAR_PRISM ) return CGNS_ENUMV( PENTA_6 ); ;
//    if ( etype == csmp::QUADRATIC_PRISM15 ) return CGNS_ENUMV( PENTA_15 ); ;
//    if ( etype == csmp::QUADRATIC_PRISM18 ) return CGNS_ENUMV( PENTA_18 ); ;

    std::cout <<"\nCGNS_ElementSpecifications::CGNS_TypeFromCSMPType: ";
    std::cout <<" Element type specifier cannot be parsed ";
    std::cout <<"\nType specifier: "<< etype << std::endl;
    return csmp::UNKNOWN;
}

int  CGNS_ElementSpecifications::CGNS_TypeFrom_CSMP_TypeName( const std::string& CSMP_element_type ) const
{
    int etype( csmp::parseFiniteElementType( CSMP_element_type ) );
    return CGNS_TypeFrom_CSMP_Type( etype );
}

std::string  CGNS_ElementSpecifications::CGNS_TypeNameFrom_CSMP_Type( int etype ) const
{
    return CGNS_TypeName( CGNS_TypeFrom_CSMP_Type( etype ) );
}

std::string  CGNS_ElementSpecifications::CGNS_TypeNameFrom_CSMP_TypeName( const std::string& CSMP_element_type ) const
{
    return CGNS_TypeName( CGNS_TypeFrom_CSMP_TypeName( CSMP_element_type ) );
}

} // end namespace csmp
