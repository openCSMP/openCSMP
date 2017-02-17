#ifndef CSMP_ELEMENT_SPECIFICATIONS_H
#define CSMP_ELEMENT_SPECIFICATIONS_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"

namespace csmp {

/**

@brief Data class to deduce properties of the CSMP finite-element types 
prior to the establishment of a csmp::Model from the
type information that is provided by @enum CSMP_FEM_TYPE
declared in file "FiniteElment.h"
 
@author S.K. Matthai
@author R. Mansipov
@date 2003,2015
@note SKM 9/2/17 made static singleton class

Used to deduce element characteristics prior to the construction of CSMP Element classes.

*/
class CSMP_ElementSpecifications {
  public:
    CSMP_ElementSpecifications() = delete;
    ~CSMP_ElementSpecifications() = delete;

    static CSMP_FEM_TYPE CSMP_Type( const std::string& FEtype );
    static std::string   CSMP_TypeName( int32 CSMP_finite_element_type );

    static size_t        InterpolationOrder( const std::string& etype );
    static size_t        InterpolationOrder( int32 etype );
    static bool          LinearElement( int32 etype );
    static bool          QuadraticElement( int32 etype );
    static bool          CubicElement( int32 etype );

    static size_t        MinimumSpatialDimension( const std::string& CSMP_finite_element_type );
    static size_t        MinimumSpatialDimension( int32 CSMP_finite_element_type );

    static bool          LineElement( int32 CSMP_finite_element_type );
    static bool          SurfaceElement( int32 CSMP_finite_element_type );
    static bool          VolumeElement( int32 CSMP_finite_element_type );

    static bool          LineElement( const std::string& CSMP_finite_element_type );
    static bool          SurfaceElement( const std::string& CSMP_finite_element_type );
    static bool          VolumeElement( const std::string& CSMP_finite_element_type );

    static void          LineElements( std::list<std::string>& line_elements );
    static void          SurfaceElements( std::list<std::string>& surf_elements );
    static void          VolumeElements( std::list<std::string>& vol_elements );

    static size_t        NodesPerElementOfType( int32 CSMP_finite_element_type );
    static size_t        FacesPerElementOfType( int32 CSMP_finite_element_type );
    static size_t        NeighborsPerElementOfType( int32 CSMP_finite_element_type );
    static size_t        NodesPerFaceForElementOfType( int32 CSMP_finite_element_type, size_t face );
    static size_t        FaceNodeForElementOfType( int32 CSMP_finite_element_type, size_t face, size_t face_node );
};

} // csmp


#endif
