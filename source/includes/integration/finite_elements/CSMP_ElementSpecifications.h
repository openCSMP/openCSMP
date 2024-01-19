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
@note SKM 9/12/21 extras for more independent VData construction

Used to deduce element characteristics prior to the construction of CSMP Element classes.

*/
class CSMP_ElementSpecifications {
  public:
    CSMP_ElementSpecifications() = delete;
    ~CSMP_ElementSpecifications() = delete;

    static CSMP_FEM_TYPE CSMP_Type( const std::string& FEtype );
    /// finds the parametric finite element type that matches the argument one that uses a global coordinate system
    static CSMP_FEM_TYPE CSMP_TypeUsingLocalCoordinates( int8_t CSMP_finite_element_type );
    static std::string   CSMP_TypeName( int8_t CSMP_finite_element_type );

    static uint32_t      InterpolationOrder( const std::string& etype );
    static uint32_t      InterpolationOrder( int8_t etype );
    /// reference element is stored in parametric space and calculation outcomes are transformed into physical space by Jacobian transformation
    static bool          UsesLocalCoordinates( int8_t etype );
    /// calls previous method because CSMP++ currently has no sub- or super-parametric elements
    static bool          IsIsoparametric( int8_t etype ) { return UsesLocalCoordinates(etype); }
    static bool          LinearElement( int8_t etype );
    static bool          QuadraticElement( int8_t etype );
    static bool          CubicElement( int8_t etype );

    static uint32_t      MinimumSpatialDimension( const std::string& CSMP_finite_element_type );
    static uint32_t      MinimumSpatialDimension( int8_t CSMP_finite_element_type );

    static bool          LineElement( int8_t CSMP_finite_element_type );
    static bool          SurfaceElement( int8_t CSMP_finite_element_type );
    static bool          VolumeElement( int8_t CSMP_finite_element_type );

    static bool          LineElement( const std::string& CSMP_finite_element_type );
    static bool          SurfaceElement( const std::string& CSMP_finite_element_type );
    static bool          VolumeElement( const std::string& CSMP_finite_element_type );

    static void          LineElements( std::list<std::string>& line_elements );
    static void          SurfaceElements( std::list<std::string>& surf_elements );
    static void          VolumeElements( std::list<std::string>& vol_elements );

    static uint32_t      NodesPerElementOfType( int8_t CSMP_finite_element_type );
    static uint32_t      SegmentsPerElementOfType( int8_t CSMP_finite_element_type );
    static uint32_t      FacesPerElementOfType( int8_t CSMP_finite_element_type );
    static uint32_t      NeighborsPerElementOfType( int8_t CSMP_finite_element_type );
    static uint32_t      NodesPerFaceForElementOfType( int8_t CSMP_finite_element_type, uint32_t face );
    static std::pair<uint32_t,uint32_t>  CornerNodesPerSegmentForElementOfType( int8_t CSMP_finite_element_type, uint32_t segm );
    static uint32_t      FaceNodeForElementOfType( int8_t CSMP_finite_element_type, uint32_t face, uint32_t face_node );
};

} // csmp


#endif
