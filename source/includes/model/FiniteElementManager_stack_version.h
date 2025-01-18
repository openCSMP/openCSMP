#ifndef CSMP_FINITE_ELEMENT_MANAGER1_H
#define CSMP_FINITE_ELEMENT_MANAGER1_H

#include "LinearLineElement.h"
#include "IsoparametricQuadraticLineElement.h"
#include "LinearTriangle.h"
#include "LinearTriangle3D.h"
#include "LinearTetrahedron.h"
#include "LinearCuboid.h"
#include "LinearRectangle.h"
#include "IsoparametricQuadraticTriangle.h"
#include "IsoparametricQuadraticTetrahedron.h"

#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearQuadrilateral.h"

#include "IsoparametricQuadraticHexahedron.h"
#include "IsoparametricQuadraticPyramid.h"
#include "IsoparametricQuadraticPrism.h"
#include "IsoparametricQuadraticTetrahedron.h"
#include "IsoparametricQuadraticQuadrilateral.h"


namespace csmp {

/**
@brief Alternative version of FiniteElementManager where all elements are stored on the stack (initialised at compile time in the body of the class)
   
   @attention IS SLOWER THAN MANAGER WHERE ELEMENTS ARE CREATED ON THE HEAP!!!! - and less flexible! (SKM30/6/2024)

 @atttention needs specialisations for
   - 1, 2, and 3D
   - with and without local coordinate system
   - with linear and with quadratic interpolation function
   
   @note USING_LOCAL_COORDS refers to whether the element is defined in parametric space coordinates (r,s,t) like all the isoparametric elements in CSMP

@author S.K. Matthai
@date 2024

*/
template<uint32_t dim, int INTPOL_ORDER, bool USING_LOCAL_COORDS>
class FiniteElementManager1 {
  public:
    FiniteElementManager1();
    
    // dummy function does nothing
    void InitializeElements( uint32_t /* dimensions */, uint32_t /* interpolation_order */, bool /* isoparametric */ ) {}

    bool ContainsElementType( CSMP_FEM_TYPE etype ) const;
    void CurrentElementTypes( std::list<CSMP_FEM_TYPE>& etypes ) const;
    int  InterpolationOrder() const { return INTPOL_ORDER; }
    bool UsesElementsWithLocalCoordinateSystem() const { return USING_LOCAL_COORDS; }
    
    /// get finite element pointer from enumeration type (TODO: maybe a reference should be returned here?)
    FiniteElement* const E( int8_t csmp_etype );
    FiniteElement* const E( CSMP_FEM_TYPE csmp_etype );
    
    /// const version needed
    const FiniteElement* const E( CSMP_FEM_TYPE csmp_etype ) const;
    
    uint32_t NodesOfElementType( CSMP_FEM_TYPE etype );

    // standard types
    FiniteElement* const LinearBarElement();
    FiniteElement* const LinearTriangleElement();
    FiniteElement* const LinearTetrahedronElement();

//    FiniteElement*    LineElementPointer();

    void              Out() const;
  
  private:
    // line elements
    LinearLineElement                   line_elmt1; // (DIM);
    IsoparametricLinearLineElement      line_elmt2;
    IsoparametricQuadraticLineElement   line_elmt3;
    // triangles
    LinearTriangle                      tria_elmt1a;
    LinearTriangle3D                    tria_elmt1b;
    IsoparametricLinearTriangle         tria_elmt2;
    IsoparametricLinearLineElement      tria_elmt3;
    // quadrilaterals
    LinearRectangle                     quad_elmt1;
    IsoparametricLinearQuadrilateral    quad_elmt2;
    IsoparametricQuadraticQuadrilateral quad_elmt3;
    // tetrahedra
    LinearTetrahedron                   tet_elmt1;
    IsoparametricLinearTetrahedron      tet_elmt2;
    IsoparametricQuadraticTetrahedron   tet_elmt3;
    // hexahedra
    LinearCuboid                        hexa_elmt1;
    IsoparametricLinearHexahedron       hexa_elmt2;
    IsoparametricQuadraticHexahedron    hexa_elmt3;
    // prisms
    IsoparametricLinearPrism            prism_elmt2;
    IsoparametricQuadraticPrism         prism_elmt3;
    // pyramids
    IsoparametricLinearPyramid          pyra_elmt2;
    IsoparametricQuadraticPyramid       pyra_elmt3;
};

} // csmp

#endif
