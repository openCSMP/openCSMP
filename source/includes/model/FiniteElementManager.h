#ifndef CSMP_FINITE_ELEMENT_MANAGER_H
#define CSMP_FINITE_ELEMENT_MANAGER_H

#include <iostream>
#include <list>
#include "FiniteElement.h"

namespace csmp {

enum { ELEMENT_POLYTYPES=7 }; // maximum number of types connected in any one model

/**
@brief Container storing instances of the finite element types used in any particular instance of the class Model.

@author S.K. Matthai
@author Stephen G. Roberts
@author S. Geiger
@date 2001

@section motivation Motivation

@todo SKM: currently all the existing element types are hardwired into this class
and must always be present. Make this dynamic, loading elements only on demand.
 
*/
class FiniteElementManager {
  public:
    FiniteElementManager();
    FiniteElementManager( uint32_t dimensions, 
                          uint32_t interpolation_order,
                          bool isoparametric );

    FiniteElementManager( const FiniteElementManager& mgr );
    FiniteElementManager& operator=( const FiniteElementManager& mgr );
    ~FiniteElementManager();

    void InitializeElements( uint32_t dim, uint32_t interpolation_order, bool isoparametric );
    
    uint32_t          Dimensions() const;
    bool              ContainsElementType( CSMP_FEM_TYPE etype ) const;
    void              CurrentElementTypes( std::list<CSMP_FEM_TYPE>& etypes ) const;
    FiniteElement*    E( int8_t csmp_etype ) const;
    FiniteElement*    E( CSMP_FEM_TYPE csmp_etype ) const;
    uint32_t          NodesOfElementType( CSMP_FEM_TYPE etype ) const;
    // standard types
    FiniteElement*    LinearBarElement() const;
    FiniteElement*    LinearTriangleElement() const;
    FiniteElement*    LinearTetrahedronElement() const;
    
    // can be changed, but intermediate nodes will become defunct
    void              InterpolationOrder( uint32_t interpolation_order ); // TODO: does not allocate new suite of elements!
    uint32_t          InterpolationOrder() const;
    bool              UsesElementsWithLocalCoordinateSystem() const;
    
    void              Out() const;
  
  private:
    bool      mixed_element_formulation; ///< different-order elements at same time
    uint32_t  dimensions, interpolation;
    // volume elements
    FiniteElement*  hexa_ptr, *pyra_ptr, *pris_ptr, *tetr_ptr;
    // surface elements
    FiniteElement*  quad_ptr, *tria_ptr;
    // line elements
    FiniteElement*  line_ptr;
};

} // csmp

#endif
