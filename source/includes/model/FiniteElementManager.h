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
    FiniteElementManager( size_t dimensions, 
                          size_t interpolation_order,
                          bool isoparametric );

    FiniteElementManager( const FiniteElementManager& mgr );
    FiniteElementManager& operator=( const FiniteElementManager& mgr );
    ~FiniteElementManager();

    void InitializeElements( size_t dim, size_t interpolation_order, bool isoparametric );
    
    size_t            Dimensions() const;
    bool              ContainsElementType( CSMP_FEM_TYPE etype ) const;
    void              CurrentElementTypes( std::list<CSMP_FEM_TYPE>& etypes ) const;
    FiniteElement*    E( int32 csmp_etype ) const;
    FiniteElement*    E( CSMP_FEM_TYPE csmp_etype ) const;
    size_t            NodesOfElementType( CSMP_FEM_TYPE etype ) const;
    // standard types
    FiniteElement*    LinearBarElement() const;
    FiniteElement*    LinearTriangleElement() const;
    FiniteElement*    LinearTetrahedronElement() const;
    
    // can be changed, but intermediate nodes will become defunct
    void              InterpolationOrder( size_t interpolation_order );
    size_t            InterpolationOrder() const;
    
    void              Out() const;
  
  private:
    bool    mixed_element_formulation; // different-order elements at same time
    size_t  dimensions, interpolation;
    // volume elements
    FiniteElement*  hexa_ptr, *pyra_ptr, *pris_ptr, *tetr_ptr;
    // surface elements
    FiniteElement*  quad_ptr, *tria_ptr;
    // line elements
    FiniteElement*  line_ptr;
};



inline  FiniteElement*  FiniteElementManager::LinearBarElement() const
 {
    if ( line_ptr != NULL ) return line_ptr;
    std::cerr <<"\nFiniteElementManager::LinearBarElement: Not available."<< std::endl;
    return NULL;
 }
 
 
inline  FiniteElement*  FiniteElementManager::LinearTriangleElement() const
 {
    if ( tria_ptr != NULL ) return tria_ptr;
    std::cerr <<"\nFiniteElementManager::LinearTriangleElement: Not available."<< std::endl;
    return NULL;
 }
 
 
inline  FiniteElement*  FiniteElementManager::LinearTetrahedronElement() const
 {
    if ( tetr_ptr != NULL ) return tetr_ptr;
    std::cerr <<"\nFiniteElementManager::LinearTetrahedronElement: Not available."<< std::endl;
    return NULL;
 }


inline  FiniteElement*  FiniteElementManager::E( int32 e_type ) const
 {
    if ( dimensions == 3 ) {
	     if ( hexa_ptr != NULL && e_type == hexa_ptr->ElementType() ) return hexa_ptr;
	     if ( pyra_ptr != NULL && e_type == pyra_ptr->ElementType() ) return pyra_ptr;
	     if ( pris_ptr != NULL && e_type == pris_ptr->ElementType() ) return pris_ptr;
	     if ( tetr_ptr != NULL && e_type == tetr_ptr->ElementType() ) return tetr_ptr;
      }
    if ( quad_ptr != NULL && e_type == quad_ptr->ElementType() ) return quad_ptr;
    if ( tria_ptr != NULL && e_type == tria_ptr->ElementType() ) return tria_ptr;
    if ( line_ptr != NULL && e_type == line_ptr->ElementType() ) return line_ptr;

    std::cerr <<"\nFiniteElementManager::Element: Requested element is not available: ";
    std::cerr << e_type <<" = "<< parseFiniteElementType(e_type) <<" returning NULL pointer."<< std::endl;
      
    return NULL;
 }


inline  FiniteElement*  FiniteElementManager::E( CSMP_FEM_TYPE e_type ) const
 {
    if ( dimensions == 3 ) {
	     if ( hexa_ptr != NULL && e_type == hexa_ptr->ElementType() ) return hexa_ptr;
	     if ( pyra_ptr != NULL && e_type == pyra_ptr->ElementType() ) return pyra_ptr;
	     if ( pris_ptr != NULL && e_type == pris_ptr->ElementType() ) return pris_ptr;
	     if ( tetr_ptr != NULL && e_type == tetr_ptr->ElementType() ) return tetr_ptr;
      }
    if ( quad_ptr != NULL && e_type == quad_ptr->ElementType() ) return quad_ptr;
    if ( tria_ptr != NULL && e_type == tria_ptr->ElementType() ) return tria_ptr;
    if ( line_ptr != NULL && e_type == line_ptr->ElementType() ) return line_ptr;

    std::cerr <<"\nFiniteElementManager::Element: Requested element is not available: ";
    std::cerr << e_type <<" = "<< parseFiniteElementType(e_type) <<" return NULL pointer."<< std::endl;
      
    return NULL;
 }


inline  size_t  FiniteElementManager::NodesOfElementType( CSMP_FEM_TYPE etype ) const
 {
    return E( etype )->Nodes();
 }



inline  size_t  FiniteElementManager::Dimensions() const
 { return dimensions; }

} // csmp

#endif
