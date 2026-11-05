// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_dNT_op_dV.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<uint32_t dim, template<uint32_t> class CELL>
Integral_dNT_op_dV<dim,CELL>::Integral_dNT_op_dV( const PropertyDatabase<dim>& pref,
                                                  const char*             oper,
                                                  const char*             test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test)
{
    MathOperatorRHS<dim,CELL>::Name("Integral_dNT_op_dV", oper, test );
    
    if ( !(MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
           MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
           MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == FACE) ||
           MathOperatorRHS<dim,CELL>::MaterialOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "Integral_dNT_op_dV<dim>::(constructor)", 
                      oper, "Operand must be a vector property placed on the element, face or element integration point." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_dNT_op_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}




template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_op_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.UsesLocalCoordinates() == false );

    // initialize output matrix
    const uint32_t n_nodes{ e.Nodes() };
    MathOperatorRHS<dim,CELL>::RHS.resize(n_nodes);
    
    //  Since the material property is an element or face property
    // -----------------------------------------------------------
    e.dN( DN_ );
    // return the transpose the shape function derivative matrix
    DN_.Transposed( DNT_ );
    
    // calculate the element contribution to RHS (lumping into vector format)
    DNT_ *= grad_prop_;

    // assigning element contribution & integrating the matrix
    const double volume = e.Volume();

    for ( uint32_t i{0U}; i<n_nodes; i++ )
      MathOperatorRHS<dim,CELL>::RHS[i] = DNT_(i,0) * volume;

} // end ComputeContribution


template class Integral_dNT_op_dV<1U,Element>;
template class Integral_dNT_op_dV<2U,Element>;
template class Integral_dNT_op_dV<3U,Element>;

template class Integral_dNT_op_dV<1U,Face>;
template class Integral_dNT_op_dV<2U,Face>;
template class Integral_dNT_op_dV<3U,Face>;

} // csmp











