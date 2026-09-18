// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_dNT_rhsop_dN_v_dV.h"

#include "PropertyDatabase.h"
#include "Exception.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_dNT_rhsop_dN_v_dV<dim,CELL>::Integral_dNT_rhsop_dN_v_dV( const PropertyDatabase<dim>& pref,
                                                              const char*             oper,
                                                              const char*             test,
                                                              const char*             grad_var )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    grad_key(pref.StorageKey( grad_var )),
    DN(2,3), DNT(3,2), VAR(3,1)
{
    MathOperatorRHS<dim,CELL>::Name("Integral_dNT_rhsop_dN_v_dV", oper, test );

    // testing the Operands 
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_v_dV::(constructor)", 
                    oper, "Operand must be placed on the element or group.");

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_v_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_v_dV::(constructor)", 
                    grad_var, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}








template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_rhsop_dN_v_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    // only if the property is an element property  something is done here
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == REGION )
      {
         MathOperatorRHS<dim,CELL>::MTRL[0].Resize(dim,dim);
         MathOperatorRHS<dim,CELL>::MTRL[0].Zero();
      
         if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
              MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonalAndZeroOffDiagonal( dim, e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey() ) ); 
           }
         if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), vc );
              MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonal( vc ); 
           }
         if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), ts );
              MathOperatorRHS<dim,CELL>::MTRL[0] = ts;
           }
      }
    else throw csmp::Exception( FATAL_ERROR, "Integral_dNT_rhsop_dN_v_dV::GetOperands", 
                               "The current finite element is analytically integrated and therefore has no integration points",
                                        "Hence, nodal properties cannot be integrated. Use numerically integrated elements in stead.");
   VAR.Resize(e.Nodes(),1);   
   // the gradient variable
   for ( uint32_t i{0}; i<e.Nodes(); i++ )
     VAR(i,0) = e.N(i)->Read( grad_key );   
          
} // end GetOperands




template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_rhsop_dN_v_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // calculate the element contribution to RHS (lumping into vector format)
    DNT *= MathOperatorRHS<dim,CELL>::MTRL[0];
    DNT *= DN;
    DNT *= VAR;
    
    // assigning element contribution & integrating the matrix
    double volume = e.Volume();
    
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    for ( uint32_t i{0U}; i<e.Nodes(); i++ )
      MathOperatorRHS<dim,CELL>::RHS[i] = DNT(i,0) * volume;

} // end ComputeContribution



template class Integral_dNT_rhsop_dN_v_dV<1U,Element>;
template class Integral_dNT_rhsop_dN_v_dV<2U,Element>;
template class Integral_dNT_rhsop_dN_v_dV<3U,Element>;

template class Integral_dNT_rhsop_dN_v_dV<1U,Face>;
template class Integral_dNT_rhsop_dN_v_dV<2U,Face>;
template class Integral_dNT_rhsop_dN_v_dV<3U,Face>;

} // csmp
