// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_NT_lhsop_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_NT_lhsop_N_dV<dim,CELL>::Integral_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                          const char* oper,
                                                          const char* basic,
                                                          const char* test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test)
 {
    MathOperatorLHS<dim,CELL>::Name("Integral_NT_lhsop_N_dV", oper, basic, test );
    
        // testing the Operands 
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_lhsop_N_dV::(constructor)", 
                   oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_NT_lhsop_N_dV::(constructor)", 
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_lhsop_N_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }






/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_lhsop_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   scalar_value_ = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
}





/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_lhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // consistent formulation
    if ( !MathOperatorLHS<dim,CELL>::LumpedFormulation() )
      {
         // Integral of the testfunction products
         e.IntegralNN( MathOperatorLHS<dim,CELL>::LHS );
         MathOperatorLHS<dim,CELL>::LHS *= scalar_value_;
      }
    // lumped formulation  
    else
      {
         MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
         MathOperatorLHS<dim,CELL>::LHS.Zero();
         MathOperatorLHS<dim,CELL>::LHS(0,0) = (e.Volume() * scalar_value_) / static_cast<double>(e.Nodes());
         for ( uint32_t i{1U}; i<e.Nodes(); i++ )
           MathOperatorLHS<dim,CELL>::LHS(i,i) = MathOperatorLHS<dim,CELL>::LHS(0,0);
      }
   // MathOperatorLHS<dim,CELL>::LHS.Out();
      
} // end ComputeContribution


template class Integral_NT_lhsop_N_dV<1U,Element>;
template class Integral_NT_lhsop_N_dV<2U,Element>;
template class Integral_NT_lhsop_N_dV<3U,Element>;

template class Integral_NT_lhsop_N_dV<1U,Face>;
template class Integral_NT_lhsop_N_dV<2U,Face>;
template class Integral_NT_lhsop_N_dV<3U,Face>;

} // csmp
