// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_SetRHS_to_Zero.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_SetRHS_to_Zero<dim,CELL>::Integral_SetRHS_to_Zero( const PropertyDatabase<dim>& pref, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,test)
 {
    MathOperatorRHS<dim,CELL>::Name("Integral_SetRHS_to_Zero", "0.", test );
    
        // testing the Operands 
    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_SetRHS_to_Zero::(constructor)", 
                             test, "Dependent variable must be a scalar property placed on the nodes." );
 }






/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_SetRHS_to_Zero<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );
      
} // end ComputeContribution


template class Integral_SetRHS_to_Zero<1U,Element>;
template class Integral_SetRHS_to_Zero<2U,Element>;
template class Integral_SetRHS_to_Zero<3U,Element>;

template class Integral_SetRHS_to_Zero<1U,Face>;
template class Integral_SetRHS_to_Zero<2U,Face>;
template class Integral_SetRHS_to_Zero<3U,Face>;

} // csmp
