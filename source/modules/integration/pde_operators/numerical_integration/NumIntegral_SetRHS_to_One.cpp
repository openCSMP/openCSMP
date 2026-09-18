// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_SetRHS_to_One.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_SetRHS_to_One<dim,CELL>::~NumIntegral_SetRHS_to_One() {}




template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_SetRHS_to_One<dim,CELL>::NumIntegral_SetRHS_to_One( const PropertyDatabase<dim>& pref, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,test)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_SetRHS_to_One", test );
 
    if ( MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR ||
         MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE )
      throw csmp::Exception( ERROR, "MathOperatorRHS->NumIntegral_SetRHS_to_One<dim>::(constructor):",
                             test, "Dependent-variable must be a scalar variable placed on the nodes." );
                      
 } // end constructor





/**

Generats a RHS vector of ones that is added to the global solution matrix.

The result is returned into the MathOperatorRHS vector<fT> 'rhs'.

@attention this 'rhs' contribution is not equivalent to an integral over an operand with a value of 1.

tested: O.K.  */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_SetRHS_to_One<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
   // this integral is only for numerically integrated isoparametric finite elements
   assert( e.FE()->Isoparametric() == true );

   // create a RHS vector of zeros
   MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
   for ( uint32_t i{0U}; i<e.Nodes(); i++ ) MathOperatorRHS<dim,CELL>::RHS[i] = 1.0/static_cast<double>(e.N(i)->Parents());

} // end ComputeContribution


template class NumIntegral_SetRHS_to_One<1U,Element>;
template class NumIntegral_SetRHS_to_One<2U,Element>;
template class NumIntegral_SetRHS_to_One<3U,Element>;

template class NumIntegral_SetRHS_to_One<1U,Face>;
template class NumIntegral_SetRHS_to_One<2U,Face>;
template class NumIntegral_SetRHS_to_One<3U,Face>;

} // csmp
