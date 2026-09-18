// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_NT_rhsop_N_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


/** basic and test variable exchanged so that the line is indicated by the basic variable (as it is
    for the lhs operators)
*/
template<uint32_t dim, template<uint32_t> class CELL>
Integral_NT_rhsop_N_dV<dim,CELL>::Integral_NT_rhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                          const char* oper,
                                                          const char* basic,
                                                          const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,basic),
    INN(3,3),
    basic_(pref.Parameter(test))
 {
    string name = "Integral_NT_rhsop_N_dV ";
    name += test;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim,CELL>::Name(cname, oper, basic );
    MathOperatorRHS<dim,CELL>::AddAccumulateLater();

    // testing the Operands 
    if ( basic_.key.place != NODE || basic_.key.type != SCALAR )
        throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                   basic, "Basic variable must be a scalar property placed on the node." );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != REGION )
        throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }






/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_rhsop_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   // reading oper and basic
   e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), sc );
   e.NodePropertyVector( basic_.key, basic_var_ );
}






/**
 
Computes the volume (area) integral over the basic function products 
multiplied with the Operand and stores the result in the test function
part of the right hand side vector.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_rhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());

    // consistent formulation
    if ( !MathOperatorRHS<dim,CELL>::LumpedFormulation() ) {
      // Integral of the testfunction products
      e.IntegralNN( INN );
      fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );
      
      // the matrix is contracted into a vector by multiplying with the basis vector
      for ( auto i{0U}; i<e.Nodes(); i++ ) 
        for ( auto j{0U}; j<e.Nodes(); j++ ) 
          MathOperatorRHS<dim,CELL>::RHS[i] += INN(i,j) * basic_var_[j]() * sc();
    }
    // lumped formulation  
    else {
      double res = (e.Volume() * sc()) / static_cast<double>(e.Nodes());
      for (auto i = 0; i < e.Nodes(); ++i) {
        MathOperatorRHS<dim,CELL>::RHS[i] = basic_var_[i]() * res;
      }
    }

} // end ComputeContribution 

template class Integral_NT_rhsop_N_dV<2U,Element>;
template class Integral_NT_rhsop_N_dV<3U,Element>;

template class Integral_NT_rhsop_N_dV<2U,Face>;
template class Integral_NT_rhsop_N_dV<3U,Face>;

} // csmp
