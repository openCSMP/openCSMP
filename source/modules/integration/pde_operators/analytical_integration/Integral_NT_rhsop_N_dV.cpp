#include "Integral_NT_rhsop_N_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


// * basic and test variable exchanged so that the line is indicated by the basic variable (as it is
// * for the lhs operators)
template<size_t dim,class SIMPLEX>
Integral_NT_rhsop_N_dV<dim,SIMPLEX>::Integral_NT_rhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                             const char* oper,
                                                             const char* basic,
                                                             const char* test )
  : MathOperatorRHS<dim>(pref,oper,basic),
    INN(3,3),
    basic_(pref.Parameter(test))
 {
    string name = "Integral_NT_rhsop_N_dV ";
    name += test;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim>::Name(cname, oper, basic );
    MathOperatorRHS<dim>::AddAccumulateLater();

    // testing the Operands 
    if ( basic_.key.place != NODE || basic_.key.type != SCALAR )
        throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                   basic, "Basic variable must be a scalar property placed on the node." );

    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != ELEMENT and MathOperatorRHS<dim>::MaterialOperandPlacement() != REGION )
        throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "Integral_NT_rhsop_N_dV<dim>::(constructor)",
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }


/** Reads the Operand values from the elements.
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_rhsop_N_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   // reading oper and basic
   e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), sc );
   e.NodePropertyVector( basic_.key, basic_var_ );
}


/**
 
Computes the volume (area) integral over the basic function products 
multiplied with the Operand and stores the result in the test function
part of the right hand side vector.  
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_rhsop_N_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());

    // consistent formulation
    if ( !MathOperatorRHS<dim>::LumpedFormulation() ) {
      // Integral of the testfunction products
      e.IntegralNN( INN );
      fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0.0 );
      
      // the matrix is contracted into a vector by multiplying with the basis vector
      for ( size_t i=0; i<e.Nodes(); i++ ) 
        for ( size_t j=0; j<e.Nodes(); j++ ) 
          MathOperatorRHS<dim>::RHS[i] += INN(i,j) * basic_var_[j]() * sc();
    }
    // lumped formulation  
    else {
      double64 res = (e.Volume() * sc()) / static_cast<double64>(e.Nodes());
      for (size_t i = 0; i < e.Nodes(); ++i) {
        MathOperatorRHS<dim>::RHS[i] = basic_var_[i]() * res;
      }
    }

} // end ComputeContribution 

template class Integral_NT_rhsop_N_dV<2U,Element<2U> >;
template class Integral_NT_rhsop_N_dV<3U,Element<3U> >;

template class Integral_NT_rhsop_N_dV<2U,Face<2U> >;
template class Integral_NT_rhsop_N_dV<3U,Face<3U> >;

} // csmp
