#include "Integral_NT_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {


template<uint32_t dim,class CELL>
Integral_NT_N_dV<dim,CELL>::Integral_NT_N_dV( const PropertyDatabase<dim>& pref,
                                                 const char* test )
  : MathOperatorLHS<dim>(pref,"permeability",test,test)
 {
    MathOperatorLHS<dim>::Name("Integral_NT_N_dV", "__", test, test );
    
        // testing the Operands 
    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_NT_N_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }






/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim,class CELL>
void Integral_NT_N_dV<dim,CELL>::ComputeContribution( const CELL& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    // consistent formulation
    if ( !MathOperatorLHS<dim>::LumpedFormulation() )
      {
         // Integral of the testfunction products
         e.IntegralNN( MathOperatorLHS<dim>::LHS );
      }
    // lumped formulation  
    else
      {
         MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
         double vol_div3 = e.Volume() / static_cast<double>(e.Nodes());
         for ( auto i=0; i<e.Nodes(); i++ )
           for ( auto j=0; j<e.Nodes(); j++ )
             if ( i == j ) MathOperatorLHS<dim>::LHS(i,j) = vol_div3;
             else          MathOperatorLHS<dim>::LHS(i,j) = 0.;
      }      
      
} // end ComputeContribution


template class Integral_NT_N_dV<1U,Element<1U> >;
template class Integral_NT_N_dV<2U,Element<2U> >;
template class Integral_NT_N_dV<3U,Element<3U> >;

template class Integral_NT_N_dV<1U,Face<1U> >;
template class Integral_NT_N_dV<2U,Face<2U> >;
template class Integral_NT_N_dV<3U,Face<3U> >;

} // csmp
