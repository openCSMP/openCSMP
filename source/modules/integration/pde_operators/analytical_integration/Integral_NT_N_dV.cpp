#include "Integral_NT_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_NT_N_dV<dim,CELL>::Integral_NT_N_dV( const PropertyDatabase<dim>& pref,
                                              const char* test )
  : MathOperatorLHS<dim,CELL>(pref,"permeability",test,test)
 {
    MathOperatorLHS<dim,CELL>::Name("Integral_NT_N_dV", "__", test, test );
    
        // testing the Operands 
    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_NT_N_dV::(constructor)", 
                             test, "Dependent variable must be a scalar property placed on the nodes." );
 }






/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    // consistent formulation
    if ( !MathOperatorLHS<dim,CELL>::LumpedFormulation() )
      {
         // Integral of the testfunction products
         e.IntegralNN( MathOperatorLHS<dim,CELL>::LHS );
      }
    // lumped formulation  
    else
      {
         MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
         double vol_div3 = e.Volume() / static_cast<double>(e.Nodes());
         for ( auto i{0U}; i<e.Nodes(); i++ )
           for ( auto j{0U}; j<e.Nodes(); j++ )
             if ( i == j ) MathOperatorLHS<dim,CELL>::LHS(i,j) = vol_div3;
             else          MathOperatorLHS<dim,CELL>::LHS(i,j) = 0.;
      }      
      
} // end ComputeContribution


template class Integral_NT_N_dV<1U,Element>;
template class Integral_NT_N_dV<2U,Element>;
template class Integral_NT_N_dV<3U,Element>;

template class Integral_NT_N_dV<1U,Face>;
template class Integral_NT_N_dV<2U,Face>;
template class Integral_NT_N_dV<3U,Face>;

} // csmp
