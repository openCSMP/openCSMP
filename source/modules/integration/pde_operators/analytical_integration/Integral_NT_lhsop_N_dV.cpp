#include "Integral_NT_lhsop_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
Integral_NT_lhsop_N_dV<dim,SIMPLEX>::Integral_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                        const char* oper,
                                                        const char* basic, 
                                                        const char* test )
  : MathOperatorLHS<dim>(pref,oper,basic,test)
 {
    MathOperatorLHS<dim>::Name("Integral_NT_lhsop_N_dV", oper, basic, test );
    
        // testing the Operands 
    if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Integral_NT_lhsop_N_dV::(constructor)", 
                   oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( CSMP_ERROR, "Integral_NT_lhsop_N_dV::(constructor)", 
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Integral_NT_lhsop_N_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }






/** Reads the Operand values from the elements.
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_lhsop_N_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
   scalar_value_ = e.Read( MathOperatorLHS<dim>::MaterialOperandKey() );
}


/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_lhsop_N_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    // consistent formulation
    if ( !MathOperatorLHS<dim>::LumpedFormulation() )
      {
         // Integral of the testfunction products
         e.IntegralNN( MathOperatorLHS<dim>::LHS );
         MathOperatorLHS<dim>::LHS *= scalar_value_;         
      }
    // lumped formulation  
    else
      {
         MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
         MathOperatorLHS<dim>::LHS.Zero();
         MathOperatorLHS<dim>::LHS(0,0) = (e.Volume() * scalar_value_) / static_cast<double64>(e.Nodes());
         for ( size_t i=1; i<e.Nodes(); i++ ) 
           MathOperatorLHS<dim>::LHS(i,i) = MathOperatorLHS<dim>::LHS(0,0); 
      }
   // MathOperatorLHS<dim>::LHS.Out();     
      
} // end ComputeContribution


template class Integral_NT_lhsop_N_dV<1U,Element<1U> >;
template class Integral_NT_lhsop_N_dV<2U,Element<2U> >;
template class Integral_NT_lhsop_N_dV<3U,Element<3U> >;

template class Integral_NT_lhsop_N_dV<1U,Face<1U> >;
template class Integral_NT_lhsop_N_dV<2U,Face<2U> >;
template class Integral_NT_lhsop_N_dV<3U,Face<3U> >;

} // csmp
