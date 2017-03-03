#include "Integral_dNT_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
Integral_dNT_dN_dV<dim,SIMPLEX>::Integral_dNT_dN_dV( const PropertyDatabase<dim>& pref,
                                                     const char*             oper,
                                                     const char*             test )
  : MathOperatorRHS<dim>(pref,oper,test),
    DN(2,3), DNT(3,2), UNITY(3,1)
{
    MathOperatorRHS<dim>::Name("Integral_dNT_dN_dV", oper, test );

    // testing the Operands 
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != ELEMENT or MathOperatorRHS<dim>::MaterialOperandPlacement() != REGION )
      throw csmp::Exception( ERROR, "Integral_dNT_dN_dV::(constructor)", 
                    oper, "Operand must be placed on the element or group.");

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_dNT_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}



template<size_t dim,class SIMPLEX>
void Integral_dNT_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // setting up the unity vector
    UNITY.Resize(e.Nodes(),1);
    UNITY = 1.;

    // calculate the element contribution to RHS (lumping into vector format)
    DNT *= MathOperatorRHS<dim>::MTRL[0];
    DNT *= DN;
    DNT *= UNITY;
    
    // assigning element contribution & integrating the matrix
    double64 volume = e.Volume();
    
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    for ( size_t i=0; i<e.Nodes(); i++ )
      MathOperatorRHS<dim>::RHS[i] = DNT(i,0) * volume;

} // end ComputeContribution



template class Integral_dNT_dN_dV<1U,Element<1U> >;
template class Integral_dNT_dN_dV<2U,Element<2U> >;
template class Integral_dNT_dN_dV<3U,Element<3U> >;

template class Integral_dNT_dN_dV<1U,Face<1U> >;
template class Integral_dNT_dN_dV<2U,Face<2U> >;
template class Integral_dNT_dN_dV<3U,Face<3U> >;

} // csmp
