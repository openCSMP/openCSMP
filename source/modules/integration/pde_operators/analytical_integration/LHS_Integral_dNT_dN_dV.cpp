#include "LHS_Integral_dNT_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
LHS_Integral_dNT_dN_dV<dim,SIMPLEX>::LHS_Integral_dNT_dN_dV( const PropertyDatabase<dim>& pref,
                                                             const char*          basic,
                                                             const char*          test )
  : MathOperatorLHS<dim>(pref,basic,test),
    DN(2,3), DNT(3,2),
    unity(dim,1.)
{
    MathOperatorLHS<dim>::Name("LHS_Integral_dNT_dN_dV", basic, test );

    // testing the Operands 
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "LHS_Integral_dNT_dN_dV::(constructor)", 
                    test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "LHS_Integral_dNT_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}






/**
   Computes the product of the interpolation function derivatives.
*/
template<size_t dim,class SIMPLEX>
void LHS_Integral_dNT_dN_dV<dim,SIMPLEX>::ComputeContribution( const SIMPLEX& e )
 {
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );
   
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );

    // calculate the element contribution to LHS (E_OP is the Basic Operand)
    DNT *= unity;
    DNT *= DN;
    // assigning element contribution
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS = DNT;
    // analytical integration over area / volume for linear triangle 
    // and tetrahedron elements, respectively
    MathOperatorLHS<dim>::LHS *= e.Volume();

} // end ComputeContribution



template class LHS_Integral_dNT_dN_dV<1U,Element<1U> >;
template class LHS_Integral_dNT_dN_dV<2U,Element<2U> >;
template class LHS_Integral_dNT_dN_dV<3U,Element<3U> >;

template class LHS_Integral_dNT_dN_dV<1U,Face<1U> >;
template class LHS_Integral_dNT_dN_dV<2U,Face<2U> >;
template class LHS_Integral_dNT_dN_dV<3U,Face<3U> >;

} // csmp
