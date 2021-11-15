#include "NumIntegral_dNT_dN_dV.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

/** Laplacian build with the interpolation function derivatives.
*/
template<size_t dim,class CELL>
NumIntegral_dNT_dN_dV<dim,CELL>::NumIntegral_dNT_dN_dV( const PropertyDatabase<dim>& pref,
                                                      const char*           basic, 
                                                      const char*           test ) 
  : MathOperatorLHS<dim>(pref,basic,test),
    B(dim,3), BT(3,dim)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_dN_dV", basic, test );
    
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}










/** Laplacian operator of shape function derivatives squared.
 */
template<size_t dim,class CELL>
void NumIntegral_dNT_dN_dV<dim,CELL>::ComputeContribution( CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         double detJ = e.dN_AtIntegrationPoint( B, i, 1 );

         // transposing B -> BT  O.K.
         B.Transposed( BT );

         // multiplying BT . B 
         BT *= B;
             
         // multiplying with determinant and weights
         BT *= e.WeightAtIntegrationPoint(i) * detJ; 

         // accumulating ME Gauss point integral contributions into element 
         // contribution to global conductance matrix
         MathOperatorLHS<dim>::LHS += BT;
      }

} // end ComputeContribution



template class NumIntegral_dNT_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_dN_dV<3U,Face<3U> >;

} // csmp











