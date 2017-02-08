#include "NumIntegral_dNT_op_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

/**
 
The Operand which is used here can be both, an element or a nodal variable
which is then interpolated to the integration points to obtain the 
integral properties.  
*/
template<size_t dim,class SIMPLEX>
NumIntegral_dNT_op_dN_dV<dim,SIMPLEX>::NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                                 const char*           oper,
                                                                 const char*           basic,
                                                                 const char*           test )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    B(dim,3), BT(3,dim)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_op_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}










/** Laplacian operator of shape function derivatives squared.
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_dNT_op_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    if ( this->MaterialOperandPlacement() == ELEMENT or
         this->MaterialOperandPlacement() == REGION or
         this->MaterialOperandPlacement() == FACE)
      {
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double64 detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );

             // transposing B -> BT  O.K.
             B.Transposed( BT );

             // multiply  BT . MTRL
             BT *= MathOperatorLHS<dim>::MTRL[0];

             // multiplying BT . B 
             BT *= B;
                 
             // multiplying with determinant and weights
             BT *= e.WeightAtIntegrationPoint(i) * detJ;

             // accumulating ME Gauss point integral contributions into element 
             // contribution to global conductance matrix
             MathOperatorLHS<dim>::LHS += BT;
          }
      }
    else { // NODE or ELEMENT_INTEGRATION_POINT
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) {
             double64 detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );
             B.Transposed( BT );
             BT *= MathOperatorLHS<dim>::MTRL[i];
             BT *= B;
             BT *= e.WeightAtIntegrationPoint(i) * detJ; 
             MathOperatorLHS<dim>::LHS += BT;
          }
     }

} // end ComputeContribution


//cout <<"\nNumIntegral_dNT_op_dN_dV: on Element "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( size_t i=0U; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nNumIntegral_dNT_op_dN_dV: zero element in diagonal of element matrix."; 


template class NumIntegral_dNT_op_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_op_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_op_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_op_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_op_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_op_dN_dV<3U,Face<3U> >;

template class NumIntegral_dNT_op_dN_dV<1U,InterFace<1U> >;
template class NumIntegral_dNT_op_dN_dV<2U,InterFace<2U> >;
template class NumIntegral_dNT_op_dN_dV<3U,InterFace<3U> >;

} // csmp











