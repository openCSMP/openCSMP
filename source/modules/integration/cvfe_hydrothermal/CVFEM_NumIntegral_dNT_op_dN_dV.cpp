#include "CVFEM_NumIntegral_dNT_op_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

/** default destructor */
template<size_t dim, class SIMPLEX>
CVFEM_NumIntegral_dNT_op_dN_dV<dim,SIMPLEX>::~CVFEM_NumIntegral_dNT_op_dN_dV() {}

/** custom constructor */
template<size_t dim, class SIMPLEX>
CVFEM_NumIntegral_dNT_op_dN_dV<dim,SIMPLEX>::CVFEM_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                                                            const char*           oper, 
                                                            const char*           basic, 
                                                            const char*           test ) 
  : CVFEM_MathOperatorLHS<dim>(pref,oper,basic,test),
    B(dim,3), BT(3,dim)
{
    MathOperatorLHS<dim>::Name("CVFEM_NumIntegral_dNT_op_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      basic, "Weighting variable must be a scalar property placed on the nodes." );


}

/** compute contribution */
template<size_t dim, class SIMPLEX>
void CVFEM_NumIntegral_dNT_op_dN_dV<dim,SIMPLEX>::ComputeContribution( Element<dim>& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    if ( this->MaterialOperandPlacement() == ELEMENT or this->MaterialOperandPlacement() == REGION or this->MaterialOperandPlacement() == FACE)
      {
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double64 detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );

             // transposing B -> BT  O.K.
             B.Transposed( BT );

             // multiply  BT . MTRL
             BT *= this->MTRL[0];

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


template class CVFEM_NumIntegral_dNT_op_dN_dV<1U,Element<1U> >;
template class CVFEM_NumIntegral_dNT_op_dN_dV<2U,Element<2U> >;
template class CVFEM_NumIntegral_dNT_op_dN_dV<3U,Element<3U> >;

template class CVFEM_NumIntegral_dNT_op_dN_dV<1U,Face<1U> >;
template class CVFEM_NumIntegral_dNT_op_dN_dV<2U,Face<2U> >;
template class CVFEM_NumIntegral_dNT_op_dN_dV<3U,Face<3U> >;

} // csmp











