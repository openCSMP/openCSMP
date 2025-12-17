#include "CVFEM_NumIntegral_dNT_op_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_NumIntegral_dNT_op_dN_dV<dim,CELL>::CVFEM_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                                          const char*           oper,
                                                                          const char*           basic,
                                                                          const char*           test )
  : CVFEM_MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    B(dim,3), BT(3,dim)
{
    MathOperatorLHS<dim,CELL>::Name("CVFEM_NumIntegral_dNT_op_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      basic, "Weighting variable must be a scalar property placed on the nodes." );


}

/** compute contribution */
template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_NumIntegral_dNT_op_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    if ( this->MaterialOperandPlacement() == ELEMENT or
         this->MaterialOperandPlacement() == REGION or this->MaterialOperandPlacement() == FACE)
      {
        for ( auto i{0}; i<e.FE()->IntegrationPoints(); i++ ) {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );

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
             MathOperatorLHS<dim,CELL>::LHS += BT;
          }
      }
    else { // NODE or ELEMENT_INTEGRATION_POINT
        for ( auto i{0}; i<e.FE()->IntegrationPoints(); i++ ) {
             double detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );
             B.Transposed( BT );
             BT *= MathOperatorLHS<dim,CELL>::MTRL[i];
             BT *= B;
             BT *= e.WeightAtIntegrationPoint(i) * detJ; 
             MathOperatorLHS<dim,CELL>::LHS += BT;
          }
     }

} // end ComputeContribution


template class CVFEM_NumIntegral_dNT_op_dN_dV<1U,Element>;
template class CVFEM_NumIntegral_dNT_op_dN_dV<2U,Element>;
template class CVFEM_NumIntegral_dNT_op_dN_dV<3U,Element>;

template class CVFEM_NumIntegral_dNT_op_dN_dV<1U,Face>;
template class CVFEM_NumIntegral_dNT_op_dN_dV<2U,Face>;
template class CVFEM_NumIntegral_dNT_op_dN_dV<3U,Face>;

} // csmp











