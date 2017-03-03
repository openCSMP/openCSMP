#include "NumIntegral_dNT_op_dV.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
NumIntegral_dNT_op_dV<dim,SIMPLEX>::NumIntegral_dNT_op_dV( const PropertyDatabase<dim>& pref,
                                                           const char*             oper,
                                                           const char*             test )
  : MathOperatorRHS<dim>(pref,oper,test),
    B(dim,3), BT(3,dim)
{
    MathOperatorRHS<dim>::Name("NumIntegral_dNT_op_dV", oper, test );
    
    if ( !(MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ||
           MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
           MathOperatorRHS<dim>::MaterialOperandPlacement() == FACE) ||
           MathOperatorRHS<dim>::MaterialOperandType() != VECTOR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_dNT_op_dV<dim>::(constructor)", 
                      oper, "Operand must be a vector property placed on the element, face or element integration point." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_dNT_op_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}



template<size_t dim,class SIMPLEX>
void NumIntegral_dNT_op_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    // initialize output matrix
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    //  There is only one case, where matrial property is an element property
    // ------------------------------------------------------------------
    if ( this->MaterialOperandPlacement() == ELEMENT || this->MaterialOperandPlacement() == FACE )
      {
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) 
          {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double64 detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );

             // transposing B -> BT 
             B.Transposed( BT );

             // multiply  BT . MTRL
             BT *= this->MTRL[0];
                 
             // multiplying with determinant and weights
             BT *= e.WeightAtIntegrationPoint(i) * detJ; 

             // accumulating ME Gauss point integral contributions into element 
             // contribution to global conductance matrix
             
             for ( size_t k=0; k < dim; k++ )
                 for ( size_t j=0; j<e.Nodes(); j++ ) MathOperatorRHS<dim>::RHS[j] += BT(j,k);
          }
      
      }
    if ( this->MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
      {
         for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ )
           {
              double64 detJ = e.dN_AtIntegrationPoint( B, i, 1 );
              B.Transposed( BT );
              BT *= MathOperatorRHS<dim>::MTRL[i];
              BT *= e.WeightAtIntegrationPoint(i) * detJ;
              for ( size_t k=0; k < dim; k++ )
                for ( size_t j=0; j<e.Nodes(); j++ ) MathOperatorRHS<dim>::RHS[j] += BT(j,k);
           }
       }

} // end ComputeContribution


template class NumIntegral_dNT_op_dV<1U,Element<1U> >;
template class NumIntegral_dNT_op_dV<2U,Element<2U> >;
template class NumIntegral_dNT_op_dV<3U,Element<3U> >;

template class NumIntegral_dNT_op_dV<1U,Face<1U> >;
template class NumIntegral_dNT_op_dV<2U,Face<2U> >;
template class NumIntegral_dNT_op_dV<3U,Face<3U> >;

} // csmp











