// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_dNT_op_dV.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dV<dim,CELL>::NumIntegral_dNT_op_dV( const PropertyDatabase<dim>& pref,
                                                        const char*             oper,
                                                        const char*             test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    B_(dim,3), BT_(3,dim)
{
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_dNT_op_dV", oper, test );
    
    if ( !(MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
           MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
           MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == FACE ||
           MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == FACE_INTEGRATION_POINT ) ||
           MathOperatorRHS<dim,CELL>::MaterialOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dV<dim>::(constructor)", 
                      oper, "Operand must be a vector property placed on the element/face, or element/face integration point." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}




/**
      Computes:  f_j​ = i ∑​ w_i​ |∣J_i​ | k ∑​  ∂x k​ ∂N_j​ ​/ ∂x_k ​  f_k
      
            This is $\frac{\partial N_j}{\partial x_k}$ at each integration point.
 */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.UsesLocalCoordinates() == true );

    // initialize output matrix
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );
    
    //  When the material property is an element or face property
    // ----------------------------------------------------------
    const uint32_t n_nodes{ e.Nodes() };
    if ( this->MaterialOperandPlacement() == ELEMENT ||
         this->MaterialOperandPlacement() == FACE )
      {
        for ( uint32_t i{0U}; i<e.IntegrationPoints(); i++ )
          {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double detJ = e.dN_AtIntegrationPoint( B_, i, SCALAR );

             // transposing B -> BT 
             B_.Transposed( BT_ );

             // multiply  BT . MTRL
             BT_ *= this->MTRL[0];
                 
             // multiplying with determinant and weights
             BT_ *= e.WeightAtIntegrationPoint(i) * detJ;

             // row sum diagonalisation of matrix into right-hand vector
             for ( uint32_t k{0U}; k < dim; k++ )
               for ( uint32_t j{0U}; j<n_nodes; j++ )
                 MathOperatorRHS<dim,CELL>::RHS[j] += BT_(j,k);
          }
        return;
      }
      
    if ( this->MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
         this->MaterialOperandPlacement() == FACE_INTEGRATION_POINT )
      {
         for ( uint32_t i{0U}; i<e.IntegrationPoints(); i++ )
           {
              double detJ = e.dN_AtIntegrationPoint( B_, i, 1 );
              B_.Transposed( BT_ );
              BT_ *= MathOperatorRHS<dim,CELL>::MTRL[i];
              BT_ *= e.WeightAtIntegrationPoint(i) * detJ;
              for ( uint32_t k{0U}; k < dim; k++ )
                for ( uint32_t j{0U}; j<n_nodes; j++ )
                  MathOperatorRHS<dim,CELL>::RHS[j] += BT_(j,k);
           }
       }

} // end ComputeContribution


template class NumIntegral_dNT_op_dV<1U,Element>;
template class NumIntegral_dNT_op_dV<2U,Element>;
template class NumIntegral_dNT_op_dV<3U,Element>;

template class NumIntegral_dNT_op_dV<1U,Face>;
template class NumIntegral_dNT_op_dV<2U,Face>;
template class NumIntegral_dNT_op_dV<3U,Face>;

} // csmp











