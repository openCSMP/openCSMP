// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  NumJacobianIntegral_BT_mN_dV.cpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 22/6/2026.
//

#include "NumJacobianIntegral_BT_m_N_dV.h"
#include "Element.h"
#include "Face.h"

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumJacobianIntegral_BT_m_N_dV<dim,CELL>::NumJacobianIntegral_BT_m_N_dV( const PropertyDatabase<dim>& pref,
                                                                        const char* oper,
                                                                        const char* basic,
                                                                        const char* test )
  : MathOperatorLHS<dim,CELL>(pref, oper, basic, test),
    alpha_(1U)
{
    MathOperatorLHS<dim,CELL>::Name("NumJacobianIntegral_BT_m_Np_dV", oper, basic, test);

    // --- Validate Biot coefficient ---
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumJacobianIntegral_BT_m_N_dV::(constructor)",
                        oper, "Biot coefficient must be a scalar property." );

    // --- Validate displacement (basic) ---
    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType()      != VECTOR )
        throw csmp::Exception( ERROR, "NumJacobianIntegral_BT_m_N_dV::(constructor)",
                        test, "operand must be a vector property discretised on the nodes." );

    // --- Validate pressure (test) ---
    // NOTE: test is scalar (pressure), unlike the stiffness matrix
    // where both basic and test are vector (displacement).
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType()      != SCALAR )
        throw csmp::Exception( ERROR, "NumJacobianIntegral_BT_m_N_dV::(constructor)",
                        basic, "operand must be a scalar property discretised on the nodes." );
}


// -----------------------------------------------------------------------
// GetOperands — mirrors your existing pattern exactly
// -----------------------------------------------------------------------
template<uint32_t dim, template<uint32_t> class CELL>
void NumJacobianIntegral_BT_m_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    assert( e.FE()->Isoparametric() == true );

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
    {
        alpha_[0] = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
    }
    else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
    {
        const size_t ipoints( e.IntegrationPoints() );
        alpha_.resize( ipoints );
        for ( uint32_t i{0}; i < ipoints; ++i )
            alpha_[i] = e.Read( i, MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
    }
    else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE )
    {
        // Interpolate nodal alpha to integration points
        const size_t ipoints( e.IntegrationPoints() );
        alpha_.resize( ipoints );
        for ( uint32_t i{0}; i < ipoints; ++i )
        {
           alpha_[i] = e.PropertyValueAtIntegrationPoint( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), i );
        }
    }
    else
        throw csmp::Exception( WARNING, "NumJacobianIntegral_BT_m_Np_dV::GetOperands",
                        MathOperatorLHS<dim,CELL>::MaterialOperandName().c_str(),
                        "Placement of Biot coefficient not handled." );
}



// -----------------------------------------------------------------------
// ComputeContribution
//
// Computes:  Q^e = sum_i [ alpha_i * B^T(xi) * m * Np(xi) * w_i * detJ_i ]
//
// Output matrix is (dim*nodes x nodes):
//   rows -> displacement DOFs (interleaved)
//   cols -> pressure DOFs
// -----------------------------------------------------------------------
template<uint32_t dim, template<uint32_t> class CELL>
void NumJacobianIntegral_BT_m_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    const uint32_t n_nodes = e.Nodes();
    const uint32_t BT_dof  = dim * n_nodes;   // rows: displacement DOFs
    const uint32_t N_dof   = n_nodes;         // cols: pressure DOFs

    // Output: (BT_dof x N_dof)
    // Fixed: Replaced undeclared disp_dof/pres_dof with BT_dof/N_dof
    MathOperatorLHS<dim,CELL>::LHS.Resize( BT_dof, N_dof );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    for ( uint32_t i{0}; i < e.IntegrationPoints(); ++i )
    {
        // --- B matrix and Jacobian determinant ---
        double detJ = e.dN_AtIntegrationPoint( B_, i, dim );

        if ( detJ <= 0. ) {
            std::cerr << "\n\tElement " << e.Idx()
                      << ": detJ at Gauss point " << i << ": " << detJ << "\n";
            throw csmp::Exception( FATAL_ERROR,
                                   "NumJacobianIntegral_BT_m_N_dV::ComputeContribution",
                                   "Jacobian transformation failed." );
        }

        // --- Pressure shape functions N_p at this integration point ---
        e.N_AtIntegrationPoint( i, Np_ );

        // Scale by alpha * weight * detJ at this integration point
        const double alpha_i = ( alpha_.size() == 1 ) ? alpha_[0] : alpha_[i];
        const double scale   = alpha_i * e.WeightAtIntegrationPoint(i) * detJ;

        // --- Direct Evaluation of B^T * m * Np * scale ---
        // B^T * m produces a column vector of size (BT_dof x 1).
        // Multiplying that by the row vector Np produces the (BT_dof x N_dof) block.
        for ( uint32_t r = 0; r < BT_dof; ++r )
        {
            // 1. Calculate the r-th component of (B^T * m)
            // (σ′ =σ + α p m (compression positive) sign convention
            double BTm_r = 0.0;
            for ( uint32_t s = 0; s < n_strains_; ++s ) {
                // Note that B^T(r, s) is identical to B_(s, r)
                BTm_r += B_(s, r) * m_[s]; 
            }
            
            // 2. Fold in the integration scale for this row
            BTm_r *= scale;

            // 3. Compute outer product against the Np row vector and accumulate
            for ( uint32_t c = 0; c < N_dof; ++c ) {
                MathOperatorLHS<dim,CELL>::LHS(r, c) += BTm_r * Np_[c];
            }
        }
        
    } // for gp

} // end ComputeContribution


// NOTES
// Q has the dimensions 6 x 3 for the linear triangle and occupies the off-diagonal block:
//
//   global row offset = 0          (displacement block starts at 0)
//   global col offset = n_disp     (pressure block starts after displacement)
//
// For each element e:
//   u_dofs[i] = interleaved displacement DOFs  (size 2*nodes)
//   p_dofs[j] = pressure DOFs = node indices   (size nodes)
//
//   global_matrix( u_dofs[i], n_disp + p_dofs[j] ) += Q^e(i,j)
//   global_matrix( n_disp + p_dofs[j], u_dofs[i] ) += Q^e(i,j)  // Q^T block

//template class NumJacobianIntegral_BT_m_N_dV<1U,Element>;
template class NumJacobianIntegral_BT_m_N_dV<2U,Element>;
template class NumJacobianIntegral_BT_m_N_dV<3U,Element>;

//template class NumJacobianIntegral_BT_m_N_dV<1U,Face>;
template class NumJacobianIntegral_BT_m_N_dV<2U,Face>;
template class NumJacobianIntegral_BT_m_N_dV<3U,Face>;

} // end csmp
