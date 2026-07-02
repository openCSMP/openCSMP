//
//  NumJacobianIntegral_N_mT_B_dV.cpp - transposed version of NumJacobianIntegral_BT_m_N_dV
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 22/6/2026.
//  Copyright © 2026 Stephan Matthai. All rights reserved.
//

#include "NumJacobianIntegral_N_mT_B_dV.h"
#include "Element.h"
#include "Face.h"

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumJacobianIntegral_N_mT_B_dV<dim,CELL>::NumJacobianIntegral_N_mT_B_dV( const PropertyDatabase<dim>& pref,
                                                                        const char* oper,
                                                                        const char* basic,
                                                                        const char* test )
  : MathOperatorLHS<dim,CELL>(pref, oper, basic, test),
    alpha_(1U)
{
    MathOperatorLHS<dim,CELL>::Name("NumJacobianIntegral_N_mT_B_dV", oper, basic, test);

    // --- Validate Biot coefficient ---
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumJacobianIntegral_N_mT_B_dV::(constructor)",
                        oper, "Biot coefficient must be a scalar property." );

    // --- Validate pressure (test) ---
    // NOTE: test is scalar (pressure), unlike the stiffness matrix
    // where both basic and test are vector (displacement).
    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType()      != SCALAR )
        throw csmp::Exception( ERROR, "NumJacobianIntegral_N_mT_B_dV::(constructor)",
                        test, "operand must be a scalar property discretised on the nodes." );

    // --- Validate displacement (basic) ---
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType()      != VECTOR )
        throw csmp::Exception( ERROR, "NumJacobianIntegral_N_mT_B_dV::(constructor)",
                        basic, "operand must be a vector property discretised on the nodes." );
}


// -----------------------------------------------------------------------
// GetOperands — mirrors your existing pattern exactly
// -----------------------------------------------------------------------
template<uint32_t dim, template<uint32_t> class CELL>
void NumJacobianIntegral_N_mT_B_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
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
// Output matrix is (2*nodes x nodes):
//   rows -> displacement DOFs (interleaved)
//   cols -> pressure DOFs
// -----------------------------------------------------------------------
template<uint32_t dim, template<uint32_t> class CELL>
void NumJacobianIntegral_N_mT_B_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    const uint32_t n_nodes   = e.Nodes();
    const uint32_t B_dof     = dim * n_nodes;   // rows: 2*nodes
    const uint32_t N_dof     = n_nodes;         // cols: nodes

    // Output: (2*nodes x nodes)
    MathOperatorLHS<dim,CELL>::LHS.Resize( N_dof, B_dof );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    // Resize working matrices to actual element size
    LHSi_.Resize(n_nodes,n_strains_);
      
    for ( uint32_t i{0}; i < e.IntegrationPoints(); ++i )
    {
        // --- B matrix and Jacobian determinant ---
        // dN_AtIntegrationPoint fills B in global coordinates,
        // exactly as in your stiffness matrix code.
        double detJ = e.dN_AtIntegrationPoint( B_, i, dim );

        if ( detJ <= 0. ) {
            std::cerr << "\n\tElement " << e.Idx()
                      << ": detJ at Gauss point " << i << ": " << detJ << "\n";
            throw csmp::Exception( FATAL_ERROR,
                           "NumJacobianIntegral_BT_m_Np_dV::ComputeContribution",
                           "Jacobian transformation failed." );
        }

        // --- Pressure shape functions N_p at this integration point ---
        // N_p is a row vector (1 x n_nodes)
        e.N_AtIntegrationPoint( i, Np_ );

        // Scale by alpha * weight * detJ at this integration point
        const double alpha_i = ( alpha_.size() == 1 ) ? alpha_[0] : alpha_[i];
        const double scale   = alpha_i * e.WeightAtIntegrationPoint(i) * detJ;

        // --- Outer product: Np * mT

        // Initialize an (N_n x 6) matrix with zeros
        LHSi_.Zero();

        for ( uint32_t m = 0; m < n_nodes; ++m )
          for ( uint32_t n = 0; n < n_strains_; ++n )
            LHSi_(m,n) = scale *  Np_[m] * m_[n];

        // Accumulate
        MathOperatorLHS<dim,CELL>::LHS += LHSi_ * B_;
    }
}

//template class NumJacobianIntegral_N_mT_B_dV<1U,Element>;
template class NumJacobianIntegral_N_mT_B_dV<2U,Element>;
template class NumJacobianIntegral_N_mT_B_dV<3U,Element>;

//template class NumJacobianIntegral_N_mT_B_dV<1U,Face>;
template class NumJacobianIntegral_N_mT_B_dV<2U,Face>;
template class NumJacobianIntegral_N_mT_B_dV<3U,Face>;


// TESTING (correct)
// In NumJacobianIntegral_N_mT_B_dV::ComputeContribution, after integration loop:
/*
if ( e.Nodes() == 4 ) {
    std::cout << "\nlhs_QT quad LHS (" 
              << MathOperatorLHS<dim,CELL>::LHS.Rows() << "x"
              << MathOperatorLHS<dim,CELL>::LHS.Cols() << "):\n";
    for ( uint32_t r=0; r<MathOperatorLHS<dim,CELL>::LHS.Rows(); ++r ) {
        std::cout << "  row " << r << ": ";
        for ( uint32_t c=0; c<MathOperatorLHS<dim,CELL>::LHS.Cols(); ++c )
            std::cout << std::setw(10) << std::fixed << std::setprecision(6)
                      << MathOperatorLHS<dim,CELL>::LHS(r,c);
        std::cout << "\n";
    }
}
*/

} // end csmp
