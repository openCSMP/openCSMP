#include "NumIntegral_dNT_lhsop_dN_dV.h"

#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"
#include "denseMatrixMethods.h"
#include <Eigen/Dense>

using namespace std;

namespace csmp {

/**
 
The Operand which is used here can be both, an element or a nodal variable
which is then interpolated to the integration points to obtain the 
integral properties.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_lhsop_dN_dV<dim,CELL>::NumIntegral_dNT_lhsop_dN_dV( const PropertyDatabase<dim>& pref,
                                                                    const char*           oper,
                                                                    const char*           basic,
                                                                    const char*           test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_lhsop_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_lhsop_dN_dV<dim>::(constructor)", 
                             basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_lhsop_dN_dV<dim>::(constructor)", 
                             test, "Operand (test) must be a scalar property placed on the nodes." );
}








/** Laplacian operator of shape function derivatives squared.

    @attention since the interpolation functions derivatives are constant across simplices (linear line, triangle and tetrahedral elements),
    the computation of DN for these elements simplifies greatly.
    
    @test works for poly-element meshes
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_lhsop_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    assert( e.UsesLocalCoordinates() );

    // Output matrix: nodes x nodes, zeroed
    auto& LHS  = MathOperatorLHS<dim,CELL>::LHS;
    auto& MTRL = MathOperatorLHS<dim,CELL>::MTRL;

    LHS.Resize( e.Nodes(), e.Nodes() );
    LHS.Zero();

    // Local scratch matrices — no shared mutable state, safe for parallel assembly
    DenseMatrix<DM_MIN> B, BT;

    const bool is_simplex( e.FE()->IsSimplex() && e.Interpolation() == 1 );

    const bool piecewise_constant_material( this->MaterialOperandPlacement() == ELEMENT ||
                                            this->MaterialOperandPlacement() == REGION  ||
                                            this->MaterialOperandPlacement() == FACE );

    // ==================================================================
    // FAST PATH: linear simplex
    //
    // For a linear simplex, grad(N) is constant throughout the element.
    // The integral therefore reduces to:
    //
    //   K_e = B^T * K * B * V_e
    //
    // where V_e = |J| * V_ref is the physical element volume, and
    // V_ref is the reference simplex volume in parametric space:
    //   triangle:     V_ref = 1/2
    //   tetrahedron:  V_ref = 1/6
    //
    // |J| is returned by dN_AtBaryCenter. The product |J| * V_ref is
    // recovered as detJ * sum(w_i), which equals detJ * V_ref exactly
    // for any quadrature rule that integrates constants exactly (all do).
    // ==================================================================
    if ( is_simplex )
    {
        // B and detJ are constant for a linear simplex — compute once
        const double detJ = e.dN_AtBaryCenter( B );
        B.Transposed( BT );

        // Physical element volume = |J| * V_ref
        // V_ref = sum of all integration weights for this element
        double V_ref = 0.;
        for ( uint32_t i{0U}; i < e.IntegrationPoints(); ++i )
            V_ref += e.WeightAtIntegrationPoint( i );
        const double V_e = detJ * V_ref;

        if ( piecewise_constant_material )
        {
            // Single material matrix — one triple product scaled by V_e
            BT *= MTRL[0];
            BT *= B;
            BT *= V_e;
            LHS += BT;
        }
        else
        {
            // NODE or ELEMENT_INTEGRATION_POINT material:
            // material varies per integration point but B does not.
            // Accumulate contribution of each integration point separately,
            // each weighted by its own w_i * |J|.
            for ( uint32_t i{0U}; i < e.IntegrationPoints(); ++i )
            {
                B.Transposed( BT );          // reset BT from the constant B
                BT *= MTRL[i];                                     // apply material at ip i
                BT *= B;                                           // apply constant B
                BT *= e.WeightAtIntegrationPoint(i) * detJ;       // scale by w_i * |J|
                LHS += BT;
            }
        }

        applyLumping( LHS,
                      MathOperatorLHS<dim,CELL>::LumpedFormulation(),
                      true, /* lump_by_volume, else you get zeroes in diagonal */
                      V_e,
                      e.Nodes() );
        return;
    }

    // ==================================================================
    // GENERAL PATH: non-simplex or higher-order elements
    //
    // B and detJ vary per integration point — full Gauss loop required.
    // ==================================================================
    double V_ref = 0.;
    for ( uint32_t i{0U}; i < e.IntegrationPoints(); ++i )
    {
        const double detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );
        B.Transposed( BT );

        // Select material matrix: constant for piecewise-constant placement,
        // per-integration-point for NODE or ELEMENT_INTEGRATION_POINT placement
        const uint32_t mtrl_idx = piecewise_constant_material ? 0U : i;

        BT *= MTRL[mtrl_idx];
        BT *= B;
        BT *= e.WeightAtIntegrationPoint(i) * detJ;
        LHS += BT;
        
        V_ref += e.WeightAtIntegrationPoint(i) * detJ;
    }

    applyLumping( LHS,
                  MathOperatorLHS<dim,CELL>::LumpedFormulation(),
                  true, /* lump_by_volume, else you get zeroes in diagonal */ 
                  V_ref,
                  e.Nodes() );

} // end ComputeContribution


template class NumIntegral_dNT_lhsop_dN_dV<1U>;
template class NumIntegral_dNT_lhsop_dN_dV<2U>;
template class NumIntegral_dNT_lhsop_dN_dV<3U>;

template class NumIntegral_dNT_lhsop_dN_dV<1U,Face>;
template class NumIntegral_dNT_lhsop_dN_dV<2U,Face>;
template class NumIntegral_dNT_lhsop_dN_dV<3U,Face>;

template class NumIntegral_dNT_lhsop_dN_dV<1U,InterFace>;
template class NumIntegral_dNT_lhsop_dN_dV<2U,InterFace>;
template class NumIntegral_dNT_lhsop_dN_dV<3U,InterFace>;

} // csmp











