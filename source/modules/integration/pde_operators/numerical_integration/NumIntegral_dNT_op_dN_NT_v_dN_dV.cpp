// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_dNT_op_dN_NT_v_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>& pref,
                                                                         const char* diffusion_oper,   // element prop, for instance thermal conductivity
                                                                         const char* advection_oper,   // element prop, for instance heat transport velocity
                                                                         const char* basic,            // e.g., fluid pressure
                                                                         const char* test ) 
  : MathOperatorLHS<dim,CELL>(pref,diffusion_oper,basic,test),
    adv_key(pref.StorageKey(advection_oper)),
    DN(dim,3), 
    DNT(3,dim),
    IPOL(3),
    VIP(dim,dim),
    NT3(3,dim)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_op_dN_NT_v_dN_dV", diffusion_oper, basic, test );
    
    for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
          it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ )
      (*it).Resize(dim,dim);
      
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      diffusion_oper, "must be an element-based variable." );

    if ( (adv_key.place != ELEMENT or adv_key.place != ELEMENT_INTEGRATION_POINT) and adv_key.type != VECTOR )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      advection_oper, "must be an element-based vector variable." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}









/**
 
The diffusion (op) and advection (adv) coefficients are read from the storage in the model. 

*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // 1. read diffusion Operand
    // -------------------------
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      {
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
              double sc = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = sc;
           }
         else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = vc[i];
           }
         else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                for ( uint32_t j{0U}; j<dim; j++ ) 
                  MathOperatorLHS<dim,CELL>::MTRL[0](i,j) = ts(i,j);
           }
      }
    else { // or an integration point variable
         for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ )
           MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim,CELL>::MaterialOperandKey(),
                                                                  i, MathOperatorLHS<dim,CELL>::MTRL[i] );
      }
   
    // 2. read element advection variable
    // -----------------------------------
    if ( adv_key.place == ELEMENT )
      e.Read( adv_key, velo_ );
        
} // end GetOperands






/**
    @brief Computes the element contribution to the advection-diffusion stiffness matrix.

    Evaluates the combined integral:

        K_jk = ∫_Ω (∇N_j)ᵀ [σ] ∇N_k dV  +  ∫_Ω N_j (v · ∇N_k) dV

    The first term is the diffusive contribution — the standard conductance
    matrix weighted by the material operand [σ] (scalar, vector or tensor
    diffusivity / conductivity / permeability).

    The second term is the advective contribution — the Petrov-Galerkin
    upwind term weighted by the advection velocity v. It is assembled as
    the outer product:

        A_jk = N_j · (v · ∇N_k) = N_j · Σ_d v_d · ∂N_k/∂x_d

    Both terms are integrated using the same full Gauss quadrature loop
    over all element integration points. This is essential for bilinear
    quadrilateral elements where the product N_j · (v · ∇N_k) is quadratic
    and requires at least a 2×2 rule for exact integration.

    @note NT3 is zeroed at the start of each integration point iteration.
    Without this reset, stale values from the previous integration point
    accumulate into the current one, producing a systematic error of
    exactly 1/8 per entry for Q4 elements with a 2×2 Gauss rule.

    @note The material operand [σ] may be element-placed (constant over
    the element, using MTRL[0]) or integration-point-placed (varying,
    using MTRL[i]). The advection velocity v is always element-placed.

    @param e  The cell over which the contribution is computed.
              Must be an isoparametric element (UsesLocalCoordinates() == true).
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // initialise output matrix to zero
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    for ( uint32_t i{0U}; i < e.IntegrationPoints(); ++i )
      {
        // ----------------------------------------------------------------
        // 1. Diffusive contribution: ∫ (∇N)ᵀ [σ] ∇N dV
        //
        // Compute the global shape function derivative matrix DN and the
        // Jacobian determinant detJ at integration point i.
        // DN is already expressed in global (physical) coordinates.
        // ----------------------------------------------------------------
        double detJ = e.dN_AtIntegrationPoint( DN, i, SCALAR );

        // transpose DN → DNT
        DN.Transposed( DNT );

        // multiply DNT · [σ]
        // use MTRL[0] for element-placed operand, MTRL[i] for ip-placed
        if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
            DNT *= MathOperatorLHS<dim,CELL>::MTRL[0];
        else
            DNT *= MathOperatorLHS<dim,CELL>::MTRL[i];

        // multiply DNT · DN  →  (∇N)ᵀ [σ] ∇N
        DNT *= DN;

        // scale by quadrature weight and Jacobian determinant
        DNT *= e.WeightAtIntegrationPoint(i) * detJ;

        // accumulate diffusive contribution into element stiffness
        MathOperatorLHS<dim,CELL>::LHS += DNT;

        // ----------------------------------------------------------------
        // 2. Advective contribution: ∫ N_j (v · ∇N_k) dV
        //
        // Assemble as the outer product:
        //   A_jk = N_j · (v · ∇N_k)
        //        = N_j · Σ_d v_d · DN(d,k)
        //
        // ----------------------------------------------------------------
        // build diagonal velocity matrix VIP = diag(v_x, v_y [, v_z])
        VIP.Zero();
        if  ( adv_key.place == ELEMENT_INTEGRATION_POINT ||
              adv_key.place == FACE_INTEGRATION_POINT ) {
             e.Read( adv_key, velo_ );
          }
        for ( uint32_t j{0U}; j < dim; ++j ) VIP(j,j) = velo_[j];

        // evaluate shape functions at integration point i
        e.N_AtIntegrationPoint( i, IPOL );

        // NT3 (n_nodes × dim) is zeroed at the start of each integration
        // point to prevent stale values from the previous point accumulating
        // into the current one — this was the source of the ±1/8 error.
        NT3.Resize( e.Nodes(), dim );
        NT3.Zero();
        // fill NT3: column d of NT3 = N_j for all nodes j
        // after multiplication by VIP this gives NT3(j,d) = N_j · v_d
        for ( uint32_t j{0U}; j < e.Nodes(); ++j ) {
            NT3(j,0) = IPOL[j];
            if constexpr ( dim >= 2U ) NT3(j,1) = IPOL[j];
            if constexpr ( dim == 3U ) NT3(j,2) = IPOL[j];
        }

        // NT3 (n_nodes × dim) · VIP (dim × dim) → NT3(j,d) = N_j · v_d
        NT3 *= VIP;

        // NT3 (n_nodes × dim) · DN (dim × n_nodes) → NT3(j,k) = N_j · (v · ∇N_k)
        // DN was already computed at this integration point in step 1
        NT3 *= DN;

        // scale by quadrature weight and Jacobian determinant
        NT3 *= e.WeightAtIntegrationPoint(i) * detJ;

        // accumulate advective contribution into element stiffness
        MathOperatorLHS<dim,CELL>::LHS += NT3;

      } // end integration point loop

} // end ComputeContribution




template class NumIntegral_dNT_op_dN_NT_v_dN_dV<1U,Element>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<2U,Element>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<3U,Element>;


template class NumIntegral_dNT_op_dN_NT_v_dN_dV<1U,Face>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<2U,Face>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<3U,Face>;

} // csmp





