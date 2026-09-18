// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_NT_rhsop_N_dV.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_rhsop_N_dV<dim,CELL>::NumIntegral_NT_rhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                           const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>( pref, oper, test ),
    op_( 1U, 1.0 )
 {
    MathOperatorRHS<dim,CELL>::Name( "NumIntegral_NT_rhsop_N_dV", oper, test );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_rhsop_N_dV<dim>::(constructor)",
                             oper, "Material operand must be a scalar." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType()      != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_rhsop_N_dV<dim>::(constructor)",
                             test, "Test operand must be a scalar property placed on the nodes." );

 } // end constructor


/**
Reads the scalar material operand into the local vector op_.
For element/face/region placement a single value is stored.
For node or integration point placement one value per integration
point is stored, interpolated to the integration points via
PropertyAtIntegrationPoint.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_rhsop_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    switch ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() )
      {
        case ELEMENT:
        case FACE:
        case REGION:
          op_.resize( 1U );
          op_[0U] = e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey() );
          break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
          {
            const auto n_ip{ e.IntegrationPoints() };
            op_.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              op_[i] = e.Read( i, MathOperatorRHS<dim,CELL>::MaterialOperandKey() );
          }
          break;

        case NODE:
          {
            const auto n_ip{ e.IntegrationPoints() };
            op_.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              op_[i] = e.PropertyValueAtIntegrationPoint( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), i );
          }
          break;

        default:
          throw csmp::Exception( FATAL_ERROR, "NumIntegral_NT_rhsop_N_dV::GetOperands",
                                 "Unsupported material operand placement." );
      }

 } // end GetOperands


/**
Computes the volume integral:

  {RHS}[j] = integral( N_j * op * N_k ) dV

Row-sum lumping is applied in the lumped formulation.
*/
/**
Computes the volume integral for the RHS vector:
  Consistent: {RHS}[j] = integral( N_j * op ) dV
  Lumped:     {RHS}[j] = op * Volume / n_nodes
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_rhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    assert( e.UsesLocalCoordinates() == true );

    const uint32_t n_nodes{    e.Nodes()             };
    const uint32_t n_ipoints{  e.IntegrationPoints() };

    // Ensure RHS vector is clean before accumulation
    MathOperatorRHS<dim,CELL>::RHS.assign( n_nodes, 0.0 );

    const bool piecewise_constant( op_.size() == 1U );
    const bool node_placed( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE );

    // ------------------------------------------------------------------
    // 1. LUMPED FORMULATION (Equal distribution)
    // ------------------------------------------------------------------
    if ( MathOperatorRHS<dim,CELL>::LumpedFormulation() )
    {
        const double share( e.Volume() / static_cast<double>(n_nodes) );

        if ( piecewise_constant )
        {
            const double val( op_[0U] * share );
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                MathOperatorRHS<dim,CELL>::RHS[j] = val;
        }
        else if ( node_placed )
        {
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                MathOperatorRHS<dim,CELL>::RHS[j] = op_[j] * share;
        }
        else
        {
            double op_avg{0.0};
            for ( uint32_t i{0U}; i < n_ipoints; ++i )
                op_avg += op_[i];
            op_avg /= static_cast<double>(n_ipoints);
            
            const double val( op_avg * share );
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                MathOperatorRHS<dim,CELL>::RHS[j] = val;
        }
        return; 
    }

    // ------------------------------------------------------------------
    // 2. CONSISTENT FORMULATION (Direct Vector Evaluation)
    // ------------------------------------------------------------------
    // We evaluate: RHS[j] = integral( N_j * op ) dV directly.
    for ( uint32_t i{0U}; i < n_ipoints; ++i )
    {
        e.N_AtIntegrationPoint( i, e.FE()->NRST );
        const double detJ( e.det_J_AtIntegrationPoint(i) );
        const double w(    e.WeightAtIntegrationPoint(i)  );

        double op_val{0.0}; 
        
        if ( piecewise_constant )
        {
            op_val = op_[0U];
        }
        else if ( node_placed )
        {
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                op_val += e.FE()->NRST[j] * op_[j];
        }
        else
        {
            op_val = op_[i]; 
        }

        // The combined differential volume scalar for this Gauss point
        const double dV = detJ * w * op_val;

        // Directly accumulate into the RHS vector (O(n) complexity)
        for ( uint32_t j{0U}; j < n_nodes; ++j )
        {
            MathOperatorRHS<dim,CELL>::RHS[j] += e.FE()->NRST[j] * dV;
        }
    }

} // end ComputeContribution

template class NumIntegral_NT_rhsop_N_dV<1U,Element>;
template class NumIntegral_NT_rhsop_N_dV<2U,Element>;
template class NumIntegral_NT_rhsop_N_dV<3U,Element>;

template class NumIntegral_NT_rhsop_N_dV<1U,Face>;
template class NumIntegral_NT_rhsop_N_dV<2U,Face>;
template class NumIntegral_NT_rhsop_N_dV<3U,Face>;

} // namespace csmp

