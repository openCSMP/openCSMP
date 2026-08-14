#include "NumIntegral_NT_lhsop_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_lhsop_N_dV<dim,CELL>::NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                                const char* oper,
                                                                const char* basic,
                                                                const char* test )
  : MathOperatorLHS<dim,CELL>( pref, oper, basic, test ),
    op_( 1U, 1.0 )
 {
    MathOperatorLHS<dim,CELL>::Name( "NumIntegral_NT_lhsop_N_dV", oper, basic, test );

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)",
                             oper, "Material operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType()      != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)",
                             test, "Test operand must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType()      != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)",
                             basic, "Basic operand must be a scalar property placed on the nodes." );

 } // end constructor


/**
Reads the scalar material operand into the local vector op_.
For element/face/region placement a single value is stored.
For node or integration point placement one value per integration
point is stored, interpolated to the integration points via
PropertyAtIntegrationPoint.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_lhsop_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    switch ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() )
      {
        case ELEMENT:
        case FACE:
        case REGION:
          op_.resize( 1U );
          op_[0U] = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
          break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
          {
            const auto n_ip{ e.IntegrationPoints() };
            op_.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              op_[i] = e.Read( i, MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
          }
          break;

        case NODE:
          {
            const auto n_ip{ e.IntegrationPoints() };
            op_.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              {
                op_[i] = e.PropertyValueAtIntegrationPoint( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), i );
              }
          }
          break;

        default:
          throw csmp::Exception( FATAL_ERROR, "NumIntegral_NT_lhsop_N_dV::GetOperands",
                                 "Unsupported material operand placement." );
      }

 } // end GetOperands


/**
Computes the volume integral of the capacitance matrix:

  C[n x n] = integral( N^T * op * N ) dV

In the lumped formulation a diagonal matrix is produced.
In the consistent formulation the full mass matrix is assembled.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_lhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    assert( e.FE()->Isoparametric() == true );

    const uint32_t n_nodes{    e.Nodes()             };
    const uint32_t n_ipoints{  e.IntegrationPoints() };

    // Initialize and clear the LHS matrix container
    MathOperatorLHS<dim,CELL>::LHS.Resize( n_nodes, n_nodes );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    const bool piecewise_constant( op_.size() == 1U );
    const bool node_placed( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE );

    // ------------------------------------------------------------------
    // 1. LUMPED FORMULATION — Diagonal Matrix
    // ------------------------------------------------------------------
    if ( MathOperatorLHS<dim,CELL>::LumpedFormulation() )
    {
        const double share( e.Volume() / static_cast<double>(n_nodes) );

        if ( piecewise_constant )
        {
            const double val( op_[0U] * share );
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                MathOperatorLHS<dim,CELL>::LHS(j, j) = val;
        }
        else if ( node_placed )
        {
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                MathOperatorLHS<dim,CELL>::LHS(j, j) = op_[j] * share;
        }
        else 
        {
            // Integration point placement — average over integration points
            double op_avg{0.0};
            for ( uint32_t i{0U}; i < n_ipoints; ++i )
                op_avg += op_[i];
            op_avg /= static_cast<double>(n_ipoints);
            
            const double val( op_avg * share );
            for ( uint32_t j{0U}; j < n_nodes; ++j )
                MathOperatorLHS<dim,CELL>::LHS(j, j) = val;
        }
        return;
    }

    // ------------------------------------------------------------------
    // 2. CONSISTENT FORMULATION — Full Matrix (Direct O(n²) Accumulation)
    // ------------------------------------------------------------------
    // Note: `is_simplex` shortcut is removed. Triangles and Tetrahedra
    // loop through their full quadrature rules to accurately integrate N_j * N_k.
    
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

        // Combine the integration scaling factors outside the inner loops
        const double dV_scaled = detJ * w * op_val;

        // Compute the outer product directly into the LHS matrix.
        // This avoids creating any temporary matrix objects or heap allocations.
        for ( uint32_t j{0U}; j < n_nodes; ++j )
        {
            const double NJ_dV = e.FE()->NRST[j] * dV_scaled;
            for ( uint32_t k{0U}; k < n_nodes; ++k )
            {
                MathOperatorLHS<dim,CELL>::LHS(j, k) += NJ_dV * e.FE()->NRST[k];
            }
        }
    }

} // end ComputeContribution

template class NumIntegral_NT_lhsop_N_dV<1U,Element>;
template class NumIntegral_NT_lhsop_N_dV<2U,Element>;
template class NumIntegral_NT_lhsop_N_dV<3U,Element>;

template class NumIntegral_NT_lhsop_N_dV<1U,Face>;
template class NumIntegral_NT_lhsop_N_dV<2U,Face>;
template class NumIntegral_NT_lhsop_N_dV<3U,Face>;

template class NumIntegral_NT_lhsop_N_dV<1U,InterFace>;
template class NumIntegral_NT_lhsop_N_dV<2U,InterFace>;
template class NumIntegral_NT_lhsop_N_dV<3U,InterFace>;

} // namespace csmp

