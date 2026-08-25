#include "NumIntegral_PT_op_dS.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {
  

template<uint32_t dim>
NumIntegral_PT_op_dS<dim>::NumIntegral_PT_op_dS( const PropertyDatabase<dim>& pref,
                                                 const char* oper, const char* test )
  : MathOperatorRHS<dim,Face>(pref,oper,test)
 {
    MathOperatorRHS<dim,Face>::Name("NumIntegral_PT_op_dS", oper, test );

    if ( MathOperatorRHS<dim,Face>::MaterialOperandPlacement() != FACE ||
        ( MathOperatorRHS<dim,Face>::MaterialOperandType() != VECTOR && MathOperatorRHS<dim,Face>::MaterialOperandType() != SCALAR ) )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_dS<dim>::(constructor):",
                             oper, "Operand must be a vector or scalar variable on the face.");


    if ( MathOperatorRHS<dim,Face>::TestOperandType() != VECTOR ||
         MathOperatorRHS<dim,Face>::TestOperandPlacement() != NODE )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_dS<dim>::(constructor):",
                             test,  "Dependent-variable Operand must be a vector variable placed on the nodes.");   
 }



/**
Reads the VectorVariable<dim> boundary stress (op) which is specified as face 
variable and flagged as Neumann condition for further processing
by ComputeContribution() method. If (op) is a ScalarVariable, it is multiplied
by the face unit normal and assumed to act inward if positive.

@param f The element from which accumulation into the righthand vector takes place.

*/
template<uint32_t dim>
void NumIntegral_PT_op_dS<dim>::GetOperands( const Face<dim>& f )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( f.FE()->Isoparametric() == true );

    if( MathOperatorRHS<dim,Face>::MaterialOperandType() == VECTOR )
      {
        f.Read( this->MaterialOperandKey(), oper_ );
        for ( uint32_t i{0U}; i < dim; ++i )
          if ( oper_.Flag(i) != NEUMANN )
            oper_(i) = 0.;
        return;
      }

    // scalar, multiply negative since positive assumed compressive and unit normal of faces pointing out
    f.Read( this->MaterialOperandKey(), scalar_ );
    if ( scalar_.Flag() != NEUMANN )
      scalar_() = 0.;
    f.UnitNormal(oper_);
    oper_ *= -scalar_();
}




/**
    @brief Computes the Neumann traction surface integral over a boundary face
    and accumulates the result into the element right-hand vector.

    Evaluates the weak-form surface integral:

        f_j^{df} = ∫_Γ N_j · t_{df} dS
                 ≈ Σ_i w_i |J_i| N_j(ξ_i) · t_{df}

    where:
    - N_j(ξ_i)  is the shape function of node j evaluated at integration point i
    - w_i       is the quadrature weight at integration point i
    - |J_i|     is the surface Jacobian determinant at integration point i
    - t_{df}    is the df-th component of the traction vector (material operand)

    The result is stored in RHS with node-major DOF ordering:
        RHS[ j * dim + df ]
    where j is the local node index and df is the spatial degree of freedom.

    @note The NEUMANN status check is performed by the PDE_Integrator framework
    in AccumulateBoundaryIntegrals before this method is called. By the time
    ComputeContribution is reached, the traction has already been confirmed as
    a Neumann condition and read into the material operand by GetOperands.

    @note RHS is zero-initialised at the start of each call. Components for
    which the traction is zero contribute nothing to the integral.

    @param f  The boundary face over which the surface integral is evaluated.

    @par Application
    Used within the PDE_Integrator framework to apply Neumann traction boundary
    conditions to the right-hand side of the linear system. Typical uses include:
    - Far-field in-situ stress conditions on outer model boundaries
    - Mud pressure on the borehole wall (BOUNDARY_WELL)
    - Any distributed surface load expressed as a traction vector

    @par References
    Zienkiewicz & Taylor, The Finite Element Method, Vol. 1, Ch. 3
*/
template<uint32_t dim>
void NumIntegral_PT_op_dS<dim>::ComputeContribution( const Face<dim>& f )
{
    // initialise RHS to zero — one entry per node per spatial DOF
    MathOperatorRHS<dim,Face>::RHS.resize( f.Nodes() * dim );
    fill( MathOperatorRHS<dim,Face>::RHS.begin(),
          MathOperatorRHS<dim,Face>::RHS.end(), 0. );

    vector<double> N;
    for ( uint32_t i{ 0U }; i < f.IntegrationPoints(); ++i ) {
        N.clear();
        f.N_AtIntegrationPoint( i, N );
        const double detJ   = f.det_J_AtIntegrationPoint( i );
        const double weight = f.WeightAtIntegrationPoint( i );

        // accumulate weighted traction contribution for each node and DOF:
        //   RHS[j*dim+df] += N_j(ξ_i) · w_i · |J_i| · t_{df}
        for ( uint32_t j{ 0U }; j < f.Nodes(); ++j )
            for ( uint32_t df{ 0U }; df < dim; ++df )
                MathOperatorRHS<dim,Face>::RHS[j*dim+df] +=
                    N[j] * weight * detJ * oper_(df);
    }
}



template class NumIntegral_PT_op_dS<1U>;
template class NumIntegral_PT_op_dS<2U>;
template class NumIntegral_PT_op_dS<3U>;

} //csmp
