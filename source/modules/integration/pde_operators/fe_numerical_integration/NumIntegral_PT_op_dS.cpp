#include "NumIntegral_PT_op_dS.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {
  

template<uint32_t dim>
NumIntegral_PT_op_dS<dim>::NumIntegral_PT_op_dS( const PropertyDatabase<dim>& pref,
                                                 const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_PT_op_dS", oper, test );

    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != FACE || ( MathOperatorRHS<dim>::MaterialOperandType() != VECTOR && MathOperatorRHS<dim>::MaterialOperandType() != SCALAR ) ) 
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_dS<dim>::(constructor):",
                             oper, "Operand must be a vector or scalar variable on the face.");


    if ( MathOperatorRHS<dim>::TestOperandType() != VECTOR || 
         MathOperatorRHS<dim>::TestOperandPlacement() != NODE ) 
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

    if( MathOperatorRHS<dim>::MaterialOperandType() == VECTOR )
      {
        f.Read( this->MaterialOperandKey(), oper_ );
        for( auto i{0U}; i < dim; ++i )
          if( oper_.Flag(i) != NEUMANN )
            oper_(i) = 0.;
        return;
      }

    // scalar, multiply negative since positive assumed compressive and unit normal of faces pointing out
    f.Read( this->MaterialOperandKey(), scalar_ );
    if( scalar_.Flag() != NEUMANN )
      scalar_() = 0.;
    f.UnitNormal(oper_);
    oper_ *= -scalar_();
}

/**
Calculates the surface integral over the boundary face of the element
and assigns this contribution in a weighted fashion to the boundary
nodes in the righthand vector. Since we checked for NEUMANN in GetOperand()
we don't need to do it here. Vector components are set to zero otherwise.

@param f The element from which accumulation into the global matrix takes placed.  

The result is returned into the MathOperatorRHS vector.

@section application Application
Within the PDE_Integrator framework to assign stress boundary conditions to a model.  
*/
template<uint32_t dim>
void NumIntegral_PT_op_dS<dim>::ComputeContribution( const Face<dim>& f )
{
   MathOperatorRHS<dim>::RHS.resize( f.Nodes() * dim );
   fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

   std::vector<double> N(0);
   for ( auto i{0U}; i < f.IntegrationPoints(); ++i )
     {
       N.clear();
       f.N_AtIntegrationPoint( i, N );
       double const detJ = f.det_JINV_AtIntegrationPoint(i);
       double const weight = f.WeightAtIntegrationPoint(i);
       for( size_t j(0); j < f.Nodes(); ++j )
         for( int df(0); df < dim; ++df )
           MathOperatorRHS<dim>::RHS[j*dim+df] += N[j] * weight * detJ * oper_(df) ;
     }
} 


template class NumIntegral_PT_op_dS<1U>;
template class NumIntegral_PT_op_dS<2U>;
template class NumIntegral_PT_op_dS<3U>;

} //csmp
