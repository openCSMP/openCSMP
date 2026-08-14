#include "NumIntegral_DNi_rhsop_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_DNi_rhsop_dV<dim,CELL>::~NumIntegral_DNi_rhsop_dV() {}

/// custom constructor that should be used
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_DNi_rhsop_dV<dim,CELL>::NumIntegral_DNi_rhsop_dV( const PropertyDatabase<dim>& pref,
                                                                 const char*                  oper,
                                                                 const char*                  test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    xyz_(Y_DIRECTION), op_vec_(dim)
{
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_DNi_rhsop_dV",oper, test );
    

    if ( MathOperatorRHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_DNi_rhsop_dV<dim>::(constructor)",
                             oper, "Operand must be a scalar property placed on the nodes." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR)
      throw csmp::Exception( ERROR, "NumIntegral_DNi_rhsop_dV<dim>::(constructor)",
                             test, "Operand (test) must be a scalar property placed on the nodes." );
}




string parse( SPATIAL_DERIVATIVE deriv ) {
    if ( deriv == X_DIRECTION ) return "X-direction";
    if ( deriv == Y_DIRECTION ) return "Y-direction";
    if ( deriv == Z_DIRECTION ) return "Z-direction";
    return string("undefined");
 }





template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNi_rhsop_dV<dim,CELL>::SpatialDerivative( SPATIAL_DERIVATIVE num_xyz )
 {
    xyz_ = num_xyz;
 }
 
 
 
 
 
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNi_rhsop_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
   e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), op_vec_ );

} // end GetOperands

 
 
 
 
//element contribution
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNi_rhsop_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // initialize output matrix
    MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() );
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );

    // if the agregated finite element is a simplex, the Jacobian and element-interpolation derivative matrix is constant throughout it
    const bool is_simplex_element_type(e.FE()->IsSimplex() && e.Interpolation() == 1);
    double det = (is_simplex_element_type) ? e.dN_AtBaryCenter( DERIV_ ) : 0.;

    vector<double> IPOL;
    for ( auto i{0U}; i<e.IntegrationPoints(); i++ )
      {
         // computing gradient of operand
         if ( !is_simplex_element_type ) det = e.dN_AtIntegrationPoint( DERIV_, i );
         double grad_op(0.);
         const size_t nodes(e.Nodes());
         for ( uint32_t j{0U}; j<nodes; ++j )
           grad_op += DERIV_(xyz_,j) * op_vec_[j]();
        
         // integration
         e.N_AtIntegrationPoint( i, IPOL );
         for ( uint32_t k=0; k<nodes; ++k )
           IPOL[k] *= grad_op * det * e.WeightAtIntegrationPoint(i);

         // RHS vector for accumulation
         for ( uint32_t k=0; k<nodes; ++k )
           MathOperatorRHS<dim,CELL>::RHS[k] += IPOL[k];
      }

} // end ComputeContribution

// TESTING
//  cerr <<"\nelement "<< e.Idx() <<": rhs: ";
//  out( MathOperatorRHS<dim,CELL>::RHS );



template class NumIntegral_DNi_rhsop_dV<1U,Element>;
template class NumIntegral_DNi_rhsop_dV<2U,Element>;
template class NumIntegral_DNi_rhsop_dV<3U,Element>;

template class NumIntegral_DNi_rhsop_dV<1U,Face>;
template class NumIntegral_DNi_rhsop_dV<2U,Face>;
template class NumIntegral_DNi_rhsop_dV<3U,Face>;

} // namespace csmp
