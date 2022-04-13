#include "NumIntegral_DNi_rhsop_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim,class CELL>
NumIntegral_DNi_rhsop_dV<dim,CELL>::~NumIntegral_DNi_rhsop_dV() {}

/// custom constructor that should be used
template<uint32_t dim,class CELL>
NumIntegral_DNi_rhsop_dV<dim,CELL>::NumIntegral_DNi_rhsop_dV( const PropertyDatabase<dim>& pref,
                                                                 const char*                  oper,
                                                                 const char*                  test )
  : MathOperatorRHS<dim>(pref,oper,test),
    xyz_(Y_DIRECTION), op_vec_(dim)
{
    MathOperatorRHS<dim>::Name("NumIntegral_DNi_rhsop_dV",oper, test );
    

    if ( MathOperatorRHS<dim>::BasicOperandPlacement() != NODE ||
         MathOperatorRHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_DNi_rhsop_dV<dim>::(constructor)",
                             oper, "Operand must be a scalar property placed on the nodes." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim>::TestOperandType() != SCALAR)
      throw csmp::Exception( ERROR, "NumIntegral_DNi_rhsop_dV<dim>::(constructor)",
                             test, "Operand (test) must be a scalar property placed on the nodes." );
}




string parse( SPATIAL_DERIVATIVE deriv ) {
    if ( deriv == X_DIRECTION ) return "X-direction";
    if ( deriv == Y_DIRECTION ) return "Y-direction";
    if ( deriv == Z_DIRECTION ) return "Z-direction";
    return string("undefined");
 }





template<uint32_t dim,class CELL>
void NumIntegral_DNi_rhsop_dV<dim,CELL>::SpatialDerivative( SPATIAL_DERIVATIVE num_xyz )
 {
    xyz_ = num_xyz;
 }
 
 
 
 
 
template<uint32_t dim,class CELL>
void NumIntegral_DNi_rhsop_dV<dim,CELL>::GetOperands( const CELL& e )
{
   e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), op_vec_ );

} // end GetOperands

 
 
//element contribution
template<uint32_t dim,class CELL>
void NumIntegral_DNi_rhsop_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // initialize output matrix
    MathOperatorRHS<dim>::RHS.resize( e.Nodes() );
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    // if the agregated finite element is a simplex, the Jacobian and element-interpolation derivative matrix is constant throughout it
    const bool is_simplex_element_type(e.FE()->IsSimplex() && e.Interpolation() == 1);
    double det = (is_simplex_element_type) ? e.dN_AtBaryCenter( MathOperatorRHS<dim>::DERIV ) : 0.;

    for ( auto i{0U}; i<e.IntegrationPoints(); i++ )
      {
         // computing gradient of operand
         if ( !is_simplex_element_type ) det = e.dN_AtIntegrationPoint( MathOperatorRHS<dim>::DERIV, i );
         double grad_op(0.);
         const size_t nodes(e.Nodes());
         for ( auto j{0U}; j<nodes; ++j )
           grad_op += MathOperatorRHS<dim>::DERIV(xyz_,j) * op_vec_[j]();
        
         // integration
         e.N_AtIntegrationPoint( i, MathOperatorRHS<dim>::IPOL );
         for ( auto k=0; k<nodes; ++k )
           MathOperatorRHS<dim>::IPOL[k] *= grad_op * det * e.WeightAtIntegrationPoint(i);

         // RHS vector for accumulation
         for ( auto k=0; k<nodes; ++k )
           MathOperatorRHS<dim>::RHS[k] += MathOperatorRHS<dim>::IPOL[k];
      }

} // end ComputeContribution

// TESTING
//  cerr <<"\nelement "<< e.Idx() <<": rhs: ";
//  out( MathOperatorRHS<dim>::RHS );



template class NumIntegral_DNi_rhsop_dV<1U,Element<1U> >;
template class NumIntegral_DNi_rhsop_dV<2U,Element<2U> >;
template class NumIntegral_DNi_rhsop_dV<3U,Element<3U> >;

template class NumIntegral_DNi_rhsop_dV<1U,Face<1U> >;
template class NumIntegral_DNi_rhsop_dV<2U,Face<2U> >;
template class NumIntegral_DNi_rhsop_dV<3U,Face<3U> >;

} // namespace csmp
