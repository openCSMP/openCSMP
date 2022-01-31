#include "NumIntegral_NT_dNi_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class CELL>
NumIntegral_NT_dNi_dV<dim,CELL>::~NumIntegral_NT_dNi_dV() {}

template<size_t dim,class CELL>
NumIntegral_NT_dNi_dV<dim,CELL>::NumIntegral_NT_dNi_dV( const PropertyDatabase<dim>& pref,
                                                                 const char*            basic,
                                                                 const char*            test )
  : MathOperatorLHS<dim>(pref,basic,test),
    TEMP(3,3), xyz_(Y_DIRECTION), transp_(false)
{
    MathOperatorLHS<dim>::Name("NumIntegral_NT_dNi_dV", basic, test );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV<dim>::(constructor)",
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim>::TestOperandType() != SCALAR)
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV<dim>::(constructor)",
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}

/** 
     specify the Cartesian direction in which the derivative shall be taken.
*/
template<size_t dim,class CELL>
void NumIntegral_NT_dNi_dV<dim,CELL>::SpatialDerivative( SPATIAL_DERIVATIVE num_xyz )
 {
    xyz_ = num_xyz;
 }

/**
    Will cause the element matrix to be transposed; default is false.
*/
template<size_t dim,class CELL>
void NumIntegral_NT_dNi_dV<dim,CELL>::Transposed()
 {
    transp_ = true;
 }
 
 
 
/** 
    Computes IPOL * DN_i product weighted by the determinant of Jacobian matrix.
*/
template<size_t dim,class CELL>
void NumIntegral_NT_dNi_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();
    const size_t nodes(e.Nodes());
    this->IPOL.resize( nodes );
    TEMP.Resize( nodes, nodes );
   
    const size_t ipoints(e.IntegrationPoints());
    for ( size_t i=0; i<ipoints; i++ )
      {
         e.N_AtIntegrationPoint( i, this->IPOL );
         double det = e.dN_AtIntegrationPoint( this->DERIV, i );
         for ( size_t j=0; j<nodes; j++ )
            for ( size_t k=0; k<nodes; k++ )
              transp_ ? TEMP(j,k) = this->IPOL[j] * this->DERIV(xyz_,k) :
                        TEMP(j,k) = this->IPOL[k] * this->DERIV(xyz_,j);
        
         TEMP *= (det * e.WeightAtIntegrationPoint(i));
        
         MathOperatorLHS<dim>::LHS += TEMP;
      }

} // end ComputeContribution

template class NumIntegral_NT_dNi_dV<1U,Element<1U> >;
template class NumIntegral_NT_dNi_dV<2U,Element<2U> >;
template class NumIntegral_NT_dNi_dV<3U,Element<3U> >;

template class NumIntegral_NT_dNi_dV<1U,Face<1U> >;
template class NumIntegral_NT_dNi_dV<2U,Face<2U> >;
template class NumIntegral_NT_dNi_dV<3U,Face<3U> >;

} //namespace csmp
