// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_NT_dNi_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_dNi_dV<dim,CELL>::NumIntegral_NT_dNi_dV( const PropertyDatabase<dim>& pref,
                                                         const char*            basic,
                                                         const char*            test )
  : MathOperatorLHS<dim,CELL>(pref,basic,test),
    TEMP(3,3), xyz_(Y_DIRECTION), transp_(false)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_NT_dNi_dV", basic, test );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV<dim>::(constructor)",
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR)
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV<dim>::(constructor)",
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}




/** 
     specify the Cartesian direction in which the derivative shall be taken.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_dNi_dV<dim,CELL>::SpatialDerivative( SPATIAL_DERIVATIVE num_xyz )
 {
    xyz_ = num_xyz;
 }




/**
    Will cause the element matrix to be transposed; default is false.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_dNi_dV<dim,CELL>::Transposed()
 {
    transp_ = true;
 }
 
 
 
/** 
    Computes N_transposed * DN_i product that is weighted by the determinant of Jacobian matrix.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_dNi_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );
    // making sure that is operator is only applied to volume elmts in 3D, surf in 2D and line elmts in 1D
    assert( e.IsEquidimensional() );

    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();
    const uint32_t nodes(e.Nodes());
    this->IPOL.resize( nodes );
    TEMP.Resize( nodes, nodes );
   
    const size_t ipoints(e.IntegrationPoints());
    for ( uint32_t i{0U}; i<ipoints; i++ )
      {
         e.N_AtIntegrationPoint( i, this->IPOL );
         double det = e.dN_AtIntegrationPoint( this->DERIV, i );
         for ( uint32_t j{0U}; j<nodes; j++ )
            for ( uint32_t k=0; k<nodes; k++ )
              transp_ ? TEMP(j,k) = this->IPOL[j] * this->DERIV(xyz_,k) :
                        TEMP(j,k) = this->IPOL[k] * this->DERIV(xyz_,j);
        
         TEMP *= (det * e.WeightAtIntegrationPoint(i));
        
         MathOperatorLHS<dim,CELL>::LHS += TEMP;
      }

} // end ComputeContribution

template class NumIntegral_NT_dNi_dV<1U,Element>;
template class NumIntegral_NT_dNi_dV<2U,Element>;
template class NumIntegral_NT_dNi_dV<3U,Element>;

template class NumIntegral_NT_dNi_dV<1U,Face>;
template class NumIntegral_NT_dNi_dV<2U,Face>;
template class NumIntegral_NT_dNi_dV<3U,Face>;

} //namespace csmp
