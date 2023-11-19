//
//  DiffusionLHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 5/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "DiffusionLHS.h"
#include "Element.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
DiffusionLHS<dim>::DiffusionLHS( const INDEX<SCALAR,ELEMENT>& diffusivity_key, 
                                 const INDEX<SCALAR,NODE>& transported_variable_key  )
  : diff_key_(diffusivity_key),
    adv_key_(transported_variable_key),
    dof_(1U) // components of the diffusivity
 {
 }

// NOT AVAILABLE
template<uint32_t dim>
void DiffusionLHS<dim>::AccumulateFiniteVolume( const Node<dim>&, SparseMatrix& ) const
 {
     throw csmp::Exception( ERROR, "DiffusionLHS<dim>::AccumulateStencil", "Use FEM diffusion operator instead!" );
 }


/**
     Uses the finite element framework of CSMP to compute the FE diffusion integral term and
     accumulate it into the solution matrix A of Ax=b.
*/
template<uint32_t dim>
void DiffusionLHS<dim>::AccumulateStencil( const Element<dim>& fe, SparseMatrix& A ) const
 {
    // assuming a scalar diffusivity
    const double  diffusion_coeff = fe.Read( diff_key_ );
   
    // if there is only a single Jacobian needed because the element is a simplex
    double detJ(0U);
    if ( fe.FE()->IsSimplex() ) {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         detJ = fe.dN_AtBaryCenter( DN_, dof_ );
         DN_.Transposed( DNT_ );
         DN_  *= diffusion_coeff;
         DNT_ *= DN_;
      }
 
    const size_t integration_points(fe.IntegrationPoints());
    for ( auto i{0U}; i<integration_points; ++i )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         if ( !fe.FE()->IsSimplex() ) {
              detJ = fe.dN_AtIntegrationPoint( DN_, i, dof_ );
              DN_.Transposed( DNT_ );
              DN_  *= diffusion_coeff;
              DNT_ *= DN_;
              RESULT_ = DNT_;
           }
         RESULT_  = DNT_;
         RESULT_ *= fe.WeightAtIntegrationPoint(i) * detJ;

         // assigning the matrix contribution to the solution matrix
         for ( uint32_t j{0U}; j<DNT_.Rows(); j++ )
           for ( uint32_t k{0U}; k<DNT_.Cols(); k++ )
             A.Add( fe.N(j)->Idx(), fe.N(k)->Idx(), RESULT_(j,k) );
      }

 } // end AccumulateStencil
 

template class DiffusionLHS<1U>;
template class DiffusionLHS<2U>;
template class DiffusionLHS<3U>;

} // end csmp
