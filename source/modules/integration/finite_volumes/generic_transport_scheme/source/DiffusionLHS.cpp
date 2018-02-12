//
//  DiffusionLHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 5/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "DiffusionLHS.h"
#include "Element.h"
#include "Model.h"

using namespace std;

namespace csmp {

template<size_t dim>
DiffusionLHS<dim>::DiffusionLHS( const Model<dim>& model, const char* diffusivity )
  : MatrixOperator<dim>(0),
    diff_key_(model.Database().StorageKey(diffusivity)),
    dof_(model.Database().Components(diffusivity))
 {
 }



/**
     Uses the finite element framework of CSMP to compute the FE diffusion integral term and
     accumulate it into the solution matrix A of Ax=b.
*/
template<size_t dim>
void DiffusionLHS<dim>::AccumulateStencil( Element<dim>& fe, SparseMatrix& A ) const
 {
    // assuming a scalar diffusivity
    double64  diffusion_coeff = fe.Read( diff_key_ );
   
    // if there is only a single Jacobian needed because the element is a simplex
    double64 detJ(0U);
    if ( fe.FE()->IsSimplex() ) {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         detJ = fe.dN_AtBaryCenter( DN_, dof_ );
         DN_.Transposed( DNT_ );
         DN_  *= diffusion_coeff;
         DNT_ *= DN_;
      }
 
    const size_t integration_points(fe.IntegrationPoints());
    for ( size_t i=0U; i<integration_points; ++i )
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
         for ( size_t j=0U; j<DNT_.Rows(); j++ )
           for ( size_t k=0U; k<DNT_.Cols(); k++ )
             A.Add( fe.N(j)->Idx(), fe.N(k)->Idx(), RESULT_(j,k) );
      }

 }
 

template class DiffusionLHS<1U>;
template class DiffusionLHS<2U>;
template class DiffusionLHS<3U>;

} // end csmp
