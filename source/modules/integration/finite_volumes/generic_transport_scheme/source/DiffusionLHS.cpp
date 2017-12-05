//
//  DiffusionLHS.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 5/12/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "DiffusionLHS.h"

namespace csmp {

/**
     Uses finite element framework to compute the diffusion term and 
     accumulate late it into the solution matrix A of Ax=b.
*/
template<size_t dim>
void DiffusionLHS<dim>::AccumulateStencil( VariableSet_TracerTransferImplicit& props, Element<dim>& e, SparseMatrix& A )
 {
    // assuming a scalar diffusivity
    double64  diffusion_coeff = e.Read( props.diff_key );
 
    const size_T integration_points(e.IntegrationPoints());
    for ( size_t i=0U; i<integration_points; ++i )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         double64 detJ = e.dN_AtIntegrationPoint( DN_, i, SCALAR );
         DN_.Transposed( DNT_ );
         DN_  *= diffusion_coeff;
         DNT_ *= DN;
         DNT_ *= e.WeightAtIntegrationPoint(i) * detJ;

         // assigning the matrix contribution to the solution matrix
         for ( size_t j=0U; j<DNT.Rows(); j++ )
           for ( size_t k=0U; k<DNT.Cols(); k++ )
             LHS.Add( e->N(j)->Idx(),
                      e->N(k)->Idx(), DNT(j,k) );
      }

 }
 
 
 
 template class DiffusionLHS<3U>;




} // end csmp
