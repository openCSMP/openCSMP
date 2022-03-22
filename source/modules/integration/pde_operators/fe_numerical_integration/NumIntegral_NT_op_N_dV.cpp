#include "NumIntegral_NT_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

template<uint32_t dim,class CELL>
NumIntegral_NT_op_N_dV<dim,CELL>::NumIntegral_NT_op_N_dV( const PropertyDatabase<dim>& pref,
                                                          const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    nodal_degrees_of_freedom(1)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_NT_op_N_dV", oper, test );

    // anisotropy can only be considered if there are multiple degrees of freedom per node
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR ) nodal_degrees_of_freedom = dim;
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR && 
         MathOperatorRHS<dim>::TestOperandType() == SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op_N_dV<dim>::(constructor)", 
                             test, "Non-scalar Operands require multiple degrees of freedom per node." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op_N_dV<dim>::(constructor)", 
                             test, "Operand (test) must be a scalar property placed on the nodes." );

    // resize material property matrix 
    for ( typename vector<DenseMatrix<DM_MIN> >::iterator
          it=MathOperatorRHS<dim>::MTRL.begin(); it!=MathOperatorRHS<dim>::MTRL.end(); it++ ) 
      (*it).Resize(dim,dim);
}





/**
 
Computes the volume integral over the element interpolation functions 
times the Operand. If the Operand is 1 over the element, then the volume
integral is naturally 1 as well.  
*/
template<uint32_t dim,class CELL>
void NumIntegral_NT_op_N_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.UsesLocalCoordinates() == true );

    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    const bool piecewise_constant_material( this->MaterialOperandPlacement() == ELEMENT ||
                                            this->MaterialOperandPlacement() == FACE    ||
                                            this->MaterialOperandPlacement() == REGION );

    const uint32_t n_nodes{ e.Nodes() }, n_ipoints{ e.IntegrationPoints() };

    // lumped formulation: only the midside nodes are used in the lumped approach
    if ( MathOperatorRHS<dim>::LumpedFormulation() ) 
      {
         const double volume(e.Volume()); // NT . N
         
         if ( piecewise_constant_material )
           {
             for ( auto j=0U; j<n_nodes; j++ )
               MathOperatorRHS<dim>::RHS[j] = 
                 (MathOperatorRHS<dim>::MTRL[0](0,0)*volume) / static_cast<double>(e.Nodes());
           }
         else if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE ||  
                   MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
           {
             for ( auto j=0U; j<n_nodes; j++ )
               MathOperatorRHS<dim>::RHS[j] = 
                 (MathOperatorRHS<dim>::MTRL[j](0,0)*volume) / static_cast<double>(e.Nodes());
           }
      }  

    // consistent formulation
    else
      {
         RHS_TEMP_.Resize(e.Nodes(), e.Nodes());
         RHS_TEMP_.Zero();
            
         // if the finite element is a linear simplex element, its Jacobian and element-interpolation derivative matrix
         // are constant throughout it
         const bool is_simplex_element_type(e.FE()->IsSimplex() && e.Interpolation() == 1 );
       
         if ( is_simplex_element_type ) {
              // initialising the interpolation function matrix
              e.N_AtBaryCenter( e.FE()->NRST );
              const double detJ = e.det_JINV_AtIntegrationPoint(0);
              // forming NT * mtrl
              NT_.Resize(n_nodes,1U);
              if ( piecewise_constant_material )
                for ( auto j=0U; j<n_nodes; j++ ) NT_(j,0U) = this->MTRL[0](0,0) * e.FE()->NRST[j];
              else {// node or integration point
                  for ( auto i{0}; i < n_ipoints; i++ )
                    for ( auto j=0U; j<n_nodes; j++ ) NT_(j,0U) = this->MTRL[i](0,0) * e.FE()->NRST[j];
                }
              // forming Wj * detJ * N
              N_.Resize(1U,n_nodes);
              for ( auto j=0U; j<n_nodes; j++ )
                N_(0U,j) = e.FE()->NRST[j] * detJ * e.WeightAtIntegrationPoint(0);
              
              // forming the mass matrix NT op N
              RHS_TEMP_ = NT_ * N_;
              return;
           }


         // if the element is not a simplex
         for ( auto i{0}; i < n_ipoints; i++ )
           {
              e.N_AtIntegrationPoint( i, e.FE()->NRST );
              const double det(e.det_JINV_AtIntegrationPoint( i ));
              // forming NT * mtrl
              NT_.Resize(n_nodes,1U);
              if ( piecewise_constant_material )
                for ( auto j=0U; j<n_nodes; j++ ) NT_(j,0U) = this->MTRL[0](0,0) * e.FE()->NRST[j];
              else // node or integration point                            ^^^
                for ( auto j=0U; j<n_nodes; j++ ) NT_(j,0U) = this->MTRL[i](0,0) * e.FE()->NRST[j];
              // forming Wj * detJ * N                                                     ^^^
              N_.Resize(1U,n_nodes);
              for ( auto j=0U; j<n_nodes; j++ ) N_(0U,j) = e.FE()->NRST[j] * det *
                                                           e.WeightAtIntegrationPoint(i);
              
              // forming the mass matrix NT op N
              RHS_TEMP_ += NT_ * N_;
           }

         // row-sum diagonalisation of matrix RHS_TEMP and addition to righthand vector
         fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
         for ( auto j=0; j<n_nodes; j++ )
           for ( auto k=0; k<n_nodes; k++ )
             MathOperatorRHS<dim>::RHS[j] += RHS_TEMP_(j,k);
      }
   
} // end ComputeContribution


//cout <<"\nElement "<< e.Idx() <<": NumIntegral_NT_op_N_dV="<< endl;
//out( MathOperatorRHS<dim>::RHS );

    

template class NumIntegral_NT_op_N_dV<1U,Element<1U> >;
template class NumIntegral_NT_op_N_dV<2U,Element<2U> >;
template class NumIntegral_NT_op_N_dV<3U,Element<3U> >;

template class NumIntegral_NT_op_N_dV<1U,Face<1U> >;
template class NumIntegral_NT_op_N_dV<2U,Face<2U> >;
template class NumIntegral_NT_op_N_dV<3U,Face<3U> >;

template class NumIntegral_NT_op_N_dV<1U,InterFace<1U> >;
template class NumIntegral_NT_op_N_dV<2U,InterFace<2U> >;
template class NumIntegral_NT_op_N_dV<3U,InterFace<3U> >;

} // csmp
