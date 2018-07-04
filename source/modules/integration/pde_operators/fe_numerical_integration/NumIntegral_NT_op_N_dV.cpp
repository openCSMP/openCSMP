#include "NumIntegral_NT_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

template<size_t dim,class CELL>
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
template<size_t dim,class CELL>
void NumIntegral_NT_op_N_dV<dim,CELL>::ComputeContribution( CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    // lumped formulation: only the midside nodes are used in the lumped approach
    if ( MathOperatorRHS<dim>::LumpedFormulation() ) 
      {
         const double64 volume(e.Volume()); // NT . N
         
         if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ||
              MathOperatorRHS<dim>::MaterialOperandPlacement() == FACE ||
              MathOperatorRHS<dim>::MaterialOperandPlacement() == INTER_FACE)
           {
             for ( size_t j=0U; j<e.Nodes(); j++ )
               MathOperatorRHS<dim>::RHS[j] = 
                 (MathOperatorRHS<dim>::MTRL[0](0,0)*volume) / static_cast<double64>(e.Nodes());
           }
         else if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE ||  
                   MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT)
           {
             for ( size_t j=0U; j<e.Nodes(); j++ )
               MathOperatorRHS<dim>::RHS[j] = 
                 (MathOperatorRHS<dim>::MTRL[j](0,0)*volume) / static_cast<double64>(e.Nodes());
           }
      }  

    // consistent formulation
    else
      {
         RHS_TEMP.Resize(e.Nodes(), e.Nodes());
         RHS_TEMP.Zero();
            
         for ( size_t i=0U; i < e.FE()->IntegrationPoints(); i++ )
           {
              e.N_AtIntegrationPoint( i, e.FE()->NRST );
              const double64 det(e.det_JINV_AtIntegrationPoint( i ));
              // forming NT * mtrl
              NT.Resize(e.Nodes(),1U);
              if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT or
                   MathOperatorRHS<dim>::MaterialOperandPlacement() == REGION  or
                   MathOperatorRHS<dim>::MaterialOperandPlacement() == FACE ) 
                for ( size_t j=0U; j<e.Nodes(); j++ ) NT(j,0U) = this->MTRL[0](0,0) * e.FE()->NRST[j];
              else // node or integration point                            ^^^
                for ( size_t j=0U; j<e.Nodes(); j++ ) NT(j,0U) = this->MTRL[i](0,0) * e.FE()->NRST[j];
              // forming Wj * detJ * N                                                 ^^^
              N.Resize(1U,e.Nodes());
              for ( size_t j=0U; j<e.Nodes(); j++ ) N(0U,j) = e.FE()->NRST[j] * det * 
                                               e.WeightAtIntegrationPoint(i);;
              
              // forming the mass matrix NT op N
              RHS_TEMP += NT * N;
           }

         // row-sum diagonalisation of matrix RHS_TEMP and addition to righthand vector
         fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
         for ( size_t j=0; j<e.Nodes(); j++ )
           for ( size_t k=0; k<e.Nodes(); k++ ) 
             MathOperatorRHS<dim>::RHS[j] += RHS_TEMP(j,k);
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
