#include "NumIntegral_NT_lhsop_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
NumIntegral_NT_lhsop_N_dV<dim,SIMPLEX>::NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                           const char* oper,
                                                           const char* basic,
                                                           const char* test )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    TEMP(3,3)
 {
    MathOperatorLHS<dim>::Name("NumIntegral_NT_lhsop_N_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)", 
                      oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)", 
                      basic, "Weighting variable must be a scalar property placed on the nodes." );
 }







/**
 
Computes the volume (area) integral over the testfunction products 
multiplied with the Operand. As an example the capacitance matrix, C,
may be used which is usually based on the storativity, S:  

C[n x n] = w(i) * { N[n x 1] * S[1 x 1] * N[1 x n] * |J[dim x dim]| } 

see, for instance, J.Istock p. 132.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_NT_lhsop_N_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS.Zero();

    if ( MathOperatorLHS<dim>::LumpedFormulation() )
      {
         const double64 volume = e.Volume();

         if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ||
              MathOperatorLHS<dim>::MaterialOperandPlacement() == FACE ||
              MathOperatorLHS<dim>::MaterialOperandPlacement() == INTER_FACE)
           {
            for ( size_t j=0; j<e.Nodes(); j++ )
                this->LHS(j,j) = (this->MTRL[0](0,0)*volume) / static_cast<double64>(e.Nodes());
           }
         else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == NODE ||  
                   MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT)
           {
            for ( size_t j=0; j<e.Nodes(); j++ )
                this->LHS(j,j) = (this->MTRL[j](0,0)*volume) / static_cast<double64>(e.Nodes());
           }
      }  

    // consistent formulation
    else 
      {
         TEMP.Resize( e.Nodes(), e.Nodes() );    

         // consistent formulation
         for ( size_t i=0U; i < e.FE()->IntegrationPoints(); i++ )
           {
              e.N_AtIntegrationPoint( i, e.FE()->NRST );
              const double64 det = e.det_JINV_AtIntegrationPoint( i );
                   
              for ( size_t j=0U; j<e.Nodes(); j++ )
                for ( size_t k=0U; k<e.Nodes(); k++ )
                  {
                    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ||
                         MathOperatorLHS<dim>::MaterialOperandPlacement() == FACE ||
                         MathOperatorLHS<dim>::MaterialOperandPlacement() == INTER_FACE)
                      TEMP(j,k) = e.FE()->NRST[j] * MathOperatorLHS<dim>::MTRL[0](0,0) * e.FE()->NRST[k];
                    else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == NODE ||
                              MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT)
                      TEMP(j,k) = e.FE()->NRST[j] * MathOperatorLHS<dim>::MTRL[i](0,0) * e.FE()->NRST[k];
                    //else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
                      //TEMP(j,k) = e.FE()->NRST[j] * MathOperatorLHS<dim>::MTRL[i](0,0) * e.FE()->NRST[k];

                  }
             
              TEMP *= (det * e.WeightAtIntegrationPoint(i));
             
              MathOperatorLHS<dim>::LHS += TEMP;
           }
      } 

} // end ComputeContribution

//cout <<"\nElement "<< e.Idx() <<": NumIntegral_NT_lhsop_N_dV="<< endl;
//MathOperatorLHS<dim>::LHS.Out();



template class NumIntegral_NT_lhsop_N_dV<1U,Element<1U> >;
template class NumIntegral_NT_lhsop_N_dV<2U,Element<2U> >;
template class NumIntegral_NT_lhsop_N_dV<3U,Element<3U> >;

template class NumIntegral_NT_lhsop_N_dV<1U,Face<1U> >;
template class NumIntegral_NT_lhsop_N_dV<2U,Face<2U> >;
template class NumIntegral_NT_lhsop_N_dV<3U,Face<3U> >;

template class NumIntegral_NT_lhsop_N_dV<1U,InterFace<1U> >;
template class NumIntegral_NT_lhsop_N_dV<2U,InterFace<2U> >;
template class NumIntegral_NT_lhsop_N_dV<3U,InterFace<3U> >;

} // csmp
