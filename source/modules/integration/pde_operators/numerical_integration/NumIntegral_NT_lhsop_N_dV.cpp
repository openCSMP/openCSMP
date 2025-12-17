#include "NumIntegral_NT_lhsop_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"

using namespace std;

namespace csmp {


template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_lhsop_N_dV<dim,CELL>::NumIntegral_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                               const char* oper,
                                                               const char* basic,
                                                               const char* test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    TEMP(3,3)
 {
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_NT_lhsop_N_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)", 
                      oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhsop_N_dV<dim>::(constructor)", 
                      basic, "Weighting variable must be a scalar property placed on the nodes." );
 }







/**
 
Computes the volume (area) integral over the testfunction products 
multiplied with the Operand. As an example the capacitance matrix, C,
may be used which is usually based on the storativity, S:  

C[n x n] = w(i) * { N[n x 1] * S[1 x 1] * N[1 x n] * |J[dim x dim]| } 

see, for instance, J.Istock p. 132.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_lhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    if ( MathOperatorLHS<dim,CELL>::LumpedFormulation() )
      {
         const double volume = e.Volume();

         if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
              MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == FACE ||
              MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == INTER_FACE)
           {
            for ( auto j{0U}; j<e.Nodes(); j++ )
                this->LHS(j,j) = (this->MTRL[0](0,0)*volume) / static_cast<double>(e.Nodes());
           }
         else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE ||
                   MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT)
           {
            for ( auto j{0U}; j<e.Nodes(); j++ )
                this->LHS(j,j) = (this->MTRL[j](0,0)*volume) / static_cast<double>(e.Nodes());
           }
      }  

    // consistent formulation
    else 
      {
         TEMP.Resize( e.Nodes(), e.Nodes() );    

         // consistent formulation
         for ( uint32_t i{0}; i < e.FE()->IntegrationPoints(); i++ )
           {
              e.N_AtIntegrationPoint( i, e.FE()->NRST );
              const double det = e.det_JINV_AtIntegrationPoint( i );
                   
              for ( uint32_t j{0U}; j<e.Nodes(); j++ )
                for ( uint32_t k{0U}; k<e.Nodes(); k++ )
                  {
                    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
                         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == FACE ||
                         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == INTER_FACE)
                      TEMP(j,k) = e.FE()->NRST[j] * MathOperatorLHS<dim,CELL>::MTRL[0](0,0) * e.FE()->NRST[k];
                    else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE ||
                              MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT)
                      TEMP(j,k) = e.FE()->NRST[j] * MathOperatorLHS<dim,CELL>::MTRL[i](0,0) * e.FE()->NRST[k];
                    //else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
                      //TEMP(j,k) = e.FE()->NRST[j] * MathOperatorLHS<dim>::MTRL[i](0,0) * e.FE()->NRST[k];

                  }
             
              TEMP *= (det * e.WeightAtIntegrationPoint(i));
             
              MathOperatorLHS<dim,CELL>::LHS += TEMP;
           }
      } 

} // end ComputeContribution

//cout <<"\nElement "<< e.Idx() <<": NumIntegral_NT_lhsop_N_dV="<< endl;
//MathOperatorLHS<dim>::LHS.Out();



template class NumIntegral_NT_lhsop_N_dV<1U,Element>;
template class NumIntegral_NT_lhsop_N_dV<2U,Element>;
template class NumIntegral_NT_lhsop_N_dV<3U,Element>;

template class NumIntegral_NT_lhsop_N_dV<1U,Face>;
template class NumIntegral_NT_lhsop_N_dV<2U,Face>;
template class NumIntegral_NT_lhsop_N_dV<3U,Face>;

template class NumIntegral_NT_lhsop_N_dV<1U,InterFace>;
template class NumIntegral_NT_lhsop_N_dV<2U,InterFace>;
template class NumIntegral_NT_lhsop_N_dV<3U,InterFace>;

} // csmp
