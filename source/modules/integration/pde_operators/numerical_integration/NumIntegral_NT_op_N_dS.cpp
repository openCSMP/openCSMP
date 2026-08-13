#include "NumIntegral_NT_op_N_dS.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
NumIntegral_NT_op_N_dS<dim>::NumIntegral_NT_op_N_dS( const PropertyDatabase<dim>& pref,
                                                     const char* oper, const char* test )
  : MathOperatorRHS<dim,Face>(pref,oper,test),
    nodal_degrees_of_freedom(1)
 {
    MathOperatorRHS<dim,Face>::Name("NumIntegral_NT_op_N_dS", oper, test );

    // anisotropy can only be considered if there are multiple degrees of freedom per node
    if ( MathOperatorRHS<dim,Face>::MaterialOperandType() != SCALAR ) nodal_degrees_of_freedom = dim;
    if ( MathOperatorRHS<dim,Face>::MaterialOperandType() != SCALAR &&
         MathOperatorRHS<dim,Face>::TestOperandType() == SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op_N_dS<dim>::(constructor)", 
                             test, "Non-scalar Operands require multiple degrees of freedom per node." );

    if ( MathOperatorRHS<dim,Face>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,Face>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op_N_dS<dim>::(constructor)", 
                             test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
Computes the volume integral over the element interpolation functions 
times the Operand. If the Operand is 1 over the element, then the volume
integral is naturally 1 as well.  

@attention accumulation takes place only for faces where the Material Operand is flagged Neumann.

*/
template<uint32_t dim>
void NumIntegral_NT_op_N_dS<dim>::ComputeContribution( const Face<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorRHS<dim,Face>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,Face>::RHS.begin(), MathOperatorRHS<dim,Face>::RHS.end(), 0. );

    // lumped formulation: only the midside nodes are used in the lumped approach
    if ( MathOperatorRHS<dim,Face>::LumpedFormulation() )
      {
         const double area(e.Area()); // NT . N
         
         if ( MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == ELEMENT ||
              MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == FACE ||
              MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == INTER_FACE )
           {
             for ( uint32_t j{0U}; j<e.Nodes(); j++ )
               MathOperatorRHS<dim,Face>::RHS[j] =
                 // after having ascertained that the material property is a scalar
                 (MathOperatorRHS<dim,Face>::MTRL[0](0,0) * area) / static_cast<double>(e.Nodes());
           }
         else if ( MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == NODE ||
                   MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
           {
             for ( uint32_t j{0U}; j<e.Nodes(); j++ )
               MathOperatorRHS<dim,Face>::RHS[j] =
                 (MathOperatorRHS<dim,Face>::MTRL[j](0,0) * area) / static_cast<double>(e.Nodes());
           }
      }  

    // consistent formulation
    else
      {
         RHS_TEMP.Resize(e.Nodes(), e.Nodes());
         RHS_TEMP.Zero();
            
         for ( uint32_t i{0U}; i < e.IntegrationPoints(); i++ )
           {
              e.N_AtIntegrationPoint( i, e.FE()->NRST );
              const double det(e.det_J_AtIntegrationPoint( i ));
              // forming NT * mtrl
              NT.Resize(e.Nodes(),1U);
              if ( MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == ELEMENT or
                   MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == REGION  or
                   MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == INTER_FACE or
                   MathOperatorRHS<dim,Face>::MaterialOperandPlacement() == FACE )
                for ( uint32_t j{0U}; j<e.Nodes(); j++ ) NT(j,0U) = this->MTRL[0](0,0) * e.FE()->NRST[j];
              else // node or integration point                               ^^^
                for ( uint32_t j{0U}; j<e.Nodes(); j++ ) NT(j,0U) = this->MTRL[i](0,0) * e.FE()->NRST[j];
              // forming Wj * detJ * N                                        ^^^
              N.Resize(1U,e.Nodes());
              for ( uint32_t j{0U}; j<e.Nodes(); j++ )
                N(0U,j) = e.FE()->NRST[j] * det * e.WeightAtIntegrationPoint(i);
              
              // forming the mass matrix NT op N
              RHS_TEMP += NT * N;
           }

         // row-sum diagonalisation of matrix RHS_TEMP and addition to righthand vector
         fill( MathOperatorRHS<dim,Face>::RHS.begin(), MathOperatorRHS<dim,Face>::RHS.end(), 0. );
         for ( uint32_t j{0U}; j<e.Nodes(); j++ )
           for ( uint32_t k{0u}; k<e.Nodes(); k++ )
             MathOperatorRHS<dim,Face>::RHS[j] += RHS_TEMP(j,k);
      }
 
   
} // end ComputeContribution

//cout <<"\nFace "<< e.Idx() <<": NumIntegral_NT_op_N_dS="<< endl;
//out( MathOperatorRHS<dim,CELL>::RHS );
    

template class NumIntegral_NT_op_N_dS<1U>;
template class NumIntegral_NT_op_N_dS<2U>;
template class NumIntegral_NT_op_N_dS<3U>;

} // csmp
