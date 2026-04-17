#include "NumIntegral_PT_op_P_dV.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_PT_op_P_dV<dim,CELL>::NumIntegral_PT_op_P_dV( const PropertyDatabase<dim>& pref,
                                                          const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    nodal_degrees_of_freedom(1)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_PT_op_P_dV", oper, test );

    // resize material property matrix 
    if constexpr ( dim == 3U )
      {
         for ( auto it=MathOperatorRHS<dim,CELL>::MTRL.begin();
               it!=MathOperatorRHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    // anisotropy can only be considered if there are multiple degrees of freedom per node
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR ) nodal_degrees_of_freedom = dim;
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR && 
         MathOperatorRHS<dim,CELL>::TestOperandType() == SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_P_dV<dim>::(constructor)", 
                      test, "Non-scalar Operands require multiple degrees of freedom per node." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_P_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
Computes the volume integral over the element interpolation functions 
times the Operand. If the Operand is 1 over the element, then the volume
integral is naturally 1 as well.  

in the special case, the integral of the testfunction products must be used


@section application Application

Use this operator to compute capacitance, storage capacity etc. matrices.
 */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_PT_op_P_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );
    
    DenseMatrix<DM_MIN> RHS_TEMP(e.Nodes(), e.Nodes()), UNITY(e.Nodes(), 1);
    
    vector<double>  N( e.Nodes() );
    double          det( 0.0 );
    double          volume( 0.0 ); // NT.N = element volume
    
    UNITY = 1.0;
    
    // if the property is constant over the element the lumped formulation is the
    // same as the consistent formulation, where
    // the element contribution which is uniform over the entire element is just distributed
    // equally over all the nodes.
       // lumped formulation
       if ( MathOperatorRHS<dim,CELL>::LumpedFormulation() )
         {  
           volume = e.Volume();
             
           uint32_t k{0};
             if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
               {   
                 for ( uint32_t n=0; n<e.Nodes(); n++ ) 
                   for ( uint32_t i{0U}; i<nodal_degrees_of_freedom; i++ )
                     MathOperatorRHS<dim,CELL>::RHS[k++] = (MathOperatorRHS<dim,CELL>::MTRL[0](i,i) * volume) / e.Nodes();
               }
             else
               {
                 cout << "\nNumIntegral_PT_op_P_dV<dim>::ComputeContribution: Lumped formulation for property placement ";
                 cout << "\nother than elements not implemented yet!" << endl;
               }

         }
         
        // TODO: consistent formulation (not tested thus far)
       else 
         {
            for ( uint32_t i = 0; i<e.FE()->IntegrationPoints(); i++ )
              {
                e.N_AtIntegrationPoint( i, N );
                det = e.det_JINV_AtIntegrationPoint( i );

                for ( uint32_t j = 0; j < e.Nodes(); j++ )
                  for ( uint32_t k = 0; k < e.Nodes(); k++ )
                    {
                      if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
                        RHS_TEMP(j,k) = N[j] * MathOperatorRHS<dim,CELL>::MTRL[0](0,0) * N[k] * det;
                      else if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE )
                        RHS_TEMP(j,k) = N[j] * MathOperatorRHS<dim,CELL>::MTRL[i](0,0) * N[k] * det;
                    }
                for ( uint32_t j = 0; j < e.Nodes(); j++ )   
                  RHS_TEMP(j,i) *= e.WeightAtIntegrationPoint(i);
                 
                RHS_TEMP *= UNITY; 
                 
                for ( uint32_t l = 0; l < e.Nodes(); l++ )
                   MathOperatorRHS<dim,CELL>::RHS[l] += RHS_TEMP(l,1);
                
                RHS_TEMP.Resize( e.Nodes(), e.Nodes() );
              }
         }
      
} // end ComputeContribution


template class NumIntegral_PT_op_P_dV<1U,Element>;
template class NumIntegral_PT_op_P_dV<2U,Element>;
template class NumIntegral_PT_op_P_dV<3U,Element>;

template class NumIntegral_PT_op_P_dV<1U,Face>;
template class NumIntegral_PT_op_P_dV<2U,Face>;
template class NumIntegral_PT_op_P_dV<3U,Face>;

} // csmp
