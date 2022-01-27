#include "NumIntegral_PT_op_P_dV.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<size_t dim,class CELL>
NumIntegral_PT_op_P_dV<dim,CELL>::NumIntegral_PT_op_P_dV( const PropertyDatabase<dim>& pref,
                                                     const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    nodal_degrees_of_freedom(1)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_PT_op_P_dV", oper, test );

    // resize material property matrix 
    if ( dim == 3 ) 
      {
         typename vector<DenseMatrix<DM_MIN> >::iterator  it;
         for ( it=MathOperatorRHS<dim>::MTRL.begin(); 
               it!=MathOperatorRHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    // anisotropy can only be considered if there are multiple degrees of freedom per node
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR ) nodal_degrees_of_freedom = dim;
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR && 
         MathOperatorRHS<dim>::TestOperandType() == SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_P_dV<dim>::(constructor)", 
                      test, "Non-scalar Operands require multiple degrees of freedom per node." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
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
template<size_t dim,class CELL>
void NumIntegral_PT_op_P_dV<dim,CELL>::ComputeContribution( const CELL& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0.0 );
    
    DenseMatrix<DM_MIN> RHS_TEMP(e.Nodes(), e.Nodes()), UNITY(e.Nodes(), 1);
    
    vector<double>  N( e.Nodes() );
    double          det( 0.0 );
    double          volume( 0.0 ); // NT.N = element volume
    size_t   k(0);
    
    UNITY = 1.0;
    
    // if the property is constant over the element the lumped formulation is the
    // same as the consistent formulation, where
    // the element contribution which is uniform over the entire element is just distributed
    // equally over all the nodes.
       // lumped formulation
       if ( MathOperatorRHS<dim>::LumpedFormulation() )
         {  
             volume = e.Volume();
             
             if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT )
               {   
                 for ( size_t n=0; n<e.Nodes(); n++ ) 
                   for ( size_t i=0; i<nodal_degrees_of_freedom; i++ )
                     MathOperatorRHS<dim>::RHS[k++] = (MathOperatorRHS<dim>::MTRL[0](i,i) * volume) / e.Nodes();
               }
             else
               {
                 cout << "\nNumIntegral_PT_op_P_dV<dim>::ComputeContribution: Lumped formulation for property placement ";
                 cout << "\nother than elements not implemented yet!" << endl;
               }

         }
         
        // consistent formulation (not tested thus far)  
       else 
         {
            for ( size_t i = 0; i<e.FE()->IntegrationPoints(); i++ )
              {
                e.N_AtIntegrationPoint( i, N );
                det = e.det_JINV_AtIntegrationPoint( i );

                for ( size_t j = 0; j < e.Nodes(); j++ )
                  for ( size_t k = 0; k < e.Nodes(); k++ )
                    {
                      if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT )
                        RHS_TEMP(j,k) = N[j] * MathOperatorRHS<dim>::MTRL[0](0,0) * N[k] * det;
                      else if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE )
                        RHS_TEMP(j,k) = N[j] * MathOperatorRHS<dim>::MTRL[i](0,0) * N[k] * det;
                    }
                for ( size_t j = 0; j < e.Nodes(); j++ )   
                  RHS_TEMP(j,k) *= e.WeightAtIntegrationPoint(i);
                 
                RHS_TEMP *= UNITY; 
                 
                for ( size_t l = 0; l < e.Nodes(); l++ )  
                   MathOperatorRHS<dim>::RHS[l] += RHS_TEMP(l,1);
                
                RHS_TEMP.Resize( e.Nodes(), e.Nodes() );
              }
         }

// nicePrint( RHS );
      
} // end ComputeContribution


template class NumIntegral_PT_op_P_dV<1U,Element<1U> >;
template class NumIntegral_PT_op_P_dV<2U,Element<2U> >;
template class NumIntegral_PT_op_P_dV<3U,Element<3U> >;

template class NumIntegral_PT_op_P_dV<1U,Face<1U> >;
template class NumIntegral_PT_op_P_dV<2U,Face<2U> >;
template class NumIntegral_PT_op_P_dV<3U,Face<3U> >;

} // csmp
