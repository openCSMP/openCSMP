#include "NumIntegral_op_PT_P_dV.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<size_t dim,class CELL>
NumIntegral_op_PT_P_dV<dim,CELL>::NumIntegral_op_PT_P_dV( const PropertyDatabase<dim>& pref,
                                                        const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    nodal_degrees_of_freedom(1)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_op_PT_P_dV", oper, test );
    
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
      throw csmp::Exception( ERROR, "NumIntegral_op_PT_P_dV<dim>::(constructor)", 
                      test, "Non-scalar Operands require multiple degrees of freedom per node." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_op_PT_P_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}



/// in the special case, the integral of the testfunction products must be used
template<size_t dim,class CELL>
void NumIntegral_op_PT_P_dV<dim,CELL>::ComputeContribution( CELL& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0.0 );
    double volume = e.Volume();
    size_t  k(0);

    if ( MathOperatorRHS<dim>::LumpedFormulation() )
      {
         // if the property is constant over the element a lumped formulation is suitable, where
         // the element contribution which is uniform over the entire element is just distributed
         // equally over all the nodes.
         if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT )
           {    
              for ( size_t n=0; n<e.Nodes(); n++ ) 
                for ( size_t i=0; i<nodal_degrees_of_freedom; i++ )
                  MathOperatorRHS<dim>::RHS[k++] = (MathOperatorRHS<dim>::MTRL[0](i,i) * volume) / e.Nodes();
           }
         // if a nodal property is accumulated, the integration must be performed 
         // explicitly  
         else 
          {
             // the property value which has already been interpolated to each integration
             // point is assembled
             for ( size_t n=0; n<e.FE()->IntegrationPoints(); n++ )
               { 
                  if ( MathOperatorRHS<dim>::MaterialOperandType() == VECTOR )
                    for ( size_t i=0; i<nodal_degrees_of_freedom; i++ )
                      MathOperatorRHS<dim>::RHS[k++] = (e.WeightAtIntegrationPoint(i) * 
                                                                    MathOperatorRHS<dim>::MTRL[n](i,i) * volume) / e.Nodes();
                  else if ( MathOperatorRHS<dim>::MaterialOperandType() == TENSOR )
                     for ( size_t i=0; i<nodal_degrees_of_freedom; i++ )
                       for ( size_t j=0; j<nodal_degrees_of_freedom; j++ )
                         MathOperatorRHS<dim>::RHS[k++] = (e.WeightAtIntegrationPoint(i) * 
                                                                      MathOperatorRHS<dim>::MTRL[n](i,j) * volume) / e.Nodes();
              }
          }
      }
    // consistent formulation (not tested thus far)  
    else
    cout <<"\nNumIntegral_op_PT_P_dV<dim>::ComputeContribution: Consistent formulation not implemented yet !"<< endl;

}


template class NumIntegral_op_PT_P_dV<1U,Element<1U> >;
template class NumIntegral_op_PT_P_dV<2U,Element<2U> >;
template class NumIntegral_op_PT_P_dV<3U,Element<3U> >;

template class NumIntegral_op_PT_P_dV<1U,Face<1U> >;
template class NumIntegral_op_PT_P_dV<2U,Face<2U> >;
template class NumIntegral_op_PT_P_dV<3U,Face<3U> >;

} // csmp
