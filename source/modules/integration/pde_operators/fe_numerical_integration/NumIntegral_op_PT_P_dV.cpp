#include "NumIntegral_op_PT_P_dV.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_op_PT_P_dV<dim,CELL>::NumIntegral_op_PT_P_dV( const PropertyDatabase<dim>& pref,
                                                          const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    nodal_degrees_of_freedom(1)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_op_PT_P_dV", oper, test );
    
    // resize material property matrix 
    if constexpr ( dim == 3U )
      for ( auto it=MathOperatorRHS<dim,CELL>::MTRL.begin();
            it!=MathOperatorRHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);

    // anisotropy can only be considered if there are multiple degrees of freedom per node
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR ) nodal_degrees_of_freedom = dim;
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR && 
         MathOperatorRHS<dim,CELL>::TestOperandType() == SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_op_PT_P_dV<dim>::(constructor)", 
                      test, "Non-scalar Operands require multiple degrees of freedom per node." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_op_PT_P_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}



/// in the special case, the integral of the testfunction products must be used
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_op_PT_P_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );
    double volume = e.Volume();
    size_t  k(0);

    if ( MathOperatorRHS<dim,CELL>::LumpedFormulation() )
      {
         // if the property is constant over the element a lumped formulation is suitable, where
         // the element contribution which is uniform over the entire element is just distributed
         // equally over all the nodes.
         if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
           {    
              for ( auto n=0; n<e.Nodes(); n++ ) 
                for ( auto i{0U}; i<nodal_degrees_of_freedom; i++ )
                  MathOperatorRHS<dim,CELL>::RHS[k++] = (MathOperatorRHS<dim,CELL>::MTRL[0](i,i) * volume) / e.Nodes();
           }
         // if a nodal property is accumulated, the integration must be performed 
         // explicitly  
         else 
          {
             // the property value which has already been interpolated to each integration
             // point is assembled
             for ( auto n=0; n<e.FE()->IntegrationPoints(); n++ )
               { 
                  if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() == VECTOR )
                    for ( auto i{0U}; i<nodal_degrees_of_freedom; i++ )
                      MathOperatorRHS<dim,CELL>::RHS[k++] = (e.WeightAtIntegrationPoint(i) * 
                                                                    MathOperatorRHS<dim,CELL>::MTRL[n](i,i) * volume) / e.Nodes();
                  else if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() == TENSOR )
                     for ( auto i{0U}; i<nodal_degrees_of_freedom; i++ )
                       for ( auto j{0U}; j<nodal_degrees_of_freedom; j++ )
                         MathOperatorRHS<dim,CELL>::RHS[k++] = (e.WeightAtIntegrationPoint(i) * 
                                                                      MathOperatorRHS<dim,CELL>::MTRL[n](i,j) * volume) / e.Nodes();
              }
          }
      }
    // consistent formulation (not tested thus far)  
    else
    cout <<"\nNumIntegral_op_PT_P_dV<dim>::ComputeContribution: Consistent formulation not implemented yet !"<< endl;

}


template class NumIntegral_op_PT_P_dV<1U,Element>;
template class NumIntegral_op_PT_P_dV<2U,Element>;
template class NumIntegral_op_PT_P_dV<3U,Element>;

template class NumIntegral_op_PT_P_dV<1U,Face>;
template class NumIntegral_op_PT_P_dV<2U,Face>;
template class NumIntegral_op_PT_P_dV<3U,Face>;

} // csmp
