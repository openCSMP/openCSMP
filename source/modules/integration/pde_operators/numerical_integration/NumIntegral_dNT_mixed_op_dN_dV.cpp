#include "NumIntegral_dNT_mixed_op_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

/**
 
The Operand which is used here can be both, an element or a nodal variable
which is then interpolated to the integration points to obtain the 
integral properties.  
 */
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_mixed_op_dN_dV<dim,CELL>::NumIntegral_dNT_mixed_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                                          const char*           oper,
                                                                          const char*           nodal_oper_multiplier,
                                                                          const char*           basic,
                                                                          const char*           test ) 
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    B(dim,3), BT(3,dim),
    ip_nmult(3),
    nkey(pref.StorageKey(nodal_oper_multiplier))
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_mixed_op_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if constexpr ( dim == 3U )
     for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
           it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
      
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      oper, "Operand (basic) must be placed on the element or group." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );

    if ( nkey.place != NODE || nkey.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                             nodal_oper_multiplier, "Operand multiplier must be a scalar property placed on the nodes." );
}





template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_mixed_op_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // only if the property is an element property  something is done here
     MathOperatorLHS<dim,CELL>::MTRL[0].Resize(dim,dim);
     MathOperatorLHS<dim,CELL>::MTRL[0].Zero();
  
     if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
          double sc = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
          MathOperatorLHS<dim,CELL>::MTRL[0](0,0) =  sc;
          if constexpr ( dim != 1U ) MathOperatorLHS<dim,CELL>::MTRL[0](1,1) = sc;
          if constexpr ( dim == 3U ) MathOperatorLHS<dim,CELL>::MTRL[0](2,2) = sc;
       }
     if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
          VectorVariable<dim>  vc;
          e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
          MathOperatorLHS<dim,CELL>::MTRL[0](0,0) = vc[0];
          if constexpr ( dim != 1U ) MathOperatorLHS<dim,CELL>::MTRL[0](1,1) = vc[1];
          if constexpr ( dim == 3U ) MathOperatorLHS<dim,CELL>::MTRL[0](2,2) = vc[2];
       }
     if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
          TensorVariable<dim>  ts;
          e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
          MathOperatorLHS<dim,CELL>::MTRL[0] = ts;
       }

    // 2. read the nodal operand multiplier and interpolate it to integration points
    // ------------------------------------
    ip_nmult.resize( e.FE()->IntegrationPoints() );
    
    // getting property values at the integration points assuming we have a scalar (which was verified earlier)
    for ( uint32_t i{0U}; i<e.FE()->IntegrationPoints(); i++ )
      ip_nmult[i] = e.PropertyValueAtIntegrationPoint( nkey, i );

 } // end GetOperands








/**
 
@section application Application

In linear elasticity computations.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_mixed_op_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    double detJ;

    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    for ( uint32_t i{0U}; i<e.IntegrationPoints(); i++ )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         detJ = e.dN_AtIntegrationPoint( B, i, 1 );

         // transposing B -> BT  O.K.
         B.Transposed( BT );

         // multiply  BT . MTRL x nodal multiplier
         BT *= (MathOperatorLHS<dim,CELL>::MTRL[i] *= ip_nmult[i]);

         // multiplying BT . B 
         BT *= B;
             
         // multiplying with determinant and weights
         BT *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // accumulating ME Gauss point integral contributions into element 
         // contribution to global conductance matrix
         MathOperatorLHS<dim,CELL>::LHS += BT;
      }

//   nicePrint( MathOperatorLHS<dim>::LHS );

} // end ComputeContribution

template class NumIntegral_dNT_mixed_op_dN_dV<1U,Element>;
template class NumIntegral_dNT_mixed_op_dN_dV<2U,Element>;
template class NumIntegral_dNT_mixed_op_dN_dV<3U,Element>;

template class NumIntegral_dNT_mixed_op_dN_dV<1U,Face>;
template class NumIntegral_dNT_mixed_op_dN_dV<2U,Face>;
template class NumIntegral_dNT_mixed_op_dN_dV<3U,Face>;

} // csmp











