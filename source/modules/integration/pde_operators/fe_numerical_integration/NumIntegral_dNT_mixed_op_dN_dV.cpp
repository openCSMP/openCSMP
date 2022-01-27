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
template<size_t dim,class CELL>
NumIntegral_dNT_mixed_op_dN_dV<dim,CELL>::NumIntegral_dNT_mixed_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                            const char*           oper, 
                                                            const char*           nodal_oper_multiplier,
                                                            const char*           basic, 
                                                            const char*           test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    B(dim,3), BT(3,dim),
    ip_nmult(3),
    nkey(pref.StorageKey(nodal_oper_multiplier))
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_mixed_op_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3 ) {
         typename vector<DenseMatrix<DM_MIN> >::iterator  it;
         for ( it=MathOperatorLHS<dim>::MTRL.begin(); 
               it!=MathOperatorLHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }
      
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      oper, "Operand (basic) must be placed on the element or group." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );

    if ( nkey.place != NODE || nkey.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_dV<dim>::(constructor)", 
                      nodal_oper_multiplier, "Operand multiplier must be a scalar property placed on the nodes." );
}



template<size_t dim,class CELL>
void NumIntegral_dNT_mixed_op_dN_dV<dim,CELL>::GetOperands( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // only if the property is an element property  something is done here
     MathOperatorLHS<dim>::MTRL[0].Resize(dim,dim);
     MathOperatorLHS<dim>::MTRL[0].Zero();
  
     if ( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR ) {
          ScalarVariable  sc;
          e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), sc );
          MathOperatorLHS<dim>::MTRL[0](0,0) =  sc();
          if ( dim != 1U ) MathOperatorLHS<dim>::MTRL[0](1,1) = sc();
          if ( dim == 3U ) MathOperatorLHS<dim>::MTRL[0](2,2) = sc();
       }
     if ( MathOperatorLHS<dim>::MaterialOperandType() == VECTOR ) {
          VectorVariable<dim>  vc;
          e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), vc );
          MathOperatorLHS<dim>::MTRL[0](0,0) = vc[0];
          if ( dim != 1U ) MathOperatorLHS<dim>::MTRL[0](1,1) = vc[1];
          if ( dim == 3U ) MathOperatorLHS<dim>::MTRL[0](2,2) = vc[2];
       }
     if ( MathOperatorLHS<dim>::MaterialOperandType() == TENSOR ) {
          TensorVariable<dim>  ts;
          e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), ts );
          MathOperatorLHS<dim>::MTRL[0] = ts;
       }

    // 2. read the nodal operand multiplier and interpolate it to integration points
    // ------------------------------------
    ip_nmult.resize( e.FE()->IntegrationPoints() );
    
    // getting property values at the integration points assuming we have a scalar (which was verified earlier)
    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      ip_nmult[i] = e.PropertyValueAtIntegrationPoint( nkey, i );

 } // end GetOperands








/**
 
@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class CELL>
void NumIntegral_dNT_mixed_op_dN_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    double detJ;

    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         detJ = e.dN_AtIntegrationPoint( B, i, 1 );

         // transposing B -> BT  O.K.
         B.Transposed( BT );

         // multiply  BT . MTRL x nodal multiplier
         NVAL.Resize(dim,dim);
         NVAL  = MathOperatorLHS<dim>::MTRL[0];
         NVAL *= ip_nmult[i];
         BT   *= NVAL; 

         // multiplying BT . B 
         BT *= B;
             
         // multiplying with determinant and weights
         BT *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // accumulating ME Gauss point integral contributions into element 
         // contribution to global conductance matrix
         MathOperatorLHS<dim>::LHS += BT;
      }

//   nicePrint( MathOperatorLHS<dim>::LHS );

} // end ComputeContribution

template class NumIntegral_dNT_mixed_op_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_mixed_op_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_mixed_op_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_mixed_op_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_mixed_op_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_mixed_op_dN_dV<3U,Face<3U> >;

} // csmp











