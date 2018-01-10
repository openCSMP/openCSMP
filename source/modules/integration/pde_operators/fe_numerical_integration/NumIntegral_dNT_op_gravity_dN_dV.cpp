#include "NumIntegral_dNT_op_gravity_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

/**
 
The Operand which is used here can be both, an element or a nodal variable
which is then interpolated to the integration points to obtain the 
integral properties.  
*/
template<size_t dim,class SIMPLEX>
NumIntegral_dNT_op_gravity_dN_dV<dim,SIMPLEX>::NumIntegral_dNT_op_gravity_dN_dV( const PropertyDatabase<dim>& pref,
                                                                                 const char*           oper,
                                                                                 const char* body_force_variable, 
                                                                                 const char*           basic,
                                                                                 const char*           test,
                                                                                 double64 acc_gravity )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    B(dim,3), BT(3,dim),
    acc_gravity_(acc_gravity),
    g_key_(pref.StorageKey(body_force_variable))
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_op_gravity_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_gravity_dN_dV<dim>::(constructor)", 
                             basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_gravity_dN_dV<dim>::(constructor)", 
                             test, "Operand (test) must be a scalar property placed on the nodes." );

    if ( g_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_gravity_dN_dV<dim>::(constructor)", 
                             body_force_variable, "Operand must be a scalar property." );
}




/** Reads the Operand from either the element or the nodes. Reads the material
property from the element, and its multipliers from the nodes for later
interpolation to the integration points.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_dNT_op_gravity_dN_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
    // if operand property is an element property
   if ( this->MaterialOperandPlacement() == ELEMENT || this->MaterialOperandPlacement() == REGION )
      {
         this->MTRL.resize(1U);
         if ( this->MaterialOperandType() == SCALAR )
           this->MTRL[0].AssignToDiagonal( dim, e.Read(  this->MaterialOperandKey() ) );
         else if ( this->MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( this->MaterialOperandKey(), vc );
              this->MTRL[0].AssignToDiagonal( vc );
           }
         else if ( this->MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read(  this->MaterialOperandKey(), ts );
              this->MTRL[0] = ts;
           }
         else if ( this->MaterialOperandType() == ARRAY )
           {
              ArrayVariable ar( this->MaterialOperandDataDepth() );
              e.Read( this->MaterialOperandKey(), ar );
              this->MTRL[0].AssignToDiagonal( ar );
           }
         else if ( this->MaterialOperandType() == FLAGGEDARRAY )
           {
              FlaggedArrayVariable fr( this->MaterialOperandDataDepth() );
              e.Read( this->MaterialOperandKey(), fr );
              this->MTRL[0].AssignToDiagonal( fr );
           }
      }
    else // if the operand is placed on the constraint-points
      {
         if ( this->MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) {
               this->MTRL.resize( e.IntegrationPoints() );
               const size_t n_integration_points( e.IntegrationPoints() );
               for ( size_t i=0U; i<n_integration_points; i++ )
                  {   
                     if ( this->MaterialOperandType() == SCALAR )
                       this->MTRL[i].AssignToDiagonal( dim, e.Read( i, this->MaterialOperandKey() ) );
                     else if ( this->MaterialOperandType() == VECTOR ) {
                          VectorVariable<dim>  vc;
                          e.Read( i, this->MaterialOperandKey(), vc );
                          this->MTRL[i].AssignToDiagonal( vc );
                       }
                     else if ( this->MaterialOperandType() == TENSOR ) {
                          TensorVariable<dim>  ts;
                          e.Read( i, this->MaterialOperandKey(), ts );
                          this->MTRL[i] = ts;
                       }
                     else if ( this->MaterialOperandType() == ARRAY )
                       {
                          ArrayVariable ar( this->MaterialOperandDataDepth() );
                          e.Read( i, this->MaterialOperandKey(), ar );
                          this->MTRL[0].AssignToDiagonal( ar );
                       }
                     else if ( this->MaterialOperandType() == FLAGGEDARRAY )
                       {
                          FlaggedArrayVariable fr( this->MaterialOperandDataDepth() );
                          e.Read( i, this->MaterialOperandKey(), fr );
                          this->MTRL[0].AssignToDiagonal( fr );
                       }
                 }
           }
           
         // if the operand is placed on the node
         else if ( this->MaterialOperandPlacement() == NODE ) {
             const size_t n_integration_points( e.IntegrationPoints() );
             for ( size_t i=0U; i<n_integration_points; i++ )
               this->PropertyAtIntegrationPoint( e, this->MaterialOperandKey(), i, this->MTRL[i] );
           }
         else
             throw csmp::Exception( FATAL_ERROR,
                                    "MathOperatorLHS<dim>::GetOperands",
                                    "Face based operands cannot be accumulated with this method");
      }
      
   // 2. reading Operand (fluid density or something like that)
   if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) {
        g_property_.resize(1U);
        e.Read( this->MaterialOperandKey(), g_property_[0U] );
     }
   else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) {
        g_property_.resize(e.IntegrationPoints());
        e.IntegrationPointPropertyVector( MathOperatorLHS<dim>::MaterialOperandKey(), g_property_ );
     }
   else {
        g_property_.resize(e.Nodes());
        e.NodePropertyVector( MathOperatorLHS<dim>::MaterialOperandKey(), g_property_ );
     }

} // end GetOperands






/** Laplacian operator of shape function derivatives squared.
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_dNT_op_gravity_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    if ( this->MaterialOperandPlacement() == ELEMENT or
         this->MaterialOperandPlacement() == REGION or
         this->MaterialOperandPlacement() == FACE)
      {
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double64 detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );

             // transposing B -> BT  O.K.
             B.Transposed( BT );

             // multiply  BT . MTRL
             BT *= MathOperatorLHS<dim>::MTRL[0];
          
             // interpolating the g property if necessary
             double64 g_prop_value(0.);
             if ( g_key_.place == NODE ) {
                  e.N_AtIntegrationPoint( i, this->IPOL );
                  for ( size_t j=0; j<e.Nodes(); j++ ) g_prop_value += IPOL[j] * g_property[j]();
               }
             else if ( g_key_.place == ELEMENT )
               g_prop_value = g_property[0]();
          
             // add the gravity component
             for ( size_t j=0; j<e.Nodes(); j++ )
               B(j,0) += B

             // multiplying BT . B 
             BT *= B;
                 
             // multiplying with determinant and weights
             BT *= e.WeightAtIntegrationPoint(i) * detJ;

             // accumulating ME Gauss point integral contributions into element 
             // contribution to global conductance matrix
             MathOperatorLHS<dim>::LHS += BT;
          }
      }
    else { // NODE or ELEMENT_INTEGRATION_POINT
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) {
             double64 detJ = e.dN_AtIntegrationPoint( B, i, SCALAR );
             B.Transposed( BT );
             BT *= MathOperatorLHS<dim>::MTRL[i];
             BT *= B;
             BT *= e.WeightAtIntegrationPoint(i) * detJ; 
             MathOperatorLHS<dim>::LHS += BT;
          }
     }

} // end ComputeContribution


//cout <<"\nNumIntegral_dNT_op_gravity_dN_dV: on Element "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( size_t i=0U; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nNumIntegral_dNT_op_gravity_dN_dV: zero element in diagonal of element matrix."; 


template class NumIntegral_dNT_op_gravity_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_op_gravity_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_op_gravity_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_op_gravity_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_op_gravity_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_op_gravity_dN_dV<3U,Face<3U> >;

//template class NumIntegral_dNT_op_gravity_dN_dV<1U,InterFace<1U> >;
//template class NumIntegral_dNT_op_gravity_dN_dV<2U,InterFace<2U> >;
//template class NumIntegral_dNT_op_gravity_dN_dV<3U,InterFace<3U> >;

} // csmp











