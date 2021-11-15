#include "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<size_t dim,class CELL>
NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_NT_op_dop_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  grad_prop,
                                                          const char*  eprop,           // e.g., property for multiplication with grad
                                                          const char* emultiplier,      // e.g., property for multiplication with eprop
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    grad_key(pref.StorageKey(grad_prop)),
    cond_key(pref.StorageKey(eprop)),
    mult_key(pref.StorageKey(emultiplier)),
    DN(dim,3), 
    DNT(3,dim),
    IPOL(3),
    VIP(dim,dim),
    NT3(3,DenseMatrix<DM_MIN>(3,dim)),
    NGRAD(3),
    with_gravity(false),
    gravity(9.80665),
    xyz(2)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_op_dN_NT_op_dop_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3U ) {
         for ( typename vector<DenseMatrix<DM_MIN> >::iterator
               it=MathOperatorLHS<dim>::MTRL.begin(); 
           it!=MathOperatorLHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }
      
    if ( grad_key.place != NODE || grad_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      grad_prop, "Operand to calculate 'v' from must be a scalar property placed on the nodes." );

    if ( (cond_key.place != ELEMENT && cond_key.place != REGION) || cond_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      eprop, "Operand to calculate 'op dop' product must be a scalar property placed on the element or group." );

    if ( mult_key.place != ELEMENT || mult_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      emultiplier, "must be a scalar property placed on th element." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





template<size_t dim,class CELL>
NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_NT_op_dop_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  grad_prop,
                                                          const char*  eprop,           // e.g., property for multiplication with grad
                                                          const char*  emultiplier,      // e.g., property for multiplication with eprop
                                                          const char*  rrho_prop,
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    grad_key(pref.StorageKey(grad_prop)),
    cond_key(pref.StorageKey(eprop)),
    mult_key(pref.StorageKey(emultiplier)),
    rrho_key(pref.StorageKey(rrho_prop)),
    DN(dim,3), 
    DNT(3,dim),
    IPOL(3),
    VIP(dim,dim),
    NT3(3,DenseMatrix<DM_MIN>(3,dim)),
    NGRAD(3),
    with_gravity(true),
    gravity(9.80665),
    xyz(2)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_op_dN_NT_op_dop_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3U )  {
         for ( typename vector<DenseMatrix<DM_MIN> >::iterator
               it=MathOperatorLHS<dim>::MTRL.begin(); 
               it!=MathOperatorLHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      grad_prop, "Operand to calculate 'v' from must be a scalar property placed on the nodes." );

    if ( (cond_key.place != ELEMENT && cond_key.place != REGION) || cond_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      eprop, "Operand to calculate 'op dop' product must be a scalar property placed on the element or group." );

    if ( mult_key.place != ELEMENT || mult_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      emultiplier, "must be a scalar property placed on th element." );

    if ( rrho_key.place != NODE || rrho_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      rrho_prop, "Operand to calculate 'rho g' term from must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
The gradients of the dependent variable are computed at the integration
points and multiplied with the n- and e-multipliers. These must be 
nodal and element variables, respectively.  

The operand is read
*/
template<size_t dim,class CELL>
void NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::GetOperands( CELL& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // 1. read Operand
    // ---------------
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
      {
         MathOperatorLHS<dim>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim>::MTRL[0].Zero();
      
         if ( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR ) {
              ScalarVariable  sc;
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), sc );
              for ( size_t i=0; i<dim; i++ )
                MathOperatorLHS<dim>::MTRL[0](i,i) = sc();
           }
         if ( MathOperatorLHS<dim>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), vc );
              for ( size_t i=0; i<dim; i++ ) 
                MathOperatorLHS<dim>::MTRL[0](i,i) = vc[i];
           }
         if ( MathOperatorLHS<dim>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), ts );
              for ( size_t i=0; i<dim; i++ ) 
                for ( size_t j=0; j<dim; j++ ) 
                  MathOperatorLHS<dim>::MTRL[0](i,j) = ts(i,j);
           }
      }
    else // if a nodal variable is dealt with
      {
         if ( e.FE()->IntegrationPoints() == 0 ) 
           throw csmp::Exception( FATAL_ERROR, "MathOperatorLHS<dim>::GetOperands", 
                                        "The current finite element has no integration points",
                                        "Therefore nodal properties cannot be integrated.");
      
         for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
           MathOperatorLHS<dim>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim>::MaterialOperandKey(), 
                                                                i, MathOperatorLHS<dim>::MTRL[i] );
      }
   
    // 2. if body force operand 'gravity' was specified
    // -------------------------------------------------
    size_t  j;
    if ( with_gravity ) {
         e.NodePropertyVector( rrho_key, rrho_vec );
         RDENS.resize(e.FE()->IntegrationPoints());
         for ( j=0; j<e.FE()->IntegrationPoints(); j++ ) RDENS[j] = 0.0;
      }

    // 3. read node variable which must be a scalar
    // --------------------------------------------
    NT3.resize(e.FE()->IntegrationPoints());

    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) {
         e.N_AtIntegrationPoint( i, IPOL );
         NT3[i].Resize(e.Nodes(),dim);
         for ( j=0; j<e.Nodes(); j++ ) 
           {
              // initializing vector of NT3 matrices at integration points
              for ( size_t k=0; k<dim; k++ ) NT3[i](j,k) = IPOL[j];
              // gravity operand if so specified
              if ( with_gravity ) RDENS[i] += rrho_vec[j]() * IPOL[j]; 
           }
      }

    // 4. read gradient variable  which must be a scalar
    // -------------------------------------------------
    NGRAD.resize(e.Nodes());
    e.NodePropertyVector( grad_key, sc_prop_vec );
    for ( size_t i=0; i<e.Nodes(); i++ ) NGRAD[i] = sc_prop_vec[i]();
    
   // 5. get conductivity multiplier for gradient property
   // ----------------------------------------------------
   e.Read( cond_key, econd );
   e.Read( mult_key, emult );

} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class CELL>
void NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::ComputeContribution( CELL& e )
 {
    double detJ;

    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();
    DNT.Resize( e.Nodes(), dim );

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      {
         // 1. Compute "diffusion" matrix DNT_K_DN_DV
         // -----------------------------------------
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (DN is already in global coordinates)
         detJ = e.dN_AtIntegrationPoint( DN, i, 1 );

         // transposing DN -> DNT and saving it in DIFF (diffusion matrix)
         DN.Transposed( DNT );

         // multiply  DNT . MTRL
         if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
           DNT *= MathOperatorLHS<dim>::MTRL[0];
         else                             
           DNT *= MathOperatorLHS<dim>::MTRL[i];

         // multiplying DNT . DN 
         DNT *= DN;
             
         // multiplying with determinant and weights
         DNT *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // accumulating ME Gauss point integral contributions into element 
         // contribution to global conductance matrix
         MathOperatorLHS<dim>::LHS += DNT;

         // 2. Compute "velocity" 'v' matrix NT3 NTNTNT_V at integration point
         // ------------------------------------------------------------------
         VIP.Resize(dim,dim);
         VIP.Zero();
         econd() *= emult();
         // establish pressure gradients
         for ( size_t j=0; j<dim; j++ ) {
              //                                x,y,z-component  d/dx,y,z  grad_prop at node
              for ( size_t k=0; k<e.Nodes(); k++ ) VIP(j,j)         += -DN(j,k) * NGRAD[k] * econd();
              // if gravity is turned on velocities are corrected correspondingly
              if ( with_gravity )                     VIP(xyz-1,xyz-1) -= gravity * RDENS[i] * econd(); 
           }

// cout <<"\nEMULT at ip "<< i << endl;
//   nicePrint( EMULT );

// cout <<"\nVIP at ip "<< i << endl;
//   nicePrint( VIP );
         
         // multiply with transposed basis function matrix at integration point
         NT3[i] *= VIP;
        
// cout <<"\nNT3 at same ip." << endl;
//   nicePrint( NT3[i] );

         // 3. Multiply with DN, integrate with weight and detJ, and add to element contribution
         // ------------------------------------------------------------------------------------
         NT3[i] *= DN;
         NT3[i] *= e.WeightAtIntegrationPoint(i) * detJ;

         MathOperatorLHS<dim>::LHS += NT3[i];
      }

// cout <<"\nLHS for element: "<< e.Idx() << endl;
//   nicePrint( MathOperatorLHS<dim>::LHS );

} // end ComputeContribution








template<size_t dim,class CELL>
void NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::SpatialDerivative( size_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }


template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<3U,Face<3U> >;

} // csmp

