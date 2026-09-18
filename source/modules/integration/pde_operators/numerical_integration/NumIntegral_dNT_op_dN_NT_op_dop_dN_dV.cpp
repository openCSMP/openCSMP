// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_NT_op_dop_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  grad_prop,
                                                          const char*  eprop,           // e.g., property for multiplication with grad
                                                          const char* emultiplier,      // e.g., property for multiplication with eprop
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
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
    gravity(ACC_GRAVITY),
    xyz(2)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_op_dN_NT_op_dop_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3U ) {
         for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
           it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
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

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_NT_op_dop_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  grad_prop,
                                                          const char*  eprop,           // e.g., property for multiplication with grad
                                                          const char*  emultiplier,      // e.g., property for multiplication with eprop
                                                          const char*  rrho_prop,
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
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
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_op_dN_NT_op_dop_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3U )  {
         for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
               it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
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

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
The gradients of the dependent variable are computed at the integration
points and multiplied with the n- and e-multipliers. These must be 
nodal and element variables, respectively.  

The operand is read
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // 1. read Operand
    // ---------------
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      {
         MathOperatorLHS<dim,CELL>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim,CELL>::MTRL[0].Zero();
      
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
              ScalarVariable  sc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), sc );
              for ( auto i{0U}; i<dim; i++ )
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = sc();
           }
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
              for ( auto i{0U}; i<dim; i++ )
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = vc[i];
           }
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
              for ( auto i{0U}; i<dim; i++ )
                for ( auto j{0U}; j<dim; j++ )
                  MathOperatorLHS<dim,CELL>::MTRL[0](i,j) = ts(i,j);
           }
      }
    else // if a nodal variable is dealt with
      {
         if ( e.FE()->IntegrationPoints() == 0 ) 
           throw csmp::Exception( FATAL_ERROR, "MathOperatorLHS<dim>::GetOperands", 
                                        "The current finite element has no integration points",
                                        "Therefore nodal properties cannot be integrated.");
      
         for ( uint32_t i{0U}; i<e.FE()->IntegrationPoints(); i++ )
           MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim,CELL>::MaterialOperandKey(),
                                                                  i, MathOperatorLHS<dim,CELL>::MTRL[i] );
      }
   
    // 2. if body force operand 'gravity' was specified
    // -------------------------------------------------
    if ( with_gravity ) {
         e.NodePropertyVector( rrho_key, rrho_vec );
         RDENS.resize(e.FE()->IntegrationPoints());
         for ( auto j{0U}; j<e.FE()->IntegrationPoints(); j++ ) RDENS[j] = 0.0;
      }

    // 3. read node variable which must be a scalar
    // --------------------------------------------
    NT3.resize(e.FE()->IntegrationPoints());

    for ( uint32_t i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
         e.N_AtIntegrationPoint( i, IPOL );
         NT3[i].Resize(e.Nodes(),dim);
         for ( uint32_t j{0U}; j<e.Nodes(); j++ )
           {
              // initializing vector of NT3 matrices at integration points
              for ( uint32_t k=0; k<dim; k++ ) NT3[i](j,k) = IPOL[j];
              // gravity operand if so specified
              if ( with_gravity ) RDENS[i] += rrho_vec[j]() * IPOL[j]; 
           }
      }

    // 4. read gradient variable  which must be a scalar
    // -------------------------------------------------
    NGRAD.resize(e.Nodes());
    e.NodePropertyVector( grad_key, sc_prop_vec );
    for ( uint32_t i{0U}; i<e.Nodes(); i++ ) NGRAD[i] = sc_prop_vec[i]();
    
   // 5. get conductivity multiplier for gradient property
   // ----------------------------------------------------
   e.Read( cond_key, econd );
   e.Read( mult_key, emult );

} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    double detJ;

    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();
    DNT.Resize( e.Nodes(), dim );

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ )
      {
         // 1. Compute "diffusion" matrix DNT_K_DN_DV
         // -----------------------------------------
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (DN is already in global coordinates)
         detJ = e.dN_AtIntegrationPoint( DN, i, 1 );

         // transposing DN -> DNT and saving it in DIFF (diffusion matrix)
         DN.Transposed( DNT );

         // multiply  DNT . MTRL
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
           DNT *= MathOperatorLHS<dim,CELL>::MTRL[0];
         else
           DNT *= MathOperatorLHS<dim,CELL>::MTRL[i];

         // multiplying DNT . DN 
         DNT *= DN;
             
         // multiplying with determinant and weights
         DNT *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // accumulating ME Gauss point integral contributions into element 
         // contribution to global conductance matrix
         MathOperatorLHS<dim,CELL>::LHS += DNT;

         // 2. Compute "velocity" 'v' matrix NT3 NTNTNT_V at integration point
         // ------------------------------------------------------------------
         VIP.Resize(dim,dim);
         VIP.Zero();
         econd() *= emult();
         // establish pressure gradients
         for ( uint32_t j{0U}; j<dim; j++ ) {
              //                                x,y,z-component  d/dx,y,z  grad_prop at node
              for ( uint32_t k=0; k<e.Nodes(); k++ ) VIP(j,j) += -DN(j,k) * NGRAD[k] * econd();
              // if gravity is turned on velocities are corrected correspondingly
              if ( with_gravity )  VIP(xyz-1,xyz-1) -= gravity * RDENS[i] * econd(); 
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

         MathOperatorLHS<dim,CELL>::LHS += NT3[i];
      }

// cout <<"\nLHS for element: "<< e.Idx() << endl;
//   nicePrint( MathOperatorLHS<dim>::LHS );

} // end ComputeContribution








template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<dim,CELL>::SpatialDerivative( uint32_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }


template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<1U,Element>;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<2U,Element>;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<3U,Element>;

template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<1U,Face>;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<2U,Face>;
template class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV<3U,Face>;

} // csmp

