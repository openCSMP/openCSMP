#include "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_dV_NT_v_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  emultiplier, 
                                                          const char*  nmultiplier, 
                                                          const char*  grad_prop,
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    emulti_key(pref.StorageKey(emultiplier)),
    nmulti_key(pref.StorageKey(nmultiplier)),
    grad_key(pref.StorageKey(grad_prop)),
    DN(dim,3), 
    DNT(3,dim),
    IPOL(3),
    VIP(dim,dim),
    NT3(3,DenseMatrix<DM_MIN>(3,dim)),
    NMULT(3),
    NGRAD(3),
    with_gravity(false),
    gravity(ACC_GRAVITY),
    xyz(2)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_op_dN_dV_NT_v_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if constexpr ( dim == 3U )
      {
         for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
               it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }
      
    if ( emulti_key.place != ELEMENT and emulti_key.place != REGION )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      emultiplier, "must be an element- or group-based variable." );

    if ( nmulti_key.place != NODE || nmulti_key.type != SCALAR )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      nmultiplier, "must be an node-based scalar variable." );

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      grad_prop, "Operand to calculate 'v' from must be a scalar property placed on the nodes." );


    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_dV_NT_v_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  emultiplier, 
                                                          const char*  nmultiplier, 
                                                          const char*  grad_prop,
                                                          const char*  rrho_prop,
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    emulti_key(pref.StorageKey(emultiplier)),
    nmulti_key(pref.StorageKey(nmultiplier)),
    grad_key(pref.StorageKey(grad_prop)),
    rrho_key(pref.StorageKey(rrho_prop)),
    DN(dim,3), 
    DNT(3,dim),
    IPOL(3),
    VIP(dim,dim),
    NT3(3,DenseMatrix<DM_MIN>(3,dim)),
    NMULT(3),
    NGRAD(3),
    with_gravity(true),
    gravity(ACC_GRAVITY),
    xyz(2)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_op_dN_dV_NT_v_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if constexpr ( dim == 3U )
      {
         for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
               it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    if ( emulti_key.place != ELEMENT )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      emultiplier, "must be an element-based variable." );

    if ( nmulti_key.place != NODE || nmulti_key.type != SCALAR )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                      nmultiplier, "must be an node-based scalar variable." );

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                             grad_prop, "Operand to calculate 'v' from must be a scalar property placed on the nodes." );

    if ( rrho_key.place != NODE || rrho_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                             rrho_prop, "Operand to calculate 'rho g' term from must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                             basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim>::(constructor)", 
                             test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
The gradients of the dependent variable are computed at the integration
points and multiplied with the n- and e-multipliers. These must be 
nodal and element variables, respectively.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
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
              const double sc = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = sc;
           }
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = vc[i];
           }
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                for ( uint32_t j{0U}; j<dim; j++ ) 
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
   
    // 2. read element multiplier variable
    // -----------------------------------
    ReadElementMultiplier( e, EMULT );
    
    // 3. if body force operand 'gravity' was specified
    // -------------------------------------------------
    uint32_t  j;
    if ( with_gravity ) {
         e.NodePropertyVector( rrho_key, rrho_vec );
         RDENS.resize(e.FE()->IntegrationPoints());
         fill( RDENS.begin(), RDENS.end(), 0. );
      }

    // 4. read node multiplier variable which must be a scalar
    // -------------------------------------------------------
    e.NodePropertyVector( nmulti_key, sc_prop_vec );
    NMULT.resize(e.FE()->IntegrationPoints());
    NT3.resize(e.FE()->IntegrationPoints());

    for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
         e.N_AtIntegrationPoint( i, IPOL );
         NT3[i].Resize(e.Nodes(),dim);
         for ( NMULT[i]=0., j=0u; j<e.Nodes(); j++ )
           {
              // initializing vector of NT3 matrices at integration points
              for ( auto k=0u; k<dim; k++ ) NT3[i](j,k) = IPOL[j];
              // interpolating NMULT to the integration points
              NMULT[i] += sc_prop_vec[j]() * IPOL[j];
              // gravity operand if so specified
              if ( with_gravity ) RDENS[i] += rrho_vec[j]() * IPOL[j]; 
           }
      }

    // 5. read gradient variable  which must be a scalar
    // -------------------------------------------------
    NGRAD.resize(e.Nodes());
    e.NodePropertyVector( grad_key, sc_prop_vec );
    for ( auto i{0U}; i<e.Nodes(); i++ ) NGRAD[i] = sc_prop_vec[i]();
    
} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

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
         double detJ = e.dN_AtIntegrationPoint( DN, i, 1 );

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
         VIP.Zero();
         for ( auto j{0U}; j<dim; j++ ) {
              //                                x,y,z-component  d/dx,y,z  grad_prop at node
              for ( auto k{0U}; k<e.Nodes(); k++ ) VIP(j,j) += -DN(j,k) * NGRAD[k];
              // adjusting vertical flow component if gravity is acting
              if ( with_gravity && j == xyz-1 )       VIP(j,j) -=  gravity * RDENS[i]; 

              // multiply with interpolated scalar nodal multiplier
              VIP(j,j) *= NMULT[i];
           }
         // multiply with element multiplier
         VIP *= EMULT;
                  
         // multiply with transposed basis function matrix at integration point
         NT3[i] *= VIP;
        
         // 3. Multiply with DN, integrate with weight and detJ, and add to element contribution
         // ------------------------------------------------------------------------------------
         NT3[i] *= DN;
         NT3[i] *= e.WeightAtIntegrationPoint(i) * detJ;

         MathOperatorLHS<dim,CELL>::LHS += NT3[i];
      }

} // end ComputeContribution








template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim,CELL>::ReadElementMultiplier( const CELL<dim>& e,
                                                                           DenseMatrix<DM_MIN>& MULT )
 {
    MULT.Resize(dim,dim);
 
    if ( emulti_key.type == SCALAR ) {
         ScalarVariable  sc;
         e.Read( emulti_key, sc );
         MULT.Zero();
         for ( auto i{0U}; i<dim; i++ ) MULT(i,i) = sc();
      }
    else if ( emulti_key.type == VECTOR ) {
         VectorVariable<dim>  vc;
         e.Read( emulti_key, vc );
         MULT.Zero();
         for ( auto i{0U}; i<dim; i++ ) MULT(i,i) = vc[i];
      }
    else // TENSOR
      {
         TensorVariable<dim>  ts;
         e.Read( emulti_key, ts );
         for ( auto i{0U}; i<dim; i++ )
           for ( auto j{0U}; j<dim; j++ ) MULT(i,j) = ts(i,j);
      }
      
 } // end ReadElementMultiplier




template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<dim,CELL>::SpatialDerivative( uint32_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }


template class NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<1U,Element>;
template class NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<2U,Element>;
template class NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<3U,Element>;

template class NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<1U,Face>;
template class NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<2U,Face>;
template class NumIntegral_dNT_op_dN_dV_NT_v_dN_dV<3U,Face>;

} // csmp





