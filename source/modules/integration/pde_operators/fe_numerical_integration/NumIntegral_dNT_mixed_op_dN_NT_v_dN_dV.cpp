#include "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<size_t dim,class CELL>
NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>::NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  emultiplier, 
                                                          const char*  nmultiplier, 
                                                          const char*  grad_prop,
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
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
    gravity(9.80665),
    xyz(2)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3 ) {
         typename vector<DenseMatrix<DM_MIN> >::iterator  it;
         for ( it=MathOperatorLHS<dim>::MTRL.begin(); 
               it!=MathOperatorLHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }
      
    if ( emulti_key.place != ELEMENT and emulti_key.place != REGION )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      emultiplier, "must be an element- or region based variable." );

    if ( nmulti_key.place != NODE || nmulti_key.type != SCALAR )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      nmultiplier, "must be an node-based scalar variable." );

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      grad_prop, "Operand to calculate 'v' from must be a scalar property placed on the nodes." );


    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





template<size_t dim,class CELL>
NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>::NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  emultiplier, 
                                                          const char*  nmultiplier, 
                                                          const char*  grad_prop,
                                                          const char*  rrho_prop,
                                                          const char*  oper, 
                                                          const char*  basic, 
                                                          const char*  test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
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
    gravity(9.80665),
    xyz(2)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3 ) 
      {
         typename vector<DenseMatrix<DM_MIN> >::iterator  it;
         for ( it=MathOperatorLHS<dim>::MTRL.begin(); 
               it!=MathOperatorLHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    if ( emulti_key.place != ELEMENT )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      emultiplier, "must be an element-based variable." );

    if ( nmulti_key.place != NODE || nmulti_key.type != SCALAR )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      nmultiplier, "must be an node-based scalar variable." );

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      grad_prop, "Operand to calculate 'v' from must be a scalar property placed on the nodes." );

    if ( rrho_key.place != NODE || rrho_key.type != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      rrho_prop, "Operand to calculate 'rho g' term from must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
The gradients of the dependent variable are computed at the integration
points and multiplied with the n- and e-multipliers. These must be 
nodal and element variables, respectively.  

The operand is read
 */
template<size_t dim,class CELL>
void NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>::GetOperands( const CELL& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // 1. read Operand
    // ---------------
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
      {
         MathOperatorLHS<dim>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim>::MTRL[0].Zero();
      
         if (MathOperatorLHS<dim>:: MaterialOperandType() == SCALAR ) {
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
   
    // 2. read element multiplier variable
    // -----------------------------------
    ReadElementMultiplier( e, EMULT );
    
    // 3. if body force operand 'gravity' was specified
    // -------------------------------------------------
    size_t  j;
    if ( with_gravity ) {
         e.NodePropertyVector( rrho_key, rrho_vec );
         RDENS.resize(e.FE()->IntegrationPoints());
         for ( j=0; j<e.FE()->IntegrationPoints(); j++ ) RDENS[j] = 0.0;
      }

    // 4. read node multiplier variable which must be a scalar
    // -------------------------------------------------------
    e.NodePropertyVector( nmulti_key, sc_prop_vec );
    NMULT.resize(e.FE()->IntegrationPoints());
    NT3.resize(e.FE()->IntegrationPoints());

    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) {
         e.N_AtIntegrationPoint( i, IPOL );
         NT3[i].Resize(e.Nodes(),dim);
         for ( NMULT[i]=0.0, j=0; j<e.Nodes(); j++ ) 
           {
              // interpolating NMULT to the integration points
              NMULT[i] += sc_prop_vec[j]() * IPOL[j];
              // initializing vector of NT3 matrices at integration points
              for ( size_t k=0; k<dim; k++ ) NT3[i](j,k) = IPOL[j];
              // gravity operand if so specified
              if ( with_gravity ) RDENS[i] += rrho_vec[j]() * IPOL[j]; 
           }
      }

    // 5. read gradient variable  which must be a scalar
    // -------------------------------------------------
    NGRAD.resize(e.Nodes());
    e.NodePropertyVector( grad_key, sc_prop_vec );
    for ( size_t i=0; i<e.Nodes(); i++ ) NGRAD[i] = sc_prop_vec[i]();
    
} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class CELL>
void NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

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
         double detJ = e.dN_AtIntegrationPoint( DN, i, 1 );

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
         VIP.Zero();
         for ( size_t j=0; j<dim; j++ ) {
              //                                x,y,z-component  d/dx,y,z  grad_prop at node
              for ( size_t k=0; k<e.Nodes(); k++ ) VIP(j,j) += -DN(j,k) * NGRAD[k];
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

         MathOperatorLHS<dim>::LHS += NT3[i];
      }

} // end ComputeContribution





template<size_t dim,class CELL>
void NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>::ReadElementMultiplier( const CELL& e,
                                                                              DenseMatrix<DM_MIN>& MULT )
 {
    MULT.Resize(dim,dim);
 
    if ( emulti_key.type == SCALAR ) {
         ScalarVariable  sc;
         e.Read( emulti_key, sc );
         MULT.Zero();
         for ( size_t i=0; i<dim; i++ ) MULT(i,i) = sc();
      }
    else if ( emulti_key.type == VECTOR ) {
         VectorVariable<dim>  vc;
         e.Read( emulti_key, vc );
         MULT.Zero();
         for ( size_t i=0; i<dim; i++ ) MULT(i,i) = vc[i];
      }
    else // TENSOR
      {
         TensorVariable<dim>  ts;
         e.Read( emulti_key, ts );
         for ( size_t i=0; i<dim; i++ )
           for ( size_t j=0; j<dim; j++ ) MULT(i,j) = ts(i,j);
      }
      
 } // end ReadElementMultiplier




template<size_t dim,class CELL>
void NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>::SpatialDerivative( size_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }

template class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<3U,Face<3U> >;

}




