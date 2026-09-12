#include "NumIntegral_DNT_op_DN_NT_v_DN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>& pref,
                                                                         const char* diffusion_oper,   // element prop, for instance thermal conductivity
                                                                         const char* advection_oper,   // element prop, for instance heat transport velocity
                                                                         const char* basic,            // e.g., fluid pressure
                                                                         const char* test ) 
  : MathOperatorLHS<dim,CELL>(pref,diffusion_oper,basic,test),
    adv_key(pref.StorageKey(advection_oper)),
    DN(dim,3), 
    DNT(3,dim),
    IPOL(3),
    VIP(dim,dim),
    NT3(3,dim)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_dNT_op_dN_NT_v_dN_dV", diffusion_oper, basic, test );
    
    for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
          it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ )
      (*it).Resize(dim,dim);
      
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT and
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      diffusion_oper, "must be an element-based variable." );

    if ( (adv_key.place != ELEMENT or adv_key.place != ELEMENT_INTEGRATION_POINT) and adv_key.type != VECTOR )
      throw csmp::Exception( ERROR,  "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      advection_oper, "must be an element-based vector variable." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_NT_v_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}









/**
 
The diffusion (op) and advection (adv) coefficients are read from the storage in the model. 

*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // 1. read diffusion Operand
    // -------------------------
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      {
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
              double sc = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = sc;
           }
         else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = vc[i];
           }
         else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
              for ( uint32_t i{0U}; i<dim; i++ ) 
                for ( uint32_t j{0U}; j<dim; j++ ) 
                  MathOperatorLHS<dim,CELL>::MTRL[0](i,j) = ts(i,j);
           }
      }
    else { // or an integration point variable
         for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ )
           MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim,CELL>::MaterialOperandKey(),
                                                                  i, MathOperatorLHS<dim,CELL>::MTRL[i] );
      }
   
    // 2. read element advection variable
    // -----------------------------------
    e.Read( adv_key, velo_ );
    
    
} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();
    
    NT3.Resize(e.Nodes(),dim);

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ )
      {
        // 1. Compute "diffusion" matrix DNT_K_DN_DV
        // -----------------------------------------
        // getting global intpol. function derivative matrix and determinant of
        // byproduct Jacobian matrix (DN is already in global coordinates)
        double detJ = e.dN_AtIntegrationPoint( DN, i, SCALAR );

        // transposing DN -> DNT and saving it in DIFF (diffusion matrix)
        DN.Transposed( DNT );

        // multiply  DNT . MTRL
        if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
          DNT *= MathOperatorLHS<dim,CELL>::MTRL[0];
        else // integration point variable
          DNT *= MathOperatorLHS<dim,CELL>::MTRL[i];

        // multiplying DNT . DN
        DNT *= DN;

        // multiplying with determinant and weights
        DNT *= e.WeightAtIntegrationPoint(i) * detJ;

        // accumulating ME Gauss point integral contributions into element
        // contribution to global conductance matrix
        MathOperatorLHS<dim,CELL>::LHS += DNT;

        // 2. Compute "velocity" 'VIP' matrix NT3 NTNTNT_V at integration point
        // --------------------------------------------------------------------
        VIP.Zero();
        for ( uint32_t j{0U}; j<dim; j++ ) VIP(j,j) = velo_[j];

        e.N_AtIntegrationPoint( i, IPOL );

        NT3.Resize(e.Nodes(),dim);
        for ( uint32_t j{0U}; j<e.Nodes(); j++ ) {
             NT3(j,0) = IPOL[j];
             if ( dim == 2U ) NT3(j,1) = IPOL[j];
             if ( dim == 3U ) NT3(j,2) = IPOL[j];
          }

        NT3 *= VIP;

        // 3. Multiply with DN, integrate with weight and detJ, and add to element contribution
        // ------------------------------------------------------------------------------------
        NT3 *= DN;
        NT3 *= e.WeightAtIntegrationPoint(i) * detJ;

        MathOperatorLHS<dim,CELL>::LHS += NT3;

      }

} // end ComputeContribution




template class NumIntegral_dNT_op_dN_NT_v_dN_dV<1U,Element>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<2U,Element>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<3U,Element>;


template class NumIntegral_dNT_op_dN_NT_v_dN_dV<1U,Face>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<2U,Face>;
template class NumIntegral_dNT_op_dN_NT_v_dN_dV<3U,Face>;

} // csmp





