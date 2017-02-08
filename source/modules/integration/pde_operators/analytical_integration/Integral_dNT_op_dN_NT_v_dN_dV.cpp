#include "Integral_dNT_op_dN_NT_v_dN_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
Integral_dNT_op_dN_NT_v_dN_dV<dim,SIMPLEX>::Integral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>& pref,
                                                                      const char* oper,
                                                                      const char* velo, 
                                                                      const char* basic, 
                                                                      const char* test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    velo_key(pref.StorageKey(velo)),
    B(2,3),      
    BT(3,2)  
{
    MathOperatorLHS<dim>::Name("Integral_dNT_op_dN_NT_v_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() != REGION )
      throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      basic, "Operand must be a property placed on the element or group." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      test, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );

    if ( velo_key.place != ELEMENT || velo_key.type != VECTOR )
      throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      test, "Operand (velocity) must be a vector property placed on the element." );
}




/**
 
Reads a scalar property describing the diffusion part and a vector
variable describing the advection part of the operator.  
 */
template<size_t dim,class SIMPLEX>
void Integral_dNT_op_dN_NT_v_dN_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
 {
    // only if the property is an element property  something is done here
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT or MathOperatorLHS<dim>::MaterialOperandPlacement() == REGION ) 
      {
         MathOperatorLHS<dim>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim>::MTRL[0].Zero();
      
         if ( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR ) {
              ScalarVariable  sc;
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), sc );
              for ( size_t i=0; i<dim; i++ ) MathOperatorLHS<dim>::MTRL[0](i,i) = sc();
           }
         if ( MathOperatorLHS<dim>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), vc );
              for ( size_t i=0; i<dim; i++ ) MathOperatorLHS<dim>::MTRL[0](i,i) = vc[i];
           }
         if ( MathOperatorLHS<dim>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), ts );
              for ( size_t i=0; i<dim; i++ ) 
                for ( size_t j=0; j<dim; j++ ) MathOperatorLHS<dim>::MTRL[0](i,j) = ts(i,j);
           }
      }
      
    e.Read( velo_key, VXYZ );
 }



/**
 
@section arguments Input Arguments 

A reference to element, the contribution of which is to be aquired and
the time-increment over which the deformation shall occur. 

The result of the computation is returned into the base class protected
member matrix [C].   

@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class SIMPLEX>
void Integral_dNT_op_dN_NT_v_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    // 1. Calculating the diffusion component
    // --------------------------------------
    e.dN( B );

    // transpose B -> BT  O.K.
    B.Transposed( BT );

    // multiply BT (3x2) with DIFF(2x2) -> (3x2) 
    BT *= MathOperatorLHS<dim>::MTRL[0];
                  
    // multiply BT (3x2) with B(2x3) -> (3x3)
    BT *= B;
         
    const size_t nodes = e.Nodes();
    MathOperatorLHS<dim>::LHS.Resize(nodes,nodes);  
    // dispersion term  
    MathOperatorLHS<dim>::LHS = BT; 

    // 2. Calculating the advection component
    // --------------------------------------
    VXYZ /= nodes;
    
    for ( size_t i=0; i<nodes; i++ )
      for ( size_t j=0; j<nodes; j++ ) {
           MathOperatorLHS<dim>::LHS(i,j)   += B(0,j) * VXYZ[0];
           if ( dim != 1U ) 
             MathOperatorLHS<dim>::LHS(i,j) += B(1,j) * VXYZ[1]; 
           if ( dim == 3U ) 
             MathOperatorLHS<dim>::LHS(i,j) += B(2,j) * VXYZ[2]; 
        }    

    // 3. multiply with the area of the Element to complete integration
    // -------------------------------------------------------
    MathOperatorLHS<dim>::LHS *= e.Volume();
    
//    LHS.Out();

} // end ComputeContribution

template class Integral_dNT_op_dN_NT_v_dN_dV<1U,Element<1U> >;
template class Integral_dNT_op_dN_NT_v_dN_dV<2U,Element<2U> >;
template class Integral_dNT_op_dN_NT_v_dN_dV<3U,Element<3U> >;

template class Integral_dNT_op_dN_NT_v_dN_dV<1U,Face<1U> >;
template class Integral_dNT_op_dN_NT_v_dN_dV<2U,Face<2U> >;
template class Integral_dNT_op_dN_NT_v_dN_dV<3U,Face<3U> >;

} // csmp
