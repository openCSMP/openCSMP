#include "Integral_dNT_op_dN_NT_v_dN_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {


template<uint32_t dim, template<uint32_t> class CELL>
Integral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::Integral_dNT_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>& pref,
                                                                        const char* oper,
                                                                        const char* velo,
                                                                        const char* basic,
                                                                        const char* test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    velo_key(pref.StorageKey(velo)),
    B(2,3),      
    BT(3,2)  
{
    MathOperatorLHS<dim,CELL>::Name("Integral_dNT_op_dN_NT_v_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != REGION )
      throw csmp::Exception( ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      basic, "Operand must be a property placed on the element or group." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      test, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );

    if ( velo_key.place != ELEMENT || velo_key.type != VECTOR )
      throw csmp::Exception( ERROR, "Integral_dNT_op_dN_NT_v_dN_dV::(constructor)", 
                      test, "Operand (velocity) must be a vector property placed on the element." );
}




/**
 
Reads a scalar property describing the diffusion part and a vector
variable describing the advection part of the operator.  
 */
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    // only if the property is an element property  something is done here
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT or
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == REGION )
      {
         MathOperatorLHS<dim,CELL>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim,CELL>::MTRL[0].Zero();
      
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
              ScalarVariable  sc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), sc );
              for ( auto i{0U}; i<dim; i++ ) MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = sc();
           }
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
              for ( auto i{0U}; i<dim; i++ ) MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = vc[i];
           }
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
              for ( auto i{0U}; i<dim; i++ ) 
                for ( auto j{0U}; j<dim; j++ ) MathOperatorLHS<dim,CELL>::MTRL[0](i,j) = ts(i,j);
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
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_op_dN_NT_v_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // 1. Calculating the diffusion component
    // --------------------------------------
    e.dN( B );

    // transpose B -> BT  O.K.
    B.Transposed( BT );

    // multiply BT (3x2) with DIFF(2x2) -> (3x2) 
    BT *= MathOperatorLHS<dim,CELL>::MTRL[0];
                  
    // multiply BT (3x2) with B(2x3) -> (3x3)
    BT *= B;
         
    const auto nodes = e.Nodes();
    MathOperatorLHS<dim,CELL>::LHS.Resize(nodes,nodes);
    // dispersion term  
    MathOperatorLHS<dim,CELL>::LHS = BT;

    // 2. Calculating the advection component
    // --------------------------------------
    VXYZ /= nodes;
    
    for ( auto i{0U}; i<nodes; i++ )
      for ( auto j{0U}; j<nodes; j++ ) {
           MathOperatorLHS<dim,CELL>::LHS(i,j)   += B(0,j) * VXYZ[0];
           if ( dim != 1U ) 
             MathOperatorLHS<dim,CELL>::LHS(i,j) += B(1,j) * VXYZ[1];
           if ( dim == 3U ) 
             MathOperatorLHS<dim,CELL>::LHS(i,j) += B(2,j) * VXYZ[2];
        }    

    // 3. multiply with the area of the Element to complete integration
    // -------------------------------------------------------
    MathOperatorLHS<dim,CELL>::LHS *= e.Volume();
    
//    LHS.Out();

} // end ComputeContribution

template class Integral_dNT_op_dN_NT_v_dN_dV<1U,Element>;
template class Integral_dNT_op_dN_NT_v_dN_dV<2U,Element>;
template class Integral_dNT_op_dN_NT_v_dN_dV<3U,Element>;

template class Integral_dNT_op_dN_NT_v_dN_dV<1U,Face>;
template class Integral_dNT_op_dN_NT_v_dN_dV<2U,Face>;
template class Integral_dNT_op_dN_NT_v_dN_dV<3U,Face>;

} // csmp
