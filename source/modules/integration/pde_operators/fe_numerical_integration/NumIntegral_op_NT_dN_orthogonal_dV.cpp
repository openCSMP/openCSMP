#include "NumIntegral_op_NT_dN_orthogonal_dV.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<size_t dim,class CELL>
NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL>::NumIntegral_op_NT_dN_orthogonal_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  oper,  // fluid pressure
                                                          const char*  test ) // stream-function
  : MathOperatorRHS<dim>(pref,oper,test),
    M(dim,3), 
    DNORTHO(3,dim),
    NPROP(3),
    IPOL(3),
    NT(3,dim),
    UNITY(3),
    RES(3)
{
    MathOperatorRHS<dim>::Name("NumIntegral_op_NT_dN_orthogonal_dV", oper, test );
    
    // only 2D is possible at the moment
    if ( dim == 3 )
      throw csmp::Exception( ERROR,  "NumIntegral_op_NT_dN_orthogonal_dV<dim>::(constructor)", 
                      "At the moment this RHS operator works only in 2D" );

    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != NODE || MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR,  "NumIntegral_op_NT_dN_orthogonal_dV<dim>::(constructor)", 
                      oper, "must be an node-based scalar variable." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR,  "NumIntegral_op_NT_dN_orthogonal_dV<dim>::(constructor)", 
                      test, "must be an node-based scalar variable." );
}









/**
 
The gradients of the dependent variable are computed at the integration
points and multiplied with the n- and e-multipliers. These must be 
nodal and element variables, respectively.  

The operand is read
*/
template<size_t dim,class CELL>
void NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL>::GetOperands( CELL& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // read node variable which must be a scalar
    // --------------------------------------------
    NPROP.resize( e.Nodes() );
    for ( size_t i=0; i<e.Nodes(); ++i )
      NPROP[i] = e.N(i)->Read( MathOperatorRHS<dim>::MaterialOperandKey() );
    
} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class CELL>
void NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL>::ComputeContribution( CELL& e )
 {
    double64  detJ;

    // resizing and initializing RHS vector to zero
    MathOperatorRHS<dim>::RHS.resize( e.Nodes() );
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    // if the number of Nodes has changed since the last element
    if ( DNORTHO.Cols() != e.Nodes() ) {
         // the transformed interpolation function derivative matrix
         DNORTHO.Resize( dim, e.Nodes() );
         // node property for which the derivative shall be taken in columns
         NT.Resize( e.Nodes(), dim );
         // unity vector
         UNITY.resize( e.Nodes() );
         fill( UNITY.begin(), UNITY.end(), 1. );
         // result vector
         RES.resize( e.Nodes() );
      }

    // 1. At each integration point the nodal property is multiplied with
    //    a shape function derivative which has been modified in the 
    //    following way dNdx = -dNdy  & dNdy = dNdx
    // ------------------------------------------------------------------
    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      {
         // interpolating basic property to integration points
         e.N_AtIntegrationPoint( i, IPOL );
         for ( size_t j=0; j<e.Nodes(); j++ )
           for ( size_t k=0; k<dim; k++ ) NT(j,k) = IPOL[j] * NPROP[j];
           
         // global intpol. function derivative matrix and determinant of Jacobian matrix
         detJ = e.dN_AtIntegrationPoint( M, i, 1 );
         
         // copying scaled derivative matrix so that spatial derivatives are rotated by 90o
         for ( size_t j=0; j<e.Nodes(); j++ ) {
              DNORTHO(0,j) = -M(1,j) * detJ; // dNdx = -dNdy P
              DNORTHO(1,j) =  M(0,j) * detJ; // dNdy =  dNdx P
           }
         // multiplying with determinant and weights and condensing result into
         // single column matrix
         // multiply  NTDIM(nodes x dim) . DNORTHO(dim x nodes)
         M  = NT * DNORTHO;
             
         // multiplying with determinant and weights and condensing result into
         // single column matrix
         RES = M * UNITY;
         
         for ( size_t j=0; j<e.Nodes(); j++ ) 
           // minus since flow is always down pressure
           MathOperatorRHS<dim>::RHS[j] += -RES[j] * e.WeightAtIntegrationPoint(i);
      }

// cout <<"\nRHS for element: "<< e.Idx() << endl;
//   out( RHS );

} // end ComputeContribution




template class NumIntegral_op_NT_dN_orthogonal_dV<1U,Element<1U> >;
template class NumIntegral_op_NT_dN_orthogonal_dV<2U,Element<2U> >;
template class NumIntegral_op_NT_dN_orthogonal_dV<3U,Element<3U> >;

template class NumIntegral_op_NT_dN_orthogonal_dV<1U,Face<1U> >;
template class NumIntegral_op_NT_dN_orthogonal_dV<2U,Face<2U> >;
template class NumIntegral_op_NT_dN_orthogonal_dV<3U,Face<3U> >;

} // csmp

