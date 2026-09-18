// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_op_NT_dN_orthogonal_dV.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_op_NT_dN_orthogonal_dV<dim,CELL>::Integral_op_NT_dN_orthogonal_dV(
                                                          const PropertyDatabase<dim>& pref,
                                                          const char*  oper,  // fluid pressure
                                                          const char*  test ) // stream-function
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    M(dim,3), 
    DNORTHO(3,dim),
    NPROP(3),
    IPOL(3),
    NT(3,dim),
    UNITY(3),
    RES(3)
{
    MathOperatorRHS<dim,CELL>::Name("Integral_op_NT_dN_orthogonal_dV", oper, test );
    
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR,  "Integral_op_NT_dN_orthogonal_dV<dim>::(constructor)", 
                             oper, "must be an node-based scalar variable." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR,  "Integral_op_NT_dN_orthogonal_dV<dim>::(constructor)", 
                             test, "must be an node-based scalar variable." );
}




/**
 
The gradients of the dependent variable are computed at the integration
points and multiplied with the n- and e-multipliers. These must be 
nodal and element variables, respectively.  

The operand is read
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_op_NT_dN_orthogonal_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    // read node variable which must be a scalar
    NPROP.resize( e.Nodes() );
    for ( uint32_t i{0U}; i<e.Nodes(); i++ )
      NPROP[i] = e.N(i)->Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey() );
    
} // end GetOperands






/**
 
@section application Application

In linear elasticity computations.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_op_NT_dN_orthogonal_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // resizing and initializing RHS vector to zero
    MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() );
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), zero );

    // if the number of Nodes has changed since the last element
    if ( DNORTHO.Cols() != e.Nodes() ) {
         // the transformed interpolation function derivative matrix
         DNORTHO.Resize( dim, e.Nodes() );
         // node property for which the derivative shall be taken in columns
         NT.Resize( e.Nodes(), dim );
         // unity vector
         UNITY.resize( e.Nodes() );
         fill( UNITY.begin(), UNITY.end(), one );
         // result vector
         RES.resize( e.Nodes() );
      }

    // 1. At the element barycenter the nodal property is multiplied with
    //    a shape function derivative which has been modified in the 
    //    following way dNdx = -dNdy  & dNdy = dNdx
    // ------------------------------------------------------------------
     // interpolating basic property to integration points
     e.N_AtBaryCenter( IPOL );
     for ( uint32_t j{0U}; j<e.Nodes(); j++ )
       for ( uint32_t k=0; k<dim; k++ ) NT(j,k) = IPOL[j] * NPROP[j];
       
     // global intpol. function derivative matrix and determinant of Jacobian matrix
     const double volume(e.Volume());
     
     // copying scaled derivative matrix so that spatial derivatives are rotated by 90o
     if ( dim == 2U ) for ( auto j{0U}; j<e.Nodes(); j++ ) {
          DNORTHO(0,j) = -M(1,j) * volume; // dNdx = -dNdy P
          DNORTHO(1,j) =  M(0,j) * volume; // dNdy =  dNdx P
       }

     else if ( dim == 3U ) for ( auto j{0U}; j<e.Nodes(); j++ ) {
          DNORTHO(0,j) = -M(1,j) * volume; // dNdx = -dNdy P
          DNORTHO(1,j) =  M(0,j) * volume; // dNdy =  dNdx P
          if ( dim == 3U )
            DNORTHO(2,j) =  M(0,j) * volume; // dNdz =  dNdz P
       }
     else
     cout <<"\nIntegral_op_NT_dN_orthogonal_dV<"<< dim <<"> undefined."<< endl;

     // multiplying with determinant and weights and condensing result into
     // single column matrix
     // multiply  NTDIM(nodes x dim) . DNORTHO(dim x nodes)
     M  = NT * DNORTHO;
         
     // multiplying with determinant and weights and condensing result into
     // single column matrix
     RES = M * UNITY;
     
     for ( auto j{0U}; j<e.Nodes(); j++ ) 
       // minus since flow is always down pressure
       MathOperatorRHS<dim,CELL>::RHS[j] += -RES[j];

// cout <<"\nRHS for element: "<< e.Idx() << endl;
//   out( RHS );

} // end ComputeContribution





template class Integral_op_NT_dN_orthogonal_dV<2U,Element>;
template class Integral_op_NT_dN_orthogonal_dV<3U,Element>;

template class Integral_op_NT_dN_orthogonal_dV<2U,Face>;
template class Integral_op_NT_dN_orthogonal_dV<3U,Face>;

} // csmp





