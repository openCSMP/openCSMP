#include "NumIntegral_BT_op_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class CELL>
NumIntegral_BT_op_dV<dim,CELL>::NumIntegral_BT_op_dV( const PropertyDatabase<dim>& pref,
                                                      const char*  oper, // pore pressure
                                                      const char*  test )
  : MathOperatorRHS<dim>(pref,oper,test), 
    B(2,6),  
    BT(6,2),
    STR(3,1) 
{
    MathOperatorRHS<dim>::Name("NumIntegral_BT_op_dV", oper, test );
    
    // verify here that the operands have the correct placement and type 
   if ( MathOperatorRHS<dim>::MaterialOperandType() == TENSOR ) 
     throw csmp::Exception( WARNING, "NumIntegral_BT_op_dV::(constructor)", 
                     oper, "Only diagonal part of tensor property will be used.");
     
   // checking how many variables are needed to store the components of the 
   // scalar, vector, or symmetric tensor variable
   if ( dim == 2 ) STR.Resize( 3, 1 ); // stress tensor components 
   else            STR.Resize( 6, 1 );

} // end constructor







/**
 
Computes the spatial (volume) integral over the matrix product:

sum_over_integration_points [B]^T [JI]^T [property] detJ

as the elements righthand contribution. The integration is performed 
numerically.  

@section arguments Input Arguments 

A reference to element, the contribution of which is to be aquired and
the time-increment over which the deformation shall occur. 

The result of the computation is returned into the base class protected
member Vec {V}.   

@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class CELL>
void NumIntegral_BT_op_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // initialize output vector
    MathOperatorRHS<dim>::RHS.resize( dim * e.Nodes() );
    
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), static_cast<double>(0.) );

    // mapping element-placed operand from material matrix into 1-column matrix
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) {
         for ( size_t i=0; i<dim; i++ ) STR(i,0) = MathOperatorRHS<dim>::MTRL[0](i,i);
         if ( dim == 2U ) STR(2,0) = MathOperatorRHS<dim>::MTRL[0](0,1);
         else { // 3D, upper diagonal elements of symmetric tensor
              STR(3,0) = MathOperatorRHS<dim>::MTRL[0](0,1); // xy 
              STR(4,0) = MathOperatorRHS<dim>::MTRL[0](1,2); // yz
              STR(5,0) = MathOperatorRHS<dim>::MTRL[0](0,2); // zx
           }
      }

    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         double detJ = e.dN_AtIntegrationPoint( B, i, dim );

         if ( detJ <= 0. ) {
              cout <<"\nDeterminant of Jacobian at Gauss point: "<< i <<": "<< detJ << endl;
              throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_op_dV::ComputeContribution",
                              "Jacobian transformation failed.");
           }

         // transposing B -> BT 
         BT.Resize( B.Cols(), B.Rows() );
         B.Transposed( BT );

         // mapping node or cpoint-placed operand from material matrix into 1-column matrix
         if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE or 
              MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) { 
             if ( dim == 1U ) {
	                STR(0,0) = MathOperatorRHS<dim>::MTRL[i](0,0);
               }
	           else if ( dim == 2U ) {
	                STR(0,0) = MathOperatorRHS<dim>::MTRL[i](0,0);
	                STR(1,0) = MathOperatorRHS<dim>::MTRL[i](1,1);
	                STR(2,0) = MathOperatorRHS<dim>::MTRL[i](0,1);
	             }
	           else { // 3D, upper diagonal elements of symmetric tensor
	                STR(0,0) = MathOperatorRHS<dim>::MTRL[i](0,0); // diagonal terms
	                STR(1,0) = MathOperatorRHS<dim>::MTRL[i](1,1);
	                STR(2,0) = MathOperatorRHS<dim>::MTRL[i](2,2);
	                STR(3,0) = MathOperatorRHS<dim>::MTRL[i](0,1); // xy 
	                STR(4,0) = MathOperatorRHS<dim>::MTRL[i](1,2); // yz
	                STR(5,0) = MathOperatorRHS<dim>::MTRL[i](0,2); // zx
	             }
           }
           
         // multiply BT(12x3) TEMP(3x1) -> BT(12x1)
         BT *= STR;

         // multiplying with determinant and weights
         for ( size_t n=0; n<BT.Rows(); n++ ) 
           BT(n,0) *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // adding to result vector
         for ( size_t n=0; n<BT.Rows(); n++ ) 
           MathOperatorRHS<dim>::RHS[n] += BT(n,0);
      }

//  cout <<"\nodal forces due to pore pressure, element: "<< e.Idx() << endl;
//  MathOperatorRHS<dim>::RHS.out();
//  cout << endl;

} // end ComputeContribution



template class NumIntegral_BT_op_dV<1U,Element<1> >;
template class NumIntegral_BT_op_dV<2U,Element<2> >;
template class NumIntegral_BT_op_dV<3U,Element<3> >;

template class NumIntegral_BT_op_dV<1U,Face<1> >;
template class NumIntegral_BT_op_dV<2U,Face<2> >;
template class NumIntegral_BT_op_dV<3U,Face<3> >;

} // end namespace csmp









