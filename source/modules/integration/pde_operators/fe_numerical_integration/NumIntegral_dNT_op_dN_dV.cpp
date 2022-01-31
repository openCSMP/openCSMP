#include "NumIntegral_dNT_op_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

/**
 
The Operand which is used here can be both, an element or a nodal variable
which is then interpolated to the integration points to obtain the 
integral properties.  
*/
template<size_t dim,class CELL>
NumIntegral_dNT_op_dN_dV<dim,CELL>::NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                              const char*           oper,
                                                              const char*           basic,
                                                              const char*           test )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    B_(dim,3), BT_(3,dim)
{
    MathOperatorLHS<dim>::Name("NumIntegral_dNT_op_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}










/** Laplacian operator of shape function derivatives squared.

    @attention since the interpolation functions derivatives are constant across simplices (linear line, triangle and tetrahedral elements),
    the computation of DN for these elements simplifies greatly.
*/
template<size_t dim,class CELL>
void NumIntegral_dNT_op_dN_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    // if the agregated finite element is a simplex, the Jacobian and element-interpolation derivative matrix is constant throughout it
    const bool is_simplex_element_type(e.FE()->IsSimplex() && e.Interpolation() == 1);
    const bool piecewise_constant_material(this->MaterialOperandPlacement() == ELEMENT or
                                           this->MaterialOperandPlacement() == REGION or
                                           this->MaterialOperandPlacement() == FACE);
   
    if ( is_simplex_element_type && piecewise_constant_material ) {
         const double detJ = e.dN_AtBaryCenter( B_ );
         // transposing B -> BT  O.K.
         B_.Transposed( BT_ );
         // multiply  BT . MTRL
         BT_ *= MathOperatorLHS<dim>::MTRL[0];
         // multiplying BT . B 
         BT_ *= B_;
         // multiplying with determinant and weights (ASSUMING that these weights are all the same for simplices)
         BT_ *= e.WeightAtIntegrationPoint(0) * e.IntegrationPoints() * detJ;
         MathOperatorLHS<dim>::LHS += BT_;
         return;
      }

    // 1. Two cases exist: The first is when the material property is an
    //    element property. In this case the material property matrix can
    //    be used as is.
    // ------------------------------------------------------------------
    if ( piecewise_constant_material )
      {
        for ( size_t i=0U; i<e.IntegrationPoints(); i++ ) {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             const double detJ = e.dN_AtIntegrationPoint( B_, i, SCALAR );
             // transposing B -> BT
             B_.Transposed( BT_ );
             // multiply  BT . MTRL
             BT_ *= MathOperatorLHS<dim>::MTRL[0];
             // multiplying BT . B
             BT_ *= B_;
             // multiplying with determinant and weights
             BT_ *= e.WeightAtIntegrationPoint(i) * detJ;
             // accumulating ME Gauss point integral contributions into element
             // contribution to global conductance matrix
             MathOperatorLHS<dim>::LHS += BT_;
          }
      }
    else { // NODE or ELEMENT_INTEGRATION_POINT material placements
        double detJ = ( is_simplex_element_type ) ? e.dN_AtBaryCenter( B_ ) : 0.;
        for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ ) {
             if ( !is_simplex_element_type ) detJ = e.dN_AtIntegrationPoint( B_, i, SCALAR );
             B_.Transposed( BT_ );
             BT_ *= MathOperatorLHS<dim>::MTRL[i];
             BT_ *= B_;
             BT_ *= e.WeightAtIntegrationPoint(i) * detJ;
             MathOperatorLHS<dim>::LHS += BT_;
          }
     }

} // end ComputeContribution


//cout <<"\nNumIntegral_dNT_op_dN_dV: on Element "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( size_t i=0U; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nNumIntegral_dNT_op_dN_dV: zero element in diagonal of element matrix."; 


template class NumIntegral_dNT_op_dN_dV<1U,Element<1U> >;
template class NumIntegral_dNT_op_dN_dV<2U,Element<2U> >;
template class NumIntegral_dNT_op_dN_dV<3U,Element<3U> >;

template class NumIntegral_dNT_op_dN_dV<1U,Face<1U> >;
template class NumIntegral_dNT_op_dN_dV<2U,Face<2U> >;
template class NumIntegral_dNT_op_dN_dV<3U,Face<3U> >;

template class NumIntegral_dNT_op_dN_dV<1U,InterFace<1U> >;
template class NumIntegral_dNT_op_dN_dV<2U,InterFace<2U> >;
template class NumIntegral_dNT_op_dN_dV<3U,InterFace<3U> >;

} // csmp











