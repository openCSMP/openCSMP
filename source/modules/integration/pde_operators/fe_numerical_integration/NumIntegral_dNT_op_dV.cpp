#include "NumIntegral_dNT_op_dV.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<uint32_t dim,class CELL>
NumIntegral_dNT_op_dV<dim,CELL>::NumIntegral_dNT_op_dV( const PropertyDatabase<dim>& pref,
                                                        const char*             oper,
                                                        const char*             test )
  : MathOperatorRHS<dim>(pref,oper,test),
    B_(dim,3), BT_(3,dim)
{
    MathOperatorRHS<dim>::Name("NumIntegral_dNT_op_dV", oper, test );
    
    if ( !(MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ||
           MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
           MathOperatorRHS<dim>::MaterialOperandPlacement() == FACE) ||
           MathOperatorRHS<dim>::MaterialOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dV<dim>::(constructor)", 
                      oper, "Operand must be a vector property placed on the element, face or element integration point." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_op_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}



template<uint32_t dim,class CELL>
void NumIntegral_dNT_op_dV<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.UsesLocalCoordinates() == true );

    // initialize output matrix
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    // if the finite element is a linear simplex element, its Jacobian and element-interpolation derivative matrix
    // are constant throughout it
    /* TODO: DOES NOT WORK YET
    const bool is_simplex_element_type(e.FE()->IsSimplex() && e.Interpolation() == 1 );
    const bool piecewise_constant_material( this->MaterialOperandPlacement() == ELEMENT ||
                                            this->MaterialOperandPlacement() == FACE    ||
                                            this->MaterialOperandPlacement() == REGION );
    const uint32_t n_nodes{ e.Nodes() };
    if ( is_simplex_element_type && piecewise_constant_material ) {
         const double detJ = e.dN_AtBaryCenter( B_ );
         // transposing B -> BT  O.K.
         B_.Transposed( BT_ );
         // multiply  BT . MTRL
         BT_ *= MathOperatorRHS<dim>::MTRL[0];
         // multiplying BT . B
         BT_ *= B_;
         // multiplying with determinant and weights (ASSUMING that for simplices these weights are all the same)
         BT_ *= e.WeightAtIntegrationPoint(0) * e.IntegrationPoints() * detJ;
         // row sum diagonalisation of matrix into right-hand vector
         for ( auto k=0; k < dim; k++ )
           for ( auto j{0U}; j<n_nodes; j++ )
             MathOperatorRHS<dim>::RHS[j] += BT_(j,k);
         return;
      }
    */
    
    //  When this is not a simplex element but the material property is an element property
    // ------------------------------------------------------------------------------------
    const uint32_t n_nodes{ e.Nodes() };
    if ( this->MaterialOperandPlacement() == ELEMENT || this->MaterialOperandPlacement() == FACE )
      {
        for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ )
          {
             // getting global intpol. function derivative matrix and determinant of
             // byproduct Jacobian matrix (B is already in global coordinates)
             double detJ = e.dN_AtIntegrationPoint( B_, i, SCALAR );

             // transposing B -> BT 
             B_.Transposed( BT_ );

             // multiply  BT . MTRL
             BT_ *= this->MTRL[0];
                 
             // multiplying with determinant and weights
             BT_ *= e.WeightAtIntegrationPoint(i) * detJ;

             // row sum diagonalisation of matrix into right-hand vector
             for ( auto k{0U}; k < dim; k++ )
               for ( auto j{0U}; j<n_nodes; j++ )
                 MathOperatorRHS<dim>::RHS[j] += BT_(j,k);
          }
        return;
      }
      
    if ( this->MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
      {
         for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ )
           {
              double detJ = e.dN_AtIntegrationPoint( B_, i, 1 );
              B_.Transposed( BT_ );
              BT_ *= MathOperatorRHS<dim>::MTRL[i];
              BT_ *= e.WeightAtIntegrationPoint(i) * detJ;
              for ( auto k{0U}; k < dim; k++ )
                for ( auto j{0U}; j<n_nodes; j++ )
                  MathOperatorRHS<dim>::RHS[j] += BT_(j,k);
           }
       }

} // end ComputeContribution


template class NumIntegral_dNT_op_dV<1U,Element<1U> >;
template class NumIntegral_dNT_op_dV<2U,Element<2U> >;
template class NumIntegral_dNT_op_dV<3U,Element<3U> >;

template class NumIntegral_dNT_op_dV<1U,Face<1U> >;
template class NumIntegral_dNT_op_dV<2U,Face<2U> >;
template class NumIntegral_dNT_op_dV<3U,Face<3U> >;

} // csmp











