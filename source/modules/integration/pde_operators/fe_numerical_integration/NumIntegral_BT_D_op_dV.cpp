#include "NumIntegral_BT_D_op_dV.h"
#include "PropertyDatabase.h"
#include "mechanics.h"
#include "Face.h"
#include "CSMP_mathUtilities.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
NumIntegral_BT_D_op_dV<dim,SIMPLEX>::NumIntegral_BT_D_op_dV( const PropertyDatabase<dim>& pref,
                                                        const char*             oper,      // strain
                                                        const char*             youngs, 
                                                        const char*             poissons, 
                                                        const char*             test,
                                                        bool plane_strain )
  : MathOperatorRHS<dim>(pref,oper,test), 
    Y_key_(pref.StorageKey( youngs )),
    nu_key_(pref.StorageKey( poissons )),
    B(3U,dim),                   // stiffness matrix
    BT(3U,dim),                  // stiffness matrix transposed
    D((dim-1U)*3U,(dim-1U)*3U),  // material property matrix
    TEMP((dim-1U)*3U,(dim-1U)*3U),
    STR((dim-1U)*3U,1U),         // column matrix
    E_(1U), nu_(1U),             // vectors of dense matrices
    plane_strain_(plane_strain)
{
    MathOperatorRHS<dim>::Name("NumIntegral_BT_D_op_dV", oper, test );

    // verify here that the operands have the correct placement and type 
    if ( Y_key_.place != nu_key_.place ) 
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)", 
                                 "Young's modulus variable must have same placement as Poisson's ratio.");

   if ( (Y_key_.place != ELEMENT or Y_key_.place != ELEMENT_INTEGRATION_POINT) && Y_key_.type != SCALAR ) {
        throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)", 
                               "Modulus variable must be a scalar placed on the element.");
     }
   if ( (nu_key_.place != ELEMENT or nu_key_.place != ELEMENT_INTEGRATION_POINT) && nu_key_.type != SCALAR ) {
        throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)", 
                               "Poisson's ratio must be a scalar placed on the element.");
     }
   if ( nu_key_.place != Y_key_.place ) {
        throw csmp::Exception( ERROR, "NumIntegral_BT_D_op_dV<dim>::(constructor)", 
                               "Young's modulus and Poisson's ratio must have the same placement.");
     }
     
} // end constructor




template<size_t dim,class SIMPLEX>
void NumIntegral_BT_D_op_dV<dim,SIMPLEX>::PlaneStress() { plane_strain_ = false; }







/**
 
Reads the Operand variable 'strain', the Young's modulus and Poisson's 
ratio. Only the x, y, z components of the strain will be used. Shear
strains are not considered.  

The righthand contribution:  {F} = ElementIntegral [BT][D]{e} 

With the strain source term: {e} = strain * {1,1,0} (2D) 

is accumulated into the righthand vector. Strain is read as a scalar
or vector variable from the element and then stored into the
E_OP vector of dimension  nodes-per-element x dim.  

@section arguments Input Arguments 

A reference to element, the contribution of which is to be aquired and
the time-increment over which the deformation shall occur. 
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_BT_D_op_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );
 
    // Young's modulus and Poisson's ratio
    if ( Y_key_.place == ELEMENT ) {
		     // reading Young's modulus (must be an element variables)
		     E_[0U] = e.Read( Y_key_ );
		     // getting Poisson's ratio
	       nu_[0U] = e.Read( nu_key_ );
      }
    else { // Node or Constraint Point
         E_.resize(e.FE()->IntegrationPoints());
         nu_.resize(e.FE()->IntegrationPoints());
         for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) {
              MathOperatorRHS<dim>::PropertyAtIntegrationPoint( e, Y_key_, i, E_[i] );
              MathOperatorRHS<dim>::PropertyAtIntegrationPoint( e, nu_key_, i, nu_[i] );
           }
      }
 
    // input strain / dilatation
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) {
         // strain
         if ( MathOperatorRHS<dim>::MaterialOperandType() == SCALAR ) {
              e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), sc_ );
              MathOperatorRHS<dim>::MTRL[0].AssignToDiagonal( dim, sc_ );
           }
         else if ( MathOperatorRHS<dim>::MaterialOperandType() == VECTOR ) {
              e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), vc_ );
              MathOperatorRHS<dim>::MTRL[0].AssignToDiagonal( vc_ );
           }
         else if ( MathOperatorRHS<dim>::MaterialOperandType() == TENSOR ) {
              e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), ts_ );
              MathOperatorRHS<dim>::MTRL[0] = ts_;
           }
      }
    else if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
      {
         MathOperatorRHS<dim>::MTRL.resize(e.IntegrationPoints());
         for ( size_t i=0; i<e.IntegrationPoints(); i++ )
           for ( size_t j=0U; j<dim; ++j )
             MathOperatorRHS<dim>::MTRL[i](j,j) = e.Read( i, MathOperatorRHS<dim>::MaterialOperandKey() );
      }
    else // if a nodal variable is dealt with
      {
         for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
          MathOperatorRHS<dim>::PropertyAtIntegrationPoint( e,
                                         MathOperatorRHS<dim>::MaterialOperandKey(), 
                                         i, MathOperatorRHS<dim>::MTRL[i] );
      }
  
} // end GetOperands




/**
 
Computes the spatial integral over the matrix product:

[B]^T [JI]^T [D] {e} detJ

as the elements righthand contribution.
The integration is performed numerically.  

@section arguments Input Arguments 

A reference to element, the contribution of which is to be aquired and
the time-increment over which the deformation shall occur. 

The result of the computation is returned into the base class protected
member vector {V}.   

@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_BT_D_op_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    // compute the material property matrix
    // compute the material property matrix based on element properties
    if ( Y_key_.place == ELEMENT ) {
         if ( dim == 2U ) {
              // assuming isotropic Young's modulus 
              if ( plane_strain_ ) 
                planeStrainMatrix( E_[0](0U,0U), nu_[0](0U,0U), D );
              else                
                planeStressMatrix( E_[0](0U,0U), nu_[0](0U,0U), D );
           }
         else   stiffnessMatrix( E_[0](0U,0U), nu_[0](0U,0U), D );
      }

    // initialize output vector
    MathOperatorRHS<dim>::RHS.resize( dim*e.Nodes() );
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );

    // mapping isostatic components of strain vector into 3x1 matrix STR
    STR.Zero();
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT )
      for ( size_t i=0; i<dim; i++ )
        STR(i,0) = MathOperatorRHS<dim>::MTRL[0](i,i);

    // numerical integration: 
    // ----------------------
    // looping over the 3 Gauss points calculating matrix products
    // and applying uniform weights (1/3) before adding integrated 
    // matrices to element - contribution matrix
    for ( size_t i=0U; i<e.FE()->IntegrationPoints(); i++ )
      {
         // the material property matrix is constructed at each integration point
         if ( Y_key_.place == ELEMENT_INTEGRATION_POINT or Y_key_.place == NODE ) {
              if ( dim == 2U ) {
                   // assuming isotropic Young's modulus 
                   if ( plane_strain_ ) 
                     planeStrainMatrix( E_[i](0U,0U), nu_[i](0U,0U), D );
                   else                
                     planeStressMatrix( E_[i](0U,0U), nu_[i](0U,0U), D );
                }
              else   stiffnessMatrix( E_[i](0U,0U), nu_[i](0U,0U), D );
           }

         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         double64  detJ = e.dN_AtIntegrationPoint( B, i, dim );

         if ( detJ <= 0. ) {
              cout <<"\nDeterminant of Jacobian at Gauss point: "<< i <<": "<< detJ << endl;
              throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_op_dV<dim>::ComputeContribution",
                              "Jacobian transformation failed.");
           }
         // transposing B -> BT 
         B.Transposed( BT );

         // computing isometric volume strain term TEMP(3x3) STR(3x1) -> TEMP(3x1)
         if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT or 
              MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE ) {
              STR.Zero();
              for ( size_t j=0; j<dim; j++ ) 
                STR(j,0U) = MathOperatorRHS<dim>::MTRL[i](j,j);
           }

         TEMP  = D; 
         TEMP *= STR; 

         // multiply BT(12x3) TEMP(3x1) -> BT(12x1)
         BT *= TEMP;

         // multiplying with determinant and weights
         for ( size_t n=0U; n<BT.Rows(); n++ ) BT(n,0U) *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // adding to result vector
         for ( size_t n=0U; n<BT.Rows(); n++ ) 
           MathOperatorRHS<dim>::RHS[n] += BT(n,0U);
      }

} // end ComputeContribution



template class NumIntegral_BT_D_op_dV<1U,Element<1U> >;
template class NumIntegral_BT_D_op_dV<2U,Element<2U> >;
template class NumIntegral_BT_D_op_dV<3U,Element<3U> >;

template class NumIntegral_BT_D_op_dV<1U,Face<1U> >;
template class NumIntegral_BT_D_op_dV<2U,Face<2U> >;
template class NumIntegral_BT_D_op_dV<3U,Face<3U> >;

} // csmp



