#include "NumIntegral_BT_D_B_dV.h"
#include "PropertyDatabase.h"
#include "mechanics.h"
#include "CSMP_mathUtilities.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
NumIntegral_BT_D_B_dV<dim,SIMPLEX>::NumIntegral_BT_D_B_dV( const PropertyDatabase<dim>& pref,
                                                           const char*             oper,  // Young's modulus
                                                           const char*             oper2, // Poisson's ratio 
                                                           const char*             basic, 
                                                           const char*             test,
                                                           bool plane_strain )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    nu_key_(pref.StorageKey( oper2 )),
    D(3,3), B(2,3), BT(3,2), nu_(1U), E_(1U),
    plane_strain_(plane_strain)
{
    MathOperatorLHS<dim>::Name("NumIntegral_BT_D_B_dV", oper, basic, test );
    // get other csmp::Index keys
    
    // Poisson's ratio
    if ( nu_key_.type != SCALAR or nu_key_.place == NODE or nu_key_.place == FACE )                  
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      oper2, "must be a scalar property placed on the constraint points or elements.");
      
    // Young's modulus
    if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR )                  
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      oper, "Currently the Operand must be a scalar property.");
    // displacement
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      basic, "Operand (basic) must be a vector property placed on the nodes.");
    // displacement
    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      test, "Operand (test) must be a vector property placed on the nodes.");

    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != nu_key_.place )                  
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                     "Both Operands must have the same placement.");
}



/**
 
Allows to toggle between plane stress and plane strain computations in
two-dimensional models. Clearly, this method only has an effect in 2D
calculations and the chosen computation must match the preceding
Algorithm.  

@section arguments Input Arguments 

A boolean variable specifying whether the material property matrix
used shall be initialized for plane strain or plane stress.   

@section messages Messages 

If the method is called in a 3D calculation, the user is warned that
it will have no effect.  
 */
template<size_t dim,class SIMPLEX>
void NumIntegral_BT_D_B_dV<dim,SIMPLEX>::PlaneStress( bool yes_no )
 { 
     plane_strain_ = yes_no; 
 }









/**
 
While Poisson's ratio must be an element variable, this method allows for
Young's modulus to be a node variable. Thus, a continuous loss of strength
can be modeled.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_BT_D_B_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.Isoparametric() == true );

    // only if the property is an element property  something is done here
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) {
         // we made sure that youngs modulus is a scalar and has the same placement as Poisson's ratio
         E_[0]  = e.Read( MathOperatorLHS<dim>::MaterialOperandKey() );
         nu_[0] = e.Read( nu_key_ );
      }
    else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) {
        const size_t ipoints(e.IntegrationPoints());
        E_.resize(ipoints);
        nu_.resize(ipoints);
        for ( size_t i=0; i<ipoints; ++i )
          {
             // we made sure that youngs modulus is a scalar and has the same placement as Poisson's ratio
             E_[i]  = e.Read( i, MathOperatorLHS<dim>::MaterialOperandKey() );
             nu_[i] = e.Read( i, nu_key_ );
          }
       }
    else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == NODE )
      {
         ScalarVariable sc;
         for ( size_t i=0U; i<e.IntegrationPoints(); i++ ) {
             e.PropertyValueAtIntegrationPoint( MathOperatorLHS<dim>::MaterialOperandKey(), i, sc );
             E_[i] = sc();
             e.PropertyValueAtIntegrationPoint( nu_key_, i, sc );
             nu_[i] = sc();
          }
      }
     else throw csmp::Exception( WARNING, "umIntegral_BT_D_B_dV<dim,SIMPLEX>::GetOperands",
                                 MathOperatorLHS<dim>::MaterialOperandName().c_str(), "placement of material operand not handled yet." );
      
} // end GetOperands





/**
 
Computes the 2D spatial integral over the matrix product:

[B]^T [JI]^T [D] [B] [JI] detJ

as the elements lefthand contribution to the global stiffness matrix.
The integration is performed numerically using a 3 Gauss-Point scheme
in the case of the quadratic triangular finite element.  

@section arguments Input Arguments 

A reference to element, the contribution of which is to be aquired and
the time-increment over which the deformation shall occur. 

The result of the computation is returned into the base class protected
member matrix [C].   

@section application Application

In linear elasticity computations.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_BT_D_B_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( dim*e.Nodes(), dim*e.Nodes() );

    // compute the material property matrix based on element properties
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) {
         if ( dim == 2 ) {
              // assuming isotropic Young's modulus 
              if ( plane_strain_ ) 
                planeStrainMatrix( E_[0], nu_[0], D );
              else                
                planeStressMatrix( E_[0], nu_[0], D );
           }
         else if( dim == 3)   
           stiffnessMatrix( E_[0], nu_[0], D );
         else if( dim == 1 )
           {
              stiffnessMatrix( E_[0],  D, e.Volume() );
              MathOperatorLHS<dim>::LHS = D;
              return;
           }

      }
      
    // analytical integration
    // -----------------------------
    if ( !e.FE()->UsesLocalCoordinates() )
      {
         const double64 volume(e.Volume());
         // setting C to 1 and its diagonal to 2
         MathOperatorLHS<dim>::LHS = volume / 12.;
         for ( size_t f=0; f<(e.Nodes()*dim); f++ ) 
           MathOperatorLHS<dim>::LHS(f,f) = volume / 6.;
      
         e.dN( B );
         dN_To2DOF( e.Nodes(), B );
      
         // transposing B -> BT  O.K.
         B.Transposed( BT );

         // multiplying BT(12x3) D(3x3) B(2x6) -> (12x12)  O.K.
         BT  *= D;
         BT  *= B;
         
         // creating the integrated output matrix
         MathOperatorLHS<dim>::LHS *= BT; 
          
         return;
      } 

    // numerical integration: 
    // ----------------------
    // looping over the 3 Gauss points calculating matrix products
    // and applying uniform weights (1/3) before adding integrated 
    // matrices to element - contribution matrix
    MathOperatorLHS<dim>::LHS.Zero();

    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      {
         // the material property matrix is constructed at each integration point
         if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) {
              if ( dim == 2 ) {
                   // assuming isotropic Young's modulus 
                   if ( plane_strain_ ) 
                     planeStrainMatrix( E_[i], nu_[i], D );
                   else                
                     planeStressMatrix( E_[i], nu_[i], D );
                }
              else   stiffnessMatrix( E_[i], nu_[i], D );
           }

         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         double64 detJ = e.dN_AtIntegrationPoint( B, i, dim );

         if ( detJ <= 0. ) {
              cerr <<"\n\tElement "<< e.Idx() <<": determinant of Jacobian at Gauss point "<< i <<": "<< detJ << endl;
              throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_B_dV<dim>::ComputeContribution",
                             "Jacobian transformation failed. Element nodes are perhaps not numbered correctly.");
           }

         // transposing B -> BT  O.K. (since matrix multiplication is associative)
         B.Transposed( BT );

         // multiply BT(12x3) D(3x3) 
         BT *= D;

         // multiplying BT'(12x3) B(3x12) -> RES(12x12)  O.K.
         BT *= B;

         // multiplying with determinant and weights
         BT *= e.WeightAtIntegrationPoint(i) * detJ; 
         
         // accumulating ME Gauss point integral contributions into element 
         // contribution to global conductance matrix
         MathOperatorLHS<dim>::LHS += BT;
      }

} // end ComputeContribution


template class NumIntegral_BT_D_B_dV<1U,Element<1U> >;
template class NumIntegral_BT_D_B_dV<2U,Element<2U> >;
template class NumIntegral_BT_D_B_dV<3U,Element<3U> >;

template class NumIntegral_BT_D_B_dV<1U,Face<1U> >;
template class NumIntegral_BT_D_B_dV<2U,Face<2U> >;
template class NumIntegral_BT_D_B_dV<3U,Face<3U> >;

} // csmp
