// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_BT_D_B_dV.h"
#include "PropertyDatabase.h"
#include "mechanics.h"
#include "CSMP_mathUtilities.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_BT_D_B_dV<dim,CELL>::NumIntegral_BT_D_B_dV( const PropertyDatabase<dim>& pref,
                                                        const char*             oper,  // Young's modulus
                                                        const char*             oper2, // Poisson's ratio
                                                        const char*             basic,
                                                        const char*             test,
                                                        bool plane_strain )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    nu_key_(pref.StorageKey( oper2 )),
    D(3,3), B(2,3), BT(3,2), nu_(1U), E_(1U),
    plane_strain_(plane_strain)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_BT_D_B_dV", oper, basic, test );
    // get other csmp::Index keys
    
    // Poisson's ratio
    if ( nu_key_.type != SCALAR or nu_key_.place == NODE or nu_key_.place == FACE )                  
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      oper2, "must be a scalar property placed on the constraint points or elements.");
      
    // Young's modulus
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      oper, "Currently the Operand must be a scalar property.");
    // displacement
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      basic, "Operand (basic) must be a vector property placed on the nodes.");
    // displacement
    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_BT_D_B_dV::(constructor)", 
                      test, "Operand (test) must be a vector property placed on the nodes.");

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != nu_key_.place )
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
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_B_dV<dim,CELL>::PlaneStress( bool yes_no )
 { 
     plane_strain_ = yes_no; 
 }









/**
 
While Poisson's ratio must be an element variable, this method allows for
Young's modulus to be a node variable. Thus, a continuous loss of strength
can be modeled.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_B_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // only if the property is an element property  something is done here
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) {
         E_.resize(1);
         nu_.resize(1);
         // we made sure that youngs modulus is a scalar and has the same placement as Poisson's ratio
         E_[0]  = e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
         nu_[0] = e.Read( nu_key_ );
      }
    else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ) {
        const size_t ipoints(e.IntegrationPoints());
        E_.resize(ipoints);
        nu_.resize(ipoints);
        for ( uint32_t i{0}; i<ipoints; ++i )
          {
             // we made sure that youngs modulus is a scalar and has the same placement as Poisson's ratio
             E_[i]  = e.Read( i, MathOperatorLHS<dim,CELL>::MaterialOperandKey() );
             nu_[i] = e.Read( i, nu_key_ );
          }
       }
    else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE )
      {
        const size_t ipoints(e.IntegrationPoints());
        E_.resize(ipoints);
        nu_.resize(ipoints);
         for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ ) {
             E_[i] = e.PropertyValueAtIntegrationPoint( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), i );
             nu_[i] = e.PropertyValueAtIntegrationPoint( nu_key_, i );
          }
      }
     else throw csmp::Exception( WARNING, "NumIntegral_BT_D_B_dV<dim,CELL>::GetOperands",
                                 MathOperatorLHS<dim,CELL>::MaterialOperandName().c_str(),
                                "material operand placement not handled yet." );
      
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
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_BT_D_B_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // initialise output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( dim*e.Nodes(), dim*e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    // detect simplex elements (linear triangles/tetrahedra) where B is constant
    // and dN only needs to be evaluated once at the barycentre
    const bool is_simplex( e.FE()->IsSimplex() && e.Interpolation() == 1 );

    // build material property matrix D for element-constant properties
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      {
        if constexpr ( dim == 1U ) {
            // in 1D, the stiffness reduces to Young's modulus
            MathOperatorLHS<dim,CELL>::LHS.AssignToDiagonal( E_[0U] );
            return;
          }
        else if constexpr ( dim == 2U ) {
            if ( plane_strain_ ) planeStrainMatrix( E_[0U], nu_[0U], D );
            else                 planeStressMatrix( E_[0U], nu_[0U], D );
          }
        else stiffnessMatrix( E_[0U], nu_[0U], D );
      }

    // for simplex elements B is constant — evaluate dN once at the barycentre
    if ( is_simplex )
      {
        // build D at integration point 0 for spatially varying properties
        if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
             MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE )
          {
            if constexpr ( dim == 1U ) {
                MathOperatorLHS<dim,CELL>::LHS.AssignToDiagonal( E_[0U] );
                return;
              }
            else if constexpr ( dim == 2U ) {
                if ( plane_strain_ ) planeStrainMatrix( E_[0U], nu_[0U], D );
                else                 planeStressMatrix( E_[0U], nu_[0U], D );
              }
            else stiffnessMatrix( E_[0U], nu_[0U], D );
          }

        const double detJ = e.dN_AtBaryCenter( B, dim );
        if ( detJ <= 0.0 )
          {
            cerr <<"\n\tElement "<< e.Idx() <<": determinant of Jacobian at barycentre: "<< detJ << endl;
            throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_B_dV<dim>::ComputeContribution",
                                   "Jacobian transformation failed. Element nodes are perhaps not numbered correctly." );
          }

        B.Transposed( BT );
        BT *= D;
        BT *= B;
        BT *= e.Volume();

        MathOperatorLHS<dim,CELL>::LHS = BT;
        return;
      }

    // numerical integration over Gauss points for non-simplex isoparametric elements
    for ( uint32_t i{0U}; i < e.IntegrationPoints(); ++i )
      {
        // rebuild D at each integration point if properties vary spatially
        if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
             MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE )
          {
            if constexpr ( dim == 1U ) {
                MathOperatorLHS<dim,CELL>::LHS.AssignToDiagonal( E_[i] );
                return;
              }
            else if constexpr ( dim == 2U ) {
                if ( plane_strain_ ) planeStrainMatrix( E_[i], nu_[i], D );
                else                 planeStressMatrix( E_[i], nu_[i], D );
              }
            else stiffnessMatrix( E_[i], nu_[i], D );
          }

        const double detJ = e.dN_AtIntegrationPoint( B, i, dim );
        if ( detJ <= 0.0 )
          {
            cerr <<"\n\tElement "<< e.Idx() <<": determinant of Jacobian at Gauss point "<< i
                 <<": "<< detJ << endl;
            throw csmp::Exception( FATAL_ERROR, "NumIntegral_BT_D_B_dV<dim>::ComputeContribution",
                                   "Jacobian transformation failed. Element nodes are perhaps not numbered correctly." );
          }

        // BT(dim*nodes x 3) = B^T
        B.Transposed( BT );

        // BT(dim*nodes x 3) * D(3x3) * B(3 x dim*nodes) -> BT(dim*nodes x dim*nodes)
        BT *= D;
        BT *= B;
        BT *= e.WeightAtIntegrationPoint(i) * detJ;

        // accumulate Gauss point contribution into element stiffness matrix
        MathOperatorLHS<dim,CELL>::LHS += BT;
      }

 } // end ComputeContribution

template class NumIntegral_BT_D_B_dV<1U,Element>;
template class NumIntegral_BT_D_B_dV<2U,Element>;
template class NumIntegral_BT_D_B_dV<3U,Element>;

template class NumIntegral_BT_D_B_dV<1U,Face>;
template class NumIntegral_BT_D_B_dV<2U,Face>;
template class NumIntegral_BT_D_B_dV<3U,Face>;

} // csmp
