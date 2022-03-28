#include "NumIntegral_NT_op_dNi_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    @attention the acceleration of gravity varies considerably around our planet;
use the local value
*/

template<uint32_t dim,class CELL>
NumIntegral_NT_op_dNi_dV<dim,CELL>::NumIntegral_NT_op_dNi_dV( const PropertyDatabase<dim>& pref,
                                                                 const char* oper, 
                                                                 const char* test,
                                                                 double acc_gravity )
  : MathOperatorRHS<dim>(pref,oper,test),
    IPOL(3),
    oper_nprop(3), 
    gravity(acc_gravity),
    xyz( (dim==1u) ? 0 : 1 ) // X in 1D, else Y-direction
 {
    MathOperatorRHS<dim>::Name("NumIntegral_NT_op_dNi_dV", oper, test );
    
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_op_dNi_dV<dim>::(constructor)", 
                           oper, "Operand must be a scalar property." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_op_dNi_dV<dim>::(constructor)", 
                           test, "Dependent variable must be a scalar property placed on the nodes." );
 }





/** The operand defines the fluid density, and the mtrl variable would for
instance be the hydraulic conductivity.  

@attention the acceleration of gravity varies considerably around our planet;
use the local value

*/
template<uint32_t dim,class CELL>
NumIntegral_NT_op_dNi_dV<dim,CELL>::NumIntegral_NT_op_dNi_dV( const PropertyDatabase<dim>& pref,
                                                         const char* oper, 
                                                         const char* mtrl, 
                                                         const char* test,
                                                         double acc_gravity )
  : MathOperatorRHS<dim>(pref,oper,test),    
    mtrl_key(pref.StorageKey(mtrl)),
    IPOL(3),
    oper_nprop(3), 
    gravity(acc_gravity),
    xyz( (dim==1u) ? 0 : 1 ) // X in 1D, else Y-direction
 {
    MathOperatorRHS<dim>::Name("NumIntegral_NT_op_dNi_dV", oper, test );
    
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_op_dNi_dV<dim>::(constructor)", 
                           oper, "Operand must be a scalar property." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_op_dNi_dV<dim>::(constructor)", 
                           test, "Dependent variable must be a scalar property placed on the nodes." );

    if ( mtrl_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_NT_op_dNi_dV<dim>::(constructor)", 
                           mtrl, "Material operand must be placed on the element." );
 }


template<uint32_t dim,class CELL>
void NumIntegral_NT_op_dNi_dV<dim,CELL>::SpatialDerivative( uint32_t num_xyz )
 {
    assert( num_xyz >= 0 && num_xyz <3U );
    xyz = num_xyz;
 }




/** Reads the Operand from either the element or the nodes. Reads the material
property from the element, and its multipliers from the nodes for later
interpolation to the integration points.  
*/
template<uint32_t dim,class CELL>
void NumIntegral_NT_op_dNi_dV<dim,CELL>::GetOperands( const CELL& e )
{
   // 1. reading Operand (fluid density or something like that)
   if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
     oper_eprop = e.Read( MathOperatorRHS<dim>::MaterialOperandKey() );
   else if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
     e.IntegrationPointPropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), oper_nprop );
   else 
     e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), oper_nprop );
     
   // 2. reading material property variable
   if ( mtrl_key.IsDefined() ) eprop = e.Read( mtrl_key );
   else eprop = 1.;

} // end GetOperands






/** Computes, for example the body force term ensuing from the density field
within the element.   

This involves numerical integration over the element volume of the 
product:
 
mtrl * volume integral (NT nmult NT op) dNi
 
Here, the op (Operand) is the fluid density and 'nmult' is the nodal
property which is used as material multiplier.  

@section arguments Input Arguments 

A reference to the finite-element from which the contribution is 
computed.  
*/
template<uint32_t dim,class CELL>
void NumIntegral_NT_op_dNi_dV<dim,CELL>::ComputeContribution( const CELL& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    double  ip_value;
    
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
    
    for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
         if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) {
               // multiply property value it with multipliers
               ip_value = oper_eprop * eprop * -gravity;
           }
         else if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
           {
               // multiply property value it with multipliers
               ip_value = oper_nprop[i]() * eprop * -gravity;
           }
         else // if Operand is placed on the NODE
           {
               // interpolating Operand value and multipliers to integration point
               e.N_AtIntegrationPoint( i, IPOL );
               double op_value(0.);
               for ( auto j{0U}; j<e.Nodes(); j++ ) op_value += IPOL[j] * oper_nprop[j]();

               // multiply property value it with multipliers
               ip_value = op_value * eprop * -gravity;
           }

         // assembling contribution to right-hand vector
         ip_value *= e.WeightAtIntegrationPoint(i) 
                   * e.dN_AtIntegrationPoint( DN, i );
         //                                                                              Y-derivative
         for ( auto j{0U}; j<e.Nodes(); j++ ) MathOperatorRHS<dim>::RHS[j] += ip_value * DN((xyz),j);
      }

} // end ComputeContribution


template class NumIntegral_NT_op_dNi_dV<1U,Element<1U> >;
template class NumIntegral_NT_op_dNi_dV<2U,Element<2U> >;
template class NumIntegral_NT_op_dNi_dV<3U,Element<3U> >;

template class NumIntegral_NT_op_dNi_dV<1U,Face<1U> >;
template class NumIntegral_NT_op_dNi_dV<2U,Face<2U> >;
template class NumIntegral_NT_op_dNi_dV<3U,Face<3U> >;

} // csmp
