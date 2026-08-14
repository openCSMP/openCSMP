#include "NumIntegral_NT_mixed_op_dNi_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {



/** The operand defines the fluid density, and the mtrl variable would for
instance be the hydraulic conductivity.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_mixed_op_dNi_dV<dim,CELL>::NumIntegral_NT_mixed_op_dNi_dV( const PropertyDatabase<dim>& pref,
                                                                          const char* nodal_mtrl_multiplier,
                                                                          const char* oper,
                                                                          const char* mtrl,
                                                                          const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    mtrl_key(pref.StorageKey(mtrl)),
    nmult_key(pref.StorageKey(nodal_mtrl_multiplier)),
    IPOL(3),
    DNI(3),
    oper_nprop(3), 
    mtrl_nprop(3),
    gravity(ACC_GRAVITY), // scalar acts to increase the pressure
    mtrl_time_multiplier(1.0),
    xyz(2) // Y-direction
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_NT_mixed_op_dNi_dV", oper, test );
    
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_mixed_op_dNi_dV<dim>::(constructor)", 
                    oper,      "Operand must be a scalar property." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_mixed_op_dNi_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );

    if ( mtrl_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_NT_mixed_op_dNi_dV<dim>::(constructor)", 
                   mtrl, "Material operand must be placed on the element." );

    if ( nmult_key.place != NODE || nmult_key.type != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_mixed_op_dNi_dV<dim>::(constructor)", 
                   nodal_mtrl_multiplier, "Nodal multiplier must be a scalar placed on the nodes." );
 }






template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_mixed_op_dNi_dV<dim,CELL>::SpatialDerivative( uint32_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }





template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_mixed_op_dNi_dV<dim,CELL>::MaterialPropertyTimeMultiplier( double time_increment )
 {
    mtrl_time_multiplier = time_increment;
 }



/**
 
Reads the Operand from either the element or the nodes. Reads the material
property from the element, and its multipliers from the nodes for later
interpolation to the integration points.  

If the MathOperator is used in a transient calculation, the function
MaterialPropertyTimeMultiplier() allows to multiply the material
parameter with the time-icrement to achieve symmetry with the lefthand
side of the equation. This is done, for instance, if the MathOperatorRHS
is used to compute the hydrostatic pressure contribution to a transient
system. Do not use MultiplyWithTimeIncrement() in this case, since
the acceleration of gravity must not be multiplied with delta t.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_mixed_op_dNi_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

   if ( MathOperatorRHS<dim,CELL>::MultiplyWithTimeIncrement() ) {
        throw csmp::Exception( FATAL_ERROR, "NumIntegral_NT_mixed_op_dNi_dV<dim>::GetOperands", 
           "Do not multiply this operator with time increment since it uses the acceleration of gravity");
     }

   // 1. reading Operand (fluid density or something like that)
   if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) 
     e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), oper_eprop );
   else 
     e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), oper_nprop );
     
   // 2. reading material property variable
   e.Read( mtrl_key, eprop );
   // multiplying material with time multiplier (the default of which is 1.0)
   eprop *= mtrl_time_multiplier;

   // 3. reading nodal material property multiplier
   e.NodePropertyVector( nmult_key, mtrl_nprop );

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
computed.  */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_mixed_op_dNi_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    double    ip_value, op_value, nmult_fac, detJ;
    uint32_t  j;
    
    DNI.resize(e.Nodes());
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );
    
    for ( uint32_t i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
         e.N_AtIntegrationPoint( i, IPOL );
         detJ = e.dN_AtIntegrationPoint( DN, i );
         
         if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) {
               // interpolating nodal multipliers to integration point
               for ( nmult_fac=0.0, j=0; j<e.Nodes(); j++ )
                 nmult_fac += IPOL[j] * mtrl_nprop[j]();

               // multiply property value it with multipliers
               ip_value = oper_eprop() * eprop() * nmult_fac * -gravity;
           }
         else // ifOperand is placed on the NODE
           {
               // interpolating Operand value and nodal multipliers to integration point
               for ( op_value=nmult_fac=0.0, j=0; j<e.Nodes(); j++ ) {
                    nmult_fac += IPOL[j] * mtrl_nprop[j]();
                    op_value  += IPOL[j] * oper_nprop[j]();
                 }
               // multiply property value it with multipliers
               ip_value = op_value * eprop() * nmult_fac * -gravity;
           }

         // assembling contribution to right-hand vector
         ip_value *= e.WeightAtIntegrationPoint(i) * detJ;
         //                                                2-1=Y-derivative
         for ( j=0U; j<e.Nodes(); j++ ) DNI[j] = ip_value * DN((xyz-1),j);
           
         for ( j=0U; j<e.Nodes(); j++ ) MathOperatorRHS<dim,CELL>::RHS[j] += DNI[j];
      }

// cout <<"\nNumIntegral_NT_mixed_op_dNi_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

} // end ComputeContribution


template class NumIntegral_NT_mixed_op_dNi_dV<1U,Element>;
template class NumIntegral_NT_mixed_op_dNi_dV<2U,Element>;
template class NumIntegral_NT_mixed_op_dNi_dV<3U,Element>;

template class NumIntegral_NT_mixed_op_dNi_dV<1U,Face>;
template class NumIntegral_NT_mixed_op_dNi_dV<2U,Face>;
template class NumIntegral_NT_mixed_op_dNi_dV<3U,Face>;

} // end csmp
