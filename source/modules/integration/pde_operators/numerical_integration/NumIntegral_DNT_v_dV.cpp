#include "NumIntegral_DNT_v_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

/**
 
The operand defines the fluid density, and the mtrl variable would for 
instance be the hydraulic conductivity.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_DNT_v_dV<dim,CELL>::NumIntegral_DNT_v_dV( const PropertyDatabase<dim>& pref,
                                                      const char* oper,       // e.g., Darcy velocity
                                                      const char* r_factor,
                                                      const char* dens,       // e.g., fluid density
                                                      const char* test )      // streaming potential
                          
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    rfac_key(pref.StorageKey(r_factor)),
    rho_key(pref.StorageKey(dens)),
    IPOL(3),
    DN(dim,3),
    DNT(3,dim),
    OPMAT(dim,1),
    oper_nprop(3), 
    dens_nprop(3)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_DNT_v_dV", oper, test );
    
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT || 
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != VECTOR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_v_dV<dim>::(constructor)", 
                   oper, "Dependent variable must be a vector property placed on the elements." );

    if ( rfac_key.type != SCALAR || rfac_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_v_dV<dim>::(constructor)", 
                    r_factor,  " must be a scalar property on the element." );

    if ( rho_key.place != NODE || rho_key.type != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_v_dV<dim>::(constructor)", 
                   dens, " operand must be a scalar placed on the nodes." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_v_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
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
void NumIntegral_DNT_v_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

   if ( MathOperatorRHS<dim,CELL>::MultiplyWithTimeIncrement() ) {
        throw csmp::Exception( FATAL_ERROR, "NumIntegral_DNT_v_dV<dim>::GetOperands", 
           "Do not multiply this operator with time increment since it uses the acceleration of gravity");
        throw invalid_argument("NumIntegral_DNT_v_dV<dim>::GetOperands");
     }


   // 0. reading the r-factor (i.e. some material property which is constant 
   //    on the element
   e.Read( rfac_key, rfac );

   // 1. reading Operand (velocity)
   if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) 
     e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), oper_eprop );
   else 
     e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), oper_nprop );
   
   OPMAT(0,0) = oper_eprop[0];
   if ( dim != 1U ) OPMAT(1,0) = oper_eprop[1];
   if ( dim == 3U ) OPMAT(2,0) = oper_eprop[2];
     
   // 2. reading fluid density variable placed on the node
   e.NodePropertyVector( rho_key, dens_nprop );

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
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNT_v_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    double detJ, fdensity;
    
    for ( uint32_t i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
         // interpolation functions to interpolate density and velocity
         e.N_AtIntegrationPoint( i, IPOL );
         detJ = e.dN_AtIntegrationPoint( DN, i );
         DN.Transposed( DNT );
         fdensity = 0.;
         
         if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) {
               // interpolating density multiplier to integration point
               for ( auto j{0U}; j<e.Nodes(); j++ ) 
                 fdensity += IPOL[j] * dens_nprop[j]();

               // multiply property value it with multipliers
               OPMAT(0,0) *= fdensity;
               if ( dim != 1U ) OPMAT(1,0) *= fdensity;
               if ( dim == 3U ) OPMAT(2,0) *= fdensity;
           }
         else // ifOperand is placed on the NODE
           {
               // interpolating Operand value and nodal multipliers to integration point
           }

         // matrix multiplication
         DNT *= OPMAT;

         // assembling contribution to right-hand vector
         for ( uint32_t j{0U}; j<e.Nodes(); j++ )
           MathOperatorRHS<dim,CELL>::RHS[j] = rfac() * DNT(j,0) * e.WeightAtIntegrationPoint(i) * detJ;
      }

// cout <<"\nNumIntegral_DNT_v_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

} // end ComputeContribution


template class NumIntegral_DNT_v_dV<1U,Element>;
template class NumIntegral_DNT_v_dV<2U,Element>;
template class NumIntegral_DNT_v_dV<3U,Element>;

template class NumIntegral_DNT_v_dV<1U,Face>;
template class NumIntegral_DNT_v_dV<2U,Face>;
template class NumIntegral_DNT_v_dV<3U,Face>;

} // csmp
