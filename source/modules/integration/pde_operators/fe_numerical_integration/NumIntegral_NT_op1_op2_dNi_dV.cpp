#include "NumIntegral_NT_op1_op2_dNi_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


/**
 
The operand defines the fluid density, and the mtrl variable would for
instance be the hydraulic conductivity.  
*/
template<size_t dim,class SIMPLEX>
NumIntegral_NT_op1_op2_dNi_dV<dim,SIMPLEX>::NumIntegral_NT_op1_op2_dNi_dV( const PropertyDatabase<dim>& pref,
                                                                          const char* oper, 
                                                                          const char* mtrl1, 
                                                                          const char* mtrl2, 
                                                                          const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    IPOL(3),
    oper_nprop(3), 
    gravity(9.80665), // acceleration of gravity
    xyz( (dim==1u) ? 0u : 1u ), // X in 1D, else Y-direction
    mtrl1_key(pref.StorageKey(mtrl1)),
    mtrl2_key(pref.StorageKey(mtrl2))
 {
    MathOperatorRHS<dim>::Name("NumIntegral_NT_op1_op2_dNi_dV", oper, test );
    
    
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_NT_op1_op2_dNi_dV<dim>::(constructor)", 
                    oper,      "Operand must be a scalar property." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op1_op2_dNi_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );

    if ( mtrl1_key.place != ELEMENT )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op1_op2_dNi_dV<dim>::(constructor)", 
                      mtrl1, "Material operand must be placed on the element." );

    if ( mtrl2_key.place != ELEMENT )
      throw csmp::Exception( ERROR, "NumIntegral_NT_op1_op2_dNi_dV<dim>::(constructor)", 
                      mtrl2, "Material operand must be placed on the element." );
 }




// true is the default
//template<size_t dim,class SIMPLEX>
//void NumIntegral_NT_op1_op2_dNi_dV<dim>::MultiplyMaterialPropertyWithTime( int mtrl_prop )
//  {
//     multiply_mtrl1_with_time_increment = (mtrl_prop==1) ? true : false;
//  }



/**
 
Reads the Operand from either the element or the nodes. Reads the material
property from the element, and its multipliers from the nodes for later
interpolation to the integration points.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_NT_op1_op2_dNi_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{

   // 1. reading Operand (fluid density or something like that)
   if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
     e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), oper_eprop );
   else 
     e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), oper_nprop );
     
   // 2. reading material property variables
   e.Read( mtrl1_key, mtrl1_prop );
   e.Read( mtrl2_key, mtrl2_prop );

} // end GetOperands






/** Computes, for example the body force term ensuing from the density field
within the element.   

This involves numerical integration over the element volume of the 
product:
 
(mtrl1 / dt + mtrl2) * volume integral (NT nmult NT op) dNi
 

Here, the op (Operand) is the fluid density and 'nmult' is the nodal
property which is used as material multiplier.  

@section arguments Input Arguments 

A reference to the finite-element from which the contribution is 
computed.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_NT_op1_op2_dNi_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
    
    // get the vertical density gradient 
    // if fluid density is an element property
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT )
      for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) 
        {
           // get dN = interpolation function derivate value at integration point
           double64 ip_value =  oper_eprop.Value(); // fluid density
           ip_value   *= -gravity;
           ip_value   *=  e.WeightAtIntegrationPoint(i);
           ip_value   *=  e.dN_AtIntegrationPoint( DN, i ); // det_J

           for ( size_t j=0; j<e.Nodes(); j++ )
             MathOperatorRHS<dim>::RHS[j] += ip_value * DN(xyz,j);
        }

    // if fluid density is a node property
    else
      for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) 
        {
           e.N_AtIntegrationPoint( i, IPOL );
           const double64 detJ = e.dN_AtIntegrationPoint( DN, i );
         
           // interpolating Operand value to integration point
           double64  ip_value(IPOL[0] * oper_nprop[0].Value());
           for ( size_t j=1; j<e.Nodes(); j++ )
             ip_value  += IPOL[j] * oper_nprop[j].Value();

           // assembling contribution to right-hand vector
           ip_value *= -gravity;
           ip_value *=  e.WeightAtIntegrationPoint(i);
           ip_value *=  detJ;
           
           for ( size_t j=0; j<e.Nodes(); j++ ) 
             MathOperatorRHS<dim>::RHS[j] += ip_value * DN(xyz,j);
      }

} // end ComputeContribution





// now everything is assembled (using the  hydraulic diffusivity kappa = k / (S mu)
template<size_t dim,class SIMPLEX>
void NumIntegral_NT_op1_op2_dNi_dV<dim,SIMPLEX>::MultiplyWithTimeFactor( double64 dt )
 {
    for ( typename vector<double64>::iterator 
          it=MathOperatorRHS<dim>::RHS.begin(); it!=MathOperatorRHS<dim>::RHS.end(); it++ ) {
          // not bad (no overshoot) 
         (*it) *= mtrl2_prop.Value();
         (*it) -= mtrl1_prop.Value() * dt;
      }

 } // end 


// cout <<"\nNumIntegral_NT_op1_op2_dNi_dV<dim>::MultiplyWithTimeFactor: ";
// out( MathOperatorRHS<dim>::RHS );
// cout << endl;


template class NumIntegral_NT_op1_op2_dNi_dV<1U,Element<1U> >;
template class NumIntegral_NT_op1_op2_dNi_dV<2U,Element<2U> >;
template class NumIntegral_NT_op1_op2_dNi_dV<3U,Element<3U> >;

template class NumIntegral_NT_op1_op2_dNi_dV<1U,Face<1U> >;
template class NumIntegral_NT_op1_op2_dNi_dV<2U,Face<2U> >;
template class NumIntegral_NT_op1_op2_dNi_dV<3U,Face<3U> >;

} // csmp
