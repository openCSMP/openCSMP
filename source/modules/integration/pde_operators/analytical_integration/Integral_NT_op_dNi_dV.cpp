#include "Integral_NT_op_dNi_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** The operand defines the fluid density, and the mtrl variable would for
instance be the hydraulic conductivity.  
*/
template<size_t dim,class SIMPLEX>
Integral_NT_op_dNi_dV<dim,SIMPLEX>::Integral_NT_op_dNi_dV( const PropertyDatabase<dim>& pref,
                                                      const char* oper, 
                                                      const char* mtrl, 
                                                      const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    prop_key(pref.StorageKey(mtrl)),
    OP(3),
    IPOL(3),
    DN(dim,3),
    coord(dim),
    gravity(9.80665), // scalar acts to increase the pressure
    prop2_time_multiplier(1.0),
    xyz(2)
 {
    MathOperatorRHS<dim>::Name("Integral_NT_op_dNi_dV", oper, test );
    
        // testing the Operands 
    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Integral_NT_op_dNi_dV::(constructor)", 
                    oper,      "Operand must be a scalar property." );

    
    if ( prop_key.place != ELEMENT and prop_key.place != REGION )
    throw csmp::Exception( CSMP_ERROR, "Integral_NT_op_dNi_dV::(constructor)", 
                   mtrl, "Material operand must be placed on the element or group." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Integral_NT_op_dNi_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }




template<size_t dim,class SIMPLEX>
void Integral_NT_op_dNi_dV<dim,SIMPLEX>::SpatialDerivative( size_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }


template<size_t dim,class SIMPLEX>
void Integral_NT_op_dNi_dV<dim,SIMPLEX>::MaterialPropertyTimeMultiplier( double64 time_increment )
 {
    prop2_time_multiplier = time_increment;
 }



/** Reads the Operand values from the elements.
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_op_dNi_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
   // reading fluid density or something like that
   if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT or MathOperatorRHS<dim>::MaterialOperandPlacement() == REGION ) 
     e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), prop1 );
   else // property is interpolated to element center
     {
        e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), OP );
        e.N_AtBaryCenter( IPOL );
        prop1() = 0.0;
        for ( size_t i=0; i<e.Nodes(); i++ ) prop1 += IPOL[i] * OP[i];
     }
     
   // reading the material property variable (for instance conductivity)
   e.Read( prop_key, prop2 );
   
   // multiplying operand with time multiplier
   prop2() *= prop2_time_multiplier;
   
//   cout <<"\ngravity: "<< gravity << endl;    
//   cout << op.Name() <<": "<< prop1() << endl;
//   cout <<"Total mobility: "<< prop2() << endl;
}




/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  

@section implementation Implementation

gravity is negative since it acts in the opposite direction of the 
coordinate axis.
*/
template<size_t dim,class SIMPLEX>
void Integral_NT_op_dNi_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    e.dN( DN );
    
    double64 vol = e.Volume();

    for ( size_t i=0; i<e.Nodes(); i++ ) 
      //                             gradZ        density    K          acc.gravity    element volume
      MathOperatorRHS<dim>::RHS[i] = DN(xyz-1,i) * prop1() * prop2() * -gravity  * vol;

// cout <<"\nIntegral_NT_op_dNi_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

} // end ComputeContribution


template class Integral_NT_op_dNi_dV<1U,Element<1U> >;
template class Integral_NT_op_dNi_dV<2U,Element<2U> >;
template class Integral_NT_op_dNi_dV<3U,Element<3U> >;

template class Integral_NT_op_dNi_dV<1U,Face<1U> >;
template class Integral_NT_op_dNi_dV<2U,Face<2U> >;
template class Integral_NT_op_dNi_dV<3U,Face<3U> >;

} // csmp
