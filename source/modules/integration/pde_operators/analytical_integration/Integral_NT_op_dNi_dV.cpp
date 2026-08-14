#include "Integral_NT_op_dNi_dV.h"
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
Integral_NT_op_dNi_dV<dim,CELL>::Integral_NT_op_dNi_dV( const PropertyDatabase<dim>& pref,
                                                        const char* oper, 
                                                        const char* mtrl, 
                                                        const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    prop_key(pref.StorageKey(mtrl)),
    OP(3),
    IPOL(3),
    DN(dim,3),
    coord(dim),
    gravity(ACC_GRAVITY), // scalar acts to increase the pressure
    prop2_time_multiplier(1.0),
    xyz(2)
 {
    MathOperatorRHS<dim,CELL>::Name("Integral_NT_op_dNi_dV", oper, test );
    
        // testing the Operands 
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_op_dNi_dV::(constructor)", 
                    oper,      "Operand must be a scalar property." );

    
    if ( prop_key.place != ELEMENT and prop_key.place != REGION )
    throw csmp::Exception( ERROR, "Integral_NT_op_dNi_dV::(constructor)", 
                   mtrl, "Material operand must be placed on the element or group." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_NT_op_dNi_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }




template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_op_dNi_dV<dim,CELL>::SpatialDerivative( uint32_t num_xyz )
 {
    assert( num_xyz > 0 && num_xyz <=3 );
    xyz = num_xyz;
 }






template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_op_dNi_dV<dim,CELL>::MaterialPropertyTimeMultiplier( double time_increment )
 {
    prop2_time_multiplier = time_increment;
 }



/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_op_dNi_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
   // reading fluid density or something like that
   if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT or
        MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == REGION )
     e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), prop1 );
   else // property is interpolated to element center
     {
        e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), OP );
        e.N_AtBaryCenter( IPOL );
        prop1() = 0.0;
        for ( auto i{0U}; i<e.Nodes(); i++ ) prop1 += IPOL[i] * OP[i];
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
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_op_dNi_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    e.dN( DN );
    
    double vol = e.Volume();

    for ( auto i{0U}; i<e.Nodes(); i++ ) 
      //                             gradZ        density    K          acc.gravity    element volume
      MathOperatorRHS<dim,CELL>::RHS[i] = DN(xyz-1,i) * prop1() * prop2() * -gravity  * vol;

// cout <<"\nIntegral_NT_op_dNi_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

} // end ComputeContribution


template class Integral_NT_op_dNi_dV<1U,Element>;
template class Integral_NT_op_dNi_dV<2U,Element>;
template class Integral_NT_op_dNi_dV<3U,Element>;

template class Integral_NT_op_dNi_dV<1U,Face>;
template class Integral_NT_op_dNi_dV<2U,Face>;
template class Integral_NT_op_dNi_dV<3U,Face>;

} // csmp
