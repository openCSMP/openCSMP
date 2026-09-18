// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NumIntegral_PT_op_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
      To consider body forces
      Gravity vector must be a vector property (involving g and density of rock)
 */
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_PT_op_dV<dim,CELL>::NumIntegral_PT_op_dV( const PropertyDatabase<dim>& pref,
                                                      const char* oper, const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    BFORCE(6) 
{
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_PT_op_dV", oper, test );
    
    // testing the specified Operands 
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_dV<dim>::(constructor", 
                             oper,  "Operand must be a vector variable.");

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_dV<dim>::(constructor", 
                             test,  "Test function Operand must be a vector variable placed on the nodes");
}





/** The nodal forces are read from each element and stored into a vector
which is of dimensions nodes x dimensions. 

If the option AdaptToQuadraticTriangle() is set, the forces
are distributed only on the midside nodes as this proves to
be the best procedure available (see Cook et al.).  

@section arguments Input Arguments 

A reference to the variable storage inside of the Model<dim>  and a 
reference to the Element from which the Operand shall be read.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_PT_op_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

   BFORCE.resize( e.Nodes() * dim );

   // reading "body" forces from the Nodes
   if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE )
     {
        vector<VectorVariable<dim> >  forces;
        e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), forces );
   
        // remapping the forces into the vector E_OP
        size_t k(0);
        for ( uint32_t i{0U}; i<e.Nodes(); i++ )
          for ( uint32_t j{0U}; j<dim; j++ ) BFORCE[k++] = forces[i][j];
     }
   else // Element property
     {
        const PLACEMENT place = MathOperatorRHS<dim,CELL>::MaterialOperandKey().place;
        assert( place == ELEMENT || place == FACE );
        VectorVariable<dim>  vc;
        e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), vc );
   
        // remapping the forces into the vector E_OP
        size_t k(0);
        for ( uint32_t i{0U}; i<e.Nodes(); i++ )
          for ( uint32_t j{0U}; j<dim; j++ ) BFORCE[k++] = vc[j];
    }

   // zeroing the corner nodes again (only for quadratic triangle in 2D
   // where number of midside nodes is equal to corner nodes)
   if ( e.FE_Type() == QUADRATIC_TRIANGLE ||
        e.FE_Type() == ISOPARAMETRIC_QUADRATIC_TRIANGLE )
     for ( uint32_t i{0U}; i<e.Nodes(); i++ ) BFORCE[i] = 0.;

} // end GetOperands









/** Computes the forces which act on the finite element due to the body
force source term and the area / volume of the element.  

The resulting contribution V to the righthand vector is of the form:  

{F} = (volume/nodes) { fx1, fy1, fz1, .... } 

The forces are assigned only to the midside nodes ! - Thus for the 
QuadraticTriangle, the force vector has the form:  

{F} = F * (volume/midside-nodes) { 0 0 0 1 1 1 } 

@section arguments Input Arguments 

A reference to element, the contribution of which is to be aquired and
the time-increment over which the deformation shall occur. 

The result of the computation is returned into the base class protected
member vector {V}.   

@section application Application

In linear elasticity computations.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_PT_op_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    MathOperatorRHS<dim,CELL>::RHS.resize( BFORCE.size() );
    // if a quadratic triangle element is used, the body forces are assigned only 
    // to the midside nodes
   if ( e.FE_Type() == QUADRATIC_TRIANGLE ||
        e.FE_Type() == ISOPARAMETRIC_QUADRATIC_TRIANGLE )
      {
         // watch out, this is not generic but restricted to 6-noded triangle in 2D
         const double volume_div_n(e.Volume() / static_cast<double>(e.FE()->MidSideNodes()));
         for ( uint32_t i{0U}; i<MathOperatorRHS<dim,CELL>::RHS.size(); i++ )
           MathOperatorRHS<dim,CELL>::RHS[i] = BFORCE[i] * volume_div_n;
     }
   else
     {
         const double volume_div_n(e.Volume() / static_cast<double>(e.Nodes()));
         for ( uint32_t i{0U}; i<MathOperatorRHS<dim,CELL>::RHS.size(); i++ ) 
           MathOperatorRHS<dim,CELL>::RHS[i] = BFORCE[i] * volume_div_n;
     }
    
} // end ComputeContribution



template class NumIntegral_PT_op_dV<1U,Element>;
template class NumIntegral_PT_op_dV<2U,Element>;
template class NumIntegral_PT_op_dV<3U,Element>;

template class NumIntegral_PT_op_dV<1U,Face>;
template class NumIntegral_PT_op_dV<2U,Face>;
template class NumIntegral_PT_op_dV<3U,Face>;

} // csmp


