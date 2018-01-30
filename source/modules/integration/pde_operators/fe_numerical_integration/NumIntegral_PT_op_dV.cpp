#include "NumIntegral_PT_op_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
NumIntegral_PT_op_dV<dim,SIMPLEX>::NumIntegral_PT_op_dV( const PropertyDatabase<dim>& pref,
                                                    const char* oper, const char* test ) 
  : MathOperatorRHS<dim>(pref,oper,test),
    BFORCE(6) 
{
    MathOperatorRHS<dim>::Name("NumIntegral_PT_op_dV", oper, test );
    
    // testing the specified Operands 
    if ( MathOperatorRHS<dim>::MaterialOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_op_dV<dim>::(constructor", 
                             oper,  "Operand must be a vector variable.");

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != VECTOR )
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
template<size_t dim,class SIMPLEX>
void NumIntegral_PT_op_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.Isoparametric() == true );

   BFORCE.resize( e.Nodes() * dim );

   // reading the nodal "body" forces
   if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE )
     {
        vector<VectorVariable<dim> >  forces;
        e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), forces );
   
        // remapping the forces into the vector E_OP
        size_t k(0);
        for ( size_t i=0; i<e.Nodes(); i++ )
          for ( size_t j=0; j<dim; j++ ) BFORCE[k++] = forces[i][j];
     }
   else // Element property
     {
        assert( MathOperatorRHS<dim>::MaterialOperandKey().place == ELEMENT );
        VectorVariable<dim>  vc;
        e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), vc );
   
        // remapping the forces into the vector E_OP
        size_t k(0);
        for ( size_t i=0; i<e.Nodes(); i++ )
          for ( size_t j=0; j<dim; j++ ) BFORCE[k++] = vc[j];
    }

   // zeroing the corner nodes again (only for quadratic triangle in 2D
   // where number of midside nodes is equal to corner nodes)
   if ( e.FE_Type() == QUADRATIC_TRIANGLE ||
        e.FE_Type() == ISOPARAMETRIC_QUADRATIC_TRIANGLE )
     for ( size_t i=0; i<e.Nodes(); i++ ) BFORCE[i] = 0.;

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
template<size_t dim,class SIMPLEX>
void NumIntegral_PT_op_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    MathOperatorRHS<dim>::RHS.resize( BFORCE.size() );
    // if a quadratic triangle element is used, the body forces are assigned only 
    // to the midside nodes
   if ( e.FE_Type() == QUADRATIC_TRIANGLE ||
        e.FE_Type() == ISOPARAMETRIC_QUADRATIC_TRIANGLE )
      {
         // watch out, this is not generic but restricted to 6-noded triangle in 2D
         const double64 volume_div_n(e.Volume() / static_cast<double64>(e.FE()->MidSideNodes()));
         for ( size_t i=0; i<MathOperatorRHS<dim>::RHS.size(); i++ ) 
           MathOperatorRHS<dim>::RHS[i] = BFORCE[i] * volume_div_n;
     }
   else
     {
         const double64 volume_div_n(e.Volume() / static_cast<double64>(e.Nodes()));
         for ( size_t i=0; i<MathOperatorRHS<dim>::RHS.size(); i++ ) 
           MathOperatorRHS<dim>::RHS[i] = BFORCE[i] * volume_div_n;
     }
    
} // end ComputeContribution



template class NumIntegral_PT_op_dV<1U,Element<1U> >;
template class NumIntegral_PT_op_dV<2U,Element<2U> >;
template class NumIntegral_PT_op_dV<3U,Element<3U> >;

template class NumIntegral_PT_op_dV<1U,Face<1U> >;
template class NumIntegral_PT_op_dV<2U,Face<2U> >;
template class NumIntegral_PT_op_dV<3U,Face<3U> >;

} // csmp


