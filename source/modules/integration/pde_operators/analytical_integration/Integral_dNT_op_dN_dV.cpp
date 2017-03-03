#include "Integral_dNT_op_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** Accumulates the conductance matrix of the interpolation function
derivatives.  
*/
template<size_t dim,class SIMPLEX>
Integral_dNT_op_dN_dV<dim,SIMPLEX>::Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                      const char*             oper, 
                                                      const char*             basic, 
                                                      const char*             test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    DN(2,3), DNT(3,2)
{
    MathOperatorLHS<dim>::Name("Integral_dNT_op_dN_dV", oper, basic, test );

    // testing the Operands 
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_dV::(constructor)", 
                    oper, "Operand must be placed on the element or group.");

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_dV::(constructor)", 
                    test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Integral_dNT_op_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}






/**
 
@section arguments Input Arguments 

A reference to the current Element.  
*/
template<size_t dim,class SIMPLEX>
void Integral_dNT_op_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );

    // calculate the element contribution to LHS (E_OP is the Basic Operand)
    DNT *= MathOperatorLHS<dim>::MTRL[0];
    DNT *= DN;
    // assigning element contribution
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS = DNT;
    // analytical integration over area / volume for linear triangle 
    // and tetrahedron elements, respectively
    MathOperatorLHS<dim>::LHS *= e.Volume();

} // end ComputeContribution



template class Integral_dNT_op_dN_dV<1U,Element<1U> >;
template class Integral_dNT_op_dN_dV<2U,Element<2U> >;
template class Integral_dNT_op_dN_dV<3U,Element<3U> >;

template class Integral_dNT_op_dN_dV<1U,Face<1U> >;
template class Integral_dNT_op_dN_dV<2U,Face<2U> >;
template class Integral_dNT_op_dN_dV<3U,Face<3U> >;

} // csmp
