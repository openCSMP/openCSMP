#include "Integral_dNT_op_dN_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** Accumulates the conductance matrix of the interpolation function
    derivatives.
*/
template<uint32_t dim, template<uint32_t> class CELL>
Integral_dNT_op_dN_dV<dim,CELL>::Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                        const char*             oper,
                                                        const char*             basic,
                                                        const char*             test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    DN(2,3), DNT(3,2)
{
    MathOperatorLHS<dim,CELL>::Name("Integral_dNT_op_dN_dV", oper, basic, test );

    // testing the Operands 
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_dNT_op_dN_dV::(constructor)", 
                           oper, "Operand must be placed on the element or group.");

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_op_dN_dV::(constructor)", 
                           test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_op_dN_dV::(constructor)", 
                           test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}






/**
 
@section arguments Input Arguments 

A reference to the current Element.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_op_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );

    // calculate the element contribution to LHS (E_OP is the Basic Operand)
    DNT *= MathOperatorLHS<dim,CELL>::MTRL[0];
    DNT *= DN;
    // assigning element contribution
    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim,CELL>::LHS = DNT;
    // analytical integration over area / volume for linear triangle 
    // and tetrahedron elements, respectively
    MathOperatorLHS<dim,CELL>::LHS *= e.Volume();

} // end ComputeContribution



template class Integral_dNT_op_dN_dV<1U,Element>;
template class Integral_dNT_op_dN_dV<2U,Element>;
template class Integral_dNT_op_dN_dV<3U,Element>;

template class Integral_dNT_op_dN_dV<1U,Face>;
template class Integral_dNT_op_dN_dV<2U,Face>;
template class Integral_dNT_op_dN_dV<3U,Face>;

} // csmp
