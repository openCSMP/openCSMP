#include "Integral_dNT_op_dN_dV_Analytical.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

	/** Accumulates the conductance matrix of the interpolation function
	derivatives.
	*/
template<uint32_t dim, template<uint32_t> class CELL>
Integral_dNT_op_dN_dV_Analytical<dim, CELL>::Integral_dNT_op_dN_dV_Analytical( const PropertyDatabase<dim>& pref,
		                                                                           const char*  oper,
                                                                               const char*  basic,
                                                                               const char*  test )
  : MathOperatorLHS<dim,CELL>(pref, oper, basic, test)
	{
		MathOperatorLHS<dim,CELL>::Name("Integral_dNT_op_dN_dV_Analytical", oper, basic, test);

		// testing the Operands 
		if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != REGION )
			throw csmp::Exception(ERROR, "Integral_dNT_op_dN_dV_Analytical::(constructor)",
				oper, "Operand must be placed on the element or group.");

		if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR)
			throw csmp::Exception(ERROR, "Integral_dNT_op_dN_dV_Analytical::(constructor)",
				test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

		if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
			throw csmp::Exception(ERROR, "Integral_dNT_op_dN_dV_Analytical::(constructor)",
				test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
	}





template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_op_dN_dV_Analytical<dim, CELL>::ComputeContribution( const CELL<dim>& e )
  {
    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(), e.Nodes());
    DenseMatrix<DM3> TMP(dim,dim);
    for ( uint32_t i{0}; i<dim; ++i )
      for ( uint32_t j{0}; j<dim; ++j )
        TMP(i,j) = MathOperatorLHS<dim,CELL>::MTRL[0](i,j);
    e.Integral_dNT_K_dN( MathOperatorLHS<dim,CELL>::LHS, TMP );
  }


template class Integral_dNT_op_dN_dV_Analytical<1U, Element>;
template class Integral_dNT_op_dN_dV_Analytical<2U, Element>;
template class Integral_dNT_op_dN_dV_Analytical<3U, Element>;

template class Integral_dNT_op_dN_dV_Analytical<1U, Face>;
template class Integral_dNT_op_dN_dV_Analytical<2U, Face>;
template class Integral_dNT_op_dN_dV_Analytical<3U, Face>;

} // csmp
