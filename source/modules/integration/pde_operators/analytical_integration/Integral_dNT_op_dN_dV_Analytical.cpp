#include "Integral_dNT_op_dN_dV_Analytical.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

	/** Accumulates the conductance matrix of the interpolation function
	derivatives.
	*/
	template<uint32_t dim, class CELL>
	Integral_dNT_op_dN_dV_Analytical<dim, CELL>::Integral_dNT_op_dN_dV_Analytical(const PropertyDatabase<dim>& pref,
		const char*  oper, const char*  basic, const char*  test)
		: MathOperatorLHS<dim>(pref, oper, basic, test)
	{
		MathOperatorLHS<dim>::Name("Integral_dNT_op_dN_dV_Analytical", oper, basic, test);

		// testing the Operands 
		if (MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() != REGION)
			throw csmp::Exception(ERROR, "Integral_dNT_op_dN_dV_Analytical::(constructor)",
				oper, "Operand must be placed on the element or group.");

		if (MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR)
			throw csmp::Exception(ERROR, "Integral_dNT_op_dN_dV_Analytical::(constructor)",
				test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

		if (MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR)
			throw csmp::Exception(ERROR, "Integral_dNT_op_dN_dV_Analytical::(constructor)",
				test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
	}

	template<uint32_t dim, class CELL>
	void Integral_dNT_op_dN_dV_Analytical<dim, CELL>::ComputeContribution( const CELL& e )
    { 
      MathOperatorLHS<dim>::LHS.Resize(e.Nodes(), e.Nodes());
      e.Integral_dNT_K_dN(MathOperatorLHS<dim>::LHS , MathOperatorLHS<dim>::MTRL[0]); 
    }


	template class Integral_dNT_op_dN_dV_Analytical<1U, Element<1U> >;
	template class Integral_dNT_op_dN_dV_Analytical<2U, Element<2U> >;
	template class Integral_dNT_op_dN_dV_Analytical<3U, Element<3U> >;

	template class Integral_dNT_op_dN_dV_Analytical<1U, Face<1U> >;
	template class Integral_dNT_op_dN_dV_Analytical<2U, Face<2U> >;
	template class Integral_dNT_op_dN_dV_Analytical<3U, Face<3U> >;

} // csmp
