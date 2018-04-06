#ifndef Integral_dNT_op_dN_dV_ANALYTICAL_h
#define Integral_dNT_op_dN_dV_ANALYTICAL_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {

	/// Known as: element conductance matrix or K div^2 P
	template<size_t dim, class SIMPLEX = Element<dim> >
	class Integral_dNT_op_dN_dV_Analytical : public MathOperatorLHS<dim> {
	public:
		Integral_dNT_op_dN_dV_Analytical(const PropertyDatabase<dim>& pref,
			const char* oper, const char* basic, const char* test);

		virtual void ComputeContribution(SIMPLEX& e);
		virtual Integral_dNT_op_dN_dV_Analytical<dim, SIMPLEX>* clone() const { return new Integral_dNT_op_dN_dV_Analytical<dim, SIMPLEX>(*this); }
	};


	/**

	@class Integral_dNT_op_dN_dV Integral_dNT_op_dN_dV "pde_operators/Integral_dNT_op_dN_dV.h"
	@author S.K. Matthaei
	@author S. Roberts
	@date 1999

	@section motivation Motivation

	PDE operator representing the divergence squared of the dependent
	variable.*/

} // csmp

#endif
#pragma once
