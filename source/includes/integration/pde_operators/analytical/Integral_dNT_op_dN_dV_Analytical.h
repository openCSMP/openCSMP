#ifndef Integral_dNT_op_dN_dV_ANALYTICAL_h
#define Integral_dNT_op_dN_dV_ANALYTICAL_h

#include "MathOperatorLHS.h"

namespace csmp {

template<uint32_t> class Element;

	/// Known as: element conductance matrix or K div^2 P
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_dNT_op_dN_dV_Analytical : public MathOperatorLHS<dim,CELL> {
	public:
		Integral_dNT_op_dN_dV_Analytical( const PropertyDatabase<dim>&,
			                                const char* oper, const char* basic, const char* test);

		void ComputeContribution( const CELL<dim>& )  override final;
    
		Integral_dNT_op_dN_dV_Analytical<dim, CELL>* clone() const override final { return new Integral_dNT_op_dN_dV_Analytical<dim, CELL>(*this); }
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
