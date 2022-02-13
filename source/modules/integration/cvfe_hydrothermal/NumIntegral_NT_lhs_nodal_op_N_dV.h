#ifndef NUMINTEGRAL_NT_LHS_NODAL_OP_N_DV_H
#define NUMINTEGRAL_NT_LHS_NODAL_OP_N_DV_H

/*   Changelog
     November 2014, Philipp Weis:
	 - initial port to CSMP++.
*/

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

namespace csmp {

/// Capacitance matrix (LHS) with the operand values applied strictly to the node

template<uint32_t dim, class CELL>
class NumIntegral_NT_lhs_nodal_op_N_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_NT_lhs_nodal_op_N_dV( const PropertyDatabase<dim>& p, 
                                      const char* oper, const char* basic, const char* test );
    
    ~NumIntegral_NT_lhs_nodal_op_N_dV();

    virtual void GetOperands( CELL& e );
    virtual void ComputeContribution( CELL& e );
    virtual NumIntegral_NT_lhs_nodal_op_N_dV<dim,CELL >* clone() const { return new NumIntegral_NT_lhs_nodal_op_N_dV<dim,Element<dim> >(*this); }
  private:

    NumIntegral_NT_lhs_nodal_op_N_dV();    
    std::vector<ScalarVariable> VAR;    
    double volume;

};

  /**
     @class NumIntegral_NT_lhs_nodal_op_N_dV NumIntegral_NT_lhs_nodal_op_N_dV.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Calculations with large point source terms and values in the capacitance matrix that can vary over several orders of magnitude require that the material operand has to be evaluated at the node.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 The material operand is read in as a nodal vector and applied to the entries of that respective node.
	 No averaging at the constraint points is performed.
          
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
} // csmp

#endif
















