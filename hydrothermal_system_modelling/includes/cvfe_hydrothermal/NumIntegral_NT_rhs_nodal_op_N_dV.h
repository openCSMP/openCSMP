// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUMINTEGRAL_NT_RHS_NODAL_OP_N_DV_H
#define NUMINTEGRAL_NT_RHS_NODAL_OP_N_DV_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

/*   Changelog
     November 2014, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

//template<uint32_t> class Element; // JK: commented out

/**
   @class NumIntegral_NT_rhs_nodal_op_N_dV NumIntegral_NT_rhs_nodal_op_N_dV.h

   Capacitance matrix (RHS) with the operand values applied strictly to the node.

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


template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_rhs_nodal_op_N_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_NT_rhs_nodal_op_N_dV( const PropertyDatabase<dim>& pref,
                                      const char* oper,
                                      const char* test,
                                      const char* thickness );//Benoit 2025 add

    ~NumIntegral_NT_rhs_nodal_op_N_dV();

    virtual void GetOperands( const CELL<dim>& );
    virtual void ComputeContribution( const CELL<dim>& );
    virtual NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL>* clone() const { return new NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL> (*this); }
  private:

    NumIntegral_NT_rhs_nodal_op_N_dV();
    std::vector<ScalarVariable> VAR;
    double volume;
    csmp::Index thickness_;//Benoit 2025 add

    
};

} // csmp

#endif
