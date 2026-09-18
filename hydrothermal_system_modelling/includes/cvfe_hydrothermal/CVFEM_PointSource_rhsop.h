// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_POINTSOURCE_RHSOP_H
#define CVFEM_POINTSOURCE_RHSOP_H

#include "CSMP_definitions.h"
#include "CVFEM_MathOperatorRHS.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<uint32_t> class Element;

  /**
     @class CVFEM_PointSource_rhsop CVFEM_PointSource_rhsop.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes
  
     @section motivation Motivation
      Equal to PointSource_rhsop, but can be applied as CVFEM_Visitor.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Applies source term at the node.
          
     @endcode
     
     @section dependencies Dependencies
	 CVFEM_MathOperatorRHS
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
template<uint32_t dim, template<uint32_t> class CELL=Element>
class CVFEM_PointSource_rhsop : public CVFEM_MathOperatorRHS<dim> {
  public:
    CVFEM_PointSource_rhsop( const PropertyDatabase<dim>& pref,
                              const char* nodal_src,
                              const char* test);

    ~CVFEM_PointSource_rhsop();
    
    virtual void GetOperands( const CELL<dim>& );
    virtual void GetOperandsCVFEM( const CELL<dim>&, csmp::Index upwind_var_key );
    virtual void ComputeContribution( const CELL<dim>& );
    
    virtual CVFEM_PointSource_rhsop<dim,CELL>* clone() const { return new CVFEM_PointSource_rhsop<dim,CELL>(*this); }
    
  private:
  
    std::vector<ScalarVariable> SRC_;
    std::vector<ScalarVariable> upwind_var_;

};


} // csmp

#endif
















