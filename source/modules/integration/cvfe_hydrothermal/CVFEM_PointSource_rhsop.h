#ifndef CVFEM_POINTSOURCE_RHSOP_H
#define CVFEM_POINTSOURCE_RHSOP_H

#include "CSMP_definitions.h"
#include "CVFEM_MathOperatorRHS.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<uint32_t dim,class CELL>
class CVFEM_PointSource_rhsop : public CVFEM_MathOperatorRHS<dim> {
  public:
    CVFEM_PointSource_rhsop( const PropertyDatabase<dim>& p, const char* nodal_src, const char* test );

    ~CVFEM_PointSource_rhsop();
    
    virtual void GetOperands( CELL& e );
    virtual void GetOperandsCVFEM( CELL& e, csmp::Index upwind_var_key );
    virtual void ComputeContribution( CELL& e );
    virtual CVFEM_PointSource_rhsop<dim,CELL>* clone() const { return new CVFEM_PointSource_rhsop<dim,CELL> (*this); }
  private:
  
    std::vector<ScalarVariable> SRC_;
    std::vector<ScalarVariable> upwind_var_;

};

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

} // csmp

#endif
















