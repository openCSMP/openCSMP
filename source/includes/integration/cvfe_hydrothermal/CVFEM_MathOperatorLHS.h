#ifndef CVFEM_MATHOPERATOR_LHS_H
#define CVFEM_MATHOPERATOR_LHS_H

#include "MathOperatorLHS.h"

/*   Changelog
     November 2014, Philipp Weis:
	 - initial port to CSMP++.
     Oct 2022, SKM
*/

namespace csmp {

template<uint32_t> class Element;

/**
     @class CVFEM_MathOperatorLHS CVFEM_MathOperatorLHS.h
     
     Base class for CVFEM Math Operators (LHS) that enable the use for control volume calculations as a post-processing step.

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes
  
     @section motivation Motivation
      In the development for the CVFEM scheme, we tested the consistency of the pressure equation and mass advection.

     @section usage Usage
      Can be used within the CVFEM scheme (Weis et al., Geofluids, 2014).
	  PDE Operators that use these Math Operators as base classes can be used with the CVFEM_Visitor.


     @code
	 Usage and functionality as the MathOperatorLHS with the addition of access functions that are used by the CVFEM_Visitor.
          
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues
	 This construct is redundant if CVFEM_Visitors are not used.
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
template<uint32_t dim, template<uint32_t> class CELL=Element>
class CVFEM_MathOperatorLHS : public MathOperatorLHS<dim,CELL> {
  public:
  
    CVFEM_MathOperatorLHS( const PropertyDatabase<dim>&,
                           const char* basic,
                           const char* test );

    CVFEM_MathOperatorLHS( const PropertyDatabase<dim>&, 
                           const char* oper, const char* basic, const char* test );
    
    virtual ~CVFEM_MathOperatorLHS();

    virtual void GetOperandsCVFEM( const CELL<dim>&, csmp::Index upwind_var_key );
    
    virtual DenseMatrix<DM_MIN> GetContribution( );
    
    virtual CVFEM_MathOperatorLHS<dim,CELL>* clone() const override { return new CVFEM_MathOperatorLHS<dim,CELL>(*this); }

  private:
    CVFEM_MathOperatorLHS();
};

} // csmp

#endif
















