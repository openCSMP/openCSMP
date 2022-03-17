#ifndef CVFEM_MATHOPERATORLHS_H
#define CVFEM_MATHOPERATORLHS_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"

/*   Changelog
     November 2014, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

/// Base class for CVFEM Math Operators (LHS) that enable the use for control volume calculations as a post-processing step.


template<uint32_t dim>
class CVFEM_MathOperatorLHS : public MathOperatorLHS<dim> {
  public:
  
    CVFEM_MathOperatorLHS( const PropertyDatabase<dim>& p, 
                     const char* basic, 
                     const char* test );

    CVFEM_MathOperatorLHS( const PropertyDatabase<dim>& p, 
                            const char* oper, const char* basic, const char* test );
    
    ~CVFEM_MathOperatorLHS();

    virtual void GetOperandsCVFEM( Element<dim>& e, csmp::Index upwind_var_key );
    virtual DenseMatrix<DM_MIN> GetContribution( );
    virtual CVFEM_MathOperatorLHS <dim>* clone() const { return new CVFEM_MathOperatorLHS <dim>(*this); }

  private:
    CVFEM_MathOperatorLHS();
    
};

  /**
     @class CVFEM_MathOperatorLHS CVFEM_MathOperatorLHS.h

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
} // csmp

#endif
















