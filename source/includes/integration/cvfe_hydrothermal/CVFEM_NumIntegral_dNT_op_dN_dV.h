#ifndef CVFEM_NUM_INTEGRAL_DNT_OP_DN_DV_H
#define CVFEM_NUM_INTEGRAL_DNT_OP_DN_DV_H

#include "CVFEM_MathOperatorLHS.h"

/**   
     Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<uint32_t> class Element;

template<uint32_t dim, template<uint32_t> class CELL=Element>
class CVFEM_NumIntegral_dNT_op_dN_dV : public CVFEM_MathOperatorLHS<dim,CELL> {
  public:
    CVFEM_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>&,
                                    const char* oper,
                                    const char* basic,
                                    const char* test );
     
    void ComputeContribution( const CELL<dim>& ) override final;
    
    CVFEM_NumIntegral_dNT_op_dN_dV<dim,CELL>* clone() const override final { return new CVFEM_NumIntegral_dNT_op_dN_dV<dim,CELL>(*this); }

private:
    DenseMatrix<DM_MIN>  B, BT;
};

  /**
     @class CVFEM_NumIntegral_dNT_op_dN_dV CVFEM_NumIntegral_dNT_op_dN_dV.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Equal to NumIntegral_dNT_op_dN_dV, but can be applied as CVFEM_Visitor.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Calculates conductance matrix.
          
     @endcode
     
     @section dependencies Dependencies
	 CVFEM_MathOperatorLHS
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */

} // csmp

#endif
