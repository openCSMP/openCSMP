#ifndef CVFEM_UPWIND_NUMINTEGRAL_DNT_OP_DN_DV_H
#define CVFEM_UPWIND_NUMINTEGRAL_DNT_OP_DN_DV_H

#include "CSMP_definitions.h"
#include "CVFEM_MathOperatorLHS.h"
#include "UpwindControlVisitor.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<uint32_t dim>
class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV : public CVFEM_MathOperatorLHS<dim> {
  public:
    CVFEM_Upwind_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                                            UpwindControlVisitor<dim>& upwind_visitor,
                                            ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                          const char* oper,
                                          const char* basic,
                                          const char* test,
                                          const char* upwind,
                                          const char* grav_trigger);

    ~CVFEM_Upwind_NumIntegral_dNT_op_dN_dV();
    
    void ComputeContribution( Element<dim>& e );
    void GetOperands( Element<dim>& e );
    void GetOperandsCVFEM( Element<dim>& e, csmp::Index upwind_var_key );

    void GetUpwindMatrix( Element<dim>& e );
    virtual CVFEM_Upwind_NumIntegral_dNT_op_dN_dV <dim>* clone() const { return new CVFEM_Upwind_NumIntegral_dNT_op_dN_dV <dim>(*this); }

  private:

    UpwindControlVisitor<dim>& UpwindVisitor;
    ExplicitFiniteVolumeTransportPHX<dim>& finite_volume;
    
    DenseMatrix<DM_MIN>  B, upwind;

    double area, normal_component;
    double contribution, operand;
    uint32_t inside_node_, outside_node_;
    
    std::vector<ScalarVariable> el_uvar;
    std::vector<ScalarVariable> upwind_var_multiplier;
    
    csmp::Index uvar_;
    csmp::Index gtvar_;
    
};

  /**
     @class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Conductance matrix for fluid pressure calculations using pre-defined upwind nodes.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Calculates conductance matrix with pre-defined upwind nodes.
	 All fluid properties are used at the upstream nodes.
          
     @endcode
     
     @section dependencies Dependencies
	 CVFEM_MathOperatorLHS
	 UpwindControlVisitor
	 ExplicitFiniteVolumeTransportPHX
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */

} // csmp

#endif
