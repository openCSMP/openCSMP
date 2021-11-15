#ifndef CVFEM_UPWIND_NUMIINTEGRAL_DNT_RHSOP_G_DV_H
#define CVFEM_UPWIND_NUMIINTEGRAL_DNT_RHSOP_G_DV_H

#include "CSMP_definitions.h"
#include "CVFEM_MathOperatorRHS.h"
#include "UpwindControlVisitor.h"
#include "Operand.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<size_t dim>
class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV : public CVFEM_MathOperatorRHS<dim> {
  public:
    CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV(const PropertyDatabase<dim>& p, 
                                   UpwindControlVisitor<dim>& upwind_visitor,
                                   ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                   const char* oper,
                                   const char* test,
                                   const char* upwind,
                                   const char* grav_trigger);
    
    ~CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV();
    
    virtual void GetOperands( Element<dim>& e );
    virtual void GetOperandsCVFEM( Element<dim>& e, csmp::Index upwind_var_key );
    virtual void ComputeContribution( Element<dim>& e );
    void GetUpwindMatrix( Element<dim>& e );
    virtual CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>* clone() const { return new CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim> (*this); }

  private:

    UpwindControlVisitor<dim>& UpwindVisitor;
    ExplicitFiniteVolumeTransportPHX<dim>& finite_volume;

    DenseMatrix<DM_MIN>  B, upwind;

    double area, normal_component;
    double contribution, operand, grav;
    size_t inside_node_, outside_node_;
        
    double                          gravity;// acceleration of gravity
    size_t                   xyz;// 1=x, 2=y, 3=z
    
    csmp::Index rho_key;
    
    Operand<dim> upwind_;
    std::vector<ScalarVariable> upwind_var_, upwind_var_multiplier;
     
};

  /**
     @class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Gravity component for fluid pressure calculations using pre-defined upwind nodes.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Calculates gravity component with pre-defined upwind nodes.
	 All fluid properties are used at the upstream nodes.
          
     @endcode
     
     @section dependencies Dependencies
	 CVFEM_MathOperatorRHS
	 UpwindControlVisitor
	 ExplicitFiniteVolumeTransportPHX
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */

} // csmp

#endif
















