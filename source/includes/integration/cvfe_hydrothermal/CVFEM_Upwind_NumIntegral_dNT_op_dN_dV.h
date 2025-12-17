#ifndef CVFEM_UPWIND_NUMINTEGRAL_DNT_OP_DN_DV_H
#define CVFEM_UPWIND_NUMINTEGRAL_DNT_OP_DN_DV_H

#include "CVFEM_MathOperatorLHS.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class UpwindControlVisitor;
template<uint32_t> class ExplicitFiniteVolumeTransportPHX;

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
template<uint32_t dim, template<uint32_t> class CELL=Element>
class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV : public CVFEM_MathOperatorLHS<dim,CELL> {
  public:
    CVFEM_Upwind_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>&, 
                                           UpwindControlVisitor<dim>& upwind_visitor,
                                           ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                           const char* oper,
                                           const char* basic,
                                           const char* test,
                                           const char* upwind,
                                           const char* grav_trigger );
    
    void ComputeContribution( const CELL<dim>& ) override final;
    void GetOperands( const CELL<dim>& ) override final;
    void GetOperandsCVFEM( const CELL<dim>&, csmp::Index upwind_var_key ) override final;

    void GetUpwindMatrix( const CELL<dim>& );
    
    CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>* clone() const override final { return new CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>(*this); }

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


} // csmp

#endif
