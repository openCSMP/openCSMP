// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_UPWIND_NUMINTEGRAL_DNT_RHSOP_G_DV_H
#define CVFEM_UPWIND_NUMINTEGRAL_DNT_RHSOP_G_DV_H

#include "CSMP_definitions.h"
#include "CVFEM_MathOperatorRHS.h"
#include "UpwindControlVisitor.h"
#include "Operand.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

namespace csmp {

template<uint32_t> class Element; //JK: needed?
template<uint32_t> class UpwindControlVisitor; //JK: needed?
template<uint32_t> class ExplicitFiniteVolumeTransportPHX; //JK: needed?


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
template<uint32_t dim, template<uint32_t> class CELL=Element>
class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV : public CVFEM_MathOperatorRHS<dim> {
  public:
    CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV(const PropertyDatabase<dim>& p, 
                                   UpwindControlVisitor<dim>& upwind_visitor,
                                   ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                   const char* oper,
                                   const char* test,
                                   const char* upwind,
                                   const char* grav_trigger,
                                   const char* thickness);//Benoit 2025 add
    
    ~CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV();
    
    virtual void GetOperands( const CELL<dim>& e );
    virtual void GetOperandsCVFEM( const CELL<dim>& e, csmp::Index upwind_var_key );
    virtual void ComputeContribution( const CELL<dim>& e );
    void GetUpwindMatrix( const CELL<dim>& e );
    virtual CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>* clone() const { return new CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim> (*this); }

  private:

    UpwindControlVisitor<dim>& UpwindVisitor;
    ExplicitFiniteVolumeTransportPHX<dim>& finite_volume;

    DenseMatrix<DM_MIN>  upwind, B;

    double area, normal_component;
    double contribution, operand, grav;
    uint32_t inside_node_, outside_node_;
        
    double gravity;// acceleration of gravity
    int    xyz;// 1=x, 2=y, 3=z
    
    csmp::Index rho_key;
    csmp::Index thickness_;//Benoit 2025 add
    
    Operand<dim> upwind_;
    std::vector<ScalarVariable> upwind_var_, upwind_var_multiplier;
     
};


} // csmp

#endif
















