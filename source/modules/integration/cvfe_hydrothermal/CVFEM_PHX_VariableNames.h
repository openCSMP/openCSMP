#ifndef CVFEM_PHX_VARIABLENAMES_H
#define CVFEM_PHX_VARIABLENAMES_H

#include "CSMP_definitions.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++.
*/

using namespace std;

namespace csmp {

class CVFEM_PHX_VariableNames {
  public:

    CVFEM_PHX_VariableNames( );
    
   std::vector<std::string> pres_grad_variables;
   std::vector<std::string> pore_visitor_variables;

   std::vector<std::string> lhs_liquid_vector, rhs_liquid_vector,
                            lhs_vapor_vector, rhs_vapor_vector;
                            
   std::vector<std::string> fv_transport_vapor_variables, fv_transport_liquid_variables;
   
   std::vector<std::string> upwind_control_variables;
   std::vector<std::string> relperm_visc, densities, saturations, velocities, pore_velocities, time_control;

   std::vector<std::string> capacitance_lhs_variables, conductance_variables,
                            capacitance_rhs_variables, heat_bottom_variables;

   std::vector<std::string> capacitance_lhs_p_variables,
                            conductance_p_upwind_liquid_variables, conductance_p_upwind_vapor_variables,
                            capacitance_rhs_p_variables,
                            grav_liq_variables, grav_vap_variables,
                            source_p_variables, source_p2_variables;

   std::vector<std::string> mass_visitor_variables,
                            enthalpy_visitor_variables,
                            enthalpy_visitor_liquid_variables,enthalpy_visitor_vapor_variables,
                            enthalpy_visitor_source_variables,
                            conduction_visitor_variables;

   std::vector<std::string> transient_variables, full_reset_variables;

   std::vector<std::string> diff_mass_variables, diff_enthalpy_variables;

   std::vector<std::string> boundary_corrections;

   void AddAdvectionVariable( std::string balanced_variable,
                              std::string new_lhs_liquid, std::string new_rhs_liquid,
                              std::string new_lhs_vapor,  std::string new_rhs_vapor  );

};

inline CVFEM_PHX_VariableNames::CVFEM_PHX_VariableNames( )
 {
  
  // variable names for contructor arguments
  // please see individual constructors and
  // the CVFEM_PHX_Scheme-Class
  
  // pressure gradient
  pres_grad_variables.push_back("fluid pressure");
  pres_grad_variables.push_back("reference pressure");
  pres_grad_variables.push_back("permeability");
  pres_grad_variables.push_back("KgradP");
  pres_grad_variables.push_back("gradP scaling");

  // pore volume
  pore_visitor_variables.push_back("nodal porosity");
  pore_visitor_variables.push_back("volume");
  pore_visitor_variables.push_back("pore volume");

  // mass advection
  lhs_liquid_vector.push_back("fluid mass liquid");
  rhs_liquid_vector.push_back("liquid mass mobility");
  lhs_vapor_vector.push_back("fluid mass vapor");
  rhs_vapor_vector.push_back("vapor mass mobility");

  // energy advection
  lhs_liquid_vector.push_back("enthalpy content liquid");
  rhs_liquid_vector.push_back("liquid enthalpy mobility");
  lhs_vapor_vector.push_back("enthalpy content vapor");
  rhs_vapor_vector.push_back("vapor enthalpy mobility");

  // salt advection
  lhs_liquid_vector.push_back("salt content liquid");
  rhs_liquid_vector.push_back("liquid salt mobility");
  lhs_vapor_vector.push_back("salt content vapor");
  rhs_vapor_vector.push_back("vapor salt mobility");

  // magmatic mass advection
  lhs_liquid_vector.push_back("magmatic fluid mass liquid");
  rhs_liquid_vector.push_back("magmatic liquid mass mobility");
  lhs_vapor_vector.push_back("magmatic fluid mass vapor");
  rhs_vapor_vector.push_back("magmatic vapor mass mobility");

  // magmatic salt advection
  lhs_liquid_vector.push_back("magmatic salt content liquid");
  rhs_liquid_vector.push_back("magmatic liquid salt mobility");
  lhs_vapor_vector.push_back("magmatic salt content vapor");
  rhs_vapor_vector.push_back("magmatic vapor salt mobility");

  // transports
  fv_transport_vapor_variables.push_back("KgradP");
  fv_transport_vapor_variables.push_back("nodal source vapor");
//  fv_transport_vapor_variables.push_back("density vapor");
  fv_transport_vapor_variables.push_back("density vapor transport");
  fv_transport_liquid_variables.push_back("KgradP");
  fv_transport_liquid_variables.push_back("nodal source liquid");
//  fv_transport_liquid_variables.push_back("density liquid");
  fv_transport_liquid_variables.push_back("density liquid transport");

  // upwind control
  upwind_control_variables.push_back("permeability");
  upwind_control_variables.push_back("porosity");
//  densities.push_back("density liquid");
//  densities.push_back("density vapor");
  densities.push_back("density liquid transport");
  densities.push_back("density vapor transport");
  relperm_visc.push_back("relperm viscosity liquid");
  relperm_visc.push_back("relperm viscosity vapor");
  saturations.push_back("saturation liquid");
  saturations.push_back("saturation vapor");
  time_control.push_back("courant liquid");
  time_control.push_back("courant vapor");
  velocities.push_back("velocity liquid");
  velocities.push_back("velocity vapor");
  pore_velocities.push_back("pore velocity liquid");
  pore_velocities.push_back("pore velocity vapor");

  // pde_operators heat diffusion
  capacitance_lhs_variables.push_back("nodal heat capacity");
  capacitance_lhs_variables.push_back("temperature");
  capacitance_lhs_variables.push_back("temperature");

  conductance_variables.push_back("thermal conductivity");
  conductance_variables.push_back("temperature");
  conductance_variables.push_back("temperature");

  capacitance_rhs_variables.push_back("nodal heat capacity");
  capacitance_rhs_variables.push_back("temperature");

  heat_bottom_variables.push_back("nodal heat flux bottom");
  heat_bottom_variables.push_back("temperature");

  // pde_operators pressure equation
  capacitance_lhs_p_variables.push_back("nodal total compressibility");
  capacitance_lhs_p_variables.push_back("fluid pressure");
  capacitance_lhs_p_variables.push_back("fluid pressure");

  conductance_p_upwind_liquid_variables.push_back("permeability");
  conductance_p_upwind_liquid_variables.push_back("fluid pressure");
  conductance_p_upwind_liquid_variables.push_back("fluid pressure");
  conductance_p_upwind_liquid_variables.push_back("liquid mass mobility");
//  conductance_p_upwind_liquid_variables.push_back("density liquid");
  conductance_p_upwind_liquid_variables.push_back("density liquid transport");

  conductance_p_upwind_vapor_variables.push_back("permeability");
  conductance_p_upwind_vapor_variables.push_back("fluid pressure");
  conductance_p_upwind_vapor_variables.push_back("fluid pressure");
  conductance_p_upwind_vapor_variables.push_back("vapor mass mobility");
//  conductance_p_upwind_vapor_variables.push_back("density vapor");
  conductance_p_upwind_vapor_variables.push_back("density vapor transport");

  capacitance_rhs_p_variables.push_back("nodal total compressibility");
  capacitance_rhs_p_variables.push_back("fluid pressure");

  grav_liq_variables.push_back("permeability");
  grav_liq_variables.push_back("fluid pressure");
  grav_liq_variables.push_back("liquid mass mobility density");
//  grav_liq_variables.push_back("density liquid");
  grav_liq_variables.push_back("density liquid transport");

  grav_vap_variables.push_back("permeability");
  grav_vap_variables.push_back("fluid pressure");
  grav_vap_variables.push_back("vapor mass mobility density");
//  grav_vap_variables.push_back("density vapor");
  grav_vap_variables.push_back("density vapor transport");

  source_p_variables.push_back("nodal fluid volume source");
  source_p_variables.push_back("fluid pressure");

  source_p2_variables.push_back("fluid source rate");
  source_p2_variables.push_back("fluid pressure");

  // CVFEM_Visitors
  // not needed at the moment
  // unless consistency check is done
  mass_visitor_variables.push_back("pore volume");
  mass_visitor_variables.push_back("mass");

  enthalpy_visitor_variables.push_back("volume");
  enthalpy_visitor_variables.push_back("enthalpy");
  
  enthalpy_visitor_liquid_variables.push_back("enthalpy liquid");
  enthalpy_visitor_vapor_variables.push_back("enthalpy vapor");
  enthalpy_visitor_source_variables.push_back("fluid source h");

  conduction_visitor_variables.push_back("volume");
  conduction_visitor_variables.push_back("delta h conduction");

  // transient variables
  transient_variables.push_back("temperature");
  transient_variables.push_back("fluid pressure");
  transient_variables.push_back("salinity");
  transient_variables.push_back("fluid mass liquid");
  transient_variables.push_back("fluid mass vapor");
  transient_variables.push_back("fluid density");
  transient_variables.push_back("enthalpy content liquid");
  transient_variables.push_back("enthalpy content vapor");
  transient_variables.push_back("salt content liquid");
  transient_variables.push_back("salt content vapor");
  transient_variables.push_back("mass");
  transient_variables.push_back("enthalpy");
  transient_variables.push_back("magmatic fluid mass");
  transient_variables.push_back("magmatic fluid mass liquid");
  transient_variables.push_back("magmatic fluid mass vapor");
  transient_variables.push_back("magmatic mass salt");
  transient_variables.push_back("magmatic salt content liquid");
  transient_variables.push_back("magmatic salt content vapor");
  transient_variables.push_back( "total fluid flux in" );
  transient_variables.push_back( "total fluid flux out" );
  transient_variables.push_back( "magmatic liquid flux in" );
  transient_variables.push_back( "magmatic liquid flux out" );
  transient_variables.push_back( "magmatic vapor flux in" );
  transient_variables.push_back( "magmatic vapor flux out" );
  transient_variables.push_back( "shell liquid flux in" );
  transient_variables.push_back( "shell liquid flux out" );
  transient_variables.push_back( "shell vapor flux in" );
  transient_variables.push_back( "shell vapor flux out" );

  // full reset variables
  full_reset_variables.push_back("liquid mass mobility");
  full_reset_variables.push_back("vapor mass mobility");
  full_reset_variables.push_back("liquid enthalpy mobility");
  full_reset_variables.push_back("vapor enthalpy mobility");
  full_reset_variables.push_back("liquid salt mobility");
  full_reset_variables.push_back("vapor salt mobility");
//  full_reset_variables.push_back("density liquid");
//  full_reset_variables.push_back("density vapor");
  full_reset_variables.push_back("density liquid transport");
  full_reset_variables.push_back("density vapor transport");
  full_reset_variables.push_back("relperm viscosity liquid");
  full_reset_variables.push_back("relperm viscosity vapor");
  full_reset_variables.push_back("liquid mass mobility density");
  full_reset_variables.push_back("vapor mass mobility density");
  full_reset_variables.push_back("enthalpy liquid");
  full_reset_variables.push_back("enthalpy vapor");
  full_reset_variables.push_back("nodal heat capacity");
  full_reset_variables.push_back("nodal total compressibility");
  full_reset_variables.push_back("fluid source rate");
  full_reset_variables.push_back("fluid source h");
  full_reset_variables.push_back("fluid source wt");
  full_reset_variables.push_back("previous total enthalpy");
  full_reset_variables.push_back("previous mass salt");
  full_reset_variables.push_back("previous fluid density");
  full_reset_variables.push_back("nodal fluid volume source");
  full_reset_variables.push_back("magmatic liquid mass mobility");
  full_reset_variables.push_back("magmatic vapor mass mobility");
  full_reset_variables.push_back("magmatic liquid salt mobility");
  full_reset_variables.push_back("magmatic vapor salt mobility");

  // diff variables
  diff_mass_variables.push_back("diff mass");
  diff_mass_variables.push_back("mass");
  diff_mass_variables.push_back("fluid density");

  diff_enthalpy_variables.push_back("diff enthalpy");
  diff_enthalpy_variables.push_back("enthalpy");
  diff_enthalpy_variables.push_back("previous total enthalpy");
  
  // boundary corrections
  boundary_corrections.push_back("boundary flow mass");
  boundary_corrections.push_back("boundary flow enthalpy");
  boundary_corrections.push_back("boundary flow salt");
  
 };

inline void CVFEM_PHX_VariableNames::AddAdvectionVariable( std::string balanced_variable,
                                                           std::string new_lhs_liquid, std::string new_rhs_liquid,
                                                           std::string new_lhs_vapor,  std::string new_rhs_vapor  )
 {

  transient_variables.push_back(balanced_variable);
  transient_variables.push_back(new_lhs_liquid);
  transient_variables.push_back(new_lhs_vapor);

  lhs_liquid_vector.push_back(new_lhs_liquid);
  rhs_liquid_vector.push_back(new_rhs_liquid);
  lhs_vapor_vector.push_back(new_lhs_vapor);
  rhs_vapor_vector.push_back(new_rhs_vapor);

 };

  /**
     @class CVFEM_PHX_VariableNames CVFEM_PHX_VariableNames.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Class to bundle all variable names needed for CVFEM scheme.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Applies source term at the node.
          
     @endcode
     
     @section dependencies Dependencies
	 Tailored for CVFEM_PHX_Schem
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */

} // end namespace csmp

#endif
