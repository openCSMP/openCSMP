#ifndef THERMAL_VISITOR_H
#define THERMAL_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#include "H2OLookup.h"

namespace csmp {

struct Index;
template<uint32_t> class Model;

/**

@brief Computes PT-dependent fluid properties for flow simulations. 

@author Thomas Driesner, ETH Zuerich
@section contact Contact
thomas.driesner@erdw.ethz.ch

*/
template<uint32_t dim>
class ThermalVisitor : public Visitor<dim>
{
  public:
      ThermalVisitor(Model<dim>& model , std::vector<std::string> *to_initialize_keys =NULL);

      ThermalVisitor( Model<dim>& model,
                      Index nCPT_key,
                      Index cp_fluid_key,
                      Index cp_rock_key,
                      Index rho_rock_key,
                      Index t_key,
                      Index p_key,
                      Index hl_key,
                      Index hCl_key,
                      Index hClp_key,
                      Index hVl_key,
                      Index mt_key,
                      Index rl_key,
                      Index mu_key,
                      Index nphi_key,
                      Index nQ_key,
                      Index dRho_key,
                      Index volume_factor_key,
                      Index beta_rock_key,
                      Index beta_tot_key,
                      Index beta_fluid_key);
          
      ~ThermalVisitor();
      void SetInitialProperties(Model<dim> *model);
	  void ApplyTemperatureBoundaryConditionsToTransportedVariables(); 
      virtual void Visit(Node<dim>* n);  
      virtual void Visit(Model<dim>* n);  

 private:
      ThermalVisitor();

      //! variables
      
      Model<dim>& model_;

      std::vector<std::string>* to_initialize_keys_;
      
      ScalarVariable  nCPT,
                      cp_fluid,
                      cp_rock,
                      rho_rock,
                      t,
                      p,
                      hl,
                      hCl,
                      hClp,
                      hVl,
                      mt,
                      rl,
                      volume_factor,
                      mu,
                      phi,
                      nQ,
                      dRho,
                      beta_rock,
					            beta_fluid,
                      beta_tot;

      csmp::Index     nCPT_key_,
                      cp_fluid_key_,
                      cp_rock_key_,
                      rho_rock_key_,
                      t_key_,
                      p_key_,
                      hl_key_,
                      hCl_key_,
                      hClp_key_,
                      hVl_key_,
                      mt_key_,
                      rl_key_,
                      mu_key_,
                      volume_factor_key_,
                      nphi_key_,
                      nQ_key_,
                      dRho_key_,
                      beta_rock_key_,
                      beta_fluid_key_,
                      beta_tot_key_;

      double        pore_volume,
                      rock_volume,
                      dhCl_,
                      dT; 
      
      H2OLookup       water;

};

}// csmp

#endif //ThermalVisitor_h

