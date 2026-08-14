#ifndef INCOMPRESSIBLE_TWOPHASEFLOWFRACTURES_VISCOUS_VVCASE
#define INCOMPRESSIBLE_TWOPHASEFLOWFRACTURES_VISCOUS_VVCASE

#include "Test.h"
#include "Index.h"
#include "PDE_Integrator.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;
template<uint32_t> class NodeCenteredFiniteVolumeTransport;
template<uint32_t> class TwoPhaseModel;

#ifdef CSMP_WITH_SAMG_SOLVER
class SAMG_Settings;
class SAMG_Solver;
#endif

template<uint32_t dim>
class IncompressibleTwoPhaseFlowFractures_Viscous_VVCase : public Test {
  public:
    
    IncompressibleTwoPhaseFlowFractures_Viscous_VVCase(	const char* prefix,
                                        const char* twophase_model = "Brooks Corey",
                                        const char* explicit_implicit = "explicit",
                                        const char* first_second_order = "1st order");


    ~IncompressibleTwoPhaseFlowFractures_Viscous_VVCase();
	
    virtual void run();

  private:

    // Debug Mode
    bool debug_;

    //*****************************************************************
    // Model Construction
    //*****************************************************************

    // Mesh File
    std::string input_file_name_;
    std::string variables_file_prefix_;
    std::string variables_file_cutted_prefix_;

    // Load Model (different dimensions)
    Model<dim>* model_;
    void LoadModel();

    //*****************************************************************
    // Variables setup
    //*****************************************************************

    // Test variables

    // model boundaries
    std::string outlet_boundary_name_,inlet_boundary_name_;

    // fluid and rock properties
    std::string porosity_,rho_w_,rho_n_,visc_n_,visc_w_;
    std::string total_mobility_,permeability_,entry_pressure_,brooks_corey_parameter_, frac_apperture_;

    // flow properties
    std::string fluid_pressure_, cap_diffusivity_, sw_,swr_,sn_,snr_, nodal_fluid_volume_source_, fluid_volume_source_;
    std::string sw_limiter,sn_limiter, mass_center;
    std::string total_velocity_,velocity_n_,velocity_w_, total_volume_flux_, total_pore_velocity_, thickness_;

    // Variable indexes
    Index porosity_idx_,rho_w_idx_,rho_n_idx_;
    Index total_mobility_idx_, permeability_idx_;
    Index fluid_pressure_idx_,cap_diffusivity_idx_,sw_idx_,sn_idx_, fluid_volume_source_idx_,nodal_fluid_volume_source_idx_;
    Index total_velocity_idx_,velocity_w_idx_,velocity_n_idx_,thickness_idx_;

	void VariablesSetup();
	void ModelSetup();

    //*****************************************************************
    // Simulation setup
    //*****************************************************************

    // Simulation Data
    double time_increment_, max_time_;
    double tolerance_l2_,tolerance_linf_,tolerance_mean_;
    double cfl_multiplier_,courant_increment_;

    // RelPerm Models
    size_t relperm_model_type_;
    std::string relperm_model_name_;
    TwoPhaseModel<dim>* relperm_model_;
    enum{LINEAR_TWOPHASE_MODEL, BROOKS_COREY_MODEL, COREY_MODEL, VAN_GENUCHTEN_MODEL};

    // Simulation Type
    bool second_order_in_space_;
    bool second_order_in_time_;
    bool with_lsm_grad_limiter_;
    bool with_capillary_forces_;
    bool with_gravity_forces_;
    bool nonlinear_scheme_;
    bool explicit_scheme_;
    std::string no_flow_bc_;

    bool assign_dirichlet_saturation_on_the_right_boundaries_;
    bool lu_solver_;

    //*****************************************************************
    // Simulation tools
    //*****************************************************************

    // Simulation Objects
    NodeCenteredFiniteVolumeTransport<dim>*  tpncfvt_;
    
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings* steady_state_pressure_solver_settings_ = nullptr;
    SAMG_Solver*   solver_ = nullptr;
#endif
    PDE_Integrator<dim,Element>* steady_state_pressure_solver_ = nullptr;

    // Velocity Computations
    void UpdateSaturations(TwoPhaseModel<dim>& saturationFunctions );
    void ComputeTotalMobility(TwoPhaseModel<dim>& saturationFunctions );
  };

}


#endif // INCOMPRESSIBLETWOPHASEFLOWFRACTURES_VISCOUS_VVCASE
