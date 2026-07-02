#ifndef CVFEM_1D_VVCASE
#define CVFEM_1D_VVCASE

#include "CSMP_definitions.h"
#include "Test.h"
#include "ComputationalSettings.h"
#include "ScalarVariable.h"
#include "PDE_Integrator.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"

using namespace std;
using namespace csmp;


namespace csmp  {

template<uint32_t> class Element;
template<uint32_t> class Model;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class VTU_Interface;
template<uint32_t> class CVFEM_PHX_Scheme;

  /** CVFEM_1D_VVCase : 1D Test Case
  =================================
  Mesh:       Not needed
  Model:      1D
  Test:       PHX_Scheme
  BC:         TBC
  Criterion:
  =================================
  */
class CVFEM_1D_VVCase : public Test
  {
    public:
          CVFEM_1D_VVCase ( const char* prefix );
          ~CVFEM_1D_VVCase();
          void run();
    private:

          const char* geometry_name_;
          const char* config_file_name_;
          const char* vars_name_;
          const char* system_file_list_name_;
          const char* recipes_file_list_name_;
          list<string> output_props_;

          Model<1U>*             model = nullptr;
          PropertyDatabase<1U>*  pd_ref = nullptr;
          VTU_Interface<1U>*     vtu = nullptr;
          ComputationalSettings  run_settings;

          //! CVFEM instance
          CVFEM_PHX_Scheme<1U>* CVFEM_PHX = nullptr;

          //! Chemical diffusion operators and variables
          std::vector<PDE_Integrator<1U>*>             diffusion;
          std::vector<NumIntegral_NT_lhsop_N_dV<1U>*>  diff_cap_lhs;
          std::vector<NumIntegral_NT_op_N_dV<1U>*>     diff_cap_rhs;
          std::vector<NumIntegral_dNT_op_dN_dV<1U>*>   diffusivity;
          Index pore_diff_key, eff_diff_key, porosity_key;
          ScalarVariable diff_coeff, porosity;

          //! time variables
          double model_time, dt, max_time, vtu_dt;
          size_t timestep, output_counter;

          int CalculateMassConductivity();
          int CalculateMassGravityTerm();
          int CalculateSteadyStatePressure();
          int CalculateSteadyStateTemperature();
          int InitializeFluidPropertiesLinearPressure();
          int InitializeFluidProperties();

  };

  } // csmp

#endif // CVFEM_1D_VVCase
