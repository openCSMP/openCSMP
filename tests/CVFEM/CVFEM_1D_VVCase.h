#ifndef CVFEM_1D_VVCASE_H
#define CVFEM_1D_VVCASE_H


#include <chrono>
#include "iostream"

// model
#include "Model1D.h"
#include "ModelTime.h"
#include "ModelBuilder.h"

// testing
#include "Test.h"

// File I/O and Initialization
#include "ComputationalSettings.h"
#include "InputDataManager.h"

// finite elements
#include "NumIntegral_dNT_op_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "NumIntegral_SetRHS_to_Zero.h"

// finite volumes
// CVFEM PHX
#include "CVFEM_MathOperatorRHS.h"
#include "CVFEM_PHX_Scheme.h"
#include "WellConfigurationFile.h"

// visitors
#include "ConductivityVisitor.h"
#include "ComputeGravityTermVisitor.h"
#include "PermeabilityVisitor.h"
#include "PoreVolumeVisitor.h"

// eos

// Solver
#include "LUdcmp_Solver.h"
#include "EigenSolver.h"

// Output
#include "VTU_Interface.h"
#include <time.h>

// compare values
#include "compareFloats.h"


using namespace std;
using namespace csmp;


namespace csmp
{

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
    CVFEM_1D_VVCase ( const char* argv , const uint32_t test_type);
    ~CVFEM_1D_VVCase();
    void run();
private:

    const char* geometry_name_;
    const char* config_file_name_;
    const char* vars_name_;
    const char* system_file_list_name_;
    const char* recipes_file_list_name_;
    list<string> output_props_;

    const uint32_t test_type;

    enum{dim=1};

    Model<dim>* model;
    PropertyDatabase<dim>*  pd_ref;
    ModelBuilder<dim> builder;
    VTU_Interface<dim>* vtu;
    ComputationalSettings  run_settings;

    //! CVFEM instance
    CVFEM_PHX_Scheme<dim>* CVFEM_PHX;

    //! Permeability Visitor
    PermeabilityVisitor<dim>* permeability_visitor;

    //! Chemical diffusion operators and variables
    PDE_Integrator<dim>              diffusion;

    Index pore_diff_key, eff_diff_key, porosity_key;
    ScalarVariable diff_coeff, porosity;

    //!time variables
    double dt, max_time, vtu_dt;
    size_t timestep, output_counter;

    int CalculateMassConductivity();
    int CalculateMassGravityTerm();
    int InitializeFluidPropertiesLinearPressure();
    void InitializeAllNaNsToZero();
    void ScalePermeabilityWithHalite();
    void PrintProgressToScreen();

    // functions to compare data from two vset files
    // used for Benchmarks
    std::map<Point<1U>, std::vector<double>> BuildPropertyMap( const Model<1U>& sg, const std::string& property_name );
    bool ComparePropertyMaps( const std::map<Point<1U>, std::vector<double>>& map_test,
                              const std::map<Point<1U>, std::vector<double>>& map_benchmark,
                              const std::string& property_name,
                              std::string& error_message );



    const double  year   = 31536000.0;     //  in sec
    const double  day    = 86400.0;        //  in sec
    const double  hour   = 3600.0;         //  in sec
    const double  minute = 60.0;           //  in sec
};

} // csmp

#endif // CVFEM_1D_VVCase
