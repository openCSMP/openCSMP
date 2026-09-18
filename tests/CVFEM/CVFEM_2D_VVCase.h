#ifndef CVFEM_2D_VVCase_H
#define CVFEM_2D_VVCase_H


#include <chrono>
#include "iostream"

// model
#include "ANSYS_Model2D.h"
#include "ModelTime.h"
#include "ModelBuilder.h"

// testing
#include "Test.h"

// File I/O and Initialization
#include "ComputationalSettings.h"
#include "InputDataManager.h"

#include "CSMP_mathUtilities.h"

// finite elements
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_NT_op_dNi_dV.h"

// finite volumes
// CVFEM PHX
#include "CVFEM_PHX_Scheme.h"
#include "WellConfigurationFile.h"


// visitors
#include "ConductivityVisitor.h"
#include "ComputeGravityTermVisitor.h"
#include "ComputeSinglePhaseGravityTermVisitor.h"
#include "PermeabilityVisitor.h"
#include "PoreVolumeVisitor.h"

#include "Boundary.h"

// eos

// Solver
#include "LUdcmp_Solver.h"

// Output
#include "VTU_Interface.h"
#include <time.h>

// compare values
#include "compareFloats.h"

using namespace std;
using namespace csmp;


namespace csmp
{

/** CVFEM_2D_VVCase : 2D Test Case
  =================================
  Mesh:       rectangle_coarse
  Model:      2D
  Test:       PHX_Scheme
  BC:         TBC
  Criterion:
  =================================
  */

class CVFEM_2D_VVCase : public Test
{
public:
    CVFEM_2D_VVCase (int argc, char **argv, const uint32_t test_type);
    ~CVFEM_2D_VVCase();
    void run();
private:

    const char* geometry_name_;
    const char* config_file_name_;
    const char* vars_name_;
    const char* system_file_list_name_;
    const char* recipes_file_list_name_;
    list<string> output_props_;

    const uint32_t test_type;

    enum{dim=2};

    Model<dim>* model;
    PropertyDatabase<dim>*  pd_ref;
    ModelBuilder<dim> builder;
    ComputationalSettings  run_settings;
    VTU_Interface<dim>* vtu;

    //! CVFEM instance
    CVFEM_PHX_Scheme<dim>* CVFEM_PHX;

    //! Permeability Visitor
    PermeabilityVisitor<dim>* permeability_visitor;



    string
        output_prefix_,
        output_name_,
        restart_file_,
        model_ID_path_,
        output_path,
        save_file_;

    double
        dt,
        max_time,
        time_multiplier;

    long
        output_increment,
        model_output,
        time_label;

    size_t
        timestep;

    bool
        with_gravity = false,
        stop = false;


    // used for time conversion
    enum TimeUnit
    {
        SECONDS,
        MINUTES,
        HOURS,
        DAYS,
        YEARS
    };

    const double  year   = 31536000.0;     //  in sec
    const double  day    = 86400.0;        //  in sec
    const double  hour   = 3600.0;         //  in sec
    const double  minute = 60.0;           //  in sec

    TimeUnit time_unit  = YEARS; // Time unit used in main-file and config-file

    void InitializeModel();
    void InitializeHydrostaticPressure();
    void InitializeTotalHeatCapacity();
    void InitializeParameters();
    void InitializeAllNaNsToZero();
    void ApplyBottomHeatFlux();
    void OutputAndBookKeeping();
    void PrintProgressToScreen();
    void OutputToVTU(Model<dim>&, std::string model_name, const std::list<std::string>& props, size_t timestep, bool with_fracture) const;
    double maxDifferenceScalarNodeProperty( Model<dim>& mdl, const char* snp1, const char* snp2 );

    // functions to compare data from two vset files
    // used for Benchmarks
    std::map<Point<dim>, std::vector<double>> BuildPropertyMap( const Model<dim>& sg, const std::string& property_name );
    bool ComparePropertyMaps( const std::map<Point<dim>, std::vector<double>>& map_test,
                             const std::map<Point<dim>, std::vector<double>>& map_benchmark,
                             const std::string& property_name,
                             std::string& error_message );
};

} // csmp

#endif // CVFEM_2D_VVCase
