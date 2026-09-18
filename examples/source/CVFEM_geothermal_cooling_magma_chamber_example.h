// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_geothermal_cooling_magma_chamber_example_H
#define CVFEM_geothermal_cooling_magma_chamber_example_H

/** @file CVFEM_geothermal_cooling_magma_chamber_example.h
 *  @author Jonas Köpping
 *  @brief
 *  @details Uses CVFEM for pressure-enthalpy-salinity (PHX) formulation for flow and transport.
 *  @date 15.09.2026
 */

// Example
#include "Example.h"

// Model
#include "ModelBuilder.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "InputDataManager.h"

// File I/O and Initialization
#include "ComputationalSettings.h"
#include "InputDataManager.h"

#include "CSMP_mathUtilities.h"

// Finite elements
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "PointSource_rhsop.h"
#include "NumIntegral_NT_op_dNi_dV.h"

// Finite volumes
// CVFEM PHX
#include "CVFEM_PHX_Scheme.h"
#include "WellConfigurationFile.h"

// Visitors
#include "ConductivityVisitor.h"
#include "PermeabilityVisitor.h"
#include "PoreVolumeVisitor.h"
#include "FailureModeVisitor.h"


//HeatFluxBottom
#include "VectorVariable.h"
#include "Region.h"
#include "Element.h"
#include "ScalarVariable.h"
#include "CSMP_definitions.h"

//
#include "VTU_Interface.h"
#include <time.h>
#include "ModelTime.h"
#include "iostream"
#include <chrono>

// Solver
#include "CVFEM_SolverChoice.h"

#ifdef CSMP_WITH_PETSC_SOLVER
#include <petscsys.h>
#endif


namespace csmp { template<uint32_t> class PoreVolumeVisitor; }

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class CVFEM_PHX_Scheme;

template<uint32_t dim>
class CVFEM_geothermal_cooling_magma_chamber_example : public Example {
public:
    CVFEM_geothermal_cooling_magma_chamber_example(int argc = 0, char** argv = nullptr);
    ~CVFEM_geothermal_cooling_magma_chamber_example ();

    virtual void Run();
    virtual void Specifications();

private:

    Model<dim>* model;
    PropertyDatabase<dim>*  pd_ref;
    CVFEM_PHX_Scheme<dim>* CVFEM_PHX;
    ComputationalSettings  run_settings;
    VTU_Interface<dim>* vtu;


    //! Visitor
    PermeabilityVisitor<dim>                *permeability_visitor;
    PoreVolumeVisitor<dim>                  *Pore_Volume_Visitor;
    FailureModeVisitor<dim>                 *failure_mode;


    SolverKind solver_kind, well_solver_kind;
    // SolverBundle init_bundle;

    std::vector<std::string> wells_list = {};

    std::string
        geometry_name_, config_file_name_, output_path_, output_name_, vars_name_;


    double
        dt, time_multiplier, max_time, largest_time_step, initial_time_step, output_increment, model_output;

    bool
        open_top,
        with_magma_chamber,
        depth_dependent,
        temperature_dependent,
        pore_fluid_factor_dependent,
        hydrofracturing,
        anisotropic_k,
        T_dependent_differential_stress;

    std::list<string> output_variables;
    std::vector<std::string> output_regions_list;

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


    static Model<dim>* CreateANSYSModel(const std::string& mesh,
                                        const std::string& region,
                                        const std::string& vars);

    bool GetBoolFromConfigFile(const std::string& var);
    void SetTimeVariables();
    void InstatiateCVFEM();
    void InstantiatePermeabilityVisitor();
    void InstantiateFailureModeVisitor();
    void ActivateFracturing(bool hydrofracturing);
    void InitializeTemperature();
    void AddMagmaticIntrusion();
    void InitializeHydrostaticPressure();
    void InitializePressure();
    void InitializeTotalHeatCapacity();
    void InitializeAllNaNsToZero();
    void InitializeParameters();

    double maxDifferenceScalarNodeProperty(const char *snp1, const char *snp2);

    void ApplyBottomHeatFlux();
    void FluxPerSecond (double dt);
    void PrintProgressToScreen();
    void DefineOutputVariables ();
    void CreateOutputFile(std::vector<std::string> output_regions_list);
    void OutputToVTU(string model_name,
                     const list<string> &props,
                     size_t timestep,
                     std::vector<std::string> additional_region) const;


};

} // csmp

#endif