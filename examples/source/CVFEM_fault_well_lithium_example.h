// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_fault_well_lithium_example_H
#define CVFEM_fault_well_lithium_example_H

#include "Example.h"
#include "ModelTimeToInteger.h"   // class template, needed by WriteTimestepOutput

namespace csmp { template<uint32_t> class PoreVolumeVisitor; }

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class CVFEM_PHX_Scheme;

class  CVFEM_fault_well_lithium_example : public Example {

public:
    virtual void Run();
    virtual void Specifications();

    std::vector<std::string> wells_list = {};

private:

    static constexpr uint32_t MODEL_DIM = 3U;

    /** Scalars that define the case, gathered so the setup steps extracted from
        Run() do not need long parameter lists. Populated once in Run(), after the
        time-unit scaling, and read-only thereafter. */
    struct RunParameters {
        double thickness;             ///< out-of-plane thickness [m], normal elements
        double thickness_LD;          ///< fault aperture [m], lower-dimensional elements
        double pressure_top;          ///< TOP boundary fluid pressure [Pa]
        double temperature_top;       ///< TOP boundary temperature [oC]
        double thermal_conductivity;  ///< uniform value used by the transient run [W/m/K]
        double heat_capacity_rock;    ///< [J/kg/K]
        double density_rock;          ///< [kg/m3]
        double compressibility_rock;  ///< [1/Pa]
        double porosity;              ///< background porosity [-]
        double total_heat_capacity;   ///< cp_rock * rho_rock * (1 - porosity)
        double salinity;              ///< background salinity [-]
        double largest_time_step;     ///< [s] — already scaled by the time unit
        bool   with_well;
        bool   with_lithium;
    };
    enum class TimeUnit { SECONDS, MINUTES, HOURS, DAYS, YEARS };
    double maxDifferenceScalarNodeProperty(Model<MODEL_DIM> &mdl, const char *snp1, const char *snp2);
    void OutputToVTU(Model<MODEL_DIM> &, const std::string &model_name, const std::list<std::string> &props, size_t timestep) const;
    void OutputToVTU_wells(Model<MODEL_DIM> &, const std::string &model_name, const std::list<std::string> &props, size_t timestep) const;
    void assignKtensor(Model<MODEL_DIM> &model);

    void SpecialSalinity(Model<MODEL_DIM> &model) const;
    void CalculateNodalPermeability(Model<MODEL_DIM> &model);
    void CopyInsideSplitBtoMiddle(Model<MODEL_DIM> &model);
    void HeatFluxBottom(Model<MODEL_DIM> &model, double heat_flux_bottom);

    /** Creates the properties PhysicalVariablesBenchmarks.txt does not declare and
        assigns every field its starting value. Must run AFTER the split boundary
        exists, so split-region and well properties can be created for the regions
        the split produced. */
    void InitialiseProperties(Model<MODEL_DIM> &model, const RunParameters &p);

    /** Per-region values and geometry-derived fields: porosity, the lithium
        starting distribution, depth/coordinate assignments, the salinity profile,
        the basal heat flux, nodal permeability and pore/bulk volumes. Runs after
        InitialiseProperties(), whose uniform values it partly overwrites. */
    void SetUpRegionsAndGeometry(Model<MODEL_DIM> &model, const RunParameters &p,
                                 double heat_flux_bottom,
                                 PoreVolumeVisitor<MODEL_DIM> &Pore_Volume_Visitor);

    /** Fills the two VTU field lists. The well list starts as a copy of the
        reservoir list, so fields added before that point appear in both. */
    void DefineOutputVariables(const RunParameters &p,
                               std::list<std::string> &output_variables,
                               std::list<std::string> &output_variables_well) const;

    /** Builds ONE well and sets its starting operating conditions: wellhead
        pressure and temperature, water table off, wellhead-pressure control. Call
        it once per well; the time loop switches them to target-rate control
        afterwards.

        Does NOT read the configuration file. The well description — geometry,
        completions, skin, physics switches — was read in Run() and handed to the
        scheme's constructor. Initialize_well() inside turns that description into
        a discretised well, which is why this runs only once the reservoir fields
        exist. */
    void StartWell(CVFEM_PHX_Scheme<MODEL_DIM> &scheme,
                   const std::string &well_name,
                   double wellhead_pressure,
                   double wellhead_temperature) const;

    /** Creates the fault split boundary from FAULT2D_HOST, inserts the middle
        (fault) region and pulls the two sides apart by the aperture, then repairs
        the well regions whose line elements referenced the duplicated nodes. */
    void CreateFaultSplitBoundary(Model<MODEL_DIM> &model, const RunParameters &p);

    /** Writes the reservoir and well VTUs if this step is due for them, and
        advances the two schedule counters (hence the references). The schedule is
        "first step at or past the threshold", so a long step can skip a slot. */
    void WriteTimestepOutput(Model<MODEL_DIM> &model, const std::string &output_name,
                             TimeUnit time_unit, ModelTimeToInteger<double> &time_conversion,
                             double model_time,
                             const std::list<std::string> &output_variables,
                             const std::list<std::string> &output_variables_well,
                             bool with_well,
                             long &model_output, long vtk_increment,
                             long &detailed_output, long detailed_vtk_increment) const;

    void SetPermeability(Model<MODEL_DIM> &model) const;

    /** Rebuilds every permeability-derived field from the per-region values:
        background value -> SetPermeability -> well override -> vertical/horizontal
        copies -> tensor -> nodal permeability. Called once during setup and again
        after each timestep, so the two paths cannot drift apart. */
    void RebuildPermeabilityFields(Model<MODEL_DIM> &model);

    /** Writes the min/max of every property in the database to a text file.
        Run after initialisation: it is the cheapest way to compare two runs
        (different solver backend, different branch) field by field. */
    void OutputAllVariableRanges(Model<MODEL_DIM> &model,
                                 const std::string &output_name,
                                 const std::string &tag) const;
    void SetPorosity(Model<MODEL_DIM> &model) const;

    void PAUSE();
};

} // csmp

#endif