// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_topo_magma_air_example_H
#define CVFEM_topo_magma_air_example_H

#include "Example.h"
#include "ModelTimeToInteger.h"   // class template, used in the declarations below

namespace csmp {

template<uint32_t> class Model;

class  CVFEM_topo_magma_air_example : public Example {

public:
    virtual void Run();
    virtual void Specifications();

private:
    // NOT named DIM: SlopeMechanics_Example.h has `#define DIM 2`, and Examples.h
    // includes it before this header, so a member called DIM would be rewritten by
    // the preprocessor into `static ... 2 = 2U;`.
    static constexpr uint32_t MODEL_DIM = 2U;

    /// Time unit of the run parameters below; they are converted to seconds
    /// by time_multiplier. Scoped so it cannot collide with the enum of the
    /// same name in aux_2D_CVFEM.h.
    enum class TimeUnit { SECONDS, MINUTES, HOURS, DAYS, YEARS };
    double maxDifferenceScalarNodeProperty(Model<MODEL_DIM> &mdl, const char *snp1, const char *snp2);
    void OutputToVTU(Model<MODEL_DIM> &, const std::string &model_name, const std::list<std::string> &props, size_t timestep) const;
    void ComputeLiquidusSolidus(Model<MODEL_DIM> &model, const char *liquidus, const char *solidus, const char *depth, double TL, double TS);

    void HeatFluxBottom(Model<MODEL_DIM> &model, double heat_flux_bottom);

    /** The feature switches that decide which fields are written. Passed as a
        struct so DefineOutputVariables() keeps one argument as more are added. */
    struct OutputFlags {
        bool with_air;
        bool with_magma_chamber;
        bool with_magma_model;
        /// Either of these enables the gas-phase mixture diagnostics.
        bool treat_air_and_vapor_as_mixture;
        bool treat_air_and_vapor_as_hydraulic_mixture;
    };

    /** Fills the VTU field list. Which fields appear depends on the features in
        play — see OutputFlags. */
    void DefineOutputVariables(const OutputFlags &p,
                               std::list<std::string> &output_variables) const;

    /** Model time expressed in the run's time unit, as the integer used to label
        output files. Was duplicated as a switch at two points in Run(). */
    long TimeLabel(TimeUnit time_unit, ModelTimeToInteger<double> &time_conversion,
                   double model_time) const;

    /** Writes the VTU if this step is due for it, and advances the schedule
        counter (hence the reference). The schedule is "first step at or past the
        threshold", so a long step can skip a slot entirely. */
    void WriteTimestepOutput(Model<MODEL_DIM> &model, const std::string &output_name,
                             TimeUnit time_unit, ModelTimeToInteger<double> &time_conversion,
                             double model_time,
                             const std::list<std::string> &output_variables,
                             long &model_output, long vtk_increment) const;

    /** Writes the min/max of every property in the database to a text file. Run
        after initialisation: diffing two of these is the cheapest way to compare
        two runs (different solver, different branch) field by field. */
    void OutputAllVariableRanges(Model<MODEL_DIM> &model,
                                 const std::string &output_name,
                                 const std::string &tag) const;


};

} // csmp

#endif