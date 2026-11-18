// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_HYDROSTATIC_COLUMN_TEST_H
#define CVFEM_HYDROSTATIC_COLUMN_TEST_H

/** @file CVFEM_HydrostaticColumn_Test.h
 *  @author Samuel Scott
 *  @brief Well-balancedness test of the CVFEM PHX gravity term: a stably conducting water column must stay at rest.
 *  @details Pure-water column, top temperature and pressure fixed, base temperature fixed, insulated no-flow side
 *  walls, initialised with the exact conductive (linear) temperature profile and the hydrostatic pressure. The
 *  Rayleigh number of the configuration is far below the onset of convection (4 pi^2, Horton and Rogers 1945,
 *  Lapwood 1948), so the exact solution is no flow at all times. Any Darcy velocity is discretisation error of the
 *  buoyancy term. The test runs the configured duration and requires the mass transported between interior control
 *  volumes (the transport scheme's "fluid flux integral", as a mean Darcy speed) to stay below max_velocity_fraction_
 *  times the column buoyancy velocity scale k g (rho_max - rho_min) / mu_min. The stored element "velocity liquid"
 *  (a visualisation quantity) and the transport at Dirichlet nodes are reported for information.
 *
 *  Initialisation follows tests/CVFEM/CVFEM_2D_VVCase.
 *  @date 17.11.2026
 */

#include "Test.h"

#include "ModelBuilder.h"
#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "ComputationalSettings.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "CVFEM_PHX_Scheme.h"
#include "WellConfigurationFile.h"
#include "PermeabilityVisitor.h"
#include "PoreVolumeVisitor.h"
#include "ModelTime.h"

#include <array>
#include <memory>
#include <string>

namespace csmp {

/// Hydrostatic water column with a conductive temperature gradient: the CVFEM PHX scheme must keep it at rest.
class CVFEM_HydrostaticColumn_Test : public Test {
  public:
    /// mesh_name and config_name are file base names (ANSYS mesh, <config_name>-configuration.txt);
    /// facet_mean_gravity_density switches CVFEM_PHX_Scheme::UseFacetMeanGravityDensity.
    CVFEM_HydrostaticColumn_Test(const std::string& mesh_name, const std::string& config_name, bool facet_mean_gravity_density);
    ~CVFEM_HydrostaticColumn_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, builds the CVFEM PHX scheme, permeability tensor and pore volumes.
    void CreateModel();
    /// Sets the linear conductive temperature profile between the top and base values on all free nodes.
    void InitializeConductiveTemperature();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Steady hydrostatic pressure with the density from PTX, iterated until the change is below 1 Pa.
    void InitializeHydrostaticPressure();
    /// Largest absolute nodal difference between two scalar properties.
    double MaxDifference(const std::string& first, const std::string& second) const;
    /// Largest liquid Darcy speed over all elements; stores its element centroid in location_of_max_speed_.
    double MaxLiquidSpeed();
    /// Largest time-averaged transport speed (fluid flux integral / run_time as a Darcy speed through half the
    /// control-volume perimeter) over interior nodes [0] and over nodes with a Dirichlet temperature or pressure [1].
    std::array<double, 2> MaxTransportSpeed(double run_time) const;
    /// Column buoyancy velocity scale k g (rho_max - rho_min) / mu_min from the current nodal state.
    double BuoyancyVelocityScale() const;
    /// Conductive temperature and hydrostatic pressure at rest; prepares the transient (first step a tenth of the largest).
    void InitializeRestingColumn(double largest_time_step);
    /// Steps to max_time; returns the largest stored element liquid speed over the run (counts steps_).
    double AdvanceAndTrackStoredVelocity(double max_time);
    /// Writes the transported-mass speeds, the stored-velocity maximum and the limit.
    void ReportResults(double velocity_scale, const std::array<double, 2>& transport_speed, double max_stored_speed) const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    bool facet_mean_gravity_density_;
    uint64_t steps_;
    double max_velocity_fraction_;
    std::array<double, 2> location_of_max_speed_;

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    ComputationalSettings run_settings_;

    static constexpr double year_{31536000.};
    static constexpr double gravity_{9.81};
};

} // namespace csmp

#endif // CVFEM_HYDROSTATIC_COLUMN_TEST_H
