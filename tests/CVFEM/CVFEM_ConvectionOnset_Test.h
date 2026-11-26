// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_CONVECTION_ONSET_TEST_H
#define CVFEM_CONVECTION_ONSET_TEST_H

/** @file CVFEM_ConvectionOnset_Test.h
 *  @author Samuel Scott
 *  @brief Growth or decay rate of a convective mode in a closed porous box heated from below, against linear
 *  stability theory (Horton-Rogers-Lapwood problem).
 *  @details Purpose: a transient, analytical check of the buoyancy (gravity) term together with conduction and the
 *  bulk heat capacity. The Horton-Rogers-Lapwood problem (Horton & Rogers 1945, Lapwood 1948) is proposed as a
 *  benchmark for density-driven flow codes by Weatherill, Simmons, Voss & Robinson (2004, "Testing density-dependent
 *  groundwater models: two-dimensional steady state unstable convection in infinite, finite and inclined porous
 *  layers", Advances in Water Resources 27, 547-562), who test the onset at Ra_c = 4 pi^2. This test measures the rate
 *  itself, below and above onset.
 *
 *  Set-up: box of length L = 2 H and height H, impermeable on all sides, fixed temperature at the base (hot) and top
 *  (cold), insulated side walls; fluid pressure fixed at one top corner node only as the pressure reference (closed
 *  top, as tests/CVFEM/CVFEM_SquareCavityConvection_Test). Initial state: the conductive (linear) temperature plus
 *  the mode theta = epsilon dT cos(pi x / H) sin(pi z / H), and hydrostatic pressure. This mode satisfies every
 *  boundary condition and is an exact eigenmode of the linearised Darcy-Boussinesq equations, with stream function
 *  proportional to sin(pi x / H) sin(pi z / H). Its amplitude evolves as A0 exp(sigma tau) with
 *      sigma = Ra a^2 / (a^2 + pi^2) - (a^2 + pi^2) = Ra / 2 - 2 pi^2   (a = pi),
 *      tau = lambda t / (C_bulk H^2),   C_bulk = (1 - phi) rho_r c_r + phi rho_f c_f,
 *      Ra = rho_f beta_T dT g H k (rho c)_f / (mu lambda),   beta_T = -(1 / rho) (d rho / d T) at constant pressure,
 *  with lambda the bulk thermal conductivity (e.g. Nield & Bejan, Convection in Porous Media). The mode decays below
 *  Ra_c = 4 pi^2 and grows above it.
 *
 *  Measurement: every case runs twice, with the mode added and subtracted. Half the difference of the two mode
 *  amplitudes is the linear response; half the sum, which does not depend on the sign, collects any source of the
 *  mode (e.g. spurious buoyancy flow of the discretisation) and the quadratic nonlinearity, and is reported. The
 *  amplitude is the least-squares projection of T - T_conductive onto the mode after every time step. Each run stops
 *  when |sigma| tau = 3; the initial amplitude is 1 % of dT for a decaying mode and 1 % of dT times exp(-3) for a
 *  growing one, so the amplitude never exceeds 1 % of dT (linear regime). ln A is fitted with a straight line against
 *  tau over the last 75 % of the run. The fitted rate is expressed as the Rayleigh number it implies,
 *  Ra_sigma = 2 (sigma + 2 pi^2), and the test passes if |Ra_sigma - Ra| < rayleigh_tolerance_ * Ra_c (3 % of
 *  Ra_c = 4 pi^2, i.e. |sigma - sigma_theory| < 0.59). The tolerance is relative to Ra_c, not Ra, because near onset
 *  sigma is a small difference of two large terms. The time step is capped so that |sigma| dtau <= 0.005: near
 *  onset sigma is a small difference of the conduction and buoyancy rates, each of which carries its own
 *  (first-order) time-stepping error. Fields of the first run are written at 40 checkpoints for inspection.
 *
 *  Finite amplitude (finite_amplitude = true, test-207e): below onset the Darcy-Boussinesq layer is also globally
 *  stable: the energy method gives the same threshold Ra_c = 4 pi^2 for all initial data (Straughan 2004, The Energy
 *  Method, Stability, and Nonlinear Convection, 2nd ed., Springer, Applied Mathematical Sciences 91, sec. 4.2.1,
 *  pp. 64-68: R_L = R_E for the Darcy equations; shown there for a periodic layer, and the impermeable, insulated side
 *  walls of the box leave the same boundary terms zero), so a perturbation of ANY size must decay. The case runs
 *  with an initial amplitude of 10 % of dT (both signs) and passes if each run has decayed by at least
 *  decay_factor_ when |sigma| tau = 3 (linear theory: a factor exp(3) = 20). The rate check is not applied (not
 *  linear).
 *
 *  Differences from the analytical problem: openCSMP uses the full pure-water equation of state, not the Boussinesq
 *  approximation. As in the square-cavity test, 140-150 C at 5 MPa keeps it close to that limit (viscosity +-3.6 %).
 *  Ra uses the THERMAL density contrast, from the model's own equation of state at constant pressure (the hydrostatic
 *  compression of the water with depth does not drive flow); the top-to-base density difference, which includes it,
 *  is reported for information.
 *
 *  Requires a structured mesh of a box whose length is an even multiple of its height. Initialisation follows
 *  tests/CVFEM/CVFEM_2D_VVCase.
 *  @date 18.11.2026
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

#include <memory>
#include <string>
#include <vector>

namespace csmp {

/// Linear growth or decay of the first convective mode in a closed porous box versus Horton-Rogers-Lapwood theory.
class CVFEM_ConvectionOnset_Test : public Test {
  public:
    /// mesh_name and config_name are file base names; facet_mean_gravity_density switches
    /// CVFEM_PHX_Scheme::UseFacetMeanGravityDensity; finite_amplitude: 10 % perturbation, decay check (below onset).
    CVFEM_ConvectionOnset_Test(const std::string& mesh_name, const std::string& config_name,
                               bool facet_mean_gravity_density, bool finite_amplitude = false);
    ~CVFEM_ConvectionOnset_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, builds the CVFEM PHX scheme, permeability tensor and pore volumes.
    void CreateModel();
    /// Box extent from the node coordinates and the fixed base and top temperatures.
    void ReadGeometryAndBoundaryTemperatures();
    /// Fixes the fluid pressure at the top-left corner node only (closed top), at its current value.
    void AnchorPressureAtOneNode();
    /// Conductive (linear) temperature plus amplitude * cos(pi x / H) sin(pi z / H) on all free nodes.
    void InitializeTemperature(double amplitude);
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Steady hydrostatic pressure with the density from PTX, iterated until the change is below 1 Pa.
    void InitializeHydrostaticPressure();
    /// Largest absolute nodal difference between two scalar properties.
    double MaxDifference(const std::string& first, const std::string& second) const;
    /// Mean of rho_f beta_T over the nodes, from the equation of state at constant pressure (temperature raised by
    /// 0.1 K and restored; the fluid properties are recomputed both times).
    double ThermalDensityDerivative();
    /// Bulk properties for Ra and tau: thermal conductivity, C_bulk, fluid (rho c), viscosity (node and cell means).
    void ReadBulkProperties();
    /// Amplitude of the mode after every time step of one run, against tau.
    struct ModeHistory {
      std::vector<double> tau;
      std::vector<double> amplitude;   ///< [K]
    };

    /// One run with the mode added (sign 1) or subtracted (sign -1); writes the fields at 40 checkpoints if
    /// write_fields.
    ModeHistory RunCase(double sign, bool write_fields);
    /// Discards any previous model and builds a fresh one (CreateModel), model time 0.
    void RebuildModel();
    /// Conductive temperature, pressure anchor, hydrostatic pressure (no perturbation yet).
    void InitializeConductiveState();
    /// Ra from the thermal density contrast at constant pressure, sigma from linear theory, and the time unit.
    void ComputeRayleighNumberAndTimeUnit();
    /// Run length (|sigma| tau = growth_span_), largest time step (|sigma| dtau = time_step_fraction_), initial
    /// amplitude.
    void ChooseRunLengthAndTimeStep();
    /// Adds the mode with the given amplitude to the conductive temperature and re-equilibrates the pressure.
    void ImposeMode(double amplitude);
    /// Top-to-base liquid density difference of the model (includes hydrostatic compression), for information.
    double TopToBaseDensityDifference() const;
    /// Writes Ra, the density contrasts, sigma (theory), the time unit, run length and largest time step.
    void ReportSetUp() const;
    /// Largest time step, transient preparation and a first time step of a tenth of the largest.
    void StartTransient();
    /// Steps to the run length, recording the mode amplitude after every step.
    ModeHistory AdvanceAndRecordMode(bool write_fields);
    /// Straight-line fit of ln((A+ - A-) / 2) against tau over the last 75 % of the run.
    double FitLinearGrowthRate(const ModeHistory& plus, const ModeHistory& minus) const;
    /// Largest |A+ + A-| / 2 during the run, relative to the initial amplitude (spurious excitation of the mode).
    double LargestSignIndependentPart(const ModeHistory& plus, const ModeHistory& minus) const;
    /// True if A+ - A- stays positive throughout the run.
    bool LinearPartStaysPositive(const ModeHistory& plus, const ModeHistory& minus) const;
    /// Rate cases: fitted sigma, as the Ra it implies, within rayleigh_tolerance_ * Ra_c of the model's Ra.
    void CheckGrowthRate(const ModeHistory& plus, const ModeHistory& minus);
    /// Finite-amplitude case: both runs decayed by at least decay_factor_ (below onset).
    void CheckFiniteAmplitudeDecay(const ModeHistory& plus, const ModeHistory& minus);
    /// Least-squares amplitude of T - T_conductive on the mode cos(pi x / H) sin(pi z / H), weighted by bulk volume.
    double ModeAmplitude() const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    bool facet_mean_gravity_density_;
    bool finite_amplitude_;
    double rayleigh_tolerance_;
    double amplitude_;              ///< initial mode amplitude [K]
    double rayleigh_;
    double sigma_theory_;
    double time_unit_;              ///< [s] per unit tau
    double thermal_density_contrast_;   ///< rho_f beta_T dT [kg m-3]
    double max_time_;               ///< run length [s]
    double largest_time_step_;      ///< [s]
    bool ran_planned_duration_;
    double temperature_base_;
    double temperature_top_;
    double x_min_;
    double y_base_;
    double y_top_;
    double conductivity_;           ///< bulk thermal conductivity [W m-1 K-1]
    double bulk_heat_capacity_;     ///< C_bulk [J m-3 K-1]
    double fluid_heat_capacity_;    ///< (rho c)_f [J m-3 K-1]
    double viscosity_;              ///< [Pa s]
    double permeability_;           ///< [m2]

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    ComputationalSettings run_settings_;

    static constexpr double year_{31536000.};
    static constexpr double gravity_{9.80665};             ///< as UpwindControlVisitor
    static constexpr double peak_relative_amplitude_{1.e-2};   ///< larger of initial and final amplitude / dT
    static constexpr double finite_relative_amplitude_{0.1};   ///< initial amplitude / dT with finite_amplitude
    static constexpr double decay_factor_{5.};                 ///< required decay with finite_amplitude
    static constexpr double growth_span_{3.};             ///< run until |sigma| tau = growth_span_
    static constexpr double time_step_fraction_{0.005};   ///< largest time step: |sigma| dtau = time_step_fraction_
};

} // namespace csmp

#endif // CVFEM_CONVECTION_ONSET_TEST_H
