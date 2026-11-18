// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVFEM_SQUARE_CAVITY_CONVECTION_TEST_H
#define CVFEM_SQUARE_CAVITY_CONVECTION_TEST_H

/** @file CVFEM_SquareCavityConvection_Test.h
 *  @author Samuel Scott
 *  @brief Thermal convection in a closed porous square: openCSMP against the CSMP benchmark of Geiger et al. (2006).
 *  @details Purpose: assess how the current CVFEM PHX scheme performs on the closed-square thermal convection
 *  benchmark that Geiger, Driesner, Heinrich & Matthai (2006, "Multiphase thermohaline convection in the Earth's
 *  crust: II. Benchmarking and application of a finite element - finite volume solution technique with a NaCl-H2O
 *  equation of state", Transport in Porous Media 63, 435-461, section 3.2.2 and Table II) used to verify the
 *  predecessor of this code. In particular it tests the treatment of the buoyancy (gravity) term, which sets the
 *  steady convective heat transport.
 *
 *  Reference values (Geiger et al. 2006, Table II; steady single-cell convection, average Nusselt number at the top):
 *    Ra_th  Caltagirone 1975  Steen & Aidun 1988  Schoofs 1999  Geiger et al. 2006 (CSMP predecessor)
 *    100    2.65              2.651               2.643         2.560
 *    200    3.81              3.810               3.806         3.818
 *    300    4.52              4.523               4.511         4.503
 *  The test compares with Caltagirone (1975) and passes within nusselt_tolerance_ (5 %); Geiger et al. (2006) were
 *  3.4 % low at Ra 100.
 *
 *  Set-up as in Geiger et al. (2006): square domain, impermeable on all four sides, fixed temperature at the base
 *  (hot) and top (cold), insulated side walls. "Closed top" means no flow through the top: the fluid pressure is
 *  fixed at one node only (top-left corner) as the pressure reference, not along the whole top, which would make
 *  the top an open (fixed-pressure) boundary with a different onset (Ra_c ~ 27 instead of 4 pi^2) and different
 *  Nusselt numbers. Nusselt number Nu = -dTbar/dzbar at the top (eq. 11), Rayleigh number
 *  Ra = (rho_top - rho_base) g H k / (mu kappa), kappa = lambda / (rho c_p)_fluid (eq. 10).
 *
 *  Differences from Geiger et al. (2006): they solved the Boussinesq problem with a density linear in T and constant
 *  viscosity, and give no temperatures (only Ra matters there). openCSMP uses the full pure-water equation of state
 *  and is not Boussinesq, so the temperatures matter. They are chosen to stay close to the Boussinesq limit:
 *  140-150 C at 5 MPa, where the viscosity varies by +-3.6 % across the cell and the density nearly linearly
 *  (rho_top - rho_base = 9.07 kg m-3); at 10-20 C the viscosity would vary by ~30 %. The permeability of each case
 *  is set for the target Ra from IAPWS water at 145 C and 5 MPa and scaled by 1.056 so that the Ra realised by the
 *  model (whose density contrast includes the hydrostatic pressure difference) is the target; the test reports it.
 *  Rock heat capacity is small (100 J kg-1 K-1): steady Nu does not depend on it, and a small value shortens the
 *  transient. The test calls AddFluidContributionToHeatCapacity(true), as the CVFEM examples do; the equilibration
 *  (NaClH2OPropertiesVisitorPHX) credits the fluid's share of the conducted heat as well as the rock's.
 *  Initial state: linear conductive temperature with a 1 % single-cell perturbation, hydrostatic pressure.
 *
 *  Requires a structured mesh (node rows of equal size) for the second-order top-boundary temperature gradient.
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

#include <memory>
#include <string>

namespace csmp {

/// Thermal convection in a closed porous square heated from below: steady Nusselt number versus published values.
class CVFEM_SquareCavityConvection_Test : public Test {
  public:
    /// mesh_name and config_name are file base names; facet_mean_gravity_density switches
    /// CVFEM_PHX_Scheme::UseFacetMeanGravityDensity; reference_nusselt is the published Nusselt number.
    CVFEM_SquareCavityConvection_Test(const std::string& mesh_name, const std::string& config_name, bool facet_mean_gravity_density,
                                      double reference_nusselt);
    ~CVFEM_SquareCavityConvection_Test() override;

    void run() override;

  private:
    /// Reads mesh and configuration, builds the CVFEM PHX scheme, permeability tensor and pore volumes.
    void CreateModel();
    /// Fixes the fluid pressure at the top-left corner node only (closed top), at its current value.
    void AnchorPressureAtOneNode();
    /// Linear conductive temperature between the base and top values plus a small single-cell perturbation.
    void InitializeTemperature();
    /// Sets every property that is NaN after configuration to zero (as tests/CVFEM/CVFEM_2D_VVCase).
    void InitializeAllNaNsToZero();
    /// Steady hydrostatic pressure with the density from PTX, iterated until the change is below 1 Pa.
    void InitializeHydrostaticPressure();
    /// Largest absolute nodal difference between two scalar properties.
    double MaxDifference(const std::string& first, const std::string& second) const;
    /// Horizontally averaged Nusselt number at the top: second-order one-sided gradient over the top three rows.
    double TopNusseltNumber() const;
    /// Rayleigh number (rho_top - rho_base) g H k / (mu kappa) from the model's current properties.
    double RayleighNumber() const;

    std::string mesh_name_;
    std::string config_name_;
    std::string vars_name_;
    bool facet_mean_gravity_density_;
    double reference_nusselt_;
    double nusselt_tolerance_;
    double steady_tolerance_;
    double temperature_base_;
    double temperature_top_;
    double y_base_;
    double y_top_;

    std::unique_ptr<Model<2U>> model_;
    PropertyDatabase<2U>* property_database_;
    std::unique_ptr<CVFEM_PHX_Scheme<2U>> cvfem_phx_;
    std::unique_ptr<PermeabilityVisitor<2U>> permeability_visitor_;
    ComputationalSettings run_settings_;

    static constexpr double year_{31536000.};
    static constexpr double gravity_{9.81};
};

} // namespace csmp

#endif // CVFEM_SQUARE_CAVITY_CONVECTION_TEST_H
