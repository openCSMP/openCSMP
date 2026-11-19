// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_HydrostaticColumn_Test.cpp
 *  @author Samuel Scott
 *  @brief Well-balancedness test of the CVFEM PHX gravity term, see CVFEM_HydrostaticColumn_Test.h.
 *  @date 17.11.2026
 */

#include "CVFEM_HydrostaticColumn_Test.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace csmp {

CVFEM_HydrostaticColumn_Test::CVFEM_HydrostaticColumn_Test(const std::string& mesh_name, const std::string& config_name, bool facet_mean_gravity_density)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    facet_mean_gravity_density_(facet_mean_gravity_density),
    steps_(0U),
    max_velocity_fraction_(1.e-3),
    location_of_max_speed_{0., 0.},
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    run_settings_()
{
}


CVFEM_HydrostaticColumn_Test::~CVFEM_HydrostaticColumn_Test() = default;


void CVFEM_HydrostaticColumn_Test::CreateModel()
{
  ModelBuilder<2U> builder;
  builder.CreateVariablesFile(false, false, false, false, false, false);
  builder.CreateLimiterLogFile();

  model_ = std::make_unique<ANSYS_Model2D>(mesh_name_.c_str(), mesh_name_.c_str(), vars_name_.c_str());
  model_->InstantiateFiniteVolumes();
  property_database_ = &(model_->Database());

  InputDataManager<2U> model_configuration;
  model_configuration.ConfigureFromFile(*model_, config_name_.c_str(),
                                        false,   ///< region name from parameter range
                                        true,    ///< default property values
                                        true,    ///< regional property values
                                        true,    ///< boundary conditions for box-shaped model
                                        true,    ///< regional property conditions
                                        true,    ///< boundary conditions for arbitrary-shaped model
                                        run_settings_);

  const bool with_gravity{static_cast<bool>(model_->Read(property_database_->StorageKey("gravity flag")))};
  const std::vector<std::string> wells_list{};
  const std::map<std::string, WellConfiguration> well_configurations{};
  cvfem_phx_ = std::make_unique<CVFEM_PHX_Scheme<2U>>(*model_, with_gravity, true, wells_list, false, false,
                                                      SolverKind::PETSc, SolverKind::PETSc, well_configurations);
  const double heat_capacity_rock{model_->Region("Model").NodeVector().front()->Read(
      property_database_->StorageKey("nodal heat capacity rock"))};
  cvfem_phx_->SetRockHeatCapacity(heat_capacity_rock);
  cvfem_phx_->AddFluidContributionToHeatCapacity(true);
  cvfem_phx_->UseFacetMeanGravityDensity(facet_mean_gravity_density_);

  permeability_visitor_ = std::make_unique<PermeabilityVisitor<2U>>(*model_);
  permeability_visitor_->AssignPermeabilityTensor(*model_);

  PoreVolumeVisitor<2U> pore_volume_visitor(*model_, "nodal porosity", "bulk volume", "pore volume", "thickness");
  model_->Accept(pore_volume_visitor);
}


void CVFEM_HydrostaticColumn_Test::InitializeConductiveTemperature()
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  double y_bottom{std::numeric_limits<double>::max()};
  double y_top{std::numeric_limits<double>::lowest()};
  double t_bottom{0.};
  double t_top{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    y_bottom = std::min(y_bottom, (*node_it)->y());
    y_top = std::max(y_top, (*node_it)->y());
  }
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    if (std::fabs((*node_it)->y() - y_bottom) < 1.e-6) t_bottom = (*node_it)->Read(temperature_key);
    if (std::fabs((*node_it)->y() - y_top) < 1.e-6) t_top = (*node_it)->Read(temperature_key);
  }
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    if (temperature.Flag() == DIRICH)
      continue;
    const double fraction{((*node_it)->y() - y_bottom) / (y_top - y_bottom)};
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), t_bottom + fraction * (t_top - t_bottom)));
  }
  std::cerr << "CVFEM_HydrostaticColumn_Test::InitializeConductiveTemperature: linear profile from " << t_bottom
            << " C (base) to " << t_top << " C (top) over " << y_top - y_bottom << " m" << std::endl;
}


void CVFEM_HydrostaticColumn_Test::InitializeAllNaNsToZero()
{
  double minimum{0.};
  double maximum{0.};
  for (auto it = model_->Database().Begin(); it != model_->Database().End(); ++it) {
    model_->MinMaxOf((*it).first.c_str(), minimum, maximum);
    if (!std::isnan(minimum) && !std::isnan(maximum))
      continue;
    if ((*it).second.key.type == SCALAR)
      model_->InputPropertyValue((*it).first.c_str(), ScalarVariable(ANY, 0.));
    else if ((*it).second.key.type == VECTOR)
      model_->InputPropertyValue((*it).first.c_str(), VectorVariable<2U>(ANY, 0.));
    else if ((*it).second.key.type == TENSOR)
      model_->InputPropertyValue((*it).first.c_str(), TensorVariable<2U>(ANY, 0.));
  }
}


double CVFEM_HydrostaticColumn_Test::MaxDifference(const std::string& first, const std::string& second) const
{
  const Index first_key(property_database_->StorageKey(first.c_str()));
  const Index second_key(property_database_->StorageKey(second.c_str()));
  double maximum{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    maximum = std::max(maximum, std::fabs((*node_it)->Read(first_key) - (*node_it)->Read(second_key)));
  return maximum;
}


void CVFEM_HydrostaticColumn_Test::InitializeHydrostaticPressure()
{
  constexpr double tolerance{1.};   // [Pa]
  constexpr uint32_t max_iterations{100U};

  SolverBundle init_bundle(SolverKind::PETSc);
#if defined(CSMP_WITH_PETSC_SOLVER)
  init_bundle.PETScSettings().SetKSPType("gmres");
  init_bundle.PETScSettings().SetPCType("gamg");
#endif
  PDE_Integrator<2U> initial_pressure(init_bundle.Get());
  NumIntegral_dNT_op_dN_dV<2U> conductance(*property_database_, "permeability", "fluid pressure", "fluid pressure");
  NumIntegral_NT_op_dNi_dV<2U> buoyancy(*property_database_, "fluid density", "permeability", "fluid pressure");
  initial_pressure.Add(&conductance);
  initial_pressure.Add(&buoyancy);
  model_->Apply(initial_pressure);

  for (uint32_t iteration = 0U; iteration < max_iterations; ++iteration) {
    model_->CopyReplace("fluid pressure", "test fluid pressure");
    cvfem_phx_->InitialFluidPropertiesFromPTX();
    model_->Apply(initial_pressure);
    const double change{MaxDifference("fluid pressure", "test fluid pressure")};
    if (change < tolerance)
      break;
  }
}


double CVFEM_HydrostaticColumn_Test::MaxLiquidSpeed()
{
  const Index velocity_key(property_database_->StorageKey("velocity liquid"));
  double maximum{0.};
  double x_at_maximum{0.};
  double y_at_maximum{0.};
  for (auto element_it = model_->Region("Model").CellsBegin(); element_it != model_->Region("Model").CellsEnd();
       ++element_it) {
    VectorVariable<2U> velocity;
    (*element_it)->Read(velocity_key, velocity);
    if (velocity.Length() > maximum) {
      maximum = velocity.Length();
      x_at_maximum = (*element_it)->BaryCenter()[0];
      y_at_maximum = (*element_it)->BaryCenter()[1];
    }
  }
  location_of_max_speed_ = {x_at_maximum, y_at_maximum};
  return maximum;
}


std::array<double, 2> CVFEM_HydrostaticColumn_Test::MaxTransportSpeed(double run_time) const
{
  const Index flux_integral_key(property_database_->StorageKey("fluid flux integral"));
  const Index density_key(property_database_->StorageKey("fluid density"));
  const Index bulk_volume_key(property_database_->StorageKey("bulk volume"));
  const Index thickness_key(property_database_->StorageKey("thickness"));
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  std::array<double, 2> maximum{0., 0.};   // [interior nodes, nodes with a Dirichlet temperature or pressure]
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double thickness{(*node_it)->Read(thickness_key)};
    const double area{(*node_it)->Read(bulk_volume_key) / thickness};
    // mass that left the control volume, as a mean Darcy speed through half its perimeter (~ 2 sqrt(area))
    const double speed{(*node_it)->Read(flux_integral_key) / run_time
                       / ((*node_it)->Read(density_key) * 2. * std::sqrt(area) * thickness)};
    ScalarVariable temperature;
    ScalarVariable pressure;
    (*node_it)->Read(temperature_key, temperature);
    (*node_it)->Read(pressure_key, pressure);
    const bool dirichlet{temperature.Flag() == DIRICH || pressure.Flag() == DIRICH};
    maximum[dirichlet ? 1U : 0U] = std::max(maximum[dirichlet ? 1U : 0U], speed);
  }
  return maximum;
}


double CVFEM_HydrostaticColumn_Test::BuoyancyVelocityScale() const
{
  // bulk fluid density and the smallest viscosity of the phases present: near the critical point the equation of
  // state labels low-density supercritical fluid as vapour
  const Index density_key(property_database_->StorageKey("fluid density"));
  const Index viscosity_liquid_key(property_database_->StorageKey("viscosity liquid"));
  const Index viscosity_vapor_key(property_database_->StorageKey("viscosity vapor"));
  const Index permeability_key(property_database_->StorageKey("permeability"));
  const Region<2U>& region(model_->Region("Model"));
  double density_min{std::numeric_limits<double>::max()};
  double density_max{0.};
  double viscosity_min{std::numeric_limits<double>::max()};
  double permeability_max{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    density_min = std::min(density_min, (*node_it)->Read(density_key));
    density_max = std::max(density_max, (*node_it)->Read(density_key));
    for (const Index& key : {viscosity_liquid_key, viscosity_vapor_key})
      if ((*node_it)->Read(key) > 0.)
        viscosity_min = std::min(viscosity_min, (*node_it)->Read(key));
  }
  for (auto element_it = region.CellsBegin(); element_it != region.CellsEnd(); ++element_it)
    permeability_max = std::max(permeability_max, (*element_it)->Read(permeability_key));
  return permeability_max * gravity_ * (density_max - density_min) / viscosity_min;
}


void CVFEM_HydrostaticColumn_Test::InitializeRestingColumn(double largest_time_step)
{
  InitializeConductiveTemperature();
  InitializeAllNaNsToZero();
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  cvfem_phx_->SetLargestTimeStep(largest_time_step);
  InitializeHydrostaticPressure();
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(largest_time_step / 10.);
}


double CVFEM_HydrostaticColumn_Test::AdvanceAndTrackStoredVelocity(double max_time)
{
  double& model_time(ModelTime::Instance().modelTime);
  double max_speed{0.};
  steps_ = 0U;
  while (model_time < max_time) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    ++steps_;
    const double speed{MaxLiquidSpeed()};
    max_speed = std::max(max_speed, speed);
    if (steps_ <= 5U || steps_ % 50U == 0U)
      std::cerr << "CVFEM_HydrostaticColumn_Test::run: step " << steps_ << ", " << model_time / year_
                << " yr, max liquid Darcy speed " << speed << " m/s at (" << location_of_max_speed_[0] << ", "
                << location_of_max_speed_[1] << ")" << std::endl;
    if (dt < 1.)
      break;
  }
  return max_speed;
}


void CVFEM_HydrostaticColumn_Test::ReportResults(double velocity_scale, const std::array<double, 2>& transport_speed,
                                                 double max_stored_speed) const
{
  std::cerr << "CVFEM_HydrostaticColumn_Test::run: " << mesh_name_ << " (facet-mean gravity density "
            << (facet_mean_gravity_density_ ? "on" : "off") << "): " << steps_ << " steps to "
            << ModelTime::Instance().modelTime / year_ << " yr, buoyancy velocity scale " << velocity_scale << " m/s"
            << std::endl
            << "  transported mass as mean Darcy speed: interior nodes " << transport_speed[0] << " m/s (ratio "
            << transport_speed[0] / velocity_scale << ", limit " << max_velocity_fraction_ << "), Dirichlet nodes "
            << transport_speed[1] << " m/s (ratio " << transport_speed[1] / velocity_scale << ")" << std::endl
            << "  stored element 'velocity liquid' (information only): max " << max_stored_speed << " m/s (ratio "
            << max_stored_speed / velocity_scale << ")" << std::endl;
  // fluid state of the nodes as classified by the equation of state (States.h); near the critical point the
  // low-density part of a single supercritical fluid is classified V, the rest F
  constexpr std::array<const char*, 11> state_names{"none", "L", "V", "F", "VL", "VH", "LH", "VLH", "H", "M", "HM"};
  const Index state_key(property_database_->StorageKey("fluid state"));
  std::map<int, size_t> states;
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    ++states[static_cast<int>(std::lround((*node_it)->Read(state_key)))];
  std::cerr << "  fluid states (information only):";
  for (const auto& [state, count] : states)
    std::cerr << " " << (state >= 0 && state < static_cast<int>(state_names.size()) ? state_names[state] : "?") << " "
              << count;
  std::cerr << std::endl;
}


void CVFEM_HydrostaticColumn_Test::run()
{
  CreateModel();
  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};

  InitializeRestingColumn(run_settings_.TimeIncrement() * year_);
  const double velocity_scale{BuoyancyVelocityScale()};
  const double max_stored_speed{AdvanceAndTrackStoredVelocity(max_time)};
  const std::array<double, 2> transport_speed{MaxTransportSpeed(model_time)};
  ReportResults(velocity_scale, transport_speed, max_stored_speed);

  do_test(model_time >= max_time, "Hydrostatic column: ran the configured duration", __FILE__, __LINE__);
  do_test(transport_speed[0] < max_velocity_fraction_ * velocity_scale,
          "Hydrostatic column: mass transported between interior control volumes below the limit fraction of the "
          "buoyancy velocity scale", __FILE__, __LINE__);
}

} // namespace csmp
