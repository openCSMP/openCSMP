// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_SquareCavityConvection_Test.cpp
 *  @author Samuel Scott
 *  @brief Steady thermal convection in a closed porous square, see CVFEM_SquareCavityConvection_Test.h.
 *  @date 17.11.2026
 */

#include "CVFEM_SquareCavityConvection_Test.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <vector>

#include "VTU_Interface.h"

namespace csmp {

CVFEM_SquareCavityConvection_Test::CVFEM_SquareCavityConvection_Test(const std::string& mesh_name,
                                                                     const std::string& config_name, bool facet_mean_gravity_density,
                                                                     double reference_nusselt)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    facet_mean_gravity_density_(facet_mean_gravity_density),
    reference_nusselt_(reference_nusselt),
    nusselt_tolerance_(0.05),
    steady_tolerance_(2.e-3),
    temperature_base_(0.),
    temperature_top_(0.),
    y_base_(0.),
    y_top_(0.),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    run_settings_()
{
}


CVFEM_SquareCavityConvection_Test::~CVFEM_SquareCavityConvection_Test() = default;


void CVFEM_SquareCavityConvection_Test::CreateModel()
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


void CVFEM_SquareCavityConvection_Test::AnchorPressureAtOneNode()
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Region<2U>& region(model_->Region("Model"));
  Node<2U>* anchor{nullptr};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    if (std::fabs((*node_it)->y() - y_top_) < 1.e-6 && (anchor == nullptr || (*node_it)->x() < anchor->x()))
      anchor = *node_it;
  const double pressure{anchor->Read(pressure_key)};
  anchor->Store(pressure_key, ScalarVariable(DIRICH, pressure));
  std::cerr << "CVFEM_SquareCavityConvection_Test::AnchorPressureAtOneNode: fluid pressure fixed at " << pressure
            << " Pa on the node (" << anchor->x() << ", " << anchor->y() << ") only (closed top)" << std::endl;
}


void CVFEM_SquareCavityConvection_Test::InitializeTemperature()
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  double x_min{std::numeric_limits<double>::max()};
  double x_max{std::numeric_limits<double>::lowest()};
  y_base_ = std::numeric_limits<double>::max();
  y_top_ = std::numeric_limits<double>::lowest();
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    x_min = std::min(x_min, (*node_it)->x());
    x_max = std::max(x_max, (*node_it)->x());
    y_base_ = std::min(y_base_, (*node_it)->y());
    y_top_ = std::max(y_top_, (*node_it)->y());
  }
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    if (std::fabs((*node_it)->y() - y_base_) < 1.e-6) temperature_base_ = (*node_it)->Read(temperature_key);
    if (std::fabs((*node_it)->y() - y_top_) < 1.e-6) temperature_top_ = (*node_it)->Read(temperature_key);
  }
  constexpr double perturbation{0.01};   // fraction of the temperature difference
  const double height{y_top_ - y_base_};
  const double width{x_max - x_min};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    if (temperature.Flag() == DIRICH)
      continue;
    const double z{((*node_it)->y() - y_base_) / height};
    const double x{((*node_it)->x() - x_min) / width};
    const double value{temperature_base_ + z * (temperature_top_ - temperature_base_)
                       + perturbation * (temperature_base_ - temperature_top_)
                         * std::sin(std::numbers::pi * z) * std::cos(std::numbers::pi * x)};
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), value));
  }
  std::cerr << "CVFEM_SquareCavityConvection_Test::InitializeTemperature: " << temperature_base_ << " C (base) to "
            << temperature_top_ << " C (top) over " << height << " m, single-cell perturbation" << std::endl;
}


void CVFEM_SquareCavityConvection_Test::InitializeAllNaNsToZero()
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


double CVFEM_SquareCavityConvection_Test::MaxDifference(const std::string& first, const std::string& second) const
{
  const Index first_key(property_database_->StorageKey(first.c_str()));
  const Index second_key(property_database_->StorageKey(second.c_str()));
  double maximum{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    maximum = std::max(maximum, std::fabs((*node_it)->Read(first_key) - (*node_it)->Read(second_key)));
  return maximum;
}


void CVFEM_SquareCavityConvection_Test::InitializeHydrostaticPressure()
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
    if (MaxDifference("fluid pressure", "test fluid pressure") < tolerance)
      break;
  }
}


double CVFEM_SquareCavityConvection_Test::TopNusseltNumber() const
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  // node rows keyed by height (rounded to mm), each a map x -> T
  std::map<long, std::map<long, double>> rows;
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    rows[std::lround((*node_it)->y() * 1000.)][std::lround((*node_it)->x() * 1000.)] = (*node_it)->Read(temperature_key);
  const auto top{rows.rbegin()};
  const auto second{std::next(top)};
  const auto third{std::next(second)};
  const double dz{(top->first - second->first) / 1000.};
  if (second->second.size() != top->second.size() || third->second.size() != top->second.size()
      || std::fabs((second->first - third->first) / 1000. - dz) > 1.e-6)
    std::cerr << "CVFEM_SquareCavityConvection_Test::TopNusseltNumber: WARNING: mesh is not structured" << std::endl;

  // trapezoidal average over x of -dT/dz at the top, second-order one-sided difference
  double integral{0.};
  double width{0.};
  double previous_x{0.};
  double previous_gradient{0.};
  bool first{true};
  for (const auto& [x_key, t_top] : top->second) {
    const double gradient{-(3. * t_top - 4. * second->second.at(x_key) + third->second.at(x_key)) / (2. * dz)};
    const double x{x_key / 1000.};
    if (!first) {
      integral += 0.5 * (gradient + previous_gradient) * (x - previous_x);
      width += x - previous_x;
    }
    previous_x = x;
    previous_gradient = gradient;
    first = false;
  }
  return integral / width * (y_top_ - y_base_) / (temperature_base_ - temperature_top_);
}


double CVFEM_SquareCavityConvection_Test::RayleighNumber()
{
  // S. Scott, 2026-11-18: drho is the thermal density contrast at constant pressure (before: the top-to-base
  // difference, which includes the hydrostatic compression of the water and does not drive flow)
  constexpr double increment{0.1};   // [K]
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Index viscosity_key(property_database_->StorageKey("viscosity liquid"));
  const Index heat_capacity_key(property_database_->StorageKey("fluid heat capacity"));
  const Index permeability_key(property_database_->StorageKey("permeability"));
  const Index conductivity_key(property_database_->StorageKey("thermal conductivity"));
  const Region<2U>& region(model_->Region("Model"));
  std::vector<double> density_before;
  double viscosity{0.};
  double density{0.};
  double heat_capacity{0.};
  double nodes{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    density_before.push_back((*node_it)->Read(density_key));
    viscosity += (*node_it)->Read(viscosity_key);
    density += (*node_it)->Read(density_key);
    heat_capacity += (*node_it)->Read(heat_capacity_key);
    nodes += 1.;
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), temperature() + increment));
  }
  cvfem_phx_->InitialFluidPropertiesFromPTX();
  double density_derivative{0.};   // -d rho / d T at constant pressure, node mean [kg m-3 K-1]
  std::size_t n{0U};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it, ++n) {
    density_derivative += (density_before[n] - (*node_it)->Read(density_key)) / increment / nodes;
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), temperature() - increment));
  }
  cvfem_phx_->InitialFluidPropertiesFromPTX();
  double permeability{0.};
  double conductivity{0.};
  double cells{0.};
  for (auto element_it = region.CellsBegin(); element_it != region.CellsEnd(); ++element_it) {
    permeability += (*element_it)->Read(permeability_key);
    conductivity += (*element_it)->Read(conductivity_key);
    cells += 1.;
  }
  const double kappa{(conductivity / cells) / ((density / nodes) * (heat_capacity / nodes))};
  return density_derivative * (temperature_base_ - temperature_top_) * gravity_ * (y_top_ - y_base_)
         * (permeability / cells) / ((viscosity / nodes) * kappa);
}


void CVFEM_SquareCavityConvection_Test::run()
{
  CreateModel();

  double& model_time(ModelTime::Instance().modelTime);
  model_time = 0.;
  const double max_time{run_settings_.Duration() * year_};
  const double largest_time_step{run_settings_.TimeIncrement() * year_};

  InitializeTemperature();
  AnchorPressureAtOneNode();
  InitializeAllNaNsToZero();
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  cvfem_phx_->SetLargestTimeStep(largest_time_step);
  InitializeHydrostaticPressure();
  const double rayleigh{RayleighNumber()};   // conductive initial state
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(largest_time_step / 10.);

  // fields at every checkpoint, for inspection: <config base name>_Model_(REGION)_<checkpoint>.vtu in the working directory
  VTU_Interface<2U> vtu(*model_);
  vtu.OmitZeroInFileName(false);
  const std::vector<std::string> output_variables{"temperature", "fluid pressure", "density liquid", "fluid density",
                                                  "velocity liquid", "fluid flux"};
  const std::string output_name{config_name_.substr(config_name_.find_last_of('/') + 1U)};
  vtu.OutputDataToVTU(output_name.c_str(), output_variables, model_->Region("Model"), 0);

  // Nusselt number at 50 checkpoints; steady if it changed by less than steady_tolerance_ over the last 5
  constexpr uint32_t checkpoints{50U};
  std::vector<double> nusselt;
  double next_checkpoint{max_time / checkpoints};
  uint64_t steps{0U};
  while (model_time < max_time) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    ++steps;
    if (model_time >= next_checkpoint) {
      nusselt.push_back(TopNusseltNumber());
      vtu.OutputDataToVTU(output_name.c_str(), output_variables, model_->Region("Model"),
                          static_cast<long>(nusselt.size()));
      std::cerr << "CVFEM_SquareCavityConvection_Test::run: " << model_time / year_ << " yr, step " << steps
                << ", Nu " << nusselt.back() << std::endl;
      next_checkpoint += max_time / checkpoints;
    }
    if (dt < 1.)
      break;
  }

  const double nusselt_final{nusselt.empty() ? 0. : nusselt.back()};
  const double change{nusselt.size() > 5U ? std::fabs(nusselt.back() - nusselt[nusselt.size() - 6U]) / nusselt.back()
                                          : 1.};
  std::cerr << "CVFEM_SquareCavityConvection_Test::run: " << config_name_ << " (facet-mean gravity density "
            << (facet_mean_gravity_density_ ? "on" : "off") << "): " << steps << " steps to "
            << model_time / year_ << " yr, Ra " << rayleigh << ", Nu " << nusselt_final << " (reference "
            << reference_nusselt_ << ", deviation " << (nusselt_final - reference_nusselt_) / reference_nusselt_
            << ", limit " << nusselt_tolerance_ << "), relative change over the last 10 % of the run " << change
            << std::endl;
  do_test(model_time >= max_time, "Square cavity: ran the configured duration", __FILE__, __LINE__);
  do_test(change < steady_tolerance_, "Square cavity: Nusselt number steady", __FILE__, __LINE__);
  do_test(std::fabs(nusselt_final - reference_nusselt_) < nusselt_tolerance_ * reference_nusselt_,
          "Square cavity: Nusselt number within the tolerance of the published value", __FILE__, __LINE__);
}

} // namespace csmp
