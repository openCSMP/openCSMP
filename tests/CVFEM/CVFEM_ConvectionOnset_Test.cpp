// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

/** @file CVFEM_ConvectionOnset_Test.cpp
 *  @author Samuel Scott
 *  @brief Growth or decay rate of a convective mode in a closed porous box, see CVFEM_ConvectionOnset_Test.h.
 *  @date 18.11.2026
 */

#include "CVFEM_ConvectionOnset_Test.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <vector>

#include "VTU_Interface.h"

namespace csmp {

CVFEM_ConvectionOnset_Test::CVFEM_ConvectionOnset_Test(const std::string& mesh_name, const std::string& config_name,
                                                       bool facet_mean_gravity_density, bool finite_amplitude)
  : mesh_name_(mesh_name),
    config_name_(config_name),
    vars_name_("PhysicalVariables.txt"),   // created by the ModelBuilder
    facet_mean_gravity_density_(facet_mean_gravity_density),
    finite_amplitude_(finite_amplitude),
    rayleigh_tolerance_(0.03),
    amplitude_(0.),
    rayleigh_(0.),
    sigma_theory_(0.),
    time_unit_(0.),
    thermal_density_contrast_(0.),
    max_time_(0.),
    largest_time_step_(0.),
    ran_planned_duration_(true),
    temperature_base_(0.),
    temperature_top_(0.),
    x_min_(0.),
    y_base_(0.),
    y_top_(0.),
    conductivity_(0.),
    bulk_heat_capacity_(0.),
    fluid_heat_capacity_(0.),
    viscosity_(0.),
    permeability_(0.),
    model_(nullptr),
    property_database_(nullptr),
    cvfem_phx_(nullptr),
    permeability_visitor_(nullptr),
    run_settings_()
{
}


CVFEM_ConvectionOnset_Test::~CVFEM_ConvectionOnset_Test() = default;


void CVFEM_ConvectionOnset_Test::CreateModel()
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


void CVFEM_ConvectionOnset_Test::ReadGeometryAndBoundaryTemperatures()
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  x_min_ = std::numeric_limits<double>::max();
  y_base_ = std::numeric_limits<double>::max();
  y_top_ = std::numeric_limits<double>::lowest();
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    x_min_ = std::min(x_min_, (*node_it)->x());
    y_base_ = std::min(y_base_, (*node_it)->y());
    y_top_ = std::max(y_top_, (*node_it)->y());
  }
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    if (std::fabs((*node_it)->y() - y_base_) < 1.e-6) temperature_base_ = (*node_it)->Read(temperature_key);
    if (std::fabs((*node_it)->y() - y_top_) < 1.e-6) temperature_top_ = (*node_it)->Read(temperature_key);
  }
}


void CVFEM_ConvectionOnset_Test::AnchorPressureAtOneNode()
{
  const Index pressure_key(property_database_->StorageKey("fluid pressure"));
  const Region<2U>& region(model_->Region("Model"));
  Node<2U>* anchor{nullptr};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    if (std::fabs((*node_it)->y() - y_top_) < 1.e-6 && (anchor == nullptr || (*node_it)->x() < anchor->x()))
      anchor = *node_it;
  const double pressure{anchor->Read(pressure_key)};
  anchor->Store(pressure_key, ScalarVariable(DIRICH, pressure));
  std::cerr << "CVFEM_ConvectionOnset_Test::AnchorPressureAtOneNode: fluid pressure fixed at " << pressure
            << " Pa on the node (" << anchor->x() << ", " << anchor->y() << ") only (closed top)" << std::endl;
}


void CVFEM_ConvectionOnset_Test::InitializeTemperature(double amplitude)
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Region<2U>& region(model_->Region("Model"));
  const double height{y_top_ - y_base_};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    if (temperature.Flag() == DIRICH)
      continue;
    const double z{((*node_it)->y() - y_base_) / height};
    const double x{((*node_it)->x() - x_min_) / height};
    const double value{temperature_base_ + z * (temperature_top_ - temperature_base_)
                       + amplitude * std::cos(std::numbers::pi * x) * std::sin(std::numbers::pi * z)};
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), value));
  }
  std::cerr << "CVFEM_ConvectionOnset_Test::InitializeTemperature: " << temperature_base_ << " C (base) to "
            << temperature_top_ << " C (top) over " << height << " m, mode amplitude " << amplitude << " K"
            << std::endl;
}


void CVFEM_ConvectionOnset_Test::InitializeAllNaNsToZero()
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


double CVFEM_ConvectionOnset_Test::MaxDifference(const std::string& first, const std::string& second) const
{
  const Index first_key(property_database_->StorageKey(first.c_str()));
  const Index second_key(property_database_->StorageKey(second.c_str()));
  double maximum{0.};
  const Region<2U>& region(model_->Region("Model"));
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it)
    maximum = std::max(maximum, std::fabs((*node_it)->Read(first_key) - (*node_it)->Read(second_key)));
  return maximum;
}


void CVFEM_ConvectionOnset_Test::InitializeHydrostaticPressure()
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


double CVFEM_ConvectionOnset_Test::ThermalDensityDerivative()
{
  constexpr double increment{0.1};   // [K]
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Region<2U>& region(model_->Region("Model"));
  std::vector<double> density;
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    density.push_back((*node_it)->Read(density_key));
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), temperature() + increment));
  }
  cvfem_phx_->InitialFluidPropertiesFromPTX();
  double sum{0.};
  std::size_t n{0U};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it, ++n) {
    sum += (density[n] - (*node_it)->Read(density_key)) / increment;
    ScalarVariable temperature;
    (*node_it)->Read(temperature_key, temperature);
    (*node_it)->Store(temperature_key, ScalarVariable(temperature.Flag(), temperature() - increment));
  }
  cvfem_phx_->InitialFluidPropertiesFromPTX();
  return sum / static_cast<double>(n);
}


void CVFEM_ConvectionOnset_Test::ReadBulkProperties()
{
  const Index porosity_key(property_database_->StorageKey("nodal porosity"));
  const Index rock_density_key(property_database_->StorageKey("nodal density rock"));
  const Index rock_heat_capacity_key(property_database_->StorageKey("nodal heat capacity rock"));
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Index heat_capacity_key(property_database_->StorageKey("fluid heat capacity"));
  const Index viscosity_key(property_database_->StorageKey("viscosity liquid"));
  const Index permeability_key(property_database_->StorageKey("permeability"));
  const Index conductivity_key(property_database_->StorageKey("thermal conductivity"));
  const Region<2U>& region(model_->Region("Model"));
  fluid_heat_capacity_ = bulk_heat_capacity_ = viscosity_ = permeability_ = conductivity_ = 0.;
  double nodes{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double porosity{(*node_it)->Read(porosity_key)};
    const double fluid{(*node_it)->Read(density_key) * (*node_it)->Read(heat_capacity_key)};
    fluid_heat_capacity_ += fluid;
    bulk_heat_capacity_ += (1. - porosity) * (*node_it)->Read(rock_density_key)
                               * (*node_it)->Read(rock_heat_capacity_key)
                           + porosity * fluid;
    viscosity_ += (*node_it)->Read(viscosity_key);
    nodes += 1.;
  }
  fluid_heat_capacity_ /= nodes;
  bulk_heat_capacity_ /= nodes;
  viscosity_ /= nodes;
  double cells{0.};
  for (auto element_it = region.CellsBegin(); element_it != region.CellsEnd(); ++element_it) {
    permeability_ += (*element_it)->Read(permeability_key);
    conductivity_ += (*element_it)->Read(conductivity_key);
    cells += 1.;
  }
  permeability_ /= cells;
  conductivity_ /= cells;
}


double CVFEM_ConvectionOnset_Test::ModeAmplitude() const
{
  const Index temperature_key(property_database_->StorageKey("temperature"));
  const Index volume_key(property_database_->StorageKey("bulk volume"));
  const Region<2U>& region(model_->Region("Model"));
  const double height{y_top_ - y_base_};
  double projection{0.};
  double norm{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    const double z{((*node_it)->y() - y_base_) / height};
    const double x{((*node_it)->x() - x_min_) / height};
    const double mode{std::cos(std::numbers::pi * x) * std::sin(std::numbers::pi * z)};
    const double anomaly{(*node_it)->Read(temperature_key) - temperature_base_
                         - z * (temperature_top_ - temperature_base_)};
    const double volume{(*node_it)->Read(volume_key)};
    projection += volume * anomaly * mode;
    norm += volume * mode * mode;
  }
  return projection / norm;
}


void CVFEM_ConvectionOnset_Test::RebuildModel()
{
  cvfem_phx_.reset();
  permeability_visitor_.reset();
  model_.reset();
  CreateModel();
  ModelTime::Instance().modelTime = 0.;
}


void CVFEM_ConvectionOnset_Test::InitializeConductiveState()
{
  ReadGeometryAndBoundaryTemperatures();
  InitializeTemperature(0.);
  AnchorPressureAtOneNode();
  InitializeAllNaNsToZero();
  cvfem_phx_->Adjust_CFL_Criterion(model_->Read(property_database_->StorageKey("cfl scaling")), true);
  InitializeHydrostaticPressure();
}


void CVFEM_ConvectionOnset_Test::ComputeRayleighNumberAndTimeUnit()
{
  ReadBulkProperties();
  const double height{y_top_ - y_base_};
  thermal_density_contrast_ = ThermalDensityDerivative() * (temperature_base_ - temperature_top_);
  rayleigh_ = thermal_density_contrast_ * gravity_ * height * permeability_ * fluid_heat_capacity_
              / (viscosity_ * conductivity_);
  sigma_theory_ = 0.5 * rayleigh_ - 2. * std::numbers::pi * std::numbers::pi;
  time_unit_ = bulk_heat_capacity_ * height * height / conductivity_;
}


void CVFEM_ConvectionOnset_Test::ChooseRunLengthAndTimeStep()
{
  const double rate{std::max(std::fabs(sigma_theory_), 1.)};
  max_time_ = std::min(run_settings_.Duration() * year_, growth_span_ / rate * time_unit_);
  largest_time_step_ = std::min(run_settings_.TimeIncrement() * year_, time_step_fraction_ / rate * time_unit_);
  // initial amplitude: the larger of the initial and final amplitude is peak_relative_amplitude_ * dT (linear regime)
  const double temperature_difference{temperature_base_ - temperature_top_};
  amplitude_ = finite_amplitude_ ? finite_relative_amplitude_ * temperature_difference
                                 : peak_relative_amplitude_ * temperature_difference
                                   * (sigma_theory_ > 0. ? std::exp(-growth_span_) : 1.);
}


void CVFEM_ConvectionOnset_Test::ImposeMode(double amplitude)
{
  InitializeTemperature(amplitude);
  InitializeHydrostaticPressure();
}


double CVFEM_ConvectionOnset_Test::TopToBaseDensityDifference() const
{
  const Index density_key(property_database_->StorageKey("density liquid"));
  const Region<2U>& region(model_->Region("Model"));
  double density_top{0.};
  double density_base{0.};
  double nodes_top{0.};
  double nodes_base{0.};
  for (auto node_it = region.NodesBegin(); node_it != region.NodesEnd(); ++node_it) {
    if (std::fabs((*node_it)->y() - y_top_) < 1.e-6) {
      density_top += (*node_it)->Read(density_key);
      nodes_top += 1.;
    }
    if (std::fabs((*node_it)->y() - y_base_) < 1.e-6) {
      density_base += (*node_it)->Read(density_key);
      nodes_base += 1.;
    }
  }
  return density_top / nodes_top - density_base / nodes_base;
}


void CVFEM_ConvectionOnset_Test::ReportSetUp() const
{
  std::cerr << "CVFEM_ConvectionOnset_Test::RunCase: Ra " << rayleigh_ << " (thermal density contrast "
            << thermal_density_contrast_ << " kg m-3; top-to-base " << TopToBaseDensityDifference()
            << " kg m-3), sigma (theory) " << sigma_theory_ << ", time unit " << time_unit_ / year_ << " yr, run to "
            << max_time_ / year_ << " yr, largest time step " << largest_time_step_ / year_ << " yr" << std::endl;
}


void CVFEM_ConvectionOnset_Test::StartTransient()
{
  cvfem_phx_->SetLargestTimeStep(largest_time_step_);
  cvfem_phx_->PrepareTransientCalculations();
  cvfem_phx_->ChangeTimeStepTo(largest_time_step_ / 10.);
}


CVFEM_ConvectionOnset_Test::ModeHistory CVFEM_ConvectionOnset_Test::AdvanceAndRecordMode(bool write_fields)
{
  // fields at 40 checkpoints, for inspection: <config base name>_Model_(REGION)_<checkpoint>.vtu in the working
  // directory
  VTU_Interface<2U> vtu(*model_);
  vtu.OmitZeroInFileName(false);
  const std::vector<std::string> output_variables{"temperature", "fluid pressure", "density liquid", "velocity liquid",
                                                  "fluid flux"};
  const std::string output_name{config_name_.substr(config_name_.find_last_of('/') + 1U)};
  if (write_fields)
    vtu.OutputDataToVTU(output_name.c_str(), output_variables, model_->Region("Model"), 0);

  constexpr uint32_t checkpoints{40U};
  uint32_t checkpoint{0U};
  double next_checkpoint{max_time_ / checkpoints};
  double& model_time(ModelTime::Instance().modelTime);
  ModeHistory history{{0.}, {ModeAmplitude()}};
  while (model_time < max_time_) {
    const double dt{cvfem_phx_->Apply()};
    model_time += dt;
    model_->InputPropertyValue("dt", ScalarVariable(PLAIN, dt));
    history.tau.push_back(model_time / time_unit_);
    history.amplitude.push_back(ModeAmplitude());
    if (write_fields && model_time >= next_checkpoint) {
      vtu.OutputDataToVTU(output_name.c_str(), output_variables, model_->Region("Model"),
                          static_cast<long>(++checkpoint));
      next_checkpoint += max_time_ / checkpoints;
    }
    if (dt < 1.)
      break;
  }
  ran_planned_duration_ = ran_planned_duration_ && model_time >= max_time_;
  std::cerr << "CVFEM_ConvectionOnset_Test::RunCase: initial mode amplitude " << history.amplitude.front() << " K, "
            << history.tau.size() - 1U << " steps to " << model_time / year_ << " yr, final amplitude "
            << history.amplitude.back() << " K" << std::endl;
  return history;
}


CVFEM_ConvectionOnset_Test::ModeHistory CVFEM_ConvectionOnset_Test::RunCase(double sign, bool write_fields)
{
  RebuildModel();
  InitializeConductiveState();
  ComputeRayleighNumberAndTimeUnit();
  ChooseRunLengthAndTimeStep();
  ImposeMode(sign * amplitude_);
  if (write_fields)
    ReportSetUp();
  StartTransient();
  return AdvanceAndRecordMode(write_fields);
}


double CVFEM_ConvectionOnset_Test::FitLinearGrowthRate(const ModeHistory& plus, const ModeHistory& minus) const
{
  // straight-line fit of ln A_lin against tau over the last 75 % of the run, A_lin = (A+ - A-) / 2
  const std::size_t samples_common{std::min(plus.tau.size(), minus.tau.size())};
  const double tau_end{plus.tau[samples_common - 1U]};
  double sum_tau{0.};
  double sum_log{0.};
  double sum_tau_tau{0.};
  double sum_tau_log{0.};
  double samples{0.};
  for (std::size_t i = 0U; i < samples_common; ++i) {
    const double linear_part{0.5 * (plus.amplitude[i] - minus.amplitude[i])};
    if (plus.tau[i] < 0.25 * tau_end || linear_part <= 0.)
      continue;
    const double log_amplitude{std::log(linear_part)};
    sum_tau += plus.tau[i];
    sum_log += log_amplitude;
    sum_tau_tau += plus.tau[i] * plus.tau[i];
    sum_tau_log += plus.tau[i] * log_amplitude;
    samples += 1.;
  }
  return (samples * sum_tau_log - sum_tau * sum_log) / (samples * sum_tau_tau - sum_tau * sum_tau);
}


double CVFEM_ConvectionOnset_Test::LargestSignIndependentPart(const ModeHistory& plus, const ModeHistory& minus) const
{
  const std::size_t samples_common{std::min(plus.amplitude.size(), minus.amplitude.size())};
  double largest{0.};
  for (std::size_t i = 0U; i < samples_common; ++i)
    largest = std::max(largest, std::fabs(0.5 * (plus.amplitude[i] + minus.amplitude[i])) / amplitude_);
  return largest;
}


bool CVFEM_ConvectionOnset_Test::LinearPartStaysPositive(const ModeHistory& plus, const ModeHistory& minus) const
{
  const std::size_t samples_common{std::min(plus.amplitude.size(), minus.amplitude.size())};
  for (std::size_t i = 0U; i < samples_common; ++i)
    if (plus.amplitude[i] - minus.amplitude[i] <= 0.)
      return false;
  return true;
}


void CVFEM_ConvectionOnset_Test::CheckGrowthRate(const ModeHistory& plus, const ModeHistory& minus)
{
  const double sigma{FitLinearGrowthRate(plus, minus)};
  const double rayleigh_critical{4. * std::numbers::pi * std::numbers::pi};
  const double rayleigh_from_rate{2. * (sigma + 2. * std::numbers::pi * std::numbers::pi)};
  const double deviation{(rayleigh_from_rate - rayleigh_) / rayleigh_critical};
  std::cerr << "CVFEM_ConvectionOnset_Test::run: " << config_name_ << " (facet-mean gravity density "
            << (facet_mean_gravity_density_ ? "on" : "off") << "): Ra " << rayleigh_ << ", sigma " << sigma
            << " (theory " << sigma_theory_ << "), Ra from the rate " << rayleigh_from_rate << " (deviation "
            << deviation << " of Ra_c, limit " << rayleigh_tolerance_ << "); sign-independent part of the mode up to "
            << LargestSignIndependentPart(plus, minus) << " of the initial amplitude" << std::endl;
  do_test(LinearPartStaysPositive(plus, minus), "Convection onset: linear mode amplitude stays positive", __FILE__,
          __LINE__);
  do_test(std::fabs(deviation) < rayleigh_tolerance_,
          "Convection onset: growth rate within the tolerance of linear stability theory", __FILE__, __LINE__);
}


void CVFEM_ConvectionOnset_Test::CheckFiniteAmplitudeDecay(const ModeHistory& plus, const ModeHistory& minus)
{
  // below onset a perturbation of any size must decay (energy stability = linear stability for this problem)
  const double decay_plus{plus.amplitude.front() / std::fabs(plus.amplitude.back())};
  const double decay_minus{-minus.amplitude.front() / std::fabs(minus.amplitude.back())};
  std::cerr << "CVFEM_ConvectionOnset_Test::run: " << config_name_ << " (facet-mean gravity density "
            << (facet_mean_gravity_density_ ? "on" : "off") << "): Ra " << rayleigh_ << ", finite amplitude "
            << amplitude_ << " K decayed by factors " << decay_plus << " (+) and " << decay_minus << " (-), required "
            << decay_factor_ << " (linear theory " << std::exp(growth_span_)
            << "); sign-independent part of the mode up to " << LargestSignIndependentPart(plus, minus)
            << " of the initial amplitude" << std::endl;
  do_test(sigma_theory_ < 0. && decay_plus > decay_factor_ && decay_minus > decay_factor_,
          "Convection onset: finite perturbation below onset decays", __FILE__, __LINE__);
}


void CVFEM_ConvectionOnset_Test::run()
{
  // the same case with the mode added and subtracted: half the difference is the linear response of the mode, half
  // the sum the part that does not depend on the sign (any source of the mode, e.g. spurious buoyancy flow, and the
  // quadratic nonlinearity)
  const ModeHistory plus{RunCase(1., true)};
  const ModeHistory minus{RunCase(-1., false)};
  do_test(ran_planned_duration_, "Convection onset: ran the planned duration", __FILE__, __LINE__);
  if (finite_amplitude_)
    CheckFiniteAmplitudeDecay(plus, minus);
  else
    CheckGrowthRate(plus, minus);
}

} // namespace csmp
