#ifndef CVFEM_1D_VVCASE_H
#define CVFEM_1D_VVCASE_H


#include <chrono>
#include <sstream>
#include <vector>
#include "iostream"

// model
#include "Model1D.h"
#include "ModelTime.h"
#include "ModelBuilder.h"

// testing
#include "Test.h"

// File I/O and Initialization
#include "ComputationalSettings.h"
#include "InputDataManager.h"

// finite elements
#include "NumIntegral_dNT_op_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "NumIntegral_SetRHS_to_Zero.h"

// finite volumes
// CVFEM PHX
#include "CVFEM_MathOperatorRHS.h"
#include "CVFEM_PHX_Scheme.h"
#include "WellConfigurationFile.h"

// visitors
#include "ConductivityVisitor.h"
#include "ComputeGravityTermVisitor.h"
#include "PermeabilityVisitor.h"
#include "PoreVolumeVisitor.h"

// eos

// Solver
#include "LUdcmp_Solver.h"
#include "EigenSolver.h"

// Output
#include "VTU_Interface.h"
#include <time.h>

// compare values
#include "compareFloats.h"


using namespace std;
using namespace csmp;


namespace csmp
{

/** CVFEM_1D_VVCase : 1D Test Case
  =================================
  Mesh:       Not needed
  Model:      1D
  Test:       PHX_Scheme
  BC:         TBC
  Criterion:
  =================================
  */
class CVFEM_1D_VVCase : public Test
{
public:
    CVFEM_1D_VVCase ( const char* argv , const uint32_t test_type);
    ~CVFEM_1D_VVCase();
    void run();
private:

    const char* geometry_name_;
    const char* config_file_name_;
    const char* vars_name_;
    const char* system_file_list_name_;
    const char* recipes_file_list_name_;
    list<string> output_props_;

    const uint32_t test_type;

    enum{dim=1};

    Model<dim>* model;
    PropertyDatabase<dim>*  pd_ref;
    ModelBuilder<dim> builder;
    VTU_Interface<dim>* vtu;
    ComputationalSettings  run_settings;

    //! CVFEM instance
    CVFEM_PHX_Scheme<dim>* CVFEM_PHX;

    //! Permeability Visitor
    PermeabilityVisitor<dim>* permeability_visitor;

    //! Chemical diffusion operators and variables
    PDE_Integrator<dim>              diffusion;

    Index pore_diff_key, eff_diff_key, porosity_key;
    ScalarVariable diff_coeff, porosity;

    //!time variables
    double dt, max_time, vtu_dt;
    size_t timestep, output_counter;
    //! S. Scott, 2026-11-18: tests 105a-c compare with the error-function solution of heat conduction from the heated
    //! LEFT end instead of benchmark data, and always run to the end; 105b/c put the fluid in the conduction heat capacity
    bool conduction_erfc_test_;
    bool fluid_heat_capacity_;
    //! largest deviation of the temperature from T_right + (T_left - T_right) erfc(x / (2 sqrt(kappa t))) [K]
    double ConductionErfcDeviation() const;

    /// S. Scott, 2026-11-20: tests 106, 107 and 108, latent heat of crystallisation in the conduction-equilibration
    /// coupling. All three are test 105b (conduction only, rock + fluid conduction heat capacity) with rock melt that
    /// crystallises: rock liquidus/solidus on, power-law crystallisation curve, melt heat capacity = rock heat capacity
    /// (as Turcotte & Schubert, Geodynamics, 2nd ed., sec. 4.18-4.19). They always run to the end.
    /// test-106 (T&S sec. 4.18, the Stefan problem): melt at its liquidus T_m, LEFT end held at T_0; front
    ///   y_m = 2 lambda1 sqrt(kappa t) (eq. 4-136), solid T = T_0 + (T_m - T_0) erf(eta) / erf(lambda1) (eq. 4-137),
    ///   lambda1 from eq. 4-141.
    /// test-107 (T&S sec. 4.19, solidification of a dike): melt at T_m for distances below dike_half_width_ from the
    ///   (insulated) LEFT end, country rock at T_0 beyond, RIGHT end held at T_0; front moves into the melt by
    ///   2 lambda2 sqrt(kappa t) (eq. 4-142), country rock T = T_0 + (T_m - T_0) erfc(eta) / (1 + erf lambda2)
    ///   (eq. 4-145), constant contact temperature (eq. 4-146), lambda2 from eq. 4-149.
    /// test-108: wide (basalt-like) interval, set-up of 106; energy balance and temperatures against an enthalpy-method
    ///   reference computed in the test.
    /// Each also writes plot-ready CSV files next to its configuration file: <test>_front.csv (106/107, every 10 years:
    /// front position and, 107, contact temperature, model and analytical with and without latent heat),
    /// <test>_energy.csv (108: heat conducted out and change of stored energy) and <test>_profiles.csv (temperature
    /// profiles at a quarter, half, three quarters and all of the run, model and analytical or reference, with and
    /// without latent heat).
    /// Pass limits (end of the run, 400 years; 10 m mesh, 1-year steps). A 1-D emulation of the same lumped, implicit
    /// conduction step with the equilibration crediting the conducted heat gives 3.1 % (106 front), 0.07 % (107
    /// travel), 0.3 K / 0.02 K (106 solid layer / 107 contact), 1.2 K and 0.08 % (108); without latent heat in the
    /// coupling the errors are 186 %, 296 %, 170 K / 114 K, 184 K and 104 %.
    bool stefan_test_;
    bool dike_test_;
    bool wide_interval_test_;
    double solidus_;              ///< [C]
    double liquidus_;             ///< [C] = initial melt temperature T_m
    double latent_heat_;          ///< [J kg-1]
    double exponent_;             ///< melt fraction X = ((T - T_s) / (T_l - T_s))^exponent_ in the interval
    double cold_temperature_;     ///< T_0 [C]
    double diffusivity_;          ///< lambda / ((1 - phi) rho_r c_r + phi rho_f c_f) at t = 0 [m2 s-1]
    double front_parameter_;      ///< lambda1 (106) or lambda2 (107)
    double boundary_heat_;        ///< heat conducted into the fixed-temperature end nodes since t = 0 [J m-2]
    double initial_energy_;       ///< stored energy at t = 0 [J m-2]
    double next_record_time_;     ///< [s]
    double next_profile_time_;    ///< [s]
    std::ostringstream front_rows_;
    std::ostringstream profile_rows_;
    std::vector<double> profile_times_;                      ///< test-108: times of the stored model profiles [s]
    std::vector<std::vector<double>> profile_temperatures_;  ///< test-108: model profiles, sorted by distance

    /// Temperatures and distances from the LEFT end, sorted by distance.
    struct Profile {
      std::vector<double> distance;
      std::vector<double> temperature;
    };

    /// Rock liquidus/solidus on, power-law crystallisation curve with the configuration's latent heat and exponent.
    void ConfigureLatentHeat();
    /// test-107: melt at the liquidus for distances below dike_half_width_ (initial state, before the fluid properties).
    void InitializeDikeTemperature();
    /// T_0, kappa and lambda1 / lambda2 from the initial state.
    void ComputeAnalyticalParameters();
    /// Melt fraction of the power-law curve, written out here independently of Rock.cpp.
    double MeltFraction(double temperature) const;
    /// Specific rock enthalpy, c_r T + L X(T), independently of Rock.cpp [J kg-1].
    double RockEnthalpy(double temperature) const;
    /// (1 - phi) rho_r c_r + phi rho_f c_f, node mean [J m-3 K-1].
    double SensibleHeatCapacity() const;
    /// Mean thermal conductivity of the cells [W m-1 K-1].
    double MeanConductivity() const;
    /// (1 - phi) rho_r, node mean [kg m-3].
    double RockMassPerVolume() const;
    /// Sum over the control volumes of (1 - phi) rho_r h_r(T) + phi rho_f h_f [J per m2 of cross-section].
    double StoredEnergy() const;
    /// Heat flow from the free nodes into the fixed-temperature end nodes, lambda (T_1 - T_0) / dx at the end of the
    /// step (implicit conduction) [W m-2].
    double BoundaryHeatFlow() const;
    Profile SortedProfile() const;
    /// Linear interpolation of the profile at a distance.
    static double Interpolate(const Profile& profile, double distance);
    /// Distance of the solidus crossing, searched from the cold end (106: LEFT, 107: RIGHT) [m].
    double FrontPosition(const Profile& profile) const;
    /// Root of the decreasing function f(lambda) = target on (0, 5], by bisection.
    template <typename Function> static double Root(Function f, double target);
    /// Analytical temperature (T&S eqs. 4-137 / 4-145, or the error-function solution without latent heat) [C].
    double AnalyticalTemperature(double distance, double time, bool with_latent_heat) const;
    /// Analytical solidus crossing with or without latent heat [m from the LEFT end].
    double AnalyticalFront(double time, bool with_latent_heat) const;
    /// Appends a row to front_rows_ and, at four times, a profile to profile_rows_ (test-108: stores the profile).
    void RecordLatentHeatState();
    /// Explicit finite-difference enthalpy method on a grid refinement_ times finer than the model, same ends and
    /// properties, T(H) by bisection; temperatures at the given distances at each of the given times.
    std::vector<std::vector<double>> EnthalpyMethodReference(const std::vector<double>& distances,
                                                             const std::vector<double>& times) const;
    /// Writes <test>_front.csv (108: <test>_energy.csv) and <test>_profiles.csv next to the configuration file.
    void WriteLatentHeatOutput();
    /// test-106: front and solid-layer temperatures against T&S sec. 4.18.
    void CheckStefanSolidification();
    /// test-107: front travel and contact temperature against T&S sec. 4.19.
    void CheckDikeSolidification();
    /// test-108: energy balance and temperatures against the enthalpy-method reference.
    void CheckWideIntervalCooling();

    static constexpr double dike_half_width_{1005.};  ///< [m], on a control-volume face of the 10 m mesh
    static constexpr double record_interval_{10. * 31536000.};   ///< front_rows_ interval [s]
    static constexpr unsigned refinement_{4U};        ///< reference grid spacing = model spacing / refinement_
    /// 106/107: |front travel / analytical - 1|; the 1 K interval sits inside one 10 m cell, so the solidus crossing is
    /// resolved to about half a cell (5 m of 155 m in 106)
    static constexpr double front_tolerance_{0.05};
    /// |T - T_analytical or reference| / (T_m - T_0); first-order time stepping and the 10 m mesh give < 0.2 %
    static constexpr double temperature_tolerance_{0.01};
    /// 108: |energy residual| / heat conducted out; the end heat flow is estimated from the end-of-step temperatures
    static constexpr double energy_tolerance_{0.02};

    int CalculateMassConductivity();
    int CalculateMassGravityTerm();
    int InitializeFluidPropertiesLinearPressure();
    void InitializeAllNaNsToZero();
    void ScalePermeabilityWithHalite();
    void PrintProgressToScreen();

    // functions to compare data from two vset files
    // used for Benchmarks
    std::map<Point<1U>, std::vector<double>> BuildPropertyMap( const Model<1U>& sg, const std::string& property_name );
    bool ComparePropertyMaps( const std::map<Point<1U>, std::vector<double>>& map_test,
                              const std::map<Point<1U>, std::vector<double>>& map_benchmark,
                              const std::string& property_name,
                              std::string& error_message );



    static constexpr double  year   = 31536000.0;     //  in sec
    static constexpr double  day    = 86400.0;        //  in sec
    static constexpr double  hour   = 3600.0;         //  in sec
    static constexpr double  minute = 60.0;           //  in sec
};

} // csmp

#endif // CVFEM_1D_VVCase
