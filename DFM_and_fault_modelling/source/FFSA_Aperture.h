#ifndef CSMP_FAR_FIELD_STRESS_APERTURE_APPROXIMATION_H
#define CSMP_FAR_FIELD_STRESS_APERTURE_APPROXIMATION_H

#include "CSMP_definitions.h"

#include <Eigen/Dense> // Eigen library for vectorized math

namespace csmp {

/**
 * @brief Input parameters for the FFSA Aperture calculation.
 *
 * All physical parameters are expected to be Eigen::VectorXd, allowing
 * for simultaneous calculation across multiple fractures (N_frac).
 * Method parameters are strings applied to all fractures.
 *
    For example application see    'Implementation Barton&Bandis-Liem.pdf'
 */
struct FFSA_InputParameters {
    // Physical Parameters (must be of size N_frac)
    Eigen::VectorXd L;           ///< Fracture length [m]
    Eigen::VectorXd alpha;       ///< Fracture angle (between fracture and x-axis) [°]
    Eigen::VectorXd sigma_H;     ///< Max. principal far field stress, SHmax [MPa]
    Eigen::VectorXd sigma_h;     ///< Min. principal far field stress, Shmin [MPa]
    Eigen::VectorXd beta;        ///< Orientation of sigma_H [°]
    Eigen::VectorXd p_f;         ///< Ambient fluid pressure [MPa]
    Eigen::VectorXd JRC;         ///< Joint roughness coefficient, JRC 0-20 [-]
    Eigen::VectorXd sigma_c;     ///< Unconfined compression strength, UCS [MPa]
    Eigen::VectorXd JCS;         ///< Joint wall compressive strength, JCS [MPa]
    Eigen::VectorXd K_ni;        ///< Initial normal stiffness [MPa/mm]
    Eigen::VectorXd vm_factor;   ///< Factor for maximum possible joint closure [-]
    Eigen::VectorXd phi_r;       ///< Residual friction angle [°]
    Eigen::VectorXd E_mod;       ///< E-modulus of fractured rock [MPa]
    Eigen::VectorXd nu;          ///< Poisson's ratio of fractured rock [-]
    Eigen::VectorXd C_g;         ///< Proportionality between displacement and shear stress [-]
    
    // Optional parameter, may be empty or size N_frac
    Eigen::VectorXd M;           ///< Damage coefficient (optional) [-]
    Eigen::VectorXd sigma_EFF;   ///< sigma_eff for mobilization calculations (optional) [MPa]

    // Method Parameters (strings)
    std::string method_displacement;    ///< Method for shear displacement ('mob', 'old')
    std::string method_dilation;        ///< Method for shear dilation ('integrate', 'integrate_pos', 'last')
    std::string method_peak_displacement; ///< Method for peak displacement ('Barton', 'Asadollahi')
    std::string method_sigma_eff;       ///< Method for sigma_EFF calculation ('end', 'onset', 'mean', 'value')

    long N_frac = 0; ///< Number of fractures (derived from vector size)
};



/**
 * @brief Structure holding the detailed debug and intermediate calculation results.
 *
 * All members are Eigen::VectorXd of size N_frac.
 */
struct FFSA_DebugResults {
    // Calculated Stress State
    Eigen::VectorXd sigma_H;
    Eigen::VectorXd sigma_h;
    Eigen::VectorXd beta;
    Eigen::VectorXd p_f;
    Eigen::VectorXd sigma_n;
    Eigen::VectorXd sigma_s;
    Eigen::VectorXd sigma_eff;
    Eigen::VectorXd sigma_EFF;

    // Normal Closure
    Eigen::VectorXd JRC;
    Eigen::VectorXd JCS;
    Eigen::VectorXd sigma_c;
    Eigen::VectorXd a_0;         ///< Initial aperture [mm]
    Eigen::VectorXd vm_factor;
    Eigen::VectorXd v_m;
    Eigen::VectorXd K_ni;
    Eigen::VectorXd delta_n;     ///< Normal closure [mm]

    // Shear Displacement
    Eigen::VectorXd phi_r;
    Eigen::VectorXd E_mod;
    Eigen::VectorXd nu;
    Eigen::VectorXd C_g;
    Eigen::VectorXd G;
    Eigen::VectorXd K_s;
    Eigen::VectorXd delta_s;     ///< Shear displacement [mm]
    Eigen::VectorXd phi_s_mob;   ///< Mobilized friction angle [°]

    // Shear Dilation & Aperture
    Eigen::VectorXd delta_peak;
    Eigen::VectorXd JRC_mob;     ///< Mobilized JRC [-]
    Eigen::VectorXd M;
    Eigen::VectorXd phi_d_mob;   ///< Mobilized dilation angle [°]
    Eigen::VectorXd delta_d;     ///< Shear dilation [mm]
    Eigen::VectorXd a;           ///< Final aperture [mm]
    Eigen::VectorXd k;           ///< Permeability factor [m^2]

    /**
     * @brief Constructor to initialize all Eigen vectors to a specific size.
     * @param size The number of fractures (N_frac).
     */
    FFSA_DebugResults(long size) {
        sigma_H.resize(size); sigma_h.resize(size); beta.resize(size); p_f.resize(size);
        sigma_n.resize(size); sigma_s.resize(size); sigma_eff.resize(size); sigma_EFF.resize(size);
        JRC.resize(size); JCS.resize(size); sigma_c.resize(size); a_0.resize(size);
        vm_factor.resize(size); v_m.resize(size); K_ni.resize(size); delta_n.resize(size);
        phi_r.resize(size); E_mod.resize(size); nu.resize(size); C_g.resize(size);
        G.resize(size); K_s.resize(size); delta_s.resize(size); phi_s_mob.resize(size);
        delta_peak.resize(size); JRC_mob.resize(size); M.resize(size); phi_d_mob.resize(size);
        delta_d.resize(size); a.resize(size); k.resize(size);
    }
};




/**
 * @brief Main class for calculating fracture aperture based on
 * Far-Field Stress Approximation (FFSA) methodology.
 *
 * This class encapsulates the logic from FFSA_Aperture.m and its helper functions.
 */
class FFSA_FractureAperture {
public:
    /**
     * @brief Calculates the shear dilation, displacement, and final aperture for a set of fractures.
     *
     * @param params Input structure containing all required fracture and method parameters.
     * @return std::pair<Eigen::MatrixXd, FFSA_DebugResults> A pair containing:
     * - Matrix of results [a_0, delta_n, delta_s, delta_d, a] (N_frac x 5) [mm].
     * - Detailed debug structure with all intermediate variables.
     */
    std::pair<Eigen::MatrixXd, FFSA_DebugResults> ApertureFromFarFieldStress( const FFSA_InputParameters& params );

private:
    // Constants
    static constexpr double PI = 3.14159265358979323846; ///< Pi constant
    static constexpr double DEG_TO_RAD = PI / 180.0;     ///< Degrees to Radians conversion

    /**
     * @brief C++ implementation of MATLAB's `tand(x)` - tangent of x in degrees.
     * @param x Angle in degrees.
     * @return Tangent of x.
     */
    double tand(double x) const {
        return std::tan(x * DEG_TO_RAD);
    }

    /**
     * @brief C++ implementation of MATLAB's `cosd(x)` - cosine of x in degrees.
     * @param x Angle in degrees.
     * @return Cosine of x.
     */
    double cosd(double x) const {
        return std::cos(x * DEG_TO_RAD);
    }

    /**
     * @brief C++ implementation of MATLAB's `sind(x)` - sine of x in degrees.
     * @param x Angle in degrees.
     * @return Sine of x.
     */
    double sind(double x) const {
        return std::sin(x * DEG_TO_RAD);
    }

    /**
     * @brief Performs 1D linear interpolation (similar to MATLAB's `interp1`).
     *
     * @param x Vector of query points.
     * @param x_data Fixed known x-coordinates of data points.
     * @param y_data Fixed known y-coordinates of data points.
     * @return Eigen::VectorXd Vector of interpolated values at query points x.
     */
    Eigen::VectorXd linear_interp1(const Eigen::VectorXd& x,
                                   const std::vector<double>& x_data,
                                   const std::vector<double>& y_data) const;


    /**
     * @brief Calculates the Mobilized Joint Roughness Coefficient (JRC_mob).
     *
     * Implements logic from FFSA_JRCmob.m.
     *
     * @param d_by_dpeak Ratio shear displacement to peak shear displacement [-].
     * @param sigma_n Normal stress [MPa].
     * @param JRC (Peak) Joint roughness coefficient [-].
     * @param JCS Joint wall compressive strength [MPa].
     * @param phi_r Residual friction angle [°].
     * @return Eigen::VectorXd Vector of Mobilized JRC values [-].
     */
    Eigen::VectorXd FFSA_JRCmob_calc(const Eigen::VectorXd& d_by_dpeak, double sigma_n, double JRC, double JCS, double phi_r) const;


    /**
     * @brief Calculates the shear displacement (delta_s) and mobilized friction angle.
     *
     * Implements logic from FFSA_ShearDisplacement.m.
     *
     * @param sigma_s Shear stress [MPa].
     * @param K_s Shear stiffness of fracture [MPa/mm].
     * @param delta_peak Peak shear displacement [mm].
     * @param JRC Joint roughness coefficient [-].
     * @param JCS Joint wall compressive strength [MPa].
     * @param sigma_eff Effective normal stress (for tau) [MPa].
     * @param sigma_EFF Effective normal stress (for JRC_mob & phi_s_mob) [MPa].
     * @param phi_r Residual friction angle [°].
     * @param method Method to calculate shear displacement ('mob' or 'old').
     * @return std::tuple<double, double, double> A tuple containing:
     * - delta_s: Shear displacement [mm].
     * - JRC_mob: Mobilized joint roughness coefficient [-].
     * - phi_s_mob: Mobilized friction angle [°].
     */
    std::tuple<double, double, double> FFSA_ShearDisplacement_calc(
              double sigma_s, double K_s, double delta_peak, double JRC, double JCS,
              double sigma_eff, double sigma_EFF, double phi_r, const std::string& method) const;


    /**
     * @brief Calculates the shear dilation (delta_d) and mobilized dilation angle.
     *
     * Implements logic from FFSA_ShearDilation.m.
     *
     * @param delta_s Shear displacement [mm].
     * @param delta_peak Peak shear displacement [mm].
     * @param JRC Joint roughness coefficient [-].
     * @param JCS Joint wall compressive strength [MPa].
     * @param sigma_eff Effective normal stress (for dilation) [MPa].
     * @param phi_r Residual friction angle [°].
     * @param M Damage coefficient [-].
     * @param method Method how to calculate shear dilation ('integrate', 'integrate_pos', 'last').
     * @return std::tuple<double, double, double> A tuple containing:
     * - delta_d: Shear dilation [mm].
     * - JRC_mob: Mobilized JRC [-].
     * - phi_d_mob: Mobilized dilation angle [°].
     */
    std::tuple<double, double, double> FFSA_ShearDilation_calc(
              double delta_s, double delta_peak, double JRC, double JCS, double sigma_eff,
              double phi_r, double M, const std::string& method) const;
};

}

#endif // CSMP_FAR_FIELD_STRESS_APERTURE_APPROXIMATION_H
