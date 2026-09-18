// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FFSA_Aperture.h"

using namespace std;

namespace csmp {

// Define a placeholder for the global variable set in MATLAB's FFSA_Aperture.m for the current fracture iteration
// Since we are iterating fracture by fracture, these are single values in the C++ implementation.
// While, the FFSA_JRCmob calculation is designed to handle vectors, in the main loop,
// it will be called with a vector of size 1, or a vector of intermediate steps.

// --- Private Helper Functions ---

/**
 * @brief Helper performs 1D linear interpolation (similar to MATLAB's `interp1`).
 *
 * @param x Vector of query points.
 * @param x_data Fixed known x-coordinates of data points.
 * @param y_data Fixed known y-coordinates of data points.
 * @return Eigen::VectorXd Vector of interpolated values at query points x.
 */
Eigen::VectorXd FFSA_FractureAperture::linear_interp1(const Eigen::VectorXd& x,
                                              const vector<double>& x_data,
                                              const vector<double>& y_data) const {
    long N = x.size();
    Eigen::VectorXd y_interp(N);

    // Find the min/max of the data for clamping
    double x_min = x_data.front();
    double x_max = x_data.back();

    for (long i = 0; i < N; ++i) {
        double xi = x(i);

        // Clamping (extrapolation in MATLAB's interp1 is disabled by default,
        // but here we clamp to the nearest data point for robustness)
        if (xi <= x_min) {
            y_interp(i) = y_data.front();
            continue;
        }
        if (xi >= x_max) {
            y_interp(i) = y_data.back();
            continue;
        }

        // Linear search for interval (x_data is sorted)
        auto it_upper = upper_bound(x_data.begin(), x_data.end(), xi);
        
        // it_upper points to x_data[k] where x_data[k-1] < xi <= x_data[k]
        long k = distance(x_data.begin(), it_upper);
        
        // k_lower must be at least 1 (since k starts at 1)
        size_t k_upper = static_cast<size_t>(k);
        size_t k_lower = static_cast<size_t>(k - 1);

        // Linear interpolation formula: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
        double x1 = x_data[k_lower];
        double y1 = y_data[k_lower];
        double x2 = x_data[k_upper];
        double y2 = y_data[k_upper];
        
        double slope = (y2 - y1) / (x2 - x1);
        y_interp(i) = y1 + (xi - x1) * slope;
    }
    return y_interp;
}




/**
 * @brief Calculates the Mobilized Joint Roughness Coefficient (JRC_mob).
 *
 * Implements logic from FFSA_JRCmob.m using linear interpolation.
 *
 * @param d_by_dpeak Ratio shear displacement to peak shear displacement [-].
 * @param sigma_n       Joint-surface normal stress [MPa].
 * @param JRC (Peak)     Joint roughness coefficient [-].
 * @param JCS                Joint wall compressive strength [MPa].
 * @param phi_r            Residual friction angle [°].
 * @return Eigen::VectorXd Vector of Mobilized JRC values [-].
 */
Eigen::VectorXd FFSA_FractureAperture::FFSA_JRCmob(const Eigen::VectorXd& d_by_dpeak, double sigma_n, double JRC, double JCS, double phi_r) const {
    // i = JRC*log10(JCS/sigma_n); % [°]
    double i_val = JRC * log10(JCS / sigma_n);
    
    // Check for division by zero or log(0) which implies tension or zero JRC
    if (i_val <= 0.0) {
        // If i is zero or negative (JCS <= sigma_n), the denominator is zero/negative.
        // The mobilization factor is based on the residual friction angle.
        // If i_val is non-positive, the JRCmob_by_JRCpeak list is likely compromised.
        // We'll proceed with JRC_mob = 0 for high displacement, or 1 for low displacement.
        // For robustness, we will assume the ratio is 1 (full JRC) if i_val is non-positive,
        // unless d_by_dpeak is large, in which case JRC_mob should approach 0.
        
        // This is a simplification based on the structure of the input data points:
        // list_JRCmob_by_JRCpeak= [-phi_r/i  0   0.75  1  0.85 0.7 0.5 0]
        // If i_val <= 0, the first element is undefined or large.
        // We will default the JRC_mob ratio list to {0, 0, 0.75, 1, 0.85, 0.7, 0.5, 0}
        
        // However, since the first data point corresponds to d_by_dpeak = 0, JRC_mob must be > 0.
        // Since this part of the original code is likely for the case when phi_r is mobilized immediately,
        // we'll enforce a minimal positive value for i_val for the first point calculation, 
        // but otherwise use a standard 0-to-1 ratio lookup table.
        
        // Use the values {0, 0.75, 1, 0.85, 0.7, 0.5, 0} for d_by_dpeak >= 0.3
        // For d_by_dpeak = 0, JRC_mob is generally expected to be low or 0 for the mobilization model.
        // Let's stick to the core interpolation points and handle the edge case.
      }

    // Fixed interpolation data points
    vector<double> list_d_by_dpeak = {0.0, 0.3, 0.6, 1.0, 2.0, 4.0, 10.0, 100.0};
    vector<double> list_JRCmob_by_JRCpeak = {0.0, 0.0, 0.75, 1.0, 0.85, 0.7, 0.5, 0.0};
    
    // Origina MATLAB code includes: list_JRCmob_by_JRCpeak= [-phi_r/i 0 0.75 1 0.85 0.7 0.5 0]
    // The entry at d_by_dpeak = 0 is non-physical if i_val is not positive.
    // Assuming phi_r is positive, we substitute -phi_r/i_val for the first element:
    if (i_val > 1e-6) { // Avoid division by zero
        list_JRCmob_by_JRCpeak[0] = -phi_r / i_val;
    } else {
        // If i_val is close to zero, the mobilization factor is extremely large/small.
        // We cap the initial ratio to 0 to prevent non-physical initial JRC_mob values, 
        // as the function is mainly used for d_by_dpeak > 0.
        list_JRCmob_by_JRCpeak[0] = 0.0; 
    }

    // Handling warnings and clamping (MATLAB code: if any(d_by_dpeak<0) -> d_by_dpeak= abs(d_by_dpeak);)
    // We assume the input vector d_by_dpeak only contains non-negative values from the caller (delta_s >= 0)
    // d_by_dpeak(d_by_dpeak>100)= 100;
    Eigen::VectorXd d_clamped = d_by_dpeak.array().min(100.0).matrix();

    // Perform interpolation
    Eigen::VectorXd JRCmob_by_JRCpeak = linear_interp1(d_clamped, list_d_by_dpeak, list_JRCmob_by_JRCpeak);

    // JRC_mob = JRCmob_by_JRCpeak*JRC;
    return JRCmob_by_JRCpeak * JRC;
}




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
 * @return tuple<double, double, double> A tuple containing:
 * - delta_d: Shear dilation [mm].
 * - JRC_mob: Mobilized JRC [-].
 * - phi_d_mob: Mobilized dilation angle [°].
 */
tuple<double, double, double> FFSA_FractureAperture::FFSA_ShearDilation(
                              double delta_s, double delta_peak, double JRC, double JCS, double sigma_eff,
                              double phi_r, double M, const string& method) const
{
    // Check if effective normal stress is positive
    if (sigma_eff < 0.0) {
        throw runtime_error("Negative sigma_eff in ShearDilation: Tensile opening is not considered.");
    }

    // Check if shear displacement is zero (or even negative?)
    if (delta_s <= 0.0) {
        // Mobilized JRC (call with d_by_dpeak = 0)
        Eigen::VectorXd d_zero(1); d_zero << 0.0;
        double JRC_mob = FFSA_JRCmob(d_zero, sigma_eff, JRC, JCS, phi_r)(0);
        // Mobilized dilation angle
        double phi_d_mob = numeric_limits<double>::quiet_NaN();
        // Shear dilation
        double delta_d = 0.0;
        return {delta_d, JRC_mob, phi_d_mob};
    }
    
    // Calculation factors
    const double log_factor = log10(JCS / sigma_eff);

    if (method == "last") {
        // Mobilized JRC
        Eigen::VectorXd d_by_dpeak_vec(1); d_by_dpeak_vec << delta_s / delta_peak;
        double JRC_mob = FFSA_JRCmob(d_by_dpeak_vec, sigma_eff, JRC, JCS, phi_r)(0);

        // Mobilized dilation angle
        double phi_d_mob = (1.0 / M) * JRC_mob * log_factor;

        // Shear dilation
        double delta_d = delta_s * tand(phi_d_mob);

        return {delta_d, JRC_mob, phi_d_mob};
    } 
    
    // --- INTEGRATION METHODS ---
    else if (method == "integrate" || method == "integrate_pos") {
        
        // Numerical integration width
        double d_delta_s = (method == "integrate") ? 0.1 : 0.01;
        
        Eigen::VectorXd all_delta_s;
        double integration_step;

        if (delta_s >= d_delta_s) {
            // MATLAB: all_delta_s = d_delta_s/2:d_delta_s:delta_s;
            // The MATLAB step size is (delta_s / d_delta_s), but here we use a fixed step.
            long num_steps = static_cast<long>(floor(delta_s / d_delta_s));
            integration_step = d_delta_s;

            // Generate sequence: [0.05, 0.15, 0.25, ..., num_steps * d_delta_s - 0.05]
            all_delta_s.resize(num_steps);
            for (long j = 0; j < num_steps; ++j) {
                all_delta_s(j) = d_delta_s / 2.0 + j * d_delta_s;
            }
            
            // Check if we need to include the final delta_s as the last point if it doesn't align
            if (abs(delta_s - all_delta_s.tail(1)(0) - d_delta_s/2.0) > 1e-9) {
                // Not exactly matching the MATLAB range behavior, but covering the range
                // We ensure the last point is less than or equal to delta_s
            }

        } else {
            // delta_s is smaller than the default step size, integrate over just one point at delta_s
            all_delta_s.resize(1);
            all_delta_s << delta_s;
            integration_step = delta_s; // Use delta_s as the integration width
        }

        // Mobilized JRC
        Eigen::VectorXd d_by_dpeak = all_delta_s / delta_peak;
        Eigen::VectorXd all_JRC_mob = FFSA_JRCmob(d_by_dpeak, sigma_eff, JRC, JCS, phi_r);

        // Mobilized dilation angle [°]
        // all_phi_d_mob = 1/M*all_JRC_mob*log10(JCS/sigma_eff);
        Eigen::VectorXd all_phi_d_mob = (1.0 / M) * all_JRC_mob.array() * log_factor;

        // Incremental shear dilation [mm]
        // d_delta_d = d_delta_s * tand(all_phi_d_mob);
        // We use a custom loop for tand to work on the Eigen VectorXd array
        Eigen::VectorXd d_delta_d = all_phi_d_mob.unaryExpr([this](double phi) {
            return this->tand(phi);
        }).array() * integration_step;
        
        // Total shear dilation
        double delta_d = 0.0;
        
        if (method == "integrate") {
            // delta_d = sum(d_delta_d);
            delta_d = d_delta_d.sum();
        } else { // 'integrate_pos'
            // delta_d = sum(d_delta_d(d_delta_d>=0));
            // Custom reduction for positive increments
            for (long j = 0; j < d_delta_d.size(); ++j) {
                if (d_delta_d(j) >= 0.0) {
                    delta_d += d_delta_d(j);
                }
            }
        }

        // Final values (last element in the vector)
        double JRC_mob = all_JRC_mob.tail(1)(0);
        double phi_d_mob = all_phi_d_mob.tail(1)(0);

        return {delta_d, JRC_mob, phi_d_mob};

    } else {
        throw runtime_error("FFSA_ShearDilation: Unknown method " + method);
    }
}






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
 * @return tuple<double, double, double> A tuple containing:
 * - delta_s: Shear displacement [mm].
 * - JRC_mob: Mobilized joint roughness coefficient [-].
 * - phi_s_mob: Mobilized friction angle [°].
 */
tuple<double, double, double> FFSA_FractureAperture::FFSA_ShearDisplacement(
                              double sigma_s, double K_s, double delta_peak, double JRC, double JCS,
                              double sigma_eff, double sigma_EFF, double phi_r, const string& method) const
{
    double delta_s = 0.0;
    double JRC_mob = 0.0;
    double phi_s_mob = 0.0;
    
    // Check for negative effective normal stress
    if (sigma_eff < 0.0 || sigma_EFF < 0.0) {
        throw runtime_error("Negative effective normal stress detected in ShearDisplacement. Tensile opening is not considered.");
    }

    // --- Legacy Code: option to calculate delta_s as in ECMOR22 paper ---
    if (method == "old") {
        // Excess shear stress [MPa]
        double delta_tau = abs(sigma_s) - sigma_eff * tand(phi_r);
        
        // Shear displacement [mm] (must be non-negative)
        delta_s = max(delta_tau / K_s, 0.0);
        
        // Ratio of shear displacement over peak shear displacement [-]
        double d_by_dpeak = delta_s / delta_peak;
        
        // Mobilized JRC [-]
        Eigen::VectorXd d_by_dpeak_vec(1);
        d_by_dpeak_vec(0) = d_by_dpeak;
        JRC_mob = FFSA_JRCmob(d_by_dpeak_vec, sigma_EFF, JRC, JCS, phi_r)(0);
        
        // Mobilized friction angle [°]
        phi_s_mob = JRC_mob * log10(JCS / sigma_EFF) + phi_r;
        
        return {delta_s, JRC_mob, phi_s_mob};
    } 
    
    // --- ITERATIVE SEARCH CODE: method='mob' ---
    else if (method != "mob") {
        throw runtime_error("Unknown method: " + method);
    }
    
    // --- Initial Candidate List Generation ---
    vector<double> temp_delta_s_cand;
    
    // [0:0.05:1.9]
    for (double d = 0.0; d <= 1.9 + 1e-9; d += 0.05) {
        temp_delta_s_cand.push_back(d * delta_peak);
    }
    // [2:10]
    for (double d = 2.0; d <= 10.0 + 1e-9; d += 1.0) {
        temp_delta_s_cand.push_back(d * delta_peak);
    }
    // [20:10:100]
    for (double d = 20.0; d <= 100.0 + 1e-9; d += 10.0) {
        temp_delta_s_cand.push_back(d * delta_peak);
    }
    
    // Convert to Eigen VectorXd (and remove duplicates)
    sort(temp_delta_s_cand.begin(), temp_delta_s_cand.end());
    temp_delta_s_cand.erase(unique(temp_delta_s_cand.begin(), temp_delta_s_cand.end(), [](double a, double b) {
        return abs(a - b) < 1e-9;
    }), temp_delta_s_cand.end());

    Eigen::VectorXd delta_s_cand = Eigen::Map<Eigen::VectorXd>(temp_delta_s_cand.data(), static_cast<long>(temp_delta_s_cand.size()) );
    
    // Shear displacement for which relaxed shear stress is zero [mm]
    double delta_s_max = abs(sigma_s) / K_s;

    // Limit candidates of delta_s so that only 1 entry is >= delta_s_max
    long count_below_max = (delta_s_cand.array() < delta_s_max).count();
    long desired_length = min(max(count_below_max + 1, 2L), (long)delta_s_cand.size());
    
    // Resize the Eigen vector conservatively
    delta_s_cand.conservativeResize(desired_length);

    // Iterative Search
    int i = 0;
    const int max_iter = 10;
    const double tolerance = 1e-3; // [MPa]
    long idx = 0;
    Eigen::VectorXd JRC_mob_vec_eigen, phi_s_mob_vec_eigen;

    while (true) {
        i++;
        
        // 1. Calculate relaxed shear stress (sigma_s_rel) [MPa]
        // sigma_s_rel = max( abs(sigma_s) - K_s*delta_s_cand , 0)
        Eigen::VectorXd sigma_s_rel = (abs(sigma_s) - K_s * delta_s_cand.array()).cwiseMax(0.0);
        
        // 2. Calculate mobilized shear strength (tau) [MPa]
        
        // Ratio of shear displacement over peak shear displacement [-]
        Eigen::VectorXd d_by_dpeak = delta_s_cand / delta_peak;
        
        // Mobilized JRC [-]
        JRC_mob_vec_eigen = FFSA_JRCmob(d_by_dpeak, sigma_EFF, JRC, JCS, phi_r);
        
        // Mobilized friction angle [°]
        // phi_s_mob = JRC_mob * log10(JCS/sigma_EFF) + phi_r
        double log_factor = log10(JCS / sigma_EFF);
        phi_s_mob_vec_eigen = JRC_mob_vec_eigen.array() * log_factor + phi_r;

        // Shear strength (tau) [MPa]
        // tau = sigma_eff * tand(phi_s_mob)
        Eigen::VectorXd tau = phi_s_mob_vec_eigen.unaryExpr([this, sigma_eff](double phi) {
            return sigma_eff * this->tand(phi);
        });
        
        // 3. Calculate residual [MPa]
        // residual = abs(sigma_s_rel - tau)
        Eigen::VectorXd residual = (sigma_s_rel - tau).cwiseAbs();
        
        // 4. Find index where residual is minimal
        double min_residual_value = residual.minCoeff(&idx);

        // 5. Break loop if residual is small enough or after max_iter iterations
        if (min_residual_value <= tolerance || i >= max_iter) {
            break;
        }
        
        // 6. Choose one index lower and higher as bounds for next iteration
        long idx1 = max(idx - 1, 0L);
        long idx2 = min(idx + 1, (long)delta_s_cand.size() - 1);
        
        // Create new candidate list (Refined grid search)
        double delta_s_lower = delta_s_cand(idx1);
        double delta_s_upper = delta_s_cand(idx2);
        double delta_s_range = delta_s_upper - delta_s_lower;

        // MATLAB: [0 .25 .4 .45 .48 .5 .52 .55 .6 .75 1]*range + delta_s_lower
        Eigen::VectorXd factors(11);
        factors << 0.0, 0.25, 0.4, 0.45, 0.48, 0.5, 0.52, 0.55, 0.6, 0.75, 1.0;
        
        // Eigen vectorized multiplication and addition
        Eigen::VectorXd new_delta_s_cand = factors.array() * delta_s_range + delta_s_lower;
        
        // Limit candidates of delta_s so that only 1 entry is >= delta_s_max
        count_below_max = (new_delta_s_cand.array() < delta_s_max).count();
        desired_length = min(max(count_below_max + 1, 2L), (long)new_delta_s_cand.size());
        
        delta_s_cand = new_delta_s_cand.head(desired_length);
    }

    // Final results are the values corresponding to the minimal residual index (idx)
    delta_s = delta_s_cand(idx);
    JRC_mob = JRC_mob_vec_eigen(idx);
    phi_s_mob = phi_s_mob_vec_eigen(idx);
    
    return {delta_s, JRC_mob, phi_s_mob};
    
} // end FFSA_ShearDisplacement






// =====================================================================
//
// --- Public Main Calculation Function ---
//
// =====================================================================

/**
 * @brief Calculates the shear dilation, displacement, and final aperture for a set of fractures.
 *
 * Implements logic from M. Liem's FFSA_Aperture.m.
 *
 * @param params Input structure containing all required fracture and method parameters.
 *
 * @return pair<Eigen::MatrixXd, FFSA_DebugResults> A pair containing:
 * - Matrix of results [a_0, delta_n, delta_s, delta_d, a] (N_frac x 5) [mm].
 * - Detailed debug structure with all intermediate variables.
 *
 * @attention While (by constract with the Cruikshank model) shear dilatation can occur at a positive (compressive) effective stress,
 * at appreciable burial / hydrostatic reservoir conditions, the effective stress is likely to be so high that Barton's model will predict asperity crushing.
 * This will manifest as a negative mobilised JRC and zero dilation and normal compression = aperture reduction during shear.
 *  Is this physically correct? - search the literature!
 */
pair<Eigen::MatrixXd, FFSA_DebugResults> FFSA_FractureAperture::ApertureFromFarFieldStress(const FFSA_InputParameters& params) {
    long N_frac = params.N_fractures;
    if (N_frac == 0) {
        throw runtime_error("Input vector size (N_frac) is zero.");
    }

    // Initialize result matrix and debug structure
    Eigen::MatrixXd result(N_frac, 5);
    FFSA_DebugResults debug(N_frac);

    // Get method strings
    const string method_disp = params.method_displacement;
    const string method_dil = params.method_dilation;
    const string method_peak = params.method_peak_displacement;
    const string method_sigma_eff = params.method_sigma_eff;

    for (long frac = 0; frac < N_frac; ++frac) {
        // --- Extract parameters for current fracture ---
        double L = params.L(frac);
        double alpha = params.alpha(frac);
        double sigma_H = params.sigma_H;
        double sigma_h = params.sigma_h;
        double beta = params.beta;
        double p_f = params.p_f(frac);
        double JRC = params.JRC(frac);
        double JCS = params.JCS(frac);
        double sigma_c = params.sigma_c(frac);
        double K_ni = params.K_ni(frac);
        double vm_factor = params.vm_factor(frac);
        double phi_r = params.phi_r(frac);
        double E_mod = params.E_mod(frac);
        double nu = params.nu(frac);
        double C_g = params.C_g(frac);
        
        // Angle [°] between sigma_H (beta) and fracture normal (alpha)
        /*
            In standard geomechanics, if sigma_H (Maximum Horizontal Stress) is aligned with the x-axis and
            sigma_h (Minimum Horizontal Stress) with the y-axis, theta represents the angle between the maximum
            principal stress direction sigma_H and the fracture plane.
        */
        double theta = 90.0 + alpha - beta;
        
        // ---------------------------------------------------------------------
        // Part 1: Calculate effective stress
        double sigma_n = sigma_H * cosd(theta) * cosd(theta) + sigma_h * sind(theta) * sind(theta);
        double sigma_s = (sigma_h - sigma_H) * cosd(theta) * sind(theta);
        double sigma_eff = sigma_n - p_f;
        
        if (sigma_eff < 0.0) {
            string err_msg = "Negative sigma_eff = " + to_string(sigma_eff) + " MPa. Tensile opening!";
            throw runtime_error(err_msg);
        }
        
        // ---------------------------------------------------------------------
        // Part 2: Normal closure
        
        // Initial fracture aperture [mm]
        double a_0 = JRC / 5.0 * (0.2 * sigma_c / JCS - 0.1);
        
        // Maximum closure [mm]
        double v_m = vm_factor * a_0;
        
        // Closure due to normal stress [mm]
        double delta_n = sigma_eff * v_m / (K_ni * v_m + sigma_eff);
        
        // ---------------------------------------------------------------------
        // Part 3: Shear displacement
        
        // Shear modulus [MPa]
        double G = E_mod / (2.0 * (1.0 + nu));
        
        // Shear stiffness [MPa/mm]. L is in [m], K_s in [MPa/mm], so L*1e3 converts m to mm
        double K_s = C_g * G / (L * 1e3);
        
        // Calculate sigma_EFF (Effective normal stress for mobilization)
        double sigma_EFF;
        
        if (method_sigma_eff == "onset") {
            sigma_EFF = sigma_n - sigma_s / tand(phi_r);
        } else if (method_sigma_eff == "mean") {
            sigma_EFF = sigma_n - 0.5 * (sigma_s / tand(phi_r) + p_f);
        } else if (method_sigma_eff == "value") {
            if (params.sigma_EFF.size() != N_frac) {
                throw runtime_error("method.sigma_eff='value' but sigma_EFF vector size is inconsistent.");
            }
            sigma_EFF = params.sigma_EFF(frac);
        } else { // 'end' (default) or unknown
            sigma_EFF = sigma_eff;
        }

        if (sigma_EFF < 0.0) {
            string err_msg = "Negative sigma_EFF = " + to_string(sigma_EFF) + " MPa. Tensile opening in mobilization calculation!";
            throw runtime_error(err_msg);
        }

        // Peak shear displacement [mm]
        double delta_peak;
        if (method_peak == "Barton") {
            delta_peak = L / 500.0 * pow(JRC / L, 0.33) * 1e3;
        } else if (method_peak == "Asadollahi") {
            double angle_term = JRC * log10(JCS / sigma_EFF);
            delta_peak = 0.0077 * pow(L, 0.45) * pow(sigma_EFF / JCS, 0.34) * cosd(angle_term) * 1e3;
        } else {
            throw runtime_error("Unknown method for calculating peak shear displacement: " + method_peak);
        }
        
        // Shear displacement, Mobilized JRC, Mobilized friction angle
        double delta_s, JRC_mob_disp, phi_s_mob;
        tie(delta_s, JRC_mob_disp, phi_s_mob) = FFSA_ShearDisplacement(
            sigma_s, K_s, delta_peak, JRC, JCS, sigma_eff, sigma_EFF, phi_r, method_disp);
        
        // ---------------------------------------------------------------------
        // Part 4: Shear dilation

        // Damage coefficient [-]
        double M;
        bool M_provided = params.M.size() == N_frac;
        if (M_provided && !isnan(params.M(frac))) {
            M = params.M(frac);
        } else {
            M = 0.7 + JRC / (12.0 * log10(JCS / sigma_EFF));
        }
        
        // Shear dilation, Mobilized JRC (for dilation), Mobilized dilation angle
        double delta_d, JRC_mob_dil, phi_d_mob;
        tie(delta_d, JRC_mob_dil, phi_d_mob) = FFSA_ShearDilation(
            delta_s, delta_peak, JRC, JCS, sigma_EFF, phi_r, M, method_dil);
        
        // We use JRC_mob_dil for debug output as it comes from the dilation calculation
        double JRC_mob = JRC_mob_dil; 

        // ---------------------------------------------------------------------
        // Part 5: Aperture
        double a = a_0 - delta_n + delta_d;
        
        // Permeability factor k [m^2]. Aperture 'a' is in [mm], so a/1e3 is in [m].
        double k = pow(a / 1e3, 2) / 12.0;
        
        // ---------------------------------------------------------------------
        // Collecting results
        result(frac, 0) = a_0;
        result(frac, 1) = delta_n;
        result(frac, 2) = delta_s;
        result(frac, 3) = delta_d;
        result(frac, 4) = a;
        
        debug.sigma_H(frac) = sigma_H; debug.sigma_h(frac) = sigma_h; debug.beta(frac) = beta;
        debug.p_f(frac) = p_f; debug.sigma_n(frac) = sigma_n; debug.sigma_s(frac) = sigma_s;
        debug.sigma_eff(frac) = sigma_eff; debug.sigma_EFF(frac) = sigma_EFF;
        
        debug.JRC(frac) = JRC; debug.JCS(frac) = JCS; debug.sigma_c(frac) = sigma_c;
        debug.a_0(frac) = a_0; debug.vm_factor(frac) = vm_factor; debug.v_m(frac) = v_m;
        debug.K_ni(frac) = K_ni; debug.delta_n(frac) = delta_n;
        
        debug.phi_r(frac) = phi_r; debug.E_mod(frac) = E_mod; debug.nu(frac) = nu;
        debug.C_g(frac) = C_g; debug.G(frac) = G; debug.K_s(frac) = K_s;
        debug.delta_s(frac) = delta_s; debug.phi_s_mob(frac) = phi_s_mob;
        
        debug.delta_peak(frac) = delta_peak; debug.JRC_mob(frac) = JRC_mob; debug.M(frac) = M;
        debug.phi_d_mob(frac) = phi_d_mob; debug.delta_d(frac) = delta_d;
        debug.a(frac) = a; debug.k(frac) = k;
    }

    return {result, debug};
}




/**
 * @brief demo_FFSA:  function for demonstration and testing.
 *
 *1. What alpha represents in Barton-style models

    In Barton / Bandis / Choubey–type joint models, alpha (degrees) is used to:
      •	resolve stresses onto the joint plane
      •	relate global stresses to joint-normal and joint-shear components
      •	define directional permeability / stiffness effects

    In 2D, alpha is,  the angle between the joint normal and the global x-axis

    $$\alpha = \angle(\mathbf n, \mathbf e_x)$$

    where:
      •	\mathbf n = unit normal to the fracture plane
      •	\mathbf e_x = (1,0)

    2. Why the normal (not the plane or trace)

    Using the normal instead of the plane itself:
      •	avoids a 90° offset everywhere
      •	gives a unique orientation (up to sign)
      •	simplifies stress projection:
   
   $$ \sigma_n = \mathbf n \cdot \boldsymbol\sigma \cdot \mathbf n$$

    This is why Barton, Bandis, and later implementations nearly always parameterise orientation via the normal.

    3. Range and symmetry

    Because a fracture plane has no direction:
      •	\mathbf n and -\mathbf n represent the same plane
      •	orientations differing by 180° are identical

    Therefore, alpha is usually reduced to:

    \alpha \in [0^\circ, 90^\circ]

    This is exactly the folding step  implemented here.

    4. What Barton does not mean by alpha

    Angle alpha is not:
      •	the dip angle (as it was originally, when sliding tests were performed rotating the sample)
      •	the strike angle
      •	the angle between the fracture trace and x
      •	a signed rotation angle

    Those are different geological conventions.

    5. Typical use of alpha in equations

    You will see \alpha appear in expressions like:
      •	directional stiffness:
    k_n(\alpha),\; k_s(\alpha)
      •	stress resolution:
    \sigma_n = \sigma_x \cos^2\alpha + \sigma_y \sin^2\alpha + 2\tau_{xy}\sin\alpha\cos\alpha
      •	permeability tensors aligned with fracture sets

    All of these assume \alpha is the normal – x-axis angle.

 */
int demo_FFSA() {
    FFSA_FractureAperture calculator;

    try {
        // --- Setup Parameters for 2 Fractures ---
        long N_frac = 2;
        FFSA_InputParameters params;
        params.N_fractures = N_frac;

        int frac_number{1};
        // Fracture 1: (strong shear, low normal effective stress)
        // Fracture 2: (weak shear, high normal stress)
        params.sigma_H = (frac_number==1) ? 15.0 : 20.0; // SH max (horizontal XZ plane is assumed) [MPa]
        params.sigma_h = (frac_number==1) ?  5.0 : 10.0; // Sh min (least principal stress) [MPa]
        params.beta    = (frac_number==1) ?  0.0 : 30.0; // Orientation of SH max (deviation from Y axis) [°]
        params.alpha.resize(N_frac); params.alpha <<        45.0,    60.0; // angle of fracture to x-axis [°]
        params.L.resize(N_frac); params.L                << 10.0,     5.0; // fracture length [m]
        params.p_f.resize(N_frac); params.p_f         <<     9.0,     2.0; // fluid pressure [MPa]
        params.JRC.resize(N_frac); params.JRC         <<    10.0,    15.0; // JRC [-] (rough fracture)
        params.sigma_c.resize(N_frac); params.sigma_c <<   100.0,   150.0; // UCS [MPa]
        params.JCS.resize(N_frac); params.JCS <<            80.0,   120.0; // [MPa]
        params.K_ni.resize(N_frac); params.K_ni <<        1000.0,  2000.0; // [MPa/mm]
        params.vm_factor.resize(N_frac); params.vm_factor << 0.5,     0.6; // maximal joint closure [-]
        params.phi_r.resize(N_frac); params.phi_r <<        30.0,    25.0; // residual friction angle [°]
        params.E_mod.resize(N_frac); params.E_mod <<     30000.0, 50000.0; // Young's modulus [MPa] (limestone ~40GPa)
        params.nu.resize(N_frac); params.nu <<               0.25,    0.2; // Poisson's ratio [-]
        params.C_g.resize(N_frac); params.C_g <<             1.0,     1.2; // shear vs. dilation proportionality [-]
        
        // Optional M parameter (testing the optional field logic)
        params.M.resize(N_frac); params.M.fill(numeric_limits<double>::quiet_NaN());
        params.M(0) = 0.8; // Specify M for frac 1

        // Method strings
        params.method_displacement = "mob"; // Iterative search
        params.method_dilation = "integrate_pos"; // Integration, positive increments only
        params.method_peak_displacement = "Asadollahi";
        params.method_sigma_eff = "end"; // sigma_EFF = sigma_eff

        // --- Run Calculation ---
        auto [results_matrix, debug_data] = calculator.ApertureFromFarFieldStress(params);

        // --- Print Results ---
        cout << "=================================================================\n";
        cout << "FFSA Aperture Calculation Results (N_frac = " << N_frac << ")\n";
        cout << "=================================================================\n";
        cout << "Units: L [m], Stress [MPa], Displacement [mm]\n\n";

        for (long i = 0; i < N_frac; ++i) {
            cout << "--- Test Fracture " << i + 1 << " ---\n";
            cout << "  Calculated Stresses: sigma_n = " << debug_data.sigma_n(i)
                      << " MPa, sigma_s = " << debug_data.sigma_s(i)
                      << " MPa, sigma_eff = " << debug_data.sigma_eff(i) << " MPa\n";
            cout << "  Initial Aperture (a_0): " << debug_data.a_0(i) << " mm\n";
            cout << "  Normal Closure (delta_n): " << debug_data.delta_n(i) << " mm\n";
            cout << "  Shear Displacement (delta_s): " << debug_data.delta_s(i) << " mm\n";
            cout << "  Shear Dilation (delta_d): " << debug_data.delta_d(i) << " mm\n";
            cout << "  Final Aperture (a): " << debug_data.a(i) << " mm\n";
            cout << "  Mobilized JRC: " << debug_data.JRC_mob(i) << " (-)\n";
            cout << "  Mobilized Dilation Angle: " << debug_data.phi_d_mob(i) << " deg\n";
            cout << "  Permeability Factor (k): " << debug_data.k(i) << " m^2\n\n";
        }

        cout << "Results Matrix [a_0, delta_n, delta_s, delta_d, a]:\n" << results_matrix << "\n";

    } catch (const exception& e) {
        cerr << "Calculation Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

} // end csmp
