// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DeformationInducedPermeability.h
//  CSMP_ReservoirSimulator
//
//

#ifndef DEFORMATION_INDUCED_PERMEABILITY_H
#define DEFORMATION_INDUCED_PERMEABILITY_H

#include <iostream>
#include <cmath>

namespace csmp {

/// collection of methods for element-wise calculations of permeabilities of fractures and fault breccia
class DeformationInducedPermeability {
  public:
    explicit DeformationInducedPermeability( double max_aperture );
    ~DeformationInducedPermeability();
  
    /// Snow (1969) formula for parallel flow-aligned fractures
    double keff_UniformAperture( double phi, double matrix_k, double frac_intensity ) const;
  
    /// Kozeny-Carman relationship, version for mean fragment size
    double kKC_BrecciaFromFragSize( double phi, double matrix_k, double meanFragSize, double tau ) const;
  
    /// negative exponential distribution with a flat center part the width of which can be controlled with the pf 
    double k_FaultCoreBreccia( double halfThickness, double x, double plateau_factor, double kmax ) const;
  
    /// skewed Gaussian distribution from combined sine and exponential functions
    double k_DamageZoneFracture( double halfThickness, double x, double skewness_factor, double kmax ) const;
  
  private:
    double AverageAperture (double phi, double fracIntensity) const;
    double max_aperture_;
};

/**
    Idea: call FaultFlowPropertyCalculator first to establish failure criteria,
    stores these in terms of a failure variable on the model.
    
    Use the combined criteria to compute, using this collection of methods
    the dilatation permeability etc. of faults, fractures or deformation bands.
*/

// To do:
//    double keff_LogNormalAperture (double phi, double max_aperture, double frac_intensity);
//    double keff_ExponentialAperture (double phi, double max_aperture, double frac_intensity );
//    double keff_PowerLawAperture (double phi, , double max_aperture, double frac_intensity );


} // end csmp

#endif /* defined(__DEFORMATION_INDUCED_PERMEABILITY_H) */
