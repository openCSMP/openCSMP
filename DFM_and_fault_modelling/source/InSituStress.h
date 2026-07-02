//
//  InSituStress.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/18/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//
#ifndef IN_SITU_STRESS_H
#define IN_SITU_STRESS_H

#include "TensorVariable.h"
#include "StressRotate.h"
#include "Point.h"

namespace csmp {

enum STRESS_REGIME { ISOSTATIC, NORMAL_FAULTING, TRANSTENSIONAL, TRANSPRESSIONAL, REVERSE_FAULTING };

std::string  parseStressRegime( STRESS_REGIME );

/// Andersonian stress measurement and suitable utility conversion functions
class InSituStress {
  public:
      /// initialization
      InSituStress( double depth, double rdensity, double overburden, double Sv, double SH, double Sh, double trend );
      InSituStress( const InSituStress& );
//      ~InSituStress() = default;
  
      /// should normally be the location of the stress measurement, but can be shifted for display purposes
      void SampleLocation( double model_X_is_east, double model_Y_is_elevation, double model_Z_is_south );

      /// the location where the stress measurement was taken
      Point<3U> const& SampleLocation() const;
      
      /// the TVD of the stress measurement
      void Depth( double tvd );
  
      /// rotated full stress tensor with principal stresses augmented by offset 
      void CartesianStress( TensorVariable<3U>&, double isotatic_stress_offset=0. ) const;
  
      /// output of compact representation of stress state
      StressRotate& StressState() const;
      /// confining pressure = overburden stress
      double P_conf() const;
      /// distance from the seafloor / earth surface (that can be different from the y-coordinate)
      double Depth() const;
      /// vertical stress normalized by overburden pressure (either of sigma1 to 3 depending on stress regime)
      double Sv() const;
      /// maximum horizontal stress, normalized by overburden stress
      double SH() const;
      /// minimum horizontal stress, normalized by overburden stress
      double Sh() const;
      /// azimuth = XZ-plane counter-clockwise deviation of sigma1 or SH expressed in degrees (N-S=0=180^o), assuming Andersonian stresses 
      double Trend() const;
      void   Trend( double new_XZ_plane_max_stress_azimuth ) { trend_azimuth_ = new_XZ_plane_max_stress_azimuth; }
      /// density of the rock around the site of stress measurement (kg/m3)
      double Density() const;
      /// the stress regime inferred from the magnitudes of the principal stresses
      STRESS_REGIME StressRegime() const;
  
      /// reports object state
      void Out() const;
  
      /// output the current stress state in a vector representation
      void OutStressTensorToVTK( const std::string& filename ) const;
      /// output for a specific clockwise rotation
      void OutStressTensorToVTK( const std::string& filename, double trend, double scale_factor ) const;
  
  private:
     InSituStress();
     /// from the internally stored stresses
     STRESS_REGIME  CalculateStressRegime() const;
     double  subsurface_depth_;
     double  overburden_pressure_; ///< vertical integral over the dry density of the rock
     double  rock_density_;        ///< density at location of sample, for interpolation in vicinity
     double  Sv_, SH_, Sh_,        ///< vertical and horizontal (Andersonian) stresses; SH_ > Sh_
             trend_azimuth_;       ///< trend of sigma1 or SHmax recorded from 0..180 (clockwise on geological compass)
     Point<3U> model_xyz_;
     mutable StressRotate  stress_state_;
     STRESS_REGIME stress_regime_;
 };




// inline functions

// confining pressure = overburden stress
inline double InSituStress::P_conf() const
 { return overburden_pressure_; }
 
// distance from the seafloor / earth surface (that can be different from the y-coordinate)
inline double InSituStress::Depth() const
 { return subsurface_depth_; }
 
// setting the distance from the seafloor / earth surface (that can be different from the y-coordinate)
inline void InSituStress::Depth( double tvd )
 { subsurface_depth_ = tvd; }

// vertical stress normalized by overburden pressure (either of sigma1 to 3 depending on stress regime)
inline double InSituStress::Sv() const
 { return Sv_; }
 
// maximum horizontal stress, normalized by overburden stress
inline double InSituStress::SH() const
 { return SH_; }
 
// minimum horizontal stress, normalized by overburden stress
inline double InSituStress::Sh() const
 { return Sh_; }
 
// azimuth (horizontal XZ plane) of sigma1 or SH expressed in degrees (N-S=0o), assuming Andersonian stresses
inline double InSituStress::Trend() const
 { return trend_azimuth_; }

// density of the rock around the site of stress measurement (kg/m3)
inline double InSituStress::Density() const
 { return rock_density_; }

// output of compact representation of stress state
inline StressRotate& InSituStress::StressState() const
 { return stress_state_; }

// the stress regime inferred from the magnitudes of the principal stresses
inline STRESS_REGIME InSituStress::StressRegime() const
 { return stress_regime_; }


} // end csmp

#endif /* IN_SITU_STRESS_H */
