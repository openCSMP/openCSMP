//
//  InSituStress.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/18/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#include "InSituStress.h"

namespace csmp {

/// converts the enum to a string that can be printed to stdio or file
std::string  parseStressRegime( STRESS_REGIME stress_regime )
  {
     if ( stress_regime == ISOSTATIC )        return std::string("ISOSTATIC");
     if ( stress_regime == NORMAL_FAULTING )  return std::string("NORMAL_FAULTING");
     if ( stress_regime == TRANSTENSIONAL )   return std::string("TRANSTENSIONAL");
     if ( stress_regime == TRANSPRESSIONAL )  return std::string("TRANSPRESSIONAL");
     if ( stress_regime == REVERSE_FAULTING ) return std::string("REVERSE_FAULTING");
     return std::string("undefined");
  }


/**  Constructing the stress measurement for a normal or strike-slip faulting regime
    
    StressRotate takes the vectors y_axis(0,1,0) and x_axis(1,0,0) for initialization of sigma1 and sigma3.
    Later this can be modified by rotation. Here the trend is the clockwise principal stress azimuth 
    from 1-180o (Clar compass).
    
    @parameter depth should match the location where the stresses were measured in the reservoir represented by the model

*/
InSituStress::InSituStress(  double depth, double rdensity, double lload, double Sv, double SH, double Sh, double trend )
 : subsurface_depth_(depth), rock_density_(rdensity),
   // since the rotation in StressRotate is counter-clockwise we add 180 and then subtract the trend
   overburden_pressure_(lload), Sv_(Sv), SH_(SH), Sh_(Sh), trend_azimuth_(180. - trend),
   //             should be sigma1  and      sigma3
   stress_state_( Point<3U>(0,1,0), Point<3U>(1,0,0), Sv*lload, SH*lload, Sh*lload ),
   // assuming the sample location is at the origin of the model at the user-defined depth
   model_xyz_(0.,depth,0.),
   stress_regime_(NORMAL_FAULTING)
 {
    stress_regime_ = CalculateStressRegime();
    assert( trend_azimuth_ <= 180. );
    assert( stress_regime_ == NORMAL_FAULTING  or  stress_regime_ == TRANSTENSIONAL );
    // TODO: implement case of a reverse faulting regime
    Out();
 }

// copy constructor
InSituStress::InSituStress( const InSituStress& iss )
 : subsurface_depth_(iss.subsurface_depth_),
   rock_density_(iss.rock_density_),
   overburden_pressure_(iss.overburden_pressure_),
   Sv_(iss.Sv_), SH_(iss.SH_), Sh_(iss.Sh_),
   trend_azimuth_(iss.trend_azimuth_),
   model_xyz_(iss.model_xyz_),
   stress_state_(iss.stress_state_),
   stress_regime_(iss.stress_regime_)
 {
 }





STRESS_REGIME  InSituStress::CalculateStressRegime() const
 {
    // starting with te most common case
    if ( Sv_ >= SH_ and SH_ >= Sh_ ) return NORMAL_FAULTING;
    if ( SH_ > Sv_ and Sh_ > Sv_ ) return REVERSE_FAULTING;
    if ( SH_ > Sv_ and Sv_ > Sh_ ) {
         if ( SH_ > Sv_ and Sv_ > 1. ) return TRANSPRESSIONAL;
         if ( SH_ > Sv_ and Sv_ > 1. ) return TRANSTENSIONAL;
      }
    return ISOSTATIC;
   
 }  // end CalculateStressRegime





/**
   sets the grid coordinate to the value where the stress measurement was taken.
   
   @attention this may create a mismatch with the subsurface_depth_  stored in InSituStress.
*/
void InSituStress::SampleLocation( double model_X_is_east, double model_Y_is_elevation, double model_Z_is_south )
 {
    model_xyz_.Set( model_X_is_east, model_Y_is_elevation, model_Z_is_south );
 }
 
 
/// returns sample location
Point<3U> const&  InSituStress::SampleLocation() const
 { return model_xyz_; }
 


/**  CartesianStress

   Computation of full stress tensor rotated clockwise according to match its compass trend 
   while assuming that its plunge is vertical (Andersonian stresses).
   
   Before the calculation, the principal stresses are adjusted to take into account
   the vertical offset of the location of interest from that of the stress
   measurement.
   
   @param stress_offset  difference in overburden stress as compared with sample location.
   @param ts  cartesian stress tensor into which the resulting stress state is returned.
   
*/
void InSituStress::CartesianStress( TensorVariable<3U>& ts, double stress_offset ) const
 {
     assert( trend_azimuth_ >= 0. );
     assert( trend_azimuth_ <= 180. );
     stress_state_.Rotate( 'y', 180. - trend_azimuth_ );
     stress_state_.CartesianStressTensor( ts, stress_offset );
     stress_state_.Reset();
 }



/// outputs stress tensor with the supplied rotation as vectors placed at origin of grid 
void InSituStress::OutStressTensorToVTK(  const std::string& filename ) const
 {
     stress_state_.OutputToVTK( filename, model_xyz_ );
 }


/// outputs stress tensor at grid origin, rotating it clockwise (like in compass) by a maximum of 180 degrees aroud vertical axis 
void InSituStress::OutStressTensorToVTK(  const std::string& filename, double trend, double scale_factor ) const
 {
     assert( trend >= 0. );
     assert( trend <= 180. );
     stress_state_.Reset();
     stress_state_.Rotate( 'y', 180. - trend );
     stress_state_.OutputToVTK( filename, model_xyz_, scale_factor );
     stress_state_.Reset();
 }



void InSituStress::Out() const
 {
    std::cout <<"\n\nInSituStress::Out:  current state (at the site of the stress measurement): ";
    /// distance from the seafloor / earth surface (that can be different from the y-coordinate)
    std::cout <<"\n\tsub-surface depth:\t\t\t"<< Depth();
    /// confining pressure = overburden stress
    std::cout <<"\n\toverburden pressure:\t\t\t"<< P_conf();
    /// azimuth of sigma1 or SH expressed in degrees (N-S=0o), assuming Andersonian stresses 
    std::cout <<"\n\tstress azimuth:\t\t\t"<< Trend();
    /// vertical stress normalized by overburden pressure (either of sigma1 to 3 depending on stress regime)
    std::cout <<"\n\tnormalized vertical stress:\t\t\t"<< Sv();
    /// maximum horizontal stress, normalized by overburden stress
    std::cout <<"\n\tSH-max, maximum horizontal stress:\t"<< SH();
    /// minimum horizontal stress, normalized by overburden stress
    std::cout <<"\n\tSh-min, minimum horizontal stress:\t"<< Sh();
    /// density of the rock around the site of stress measurement (kg/m3)
    std::cout <<"\n\tambient rock density:\t\t\t"<< Density();
    /// the stress regime inferred from the magnitudes of the principal stresses
    std::cout <<"\n\tstress regime:\t\t\t"<< parseStressRegime(StressRegime()) << std::endl;
 }


} // end csmp
