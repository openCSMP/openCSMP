//
//  OtwayCRC3_RockTypes.h
//  CSMP_CO2GeoSequestrationSimulator
//
//  Created by Stephan Matthai on 23/9/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_OTWAY_CRC3_ROCK_TYPES_H
#define CSMP_OTWAY_CRC3_ROCK_TYPES_H

#include <tuple>
#include "CSMP_definitions.h"

// HARDCODED ROCK PROPERTIES FOR SIMULATION OF CRC3-CRC2 CROSS SECTION
// based on data compiled by Maartje Boon (23/9/2019)
// SKM 24/9/19 - replaced layer thicknesses with volume fractions

namespace csmp {

/**
  Using the "tuple" container of different (rock) types,  a mapping is created between the rocktype numbers sampled from Otway CRC3 and
  the names of the rock so that they can be seen during debugging.
  
  The properties of the rocktypes are hard-coded into the corresponding classes CRC3_RockType1 to 15.
  
  The HeterogeneityAndRateAwareModel gets the rock types from the 'rocktype' variable and choses the corresponding sets of curve fits.
  
  For any of the composite rocktypes there are lamination parallel and perpendicular saturation functions. These are averaged dependent on the flow direction.
  
  For composites, only the drainage curves are dependend on flow rate, the imbibition ones are not.
  
  The imbibition curves are at the capillary limit.

  Some of the composites are treated as mixtures between impermeable nonporous carbonate cemented rocks
  and single rock types. In these cases conventional relperms are used for the remaining porous subtypes, 
  together with a downscaled porosity, and an anisotropic permeability.
  These are:

    M_CbSst_Mst      Massive bedding (Carbonate-cemented sandstone lamina and mudstone matrix) - ONLY ONE SUBTYPE FLOWS
    M_CbSst_Slt        Massive bedding (Carbonate-cemented sandstone lamina and siltstone matrix) - ONLY ONE SUBTYPE FLOWS
    M_CbSst_FSs      Massive bedding (Carbonate-cemented sandstone lamina and fine sandstone matrix) - ONLY ONE SUBTYPE FLOWS
    M_CbSst_CSst     Massive bedding (Carbonate-cemented sandstone lamina and coarse sandstone matrix) - ONLY ONE SUBTYPE FLOWS

 No flow will occur through:

    H_CbSst,      // -- Homogeneous (Carbonate-cemented sandstone) -- NO_FLOW

 which leaves the 6 composites:

    P_Mst_Slt,         4 COMPOSITE -- Planar bedding (Mudstone lamina and siltstone matrix)
    P_Mst_FSst,      6 COMPOSITE: Planar bedding (Mudstone lamina and fine sandstone matrix)
    P_Mst_CSst,   // 8 COMPOSITE: Planar bedding (Mudstone lamina and coarse sandstone matrix)
    P_Slt_FSst,   // 11 COMPOSITE: MAARTJE SAMPLE# - planar bedding (Siltstone lamina with fine sandstone matrix)
    P_Slt_CSst,   // 12 COMPOSITE: Planar bedding (Siltstone lamina with coarse sandstone matrix)
    X_CSst_FSst,  // 14 COMPOSITE: Cross bedding (Fine sandstone lamina and coarse sandstone matrix)

 Then there are 4 rocktypes for which we can use standard VG / BC relperms:

    H_Mst,        // -- Homogeneous Mudstone
    H_Slt,        // -- Homogeneous Siltstone
    H_FSst,       // -- Homogeneous (Fine Sandstone)
    H_CSst        // -- 15 Homogeneous (Coarse sandstone)

  Those composite rocktypes where one of the components is carbonate cemented, can treated as single types with anisotropic perm and high irreducible sw
  
    @attention the averaging is not done inside of this class.

  Extra rock types to define the laminations within the composites?

*/
enum class RockName : int8_t
  {
    // from Achyut & Kuncho's file (July 2019), mostly composites
    WELL,         // 0 by default
    H_Mst,        // 1 Homogeneous Mudstone
    M_CbSst_Mst,  // 2 Massive bedding (Carbonate-cemented sandstone lamina and mudstone matrix) - ONLY ONE SUBTYPE FLOWS
    H_CbSst,      // 3 Homogeneous (Carbonate-cemented sandstone) -- NO_FLOW
    P_Mst_Slt,    // 4 COMPOSITE -- Planar bedding (Mudstone lamina and siltstone matrix)
    M_CbSst_Slt,  // 5 Massive bedding (Carbonate-cemented sandstone lamina and siltstone matrix) - ONLY ONE SUBTYPE FLOWS
    P_Mst_FSst,   // 6 COMPOSITE: Planar bedding (Mudstone lamina and fine sandstone matrix)
    M_CbSst_FSst, // 7 Massive bedding (Carbonate-cemented sandstone lamina and fine sandstone matrix) - ONLY ONE SUBTYPE FLOWS
    P_Mst_CSst,   // 8 COMPOSITE: Planar bedding (Mudstone lamina and coarse sandstone matrix)
    M_CbSst_CSst, // 9 Massive bedding (Carbonate-cemented sandstone lamina and coarse sandstone matrix) - ONLY ONE SUBTYPE FLOWS
    H_Slt,        // 10 Homogeneous Siltstone
    P_Slt_FSst,   // 11 COMPOSITE: MAARTJE SAMPLE# - planar bedding (Siltstone lamina with fine sandstone matrix)
    P_Slt_CSst,   // 12 COMPOSITE: Planar bedding (Siltstone lamina with coarse sandstone matrix)
    H_FSst,       // 11 Homogeneous (Fine Sandstone)
    X_CSst_FSst,  // 14 COMPOSITE: Cross bedding (Fine sandstone lamina and coarse sandstone matrix)
    H_CSst,       // 15 Homogeneous (Coarse sandstone)
    Baffle
  };
  
inline std::string parseRockType( long type ) {
   if ( type == 0 ) return "WELL";
   if ( type == 1 ) return "H_Mst";
   if ( type == 2 ) return "M_CbSst_Mst";
   if ( type == 3 ) return "H_CbSst";
   if ( type == 4 ) return "P_Mst_Slt";
   if ( type == 5 ) return "M_CbSst_Slt";
   if ( type == 6 ) return "P_Mst_FSst";
   if ( type == 7 ) return "M_CbSst_FSst";
   if ( type == 8 ) return "P_Mst_CSst";
   if ( type == 9 ) return "M_CbSst_CSst";
   if ( type == 10 ) return "H_Slt";
   if ( type == 11 ) return "P_Slt_FSst";
   if ( type == 12 ) return "P_Slt_CSst";
   if ( type == 13 ) return "H_FSst";
   if ( type == 14 ) return "X_CSst_FSst";
   if ( type == 15 ) return "H_CSst";
   if ( type == 16 ) return "BAFFLE";
   return "UNDEFINED";
 }

// UTILITY FUNCTIONS
// functions used repeatedly in the curve fitting for layer-parallel and perpendicular flows

/// effective wetting phase saturation
inline double seff( double Sw, double Swr ) {
     return (Sw - Swr) / (1. - Swr);
  } 

/// effective wetting phase saturation taking into account residual saturation of the gas phase
inline double seff( double Sw, double Swr, double Sgr ) {
     return (Sw - Swr) / (1. - Swr - Sgr);
  } 

/// effective wetting phase saturation, limited to range 0..1
inline double seffL( double Sw, double Swr ) {
     return std::max( std::min( (Sw - Swr) / (1. - Swr), 1. ), 0. );
  } 

/// effective wetting phase saturation taking into account residual saturation of the gas phase, limited to range 0..1
inline double seffL( double Sw, double Swr, double Sgr ) {
     return std::max( std::min( (Sw - Swr) / (1. - Swr - Sgr), 1. ), 0. );
  } 

/// repetitive square product involving effective water saturation
inline double seffL_Product( double Sw, double Swr ) {
     const double Seff = seffL(Sw,Swr);
     return (1. - Seff * Seff) * (1 - Seff) * (1 - Seff);
  }

/// repetitive square product involving effective water saturation
inline double seffL_Product( double Sw, double Swr, double Sgr ) {
     const double Seff = seffL(Sw,Swr,Sgr);
     return (1. - Seff * Seff) * (1 - Seff) * (1 - Seff);
  }

/// capillary limit Van Genuchten for wetting phase; the exponents are averaged
inline double krw_VG( double Sw, double m_ave ) {
    const double t1 = std::sqrt(Sw);
    const double t3 = std::pow(Sw, 0.1e1 / m_ave);
    const double t5 = std::pow(0.1e1 - t3, m_ave);
    const double t7 = std::pow(0.1e1 - t5, 0.2e1);
    return t7 * t1;
 }
  
/// Brooks-Corey, capillary limit approximation
inline double krn_BC( double Sw ) {
     const double temp = 1. - Sw * Sw;
     return temp * temp;
  }
  
/// pc(sw) model, Maartje 20/12/19
inline double pc_VG( double Sw, double pd, double m, double Swi_pc ) {
     const double Sw_star = seff(Sw,Swi_pc), Pc_MAX(1.0e+7);
     const double t2 = std::pow(Sw_star, -0.1e1 / m);
     const double t5 = std::pow(t2 - 0.1e1, 0.1e1 - m);
     return std::min( t5 * pd, Pc_MAX );
  }

/// Brooks-Corey lambda parameter from VG m parameter, Lenhard et al. (1989)
inline double lambdaFrom_VG( double m ) {
     assert( m > 0. );
     // Lenhard model does not work; used curve fit instead
     // const double Sweff_ref(0.8);
     // return  (m / (1. - m)) * (1. - std::pow( std::min(Sweff_ref,0.99), 1./m ));
     return 0.5e-1 * std::exp(5.8 * m) + 1.;
  }


/// log10 slope for logarithmic extension of pc following Webb (2000), Pc0(sw=0)=PC_MAX, Pc_star is Pc(Sw_star) = tangent poin where extension meets standard pc curve
inline double pc_Slope( double Pc0, double Pc_star, double Sw_star ) {
     const double t1 = log10(Pc0);
     const double t2 = std::log10(Pc_star);
     return -0.1e1 / Sw_star * (t1 - t2);
  }


/// logarithmic extension of PC following Webb (2000); Sw=actual water saturation, Pc0(sw=0)=PC_MAX, Pc_star is Pc(Sw_star) = tangent poin where log extension starts
inline double pc_LogExtension( double Sw, double Sw_star, double Pc_star ) {
     const double MAX_CAPILLARY_PRESSURE(4e7),
                    slope_pc(pc_Slope(MAX_CAPILLARY_PRESSURE,Pc_star,Sw_star)),
                    t3 = std::log10(Pc_star);
     return std::pow(0.10e2, slope_pc * (Sw - Sw_star) + t3);
  }

/**
    Brooks-Corey capillary pressure correlation for drainage of a water wet medium.
    
    Piecewise definition over entire saturation range:
    BC function down to effective water saturation, Seff= 1%, linear slope below Seff=0.01.
*/
inline double pc_BC( double Sw, double swr, double pd, double bcp ) {
     assert( bcp > 0. );
     assert( Sw >= 0. );
     const double Seff = seff(Sw,swr);
     if ( Seff >= 0.01 )
       return pd * std::pow( Seff, -1. / bcp ); 
       
     // linear extension from Sw to Sw=0 at pc_max
     const double pc01     = pd * std::pow( 0.01, -1. / bcp ),
                    sw_star  = -0.01 * swr + 0.01 + swr;
                    
     return pc_LogExtension( Sw, sw_star, pc01 );
  }


/// curve fitting polynomial with 2 coefficients used by Maartje for water relperms
inline double polyC2( double Sw, double C1, double C2 ) {
     if ( Sw > 1. ) return 1.;
     return std::sqrt(Sw) * (1. - std::pow( 1. - std::pow(Sw,C1), C2));
  }

/// curve fitting polynomial with 3 coefficients used by Maartje for CO2 relperms
inline double polyC3( double Sw, double C1, double C2, double C3 ) {
     if ( Sw > 1. ) return 0.;
     return C1 * std::pow(1. - Sw, C2) * (1. - std::pow(Sw, C3));
  }

inline double poly_abc( double ux, double a, double b, double c ) {
     return a * std::pow( ux, b ) + c;
  }
  
  

// ROCK-TYPE SPECIFICATIONS

/// WELL - rocktype representing a sand-face completion of the well
struct CRC3_RockType0 {
  /// this function establishes whether kri is rate- and flow-direction dependent
  bool IsComposite() const { return false; }

  /// Linear relative permeability: water
  double Krw( double Sw ) const {
      Sw = std::max( std::min(Sw,1.), 0. );
      return Sw;
   }
  /// Linear relative permeability: CO2 
  double Krn( double Sw ) const {
       Sw = std::max( std::min(Sw,1.), 0. );
       return 1. - Sw;
    }
  /// capillary pressure that does not depend on Sw but on radius of well completion  
  double Pc( double /* Sw */ ) const {
       const double IFT = 0.035; // interfacial tension water/CO2 (N/m) 
std::cerr << IFT / diameter_ << "\n";
       return IFT / diameter_;
    }

  const std::string name = "well";        
  const int      rocktype_ = 0;
  const int      subtypes_ = 1;
  const double diameter_ = 0.2,
                 k_        = 1.0e-7,  
                 phi_      = 1.,
                 Swi_pc_   = 0.,
                 Swi_      = 0.,
                 Sgr_      = 0.,
                 m_        = 0.1,
                 pd_       = Pc(1.);  // of the sandface section
};





/** 
     homogeneous mudstone - no special provisions
 */
struct CRC3_RockType1 {
  bool IsComposite() const { return false; }
  /// Water relative permeability
  double Krw( double Sw ) const {
      if ( Sw < Swi_ ) return 0.;
      const double C1(11.73), C2(0.3316);
      return polyC2( Sw, C1, C2 );
   }
  /// CO2 relative permeability:
  double Krn( double Sw ) const {
       if ( Sw < Swi_ ) return 1.; // otherwise Krn will have a value sighltly above zero 
       const double C1(2.848), C2(2.042), C3(3.892);
       return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
    }
    
  const std::string name = "H-Mst";        
  const int      rocktype_ = 1;       // integer code
  const int      subtypes_ = 1;       // out of how many petrotypes the rock consists
  const double k_      = 1.73e-15,  // permeability
                 phi_    = 0.176,     // porosity
                 Swi_    = 0.39,      // irreducible water saturation during drainage
                 Swi_pc_ = 0.39,      // irreducible saturation for pc calc; same as Swi for noncomposites
                 Sgr_    = 0.341,     // max residual gas saturation during layer-parallel imbibition (Maartje MS Word file)
                 m_      = 0.4,       // van Genuchten exponent
                 pd_     = 5000.;     // capillary (drainage) entry pressure of low and high
};








///  homogeneous carbonate-cemented sandstone
struct CRC3_RockType3 {
  bool IsComposite() const { return false; }
  /// Water relative permeability
  double Krw( double Sw ) const {
      if ( Sw < Swi_ ) return 0.;
      const double C1(16.6), C2(0.3374);
      return polyC2( Sw, C1, C2 );
   }
  /// CO2 relative permeability:
  double Krn( double Sw ) const {
       if ( Sw < Swi_ ) return 1.;
       const double C1(5.458), C2(2.051), C3(5.558);
       return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
    }
  const std::string name = "H-CbSst";  
  const int      rocktype_ = 3;
  const int      subtypes_ = 1;
  const double k_        = 1.54e-15, // permeability
                 phi_      = 0.176,
                 Swi_      = 0.552,    // irreducible saturation
                 Swi_pc_   = 0.552,    // irreducible saturation for pc calc; same as Swi for noncomposites
                 Sgr_      = 0.25,
                 m_        = 0.4,      // van Genuchten exponent
                 pd_       = 10000.;   // capillary (drainage) entry pressure of low and high
};






/// Homogeneous (Carbonate-cemented sandstone):
struct CRC3_RockType10 {
  bool IsComposite() const { return false; }
  /// Water relative permeability
     double Krw( double Sw ) const {
     if ( Sw < Swi_ ) return 0.;
     const double C1(6.337), C2(0.4387);
     return polyC2( Sw, C1, C2 );
  }
  /// CO2 relative permeability:
  double Krn( double Sw ) const {
      if ( Sw < Swi_ ) return 1.;
      const double C1(1.522), C2(2.027), C3(2.695);
      return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
   }
  const std::string name = "H-Slt";  
  const int      rocktype_ = 10;
  const int      subtypes_ = 1;
  const double k_        = 3.48e-14, // permeability
                 phi_      = 0.19, 
                 Swi_      = 0.18,     // irreducible saturation
                 Swi_pc_   = 0.18,     // irreducible saturation for pc calc; same as Swi for noncomposites
                 Sgr_      = 0.458,
                 m_        = 0.5,      // van Genuchten exponent
                 pd_       = 3000.;    // capillary (drainage) entry pressure of low and high
};




struct CRC3_RockType13 {
   bool IsComposite() const { return false; }
   /// Water relative permeability
   double Krw( double Sw ) const {
       if ( Sw < Swi_ ) return 0.;
       const double C1(4.936), C2(0.5563);
       return polyC2( Sw, C1, C2 );
    }
   /// CO2 relative permeability:
   double Krn( double Sw ) const {
        if ( Sw < Swi_ ) return 1.;
        const double C1(1.442), C2(2.022), C3(2.594);
        return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
     }
  const std::string name = "H-FSlt";
  const int      rocktype_ = 13;
  const int      subtypes_ = 1;
  const double k_        = 3.61e-13, // permeability
                 phi_      = 0.28, 
                 Swi_      = 0.159,    // irreducible saturation
                 Swi_pc_   = 0.159,    // irreducible saturation for pc calc; same as Swi for noncomposites
                 Sgr_      = 0.341,
                 m_        = 0.6,      // van Genuchten exponent
                 pd_       = 1000.;    // capillary (drainage) entry pressure of low and high
};




/// coarse sandstone
struct CRC3_RockType15 {
  bool IsComposite() const { return false; }
  /// Water relative permeability
  double Krw( double Sw ) const {
     if ( Sw < Swi_ ) return 0.;
     const double C1(3.755), C2(0.6705);
     return polyC2( Sw, C1, C2 );
  }
  /// CO2 relative permeability
  double Krn( double Sw ) const {
      if ( Sw < Swi_ ) return 1.;
      const double C1(1.26), C2(2.012), C3(2.362);
      return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
   }
  const std::string name = "H-CSst"; ///< carbonate-cemented sandstones
  const int      rocktype_ = 15;
  const int      subtypes_ = 1;
  const double k_        = 2.50e-12,  // permeability
                 phi_      = 0.286,     // porosity 
                 Swi_      = 0.104,     // irreducible saturation
                 Swi_pc_   = 0.104,     // irreducible saturation for pc calc; same as Swi for noncomposites
                 Sgr_      = 0.353,
                 m_        = 0.7,       // van Genuchten exponent
                 pd_       = 750.;      // capillary (drainage) entry pressure of low and high
};



/**
    Essentially impermeable rock with an entry pressure that is so high that it will not be exceeded
    by a plume with 100-m column height.
*/
struct CRC3_RockType16 {
   bool IsComposite() const { return false; }
   /// Linear water relative permeability for drainage
   double Krw( double Sw ) const { 
        if ( Sw < Swi_ ) return 0.;
        return seffL(Sw,Swi_);       
     }
   /// Linear CO2 relative permeability for drainage
   double Krn( double Sw ) const { 
        if ( Sw < Swi_ ) return 0.8;
        return 1. - seffL(Sw,Swi_,Sgr_);  
     }

   const std::string name = "baffle"; 
   const int      rocktype_ = 16;
   const int      subtypes_ = 1;
   const double k_        = 1.0e-16,  // permeability
                  phi_      = 0.09,     // porosity
                  Swi_      = 0.3,      // irreducible saturation
                  Swi_pc_   = 0.3,      // irreducible saturation for pc calc; same as Swi for noncomposites
                  Sgr_      = 0.1,      // guestimate
                  m_        = 0.7,      // van Genuchten exponent
                  pd_       = 1.0e+6;   // capillary (drainage) entry pressure
};



// ==============================================

// COMPOSITE ROCK TYPES

// ==============================================

const double sw_for_krw_eq1(0.99);


/** 
     Massive bedding: Carbonate-cemented Sandstone - Mudstone
     
     @attention SKM: difference of relperms for vertical flow is ignored here because it is so small.
     
           @TODO this should be made a composite to honour different end-point krn for horizontal and vertical
 */
struct CRC3_RockType2 {
  bool IsComposite() const { return true; }
  
  /// Water relative permeability (CL)
  double Krw_ParallelDrainage( double Sw ) const {
       if ( Sw < Swi_ ) return 0.;
       return krw_VG( seffL(Sw,Swi_), m_ave_ );
    }
  /// CO2 relative permeability:
  double Krn_ParallelDrainage( double Sw ) const {
       if ( Sw < Swi_ ) return 0.37;  
       return std::min( seffL_Product(Sw,Swi_), 0.37 );
    }

  /// Water relative permeability (CL)
  double Krw_CrossDrainage( double Sw ) const {
       if ( Sw < Swi_ ) return 0.;
       return krw_VG( seffL(Sw,Swi_), 0.4054 );
    }
  /// CO2 relative permeability:
  double Krn_CrossDrainage( double Sw ) const {
       if ( Sw < Swi_ ) return 0.22;  
       return std::min( seffL_Product(Sw,Swi_), 0.22 );
    }

  const std::string name = "M-CbSst-Mst";        
  const int      rocktype_ = 2;
  const int      subtypes_ = 1;
  const double k_low_    = 1.539564e-15, k_high_ = 1.727075e-15, // permeability
                 k_        = k_low_ + k_high_ / 2.,
                 phi_     = 0.176,                           // porosity
                 LY_low_  = 0.025,  LY_high_  = 0.025,        // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.552,  Swi_high_ = 0.39,        // irreducible saturations of the 2 different layers
                 Swi_pc_  = 0.472,      // irreducible saturation for pc calc; same as Swi for noncomposites
                 Sgr_     = 0.288,
                 m_low_   = 0.4,    m_high_ = 0.4,           // van Genuchten exponents for the 2 different layers
                 m_       = 0.4,   // average value for c(sw) calculation 
                 pd_low_  = 10000., pd_high_ = 5000.,        // capillary (drainage) entry pressure of low and high
                 pd_      = 5000.,
                 dPc_     = 5000., // UNSPECIFIED: guess of pc-difference between high_k and low_k layer at swc
                 Swi_     = 0.472,
                 m_ave_   = 0.427;
};



///< Planar bedding mudstone - silt
struct CRC3_RockType4 {
   bool IsComposite() const { return true; }

   /// Water relative permeability
   double Krw_ParallelDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= 6e-3 ) return 0.61 * pow( seffL(Sw,Swi_), 4.06 );
       if ( ux <= 1e-5 ) return 0.626 * pow( seffL(Sw,Swi_), 5.11 );
       const double C1_krw = poly_abc(ux,0.2485,0.01877,0.3717);
       const double a2 =  0.1286,
                      b2 = -0.2088,
                      c2 =  3.691,
                      C2_krw = a2 * std::pow(ux,b2) + c2;
       
       // Krw_ave (for averaged velocity?)
       return std::min( C1_krw * std::pow( seffL(Sw,Swi_),C2_krw), 1. );
    }
 
   /// CO2 relative permeability:
   double Krn_ParallelDrainage( double Sw, double ux ) const {
      if ( ux >= 6e-3 ) return std::min( 0.323 * seffL_Product(Sw,Swi_), 0.47 );
      if ( ux <= 1e-5 ) {
           if ( Sw >= 1. - 0.0206 ) return 0.; // to avoid capillary breakthrough
           return std::min( 0.488 * seffL_Product(Sw,Swi_), 0.47 );
        }
      const double a_c1 =  0.04757,
                     b_c1 =  -0.1498,
                     c_c1 =   0.2208,
                     C1_krnw = a_c1 * std::pow(ux,b_c1) + c_c1;
        // Krnw_ave=
        return std::min( C1_krnw * seffL_Product(Sw,Swi_), 0.47 );
     }
   
   /// Water relative permeability
   double Krw_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= ux_VL_cross_ ) return 0.581 * std::pow( seffL(Sw,Swi_), 5.26 );
       if ( ux <= ux_CL_cross_ ) return 0.559 * std::pow( seffL(Sw,Swi_), 3.15 );
       const double a1 =  0.2485,
                      b1 =  0.01877,
                      c1 =  0.3717,
                      C1_krw = a1 * std::pow(ux,b1) + c1;

       const double a2 =  0.2485,
                      b2 =  0.01877,
                      c2 =  0.3717,
                      C2_krw = a2 * std::pow(ux,b2) + c2;
       // Krw_ave=
       return C1_krw * std::pow( seffL(Sw,Swi_),C2_krw);
    }
   
   /// CO2 relative permeability:
   double Krn_CrossDrainage( double Sw, double ux ) const {
        // making sure that for the VL case Krn>Sgr = 0
        if ( Sw >= 1. - 0.0153 ) return 0.;
        if ( ux > ux_VL_cross_ ) ux = ux_VL_cross_; // return std::max( 0.267 * seffL_Product(Sw,Swi_,0.0153), 0. );
        if ( ux < ux_CL_cross_ ) ux = ux_CL_cross_; // return std::max( 0.0794 * seffL_Product(Sw,Swi_,0.0206), 0. );
        const double Sgr_ux = poly_abc(ux,0.0075237,-0.05069,0.002388); // between 1-2%, not sure what this means for drainage
        const double C1_krn = poly_abc(ux,3.612,0.3236,0.06067);
        // Krnw_ave=
        return std::min( std::max( C1_krn * seffL_Product(Sw,Swi_,Sgr_ux), 0. ), 0.34 );
     }
   
  // ROCK PROPERTIES 
  const std::string name = "P-Mst-Slt";  
  const int      rocktype_ = 4;
  const int      subtypes_ = 2;
  const double k_low_   = 1.73e-15, k_high_ = 3.48e-14,      // layer permeabilities
                 phi_     = 0.183, 
                 LY_low_  = 0.25,  LY_high_ = 0.25,   // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.39,  Swi_high_ = 0.18,  // irreducible saturations of the 2 different layers
                 Swi_     = 0.562,
                 Swi_pc_  = 0.281, // for pc calculation
                 Sgr_     = 0.41,  // max Sgr for layer parallel imbibition (Maartje curves)
                 m_low_   = 0.4,   m_high_ = 0.5,     // van Genuchten exponents for the 2 different layers
                 m_       = 0.431, // for pc calculation
                 pd_low_  = 5000., pd_high_ = 3000.,  // capillary (drainage) entry pressure of low and high
                 pd_      = 3000.,
                 ux_VL_cross_ = 1.0e-4,
                 ux_CL_cross_ = 1.0e-7,
                 dPc_     = 2000.; // UNSPECIFIED: guess of pc difference between high_k and low_k layer at swc

}; // end CRC3_RockType4







/**  
    homogeneous carbonate-cemented sandstone (originally perceived as a composite, but treated as single rocktype)
    
    @attention although this is a composite with a layered K structure and a directionally dependent permeability, rate dependence is ignored
    
    @note we have a different residual saturation dependent on the flow direction
*/
struct CRC3_RockType5 {
  bool IsComposite() const { return true; } 
  /// Water relative permeability
  double Krw_ParallelDrainage( double Sw ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       const double m_ave(0.493), Swi(0.5);
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability:
  double Krn_ParallelDrainage( double Sw ) const {
       const double Swi(0.5);       
       return std::min( seffL_Product(Sw,Swi), 0.56 );
    }

  /// Water relative permeability (capillary limit)
  double Krw_CrossDrainage( double Sw ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       const double m_ave(0.636), Swi(0.359);
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability (capillary limit)
  double Krn_CrossDrainage( double Sw ) const {
       const double Swi(0.359);       
       return std::min( seffL_Product(Sw,Swi), 0.22 );
    }

  const std::string name = "M-CbSst-Slt";  
  const int      rocktype_ = 5;
  const int      subtypes_ = 1;
  const double k_low_= 1.54e-15, k_high_=3.48e-14, 
                 phi_      = 0.183, 
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.18,
                 Swi_pc_   = 0.359, 
                 Sgr_      = 0.353,       
                 m_low_    = 0.4,      
                 m_high_   = 0.5,
                 m_        = 0.415,      
                 pd_low_   = 10000.,    
                 pd_high_  = 3000.,
                 pd_       = 3000.,
                 dPc_      = 0.;    
};






/// composite - planar bedding mudstone - fine sandstone
struct CRC3_RockType6 {
  bool IsComposite() const { return true; }
  
  double Krw_ParallelDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       // viscous limit approximation
       if ( ux >= 1.0e-2 ) return 0.728 * pow( seffL(Sw,Swi_), 3.47 );
       // capillary limit approximation
       if ( ux <= 1.0e-5 ) return 0.75  * pow( seffL(Sw,Swi_), 5.08 );
      // between CL and VL
       const double a1(0.0002717), b1(-0.3894), c1(0.7265);
       const double C1_krw = a1 * std::pow( ux, b1 ) + c1; 
       const double a2(2.2), b2(-0.06136), c2(0.6185);
       const double C2_krw = a2 * std::pow( ux, b2) + c2;
       // Krw_ave 
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }
    
   double Krn_ParallelDrainage( double Sw, double ux ) const {
        if ( ux >= 1.0e-2 ) return std::min( 0.4  * seffL_Product(Sw,Swi_), 0.61 );
        if ( ux <= 1.0e-5 ) return std::min( 0.87 * seffL_Product(Sw,Swi_), 0.61 );
        const double C1_krn   = poly_abc( ux, 3.904, -0.01475, -3.756 ),
                       Krnw_ave = C1_krn * seffL_Product(Sw,Swi_);
        // limiting value to >=0
        return std::min( std::max( Krnw_ave, 0. ), 0.61 );
     }

  double Krw_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       // viscous limit approximation
       if ( ux >= 7.0e-4 ) return 0.557 * pow( seffL(Sw,Swi_), 5.74 );
       // capillary limit approximation
       if ( ux <= 5.0e-7 ) return 0.834 * pow( seffL(Sw,Swi_), 1.57 );
      // between CL and VL
       const double a1(0.2151), b1(-0.07629), c1(0.1834);
       const double C1_krw = a1 * std::pow( ux, b1 ) + c1; 
       const double a2(-7.324), b2(-0.04686), c2(16.02);
       const double C2_krw = a2 * std::pow( ux, b2) + c2;
       // Krw_ave 
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }
    
   double Krn_CrossDrainage( double Sw, double ux ) const {
        if ( Sw >= 1. - 0.00397 ) return 0.;
        if ( ux >= 7.0e-4 ) return std::min( std::max( 0.2014 * seffL_Product(Sw,Swi_,0.00397), 0. ), 0.33 );
        if ( ux <= 5.0e-7 ) return std::min( std::max( 0.018  * seffL_Product(Sw,Swi_,0.1154), 0. ), 0.33 );
        const double Sgr_ux = poly_abc( ux, 0.08238, -0.0783, -0.1412 ),
                       C1_krn = poly_abc( ux, 0.7342, 0.1051, -0.1418 ),
                       Krnw_ave = C1_krn * seffL_Product(Sw,Swi_,Sgr_ux);
        // limiting value to >=0
        return std::min( std::max( Krnw_ave, 0. ), 0.33 );
     }

  const std::string name = "P-Mst-FSst";  
  const int      rocktype_ = 6;
  const int      subtypes_ = 2;
  const double k_low_=1.73e-15, k_high_=3.61e-13, 
                 phi_  = 0.228,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.159, 
                 m_low_=0.4, m_high_=0.6,        // van Genuchten exponents for the 2 different layers
                 m_      = 0.408,
                 pd_low_=5000., pd_high_=1000.,  // capillary (drainage) entry pressure of low and high
                 pd_     = 1000.,
                 dPc_    = 4000., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 Swi_    = 0.496,
                 Swi_pc_ = 0.248,
                 Sgr_     = 0.327;    
};







///  massive carbonate cemented Sandstone – Fine Sandstone (composite, but not rate dependent)
struct CRC3_RockType7 {
  bool IsComposite() const { return true; }
  /// Water relative permeability
  double Krw_ParallelDrainage( double Sw ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       const double m_ave(0.6), Swi(0.47);
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability:
  double Krn_ParallelDrainage( double Sw ) const {
       const double Swi(0.47);
       return std::min( seffL_Product(Sw,Swi), 0.68 );
    }

  /// Water relative permeability
  double Krw_CrossDrainage( double Sw ) const {
      if ( Sw >= sw_for_krw_eq1 ) return 1.;
      const double m_ave(0.593), Swi(0.46);
      return krw_VG( seffL(Sw,Swi), m_ave );
   }
  
  /// CO2 relative permeability:
  double Krn_CrossDrainage( double Sw ) const {
       const double Swi(0.46);
       return std::min( seffL_Product(Sw,Swi), 0.33 );
    }

  const std::string name = "M-CbSst-FSst";  
  const int      rocktype_ = 7;
  const int      subtypes_ = 1;
  const double k_low_    = 1.54e-15, 
                 k_high_   = 3.61e-13, 
                 phi_      = 0.228,
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.18,  
                 Swi_pc_   = 0.311,
                 Sgr_      = 0.298,      
                 m_low_    = 0.4,      
                 m_high_   = 0.6,  
                 m_        = 0.401,    
                 pd_low_   = 10000.,    
                 pd_high_  = 1000.,
                 pd_       = 1000.,
                 dPc_      = 9000.; // TODO: pre-compute    
};







/// composite - Planar bedding Mudstone – Coarse Sandstone - @todo: appears to be vertically layered?
struct CRC3_RockType8 {
  bool IsComposite() const { return true; }
  
  double Krw_ParallelDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux <= 5.0e-5 ) return 0.848 * std::pow( seffL(Sw,Swi_), 4.475 );
       if ( ux >= 4.0e-2 ) return 0.824 * std::pow( seffL(Sw,Swi_), 3.176 );
       // between CL and VL
       const double C1_krw = poly_abc( ux, 0.0001377, -0.5243, 0.8232 );
       const double C2_krw = poly_abc( ux, 0.775, -0.114, 2.078 );
       return C1_krw * pow( seffL(Sw,Swi_), C2_krw );
    }

  double Krn_ParallelDrainage( double Sw, double ux ) const {
       if ( ux <= 5.0e-5 ) return std::min( 1.07 * seffL_Product(Sw, Swi_), 0.67 );
       if ( ux >= 4.0e-2 ) return std::min( 0.53 * seffL_Product(Sw, Swi_), 0.67 );
       // between CL and VL
       const double C1_krn = poly_abc( ux, 0.755, -0.06682, -0.3927 );
       // limiting range between 0. and 0.43
       return std::max( std::min( C1_krn * seffL_Product(Sw,Swi_), 0.67 ), 0. );
    }

  double Krw_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= 2e-4 ) return 0.492 * std::pow( seffL(Sw,Swi_), 4.19 );
       if ( ux <= 5e-7 ) return std::min( 2.98  * std::pow( seffL(Sw,Swi_), 1.66 ), 0.4921 );
       // between CL and VL
       const double C1_krw = poly_abc( ux, 2.446, -0.0728, -4.056 );
       const double C2_krw = poly_abc( ux, -133.2, -0.003056, 140.9 );
       return std::min( C1_krw * pow( seffL(Sw,Swi_), C2_krw ), 0.4921 );
    }

  double Krn_CrossDrainage( double Sw, double ux ) const {
       if ( ux >= 2e-4 ) return std::min( std::max( 0.1361 * seffL_Product(Sw, Swi_,0.0658), 0. ), 0.33 );
       if ( ux <= 5e-7 ) return std::min( std::max( 0.0134 * seffL_Product(Sw, Swi_,0.427), 0. ), 0.33 );
       // between CL and VL
       const double Sgr_ux = poly_abc( ux, 1.687, -0.02636, -2.046 ),
                      C1_krn = poly_abc( ux, 0.9455, 0.1964, -0.04131 );
       // limiting range >= 0. 
       return std::min( std::max( C1_krn * seffL_Product(Sw,Swi_,Sgr_ux), 0. ), 0.33 );
    }

  const std::string name = "P-Mst-CSst";  
  const int      rocktype_ = 8;
  const int      subtypes_ = 2;
  const double k_low_=1.73e-15, k_high_=2.50e-12,
                 phi_  = 0.231, 
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.104, 
                 Swi_pc_=0.213,
                 Sgr_    = 0.333,    
                 m_low_=0.4, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 m_=0.399,
                 pd_low_=5000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 pd_=750.,
                 dPc_=4250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 Swi_=0.426;
};








///  Massive bedding: Carbonate cemented Sandstone – Coarse Sandstone (not rate dependent, but composite)
struct CRC3_RockType9 {
  bool IsComposite() const { return true; }
  
  /// Water relative permeability
  double Krw_ParallelDrainage( double Sw ) const {
      if ( Sw >= sw_for_krw_eq1 ) return 1.;
      const double m_ave(0.7), Swi(0.44);
      return krw_VG( seffL(Sw,Swi), m_ave );
   }
  
  /// CO2 relative permeability:
  double Krn_ParallelDrainage( double Sw ) const {
       const double Swi(0.44);
       if ( Sw < Swi ) return 0.67;
       return std::min( seffL_Product(Sw,Swi), 0.67 );
    }

  /// Water relative permeability
  double Krw_CrossDrainage( double Sw ) const {
      if ( Sw >= sw_for_krw_eq1 ) return 1.;
      const double m_ave(4.29), Swi(0.4 ); // really high m_ave, just to get curvefit
      return krw_VG( seffL(Sw,Swi), m_ave );
   }
  
  /// CO2 relative permeability:
  double Krn_CrossDrainage( double Sw ) const {
       const double Swi(0.275);
       if ( Sw < Swi ) return 0.33;
       return std::min( seffL_Product(Sw,Swi,0.55), 0.33 );
    }

  const std::string name = "M-CbSst-CSst";  
  const int      rocktype_ = 9;
  const int      subtypes_ = 1;
  const double k_low_    = 1.54e-15, 
                 k_high_   = 2.5e-12, 
                 phi_      = 0.231,   
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.104, 
                 Swi_pc_   = 0.275,
                 Sgr_      = 0.299,    // drainage    
                 m_low_    = 0.4,      
                 m_high_   = 0.7,
                 m_        = 0.388,      
                 pd_low_   = 10000.,    
                 pd_high_  = 750.,
                 pd_       = 750.,
                 dPc_      = 9250.; // TODO: pre-compute    
};









/// composite - planar bedding fine sandstone - silt
struct CRC3_RockType11 {
 bool IsComposite() const { return true; }

 /// Water relative permeability
 double Krw_ParallelDrainage( double Sw, double ux ) const {
     if ( Sw >= sw_for_krw_eq1 ) return 1.;
     if ( ux >= 3.0e-2 ) return 0.731 * std::pow( seffL(Sw,Swi_), 4.17 ); // VL;
     if ( ux <= 1.0e-5 ) return 0.739 * std::pow( seffL(Sw,Swi_), 5.95 ); // CL;     
     const double C1_krw = poly_abc( ux, 0.00273, -0.1325, 0.7263 );
     const double C2_krw = poly_abc( ux, 49.75, -0.004283, -46.31 );
     // Krw_ave
     return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
  }
 
 /// CO2 relative permeability:
 double Krn_ParallelDrainage( double Sw, double ux  ) const {
     if ( ux >= 3.0e-2 ) return std::min( 0.545 * seffL_Product(Sw,Swi_), 0.74 );
     if ( ux <= 1.0e-5 ) return std::min( 1.33 * seffL_Product(Sw,Swi_), 0.74 );    
     const double C1_krn = poly_abc( ux, -30.39, 0.003259, 30.6 );
     // Krnw_ave=
     const double Krnw_ave = C1_krn * seffL_Product(Sw,Swi_);
     return std::min( std::max( Krnw_ave, 0. ), 0.74 );
  }

  /// Water relative permeability
  double Krw_CrossDrainage( double Sw, double ux ) const {
      if ( Sw >= sw_for_krw_eq1 ) return 1.;
      if ( ux >= 3.0e-2 ) return 0.6197 * std::pow( seffL(Sw,Swi_), 4.74 ); // VL;
      if ( ux <= 1.0e-5 ) return std::min( 0.993 * std::pow( seffL(Sw,Swi_), 3.26 ), 0.6197 ); // CL;     
      const double C1_krw = poly_abc( ux, 0.0066, -0.2977, 0.5443 );
      const double C2_krw = poly_abc( ux, -0.7936, -0.1001, 6.542 );
      // Krw_ave
      return std::min( C1_krw * std::pow( seffL(Sw,Swi_), C2_krw ), 0.6197 );
   }
  
  /// CO2 relative permeability:
  double Krn_CrossDrainage( double Sw, double ux  ) const {
      if ( Sw >= 1. - 0.037 ) return 0.;
      if ( ux >= 3.0e-2 ) return std::min( std::max( 1.33 * seffL_Product(Sw,Swi_,0.037), 0. ), 0.46 ); 
      if ( ux <= 1.0e-5 ) return std::min( std::max( 0.0416 * seffL_Product(Sw,Swi_,0.159), 0. ), 0.46 );  
      const double Sgr_ux = poly_abc( ux, 0.003424, -0.2677, 0.006397 ),
                     C1_krn = poly_abc( ux, 3.153, 0.1897, -0.1728 );
      // Krnw_ave=
      const double Krnw_ave = C1_krn * seffL_Product(Sw,Swi_,Sgr_ux);
      return std::min( std::max( Krnw_ave, 0. ), 0.46 );
   }

  const std::string name = "P-FSst-Slt";  
  // Composite1 - hypothetical sample based on Achyut's rocktypes
  const int      rocktype_ = 11;
  const int      subtypes_ = 2;
  const double k_low_   = 3.48e-14,
                 k_high_  = 3.61e-13,  // layer permeabilities
                 phi_     = 0.235,
                 LY_low_  = 0.25,  LY_high_ = 0.25,    // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.18,  Swi_high_ = 0.159, // irreducible saturations of the 2 different layers
                 Swi_pc_  = 0.167,
                 Sgr_     = 0.376,    
                 m_low_   = 0.5,   m_high_ = 0.6,       // van Genuchten exponents for the 2 different layers
                 m_       = 0.49,
                 pd_low_  = 3000., pd_high_ = 1000.,  // capillary (drainage) entry pressure of low and high
                 pd_      = 1000.,
                 dPc_     = 2000., // was: 1.3061e+06, pc difference between high_k and low_k layer at swc
                 Swi_     = 0.335; // for both cross-layer and layer parallel flow
};







/// composite - Planar bedding Siltstone – Coarse Sandstone 
struct CRC3_RockType12 {
  bool IsComposite() const { return true; }

  double Krw_ParallelDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= ux_VL_ ) return 0.83  * std::pow( seffL(Sw,Swi_), 3.44 );
       if ( ux <= ux_CL_ ) return 0.869 * std::pow( seffL(Sw,Swi_), 6.14 );
       const double C1_krw = poly_abc( ux, 0.002649, -0.2458, 0.8242 );
       const double C2_krw = poly_abc( ux, -29.82, 0.01156, 32.25 );
       // Krw_ave =
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }

  double Krn_ParallelDrainage( double Sw, double ux ) const {
       if ( ux >= ux_VL_ ) return std::min( 0.637 * seffL_Product(Sw,Swi_), 0.43 );
       if ( ux <= ux_CL_ ) return std::min( 2.094 * seffL_Product(Sw,Swi_), 0.43 );
       // NB: Sgr=0 for parallel flow
       const double C1_krn = poly_abc( ux, -6.277, 0.03536, 6.239 );
       const double Krn_ave = C1_krn * seffL_Product(Sw,Swi_);
       return std::min( Krn_ave, 0.43 );
    }

  double Krw_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= ux_cross_VL_ ) return 0.544966  * std::pow( seffL(Sw,Swi_), 4.437 );
       if ( ux <= ux_cross_CL_ ) return std::min( 2.55 * std::pow( seffL(Sw,Swi_), 3.066 ), 0.545 );
       const double C1_krw = poly_abc( ux, 2.036, -0.06943, -3.025 );
       const double C2_krw = poly_abc( ux, 6.683, 0.1167, 1.837 );
       // Krw_ave =
       return std::min( C1_krw * std::pow( seffL(Sw,Swi_), C2_krw ), 0.545 );
    }

  double Krn_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= 1. - 0.0463 ) return 0.;
       if ( ux >= ux_cross_VL_ ) return std::min( std::max( 0.3541 * seffL_Product(Sw,Swi_,0.0463), 0. ), 0.45 );
       if ( ux <= ux_cross_CL_ ) return std::min( std::max( 0.0239 * seffL_Product(Sw,Swi_,0.327), 0. ), 0.45 );
       const double Sgr_ux = poly_abc( ux, -1.888, 0.03386, 1.482 ),
                      C1_krn = poly_abc( ux, 3.627, 0.2739, -0.04431 );
       const double Krn_ave = C1_krn * seffL_Product(Sw,Swi_,Sgr_ux);
       return std::min( std::max( Krn_ave, 0. ), 0.45 );
    }

  const std::string name = "P-Slt-CSst";  
  const int      rocktype_ = 12;
  const int      subtypes_ =  2;
  const double k_low_=3.48e-14, k_high_=2.50e-12, 
                 phi_  = 0.238,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.18, Swi_high_=0.104, 
                 Swi_pc_ =0.134,
                 m_low_=0.5, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 m_=0.394,
                 pd_low_=3000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 pd_=750.,
                 dPc_=2250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 5.0e-2, ux_CL_ = 1e-5,
                 ux_cross_VL_ = 3.0e-4, ux_cross_CL_ = 5e-7,
                 Swi_ = 0.269;
};






/// composite - cross-bedded coarse - fine laminated sandstone
struct CRC3_RockType14 {
  bool IsComposite() const { return true; }

  double Krw_ParallelDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= ux_VL_ ) return 0.82 * std::pow( seffL(Sw,Swi_), 3.69 );
       if ( ux <= ux_CL_ ) return 0.86 * std::pow( seffL(Sw,Swi_), 4.64 );
       const double C1_krw = poly_abc( ux, 0.008275, -0.1924, 0.811 );
       const double C2_krw = poly_abc( ux, 0.4083, -0.1338, 3.239 );
       // Krw_ave =
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }

  double Krn_ParallelDrainage( double Sw, double ux ) const {
       if ( ux >= ux_VL_ ) return std::min( 0.68 * seffL_Product(Sw,Swi_), 0.69 );
       if ( ux <= ux_CL_ ) return std::min( 1.07 * seffL_Product(Sw,Swi_), 0.69 );
       const double C1_krn = poly_abc( ux, 0.2568, -0.1039, 0.4035 );
       const double Krn_ave = C1_krn * seffL_Product(Sw,Swi_);
       return std::min( Krn_ave, 0.69 );
    }

  double Krw_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( ux >= ux_cross_VL_ ) return 0.756 * std::pow( seffL(Sw,Swi_), 4.22 );
       if ( ux <= ux_cross_CL_ ) return std::min( 0.715 * std::pow( seffL(Sw,Swi_), 2.944 ), 0.7557 );
       const double C1_krw = poly_abc( ux, 0.8557, 0.5453, 0.7148 );
       const double C2_krw = poly_abc( ux, 4.545, 0.1916, 2.662 );
       // Krw_ave =
       return std::min( C1_krw * std::pow( seffL(Sw,Swi_), C2_krw ), 0.7557 );
    }

  double Krn_CrossDrainage( double Sw, double ux ) const {
       if ( Sw >= 1. - 0.0516 ) return 0.;
       if ( ux >= ux_cross_VL_ ) return std::min( std::max( 0.793 * seffL_Product(Sw,Swi_,0.436), 0. ), 0.57 );
       if ( ux <= ux_cross_CL_ ) return std::min( std::max( 0.316 * seffL_Product(Sw,Swi_,0.0516), 0. ), 0.57 );
       const double Sgr_ux = poly_abc( ux, 0.0175, -0.03539, 0.02231 ),
                      C1_krn = poly_abc( ux, 3.329, 0.3396, 0.2929 );
       const double Krn_ave = C1_krn * seffL_Product(Sw,Swi_,Sgr_ux);
       return std::min( std::max( Krn_ave, 0. ), 0.57 );
    }

  const std::string name = "X-CSst-FSst";  
  const int      rocktype_ = 14;
  const int      subtypes_ =  2;
  const double k_low_=3.61e-13, k_high_=2.50e-12, 
                 phi_  = 0.283,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.159, Swi_high_=0.104, 
                 Swi_pc_=0.131,
                 Sgr_   = 0.33,    
                 m_low_=0.6, m_high_=0.7,    
                 m_=0.581,    
                 pd_low_=1000., pd_high_=750., 
                 pd_=750., 
                 dPc_=250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 0.5, ux_CL_ = 1.0e-4,
                 ux_cross_VL_ = 4.0e-3, ux_cross_CL_ = 5.0e-7,
                 Swi_ = 0.262;
};



/// Tuple container to return rocktype by number (std::get<1>(rocktypes_);
struct OtwayRockTypes {
    bool IsComposite( int RT ) const {
         if ( RT == 2 || RT == 4 || RT == 5 || RT == 6 || 
              RT == 7 || RT == 8 || RT == 9 || RT == 11 || RT == 12 || RT == 14 ) return true;
         return false;
      }
    
    std::tuple<CRC3_RockType0, // well
               CRC3_RockType1,CRC3_RockType2,CRC3_RockType3,
               CRC3_RockType4,CRC3_RockType5,CRC3_RockType6,
               CRC3_RockType7,CRC3_RockType8,CRC3_RockType9,
               CRC3_RockType10,CRC3_RockType11,CRC3_RockType12,
               CRC3_RockType13,CRC3_RockType14,CRC3_RockType15,
               CRC3_RockType16>  rocktype_;  
};










} // end csmp


/** ADDITIONAL DOCUMENTATION

 Anisotropy in absolute permeability
 
 In the model absolute permeability is a scalar value. We choose this scalar value to be the permeability in the horizontal direction. Anisotropy in the absolute permeability is incorporated through the relative permeability as described below.
 
 Kave_hor = the average cell permeability in the horizontal direction
 Kave_ver = the average cell permeability in the vertical direction direction
 
 The water permeability in the vertical direction is given by Kwater=Kave_ver * Krw
 The CO2 permeability in the vertical direction is  given by KCO2=Kave_ver * Krnw
 
 In the model the absolute permeability is chosen to be Kave_hor. To arrive at the correct phase permeabilities in the vertical direction we therefore have to multiply the relative permeability in the vertical direction by Kave_ver/Kave_hor:
 
 Kwater=Kave_hor*(Kave_ver/Kave_hor)* Krw = Kave_ver * Krw
 KCO2= Kave_hor*(Kave_ver/Kave_hor)* Krnw = Kave_ver * Krnw
 
 For the laminated sand/silt this is  35.22mD/365.53mD = 0.0964. So all of the above curve fittings need to be multiplied by 0.0964 to incorporate anisotropy into the model.
 
 */


#endif /* CSMP_OTWAY_CRC3_ROCK_TYPES_H */
