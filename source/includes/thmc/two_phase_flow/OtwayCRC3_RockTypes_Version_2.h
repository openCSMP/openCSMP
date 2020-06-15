//
//  OtwayCRC3_RockTypes.h
//  CSMP_CO2GeoSequestrationSimulator
//
//  Created by Stephan Matthai on 23/9/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_OTWAY_CRC3_ROCK_TYPES_VERSION_2_H
#define CSMP_OTWAY_CRC3_ROCK_TYPES_VERSION_2_H

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
// functions used repeatedly in the curve fitting for layer-parallel and perpendicular flows that are produced by Maartje

/// r_w used in Chierici model (AbdAllah addition)
inline double64 r_w( double64 Sw, double64 Swr ) {
     assert(Sw != 1.);
     return (Sw - Swr) / (1. - Sw);
  } 

/// gamma used in Chierici model (AbdAllah addition)
inline double64 gamma( double64 ux, double64 a, double64 b, double64 c, double64 d ) {
     return d + (a / (1. + b * std::pow( ux, c )));
  }   
 
/// Wetting phase relative permeability (Chierici model) (AbdAllah addition)
inline double64 k_w( double64 rw, double64 gamma_W_1, double64 gamma_W_2 ) {
     return std::exp( -1. * gamma_W_1 * std::pow( rw , -1. * gamma_W_2) );
  }
  
/// Non-wetting phase relative permeability (Chierici model) (AbdAllah addition)
inline double64 k_nw( double64 rw, double64 gamma_NW_1, double64 gamma_NW_2 ) {
     return std::exp( -1. * gamma_NW_1 * std::pow( rw , gamma_NW_2) );
  }  

  

/// effective wetting phase saturation
inline double64 seff( double64 Sw, double64 Swr ) {
     return (Sw - Swr) / (1. - Swr);
  } 

/// effective wetting phase saturation taking into account residual saturation of the gas phase
inline double64 seff( double64 Sw, double64 Swr, double64 Sgr ) {
     return (Sw - Swr) / (1. - Swr - Sgr);
  } 

/// effective wetting phase saturation, limited to range 0..1
inline double64 seffL( double64 Sw, double64 Swr ) {
     return std::max( std::min( (Sw - Swr) / (1. - Swr), 1. ), 0. );
  } 

/// effective wetting phase saturation taking into account residual saturation of the gas phase, limited to range 0..1
inline double64 seffL( double64 Sw, double64 Swr, double64 Sgr ) {
     return std::max( std::min( (Sw - Swr) / (1. - Swr - Sgr), 1. ), 0. );
  } 

/// repetitive square product involving effective water saturation
inline double64 seffL_Product( double64 Sw, double64 Swr ) {
     const double64 Seff = seffL(Sw,Swr);
     return (1. - Seff * Seff) * (1 - Seff) * (1 - Seff);
  }

/// repetitive square product involving effective water saturation
inline double64 seffL_Product( double64 Sw, double64 Swr, double64 Sgr ) {
     const double64 Seff = seffL(Sw,Swr,Sgr);
     return (1. - Seff * Seff) * (1 - Seff) * (1 - Seff);
  }

/// capillary limit Van Genuchten for wetting phase; the exponents are averaged
inline double64 krw_VG( double64 Sw, double64 m_ave ) {
    const double64 t1 = std::sqrt(Sw);
    const double64 t3 = std::pow(Sw, 0.1e1 / m_ave);
    const double64 t5 = std::pow(0.1e1 - t3, m_ave);
    const double64 t7 = std::pow(0.1e1 - t5, 0.2e1);
    return t7 * t1;
 }
  
/// Brooks-Corey, capillary limit approximation
inline double64 krn_BC( double64 Sw ) {
     const double64 temp = 1. - Sw * Sw;
     return temp * temp;
  }
    

/// Brooks-Corey lambda parameter from VG m parameter, Lenhard et al. (1989)
inline double64 lambdaFrom_VG( double64 m ) {
     assert( m > 0. );
     // Lenhard model does not work; used curve fit instead
     // const double64 Sweff_ref(0.8);
     // return  (m / (1. - m)) * (1. - std::pow( std::min(Sweff_ref,0.99), 1./m ));
     return 0.5e-1 * std::exp(5.8 * m) + 1.;
  }


/// log10 slope for logarithmic extension of pc following Webb (2000), Pc0(sw=0)=PC_MAX, Pc_star is Pc(Sw_star) = tangent poin where extension meets standard pc curve
inline double pc_Slope( double64 Pc0, double64 Pc_star, double64 Sw_star ) {
     const double64 t1 = log10(Pc0);
     const double64 t2 = std::log10(Pc_star);
     return -0.1e1 / Sw_star * (t1 - t2);
  }


/// logarithmic extension of PC following Webb (2000); Sw=actual water saturation, Pc0(sw=0)=PC_MAX, Pc_star is Pc(Sw_star) = tangent poin where log extension starts
inline double pc_LogExtension( double64 Sw, double64 Sw_star, double64 Pc_star ) {
     const double64 MAX_CAPILLARY_PRESSURE(4e7),
                    slope_pc(pc_Slope(MAX_CAPILLARY_PRESSURE,Pc_star,Sw_star)),
                    t3 = std::log10(Pc_star);
     return std::pow(0.10e2, slope_pc * (Sw - Sw_star) + t3);
  }

/**
    Brooks-Corey capillary pressure correlation for drainage of a water wet medium.
    
    Piecewise definition over entire saturation range:
    BC function down to effective water saturation, Seff= 1%, linear slope below Seff=0.01.
*/
inline double64 pc_BC( double64 Sw, double swr, double64 pd, double64 bcp ) {
     assert( bcp > 0. );
     assert( Sw >= 0. );
     const double64 Seff = seff(Sw,swr);
     if ( Seff >= 0.01 )
       return pd * std::pow( Seff, -1. / bcp ); 
       
     // linear extension from Sw to Sw=0 at pc_max
     const double64 pc01     = pd * std::pow( 0.01, -1. / bcp ),
                    sw_star  = -0.01 * swr + 0.01 + swr;
                    
     return pc_LogExtension( Sw, sw_star, pc01 );
  }
  
  
  
/// pc(sw) model, Maartje 20/12/19
inline double64 pc_VG( double64 Sw, double64 pd, double64 m, double64 Swi_pc ) {

     const double64 Sw_star = seff(Sw,Swi_pc), Pc_MAX(1.0e+7);
     const double64 t2 = std::pow(Sw_star, -0.1e1 / m);
     const double64 t5 = std::pow(t2 - 0.1e1, 0.1e1 - m);
     return std::min( t5 * pd, Pc_MAX );     
               
  }
    


/// curve fitting polynomial with 2 coefficients used by Maartje for water relperms
inline double64 polyC2( double64 Sw, double64 C1, double64 C2 ) {
     if ( Sw > 1. ) return 1.;
     return std::sqrt(Sw) * (1. - std::pow( 1. - std::pow(Sw,C1), C2));
  }
   
/// curve fitting polynomial with 3 coefficients used by Maartje for CO2 relperms
inline double64 polyC3( double64 Sw, double64 C1, double64 C2, double64 C3 ) {
     if ( Sw > 1. ) return 0.;
     return C1 * std::pow(1. - Sw, C2) * (1. - std::pow(Sw, C3));
  }

inline double64 poly_abc( double64 ux, double64 a, double64 b, double64 c ) {
     return a * std::pow( ux, b ) + c;
  }
  
  

// ROCK-TYPE SPECIFICATIONS

/// WELL - rocktype representing a sand-face completion of the well
struct CRC3_RockType0 {
  /// this function establishes whether kri is rate- and flow-direction dependent
  bool IsComposite() const { return false; }

  /// Linear relative permeability: water
  double64 Krw( double64 Sw ) const {
      Sw = std::max( std::min(Sw,1.), 0. );
      return Sw;
   }
  /// Linear relative permeability: CO2 
  double64 Krn( double64 Sw ) const {
       Sw = std::max( std::min(Sw,1.), 0. );
       return 1. - Sw;
    }
  /// capillary pressure that does not depend on Sw but on radius of well completion  
  double64 Pc( double64 /* Sw */ ) const {
       const double64 IFT = 0.035; // interfacial tension water/CO2 (N/m) 
std::cerr << IFT / diameter_ << "\n";
       return IFT / diameter_;
    }

  const std::string name = "well";        
  const int      rocktype_ = 0;
  const int      subtypes_ = 1;
  const double64 diameter_ = 0.2,
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
  double64 Krw( double64 Sw ) const {
      if ( Sw < Swi_ ) return 0.;
      const double64 C1(11.73), C2(0.3316);
      return polyC2( Sw, C1, C2 );
   }
  /// CO2 relative permeability:
  double64 Krn( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.; // otherwise Krn will have a value sighltly above zero 
       const double64 C1(2.848), C2(2.042), C3(3.892);
       return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
    }
    
  const std::string name = "H-Mst";        
  const int      rocktype_ = 1;       // integer code
  const int      subtypes_ = 1;       // out of how many petrotypes the rock consists
  const double64 k_      = 1.73e-15,  // permeability
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
  double64 Krw( double64 Sw ) const {
      if ( Sw < Swi_ ) return 0.;
      const double64 C1(16.6), C2(0.3374);
      return polyC2( Sw, C1, C2 );
   }
  /// CO2 relative permeability:
  double64 Krn( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.;
       const double64 C1(5.458), C2(2.051), C3(5.558);
       return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
    }
  const std::string name = "H-CbSst";  
  const int      rocktype_ = 3;
  const int      subtypes_ = 1;
  const double64 k_        = 1.54e-15, // permeability
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
     double64 Krw( double64 Sw ) const {
     if ( Sw < Swi_ ) return 0.;
     const double64 C1(6.337), C2(0.4387);
     return polyC2( Sw, C1, C2 );
  }
  /// CO2 relative permeability:
  double64 Krn( double64 Sw ) const {
      if ( Sw < Swi_ ) return 1.;
      const double64 C1(1.522), C2(2.027), C3(2.695);
      return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
   }
  const std::string name = "H-Slt";  
  const int      rocktype_ = 10;
  const int      subtypes_ = 1;
  const double64 k_        = 3.48e-14, // permeability
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
   double64 Krw( double64 Sw ) const {
       if ( Sw < Swi_ ) return 0.;
       const double64 C1(4.936), C2(0.5563);
       return polyC2( Sw, C1, C2 );
    }
   /// CO2 relative permeability:
   double64 Krn( double64 Sw ) const {
        if ( Sw < Swi_ ) return 1.;
        const double64 C1(1.442), C2(2.022), C3(2.594);
        return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
     }
  const std::string name = "H-FSlt";
  const int      rocktype_ = 13;
  const int      subtypes_ = 1;
  const double64 k_        = 3.61e-13, // permeability
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
  double64 Krw( double64 Sw ) const {
     if ( Sw < Swi_ ) return 0.;
     const double64 C1(3.755), C2(0.6705);
     return polyC2( Sw, C1, C2 );
  }
  /// CO2 relative permeability
  double64 Krn( double64 Sw ) const {
      if ( Sw < Swi_ ) return 1.;
      const double64 C1(1.26), C2(2.012), C3(2.362);
      return std::min( polyC3( Sw, C1, C2, C3 ), 1. );
   }
  const std::string name = "H-CSst"; ///< carbonate-cemented sandstones
  const int      rocktype_ = 15;
  const int      subtypes_ = 1;
  const double64 k_        = 2.50e-12,  // permeability
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
   double64 Krw( double64 Sw ) const { 
        if ( Sw < Swi_ ) return 0.;
        return seffL(Sw,Swi_);       
     }
   /// Linear CO2 relative permeability for drainage
   double64 Krn( double64 Sw ) const { 
        //if ( Sw < Swi_ ) return 0.8;
        if ( Sw < Swi_ ) return 1.; 
        return 1. - seffL(Sw,Swi_,Sgr_);  
     }

   const std::string name = "baffle"; 
   const int      rocktype_ = 16;
   const int      subtypes_ = 1;
   const double64 k_        = 1.0e-16,  // permeability
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

const double64 sw_for_krw_eq1(0.99);


/** 
     Massive bedding: Carbonate-cemented Sandstone - Mudstone
     
     @attention SKM: difference of relperms for vertical flow is ignored here because it is so small.
     
           @TODO this should be made a composite to honour different end-point krn for horizontal and vertical
 */
struct CRC3_RockType2 {
  bool IsComposite() const { return true; }
  
  /// Water relative permeability (CL)
  double64 Krw_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 0.;
       return krw_VG( seffL(Sw,Swi_), m_ave_ );
    }
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.;  
       return std::min( seffL_Product(Sw,Swi_), 1. );
    }

  /// Water relative permeability (CL)
  double64 Krw_CrossDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 0.;
       return krw_VG( seffL(Sw,Swi_), 0.4054 );
    }
  /// CO2 relative permeability:
  double64 Krn_CrossDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.;  
       return std::min( seffL_Product(Sw,Swi_), 1. );
    }

  const std::string name = "M-CbSst-Mst";        
  const int      rocktype_ = 2;
  const int      subtypes_ = 1;
  const double64 k_low_    = 1.539564e-15, k_high_ = 1.727075e-15, // permeability
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

  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
       const double64 C1_krw = gamma( ux, aw1_, bw1_, cw1_,dw1_ ); 
       const double64 C2_krw = 0.6;
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  );
       return std::min( std::max( Krw_ave, 0. ), 1. ); 
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.; 
       if ( ux >= ux_VL_ ) ux = ux_VL_;;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krnw = gamma( ux, anw1_, bnw1_, cnw1_,dnw1_ ); 
	   const double64 C2_krnw = 0.7;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
       if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krw = gamma( ux, aw1_cross, bw1_cross, cw1_cross,dw1_cross ); 
	   const double64 C2_krw = 0.7;
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  ); 
       return std::min( std::max( Krw_ave, 0. ), 1. );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
       if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krnw = gamma( ux, anw1_cross, bnw1_cross, cnw1_cross,dnw1_cross ); 
	   const double64 C2_krnw = 0.6;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }
	
  double64 Pd( double64 ux ) const {
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       return gamma( ux, apd_, bpd_, cpd_,dpd_ );	   
  }	
   
  // ROCK PROPERTIES 
  const std::string name = "P-Mst-Slt";  
  const int      rocktype_ = 4;
  const int      subtypes_ = 2;
  const double64 k_low_   = 1.73e-15, k_high_ = 3.48e-14,      // layer permeabilities
                 phi_     = 0.183, 
                 LY_low_  = 0.25,  LY_high_ = 0.25,   // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.39,  Swi_high_ = 0.18,  // irreducible saturations of the 2 different layers
                 Swi_pc_  = 0.285, // for pc calculation
                 Sgr_     = 0.41,  // max Sgr for layer parallel imbibition (Maartje curves)
                 m_low_   = 0.4,   m_high_ = 0.5,     // van Genuchten exponents for the 2 different layers
                 m_       = 0.445, // for pc calculation
		 apd_ = 550., bpd_ = 0.001903616, cpd_ = -0.58, dpd_ = 3300.,   // used for rate dependent pd
                 pd_low_  = 5000., pd_high_ = 3000.,  // capillary (drainage) entry pressure of low and high
                 pd_      = 3000.,
                 dPc_     = 2000., // UNSPECIFIED: guess of pc difference between high_k and low_k layer at swc
                 ux_VL_ = 9.00e-03, ux_CL_ = 2.20e-06,   // Threshold velocity for VL and CL of parallel flow 
                 ux_cross_VL_ = 3.00e-03, ux_cross_CL_ = 1.00e-07,   // Threshold velocity for VL and CL of cross flow 
		 aw1_ = 0.39, bw1_ = 12088.38073, cw1_ = 1., dw1_ = 4.29,   // used for gamma_W_1 (parallel flow)
		 anw1_ = 0.31, bnw1_ = 0.009279014, cnw1_ = -0.65, dnw1_ = 1.58, // used for gamma_NW_1 (parallel flow)	 		 
		 aw1_cross = 4.2525, bw1_cross = 3.6414e-06, cw1_cross = -1., dw1_cross = 3.45,  // used for gamma_W_1 (cross flow)
		 anw1_cross = 0.86, bnw1_cross = 11489.62807, cnw1_cross = 1., dnw1_cross = 2.07, // used for gamma_NW_1 (cross flow)		 
                 Swi_ = 0.38,           // for parallel flow for kr calculation
		 Swi_cross = 0.3458;      // for cross flow for kr calculation		 

}; // end CRC3_RockType4







/**  
    homogeneous carbonate-cemented sandstone (originally perceived as a composite, but treated as single rocktype)
    
    @attention although this is a composite with a layered K structure and a directionally dependent permeability, rate dependence is ignored
    
    @note we have a different residual saturation dependent on the flow direction
*/
struct CRC3_RockType5 {
  bool IsComposite() const { return true; } 
  /// Water relative permeability
  double64 Krw_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 0.;	  
       return krw_VG( seffL(Sw,Swi_), m_ave_ );
    }
  
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.;  
       return std::min( seffL_Product(Sw,Swi_), 1. );
    }

  /// Water relative permeability (capillary limit)
  double64 Krw_CrossDrainage( double64 Sw ) const {
	  const double64 Swi(0.359), m_ave(0.636);
      if ( Sw < Swi ) return 0.;	  
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability (capillary limit)
  double64 Krn_CrossDrainage( double64 Sw ) const {
       const double64 Swi(0.359);
	   if ( Sw < Swi ) return 1.;	       
       return std::min( seffL_Product(Sw,Swi), 1. );
    }

  const std::string name = "M-CbSst-Slt";  
  const int      rocktype_ = 5;
  const int      subtypes_ = 1;
  const double64 k_low_= 1.54e-15, k_high_=3.48e-14, 
                 phi_      = 0.183, 
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.18,
		 Swi_ = 0.5,
		 m_ave_ = 0.493,
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
  
  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krw = gamma( ux, aw1_, bw1_, cw1_,dw1_ ); 
	   const double64 C2_krw = 0.66;
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  );
       return std::min( std::max( Krw_ave, 0. ), 1. ); 
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_VL_ ) ux = ux_VL_;;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krnw = gamma( ux, anw1_, bnw1_, cnw1_,dnw1_ ); 
	   const double64 C2_krnw = 0.75;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krw = gamma( ux, aw1_cross, bw1_cross, cw1_cross,dw1_cross ); 
	   const double64 C2_krw = 0.66;
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  ); 
       return std::min( std::max( Krw_ave, 0. ), 1. );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krnw = gamma( ux, anw1_cross, bnw1_cross, cnw1_cross,dnw1_cross ); 
	   const double64 C2_krnw = 0.6;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }
	
  double64 Pd( double64 ux ) const {
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
	   return gamma( ux, apd_, bpd_, cpd_,dpd_ );	   
  }	

  const std::string name = "P-Mst-FSst";  
  const int      rocktype_ = 6;
  const int      subtypes_ = 2;
  const double64 k_low_=1.73e-15, k_high_=3.61e-13, 
                 phi_  = 0.228,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.159, 
                 m_low_=0.4, m_high_=0.6,        // van Genuchten exponents for the 2 different layers
                 m_      = 0.48,
		 apd_ = 1200., bpd_ = 0.000338834, cpd_ = -0.85, dpd_ = 1100.,   // used for rate dependent pd
                 pd_low_=5000., pd_high_=1000.,  // capillary (drainage) entry pressure of low and high
                 pd_     = 3000.,
                 dPc_    = 4000., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 Swi_pc_ = 0.28,
                 ux_VL_ = 9.00e-03, ux_CL_ = 2.20e-06,   // Threshold velocity for VL and CL of parallel flow 
                 ux_cross_VL_ = 8.50e-04, ux_cross_CL_ = 5.00e-08,   // Threshold velocity for VL and CL of cross flow 
		 aw1_ = 0.49, bw1_ = 18769.71602, cw1_ = 1.2, dw1_ = 3.49,   // used for gamma_W_1 (parallel flow)
		 anw1_ = 0.56, bnw1_ = 0.001230912, cnw1_ = -1., dnw1_ = 1.42, // used for gamma_NW_1 (parallel flow)	 		 
		 aw1_cross = 3.58, bw1_cross = 9.9295e-06, cw1_cross = -1., dw1_cross = 1.75,  // used for gamma_W_1 (cross flow)
		 anw1_cross = 4.66, bnw1_cross = 215.0348268, cnw1_cross = 0.58, dnw1_cross = 1.725, // used for gamma_NW_1 (cross flow)		 
                 Swi_ = 0.38,           // for parallel flow for kr calculation
		 Swi_cross = 0.3,      // for cross flow for kr calculation
		 Sgr_     = 0.327; 
};







///  massive carbonate cemented Sandstone – Fine Sandstone (composite, but not rate dependent)
struct CRC3_RockType7 {
  bool IsComposite() const { return true; }
  /// Water relative permeability
  double64 Krw_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 0.;	  
       return krw_VG( seffL(Sw,Swi_), m_ave_ );
    }
  
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.;  
       return std::min( seffL_Product(Sw,Swi_), 1. );
    }

  /// Water relative permeability (capillary limit)
  double64 Krw_CrossDrainage( double64 Sw ) const {
	  const double64 Swi(0.46), m_ave(0.593 );
      if ( Sw < Swi ) return 0.;	  
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability (capillary limit)
  double64 Krn_CrossDrainage( double64 Sw ) const {
       const double64 Swi(0.46);
	   if ( Sw < Swi ) return 1.;	       
       return std::min( seffL_Product(Sw,Swi), 1. );
    }

  const std::string name = "M-CbSst-FSst";  
  const int      rocktype_ = 7;
  const int      subtypes_ = 1;
  const double64 k_low_    = 1.54e-15, 
                 k_high_   = 3.61e-13, 
                 phi_      = 0.228,
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.18,  
		 Swi_ = 0.47,
		 m_ave_ = 0.6,
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
  
  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw >= sw_for_krw_eq1 ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
	   if ( Sw <= Swi_ ) return 0.;
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krw = gamma( ux, aw1_, bw1_, cw1_,dw1_ ); 
	   const double64 C2_krw = 0.68;
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  );
       return std::min( std::max( Krw_ave, 0. ), 1. ); 
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_VL_ ) ux = ux_VL_;;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krnw = gamma( ux, anw1_, bnw1_, cnw1_,dnw1_ ); 
	   const double64 C2_krnw = 0.8;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
       if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krw = gamma( ux, aw1_cross, bw1_cross, cw1_cross,dw1_cross ); 
	   const double64 C2_krw = gamma( ux, aw2_cross, bw2_cross, cw2_cross,dw2_cross );
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  ); 
       return std::min( std::max( Krw_ave, 0. ), 1. );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.; 
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krnw = gamma( ux, anw1_cross, bnw1_cross, cnw1_cross,dnw1_cross ); 
	   const double64 C2_krnw = 0.55;
       const double64 Krn_ave =  k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }
	
  double64 Pd( double64 ux ) const {
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
	   return gamma( ux, apd_, bpd_, cpd_,dpd_ );	   
  }	

  const std::string name = "P-Mst-CSst";  
  const int      rocktype_ = 8;
  const int      subtypes_ = 2;
  const double64 k_low_=1.73e-15, k_high_=2.50e-12,
                 phi_  = 0.231, 
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.104, 
                 Swi_pc_=0.29,
                 Sgr_    = 0.333,    
                 m_low_=0.4, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 m_=0.55,
		 apd_ = 1550., bpd_ = 0.010051836, cpd_ = -0.5, dpd_ = 600.,   // used for rate dependent pd
                 pd_low_=5000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 pd_=750.,
                 dPc_=4250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 1.50e-01, ux_CL_ = 2.20e-05,   // Threshold velocity for VL and CL of parallel flow 
                 ux_cross_VL_ = 5.00e-03, ux_cross_CL_ = 5.00e-08,   // Threshold velocity for VL and CL of cross flow 
		 aw1_ = 0.46, bw1_ = 271.7819254, cw1_ = 0.95, dw1_ = 3.,   // used for gamma_W_1 (parallel flow)
		 anw1_ = 0.52, bnw1_ = 0.002029431, cnw1_ = -1., dnw1_ = 1.28, // used for gamma_NW_1 (parallel flow)	 		 
		 aw1_cross = 2.96, bw1_cross = 1.14055e-06, cw1_cross = -1.244, dw1_cross = 1.2, // used for gamma_W_1 (cross flow)
		 aw2_cross =-0.59, bw2_cross = 8.63605e-07, cw2_cross = -1.311, dw2_cross = 1.11, // used for gamma_W_1 (cross flow)				 
		 anw1_cross = 9.79, bnw1_cross = 371.6676348, cnw1_cross = 0.44, dnw1_cross = 2.27, // used for gamma_NW_1 (cross flow)		 
                 Swi_ = 0.3,           // for parallel flow for kr calculation
		 Swi_cross = 0.27;      // for cross flow for kr calculation
};








///  Massive bedding: Carbonate cemented Sandstone – Coarse Sandstone (not rate dependent, but composite)
struct CRC3_RockType9 {
  bool IsComposite() const { return true; }
  
  /// Water relative permeability
  double64 Krw_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 0.;	  
       return krw_VG( seffL(Sw,Swi_), m_ave_ );
    }
  
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       if ( Sw < Swi_ ) return 1.;  
       return std::min( seffL_Product(Sw,Swi_), 1. );
    }

  /// Water relative permeability (capillary limit)
  double64 Krw_CrossDrainage( double64 Sw ) const {
	  const double64 Swi(0.4), m_ave(4.29 );
      if ( Sw < Swi ) return 0.;	  
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability:
  double64 Krn_CrossDrainage( double64 Sw ) const {
       const double64 Swi(0.275), Sgr(0.55);
       if ( Sw < Swi ) return 1.;
       //if ( Sw > (1-Sgr) ) return 1.;
       if ( Sw > (1-Sgr) ) return 0.;
       return std::min( seffL_Product(Sw,Swi,Sgr), 1. );
    }

  const std::string name = "M-CbSst-CSst";  
  const int      rocktype_ = 9;
  const int      subtypes_ = 1;
  const double64 k_low_    = 1.54e-15, 
                 k_high_   = 2.5e-12, 
                 phi_      = 0.231,   
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.104, 
		 Swi_ = 0.44,
		 m_ave_ = 0.7,
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

  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krw = gamma( ux, aw1_, bw1_, cw1_,dw1_ ); 
	   const double64 C2_krw = gamma( ux, aw2_, bw2_, cw2_,dw2_ );
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  );
       return std::min( std::max( Krw_ave, 0. ), 1. ); 
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.; 
	   if ( ux >= ux_VL_ ) ux = ux_VL_;;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krnw = gamma( ux, anw1_, bnw1_, cnw1_,dnw1_ ); 
	   const double64 C2_krnw = 0.78;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krw = gamma( ux, aw1_cross, bw1_cross, cw1_cross,dw1_cross ); 
	   const double64 C2_krw = gamma( ux, aw2_cross, bw2_cross, cw2_cross,dw2_cross );
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  ); 
       return std::min( std::max( Krw_ave, 0. ), 1. );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krnw = gamma( ux, anw1_cross, bnw1_cross, cnw1_cross,dnw1_cross ); 
	   const double64 C2_krnw = 0.75; 
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }
	
  double64 Pd( double64 ux ) const {
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
	   return gamma( ux, apd_, bpd_, cpd_,dpd_ );	   
  }	

  const std::string name = "P-FSst-Slt";  
  // Composite1 - hypothetical sample based on Achyut's rocktypes
  const int      rocktype_ = 11;
  const int      subtypes_ = 2;
  const double64 k_low_   = 3.48e-14,
                 k_high_  = 3.61e-13,  // layer permeabilities
                 phi_     = 0.235,
                 LY_low_  = 0.25,  LY_high_ = 0.25,    // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.18,  Swi_high_ = 0.159, // irreducible saturations of the 2 different layers
                 Swi_pc_  = 0.18,
                 Sgr_     = 0.376,    
                 m_low_   = 0.5,   m_high_ = 0.6,       // van Genuchten exponents for the 2 different layers
                 m_       = 0.53,
		 apd_ = 495., bpd_ = 0.003388053, cpd_ = -0.65, dpd_ = 1210.,   // used for rate dependent pd
                 pd_low_  = 3000., pd_high_ = 1000.,  // capillary (drainage) entry pressure of low and high
                 pd_      = 1000.,
                 dPc_     = 2000., // was: 1.3061e+06, pc difference between high_k and low_k layer at swc
                 ux_VL_ = 1.00e-01, ux_CL_ = 2.20e-06,   // Threshold velocity for VL and CL of parallel flow 
                 ux_cross_VL_ = 5.00e-03, ux_cross_CL_ = 5.00e-08,   // Threshold velocity for VL and CL of cross flow 
		 aw1_ = 0.38, bw1_ = 270.4264074, cw1_ = 1., dw1_ = 3.6,   // used for gamma_W_1 (parallel flow)
		 aw2_ = 0.043, bw2_ = 1.69541e-05, cw2_ = -1.3, dw2_ = 0.635,  // used for gamma_W_2 (parallel flow)
		 anw1_ = 0.55, bnw1_ = 0.01227734, cnw1_ = -0.8, dnw1_ = 1.25, // used for gamma_NW_1 (parallel flow)	 		 
		 aw1_cross = 0.46, bw1_cross = 7.83504e-13, cw1_cross = -2.5, dw1_cross = 3.69,  // used for gamma_W_1 (cross flow)
		 aw2_cross = -0.63, bw2_cross = 6.43602e-06, cw2_cross = -0.62, dw2_cross = 1.23,  // used for gamma_W_2 (cross flow)
		 anw1_cross = 2., bnw1_cross = 1373.063975, cnw1_cross = 0.88, dnw1_cross = 1.75, // used for gamma_NW_1 (cross flow)		 
                 Swi_ = 0.27,           // for parallel flow for kr calculation
		 Swi_cross = 0.17;      // for cross flow for kr calculation
};







/// composite - Planar bedding Siltstone – Coarse Sandstone 
struct CRC3_RockType12 {
  bool IsComposite() const { return true; }

  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.; 
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krw = gamma( ux, aw1_, bw1_, cw1_,dw1_ ); 
	   const double64 C2_krw = gamma( ux, aw2_, bw2_, cw2_,dw2_ );
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  );
       return std::min( std::max( Krw_ave, 0. ), 1. ); 
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_VL_ ) ux = ux_VL_;;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krnw = gamma( ux, anw1_, bnw1_, cnw1_,dnw1_ ); 
	   const double64 C2_krnw = 0.75;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krw = gamma( ux, aw1_cross, bw1_cross, cw1_cross,dw1_cross ); 
	   const double64 C2_krw = gamma( ux, aw2_cross, bw2_cross, cw2_cross,dw2_cross );
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  ); 
       return std::min( std::max( Krw_ave, 0. ), 1. );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krnw = gamma( ux, anw1_cross, bnw1_cross, cnw1_cross,dnw1_cross ); 
	   const double64 C2_krnw = 0.71; 
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Pd( double64 ux ) const {
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       return gamma( ux, apd_, bpd_, cpd_,dpd_ );	   
  }

  const std::string name = "P-Slt-CSst";  
  const int      rocktype_ = 12;
  const int      subtypes_ =  2;
  const double64 k_low_=3.48e-14, k_high_=2.50e-12, 
                 phi_  = 0.238,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.18, Swi_high_=0.104, 
                 Swi_pc_ =0.165,
                 m_low_=0.5, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 m_=0.58,
		 apd_ = 580, bpd_ = 0.004190229, cpd_ = -0.75, dpd_ = 950.,   // used for rate dependent pd
                 pd_low_=3000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 pd_=750.,
                 dPc_=2250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 1.00e-01, ux_CL_ = 1.70e-05,   // Threshold velocity for VL and CL of parallel flow 
                 ux_cross_VL_ = 5.00e-03, ux_cross_CL_ = 5.00e-08,   // Threshold velocity for VL and CL of cross flow 
		 aw1_ = 0.95, bw1_ = 200.33681, cw1_ = 1., dw1_ = 2.8,   // used for gamma_W_1 (parallel flow)
		 aw2_ = 0.04, bw2_ = 1.69541e-05, cw2_ = -1.3, dw2_ = 0.7,  // used for gamma_W_2 (parallel flow)
		 anw1_ = 0.64, bnw1_ = 0.006737947, cnw1_ = -1., dnw1_ = 1.26, // used for gamma_NW_1 (parallel flow)	 		 		 
		 aw1_cross = 1.819, bw1_cross = 6.9598e-08, cw1_cross = -1.485, dw1_cross = 1.945,  // used for gamma_W_1 (cross flow)
		 aw2_cross = -0.384, bw2_cross = 1.2813e-06, cw2_cross = -1.387, dw2_cross = 0.955,  // used for gamma_W_2 (cross flow)
		 anw1_cross = 3.88, bnw1_cross = 617.3894351, cnw1_cross = 0.71, dnw1_cross = 2.09, // used for gamma_NW_1 (cross flow)		 
                 Swi_ = 0.27,          // for parallel flow for kr calculation
		 Swi_cross = 0.15;      // for cross flow for kr calculation
};






/// composite - cross-bedded coarse - fine laminated sandstone
struct CRC3_RockType14 {
  bool IsComposite() const { return true; }

  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.;
       if ( ux >= ux_VL_ ) ux = ux_VL_;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krw = gamma( ux, aw1_, bw1_, cw1_,dw1_ ); 
	   const double64 C2_krw = gamma( ux, aw2_, bw2_, cw2_,dw2_ );
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  );
       return std::min( std::max( Krw_ave, 0. ), 1. ); 
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_ ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.; 
	   if ( ux >= ux_VL_ ) ux = ux_VL_;;
       if ( ux <= ux_CL_ ) ux = ux_CL_;
       const double64 rw = r_w( Sw, Swi_ );
	   const double64 C1_krnw = gamma( ux, anw1_, bnw1_, cnw1_,dnw1_ ); 
	   const double64 C2_krnw = 0.76;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 0.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 1.; 
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krw = gamma( ux, aw1_cross, bw1_cross, cw1_cross,dw1_cross ); 
	   const double64 C2_krw = 0.78;
       const double64 Krw_ave = k_w( rw, C1_krw, C2_krw  ); 
       return std::min( std::max( Krw_ave, 0. ), 1. );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( Sw <= Swi_cross ) return 1.;
       if ( fabs(Sw - 1.) < std::numeric_limits<double64>::epsilon() * 100.) return 0.;
	   if ( ux >= ux_cross_VL_ ) ux = ux_cross_VL_;
       if ( ux <= ux_cross_CL_ ) ux = ux_cross_CL_;
	   const double64 rw = r_w( Sw, Swi_cross );
	   const double64 C1_krnw = gamma( ux, anw1_cross, bnw1_cross, cnw1_cross,dnw1_cross ); 
	   const double64 C2_krnw = 0.78;
       const double64 Krn_ave = k_nw( rw, C1_krnw, C2_krnw  );
       return std::min( std::max( Krn_ave, 0. ), 1.);
    }

  const std::string name = "X-CSst-FSst";  
  const int      rocktype_ = 14;
  const int      subtypes_ =  2;
  const double64 k_low_=3.61e-13, k_high_=2.50e-12, 
                 phi_  = 0.283,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.159, Swi_high_=0.104, 
                 Swi_pc_=0.14,
                 Sgr_   = 0.33,    
                 m_low_=0.6, m_high_=0.7,    
                 m_=0.646422242,
                 pd_low_=1000., pd_high_=750., 
                 pd_=8.545196768087378e+02,    // This value is the one used in VG model, is this correct?
                 dPc_=250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 5.00e-01, ux_CL_ = 5.00e-08,   // Threshold velocity for VL and CL of parallel flow 
                 ux_cross_VL_ = 5.00e-02, ux_cross_CL_ = 5.00e-08,   // Threshold velocity for VL and CL of cross flow 
		 aw1_ = 0.85, bw1_ = 897.8472917, cw1_ = 1., dw1_ = 2.55,   // used for gamma_W_1 (parallel flow)
		 aw2_ = 0.13, bw2_ = 601.8450379, cw2_ = 1., dw2_ = 0.605,  // used for gamma_W_2 (parallel flow)
		 anw1_ = 0.43, bnw1_ = 0.004086771, cnw1_ = -1., dnw1_ = 1.5, // used for gamma_NW_1 (parallel flow)	 		 
		 aw1_cross = 1.2, bw1_cross = 6.9598e-08, cw1_cross = -1.5, dw1_cross = 2.8,  // used for gamma_W_1 (cross flow)
		 anw1_cross = 0.75, bnw1_cross = 403.4287935, cnw1_cross = 1., dnw1_cross = 1.56, // used for gamma_NW_1 (cross flow)		 
                 Swi_ = 0.22,           // for parallel flow for kr calculation
		 Swi_cross = 0.13;      // for cross flow for kr calculation
};



/// Tuple container to return rocktype by number (std::get<1>(rocktypes_);
struct OtwayRockTypes {
    bool IsComposite( int RT ) const {
         // TODO: update if changes are made
         if ( RT <= 3 || RT == 10 || RT == 13 || RT == 15 ) return false;
         return true;
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


#endif /* CSMP_OTWAY_CRC3_ROCK_TYPES_VERSION_2_H */
