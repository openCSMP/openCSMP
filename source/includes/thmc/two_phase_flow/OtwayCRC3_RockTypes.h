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
  
  

// UTILITY FUNCTIONS
// functions used repeatedly in the curve fitting for layer-parallel and perpendicular flows

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

inline double64 polyC2( double64 Sw, double64 C1, double64 C2 ) {
     return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
  }

inline double64 polyC3( double64 Sw, double64 C1, double64 C2, double64 C3 ) {

     return C1 * std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3);
  }

inline double64 poly_abc( double64 a, double64 b, double64 c ) {
     return std::pow( a, b ) + c;
  }
  
  

// ROCK-TYPE SPECIFICATIONS

/// WELL - rocktype representing a sand-face completion of the well
struct CRC3_RockType0 {
  /// this function establishes whether kri is rate- and flow-direction dependent
  bool IsComposite() const { return false; }

  /// Linear relative permeability: water
  double64 Krw( double64 Sw ) const {
      Sw = std::max( std::max(Sw,1.), 0. );
      return Sw;
   }
  
  /// Linear relative permeability: CO2 
  double64 Krn( double64 Sw ) const {
       Sw = std::max( std::min(Sw,1.), 0. );
       return 1. - Sw;
    }
    
  /// capillary pressure that does not depend on Sw but on radius of well completion  
  double64 Pc( double64 /* Sw */ ) const {
       const double64 IFT = 0.05; // interfacial tension (N/m) 
       return IFT / diameter_;
    }

  const std::string name = "well";        
  const int      rocktype_ = 0;
  const int      subtypes_ = 1;
  const double64 k_        = 1.0e-7,  // permeability
                 diameter_ = 0.2;     // of the sandface section
};



/** 
     homogeneous mudstone - no special provisions
 */
struct CRC3_RockType1 {
  bool IsComposite() const { return false; }
  
  /// Water relative permeability
  double64 Krw( double64 Sw ) const {
      const double64 C1(11.73), C2(0.3316);
      return polyC2( Sw, C1, C2 );
   }
  
  /// CO2 relative permeability:
  double64 Krnw( double64 Sw ) const {
       const double64 C1(2.848), C2(2.042), C3(3.892);
       return polyC3( Sw, C1, C2, C3 );
    }
    
  const std::string name = "H-Mst";        
  const int      rocktype_ = 1;
  const int      subtypes_ = 1;
  const double64 k_    = 1.73e-15,  // permeability
                 LY_   = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                 Swi_  = 0.39,    // irreducible saturation
                 m_    = 0.4,       // van Genuchten exponent
                 pd_   = 5000.;     // capillary (drainage) entry pressure of low and high
};




/** 
     Massive bedding: Carbonate-cemented Sandstone - Mudstone
 */
struct CRC3_RockType2 {
  bool IsComposite() const { return false; }
  
  /// Water relative permeability (CL)
  double64 Krw( double64 Sw ) const {
       return krw_VG( seffL(Sw,Swi_), m_ave_ );
    }
  
  /// CO2 relative permeability:
  double64 Krn( double64 Sw ) const {
       return krn_BC( seffL(Sw,Swi_) );
    }

  const std::string name = "M-CbSst-Mst";        
  const int      rocktype_ = 1;
  const int      subtypes_ = 1;
  const double64 k_    = (1.539564e-15 + 1.727075e-15) / 2., // permeability
                 LY_low_  = 0.025,  LY_high_ = 0.025,        // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.552,  Swi_high_ = 0.39,        // irreducible saturations of the 2 different layers
                 m_low_   = 0.4,    m_high_ = 0.4,           // van Genuchten exponents for the 2 different layers
                 pd_low_  = 10000., pd_high_ = 5000.,        // capillary (drainage) entry pressure of low and high
                 dPc_     = 5000., // UNSPECIFIED: guess of pc-difference between high_k and low_k layer at swc
                 Swi_     = 0.472,
                 m_ave_   = 0.427;
};




///  homogeneous carbonate-cemented sandstone
struct CRC3_RockType3 {
  bool IsComposite() const { return false; }
  /// Water relative permeability
  double64 Krw( double64 Sw ) const {
      const double64 C1(16.6), C2(0.3374);
      return polyC2( Sw, C1, C2 );
   }
  
  /// CO2 relative permeability:
  double64 Krnw( double64 Sw ) const {
       const double64 C1(5.458), C2(2.051), C3(5.558);
       return polyC3( Sw, C1, C2, C3 );
    }

  const std::string name = "H-CbSst";  
  const int      rocktype_ = 3;
  const int      subtypes_ = 1;
  const double64 k_        = 1.54e-15, // permeability
                 LY_       = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                 Swi_      = 0.552,    // irreducible saturation
                 m_        = 0.4,      // van Genuchten exponent
                 pd_       = 10000.;   // capillary (drainage) entry pressure of low and high
};






/// Homogeneous (Carbonate-cemented sandstone):
struct CRC3_RockType10 {
  bool IsComposite() const { return false; }
 
  /// Water relative permeability
  double64 Krw( double64 Sw ) const {
     const double64 C1(6.337), C2(0.4387);
     return polyC2( Sw, C1, C2 );
  }
 
  /// CO2 relative permeability:
  double64 Krnw( double64 Sw ) const {
      const double64 C1(1.522), C2(2.027), C3(2.695);
      return polyC3( Sw, C1, C2, C3 );
   }

  const std::string name = "H-Slt";  
  const int      rocktype_ = 10;
  const int      subtypes_ = 1;
  const double64 k_        = 3.48e-14, // permeability
                 LY_1_     = 0.05,     //  permeabilities individual layers
                 LY_       = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                 Swi_      = 0.18,     // irreducible saturation
                 m_        = 0.5,      // van Genuchten exponent
                 pd_       = 3000.;    // capillary (drainage) entry pressure of low and high

};




struct CRC3_RockType13 {
   bool IsComposite() const { return false; }

   /// Water relative permeability
   double64 Krw( double64 Sw ) const {
       const double64 C1(4.936), C2(0.5563);
       return polyC2( Sw, C1, C2 );
    }
   
   /// CO2 relative permeability:
   double64 Krnw( double64 Sw ) const {
        const double64 C1(1.442), C2(2.022), C3(2.594);
        return polyC3( Sw, C1, C2, C3 );
     }

  const std::string name = "H-FSlt";
  const int      rocktype_ = 13;
  const int      subtypes_ = 1;
  const double64 k_        = 3.61e-13, // permeability
                 LY_       = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                 Swi_      = 0.159,    // irreducible saturation
                 m_        = 0.6,      // van Genuchten exponent
                 pd_       = 1000.;    // capillary (drainage) entry pressure of low and high
};




/// coarse sandstone
struct CRC3_RockType15 {
  bool IsComposite() const { return false; }

  /// Water relative permeability
  double64 Krw( double64 Sw ) const {
     const double64 C1(3.755), C2(0.6705);
     return polyC2( Sw, C1, C2 );
  }
 
  /// CO2 relative permeability
  double64 Krnw( double64 Sw ) const {
      const double64 C1(1.26), C2(2.012), C3(2.362);
      return polyC3( Sw, C1, C2, C3 );
   }

  const std::string name = "H-CSst"; ///< carbonate-cemented sandstones
  const int      rocktype_ = 15;
  const int      subtypes_ = 1;
  const double64 k_        = 2.50e-12,  // permeability
                 LY_       = 0.05,      // cumulative layer thickness in the vertical direction (Y)
                 Swi_      = 0.104,     // irreducible saturation
                 m_        = 0.7,       // van Genuchten exponent
                 pd_       = 750.;      // capillary (drainage) entry pressure of low and high

};



/**
    Essentially impermeable rock with an entry pressure that is so high that it will not be exceeded
    by a plume with 100-m column height.
*/
struct CRC3_RockType16 {
   bool IsComposite() const { return false; }
   const std::string name = "baffle"; 
   const int      rocktype_ = 16;
   const int      subtypes_ = 1;
   const double64 k_        = 1.0e-16,  // permeability
                  phi_      = 0.09,     // porosity
                  LY_1_     = 0.05,      //  permeabilities individual layers
                  LY_       = 0.05,      // cumulative layer thickness in the vertical direction (Y)
                  Swi_      = 0.4,     // irreducible saturation
                  m_        = 0.7,       // van Genuchten exponent
                  pd_       = 1.0e-6;      // capillary (drainage) entry pressure of low and high

   /// Linear water relative permeability for drainage
   double64 Krw( double64 Sw ) const { return (Sw - Swi_) / (1. - Swi_); }
   
   /// Linear CO2 relative permeability for drainage
   double64 Krnw( double64 Sw ) const { return ((1.-Sw) - Swi_) / (1. - Swi_); }
};



// ==============================================

// COMPOSITE ROCK TYPES

// ==============================================

///< Planar bedding mudstone - silt
struct CRC3_RockType4 {
 bool IsComposite() const { return true; }

 /// Water relative permeability
 double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
     if ( ux >= 6e-3 ) return 0.61 * pow( seff(Sw,Swi_), 4.06 );
     if ( ux <= 1e-5 ) return 0.626 * pow( seff(Sw,Swi_), 5.11 );
     const double64 C1_krw = 10.;
   
     const double64 a2 =  0.1286,
                    b2 = -0.2088,
                    c2 =  3.691,
                    C2_krw = a2 * std::pow(ux,b2) + c2;
     
     // Krw_ave (for averaged velocity?)
     return C1_krw * std::pow( seffL(Sw,Swi_),C2_krw);
  }
 
 /// CO2 relative permeability:
 double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
    if ( ux >= 1e-5 ) return 0.323 * krn_BC( seff(Sw,Swi_) );
    if ( ux <= 5e-8 ) return 0.488 * krn_BC( seff(Sw,Swi_) );
    const double64 a_c1 =  0.04757,
                   b_c1 =  -0.1498,
                   c_c1 =   0.2208,
                   C1_krnw = a_c1 * std::pow(ux,b_c1) + c_c1;
      // Krnw_ave=
      return C1_krnw * krn_BC( seffL(Sw,Swi_) );
   }
   
   /// Water relative permeability
   double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= 1.0e-4 ) return 0.581 * std::pow( seffL(Sw,Swi_), 5.26 );
       if ( ux <= 1.0e-7 ) return 0.559 * std::pow( seffL(Sw,Swi_), 3.15 );
       const double64 a1 =  0.2485,
                      b1 =  0.01877,
                      c1 =  0.3717,
                      C1_krw = a1 * std::pow(ux,b1) + c1;

       const double64 a2 =  0.2485,
                      b2 =  0.01877,
                      c2 =  0.3717,
                      C2_krw = a2 * std::pow(ux,b2) + c2;
       // Krw_ave=
       return C1_krw * std::pow( seffL(Sw,Swi_),C2_krw);
    }
   
   /// CO2 relative permeability:
   double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
        if ( ux >= 1.0e-4 ) return 0.267 * krn_BC( seff(Sw,Swi_,0.0153) );
        if ( ux <= 1.0e-7 ) return 0.0794 * krn_BC( seff(Sw,Swi_,0.0206) );
        const double64 a =  0.0075237,
                       b = -0.05069,
                       c =  0.002388,
                       Sgr_ux = a * std::pow(ux,b) + c;

        const double64 a1 =  3.612,
                       b1 =  0.3236,
                       c1 =  0.06067,
                       C1_krn = a1 * std::pow(ux,b1) + c1;
        // Krnw_ave=
        return C1_krn * krn_BC( seffL(Sw,Swi_,Sgr_ux) );
     }
   
  // ROCK PROPERTIES 
  const std::string name = "P-Mst-Slt";  
  const int      rocktype_ = 4;
  const int      subtypes_ = 2;
  const double64 k_low_   = 1.73e-15, k_high_ = 3.48e-14,      // layer permeabilities
                 LY_1     = 0.015,  LY_2 = 0.01, LY_3 = 0.005, // layer thicknesses
                 LY_4     = 0.015, LY_5 = 0.005,      // (layer 1,3,5 are high_k layers, layer 2, 4 are low_k layers)
                 LY_low_  = 0.25,  LY_high_ = 0.25,   // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.39,  Swi_high_ = 0.18,  // irreducible saturations of the 2 different layers
                 Swi_     = 0.562,
                 m_low_   = 0.4,   m_high_ = 0.5,     // van Genuchten exponents for the 2 different layers
                 pd_low_  = 5000., pd_high_ = 3000.,  // capillary (drainage) entry pressure of low and high
                 dPc      = 2000.; // UNSPECIFIED: guess of pc difference between high_k and low_k layer at swc

}; // end CRC3_RockType4







/**  
    homogeneous carbonate-cemented sandstone (originally perceived as a composite, but treated as single rocktype)
    
    @attention although this is a composite with a layered K structure and a directionally dependent permeability, rate dependence is ignored
*/
struct CRC3_RockType5 {
  bool IsComposite() const { return true; } // TODO: will this work like this?
  /// Water relative permeability
  double64 Krw_ParallelDrainage( double64 Sw ) const {
       const double64 m_ave(0.493), Swi(0.5);
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       const double64 Swi(0.5);       
       return krn_BC( seffL(Sw,Swi) );
    }

  /// Water relative permeability (capillary limit)
  double64 Krw_CrossDrainage( double64 Sw ) const {
       const double64 m_ave(0.636), Swi(0.359);
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability (capillary limit)
  double64 Krn_CrossDrainage( double64 Sw ) const {
       const double64 Swi(0.359);       
       return krn_BC( seffL(Sw,Swi) );
    }

  const std::string name = "H-CbSst-Slt";  
  const int      rocktype_ = 5;
  const int      subtypes_ = 1;
  const double64 k_        = 1.54e-15, 
                 LY_1_     = 0.025,    
                 LY_2_     = 0.025,    
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.18,     
                 m_low_    = 0.4,      
                 m_high_   = 0.5,      
                 pd_low_   = 10000.,    
                 pd_high_  = 3000.;    
};






/// composite - planar bedding mudstone - fine sandstone
struct CRC3_RockType6 {
  bool IsComposite() const { return true; }
  
  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       // viscous limit approximation
       if ( ux >= 1.0e-2 ) return 0.728 * pow( seff(Sw,Swi_), 3.47 );
       // capillary limit approximation
       if ( ux <= 1.0e-5 ) return 0.75  * pow( seff(Sw,Swi_), 5.08 );
      // between CL and VL
       const double64 a1(0.0002717), b1(-0.3894), c1(0.7265);
       const double64 C1_krw = a1 * std::pow( ux, b1 ) + c1; 
       const double64 a2(2.2), b2(-0.06136), c2(0.6185);
       const double64 C2_krw = a2 * std::pow( ux, b2) + c2;
       // Krw_ave 
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }
    
   double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
        if ( ux >= 1.0e-2 ) return 0.4 * krn_BC( seffL(Sw,Swi_));
        if ( ux <= 1.0e-5 ) return 0.87 * krn_BC( seffL(Sw,Swi_));
        const double64 Krnw_ave = 3.904 * pow( ux, -0.01475 ) - 3.756;
        // limiting value to >=0
        return std::max( Krnw_ave, 0. );
     }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       // viscous limit approximation
       if ( ux >= 7.0e-4 ) return 0.557 * pow( seff(Sw,Swi_), 5.74 );
       // capillary limit approximation
       if ( ux <= 5.0e-7 ) return 0.834 * pow( seff(Sw,Swi_), 1.57 );
      // between CL and VL
       const double64 a1(0.2151), b1(-0.07629), c1(0.1834);
       const double64 C1_krw = a1 * std::pow( ux, b1 ) + c1; 
       const double64 a2(-7.324), b2(-0.04686), c2(16.2);
       const double64 C2_krw = a2 * std::pow( ux, b2) + c2;
       // Krw_ave 
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }
    
   double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
        if ( ux >= 7.0e-4 ) return 0.2014 * krn_BC( seffL(Sw,Swi_,0.00397));
        if ( ux <= 5.0e-7 ) return 0.018  * krn_BC( seffL(Sw,Swi_,0.1154));
        const double64 Sgr_ux = poly_abc( 0.08238, -0.0783, -0.1412 ),
                       C1_krn = poly_abc( 0.7342, 0.1051, -0.1418 ),
                       Krnw_ave = C1_krn * krn_BC( seffL(Sw,Swi_,Sgr_ux) );
        // limiting value to >=0
        return std::max( Krnw_ave, 0. );
     }

  const std::string name = "P-Mst-FSst";  
  const int      rocktype_ = 6;
  const int      subtypes_ = 2;
  const double64 k_low_=1.73e-15, k_high_=3.61e-13, 
                 LY_1_=0.02, LY_2_=0.01, LY_3_=0.005, LY_4_=0.015,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.159, 
                 m_low_=0.4, m_high_=0.6,        // van Genuchten exponents for the 2 different layers
                 pd_low_=5000., pd_high_=1000.,  // capillary (drainage) entry pressure of low and high
                 dPc=4000., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 Swi_ = 0.496;
};







///  massive carbonate cemented Sandstone – Fine Sandstone (composite, but not rate dependent)
struct CRC3_RockType7 {
  bool IsComposite() const { return true; }
  /// Water relative permeability
  double64 Krw_ParallelDrainage( double64 Sw ) const {
       const double64 m_ave(0.6), Swi(0.47);
       return krw_VG( seffL(Sw,Swi), m_ave );
    }
  
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       const double64 Swi(0.47);
       return krn_BC( seffL(Sw,Swi) );
    }

  /// Water relative permeability
  double64 Krw_CrossDrainage( double64 Sw ) const {
      const double64 m_ave(0.593), Swi(0.46);
      return krw_VG( seffL(Sw,Swi), m_ave );
   }
  
  /// CO2 relative permeability:
  double64 Krn_CrossDrainage( double64 Sw ) const {
       const double64 Swi(0.46);
       return krn_BC( seffL(Sw,Swi) );
    }

  const std::string name = "M-CbSst-FSst";  
  const int      rocktype_ = 7;
  const int      subtypes_ = 1;
  const double64 k_low_    = 1.54e-15, 
                 k_high_   = 3.61e-13,    
                 LY_1_     = 0.025,    
                 LY_2_     = 0.025,    
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.18,     
                 m_low_    = 0.4,      
                 m_high_   = 0.6,      
                 pd_low_   = 10000.,    
                 pd_high_  = 3000.;    
};







/// composite - Planar bedding Mudstone – Coarse Sandstone - @todo: appears to be vertically layered?
struct CRC3_RockType8 {
  bool IsComposite() const { return true; }
  
  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( ux <= 5.0e-5 ) return 0.75 * std::pow( seffL(Sw,Swi_), 5.08 );
       if ( ux >= 4.0e-2 ) return 0.728 * std::pow( seffL(Sw,Swi_), 3.47 );
       // between CL and VL
       const double64 C1_krw = poly_abc( 0.0001377, -0.5243, 0.8232 );
       const double64 C2_krw = poly_abc( 0.775, -0.114, 2.078 );
       return C1_krw * pow( seff(Sw,Swi_), C2_krw );
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( ux <= 5.0e-5 ) return 0.87 * krn_BC( seff(Sw, Swi_) );
       if ( ux >= 4.0e-2 ) return 0.4 * krn_BC( seff(Sw, Swi_) );
       // between CL and VL
       const double64 C1_krn = poly_abc( 0.755, -0.06682, -0.3927 );
       // limiting range between 0. and 0.43
       return std::max( std::min( C1_krn * krn_BC( seff(Sw,Swi_) ), 0.43 ), 0. );
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux <= 2e-4 ) return 0.492 * std::pow( seffL(Sw,Swi_), 4.19 );
       if ( ux >= 5e-7 ) return std::min( 2.98  * std::pow( seffL(Sw,Swi_), 1.66 ), 0.4921 );
       // between CL and VL
       const double64 C1_krw = poly_abc( 2.446, -0.0728, -4.056 );
       const double64 C2_krw = poly_abc( -133.2, -0.003056, 140.9 );
       return C1_krw * pow( seff(Sw,Swi_), C2_krw );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux <= 2e-4 ) return std::max( 0.1361 * krn_BC( seff(Sw, Swi_,0.0658) ), 0. );
       if ( ux >= 5e-7 ) return std::max( 0.0134 * krn_BC( seff(Sw, Swi_,0.427) ), 0. );
       // between CL and VL
       const double64 Sgr_ux = poly_abc( 1.687, -0.02636, -2.046 ),
                      C1_krn = poly_abc( 0.9455, 0.1964, -0.04131 );
       // limiting range >= 0. 
       return std::max( C1_krn * krn_BC( seff(Sw,Swi_,Sgr_ux) ), 0. );
    }

  const std::string name = "P-Mst-CSst";  
  const int      rocktype_ = 8;
  const int      subtypes_ = 2;
  const double64 k_low_=1.73e-15, k_high_=2.50e-12, 
                 LY_1_=0.02, LY_2_=0.01, LY_3_=0.005, LY_4_=0.015,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.104, 
                 m_low_=0.4, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 pd_low_=5000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 dPc=4250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 Swi_=0.426;
};








///  Massive bedding: Carbonate cemented Sandstone – Coarse Sandstone (not rate dependent, but composite)
struct CRC3_RockType9 {
  bool IsComposite() const { return true; }
  
  /// Water relative permeability
  double64 Krw_ParallelDrainage( double64 Sw ) const {
      const double64 m_ave(0.7), Swi(0.44);
      return krw_VG( seffL(Sw,Swi), m_ave );
   }
  
  /// CO2 relative permeability:
  double64 Krn_ParallelDrainage( double64 Sw ) const {
       const double64 Swi(0.44);
       return krn_BC( seffL(Sw,Swi) );
    }

  /// Water relative permeability
  double64 Krw_CrossDrainage( double64 Sw ) const {
      const double64 m_ave(0.469), Swi(0.45);
      return krw_VG( seffL(Sw,Swi), m_ave );
   }
  
  /// CO2 relative permeability:
  double64 Krn_CrossDrainage( double64 Sw ) const {
       const double64 Swi(0.275);
       return krn_BC( seffL(Sw,Swi) );
    }

  const std::string name = "M-CbSst-CSst";  
  const int      rocktype_ = 9;
  const int      subtypes_ = 1;
  const double64 k_low_    = 1.54e-15, 
                 k_high_   = 2.5e-12,    
                 LY_1_     = 0.025,    
                 LY_2_     = 0.025,    
                 LY_low_   = 0.025,    
                 LY_high_  = 0.025,    
                 Swi_low_  = 0.552,    
                 Swi_high_ = 0.104,     
                 m_low_    = 0.4,      
                 m_high_   = 0.7,      
                 pd_low_   = 10000.,    
                 pd_high_  = 750.;    
};









/// composite - planar bedding fine sandstone - silt
struct CRC3_RockType11 {
 bool IsComposite() const { return true; }

 /// Water relative permeability
 double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
     if ( ux >= 3.0e-2 ) return 0.731 * std::pow( seffL(Sw,Swi_), 4.17 ); // VL;
     if ( ux <= 1.0e-5 ) return 0.739 * std::pow( seffL(Sw,Swi_), 5.95 ); // CL;     
     const double64 C1_krw = poly_abc( 0.00273, -0.1325, 0.7263 );
     const double64 C2_krw = poly_abc( 49.75, -0.004283, -46.31 );
     // Krw_ave
     return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
  }
 
 /// CO2 relative permeability:
 double64 Krn_ParallelDrainage( double64 Sw, double64 ux  ) const {
     if ( ux >= 3.0e-2 ) return 0.545 * krn_BC( seffL(Sw,Swi_) );
     if ( ux <= 1.0e-5 ) return std::min( 1.33 * krn_BC( seffL(Sw,Swi_) ), 0.49 );    
     const double64 C1_krn = poly_abc( -30.39, 0.003259, 30.6 );
     // Krnw_ave=
     const double64 Krnw_ave = C1_krn * krn_BC( seffL(Sw,Swi_) );
     return std::min( std::max( Krnw_ave, 0. ), 0.43 );
  }

  /// Water relative permeability
  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
      if ( ux >= 3.0e-2 ) return 0.6197 * std::pow( seffL(Sw,Swi_), 4.74 ); // VL;
      if ( ux <= 1.0e-5 ) return std::min( 0.993 * std::pow( seffL(Sw,Swi_), 3.26 ), 0.6197 ); // CL;     
      const double64 C1_krw = poly_abc( 0.0066, -0.2977, 0.5443 );
      const double64 C2_krw = poly_abc( -0.7936, -0.1001, 6.542 );
      // Krw_ave
      return std::min( C1_krw * std::pow( seffL(Sw,Swi_), C2_krw ), 0.6197 );
   }
  
  /// CO2 relative permeability:
  double64 Krn_CrossDrainage( double64 Sw, double64 ux  ) const {
      if ( ux >= 3.0e-2 ) return std::max( 1.33 * krn_BC( seffL(Sw,Swi_,0.037) ), 0. ); 
      if ( ux <= 1.0e-5 ) return std::max( 0.0416 * krn_BC( seffL(Sw,Swi_,0.159) ), 0. );    
      const double64 Sgr_ux = poly_abc( 0.003424, -0.2677, 0.006397 ),
                     C1_krn = poly_abc( 3.153, 0.1897, -0.1728 );
      // Krnw_ave=
      const double64 Krnw_ave = C1_krn * krn_BC( seffL(Sw,Swi_,Sgr_ux) );
      return std::max( Krnw_ave, 0. );
   }



  const std::string name = "P-FSst-Slt";  
  // Composite1 - hypothetical sample based on Achyut's rocktypes
  const int      rocktype_ = 11;
  const int      subtypes_ = 2;
  const double64 k_low_   = 3.48e-14,
                 k_high_  = 3.61e-13,  // layer permeabilities
                 LY_1     = 0.015, LY_2 = 0.01, LY_3 = 0.005,
                 LY_4     = 0.015, LY_5 = 0.005,
                 LY_low_  = 0.25,  LY_high_ = 0.25,    // cumulative layer thickness in the vertical direction (Y)
                 Swi_low_ = 0.18,  Swi_high_ = 0.159, // irreducible saturations of the 2 different layers
                 m_low_   = 0.5,   m_high_ = 0.6,       // van Genuchten exponents for the 2 different layers
                 pd_low_  = 3000., pd_high_ = 1000.,  // capillary (drainage) entry pressure of low and high
                 dPc      = 1.3061e+06, // pc difference between high_k and low_k layer at swc
                 Swi_     = 0.335; // for both cross-layer and layer parallel flow
};







/// composite - Planar bedding Siltstone – Coarse Sandstone 
struct CRC3_RockType12 {
  bool IsComposite() const { return true; }

  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_VL_ ) return 0.83  * std::pow( seffL(Sw,Swi_), 3.44 );
       if ( ux <= ux_CL_ ) return 0.869 * std::pow( seffL(Sw,Swi_), 6.14 );
       const double64 C1_krw = poly_abc( 0.002649, -0.2458, 0.8242 );
       const double64 C2_krw = poly_abc( -29.82, 0.01156, 32.25 );
       // Krw_ave =
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_VL_ ) return std::min( 0.637 * krn_BC( seffL(Sw,Swi_) ), 0.43 );
       if ( ux <= ux_CL_ ) return std::min( 2.094 * krn_BC( seffL(Sw,Swi_) ), 0.43 );
       const double64 C1_krn = poly_abc( -6.277, 0.03536, 6.239 );
       const double64 Krn_ave = C1_krn * krn_BC( seffL(Sw,Swi_) );
       return std::min( Krn_ave, 0.43 );
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_cross_VL_ ) return 0.544966  * std::pow( seffL(Sw,Swi_), 4.437 );
       if ( ux <= ux_cross_CL_ ) return std::min( 2.55 * std::pow( seffL(Sw,Swi_), 3.066 ), 0.545 );
       const double64 C1_krw = poly_abc( 2.036, -0.06943, -3.025 );
       const double64 C2_krw = poly_abc( 6.683, 0.1167, 1.837 );
       // Krw_ave =
       return std::min( C1_krw * std::pow( seffL(Sw,Swi_), C2_krw ), 0.545 );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_cross_VL_ ) return std::max( 0.3541 * krn_BC( seffL(Sw,Swi_,0.0463) ), 0. );
       if ( ux <= ux_cross_CL_ ) return std::max( 0.0239 * krn_BC( seffL(Sw,Swi_,0.327) ), 0. );
       const double64 Sgr_ux = poly_abc( -1.888, 0.03386, 1.482 ),
                      C1_krn = poly_abc( 3.627, 0.2739, -0.04431 );
       const double64 Krn_ave = C1_krn * krn_BC( seffL(Sw,Swi_,Sgr_ux) );
       return std::max( Krn_ave, 0. );
    }


  const std::string name = "P-Slt-CSst";  
  const int      rocktype_ = 12;
  const int      subtypes_ =  2;
  const double64 k_low_=3.48e-14, k_high_=2.50e-12, 
                 LY_1_=0.02, LY_2_=0.01, LY_3_=0.005, LY_4_=0.015,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.18, Swi_high_=0.104, 
                 m_low_=0.5, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 pd_low_=3000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 dPc=4250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 5.0e-2, ux_CL_ = 1e-5,
                 ux_cross_VL_ = 3.0e-4, ux_cross_CL_ = 5e-7,
                 Swi_ = 0.269;
};






/// composite - cross-bedded coarse - fine laminated sandstone
struct CRC3_RockType14 {
  bool IsComposite() const { return true; }

  double64 Krw_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_VL_ ) return 0.82 * std::pow( seffL(Sw,Swi_), 3.69 );
       if ( ux <= ux_CL_ ) return 0.86 * std::pow( seffL(Sw,Swi_), 4.64 );
       const double64 C1_krw = poly_abc( 0.008275, -0.1924, 0.811 );
       const double64 C2_krw = poly_abc( 0.4083, -0.1338, 3.239 );
       // Krw_ave =
       return C1_krw * std::pow( seffL(Sw,Swi_), C2_krw );
    }

  double64 Krn_ParallelDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_VL_ ) return std::min( 0.68 * krn_BC( seffL(Sw,Swi_) ), 0.45 );
       if ( ux <= ux_CL_ ) return std::min( 1.07 * krn_BC( seffL(Sw,Swi_) ), 0.45 );
       const double64 C1_krn = poly_abc( 0.2568, -0.1039, 0.4035 );
       const double64 Krn_ave = C1_krn * krn_BC( seffL(Sw,Swi_) );
       return std::min( Krn_ave, 0.45 );
    }

  double64 Krw_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_VL_ ) return 0.756 * std::pow( seffL(Sw,Swi_), 4.22 );
       if ( ux <= ux_CL_ ) return std::min( 0.715 * std::pow( seffL(Sw,Swi_), 2.944 ), 0.7557 );
       const double64 C1_krw = poly_abc( 0.8557, 0.5453, 0.7148 );
       const double64 C2_krw = poly_abc( 4.545, 0.1916, 2.662 );
       // Krw_ave =
       return std::min( C1_krw * std::pow( seffL(Sw,Swi_), C2_krw ), 0.7557 );
    }

  double64 Krn_CrossDrainage( double64 Sw, double64 ux ) const {
       if ( ux >= ux_VL_ ) return std::max( 0.793 * krn_BC( seffL(Sw,Swi_,0.436) ), 0. );
       if ( ux <= ux_CL_ ) return std::max( 0.316 * krn_BC( seffL(Sw,Swi_,0.0516) ), 0. );
       const double64 Sgr_ux = poly_abc( 0.0175, -0.03539, 0.02231 ),
                      C1_krn = poly_abc( 3.329, 0.3396, 0.2929 );
       const double64 Krn_ave = C1_krn * krn_BC( seffL(Sw,Swi_,Sgr_ux) );
       return std::max( Krn_ave, 0. );
    }

  const std::string name = "X-CSst-FSst";  
  const int      rocktype_ = 14;
  const int      subtypes_ =  2;
  const double64 k_low_=3.61e-13, k_high_=2.50e-12, 
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.159, Swi_high_=0.104, 
                 m_low_=0.6, m_high_=0.7,        
                 pd_low_=1000., pd_high_=750.,  
                 dPc=250., // NOT SPECIFIED: pc difference between high_k and low_k layer at swc
                 ux_VL_ = 5.0e-1, ux_CL_ = 1e-4,
                 ux_cross_VL_ = 4.0e-3, ux_cross_CL_ = 5e-7,
                 Swi_ = 0.262;
};



/// Container that returns rocktype by number
struct OtwayRockTypes {
//    bool     IsComposite( ROCKTYPE rocktype ) const { return std::get<rocktype>(rocktypes_).IsComposite(); }
    double64 Krw_Parallel( double64 sw, double64 v ) const;
  
    std::tuple<CRC3_RockType0, // well
               CRC3_RockType1,CRC3_RockType2,CRC3_RockType3,
               CRC3_RockType4,CRC3_RockType5,CRC3_RockType6,
               CRC3_RockType7,CRC3_RockType8,CRC3_RockType9,
               CRC3_RockType10,CRC3_RockType11,CRC3_RockType12,
               CRC3_RockType13,CRC3_RockType14,CRC3_RockType15,
               CRC3_RockType16>  rocktypes_;  
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
