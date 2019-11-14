//
//  OtwayCRC3_RockTypes.h
//  CSMP_CO2GeoSequestrationSimulator
//
//  Created by Stephan Matthai on 23/9/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_OTWAY_CRC3_ROCK_TYPES_H
#define CSMP_OTWAY_CRC3_ROCK_TYPES_H

#include "CSMP_definitions.h"

// HARDCODED ROCK PROPERTIES FOR SIMULATION OF CRC3-CRC2 CROSS SECTION
// based on data compiled by Maartje Boon (23/9/2019)
// SKM 24/9/19 - replaced layer thicknesses with volume fractions

namespace csmp {

/**
  Creating a mapping which associates the rocktype numbers in Otway cross-section CRC3-CR2 with
  the names of the rock so that they can be seen during debugging.

  many of the composites can be treated as mixtures between impermeable nonporous carbonate cemented rocks
  and single rock types.
  In these cases we can run with the conventional relperms for the remaining porous subtypes, a downscaled porosity, and an anisotropic permeability.
  These are:

    M_CbSst_Mst      Massive bedding (Carbonate-cemented sandstone lamina and mudstone matrix) - ONLY ONE SUBTYPE FLOWS
    M_CbSst_Slt        Massive bedding (Carbonate-cemented sandstone lamina and siltstone matrix) - ONLY ONE SUBTYPE FLOWS
    M_CbSst_FSs      Massive bedding (Carbonate-cemented sandstone lamina and fine sandstone matrix) - ONLY ONE SUBTYPE FLOWS
    M_CbSst_CSst     Massive bedding (Carbonate-cemented sandstone lamina and coarse sandstone matrix) - ONLY ONE SUBTYPE FLOWS

 No flow will occur through:

    H_CbSst,      // -- Homogeneous (Carbonate-cemented sandstone) -- NO_FLOW

 which leaves us with the 6 composites:

    P_Mst_Slt,    // 4 COMPOSITE -- Planar bedding (Mudstone lamina and siltstone matrix)
    P_Mst_FSst,   // 6 COMPOSITE: Planar bedding (Mudstone lamina and fine sandstone matrix)
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
    H_CSst        // 15 Homogeneous (Coarse sandstone)
  };



/// homogeneous mudstone
struct CRC3_RockType1 {
  const std::string name = "H-Mst";        
  const int      rocktype_ = 1;
  const int      subtypes_ = 1;
  const double64 k_    = 1.73e-15,  // permeability
                 LY_1_ = 0.05,   //  permeabilities individual layers
                 LY_   = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                 Swi_  = 0.39,    // irreducible saturation
                 m_    = 0.4,       // van Genuchten exponent
                 pd_   = 5000.;     // capillary (drainage) entry pressure of low and high
  
  /// Water relative permeability
  double64 Krw( double64 Sw ) const {
      const double64 C1(11.73), C2(0.3316);
      return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
   }
  
  /// CO2 relative permeability:
  double64 Krnw( double64 Sw ) const {
       const double64 C1(2.848), C2(2.042), C3(3.892);
       return C1 * std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3);
    }
};




///  homogeneous carbonate-cemented sandstone
struct CRC3_RockType3 {
  const std::string name = "H-CbSst";  
  const int      rocktype_ = 3;
  const int      subtypes_ = 1;
  const double64 k_        = 1.54e-15, // permeability
                 LY_1_     = 0.05,     //  permeabilities individual layers
                 LY_       = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                 Swi_      = 0.552,    // irreducible saturation
                 m_        = 0.4,      // van Genuchten exponent
                 pd_       = 10000.;   // capillary (drainage) entry pressure of low and high

  /// Water relative permeability
  double64 Krw( double64 Sw ) const {
      const double64 C1(16.6), C2(0.3374);
      return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
   }
  
  /// CO2 relative permeability:
  double64 Krnw( double64 Sw ) const {
       const double64 C1(5.458), C2(2.051), C3(5.558);
       return C1 * std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3);
    }
};




/// Homogeneous (Carbonate-cemented sandstone):
struct CRC3_RockType10 {
 const std::string name = "H-Slt";  
 const int      rocktype_ = 10;
 const int      subtypes_ = 1;
 const double64 k_        = 3.48e-14, // permeability
                LY_1_     = 0.05,     //  permeabilities individual layers
                LY_       = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                Swi_      = 0.18,     // irreducible saturation
                m_        = 0.5,      // van Genuchten exponent
                pd_       = 3000.;    // capillary (drainage) entry pressure of low and high

 /// Water relative permeability
 double64 Krw( double64 Sw ) const {
     const double64 C1(6.337), C2(0.4387);
     return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
  }
 
 /// CO2 relative permeability:
 double64 Krnw( double64 Sw ) const {
      const double64 C1(1.522), C2(2.027), C3(2.695);
      return C1 * std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3);
   }
};




struct CRC3_RockType13 {
 const std::string name = "H-FSlt";
 const int      rocktype_ = 13;
 const int      subtypes_ = 1;
 const double64 k_        = 3.61e-13, // permeability
                LY_1_     = 0.05,     //  permeabilities individual layers
                LY_       = 0.05,     // cumulative layer thickness in the vertical direction (Y)
                Swi_      = 0.159,    // irreducible saturation
                m_        = 0.6,      // van Genuchten exponent
                pd_       = 1000.;    // capillary (drainage) entry pressure of low and high

 /// Water relative permeability
 double64 Krw( double64 Sw ) const {
     const double64 C1(4.936), C2(0.5563);
     return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
  }
 
 /// CO2 relative permeability:
 double64 Krnw( double64 Sw ) const {
      const double64 C1(1.442), C2(2.022), C3(2.594);
      return C1 * std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3);
   }
};




struct CRC3_RockType15 {
 const std::string name = "H-Sst"; ///< carbonate-cemented sandstones
 const int      rocktype_ = 15;
 const int      subtypes_ = 1;
 const double64 k_        = 2.50e-12,  // permeability
                LY_1_     = 0.05,      // permeabilities individual layers
                LY_       = 0.05,      // cumulative layer thickness in the vertical direction (Y)
                Swi_      = 0.104,     // irreducible saturation
                m_        = 0.7,       // van Genuchten exponent
                pd_       = 750.;      // capillary (drainage) entry pressure of low and high

 /// Water relative permeability
 double64 Krw( double64 Sw ) const {
     const double64 C1(3.755), C2(0.6705);
     return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
  }
 
 /// CO2 relative permeability
 double64 Krnw( double64 Sw ) const {
      const double64 C1(1.26), C2(2.012), C3(2.362);
      return C1 * std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3);
   }
};



/**
    Essentially impermeable rock with an entry pressure that is so high that it will not be exceeded
    by a plume with 100-m column height.
*/
struct CRC3_RockType16 {
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



// COMPOSITE ROCK TYPES

///< Planar bedding mudstone - silt
struct CRC3_RockType4 {
 const std::string name = "P-Mst-Slt";  
 const int      rocktype_ = 4;
 const int      subtypes_ = 2;
 const double64 k_low_   = 1.73e-15, k_high_ = 3.48e-14,  // layer permeabilities
                LY_1     = 0.015,  LY_2 = 0.01, LY_3 = 0.005, 
                LY_4     = 0.015, LY_5 = 0.005,      // layers thicknesses (layer 1,3,5 are high_k layers, layer 2, 4 are low_k layers)
                LY_low_  = 0.25,  LY_high_ = 0.25,   // cumulative layer thickness in the vertical direction (Y)
                Swi_low_ = 0.39,  Swi_high_ = 0.18,  // irreducible saturations of the 2 different layers
                m_low_   = 0.4,   m_high_ = 0.5,     // van Genuchten exponents for the 2 different layers
                pd_low_  = 5000., pd_high_ = 3000.,  // capillary (drainage) entry pressure of low and high
                dPc      = 2000.; // UNSPECIFIED: capillary pressure difference between high_k and low_k layer at connate water saturation

 /// Water relative permeability
 double64 Krw( double64 Sw, double64 ux ) const {
 // ==============================================
     if ( ux >= 1e-5 ) return Krw_VL( Sw, ux );
     if ( ux <= 5e-8 ) return Krw_CL( Sw, ux );
     const double64 a_c1 =  -0.0005899,
                    b_c1 =     -0.3025,
                    c_c1 =      0.4732,
                    C1_krw = a_c1 * std::pow(ux,b_c1) + c_c1;
   
     const double64 a_c2 =   -0.006556,
                    b_c2 =     -0.4139,
                    c_c2 =       14.96,
                    C2_krw = a_c2 * std::pow(ux,b_c2) + c_c2;
     
     const double64 a_c3 =   3.857e-05,
                    b_c3 =     -0.3863,
                    c_c3 =    0.003372,
                    C3_krw = a_c3 * std::pow(ux,b_c3) + c_c3;
     // Krw_ave
     return C1_krw * std::pow(Sw,C2_krw) + C3_krw;
  }
 
 /// CO2 relative permeability:
 double64 Krnw( double64 Sw, double64 ux ) const {
 // ==============================================
    if ( ux >= 1e-5 ) return Krnw_VL( Sw, ux );
    if ( ux <= 5e-8 ) return Krw_CL( Sw, ux );
    const double64 a_c1 =   3.698e+04,
                   b_c1 =  -3.692e+07,
                   c_c1 =       585.6,
                   d_c1 =    3.72e+04,
                   C1_krnw = a_c1 * std::exp(b_c1 * ux) + c_c1 * exp(d_c1 * ux);
      
    const double64 a_c2 = -10.5,
                   b_c2 =  -1.278e+07,
                   c_c2 =      -11.44,
                   d_c2 =        5454,
                   C2_krnw = a_c2 * std::exp(b_c2 * ux) + c_c2 * exp(d_c2 * ux);
      // Krnw_ave=
      return C1_krnw * std::exp(C2_krnw * Sw);
   }
   
   // viscous limit aproximation
   
   // Krw_VL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
   double64 Krw_VL( double64 Sw, double64 ux ) const {
        const double64 C1(11.97), C2(0.3471);
        return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
     }     
   // Krnw_VL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
   double64 Krnw_VL( double64 Sw, double64 ux ) const {
        const double64 C1(2.665), C2(1.992), C3(3.623);
        return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
     }
   // Krw_CL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
   double64 Krw_CL( double64 Sw, double64 ux ) const {
        const double64 C1(5.776), C2(0.3224);
        return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
     }
   // Krnw_CL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
   double64 Krnw_CL( double64 Sw, double64 ux ) const {
        const double64 C1(3.632), C2(2.896), C3(1.052);
        return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
     }
   
}; // end CRC3_RockType4





/// composite - planar bedding mudstone - fine sandstone
struct CRC3_RockType6 {
  const std::string name = "P-Mst-FSst";  
  const int      rocktype_ = 6;
  const int      subtypes_ = 2;
  const double64 k_low_=1.73e-15, k_high_=3.61e-13, 
                 LY_1_=0.015, LY_2_=0.01, LY_3_=0.005, LY_4_=0.015, LY_5_=0.005,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.159, 
                 m_low_=0.4, m_high_=0.6,        // van Genuchten exponents for the 2 different layers
                 pd_low_=5000., pd_high_=1000.,  // capillary (drainage) entry pressure of low and high
                 dPc=4000., // NOT SPECIFIED: capillary pressure difference between high_k and low_k layer at connate water saturation
                 VL_ux = 1.0e-5, CL_ux = 8e-8;
  
  // Krw_ave= C1_krw * Sw(ux)^C2_krw 
  double64 Krw( double64 Sw, double64 ux ) const {
  // =============================================
      if ( ux >= 1.0e-5 ) return Krw_VL( Sw, ux );
       const double64 a1(3.203), b1(-2.415e+06), c1(0.4597), d1(-351.2);
       const double64 C1_krw = a1 * std::exp(b1*ux) + c1 * std::exp(d1*ux); 
       const double64 a2(-6.609e-05), b2(-0.6526), c2(13.28);
       const double64 C2_krw = a2 * std::pow( ux, b2) + c2;
       const double64 Krw_ave = C1_krw * std::pow( Sw, C2_krw );
       return std::min( Krw_ave, Krw_CL( Sw, ux ) );
    }
    
   double64 Krnw( double64 Sw, double64 ux ) const {
   // ==============================================
        if ( ux < 1.0e-5 ) {
             // capillary limit approximation is used
             const double64 C1(11.49), C2(6.957), C3(2.826);
             return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
          }
        // Krnw_VL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
        const double64 C1(2.736), C2(1.964), C3(3.595);
        return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
     }

   // Krw_VL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2) 
   double64 Krw_VL( double64 Sw, double64 ux ) const {
        const double64 C1(12.23), C2(0.3462);
        return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
     }
   // Krw_CL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2) - not a very good fit
   double64 Krw_CL( double64 Sw, double64 ux ) const {
        const double64 C1(4.211), C2(0.6845);
        return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
     }
};




/// composite - Planar bedding Mudstone – Coarse Sandstone - @todo: appears to be vertically layered?
struct CRC3_RockType8 {
  const std::string name = "P-Mst-CSst";  
  const int      rocktype_ = 8;
  const int      subtypes_ = 2;
  const double64 k_low_=1.73e-15, k_high_=2.50e-12, 
                 LY_1_=0.015, LY_2_=0.01, LY_3_=0.005, LY_4_=0.015, LY_5_=0.005,
                 LY_low_=0.25, LY_high_=0.25,
                 Swi_low_=0.39, Swi_high_=0.104, 
                 m_low_=0.4, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
                 pd_low_=5000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
                 dPc=4250., // NOT SPECIFIED: capillary pressure difference between high_k and low_k layer at connate water saturation
                 VL_ux = 1.0e-5, CL_ux = 8e-8;
  
  // Krw_ave= C1_krw * Sw(ux)^C2_krw 
  double64 Krw( double64 Sw, double64 ux ) const {
  // =============================================
    }

  double64 Krnw( double64 Sw, double64 ux ) const {
  // ==============================================
    }
    
  // Krw_VL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
  double64 Krw_VL( double64 Sw, double64 ux ) const {
       const double64 C1(12.46), C2(0.3495);
       return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
    }
  // Krw_CL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
  double64 Krw_CL( double64 Sw, double64 ux ) const {
       const double64 C1(4.271), C2(1.241);
       return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
    }
  // Krnw_VL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
  double64 Krnw_VL( double64 Sw, double64 ux ) const {
       const double64 C1(2.736), C2(1.955), C3(3.586);
       return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
    }
  // Krnw_CL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
  double64 Krnw_CL( double64 Sw, double64 ux ) const {
       const double64 C1(13.15), C2(8.095), C3(2.653);
       return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
    }
};







/// composite - planar bedding fine sandstone - silt
struct CRC3_RockType11 {
 const std::string name = "P-FSst-Slt";  
 // Composite1 - hypothetical sample based on Achyut's rocktypes
 const int      rocktype_ = 11;
 const int      subtypes_ = 2;
 const double64 k_low_   = 3.4759e-14,
                k_high_  = 3.6075e-13,  // layer permeabilities
                LY_1     = 0.015, LY_2 = 0.01, LY_3 = 0.005,
                LY_4     = 0.015, LY_5 = 0.005,
                LY_low_  = 0.25,  LY_high_ = 0.25,    // cumulative layer thickness in the vertical direction (Y)
                Swi_low_ = 0.18,  Swi_high_ = 0.159, // irreducible saturations of the 2 different layers
                m_low_   = 0.5,   m_high_ = 0.6,       // van Genuchten exponents for the 2 different layers
                pd_low_  = 3000., pd_high_ = 1000.,  // capillary (drainage) entry pressure of low and high
                dPc      = 1.3061e+06, // capillary pressure difference between high_k and low_k layer at connate water saturation
                VL_ux    = 4.0e-4, CL_ux = 5.0e-7;

 /// Water relative permeability
 double64 Krw( double64 Sw, double64 ux ) const {
     if ( ux >= VL_ux ) return Krw_VL( Sw );
     if ( ux <= CL_ux ) return Krw_CL( Sw );
     
     const double64 a_c1 = 0.4499,
                    b_c1 = 4.609e+05,
                    c_c1 = 0.6106,
                    C1_krw = a_c1 * std::exp(-b_c1 * ux) + c_c1;
   
     const double64 p1_2 =  -2.837e+32,
                    p2_2 =   5.895e+28,
                    p3_2 =  -4.967e+24,
                    p4_2 =   2.174e+20,
                    p5_2 =  -5.235e+15,
                    p6_2 =    6.63e+10,
                    p7_2 =  -3.104e+05,
                    p8_2 =       5.662,
                    // efficient polynomial evaluation
                    // p=(((c[4]*x+c[3])*x+c[2])*x+c[1])*x+c[0]
                    // C2_krw= p1_2*ux^7 + p2_2*ux^6 + p3_2*ux^5 + p4_2*ux^4 + p5_2*ux^3 + p6_2*ux^2 + p7_2*ux + p8_2
                    C2_krw = (((((((p1_2 * ux + p2_2) * ux + p3_2) * ux + p4_2) * ux) + p5_2) * ux + p6_2) * ux + p7_2) * ux + p8_2;
     
     const double64 p1_3 =  -1.936e+30,
                    p2_3 =   3.925e+26,
                    p3_3 =  -3.169e+22,
                    p4_3 =   1.281e+18,
                    p5_3 =    -2.6e+13,
                    p6_3 =   1.913e+08,
                    p7_3 =        1733,
                    p8_3 =    -0.03431,
                    C3_krw = (((((((p1_3 * ux + p2_3) * ux + p3_3) * ux + p4_3) * ux) + p5_3) * ux + p6_3) * ux + p7_3) * ux + p8_3;
     // Krw_ave
     return C1_krw * std::pow(Sw,C2_krw) + C3_krw;
  }
 
 /// CO2 relative permeability:
 double64 Krnw( double64 Sw, double64 ux  ) const {
    if ( ux >= VL_ux ) return Krnw_VL( Sw );
    if ( ux <= CL_ux ) return Krnw_CL( Sw );
 
    const double64 a_c1 = 582.9,
                   b_c1 =  -1.416e+06,
                   c_c1 =  18.09,
                   C1_krnw = a_c1 * std::exp(b_c1 * ux) + c_c1;
      
    const double64 a_c2 = -0.08602,
                   b_c2 = -0.328,
                   c_c2 = -6.47,
                   C2_krnw = a_c2 * std::exp(b_c2 * ux) + c_c2;
      // Krnw_ave=
      return C1_krnw * exp(C2_krnw * Sw);
   }
     double64 Krw_VL( double64 Sw ) const {
    const double64 C1 = 6.342,
                   C2 = 0.4549;
    return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
 }
   
 double64 Krnw_VL( double64 Sw ) const {
    const double64 C1 = 1.491,
                   C2 = 1.945,
                   C3 = 2.472;
    // Krnw_VL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))               
    return C1 * (1. - std::pow(Sw,C2)) * (1. - std::pow(Sw,C3));
 }
   
double64 Krw_CL( double64 Sw ) const {
  double64 C1 = 2.717,
           C2 = 0.5099;
  if ( Sw < 0.83 ) {
       C1 = 7.39;
       C2 = 1.557;       
    }
  return std::sqrt(Sw) * (1. - std::pow(1. - std::pow(Sw,C1), C2));
}

double64 Krnw_CL( double64 Sw ) const {
  double64 C1 = 3.1, 
           C2 = 5.033, 
           C3 = 2.758; 
      if ( Sw < 0.3 ) {
           C1 =  2.199;
           C2 =  4.148;
           C3 = 30.1;       
        }
  // Krnw_CL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
  return C1 * (1. - std::pow(Sw,C2)) * (1. - std::pow(Sw,C3)); 
}

};



/// composite - Planar bedding Siltstone – Coarse Sandstone 
struct CRC3_RockType12 {
const std::string name = "P-Slt-CSst";  
const int      rocktype_ = 12;
const int      subtypes_ =  2;
const double64 k_low_=3.48e-14, k_high_=2.50e-12, 
               LY_1_=0.015, LY_2_=0.01, LY_3_=0.005, LY_4_=0.015, LY_5_=0.005,
               LY_low_=0.25, LY_high_=0.25,
               Swi_low_=0.18, Swi_high_=0.104, 
               m_low_=0.4, m_high_=0.7,        // van Genuchten exponents for the 2 different layers
               pd_low_=5000., pd_high_=750.,  // capillary (drainage) entry pressure of low and high
               dPc=4250., // NOT SPECIFIED: capillary pressure difference between high_k and low_k layer at connate water saturation
               VL_ux = 1.0e-5, CL_ux = 8e-8;

// Krw_ave= C1_krw * Sw(ux)^C2_krw 
double64 Krw( double64 Sw, double64 ux ) const {
// =============================================
  }

double64 Krnw( double64 Sw, double64 ux ) const {
// ==============================================
  }

  // Krw_VL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
  double64 Krw_VL( double64 Sw, double64 ux ) const {
       const double64 C1(6.459), C2(0.4515);
       return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
    }
  // Krw_CL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
  double64 Krw_CL( double64 Sw, double64 ux ) const {
       const double64 C1(3.962), C2(1.183);
       return std::sqrt(ux) * std::pow(1. - std::pow( 1.-Sw, C1 ), C2 );
    }
  // Krnw_VL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
  double64 Krnw_VL( double64 Sw, double64 ux ) const {
       const double64 C1(1.484), C2(1.909), C3(2.383);
       return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
    }
  // Krnw_CL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
  double64 Krnw_CL( double64 Sw, double64 ux ) const {
       const double64 C1(2.373), C2(5.23), C3(2.133);
       return C1 * ( std::pow(1. - Sw, C2) * std::pow(1. - Sw, C3) );
    }
};




// COMPOSITE ROCK-TYPES WHERE ONLY OF THE TWO TYPES IS FLOWING

/**
    M_CbSst_Mst      Massive bedding (Carbonate-cemented sandstone lamina and mudstone matrix) 
    M_CbSst_Slt        Massive bedding (Carbonate-cemented sandstone lamina and siltstone matrix) 
    M_CbSst_FSs      Massive bedding (Carbonate-cemented sandstone lamina and fine sandstone matrix) 
    M_CbSst_CSst     Massive bedding (Carbonate-cemented sandstone lamina and coarse sandstone matrix)
    
    In these rocks - in the best case - - ONLY ONE SUBTYPE FLOWS.

    We treat the carbonate-cemented rocks as singlets  because their entry pressure is greater than 10kP,
    which amounts to 1.5-m CO2 column height and our grid cells are only 5-cm tall.
        
    When k_low < 0.1mD or pd_low >=100kPa, we do not need the treatment as a composite anymore
    because non-wetting phase will never get in. We just model the permeability anisotropy. Note that
    100kPa is roughtly equivalent to a 15m CO2 column height which is >> than the height of our grid cells for which we need the composite rock types.
    
    Due to the presence of carbonate cement, these all are low permeability rocks with negligible porosity.
    Their relative permeability curves are treated as linear.
    Also, permeability anisotropy is ignored because even their maximum values are very low.
*/


/// M_CbSst_Mst   - 2 Massive bedding (Carbonate-cemented sandstone lamina and mudstone matrix) - ONLY ONE SUBTYPE FLOWS
/// M_CbSst_Slt     - 5 Massive bedding (Carbonate-cemented sandstone lamina and siltstone matrix) - ONLY ONE SUBTYPE FLOWS
/// M_CbSst_FSst -  7 Massive bedding (Carbonate-cemented sandstone lamina and fine sandstone matrix) - ONLY ONE SUBTYPE FLOWS
/// M_CbSst_CSst - 9 Massive bedding (Carbonate-cemented sandstone lamina and coarse sandstone matrix) - ONLY ONE SUBTYPE FLOWS


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
