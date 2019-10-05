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

struct CRC3_RockType1 {
  const std::string name = "H-Mst";        ///< homogeneous mudstone
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



struct CRC3_RockType3 {
  const std::string name = "H-CbSst";  ///< homogeneous carbonate-cemented sandstone
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




struct CRC3_RockType10 {
 const std::string name = "H-Slt";  ///< Homogeneous (Carbonate-cemented sandstone):
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
                LY_1_     = 0.05,      //  permeabilities individual layers
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









// COMPOSITE ROCK TYPES

struct CRC3_RockType4 {
 const std::string name = "P-Mst-Slt";  ///< Planar bedding mudstone - silt
 // Composite1 - hypothetical sample based on Achyut's rocktypes
 const int      rocktype_ = 4;
 const int      subtypes_ = 2;
 const double64 k_low_   = 1.73e-15, k_high_ = 3.48e-14,  // layer permeabilities
                LY_1     = 0.015,  LY_2 = 0.01, LY_3 = 0.005, 
                LY_4     = 0.015, LY_5 = 0.005,      // layers thicknesses (layer 1,3,5 are high_k layers, layer 2, 4 are low_k layers)
                LY_low_  = 0.25,  LY_high_ = 0.25,   // cumulative layer thickness in the vertical direction (Y)
                Swi_low_ = 0.39,  Swi_high_ = 0.18,  // irreducible saturations of the 2 different layers
                m_low_   = 0.4,   m_high_ = 0.5,     // van Genuchten exponents for the 2 different layers
                pd_low_  = 5000., pd_high_ = 3000.,  // capillary (drainage) entry pressure of low and high
                dPc      = 2000.; // capillary pressure difference between high_k and low_k layer at connate water saturation

 /// Water relative permeability
 double64 Krw( double64 Sw, double64 ux ) const {
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
 double64 Krnw( double64 Sw, double64 ux  ) const {
    const double64 a_c1 =   3.698e+04,
                   b_c1 =  -3.692e+07,
                   c_c1 =       585.6,
                   d_c1 =    3.72e+04,
                   C1_krnw = a_c1 * exp(b_c1 * ux) + c_c1 * exp(d_c1 * ux);
      
    const double64 a_c2 = -10.5,
                   b_c2 =  -1.278e+07,
                   c_c2 =      -11.44,
                   d_c2 =        5454,
                   C2_krnw = a_c2 * exp(b_c2 * ux) + c_c2 * exp(d_c2 * ux);
      // Krnw_ave=
      return C1_krnw * exp(C2_krnw * Sw);
   }
};






struct CRC3_RockType11 {
 const std::string name = "P-FSst-Slt";  ///< Planar bedding fine sandstone - silt
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
