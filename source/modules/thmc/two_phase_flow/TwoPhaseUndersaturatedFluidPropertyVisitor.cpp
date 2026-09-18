// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  TwoPhaseUndersaturatedFluidPropertyVisitor.cpp
//  CSMP-library
//
//  Created by Lukas Mosser on 10/30/13.
//

#include "TwoPhaseUndersaturatedFluidPropertyVisitor.h"
#include "Exception.h"
#include <algorithm>
#include "Node.h"

namespace csmp {
  
  //API_in = API of Oil
  //spGr_sp_in = Specific Gravity of Oil measured at the separator
  //Pb_in = Bubblepoint Pressure psia
  //Rsb_in = Solution gas ratio at the bubble point scf/ft3
  //my_o_b_in = Viscosity of oil at bubble point
  //rho_o_b_in = density of oil at bubble point psia
  //salinity_in = Salinity in  wt%
  template<uint32_t dim>
  TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::TwoPhaseUndersaturatedFluidPropertyVisitor( Model<dim>& model,
                                                                                              const FluidPropertyVisitorConfiguration& config )
  : Visitor<dim>(MODEL, NODE),
  model_ref_(model),
  prop_ref_(model.Database()),
  API(config.API),
  spGr_sp(config.spGr_sp),
  Pb(config.Pb),
  Rsb(config.Rsb),
  my_o_b(config.my_o_b),
  rho_o_b(config.rho_o_b),
  salinity(config.salinity),
  water_pressure_key_(model.Database().StorageKey(config.pres_w)),
  oil_pressure_key_(model.Database().StorageKey(config.pres_o)),
  temperature_key_(model.Database().StorageKey(config.temp)),
  water_visc_key_(model.Database().StorageKey(config.visc_w)),
  oil_visc_key_(model.Database().StorageKey(config.visc_o)),
  water_comp_key_(model.Database().StorageKey(config.comp_w)),
  oil_comp_key_(model.Database().StorageKey(config.comp_o)),
  water_density_key_(model.Database().StorageKey(config.dens_w)),
  oil_density_key_(model.Database().StorageKey(config.dens_o))
  {

    if ( water_pressure_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.pres_w, " must be a node property." );
    if ( oil_pressure_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.pres_o, " must be a node property." );
    if ( temperature_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.temp, " must be a node property." );
    if ( water_visc_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.visc_w, " must be a node property." );
    if ( oil_visc_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.visc_o, " must be a node property." );
    if ( water_comp_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.comp_w, " must be a node property." );
    if ( oil_comp_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.comp_o, " must be a node property." );
    if ( water_density_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.dens_w, " must be a node property." );
    if ( oil_density_key_.place != NODE )
      throw csmp::Exception( ERROR, "TwoPhaseUndersaturatedFluidPropertyVisitor", config.dens_o, " must be a node property." );
    
    assert(API>0.0);
    assert(spGr_sp>0.0);
    assert(Pb>0.0);
    assert(Rsb>=0.0);
    assert(my_o_b>=0.0);
    assert(rho_o_b>=0.0);
    
    //Initialization of coefficients
    double rho_w_pure_coeffs_c[5] = {-0.127213,   0.645486,   1.03265,        -0.070291,  0.639589    };
    double Ew_coeffs_c[5]         = {4.221,       -3.478,     6.221,          0.5182,     -0.4405     };
    double Fw_coeffs_c[5]         = {-11.403,     29.932,     27.952,         0.20684,    0.3768      };

    double dm_2_coeffs_c[5]       = {-1.1149E-4,  1.7105E-4,  -4.3766E-4,     0.0,        0.0         };
    double dm_3_2_coeffs_c[5]     = {-8.878E-4,   -1.388E-4,  -2.96318E-3,    0.0,        0.51103     };
    double dm_1_coeffs_c[5]       = {2.1466E-3,   1.2427E-2,  4.2648E-2,      -8.1009E-2, 0.525417    };
    double dm_1_2_coeffs_c[5]     = {2.356E-4, -  3.636E-4,   -2.278E-4,      0.0,        0.0         };

    double em_coeffs_c[5]         = {0.0,         0.0,        0.1249,         0.0,        0.0         };
    double fm_3_2_coeffs_c[5]     = {-0.617,      -0.747,     -0.4339,        0.0,        10.26       };
    double fm_1_coeffs_c[5]       = {0.0,         9.917,      5.1128,         0.0,        3.892       };
    double fm_1_2_coeffs_c[5]     = {0.0365,      -0.0369,    0.0,            0.0,        0.0         };

    double di_c[10] = {0.28853170E7, -0.11072577E5, -0.90834095E1, 0.30925651E-1, -0.27407100E-4, -0.192827250E7, 0.56216046E4, 0.13827250E2, -0.47609523E-1, 0.35545041E-4};


    std::copy(rho_w_pure_coeffs,    rho_w_pure_coeffs + 5,  rho_w_pure_coeffs_c );
    std::copy(Ew_coeffs,            Ew_coeffs + 5,          Ew_coeffs_c         );
    std::copy(Fw_coeffs,            Fw_coeffs + 5,          Fw_coeffs_c         );

    std::copy(dm_2_coeffs,          dm_2_coeffs+5,          dm_2_coeffs_c       );
    std::copy(dm_3_2_coeffs,        dm_3_2_coeffs+5,        dm_3_2_coeffs_c     );
    std::copy(dm_1_coeffs,          dm_1_coeffs+5,          dm_1_coeffs_c       );
    std::copy(dm_1_2_coeffs,        dm_1_2_coeffs+5,        dm_1_2_coeffs_c     );

    std::copy(em_coeffs,            em_coeffs+5,            em_coeffs_c         );
    std::copy(fm_3_2_coeffs,        fm_3_2_coeffs+5,        fm_3_2_coeffs_c     );
    std::copy(fm_1_coeffs,          fm_1_coeffs+5,          fm_1_coeffs_c       );
    std::copy(fm_1_2_coeffs,        fm_1_2_coeffs+5,        fm_1_2_coeffs_c     );

    std::copy(di, di+10, di_c);

    //Convert Salinity to Molality
    M = SalinityToConcentration(salinity);
    //Calculate the Variables that only need to be calculated once i.e. independent on P,T,S
    CalculateConstantVariables();

  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::CalculateConstantVariables()
  {
    double log_API, log_spGr_SP, log_pB, log_Rsb;
    log_API = log(API);
    log_spGr_SP = log(spGr_sp);
    log_pB = log(Pb);
    log_Rsb = log(Rsb);
    
    //cout << log_API << " " <<log_spGr_SP << " " << log_pB << " " << log_Rsb << endl;
    Z1 = 3.011-2.6254*log_API+0.497*pow(log_API,2);
    Z2 = -0.0835-0.259*log_spGr_SP+0.382*pow(log_spGr_SP,2);
    Z3 = 3.51-0.0289*log_pB-0.0584*pow(log_pB,2);
    Z5 = -1.918-0.642*log_Rsb+0.154*pow(log_Rsb, 2);
    
    A = -1.0146+1.3322*log(my_o_b)-0.4876*pow(log(my_o_b), 2)-1.15036*pow(log(my_o_b), 3);
    //cout << log(my_o_b) << " " << my_o_b << " " << Z1 << " " << Z2 << " " << Z3 << " " << Z5 << " " << A << endl;
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyOilCompressibilityCorrelation(double P, double T)
  {
    
    double log_P_Pb, log_Tr, dz_dp, dcofb_dp;
    Tr = T;
    log_P_Pb = log(P/Pb);
    log_Tr = log(Tr);
    
    Z4 =0.327-0.608*log_P_Pb+0.0911*pow(log_P_Pb, 2);
    Z6 = 2.52-2.73*log_Tr+0.429*pow(log_Tr, 2);
    
    /*Z = 3.011-2.6254*log_API+0.497*pow(log_API,2)
     -0.0835-0.259*log_spGrSP+0.382*pow(log_spGrSP,2)
     +3.51-0.0289*log_pB-0.0584*pow(log_pB,2)
     +0.327-0.608*log_P_Pb+0.0911*pow(log_P_Pb, 2)
     -1.918-0.642*log_Rsb+0.154*pow(log_Rsb, 2)
     +2.52-2.73*log_Tr+0.429*pow(log_Tr, 2);*/
    
    Z = Z1 + Z2 + Z3 + Z4 + Z5 + Z6;
    //cout << " " << Z1 << " " << Z2 << " " << Z3 << " " << Z4 << " " << Z5 << " " << Z6 << " " << endl;
    cofb = exp(2.434+0.475*Z+0.048*pow(Z, 2)-log(1E6));
    
    dz_dp = (-0.608+0.1822*log_P_Pb)/P;
    dcofb_dp = cofb*(0.475+0.096*Z)*dz_dp;
    
    co = cofb + (P-Pb)*dcofb_dp;
  }
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::DensityOilGreaterPb(double P)
  {
    rho_o = rho_o_b*exp(cofb*(P-Pb));
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::PetroskyFarshadOilViscosityGreaterPb(double P)
  {
    my_o = my_o_b + 1.3449E-3*(P-Pb)*pow(10, A);
  }
  
  template<uint32_t dim>
  double TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyA(const double ai[], double T) const
  {
    return (ai[0]*pow(T/100,2)+ai[1]*(T/100)+ai[2])/(ai[3]*pow(T/100,2)+ai[4]*(T/100)+1);
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyBrineAFactors(double T)
  {
    rho_w_pure_p_ref = SpiveyA(Ew_coeffs, T);
    Ew = SpiveyA(Ew_coeffs, T);
    Fw = SpiveyA(Fw_coeffs, T);
    
    dm_2 = SpiveyA(dm_2_coeffs, T);
    dm_3_2 = SpiveyA(dm_3_2_coeffs, T);
    dm_1 = SpiveyA(dm_1_coeffs, T);
    dm_1_2 = SpiveyA(dm_1_2_coeffs, T);
    
    Em = SpiveyA(em_coeffs, T);
    Fm_3_2 = SpiveyA(fm_3_2_coeffs, T);
    Fm_1 = SpiveyA(fm_1_coeffs, T);
    Fm_1_2 = SpiveyA(fm_1_2_coeffs, T);
    
    //cout << rho_w_pure_p_ref << " " << Ew << " " << Fw << " " <<  dm_2 << " " <<  dm_3_2 << " " <<  dm_1 << " " <<  dm_1_2 << " " <<  Em << " " <<  Fm_3_2 << " " <<  Fm_1 << " " <<  Fm_1_2 <<endl;
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyBrineCompressibilityPure(double P)
  {
    c_w_pure = (1.0/70.0)*1.0/(Ew*(P/70.0)+Fw);
    //cout << c_w_pure << endl;
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyBrineDensityPure(double P)
  {
    Iw_p_ref = (1.0/Ew)*log(fabs(Ew+Fw));
    Iw = (1.0/Ew)*log(fabs(Ew*(P/70.0)+Fw));
    rho_w_pure = rho_w_pure_p_ref*exp(Iw-Iw_p_ref);
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyBrindeDensityPRef(double m)
  {
    rho_b_p_ref = rho_w_pure_p_ref + dm_2*pow(m, 2) + dm_3_2*pow(m, 1.5) + dm_1*m + dm_1_2*sqrt(m);
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyBrineCompressibility(double P, double m)
  {
    Eb = Ew + Em;
    Fb = Fw + Fm_3_2*pow(m, 1.5) + Fm_1 * m + Fm_1_2 * sqrt(m);
    c_b = (1.0/70.0)*(1.0/(Eb*(P/70.0) + Fb));
    //cout << c_b << endl;
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SpiveyBrineDensity(double P)
  {
    Ib_p_ref = 1.0/Eb * log(fabs(Eb+Fb));
    Ib = (1.0/Eb)*log(fabs(Eb*(P/70.0)+Fb));
    rho_b = rho_b_p_ref*exp(Ib - Ib_p_ref);
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::MaoDuanWaterViscosityPure(double T)
  {
    double a = 0;
    double b = 0;
    for(int i = 1; i <= 5; ++i) {
      a += di[i-1]*pow(T, i-3);
    }
    
    for(int i = 6; i <= 10; ++i) {
      b += di[i-1]*pow(T, i-8);
    }
    
    log_my_w_pure = a+b;
    my_w_pure = exp(log_my_w_pure);
  }
  
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::MaoDuanBrineRelativeViscosity( double T, double m ) {
    const double Tsquared(T*T);
    const double Ahere = -0.21319213+0.13651589E-2*T - 0.12191756E-5*Tsquared;
    const double B = 0.69161945E-1 - 0.27292263E-3*T + 0.20852448E-6*Tsquared;
    const double C = -0.25988855E-2 + 0.77989227E-5*T;
    
    log_my_b_rel = Ahere*m + B*pow(m, 2) + C*pow(m, 3);
    my_b_rel = exp(log_my_b_rel);
    my_b = my_b_rel * my_w_pure;
  }
  
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::Visit( Node<dim>* node )
  {
    const double Pw = node->Read(water_pressure_key_ );
    const double Po = node->Read(oil_pressure_key_ );
    const double T = node->Read(temperature_key_ );
    
// NOT USED   double P_w_field = Pw * 1E-6 * 1.45037738E2; // Pa -> MPa -> PSIA
    double P_w_SI = Pw * 1E-6; //MPa
    double P_o_field = Po * 1E-6 * 1.45037738E2; // Pa -> MPa -> PSIA;
// NOT USED   double P_o_SI = Po * 1E-6; //MPa
    double T_field = T * 1.8 + 32; //T [F]
    double T_cels = T; //T [Celsius]
    double T_SI = T + 273.15; // T [Kelvin]
    
    //Calculate Oil Properties
    SpiveyOilCompressibilityCorrelation(P_o_field, T_field);
    DensityOilGreaterPb(P_o_field);
    PetroskyFarshadOilViscosityGreaterPb(P_o_field);
    
    SpiveyBrineAFactors(T_cels);
    SpiveyBrineCompressibilityPure(P_w_SI);
    SpiveyBrineDensityPure(P_w_SI);
    SpiveyBrindeDensityPRef(M);
    SpiveyBrineCompressibility(P_w_SI,M);
    SpiveyBrineDensity(P_w_SI);
    MaoDuanWaterViscosityPure(T_SI);
    MaoDuanBrineRelativeViscosity(T_SI, M);
    MaoDuanBrineRelativeViscosity(T_SI, M);

// to fix:  use makeScalar() if scalar is to be used only once
    const ScalarVariable oil_viscosity( PLAIN, my_o*1E-3); //cP->Pa*s
    const ScalarVariable brine_viscosity( PLAIN, my_b); //Pa*s
    const ScalarVariable oil_compressibility( PLAIN, co/6897.0); //PSIA-1 -> Pa-1
    const ScalarVariable brine_compressibility( PLAIN, c_b*1.0E-6); //MPa-1 -> Pa-1 //IS THIS CORRECT!!!
    const ScalarVariable oil_density(PLAIN, rho_o*1.60184634E-2*1000.); //ft3/lbm->g/cm3->kg/m3
    const ScalarVariable brine_density(PLAIN, rho_b*1000.); // g/cm3 -> kg/m3
    
    node->Store(water_visc_key_, brine_viscosity);
    node->Store(oil_visc_key_, oil_viscosity);
    node->Store(oil_comp_key_,oil_compressibility);
    node->Store(water_comp_key_, brine_compressibility);
    node->Store(oil_density_key_, oil_density);
    node->Store(water_density_key_, brine_density);
  }
  
  template<uint32_t dim>
  void TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::Visit(Model<dim>* model)
  {
    std::cout << "\nTwoPhaseUndersaturatedFluidPropertyVisitor::Visit(Model): Visiting model.\n";
    model->Accept( *this );
    std::cout << "\nTwoPhaseUndersaturatedFluidPropertyVisitor::Visit(Model): Values Calculated Successfuly!\n";
  }
  
  
  template class TwoPhaseUndersaturatedFluidPropertyVisitor<1U>;
  template class TwoPhaseUndersaturatedFluidPropertyVisitor<2U>;
  template class TwoPhaseUndersaturatedFluidPropertyVisitor<3U>;
}
