//
//  TwoPhaseUndersaturatedFluidPropertyVisitor.h
//  CSMP-library
//
//  Created by Lukas Mosser on 10/30/13.
//  Copyright (c) 2013 Lukas Mosser. All rights reserved.
//

#ifndef CSMP_library_TwoPhaseUndersaturatedFluidPropertyVisitor_h
#define CSMP_library_TwoPhaseUndersaturatedFluidPropertyVisitor_h

#include "Model.h"
#include "Visitor.h"
#include <vector>

namespace csmp {
 
  /**
   @author Lukas Mosser
   @date 2013
   */
  
template<uint32_t> class Model;

  template<uint32_t> class Model;
  
  typedef struct {
    const char* pres_w;
    const char* pres_o;
    const char* temp;
    const char* dens_w;
    const char* dens_o;
    const char* comp_w;
    const char* comp_o;
    const char* visc_w;
    const char* visc_o;
    
    double API;
    double spGr_sp;
    double Pb;
    double Rsb;
    double my_o_b;
    double rho_o_b;
    double salinity;
  } FluidPropertyVisitorConfiguration;
  
  
  template<uint32_t dim>
  class TwoPhaseUndersaturatedFluidPropertyVisitor : public Visitor<dim> {
  public:
    
    TwoPhaseUndersaturatedFluidPropertyVisitor(
                                               Model<dim>& model,
                                               const FluidPropertyVisitorConfiguration& config);
    
    virtual ~TwoPhaseUndersaturatedFluidPropertyVisitor() {}
    
    virtual void Visit(Node<dim>* node);
    virtual void Visit(Model<dim>* model);
    
    void CalculateConstantVariables();
    
    void SpiveyOilCompressibilityCorrelation(double P, double T);
    
    void DensityOilGreaterPb(double P);
    
    void PetroskyFarshadOilViscosityGreaterPb(double P) ;
    
    //T [C]
    void SpiveyBrineAFactors(double T);
    
    //T [C]
    double SpiveyA(const double ai[], double T) const;
    
    //P[MPa], T [C]
    void SpiveyBrineCompressibilityPure(double P);
    
    //P[MPa], T [C]
    void SpiveyBrineDensityPure(double P);
    
    //P[MPa], T [C], m [g-mol/kg H20 NaCl]
    void SpiveyBrindeDensityPRef(double m);
    
    //P[MPa], T [C], m [g-mol/kg H20 NaCl]
    void SpiveyBrineCompressibility(double P, double m);
    
    //P[MPa], T [C], m [g-mol/kg H20 NaCl]
    void SpiveyBrineDensity(double P);
    
    void MaoDuanWaterViscosityPure(double T);
    
    //T [Kelvin], m [g-mol/kg H20 NaCl]
    //Brine Rel. Viscosity Pa*s
    void MaoDuanBrineRelativeViscosity(double T, double m);
    
    double SalinityToConcentration(double S) const;
  
  private:
    const PropertyDatabase<dim>&  prop_ref_;
    Model< dim>&                  model_ref_;
  
    csmp::Index   water_pressure_key_,oil_pressure_key_, temperature_key_, water_density_key_,
                  oil_density_key_, water_visc_key_, oil_visc_key_, water_comp_key_, oil_comp_key_;
  
    //Input properties
    double API, Rs, spGr_gas, spGr_STO, Pb, Tr, Rsb, spGr_sp, rho_o_b, my_o_b, salinity;
  
    //Variables only calculated once
    //Variables calculated once for spivey oil compressibility
    double Z1, Z2, Z3, Z5;
    //Variables calculated once for Petrosky oil viscosity
    double A;
  
    //Variables calculated for each pressure/temperature for each element/node
    double Z,Z4,Z6, cofb, co, rho_o,my_o;
 
/// SKM: MOVE COEFFICIENTS INTO *.cpp FILE - CONSTRUCTOR ARGUMENT LIST
/// Roman, 2013: Initialization of const arrays was not allowed in C++ before C++11 standard due to unique definition rule
///              keep non-const arrays until the constructor argument list will be implemented

//    const double rho_w_pure_coeffs[5] = {-0.127213, 0.645486, 1.03265, -0.070291, 0.639589};
//    const double Ew_coeffs[5] = {4.221, -3.478, 6.221, 0.5182, -0.4405};
//    const double Fw_coeffs[5] = {-11.403, 29.932, 27.952, 0.20684, 0.3768};
    double rho_w_pure_coeffs[5];
    double Ew_coeffs[5];
    double Fw_coeffs[5];

    double rho_w_pure_p_ref, Ew, Fw, Iw_p_ref, Iw;
    double c_w_pure, rho_w_pure;
  
//    const double dm_2_coeffs[5] = {-1.1149E-4, 1.7105E-4, -4.3766E-4, 0.0, 0.0};
//    const double dm_3_2_coeffs[5] = {-8.878E-4, -1.388E-4, -2.96318E-3, 0.0, 0.51103};
//    const double dm_1_coeffs[5] = {2.1466E-3, 1.2427E-2, 4.2648E-2, -8.1009E-2, 0.525417};
//    const double dm_1_2_coeffs[5] = {2.356E-4, -3.636E-4, -2.278E-4, 0.0, 0.0};
    double dm_2_coeffs[5];
    double dm_3_2_coeffs[5];
    double dm_1_coeffs[5];
    double dm_1_2_coeffs[5];

    double M, dm_2, dm_3_2, dm_1, dm_1_2;
    double Eb, Fb;
  
    double rho_b_p_ref;
  
//    const double em_coeffs[5] = {0.0, 0.0, 0.1249, 0.0, 0.0};
//    const double fm_3_2_coeffs[5] = {-0.617, -0.747, -0.4339, 0.0, 10.26};
//    const double fm_1_coeffs[5] = {0.0, 9.917, 5.1128, 0.0, 3.892};
//    const double fm_1_2_coeffs[5] = {0.0365, -0.0369, 0.0, 0.0, 0.0};
    double em_coeffs[5];
    double fm_3_2_coeffs[5];
    double fm_1_coeffs[5];
    double fm_1_2_coeffs[5];

    double Em, Fm_3_2, Fm_1, Fm_1_2;
    double Ib_p_ref, Ib;
    double c_b, rho_b;
  
//    const double di[10] = {0.28853170E7, -0.11072577E5, -0.90834095E1, 0.30925651E-1, -0.27407100E-4, -0.192827250E7, 0.56216046E4, 0.13827250E2, -0.47609523E-1, 0.35545041E-4};
    double di[10];

    double log_my_w_pure, my_w_pure, log_my_b_rel, my_b_rel, my_b;
};

  
  ///Spivey et. al 1999 Compressibility of oil at pressures higher than Pb
  
template<uint32_t dim>
inline double TwoPhaseUndersaturatedFluidPropertyVisitor<dim>::SalinityToConcentration( double S ) const {
    return 1000.0*S/(58.4428*(1-S));
 }
  

  
}


#endif
