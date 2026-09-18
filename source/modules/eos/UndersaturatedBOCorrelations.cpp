// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  UndersaturatedBOCorrelations.cpp
//  CSMPLibraryTestProject
//
//  Created by Lukas Mosser on 7/9/14.
//

#include "UndersaturatedBOCorrelations.h"

namespace csmp {

/** StandingRSCorrelation constructor
 
Standard constructor for correlation type class.
Precomputes set of parameters for later reuse.
 
@param b An input struct of type BOInputProperties containing the required input properties.

@test Tested OK
*/
StandingRS_Correlation::StandingRS_Correlation(BO_InputProperties& b)
 {
    assert(b.gamma_API > 0);
    psi_ = pow(10, 0.0125*b.gamma_API);
 }

/** StandingRSCorrelation compute method
 
 Implementation of CorrelationInterface compute method.
 Calculates pressure and temperature dependent Solution gas oil ratio.
 Converts input pressure and temperature from SI unit to field units.
 Writes results to a struct of type BODynamicProperties.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 @param P Oil phase pressure in SI unit Pa.
 @param T Oil phase temperature in degree C.
 @param d An output struct of type BOInputProperties containing the calculated output property.
 
 @test Tested OK
 */
void StandingRS_Correlation::Compute(BO_InputProperties& b, double P, double T, BO_DynamicProperties& d)
 {
    assert(P >= 100325.);
    assert(T > 0.);
    assert(b.gamma_g > 0.);
    
    P *=1.4504e-04;
    T *= 1.8;
    T += 32;
    d.Rs = b.gamma_g*pow(((0.055*P + 1.4)*psi_)*pow(10, -0.00091*T), 1.205);
 }



/** StandingOilDensityCorrelation constructor
 
 Standard constructor for correlation type class.
 Precomputes set of parameters for later reuse.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 
 @test Tested OK
 */
StandingOilDensityCorrelation::StandingOilDensityCorrelation(BO_InputProperties& b)
 {
    assert(b.gamma_o > 0.);
    psi_ = 62.4*b.gamma_o;
 }



/** StandingRSCorrelation compute method
 
 Implementation of CorrelationInterface compute method.
 Calculates pressure and temperature dependent oil density.
 Converts input pressure and temperature from SI unit to field units.
 Writes results to a struct of type BODynamicProperties.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 @param P Oil phase pressure in SI unit Pa.
 @param T Oil phase temperature in degree C.
 @param d An output struct of type BOInputProperties containing the calculated output property.
 
 @attention Requires an updates Rs value for current temperature and pressure conditions.
 Always compute Rs before the oil density as is seen in BOPropertyCalculatorImplementation.
 There calling update will always call the Rs correlation first to update this property to new P-T conditions.
 
 @test Tested OK
 */
void StandingOilDensityCorrelation::Compute(BO_InputProperties &b, double, double, BO_DynamicProperties& d )
 {
    d.rho_o = ((psi_ + 0.0136*b.gamma_g*d.Rs)/d.Bo);
    d.rho_o_SI = d.rho_o*16.018463;
 }

/** VasquezBeggsUndersaturatedOilCompressibilityCorrelationImplementation constructor
 
 Standard constructor for correlation type class.
 Precomputes set of parameters for later reuse.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 
 @test Tested OK
 */
VasquezBeggsOilCompressibilityCorrelation::VasquezBeggsOilCompressibilityCorrelation(BO_InputProperties& b)
 {
    assert(b.gamma_g > 0.);
    assert(b.gamma_API > 0.);
    assert(b.T_sep > 0.);
    assert(b.P_sep >= 100325.);
    b.gamma_gc = b.gamma_g*(1+0.5912e-04*b.gamma_API*b.T_sep*log10(b.P_sep/114.7));
 }



/** VasquezBeggsOilCompressibilityCorrelation compute method
 
 Implementation of CorrelationInterface compute method.
 Calculates pressure and temperature dependent oil compressibility.
 Converts input pressure and temperature from SI unit to field units.
 Writes results to a struct of type BODynamicProperties.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 @param P Oil phase pressure in SI unit Pa.
 @param T Oil phase temperature in degree C.
 @param d An output struct of type BOInputProperties containing the calculated output property.
 
 @test Tested OK
 */
void VasquezBeggsOilCompressibilityCorrelation::Compute(BO_InputProperties& b, double P, double T, BO_DynamicProperties& d)
 {
    assert(b.Rsb > 0.);
    assert(b.Bob > 0.);
    assert(b.Pb > 0.);
    assert(b.gamma_API > 0.);
    assert(b.gamma_gc > 0.);
    assert(P >= 100325.);
    assert(T > 0.);
    
    P *=1.4504e-04;
    
    assert(P > b.Pb); //check for undersaturated condition
    
    T *= 1.8;
    T += 32;
    a_ =(1e-5*(5*b.Rsb + 17.2*T - 1180.0*b.gamma_gc + 12.61*b.gamma_API -1433));
    d.c_o = a_/P;
    d.c_o_SI = d.c_o/6897.0;
    d.Bo = b.Bob*pow((b.Pb/P),a_);
 }



/** AbdulMajeedOilViscosityCorrelation constructor
 
 Standard constructor for correlation type class.
 Precomputes set of parameters for later reuse.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 
 @test Tested OK
 */

AbdulMajeedOilViscosityCorrelation::AbdulMajeedOilViscosityCorrelation(const BO_InputProperties& b)
 {
    assert(b.gamma_API > 0.);
    alpha_ = 1.9311 - 0.001104*b.gamma_API*b.gamma_API - 5.2106;
 }

/** AbdulMajeedOilViscosityCorrelation compute method
 
 Implementation of CorrelationInterface compute method.
 Calculates pressure and temperature dependent oil viscosity.
 Converts input pressure and temperature from SI unit to field units.
 Writes results to a struct of type BODynamicProperties.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 @param P Oil phase pressure in SI unit Pa.
 @param T Oil phase temperature in degree C.
 @param d An output struct of type BOInputProperties containing the calculated output property.
 
 @attention Requires an updates Rs value for current temperature and pressure conditions.
 Always compute Rs before the oil density as is seen in BOPropertyCalculator.
 There calling update will always call the Rs correlation first to update this property to new P-T conditions.
 
 @test Tested OK
 */
void AbdulMajeedOilViscosityCorrelation::Compute(BO_InputProperties &b, double P, double, BO_DynamicProperties&d )
 {
    assert(d.Rs > 0.);
    assert(b.gamma_API > 0.);
    assert(b.my_ob > 0.);
    assert(b.Pb > 0.);
    
    P *=1.4504e-04;
    
    assert(P > b.Pb); //check for undersaturated
    
    ln_Rs_ = log(d.Rs);
    A_ = alpha_ - 0.89941*ln_Rs_ + 0.0092545*b.gamma_API*ln_Rs_;
    d.my_o = b.my_ob + pow(10.0, A_ + 1.11*log10(P-b.Pb));
    d.my_o_SI = d.my_o*0.001;
 }


/** BOPropertyCalculator Constructor
 
 Constructor for BOPropertyCalculator Implementation.
 Aggregates a number of correlations and initializes them with 
 a set of predefined oil fluid properties wrapped in the input 
 struct of BOInputProperties.
 
 @param b An input struct of type BOInputProperties containing the required input properties.
 
 @test Tested OK
 */
BO_PropertyCalculator::BO_PropertyCalculator( const BO_InputProperties& b )
: b_(b)
 {
    correlations_.push_back(new StandingRS_Correlation(b_));
    correlations_.push_back(new VasquezBeggsOilCompressibilityCorrelation(b_));
    correlations_.push_back(new StandingOilDensityCorrelation(b_));
    correlations_.push_back(new AbdulMajeedOilViscosityCorrelation(b_));
 }
 
 
 BO_PropertyCalculator::~BO_PropertyCalculator()
  {
     for ( auto& it : correlations_ ) {
          delete it;
          it = nullptr;
       }
 }

 

/** BOPropertyCalculator compute method
 
 Implementation of the BOPropertyCalculatorInterface compute method.
 Calculates pressure and temperature dependent variables from aggregated correlations.
 Properties can then be acceses via the accessor methods exposed in the BOPropertyCalculatorInterface.

 @param P Oil phase pressure in SI unit Pa.
 @param T Oil phase temperature in degree C.
 
 @test Tested OK
 */
void BO_PropertyCalculator::Update(double P, double T)
 {
    assert(P*1.4505e-4 > b_.Pb);
    assert(P >= 100325.);
    assert(T > 0.);
    
    for (std::vector<CorrelationInterface*>::iterator it = correlations_.begin();
         it != correlations_.end();
         ++it)
    {
      (*it)->Compute(b_, P, T, d_);
    }
 }

/** BOPropertyCalculator out method
 
 Implementation of the BOPropertyCalculatorInterface out method.
 Prints all calculated parameters exposed by the BOPropertyCalculatorInterface
 to console.
 
 @test Tested OK
 */

void BO_PropertyCalculator::Out() const
 {
    std::cout<<
    getUndersaturatedOilDensity()<< " " <<
    getUndersaturatedOilCompressibility() << " " <<
    getUndersaturatedOilViscosity()
    << std::endl;
 }

} //end csmp
