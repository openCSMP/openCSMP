//
//  UndersaturatedBOCorrelations.h
//  CSMPLibraryTestProject
//
//  Created by Lukas Mosser on 7/9/14.
//  Copyright (c) 2014 Lukas Mosser. All rights reserved.
//

#ifndef CSMP_UNDERSATURATED_BLACK_OIL_CORRELATIONS_H
#define CSMP_UNDERSATURATED_BLACK_OIL_CORRELATIONS_H

#include "CSMP_definitions.h"

namespace csmp {

/**
 @struct BO_InputProperties UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Struct handles all properties required for current implementation of 
 Black Oil input fluid properties for undersaturated oil reservoirs.
 Properties are in field units as specified in the member descriptions.
 
 @section design Design Intent
 Wraps all required data for use with current implementation of 
 Black Oil fluid propery correlation calculations.
 */
struct BO_InputProperties
{
  double gamma_g;     ///<Specific reservoir gas gravity
  double gamma_API;   ///<API gravity of reservoir oil
  double gamma_o;     ///<specific gravity of reservoir oil
  double gamma_gc;    ///<Vasquez Beggs corrected gas gravity
  double Pb;          ///<Bubblepoint pressure of the reservoir fluid. Unit: psia
  double Rsb;         ///<Solution gas oil ratio at bubble point pressure. Unit scf/STB
  double Bob;         ///<Formation Volume Factor of oil phase at bubble point pressure bbl/STB
  double my_ob;       ///<Oil viscosity at bubble point pressure. Unit: centipoise
  double P_sep;       ///<Seperator Pressure Unit: psia
  double T_sep;       ///<Seperator Temperature Unit: degree F
};


/**
 @struct BO_DynamicProperties UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Struct handles all properties calculated from current implementation of
 Black Oil input fluid properties for undersaturated oil reservoirs.
 Properties are in field units as specified in the member descriptions.
 Simulation relevant properties are stored in SI units accordingly.
 
 @section design Design Intent
 Wraps all provided data from current implementation of
 Black Oil fluid property correlation calculations.
 */
struct BO_DynamicProperties
{
  double Rs;        ///<Solution gas oil ratio at reservoir pressure and temperature Unit: scf/STB
  double Bo;        ///<Formation Volume Factor of oil phase at reservoir pressure and temperature Unit: bbl/STB
  double rho_o;     ///<Density of oil phase at reservoir temperature and pressure Unit: lbm/ft3
  double rho_o_SI;  ///<Density of oil phase at reservoir temperature and pressure. Unit: kg/m3
  double c_o;       ///<Oil phase compressibility at reservoir temperature and pressure. Unit: 1/psia
  double c_o_SI;    ///<Oil phase compressibility at reservoir temperature and pressure. Unit: 1/Pa
  double my_o;      ///<Reservoir oil viscosity at reservoir temperature and pressure. Unit: centipoise
  double my_o_SI;   ///<Reservoir oil viscosity at reservoir temperature and pressure. Unit: Pa*s
};


/**
 @class CorrelationInterface UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
Interface describing the basic methods required by all correlations
type classes. Correlation classes should inherit and implement the CorrelationInterface
class to allow interoperability with preexisting code.
 
 @section design Design Intent
Allows interoperability with future implementations of correlations.
Exposes a minimal interface to allow implementation of interoperable code
without exposing existing implementation information.
 */
class CorrelationInterface
{
public:
  virtual ~CorrelationInterface(){}
  virtual void Compute( BO_InputProperties&, double P, double T, BO_DynamicProperties& ) = 0;
};


/**
 @class StandingRSCorrelation UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
Implementation of Solution gas oil ratio based on a correlation published by Standing et. al.
 
 @section design Design Intent
Inherits from CorrelationInterface to allow interoperability.
Converts input SI units for pressure temperature.
Precomputes values in the constructor save computational cost.
 */
class StandingRS_Correlation : public CorrelationInterface
{
public:
  StandingRS_Correlation(BO_InputProperties& );
  void Compute(BO_InputProperties&, double P, double T, BO_DynamicProperties& ) override final;

private:
  double psi_;
};



/**
@class StandingOilDensityCorrelation UndersaturatedBOCorrelations.h

@author L.J. Mosser
@section motivation Motivation

Implementation of pressure and temperature dependent oil density
based on a correlation published by Standing et. al.
Valid for undersaturated regimes P > Pb.

@section design Design Intent
Inherits from CorrelationInterface to allow interoperability.
Converts input SI units for pressure temperature.
Precomputes values in the constructor save computational cost.
*/
class StandingOilDensityCorrelation : public CorrelationInterface
{
public:
  StandingOilDensityCorrelation(BO_InputProperties& b);
  void Compute( BO_InputProperties&, double P, double T, BO_DynamicProperties& ) override final;
  
private:
  double psi_;
};



/**
 @class VasquezBeggsOilCompressibilityCorrelation UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Implementation of pressure and temperature dependent oil compressibility
 based on a correlation published by Vasquey and Beggs.
 Valid for undersaturated regimes P > Pb.
 
 @section design Design Intent
 Inherits from CorrelationInterface to allow interoperability.
 Converts input SI units for pressure temperature.
 Precomputes values in the constructor save computational cost.
 */
class VasquezBeggsOilCompressibilityCorrelation : public CorrelationInterface
{
public:
  VasquezBeggsOilCompressibilityCorrelation(BO_InputProperties& b);
  void Compute( BO_InputProperties&, double P, double T, BO_DynamicProperties& ) override final;

private:
  double a_;
};



/**
 @class AbdulMajeedOilViscosityCorrelation UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Implementation of pressure and temperature dependent oil viscosity
 based on a correlation published by Abdul and Majeed.
 Only valid for undersaturated regime.
 
 @section design Design Intent
 Inherits from CorrelationInterface to allow interoperability.
 Converts input SI units for pressure temperature.
 Precomputes values in the constructor save computational cost.
 */
class AbdulMajeedOilViscosityCorrelation : public CorrelationInterface
{
public:
  AbdulMajeedOilViscosityCorrelation(const BO_InputProperties& b);
  void Compute( BO_InputProperties&, double P, double T, BO_DynamicProperties& ) override final;
  
private:
  double alpha_, A_, ln_Rs_;
};



/**
 @class PropertyCalculatorInterface UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Interface describing the basic methods required by all PropertyCalculators.
 PropertyCalculator classes should inherit and implement the PropertyCalculatorInterface
 class to allow interoperability with preexisting and future code implementations
 Allows for at runtime exchange of different black oil property calculators
 implementing a different set of correlations.
 
 @section design Design Intent
 Allows interoperability with future implementations of PropertyCalculators.
 Exposes a minimal interface to allow implementation of interoperable code
 without exposing existing implementation information.
 */
class PropertyCalculatorInterface
{
public:
  virtual ~PropertyCalculatorInterface(){}
  virtual void Update( double P, double T ) = 0;
  virtual void Out() const = 0;
};


/**
 @class ReservoirWaterPropertyCalculatorInterface UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Interface describing the basic methods required by all ReservoirWaterProperty calculators.
 Correlation classes should inherit and implement the PropertyCalculatorInterface
 class to allow interoperability with preexisting and future code implementations
 Allows for at runtime exchange of different black oil property calculators
 implementing a different set of correlations.
 
 @section design Design Intent
 Allows interoperability with future implementations of ReservoirWaterPropertyCalculatorInterface.
 Exposes a minimal interface to allow implementation of interoperable code
 without exposing existing implementation information.
 */
class ReservoirWaterPropertyCalculatorInterface
: public PropertyCalculatorInterface
{
public:
  virtual ~ReservoirWaterPropertyCalculatorInterface(){}
  virtual double getReservoirWaterDensity() const = 0;
  virtual double getReservoirWaterCompressibility() const = 0;
  virtual double getReservoirWaterViscosity() const = 0;
};



/**
 @class BO_PropertyCalculatorInterface UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Interface describing the basic methods required by all BOProperty calculators.
 Correlation classes should inherit and implement the CorrelationInterface
 class to allow interoperability with preexisting and future code implementations
 Allows for at runtime exchange of different black oil property calculators
 implementing a different set of correlations.
 
 @section design Design Intent
 
 Allows interoperability with future implementations of BOPropertyCalculatorImplementations.
 Exposes a minimal interface to allow implementation of interoperable code
 without exposing existing implementation information.
 */
class BO_PropertyCalculatorInterface
: public PropertyCalculatorInterface
{
public:
  virtual ~BO_PropertyCalculatorInterface(){}
  virtual double getUndersaturatedOilDensity() const = 0;
  virtual double getUndersaturatedOilCompressibility() const = 0;
  virtual double getUndersaturatedOilViscosity() const = 0;
};



/**
 @class BOPropertyCalculator UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 @section motivation Motivation
 
 Basic implementation of BOPropertyCalculatorInterface for use with temperature and pressure
 dependent formulation of black oil properties. Provides oil properties to 
 classes requiring fluid properties for undersaturated oil phases.
 Based on Standings correlations for Solution gas oil ratio and oil density.
 Implements Vasquez-Beggs correlation for compressibility.
 Implements Abdul-Majeed correlation for undersaturated oil viscosity.
 
 Has been tested for speed. 
 Current implementation allows for 1e7 update calls per second
 on Core i7 2.4 Ghz with 8 Gb ram (Macbook Pro 2012).
 
 @section design Design Intent
 Inherits from BOPropertyCalculatorInterface to allow interoperability.
 Implements update function from the PropertyCalculatorInterface to calculate properties at new P-T
 conditions.
 Implements getter functions defined in the BOPropertyCalculatorInterface
 to allow access to oil fluid properties required in an extended non-isothermal black oil formulation.
 
 @section examples Application Examples
 
 Example usage looks as following:
 
 @code
 BO_InputProperties bo_in = {{...}};
 double P = 3.e7;  // 300 bars
 double T = 120.0; // 120 degrees C
 BO_PropertyCalculatorInterface* calc_ptr = new BO_PropertyCalculator(bo_in);
 calc_ptr->Update(P, T);
 std::cout << calc_ptr->getUndersaturatedOilViscosity() << std::endl;
 P = 2.e7;
 T = 100.0;
 calc_ptr->Update(P,T);
 std::cout << calc_ptr->getUndersaturatedOilViscosity() << std::endl;
 @endcode
 
 */
class BO_PropertyCalculator : public BO_PropertyCalculatorInterface
{
public:
  explicit BO_PropertyCalculator( const BO_InputProperties& );
  virtual ~BO_PropertyCalculator();
  
  virtual void Update( double P, double T ) override final;
  
  double getUndersaturatedOilDensity() const override final {
    return d_.rho_o_SI;
  }
  
  double getUndersaturatedOilCompressibility() const override final {
    return d_.c_o_SI;
  }
  
  double getUndersaturatedOilViscosity() const override final {
    return d_.my_o_SI;
  }
  
  void Out() const override final;
  
private:
  BO_InputProperties   b_;
  BO_DynamicProperties d_;
  std::vector<CorrelationInterface*> correlations_;
};




/**
 @class ConstantReservoirWaterProperties UndersaturatedBOCorrelations.h
 
 @author L.J. Mosser
 
 @section motivation Motivation
 
 Basic implementation of ReservoirWaterPropertyCalculatorInterface for use with temperature and pressure
 dependent formulation of black oil properties. Provides constant water properties to
 classes requiring fluid properties for undersaturated oil water two phase simulations.
 All values are assumed to be constant and P-T independent.
 
 
 @section design Design Intent
 
 Inherits from ReservoirWaterPropertyCalculatorInterface to allow interoperability.
 Implements update function from the PropertyCalculatorInterface to calculate properties at new P-T
 conditions. Implements getter functions defined in the ReservoirWaterPropertyCalculatorInterface
 to allow access to water properties required in an extended non-isothermal black oil formulation.
 
 @section examples Application Examples
 
 Example usage looks as following:
 
 @code
 ReservoirWaterPropertyCalculatorInterface* calc_ptr = new ConstantReservoirWaterProperties(1.e-09, 1000.0, 0.0001);
 calc_ptr->Update(P, T);
 std::cout << calc_ptr->getReservoirWaterDensity() << std::endl;
 P = 2.e7;
 T = 100.0;
 calc_ptr->Update(P,T);
 std::cout << calc_ptr->getReservoirWaterDensity() << std::endl; //nothing will have changed becouse P-T independent
 @endcode
 */
class ConstantReservoirWaterProperties : public ReservoirWaterPropertyCalculatorInterface
{
public:
  ConstantReservoirWaterProperties( double cw, double rho_w, double visc_w )
  : cw_(cw), rho_w_(rho_w), visc_(visc_w) {}
  
  void Update( double, double ) override final {}
  
  double getReservoirWaterDensity() const override final {
    return rho_w_;
  }
  double getReservoirWaterCompressibility() const override final {
    return cw_;
  }
  double getReservoirWaterViscosity() const override final {
    return visc_;
  }
  
  void Out() const override final {}
  
private:
  const double visc_;
  const double rho_w_;
  const double cw_;
};

} //end csmp

#endif /* CSMP_UNDERSATURATED_BLACK_OIL_CORRELATIONS_H */
