// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef EOS_CO2_H2O_NACL_SPYCHER_2004_H
#define EOS_CO2_H2O_NACL_SPYCHER_2004_H

#include <valarray>
#include <complex>
#include <cmath>
#include <iostream>
#include <vector>
#include"CSMP_definitions.h"

using namespace std;

namespace csmp {

    double  psiToPa( double pressureInPsi );

    double  paToPsi( double pressureInPa );

    double  paTobar( double pressureInPa );

    double  barTopa( double pressureInbar );

    double  degreeCToKelvin( double temperatureInC );

    double  KelvinTodegreeC( double temperatureInK );




class EOS_CO2H2ONaCl_Spycher04 {
  public:
    EOS_CO2H2ONaCl_Spycher04();
    ~EOS_CO2H2ONaCl_Spycher04();

// thermodynamic properties as a function of P, T and salt molality

    // compressed volume molar CO2
    double  V_Co2(double pressure,double temperature );

// molar fractions

    /// mole fraction (0..1) of CO2 in carbonic phase (the rest is H2O)
    double  y_Co2(double pressure,double temperature, double msalt );

    /// mole fraction (0..1) of H20 in carbonic phase
    double  y_H2o(double pressure,double temperature, double msalt );

    /// mole fraction (0..1) of CO2 in aqueous phase
    double  x_Co2(double pressure,double temperature, double msalt );

    /// mole fraction (0..1) of H20 in aqueous phase (the rest is CO2)
    double  x_H2o(double pressure,double temperature, double msalt );

    /// mole fraction of salt from salt molality in aqueous phase (may contain CO2 : molality)
    double  x_salt( double molalityCO2, double msalt );

    /// mass fraction (0..1) of salt in aqueous phase (weight percent)
    double  X_salt(double pressure,double temperature, double msalt );

    /// mass fraction of CO2 in carbonic phase (weight percent)
    double  Y_Co2(double pressure,double temperature, double msalt );

    /// mass fraction of H2O in carbonic phase (weight percent)
    double  Y_H2o(double pressure,double temperature, double msalt );

    /// mass fraction of CO2 in aqueous phase (weight percent)
    double  X_Co2(double pressure,double temperature, double msalt );

    /// mass fraction H2O in aqueous phase (weight percent)
    double  X_H2o(double pressure,double temperature, double msalt );

    /// mass fraction total dissolved salts in aqueous phase (weight percent)
    double  X_s( double pressure, double temperature, double msalt );


    /// molality of CO2 in aqueous phase
    double  m_Co2(double pressure,double temperature, double msalt );

// transport properties as a function of P, T and msalt

    /// molar volume of dissolved CO2 in brine
    double  Vdiss_Co2(double temperature );

    /// density brine
    double  Rho_brine(double pressure,double temperature, double msalt );

    /// viscosity brine
    double  mu_brine(double pressure,double temperature, double msalt );

    /// compressibility of brine
    double  C_brine(double pressure,double temperature, double msalt );

    /// density of aqueous phase with contained dissolved CO2
    double  Rho_AqueousPhase(double pressure,double temperature, double msalt );

    /// viscosity of aqueous phase with contain dissolved CO2
    double  mu_AqueousPhase(double pressure,double temperature, double msalt );

    /// compressibility of aqueous phase with dissolved CO2
    double  C_AqueousPhase(double pressure,double temperature, double msalt );

    /// density of the carbonic phase
    double  Rho_CarbonicPhase(double pressure,double temperature);

    /// viscosity of the carbonic phase
    double  mu_CarbonicPhase(double pressure,double temperature );

    /// compressibility of the carbonic phase
    double  C_CarbonicPhase(double pressure,double temperature );

    /// Z compressibility factor for carbonic phase
    double  Z_CarbonicPhase(double pressure,double temperature);

    /// molecular diffusion coefficient of CO2 in brine
    double  D_Co2(double pressure,double temperature, double msalt );

    /// equilibrium const KCO2
    double  KCo2(double pressure,double temperature );

    /// equilibrium const KH2O
    double  KH2o(double pressure,double temperature );

//  Reservoir properties as a function of T, P and msalt

    /// reservoir solution water-CO2 ratio
    double  Rs(double pressure, double temperature, double msalt );

    /// gas (CO2) formation volume factor, Bg
    double  Bg(double pressure, double temperature );

    /// water formation volume factor, Bw
    double  Bw(double pressure, double temperature ,double msalt);


//  unit conversions

    double  molalNaClToMassFracNaClInAqueousPhase( double mSalt); // input NaCl molality

    double  massFracNaClToMolalNaClInAqueousPhase( double massFracSalt); // input XNaCl

    double  massFracNaClToMolarFracNaClInAqueousPhase( double massFracSalt);// input XNaCl

    double  molalNaClToMolarFracNaClInAqueousPhase( double mSalt); // input NaCl molality

    double  ppmNaClToMolalNaClInAqueousPhase( double ppmSalt );

    double  molalNaClToPpmInAqueousPhase( double mSalt );

    double           CompressedVolumeCo2( double pressure, double temperature );

    double           FugacityCo2( double pressure,
                                    double temperature,
                                    double phaseVolumeCo2 );

    double           FugacityH2o( double pressure,
                                    double temperature,
                                    double phaseVolumeH2o );

    double           thermEquilConstCo2L( double temperature );

    double           thermEquilConstCo2G( double temperature );

    double           thermEquilConstH2o( double temperature );


    // need conversion (iterative approach) Equation A8 Spychler 2004
    double           activityCoefficientDrummond1981( double pressure,double temperature,
                                                        double mSalt );

    double           activityCoefficientDuanSun2003( double temperature,
                                                       double pressure,
                                                       double mSalt );
    /// ready to use
    double           activityCoefficientBattistelliEtal1997( double temperature,
                                                               double mSalt );

    /// needs conversion (iterative approach) Equation A8 Spychler 2004
    double           activityCoefficientRumpf1994( double pressure,double temperature,
                                                     double mSalt );

    double           calculateSpycherA( double pressure,
                                          double temperature,
                                          double kH2o,
                                          double phiH2o );

    double           calculateSpycherB( double pressure,
                                          double temperature,
                                          double phaseVolumeCo2,
                                          double  phiCo2,
                                          double kCo2L,
                                          double kCo2G,
                                          double activityCoeffCO2 );

    double           molarFracH2oCarbon( double spycherA,
                                           double spycherB,
                                           double mSalt );


    double           molarFracCO2Brine( double yH2o,
                                          double spycherB );

    double           molalNaClToMolarFracNaClInCo2SatAqueousPhase( double molalCo2,
                                                                     double mSalt );

    double           molarFracH2oBrine ( double xCo2,
                                           double xSalt );

    double           molarFracCo2Carbon( double yH2oBrine );


    double           massFracCo2InCarbonicPhase( double molarFracCo2InCarbonicPhase,
                                                   double molarFracH2oInCarbonicPhase );

    double           massFracH2oInCarbonicPhase( double molarFracCo2InCarbonicPhase,
                                                   double molarFracH2oInCarbonicPhase );

    double           massFracCo2inAqueousPhase( double molarFracCo2InAqueousPhase,
                                                  double molarFracH2oInAqueousPhase,
                                                  double molarFracNaclInAqueousPhase );

    double           massFracH2oInAqueousPhase( double molarFracCo2InAqueousPhase,
                                                  double molarFracH2oInAqueousPhase,
                                                  double molarFracNaclInAqueousPhase );

    double           massFracNaClInAqueousPhase( double molarFracCo2InAqueousPhase,
                                                   double molarFracH2oInAqueousPhase,
                                                   double molarFracNaclInAqueousPhase );

    double           molalCo2FromBrineMoleFractionCo2( double xCo2,
                                                         double mSalt );

    // I havent used them???
    double           molalCo2FromPureWaterMolarFracCo2( double MolarFracCo2 );

    double           molalCo2FromPureWaterMolalCo2( double molalCo2Pure,
                                                      double activityCoeff );

    double           molarFracCo2FromBrineMolalCo2( double molalCo2,
                                                      double mSalt);

    // caculation of a_co2 = a_mix as a function of T

    // activity of Co2 a_co2
    double  a_Co2(double temperature );

    // activity of mixture a_mix
    double  a_mix(double temperature );


    //  Auxilliary Functions
    complex<double>  complex_acos(const complex<double> & x);


    //  Transport Properties
    double           volumePartialMolarCo2( double temperature );

    double           densityBrine( double pressure,
                                     double temperature,
                                     double mSalt );


    double           viscosityBrine( double pressure,
                                       double temperature,
                                       double mSalt );

    double           compressibilityBrine( double densBrine,
                                             double pressure,
                                             double densBrineRef );

    double           densityAqueousPhase( double vPartialmolar,
                                            double densBrine,
                                            double xCo2 );

    double           densityCarbonicPhase( double phaseVolumeCo2 );

    double           compressibilityCarbonicPhase( double temperature,
                                                     double phaseVolumeCo2 );

    double           compressibilityCarbonicPhaseZ( double pressure,
                                                      double temperature,
                                                      double phaseVolumeCo2 );

    double           viscosityCarbonicPhase( double temperature,
                                               double phaseVolumeCo2 );

    double           molecularDiffCoeffCo2intoBrine( double temperature,
                                                       double muBrine );

    double           equilKH2o( double temperature,
                                  double pressure,
                                  double kH2o );

    double           equilKCo2( double temperature,
                                  double pressure,
                                  double kCo2G );

    //  Reservoir Properties
    double           solutionAqueousCarbonicRatio( double pressure, double temperature ,double msalt );

    double           GasFormationVolumeFactor( double pressure, double temperature );

    double           WaterFormationVolumeFactor( double pressure, double temperature ,double msalt);



    // plotting
    void plot_brine();

    void plot_AqueousPhase();

    void plot_CarbonicPhase();

    void plot_thermodynamics();

    const double R, molarMassH2o, molarMassCo2, molarMassNacl;

  private:
    // constants
    const double a_h2oco2, b_co2, b_h2o, b_mix,
                   stoichio,p0,vH2o, vCo2,
                   pSC,tSC;
};


} // end csmp

#endif /* _EOS_CO2H2ONaCl_Spycher04_SPYCHER_2004_H */
