#ifndef EOS_CO2_H2O_NACL_SPYCHER_2004_H
#define EOS_CO2_H2O_NACL_SPYCHER_2004_H

#include <valarray>
#include <complex>
#include <cmath>
#include <iostream>
#include <vector>
#include"CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

class EOS_CO2H2ONaCl_Spycher04 {

public:
    EOS_CO2H2ONaCl_Spycher04();
    ~EOS_CO2H2ONaCl_Spycher04();

// thermodynamic properties as a function of P, T and salt molality

    // compressed volume molar CO2
    double64  V_Co2(double64 pressure,double64 temperature );

// molar fractions

    /// mole fraction (0..1) of CO2 in carbonic phase (the rest is H2O)
    double64  y_Co2(double64 pressure,double64 temperature, double64 msalt );

    /// mole fraction (0..1) of H20 in carbonic phase
    double64  y_H2o(double64 pressure,double64 temperature, double64 msalt );

    /// mole fraction (0..1) of CO2 in aqueous phase
    double64  x_Co2(double64 pressure,double64 temperature, double64 msalt );

    /// mole fraction (0..1) of H20 in aqueous phase (the rest is CO2)
    double64  x_H2o(double64 pressure,double64 temperature, double64 msalt );

    /// mole fraction of salt from salt molality in aqueous phase (may contain CO2 : molality)
    double64  x_salt( double64 molalityCO2, double64 msalt );

    /// mass fraction (0..1) of salt in aqueous phase
    double64  X_salt(double64 pressure,double64 temperature, double64 msalt );

    /// mass fraction of CO2 in carbonic phase (weight percent)
    double64  Y_Co2(double64 pressure,double64 temperature, double64 msalt );

    /// mass fraction of H2O in carbonic phase (weight percent)
    double64  Y_H2o(double64 pressure,double64 temperature, double64 msalt );

    /// mass fraction of CO2 in aqueous phase (weight percent)
    double64  X_Co2(double64 pressure,double64 temperature, double64 msalt );

    /// mass fraction H2O in aqueous phase (weight percent)
    double64  X_H2o(double64 pressure,double64 temperature, double64 msalt );

    /// mass fraction total dissolved salts in aqueous phase (weight percent)
    double64  X_s( double64 pressure, double64 temperature, double64 msalt );


    /// molality of CO2 in aqueous phase
    double64  m_Co2(double64 pressure,double64 temperature, double64 msalt );

// transport properties as a function of P, T and msalt

    /// molar volume of dissolved CO2 in brine
    double64  Vdiss_Co2(double64 temperature );

    /// density brine
    double64  Rho_brine(double64 pressure,double64 temperature, double64 msalt );

    /// viscosity brine
    double64  mu_brine(double64 pressure,double64 temperature, double64 msalt );

    /// compressibility of brine
    double64  C_brine(double64 pressure,double64 temperature, double64 msalt );

    /// density of aqueous phase with contained dissolved CO2
    double64  Rho_AqueousPhase(double64 pressure,double64 temperature, double64 msalt );

    /// viscosity of aqueous phase with contain dissolved CO2
    double64  mu_AqueousPhase(double64 pressure,double64 temperature, double64 msalt );

    /// compressibility of aqueous phase with dissolved CO2
    double64  C_AqueousPhase(double64 pressure,double64 temperature, double64 msalt );

    /// density of the carbonic phase
    double64  Rho_CarbonicPhase(double64 pressure,double64 temperature);

    /// viscosity of the carbonic phase
    double64  mu_CarbonicPhase(double64 pressure,double64 temperature );

    /// compressibility of the carbonic phase
    double64  C_CarbonicPhase(double64 pressure,double64 temperature );

    /// Z compressbility factor for carbonic phase
    double64  Z_CarbonicPhase(double64 pressure,double64 temperature);

    /// molecular diffusion coefficient of CO2 in brine
    double64  D_Co2(double64 pressure,double64 temperature, double64 msalt );

    /// equilibrium const KCO2
    double64  KCo2(double64 pressure,double64 temperature );

    /// equilibrium const KH2O
    double64  KH2o(double64 pressure,double64 temperature );

//  Reservoir properties as a function of T, P and msalt

    /// reservoir solution water-CO2 ratio
    double64  Rs(double64 pressure, double64 temperature, double64 msalt );

    /// gas (CO2) formation volume factor, Bg
    double64  Bg(double64 pressure, double64 temperature );

    /// water formation volume factor, Bw
    double64  Bw(double64 pressure, double64 temperature ,double64 msalt);


//  unit conversions

    double64  molalNaClToMassFracNaClInAqueousPhase( double64 mSalt); // no CO2

    double64  massFracNaClToMolalNaClInAqueousPhase( double64 massFracSalt); // no CO2

    double64  massFracNaClToMolarFracNaClInAqueousPhase( double64 massFracSalt);// no CO2

    double64  molalNaClToMolarFracNaClInAqueousPhase( double64 mSalt);// no CO2

    double64  ppmNaClToMolalNaClInAqueousPhase( double64 ppmSalt );

    double64  molalNaClToPpmInAqueousPhase( double64 mSalt );

    double64  psiToPa( double64 pressureInPsi );

    double64  paToPsi( double64 pressureInPa );

    double64  paTobar( double64 pressureInPa );

    double64  barTopa( double64 pressureInbar );

    double64  degreeCToKelvin( double64 temperatureInC );

    double64  KelvinTodegreeC( double64 temperatureInK );


    double64 temp( double64 depth );

    double64 pres( double64 depth );


    void plot_brine();

    void plot_AqueousPhase();

    void plot_CarbonicPhase();

    void plot_thermodynamics();

    // plotting
    const double64 Tmax,Tmin,Pmax,Pmin,msaltmin,msaltmax,nTP,nmsalt;

    double64           CompressedVolumeCo2( double64 pressure, double64 temperature );

    double64           FugacityCo2( double64 pressure,
                                    double64 temperature,
                                    double64 phaseVolumeCo2 );

    double64           FugacityH2o( double64 pressure,
                                    double64 temperature,
                                    double64 phaseVolumeH2o );

    double64           thermEquilConstCo2L( double64 temperature );

    double64           thermEquilConstCo2G( double64 temperature );

    double64           thermEquilConstH2o( double64 temperature );


    // need conversion (iterative approach) Equation A8 Spychler 2004
    double64           activityCoefficientDrummond1981( double64 pressure,double64 temperature,
                                                        double64 mSalt );

    double64           activityCoefficientDuanSun2003( double64 temperature,
                                                       double64 pressure,
                                                       double64 mSalt );
    /// ready to use
    double64           activityCoefficientBattistelliEtal1997( double64 temperature,
                                                               double64 mSalt );

    /// needs conversion (iterative approach) Equation A8 Spychler 2004
    double64           activityCoefficientRumpf1994( double64 pressure,double64 temperature,
                                                     double64 mSalt );

    double64           calculateSpycherA( double64 pressure,
                                          double64 temperature,
                                          double64 kH2o,
                                          double64 phiH2o );

    double64           calculateSpycherB( double64 pressure,
                                          double64 temperature,
                                          double64 phaseVolumeCo2,
                                          double64  phiCo2,
                                          double64 kCo2L,
                                          double64 kCo2G,
                                          double64 activityCoeffCO2 );

    double64           molarFracH2oCarbon( double64 spycherA,
                                           double64 spycherB,
                                           double64 mSalt );


    double64           molarFracCO2Brine( double64 yH2o,
                                          double64 spycherB );

    double64           molalNaClToMolarFracNaClInCo2SatAqueousPhase( double64 molalCo2,
                                                                     double64 mSalt );

    double64           molarFracH2oBrine ( double64 xCo2,
                                           double64 xSalt );

    double64           molarFracCo2Carbon( double64 yH2oBrine );


    double64           massFracCo2InCarbonicPhase( double64 molarFracCo2InCarbonicPhase,
                                                   double64 molarFracH2oInCarbonicPhase );

    double64           massFracH2oInCarbonicPhase( double64 molarFracCo2InCarbonicPhase,
                                                   double64 molarFracH2oInCarbonicPhase );

    double64           massFracCo2inAqueousPhase( double64 molarFracCo2InAqueousPhase,
                                                  double64 molarFracH2oInAqueousPhase,
                                                  double64 molarFracNaclInAqueousPhase );

    double64           massFracH2oInAqueousPhase( double64 molarFracCo2InAqueousPhase,
                                                  double64 molarFracH2oInAqueousPhase,
                                                  double64 molarFracNaclInAqueousPhase );

    double64           massFracNaClInAqueousPhase( double64 molarFracCo2InAqueousPhase,
                                                   double64 molarFracH2oInAqueousPhase,
                                                   double64 molarFracNaclInAqueousPhase );

    double64           molalCo2FromBrineMoleFractionCo2( double64 xCo2,
                                                         double64 mSalt );

    // I havent used them???
    double64           molalCo2FromPureWaterMolarFracCo2( double64 MolarFracCo2 );

    double64           molalCo2FromPureWaterMolalCo2( double64 molalCo2Pure,
                                                      double64 activityCoeff );

    double64           molarFracCo2FromBrineMolalCo2( double64 molalCo2,
                                                      double64 mSalt);

    // caculation of a_co2 = a_mix as a function of T

    // activity of Co2 a_co2
    double64  a_Co2(double64 temperature );

    // activity of mixture a_mix
    double64  a_mix(double64 temperature );


    //  Auxilliary Functions
    complex<double64>  complex_acos(const complex<double64> & x);


    //  Transport Properties
    double64           volumePartialMolarCo2( double64 temperature );

    double64           densityBrine( double64 pressure,
                                     double64 temperature,
                                     double64 mSalt );


    double64           viscosityBrine( double64 pressure,
                                       double64 temperature,
                                       double64 mSalt );

    double64           compressibilityBrine( double64 densBrine,
                                             double64 pressure,
                                             double64 densBrineRef );

    double64           densityAqueousPhase( double64 vPartialmolar,
                                            double64 densBrine,
                                            double64 xCo2 );

    double64           densityCarbonicPhase( double64 phaseVolumeCo2 );

    double64           compressibilityCarbonicPhase( double64 temperature,
                                                     double64 phaseVolumeCo2 );

    double64           compressibilityCarbonicPhaseZ( double64 pressure,
                                                      double64 temperature,
                                                      double64 phaseVolumeCo2 );

    double64           viscosityCarbonicPhase( double64 temperature,
                                               double64 phaseVolumeCo2 );

    double64           molecularDiffCoeffCo2intoBrine( double64 temperature,
                                                       double64 muBrine );

    double64           equilKH2o( double64 temperature,
                                  double64 pressure,
                                  double64 kH2o );

    double64           equilKCo2( double64 temperature,
                                  double64 pressure,
                                  double64 kCo2G );

    //  Reservoir Properties
    double64           solutionAqueousCarbonicRatio( double64 pressure, double64 temperature ,double64 msalt );

    double64           GasFormationVolumeFactor( double64 pressure, double64 temperature );

    double64           WaterFormationVolumeFactor( double64 pressure, double64 temperature ,double64 msalt);

  private:
    // constants
    const double64 R,a_h2oco2, b_co2, b_h2o, b_mix,
                   stoichio,p0,vH2o, vCo2,
                   molarMassH2o, molarMassCo2, molarMassNacl,
                   pSC,tSC,therm_grad,pres_grad,temp_surface,pres_surface,max_depth;

};


} // end csmp

#endif /* _EOS_CO2H2ONaCl_Spycher04_SPYCHER_2004_H */
