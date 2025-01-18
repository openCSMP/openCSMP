#ifndef H2OFLUIDPROPERTIES_H
#define H2OFLUIDPROPERTIES_H

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

#include "States.h"
#include "Fluidproperties.h" 
#include "CriticalPointH2O.h"
#include "CriticalCurveLookup.h"
#include "H2OLookup.h"
/* #include "Rock.h" */



namespace csmp
{
  class H2OFluidProperties
  {
  public:
    H2OFluidProperties(   const double& externaltemperature_in_C,
                          const double& externalpressure_in_Pa,
                          const double& external_fluid_enthalpy_in_J_per_kg,
                          const double& external_cp_rock,
                          const double& external_rho_rock,
                          const double& external_phi);
    ~H2OFluidProperties();

    Fluidproperties       LiquidPropertiesFromTHP();
    Fluidproperties       VaporPropertiesFromTHP();
    Fluidproperties       BulkPropertiesFromTHP();

    Fluidproperties       LiquidPropertiesFromTP();
    Fluidproperties       VaporPropertiesFromTP();
    Fluidproperties       BulkPropertiesFromTP();

    Fluidproperties       ReportLiquidProperties();
    Fluidproperties       ReportVaporProperties();
    Fluidproperties       ReportBulkProperties();

    States                StateFromTHP();
    States                StateFromTP();

    // for debugging purposes only
    // int                   EqType();

    bool                  Equilibrated();

    void                  PrintProperties();
    void                  ReportAllProperties();
    void                  InitializeToBogus();

  private:

    const double&       temperature;
    const double&       pressure;
    const double&       enthalpy;
    const double&       cpr;
    const double&       rr;
    const double&       phi;

    double              tcurrent;
    double              pcurrent;
    double              hcurrent;
    double              ph2o;
    double              b;

    // for debugging purposes only
    //    int                   eqtype;

    bool                  equilibrated;
    bool                  below_pcrith2o;

    States                state;

    Fluidproperties       liq;
    Fluidproperties       vap;
    Fluidproperties       bulk;
    
    //    Rock                  rock;

    CriticalPointH2O      cp_h2o;
    CriticalCurveLookup   critcurve;
    H2OLookup             water;

    double              TwophaseCompressibility();

    void                  CheckEquilibrated();
    void                  ErrorCheck();
    void                  DumpStatus();

    void                  UpdatePropertiesFromTP();
    void                  UpdatePropertiesFromTHP();

    void                  PerformUpdate();

    void                  UpdatePropertiesVL();
    void                  UpdatePropertiesV_HighT();
    void                  UpdatePropertiesF_HighT();
    void                  UpdatePropertiesV_LowT();
    void                  UpdatePropertiesF_LowT( bool below_pcrith2o );
    void                  UpdatePropertiesToSaturatedLiquid();
    void                  UpdatePropertiesToSaturatedVapor();

    ErrorHandler&         csmp_error;
  };

  inline States           H2OFluidProperties::StateFromTHP()            { UpdatePropertiesFromTHP(); return state; } 
  inline Fluidproperties  H2OFluidProperties::LiquidPropertiesFromTHP() { UpdatePropertiesFromTHP(); return liq;   }
  inline Fluidproperties  H2OFluidProperties::VaporPropertiesFromTHP()  { UpdatePropertiesFromTHP(); return vap;   }
  inline Fluidproperties  H2OFluidProperties::BulkPropertiesFromTHP()   { UpdatePropertiesFromTHP(); return bulk;  }

  inline States           H2OFluidProperties::StateFromTP()            { UpdatePropertiesFromTP(); return state; } 
  inline Fluidproperties  H2OFluidProperties::LiquidPropertiesFromTP() { UpdatePropertiesFromTP(); return liq;   }
  inline Fluidproperties  H2OFluidProperties::VaporPropertiesFromTP()  { UpdatePropertiesFromTP(); return vap;   }
  inline Fluidproperties  H2OFluidProperties::BulkPropertiesFromTP()   { UpdatePropertiesFromTP(); return bulk;  }

  // caution with the next three. These report the current properties only without updating them before !!!
  inline Fluidproperties  H2OFluidProperties::ReportLiquidProperties() { return liq;   }
  inline Fluidproperties  H2OFluidProperties::ReportVaporProperties()  { return vap;   }
  inline Fluidproperties  H2OFluidProperties::ReportBulkProperties()   { return bulk;  }

  inline void             H2OFluidProperties::PrintProperties()        { UpdatePropertiesFromTHP(); DumpStatus(); }
  inline void             H2OFluidProperties::ReportAllProperties()    { DumpStatus(); }

  //  inline int              H2OFluidProperties::EqType()                 { return eqtype; }
  inline bool             H2OFluidProperties::Equilibrated()           { return equilibrated; }


  /**
     @class H2OFluidProperties H2OFluidProperties.h 

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch
     
     @section motivation Motivation
     H2OFluidProperties is an single interface to pure water single- and twophase properties in a porous medium via an equation of state. Information about pure water properties might also be obtained from H2OLookup but would require a reasonable knowledge of the thermodynamics involved (e.g., how to compute liquid-vapor saturations in the case of a twophase fluid) and would be incomplete because some properties vary with the exact properties of the rock matrix (e.g., the compressibility of a twophase liquid+vapor mixture is a function of the rock's heat capacity and porosity, see Grant&Sorey (1979), The Compressibility and Hydraulic Diffusivity of a Water-Steam Flow, Water Resources Research, 15, 3, 684-686).

     H2OFluidProperties eliminates the need of such a detailed implementation in terms of the interface of H2OLookup by the user. 
     
     @section usage Usage
     A H2OFluidProperties is instantiated with six constrcutor arguments. These are temperature in Celsius, pressure in bar, specific fluid enthalpy (in J/kg), the rock's heat capacity in J/kg/K, the rock's density in kg/m3, and the porosity (dimensionless). These are passed as const references to the respective existing objects in the user's code. H2OFluidProperties has an internal mechanism that makes sure that whenever it is queried (except for the Report... functions, see below) it checks and updates its internal placeholders to the current values of these variables, i.e., the user does NOT have to trigger an update. Although H2OFluidProperties can be used standalone, it is most often applied as a member of H2OPropertiesVisitor.

     The most important interfaces return objects of type Fluidproperties. These are:

     (1) LiquidPropertiesFromTHP(), VaporPropertiesFromTHP(), and BulkPropertiesFromTHP(). These determine the fluid state and properties primarily at the given temperature and pressure. If these conditions fall - within 5 x the available numerical precision of double on the given platform - very close to the boiling curve, sepcific enthalpy is used as an additional criterion to distinguish between liquid, vapor and twophase. In the case of twophase liquid+vapor, the actual enthalpy value is use to compute the saturations and mass fractions, respectively. If one of these member functions is invoked, all properties are computed. To avoid re-computation(*), properties of the other two may the also be extracted via the respective Report... functions described below. 

     (2) LiquidPropertiesFromTP(), VaporPropertiesFromTP(), and BulkPropertiesFromTP(). These determine the fluid state and properties at the given temperature and pressure. Where enthalpy is needed - again if t is within 5 x the available numerical precision of double on the given platform - the critical enthalpy is used, which always enforces a twophase state. This may be useful for initial and boundary conditions. The user should, however, be aware that this is an arbitrary value. If one of these member functions is invoked, all properties are computed. To avoid re-computation, properties of the other two may the also be extracted via the respective Report... functions described below.

     (3) ReportLiquidProperties(), ReportVaporProperties(), ReportBulkProperties(). These should be used with care. These functions bypass the internal automatic update of temperature etc. and will hence return the fluid properties as they were computed after the last internal update. This is likely to disagree with the current values of temperature, pressure etc. in the user's code. However, if one of the above ...FreomTHP() or ...FromTP() functons has been called, properties for the other two calls have already been computed and might equally well be extracted via the respective Report... functions (e.g.: if you called BulkPropertiesFromTHP(), you can safely extract vapor and liquid properties by immediately afterwards calling ReportVaporProperties() and ReportLiquidProperties().
    
     (*) This has changed from previous versions: the ...FromTHP and ...FromTP() functions used to have an internal check that avoided re-computation in case that t,p,h,cp_rock,rho_rock and phi didn't change. This has currently been disabled and may be re-activated after intense testing of this new version has been done.
     @section dependencies Dependencies
   
     
     @section issues Known issues and to do list
     - Rather than cp_rock etc. provide reference to a rock class or struct or function to allow for simple handling of t-dependent heat capacity etc.?
     - Have slightly improved the latter in TwoPhaseCompressibility but it needs re-check
     - maybe do two constructors with/without rock?
     - or re-check compressibility kieffer1977 vs. GrantSorey1979
     - implement a true FromHP formulation once requests emerge
     - avoid re-computation in FromTHP and FromTP functions

     @section testing Testing
     Testing is in progress.

     @section changes Changes
     H2OFluidProperties is an essentially completely changed and re-written version of what used to be FluidLookupWater.
     Thomas Driesner, Jan-Feb 2011

     @todo
     - check faster equilibration route for twophase conditions via enthalpy check
     - lots of assignments in UpdatePropertiesXXX not needed anymore?
  */

}//csmp

#endif
