#ifndef H2ONACLFLUIDPROPERTIES_H
#define H2ONACLFLUIDPROPERTIES_H

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

#include "States.h"
#include "Fluidproperties.h" 

#include "Rock.h"
#include "TriplePointNaCl.h"
#include "CriticalPointH2O.h"
#include "CriticalCurveLookup.h"
#include "NaClMeltingCurveLiquidLookup.h"
#include "NaClMeltingCurveHaliteLookup.h"
#include "HaliteLookup.h"
#include "HaliteLiquidusLookup.h"
#include "VLH_VaporLookup.h"
#include "VLH_LiquidLookup.h"
#include "VLH_HaliteLookup.h"
#include "NaClSaturatedVaporLookup.h"
#include "TwophaseLiquidLookup.h"
#include "TwophaseVaporLookup.h"
#include "H2ONaClLookup.h"
#include "H2OLookup.h"
#include "H2OFluidProperties.h"
#include "VH_HaliteLookup.h"
#include "LH_HaliteLookup.h"

namespace csmp
{
  class H2ONaClFluidProperties
  {
  public:
    H2ONaClFluidProperties(const   double64& externaltemperature_in_C,
                           const   double64& externalpressure_in_Pa,
                           const   double64& externalcomposition_in_mole_fraction,
                           const   double64& external_fluid_enthalpy,
                           const   double64& external_cp_rock,
                           const   double64& external_rho_rock,
                           const   double64& external_phi, 
                           const   bool&     verbose);
    ~H2ONaClFluidProperties();
    
    double64                       BulkEnthalpy();

    Fluidproperties                LiquidProperties();
    Fluidproperties                VaporProperties();
    Fluidproperties                SaltProperties();
    Fluidproperties                BulkProperties();

    Fluidproperties                ReportLiquidProperties();
    Fluidproperties                ReportVaporProperties();
    Fluidproperties                ReportSaltProperties();
    Fluidproperties                ReportBulkProperties();

    States                         State();
    // debug only
    //    int                            EqType();

    bool                           Equilibrated();
    bool                           Fatal();

    void                           PrintProperties();
    void                           ReportAllProperties();
    void                           InitializeToBogus();

  private:

    const double64&                temperature;
    const double64&                pressure;
    const double64&                composition;
    const double64&                enthalpy;
    const double64&                cpr;
    const double64&                rr;
    const double64&                phi;

    const double64                 safety_limit; 

    double64                       tcurrent; 
    double64                       pcurrent; 
    double64                       xcurrent; 
    double64                       hcurrent;
    double64                       dliqmfdp;
    double64                       dliqmfdt;
    double64                       dsaltmfdt;
    double64                       rl; 
    double64                       rv; 
    double64                       cpl; 
    double64                       cpv; 
    double64                       hl; 
    double64                       hv; 
    double64                       sl; 
    double64                       sv; 
    double64                       product; 
    double64                       b;
    double64                       liq_dsmfdt; 
    double64                       vap_dsmfdt;
    double64                       t_vlh_low; 
    double64                       t_vlh_high;
    double64                       h_vh; 
    double64                       h_lh; 
    double64                       h_vl; 
    double64                       liq_mf; 
    double64                       salt_mf; 
    double64                       vap_mf; 
    double64                       tsat;

    // debug only
    //    int                            eqtype;

    bool                           equilibrated;
    bool                           vl_boundary_encountered;
    bool                           low_x;
    bool                           fatal;
    const bool& verbose;
 
    States                         state;

    Fluidproperties                liq;
    Fluidproperties                vap;
    Fluidproperties                salt;
    Fluidproperties                bulk;
    
    Rock                           rock;
    TriplePointNaCl                tp_nacl;
    CriticalPointH2O               cp_h2o;
    CriticalCurveLookup            critcurve;
    NaClMeltingCurveLiquidLookup   naclmelt_l;
    NaClMeltingCurveHaliteLookup   naclmelt_h;
    HaliteLookup                   halite;
    HaliteLiquidusLookup           liquidus;
    VLH_VaporLookup                vlh_v;
    VLH_LiquidLookup               vlh_l;
    VLH_HaliteLookup               vlh_h;
    NaClSaturatedVaporLookup       naclsatvap;
    TwophaseLiquidLookup           twophase_l;
    TwophaseVaporLookup            twophase_v;
    H2ONaClLookup                  brine;
    H2OLookup                      water;
    H2OFluidProperties             fluid_water;
    VH_HaliteLookup                vh_halite;
    LH_HaliteLookup                lh_halite;

    double64                       TwophaseCompressibility();
    double64                       NewTwophaseCompressibility();

    void                           UpdateProperties();
    void                           UpdatePropertiesL();
    void                           UpdatePropertiesF();
    void                           UpdatePropertiesV();
    void                           UpdatePropertiesVH();
    void                           UpdatePropertiesLH();
    void                           UpdatePropertiesVL();
    void                           UpdatePropertiesAtVLH_For_P(double64 enthalpy);
    void                           UpdatePropertiesNaClMelt();
    void                           UpdatePropertiesHalite();
    void                           UpdatePropertiesNaClMeltingCurve_ForP(double64 enthalpy);

    void                           UpdateEnthalpyL();
    void                           UpdateEnthalpyF();
    void                           UpdateEnthalpyV();
    void                           UpdateEnthalpyVH();
    void                           UpdateEnthalpyLH();
    void                           UpdateEnthalpyVL();
    void                           UpdateEnthalpyAtVLH_For_P(double64 enthalpy);
    void                           UpdateEnthalpyNaClMelt();
    void                           UpdateEnthalpyHalite();
    void                           UpdateEnthalpyNaClMeltingCurve_ForP(double64 enthalpy);

    void                           PrintTwophaseCompressibilityParameters();
    void                           DumpStatus();
    void                           CheckEquilibrated();
    
    void Snapshot(const int& is);

    ErrorHandler&                  csmp_error;
  };

  inline States                    H2ONaClFluidProperties::State()            { UpdateProperties(); return state; } 
  //  inline int                       H2ONaClFluidProperties::EqType()           { UpdateProperties(); return eqtype; } 
  inline Fluidproperties           H2ONaClFluidProperties::LiquidProperties() { UpdateProperties(); return liq;   }
  inline Fluidproperties           H2ONaClFluidProperties::VaporProperties()  { UpdateProperties(); return vap;   }
  inline Fluidproperties           H2ONaClFluidProperties::SaltProperties()   { UpdateProperties(); return salt;  }
  inline Fluidproperties           H2ONaClFluidProperties::BulkProperties()   { UpdateProperties(); return bulk;  }
  // caution with the next four. These report the current properties only
  // without updating them before !!!
  inline Fluidproperties           H2ONaClFluidProperties::ReportLiquidProperties() { return liq;   }
  inline Fluidproperties           H2ONaClFluidProperties::ReportVaporProperties()  { return vap;   }
  inline Fluidproperties           H2ONaClFluidProperties::ReportSaltProperties()   { return salt;  }
  inline Fluidproperties           H2ONaClFluidProperties::ReportBulkProperties()   { return bulk;  }

  inline void                      H2ONaClFluidProperties::PrintProperties()  { UpdateProperties(); DumpStatus(); }
  inline void                      H2ONaClFluidProperties::ReportAllProperties()  { DumpStatus(); }
  inline bool                      H2ONaClFluidProperties::Equilibrated(){ return equilibrated; }
  inline bool                      H2ONaClFluidProperties::Fatal(){ return fatal; }


  /**
     @class H2ONaClFluidProperties H2ONaClFluidProperties.h 

     @author Thomas Driesner, ETH Zuerich
     @section contact Contact
     thomas.driesner@erdw.ethz.ch
     
     @section motivation Motivation
     H2ONaClFluidProperties is an single interface to H2O single-, two-, and three-phase properties in a porous medium via an equation of state. This information might also be obtained from the individual LookupTables for the H2O_naCl system but would require a detailed knowledge of the phase relations and thermodynamics involved and would be incomplete because some properties vary with the exact properties of the rock matrix (e.g., the compressibility of a twophase liquid+vapor mixture is a function of the rock's heat capacity and porosity, see Grant&Sorey (1979), The Compressibility and Hydraulic Diffusivity of a Water-Steam Flow, Water Resources Research, 15, 3, 684-686).

     H2ONaClFluidProperties eliminates the need of such a detailed implementation in terms of the interface of the various H2O-NaCl lookup tables (the list of which is apparent from the #include statements in this file) by the user. 
     
     @section usage Usage
     A H2ONaClFluidProperties is instantiated with seven constrcutor arguments. These are temperature in Celsius, pressure in bar, NaCl content in mole fraction, specific fluid enthalpy (in J/kg), the rock's heat capacity in J/kg/K, the rock's density in kg/m3, and the porosity (dimensionless). These are passed as const references to the respective existing objects in the user's code. H2ONaClFluidProperties has an internal mechanism that makes sure that whenever it is queried (except for the Report... functions, see below) it checks and updates its internal placeholders to the current values of these variables, i.e., the user does NOT have to trigger an update. Although H2ONaClFluidProperties can be used standalone, it is most often applied as a member of NaClH2OPropertiesVisitor.

     For accurate computation in the pure H2O limit, H2ONaClFluidProperties has an instance of H2OFluidProperties as a member.

     The most important interfaces return objects of type Fluidproperties. These are:

     (1) LiquidProperties(), VaporProperties(), SaltProperties(), and BulkProperties(). These determine the fluid state and properties primarily at the given temperature, pressure, and composition. Specific fluid enthalpy comes into play if the three-phase (vapor-liquid-halite) case or the pure H2O two-phase (vapor-liquid) or pure NaCl twophase (Halite+melt) case is encountered. Very low pressure cases along the respective sublimation curves are outside the range of validity as is the low-temperature cases involving ice or hydrohalite. Those conditions where enthalpy is required become active when temperature is identical to the temperature at the phase boundary within 5 x the available numerical precision of double64 on the given platform. 

     If one of these member functions is invoked, all properties are computed. To avoid re-computation(*), properties of the other two may the also be extracted via the respective Report... functions described below. This has changed from previous versions: the functions used to have an internal check that avoided re-computation in case that t,p,h,x,cp_rock,rho_rock and phi didn't change. This has currently been disabled and may be re-activated after intense testing of this new version has been done.

     (2) ReportLiquidProperties(), ReportVaporProperties(), ReportSaltProperties(), ReportBulkProperties(). These should be used with care. These functions bypass the internal automatic update of temperature etc. and will hence return the fluid properties as they were computed after the last internal update. This is likely to disagree with the current values of temperature, pressure etc. in the user's code. However, if one of the above LiquidProperties(), VaporProperties(), SaltProperties(), or BulkProperties() functions has been called, properties for the other two calls have already been computed and might equally well be extracted via the respective Report... functions (e.g.: if you called BulkProperties(), you can safely extract vapor, liquid, and salt properties by immediately afterwards calling ReportVaporProperties(), ReportLiquidProperties(), and ReportSaltProperties().
    
     @section dependencies Dependencies
        
     @section issues Known issues and to do list
     - Rather than cp_rock etc. provide reference to a rock class or struct or function to allow for simple handling of t-dependent heat capacity etc.?
     - implement a true FromHP formulation once requests emerge
     - avoid re-computation in FromTHP and FromTP functions

     @section testing Testing
     Testing is in progress.

     @section changes Changes
     H2ONaClFluidProperties is an essentially completely changed and re-written version of what used to be FluidLookup.
     Thomas Driesner, Jan-Feb 2011
  */


}//csmp


#endif
