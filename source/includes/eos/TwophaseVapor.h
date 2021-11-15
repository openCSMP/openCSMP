#ifndef TWOPHASEVAPOR_H
#define TWOPHASEVAPOR_H

#include "HaliteLiquidus.h"
#include "TwophaseLiquid.h"
#include "NaClSaturatedVapor.h"

/* Changelog
   
   October 16, 2014, Thomas Driesner
   - re-ported to CSMP++
   - cleaned up interface, eliminating functionaility that is currently not being used anymore
   - changed units in interface to [Pa] for pressure and [mass fraction NaCl] for composition
   - added basic doxygen documentation
   - changed coding style towards CSMP++ style guide recommendations
   - performed manual testing against CSP version - OK
 */
namespace csmp
{
  /// Composition of vapor phase on twophase vapor+liquid envelope
  class TwophaseVapor
  {

  public:
    // Constructor etc.
    TwophaseVapor(const double& externaltemperature, // [C]
                  const double& externalpressure,    // [Pa]
                  const double& ph2o);               // [Pa]
    ~TwophaseVapor();

    double MassFractionNaCl();                       // [XNaCl]

  private:

    TwophaseVapor();

    const double& temperature_;                      ///< internal reference to temperature [C] in flow code
    const double& pressure_;                         ///< internal reference to fluid pressure [P] in flow code
    const double& ph2o_;                             ///< internal reference to water boiling pressure (from Water object, which must have been constructed with the same temperature and fluid pressure variables as constructor arguments)

    double tcurrent_;                                ///< internal temperature variable [C]
    double pcurrent_;                                ///< internal pressure variable [Pa]
    double massfractionnacl_;                        ///< internal variable for composition 
    
    void     CheckState();
    double MoleFractionNaCl();

    HaliteLiquidus         liquidus;
    TwophaseLiquid         twophase_l;
    NaClSaturatedVapor     naclsatvap;
  };

  /**
     @class TwophaseVapor TwophaseVapor.h "eos/h2o_nacl/TwophaseVapor.h"

     @author Thomas Driesner, ETH Zuerich, thomas.driesner@erdw.ethz.ch

     @section motivation Motivation

     "TwophaseVapor" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, and computes the composition of vapor on the twophase liquid+vapor envelope as a function of temperature and pressure. It can be used standalone if only this information is needed. Otherwise, it is only used in the CSMP++ context for building accurate lookup tables for the H2O-NaCl system, which are significantly faster if the full system properties are needed. 

     Details of the underlying formulation can be found in

     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage
     Construct an instance of "TwophaseVapor" with the temperature ("externaltemperature", in C), fluid pressure ("externalpressure", in Pa) and water boiling pressure ("ph2o", in Pa, see attention remark below) as constructor arguments.

     "TwophaseVapor" has an internal mechanism to make sure that it always uses the current value of temperature, fluid pressure, and water boiling pressure as they exist in your code. Public member names should be self-explanatory, I hope.

     The NaCl concentration can be retrieved in units of mass fraction. If you want to convert these or need values for molar masses of NaCl and H2O, please use ConvertConcentrationUnitsNaCl.h to stay consistent.

     @attention The water boiling pressure as constructor argument requires that a Water object is being instantiated before twophase liquid, with the same temperature and fluid pressure as its constructor arguments. Also, before querying TwophaseVapor::MassFractionNaCl(), one has to make sure that ph2o has been updated. See code example. I apologize for this odd design but for performance reasons it made sense when created originally.

     @code
     double t(somevalue);          // temperature [C] in user's application
     double p(anothervalue);       // fluid pressure [in Pa] in user's application
     double ph2o(yetanothervalue); // water boling pressure [in Pa] in user's application
     ...
     Water          water(t,p);
     TwophaseVapor  twophase_vapor(t,p,ph2o); 
     ...
     // Updating t,p,ph2o:
     t    = my_new_t_value;
     p    = my_new_p_value;
     ph2o = water.SaturationPressure();
     // And retrieving the new value of mass fraction NaCl in TwophaseVapor:
     x = twophase_vapor.MassFractionNaCl(); // compute and assign to x the mass fraction of NaCl of TwophaseVapor at t and p
     ...
     @endcode

     @section testing Testing
     Tested manually against CSP version OK
  */

}//csmp
#endif





