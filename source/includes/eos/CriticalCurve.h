#ifndef CRITICALCURVE_H
#define CRITICALCURVE_H

#include "CSMP_definitions.h"
#include "CriticalPointH2O.h"

/*   Changelog
     March 27, 2013, Thomas Driesner: ported to CSMP++, changed pressure units in interface to Pa.
     September 9, 2014, Thomas Driesner: re-visited the port, cleaned up and adopted CSMP++ style rules
     September 25, 2014, Thomas Driesner: minor revision of documentation, basic testing against CSP version - OK
*/

namespace csmp
{

  /// Computation of temperature-pressure-compsition relations for the critical curve of the system H2O-NaCl, non-lookup
  class CriticalCurve
  {
  public:
    CriticalCurve(const double& externaltemperature);
    ~CriticalCurve();

    double         Pressure();
    double         MassFractionNaCl();
    double         MoleFractionNaCl();
    double         TfromP(const double& press_ ); 
    double         TfromMolefractionNaCl( const double& molefraction_ ); 

  private:
    CriticalCurve();

    CriticalPointH2O cp_h2o;

    const double&  temperature_;
    
    const double   c0_;   ///< Parameter c0  (= p_crit_h2o) of eqs 5, Driesner & Heinrich (2007)
    const double   c1_;   ///< Parameter c1   of eqs 5, Driesner & Heinrich (2007)
    const double   c2_;   ///< Parameter c2   of eqs 5, Driesner & Heinrich (2007)
    const double   c3_;   ///< Parameter c3   of eqs 5, Driesner & Heinrich (2007)
    const double   c4_;   ///< Parameter c4   of eqs 5, Driesner & Heinrich (2007)
    const double   c5_;   ///< Parameter c5   of eqs 5, Driesner & Heinrich (2007)
    const double   c6_;   ///< Parameter c6   of eqs 5, Driesner & Heinrich (2007)
    const double   c7_;   ///< Parameter c7   of eqs 5, Driesner & Heinrich (2007)
    const double   c8_;   ///< Parameter c8   of eqs 5, Driesner & Heinrich (2007)
    const double   c9_;   ///< Parameter c9   of eqs 5, Driesner & Heinrich (2007)
    const double   c10_;  ///< Parameter c10  of eqs 5, Driesner & Heinrich (2007)
    const double   c11_;  ///< Parameter c11  of eqs 5, Driesner & Heinrich (2007)
    double         c12_;  ///< Parameter c12  of eqs 5, Driesner & Heinrich (2007)
    double         c13_;  ///< Parameter c13  of eqs 5, Driesner & Heinrich (2007)
    const double   c14_;  ///< Parameter c14  of eqs 5, Driesner & Heinrich (2007)

    const double   c1A_;  ///< Parameter c1A  of eqs 5, Driesner & Heinrich (2007)
    const double   c2A_;  ///< Parameter c2A  of eqs 5, Driesner & Heinrich (2007)
    const double   c3A_;  ///< Parameter c3A  of eqs 5, Driesner & Heinrich (2007)
    const double   c4A_;  ///< Parameter c4A  of eqs 5, Driesner & Heinrich (2007)
    const double   c5A_;  ///< Parameter c5A  of eqs 5, Driesner & Heinrich (2007)
    const double   c6A_;  ///< Parameter c6A  of eqs 5, Driesner & Heinrich (2007)
    const double   c7A_;  ///< Parameter c7A  of eqs 5, Driesner & Heinrich (2007)
    const double   c8A_;  ///< Parameter c8A  of eqs 5, Driesner & Heinrich (2007)
    const double   c9A_;  ///< Parameter c9A  of eqs 5, Driesner & Heinrich (2007)
    const double   c10A_; ///< Parameter c10A of eqs 5, Driesner & Heinrich (2007)
    const double   c11A_; ///< Parameter c11A of eqs 5, Driesner & Heinrich (2007)

    const double   d1_;   ///< Parameter d1   of eqs 7, Driesner & Heinrich (2007)
    const double   d2_;   ///< Parameter d2   of eqs 7, Driesner & Heinrich (2007)
    const double   d3_;   ///< Parameter d3   of eqs 7, Driesner & Heinrich (2007)
    const double   d4_;   ///< Parameter d4   of eqs 7, Driesner & Heinrich (2007)
    const double   d5_;   ///< Parameter d5   of eqs 7, Driesner & Heinrich (2007)
    const double   d6_;   ///< Parameter d6   of eqs 7, Driesner & Heinrich (2007)
    const double   d7_;   ///< Parameter d7   of eqs 7, Driesner & Heinrich (2007)
    const double   d8_;   ///< Parameter d8   of eqs 7, Driesner & Heinrich (2007)
    const double   d9_;   ///< Parameter d9   of eqs 7, Driesner & Heinrich (2007)
    const double   d10_;  ///< Parameter d10  of eqs 7, Driesner & Heinrich (2007)
    const double   d11_;  ///< Parameter d11  of eqs 7, Driesner & Heinrich (2007)

    double         tcurrent_;        ///< temperature [C] for internal use
    double         molefraction_;    ///< mole fraction NaCl for internal use
    double         myt_;             ///< reduced temperature variable for internal use
    double         t_;               ///< another temperature for internal use
    double         delta_;           ///< convergence criterion
    double         dpdt_;            ///< tmeperature derivative of pressure (bar), internal use

    void             CheckState();
    double         Pressure(     const double& t_ );
    double         DPressureDT(  const double& t_ );
    double         Molefraction( const double& t_ );
    double         DMolefractionDT( const double& t_ );
  };


  /**
     @class CriticalCurve CriticalCurve.h "eos/h2o_nacl/CriticalCurve.h"                                              
                                                                                                                      
     @author Thomas Driesner, ETH Zuerich                                                                             
     @section contact Contact                                                                                         
     thomas.driesner@erdw.ethz.ch                                                                                     
    
     @section motivation Motivation                                                                                   
     Although "CriticalCurve" is legacy code from Thomas Driesner's developments for the H2O-NaCl system, it is REQUIRED in the CSMP++ context. There, it is essential to build accurate lookup tables for the H2O-NaCl system, which are then actually used by CSMP++. As long as we don't distribute the lookup tables per se for the various platforms, the "CriticalCurve" legacy will remain in CSMP++.                                                                                    
    
     "CriticalCurve" computes the temperature-pressure-molefraction [Celsius, Pa, mole fraction NaCl] coordinates of the critical curve of the system H2O-NaCl, according to the paper                                                     
                                                                                                                      
     Driesner T. and Heinrich C.A. (2007): The system H2O-NaCl. Part I: Correlation formulae for phase relations in temperature-pressure-molefraction space from 0 to 1000oC, 0 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta 71, 4880-4901.
                                                           
     @section usage Usage                                                                                             
     Construct an instance of "CriticalCurve" with the temperature (in C) as constructor variable. "CriticalCurve" has an internal mechanism to make sure that it always uses the current value of that temperature. Public member names should be self-explanatory, I hope.

     @code                                                                                                            
     double t; // temperature [C] in user's application                                                             
     CriticalCurve critcurve(t);                                                                                      
     @endcode                                                                                                         
                                                                                                                      
     @section dependencies Dependencies                                                                               
     requires "CriticalPointH2O.h"   

     @section issues Known issues                                                                                     
     Molefraction() is in mole fraction NaCl, to convert to other units, the functions in "ConvertConcentrationUnitsNaCl.h" may be useful                                                                                                   
                                                                                                                      
     @section testing Testing                                                                                         
     testing was done in the period before publication in 2007, prior to writing this documentation, no details are available anymore                                                                                                      
  */                                                                                                                  
       
}
#endif


