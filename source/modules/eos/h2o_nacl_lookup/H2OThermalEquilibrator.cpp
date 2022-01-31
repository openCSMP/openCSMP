#include "H2OThermalEquilibrator.h"
#include "compareFloats.h"
#include <limits>

namespace csmp{

  H2OThermalEquilibrator::H2OThermalEquilibrator(const double& external_mass_rock,  // [kg]
                                                 const double& external_cp_rock,    // [J/kg/K]
                                                 const double& external_rho_rock,   // [kg/m^3]
                                                 const double& external_porosity,   // [dimensionless]
                                                 const double& external_mass_fluid, // [kg]
                                                 const double& external_t_previous, // [C]
                                                 const double& external_p_current,  // [bar]
                                                 const double& external_H_current,  // [J] NOT [J/kg]
                                                 const double& external_H_previous, // [J]
                                                 const bool&      fixed_external_t,
                                                 const double& t_fix,
                                                 const double& external_t_diff)
  : mass_rock(external_mass_rock),
    cp_rock(external_cp_rock),
    rho_rock(external_rho_rock),
    phi(external_porosity),
    mass_fluid(external_mass_fluid),
    t_previous(external_t_previous),
    p_current(external_p_current),
    H_current(external_H_current),
    H_previous(external_H_previous),
    t_fix(t_fix),
    t_diff(external_t_diff),
    fixed_t(fixed_external_t),
    equilibrated(false),
    fatal(false),
    convergence_criterion(1.0e-6),
    t_eq(0.0), 
    h_fluid_eq(0.0),
  // change values of tmin and tmax only if you really know what you are doing
  //     these should be upper and lower limits of lookup table,
  //     and they MUST be initialized to the proper values before the 
  //     constructor body is reached
  // *** DEBUG info -> shall they be initialized by the "fluid" object?
    tmin(5.0), \
    tmax(1000.0), 
    resid(0.0), 
    Hmax(0.0), 
    Hmin(0.0),
    icrit(0),
    icrit_max(0),
    mass_rock_eq(mass_rock-1.0e-10),
    cp_rock_eq(cp_rock-1.0e-10),
    mass_fluid_eq(mass_fluid-1.0e-10),
    t_previous_eq(t_previous-1.0e-10),
    p_current_eq(p_current-1.0e-10),
    H_current_eq(H_current-1.0e-10),
    H_previous_eq(H_previous-1.0e-10),
    H_test(0.0),
    tdummy(0.0),
    hdummy(0.0),
    cpr_min(0.0), 
    cpr_max(0.0), 
    t_min(0.0), 
    t_max(0.0), 
    cpr_t_dep(0.0),
    t_dependent_cpr(false),
    fluid(tdummy, p_current_eq, hdummy, cp_rock, rho_rock, phi),
    csmp_error( ErrorHandler::Instance() )
  {
    // compute number of bisections (icrit_max) that would lead to reaching numerical precision
    //     following the formula (tmax-tmin) * 0.5^icrit_max <= epsilon;
    double dummy;
    dummy               = log(numeric_limits<double>::epsilon()/(tmax-tmin));
    dummy              /= log(0.5);
    icrit_max           = int(dummy)+1; // +1 for possible rounding problem
    icrit_max          += 1;            // +1 to ensure above formula works properly

    icrit_max          += 10;           // 10 = switch between weighted and true bisection
    // *** DEBUG ***
    // cout << "icrit_max = " << icrit_max << endl;
    // char mychar;
    // cin >> mychar;
    equilibrated        = false;
  }


  H2OThermalEquilibrator::~H2OThermalEquilibrator()
  {
  }


  Fluidproperties H2OThermalEquilibrator::Liquid()
  {      
    Equilibrate(); 
    return liquidprops;
  }


  Fluidproperties H2OThermalEquilibrator::Vapor()
  {
    Equilibrate();
    return vaporprops;
  }


  Fluidproperties H2OThermalEquilibrator::Bulk()
  {
    Equilibrate();
    return bulkprops;
  }


  int H2OThermalEquilibrator::n_iterations()
  { 
    return icrit;
  }


  bool H2OThermalEquilibrator::Equilibrated()
  { 
    return equilibrated;
  }


  bool H2OThermalEquilibrator::Fatal()
  {        
    return fatal;
  }


  void H2OThermalEquilibrator::Equilibrate()
  {
    // 0a. Set diagnostic bools to false
    equilibrated = false;
    fatal        = false;

    // 0b. If fixed_t, no equilibration is required (e.g. thermal DIRICH boundary condition);
    //     compute fluid properties and return.
    if(fixed_t == true)
      {
        mass_rock_eq    = mass_rock;
        cp_rock_eq      = cp_rock;
        mass_fluid_eq   = mass_fluid;
        t_previous_eq   = t_fix;
        p_current_eq    = p_current;
        H_current_eq    = H_current;
        H_previous_eq   = H_previous;
        tdummy          = t_fix;
	
        // *** DEBUG note: Here we have to decide whether we go via HP or TP
        //                 dilemma - HP
        h_fluid_eq      = H_current_eq;
        h_fluid_eq     -= mass_rock_eq*rock.Enthalpy(tdummy)/*cp_rock_eq*t_fix*/;
        h_fluid_eq     /= mass_fluid_eq;
        hdummy          = h_fluid_eq;
        bulkprops       = fluid.BulkPropertiesFromTHP();
        liquidprops     = fluid.LiquidPropertiesFromTHP();
        vaporprops      = fluid.VaporPropertiesFromTHP();
        equilibrated    = true;
        return;
      }
    // 0c. Otherwise, update internal variables to have same values as external references
    else
      {
        mass_rock_eq    = mass_rock;
        cp_rock_eq      = cp_rock;
        mass_fluid_eq   = mass_fluid;
        t_previous_eq   = t_previous;
        p_current_eq    = p_current;
        H_current_eq    = H_current;
        H_previous_eq   = H_previous;
      }  

    // 1. Set initial values for bisection limits tmin and tmax at lower and upper limit
    //    of lookup table and initialize those variables used in convergence evaluation;
    //    resid         : difference between H_current_eq (total enthalpy (rock + fluid)
    //                    in control volume) and that at the temperature of interest; is 
    //                    the prime parameter to check for convergence. Convergence is 
    //                    achieved if resid < convergence_criterion*bulkprops.h, where
    //                    bulkprops.h is the specific fluid enthalpy at the t-p-h condition
    //                    of interest, it is typically set to 1e-4, lower values showed
    //                    non-convergence in earlier version; worth re-testing
    //    icrit         : a counter for the number of bisections; maximum is set to 50,
    //                    although that may be insufficient to bisect down to numerical 
    //                    precision
    tmin                = 5.0;
    tmax                = 1000.0;
    resid               = 1.0e0;
    icrit               = 0;
    equilibrated        = false; // was already set so, just making sure ...

    Hmin = ComputeTotalEnthalpyAtTemperature( tmin );
    ErrorCheckHmin(icrit);

    Hmax = ComputeTotalEnthalpyAtTemperature( tmax );
    ErrorCheckHmax(icrit);

    // 2. Bisection loop
    while(fabs(resid) > 1.0e-10)
      {

        //	cout << "equilibrator starting bisection ...\n";
        // 2a. Compute new t_eq at mid-point between tmin and tmax, determine enthalpies,
        //     and compute residual
        if(icrit < 10)
          t_eq    = tmin + (H_current_eq-Hmin)/(Hmax-Hmin)*(tmax-tmin); //0.5*(tmin+tmax);
        else t_eq = 0.5*(tmin+tmax);

        H_test    = ComputeTotalEnthalpyAtTemperature( t_eq );
        resid     = H_current_eq-H_test;
        //	cout << icrit << "\t" << resid << endl;
        // 2b. Check for convergence
        if( definitelyLessThan( fabs(resid), convergence_criterion*bulkprops.h, numeric_limits<double>::epsilon()) )
          {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
          }

        // 2c. If not converged, check if tmin = tmax, that'd be definite non-convergence
        if( essentiallyEqual( tmin, tmax, numeric_limits<double>::epsilon()) )
          {
            cerr << "---------------------------------------------------------------------------------\n";
            cerr << "Upcoming non-convergence message for H2OThermalEquilibrator::Equilibrate(), this \n";
            cerr << "was the status:\n";
            cerr << "tmin         = " << tmin              << endl;
            cerr << "tmax         = " << tmax              << endl;
            cerr << "Hmin         = " << Hmin              << endl;
            cerr << "Hmax         = " << Hmax              << endl;
            cerr << "H_current_eq = " << H_current_eq      << endl;
            cerr << "Hmin-H_curr  = " << Hmin-H_current_eq << endl;
            cerr << "Hmax-H_curr  = " << Hmax-H_current_eq << endl;
            cerr << "t_eq         = " << t_eq              << endl;
            cerr << "tdummy       = " << tdummy            << endl;
            cerr << "hdummy       = " << hdummy            << endl;
            cerr << "p_current_eq = " << p_current_eq      << endl;

            csmp_error.notice( FATAL_ERROR, 
                            "H2OThermalEquilibrator::Equilibrate() -",
                            "tmin == tmax but no convergence based on enthalpy-criterion, see above data!\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch");
            return;
          }
	
        // 2d. Define new tmin or tmax, depending on sign of H_test - H_current_eq
        // *** DEBUG info:
        //     Convergence has already been checked, so it seems unlikely that we have problems here with
        //     H_test and H_current_eq being close to each other even near epsilon.
        //     However, there remains a vague chance that tmin ot tmax may essentially be
        //     equal to t_eq, which may eventually show up either as an error at 2c or as
        //     icrit > icrit_max. I guess this remains to be seen by just running the code.
        if( definitelyLessThan( H_test, H_current_eq, numeric_limits<double>::epsilon()) ) 
          tmin = t_eq;
        else  
          tmax = t_eq;


        Hmin = ComputeTotalEnthalpyAtTemperature( tmin );
        ErrorCheckHmin(icrit);
	
        Hmax = ComputeTotalEnthalpyAtTemperature( tmax );
        ErrorCheckHmax(icrit);
	
        // 2e. It shouldn't happen but let's nevertheless check if tmin accidentally got larger than
        //     tmax, that's be fatal
        if( definitelyGreaterThan( tmin, tmax, numeric_limits<double>::epsilon()) )
          {	    
            csmp_error.notice( FATAL_ERROR, 
                            "H2OThermalEquilibrator::Equilibrate() -",
                            "tmin > tmax!\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch");
            return;
          }
	

        // 2f. If icrit > icrit_max -> non-convergence
        icrit++;
        if(icrit > icrit_max){
          fluid.InitializeToBogus();
          equilibrated = false; // to make sure the bogus values got to CSMP's range check
          fatal = true;
          return;
        }
      }   

    return;    
  }


  void H2OThermalEquilibrator::CheckForTminTmaxOutOfRange(){
    if( (tmin < 5.0e0) || (tmin > 1000.0e0) ||
        (tmax < 5.0e0) || (tmax > 1000.0e0) ){
      // in a newer version replace these by a backup solution that sets them to reasonable values
      cerr << "*** tmin = " << tmin << ", tmax = " << tmax << " ***\n";
      cerr << "H2OThermalEquilibrator::FindInitialValues() - \n";
      cerr << "tmin and/or tmax out of range, passing fatal to PropertiesVisitor!\n\n";
      // 	csmp_error.notice( FATAL_ERROR, 
      // 			"H2OThermalEquilibrator::FindInitialValues() -",
      // 			"tmin and/or tmax out of range!!!\n");
      fluid.InitializeToBogus();
      equilibrated = false;
      fatal = true;
      return;
    }
    return;
  }

  void H2OThermalEquilibrator::AssignAllPropertiesViaReport(){
    bulkprops   = fluid.ReportBulkProperties();
    liquidprops = fluid.ReportLiquidProperties();
    vaporprops  = fluid.ReportVaporProperties();
    FluidPropertiesErrorCheck();
    return;
  }

  // JPW Nov 2010
  void H2OThermalEquilibrator::TemperatureDependentHeatCapacityRock( double cpr_min_ext, double t_min_ext,
                                                                     double cpr_max_ext, double t_max_ext )
  {
    t_dependent_cpr = true;
    cpr_min = cpr_min_ext;
    cpr_max = cpr_max_ext;
    t_min = t_min_ext;
    t_max = t_max_ext;
  }

  // JPW Nov 2010
  double H2OThermalEquilibrator::TemperatureDependent_cpr( double T )
  {

    if (T < t_min)
      return cpr_min;
    else if (T < t_max)
      return (cpr_min*(t_max-T)+cpr_max*(T-t_min))/(t_max-t_min);
    else
      return cpr_max;

  }

  double H2OThermalEquilibrator::ComputeTotalEnthalpyAtTemperature( const double& t )
  {
    tdummy    =  t;
    hdummy    = (H_current_eq - mass_rock_eq*rock.Enthalpy(tdummy)/*cp_rock_eq*t*/); 
    hdummy   /=  mass_fluid_eq;
    bulkprops =  fluid.BulkPropertiesFromTHP();
    // cout << "ComputeTotalEnthalpyAtTemperature for " << t << "is " << mass_rock_eq*cp_rock_eq*t + mass_fluid_eq*bulkprops.h << endl;
    // cout << "with bulkprops.h = " << bulkprops.h << endl;
    return       mass_rock_eq*rock.Enthalpy(tdummy)/*cp_rock_eq*t*/ + mass_fluid_eq*bulkprops.h;
  }



  void H2OThermalEquilibrator::ErrorCheckHmin(const int& i)
  {
    if( definitelyLessThan( H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
      {
        cout << "tmin = " << tmin << endl;
        cout << "tmax = " << tmax << endl;
        cout << "Hmin = " << Hmin << endl;
        cout << "Hmax = " << Hmax << endl;
        cout << "H_current_eq = " << H_current_eq << endl;
        std::stringstream istring;
        istring << i;
        std::string source( "H2OThermalEquilibrator::Equilibrate()");
        std::string message("in iteration #");
        message += istring.str();
        message += ": H_current_eq is less than Hmin at 5 C, and, hence, we are out of range of validity.\nFIRST: check if you initialized your domain correctly to a valid temperatur.\nIf that check doesn't indicate any anomalies, report this incident to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch\n";
        csmp_error.notice( FATAL_ERROR, source, message );
      }
    return;
  }

  void H2OThermalEquilibrator::ErrorCheckHmax(const int& i)
  {
    if( definitelyGreaterThan( H_current_eq, Hmax, numeric_limits<double>::epsilon()) )
      {
        std::stringstream istring;
        istring << i;
        std::string source( "H2OThermalEquilibrator::Equilibrate()");
        std::string message("in iteration #");
        message += istring.str();
        message += ": H_current_eq is larger than Hmax at 1000 C, and, hence, we are out of range of validity.\nFIRST: check if you initialized your domain correctly to a valid temperatur.\nIf that check doesn't indicate any anomalies, report this incident to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch\n";
        csmp_error.notice( FATAL_ERROR, source, message );
      }
    return;
  }


  void H2OThermalEquilibrator::FluidPropertiesErrorCheck()
  { 
    if( 
       definitelyLessThan( bulkprops.t, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.p, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.x, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.wt, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.rho, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.cp, 0.0, numeric_limits<double>::epsilon())
       ||
       // definitelyLessThan( bulkprops.dhdt, 0.0, numeric_limits<double>::epsilon())
       // ||
       definitelyLessThan( bulkprops.beta, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.s, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.mf, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( bulkprops.mu, 0.0, numeric_limits<double>::epsilon())
       ||

       definitelyLessThan( liquidprops.t, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.p, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.x, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.wt, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.rho, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.h, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.cp, 0.0, numeric_limits<double>::epsilon())
       ||
       // definitelyLessThan( liquidprops.dhdt, 0.0, numeric_limits<double>::epsilon())
       // ||
       definitelyLessThan( liquidprops.beta, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.s, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.mf, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( liquidprops.mu, 0.0, numeric_limits<double>::epsilon())
       ||

       definitelyLessThan( vaporprops.t, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.p, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.x, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.wt, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.rho, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.h, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.cp, 0.0, numeric_limits<double>::epsilon())
       ||
       // definitelyLessThan( vaporprops.dhdt, 0.0, numeric_limits<double>::epsilon())
       // ||
       definitelyLessThan( vaporprops.beta, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.s, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.mf, 0.0, numeric_limits<double>::epsilon())
       ||
       definitelyLessThan( vaporprops.mu, 0.0, numeric_limits<double>::epsilon())
        )
      {
        cerr << "bulkprops.t"      << bulkprops.t      << endl;
        cerr << "bulkprops.p"      << bulkprops.p      << endl;
        cerr << "bulkprops.x"      << bulkprops.x      << endl;
        cerr << "bulkprops.wt"     << bulkprops.wt     << endl;
        cerr << "bulkprops.rho"    << bulkprops.rho    << endl;
        cerr << "bulkprops.h"      << bulkprops.h      << endl;
        cerr << "bulkprops.cp"     << bulkprops.cp     << endl;
        // cerr << "bulkprops.dhdt"   << bulkprops.dhdt   << endl;
        cerr << "bulkprops.beta"   << bulkprops.beta   << endl;
        cerr << "bulkprops.s"      << bulkprops.s      << endl;
        cerr << "bulkprops.mf"     << bulkprops.mf     << endl;
        cerr << "bulkprops.mu"     << bulkprops.mu     << endl;
        cerr << "liquidprops.t"    << liquidprops.t    << endl;
        cerr << "liquidprops.p"    << liquidprops.p    << endl;
        cerr << "liquidprops.x"    << liquidprops.x    << endl;
        cerr << "liquidprops.wt"   << liquidprops.wt   << endl;
        cerr << "liquidprops.rho"  << liquidprops.rho  << endl;
        cerr << "liquidprops.h"    << liquidprops.h    << endl;
        cerr << "liquidprops.cp"   << liquidprops.cp   << endl;
        // cerr << "liquidprops.dhdt" << liquidprops.dhdt << endl;
        cerr << "liquidprops.beta" << liquidprops.beta << endl;
        cerr << "liquidprops.s"    << liquidprops.s    << endl;
        cerr << "liquidprops.mf"   << liquidprops.mf   << endl;
        cerr << "liquidprops.mu"   << liquidprops.mu   << endl;
        cerr << "vaporprops.t"     << vaporprops.t     << endl;
        cerr << "vaporprops.p"     << vaporprops.p     << endl;
        cerr << "vaporprops.x"     << vaporprops.x     << endl;
        cerr << "vaporprops.wt"    << vaporprops.wt    << endl;
        cerr << "vaporprops.rho"   << vaporprops.rho   << endl;
        cerr << "vaporprops.h"     << vaporprops.h     << endl;
        cerr << "vaporprops.cp"    << vaporprops.cp    << endl;
        // cerr << "vaporprops.dhdt"  << vaporprops.dhdt  << endl;
        cerr << "vaporprops.beta"  << vaporprops.beta  << endl;
        cerr << "vaporprops.s"     << vaporprops.s     << endl;
        cerr << "vaporprops.mf"    << vaporprops.mf    << endl;
        cerr << "vaporprops.mu"    << vaporprops.mu    << endl;

      }
    return;
  }




  Fluidproperties H2OThermalEquilibrator::ReportLiquidProperties(const double& t, const double& p, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.LiquidPropertiesFromTHP();
  }

  Fluidproperties H2OThermalEquilibrator::ReportVaporProperties(const double& t, const double& p, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.VaporPropertiesFromTHP();
  }

  Fluidproperties H2OThermalEquilibrator::ReportBulkProperties(const double& t, const double& p, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.BulkPropertiesFromTHP();
  }



} // csmp

