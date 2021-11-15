#include <cmath>
#include <iostream>

#include "H2ONaClThermalEquilibrator.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "compareFloats.h"

using namespace std;

namespace csmp
{
  H2ONaClThermalEquilibrator::H2ONaClThermalEquilibrator(const double& external_mass_rock,  // [kg]
                                                         const double& external_cp_rock,    // [J/kg/K]
                                                         const double& external_rho_rock,   // [kg/m^3]
                                                         const double& external_porosity,   // [dimensionless]
                                                         const double& external_mass_fluid, // [kg]
                                                         const double& external_wt_current, // [wt% NaCl]
                                                         const double& external_t_previous, // [C]
                                                         const double& external_p_current,  // [bar]
                                                         const double& external_H_current,  // [J] NOT [J/kg]
                                                         const double& external_H_previous, // [J]
                                                         const bool&      fixed_external_t,
                                                         const double& t_fix,
                                                         const double& external_t_diff, 
                                                         const bool& verbose)
  : mass_rock(external_mass_rock),
    cp_rock(external_cp_rock),
    rho_rock(external_rho_rock),
    phi(external_porosity),
    mass_fluid(external_mass_fluid),
    wt_current(external_wt_current),
    t_previous(external_t_previous),
    p_current(external_p_current),
    H_current(external_H_current),
    H_previous(external_H_previous),
    t_fix(t_fix),
    t_diff(external_t_diff),
    fixed_t(fixed_external_t),
    verbose(verbose),
    equilibrated(false),
    fatal(false),
    convergence_criterion(1.0e-4),
    minimum_dt(20.0),
    t_eq(0.0),
    x_current_eq(1.0),
    h_fluid_eq(0.0),
    h_fluid_test(0.0),
  // change values of tmin and tmax only if you really know what you are doing
  //     these should be upper and lower limits of lookup table,
  //     and they MUST be initialized to the proper values before the 
  //     constructor body is reached
  // *** DEBUG info -> shall they be initialized by the "fluid" object?
    tmin(5.0), 
    tmax(990.0), 
    resid(0.0), 
    Hmax(0.0), 
    Hmin(0.0),
    mass_rock_eq(mass_rock-1.0e-10),
    cp_rock_eq(cp_rock-1.0e-10),
    mass_fluid_eq(mass_fluid-1.0e-10),
    wt_current_eq(wt_current-1.0e-10),
    t_previous_eq(t_previous-1.0e-10),
    p_current_eq(p_current-1.0e-10),
    H_current_eq(H_current-1.0e-10),
    H_previous_eq(H_previous-1.0e-10),
    H_test(0.0),
    tdummy(0.0),
    hdummy(0.0),
    t_previous_ini(0.0),
    mass_rock_ini(0.0),
    cp_rock_ini(0.0),
    rho_rock_ini(0.0),
    porosity_ini(0.0),
    mass_fluid_ini(0.0),
    wt_current_ini(0.0),
    x_current_ini(0.0),
    p_current_ini(0.0),
    H_current_ini(0.0),
    H_previous_ini(0.0),
    icrit(0),
    icrit_max(0),
  // fluid initialization: since p and x are given coordinates for iteration,
  // only the variable t and h are initialized via "dummy", which is needed
  // for member ComputeTotalEnthalpyAtTemperature( const double& t )
    fluid(tdummy, p_current_eq, x_current_eq, hdummy, cp_rock, rho_rock, phi, verbose),
    water_equilibrator(mass_rock, cp_rock, rho_rock, phi, mass_fluid,t_previous,
                       p_current, H_current, H_previous, fixed_t, t_fix, t_diff),
    csmp_error( ErrorHandler::Instance() )
  {
    // compute number of bisections (icrit_max) that would lead to reaching numerical precision
    //     following the formula (tmax-tmin) * 0.5^icrit_max <= epsilon;
    double dummy;
    dummy               = log(numeric_limits<double>::epsilon()/(tmax-tmin));
    dummy              /= log(0.5);
    icrit_max           = int(dummy)+1; // +1 for possible rounding problem
    icrit_max          += 1;            // +1 to ensure above formula works properly
    icrit_max          += 10;           // +10 = switch between weighted and true bisection
    // *** DEBUG ***
    // cerr << "icrit_max = " << icrit_max << endl;
    // char mychar;
    // cin >> mychar;
    equilibrated        = false;
  }

  H2ONaClThermalEquilibrator::~H2ONaClThermalEquilibrator()
  {
  }


  Fluidproperties H2ONaClThermalEquilibrator::Liquid()      
  { 
    Equilibrate(); 
    return liquidprops; 
  }


  Fluidproperties H2ONaClThermalEquilibrator::Vapor()
  { 
    Equilibrate(); 
    return vaporprops;  
  }


  Fluidproperties H2ONaClThermalEquilibrator::Bulk()
  { 
    Equilibrate();
    return bulkprops;
  }


  Fluidproperties H2ONaClThermalEquilibrator::Salt()
  { 
    Equilibrate(); 
    return saltprops;   
  }


  double H2ONaClThermalEquilibrator::Resid()
  { 
    return resid;
  }


  int H2ONaClThermalEquilibrator::n_iterations()
  { 
    return icrit;
  }


  bool H2ONaClThermalEquilibrator::Equilibrated()
  { 
    return equilibrated;
  }


  bool H2ONaClThermalEquilibrator::Fatal()
  { 
    return fatal;
  }


  void H2ONaClThermalEquilibrator::Equilibrate()
  {
    if(verbose) cerr << "***** H2ONaClThermalEquilibrator::Equilibrate() *****\n";
    // 0a. Set diagnostic bools to false
    equilibrated = false;
    fatal        = false;
    char           mychar;

    mass_rock_eq    = mass_rock;
    cp_rock_eq      = cp_rock;
    mass_fluid_eq   = mass_fluid;
    wt_current_eq   = wt_current;
    x_current_eq    = Weight2XNaCl(wt_current);
    p_current_eq    = p_current;
    H_current_eq    = H_current;
    H_previous_eq   = H_previous;

    if( essentiallyEqual( x_current_eq, 0.0, numeric_limits<double>::epsilon()) )
      {
        bulkprops   = water_equilibrator.Bulk();
        liquidprops = water_equilibrator.Liquid();
        vaporprops  = water_equilibrator.Vapor();
        saltprops.InitToZero();
        FluidPropertiesErrorCheck();
        return;
      }

    // 0b. If fixed_t, no equilibration is required (e.g. thermal DIRICH boundary condition);
    //     compute fluid properties and return.
    if(fixed_t == true)
      {
        t_previous_eq   = t_fix;
        tdummy          = t_fix;
	
        // *** DEBUG note: Check whether here we have to decide to go via HP or TP
        //                 This was at least an issue in the pure water limit,
        //                 maybe not here
        h_fluid_eq      = H_current_eq;
        h_fluid_eq     -= mass_rock_eq*rock.Enthalpy(tdummy);  //cp_rock_eq*t_fix;
        h_fluid_eq     /= mass_fluid_eq;
        hdummy          = h_fluid_eq;
        bulkprops       = fluid.BulkProperties();
        liquidprops     = fluid.LiquidProperties();
        vaporprops      = fluid.VaporProperties();
        saltprops       = fluid.SaltProperties();
        equilibrated    = true;
        return;
      }
    // 0c. Otherwise, update internal variables to have same values as external references
    else
      {
        t_previous_eq   = t_previous;

        mass_rock_ini   = mass_rock_eq;
        cp_rock_ini     = cp_rock_eq;
        mass_fluid_ini  = mass_fluid_eq;
        wt_current_ini  = wt_current_eq;
        x_current_ini   = x_current_eq;
        t_previous_ini  = t_previous_eq;
        p_current_ini   = p_current_eq;
        H_current_ini   = H_current_eq;
        H_previous_ini  = H_previous_eq;
      }  

    tmin                = 5.0;
    tmax                = 990.0;

    FindInitialValues();

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

        // cerr << "equilibrator starting bisection ...\n";
	
        // 2a. Compute new t_eq at mid-point between tmin and tmax, determine enthalpies,
        //     and compute residual
        if(icrit < 10
           && !essentiallyEqual( tmin, tmax, numeric_limits<double>::epsilon())
           && !essentiallyEqual( Hmin, Hmax, numeric_limits<double>::epsilon())
           )
          t_eq    = tmin + (H_current_eq-Hmin)/(Hmax-Hmin)*(tmax-tmin); //0.5*(tmin+tmax);
        else t_eq = 0.5*(tmin+tmax);
	
        H_test    = ComputeTotalEnthalpyAtTemperature( t_eq );
        resid     = H_current_eq - H_test;

        // 2b. Check for convergence
        //     It is absolutely crucial that the comparison is between the two fabs, since halite enthalpy
        //     can formally have negative values!!!
        if( definitelyLessThan( fabs(resid), 1.0, numeric_limits<double>::epsilon()) )
          {
            //d cerr << "converged with \n";
            //d cerr << "resid         = " << resid             << endl;
            //d cerr << "H_current_eq  = " << H_current_eq      << endl;
            //d cerr << "H_test        = " << H_test            << endl;
            //d cerr << "Hmin-H_curr   = " << Hmin-H_current_eq << endl;
            //d cerr << "Hmax-H_curr   = " << Hmax-H_current_eq << endl;
            //d cerr << "tmax-t_eq     = " << tmax-t_eq         << endl;
            //d cerr << "t_eq-tmin     = " << t_eq-tmin         << endl;
            //d cerr << "t_eq          = " << t_eq              << endl;
	    
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
          }
	
        // 2c. If not converged, check if tmin == tmax, that'd be definite non-convergence
        if( essentiallyEqual( tmin, tmax, numeric_limits<double>::epsilon()) )
          {
            double new_resid;
            new_resid = (Hmax - Hmin)/mass_fluid_ini;
            cerr << "*** new_resid                  = " << new_resid << "***\n";
            cerr << "*** (1.0e-1)*fabs(h_fluid_test) = " << (1.0e-1)*fabs(h_fluid_test) << "***\n";
            if(fabs(new_resid) < (1.0e-1)*fabs(h_fluid_test)) // quick fix
              {
                AssignAllPropertiesViaReport();
                FluidPropertiesErrorCheck();
                equilibrated = true;
                //d cerr << "survived FluidPropertiesErrorCheck()\n";
                //d PrintStatusToCerr();
                return;
                // since the loose criterion can cause very small negative saturations, these will be taken care of here
                // if(vaporprops.s < 0.0 && essentiallyEqual(vaporprops.s, 0.0, 5.0*numeric_limits<double>::epsilon()))
                //   { vaporprops.s = 0.0; return; }
                // else if(liquidprops.s < 0.0 && essentiallyEqual(liquidprops.s, 0.0, 5.0*numeric_limits<double>::epsilon()))
                //   { liquidprops.s = 0.0; return; }
                // else if(vaporprops.s < 0.0 && essentiallyEqual(vaporprops.s, 0.0, 5.0*numeric_limits<double>::epsilon()))
                //   { vaporprops.s = 0.0; return; }
                // else return;
              }
            cerr << "Upcoming non-convergence message for H2ONaClThermalEquilibrator::Equilibrate(), this \n";
            cerr << "was the status:\n";
            ErrorConditionsToScreen();
            csmp_error.notice( FATAL_ERROR,
                               //csmp_error.notice( WARNING,
                               "H2ONaClThermalEquilibrator::Equilibrate() -",
                               "tmin == tmax but no convergence based on enthalpy-criterion, see above data!\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                               );
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
          {
            tmin = t_eq;
            Hmin = ComputeTotalEnthalpyAtTemperature( tmin );
            ErrorCheckHmin(icrit);
            // tmax, Hmax remain as is
          }
        else  
          {
            tmax = t_eq;
            Hmax = ComputeTotalEnthalpyAtTemperature( tmax );
            ErrorCheckHmax(icrit);
            // tmin, Hmin remains as is
          }

        //d cerr << icrit << "\t" << Hmin << "\t" << H_current_eq << "\t" << Hmax << "\t" << resid << "\t" << bulkprops.h << endl;
        // 2e. It shouldn't happen but let's nevertheless check if tmin accidentally got larger than
        //     tmax, that's be fatal
        if( definitelyGreaterThan( tmin, tmax, numeric_limits<double>::epsilon()) )
          {	    
            csmp_error.notice( FATAL_ERROR,
                               //csmp_error.notice( WARNING,
                               "H2ONaClThermalEquilibrator::Equilibrate() -",
                               "tmin > tmax!\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                               );
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


  void H2ONaClThermalEquilibrator::AssignAllPropertiesViaReport(){
    // IMPORTANT: since bulkprops.h is now computed directly, all other 
    //            properties will NOT be up to date if convergence!
    //            Hence, the first line here was changed to 
    bulkprops   = fluid.BulkProperties();//ReportBulkProperties();
    liquidprops = fluid.ReportLiquidProperties();
    vaporprops  = fluid.ReportVaporProperties();
    saltprops   = fluid.ReportSaltProperties();
    FluidPropertiesErrorCheck();
    return;
  }

  double H2ONaClThermalEquilibrator::ComputeTotalEnthalpyAtTemperature( const double& t )
  {
    //    cerr << "H2ONaClThermalEquilibrator::ComputeTotalEnthalpyAtTemperature = " << t << endl;
    tdummy    =  t;
    hdummy    = (H_current_eq - mass_rock_eq*rock.Enthalpy(tdummy)); 
    hdummy   /=  mass_fluid_eq;
    //d cerr << "hdummy      = " << hdummy << endl;
    h_fluid_test = fluid.BulkEnthalpy();
    if(fluid.Fatal())
      {
        // cerr << "H2ONaClThermalEquilibrator::ComputeTotalEnthalpyAtTemperature( const double& t ) encountered fatal from Fluid objetc\n";
        PrintStatusToCerr();
        csmp_error.notice( FATAL_ERROR,
                           //                           csmp_error.notice( WARNING,
                           "H2ONaClThermalEquilibrator::ComputeTotalEnthalpyAtTemperature( const double& t ) encountered fatal from Fluid objetc\n",
                           "send above data to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                           ); 
      }
    return  mass_rock_eq*rock.Enthalpy(tdummy) + mass_fluid_eq*h_fluid_test;
  }



  void H2ONaClThermalEquilibrator::ErrorCheckHmin(const int& i)
  {
    if( definitelyLessThan( H_current_eq, Hmin, numeric_limits<double>::epsilon()) )
      {
        ErrorConditionsToScreen();
        liquidprops  = fluid.ReportLiquidProperties();
        vaporprops   = fluid.ReportVaporProperties();
        saltprops    = fluid.ReportSaltProperties();
        PrintStatusToCerr();
        std::stringstream istring, tstring;
        istring << i;
        std::string source( "H2ONaClThermalEquilibrator::Equilibrate()");
        std::string message("in iteration #");
        message += istring.str();
        message += ": H_current_eq is less than Hmin at ";
        tstring << tmin;
        message += tstring.str();
        message += "C, and, hence, we are out of range of validity.\nFIRST: check if you initialized your domain correctly to a valid temperatur.\nIf that check doesn't indicate any anomalies, report this incident to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch\n";
        csmp_error.notice( FATAL_ERROR, source, message );
        //csmp_error.notice( WARNING, source, message );
      }
    return;
  }

  void H2ONaClThermalEquilibrator::ErrorCheckHmax(const int& i)
  {
    if( definitelyGreaterThan( H_current_eq, Hmax, numeric_limits<double>::epsilon()) )
      {
        ErrorConditionsToScreen();
        liquidprops  = fluid.ReportLiquidProperties();
        vaporprops   = fluid.ReportVaporProperties();
        saltprops    = fluid.ReportSaltProperties();
        PrintStatusToCerr();
        std::stringstream istring, tstring;
        istring << i;
        std::string source( "H2ONaClThermalEquilibrator::Equilibrate()");
        std::string message("in iteration #");
        message += istring.str();
        message += ": H_current_eq is larger than Hmax at ";
        tstring << tmin;
        message += tstring.str();
        message += "C, and, hence, we are out of range of validity.\nFIRST: check if you initialized your domain correctly to a valid temperatur.\nIf that check doesn't indicate any anomalies, report this incident to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch\n";
        csmp_error.notice( FATAL_ERROR, source, message );
        //csmp_error.notice( WARNING, source, message );
      }
    return;
  }


  void H2ONaClThermalEquilibrator::FluidPropertiesErrorCheck()
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
        cerr << "H2ONaClThermalEquilibrator::FluidPropertiesErrorCheck() found a match\n";
        PrintStatusToCerr();
        ErrorConditionsToScreen();
        return;
      }
  }

  void H2ONaClThermalEquilibrator::PrintStatusToCerr()
  {
    cerr << "tdummy           = " << tdummy << endl;
    cerr << "p_current_eq     = " << p_current_eq << endl;
    cerr << "x_current_eq     = " << x_current_eq << endl;
    cerr << "t_previous_eq    = " << t_previous_eq << endl;
    cerr << "tmin             = " << tmin << endl;
    cerr << "tmax             = " << tmax << endl;
    cerr << "Hmin             = " << Hmin << endl;
    cerr << "Hmax             = " << Hmax << endl;
    cerr << "H_current_eq     = " << H_current_eq << endl;
    cerr << "H_previous_eq    = " << H_previous_eq << endl;
    cerr << "bulkprops.t      : " << bulkprops.t      << endl;
    cerr << "bulkprops.p      : " << bulkprops.p      << endl;
    cerr << "bulkprops.x      : " << bulkprops.x      << endl;
    cerr << "bulkprops.wt     : " << bulkprops.wt     << endl;
    cerr << "bulkprops.rho    : " << bulkprops.rho    << endl;
    cerr << "bulkprops.h      : " << bulkprops.h      << endl;
    cerr << "bulkprops.cp     : " << bulkprops.cp     << endl;
    // cerr << "bulkprops.dhdt   : " << bulkprops.dhdt   << endl;
    cerr << "bulkprops.beta   : " << bulkprops.beta   << endl;
    cerr << "bulkprops.s      : " << bulkprops.s      << endl;
    cerr << "bulkprops.mf     : " << bulkprops.mf     << endl;
    cerr << "bulkprops.mu     : " << bulkprops.mu     << endl;
    cerr << "bulkprops.state  : " << bulkprops.state  << endl << endl;
    cerr << "liquidprops.t    : " << liquidprops.t    << endl;
    cerr << "liquidprops.p    : " << liquidprops.p    << endl;
    cerr << "liquidprops.x    : " << liquidprops.x    << endl;
    cerr << "liquidprops.wt   : " << liquidprops.wt   << endl;
    cerr << "liquidprops.rho  : " << liquidprops.rho  << endl;
    cerr << "liquidprops.h    : " << liquidprops.h    << endl;
    cerr << "liquidprops.cp   : " << liquidprops.cp   << endl;
    // cerr << "liquidprops.dhdt : " << liquidprops.dhdt << endl;
    cerr << "liquidprops.beta : " << liquidprops.beta << endl;
    cerr << "liquidprops.s    : " << liquidprops.s    << endl;
    cerr << "liquidprops.mf   : " << liquidprops.mf   << endl;
    cerr << "liquidprops.mu   : " << liquidprops.mu   << endl;
    cerr << "liquidprops.state: " << liquidprops.state<< endl << endl;
    cerr << "vaporprops.t     : " << vaporprops.t     << endl;
    cerr << "vaporprops.p     : " << vaporprops.p     << endl;
    cerr << "vaporprops.x     : " << vaporprops.x     << endl;
    cerr << "vaporprops.wt    : " << vaporprops.wt    << endl;
    cerr << "vaporprops.rho   : " << vaporprops.rho   << endl;
    cerr << "vaporprops.h     : " << vaporprops.h     << endl;
    cerr << "vaporprops.cp    : " << vaporprops.cp    << endl;
    // cerr << "vaporprops.dhdt  : " << vaporprops.dhdt  << endl;
    cerr << "vaporprops.beta  : " << vaporprops.beta  << endl;
    cerr << "vaporprops.s     : " << vaporprops.s     << endl;
    cerr << "vaporprops.mf    : " << vaporprops.mf    << endl;
    cerr << "vaporprops.mu    : " << vaporprops.mu    << endl;
    cerr << "vaporprops.state : " << vaporprops.state << endl << endl;
    cerr << "saltprops.t      : " << saltprops.t     << endl;
    cerr << "saltprops.p      : " << saltprops.p     << endl;
    cerr << "saltprops.x      : " << saltprops.x     << endl;
    cerr << "saltprops.wt     : " << saltprops.wt    << endl;
    cerr << "saltprops.rho    : " << saltprops.rho   << endl;
    cerr << "saltprops.h      : " << saltprops.h     << endl;
    cerr << "saltprops.cp     : " << saltprops.cp    << endl;
    // cerr << "saltprops.dhdt   : " << saltprops.dhdt  << endl;
    cerr << "saltprops.beta   : " << saltprops.beta  << endl;
    cerr << "saltprops.s      : " << saltprops.s     << endl;
    cerr << "saltprops.mf     : " << saltprops.mf    << endl;
    cerr << "saltprops.mu     : " << saltprops.mu    << endl;
    cerr << "saltprops.state  : " << saltprops.state << endl;
  }

  void H2ONaClThermalEquilibrator::ErrorConditionsToScreen()
  {
    cerr.setf(ios::scientific);
    cerr.precision(20);
    cerr << "---H2ONaClThermalEquilibrator::ErrorConditionsToScreen()---:\n";
    cerr << "tmin          = " << tmin              << endl;
    cerr << "tmax          = " << tmax              << endl;
    cerr << "Hmin          = " << Hmin              << endl;
    cerr << "Hmax          = " << Hmax              << endl;
    cerr << "H_current_eq  = " << H_current_eq      << endl;
    cerr << "Hmin-H_curr   = " << Hmin-H_current_eq << endl;
    cerr << "Hmax-H_curr   = " << Hmax-H_current_eq << endl;
    cerr << "tmax-t_eq     = " << tmax-t_eq         << endl;
    cerr << "t_eq-tmin     = " << t_eq-tmin         << endl;
    cerr << "t_eq          = " << t_eq              << endl;
    cerr << "tdummy        = " << tdummy            << endl;
    cerr << "hdummy        = " << hdummy            << endl;
    cerr << "p_current_eq  = " << p_current_eq      << endl;
    cerr << "wt_current_eq = " << wt_current_eq     << endl;
    cerr << "x_current_eq  = " << x_current_eq      << endl;
    cerr << "bulk.state    = " << bulkprops.state   << endl;
    
    cerr << "---These started from the initial conditions---:\n";
    cerr << "m_rock      = " <<  mass_rock_ini << ";\n";
    cerr << "cp_rock     = " <<  cp_rock_ini << ";\n";
    cerr << "m_fluid     = " <<  mass_fluid_ini << ";\n";
    cerr << "wt          = " <<  wt_current_ini << ";\n";
    cerr << "x           = " <<  x_current_ini << ";\n";
    cerr << "t_previous  = " <<  t_previous_ini << ";\n";
    cerr << "p           = " <<  p_current_ini << ";\n";
    cerr << "H_current   = " <<  H_current_ini << ";\n";
    cerr << "H_previous  = " <<  H_previous_ini << ";\n";
    return;
  }



  Fluidproperties H2ONaClThermalEquilibrator::ReportLiquidProperties(const double& t, const double& p, const double& x, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.LiquidProperties();
  }

  Fluidproperties H2ONaClThermalEquilibrator::ReportVaporProperties(const double& t, const double& p, const double& x, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.VaporProperties();
  }

  Fluidproperties H2ONaClThermalEquilibrator::ReportBulkProperties(const double& t, const double& p, const double& x, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.BulkProperties();
  }

  Fluidproperties H2ONaClThermalEquilibrator::ReportSaltProperties(const double& t, const double& p, const double& x, const double& h)
  {
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.SaltProperties();
  }

  // JPW 12.4.2011
  void  H2ONaClThermalEquilibrator::ThreePhaseProperties(Fluidproperties& bulk_external,
                                                         Fluidproperties& liquid_external,
                                                         Fluidproperties& vapor_external,
                                                         Fluidproperties& salt_external)
  {
    Equilibrate();
    bulk_external   = bulkprops;
    liquid_external = liquidprops;
    vapor_external  = vaporprops;
    salt_external   = saltprops;
  }

  void H2ONaClThermalEquilibrator::FindInitialValues()
  {
    H_test = ComputeTotalEnthalpyAtTemperature( t_previous_ini );

    if(H_current_eq-H_test > 0.0e0)
      {
        //d  cerr << "H2ONaClThermalEquilibrator::FindInitialValues() ::FindInitialValues() H_current_eq-H_test > 0.0e0\n";
        // t_eq will be larger than t_previous_eq
        // current estimate is tmax because we used
        // too small a bulk heat capacity (i.e. rock only)
        tmax                   = t_previous_ini+(H_current_eq-H_test)/rock.MinimumHeatCapacity()/mass_rock_eq;
        tmax                  += minimum_dt;
        if(tmax > 990.) tmax   = 990.;
        tmin                   = t_previous_ini;
        tmin                  -= minimum_dt;
        if(tmin < 5.0) tmin    = 5.0;
      }
    
    else if(H_current_eq-H_test < 0.0e0)
      {
        // t_eq will be smaller than t_previous_eq
        // current estimate is tmax because we used 
        // too small a bulk heat capacity (i.e. rock only)
        //d  cerr << "H2ONaClThermalEquilibrator::FindInitialValues() ::FindInitialValues() H_current_eq-H_test < 0.0e0\n";
        tmin                   = t_previous_eq+(H_current_eq-H_test)/rock.MinimumHeatCapacity()/mass_rock_eq;
        tmin                  -= minimum_dt;
        if(tmin < 5.0e0) tmin  = 5.0e0;
        tmax                   = t_previous_ini;
        tmax                  += minimum_dt;
        if(tmax > 990.0) tmax  = 990.0;
      }

    else
      {
        // can only be H_current_eq == H_test
        tmin                  = t_previous_ini-minimum_dt;
        if(tmin < 5.0e0) tmin = 5.0e0;
        tmax                  = t_previous_ini+minimum_dt;
        if(tmax > 990.0) tmax = 990.0;

        //  char mychar;
        //d cerr << "H2ONaClThermalEquilibrator::FindInitialValues() ::FindInitialValues() H_current_eq-H_test == 0.0e0\n";
        //  cerr << "contact Thomas Driesner and provide him with the following data:\n";
        //  ErrorConditionsToScreen();
        //  cin >> mychar;
      }

    if( tmin > tmax )
      {
        ErrorConditionsToScreen();
        csmp_error.notice( FATAL_ERROR,
                           //csmp_error.notice( WARNING,
                           "H2ONaClThermalEquilibrator::FindInitialValues() -",
                           "tmin > tmax, send above data to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                           );
      }
    
    return;
  }
  
  
} // csmp

