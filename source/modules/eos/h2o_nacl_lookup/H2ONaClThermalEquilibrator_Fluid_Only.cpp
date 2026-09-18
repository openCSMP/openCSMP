// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "H2ONaClThermalEquilibrator_Fluid_Only.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "compareFloats.h"

using namespace std;

namespace csmp
{
//extern double                   global_time;

H2ONaClThermalEquilibrator_Fluid_Only::H2ONaClThermalEquilibrator_Fluid_Only(
        const double& external_mass_fluid, // [kg]
        const double& external_wt_current, // [wt% NaCl]
        const double& external_t_previous, // [C]
        const double& external_p_current,  // [bar]
        const double& external_H_current )

    : mass_fluid(external_mass_fluid),
      wt_current(external_wt_current),
      t_previous(external_t_previous),
      p_current(external_p_current),
      H_current(external_H_current),

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
      mass_fluid_eq(mass_fluid-1.0e-10),
      wt_current_eq(wt_current-1.0e-10),
      t_previous_eq(t_previous-1.0e-10),
      p_current_eq(p_current-1.0e-10),
      H_current_eq(H_current-1.0e-10),
      H_test(0.0),
      tdummy(0.0),
      hdummy(0.0),
      icrit(0),
      icrit_max(0),

      // fluid initialization: since p and x are given coordinates for iteration,
      // only the variable t and h are initialized via "dummy", which is needed
      // for member ComputeTotalEnthalpyAtTemperature( const double& t )
      fluid(tdummy, p_current_eq, x_current_eq, hdummy, 0., 0., 1., false),


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
    // char mychar;
    // cin >> mychar;
    equilibrated        = false;
}

H2ONaClThermalEquilibrator_Fluid_Only::~H2ONaClThermalEquilibrator_Fluid_Only()
{
}


Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::Liquid()
{
    Equilibrate();
    return liquidprops;
}


Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::Vapor()
{
    Equilibrate();
    return vaporprops;
}


Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::Bulk()
{
    Equilibrate();
    return bulkprops;
}


Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::Salt()
{
    Equilibrate();
    return saltprops;
}


double H2ONaClThermalEquilibrator_Fluid_Only::Resid()
{
    return resid;
}


int H2ONaClThermalEquilibrator_Fluid_Only::n_iterations()
{
    return icrit;
}


bool H2ONaClThermalEquilibrator_Fluid_Only::Equilibrated()
{
    return equilibrated;
}


bool H2ONaClThermalEquilibrator_Fluid_Only::Fatal()
{
    return fatal;
}


void H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate()
{
    ///*if(verbose) */cout << "***** H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate() *****\n";
    // 0a. Set diagnostic bools to false, update internal variables to have same values as external references
    equilibrated = false;
    fatal        = false;

    mass_fluid_eq   = mass_fluid;
    p_current_eq    = p_current;
    H_current_eq    = H_current;
    t_previous_eq   = t_previous;
    wt_current_eq   = wt_current;
    x_current_eq    = Weight2XNaCl(wt_current);

    tmin                = 5.0;
    tmax                = 990.0;

    resid               = 1.0e0;
    icrit               = 0;
    equilibrated        = false; // was already set so, just making sure ...

    InitialValues();

    if (equilibrated)
        return;


    // 2. Bisection loop
    while(fabs(resid) > 1.0e-10)
    {
        //cerr<<endl<<"Bisection loop ";
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
        if( definitelyLessThan( fabs(resid), 0.1, numeric_limits<double>::epsilon()) )
        {
            if(attempt_to_survive_fluid_properties_error) //BenoitLC add boolean
            {
                // Under some conditions (e.g. near the critical point), altough the solution converged, the returned fluid properties can be unphysical (e.g. negative viscosity)
                // Sometimes it is worth applying a small temperature correction in an attempt to recover physically meaningful fluid properties.
                // The correction is made in increments of 0.01C away from the converged solution, -0.01C if the temperature is on a downward trend or +0.01 in case of upward trend.
                // The correction cannot exceed 0.1C.
                // This is a coarse fix...

                properties_check_fail = false;
                temp_correction=0.;
                AssignAllPropertiesViaReport();
                if(properties_check_fail) cerr<<"Properties check FAIL\n";

                while(properties_check_fail && fabs(temp_correction)<1.)
                {
                    if(resid>=0.)
                    {
                        cerr<<"resid positive, we increase t_eq and H_test(t_eq) slightly\n";
                        t_eq+=0.05;
                        temp_correction+=0.05;
                        H_test    = ComputeTotalEnthalpyAtTemperature( t_eq );
                        cerr<<", temperature correction applied: "<<temp_correction<<endl;
                    }

                    else
                    {
                        cerr<<"resid negative, we decrease t_eq and H_test(t_eq) slightly\n";
                        t_eq-=0.05;
                        temp_correction-=0.05;
                        H_test    = ComputeTotalEnthalpyAtTemperature( t_eq );
                        cerr<<", temperature correction applied: "<<temp_correction<<endl;
                    }

                    AssignAllPropertiesViaReport();
                    if(properties_check_fail) cerr<<"Properties check FAIL\n";
                    else cerr<<"Properties check PASS!\n";
                }

                if ( properties_check_fail )
                    cerr<<"H2ONaClThermalEquilibrator_Fluid_Only, Converged with fluid properties error, could not fix it\n";
            }

            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }

        // 2c. If not converged, check if tmin == tmax, that'd be definite non-convergence
        if( essentiallyEqual( tmin, tmax, numeric_limits<double>::epsilon()) )
        {
            csmp_error.Note( WARNING,
                               "H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate() -",
                               "tmin == tmax but no convergence based on enthalpy-criterion, see above data!\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                                               );

                double new_resid;
                new_resid = (Hmax - Hmin)/mass_fluid;
                cout << "*** resid                  = " << resid << "***\n";
                cout << "*** new_resid                  = " << new_resid << "***\n";
                cout << "*** (3.0e-1)*fabs(h_fluid_test) = " << (3.0e-1)*fabs(h_fluid_test) << "***\n";
                if(fabs(resid) < (3.0e-1)*fabs(h_fluid_test) or fabs(new_resid) < (3.0e-1)*fabs(h_fluid_test)) // quick fix, "new_resid" added by Benoit
                {
                    AssignAllPropertiesViaReport();
                    FluidPropertiesErrorCheck();
                    equilibrated = true;
                    return;
                }

                cerr << "Upcoming non-convergence message for H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate(), this \n";
                cerr << "was the status:\n";
                ErrorConditionsToScreen();
                cout << "*** resid                  = " << resid << "***\n";
                cerr << "*** new_resid                  = " << new_resid << "***\n";
                cerr << "*** (3.0e-1)*fabs(h_fluid_test) = " << (3.0e-1)*fabs(h_fluid_test) << "***\n";

                return;

        }

        // 2d. Define new tmin or tmax, depending on sign of H_test - H_current_eq
        // *** DEBUG info:
        //     Convergence has already been checked, so it seems unlikely that we have problems here with
        //     H_test and H_current_eq being close to each other even near epsilon.
        //     However, there remains a vague chance that tmin ot tmax may essentially be
        //     equal to t_eq, which may eventually show up either as an error at 2c or as
        //     icrit > icrit_max. I guess this remains to be seen by just running the code.


        if (definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()))
        {
            tmin = t_eq;
            Hmin = H_test;
            ErrorCheckHmin(icrit);
            // tmax, Hmax remain as is
        }

        else
        {
            tmax = t_eq;
            Hmax = H_test;
            ErrorCheckHmax(icrit);
            // tmin, Hmin remains as is
        }

        // 2e. It shouldn't happen but let's nevertheless check if tmin accidentally got larger than
        //     tmax, that's be fatal
        if( definitelyGreaterThan( tmin, tmax, numeric_limits<double>::epsilon()) )
        {
            csmp_error.Note( FATAL_ERROR,
                               //csmp_error.Note( WARNING,
                               "H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate() -",
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


void H2ONaClThermalEquilibrator_Fluid_Only::AssignAllPropertiesViaReport(){
    // IMPORTANT: since bulkprops.h is now computed directly, all other
    //            properties will NOT be up to date if convergence!
    //            Hence, the first line here was changed to
    bulkprops   = fluid.BulkProperties();//ReportBulkProperties();
    liquidprops = fluid.ReportLiquidProperties();
    vaporprops  = fluid.ReportVaporProperties();
    saltprops   = fluid.ReportSaltProperties();
    FluidPropertiesErrorCheck();
    //return;
}

double H2ONaClThermalEquilibrator_Fluid_Only::ComputeTotalEnthalpyAtTemperature( const double& t )
{
    tdummy    = t;
    hdummy    = H_current_eq;
    hdummy   /= mass_fluid_eq;

    h_fluid_test = fluid.BulkEnthalpy();

    if(fluid.Fatal())
    {
        PrintStatusToCerr();
        csmp_error.Note( FATAL_ERROR,
                           "H2ONaClThermalEquilibrator_Fluid_Only::ComputeTotalEnthalpyAtTemperature( const double& t ) encountered fatal from Fluid objetc\n",
                           "send above data to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                           );
    }
    return  mass_fluid_eq * h_fluid_test;
}

void H2ONaClThermalEquilibrator_Fluid_Only::ErrorCheckHmin(const int& i)
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
        std::string source( "H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate()");
        std::string message("in iteration #");
        message += istring.str();
        message += ": H_current_eq is less than Hmin at ";
        tstring << tmin;
        message += tstring.str();
        message += "C, and, hence, we are out of range of validity.\nFIRST: check if you initialized your domain correctly to a valid temperature.\nIf that check doesn't indicate any anomalies, report this incident to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch\n";
        csmp_error.Note( FATAL_ERROR, source, message );
    }
    return;
}

void H2ONaClThermalEquilibrator_Fluid_Only::ErrorCheckHmax(const int& i)
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
        std::string source( "H2ONaClThermalEquilibrator_Fluid_Only::Equilibrate()");
        std::string message("in iteration #");
        message += istring.str();
        message += ": H_current_eq is larger than Hmax at ";
        tstring << tmin;
        message += tstring.str();
        message += "C, and, hence, we are out of range of validity.\nFIRST: check if you initialized your domain correctly to a valid temperature.\nIf that check doesn't indicate any anomalies, report this incident to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch\n";
        csmp_error.Note( FATAL_ERROR, source, message );
    }
    return;
}


void H2ONaClThermalEquilibrator_Fluid_Only::FluidPropertiesErrorCheck()
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
            //definitelyLessThan( bulkprops.beta, 0.0, numeric_limits<double>::epsilon())REMOVED BY BENOIT 2022
            //||
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
        cerr << "H2ONaClThermalEquilibrator_Fluid_Only::FluidPropertiesErrorCheck() found an error\n";
        PrintStatusToCerr();
        ErrorConditionsToScreen();
        //csmp_error.Note(FATAL_ERROR,"H2ONaClThermalEquilibrator_Fluid_Only","Fluid prop error");//BB
        properties_check_fail = true;
    }

    else properties_check_fail = false;
    return;
}

void H2ONaClThermalEquilibrator_Fluid_Only::PrintStatusToCerr()
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

void H2ONaClThermalEquilibrator_Fluid_Only::ErrorConditionsToScreen()
{
    //cerr.setf(ios::scientific);
    //cerr.precision(20);
    cerr << "---H2ONaClThermalEquilibrator_Fluid_Only::ErrorConditionsToScreen()---:\n";
    cerr << "tmin          = " << tmin              << endl;
    cerr << "tmax          = " << tmax              << endl;
    cerr << "Hmin          = " << Hmin              << endl;
    cerr << "Hmax          = " << Hmax              << endl;
    cerr << "H_current_eq  = " << H_current_eq      << endl;
    cerr << "H_test  = " << H_test      << endl;
    cerr << "resid  = " << resid      << endl;
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
    cerr << "m_fluid     = " <<  mass_fluid << ";\n";
    cerr << "wt          = " <<  wt_current << ";\n";
    cerr << "x           = " <<  Weight2XNaCl(wt_current) << ";\n";
    cerr << "t_previous  = " <<  t_previous << ";\n";
    cerr << "p           = " <<  p_current << ";\n";
    cerr << "H_current   = " <<  H_current << ";\n";
    return;
}

Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::ReportLiquidProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.LiquidProperties();
}

Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::ReportVaporProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.VaporProperties();
}

Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::ReportBulkProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.BulkProperties();
}

Fluidproperties H2ONaClThermalEquilibrator_Fluid_Only::ReportSaltProperties(const double& t, const double& p, const double& x, const double& h)
{
    h_fluid_eq    = h;
    t_eq          = t;
    p_current_eq  = p;
    x_current_eq  = x;
    t_previous_eq = -100.0; // to make sure Equilibrate doe not assume it has been called succesfully
    return fluid.SaltProperties();
}

// JPW 12.4.2011
bool  H2ONaClThermalEquilibrator_Fluid_Only::ThreePhaseProperties(Fluidproperties& bulk_external,
                                                                  Fluidproperties& liquid_external,
                                                                  Fluidproperties& vapor_external,
                                                                  Fluidproperties& salt_external)
{
    Equilibrate();
    bulk_external   = bulkprops;
    liquid_external = liquidprops;
    vapor_external  = vaporprops;
    salt_external   = saltprops;

    return equilibrated;
}

// added function PW May 2016
void H2ONaClThermalEquilibrator_Fluid_Only::InitialValues()
{
    // check for VLH first
    if (essentiallyEqual(p_current_eq, fluid.VLH_Pmax()))
    {
        t_eq = fluid.VLH_Tmax();
        H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
        resid = H_current_eq - H_test;

        if (definitelyLessThan(fabs(resid), convergence_criterion*H_current_eq, numeric_limits<double>::epsilon()))
        {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }
    }

    else if (definitelyLessThan(p_current_eq, fluid.VLH_Pmax()))
    {
        // VLH case #1
        t_eq = fluid.VLH_T_low(p_current_eq);
        H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
        resid = H_current_eq - H_test;

        if (definitelyLessThan(fabs(resid), convergence_criterion*H_current_eq, numeric_limits<double>::epsilon()))
        {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }

        // VLH case #1
        t_eq = fluid.VLH_T_high(p_current_eq);
        H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
        resid = H_current_eq - H_test;

        if (definitelyLessThan(fabs(resid), convergence_criterion*H_current_eq, numeric_limits<double>::epsilon()))
        {
            AssignAllPropertiesViaReport();
            equilibrated = true;
            return;
        }
    }

    t_eq = t_previous_eq;

    H_test = ComputeTotalEnthalpyAtTemperature(t_eq);
    resid = H_current_eq - H_test;

    if (definitelyLessThan(fabs(resid), convergence_criterion*H_current_eq, numeric_limits<double>::epsilon()))
    {
        AssignAllPropertiesViaReport();
        equilibrated = true;
        return;
    }

    else if (definitelyLessThan(H_test, H_current_eq, numeric_limits<double>::epsilon()))
    {
        tmin = t_eq;
        Hmin = H_test;
        double dT(1.);//(fabs(H_current_eq - H_test) / (bulkprops.cp * mass_fluid));

        tmax = tmin + dT;
        Hmax = ComputeTotalEnthalpyAtTemperature(tmax);

        // to be on the safe side:
        while (definitelyLessThan(Hmax, H_current_eq, numeric_limits<double>::epsilon()))
        {
            tmax += dT;
            Hmax = ComputeTotalEnthalpyAtTemperature(tmax);
        }
    }

    else
    {

        tmax = t_eq;
        Hmax = H_test;
        double dT(1.);//(fabs(H_test - H_current_eq) / (bulkprops.cp * mass_fluid));

        tmin = tmax - dT;
        Hmin = ComputeTotalEnthalpyAtTemperature(tmin);

        // to be on the safe side:
        while (definitelyLessThan(H_current_eq, Hmin, numeric_limits<double>::epsilon()))
        {
            tmin -= dT;
            Hmin = ComputeTotalEnthalpyAtTemperature(tmin);
            //cerr<<endl<<"tmin: "<<tmin<<", Hmin: "<<Hmin;
        }
    }

    ErrorCheckHmin(icrit);
    ErrorCheckHmax(icrit);

    if (tmin > tmax)
    {
        ErrorConditionsToScreen();
        csmp_error.Note(FATAL_ERROR,
                          "H2ONaClThermalEquilibrator_Fluid_Only::FindInitialValues() -",
                          "tmin > tmax, send above data to responsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch"
                          );
    }
    return;
}

//BenoitLC add
void H2ONaClThermalEquilibrator_Fluid_Only::AttemptToSurviveFluidPropertiesError(bool attempt_to_survive_fluid_properties_error_)
{
    attempt_to_survive_fluid_properties_error = attempt_to_survive_fluid_properties_error_;
}

} // csmp
