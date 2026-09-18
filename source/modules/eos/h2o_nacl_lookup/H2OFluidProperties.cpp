// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <limits>

#include "H2OFluidProperties.h"

#include "compareFloats.h"
#include "ConvertConcentrationUnitsNaCl.h"

using namespace std;

// August 2009:
// consolidated several functions and removed redundancies, eliminated need for
// special objects for saturated vapor and liquid by strict use of H2OLookup 

namespace csmp
{
  H2OFluidProperties::H2OFluidProperties(const double& externaltemperature_in_C,
                                         const double& externalpressure_in_Pa,
                                         const double& external_fluid_enthalpy_in_J_per_kg,
                                         const double& external_cp_rock,
                                         const double& external_rho_rock,
                                         const double& external_phi)
    : temperature(externaltemperature_in_C),
      pressure(externalpressure_in_Pa),
    enthalpy(external_fluid_enthalpy_in_J_per_kg),
    cpr(external_cp_rock),
    rr(external_rho_rock),
    phi(external_phi),
    tcurrent(temperature-1.0),
    pcurrent(pressure-1.0),
    hcurrent( 2086.e3 ),
    ph2o(0.0),
    b(0.0),
  //    eqtype(0),
    equilibrated(false),
    below_pcrith2o(false),
    with_rock_liquidus_solidus(false),
    tl(-1.0),
    ts(-1.0),
    state(none),
    critcurve(tcurrent), 
    water(),
    csmp_error( ErrorHandler::Instance() )
  {
    liq.InitToZero();
    vap.InitToZero();
    bulk.InitToZero();
  }
  
  H2OFluidProperties::~H2OFluidProperties()
  {
  }
  

  void H2OFluidProperties::UpdatePropertiesFromTP()
  {
    pcurrent = pressure;
    tcurrent = temperature;
    hcurrent = cp_h2o.Enthalpy();
    ph2o     = water.SaturationPressureFromT(tcurrent);
    
    PerformUpdate();
    return;
  }


  void H2OFluidProperties::UpdatePropertiesFromTHP()
  {
    pcurrent = pressure;
    tcurrent = temperature;
    hcurrent = enthalpy;
    ph2o     = water.SaturationPressureFromT(tcurrent);
    
    PerformUpdate();
    return;
  }


  void H2OFluidProperties::PerformUpdate()
  {
   
    equilibrated = false;

    // Determination of phase state, fully P-H based.
    // However, WHICH enthalpy value is used is determined via the interfaces UpdatePropertiesFromTP() or
    // UpdatePropertiesFromTHP() that call this member.
    // Numerical uncertainties should be catched in the respective UpdatePropertiesXXX functions, NOT here!

    if( !definitelyLessThan( tcurrent, cp_h2o.Temperature()) )
      {
        // Here we only check for what was intodruced as the "artificical phase boundary" in the H2O-NaCl
        // simulations to have consistent naming conventions. Since the critical curve in H2O-NaCl is almost
        // coincident with the critical isochore of pure H2O, calling everything "V" below this pressure at
        // any temperature above Tcrit(H2O) is a safe choice.
        //
        // Notice, however, that this may lead to artificial two-phase situations with both "V" (vapor) 
        // and "F" (single phase fluid, typically treated as liquid although this is potentially incorrect
        // at some T-P-X) in a given finite element. The transport code must be designed to detect this problem
        // and avoid assigning two-phase flow. In true two-phase vapor-liquid cases, the fluid state is "VL".
      
        if( !definitelyLessThan( pcurrent, critcurve.Pressure() ) )
          { /*eqtype = 0;*/ UpdatePropertiesF_HighT(); return; }

        else 
          { /*eqtype = 1;*/ UpdatePropertiesV_HighT(); return; }
      }	

    else
      {
        // Temperature is less than critical temperature of water.
        // First check if p is less than critical pressure of water; this is essential because
        // at p equal to or higher than the critical pressure, the very sensitive enthalpy checks
        // may be in error
        if(!definitelyLessThan( pcurrent, cp_h2o.Pressure() ) )
          { /*eqtype = 2;*/ below_pcrith2o = false; UpdatePropertiesF_LowT( below_pcrith2o ); return; }

        // What follows is an isobaric enthalpy check at p < pcrit
        // this should be most sensitive and much better than temperature-based
        else
          {
            if( definitelyLessThan( tcurrent, water.SaturationTemperatureFromP( pcurrent ), 5.0*numeric_limits<double>::epsilon() ) )
              { /*eqtype = 3;*/ below_pcrith2o = true; UpdatePropertiesF_LowT( below_pcrith2o ); return; }
	
            else if( definitelyGreaterThan( tcurrent, water.SaturationTemperatureFromP( pcurrent ), 5.0*numeric_limits<double>::epsilon() ) )
              { /*eqtype = 4;*/ UpdatePropertiesV_LowT(); return; }
	
            // Now we know that we are within 5 epsilons of tsat, i.e., we can do a simple check for VL
            else if( definitelyLessThan( hcurrent, water.LiquidEnthalpyFromP( pcurrent ) ) )
              { /*eqtype = 5;*/ UpdatePropertiesToSaturatedLiquid(); return; }
	
            else if( definitelyGreaterThan( hcurrent, water.VaporEnthalpyFromP( pcurrent ) ) )
              { /*eqtype = 6;*/ UpdatePropertiesToSaturatedVapor(); return; }
	
            else
              { /*eqtype = 7;*/ UpdatePropertiesVL(); return; }
          }
      }
  }



  void H2OFluidProperties::UpdatePropertiesF_HighT()
  {
    // looks ok, not yet tested, 27-Jan-2011, TD
    vap.InitToZero();
    state      = F;
    
	// Changing back to original version -- PW March 2017
	// requires in-depth code comparison

	//bulk = water.UpdatePropertiesSinglePhase(tcurrent, pcurrent);
	bulk.rho = water.Density(tcurrent, pcurrent);
	bulk.h = water.Enthalpy(tcurrent, pcurrent);
	bulk.cp = water.HeatCapacity(tcurrent, pcurrent);
	bulk.mu = water.Viscosity(tcurrent, pcurrent);
	bulk.beta = water.Compressibility(tcurrent, pcurrent);

    bulk.state = F;
    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = 0.0e0;
    bulk.smf   = 0.0e0;
    bulk.wt    = 0.0;
    
    //// bulk.dhdt  = water.HeatCapacity(tcurrent,pcurrent);
    
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    liq        = bulk;
    //ErrorCheck();
    CheckEquilibrated();
    // cout << "H2OFluidProperties::UpdatePropertiesF_HighT() - \n";
    // DumpStatus();
    return;
  }

  void H2OFluidProperties::UpdatePropertiesF_LowT( bool below_pcrith2o )
  {
    // looks ok, not yet tested, 27-Jan-2011, TD
    vap.InitToZero();
    state      = F;

	// Changing back to original version -- PW March 2017
	// requires in-depth code comparison

	//bulk = water.UpdatePropertiesSinglePhase(tcurrent,pcurrent);
	bulk.rho = water.Density(tcurrent, pcurrent);
	bulk.h = water.Enthalpy(tcurrent, pcurrent);
	bulk.cp = water.HeatCapacity(tcurrent, pcurrent);
	bulk.mu = water.Viscosity(tcurrent, pcurrent);
	bulk.beta = water.Compressibility(tcurrent, pcurrent);

    bulk.state = F;
    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = 0.0e0;
    bulk.smf   = 0.0e0;
    bulk.wt    = 0.0;

	// bulk.dhdt  = bulk.cp;

    //************************************
    // enforce error handling!!!!!!!!!!!!!
    //************************************

    // The following errorcheck makes only sense if p < pcrith2o
    if( below_pcrith2o )
      {
        if( definitelyLessThan(    bulk.rho,  water.LiquidDensityFromP(         pcurrent ))
            ||
            definitelyGreaterThan( bulk.h,    water.LiquidEnthalpyFromP(        pcurrent ))
            ||
            definitelyLessThan(    bulk.mu,   water.LiquidViscosityFromP(       pcurrent ))
            // ||
            // definitelyGreaterThan( bulk.beta, water.LiquidCompressibilityFromP( pcurrent ))
            // did not include compressibility  as I am not sure if there is a truly unique criterion!
            )
          {
            cerr << "---------------------------------------------------------------------------\n";
			cerr << "tcurrent: " << tcurrent << "; pcurrent: " << pcurrent << endl;
			cerr << "H2OFluidProperties::UpdatePropertiesF_LowT() - encountered error condition:\n";
            cerr << "bulk.rho,   water.LiquidDensityFromP         : " << bulk.rho  << "\t" <<  water.LiquidDensityFromP(pcurrent)         << endl;
            cerr << "bulk.h,     water.LiquidEnthalpyFromP        : " << bulk.h    << "\t" <<  water.LiquidEnthalpyFromP(pcurrent)        << endl;
            cerr << "bulk.mu,    water.LiquidViscosityFromP       : " << bulk.mu   << "\t" <<  water.LiquidViscosityFromP(pcurrent)       << endl;
            cerr << "bulk.beta, water.LiquidCompressibilityFromP : " << bulk.beta << "\t" <<  water.LiquidCompressibilityFromP(pcurrent) << endl;
            //	    DumpStatus();
            cerr << "---------------------------------------------------------------------------\n\n";
            csmp_error.Note( FATAL_ERROR,
                             "H2OFluidProperties::UpdatePropertiesF_LowT() - encountered error condition:",
                             "see output above, terminating ...\n");

            // bulk.rho   = water.LiquidDensityFromP( pcurrent );
            // bulk.h     = water.LiquidEnthalpyFromP(pcurrent);
            // bulk.cp    = water.LiquidHeatCapacityFromP(pcurrent);
            // bulk.dhdt  = bulk.cp;
            // bulk.mu    = water.LiquidViscosityFromP(pcurrent);
            // bulk.beta  = water.LiquidCompressibilityFromP(pcurrent);
          }
      }
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    liq        = bulk;
    //ErrorCheck();
    CheckEquilibrated();
    // cout << "H2OFluidProperties::UpdatePropertiesF_LowT() - \n";
    // DumpStatus();
    return;
  }
  
  
  void H2OFluidProperties::UpdatePropertiesV_LowT()
  {
    // looks ok, not yet tested, 27-Jan-2011, TD
    liq.InitToZero();
    state      = V;
    
    //pcurrent   = pressure;

	// Changing back to original version -- PW March 2017
	// requires in-depth code comparison

    //bulk = water.UpdatePropertiesSinglePhase(tcurrent,pcurrent);
	bulk.rho = water.Density(tcurrent, pcurrent);
	bulk.h = water.Enthalpy(tcurrent, pcurrent);
	bulk.cp = water.HeatCapacity(tcurrent, pcurrent);
	bulk.mu = water.Viscosity(tcurrent, pcurrent);
	bulk.beta = water.Compressibility(tcurrent, pcurrent);

    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = 0.0e0;
    bulk.smf   = 0.0e0;
    bulk.wt    = 0.0e0;
    
    // bulk.dhdt  = bulk.cp;

    // here, however, the check should always make sense since it is ONLY called if p<pcrith2o anyhow
    if( definitelyGreaterThan( bulk.rho, water.VaporDensityFromP( pcurrent ) )
        || 
        definitelyLessThan(    bulk.h  , water.VaporEnthalpyFromP( pcurrent )) 
        // other props might be ambigous ...
        // how well this works, I don't know yet
        )
      {
        bulk.rho   = water.VaporDensityFromP(pcurrent);
        bulk.beta  = water.VaporCompressibilityFromP(pcurrent);
        bulk.h     = water.VaporEnthalpyFromP(pcurrent);
        bulk.cp    = water.VaporHeatCapacityFromP(pcurrent);
        // bulk.dhdt  = bulk.cp;
        bulk.mu    = water.VaporViscosityFromP(pcurrent);
      }
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    bulk.state = V;
    vap        = bulk;
    // ErrorCheck();
    CheckEquilibrated();
    //cout << "H2OFluidProperties::UpdatePropertiesV_LowT() - \n";
    // DumpStatus();
    return;
  }

  void H2OFluidProperties::UpdatePropertiesV_HighT()
  {
    // looks ok, not yet tested, 27-Jan-2011, TD
    liq.InitToZero();
    state      = V;

	// Changing back to original version -- PW March 2017
	// requires in-depth code comparison

    //bulk = water.UpdatePropertiesSinglePhase(tcurrent, pcurrent);
	bulk.rho = water.Density(tcurrent, pcurrent);
	bulk.h = water.Enthalpy(tcurrent, pcurrent);
	bulk.cp = water.HeatCapacity(tcurrent, pcurrent);
	bulk.mu = water.Viscosity(tcurrent, pcurrent);
	bulk.beta = water.Compressibility(tcurrent, pcurrent);

    // pcurrent   = pressure;
    // don't know why  that line was there
    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = 0.0e0;
    bulk.smf   = XNaCl2Massfraction(bulk.x);
    bulk.wt    = 100.0*bulk.smf;
    
    // bulk.dhdt  = bulk.cp;
    
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    bulk.state = V;
    vap        = bulk;
    //    ErrorCheck();
    CheckEquilibrated();
    //cout << "H2OFluidProperties::UpdatePropertiesV_HighT() - \n";
    //DumpStatus();
    return;
  }

  void H2OFluidProperties::UpdatePropertiesVL()
  {
    bulk.InitToZero();
    state = VL;
    
    pcurrent  = pressure;
    tcurrent  = water.SaturationTemperatureFromP(pcurrent); // to be on the safe side
    hcurrent  = enthalpy;
    
    liq.t     = tcurrent;
    liq.p     = pcurrent;
    liq.x     = 0.0;
    liq.smf   = 0.0;
    liq.rho   = water.LiquidDensityFromP(pcurrent);
    liq.h     = water.LiquidEnthalpyFromP(pcurrent);
    liq.cp    = water.LiquidHeatCapacityFromP(pcurrent);
    // liq.dhdt  = liq.cp;
    liq.mu    = water.LiquidViscosityFromP(pcurrent);
    liq.beta  = water.LiquidCompressibilityFromP(pcurrent);
    liq.state = VL;
    
    vap.t     = tcurrent;
    vap.p     = pcurrent;
    vap.x     = 0.0;
    vap.smf   = 0.0;
    vap.rho   = water.VaporDensityFromP(pcurrent);
    vap.h     = water.VaporEnthalpyFromP(pcurrent);
    vap.cp    = water.VaporHeatCapacityFromP(pcurrent);
    // vap.dhdt  = vap.cp;
    vap.mu    = water.VaporViscosityFromP(pcurrent);
    vap.beta  = water.VaporCompressibilityFromP(pcurrent); 
    vap.state = VL;
    
    bulk.t     = liq.t;
    bulk.p     = liq.p;
    
    bulk.h     = enthalpy; // can be wrong if not properly updated outside !!!!!!!!
    bulk.s     = 1.0e0;    // not really, may be larger than porosity
    bulk.mf    = 1.0e0;
    bulk.state = VL;
    
    liq.mf     = (bulk.h-vap.h)/(liq.h-vap.h);
    vap.mf     = 1.0-liq.mf;
    liq.s      = (liq.mf/liq.rho)/(liq.mf/liq.rho+vap.mf/vap.rho);
    vap.s      = 1.0e0-liq.s;
    bulk.rho   = liq.s*liq.rho + vap.s*vap.rho;
    bulk.beta  = TwophaseCompressibility(); // get it externally!
    bulk.cp    = liq.mf*liq.cp+vap.mf*vap.cp;
    // bulk.dhdt  = bulk.cp; // or very high value ???

/*    if(vap.s < 0.05)
      {
        bulk.beta = vap.s/0.05*bulk.beta + (0.05-vap.s)/0.05*liq.beta;
      }*/
    
    //    ErrorCheck();
    equilibrated = true;
    // cout << "H2OFluidProperties::UpdatePropertiesVL() - \n";
    // DumpStatus();
    
    return;
  }


  void H2OFluidProperties::UpdatePropertiesToSaturatedLiquid()
  {
    bulk.InitToZero();
    vap.InitToZero();
    state     = F;
    
    pcurrent  = pressure;
    tcurrent  = water.SaturationTemperatureFromP(pcurrent); // to be on the safe side
    hcurrent  = enthalpy;
    
    liq.t     = tcurrent;
    liq.p     = pcurrent;
    liq.x     = 0.0;
    liq.smf   = 0.0;
    liq.rho   = water.LiquidDensityFromP(pcurrent);
    liq.h     = water.LiquidEnthalpyFromP(pcurrent);
    liq.cp    = water.LiquidHeatCapacityFromP(pcurrent);
    // liq.dhdt  = liq.cp;
    liq.mu    = water.LiquidViscosityFromP(pcurrent);
    liq.beta  = water.LiquidCompressibilityFromP(pcurrent);
    liq.s     = 1.0e0;
    liq.mf    = 1.0e0;
    
    liq.state = L;
    
    bulk      = liq;
    //    ErrorCheck();
    CheckEquilibrated();
    // cout << "H2OFluidProperties::UpdatePropertiesVL() - \n";
    // DumpStatus();
    
    return;
  }

  void H2OFluidProperties::UpdatePropertiesToSaturatedVapor()
  {
    bulk.InitToZero();
    liq.InitToZero();
    state = V;
    
    pcurrent  = pressure;
    tcurrent  = water.SaturationTemperatureFromP(pcurrent); // to be on the safe side
    hcurrent  = enthalpy;
    
    vap.t     = tcurrent;
    vap.p     = pcurrent;
    vap.x     = 0.0;
    vap.smf   = 0.0;
    vap.rho   = water.VaporDensityFromP(pcurrent);
    vap.h     = water.VaporEnthalpyFromP(pcurrent);
    vap.cp    = water.VaporHeatCapacityFromP(pcurrent);
    // vap.dhdt  = vap.cp;
    vap.mu    = water.VaporViscosityFromP(pcurrent);
    vap.beta  = water.VaporCompressibilityFromP(pcurrent); 
    vap.s     = 1.0e0;
    vap.mf    = 1.0e0;
    vap.state = V;
    
    bulk      = vap;
    //    ErrorCheck();
    CheckEquilibrated();
    // cout << "H2OFluidProperties::UpdatePropertiesVL() - \n";
    // DumpStatus();
    
    return;
  }

  
  void H2OFluidProperties::DumpStatus()
  {
    cerr << "H2OFluidProperties Status summary: \n";
    cerr << "State is             ";
    if(bulk.state == L)  cerr << "L\n";
    if(bulk.state == F)  cerr << "F\n";
    if(bulk.state == V)  cerr << "V\n";
    if(bulk.state == VL) cerr << "VL\n";
      
    cerr.setf(ios::scientific);
    cerr << "H2OFluidProperties::DumpStatus() Properties for\tbulk\tliq\tvap\tsalt\n";
    cerr << "state   = " << bulk.state << "\t" << liq.state <<  "\t" << vap.state << endl;
    cerr << "t       = " << bulk.t     << "\t" << liq.t     <<  "\t" << vap.t     << endl;
    cerr << "p       = " << bulk.p     << "\t" << liq.p     <<  "\t" << vap.p     << endl;
    cerr << "x       = " << bulk.x     << "\t" << liq.x     <<  "\t" << vap.x     << endl;
    cerr << "wt      = " << bulk.wt    << "\t" << liq.wt    <<  "\t" << vap.wt    << endl;
    cerr << "rho     = " << bulk.rho   << "\t" << liq.rho   <<  "\t" << vap.rho   << endl;
    cerr << "h       = " << bulk.h     << "\t" << liq.h     <<  "\t" << vap.h     << endl;
    cerr << "cp      = " << bulk.cp    << "\t" << liq.cp    <<  "\t" << vap.cp    << endl;
    // cerr << "dhdt    = " << bulk.dhdt  << "\t" << liq.dhdt  <<  "\t" << vap.dhdt  << endl;
    cerr << "beta    = " << bulk.beta  << "\t" << liq.beta  <<  "\t" << vap.beta  << endl;
    cerr << "s       = " << bulk.s     << "\t" << liq.s     <<  "\t" << vap.s     << endl;
    cerr << "mf      = " << bulk.mf    << "\t" << liq.mf    <<  "\t" << vap.mf    << endl;
    cerr << "mu      = " << bulk.mu    << "\t" << liq.mu    <<  "\t" << vap.mu    << endl;
      
    return;
  }
  
  void H2OFluidProperties::InitializeToBogus()
  {
    bulk.InitToBogus();    
    liq.InitToBogus();    
    vap.InitToBogus();    
  }
  
  double H2OFluidProperties::TwophaseCompressibility()
  {
    // Grant & Sorey (1979) two-phase compressibility.
    // Uses rock.HeatCapacity(T) at the current trial temperature so that the
    // liquidus/solidus latent-heat model (if enabled) is evaluated at the
    // correct temperature during the equilibrator's bisection loop.
    if(with_rock_liquidus_solidus)
        b  = ( 1.0 - phi ) * rock.HeatCapacity(liq.t, tl, ts) * rr;
    else
        b  = ( 1.0 - phi ) * rock.HeatCapacity(liq.t) * rr;

    b +=   phi * liq.rho * liq.cp * liq.s;
    b +=   phi * vap.rho * vap.cp * vap.s;
    b *= ( liq.rho - vap.rho ) / ( ( vap.h - liq.h ) * liq.rho * vap.rho );
    b *= ( liq.rho - vap.rho ) / ( ( vap.h - liq.h ) * liq.rho * vap.rho );
    b *= ( liq.t + 273.15 ) * 1.0 / phi;
    return b;    
  }
  
  void H2OFluidProperties::CheckEquilibrated()
  {
    // This does a final check if - after all other criteria are matched - also the
    // specific enthalpy is matched within a safety margin
    // Then, thermal equilibrium has been reached and equilibrated = true can
    // be issued to tell H2ONaClThermalEquilibrator that equilibration has converged
    if( fabs((enthalpy-bulk.h)/enthalpy) < 1.0e-4 ) equilibrated = true;
    else equilibrated = false;
  }
  
  void H2OFluidProperties::ErrorCheck()
  { 
    if( 
       definitelyLessThan( bulk.t, 0.0)
       ||
       definitelyLessThan( bulk.p, 0.0)
       ||
       definitelyLessThan( bulk.x, 0.0)
       ||
       definitelyLessThan( bulk.wt, 0.0)
       ||
       definitelyLessThan( bulk.rho, 0.0)
       ||
       definitelyLessThan( bulk.h, 0.0)
       ||
       definitelyLessThan( bulk.cp, 0.0)
       ||
       // definitelyLessThan( bulk.dhdt, 0.0)
       // ||
       definitelyLessThan( bulk.beta, 0.0)
       ||
       definitelyLessThan( bulk.s, 0.0)
       ||
       definitelyLessThan( bulk.mf, 0.0)
       ||
       definitelyLessThan( bulk.mu, 0.0)
       ||

       definitelyLessThan( liq.t, 0.0)
       ||
       definitelyLessThan( liq.p, 0.0)
       ||
       definitelyLessThan( liq.x, 0.0)
       ||
       definitelyLessThan( liq.wt, 0.0)
       ||
       definitelyLessThan( liq.rho, 0.0)
       ||
       definitelyLessThan( liq.h, 0.0)
       ||
       definitelyLessThan( liq.cp, 0.0)
       ||
       // definitelyLessThan( liq.dhdt, 0.0)
       // ||
       definitelyLessThan( liq.beta, 0.0)
       ||
       definitelyLessThan( liq.s, 0.0)
       ||
       definitelyLessThan( liq.mf, 0.0)
       ||
       definitelyLessThan( liq.mu, 0.0)
       ||

       definitelyLessThan( vap.t, 0.0)
       ||
       definitelyLessThan( vap.p, 0.0)
       ||
       definitelyLessThan( vap.x, 0.0)
       ||
       definitelyLessThan( vap.wt, 0.0)
       ||
       definitelyLessThan( vap.rho, 0.0)
       ||
       definitelyLessThan( vap.h, 0.0)
       ||
       definitelyLessThan( vap.cp, 0.0)
       ||
       // definitelyLessThan( vap.dhdt, 0.0)
       // ||
       definitelyLessThan( vap.beta, 0.0)
       ||
       definitelyLessThan( vap.s, 0.0)
       ||
       definitelyLessThan( vap.mf, 0.0)
       ||
       definitelyLessThan( vap.mu, 0.0)
        )
      {
        cerr << "H2OFluidProperties::ErrorCheck() found invalid result:\n\n";
        DumpStatus();
      }
    return;
  }


  // PW May 2016 - added functions for convergence speed-up

  bool H2OFluidProperties::AboveCriticalPressureH2O(double& pressure_external)
  {
	  return definitelyGreaterThan(pressure_external, cp_h2o.Pressure());
  }


  double H2OFluidProperties::BoilingPointAtPressure(double& pressure_external)
  {
	  return water.SaturationTemperatureFromP(pressure_external);
  }

  double H2OFluidProperties::VaporEnthalpyVL(double& pressure_external)
  {
	  return water.VaporEnthalpyFromP(pressure_external);
  }

  double H2OFluidProperties::LiquidEnthalpyVL(double& pressure_external)
  {
	  return water.LiquidEnthalpyFromP(pressure_external);
  }
// end added PW May 2016

  double H2OFluidProperties::SaturationPressureFromT(double temperature_C)
  {
      return water.SaturationPressureFromT(temperature_C);
}

void H2OFluidProperties::WithRockLiquidusSolidus(bool with_rock_liquidus_solidus_)
{
    with_rock_liquidus_solidus = with_rock_liquidus_solidus_;
}

void H2OFluidProperties::SetRockLiquidusSolidusTemperatures(double tl_, double ts_)
{
    tl = tl_;
    ts = ts_;
}

void H2OFluidProperties::SetRockHeatCapacity( double mini_cp )
{
    rock.SetRockHeatCapacity( mini_cp );
}

void H2OFluidProperties::SetRockCrystallizationCurve(double nu_coefficient, double sigma1_coefficient,
                                                     double latent_heat_of_fusion, double b_coefficient, std::string crystallization_curve_)
{
    rock.SetRockCrystallizationCurve(nu_coefficient, sigma1_coefficient,
                                     latent_heat_of_fusion, b_coefficient, crystallization_curve_);
  }


}//csmp
