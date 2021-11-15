#include "H2ONaClFluidProperties.h"
#include "compareFloats.h"
#include "ConvertConcentrationUnitsNaCl.h"

#include <limits>

using namespace std;

// old stuff, still to be checked
// Things to remember: check mass fraction stuff

// ********** this is probably not a crucial point anymore, nevertheless check: **********
// VERY CRUCIAL (???): sequence of LookupInitializations 

namespace csmp
{
  H2ONaClFluidProperties::H2ONaClFluidProperties(const double& externaltemperature_in_C,
                                                 const double& externalpressure_in_Pa,
                                                 const double& externalcomposition_in_mole_fraction,
                                                 const double& external_fluid_enthalpy_in_J_per_kg,
                                                 const double& external_cp_rock,
                                                 const double& external_rho_rock,
                                                 const double& external_phi, const bool& verbose ) 
  : temperature(externaltemperature_in_C),
    pressure(externalpressure_in_Pa),
    composition(externalcomposition_in_mole_fraction),
    enthalpy(external_fluid_enthalpy_in_J_per_kg),
    cpr(external_cp_rock),
    rr(external_rho_rock),
    phi(external_phi),
  //------------------
    safety_limit(5.0),
  //------------------
    tcurrent(0.0), 
    pcurrent(0.0),
    xcurrent(0.0), 
    hcurrent(0.0),
    dliqmfdp(0.0),
    dliqmfdt(0.0),
    dsaltmfdt(0.0),
    rl(0.0),
    rv(0.0),
    cpl(0.0),
    cpv(0.0),
    hl(0.0),
    hv(0.0),
    sl(0.0),
    sv(0.0),
    product(0.0),
    b(0.0),
    liq_dsmfdt(0.0), 
    vap_dsmfdt(0.0),
    t_vlh_low(0.0),
    t_vlh_high(0.0),
    h_vh(0.0),
    h_lh(0.0),
    h_vl(0.0),
    liq_mf(0.0),
    salt_mf(0.0),
    vap_mf(0.0),
    tsat(0.0),
  //    eqtype(0),
    equilibrated(false),
    vl_boundary_encountered(false),
    low_x(false),
    fatal(false),
    verbose(verbose),
    state(none),
    rock(), // CAUTION: RE-CHECK THIS!
    critcurve(tcurrent), 
    naclmelt_l(tcurrent,pcurrent), 
    naclmelt_h(tcurrent,pcurrent), 
    halite(tcurrent,pcurrent), 
    liquidus(tcurrent,pcurrent),
    vlh_v(tcurrent),
    vlh_l(tcurrent),
    vlh_h(tcurrent),
    naclsatvap(tcurrent,pcurrent), 
    twophase_l(tcurrent,pcurrent), 
    twophase_v(tcurrent,pcurrent),  
    brine(tcurrent,pcurrent,xcurrent),
    water(),
    fluid_water(tcurrent, pcurrent, hcurrent, cpr, rr, phi),
    vh_halite(tcurrent,pcurrent),
    lh_halite(tcurrent,pcurrent),     
    csmp_error( ErrorHandler::Instance() )
  {
    liq.InitToZero();
    vap.InitToZero();
    bulk.InitToZero();
    salt.InitToZero();
  }


  H2ONaClFluidProperties::~H2ONaClFluidProperties()
  {
  }


  double H2ONaClFluidProperties::BulkEnthalpy()
  {
    pcurrent = pressure;
    tcurrent = temperature;
    xcurrent = composition;
    hcurrent = enthalpy;
    low_x        = false;
    fatal        = false;

    if(definitelyLessThan( xcurrent, 2.0*numeric_limits<double>::epsilon() ) )
      {
        bulk               = fluid_water.BulkPropertiesFromTHP();
        return bulk.h;
      }

    else if( essentiallyEqual( xcurrent, 1.0e0 ) )
      {
        if( definitelyLessThan(tcurrent, naclmelt_l.TfromP(pcurrent), 5.0*numeric_limits<double>::epsilon()) )
          { UpdateEnthalpyHalite(); return bulk.h;}
        else if( definitelyGreaterThan(tcurrent, naclmelt_l.TfromP(pcurrent), 5.0*numeric_limits<double>::epsilon()) )
          { UpdateEnthalpyNaClMelt(); return bulk.h;}
        else if ( !definitelyLessThan( enthalpy, naclmelt_l.Enthalpy() ) )
          { UpdateEnthalpyNaClMelt(); return bulk.h;}
        else if ( !definitelyGreaterThan( enthalpy, naclmelt_h.Enthalpy() ) )
          { UpdateEnthalpyHalite(); return bulk.h;}
        else
          { UpdateEnthalpyNaClMeltingCurve_ForP(enthalpy); return bulk.h;}
      }

    else
      {
	
        if( definitelyLessThan(tcurrent, tp_nacl.Temperature() ))
          {
            if( definitelyGreaterThan( pcurrent, vlh_v.Pmax() ) )
              {
                if( definitelyGreaterThan( tcurrent, cp_h2o.Temperature() )
                    &&
                    definitelyLessThan( pcurrent, critcurve.Pressure() ) 
                    )
                  {
                    if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                      { 
                        //d cout << "H2ONaClFluidEnthalpy went into 0\n"; 
                        UpdateEnthalpyV(); 
                        //d Snapshot(0); 
                        return bulk.h;
                      } 
                    // Check if VL 
                    else if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 1\n";
                        UpdateEnthalpyVL(); 
                        //d Snapshot(1); 
                        return bulk.h;
                      }
                    // Check if L
                    else if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 2\n";
                        UpdateEnthalpyL(); 
                        //d Snapshot(2); 
                        return bulk.h; 
                      }
                    // If all previous checks passed, it must be LH
                    else
                      { //d cout << "H2ONaClFluidEnthalpy went into 3\n";
                        UpdateEnthalpyLH(); 
                        //d Snapshot(3); 
                        return bulk.h;
                      }
                  }

                else
                  {
                    // Check if F
                    if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 4\n";
                        UpdateEnthalpyF(); 
                        //d Snapshot(4); 
                        return bulk.h;
                      }
                    // Otherwise must be LH
                    else
                      { //d cout << "H2ONaClFluidEnthalpy went into 5\n";
                        UpdateEnthalpyLH(); 
                        //d Snapshot(5); 
                        return bulk.h;
                      }
                  }
              }

            // --------------- P identical to pressure maximum of the VLH surface ---------------
            else if ( essentiallyEqual( pcurrent, vlh_v.Pmax() ) )
              {
                //d cout << "H2ONaClFluidEnthalpy::UpdateEnthalpy() - p == pmax\n";
                //d  cout << "vlh 1\n";
                if( essentiallyEqual( tcurrent, vlh_v.Tmax() ) )
                  {
                    if( !definitelyGreaterThan( xcurrent, vlh_v.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 6\n";
                        UpdateEnthalpyV(); 
                        //d Snapshot(6); 
                        return bulk.h;
                      } 
                    else 
                      // *** DEBUG info: does this need an enthalpy check?
                      //                analogously to pure water case?
                      //                Probably yes, cross-check
                      { //d cout << "H2ONaClFluidEnthalpy went into 7\n";
                        UpdateEnthalpyAtVLH_For_P(enthalpy); 
                        //d Snapshot(7); 
                        return bulk.h;
                      }
                  }
                // ******************** Here an else is missing *********************
              }
            
            // --------------- P less than pressure maximum of the VLH surface ---------------
            else
              {
                //d cout << "H2ONaClFluidEnthalpy::UpdateEnthalpy() - p < pmax\n";
                t_vlh_low  = vlh_v.TfromP(pcurrent, 100.0);
                t_vlh_high = vlh_v.TfromP(pcurrent, 700.0);

                //************************************* re-check for complete control statements !!! *************************

                // Now VLH case #1
                if( essentiallyEqual( tcurrent, t_vlh_high, 5.0*numeric_limits<double>::epsilon()) )
                  { 
                    //d cout << "H2ONaClFluidEnthalpy::UpdateEnthalpy() - VLH case 1\n";
                    // Check if V only
                    if(!definitelyGreaterThan( xcurrent, vlh_v.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 14\n";
                        UpdateEnthalpyV(); 
                        //d Snapshot(14); 
                        return bulk.h;
                      }
	      
                    // Else compositionally, VLH is possible 
                    // 
                    else
                      {  
                        UpdateEnthalpyAtVLH_For_P(enthalpy); // dummy only to have all props available
                        // Compute h_vh, i.e., the enthalpy of a VH fluid at t_vlh, for testing purposes
                        salt_mf    = (bulk.smf-vap.smf)/(salt.smf-vap.smf);
                        vap_mf     = 1.0-salt_mf;
                        h_vh       = vap_mf*vap.h+salt_mf*salt.h;

                        // case 1: just pure liquid
                        if( essentiallyEqual( xcurrent, liq.x ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 15\n";
                            UpdateEnthalpyL(); 
                            //d Snapshot(15); 
                            return bulk.h;
                          }

                        // case 2: between vap.x and liq.x 
                        else if( !definitelyGreaterThan( xcurrent, liq.x ) )
                          {
                            liq_mf     = (bulk.smf-vap.smf)/(liq.smf-vap.smf);
                            vap_mf     = 1.0-liq_mf;
                            h_vl       = liq_mf*liq.h+vap_mf*vap.h;
			  
                            if( !definitelyLessThan( hcurrent, h_vl ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 16\n";
                                UpdateEnthalpyVL(); 
                                //d Snapshot(16);
                                return bulk.h;
                              } // enthalpy too high for VLH
                            else if( !definitelyGreaterThan( hcurrent, h_vh ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 17\n";
                                UpdateEnthalpyVH(); 
                                //d Snapshot(17); 
                                return bulk.h;
                              } // enthalpy too low for VLH
                            else // enthalpy right for VLH, this also implies, tcurrent is the desired temperature 
                              { //d cout << "H2ONaClFluidEnthalpy went into 98\n";
                                UpdateEnthalpyAtVLH_For_P(enthalpy); equilibrated = true; 
                                //d Snapshot(98); 
                                return bulk.h;
                              } // redo update because mfs were changed
                          }

                        // case 3: between liq.x and halite
                        else 
                          {
                            liq_mf     = (bulk.smf-salt.smf)/(liq.smf-salt.smf);
                            salt_mf    = 1.0-liq_mf;
                            h_lh       = liq_mf*liq.h+salt_mf*salt.h;
			  
                            if( !definitelyLessThan( hcurrent, h_lh ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 18\n";
                                UpdateEnthalpyLH(); 
                                //d Snapshot(18); 
                                return bulk.h;
                              } // enthalpy too high for VLH
                            else if( !definitelyGreaterThan( hcurrent, h_vh ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 19\n";
                                UpdateEnthalpyVH(); 
                                //d Snapshot(19); 
                                return bulk.h;
                              } // enthalpy too low for VLH
                            else // MUST be on VLH between xv and xl and is the desired temperature
                              { //d cout << "H2ONaClFluidEnthalpy went into 99\n";
                                UpdateEnthalpyAtVLH_For_P(enthalpy); equilibrated = true; 
                                //d Snapshot(99); 
                                return bulk.h;
                              } // redo update because mfs were changed	   
                          }
                      }
                  }

                // Now VLH case #2
                if( essentiallyEqual( tcurrent, t_vlh_low, 5.0*numeric_limits<double>::epsilon()) )
                  { 
                    //d cout << "H2ONaClFluidEnthalpy::UpdateEnthalpy() - VLH case 2\n";
                    // Check if V only
                    if(!definitelyGreaterThan( xcurrent, vlh_v.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 20\n";
                        UpdateEnthalpyV(); 
                        //d Snapshot(20); 
                        return bulk.h;
                      }

                    // Else compositionally, VLH is possible, 
                    // error check via negative saturations in UpdateEnthalpyAtVLH_For_P(enthalpy)
                    else
                      {
                        UpdateEnthalpyAtVLH_For_P(enthalpy); // dummy only to have all props available
                        // Compute h_vh, i.e., the enthalpy of a VH fluid at t_vlh, for testing purposes
                        salt_mf    = (bulk.smf-vap.smf)/(salt.smf-vap.smf);
                        vap_mf     = 1.0-salt_mf;
                        h_vh       = vap_mf*vap.h+salt_mf*salt.h;
                        // *** MAYBE introduce check for negative saturations???

                        // case 1: just pure liquid
                        if( essentiallyEqual( xcurrent, liq.x ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 21\n";
                            UpdateEnthalpyL(); 
                            //d Snapshot(21); 
                            return bulk.h;
                          }

                        // case 2: between vap.x and liq.x 
                        else if( !definitelyGreaterThan( xcurrent, liq.x ) )
                          {
                            liq_mf     = (bulk.smf-vap.smf)/(liq.smf-vap.smf);
                            vap_mf     = 1.0-liq_mf;
                            h_vl       = liq_mf*liq.h+vap_mf*vap.h;
			  
                            // this sequence should ensure that the above mentioned inconsistencies will not appear 
                            if( !definitelyGreaterThan( hcurrent, h_vl ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 22\n";
                                UpdateEnthalpyVL(); 
                                //d Snapshot(22); 
                                return bulk.h;
                              } // enthalpy too low for VLH
                            else if( !definitelyLessThan( hcurrent, h_vh ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 23\n";
                                UpdateEnthalpyVH(); 
                                //d Snapshot(23); 
                                return bulk.h;
                              } // enthalpy too high for VLH
                            else // enthalpy right for VLH, this also implies, tcurrent is the desired temperature 
                              { //d cout << "H2ONaClFluidEnthalpy went into 100\n";
                                UpdateEnthalpyAtVLH_For_P(enthalpy); equilibrated = true; 
                                //d Snapshot(100); 
                                return bulk.h;
                              } // redo update because mfs were changed
                          }

                        // case 3: between liq.x and halite
                        else 
                          {
                            liq_mf     = (bulk.smf-salt.smf)/(liq.smf-salt.smf);
                            salt_mf    = 1.0-liq_mf;
                            h_lh       = liq_mf*liq.h+salt_mf*salt.h;
			  
                            if( !definitelyGreaterThan( hcurrent, h_lh ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 24\n";
                                UpdateEnthalpyLH(); 
                                //d Snapshot(24); 
                                return bulk.h;
                              } // enthalpy too low for VLH
                            else if( !definitelyLessThan( hcurrent, h_vh ) )
                              { //d cout << "H2ONaClFluidEnthalpy went into 25\n";
                                UpdateEnthalpyVH(); 
                                //d Snapshot(25); 
                                return bulk.h;
                              } // enthalpy too high for VLH
                            else // MUST be on VLH between xv and xl and is the desired temperature
                              { //d cout << "H2ONaClFluidEnthalpy went into 101\n";
                                UpdateEnthalpyAtVLH_For_P(enthalpy); equilibrated = true; 
                                //d Snapshot(101); 
                                return bulk.h;
                              } // redo update because mfs were changed	   
                          }
                      }
                  }


                if( definitelyGreaterThan( tcurrent, t_vlh_high ) )
                  {
                    //d cout << "H2ONaClFluidEnthalpy::UpdateEnthalpy() - t > t_vlh_high\n";
                    // Check if vapor only
                    if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 8\n";
                        UpdateEnthalpyV();
                        //d cout << "twophase_v.MassFractionNaCl() = " << twophase_v.MassFractionNaCl() << endl;
                        //d Snapshot(8); 
                        return bulk.h;
                      } 
		    
                    // Check if VL
                    if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 9\n";
                        UpdateEnthalpyVL(); 
                        //d Snapshot(9); 
                        return bulk.h;
                      }
		    
                    // Check if L
                    if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 10\n";
                        //d cout << "xcurrent, liquidus.MassFractionNaCl() : " << xcurrent << "\t" << liquidus.MassFractionNaCl() << endl; 
                        UpdateEnthalpyL(); 
                        //d Snapshot(10); 
                        return bulk.h;
                      }
		    
                    // Else must be LH
                    else
                      { //d cout << "H2ONaClFluidEnthalpy went into 11\n";
                        UpdateEnthalpyLH(); 
                        //d Snapshot(11); 
                        return bulk.h;
                      }
                  }

                if( definitelyGreaterThan( tcurrent, t_vlh_low ) 
                    &&
                    definitelyLessThan(    tcurrent, t_vlh_high ) 
                    )
                  {
                    //d cout << "H2ONaClFluidEnthalpy::UpdateEnthalpy() - t between t_vlh_low and t_vlh_high\n";

                    // Check if vapor only
                    if( !definitelyGreaterThan( xcurrent, naclsatvap.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 12\n";
                        UpdateEnthalpyV(); 
                        //d Snapshot(12); 
                        return bulk.h;
                      }
		    
                    // Else must be VH
                    else
                      { //d cout << "H2ONaClFluidEnthalpy went into 13\n";
                        UpdateEnthalpyVH(); 
                        //d Snapshot(13); 
                        return bulk.h;
                      }
                  }


                //--------------------------------------------
                // Remaining cases are always at t < t_vlh_low!

                // First, the easy ones: p is higher than critical  pressure of H2O
                if( definitelyGreaterThan( pcurrent, cp_h2o.Pressure() ) )
                  {
                    if( definitelyGreaterThan( tcurrent, cp_h2o.Temperature() )
                        &&
                        definitelyLessThan( pcurrent, critcurve.Pressure() ) )
                      {
                        // - if pcurrent is less that the pressure on the critical curve at tcurrent, then we
                        //   can encounter V, VL, L, LH. Distinction between these cases should obviously be
                        //   based on xcurrent

                        // Check if V only
                        if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                          { 
                            //d cout << "H2ONaClFluidEnthalpy went into 26\n";
                            //d cout << "twophase_v.MassFractionNaCl() was " << twophase_v.MassFractionNaCl() << endl;
                            UpdateEnthalpyV(); 
                            //d Snapshot(26); 
                            return bulk.h; 
                          }
			
                        // Check if VL
                        if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                          { 
                            //d cout << "H2ONaClFluidEnthalpy went into 27\n";
                            UpdateEnthalpyVL(); 
                            //d Snapshot(27); 
                            return bulk.h;
                          }
			
                        // Check if L
                        if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                          { 
                            //d cout << "H2ONaClFluidEnthalpy went into 28\n";
                            UpdateEnthalpyL(); 
                            //d Snapshot(28); 
                            return bulk.h; 
                          }
			
                        // Else must be LH
                        else
                          { 
                            //d cout << "H2ONaClFluidEnthalpy went into 29\n";
                            UpdateEnthalpyLH(); 
                            //d Snapshot(29); 
                            return bulk.h;
                          }
                      }
		    
                    else
                      {
                        // - if pcurrent is greater or equal than the pressure on the critical curve (which btw has NO 
                        //   metastable extension at t < tcrit(H2O) in the Lookup version), only F and LH states are 
                        //   possible
                        //
                        // - pcurrent is essentially identical to critical pressure at tcurrent
                        //   this is sort of problematic because within numerical limits we could be
                        //   either in VL field or F if xcurrent is near xcrit.
                        //   The decision has be made that this is always weighted as F. Hence, again only the two
                        //   states F and LH need to be considered

                        // Check if F
                        if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 30\n";
                            UpdateEnthalpyF(); 
                            //d Snapshot(30); 
                            return bulk.h;  
                          }

                        // Else must be LH
                        else
                          { //d cout << "H2ONaClFluidEnthalpy went into 31\n";
                            UpdateEnthalpyLH(); 
                            //d Snapshot(31); 
                            return bulk.h; 
                          }
                      }
                  }
		    
                // Now the cases where p is less than or equal to critical pressure of H2O
                else
                  {
                    //		  cout << "caught the right general if\n";
                    tsat = water.SaturationTemperatureFromP( pcurrent );
                    // If less than tsat, topology is relatively easy
                    if( definitelyLessThan( tcurrent, tsat ) )
                      {
                        // Check if F
                        if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 32\n";
                            UpdateEnthalpyF(); 
                            //d Snapshot(32); 
                            return bulk.h;
                          }
                        // *** DEBUG info: does this need an enthalpy check for case near tsat?
                        //                analogously to pure water case?
                        //                Probably yes, cross-check

                        // Else must be LH
                        else
                          { //d cout << "H2ONaClFluidEnthalpy went into 33\n";
                            UpdateEnthalpyLH(); 
                            //d Snapshot(33); 
                            return bulk.h;
                          }
                      }

                    // *** DEBUG info: double check the following conditions !(...)
                    else if( /*! (definitelyGreaterThan( twophase_v.MassFractionNaCl(), 0.0 ) 
                               &&*/					
                            definitelyLessThan( xcurrent, twophase_l.MinXResolution() )
                            &&
                            // definitelyGreaterThan( twophase_l.MassFractionNaCl(), 0.0 ) 
                            // &&
                            definitelyLessThan( twophase_l.MassFractionNaCl(), twophase_l.MinXResolution() )
                            // && 
                            // !definitelyGreaterThan( xcurrent, twophase_l.MassFractionNaCl() )
                            // && 
                            // !definitelyLessThan( xcurrent, twophase_v.MassFractionNaCl() )

                            /*&&
                              definitelyGreaterThan( twophase_l.MassFractionNaCl(), twophase_v.MassFractionNaCl(), safety_limit*numeric_limits<double>::epsilon()) )*/
                             )
                      // At least one of the above conditions is not fulfilled, which very likely means very low x
                      // probably requires enthalpy-based special treatment

                      {
                        //		      cout << "failed vs. 1\n";
                        // ******* I am a bit worried that these enthalpy checks may let us get stuck in unlucky cases *********
                        //d cout << "***\n" << enthalpy << "\t" << twophase_l.Enthalpy() << "\t" << twophase_v.Enthalpy() << "***\n" << endl;
                        if( definitelyLessThan( enthalpy, /*twophase_l.Enthalpy()*/water.LiquidEnthalpyFromP(pcurrent)  ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 34\n";
                            UpdateEnthalpyL(); 
                            //d Snapshot(34); 
                            return bulk.h;
                          }

                        if( definitelyGreaterThan( enthalpy, /*twophase_v.Enthalpy()*/water.VaporEnthalpyFromP(pcurrent)  ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 35\n";
                            UpdateEnthalpyV(); 
                            //d Snapshot(35); 
                            return bulk.h;
                          }
		      
                        else
                          { //d cout << "H2ONaClFluidEnthalpy went into 36\n";
                            low_x = true; UpdateEnthalpyVL(); 
                            //d Snapshot(36); 
                            return bulk.h;
                          }
                      }
		  
                    else 
                      {
                        //		      cout << "failed vs. 2\n";
                        // Check if vapor only
                        if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 37\n";
                            UpdateEnthalpyV(); 
                            //d Snapshot(37); 
                            return bulk.h; 
                          } 
		      
                        // Check if VL
                        if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 38\n";
                            UpdateEnthalpyVL(); 
                            //d Snapshot(38); 
                            return bulk.h;
                          }
		      
                        // Check if L
                        if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                          { //d cout << "H2ONaClFluidEnthalpy went into 39\n";
                            UpdateEnthalpyL(); 
                            //d Snapshot(39); 
                            return bulk.h; 
                          }
		      
                        // Else must be LH
                        else
                          { //d cout << "H2ONaClFluidEnthalpy went into 40\n";
                            UpdateEnthalpyLH(); 
                            //d Snapshot(40); 
                            return bulk.h;
                          }
                      }
                  }
              }
          }
	
        else
          // t > t_triple_nacl
          {
            // if p >= pcrit, only single-phase fluid or LH possible
            if( !definitelyLessThan( pcurrent, critcurve.Pressure() ) )
              {
                if( !definitelyGreaterThan( tcurrent, naclmelt_h.TfromP( pcurrent ) ) )
                  {
                    if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                      { //d cout << "H2ONaClFluidEnthalpy went into 41\n";
                        UpdateEnthalpyF(); 
                        //d Snapshot(41); 
                        return bulk.h; 
                      }
                    else
                      { //d cout << "H2ONaClFluidEnthalpy went into 42\n";
                        UpdateEnthalpyLH(); 
                        //d Snapshot(42); 
                        return bulk.h;
                      }
                  }
                else
                  { //d cout << "H2ONaClFluidEnthalpy went into 43\n";
                    UpdateEnthalpyF(); 
                    //d Snapshot(43); 
                    return bulk.h; 
                  }
              }
            // if p < pcrit, V, VL and L are possible
            else
              {
                // Check if vapor only
                if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                  { //d cout << "H2ONaClFluidEnthalpy went into 44\n";
                    UpdateEnthalpyV(); 
                    //d Snapshot(44); 
                    return bulk.h; 
                  } 
	      
                // Check if VL
                if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                  { //d cout << "H2ONaClFluidEnthalpy went into 45\n";
                    UpdateEnthalpyVL(); 
                    //d Snapshot(45); 
                    return bulk.h;
                  }
	      
                // else must be L
                else
                  { //d cout << "H2ONaClFluidEnthalpy went into 46\n";
                    UpdateEnthalpyL(); 
                    //d Snapshot(46); 
                    return bulk.h; 
                  }
              }
          }
      }
    
      // SKM FIX: else return value may not be defined (compilation problem)
      return std::numeric_limits<double>::quiet_NaN();
  }


void H2ONaClFluidProperties::UpdateEnthalpyNaClMeltingCurve_ForP(double enthalpy)
  {
    // The correct assignment of bulk.h and other properties in this function depends on the
    // assumption that prior to calling this function there was a check that enthalpy is
    // indeed between salt.h and liq.h for pcurrent
    //
    // I strongly recommend not to fiddle around in this function or functions that call it
    // if you don't stick to this assumption
    //
    bulk.h       = enthalpy;
    return;
  }
  
  void H2ONaClFluidProperties::UpdateEnthalpyHalite(){
    bulk.h     = halite.Enthalpy();
    return;
  }


  void H2ONaClFluidProperties::UpdateEnthalpyNaClMelt(){
    bulk.h     = brine.Enthalpy();
    return;
  }

  void H2ONaClFluidProperties::UpdateEnthalpyAtVLH_For_P(double enthalpy){
    // TD, 06-Feb-2006, modified 12-Sep-2006

    // tcurrent, pcurrent, xcurrent are known externally when this function is called;
    // CAUTION: for isobaric conditions, there are two possible temperatures - which one is
    // correct needs determination by the program that calls H2ONaClFluidProperties and passes it on as
    // tcurrent
    bulk.h           = enthalpy; 
    return;
  }




  void H2ONaClFluidProperties::UpdateEnthalpyVH(){
    bulk.smf   = XNaCl2Massfraction(xcurrent);
    
    vap.smf    = XNaCl2Massfraction(naclsatvap.MassFractionNaCl());
    vap.h      = naclsatvap.Enthalpy();

    salt.smf   = 1.0e0;
    salt.h     = vh_halite.Enthalpy();

    salt.mf    = (bulk.smf-vap.smf)/(salt.smf-vap.smf);
    vap.mf     = 1.0-salt.mf;

    bulk.h     = vap.mf*vap.h+salt.mf*salt.h;
    return;
  }



  void H2ONaClFluidProperties::UpdateEnthalpyLH(){
    bulk.smf   = XNaCl2Massfraction(xcurrent);

    liq.smf    = XNaCl2Massfraction(liquidus.MassFractionNaCl());
    liq.h      = liquidus.Enthalpy();
    
    salt.smf   = 1.0e0;
    salt.h     = lh_halite.Enthalpy();
    
    salt.mf    = (bulk.smf-liq.smf)/(salt.smf-liq.smf);
    liq.mf     = 1.0-salt.mf;
    
    bulk.h     = liq.mf*liq.h+salt.mf*salt.h;
    return;
  }    


  void H2ONaClFluidProperties::UpdateEnthalpyL(){
    bulk.h     = brine.Enthalpy();
    return;
  }

  void H2ONaClFluidProperties::UpdateEnthalpyF(){
    bulk.h     = brine.Enthalpy();
    return;
  }

    
  void H2ONaClFluidProperties::UpdateEnthalpyV()
  {
    if( tcurrent > cp_h2o.Temperature() )
      { 
        bulk.h     = brine.Enthalpy();
      }
    else
      {
        bulk.h     = brine.SubcriticalVaporEnthalpy();
      }      
    return;
  }


  void H2ONaClFluidProperties::UpdateEnthalpyVL(){
    if( !low_x)
      {
        //d cout << "UpdateEnthalpyVL() -> !low_x\n";
        bulk.smf      = XNaCl2Massfraction(xcurrent);

        liq.smf       = XNaCl2Massfraction(twophase_l.MassFractionNaCl());
        liq.h         = twophase_l.Enthalpy();

        vap.smf       = XNaCl2Massfraction(twophase_v.MassFractionNaCl());
        vap.h         = twophase_v.Enthalpy();
	  
        liq.mf        = (bulk.smf-vap.smf)/(liq.smf-vap.smf);
        vap.mf        = 1.0-liq.mf;

        bulk.h        = liq.mf*liq.h  + vap.mf*vap.h;
        return;
      }
    else
      {
        // if(verbose)
        //   {
        //     cout << "*****************************\n";
        //     cout << "UpdateEnthalpyVL() -> low_x\n";
        //     cout << "*****************************\n";
        //   }
        bulk.h        = enthalpy;
        return;
      }
  }



  void H2ONaClFluidProperties::UpdateProperties()
  {
    // // using this way of comparing doubles
    // if( essentiallyEqual( pcurrent, pressure    ) && 
    // 	essentiallyEqual( tcurrent, temperature ) && 
    // 	essentiallyEqual( xcurrent, composition ) && 
    // 	essentiallyEqual( hcurrent, enthalpy    ) )
    //   {
    // 	return;
    //   }
    // else
    //   {
    pcurrent = pressure;
    tcurrent = temperature;
    xcurrent = composition;
    hcurrent = enthalpy;
    // }

    equilibrated = false;
    low_x        = false;
    fatal        = false;

    //d cout << "xcurrent = " << xcurrent << " " << numeric_limits<double>::epsilon() << endl;
    // ---------------------------------------------------------------------------------
    // If pure water, let H2ONaClFluidPropertiesWater do  the job
    // if( essentiallyEqual( xcurrent, 0.0 ) ) // "== 0.0"
    if(definitelyLessThan( xcurrent, 2.0*numeric_limits<double>::epsilon() ) )
      {
        //d cout << "H2ONaClFluidProperties::UpdateProperties() - pure water case\n";
        bulk               = fluid_water.BulkPropertiesFromTHP();
        liq                = fluid_water.LiquidPropertiesFromTHP();
        vap                = fluid_water.VaporPropertiesFromTHP();
        salt.InitToZero();
        equilibrated       = fluid_water.Equilibrated();
        return;
      }
    // This should be ok as long as H2ONaClFluidPropertiesWater works properly
    

    // ---------------------------------------------------------------------------------
    // If pure NaCl, do specialized stuff
    else if( essentiallyEqual( xcurrent, 1.0e0 ) )
      {
        //d cout << "H2ONaClFluidProperties::UpdateProperties() - pure NaCl case\n";
        // t-based check if halite
        if( definitelyLessThan(tcurrent, naclmelt_l.TfromP(pcurrent), 5.0*numeric_limits<double>::epsilon()) )
          { UpdatePropertiesHalite(); return; }

        // t-based check if NaCl melt
        else if( definitelyGreaterThan(tcurrent, naclmelt_l.TfromP(pcurrent), 5.0*numeric_limits<double>::epsilon()) )
          { UpdatePropertiesNaClMelt(); return; }

        // remainig possibilities: very near true tmelt, hence, double-check via enthalpy
        else if ( !definitelyLessThan( enthalpy, naclmelt_l.Enthalpy() ) )
          { UpdatePropertiesNaClMelt(); return; }

        else if ( !definitelyGreaterThan( enthalpy, naclmelt_h.Enthalpy() ) )
          { UpdatePropertiesHalite(); return; }
	  
        else
          // true NaCl Vl-twophase
          { UpdatePropertiesNaClMeltingCurve_ForP(enthalpy); return; }
      }
    // This should be ok as long as enthalpy was computed correctly



    // ---------------------------------------------------------------------------------
    // Else, do H2O-NaCl phase determination and property assignment
    else{
	
      if( definitelyLessThan(tcurrent, tp_nacl.Temperature() ))
        {
          // If temperature is lower than triple point temperature of NaCl, we may envisage the
          // three-phase VLH equilibrium.
          //
          // Based on pressure, we will distinguish different topologies of the phase diagram in
          // T-X space (the one in which thermal equilibration is done
	  
          // --------------- P > pressure maximum of the VLH surface ---------------
          //
          if( definitelyGreaterThan( pcurrent, vlh_v.Pmax() ) )
            {
              //d cout << "H2ONaClFluidProperties::UpdateProperties() - p > pmax\n";
              //d cout << "vlh 0\n";
              // At pressure higher than the pressure maximum of the VLH surface, no complications by 
              // VLH or low pressure VL are possible.
              // Further distiction of regions in the isobaric T-X plane is again based on pressure:
              if( definitelyGreaterThan( tcurrent, cp_h2o.Temperature() )
                  &&
                  definitelyLessThan( pcurrent, critcurve.Pressure() ) 
                  )
                {
                  // - if pcurrent is less than the pressure on the critical curve AT TCURRENT, then we
                  //   can encounter V, VL, L, LH. Distinction between these cases should obviously be
                  //   based on xcurrent
		  
                  // Check if only vapor
                  if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                    { 
                      //d cout << "H2ONaClFluidProperties went into 0\n"; 
                      UpdatePropertiesV(); 
                      //d Snapshot(0); 
                      return; } 
		    
                  // Check if VL 
                  else if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 1\n";
                      UpdatePropertiesVL(); 
                      //d Snapshot(1); 
                      return; }

                  // Check if L
                  else if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 2\n";
                      UpdatePropertiesL(); 
                      //d Snapshot(2); 
                      return;  }

                  // If all previous checks passed, it must be LH
                  else
                    { //d cout << "H2ONaClFluidProperties went into 3\n";
                      UpdatePropertiesLH(); 
                      //d Snapshot(3); 
                      return; }
                }

              else
                {
                  // - if pcurrent is greater than the pressure on the critical curve (which btw has NO metastable 
                  //   extension at t < tcrit(H2O) in the Lookup version), only F and LH states are possible
                  //
                  // - pcurrent is essentially identical to critical pressure at tcurrent
                  //   this is sort of problematic because within numerical limits we could be
                  //   either in VL field or F if xcurrent is near xcrit.
                  //   The decision has been made that this is always weighted as F. Hence, again only the two
                  //   states F and LH need to be considered

                  // Check if F
                  if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 4\n";
                      UpdatePropertiesF(); 
                      //d Snapshot(4); 
                      return; }

                  // Otherwise must be LH
                  else
                    { //d cout << "H2ONaClFluidProperties went into 5\n";
                      UpdatePropertiesLH(); 
                      //d Snapshot(5); 
                      return; }
                }
            }

          // --------------- P identical to pressure maximum of the VLH surface ---------------
          else if ( essentiallyEqual( pcurrent, vlh_v.Pmax() ) )
            {
              //d cout << "H2ONaClFluidProperties::UpdateProperties() - p == pmax\n";
              //d  cout << "vlh 1\n";
              if( essentiallyEqual( tcurrent, vlh_v.Tmax() ) )
                {
                  if( !definitelyGreaterThan( xcurrent, vlh_v.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 6\n";
                      UpdatePropertiesV(); 
                      //d Snapshot(6); 
                      return; } 

                  else 
                    // *** DEBUG info: does this need an enthalpy check?
                    //                analogously to pure water case?
                    //                Probably yes, cross-check
                    { //d cout << "H2ONaClFluidProperties went into 7\n";
                      UpdatePropertiesAtVLH_For_P(enthalpy); 
                      //d Snapshot(7); 
                      return; }
                }

            }
	    
          // --------------- P less than pressure maximum of the VLH surface ---------------
          else
            {
              //d cout << "H2ONaClFluidProperties::UpdateProperties() - p < pmax\n";
	
              // Now we will have complications arising from determining precisely a VLH state PLUS
              // we can run into problems related to VL determination at low X at pressure below
              // the critical pressure of pure H2O
              // Now is the time to introduce VLH checks
              t_vlh_low  = vlh_v.TfromP(pcurrent, 100.0);
              t_vlh_high = vlh_v.TfromP(pcurrent, 700.0);
              // cout << "t_vlh_low        = " << t_vlh_low  << endl;
              // cout << "t_vlh_high       = " << t_vlh_high << endl;
              // cout << "tcurrent         = " << tcurrent   << endl;
              // cout << "current-low rel  = " << (tcurrent-t_vlh_low)/tcurrent << endl;
              // cout << "high-current rel = " << (t_vlh_high-tcurrent)/tcurrent << endl;
              // cout << "5*epsilon        = " << 5.0*numeric_limits<double>::epsilon() << endl;

              // Now VLH case #1
              if( essentiallyEqual( tcurrent, t_vlh_high, 5.0*numeric_limits<double>::epsilon()) )
                { 
                  //d cout << "H2ONaClFluidProperties::UpdateProperties() - VLH case 1\n";
                  // Check if V only
                  if(!definitelyGreaterThan( xcurrent, vlh_v.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 14\n";
                      UpdatePropertiesV(); 
                      //d Snapshot(14); 
                      return; }
	      
                  // Else compositionally, VLH is possible 
                  // 
                  else
                    {  
                      UpdatePropertiesAtVLH_For_P(enthalpy); // dummy only to have all props available
                      // Compute h_vh, i.e., the enthalpy of a VH fluid at t_vlh, for testing purposes
                      salt_mf    = (bulk.smf-vap.smf)/(salt.smf-vap.smf);
                      vap_mf     = 1.0-salt_mf;
                      h_vh       = vap_mf*vap.h+salt_mf*salt.h;

                      // case 1: just pure liquid
                      if( essentiallyEqual( xcurrent, liq.x ) )
                        { //d cout << "H2ONaClFluidProperties went into 15\n";
                          UpdatePropertiesL(); 
                          //d Snapshot(15); 
                          return; }

                      // case 2: between vap.x and liq.x 
                      else if( !definitelyGreaterThan( xcurrent, liq.x ) )
                        {
                          liq_mf     = (bulk.smf-vap.smf)/(liq.smf-vap.smf);
                          vap_mf     = 1.0-liq_mf;
                          h_vl       = liq_mf*liq.h+vap_mf*vap.h;
			  
                          if( !definitelyLessThan( hcurrent, h_vl ) )
                            { //d cout << "H2ONaClFluidProperties went into 16\n";
                              UpdatePropertiesVL(); 
                              //d Snapshot(16);
                              // cout << "hcurrent, bulk.h " << hcurrent << "\t" << bulk.h << endl;
                              // cout << "h_vh, h_vl " << h_vh << "\t" << h_vl << endl;
                              // cout << xcurrent << "\t" << liq.x << endl;
                              return; } // enthalpy too high for VLH
                          else if( !definitelyGreaterThan( hcurrent, h_vh ) )
                            { //d cout << "H2ONaClFluidProperties went into 17\n";
                              UpdatePropertiesVH(); 
                              //d Snapshot(17); 
                              return; } // enthalpy too low for VLH
                          else // enthalpy right for VLH, this also implies, tcurrent is the desired temperature 
                            { //d cout << "H2ONaClFluidProperties went into 98\n";
                              UpdatePropertiesAtVLH_For_P(enthalpy); equilibrated = true; 
                              //d Snapshot(98); 
                              return; } // redo update because mfs were changed
                        }

                      // case 3: between liq.x and halite
                      else 
                        {
                          liq_mf     = (bulk.smf-salt.smf)/(liq.smf-salt.smf);
                          salt_mf    = 1.0-liq_mf;
                          h_lh       = liq_mf*liq.h+salt_mf*salt.h;
			  
                          if( !definitelyLessThan( hcurrent, h_lh ) )
                            { //d cout << "H2ONaClFluidProperties went into 18\n";
                              UpdatePropertiesLH(); 
                              //d Snapshot(18); 
                              return; } // enthalpy too high for VLH
                          else if( !definitelyGreaterThan( hcurrent, h_vh ) )
                            { //d cout << "H2ONaClFluidProperties went into 19\n";
                              UpdatePropertiesVH(); 
                              //d Snapshot(19); 
                              return; } // enthalpy too low for VLH
                          else // MUST be on VLH between xv and xl and is the desired temperature
                            { //d cout << "H2ONaClFluidProperties went into 99\n";
                              UpdatePropertiesAtVLH_For_P(enthalpy); equilibrated = true; 
                              //d Snapshot(99); 
                              return; } // redo update because mfs were changed	   
                        }
                    }
                }

              // Now VLH case #2
              if( essentiallyEqual( tcurrent, t_vlh_low, 5.0*numeric_limits<double>::epsilon()) )
                { 
                  //d cout << "H2ONaClFluidProperties::UpdateProperties() - VLH case 2\n";
                  // Check if V only
                  if(!definitelyGreaterThan( xcurrent, vlh_v.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 20\n";
                      UpdatePropertiesV(); 
                      //d Snapshot(20); 
                      return; }

                  // Else compositionally, VLH is possible, 
                  // error check via negative saturations in UpdatePropertiesAtVLH_For_P(enthalpy)
                  else
                    {
                      UpdatePropertiesAtVLH_For_P(enthalpy); // dummy only to have all props available
                      // Compute h_vh, i.e., the enthalpy of a VH fluid at t_vlh, for testing purposes
                      salt_mf    = (bulk.smf-vap.smf)/(salt.smf-vap.smf);
                      vap_mf     = 1.0-salt_mf;
                      h_vh       = vap_mf*vap.h+salt_mf*salt.h;
                      // *** MAYBE introduce check for negative saturations???

                      // case 1: just pure liquid
                      if( essentiallyEqual( xcurrent, liq.x ) )
                        { //d cout << "H2ONaClFluidProperties went into 21\n";
                          UpdatePropertiesL(); 
                          //d Snapshot(21); 
                          return; }

                      // case 2: between vap.x and liq.x 
                      else if( !definitelyGreaterThan( xcurrent, liq.x ) )
                        {
                          liq_mf     = (bulk.smf-vap.smf)/(liq.smf-vap.smf);
                          vap_mf     = 1.0-liq_mf;
                          h_vl       = liq_mf*liq.h+vap_mf*vap.h;
			  
                          // this sequence should ensure that the above mentioned inconsistencies will not appear 
                          if( !definitelyGreaterThan( hcurrent, h_vl ) )
                            { //d cout << "H2ONaClFluidProperties went into 22\n";
                              UpdatePropertiesVL(); 
                              //d Snapshot(22); 
                              return; } // enthalpy too low for VLH
                          else if( !definitelyLessThan( hcurrent, h_vh ) )
                            { //d cout << "H2ONaClFluidProperties went into 23\n";
                              UpdatePropertiesVH(); 
                              //d Snapshot(23); 
                              return; } // enthalpy too high for VLH
                          else // enthalpy right for VLH, this also implies, tcurrent is the desired temperature 
                            { //d cout << "H2ONaClFluidProperties went into 100\n";
                              UpdatePropertiesAtVLH_For_P(enthalpy); equilibrated = true; 
                              //d Snapshot(100); 
                              return; } // redo update because mfs were changed
                        }

                      // case 3: between liq.x and halite
                      else 
                        {
                          liq_mf     = (bulk.smf-salt.smf)/(liq.smf-salt.smf);
                          salt_mf    = 1.0-liq_mf;
                          h_lh       = liq_mf*liq.h+salt_mf*salt.h;
			  
                          if( !definitelyGreaterThan( hcurrent, h_lh ) )
                            { //d cout << "H2ONaClFluidProperties went into 24\n";
                              UpdatePropertiesLH(); 
                              //d Snapshot(24); 
                              return; } // enthalpy too low for VLH
                          else if( !definitelyLessThan( hcurrent, h_vh ) )
                            { //d cout << "H2ONaClFluidProperties went into 25\n";
                              UpdatePropertiesVH(); 
                              //d Snapshot(25); 
                              return; } // enthalpy too high for VLH
                          else // MUST be on VLH between xv and xl and is the desired temperature
                            { //d cout << "H2ONaClFluidProperties went into 101\n";
                              UpdatePropertiesAtVLH_For_P(enthalpy); equilibrated = true; 
                              //d Snapshot(101); 
                              return; } // redo update because mfs were changed	   
                        }
                    }
                }


              if( definitelyGreaterThan( tcurrent, t_vlh_high ) )
                {
                  //d cout << "H2ONaClFluidProperties::UpdateProperties() - t > t_vlh_high\n";
                  // Check if vapor only
                  if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 8\n";
                      UpdatePropertiesV();
                      //d cout << "twophase_v.MassFractionNaCl() = " << twophase_v.MassFractionNaCl() << endl;
                      //d Snapshot(8); 
                      return; } 
		    
                  // Check if VL
                  if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 9\n";
                      UpdatePropertiesVL(); 
                      //d Snapshot(9); 
                      return; }
		    
                  // Check if L
                  if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 10\n";
                      //d cout << "xcurrent, liquidus.MassFractionNaCl() : " << xcurrent << "\t" << liquidus.MassFractionNaCl() << endl; 
                      UpdatePropertiesL(); 
                      //d Snapshot(10); 
                      return; }
		    
                  // Else must be LH
                  else
                    { //d cout << "H2ONaClFluidProperties went into 11\n";
                      UpdatePropertiesLH(); 
                      //d Snapshot(11); 
                      return; }
                }

              if( definitelyGreaterThan( tcurrent, t_vlh_low ) 
                  &&
                  definitelyLessThan(    tcurrent, t_vlh_high ) 
                  )
                {
                  //d cout << "H2ONaClFluidProperties::UpdateProperties() - t between t_vlh_low and t_vlh_high\n";

                  // Check if vapor only
                  if( !definitelyGreaterThan( xcurrent, naclsatvap.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 12\n";
                      UpdatePropertiesV(); 
                      //d Snapshot(12); 
                      return; }
		    
                  // Else must be VH
                  else
                    { //d cout << "H2ONaClFluidProperties went into 13\n";
                      UpdatePropertiesVH(); 
                      //d Snapshot(13); 
                      return; }
                }


              //--------------------------------------------
              // Remaining cases are always at t < t_vlh_low!

              // First, the easy ones: p is higher than critical  pressure of H2O
              if( definitelyGreaterThan( pcurrent, cp_h2o.Pressure() ) )
                {
                  if( definitelyGreaterThan( tcurrent, cp_h2o.Temperature() )
                      &&
                      definitelyLessThan( pcurrent, critcurve.Pressure() ) )
                    {
                      // - if pcurrent is less that the pressure on the critical curve at tcurrent, then we
                      //   can encounter V, VL, L, LH. Distinction between these cases should obviously be
                      //   based on xcurrent

                      // Check if V only
                      if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                        { 
                          //d cout << "H2ONaClFluidProperties went into 26\n";
                          //d cout << "twophase_v.MassFractionNaCl() was " << twophase_v.MassFractionNaCl() << endl;
                          UpdatePropertiesV(); 
                          //d Snapshot(26); 
                          return;  }
			
                      // Check if VL
                      if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                        { 
                          //d cout << "H2ONaClFluidProperties went into 27\n";
                          UpdatePropertiesVL(); 
                          //d Snapshot(27); 
                          return; }
			
                      // Check if L
                      if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                        { 
                          //d cout << "H2ONaClFluidProperties went into 28\n";
                          UpdatePropertiesL(); 
                          //d Snapshot(28); 
                          return;  }
			
                      // Else must be LH
                      else
                        { 
                          //d cout << "H2ONaClFluidProperties went into 29\n";
                          UpdatePropertiesLH(); 
                          //d Snapshot(29); 
                          return; }
                    }
		    
                  else
                    {
                      // - if pcurrent is greater or equal than the pressure on the critical curve (which btw has NO 
                      //   metastable extension at t < tcrit(H2O) in the Lookup version), only F and LH states are 
                      //   possible
                      //
                      // - pcurrent is essentially identical to critical pressure at tcurrent
                      //   this is sort of problematic because within numerical limits we could be
                      //   either in VL field or F if xcurrent is near xcrit.
                      //   The decision has be made that this is always weighted as F. Hence, again only the two
                      //   states F and LH need to be considered

                      // Check if F
                      if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                        { //d cout << "H2ONaClFluidProperties went into 30\n";
                          UpdatePropertiesF(); 
                          //d Snapshot(30); 
                          return;   }

                      // Else must be LH
                      else
                        { //d cout << "H2ONaClFluidProperties went into 31\n";
                          UpdatePropertiesLH(); 
                          //d Snapshot(31); 
                          return;  }
                    }
                }
		    
              // Now the cases where p is less than or equal to critical pressure of H2O
              else
                {
                  //		  cout << "caught the right general if\n";
                  tsat = water.SaturationTemperatureFromP( pcurrent );
                  // If less than tsat, topology is relatively easy
                  if( definitelyLessThan( tcurrent, tsat ) )
                    {
                      // Check if F
                      if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                        { //d cout << "H2ONaClFluidProperties went into 32\n";
                          UpdatePropertiesF(); 
                          //d Snapshot(32); 
                          return; }
                      // *** DEBUG info: does this need an enthalpy check for case near tsat?
                      //                analogously to pure water case?
                      //                Probably yes, cross-check

                      // Else must be LH
                      else
                        { //d cout << "H2ONaClFluidProperties went into 33\n";
                          UpdatePropertiesLH(); 
                          //d Snapshot(33); 
                          return; }
                    }

                  // *** DEBUG info: double check the following conditions !(...)
                  else if( /*! (definitelyGreaterThan( twophase_v.MassFractionNaCl(), 0.0 ) 
                             &&*/					
                          definitelyLessThan( xcurrent, twophase_l.MinXResolution() )
                          &&
                          // definitelyGreaterThan( twophase_l.MassFractionNaCl(), 0.0 ) 
                          // &&
                          definitelyLessThan( twophase_l.MassFractionNaCl(), twophase_l.MinXResolution() )
                          // && 
                          // !definitelyGreaterThan( xcurrent, twophase_l.MassFractionNaCl() )
                          // && 
                          // !definitelyLessThan( xcurrent, twophase_v.MassFractionNaCl() )

                          /*&&
                            definitelyGreaterThan( twophase_l.MassFractionNaCl(), twophase_v.MassFractionNaCl(), safety_limit*numeric_limits<double>::epsilon()) )*/
                           )
                    // At least one of the above conditions is not fulfilled, which very likely means very low x
                    // probably requires enthalpy-based special treatment

                    {
                      //		      cout << "failed vs. 1\n";
                      // ******* I am a bit worried that these enthalpy checks may let us get stuck in unlucky cases *********
                      //d cout << "***\n" << enthalpy << "\t" << twophase_l.Enthalpy() << "\t" << twophase_v.Enthalpy() << "***\n" << endl;
                      if( definitelyLessThan( enthalpy, /*twophase_l.Enthalpy()*/water.LiquidEnthalpyFromP(pcurrent)  ) )
                        { //d cout << "H2ONaClFluidProperties went into 34\n";
                          UpdatePropertiesL(); 
                          //d Snapshot(34); 
                          return; }

                      if( definitelyGreaterThan( enthalpy, /*twophase_v.Enthalpy()*/water.VaporEnthalpyFromP(pcurrent)  ) )
                        { //d cout << "H2ONaClFluidProperties went into 35\n";
                          UpdatePropertiesV(); 
                          //d Snapshot(35); 
                          return; }
		      
                      else
                        { //d cout << "H2ONaClFluidProperties went into 36\n";
                          low_x = true; UpdatePropertiesVL(); 
                          //d Snapshot(36); 
                          return; }
                    }
		  
                  else 
                    {
                      //		      cout << "failed vs. 2\n";
                      // Check if vapor only
                      if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                        { //d cout << "H2ONaClFluidProperties went into 37\n";
                          UpdatePropertiesV(); 
                          //d Snapshot(37); 
                          return;  } 
		      
                      // Check if VL
                      if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                        { //d cout << "H2ONaClFluidProperties went into 38\n";
                          UpdatePropertiesVL(); 
                          //d Snapshot(38); 
                          return; }
		      
                      // Check if L
                      if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                        { //d cout << "H2ONaClFluidProperties went into 39\n";
                          UpdatePropertiesL(); 
                          //d Snapshot(39); 
                          return;  }
		      
                      // Else must be LH
                      else
                        { //d cout << "H2ONaClFluidProperties went into 40\n";
                          UpdatePropertiesLH(); 
                          //d Snapshot(40); 
                          return; }
		      
                    }
                }
            }
        }
	
      else
        // t > t_triple_nacl
        {
          // if p >= pcrit, only single-phase fluid or LH possible
          if( !definitelyLessThan( pcurrent, critcurve.Pressure() ) )
            {
              if( !definitelyGreaterThan( tcurrent, naclmelt_h.TfromP( pcurrent ) ) )
                {
                  if( !definitelyGreaterThan( xcurrent, liquidus.MassFractionNaCl() ) )
                    { //d cout << "H2ONaClFluidProperties went into 41\n";
                      UpdatePropertiesF(); 
                      //d Snapshot(41); 
                      return;  }
                  else
                    { //d cout << "H2ONaClFluidProperties went into 42\n";
                      UpdatePropertiesLH(); 
                      //d Snapshot(42); 
                      return; }
                }
              else
                { //d cout << "H2ONaClFluidProperties went into 43\n";
                  UpdatePropertiesF(); 
                  //d Snapshot(43); 
                  return;  }
            }
          // if p < pcrit, V, VL and L are possible
          else
            {
              // Check if vapor only
              if( !definitelyGreaterThan( xcurrent, twophase_v.MassFractionNaCl() ) )
                { //d cout << "H2ONaClFluidProperties went into 44\n";
                  UpdatePropertiesV(); 
                  //d Snapshot(44); 
                  return;  } 
	      
              // Check if VL
              if( definitelyLessThan( xcurrent, twophase_l.MassFractionNaCl() ) )
                { //d cout << "H2ONaClFluidProperties went into 45\n";
                  UpdatePropertiesVL(); 
                  //d Snapshot(45); 
                  return; }
	      
              // else must be L
              else
                { //d cout << "H2ONaClFluidProperties went into 46\n";
                  UpdatePropertiesL(); 
                  //d Snapshot(46); 
                  return;  }
            }
        }
    }
  }

  void H2ONaClFluidProperties::UpdatePropertiesNaClMeltingCurve_ForP(double enthalpy)
  {
    // The correct assignment of bulk.h and other properties in this function depends on the
    // assumption that prior to calling this function there was a check that enthalpy is
    // indeed between salt.h and liq.h for pcurrent
    //
    // I strongly recommend not to fiddle around in this function or functions that call it
    // if you don't stick to this assumption
    //
    // T. Driesner, 24 January 2011
    
    vap.InitToZero();
    
    // "state"   is for internal use only; phase properties will be set to LH to ensure 
    //           that external modules in the context of FV transport etc. know the involved 
    //           phases (those don't know M as of January 2011)
    state        = HM; 
    pcurrent     = pressure;
    tcurrent     = naclmelt_l.TfromP(pcurrent);
    xcurrent     = 1.0; 
    
    bulk.t       = tcurrent;
    bulk.p       = pressure;
    bulk.x       = xcurrent; 
    // Caution:  don't modify the following line unless you really know what you are doing
    //           later on, CheckEquilibrated() is queried that compares bulk.h and enthalpy; if true
    //           H2ONaClThermalEquilibrator will think it converged.
    //           see other caution note in first line of function body a few lines above here!
    bulk.h       = enthalpy;
    bulk.smf     = 1.0e0;
    bulk.wt      = 100.0*bulk.smf;
    bulk.s       = 1.0e0; // not really, may be larger than porosity
    bulk.mf      = 1.0e0;
    bulk.mu      = 0.0e0; // because it makes no sense
    bulk.state   = LH;
    
    liq          = bulk;
    // TD: For the next few re-check if you want to revert to the commented version
    liq.rho      = brine.Density();//naclmelt_l.Density(); 
    liq.h        = brine.Enthalpy();//naclmelt_l.Enthalpy();
    liq.cp       = brine.HeatCapacity();//naclmelt_l.HeatCapacity();
    //speed    liq.dhdt     = liq.cp; // or very high?
    liq.mu       = brine.Viscosity();//naclmelt_l.Viscosity();
    liq.beta     = brine.Compressibility();//naclmelt_l.Compressibility();
    liq.state    = LH;
    
    salt.t       = bulk.t;
    salt.p       = pressure;
    salt.x       = bulk.x;
    salt.smf     = bulk.smf;
    salt.wt      = 100.0*salt.smf;
    salt.rho     = halite.Density();
    salt.h       = halite.Enthalpy();
    salt.cp      = halite.HeatCapacity();
    //speed    salt.dhdt    = salt.cp;
    salt.mu      = 1.0e20; // needs P-adjustment
    salt.state   = LH;
    salt.beta    = halite.Compressibility(); 
    
    liq.mf       = (bulk.h-salt.h)/(liq.h-salt.h);
    salt.mf      = 1.0-liq.mf;
    liq.s        = (liq.mf/liq.rho)/(liq.mf/liq.rho+salt.mf/salt.rho);
    salt.s       = 1.0e0-liq.s;
    bulk.rho     = liq.s*liq.rho + salt.s*salt.rho;
    bulk.beta    = TwophaseCompressibility(); // get it externally!
    //speed    bulk.dhdt    = 0.0e0; // or very high value ???
    
    equilibrated = true; // by definition, see above comments on enthalpy
    return;
  }
  
  
  
  void H2ONaClFluidProperties::UpdatePropertiesHalite(){
    vap.InitToZero();
    liq.InitToZero();
    //    bulk.InitToZero();
    //    salt.InitToZero();
    state      = H;
    pcurrent   = pressure;
    tcurrent   = temperature;
    xcurrent   = 1.0;

    bulk.t     = temperature;
    bulk.p     = pressure;
    bulk.x     = 1.0; 
    bulk.h     = halite.Enthalpy();
    bulk.smf   = 1.0;
    bulk.wt    = 100.0*bulk.smf;
    bulk.s     = 1.0e0; // not really, may be larger than porosity
    bulk.mf    = 1.0e0;
    bulk.mu    = 0.0e0; // because it makes no sense
    bulk.state = H;
    bulk.rho   = halite.Density();
    bulk.cp    = halite.HeatCapacity();
    //speed    bulk.dhdt  = bulk.cp; 
    bulk.mu    = 1.0e20;
    bulk.beta  = halite.Compressibility();
    salt       = bulk;
    CheckEquilibrated();
    return;
  }


  void H2ONaClFluidProperties::UpdatePropertiesNaClMelt(){
    vap.InitToZero();
    //    liq.InitToZero();
    //    bulk.InitToZero();
    salt.InitToZero();

    state      = L;
    // state should rather be M but to avoid confusion in dependent codes
    // assign L
    pcurrent   = pressure;
    tcurrent   = temperature;
    xcurrent   = 1.0;
    bulk.t     = temperature;
    bulk.p     = pressure;
    bulk.x     = 1.0; 
    bulk.h     = brine.Enthalpy();
    bulk.smf   = 1.0;
    bulk.wt    = 100.0*bulk.smf;
    bulk.s     = 1.0e0; // not really, may be larger than porosity
    bulk.mf    = 1.0e0;
    bulk.state = L;
    bulk.rho   = brine.Density();
    bulk.cp    = brine.HeatCapacity();
    //speed    bulk.dhdt  = bulk.cp; 
    bulk.mu    = brine.Viscosity();
    bulk.beta  = brine.Compressibility();
    liq        = bulk;
    CheckEquilibrated();
    return;
  }





  void H2ONaClFluidProperties::UpdatePropertiesAtVLH_For_P(double enthalpy){
    // TD, 06-Feb-2006, modified 12-Sep-2006

    // tcurrent, pcurrent, xcurrent are known externally when this function is called;
    // CAUTION: for isobaric conditions, there are two possible temperatures - which one is
    // correct needs determination by the program that calls H2ONaClFluidProperties and passes it on as
    // tcurrent

    state            = VLH;

    pcurrent         = pressure;
    tcurrent         = temperature;
    xcurrent         = composition; 
    
    bulk.state       = VLH;
    bulk.t           = tcurrent;
    bulk.p           = pcurrent;
    bulk.x           = xcurrent;
    bulk.smf         = XNaCl2Massfraction(bulk.x);
    bulk.wt          = 100.0*bulk.smf;
    // ********************** caution: "equilibrated" logic***************************
    bulk.h           = enthalpy; 
    bulk.cp          = 0.0; // rather: not defined
    //speed      bulk.dhdt        = 0.0; // rather: infinity
    bulk.s           = 1.0e0;
    bulk.mf          = 1.0e0;
    bulk.mu          = 0.0e0; // not defined
    
    vap.state        = VLH;
    vap.t            = temperature;
    vap.p            = pressure;
    vap.x            = vlh_v.MassFractionNaCl();
    // watch this: xcurrent is modified!
    xcurrent         = vap.x;
    vap.smf          = XNaCl2Massfraction(vap.x);
    vap.wt           = 100.0*vap.smf;
    vap.rho          = vlh_v.Density();
    vap.h            = vlh_v.Enthalpy();
    vap.cp           = vlh_v.HeatCapacity();
    //speed    vap.dhdt         = vlh_v.DEnthalpyDT();
    vap.beta         = vlh_v.Compressibility();
    vap.mu           = vlh_v.Viscosity();

    liq.state        = VLH;
    liq.t            = temperature;
    liq.p            = pressure;
    liq.x            = vlh_l.MassFractionNaCl();
    // watch this: xcurrent is modified!
    xcurrent         = liq.x;
    liq.smf          = XNaCl2Massfraction(liq.x);
    liq.wt           = 100.0*liq.smf;
    liq.rho          = vlh_l.Density();
    liq.h            = vlh_l.Enthalpy();
    liq.cp           = vlh_l.HeatCapacity();
    //speed    liq.dhdt         = vlh_l.DEnthalpyDT();
    liq.mu           = vlh_l.Viscosity();

    salt.state       = VLH;
    salt.t           = temperature;
    salt.p           = pressure;
    salt.x           = 1.0e0;
    salt.smf         = 1.0e0;
    salt.wt          = 100.0*salt.smf;
    // salt.rho         = halite.Density();
    // salt.h           = halite.Enthalpy();
    // salt.cp          = halite.HeatCapacity();
    // salt.dhdt        = vlh_h.DEnthalpyDT();
    // salt.beta        = halite.Compressibility()*1.0e-5;
    salt.rho         = vlh_h.Density();
    salt.h           = vlh_h.Enthalpy();
    salt.cp          = vlh_h.HeatCapacity();
    //speed    salt.dhdt        = vlh_h.DEnthalpyDT();
    salt.beta        = vlh_h.Compressibility()*1.0e-5;
    salt.mu          = 1.0e10;

    liq.mf           = bulk.h*(vap.smf-1.0)+vap.h*(1.0-bulk.smf)+salt.h*(bulk.smf-vap.smf);
    liq.mf          /= liq.h*(vap.smf-1.0)+vap.h*(1.0-liq.smf)+salt.h*(liq.smf-vap.smf);
    vap.mf           = (bulk.smf-1.0+liq.mf*(1.0-liq.smf))/(vap.smf-1.0);
    salt.mf          = 1.0-liq.mf-vap.mf;

    bulk.rho         = 1.0/(liq.mf/liq.rho+vap.mf/vap.rho+salt.mf/salt.rho);

    liq.s            = liq.mf*bulk.rho/liq.rho;
    vap.s            = vap.mf*bulk.rho/vap.rho;
    salt.s           = 1.0e0-liq.s-vap.s;

    // reset although not strictly necessary
    xcurrent         = composition;

    // new threephase compressibility
    double dvdh;
    dvdh             = (1.0-vap.smf)/liq.rho;      
    dvdh            += (liq.smf-1.0)/vap.rho;      
    dvdh            += (vap.smf-liq.smf)/salt.rho; 
    dvdh            /= liq.h*(1.0-vap.smf)+vap.h*(liq.smf-1.0)+salt.h*(vap.smf-liq.smf); // Nenner

    double cptot  = liq.mf*liq.rho*liq.cp;
    cptot           += vap.mf*vap.rho*vap.cp;
    cptot           += salt.mf*salt.rho*salt.cp;
    cptot           *= phi;
    cptot           += (1.0-phi)*rock.HeatCapacity(tcurrent)/*cpr*/*rr;
    double dpdt   = vlh_l.DPressureDT()*1.0e5;

    bulk.beta        = 1.0/phi;// the minus sign is skipped becaus Grant&Sorey forgot the added/released convention
    bulk.beta       *= dvdh;// *phi
    bulk.beta       *= cptot;
    bulk.beta       /= dpdt;
    if(bulk.beta < 0 )
      {
        cout << "bulk.beta < 0 for VLH; Threephase compressibility components: \n";
        cout << "dvdh    = " << dvdh  << endl;
        cout << "cptot   = " << cptot << endl;
        cout << "dpdt    = " << dpdt  << endl;
        cout << "t-tmax  = " << bulk.t - vlh_l.Tmax() << endl;
      }  
    // end new threephase compressibility

    return;
  }




  void H2ONaClFluidProperties::UpdatePropertiesVH(){
    //    bulk.InitToZero();
    liq.InitToZero();
    //    vap.InitToZero();
    //    salt.InitToZero();
    //d cout << "upon entering H2ONaClFluidProperties::UpdatePropertiesVH(), enthalpy is " << enthalpy << ", " << hcurrent << endl;
    state      = VH;

    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = xcurrent;
    bulk.smf   = XNaCl2Massfraction(bulk.x);
    bulk.wt    = 100.0*bulk.smf;
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    bulk.state = VH;
    
    vap.t      = tcurrent;
    vap.p      = pressure;
    vap.x      = naclsatvap.MassFractionNaCl();
    // watch this: xcurrent is modified!
    xcurrent   = vap.x;
    vap.smf    = XNaCl2Massfraction(vap.x);
    vap.wt     = 100.0*vap.smf;
    //speed    vap_dsmfdt = naclsatvap.DSaltMassFractionDT();
    //speed    vap.dhdt   = naclsatvap.DEnthalpyDT();
    vap.rho    = naclsatvap.Density();
    vap.h      = naclsatvap.Enthalpy();
    vap.cp     = naclsatvap.HeatCapacity();
    vap.beta   = naclsatvap.Compressibility();
    vap.mu     = naclsatvap.Viscosity();
    vap.state  = VH;

    // reset
    xcurrent   = composition;

    salt.t     = tcurrent;
    salt.p     = pressure;
    salt.x     = 1.0e0;
    salt.smf   = 1.0e0;
    salt.wt     = 100.0*salt.smf;
    salt.rho   = vh_halite.Density();
    salt.h     = vh_halite.Enthalpy();
    salt.cp    = vh_halite.HeatCapacity();
    //speed    salt.dhdt  = vh_halite.HeatCapacity();
    salt.beta  = vh_halite.Compressibility()*1.0e-5;
    // salt.rho   = halite.Density();
    // salt.h     = halite.Enthalpy();
    // salt.cp    = halite.HeatCapacity();
    // salt.dhdt  = halite.HeatCapacity();
    // salt.beta  = halite.Compressibility()*1.0e-5;
    salt.mu    = 1.0e10;
    salt.state = VH;
    
    salt.mf    = (bulk.smf-vap.smf)/(salt.smf-vap.smf);
    vap.mf     = 1.0-salt.mf;
    salt.s     = salt.mf/salt.rho/(salt.mf/salt.rho+vap.mf/vap.rho);
    vap.s      = 1.0e0 - salt.s;
    bulk.rho   = vap.s*vap.rho + salt.s*salt.rho;
    //speed    dsaltmfdt  = vap_dsmfdt*(-salt.smf+bulk.smf)/(salt.smf-vap.smf)/(salt.smf-vap.smf);
    bulk.beta  = TwophaseCompressibility();
    if(bulk.beta < 0.0e0)
      {
        DumpStatus();
        //	csmp_error.notice( FATAL_ERROR,
        csmp_error.notice( WARNING,
                        "H2ONaClFluidProperties::UpdatePropertiesVH() -",
                        "bulk compressibility < 0 for above condition: Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
        bulk.beta  = -bulk.beta; // just simple workaround if you downgrade fatal error to warning or less
      }
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0;
    bulk.mu    = vap.mu; // how to implement better?
    bulk.h     = vap.mf*vap.h+salt.mf*salt.h;
    //speed    bulk.dhdt  = vap.dhdt*(1.0-salt.mf)+dsaltmfdt*(salt.h-vap.h)+salt.mf*salt.dhdt;
    CheckEquilibrated();

    return;
  }



  void H2ONaClFluidProperties::UpdatePropertiesLH(){
    vap.InitToZero();

    state      = LH; // must be set before bulk.beta is computed
	    
    bulk.t     = tcurrent; 
    bulk.p     = pressure; 
    bulk.x     = xcurrent; 
    bulk.smf   = XNaCl2Massfraction(bulk.x);
    bulk.wt    = 100.0*bulk.smf;
    bulk.cp    = 0.0; // rather: not defined
    bulk.state = LH;
    
    liq.t      = tcurrent;
    liq.p      = pressure;
    liq.x      = liquidus.MassFractionNaCl();
    liq.smf    = XNaCl2Massfraction(liq.x);
    liq.wt     = 100.0*liq.smf;
    //speed    liq_dsmfdt = liquidus.DSaltMassFractionDT();
    //    liq.dwtdt  = 100.0*liquidus.DSaltMassFractionDT();
    //speed    liq.dhdt   = liquidus.DEnthalpyDT();
    // watch this: xcurrent is modified!
    xcurrent   = liq.x;
    liq.rho    = liquidus.Density();
    liq.h      = liquidus.Enthalpy();
    liq.cp     = liquidus.HeatCapacity();
    liq.beta   = liquidus.Compressibility();
    liq.mu     = liquidus.Viscosity();
    liq.state  = LH;
    //reset
    xcurrent   = composition;
    
    salt.t     = tcurrent;
    salt.p     = pressure;
    salt.x     = 1.0e0;
    salt.smf   = 1.0e0;
    salt.wt    = 100.0*salt.smf;
    // salt.rho   = halite.Density();
    // salt.h     = halite.Enthalpy();
    // salt.cp    = halite.HeatCapacity();
    // salt.dhdt  = salt.cp;
    // salt.beta  = halite.Compressibility()*1.0e-5;
    salt.rho   = lh_halite.Density();
    salt.h     = lh_halite.Enthalpy();
    salt.cp    = lh_halite.HeatCapacity();
    //speed    salt.dhdt  = salt.cp;
    salt.beta  = lh_halite.Compressibility()*1.0e-5;
    salt.mu    = 1.e10;
    salt.state = LH;
    
    salt.mf    = (bulk.smf-liq.smf)/(salt.smf-liq.smf);
    liq.mf     = 1.0-salt.mf;
    salt.s     = salt.mf/salt.rho/(salt.mf/salt.rho+liq.mf/liq.rho);
    liq.s      = 1.0e0 - salt.s;
    bulk.rho   = liq.s*liq.rho + salt.s*salt.rho;

    //speed    dsaltmfdt  = liq_dsmfdt*(-salt.smf+bulk.smf)/(salt.smf-liq.smf)/(salt.smf-liq.smf);

    bulk.s     = 1.0e0;
    bulk.mf    = 1.0;
    bulk.mu    = liq.mu; // how to implement better?
    
    bulk.h     = liq.mf*liq.h+salt.mf*salt.h;
    //speed    bulk.dhdt  = -dsaltmfdt*liq.h
    //speed      +liq.dhdt*(1.0-salt.mf)+dsaltmfdt*salt.h+salt.mf*salt.dhdt;
    //d cout << "H2ONaClFluidProperties::UpdatePropertiesLH() - hbulk,liq,salt" << bulk.h << ", " << liq.h << ", " << salt.h << endl;
    bulk.beta  = TwophaseCompressibility();
    if(bulk.beta < 0.0e0)
      {
        DumpStatus();
        //	csmp_error.notice( FATAL_ERROR,
        csmp_error.notice( WARNING,
                        "H2ONaClFluidProperties::UpdatePropertiesLH() -",
                        "bulk compressibility < 0 for above condition: Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
        bulk.beta  = -bulk.beta; // just simple workaround if you downgrade fatal error to warning or less
      }

    CheckEquilibrated();

    return;
  }    


  void H2ONaClFluidProperties::UpdatePropertiesL(){
    vap.InitToZero();
    salt.InitToZero();
    state      = L;
    bulk.state = L;

    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = xcurrent;
    bulk.smf   = XNaCl2Massfraction(bulk.x);
    bulk.wt    = 100.0*bulk.smf;
    bulk.rho   = brine.Density();
    bulk.h     = brine.Enthalpy();
    bulk.cp    = brine.HeatCapacity();
    //speed    bulk.dhdt  = bulk.cp;
    bulk.beta  = brine.Compressibility();

    if(bulk.beta < 0.0e0)
      {
        DumpStatus();
        //	csmp_error.notice( FATAL_ERROR,
        csmp_error.notice( WARNING,
                        "H2ONaClFluidProperties::UpdatePropertiesL() -",
                        "liquid compressibility < 0 for above condition: Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
        bulk.beta  = -bulk.beta; // just simple workaround if you downgrade fatal error to warning or less
      }
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    bulk.mu    = brine.Viscosity();
    
    liq        = bulk;
    return;
  }

  void H2ONaClFluidProperties::UpdatePropertiesF(){
    vap.InitToZero();
    salt.InitToZero();
    state      = F;
    bulk.state = F;
    
    bulk.t     = tcurrent;
    bulk.p     = pressure;
    bulk.x     = xcurrent;
    bulk.smf   = XNaCl2Massfraction(bulk.x);
    bulk.wt    = 100.0*bulk.smf;
    //d cout << "H2ONaClFluidProperties::UpdatePropertiesF() calling brine.Density()\n";
    bulk.rho   = brine.Density();
    //d cout << "H2ONaClFluidProperties::UpdatePropertiesF() calling brine.Enthalpy()\n";
    bulk.h     = brine.Enthalpy();
    //d cout << "H2ONaClFluidProperties::UpdatePropertiesF() calling brine.HeatCapacity()\n";
    bulk.cp    = brine.HeatCapacity();
    //speed    bulk.dhdt  = brine.HeatCapacity();
    bulk.beta  = brine.Compressibility();
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    bulk.mu    = brine.Viscosity();    
    liq        = bulk;

    if(bulk.beta < 0.0e0)
      {
        DumpStatus();
        //	csmp_error.notice( FATAL_ERROR,
        csmp_error.notice( WARNING,
                        "H2ONaClFluidProperties::UpdatePropertiesF() -",
                        "liquid compressibility < 0 for above condition: Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
        bulk.beta  = -bulk.beta; // just simple workaround if you downgrade fatal error to warning or less
      }
    CheckEquilibrated();
    return;
  }

    
  void H2ONaClFluidProperties::UpdatePropertiesV()
  {
    liq.InitToZero();
    salt.InitToZero();
    state        = V;
      
    pcurrent     = pressure; // why is this here?
      
    bulk.t       = tcurrent;
    bulk.p       = pressure;
    bulk.x       = xcurrent;
    bulk.smf     = XNaCl2Massfraction(bulk.x);
    bulk.wt      = 100.0*bulk.smf;
      
    if( tcurrent > cp_h2o.Temperature() )
      { 
        bulk.rho   = brine.Density();
        bulk.beta  = brine.Compressibility();
        bulk.h     = brine.Enthalpy();
        bulk.cp    =brine.HeatCapacity();
        //speed    bulk.dhdt  = brine.HeatCapacity();
        bulk.mu    = brine.Viscosity();
      }
    else
      {
        bulk.rho   = brine.SubcriticalVaporDensity();
        bulk.beta  = brine.SubcriticalVaporCompressibility();
        bulk.h     = brine.SubcriticalVaporEnthalpy();
        //speed    bulk.dhdt  = brine.SubcriticalVaporHeatCapacity();
        bulk.cp    = brine.SubcriticalVaporHeatCapacity();
        bulk.mu    = brine.SubcriticalVaporViscosity();
      }
      
    if(bulk.beta < 0.0e0)
      {
        if( tcurrent > cp_h2o.Temperature() && pcurrent < 50.0 ) bulk.beta = - bulk.beta; // ** quick fix
        else
          {
            DumpStatus();
            //	    csmp_error.notice( FATAL_ERROR,
            csmp_error.notice( WARNING,
                            "H2ONaClFluidProperties::UpdatePropertiesV() -",
                            "vapor compressibility < 0 for above condition: Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
            // simple workaround if you want to downgrade fatal error ...
            bulk.beta = -bulk.beta;
          }
      
      }
    bulk.s     = 1.0e0;
    bulk.mf    = 1.0e0;
    bulk.state = V;
    vap        = bulk;
    CheckEquilibrated();
    return;
  }


  void H2ONaClFluidProperties::UpdatePropertiesVL(){
    salt.InitToZero();
    state         = VL;

    bulk.t        = tcurrent;
    bulk.p        = pressure;
    
    if( !low_x)
      {
        //d cout << "UpdatePropertiesVL() -> !low_x\n";
        bulk.x        = xcurrent;
        bulk.smf      = XNaCl2Massfraction(bulk.x);
        bulk.wt       = 100.0*bulk.smf;
        bulk.s        = 1.0e0; // obvious
        bulk.mf       = 1.0e0; // maybe these can be set by default, so that the assignment doesn't need to be done ...
        bulk.mu       = 0.0e0; // because it makes no sense; but maybe there is something in Garg&Pritchett?
        bulk.state    = VL;
        //d cout << "UpdatePropertiesVL() -> did block 1\n";
        liq.t         = tcurrent;
        liq.p         = pressure;
        liq.x         = twophase_l.MassFractionNaCl();
        liq.smf       = XNaCl2Massfraction(liq.x);
        liq.wt        = 100.0*liq.smf;
        //speed    liq_dsmfdt    = twophase_l.DSaltMassFractionDT();
        //speed      liq.dhdt      = twophase_l.DEnthalpyDT();
        // watch this: xcurrent is modified!
        xcurrent      = liq.x;
        liq.rho       = twophase_l.Density();
        liq.h         = twophase_l.Enthalpy();
        liq.cp        = twophase_l.HeatCapacity();
        liq.beta      = twophase_l.Compressibility();
        liq.mu        = twophase_l.Viscosity();
        liq.state     = VL;
        //d cout << "UpdatePropertiesVL() -> did block 2\n";
	  
        vap.t         = tcurrent;
        vap.p         = pressure;
        vap.x         = twophase_v.MassFractionNaCl();
        vap.smf       = XNaCl2Massfraction(vap.x);
        vap.wt        = 100.0*vap.smf;
        //speed    vap_dsmfdt    = twophase_v.DSaltMassFractionDT();
        //speed      vap.dhdt      = twophase_v.DEnthalpyDT();
        // watch this: xcurrent is modified!
        xcurrent      = vap.x;
        vap.rho       = twophase_v.Density();
        vap.h         = twophase_v.Enthalpy();
        vap.cp        = twophase_v.HeatCapacity();
        vap.beta      = twophase_v.Compressibility();
        vap.mu        = twophase_v.Viscosity();
        vap.state     = VL;
        // reset xcurrent
        xcurrent      = composition;
        //d cout << "UpdatePropertiesVL() -> did block 3\n";
	  
        liq.mf        = (bulk.smf-vap.smf)/(liq.smf-vap.smf);
	  
        vap.mf        = 1.0-liq.mf;
        liq.s         = (liq.mf/liq.rho)/(liq.mf/liq.rho+vap.mf/vap.rho);
        vap.s         = 1.0-liq.s;
        //speed    dliqmfdt      = -(vap_dsmfdt)/(liq.smf-vap.smf)-(bulk.smf-vap.smf)*(liq_dsmfdt-vap_dsmfdt)/pow(liq.smf-vap.smf,2.0);
	  
        bulk.rho      = liq.s*liq.rho + vap.s*vap.rho;
        bulk.h        = liq.mf*liq.h  + vap.mf*vap.h;
        //speed      bulk.dhdt     = vap.dhdt + liq.mf*(liq.dhdt-vap.dhdt) + dliqmfdt*(liq.h-vap.h); //ok 13Nov2006
        bulk.cp       = 0.0; // rather: not defined
        //d cout << "UpdatePropertiesVL() -> did block 3a\n";
        bulk.beta     = TwophaseCompressibility();
        //d cout << "UpdatePropertiesVL() -> did block 4\n";
        CheckEquilibrated();
        //d cout << "UpdatePropertiesVL() -> DumpStatus\n";
        //d DumpStatus();
        //d cout << "pcrit = " << critcurve.Pressure() << endl; 
        return;
      }
    else
      {
        // if(verbose)
        //   {
        //     cout << "*****************************\n";
        //     cout << "UpdatePropertiesVL() -> low_x\n";
        //     cout << "*****************************\n";
        //   }
        bulk.x        = xcurrent;
        bulk.smf      = XNaCl2Massfraction(bulk.x);
        bulk.wt       = 100.0*bulk.smf;
        bulk.s        = 1.0e0; // obvious
        bulk.mf       = 1.0e0; // maybe these can be set by default, so that the assignment doesn't need to be done ...
        bulk.mu       = 0.0e0; // because it makes no sense; but maybe there is something in Garg&Pritchett?
        bulk.state    = VL;

        // this has been pre-checked to be twophase, -> we can use enthalpy for saturations
        // ans since we are at very low x this is done with the pure water values
        bulk.h        = enthalpy;

        // to reduce expesnive calls to the "fromp" functions in "water", we go via "fromt"
        double th2o = water.SaturationTemperatureFromP( pcurrent );
        liq.h         = water.LiquidEnthalpyFromT( th2o );
        vap.h         = water.VaporEnthalpyFromT(  th2o );
        liq.rho       = water.LiquidDensityFromT(  th2o );
        vap.rho       = water.VaporDensityFromT(   th2o );

        liq.mf        = (bulk.h-vap.h)/(liq.h-vap.h);
        vap.mf        = 1.0-liq.mf;
        liq.s         = (liq.mf/liq.rho)/(liq.mf/liq.rho+vap.mf/vap.rho);
        vap.s         = 1.0e0-liq.s;

        double   kd = -1.11605*pow(cp_h2o.Temperature()-tsat,0.15201)-0.027387*(cp_h2o.Temperature()-tsat);
        kd            = pow(10.0,kd);
        liq.x         = xcurrent/(kd+liq.mf*(1.0-kd));
        vap.x         = kd*liq.x;
        liq.smf       = XNaCl2Weight(liq.x)*0.01;
        vap.smf       = XNaCl2Weight(vap.x)*0.01;
        liq.wt        = liq.smf*100.;
        vap.wt        = vap.smf*100.;
	
        liq.t         = tcurrent;
        liq.p         = pressure;
        vap.t         = tcurrent;
        vap.p         = pressure;

        // since x is so low, there is no chance to represent vapor salinity
        // numerically, hence I put all salt into liquid
        // liq.wt        = 1.0/liq.mf*bulk.wt;
        // liq.x         = Weight2XNaCl(liq.wt);
        // liq.smf       = XNaCl2Massfraction(liq.x);

        // vap.x         = 0.0;
        // vap.smf       = 0.0;
        // vap.wt        = 0.0;

        bulk.rho      = liq.s*liq.rho + vap.s*vap.rho;
        //speed      bulk.dhdt     = vap.dhdt + liq.mf*(liq.dhdt-vap.dhdt) + dliqmfdt*(liq.h-vap.h); //ok 13Nov2006
        bulk.cp       = 0.0; // rather: not defined

        //rather undefined:	  
        //speed    liq_dsmfdt    = twophase_l.DSaltMassFractionDT();
        //speed      liq.dhdt      = twophase_l.DEnthalpyDT();
        // watch this: xcurrent is modified!
        liq.cp        = water.LiquidHeatCapacityFromT( th2o );
        liq.beta      = water.LiquidCompressibilityFromT( th2o );
        liq.mu        = water.LiquidViscosityFromT( th2o );
        liq.state     = VL;

        //rather undefined:	  
        //speed    vap_dsmfdt    = twophase_v.DSaltMassFractionDT();
        //speed      vap.dhdt      = twophase_v.DEnthalpyDT();
        // watch this: xcurrent is modified!
        vap.cp        = water.VaporHeatCapacityFromT( th2o );
        vap.beta      = water.VaporCompressibilityFromT( th2o );
        vap.mu        = water.VaporViscosityFromT( th2o );
        vap.state     = VL;
	  
        //speed    dliqmfdt      = -(vap_dsmfdt)/(liq.smf-vap.smf)-(bulk.smf-vap.smf)*(liq_dsmfdt-vap_dsmfdt)/pow(liq.smf-vap.smf,2.0);
        bulk.beta     = TwophaseCompressibility();
	  
        // if(vap.s < 0.05)
        //   {
        //     bulk.beta = vap.s/0.05*bulk.beta + (0.05-vap.s)/0.05*liq.beta;
        //   }

        CheckEquilibrated();
	  
        return;
      }

  }



  void H2ONaClFluidProperties::DumpStatus(){
    cerr << "H2ONaClFluidProperties Status summary: \n";
    cerr << "State is             ";
    if(bulk.state == L) cerr << "L\n";
    if(bulk.state == F) cerr << "F\n";
    if(bulk.state == V) cerr << "V\n";
    if(bulk.state == VL) cerr << "VL\n";
    if(bulk.state == VH) cerr << "VH\n";
    if(bulk.state == LH) cerr << "LH\n";
    if(bulk.state == VLH) cerr << "VLH\n";

    cerr.setf(ios::scientific);
    cerr.precision(20);
    cerr << "H2ONaClFluidProperties::DumpStatus() Properties for\tbulk\tliq\tvap\tsalt\n";
    cerr << "state   = " << bulk.state << "\t" << liq.state <<  "\t" << vap.state <<  "\t"  << salt.state << endl;
    cerr << "t       = " << bulk.t     << "\t" << liq.t     <<  "\t" << vap.t     <<  "\t"  << salt.t     << endl;
    cerr << "p       = " << bulk.p     << "\t" << liq.p     <<  "\t" << vap.p     <<  "\t"  << salt.p     << endl;
    cerr << "x       = " << bulk.x     << "\t" << liq.x     <<  "\t" << vap.x     <<  "\t"  << salt.x     << endl;
    cerr << "wt      = " << bulk.wt    << "\t" << liq.wt    <<  "\t" << vap.wt    <<  "\t"  << salt.wt    << endl;
    cerr << "smf     = " << bulk.smf    << "\t" << liq.smf    <<  "\t" << vap.smf    <<  "\t"  << salt.smf    << endl;
    cerr << "rho     = " << bulk.rho   << "\t" << liq.rho   <<  "\t" << vap.rho   <<  "\t"  << salt.rho   << endl;
    cerr << "h       = " << bulk.h     << "\t" << liq.h     <<  "\t" << vap.h     <<  "\t"  << salt.h     << endl;
    cerr << "cp      = " << bulk.cp    << "\t" << liq.cp    <<  "\t" << vap.cp    <<  "\t"  << salt.cp    << endl;
    //speed          cerr << "dhdt    = " << bulk.dhdt  << "\t" << liq.dhdt  <<  "\t" << vap.dhdt  <<  "\t"  << salt.dhdt  << endl;
    cerr << "beta    = " << bulk.beta  << "\t" << liq.beta  <<  "\t" << vap.beta  <<  "\t"  << salt.beta  << endl;
    cerr << "s       = " << bulk.s     << "\t" << liq.s     <<  "\t" << vap.s     <<  "\t"  << salt.s     << endl;
    cerr << "mf      = " << bulk.mf    << "\t" << liq.mf    <<  "\t" << vap.mf    <<  "\t"  << salt.mf    << endl;
    cerr << "mu      = " << bulk.mu    << "\t" << liq.mu    <<  "\t" << vap.mu    <<  "\t"  << salt.mu    << endl;

    return;
  }

  void H2ONaClFluidProperties::Snapshot( const int& is)
  {
    // if(verbose)
    //   {
    // 	cout << "H2ONaClFluidProperties::Snapshot( const int& is): ";
    // 	cout << "eqtype was " << is << endl;
    // 	cout << "tcurrent = " << tcurrent << endl;
    // 	cout << "pcurrent = " << pcurrent << endl;
    // 	cout << "xcurrent = " << xcurrent << endl;
    // 	cout << "bulk.h   = " << bulk.h   << "\tenthalpy = " << enthalpy << endl;
    // 	cout << "state    = " << bulk.state << endl;
    // 	cout << "state   = " << bulk.state << "\t" << liq.state <<  "\t" << vap.state <<  "\t"  << salt.state << endl;
    // 	cout << "t       = " << bulk.t     << "\t" << liq.t     <<  "\t" << vap.t     <<  "\t"  << salt.t     << endl;
    // 	cout << "p       = " << bulk.p     << "\t" << liq.p     <<  "\t" << vap.p     <<  "\t"  << salt.p     << endl;
    // 	cout << "x       = " << bulk.x     << "\t" << liq.x     <<  "\t" << vap.x     <<  "\t"  << salt.x     << endl;
    // 	cout << "rho     = " << bulk.rho   << "\t" << liq.rho   <<  "\t" << vap.rho   <<  "\t"  << salt.rho   << endl;
    // 	cout << "h       = " << bulk.h     << "\t" << liq.h     <<  "\t" << vap.h     <<  "\t"  << salt.h     << endl;
    // 	cout << "cp      = " << bulk.cp    << "\t" << liq.cp    <<  "\t" << vap.cp    <<  "\t"  << salt.cp    << endl;
    // 	cout << "s       = " << bulk.s     << "\t" << liq.s     <<  "\t" << vap.s     <<  "\t"  << salt.s     << endl;
    // 	cout << "mf      = " << bulk.mf    << "\t" << liq.mf    <<  "\t" << vap.mf    <<  "\t"  << salt.mf    << endl;
	
    // 	//d DumpStatus();
    //   }
  }


  void H2ONaClFluidProperties::InitializeToBogus(){
    bulk.InitToBogus();    
    liq.InitToBogus();    
    vap.InitToBogus();    
    salt.InitToBogus();
  }

  // double H2ONaClFluidProperties::NewTwophaseCompressibility(){

  //     double Clapeyron;
  //     dldt  = (1.0-liq.mf)*vap.dhdt-liq.mf*liq.dhdt;
  //     dldt += -dliqmfdt*(vap.h+liq.h);
  //     dldp  = (1.0-liq.mf)*vap.dhdp-liq.mf*liq.dhdp;
  //     dldp += -dliqmfdp*(vap.h+liq.h);
  // //    Clapeyron = -dldt/(dldp*1.0e-5);
  //     Clapeyron  = (vap.h-liq.h)*dldt*liq.mf+vap.dhdt+liq.mf*dldt;
  //     Clapeyron /= (vap.h-liq.h)*dldp*1.0e-5*liq.mf+vap.dhdp*1.0e-5+liq.mf*dldp*1.0e-5;
  //     Clapeyron *= -1.0e0;
  //     cout.setf(ios::scientific);
  // //    cout << pcurrent << "\t" << -Clapeyron - (-dldt/(dldp*1.0e-5)) << endl;
  // //    cout << pcurrent << "\t" << (vap.h-liq.h)/((tcurrent+273.15)*(vap.v-liq.v)) << endl;
  //     double dVdT;
  //     dVdT  = dliqmfdt*(1.0/liq.rho-1.0/vap.rho);
  //     dVdT += liq.mf*(liq.alfa/liq.rho-vap.alfa/vap.rho);
  //     dVdT += vap.alfa/vap.rho;

  //     double dVdP;
  //     dVdP  = dliqmfdp*(1.0/liq.rho-1.0/vap.rho);
  //     dVdP += liq.mf*(-liq.beta/liq.rho+vap.beta/vap.rho);
  //     dVdP -= vap.beta/vap.rho;

  //     b  = (1.0-phi)*rr*cpr+phi*(liq.s*liq.rho*liq.cp+vap.s*vap.rho*vap.cp);
  // //    b  = (1.0-phi)*rr*cpr+phi*(bulk.dhdt);
  // //    b *= bulk.v*(bulk.alfa/bulk.dhdt-bulk.beta/bulk.dhdp);
  //     b *= bulk.v/bulk.dhdt*(bulk.alfa //-bulk.beta*1.0e-5*Clapeyron);
  //     b /= phi*Clapeyron;

  //     return b;
  // }

  void H2ONaClFluidProperties::PrintTwophaseCompressibilityParameters(){
    cout << pcurrent << "\t" << rl << "\t" << rv << "\t" << hl << "\t" << hv << "\t" << cpl << "\t" << cpv << "\t" << sl << "\t" << sv << endl;
    return;
  }

  double H2ONaClFluidProperties::TwophaseCompressibility()
  {
    // double t        = tcurrent;
    // double p        = water.SaturationPressureFromT(t);
    // double bl       = water.LiquidCompressibilityFromP(p);
    // double bv       = water.VaporCompressibilityFromP(p);
    // double rl       = water.LiquidDensityFromP(p);
    // double rv       = water.VaporDensityFromP(p);
    // double hl       = water.LiquidEnthalpyFromP(p);
    // double hv       = water.VaporEnthalpyFromP(p);
    // double cpl      = water.LiquidHeatCapacityFromP(p);
    // double cpv      = water.VaporHeatCapacityFromP(p);
    // double L        = hv-hl;
    // double dhldpsat = water.LiquidEnthalpyFromP(p+1.0e-5)-hl;
    // double dhvdpsat = water.VaporEnthalpyFromP(p+1.0e-5)-hv;
    // double mfv      = (enthalpy-hl)/(hv-hl);     
    // double mfl      = 1.0-mfv;
    // double sl       = mfl/rl / (mfl/rl+mfv/rv);
    // double sv       = 1.0-sl;
    
    // double dmfvdp   = (-dhldpsat-mfv*(dhvdpsat-dhldpsat))/L;
    // double b        = sl*bl + sv*bv;
    // b       += -1.0/(mfl/rl+mfv/rv)*(1.0/rv-1.0/rl)*dmfvdp;
    
    if(xcurrent == 0.0e0 && state == VL){
      rl  = liq.rho;
      rv  = vap.rho;
      hl  = liq.h;
      hv  = vap.h;
      cpl = liq.cp;
      cpv = vap.cp;
      sl  = liq.s;
      sv  = vap.s;
    }
    if(state == HM){
      rl  = salt.rho;
      rv  = liq.rho;
      hl  = salt.h;
      hv  = liq.h;
      cpl = salt.cp;
      cpv = liq.cp;
      sl  = salt.s;
      sv  = liq.s;
    }
    if(xcurrent > 0.0e0 && state == VL){
      rl  = liq.rho;
      rv  = vap.rho;
      hl  = liq.h;
      hv  = vap.h;
      cpl = liq.cp;
      cpv = vap.cp;
      sl  = liq.s;
      sv  = vap.s;
      //d PrintTwophaseCompressibilityParameters();
    }
    if(state == VLH){
      if(xcurrent <= liq.x){
        // should be approximately like VL but that is probably incorrect
        rl  = liq.rho;
        rv  = vap.rho;
        hl  = liq.h;
        hv  = vap.h;
        cpl = liq.cp;
        cpv = vap.cp;
        sl  = liq.s;
        sv  = vap.s;
      }
      else{
        // should be approximately like LH but that is probably incorrect
        rl  = salt.rho;
        rv  = liq.rho;
        hl  = salt.h;
        hv  = liq.h;
        cpl = salt.cp;
        cpv = liq.cp;
        sl  = salt.s;
        sv  = liq.s;
      }
    }
    if(state == VH){
      rl  = salt.rho;
      rv  = vap.rho;
      hl  = salt.h;
      hv  = vap.h;
      cpl = salt.cp;
      cpv = vap.cp;
      sl  = salt.s;
      sv  = vap.s;
    }
    if(state == LH){
      rl  = salt.rho;
      rv  = liq.rho;
      hl  = salt.h;
      hv  = liq.h;
      cpl = salt.cp;
      cpv = liq.cp;
      sl  = salt.s;
      sv  = liq.s;
    }

    product = ( rl - rv ) / ( ( hv - hl ) * rl * rv );
    b  = ( 1.0 - phi ) * rock.HeatCapacity(tcurrent) * rr;
    b += phi * rl * cpl * sl;
    b += phi * rv * cpv * sv;
    b *= product * product;
    b *= ( liq.t + 273.15 ) / phi;
    //    cout << pcurrent << "\t" << (hv-hl)/((tcurrent+273.15)*(1.0/rv-1.0/rl)) << endl;

    return b;

    // alternatively
    //    product = ( rl - rv ) / ( ( hv - hl ) * rl * rv );
    //    b  = ( 1.0 - phi ) * cpr * rr;
    //    b += phi * bulk.rho * bulk.dhdt; // what is bulk.dhdt for pure water or at xcrit(p const)?
    //    b += phi * bulk.rho * bulk.dhdt;
    //    b *= product * product;
    //    b *= ( tcurrent + 273.15 ) / phi;
    //    return b;
  }

  void H2ONaClFluidProperties::CheckEquilibrated()
  {
    // This does a final check if - after all other criteria are matched - also the
    // specific enthalpy is matched within a safety margin
    // Then, thermal equilibrium has been reached and equilibrated = true can
    // be issued to tell H2ONaClThermalEquilibrator that equilibration has converged
    if( std::fabs( (enthalpy-bulk.h)/enthalpy ) < 1.0e-6 )  equilibrated = true;
    else equilibrated = false;
  }


}//csmp

