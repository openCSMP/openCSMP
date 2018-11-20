//
//  equilibrateH2O_CO2_NaCl.cpp
//  CSMP_CO2GeoSequestrationSimulator
//
//  Created by Stephan Matthai on 21/10/18.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "equilibrateH2O_CO2_NaCl.h"
#include "EOS_CO2H2ONaCl_Spycher2004.h"
#include "HaliteLiquidus.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "ErrorHandler.h"

using namespace std;
using namespace csmp::variables;


#define CSMP_DEBUG_aqueousPhaseDensity_H2O_CO2_NaCl

namespace csmp {
  
EOS_CO2H2ONaCl_Spycher04  eos_csmp;


/**
   Pruess (2005), compute aqueous phase density from brine density and amount of CO2 dissolved
 
   @test OK SKM 6/11/20182
*/
double64 aqueousPhaseDensity_H2O_CO2_NaCl( double64 rho_brine, double64 rhoCO2, double64 X_CO2 )
 {
    return 1. / ((1.-X_CO2) / rho_brine + X_CO2 / rhoCO2);
   
 } // end aqueousPhaseDensity_H2O_CO2_NaCl


/**
    computes phase state of pure CO2, critical point 7.39 MPa, 31.1 oC
*/
THERMODYNAMIC_STATE stateOfCO2( double64 pCO2, double64 ToC)
 {
    // liquid
    if ( ToC<= 31.1 ) return CSMP_LIQUID;
    // vapour
    if ( pCO2 <= 7.39e6 ) return GASEOUS;
    // supercritical
    return SUPERCRITICAL;
 }

  
  
  
  
/**  equilibrateH2O_CO2_NaCl()

   Computes new compositions of aqueous and carbonic phase, precipitates salt if any, updates and saturations, densities and viscosities.
   Call this after each transport step. Then merely read (do not recompute) densities etc. when the fractional flow functions etc.
   are computed.
 
 Revised approach (roughly following Islam & Carlsson, 2012, Stanford geothermal workshop):
 
0. Read input compositions and PT data
 
1. Compute fugacity of (pure) CO2 phase and vapourisation pressure of aqueous phase

2. Compute solubilities of NaCl and CO2 within the aqueous phase

3. Compute solubility of water in the CO2 = carbonic phase

4. Compare equilibrium mass fractions of components in the phases for the new composition and PT conditions

   4.1 Dissolve whatever extra fits in
   4.2 Exsolve water into CO2 phase
   4.3 Evaporate CO2 from aqueous phase as the pressure decreases
   4.4 Precipitate salt if there is more salt than the saturation level
   4.5 Recompute phase densities and saturations

5. Compute aqueous diffusivity of CO2 (TODO: not stored at moment)
 
6. compute new phase compressibilities

7. Update the phase viscosities

8. Update phase compositions and store them

9. compute the fluid mass sources or sinks that are due to the compositional changes and the  expansion or compression of the fluid that occurred
 
We need to take into account the capillary pressure, i.e. the elevated pressure of the non-wetting phase, see Reichenberger paper.

@assumption input phase saturations sum up to 1.

@attention EOS functions take and return mass fraction percentages 0..100 -> in CSMP we store mass fractions, therefore we have conversions

@attention notation:  x denotes mole fractions and X mass fractions in aqueous phase; y,Y in carbonic phase

@attention this method updates the fluid densities. @todo We have to avoid the extra workload of doing this twice

@note SKM refactored 8/11/2018: do all computations at constant mass; then recompute densities and saturations and normalise to keep mass constant while changing volume; now compute source term
 
conversions:   vf1 = 1/(1+(rho1/rho2)*(1/mf1-1)) = -(rho2 * mf1) / (rho1 * mf1 - rho2 * mf1 - rho1)

@author SKM
@date 9/11/2018

@copyright CSMP at University of Melbourne.

*/
template<size_t dim>
void equilibrateH2O_CO2_NaCl( const variables::VariableSet_CO2GeoSequestration& props, Node<dim>& n, double64 del_t )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 0. Reading fluid composition after transport and computing mass balances
    // ------------------------------------------------------------------------
    ArrayVariable  H2Ocomposition(3), CO2composition(2);
    n.Read( props.key_H2O_comp, H2Ocomposition );
    n.Read( props.key_CO2_comp, CO2composition );
    
    if (fabs(H2Ocomposition(XH2O)+H2Ocomposition(XCO2)+H2Ocomposition(XNACl_aq)-1.) > numeric_limits<double64>::epsilon() )
      throw csmp::Exception( ERROR, "equilibrateH2O_CO2_NaCl:", "mass fractions of aqueous phase do not add up to 1.");

    if (fabs(CO2composition(YH2O)+CO2composition(YCO2)-1.) > numeric_limits<double64>::epsilon() )
      throw csmp::Exception( ERROR, "equilibrateH2O_CO2_NaCl:", "mass fractions of carbonic phase do not add up to 1."); 
         
      
    // saturation
    double64 old_sCO2 = n.Read( props.key_sCO2 );
    // phase densities (0) vs. new (1)
    double64 old_rhoH2O = n.Read( props.key_rhoH2O );
    double64 old_rhoCO2 = n.Read( props.key_rhoCO2 );
    const double64 rhoNaCl(2170.); // density = 2170 kg/m3
    // water
    const double64 old_sw(1. - old_sCO2);
    const double64 mass_H2O = H2Ocomposition(XH2O) * old_rhoH2O * old_sw  +
                              CO2composition(YH2O) * old_rhoCO2 * old_sCO2;
    // carbondioxide
    const double64 mass_CO2 = H2Ocomposition(XCO2) * old_rhoH2O * old_sw  +
                              CO2composition(YCO2) * old_rhoCO2 * old_sCO2;
    // dissolved salt
    double64 mass_NaClaq = H2Ocomposition(XNACl_aq) * old_rhoH2O * old_sw;
    // total mass
    double64 total_mass = mass_H2O + mass_CO2 + mass_NaClaq;
    // bulk density
    const double64 bulk_density(old_rhoH2O * old_sw + old_rhoCO2 * old_sCO2);
    // volume
    const double64 old_fluid_volume = total_mass / bulk_density;


    // 1. From P,T and pure CO2-phase molar volume, compute CO2 state and partial pressures of phases
    // ----------------------------------------------------------------------------------------------
    // TODO: local capillary pressure not taken into account yet correctly
    const double64 Pf(n.Read( props.key_pf ));
    const double64 ToC(n.Read( props.key_T )); // temperature in centrigrade
    // compute phase state
    //THERMODYNAMIC_STATE CO2_phase_state = stateOfCO2( pCO2, ToC);
    //Qi: use total pressure for now SKM - OK
    THERMODYNAMIC_STATE CO2_phase_state = stateOfCO2( Pf, ToC );
#ifdef DEBUG
   if ( CO2_phase_state == GASEOUS )
      cerr <<"\nequilibrateH2O_CO2_NaCl: the CO2 us in the gaseous state.";
#endif

    // 2. Compute NaCl and CO2 solubilities in the brine phase as mass fractions
    //    precipitating salt if brine is oversaturated,
    //    dissolving salt, if any, into undersaturated brine
    //    updating overall mass balance.
    // ----------------------------------------------------------------------------
    // NaCl dissolved in aqueous phase is converted into molality (notice stupid mass percentages)
    double64 molalityNaCl( eos_csmp.massFracNaClToMolalNaClInAqueousPhase(H2Ocomposition(XNACl_aq)*100.) );
    // solid salt if any ('key_NaCl' gives volume fraction of salt in pore space)
    double64  salt(n.Read(props.key_NaCl));
    // 'xNACl_aq' in H2Ocomposition, taking CO2 content of brine (tested against Driesner EOS, SKM OK)
//    double64 max_XNaCl_aq = eos_csmp.X_salt( Pf, ToC, molalityNaCl ) / 100.; // limit increases with moality to unplausible values

    HaliteLiquidus liquidus( ToC, Pf );
    double64 max_XNaCl_aq = liquidus.MassFractionNaCl(); // at 25C and reasonably low P should give ca.

    double64 delta_XNaCl_aq(H2Ocomposition(XNACl_aq) - max_XNaCl_aq);
    // precipitating excess salt if brine becomes supersaturated
    if ( delta_XNaCl_aq > 0. ) {
         // RESULT
         double64 kg_NaCl_precipitated = delta_XNaCl_aq * old_rhoH2O * old_sw;
         double64 V_NaCl_precipitated = kg_NaCl_precipitated / rhoNaCl;
         // storing salt precipitate (recorded as volume fraction occupied by solid in fluid phase)
         if ( n.Status(props.key_NaCl) != DIRICH )
           n.Store( props.key_NaCl, makeScalar(n.Status(props.key_NaCl),salt+V_NaCl_precipitated) );
         // updating mass balance
         mass_NaClaq -= kg_NaCl_precipitated;
         total_mass  -= kg_NaCl_precipitated;
      }
    // conversely, (re)dissolving solid salt if aqueous fluid was undersaturated
    else if ( salt > 0. ) { //                                 kg_NaCl_precipitated
         double64 V_salt_to_dissolve = std::min( salt, -(delta_XNaCl_aq * old_rhoH2O * old_sw) / rhoNaCl );
         // removing the dissolved solid salt
         if ( n.Status(props.key_NaCl) != DIRICH )
           n.Store( props.key_NaCl, makeScalar(n.Status(props.key_NaCl),salt-V_salt_to_dissolve) );
         // RESULT - conversion into mass
         double64 kg_salt_to_dissolve = V_salt_to_dissolve * rhoNaCl;
         // updating mass balance
         mass_NaClaq += kg_salt_to_dissolve;
         total_mass  += kg_salt_to_dissolve;
      }
    // updating NaCL molality, to account for the removal of salt
    molalityNaCl = eos_csmp.massFracNaClToMolalNaClInAqueousPhase( max_XNaCl_aq * 100. );
    // updating fluid salinity is done after the density has been recomputed, see below

    // maximum aqueous CO2 mass fraction that can be dissolved
    // fugacity is calculated inside of X_CO2 function
    double64 maxXCO2_aq = eos_csmp.X_Co2( Pf, ToC, molalityNaCl ) / 100.; // OK
   
   
    // 3. Compute solubility of H2O into CO2 phase; salt precipitation considered next time
    // ---------------------------------------------------------------------------------------------
    const double64 maxYH2O_carbonic = eos_csmp.Y_H2o( Pf, ToC, molalityNaCl ) / 100.; // OK
   

    // 4. Computing mass transfer between phases
    /* ---------------------------------------------------------------------------------------------------

        - CO2 saturation in water
        - H2O evaporated into CO2
        - new densities, including brine density with CO2 (Pruess approach)
        - new saturations
        - renormalised mass fractions
        - computing aqueous diffusivity of CO2 */

    // ---------------------------------------------------------------------------------------------------
#ifdef CSMP_DEBUG_aqueousPhaseDensity_H2O_CO2_NaCl
      double64 old_XCO2(H2Ocomposition(XCO2)), old_XH2O(H2Ocomposition(XH2O)), old_XNACl_aq(H2Ocomposition(XNACl_aq)),
               old_YCO2(CO2composition(YCO2)), old_YH2O(CO2composition(YH2O));
#endif
   
    // 4.1 Calculating adjustment of CO2 content of brine
    // ---------------------------------------------------------------------
    // comparing CO2-solubility in aqueous phase with current CO2 mass fraction of
    // and calculating how much extra CO2 could be dissolved in aqueous (or degree of oversaturation)
    // (positive value indicates that more CO2 can be dissolved)
    //                   max-massfraction   current mass fraction
    const double64 delta_XCO2aq = maxXCO2_aq - H2Ocomposition(XCO2);
    // imposing limit on dissolution due to restricted amount of CO2 left in carbonic phase
    const double64 mass_CO2_carbonic = CO2composition(YCO2) * old_rhoCO2 * old_sCO2; 
    // imposing limit on degasification due to restricted amount of CO2 left in aquesous phase 
    const double64 mass_CO2_aq = H2Ocomposition(XCO2) * old_rhoH2O * old_sw; 
    // RESULT - how much CO2 to dissolve in brine
    const double64 kg_CO2_to_dissolve = (delta_XCO2aq > 0. ) ? min( delta_XCO2aq * old_rhoH2O * old_sw, mass_CO2_carbonic) : max(delta_XCO2aq * old_rhoH2O * old_sw, -mass_CO2_aq);
    //double64 XCO2aq_new = H2Ocomposition(XCO2) + kg_CO2_to_dissolve / (old_rhoH2O * old_sw);
    double64 XCO2aq_new = H2Ocomposition(XCO2);
    if (old_rhoH2O * old_sw != 0.) XCO2aq_new  += kg_CO2_to_dissolve / (old_rhoH2O * old_sw);
    XCO2aq_new = max(XCO2aq_new, 0.);

    // 4.2 Calculating adjustment of water content of CO2 phase
    // ---------------------------------------------------------------------
    // calculating how much extra H2O could be dissolved (or degree of oversaturation)
    // (positive value indicates that more H2O can be dissolved)
    //                            max-massfraction   current mass fraction
    const double64 delta_YH2O_carbonic = (maxYH2O_carbonic - CO2composition(YH2O));
    // imposing limit on condensation due to restricted amount of water left in carbonic phase
    const double64 mass_H2O_carbonic = CO2composition(YH2O) * old_rhoCO2 * old_sCO2; 
    // imposing limit on evaporation due to restricted amount of water left in aqueous phase
    const double64 mass_H2O_aq = H2Ocomposition(XH2O) * old_rhoH2O * old_sw;
    // RESULT - amount of water to evaporate into CO2 phase
    const double64 kg_H2O_to_transfer = (delta_YH2O_carbonic > 0. ) ? min( delta_YH2O_carbonic * old_rhoCO2 * old_sCO2, mass_H2O_aq ) : max(delta_YH2O_carbonic * old_rhoCO2 * old_sCO2, -mass_H2O_carbonic);
 

    // 4.3 Updating fluid densities: CO2 and brine phase
    // ---------------------------------------------------------------------
    // RESULT
    const double64 rhoCO2    = eos_csmp.Rho_CarbonicPhase( Pf, ToC);
    n.Store( props.key_rhoCO2, makeScalar(n.Status(props.key_rhoCO2),rhoCO2) );
    const double64 rho_brine = eos_csmp.Rho_brine(Pf, ToC, molalityNaCl );
    // RESULT - harmonic mean of the densities of brine and CO2
    double64 rho_brine_CO2   = aqueousPhaseDensity_H2O_CO2_NaCl( rho_brine, rhoCO2, XCO2aq_new );
rho_brine_CO2 = std::max( rho_brine_CO2, rho_brine );
    assert( rho_brine_CO2 >= rho_brine ); // Pruess, 2005
    // saving
    if ( n.Status(props.key_rhoH2O) != DIRICH )
      n.Store( props.key_rhoH2O, makeScalar(n.Status(props.key_rhoH2O),rho_brine_CO2) );


    // 4.4 Updating saturations using the new mass fractions (vf=volume fraction, mf=mass fraction)
    // --------------------------------------------------------------------------------------------
    const double64 mf_Y_phase = (old_rhoCO2 * old_sCO2 - kg_CO2_to_dissolve + kg_H2O_to_transfer) / total_mass;
    //const double64 mf_X_phase = (old_rhoH2O * old_sH2O + kg_CO2_to_dissolve - kg_H2O_to_transfer) / total_mass;
    //  s1 = vf1=1/(1+(rho1/rho2)*(1/mf1-1)) = -(rho2 * mf1) / (rho1 * mf1 - rho2 * mf1 - rho1)
    double64 sCO2new = -(rho_brine_CO2 * mf_Y_phase) / (rhoCO2 * mf_Y_phase - rho_brine_CO2 * mf_Y_phase - rhoCO2);
    assert( sCO2new >= 0. );
    assert( sCO2new < 1. + numeric_limits<double64>::epsilon() );
    sCO2new = std::max( sCO2new, 0. );
    sCO2new = std::min( sCO2new, 1. );
    // RESULT - assuming that saturations add up to 1 and the brine phase is less compressible than the CO2 phase
    if ( n.Status(props.key_sCO2) != DIRICH )
      n.Store( props.key_sCO2, makeScalar(n.Status(props.key_sCO2),sCO2new) );
    if ( n.Status(props.key_sH2O) != DIRICH )
      n.Store( props.key_sH2O, makeScalar(n.Status(props.key_sH2O),1. - sCO2new) );
   

    // 4.5 Computing and storing new compositions
    // ------------------------------------------------------------------------------------
    double64 sH2Onew = 1. - sCO2new;
    //const double64 new_mass_aqueous_phase(rho_brine_CO2 * sH2Onew);
    const double64 new_mass_aqueous_phase = old_rhoH2O * old_sw + kg_CO2_to_dissolve - kg_H2O_to_transfer;
    // aqueous phase
    if(new_mass_aqueous_phase > numeric_limits<double64>::epsilon()) { 
        H2Ocomposition(XCO2) = ( H2Ocomposition(XCO2) * old_rhoH2O * old_sw + kg_CO2_to_dissolve) / new_mass_aqueous_phase;
        H2Ocomposition(XCO2)     = std::min( H2Ocomposition(XCO2), maxXCO2_aq );
        H2Ocomposition(XCO2)     = std::max( H2Ocomposition(XCO2), 0. );
        //H2Ocomposition(XH2O)    -= kg_H2O_to_transfer / new_mass_aqueous_phase;
        H2Ocomposition(XH2O) = ( H2Ocomposition(XH2O) * old_rhoH2O * old_sw - kg_H2O_to_transfer) / new_mass_aqueous_phase;
        H2Ocomposition(XH2O)     = std::max( H2Ocomposition(XH2O), 0. );
        H2Ocomposition(XNACl_aq) = mass_NaClaq / new_mass_aqueous_phase;
        H2Ocomposition(XNACl_aq) = std::max( H2Ocomposition(XNACl_aq), 0. );
      
        // renormalisation
        H2Ocomposition /= (H2Ocomposition(XH2O) + H2Ocomposition(XCO2) + H2Ocomposition(XNACl_aq));
        assert( fabs((H2Ocomposition(XH2O) + H2Ocomposition(XCO2) + H2Ocomposition(XNACl_aq))-1.) <= numeric_limits<double64>::epsilon() );
    } else {
        H2Ocomposition(XH2O) = 1.;
        H2Ocomposition(XCO2) = 0.;
        H2Ocomposition(XNACl_aq) = 0.;
    }
    // carbonic phase: OK
    //const double64 new_mass_carbonic_phase(rhoCO2 * sCO2new);
    const double64 new_mass_carbonic_phase = old_rhoCO2 * old_sCO2 - kg_CO2_to_dissolve + kg_H2O_to_transfer;
    if(new_mass_carbonic_phase > numeric_limits<double64>::epsilon()) { 
        CO2composition(YH2O) = ( CO2composition(YH2O) * old_rhoCO2 * old_sCO2 + kg_H2O_to_transfer) / new_mass_carbonic_phase;
        CO2composition(YH2O)  = std::min( CO2composition(YH2O), maxYH2O_carbonic );
        CO2composition(YH2O)  = std::max( CO2composition(YH2O), 0. );
        //CO2composition(YCO2) -= kg_CO2_to_dissolve / new_mass_carbonic_phase;
        CO2composition(YCO2) = ( CO2composition(YCO2) * old_rhoCO2 * old_sCO2 - kg_CO2_to_dissolve) / new_mass_carbonic_phase;
        CO2composition(YCO2)  = std::max( CO2composition(YCO2), 0. );
        // renormalisation
        CO2composition /= (CO2composition(YH2O) + CO2composition(YCO2));
        assert( fabs((CO2composition(YH2O) + CO2composition(YCO2))-1.) <= numeric_limits<double64>::epsilon() );

    } else {
        CO2composition(YCO2) = 1.;
        CO2composition(YH2O) = 0.;
    }

    // storing compositions
    if ( n.Status( props.key_CO2_comp ) != DIRICH ) n.Store( props.key_CO2_comp, CO2composition );
    if ( n.Status( props.key_H2O_comp ) != DIRICH ) n.Store( props.key_H2O_comp, H2Ocomposition );
    // fluid salinity
    
    // compute and store dissolved CO2, evaporated water and salinity for visualisation
    // dissolved CO2 (kg/m3)
    if ( n.Status(props.key_CO2aq) != DIRICH ) {
      //const double64 dissolved_CO2 = H2Ocomposition(XCO2) * rho_brine_CO2 * sH2Onew; 
      const double64 dissolved_CO2 = H2Ocomposition(XCO2);
      n.Store( props.key_CO2aq, makeScalar(n.Status(props.key_CO2aq),dissolved_CO2) );
    }     
    // salinity (kg/m3)
    if ( n.Status(props.key_NaClaq) != DIRICH ) {
         // mass salt within a cubic meter of brine (kg/m3)
         double64 salinity = H2Ocomposition(XNACl_aq) * rho_brine_CO2;
         n.Store( props.key_NaClaq, makeScalar(n.Status(props.key_NaClaq),salinity) );
      }
   
    // evaporated water (kg/m3)
    if ( n.Status(props.key_H2Og) != DIRICH ) {
      //const double64 evaporated_water = CO2composition(YH2O) * rhoCO2 * sCO2new; 
      const double64 evaporated_water = CO2composition(YH2O); 
      n.Store( props.key_H2Og, makeScalar(n.Status(props.key_H2Og),evaporated_water) );
    }      
    
      
   
    // 5. computing diffusivity of CO2 in water at the given P, T, salinity
    // ---------------------------------------------------------------------
    const double64 D_Co2 = eos_csmp.D_Co2( Pf, ToC, molalityNaCl );
    if ( D_Co2 < 2.0e-9 || D_Co2 > 12.5e-9 ) {
         cerr <<"\n\t"<< D_Co2 <<" m2/s.";
         csmp_error.notice( WARNING, "equilibrateH2O_CO2_NaCl:", "CO2 diffusivity out of bounds defined by Cardogan et al. 2014.");
      }
    // store value
    // TODO: add aqueous diffusivity of CO2 to variable list and consider it in transport calculation
    // if ( n.Status(props.key_DCO2aq) != DIRICH )
    //  n.Store( props.key_DCO2aq, makeScalar(n.Status(props.key_DCO2aq),D_Co2) );   // nodal CO2 diffusivity in brine
   

    // 6. Viscosities
    // ------------------------------------------------------------------------------------
    const double64 muH2O = eos_csmp.mu_AqueousPhase( Pf, ToC, molalityNaCl );
    const double64 muCO2 = eos_csmp.mu_CarbonicPhase( Pf, ToC);
    n.Store( props.key_muH2O, makeScalar(n.Status(props.key_muH2O),muH2O) );
    n.Store( props.key_muCO2, makeScalar(n.Status(props.key_muCO2),muCO2) );



    // 7. Phase compressibilities
    // ------------------------------------------------------------------------------------
    const double64 betaH2O = eos_csmp.C_AqueousPhase( Pf, ToC, molalityNaCl );
    const double64 betaCO2 = eos_csmp.C_CarbonicPhase( Pf, ToC);
    n.Store( props.key_cH2O, makeScalar(n.Status(props.key_cH2O),betaH2O) );
    n.Store( props.key_cCO2, makeScalar(n.Status(props.key_cCO2),betaCO2) );



    // 8. Computing mass sources or sinks
    // ----------------------------------
    // do the mass fractions add up?
assert( fabs(sH2Onew - sCO2new) < 1. + numeric_limits<double64>::epsilon() );

    // water
    const double64 new_mass_H2O = H2Ocomposition(XH2O) * rho_brine_CO2 * sH2Onew  +
                                  CO2composition(YH2O) * rhoCO2 * sCO2new;
    // carbondioxide
    const double64 new_mass_CO2 = H2Ocomposition(XCO2) * rho_brine_CO2 * sH2Onew  +
                                  CO2composition(YCO2) * rhoCO2 * sCO2new;
    // dissolved salt
    const double64 new_mass_NaClaq = H2Ocomposition(XNACl_aq) * rho_brine_CO2 * sH2Onew;
    // total mass
    const double64 new_total_mass = new_mass_H2O + new_mass_CO2 + new_mass_NaClaq;
    // bulk density
    const double64 new_bulk_density(rho_brine_CO2 * sH2Onew + rhoCO2 * sCO2new);
    // new volume (mass/density)
    const double64 new_fluid_volume = new_total_mass / new_bulk_density;

    // RESULT - computing 'nodal fluid volume source' from apparent change in fluid masss
    double64  nodal_mass_source = new_total_mass - total_mass;
    if ( n.Status(props.key_nQM) != DIRICH )
      n.Store( props.key_nQM, makeScalar(n.Status(props.key_nQM),nodal_mass_source) );

    // RESULT - computing 'nodal fluid volume source' from apparent change in fluid masss
    const double64  expansion_factor  = 1. - (nodal_mass_source / new_bulk_density) / new_fluid_volume;
    
    // store nodal fluid mass source for pressure computation
    const double64 PV = n.Read(props.key_fvPV);
    
    nodal_mass_source *= (PV/del_t);
    if (n.Status(props.key_nQM) != DIRICH)
        n.Store( props.key_nQM, makeScalar(n.Status(props.key_nQM),nodal_mass_source) );
    
    // store mass
    n.Store( props.key_mH2O, makeScalar(n.Status(props.key_mH2O),new_mass_H2O*PV) );
    n.Store( props.key_mCO2, makeScalar(n.Status(props.key_mCO2),new_mass_CO2*PV) );
    n.Store( props.key_mNaClaq, makeScalar(n.Status(props.key_mNaClaq),new_mass_NaClaq*PV) );
    const double64 new_V_salt(n.Read(props.key_NaCl));
    const double64 new_mass_NaClsd = new_V_salt *PV *rhoNaCl;
    n.Store( props.key_mNaClsd, makeScalar(n.Status(props.key_mNaClsd),new_mass_NaClsd) );    

    // reporting mass fractions
#ifdef CSMP_DEBUG_aqueousPhaseDensity_H2O_CO2_NaCl
      cerr <<"\n\nNode "<< n.Idx() <<": y-coordinate (datum): "<< n.y() <<" m, pf: "<< Pf <<" Pa, T: "<< ToC <<" oC, NaCl molality: "<< molalityNaCl;
      cerr <<"\n\ninitial sum of mass fractions: ";
      cerr <<"\n\tH2O: "<< old_XCO2+old_XH2O+old_XNACl_aq;
      cerr <<"\n\tCO2: "<< old_YCO2+old_YH2O;
      cerr <<"\n";
      cerr <<"\nequilibration calculation: phase properties (BEFORE) and AFTER calculation: ";
      cerr <<"\n\tCO2 saturation:         ("<< old_sCO2     <<") "<< sCO2new;
      cerr <<"\n\tXCO2:                   ("<< old_XCO2     <<") "<< H2Ocomposition(XCO2);
      cerr <<"\n\tXH2O:                   ("<< old_XH2O     <<") "<< H2Ocomposition(XH2O);
      cerr <<"\n\tXNaCl:                  ("<< old_XNACl_aq <<") "<< H2Ocomposition(XNACl_aq);
      cerr <<"\n\tYCO2:                   ("<< old_YCO2     <<") "<< CO2composition(YCO2);
      cerr <<"\n\tYH2O:                   ("<< old_YH2O     <<") "<< CO2composition(YH2O);
      cerr <<"\n\tnew CO2 density:        ("<< old_rhoCO2   <<") "<< rhoCO2 <<" kg/m3";
      cerr <<"\n\tnew brine density:      ("<< old_rhoH2O   <<") "<< rho_brine_CO2 <<" kg/m3";
      cerr <<"\n\tnew salinity:           "<< n.Read(props.key_NaClaq) <<" kg/m3";
      cerr <<"\n\tnew salt:               "<< n.Read(props.key_NaCl) <<" X/m3";
      cerr <<"\n\tCO2 diffusivity in H2O: "<< D_Co2 * 1.0e9 <<" 10^-9 m2/s";
      cerr <<"\n\tCO2 dynamic viscosity:  "<< muCO2 <<" Pa.s";
      cerr <<"\n\tH2O dynamic viscosity:  "<< muH2O <<" Pa.s";
      cerr <<"\n\tfluid expansion caused: "<< expansion_factor <<" (fraction of original volume.)";
      const double64 dV = (old_fluid_volume - new_fluid_volume) / old_fluid_volume;
      cerr <<"\n\nfluid volume mismatch error (divergence): "<< dV;
      cerr <<"\n\n";
#endif

} // end equilibrateH2O_CO2_NaCl


template void equilibrateH2O_CO2_NaCl<1U>( const variables::VariableSet_CO2GeoSequestration&, Node<1U>&, double64 );
template void equilibrateH2O_CO2_NaCl<2U>( const variables::VariableSet_CO2GeoSequestration&, Node<2U>&, double64 );
template void equilibrateH2O_CO2_NaCl<3U>( const variables::VariableSet_CO2GeoSequestration&, Node<3U>&, double64 );



/**
    Applies  equilibrateH2O_CO2_NaCl() method to nodes of region
*/
template<size_t dim>
void equilibrateFluid( const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain, double64 delta_t )
 {
    for ( auto nit=model_subdomain.NodesBegin(); nit!=model_subdomain.NodesEnd(); ++nit )
      equilibrateH2O_CO2_NaCl( props, *(*nit), delta_t );
 }

template void equilibrateFluid<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>&, double64 );
template void equilibrateFluid<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>&, double64 );
template void equilibrateFluid<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>&, double64 );




/**
    Computes (barycentric) porosity from volume fraction of salt that is occupying the pore space.
    The new value is stored.
 
    @attention permeability change has to be calculated separately.
*/
template<size_t dim>
double64 porosityWithSalt( const variables::VariableSet_CO2GeoSequestration& props, Element<dim>& e )
 {
    // getting the salt volume fractions from the nodes
    double64 salt_volume_fraction = e.PropertyValueAtBaryCenter( props.key_NaCl );
    double64 phi = e.Read( props.key_phi );
   
    // changing the porosity (preventing it from going negative)
    assert( salt_volume_fraction >= 0. );
    phi = max( 0., phi - salt_volume_fraction );
   
    e.Store( props.key_phi, makeScalar( e.Status(props.key_phi), phi ) );
   
    return phi;
   
 } // end porosityAccoutingForSalt

template double64 porosityWithSalt<1U>( const variables::VariableSet_CO2GeoSequestration&, Element<1U>& );
template double64 porosityWithSalt<2U>( const variables::VariableSet_CO2GeoSequestration&, Element<2U>& );
template double64 porosityWithSalt<3U>( const variables::VariableSet_CO2GeoSequestration&, Element<3U>& );



/**
    Applies  porosityWithSalt() method to all model nodes
*/
template<size_t dim>
void updatePorosityAndPermeability( const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain )
 {
    for ( auto eit=model_subdomain.ElementsBegin(); eit!=model_subdomain.ElementsEnd(); ++eit ) {
        double64 phi_old = (*eit)->Read(props.key_phi);
        porosityWithSalt( props, *(*eit) );
        double64 phi_new = (*eit)->Read(props.key_phi);
            
        if (phi_new != phi_old && phi_old != 0. && (1.-phi_new) != 0.) {
            double64 k_old = (*eit)->Read(props.key_k);
            double64 phi3 = (phi_new*phi_new*phi_new)/(phi_old*phi_old*phi_old);
            double64 t1 = (1.-phi_old)*(1.-phi_old)/(1.-phi_new)/(1.-phi_new);
            double64 k_new = k_old * phi3 * t1;
            (*eit)->Store( props.key_k, makeScalar( (*eit)->Status( props.key_k), k_new ) );
        }
    }
 }
 
template void updatePorosityAndPermeability<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>& ); 
template void updatePorosityAndPermeability<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>& );
template void updatePorosityAndPermeability<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>& );
 


/**
    update pore volume to nodes of region
*/
template<size_t dim>
void updatePoreVolume( const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain )
 {
    for ( auto nit=model_subdomain.NodesBegin(); nit!=model_subdomain.NodesEnd(); ++nit ) {
        double64 pore_volume (0.); 
        for ( size_t t=0U; t<(*nit)->Parents(); t++ )
        {
            Element<dim>* const eptr((*nit)->Parent(t));
            double64 phi = eptr->Read( props.key_phi);
            const double64 thickness = eptr->Read( props.key_thi );
            if (!isnan(thickness)) phi *= thickness; //if thickness is initialised
            const size_t pnid((*nit)->ParentNodeNumber(t));
            const double64 sector_volume = eptr->SectorVolume(pnid);
            pore_volume   += phi * sector_volume;
        }
        (*nit)->Store( props.key_fvPV, makeScalar( (*nit)->Status( props.key_fvPV), pore_volume ) ); 
    }
 }

template void updatePoreVolume<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>& );
template void updatePoreVolume<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>& );
template void updatePoreVolume<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>& );



/**
    update PTX based fluid properties (density, viscosity, compressibility) for all model nodes
*/
template<size_t dim>
void updatePTXBasedFluidProperties(const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain )
{
    const double64 rhoNaCl(2170.); // density = 2170 kg/m3
    for ( auto nit=model_subdomain.NodesBegin(); nit!=model_subdomain.NodesEnd(); ++nit )
    {    
        ArrayVariable  H2Ocomposition(3), CO2composition(2);
        (*nit)->Read( props.key_H2O_comp, H2Ocomposition ); 
        (*nit)->Read( props.key_CO2_comp, CO2composition ); 
        const double64 Pf((*nit)->Read( props.key_pf ));
        const double64 TC((*nit)->Read( props.key_T )); // in centrigrade
        const double64 molalityNaCl( eos_csmp.massFracNaClToMolalNaClInAqueousPhase(H2Ocomposition(XNACl_aq)* 100.) );
        
        const double64 rhoCO2 = eos_csmp.Rho_CarbonicPhase( Pf, TC );
        if ( (*nit)->Status(props.key_rhoCO2) != DIRICH )
            (*nit)->Store( props.key_rhoCO2, makeScalar((*nit)->Status(props.key_rhoCO2),rhoCO2) );
            
        const double64 rho_brine = eos_csmp.Rho_brine(Pf, TC, molalityNaCl);
        double64 rho_brine_CO2 = rho_brine;
        if ( (*nit)->Status(props.key_rhoH2O) != DIRICH )
           (*nit)->Store( props.key_rhoH2O, makeScalar((*nit)->Status(props.key_rhoH2O),rho_brine_CO2) );
           
        const double64 muH2O = eos_csmp.mu_AqueousPhase( Pf, TC, molalityNaCl );
        const double64 muCO2 = eos_csmp.mu_CarbonicPhase( Pf, TC );
        if ( (*nit)->Status(props.key_muH2O) != DIRICH )
            (*nit)->Store( props.key_muH2O, makeScalar((*nit)->Status(props.key_muH2O),muH2O) );
        if ( (*nit)->Status(props.key_muCO2) != DIRICH )
            (*nit)->Store( props.key_muCO2, makeScalar((*nit)->Status(props.key_muCO2),muCO2) );

        const double64 betaH2O = eos_csmp.C_AqueousPhase( Pf, TC, molalityNaCl );
        const double64 betaCO2 = eos_csmp.C_CarbonicPhase( Pf, TC );
        if ( (*nit)->Status(props.key_cH2O) != DIRICH )
            (*nit)->Store( props.key_cH2O, makeScalar((*nit)->Status(props.key_cH2O),betaH2O) );
        if ( (*nit)->Status(props.key_cCO2) != DIRICH )
            (*nit)->Store( props.key_cCO2, makeScalar((*nit)->Status(props.key_cCO2),betaCO2) );

        const double64 sCO2new = (*nit)->Read(props.key_sCO2);
        const double64 sH2Onew = (*nit)->Read(props.key_sH2O);
        // water mass
        const double64 new_mass_H2O = H2Ocomposition(XH2O) * rho_brine * sH2Onew  +
                                      CO2composition(YH2O) * rhoCO2 * sCO2new;
        // carbondioxide mass
        const double64 new_mass_CO2 = H2Ocomposition(XCO2) * rho_brine * sH2Onew  +
                                      CO2composition(YCO2) * rhoCO2 * sCO2new;
        // dissolved salt mass
        const double64 new_mass_NaClaq = H2Ocomposition(XNACl_aq) * rho_brine * sH2Onew;
        const double64 PV = (*nit)->Read(props.key_fvPV);
        
        // solid/precipitated salt mass
        const double64 new_V_salt((*nit)->Read(props.key_NaCl));
        const double64 new_mass_NaClsd = new_V_salt *PV *rhoNaCl;       
    
        // store mass
        (*nit)->Store( props.key_mH2O, makeScalar((*nit)->Status(props.key_mH2O),new_mass_H2O*PV) );
        (*nit)->Store( props.key_mCO2, makeScalar((*nit)->Status(props.key_mCO2),new_mass_CO2*PV) );
        (*nit)->Store( props.key_mNaClaq, makeScalar((*nit)->Status(props.key_mNaClaq),new_mass_NaClaq*PV) );
        (*nit)->Store( props.key_mNaClsd, makeScalar((*nit)->Status(props.key_mNaClsd),new_mass_NaClsd) );
            
                       
    }
} 

template void updatePTXBasedFluidProperties<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>& ); 
template void updatePTXBasedFluidProperties<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>& );
template void updatePTXBasedFluidProperties<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>& );
 
 
} // end csmp


