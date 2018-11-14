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

using namespace std;
using namespace csmp::variables;

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
THERMODYNAMIC_STATE stateOfCO2( double64 pCO2, double64 TC )
 {
    // liquid
    if ( TC <= 31.1 ) return LIQUID;
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
void equilibrateH2O_CO2_NaCl( const variables::VariableSet_CO2GeoSequestration& props, Node<dim>& n )
 {
    const bool verbose(true);
   
    // 0. Reading fluid composition after transport and computing mass balances
    // ------------------------------------------------------------------------
    ArrayVariable  H2Ocomposition(3), CO2composition(2);
    n.Read( props.key_H2O_comp, H2Ocomposition );
    n.Read( props.key_CO2_comp, CO2composition );
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
    const double64 entry_pressure(0.); // TODO: not a node property: n.Read(props.key_pd)); // of non-wetting phase
    const double64 TC(n.Read( props.key_T )); // in centrigrade
    double64 phaseVolumeCo2_m3_per_mol = eos_csmp.CompressedVolumeCo2( Pf + entry_pressure, TC ); // TODO: check whether correct?
    //                   returns the fugacity coefficient
    double64 pCO2 = Pf * eos_csmp.FugacityCo2( Pf, TC, phaseVolumeCo2_m3_per_mol );
    assert( pCO2 <= Pf );
    // water
    double64 pH2O = Pf * eos_csmp.FugacityH2o( Pf, TC, phaseVolumeCo2_m3_per_mol );
    assert( pH2O <= Pf );
    // compute phase state
    THERMODYNAMIC_STATE CO2_phase_state = stateOfCO2( pCO2, TC );
    if ( CO2_phase_state == GASEOUS )
      throw csmp::Exception( WARNING, "equilibrateH2O_CO2_NaCl", "the CO2 us in the gaseous state.");
   
   
    // 2. Compute NaCl and CO2 solubilities in the brine phase as mass fractions
    //    precipitating salt if brine is oversaturated,
    //    dissolving salt, if any, into undersaturated brine
    //    updating overall mass balance.
    // ----------------------------------------------------------------------------
    // solid salt if any ('key_NaCl' gives volume fraction of salt in pore space)
    double64  salt(n.Read(props.key_NaCl));
    // 'xNACl_aq' inside of H2Ocomposition = advected quantity
    // H2O-NaCl equation of Thomas Drieser (06, Geochim. CosmoChim. Acta)
    HaliteLiquidus liquidus( TC, Pf ); // TODO: should pH2O be used instead of Pf?
    double64 max_XNaCl_aq = liquidus.MassFractionNaCl(); // at 25C and reasonably low P should give ca. 0.26
    double64 delta_XNaCl_aq(H2Ocomposition(XNACl_aq) - max_XNaCl_aq);
    // precipitating excess salt if brine becomes supersaturated
    if ( delta_XNaCl_aq > 0. ) {
         // RESULT
         const double64 kg_NaCl_precipitated = delta_XNaCl_aq * old_rhoH2O * old_sw;
         // converting this mass into a volume fraction of the space filled up by the fluid
         const double V_NaCl = (kg_NaCl_precipitated * rhoNaCl + salt) / (total_mass * bulk_density);
         // storing salt precipitate (recorded as volume fraction occupied by solid in fluid phase)
         if ( n.Status(props.key_NaCl) != DIRICH ) n.Store( props.key_NaCl, makeScalar(n.Status(props.key_NaCl),V_NaCl) );
         mass_NaClaq -= kg_NaCl_precipitated;
         total_mass  -= kg_NaCl_precipitated;
      }
    // conversely, (re)dissolving solid salt if aqueous fluid was undersaturated
    else if ( salt > 0. ) {
         const double64 V_salt_to_dissolve = std::min( salt, -delta_XNaCl_aq * old_rhoH2O * old_sw );
         if ( n.Status(props.key_NaCl) != DIRICH ) n.Store( props.key_NaCl, makeScalar(PLAIN,salt-V_salt_to_dissolve) );
         // RESULT - conversion into mass
         const double64 kg_salt_to_dissolve = V_salt_to_dissolve * rhoNaCl;
         mass_NaClaq += kg_salt_to_dissolve;
         total_mass  += kg_salt_to_dissolve;
      }

    // NaCl dissolved in aqueous phase is converted into molality (notice stupid mass percentages)
    const double64 molalityNaCl( eos_csmp.massFracNaClToMolalNaClInAqueousPhase(H2Ocomposition(XNACl_aq)) * 100. ); // OK
    // maximum aqueous CO2 mass fraction that can be dissolved
    double64 maxXCO2_aq = eos_csmp.X_Co2( pH2O, TC, molalityNaCl ) / 100.; // OK
   
   
    // 3. Compute solubility of H2O into CO2 phase; salt precipitation considered next time
    // ---------------------------------------------------------------------------------------------
    const double64 maxYH2O_carbonic = eos_csmp.Y_H2o( pH2O, TC, molalityNaCl ) / 100.; // OK
   

    // 4. Computing mass transfer between phases
    /* ---------------------------------------------------------------------------------------------------

        - CO2 saturation in water
        - H2O evaporated into CO2
        - new densities, including brine density with CO2 (Pruess approach)
        - new saturations
        - renormalised mass fractions
        - computing aqueous diffusivity of CO2 */

    // ---------------------------------------------------------------------------------------------------
    if ( verbose ) {
      // looking at the current fluid composition
      cerr <<"\nNode "<< n.Idx() <<": y-coordinate (datum): "<< n.y() <<" m, pf: "<< Pf <<" Pa, T: "<< TC <<" oC, NaCl molality: "<< molalityNaCl;
      cerr <<"\n\nfugacities CO2, H2O: "<< eos_csmp.FugacityCo2( Pf, TC, phaseVolumeCo2_m3_per_mol ) <<", "<< eos_csmp.FugacityH2o( Pf, TC, phaseVolumeCo2_m3_per_mol );
      cerr <<"\n\tCO2 saturation: "<< old_sCO2;
      cerr <<"\n\tdensity CO2:    "<< old_rhoCO2;
      cerr <<"\n\tdensity H2O:    "<< old_rhoH2O;
      cerr <<"\n\nmass fractions (X=aqueous, Y=carbonic phase) after transport (BEFORE re-equilibration): ";
      cerr <<"\n\tXCO2:  "<< H2Ocomposition(XCO2);
      cerr <<"\n\tXH2O:  "<< H2Ocomposition(XH2O);
      cerr <<"\n\tXNaCl: "<< H2Ocomposition(XNACl_aq);
      cerr <<"\n\tYCO2:  "<< CO2composition(YCO2);
      cerr <<"\n\tYH2O:  "<< CO2composition(YH2O);
      cerr <<"\n";
    }
   
    // 4.1 Calculating adjustment of CO2 content of brine
    // ---------------------------------------------------------------------
    // comparing CO2-solubility in aqueous phase with current CO2 mass fraction of
    // and calculating how much extra CO2 could be dissolved in aqueous (or degree of oversaturation)
    // (positive value indicates that more CO2 can be dissolved)
    //                   max-massfraction   current mass fraction
    const double64 delta_XCO2aq = maxXCO2_aq - H2Ocomposition(XCO2);
    // RESULT - how much CO2 to dissolve in brine
    const double64 kg_CO2_to_dissolve = std::min( CO2composition(YCO2) * old_rhoCO2 * old_sCO2, delta_XCO2aq * old_rhoCO2 * old_sCO2 );
    const double64 XCO2aq_new = H2Ocomposition(XCO2) - kg_CO2_to_dissolve / (old_rhoH2O * old_sw);

    // 4.2 Calculating adjustment of water content of CO2 phase
    // ---------------------------------------------------------------------
    // calculating how much extra H2O could be dissolved (or degree of oversaturation)
    // (positive value indicates that more H2O can be dissolved)
    //                            max-massfraction   current mass fraction
    const double64 delta_YH2O_carbonic = (maxYH2O_carbonic - CO2composition(YH2O));
    // imposing limit on evaporation due to restricted amount of water that can be evaoporated into the CO2 phase
    const double64 mass_CO2_carbonic = CO2composition(YCO2) * old_rhoCO2 * old_sCO2;
    const double64 mass_H2O_carbonic = CO2composition(YH2O) * old_rhoH2O * old_sCO2;
    // imposing limit on evaporation due to restricted amount of water left
    const double64 mass_H2O_aq = H2Ocomposition(XH2O) * old_rhoH2O * old_sw;
    // RESULT - amount of water to evaporate into CO2 phase
    const double64 kg_H2O_to_transfer = (delta_YH2O_carbonic > 0. ) ? min( delta_YH2O_carbonic * mass_CO2_carbonic, mass_H2O_aq ) : mass_H2O_carbonic;
 

    // 4.3 Updating fluid densities: CO2 and brine phase
    // ---------------------------------------------------------------------
    // RESULT
    const double64 rhoCO2    = eos_csmp.Rho_CarbonicPhase( Pf, TC );
    n.Store( props.key_rhoCO2, makeScalar(n.Status(props.key_rhoCO2),rhoCO2) );
    const double64 rho_brine = eos_csmp.Rho_brine(Pf, TC, molalityNaCl);
    // RESULT
    double64 rho_brine_CO2   = aqueousPhaseDensity_H2O_CO2_NaCl( rho_brine, rhoCO2, XCO2aq_new );
// TODO: why is density greater?
rho_brine_CO2 = std::max( rho_brine_CO2, rho_brine );
    assert( rho_brine_CO2 >= rho_brine ); // Pruess, 2005
    // saving
    if ( n.Status(props.key_rhoH2O) != DIRICH )
      n.Store( props.key_rhoH2O, makeScalar(n.Status(props.key_rhoH2O),rho_brine_CO2) );
    // RESULT - salinity (kg NaCl / m3)
    if ( n.Status(props.key_NaCl) != DIRICH )
      n.Store( props.key_NaCl, makeScalar(n.Status(props.key_NaCl),salinity(H2Ocomposition(XNACl_aq),rho_brine_CO2)) );


    // 4.2 Updating saturations using the new mass fractions (vf=volume fraction, mf=mass fraction)
    // --------------------------------------------------------------------------------------------
    const double64 mf_Y_phase = (old_rhoCO2 * old_sCO2 - kg_CO2_to_dissolve + kg_H2O_to_transfer) / total_mass;
    //const double64 mf_X_phase = (old_rhoH2O * old_sH2O + kg_CO2_to_dissolve - kg_H2O_to_transfer) / total_mass;
    //  s1 = vf1=1/(1+(rho1/rho2)*(1/mf1-1)) = -(rho2 * mf1) / (rho1 * mf1 - rho2 * mf1 - rho1)
    double64 sCO2new = -(rho_brine_CO2 * mf_Y_phase) / (rhoCO2 * mf_Y_phase - rho_brine_CO2 * mf_Y_phase - rho_brine_CO2);
    assert( sCO2new >= 0. );
    assert( sCO2new < 1. + numeric_limits<double64>::epsilon() );
    sCO2new = std::max( sCO2new, 0. );
    sCO2new = std::min( sCO2new, 1. );
    if ( n.Status(props.key_sCO2) != DIRICH )
      n.Store( props.key_sCO2, makeScalar(n.Status(props.key_sCO2),sCO2new) );
    if ( n.Status(props.key_sH2O) != DIRICH )
      n.Store( props.key_sH2O, makeScalar(n.Status(props.key_sH2O),1. - sCO2new) );
    // RESULT - assuming that saturations add up to 1 and the brine phase is less compressible than the CO2 phase
    double64 sH2Onew = 1. - sCO2new;
    if ( n.Status(props.key_sH2O) != DIRICH )
      n.Store( props.key_sH2O, makeScalar(n.Status(props.key_sH2O),sH2Onew) );
    if ( n.Status(props.key_sCO2) != DIRICH )
      n.Store( props.key_sCO2, makeScalar(n.Status(props.key_sCO2),sCO2new) );
   

    // 4.3 Computing and storing new compositions
    // ------------------------------------------------------------------------------------
    const double64 new_mass_aqueous_phase(rho_brine_CO2 * sH2Onew);
    // aqueous phase
    H2Ocomposition(XCO2)    += kg_CO2_to_dissolve / new_mass_aqueous_phase;
    H2Ocomposition(XH2O)    -= kg_H2O_to_transfer / new_mass_aqueous_phase;
    H2Ocomposition(XNACl_aq) = mass_NaClaq / new_mass_aqueous_phase;
    // renormalisation
    H2Ocomposition /= (H2Ocomposition(XH2O) + H2Ocomposition(XCO2) + H2Ocomposition(XNACl_aq));
    assert( fabs(H2Ocomposition(XH2O) + H2Ocomposition(XCO2) + H2Ocomposition(XNACl_aq)) < 1. + numeric_limits<double64>::epsilon() );
    // carbonic phase: OK
    const double64 new_mass_carbonic_phase(new_mass_aqueous_phase);
    CO2composition(YH2O) += kg_H2O_to_transfer / new_mass_carbonic_phase;
    CO2composition(YH2O)  = std::min( CO2composition(YH2O), maxYH2O_carbonic );
    CO2composition /= (CO2composition(YH2O) + CO2composition(YCO2));
    // storing compositions
    if ( n.Status( props.key_CO2_comp ) != DIRICH ) n.Store( props.key_CO2_comp, CO2composition );
    if ( n.Status( props.key_H2O_comp ) != DIRICH ) n.Store( props.key_H2O_comp, H2Ocomposition );

   
    // 5. computing diffusivity of CO2 in water at the given P, T, salinity
    // ---------------------------------------------------------------------
    const double64 D_Co2 = eos_csmp.D_Co2( Pf, TC, molalityNaCl );
    assert( D_Co2 > 2.0e-9 ); // m2/s, Cadogan et al., 2014
    assert( D_Co2 < 12.5e-9 );
    // store value
    // TODO: add aqueous diffusivity of CO2 to variable list and consider it in transport calculation
    // if ( n.Status(props.key_DCO2aq) != DIRICH )
    //  n.Store( props.key_DCO2aq, makeScalar(n.Status(props.key_DCO2aq),D_Co2) );   // nodal CO2 diffusivity in brine
   

    // 6. Viscosities
    // ------------------------------------------------------------------------------------
    const double64 muH2O = eos_csmp.mu_AqueousPhase( Pf, TC, molalityNaCl );
    const double64 muCO2 = eos_csmp.mu_CarbonicPhase( Pf, TC );
    n.Store( props.key_muH2O, makeScalar(n.Status(props.key_muH2O),muH2O) );
    n.Store( props.key_muCO2, makeScalar(n.Status(props.key_muCO2),muCO2) );



    // 7. Phase compressibilities
    // ------------------------------------------------------------------------------------
    const double64 betaH2O = eos_csmp.C_AqueousPhase( Pf, TC, molalityNaCl );
    const double64 betaCO2 = eos_csmp.C_CarbonicPhase( Pf, TC );
    n.Store( props.key_cH2O, makeScalar(n.Status(props.key_cH2O),betaH2O) );
    n.Store( props.key_cCO2, makeScalar(n.Status(props.key_cCO2),betaCO2) );



    // 8. Volume changes
    // -----------------
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
    const double64 dV = (old_fluid_volume - new_fluid_volume) / old_fluid_volume;
assert( fabs(dV) <= numeric_limits<double64>::epsilon() );
    // RESULT - computing 'nodal fluid volume source' from apparent change in fluid masss
    const double64  nodal_mass_source = new_total_mass - total_mass;
    const double64  expansion_factor  = 1. - (nodal_mass_source / new_bulk_density) / new_fluid_volume;
    // TODO - use this (nodal) mass source
//    n.Store( props.key_msrc, makeScalar(n.Status(props.key_msrc),nodal_mass_source) );
   
    // reporting mass fractions
    if ( verbose ) {
      cerr <<"\nequilibration calculation: phase properties AFTER calculation: ";
      cerr <<"\n\tCO2 saturation:         "<< sCO2new;
      cerr <<"\n\tXCO2:                   "<< H2Ocomposition(XCO2);
      cerr <<"\n\tXH2O:                   "<< H2Ocomposition(XH2O);
      cerr <<"\n\tXNaCl:                  "<< H2Ocomposition(XNACl_aq);
      cerr <<"\n\tYCO2:                   "<< CO2composition(YCO2);
      cerr <<"\n\tYH2O:                   "<< CO2composition(YH2O);
      cerr <<"\n\tnew CO2 density:        "<< rhoCO2 <<" kg/m3";
      cerr <<"\n\tnew brine density:      "<< rho_brine_CO2 <<" kg/m3";
      cerr <<"\n\tCO2 diffusivity in H2O: "<< D_Co2 * 1.0e9 <<" 10^-9 m2/s";
      cerr <<"\n\tCO2 dynamic viscosity:  "<< muCO2 <<" Pa.s";
      cerr <<"\n\tH2O dynamic viscosity:  "<< muH2O <<" Pa.s";
      cerr <<"\n\tfluid expansion caused: "<< expansion_factor <<" (fraction of original volume.)";
      cerr <<"\n";
    }


} // end equilibrateH2O_CO2_NaCl


template void equilibrateH2O_CO2_NaCl<1U>( const variables::VariableSet_CO2GeoSequestration&, Node<1U>&);
template void equilibrateH2O_CO2_NaCl<2U>( const variables::VariableSet_CO2GeoSequestration&, Node<2U>&);
template void equilibrateH2O_CO2_NaCl<3U>( const variables::VariableSet_CO2GeoSequestration&, Node<3U>&);



/**
    Applies  equilibrateH2O_CO2_NaCl() method to nodes of region
*/
template<size_t dim>
void equilibrateFluid( const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain )
 {
    for ( auto nit=model_subdomain.NodesBegin(); nit!=model_subdomain.NodesEnd(); ++nit )
      equilibrateH2O_CO2_NaCl( props, *(*nit) );
 }

template void equilibrateFluid<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>& );
template void equilibrateFluid<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>& );
template void equilibrateFluid<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>& );

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
void updatePorosity( const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain )
 {
    for ( auto eit=model_subdomain.ElementsBegin(); eit!=model_subdomain.ElementsEnd(); ++eit )
      porosityWithSalt( props, *(*eit) );
 }
 
template void updatePorosity<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>& ); 
template void updatePorosity<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>& );
template void updatePorosity<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>& );
 
 
 
} // end csmp


