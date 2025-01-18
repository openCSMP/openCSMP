//
//  PVTX_Calculator_H2O_CO2_NaCl.cpp
//  CSMP_ACGSS_Simulator
//
//  Created by Stephan Matthai on 3/1/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#include "PVTX_Calculator_H2O_CO2_NaCl.h"
#include "Node.h"
#include "Element.h"
#include "ErrorHandler.h"

using namespace std;

/* DEVELOPMENT LOG - reverse chronological order

4/1/2018 SKM
============






4/1/2018 Thomas Driesner
========================

Die erste volle beta-Version ist jetzt im bitbucket. Sazu einige Anmerkungen:
 
- Ich habe darauf verzichtet, die eos public member zu machen.
Da die immer wieder spezielle Werte von mnacl etc. braucht halte ich es fuer sinnvoler, dass in meiner Klasse zu kapseln.
Statt dessen kannst Du die Variablen jetzt direkt ueber  das interface bekommen, siehe header file.
ACTION: Das heisst, dass Du in Deinem PVTX_Calculator einfach phasestate_.X_h2o() aufrufen kannst, z.B.
 
- Das Ding ist nur krude getestet, fuer mehr fehlt mir jetzt die Zeit.
ATTENTION: Fuer den salty aq_carb Fall gibt es Unstimmigkeiten auf der vierten Stelle der mole fractions,
die ich mir nicht erklaeren kann, aber mehr kann ich momentan nicht bieten.

- Ich nehme mal an, dass Dein transport schemes selber weiss, wann ein- oder zwei-phasiger Zustand vorliegt?
ACTION: discern single- from 2-phase case

- In Deiner Funktion EstablishMassBalanceEtc: dort gibt es die total_mass. Diese braucht meine PhaseStateFinder als constructor Argument.
  Es muss zu jedem Zeitpunkt sichergestellt sein, dass darin auch die flash_.massHalite() inbegriffen ist!
ACTION: add solid salt as well

- Weitere Konstruktorargumente sind die bulk mass fractions der drei Komponenten - die muessten also ausserhalb ausgerechnet werden.
  Alternativ koennten wir das anpassen
 
ATTENTION special case: what if all porosity gets clogged by halite
 
Ich hoffe, dass es ansonsten halbwegs selbsterklaerend und bug-arm ist, aber das wissen wir wohl erst morgen. Bin ab vielleicht 8 oder 9 im home office.

*/


namespace csmp {

template<uint32_t dim>
PVTX_Calculator_H2O_CO2_NaCl<dim>::PVTX_Calculator_H2O_CO2_NaCl( const variables::VariableSet_CO2GeoSequestration& props )
  : props_(props),
    Xbulk_(4,numeric_limits<double>::quiet_NaN(),ANY), aq_ph_composition_(3), carb_ph_composition_(2),
    Pf_(numeric_limits<double>::quiet_NaN()), ToC_(numeric_limits<double>::quiet_NaN()),
    flash_( ToC_, Pf_, Xbulk_(TMASS), Xbulk_(XAQ), Xbulk_(XCARB), Xbulk_(XSALT) )
  {
     cout <<"\nPVTX_Calculator_H2O_CO2_NaCl(custom constructor): Make sure to build only once because this is costly!\n";
  } // end constructor
  
  
 
 
/**
    Calculates the volume fraction of water from the masses of the aqueous, carbonic and salt phase.

    @note Vi = mass_i / density_i, the volume fraction is the phase volume over the total volume.

    @attention since precipitatated salt also occupies a fraction of the pore space, the water and carbondioxide saturations
    will be reduced correspondingly.
*/
template<uint32_t dim>
double PVTX_Calculator_H2O_CO2_NaCl<dim>::WaterSaturation() const
{
  return (flash_.massAqueousPhase() / flash_.rho_aq()) /
         (flash_.massAqueousPhase() / flash_.rho_aq() + flash_.massCarbonicPhase() / flash_.rho_carb() + flash_.massHalite() / rhoNaCl_);

} // end WaterSaturation


/**
    Calculate the halite saturation.
*/
template<uint32_t dim>
double PVTX_Calculator_H2O_CO2_NaCl<dim>::HaliteVolumeFraction() const
 {
  return (flash_.massHalite() / rhoNaCl_) /
         (flash_.massAqueousPhase() / flash_.rho_aq() + flash_.massCarbonicPhase() / flash_.rho_carb() + flash_.massHalite() / rhoNaCl_);

 } // end HaliteVolumeFraction

  
/**
     MASTER FUNCTION to be executed after each transport step to establish saturations, densities etc.

    Input: (read from Model via the ArrayVariable key_Xbulk)
      - pf, T
      - composition: total_mass, Xaq, Xcarb, Xsalt
 
    Output: (everything gets stored directly on the model
      - mass_aq_phase
      - mass_carb_phase
      (-water saturation)
      - aq_ph_composition_, CO2composition - composition of the phases in terms of new mass fractions
      - the densities, viscosities and compressibilities of the phases
 */
template<uint32_t dim>
SYSTEM_STATE PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate( Node<dim>* n )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
    // 0. establish fluid pressure and temperature
    Pf_ = n->Read( props_.key_pf );
    ToC_ = n->Read( props_.key_T ); // temperature in centrigrade
    if ( Pf_ <= 100325. ) {
         cerr <<"\n\t"<< Pf_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'fluid pressure' too low for EOS-scheme (min pf=1bar); errors may occur");
      }
    else if ( Pf_ >= 6e7 ) {
         cerr <<"\n\t"<< Pf_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'fluid pressure' too high for EOS-scheme (max pf=600bar); errors may occur");
      }
    if ( ToC_ <= 12. ) {
         cerr <<"\n\t"<< ToC_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'temperature' too low for EOS-scheme (min T(oc)=12); errors may occur");
      }
    else if ( ToC_ >= 100. ) {
         cerr <<"\n\t"<< ToC_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'temperature' too high for EOS-scheme (max T(oc)=100); errors may occur");
      }

    // 1. reading current mass balance from last transport step; mass in unit volume of fluid
    n->Read( props_.key_Xbulk, Xbulk_ );
 
    // 2. flash calculation
    SYSTEM_STATE phase_state = flash_.Equilibrate();
 
    // 3. updating YCO2, YH2O, XCO2, XH2O, XNaCl mass fractions of the phases are stored in aq_ph_composition_, carb_ph_composition_
    switch(phase_state)
    {
    case SYSTEM_STATE::aq :
        aq_ph_composition_(XH2O)     = flash_.X_h2o();
        aq_ph_composition_(XCO2)     = flash_.X_co2();
        aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
        assert( !isnan(flash_.beta_aq()) );
        // saturations
        if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O), 1.) );
        if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2), 0.) );
        if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl), 0.) );
        // density and viscosity
        if ( n->Status( props_.key_H2O_comp ) != DIRICH ) n->Store( props_.key_H2O_comp, aq_ph_composition_ );
        if ( n->Status(props_.key_rhoH2O) != DIRICH ) n->Store( props_.key_rhoH2O, makeScalar(n->Status(props_.key_rhoH2O),flash_.rho_aq()) );
        if ( n->Status(props_.key_muH2O) != DIRICH ) n->Store( props_.key_muH2O, makeScalar(n->Status(props_.key_muH2O),flash_.mu_aq()) );
        if ( n->Status(props_.key_cH2O) != DIRICH ) n->Store( props_.key_cH2O, makeScalar(n->Status(props_.key_cH2O),flash_.beta_aq()) );
        // dissolved CO2 (kg/m3)
        if ( n->Status(props_.key_CO2aq) != DIRICH ) {
          //const double dissolved_CO2 = aq_ph_composition_(XCO2) * flash_.rho_aq();
          const double dissolved_CO2 = aq_ph_composition_(XCO2);
          n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),dissolved_CO2) );
        }
        // salinity (kg/m3)
        if ( n->Status(props_.key_NaClaq) != DIRICH ) {
          //double salinity = aq_ph_composition_(XNACl_aq) * flash_.rho_aq();
          double salinity = aq_ph_composition_(XNACl_aq);
          n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),salinity) );
        }
        // evaporated water (kg/m3)
        if ( n->Status(props_.key_H2Og) != DIRICH ) {
          n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),0.) );
        }                
      break;
     
    case SYSTEM_STATE::carb :
        carb_ph_composition_(YCO2) = flash_.Y_co2();
        carb_ph_composition_(YH2O) = flash_.Y_h2o();
        // saturations
        if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O), 0.) );
        if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2), 1.) );
        if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl), 0.) );
       // density and viscosity
        if ( n->Status( props_.key_CO2_comp ) != DIRICH ) n->Store( props_.key_CO2_comp, carb_ph_composition_ );
        if ( n->Status(props_.key_rhoCO2) != DIRICH ) n->Store( props_.key_rhoCO2, makeScalar(n->Status(props_.key_rhoCO2),flash_.rho_carb()) );
        if ( n->Status(props_.key_muCO2) != DIRICH ) n->Store( props_.key_muCO2, makeScalar(n->Status(props_.key_muCO2),flash_.mu_carb()) );
        if ( n->Status(props_.key_cCO2) != DIRICH ) n->Store( props_.key_cCO2, makeScalar(n->Status(props_.key_cCO2),flash_.beta_carb()) );
        // dissolved CO2 (kg/m3)
        if ( n->Status(props_.key_CO2aq) != DIRICH ) {
          n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),0.) );
        }
        // salinity (kg/m3)
        if ( n->Status(props_.key_NaClaq) != DIRICH ) {
          n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),0.) );
        }
        // evaporated water (kg/m3)
        if ( n->Status(props_.key_H2Og) != DIRICH ) {
          //double evaporated_water = carb_ph_composition_(YH2O) * flash_.rho_carb();
          double evaporated_water = carb_ph_composition_(YH2O);
          n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),evaporated_water) );
        }              
      break;
      
    case SYSTEM_STATE::salt :
        // saturations
        if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O), 0.) );
        if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2), 0.) );
        if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl), 1.) );
        // dissolved CO2 (kg/m3)
        if ( n->Status(props_.key_CO2aq) != DIRICH ) {
          n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),0.) );
        }
        // salinity (kg/m3)
        if ( n->Status(props_.key_NaClaq) != DIRICH ) {
          n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),0.) );
        }
        // evaporated water (kg/m3)
        if ( n->Status(props_.key_H2Og) != DIRICH ) {
          n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),0.) );
        }               
      break;
      
    case SYSTEM_STATE::aq_salt : {
        aq_ph_composition_(XH2O)     = flash_.X_h2o();
        aq_ph_composition_(XCO2)     = flash_.X_co2();
        aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
        assert( !isnan(flash_.beta_aq()) );
        // saturations
        const double sw = (flash_.massAqueousPhase() / flash_.rho_aq()) /
                            (flash_.massAqueousPhase() / flash_.rho_aq() + flash_.massHalite() / rhoNaCl_);
      
        if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O), sw) );
        if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2), 0.) );
        if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl), 1. - sw) );
        // density and viscosity
        if ( n->Status( props_.key_H2O_comp ) != DIRICH ) n->Store( props_.key_H2O_comp, aq_ph_composition_ );
        if ( n->Status(props_.key_rhoH2O) != DIRICH ) n->Store( props_.key_rhoH2O, makeScalar(n->Status(props_.key_rhoH2O),flash_.rho_aq()) );
        if ( n->Status(props_.key_muH2O) != DIRICH ) n->Store( props_.key_muH2O, makeScalar(n->Status(props_.key_muH2O),flash_.mu_aq()) );
        if ( n->Status(props_.key_cH2O) != DIRICH ) n->Store( props_.key_cH2O, makeScalar(n->Status(props_.key_cH2O),flash_.beta_aq()) );
        // dissolved CO2 (kg/m3)
        if ( n->Status(props_.key_CO2aq) != DIRICH ) {
          //const double dissolved_CO2 = aq_ph_composition_(XCO2) * flash_.rho_aq();
          const double dissolved_CO2 = aq_ph_composition_(XCO2);
          n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),dissolved_CO2) );
        }
        // salinity (kg/m3)
        if ( n->Status(props_.key_NaClaq) != DIRICH ) {
          //double salinity = aq_ph_composition_(XNACl_aq) * flash_.rho_aq();
          double salinity = aq_ph_composition_(XNACl_aq);
          n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),salinity) );
        }
        // evaporated water (kg/m3)
        if ( n->Status(props_.key_H2Og) != DIRICH ) {
          n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),0.) );
        }
        }
      break;
      
    case SYSTEM_STATE::aq_carb : {
          aq_ph_composition_(XH2O) = flash_.X_h2o();
          aq_ph_composition_(XCO2) = flash_.X_co2();
          aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
          assert( !isnan(flash_.beta_aq()) );
          carb_ph_composition_(YCO2) = flash_.Y_co2();
          carb_ph_composition_(YH2O) = flash_.Y_h2o();
          // saturations
          const double sw = (flash_.massAqueousPhase() / flash_.rho_aq()) /
                              (flash_.massAqueousPhase() / flash_.rho_aq() + flash_.massCarbonicPhase() / flash_.rho_carb());

          if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O),sw) );
          if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2),1.-sw) ); // since there is no salt
          if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl),0.) );
          // density and viscosity
          if ( n->Status( props_.key_H2O_comp ) != DIRICH ) n->Store( props_.key_H2O_comp, aq_ph_composition_ );
          if ( n->Status(props_.key_rhoH2O) != DIRICH ) n->Store( props_.key_rhoH2O, makeScalar(n->Status(props_.key_rhoH2O),flash_.rho_aq()) );
          if ( n->Status(props_.key_muH2O) != DIRICH ) n->Store( props_.key_muH2O, makeScalar(n->Status(props_.key_muH2O),flash_.mu_aq()) );
          if ( n->Status(props_.key_cH2O) != DIRICH ) n->Store( props_.key_cH2O, makeScalar(n->Status(props_.key_cH2O),flash_.beta_aq()) );
          if ( n->Status( props_.key_CO2_comp ) != DIRICH ) n->Store( props_.key_CO2_comp, carb_ph_composition_ );
          if ( n->Status(props_.key_rhoCO2) != DIRICH ) n->Store( props_.key_rhoCO2, makeScalar(n->Status(props_.key_rhoCO2),flash_.rho_carb()) );
          if ( n->Status(props_.key_muCO2) != DIRICH ) n->Store( props_.key_muCO2, makeScalar(n->Status(props_.key_muCO2),flash_.mu_carb()) );
          if ( n->Status(props_.key_cCO2) != DIRICH ) n->Store( props_.key_cCO2, makeScalar(n->Status(props_.key_cCO2),flash_.beta_carb()) );
          // dissolved CO2 (kg/m3)
          if ( n->Status(props_.key_CO2aq) != DIRICH ) {
            //const double dissolved_CO2 = aq_ph_composition_(XCO2) * flash_.rho_aq();
            const double dissolved_CO2 = aq_ph_composition_(XCO2);
            n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),dissolved_CO2) );
          }
          // salinity (kg/m3)
          if ( n->Status(props_.key_NaClaq) != DIRICH ) {
            //double salinity = aq_ph_composition_(XNACl_aq) * flash_.rho_aq();
            double salinity = aq_ph_composition_(XNACl_aq);
            n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),salinity) );
          }
          // evaporated water (kg/m3)
          if ( n->Status(props_.key_H2Og) != DIRICH ) {
            //double evaporated_water = carb_ph_composition_(YH2O) * flash_.rho_carb();
            double evaporated_water = carb_ph_composition_(YH2O);
            n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),evaporated_water) );
          }
        }
      break;
      
    case SYSTEM_STATE::carb_salt : {
        carb_ph_composition_(YCO2) = flash_.Y_co2();
        carb_ph_composition_(YH2O) = flash_.Y_h2o();
        // saturations
        const double halite_saturation = (flash_.massHalite() / rhoNaCl_) /
                                           (flash_.massCarbonicPhase() / flash_.rho_carb() + flash_.massHalite() / rhoNaCl_);
      
        if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O),0.) );
        if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2),1.-halite_saturation) );
        if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl),halite_saturation) );
        // density and viscosity
        if ( n->Status( props_.key_CO2_comp ) != DIRICH ) n->Store( props_.key_CO2_comp, carb_ph_composition_ );
        if ( n->Status(props_.key_rhoCO2) != DIRICH ) n->Store( props_.key_rhoCO2, makeScalar(n->Status(props_.key_rhoCO2),flash_.rho_carb()) );
        if ( n->Status(props_.key_muCO2) != DIRICH ) n->Store( props_.key_muCO2, makeScalar(n->Status(props_.key_muCO2),flash_.mu_carb()) );
        if ( n->Status(props_.key_cCO2) != DIRICH ) n->Store( props_.key_cCO2, makeScalar(n->Status(props_.key_cCO2),flash_.beta_carb()) );
        // dissolved CO2 (kg/m3)
        if ( n->Status(props_.key_CO2aq) != DIRICH ) {
          n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),0.) );
        }
        // salinity (kg/m3)
        if ( n->Status(props_.key_NaClaq) != DIRICH ) {
          n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),0.) );
        }
        // evaporated water (kg/m3)
        if ( n->Status(props_.key_H2Og) != DIRICH ) {
          //double evaporated_water = carb_ph_composition_(YH2O) * flash_.rho_carb();
          double evaporated_water = carb_ph_composition_(YH2O);
          n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),evaporated_water) );
        }
        }
      break;
      
    case SYSTEM_STATE::aq_carb_salt : {
          aq_ph_composition_(XH2O)     = flash_.X_h2o();
          aq_ph_composition_(XCO2)     = flash_.X_co2();
          aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
          assert( !isnan(flash_.beta_aq()) );
          carb_ph_composition_(YCO2)   = flash_.Y_co2();
          carb_ph_composition_(YH2O)   = flash_.Y_h2o();
          // saturations
          const double sw                = WaterSaturation();
          const double halite_saturation = HaliteVolumeFraction();
          if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O),sw) );
          if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2),1.-sw-halite_saturation) );
          if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl),halite_saturation) );
          // density and viscosity
          if ( n->Status( props_.key_H2O_comp ) != DIRICH ) n->Store( props_.key_H2O_comp, aq_ph_composition_ );
          if ( n->Status(props_.key_rhoH2O) != DIRICH ) n->Store( props_.key_rhoH2O, makeScalar(n->Status(props_.key_rhoH2O),flash_.rho_aq()) );
          if ( n->Status(props_.key_muH2O) != DIRICH ) n->Store( props_.key_muH2O, makeScalar(n->Status(props_.key_muH2O),flash_.mu_aq()) );
          if ( n->Status(props_.key_cH2O) != DIRICH ) n->Store( props_.key_cH2O, makeScalar(n->Status(props_.key_cH2O),flash_.beta_aq()) );
          if ( n->Status( props_.key_CO2_comp ) != DIRICH ) n->Store( props_.key_CO2_comp, carb_ph_composition_ );
          if ( n->Status(props_.key_rhoCO2) != DIRICH ) n->Store( props_.key_rhoCO2, makeScalar(n->Status(props_.key_rhoCO2),flash_.rho_carb()) );
          if ( n->Status(props_.key_muCO2) != DIRICH ) n->Store( props_.key_muCO2, makeScalar(n->Status(props_.key_muCO2),flash_.mu_carb()) );
          if ( n->Status(props_.key_cCO2) != DIRICH ) n->Store( props_.key_cCO2, makeScalar(n->Status(props_.key_cCO2),flash_.beta_carb()) );
          // dissolved CO2 (kg/m3)
          if ( n->Status(props_.key_CO2aq) != DIRICH ) {
            //const double dissolved_CO2 = aq_ph_composition_(XCO2) * flash_.rho_aq();
            const double dissolved_CO2 = aq_ph_composition_(XCO2);
            n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),dissolved_CO2) );
          }
          // salinity (kg/m3)
          if ( n->Status(props_.key_NaClaq) != DIRICH ) {
            //double salinity = aq_ph_composition_(XNACl_aq) * flash_.rho_aq();
            double salinity = aq_ph_composition_(XNACl_aq);
            n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),salinity) );
          }
          // evaporated water (kg/m3)
          if ( n->Status(props_.key_H2Og) != DIRICH ) {
            //double evaporated_water = carb_ph_composition_(YH2O) * flash_.rho_carb();
            double evaporated_water = carb_ph_composition_(YH2O);
            n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),evaporated_water) );
          }
        }
      break;
      
    case SYSTEM_STATE::full : {
          aq_ph_composition_(XH2O)     = flash_.X_h2o();
          aq_ph_composition_(XCO2)     = flash_.X_co2();
          aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
          assert( !isnan(flash_.beta_aq()) );
          carb_ph_composition_(YCO2)   = flash_.Y_co2();
          carb_ph_composition_(YH2O)   = flash_.Y_h2o();
          // saturations
          const double sw = WaterSaturation();
          const double halite_saturation = HaliteVolumeFraction();
          if ( n->Status( props_.key_sH2O ) != DIRICH ) n->Store( props_.key_sH2O, makeScalar(n->Status(props_.key_sH2O),sw) );
          if ( n->Status( props_.key_sCO2 ) != DIRICH ) n->Store( props_.key_sCO2, makeScalar(n->Status(props_.key_sCO2),1.-sw-halite_saturation) );
          if ( n->Status( props_.key_NaCl ) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl),halite_saturation) );
          // density and viscosity
          if ( n->Status( props_.key_H2O_comp ) != DIRICH ) n->Store( props_.key_H2O_comp, aq_ph_composition_ );
          if ( n->Status(props_.key_rhoH2O) != DIRICH ) n->Store( props_.key_rhoH2O, makeScalar(n->Status(props_.key_rhoH2O),flash_.rho_aq()) );
          if ( n->Status(props_.key_muH2O) != DIRICH ) n->Store( props_.key_muH2O, makeScalar(n->Status(props_.key_muH2O),flash_.mu_aq()) );
          if ( n->Status(props_.key_cH2O) != DIRICH ) n->Store( props_.key_cH2O, makeScalar(n->Status(props_.key_cH2O),flash_.beta_aq()) );
          if ( n->Status( props_.key_CO2_comp ) != DIRICH ) n->Store( props_.key_CO2_comp, carb_ph_composition_ );
          if ( n->Status(props_.key_rhoCO2) != DIRICH ) n->Store( props_.key_rhoCO2, makeScalar(n->Status(props_.key_rhoCO2),flash_.rho_carb()) );
          if ( n->Status(props_.key_muCO2) != DIRICH ) n->Store( props_.key_muCO2, makeScalar(n->Status(props_.key_muCO2),flash_.mu_carb()) );
          if ( n->Status(props_.key_cCO2) != DIRICH ) n->Store( props_.key_cCO2, makeScalar(n->Status(props_.key_cCO2),flash_.beta_carb()) );
          if ( n->Status(props_.key_NaCl) != DIRICH ) n->Store( props_.key_NaCl, makeScalar(n->Status(props_.key_NaCl),flash_.massHalite()/rhoNaCl_) );
          // dissolved CO2 (kg/m3)
          if ( n->Status(props_.key_CO2aq) != DIRICH ) {
            //const double dissolved_CO2 = aq_ph_composition_(XCO2) * flash_.rho_aq();
            const double dissolved_CO2 = aq_ph_composition_(XCO2);
            n->Store( props_.key_CO2aq, makeScalar(n->Status(props_.key_CO2aq),dissolved_CO2) );
          }
          // salinity (kg/m3)
          if ( n->Status(props_.key_NaClaq) != DIRICH ) {
            //double salinity = aq_ph_composition_(XNACl_aq) * flash_.rho_aq();
            double salinity = aq_ph_composition_(XNACl_aq);
            n->Store( props_.key_NaClaq, makeScalar(n->Status(props_.key_NaClaq),salinity) );
          }
          // evaporated water (kg/m3)
          if ( n->Status(props_.key_H2Og) != DIRICH ) {
            //double evaporated_water = carb_ph_composition_(YH2O) * flash_.rho_carb();
            double evaporated_water = carb_ph_composition_(YH2O);
            n->Store( props_.key_H2Og, makeScalar(n->Status(props_.key_H2Og),evaporated_water) );
          }
        }
      break;
      
    default:   
        cerr << "ERROR: PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate( Node<dim>* n, double delta_t ):\n";  
        cerr << "wrong state, state provided was " << parseState(phase_state) << endl;
    }  
 
    // 5. diffusivity of CO2 in the aqueous phase
    // ------------------------------------------
    static bool first_warning(true);
    if ( flash_.massAqueousPhase() > 0. && Xbulk_(XCARB) > 0. ) {
        const double diff_co2_aq = flash_.D_Co2();
        if ( first_warning && (diff_co2_aq < 2.0e-10 || diff_co2_aq > 12.5e-9) ) {
             cerr <<"\n\t"<< diff_co2_aq <<" m2/s.";
             csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate:",
                               "CO2 diffusivity in brine is out of bounds as defined by Cardogan et al. 2014.");
             first_warning = false;
          }
      }

    n->Store( props_.key_nPHS, makeScalar(n->Status(props_.key_nPHS),static_cast<int>(phase_state)) );
 
    return phase_state;
 
} // end Equilibrate ( Node version)








/**
     @attention does not compute source term.
*/
template<uint32_t dim>
SYSTEM_STATE PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate( Element<dim>* e )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
    // 1. computes current mass balance; apply after last transport step; returns the total mass in finite volume cell connected to Node
    InterpolateInputVariablesToBaryCenter( e );
    if ( Pf_ <= 100325. ) {
         cerr <<"\n\t"<< Pf_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'fluid pressure' too low for EOS-scheme (min pf=1bar); errors may occur");
      }
    else if ( Pf_ >= 6e7 ) {
         cerr <<"\n\t"<< Pf_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'fluid pressure' too high for EOS-scheme (max pf=600bar); errors may occur");
      }
    if ( ToC_ <= 12. ) {
         cerr <<"\n\t"<< ToC_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'temperature' too low for EOS-scheme (min T(oc)=12); errors may occur");
      }
    else if ( ToC_ >= 100. ) {
         cerr <<"\n\t"<< ToC_;
         csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate(Node):",
                           "input 'temperature' too high for EOS-scheme (max T(oc)=100); errors may occur");
      }

    // 2. flash calculation
    SYSTEM_STATE phase_state = flash_.Equilibrate();

    switch(phase_state)
      {
          case SYSTEM_STATE::aq :
              aq_ph_composition_(XH2O)     = flash_.X_h2o();
              aq_ph_composition_(XCO2)     = flash_.X_co2();
              aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
            break;
          case SYSTEM_STATE::carb :
              carb_ph_composition_(YCO2)   = flash_.Y_co2();
              carb_ph_composition_(YH2O)   = flash_.Y_h2o();
            break;
          case SYSTEM_STATE::salt :
              throw csmp::Exception( ERROR, "PVTX_Calculator::Equilibrate(Element):", "salt only case not handled yet.");
            break;
          case SYSTEM_STATE::aq_salt :
              aq_ph_composition_(XH2O)     = flash_.X_h2o();
              aq_ph_composition_(XCO2)     = flash_.X_co2();
              aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
            break;
          case SYSTEM_STATE::aq_carb :
              aq_ph_composition_(XH2O)     = flash_.X_h2o();
              aq_ph_composition_(XCO2)     = flash_.X_co2();
              aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
              carb_ph_composition_(YCO2)   = flash_.Y_co2();
              carb_ph_composition_(YH2O)   = flash_.Y_h2o();
            break;
          case SYSTEM_STATE::carb_salt :
              carb_ph_composition_(YCO2)   = flash_.Y_co2();
              carb_ph_composition_(YH2O)   = flash_.Y_h2o();
            break;
          case SYSTEM_STATE::aq_carb_salt :
             aq_ph_composition_(XH2O)     = flash_.X_h2o();
              aq_ph_composition_(XCO2)     = flash_.X_co2();
              aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
              carb_ph_composition_(YCO2)   = flash_.Y_co2();
              carb_ph_composition_(YH2O)   = flash_.Y_h2o();
            break;
          case SYSTEM_STATE::full :
              aq_ph_composition_(XH2O)     = flash_.X_h2o();
              aq_ph_composition_(XCO2)     = flash_.X_co2();
              aq_ph_composition_(XNACl_aq) = flash_.X_nacl();
              carb_ph_composition_(YCO2)   = flash_.Y_co2();
              carb_ph_composition_(YH2O)   = flash_.Y_h2o();
            break;
          default:
              cerr << "ERROR: PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate( Node<dim>* n, double delta_t ):\n";
              cerr << "wrong state, state provided was " << parseState(phase_state) << endl;
    }

    // 5. diffusivity of CO2 in the aqueous phase, if there is CO2 around
    // ------------------------------------------------------------------
    const double diff_CO2_min(2.0e-10), diff_CO2_max(12.5e-9); // bounds by Cardogan et al. 2014
    double  diff_co2_aq(0.);
    if ( flash_.massAqueousPhase() > 0. && Xbulk_(XCARB) > 0. ) {
        diff_co2_aq = flash_.D_Co2();
        if ( diff_co2_aq < diff_CO2_min || diff_co2_aq > diff_CO2_max ) {
             cerr <<"\n\t"<< diff_co2_aq <<" m2/s.";
             csmp_error.Note( WARNING, "PVTX_Calculator_H2O_CO2_NaCl<dim>::Equilibrate:",
                                         "CO2 diffusivity in brine is out of bounds as defined by Cardogan et al. 2014; setting it to minimum value.");
             diff_co2_aq = diff_CO2_min;
          }
      }
    // storing the aqueous diffusivity of CO2
    e->Store( props_.key_diff, makeScalar(e->Status(props_.key_diff),0.) );

    e->Store( props_.key_ePHS, makeScalar(e->Status(props_.key_ePHS),static_cast<int>(phase_state)) );
    return phase_state;
 
} // end Equilibrate (Element verson)







/**
    TODO: deal with the salt-only case, in terms of saturations and mass balances, or use trick that there is always some water left
*/
template<uint32_t dim>
void PVTX_Calculator_H2O_CO2_NaCl<dim>::InitialisePVTX_FromFieldData( const Node<dim>* n )
 {
    // 0. establish locai conditions
    // -----------------------------
    Pf_ = n->Read( props_.key_pf );
    ToC_ = n->Read( props_.key_T ); // temperature in centrigrade
    assert( Pf_ >= 100325. );
    assert( ToC_ <= 100. );
    const double entry_pressure( MaxEntryPressureOfParentElements( n ) );
    Pf_ += entry_pressure;

    // 1. Reading fluid composition after transport and computing mass balances
    // ------------------------------------------------------------------------
    // saturation
    const double min_sat(numeric_limits<double>::epsilon() * 2.);
    const double sw     = n->Read( props_.key_sH2O );
    const double sCO2   = n->Read( props_.key_sCO2 );
    const double halite = n->Read( props_.key_NaCl ); // volume fraction of salt in pore space
    assert( !isnan(halite) );
    assert( halite >= 0. and halite <= 1. );
    assert( sw >= 0. and sw <= 1. );
    assert( sCO2 >= 0. and sCO2 <= 1. );
    assert( fabs(sw+sCO2+halite) <= 1. + min_sat );
    const double rhoH2O = (sw > 0.) ? n->Read( props_.key_rhoH2O ) : 0.;
    const double rhoCO2 = (sCO2 > 0.) ? n->Read( props_.key_rhoCO2 ) : 0.;

    // mass totals used below for the computation of the new phase saturations and compositions - OK rsvp with Thomas
    double mass_H2O, mass_CO2, mass_NaCl_aq;
    // H2O
    if ( sw >= min_sat ) {
         n->Read( props_.key_H2O_comp, aq_ph_composition_ );
         if ( fabs(aq_ph_composition_(XH2O)+aq_ph_composition_(XCO2)+aq_ph_composition_(XNACl_aq)-1.) > numeric_limits<double>::epsilon() )
           throw csmp::Exception( ERROR, "PVTX_Calculator_H2O_CO2_NaCl<dim>::EstablishMassBalanceAnd_PT_Conditions:",
                                         "mass fractions of aqueous phase do not add up to 1.");
         mass_H2O     = aq_ph_composition_(XH2O) * rhoH2O * sw;
         mass_NaCl_aq = aq_ph_composition_(XNACl_aq) * rhoH2O * sw;
         if ( sCO2 >= min_sat )
           mass_H2O += carb_ph_composition_(YH2O) * rhoCO2 * sCO2;
      }
    else {
         mass_H2O     = 0.;
         mass_NaCl_aq = 0.;
      }
    // CO2
    if ( sCO2 >= min_sat ) {
         n->Read( props_.key_CO2_comp, carb_ph_composition_ );
         if ( fabs(carb_ph_composition_(YH2O)+carb_ph_composition_(YCO2)-1.) > numeric_limits<double>::epsilon() )
           throw csmp::Exception( ERROR, "PVTX_Calculator_H2O_CO2_NaCl<dim>::EstablishMassBalanceAnd_PT_Conditions:",
                                         "mass fractions of carbonic phase do not add up to 1.");
         mass_CO2 = carb_ph_composition_(YCO2) * rhoCO2 * sCO2;
         if ( sw >= min_sat )
           mass_CO2 += aq_ph_composition_(XCO2) * rhoH2O * sw;
      }
    else mass_CO2 = 0.;

    // NaCl - solid salt volume converted into mass, assuming a unit pore volume
    const double mass_NaCl_cr = halite * rhoNaCl_;

    // derived quantities
    Xbulk_(TMASS) = mass_H2O + mass_CO2 + mass_NaCl_aq + mass_NaCl_cr;

    // global mass fractions Z in system
    Xbulk_(XAQ)   = mass_H2O / Xbulk_(TMASS);
    Xbulk_(XCARB) = mass_CO2 / Xbulk_(TMASS);
    Xbulk_(XSALT) = (mass_NaCl_aq + mass_NaCl_cr) / Xbulk_(TMASS);
   
//cerr <<"\nNode: "<< n->Idx();
//Xbulk_.Out();
 
 } // end InitialisePVTX_FromFieldData (Node)










/**
    Establishes the input values: P,T conditions and mass fractions of the system for Equilibrate()
    interpolating them to the element barycentre.
 
    @test SKM 30/1/19
*/
template<uint32_t dim>
void PVTX_Calculator_H2O_CO2_NaCl<dim>::InterpolateInputVariablesToBaryCenter( const Element<dim>* e )
 {
    // 0. interpolating all relevant properties to element barycentre
    const auto nodes(e->Nodes());
    IPOL_.resize(nodes);
    e->N_AtBaryCenter( IPOL_ );
    const double PV = e->Volume() * e->Read(props_.key_thi) * e->Read(props_.key_phi);
    ArrayVariable  temp(4,0.,ANY);
    double       ePV_from_sectors(0.);

    // properties
    Xbulk_ = Pf_ = ToC_ = 0.;
    double NaCl(0.);
    // interpolation
    for ( auto i{0U}; i<nodes; ++i ) {
         // pressure and temperature
         Pf_  += IPOL_[i] * e->N(i)->Read( props_.key_pf );
         ToC_ += IPOL_[i] * e->N(i)->Read( props_.key_T );
         // salt as a mineral
         NaCl += IPOL_[i] * e->N(i)->Read( props_.key_NaCl );
         // mass fractions and total mass using the sector pore volumes as weights
         //const double sector_pv = e->Read( i, sector_ip, props_.key_sPV );
         const double sector_pv = e->SectorVolume(i) * e->Read(props_.key_thi) * e->Read(props_.key_phi);
         assert( !isnan(sector_pv) );
         ePV_from_sectors += sector_pv;
         e->N(i)->Read( props_.key_Xbulk, temp );
         // TODO: verify this weighting, does it need IPOL[i] ? - is the thickness taken into account in the sector PV
         Xbulk_(TMASS) += temp(TMASS) * (sector_pv / PV);
         Xbulk_(XAQ)   += IPOL_[i] * temp(XAQ);
         Xbulk_(XCARB) += IPOL_[i] * temp(XCARB);
         Xbulk_(XSALT) += IPOL_[i] * temp(XSALT);
      }
    assert( fabs(PV - ePV_from_sectors)  <= numeric_limits<double>::epsilon()*PV );
   
    // checks and taking into account capillary pressure
    const double entry_pressure = e->Read( props_.key_pd );
    assert( !isnan(entry_pressure) );
    Pf_ += entry_pressure;
    assert( Pf_ >= 100325. );
    assert( ToC_ <= 100. );

    // 1. Checking the fluid composition after the interpolation
    // ------------------------------------------------------------------------
    if ( fabs(Xbulk_(XAQ)+Xbulk_(XCARB)+Xbulk_(XSALT)-1.) > 1.0e-8 ) //numeric_limits<double>::epsilon() is too strict
      throw csmp::Exception( ERROR, "PVTX_Calculator_H2O_CO2_NaCl<dim>::EstablishMassBalanceAnd_PT_Conditions:",
                                    "mass fractions of system do not add up to 1.");
 
 } // end InterpolateInputVariablesToBaryCenter (Element)
















/**
     Loops over the parent elements of the node determining the maximum capillary entry pressure therein.
 
     @return returns the maximum pd detected.
*/
template<uint32_t dim>
double PVTX_Calculator_H2O_CO2_NaCl<dim>::MaxEntryPressureOfParentElements( const Node<dim>* n ) const
 {
    assert( n != nullptr );
    double pd(0.);
   
    for ( auto i{0U}; i<n->Parents(); ++i ) {
         const Element<dim>* eptr(n->Parent(i));
         assert( eptr != nullptr );
         pd = max( pd, eptr->Read(props_.key_pd) );
      }
   
    return pd;
   
 } // end MaxEntryPressureOfParentElements




// NON-MEMBER FUNCTIONS


/**
    Returns which phases which are present on all nodes of the element using either of the options

    NEITHER, AQUEOUS, CARBONIC, HALITE, AQUEOUS_CARBONIC, AQUEOUS_CARBONIC_HALITE.

    @test OK SKM 20/1/2019
*/
template<uint32_t dim>
PHASES_CONTINUOUS_ACROSS_ELEMENT  continuousPhases( const variables::VariableSet_CO2GeoSequestration& props, const Element<dim>* eptr )
 {
    assert( eptr != nullptr );
    int aqueous(0), carbonic(0), salt(0);

    const auto nodes{eptr->Nodes()};
    for ( auto i{0U}; i<nodes; ++i ) {
         const SYSTEM_STATE state(static_cast<SYSTEM_STATE>(static_cast<int>(eptr->N(i)->Read(props.key_nPHS))));
         // counting the phases
         if      ( state == SYSTEM_STATE::aq   ) aqueous++;
         else if ( state == SYSTEM_STATE::carb ) carbonic++;
         else if ( state == SYSTEM_STATE::aq_carb ) {
              aqueous++;
              carbonic++;
           }
         // cases with salt
         else if ( state == SYSTEM_STATE::salt ) salt++;
         else if ( state == SYSTEM_STATE::aq_salt ) {
              aqueous++;
              salt++;
           }
         // carbonic only
         else if ( state == SYSTEM_STATE::carb_salt ) {
              carbonic++;
              salt++;
           }
         // aqueus + carbonic + salt
         else if ( state == SYSTEM_STATE::aq_carb_salt || state == SYSTEM_STATE::full ) {
              aqueous++;
              carbonic++;
              salt++;
           }              
      }
    
    // diagnosing the state of the element
    if ( aqueous == nodes ) {
          if ( carbonic < nodes ) {
               if ( salt < nodes ) return AQUEOUS;
               else return AQUEOUS_HALITE;
            }
          else { // if carbonic == nodes
               if ( salt < nodes ) return AQUEOUS_CARBONIC;
               else return AQUEOUS_CARBONIC_HALITE;
            }
      }

    if ( carbonic == nodes ) {
         // if  aqueous < nodes
         if ( salt < nodes ) return CARBONIC;
         else return CARBONIC_HALITE;
      }

    // the cases of salt with brine and carbondioxide were already handled
    if ( salt == nodes ) return HALITE;    

    return NEITHER;
   
 } // phasesFlowingAcrossElement
 
template PHASES_CONTINUOUS_ACROSS_ELEMENT  continuousPhases( const variables::VariableSet_CO2GeoSequestration&, const Element<1U>* );
template PHASES_CONTINUOUS_ACROSS_ELEMENT  continuousPhases( const variables::VariableSet_CO2GeoSequestration&, const Element<2U>* );
template PHASES_CONTINUOUS_ACROSS_ELEMENT  continuousPhases( const variables::VariableSet_CO2GeoSequestration&, const Element<3U>* );



/**
    Reports which SYSTEM_STATE an element is in when the given phases are continuous across it.

    PHASES_CONTINUOUS_ACROSS_ELEMENT { AQUEOUS, CARBONIC, HALITE, AQUEOUS_HALITE, AQUEOUS_CARBONIC,
                                        CARBONIC_HALITE, AQUEOUS_CARBONIC_HALITE, NEITHER };
 
    @note needed for unit testing.
*/
SYSTEM_STATE  elementState( PHASES_CONTINUOUS_ACROSS_ELEMENT cphases )
 {
    switch( cphases ) {
         case AQUEOUS:
           return SYSTEM_STATE::aq;
         case CARBONIC:
           return SYSTEM_STATE::aq;
         case HALITE:
           return SYSTEM_STATE::aq;
         case AQUEOUS_HALITE:
           return SYSTEM_STATE::aq;
         case AQUEOUS_CARBONIC:
           return SYSTEM_STATE::aq;
         case CARBONIC_HALITE:
           return SYSTEM_STATE::aq;
         case AQUEOUS_CARBONIC_HALITE:
           return SYSTEM_STATE::aq;
         case NEITHER:
           cerr <<"\nelementState: no suitable conversion found.\n";
           return SYSTEM_STATE::full;
      }
    return SYSTEM_STATE::undefined;
    
 } //end elementState






template class PVTX_Calculator_H2O_CO2_NaCl<1U>;
template class PVTX_Calculator_H2O_CO2_NaCl<2U>;
template class PVTX_Calculator_H2O_CO2_NaCl<3U>;
 

} // end csmp
