#include "Reactant.h"

using namespace std;

namespace csmp {

Reactant::Reactant()
  :  name("unspecified"),
     phase_state(1), // aqueous 
     charge(0),
     density(0.0),
     molar_weight(0.0)
 {
 }
 
 
Reactant::Reactant( const Reactant& r )
 {
    *this = r;
 }
 
 
Reactant::Reactant( const char* n, short st, int ch, double dens, double mw )
  :  name(n),
     phase_state(st), // aqueous 
     charge(ch),
     density(dens),
     molar_weight(mw)
 {
    // solids shall not be charged 
    if ( phase_state == 1 ) charge = 0;
    // aqueous species densities are unlikely measured
    if ( phase_state == 2 ) density = 0.0;
 }
 
  
Reactant&  Reactant::operator=( const Reactant& r )
 {
    if ( &r != this )
      {
         name         = r.name;
         phase_state  = r.phase_state; 
         charge       = r.charge;
         density      = r.density;
         molar_weight = r.molar_weight;
      }
    return *this;
 }
 
 


double  Reactant::VolumeFractionToMolality( double fraction, double fluid_mass ) const
 {
    if ( phase_state != 1 )
      {
         cout <<"\nReactant::VolumeFractionToMolality: ";
         cout <<"Only volume fractions of solids are constrained. ";
         cout <<"\nNothing was done."<< endl;
         return fraction;
      }
    //                                kg -> g 
    else return (fraction * density * 1.0e+3) / (molar_weight * fluid_mass);
 }
 

double  Reactant::MolalityToVolumeFraction( double molality, double fluid_mass ) const
 {
    if ( phase_state != 1 )
      {
         cout <<"\nReactant::MolalityToVolumeFraction: ";
         cout <<"Can only convert solids to volume fractions. ";
         cout <<"\nNothing was done."<< endl;
         return molality;
      }
    //                                     g -> kg  
    else return (molality * molar_weight * 1.0e-3 * fluid_mass) / density;
 }
 
 
 
 
 
/// converts moles of aqueous species per m3 fluid into molalites (moles per kg fluid)
double  Reactant::MolesPerVolumeToMolality( double moles, double fluid_dens ) const
 {
    if ( phase_state != 2 )
      {
         cout <<"\nReactant::MolesPerVolumeToMolality: ";
         cout <<"Only aqueous species are targetted by this treatment. ";
         cout <<"\nNothing was done."<< endl;
         return moles;
      }
    else return moles / fluid_dens;
 }
 

double  Reactant::MolalityToMolesPerVolume( double molality, double fluid_dens ) const
 {
    if ( phase_state != 2 )
      {
         cout <<"\nReactant::MolalityToMolesPerVolume: ";
         cout <<"Only aqueous species are targetted by this treatment. ";
         cout <<"\nNothing was done."<< endl;
         return molality;
      }
    else return molality * fluid_dens;
 }
 


void  Reactant::Out() const
 {
    cout <<"\nReactant::Out: Printing private data:"<< endl;
    cout <<"Reactant name: "<< name << endl;
    cout <<"Phase state:   ";
    if      ( phase_state == 1 ) cout <<"solid"   << endl;
    else if ( phase_state == 2 ) cout <<"aqueous" << endl;
    else if ( phase_state == 3 ) cout <<"gaseous" << endl;
    cout <<"Charge:        "<< charge << endl;
    cout <<"Density:       "<< density << endl;
    cout <<"Molar weight:  "<< molar_weight << endl;
    cout << endl;
 }

} // csp 
 
 
 
 
 
