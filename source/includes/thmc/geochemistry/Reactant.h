#ifndef CSMP_REACTANT_H
#define CSMP_REACTANT_H

#include "CSMP_definitions.h"

namespace csmp {

class Reactant {
    std::string  name;
    short   phase_state; // 1=solid, 2=aqueous, 3=gaseous
    int     charge;
    double  density;
    double  molar_weight;
    
 public:
    Reactant();
    Reactant( const Reactant& r );
    Reactant( const char* n, short st, int ch, double dens, double mw ); 
    Reactant& operator=( const Reactant& r );
    std::string Name()                const { return name; };
    void   Name( const char* s )       { name=s; };
    void   Name( const std::string& s ){ name=s; };
    short  State()               const { return phase_state; };
    void   State( short s )            { phase_state=s; };
    int    Charge()              const { return charge; };
    void   Charge( int ch )            { charge=ch; };
    double Density()             const { return density; };
    void   Density( double d )         { density=d; };      
    double MolarWeight()         const { return molar_weight; };
    void   MolarWeight( double mw )    { molar_weight=mw; };
    // conversions
    double VolumeFractionToMolality( double fraction, double fluid_mass ) const;
    double MolalityToVolumeFraction( double molality, double fluid_mass ) const;
    double MolesPerVolumeToMolality( double moles,    double fluid_dens ) const;
    double MolalityToMolesPerVolume( double molality, double fluid_dens ) const;
    void   Out() const;
};

} // csmp

#endif
