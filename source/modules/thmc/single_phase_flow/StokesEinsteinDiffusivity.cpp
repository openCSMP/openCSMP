#include "StokesEinsteinDiffusivity.h"

using namespace std;

namespace csmp {

template<size_t dim>
StokesEinsteinDiffusivity<dim>::StokesEinsteinDiffusivity( 
                                      const PropertyDatabase<dim>& p, 
                                      const char* species_diffusivity, 
                                      double64      diffusivity, 
                                      double64      Omega ) 
      : Interrelation<dim>(p),
        T( Interrelation<dim>::GlobalProperty("temperature") ),
        E( Interrelation<dim>::GlobalProperty("nodal fluid viscosity") ),
        D( Interrelation<dim>::GlobalProperty(species_diffusivity) ),
        solute(species_diffusivity),
        Bo(1.380662e-23),
        diffus(diffusivity),
        sixPi(18.849556),
        omega(Omega)
 {
    Interrelation<dim>::Name("StokesEinsteinDiffusivity");
    Interrelation<dim>::OutputCondition( D, PLAIN );
    Interrelation<dim>::ResultProperty(species_diffusivity);
    
    // calculating R0 from the Stokes-Einstein relationship and the given values
    // assuming that the diffusivity is quoted at 25oC and that the viscosity of
    // water at this temperature is 1.6e-5 Pa s-1
    R0 = (Bo * 25.0) / (diffus * sixPi * 1.6e-5);
    cout <<"\nStokesEinsteinDiffusivity (ctor): Solute radius calculated for '"<< solute;
    cout <<"' from the Stokes-Einstein rln (nano-m): "<< (R0*1.0e+9) << endl;
 } 



/// Calculating the diffusivity of the target species from temperature,
/// fluid viscosity and R0
template<size_t dim>
void StokesEinsteinDiffusivity<dim>::Calculate()
 {
    T.AssignTo( Tc );
    E.AssignTo( eta );
    
    D  = (Bo * Tc()) / (sixPi * eta() * R0);
    D *= omega;
 } 


template class StokesEinsteinDiffusivity<1U>;
template class StokesEinsteinDiffusivity<2U>;
template class StokesEinsteinDiffusivity<3U>;

} // csmp
