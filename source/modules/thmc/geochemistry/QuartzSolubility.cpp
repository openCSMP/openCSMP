#include "QuartzSolubility.h"

using namespace std;

namespace csmp {

template<size_t dim>
QuartzSolubility<dim>::QuartzSolubility( const PropertyDatabase<dim>& p ) 
      : Interrelation<dim>(p), 
        T( Interrelation<dim>::GlobalProperty("temperature") ),
        R( Interrelation<dim>::GlobalProperty("nodal fluid density") ),
        S( Interrelation<dim>::GlobalProperty("quartz solubility") )
 {
    Interrelation<dim>::Name("QuartzSolubility");
    Interrelation<dim>::OutputCondition( S, PLAIN );
    Interrelation<dim>::ResultProperty("quartz solubility");
 }


template<size_t dim>
void QuartzSolubility<dim>::Calculate()
 {
     // 1. getting input data for the calculation 
     // -----------------------------------------
     T.AssignTo( sc_tC );
     R.AssignTo( rho );
     tC = sc_tC();
 
     // converting density: kg m-3 -> g cm-3
     v = 1.0 / (1.0e-3 * rho());
     
    // 3. Craig Mannings fit of quartz solubility
    // ------------------------------------------
    sol  = 4.262 - 5764.2/tC + 1751300.0/tC/tC - 2.2869e+8/tC/tC/tC + (2.8454 - 1006.9/tC + 356890.0/tC/tC); 
	  sol *= log10(v);
    sol  = pow( 10.0, sol ) * 0.060048;

//    cout <<"\nQuartzSolubility::Calculate: Craig Manning's quartz solubility (g silica /  kgH2O): "<< sol << endl;
    
    // 4. assigning quartz solubility as kg silica / kg fluid
    // ------------------------------------------------------
    S = sol;
 } 
 
 
 
template class QuartzSolubility<1U>; 
template class QuartzSolubility<2U>; 
template class QuartzSolubility<3U>; 
 
} // csmp
