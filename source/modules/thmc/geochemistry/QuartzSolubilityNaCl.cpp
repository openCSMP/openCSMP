// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "QuartzSolubilityNaCl.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
QuartzSolubilityNaCl<dim>::QuartzSolubilityNaCl( const PropertyDatabase<dim>& p ) 
      : Interrelation<dim>(p), 
        T( Interrelation<dim>::GlobalProperty("temperature") ),
        R( Interrelation<dim>::GlobalProperty("nodal fluid density") ),
        X( Interrelation<dim>::GlobalProperty("mass fraction NaCl") ),
        S( Interrelation<dim>::GlobalProperty("quartz solubility") ),
        h(0.0) // cation hydration number
 {
    Interrelation<dim>::Name("QuartzSolubilityNaCl");
    Interrelation<dim>::OutputCondition( S, PLAIN );
    Interrelation<dim>::ResultProperty("quartz solubility");
 }


template<uint32_t dim>
void QuartzSolubilityNaCl<dim>::Calculate()
 {
    // 1. getting input data for the calculation 
    // -----------------------------------------
    T.AssignTo( tC );
    R.AssignTo( rho );     // kg m-3
    K    = tC() + 273.15;  // conversion to temperature in Kelvin

    // 3. Acounting for the presence of NaCl and a cation hydration number, h of 0
    //    where Xaq=weight fraction water and Xsalt is weight fraction salt
    // --------------------------------------------
    X.AssignTo( x_salt );
    // calculating weight fraction water, F, from salt weight fraction
    F = 1.0 - x_salt();
    //                            conversion of molality to mole fraction
    rho_e = rho()*1.0e-3 * F * (1.0 - h*((x_salt()/58.44)/55.51));
    
    logV = log10( 1.0 / rho_e ); // density (g cm-3) -> log base 10 of specific volume (cm3 g-1)
    
     // 2. calculating quartz solubility a la Fournier and Potter GCA V.46 P.1969 1982.
    //    (g silica/ g H2O) -> to kg / kg
    // -------------------------------------------------------------------------------
    A = -4.66206 + 0.0034063*K + 2179.7/K - 1.1292e+6/(K*K) + 1.3543e+8/(K*K*K);
    B = -0.001418*K - 806.97/K;
    C =  3.9465e-4*K;   
 //       m    = A + B*logV + C*(logV*logV);  // result in log molality -> kg SiO2 kg-1 solution (m.w. 60.0843)
    m = pow( 10.0, A + B*logV + C*(logV*logV) );   
 
//    cout <<"\nQuartzSolubility::Calculate: Fournier (1983) quartz solubility (molality): "<< m << endl;
    
    // 4. assigning quartz solubility as kg silica / kg fluid
    // (converting from molality to kg kg-1 of saline fluid
    //  mole weight SiO2 = 60.0843 g). 
    // ------------------------------------------------------
    S = m * 0.0600843;
//    S.PrintValue();
} 


template class QuartzSolubilityNaCl<1U>;
template class QuartzSolubilityNaCl<2U>;
template class QuartzSolubilityNaCl<3U>;

} // csmp
