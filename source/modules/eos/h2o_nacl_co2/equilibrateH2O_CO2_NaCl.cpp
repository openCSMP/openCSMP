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
  
EOS_CO2H2ONaCl_Spycher04  csmp_eos;


/**
    Applies PVTX_Calculator_H2O_CO2_NaCl to nodes of region
*/
template<size_t dim>
void equilibrateFluid( Region<dim>& model_subdomain, PVTX_Calculator_H2O_CO2_NaCl<dim>& pvtx_calculator, double64 delta_t )
 {
    for ( auto nit=model_subdomain.NodesBegin(); nit!=model_subdomain.NodesEnd(); ++nit )
      pvtx_calculator.Equilibrate( *nit, delta_t );
 }

template void equilibrateFluid<1U>( Region<1U>&, PVTX_Calculator_H2O_CO2_NaCl<1U>&, double64 );
template void equilibrateFluid<2U>( Region<2U>&, PVTX_Calculator_H2O_CO2_NaCl<2U>&, double64 );
template void equilibrateFluid<3U>( Region<3U>&, PVTX_Calculator_H2O_CO2_NaCl<3U>&, double64 );




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
    update PTX based fluid properties (density, viscosity, compressibility, nodal mass source)
*/
template<size_t dim>
void updateFluidProperties(const variables::VariableSet_CO2GeoSequestration& props, Node<dim>& n, double64 del_t)
{    
    const double64 rhoNaCl(2170.); // density = 2170 kg/m3
     
    ArrayVariable  H2Ocomposition(3), CO2composition(2);
    n.Read( props.key_H2O_comp, H2Ocomposition ); 
    n.Read( props.key_CO2_comp, CO2composition ); 
       
    // saturation
    double64 old_sCO2 = n.Read( props.key_sCO2 );
    // phase densities
    double64 old_rhoH2O = n.Read( props.key_rhoH2O );
    double64 old_rhoCO2 = n.Read( props.key_rhoCO2 );
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
/*
#ifdef CSMP_DEBUG_aqueousPhaseDensity_H2O_CO2_NaCl
      double64 old_XCO2(H2Ocomposition(XCO2)), old_XH2O(H2Ocomposition(XH2O)), old_XNACl_aq(H2Ocomposition(XNACl_aq)),
               old_YCO2(CO2composition(YCO2)), old_YH2O(CO2composition(YH2O));
#endif
*/        
    const double64 Pf(n.Read( props.key_pf ));
    const double64 TC(n.Read( props.key_T )); // in centrigrade
    const double64 molalityNaCl( csmp_eos.massFracNaClToMolalNaClInAqueousPhase(H2Ocomposition(XNACl_aq)* 100.) );

    const double64 rhoCO2 = csmp_eos.Rho_CarbonicPhase( Pf, TC );
    if ( n.Status(props.key_rhoCO2) != DIRICH )
        n.Store( props.key_rhoCO2, makeScalar(n.Status(props.key_rhoCO2),rhoCO2) );
        

    const double64 rho_brine_CO2 = csmp_eos.Rho_AqueousPhase( Pf,TC, molalityNaCl );  
    if ( n.Status(props.key_rhoH2O) != DIRICH )
       n.Store( props.key_rhoH2O, makeScalar(n.Status(props.key_rhoH2O),rho_brine_CO2) );
           
    const double64 muH2O = csmp_eos.mu_AqueousPhase( Pf, TC, molalityNaCl );
    const double64 muCO2 = csmp_eos.mu_CarbonicPhase( Pf, TC );
    if ( n.Status(props.key_muH2O) != DIRICH )
        n.Store( props.key_muH2O, makeScalar(n.Status(props.key_muH2O),muH2O) );
    if ( n.Status(props.key_muCO2) != DIRICH )
        n.Store( props.key_muCO2, makeScalar(n.Status(props.key_muCO2),muCO2) );

    const double64 betaH2O = csmp_eos.C_AqueousPhase( Pf, TC, molalityNaCl );
    const double64 betaCO2 = csmp_eos.C_CarbonicPhase( Pf, TC );
    if ( n.Status(props.key_cH2O) != DIRICH )
        n.Store( props.key_cH2O, makeScalar(n.Status(props.key_cH2O),betaH2O) );
    if ( n.Status(props.key_cCO2) != DIRICH )
        n.Store( props.key_cCO2, makeScalar(n.Status(props.key_cCO2),betaCO2) );

    const double64 sCO2new = n.Read(props.key_sCO2);
    const double64 sH2Onew = n.Read(props.key_sH2O);
    // new water mass
    const double64 new_mass_H2O = H2Ocomposition(XH2O) * rho_brine_CO2 * sH2Onew  +
                                  CO2composition(YH2O) * rhoCO2 * sCO2new;
    // new carbondioxide mass
    const double64 new_mass_CO2 = H2Ocomposition(XCO2) * rho_brine_CO2 * sH2Onew  +
                                  CO2composition(YCO2) * rhoCO2 * sCO2new;
    // new dissolved salt mass
    const double64 new_mass_NaClaq = H2Ocomposition(XNACl_aq) * rho_brine_CO2 * sH2Onew;
    // new total mass
    const double64 new_total_mass = new_mass_H2O + new_mass_CO2 + new_mass_NaClaq;
    // new bulk density
    const double64 new_bulk_density(rho_brine_CO2 * sH2Onew + rhoCO2 * sCO2new);
    // new volume (mass/density)
    const double64 new_fluid_volume = new_total_mass / new_bulk_density;

    // total mass change
    double64 nodal_mass_source = new_total_mass - total_mass;
    // expansion factor
    const double64 expansion_factor = 1. - (nodal_mass_source / new_bulk_density) / new_fluid_volume;
    // store nodal fluid mass source for pressure computation
    const double64 PV = n.Read(props.key_fvPV);
    //if (del_t <= 0.) nodal_mass_source =0.;
    //else nodal_mass_source *= (PV/del_t);
    
    if (del_t <= 0.) {
        nodal_mass_source =0.;
    } else {
        const double64 drhof_dt = ( new_bulk_density - bulk_density ) / del_t;
        nodal_mass_source = drhof_dt * PV;
    }
    
    if (n.Status(props.key_nQM) != DIRICH)
        n.Store( props.key_nQM, makeScalar(n.Status(props.key_nQM),nodal_mass_source) );

    // solid/precipitated salt mass
    const double64 V_salt(n.Read(props.key_NaCl));
    const double64 mass_NaClsd = V_salt * PV * rhoNaCl;       
    
    // store mass
    n.Store( props.key_mH2O, makeScalar(n.Status(props.key_mH2O),new_mass_H2O*PV) );
    n.Store( props.key_mCO2, makeScalar(n.Status(props.key_mCO2),new_mass_CO2*PV) );
    n.Store( props.key_mNaClaq, makeScalar(n.Status(props.key_mNaClaq),new_mass_NaClaq*PV) );
    n.Store( props.key_mNaClsd, makeScalar(n.Status(props.key_mNaClsd),mass_NaClsd) );
        
} 

template void updateFluidProperties<1U>( const variables::VariableSet_CO2GeoSequestration&, Node<1U>&, double64 ); 
template void updateFluidProperties<2U>( const variables::VariableSet_CO2GeoSequestration&, Node<2U>&, double64 );
template void updateFluidProperties<3U>( const variables::VariableSet_CO2GeoSequestration&, Node<3U>&, double64 );



/**
    update PTX based fluid properties (density, viscosity, compressibility) for all model nodes
*/
template<size_t dim>
void updatePTXBasedFluidProperties(const variables::VariableSet_CO2GeoSequestration& props, Region<dim>& model_subdomain, double64 del_t)
{
    
    for ( auto nit=model_subdomain.NodesBegin(); nit!=model_subdomain.NodesEnd(); ++nit )
        updateFluidProperties( props, *(*nit), del_t);
}

template void updatePTXBasedFluidProperties<1U>( const variables::VariableSet_CO2GeoSequestration&, Region<1U>&, double64 ); 
template void updatePTXBasedFluidProperties<2U>( const variables::VariableSet_CO2GeoSequestration&, Region<2U>&, double64 );
template void updatePTXBasedFluidProperties<3U>( const variables::VariableSet_CO2GeoSequestration&, Region<3U>&, double64 );
 
 
} // end csmp


