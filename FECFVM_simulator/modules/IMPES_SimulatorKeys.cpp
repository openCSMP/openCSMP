// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IMPES_SimulatorKeys.h"
#include "PropertyDatabase.h"
#include "Exception.h"


namespace csmp {

template<uint32_t dim>
IMPES_SimulatorKeys<dim>::IMPES_SimulatorKeys( const PropertyDatabase<dim>& p )
 :  pf_key(p.StorageKey("fluid pressure")),            // fluid pressure
    vt_key(p.StorageKey("velocity")),                  // Darcy total velocity
    vo_key(p.StorageKey("velocity oil")),              // Darcy oil velocity
    qf_key(p.StorageKey("fluid volume source")),       // fluid volume flux
    qw_key(p.StorageKey("water volume source")),       // water volume flux
    qo_key(p.StorageKey("oil volume source")),         // oil volume flux
    sw_key(p.StorageKey("saturation water")),          // water saturation
    so_key(p.StorageKey("saturation oil")),            // oil saturation
    mw_key(p.StorageKey("viscosity water")),           // water viscosity
    mo_key(p.StorageKey("viscosity oil")),             // oil viscosity
    rw_key(p.StorageKey("density water")),             // density water
    ro_key(p.StorageKey("density oil")),               // density oil
    gt_key(p.StorageKey("gravity term")),              // gravity term
    k_key(p.StorageKey("permeability")),               // permeability
    tm_key(p.StorageKey("total mobility")),            // total mobility
    ct_key(p.StorageKey("capillary term")),            // capillary term

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING
    pc_key(p.StorageKey("capillary pressure")),        // capillary pressure
#endif

    ns_key(p.StorageKey("nodal fluid volume source")), // nodal source
    phi_key(p.StorageKey("porosity")),                 // porosity
    vm_key(p.StorageKey("volume modifier"))            // volume modifier
{
     if ( pf_key.place != NODE || pf_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'fluid pressure' variable must be a scalar placed on the node" );

     if ( vt_key.place != ELEMENT || vt_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'velocity' variable must be a vector placed on the element" );

     if ( vo_key.place != ELEMENT || vo_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'velocity oil' variable must be a vector placed on the element" );

     if ( qf_key.place != ELEMENT || qw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
         "The 'fluid volume source' must be a scalar element variable" );

     if ( qw_key.place != ELEMENT || qw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
         "The 'water volume source' must be a scalar element variable" );

     if ( qo_key.place != ELEMENT || qo_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
         "The 'oil volume source' must also be a scalar element variable" );

     if ( sw_key.place != ELEMENT || sw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'saturation water' variable must be scalar element property" );

     if ( so_key.place != ELEMENT || so_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'saturation oil' variable must be scalar element property" );

     if ( mo_key.place != ELEMENT || mo_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'viscosity oil' must be scalar element property" );

     if ( mw_key.place != ELEMENT || mw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'viscosity water' must be scalar element property" );

     if ( ro_key.place != ELEMENT || ro_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'density oil' must be scalar element property" );

     if ( rw_key.place != ELEMENT || rw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'density water' must be scalar element property" );

     if ( gt_key.place != ELEMENT || gt_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
         "The 'gravity term' must be a vector element variable" );

    if ( k_key.place != ELEMENT )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'permeability' variable must be an element property" );

    if ( tm_key.place != ELEMENT || tm_key.type != SCALAR  )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'total mobility' variable must be a scalar element property" );

    if ( ct_key.place != ELEMENT || ct_key.type != VECTOR  )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'capillary term' variable must be a vector element property" );

#ifdef FECFV_SIMULATOR_EXTRA_MONITORING
    if ( pc_key.place != ELEMENT || pc_key.type != SCALAR)
        throw csmp::Exception(FATAL_ERROR,
                              "IMPES_SimulatorKeys(constructor)",
                              "The 'capillary pressure' variable must be a scalar element property");
#endif

    if ( (ns_key.place != NODE || ns_key.type != SCALAR) )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
                      "The 'nodal fluid volume source' variable must be scalar element property" );

     if ( phi_key.place != ELEMENT || phi_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'porosity' variable must be scalar element property" );

     if ( vm_key.place != ELEMENT || vm_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'volume modifier' variable must be scalar element property" );

}

template<uint32_t dim>
IMPES_SimulatorKeys<dim>::IMPES_SimulatorKeys()
{
}


template<uint32_t dim>
IMPES_SimulatorKeys<dim>::~IMPES_SimulatorKeys()
{
}


template<uint32_t dim>
void IMPES_SimulatorKeys<dim>::ConstructKeys( const PropertyDatabase<dim>& p )
{
    pf_key = p.StorageKey("fluid pressure");            // fluid pressure
    vt_key = p.StorageKey("velocity");                  // Darcy total velocity
    vo_key = p.StorageKey("velocity oil");              // Darcy oil velocity
    qf_key = p.StorageKey("fluid volume source");       // fluid volume flux
    qw_key = p.StorageKey("water volume source");       // water volume flux
    qo_key = p.StorageKey("oil volume source");         // oil volume flux
    sw_key = p.StorageKey("saturation water");          // water saturation
    so_key = p.StorageKey("saturation oil");            // oil saturation
    mw_key = p.StorageKey("viscosity water");           // water viscosity
    mo_key = p.StorageKey("viscosity oil");             // oil viscosity
    rw_key = p.StorageKey("density water");             // density water
    ro_key = p.StorageKey("density oil");               // density oil
    gt_key = p.StorageKey("gravity term");              // gravity term
    k_key = p.StorageKey("permeability");               // permeability
    tm_key = p.StorageKey("total mobility");            // total mobility
    ct_key = p.StorageKey("capillary term");            // capillary term
	pc_key = p.StorageKey("capillary pressure");        // capillary pressure
    ns_key = p.StorageKey("nodal fluid volume source"); // nodal source
    phi_key = p.StorageKey("porosity");                 // porosity
    vm_key = p.StorageKey("volume modifier");           // volume modifier
    
    
     if ( pf_key.place != NODE || pf_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'fluid pressure' variable must be a scalar placed on the node" );

     if ( vt_key.place != ELEMENT || vt_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'velocity' variable must be a vector placed on the element" );

     if ( vo_key.place != ELEMENT || vo_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'velocity oil' variable must be a vector placed on the element" );

     if ( qf_key.place != ELEMENT || qw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
         "The 'fluid volume source' must be a scalar element variable" );

     if ( qw_key.place != ELEMENT || qw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
         "The 'water volume source' must be a scalar element variable" );

     if ( qo_key.place != ELEMENT || qo_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
         "The 'oil volume source' must also be a scalar element variable" );

     if ( sw_key.place != ELEMENT || sw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'saturation water' variable must be scalar element property" );

     if ( so_key.place != ELEMENT || so_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'saturation oil' variable must be scalar element property" );

     if ( mo_key.place != ELEMENT || mo_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'viscosity oil' must be scalar element property" );

     if ( mw_key.place != ELEMENT || mw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'viscosity water' must be scalar element property" );

     if ( ro_key.place != ELEMENT || ro_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'density oil' must be scalar element property" );

     if ( rw_key.place != ELEMENT || rw_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'density water' must be scalar element property" );

     if ( gt_key.place != ELEMENT || gt_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
         "The 'gravity term' must be a vector element variable" );

    if ( k_key.place != ELEMENT )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'permeability' variable must be an element property" );

    if ( tm_key.place != ELEMENT || tm_key.type != SCALAR  )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'total mobility' variable must be a scalar element property" );

    if ( ct_key.place != ELEMENT || ct_key.type != VECTOR  )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'capillary term' variable must be a vector element property" );

	if (pc_key.place != ELEMENT || pc_key.type != SCALAR)
		throw csmp::Exception(FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
		"The 'capillary pressure' variable must be a scalar element property");

    if ( (ns_key.place != NODE || ns_key.type != SCALAR) )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
                      "The 'nodal fluid volume source' variable must be scalar element property" );

     if ( phi_key.place != ELEMENT || phi_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys<dim>::ConstructKeys",
          "The 'porosity' variable must be scalar element property" );

     if ( vm_key.place != ELEMENT || vm_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "IMPES_SimulatorKeys(constructor)",
          "The 'volume modifier' variable must be scalar element property" );

}

template class IMPES_SimulatorKeys<2>;

} // end csmp
