#ifndef IMPES_SIMULATOR_KEYS_H
#define IMPES_SIMULATOR_KEYS_H

//#define FECFV_SIMULATOR_EXTRA_MONITORING

#include "Index.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;

template<uint32_t dim>
struct IMPES_SimulatorKeys
{
    IMPES_SimulatorKeys(const PropertyDatabase<dim>&);
    
    IMPES_SimulatorKeys();
    
    ~IMPES_SimulatorKeys();
    
    void ConstructKeys(const PropertyDatabase<dim>&);

    csmp::Index  pf_key; // fluid pressure
    csmp::Index  vt_key; // Darcy total velocity
    csmp::Index  vo_key; // Darcy oil velocity
    csmp::Index  qf_key; // fluid volume flux
    csmp::Index  qw_key; // water volume flux
    csmp::Index  qo_key; // oil volume flux
    csmp::Index  sw_key; // water saturation
    csmp::Index  so_key; // oil saturation
    csmp::Index  mw_key; // viscosity water
    csmp::Index  mo_key; // oil viscosity
    csmp::Index  rw_key; // density water
    csmp::Index  ro_key; // density oil
    csmp::Index  gt_key; // gravity term
    csmp::Index  k_key;  // permeability
    csmp::Index  tm_key; // total mobility
    csmp::Index  ct_key; // capillary term
	  csmp::Index  pc_key; // capillary pressure
    csmp::Index  ns_key; // nodal source
    csmp::Index  phi_key; // porosity
    csmp::Index  vm_key; // volume modifier
};


} // end csmp

#endif
