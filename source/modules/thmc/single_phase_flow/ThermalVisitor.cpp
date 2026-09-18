// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ThermalVisitor.h"
#include "Model.h"
#include "Region.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ThermalVisitor<dim>::ThermalVisitor( Model<dim>& model, std::vector<std::string>* to_initialize_keys):
      Visitor<dim>( MODEL, NODE ),
      model_ (model),
      to_initialize_keys_(to_initialize_keys)
{
    if (!to_initialize_keys_){
        nCPT_key_=           model.Database().StorageKey( "nodal total heat capacity"        ) ;
        cp_fluid_key_=       model.Database().StorageKey( "fluid heat capacity"              ) ;
        cp_rock_key_=        model.Database().StorageKey( "nodal heat capacity rock"         ) ;
        rho_rock_key_=       model.Database().StorageKey( "nodal density rock"               ) ;
        t_key_=              model.Database().StorageKey( "temperature"                      ) ;
        p_key_=              model.Database().StorageKey( "fluid pressure"                   ) ;
        hl_key_=             model.Database().StorageKey( "enthalpy liquid"                  ) ;
        hCl_key_=            model.Database().StorageKey( "enthalpy content liquid"          ) ;
        hClp_key_=           model.Database().StorageKey( "previous enthalpy content liquid" ) ;
        hVl_key_=            model.Database().StorageKey( "volumetric enthalpy liquid"       ) ;
        mt_key_=             model.Database().StorageKey( "fluid density"                    ) ;
        rl_key_=             model.Database().StorageKey( "density liquid"                   ) ;
        mu_key_=             model.Database().StorageKey( "fluid viscosity"                  ) ;
        nphi_key_=           model.Database().StorageKey( "nodal porosity"                   ) ;
        nQ_key_=             model.Database().StorageKey( "nodal fluid volume source"        ) ;
        dRho_key_=           model.Database().StorageKey( "density difference"               ) ;
        volume_factor_key_=  model.Database().StorageKey( "volume factor"                    ) ;
        beta_rock_key_=      model.Database().StorageKey( "nodal compressibility rock"       ) ;
        beta_tot_key_=       model.Database().StorageKey( "nodal total compressibility"      ) ;
        beta_fluid_key_=     model.Database().StorageKey( "fluid compressibility"            ) ;
    }
    else if (to_initialize_keys_->size()==20){
        nCPT_key_=           model.Database().StorageKey( to_initialize_keys_->operator [](0).c_str()) ;
        cp_fluid_key_=       model.Database().StorageKey( to_initialize_keys_->operator [](1).c_str()) ;
        cp_rock_key_=        model.Database().StorageKey( to_initialize_keys_->operator [](2).c_str()) ;
        rho_rock_key_=       model.Database().StorageKey( to_initialize_keys_->operator [](3).c_str()) ;
        t_key_=              model.Database().StorageKey( to_initialize_keys_->operator [](4).c_str()) ;
        p_key_=              model.Database().StorageKey( to_initialize_keys_->operator [](5).c_str()) ;
        hl_key_=             model.Database().StorageKey( to_initialize_keys_->operator [](6).c_str()) ;
        hCl_key_=            model.Database().StorageKey( to_initialize_keys_->operator [](7).c_str()) ;
        hClp_key_=           model.Database().StorageKey( to_initialize_keys_->operator [](8).c_str()) ;
        hVl_key_=            model.Database().StorageKey( to_initialize_keys_->operator [](9).c_str()) ;
        mt_key_=             model.Database().StorageKey( to_initialize_keys_->operator [](10).c_str()) ;
        rl_key_=             model.Database().StorageKey( to_initialize_keys_->operator [](11).c_str()) ;
        mu_key_=             model.Database().StorageKey( to_initialize_keys_->operator [](12).c_str()) ;
        nphi_key_=           model.Database().StorageKey( to_initialize_keys_->operator [](13).c_str()) ;
        nQ_key_=             model.Database().StorageKey( to_initialize_keys_->operator [](14).c_str()) ;
        dRho_key_=           model.Database().StorageKey( to_initialize_keys_->operator [](15).c_str()) ;
        volume_factor_key_=  model.Database().StorageKey( to_initialize_keys_->operator [](16).c_str()) ;
        beta_rock_key_=      model.Database().StorageKey( to_initialize_keys_->operator [](17).c_str()) ;
        beta_tot_key_=       model.Database().StorageKey( to_initialize_keys_->operator [](18).c_str()) ;
        beta_fluid_key_=     model.Database().StorageKey( to_initialize_keys_->operator [](19).c_str()) ;
    }
    else
    {
        throw csmp::Exception(FATAL_ERROR,"ThermalVisitor(constructor)"," Incorrect number of variable names provided"," Please use the correct number.");
    }
  // Variable initialization
  // 
  // "heat capacity rock" should have been defined and assigned a value before construction of this visitor.
  //  As other code parts might need it, e.g., the thermal diffusion algorithm, I thought it might be better
  //  not to do it here.

  cout <<"\nThermalVisitor Constructor: Initialised"<< endl;
}

template<uint32_t dim>
ThermalVisitor<dim>::ThermalVisitor( Model<dim>& model,
                                     Index nCPT_key,
                                     Index cp_fluid_key,
                                     Index cp_rock_key,
                                     Index rho_rock_key,
                                     Index t_key,
                                     Index p_key,
                                     Index hl_key,
                                     Index hCl_key,
                                     Index hClp_key,
                                     Index hVl_key,
                                     Index mt_key,
                                     Index rl_key,
                                     Index mu_key,
                                     Index nphi_key,
                                     Index nQ_key,
                                     Index dRho_key,
                                     Index volume_factor_key,
                                     Index beta_rock_key,
                                     Index beta_tot_key,
                                     Index beta_fluid_key):
      Visitor<dim>( MODEL, NODE ),
      model_ (model),
      nCPT_key_(nCPT_key),
      cp_fluid_key_(cp_fluid_key),
      cp_rock_key_(cp_rock_key),
      rho_rock_key_(rho_rock_key),
      t_key_(t_key),
      p_key_(p_key),
      hl_key_(hl_key),
      hCl_key_(hCl_key),
      hClp_key_(hClp_key),
      hVl_key_(hVl_key),
      mt_key_(mt_key),
      rl_key_(rl_key),
      mu_key_(mu_key),
      nphi_key_(nphi_key),
      nQ_key_(nQ_key),
      dRho_key_(dRho_key),
      volume_factor_key_(volume_factor_key),
      beta_rock_key_(beta_rock_key),
      beta_tot_key_(beta_tot_key),
      beta_fluid_key_(beta_fluid_key)
{
  // Variable initialization
  //
  // "heat capacity rock" should have been defined and assigned a value before construction of this visitor.
  //  As other code parts might need it, e.g., the thermal diffusion algorithm, I thought it might be better
  //  not to do it here.

  cout <<"\nThermalVisitor Constructor: Initialised"<< endl;
}



/**
    This is the key part of this Node visitor.
*/
template<uint32_t dim>
void ThermalVisitor<dim>::Visit( Node<dim>* n )
{
  // fluid properties need to be updated if T=DIRICH and pressure is not
  if ((n->Status(t_key_) != DIRICH) || (n->Status(p_key_) != DIRICH))
    {
      // 1. Read relevant variables
      n->Read( cp_fluid_key_, cp_fluid );
      n->Read( cp_rock_key_,  cp_rock );
      n->Read( rho_rock_key_, rho_rock ); // needs to be updated in GEMS-visitor!!!!!!!!!!!!!!!!!!!!
      n->Read( t_key_,        t        ); // Celsius, from CSMP
      n->Read( p_key_,        p        ); // Celsius, from CSMP
      n->Read( hCl_key_,      hCl      );
      n->Read( hClp_key_,     hClp     );
      n->Read( mt_key_,       mt       );
      n->Read( rl_key_,       rl       );
      n->Read( nphi_key_,     phi      ); // is updated from GEMSVisitor
      n->Read( beta_rock_key_, beta_rock);
      //n->Read( finite_volume_key, finite_volume );

      // 2. Compute values of transient properties
      pore_volume    = phi();
      rock_volume    = 1.0 - phi();

      //after temperature diffusion cp_fluid has changed, total heat capacity has changed
      cp_fluid()     = water.HeatCapacity (t(), p());
      nCPT()         = (rock_volume*cp_rock()*rho_rock() + pore_volume*cp_fluid()*rl());  
      
      // 3. Compute enthalpy change and convert to delta T
      // temperature change due to enthalpy advection
      if (n->Status(t_key_) != DIRICH)
        {
              dhCl_          = pore_volume * (hCl() - hClp());
              dT             = dhCl_/nCPT();

      // 4. Correct temperature
              t()           += dT;
        }

      // density liquid
      rl()           = water.Density(t(),p());
      // fluid density - mt - is changed only in advection
      // fluid viscosity
      mu()           = water.Viscosity(t(), p());
      
      //! for wells
      if (n->Status(mt_key_) == DIRICH)
        {
          mt() = rl();   
          n->Store( mt_key_, mt );
        }

      // simple check, should be close to 1
      volume_factor()= mt()/rl();
    
      // density difference is used later to calculate the
      // source term that goes into the pressure equation
      dRho()           = (mt()-rl()); 
      nQ = 0;
      
      // enthalpy liquid 
      hl()           = water.Enthalpy(t(),p());
      // enthalpy content liquid - advected variable
      hCl()          = hl()*mt(); 
      // volumetric enthalpy liquid
      hVl()          = hl()*rl(); 
      
      //other formulation is possible, but in general beta_rock << beta_fluid anyway
      beta_tot()     = mt()*phi()*(water.Compressibility(t(),p()) + beta_rock()); 
    
      // addition Sebastian and Alessandro - output fluid compressibility
      beta_fluid     = water.Compressibility(t(),p());
      
      n->Store( t_key_, t );
      n->Store( cp_fluid_key_, cp_fluid);
      n->Store( nCPT_key_, nCPT );
      n->Store( hCl_key_, hCl);
      n->Store( hClp_key_, hCl);
      n->Store( hVl_key_, hVl);
      n->Store( rl_key_,  rl );
      n->Store( mu_key_,  mu );
      n->Store( nQ_key_,  nQ );
      n->Store( dRho_key_,  dRho );
      n->Store( volume_factor_key_,  volume_factor );
      n->Store( beta_tot_key_,  beta_tot );
      n->Store( beta_fluid_key_,  beta_fluid );
    
    }
  
} // end Visit(node)


template<uint32_t dim>
void ThermalVisitor<dim>::ApplyTemperatureBoundaryConditionsToTransportedVariables( )
{
  const typename vector<Node<dim>*>::const_iterator modelNodesEnd( model_.Region("Model").NodesEnd() );
  VARIABLE_FLAG status;

  for( typename vector<Node<dim>*>::const_iterator it( model_.Region("Model").NodesBegin() ); it != modelNodesEnd; ++it )
  {
    status = (*it)->Status(t_key_);      // variable flag from temeprature
    (*it)->Status( hCl_key_, status );   // enthalpy content liqiuid
    (*it)->Status( mt_key_, status );    // fluid density

  }

}


/**

Loops over the nodes writing basis fluid (from H2O-EOS) and rock properties there.

@todo  if this visitor gets applied to another region than the model, all hell will break loose! - fix

*/
template<uint32_t dim>
void ThermalVisitor<dim>::SetInitialProperties(Model<dim>* model)
{
  const typename vector<Node<dim>*>::const_iterator modelNodesEnd( model_.Region("Model").NodesEnd() );
  for( typename vector<Node<dim>*>::const_iterator it( model_.Region("Model").NodesBegin() ); it != modelNodesEnd; ++it )
  {
    (*it)->Read(t_key_,         t   ); // Celsius, from CSMP
    (*it)->Read(p_key_,         p   ); // Pascal, from CSMP
    (*it)->Read(nphi_key_,      phi );
    (*it)->Read(rho_rock_key_,  rho_rock );
    (*it)->Read(cp_fluid_key_,  cp_fluid );
    (*it)->Read(cp_rock_key_,   cp_rock );
    (*it)->Read(beta_rock_key_, beta_rock);
    
    pore_volume     = phi();
    rock_volume     = 1.0 - phi();
    rl()            = water.Density(t(),p() ); 
    mt()            = rl();
    volume_factor() = 1.;
    mu()            = water.Viscosity(t(), p());

    hl()           = water.Enthalpy(t(),p() );
    hVl()          = hl()*rl();
    hCl()          = hl()*mt();
    nQ()           = 0.;
    dRho =         0.;

    beta_tot()     = mt()*phi()*(water.Compressibility(t(),p()) + beta_rock()); 
    cp_fluid()     = water.HeatCapacity (t(), p());
    nCPT()         = (rock_volume*cp_rock()*rho_rock() + pore_volume*cp_fluid()*rl());  

	  beta_fluid()   = water.Compressibility(t(),p()); // addition Sebastian and Alessandro, adding calculation of fluid compressibility
	
    (*it)->Store( beta_tot_key_, beta_tot );
    (*it)->Store( beta_fluid_key_,  beta_fluid );
    (*it)->Store( cp_fluid_key_, cp_fluid );
    (*it)->Store( nCPT_key_, nCPT );
    (*it)->Store( hl_key_,  hl );
    (*it)->Store( hVl_key_, hVl );
    (*it)->Store( hCl_key_, hCl );
    (*it)->Store( hClp_key_, hCl );
    (*it)->Store( mt_key_,  mt );
    (*it)->Store( rl_key_,  rl );
    (*it)->Store( mu_key_,  mu );
    (*it)->Store( nQ_key_,  nQ );
    (*it)->Store( dRho_key_,  dRho );
    (*it)->Store( volume_factor_key_, volume_factor );
  }
  this->ApplyTemperatureBoundaryConditionsToTransportedVariables();

} // end SetInitialProperties



template class ThermalVisitor<1>;
template class ThermalVisitor<2>;
template class ThermalVisitor<3>;

/* TODO: convert into template so that visitor can be applied to Faces
template<uint32_t dim, template<uint32_t> class CELL>
template class ThermalVisitor<1,Face>;
template class ThermalVisitor<2,Face>;
template class ThermalVisitor<3,Face>;
*/

} // end csmp

