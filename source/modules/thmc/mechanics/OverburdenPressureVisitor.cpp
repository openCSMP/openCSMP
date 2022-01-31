//
//  OverburdenPressureVisitor.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include "OverburdenPressureVisitor.h"
#include "Exception.h"
#include "Model.h"

namespace csmp {

/**
     For the computation of the nodal 'gravity force'
     
     Hardwired variables:
     'fluid density'
     'dry rock density'
     'porosity'
     'permeability'
     'conductivity'
     'gravity term'  for hydrostatic fluid pressure gradient
     'gravity force' for body forces acting on the rock skeleton
*/
template<size_t dim>
OverburdenPressureVisitor<dim>::OverburdenPressureVisitor( Model<dim>& model, // Pa.s
                                                           double acc_gravity )
    : Visitor<dim>( MODEL, ELEMENT ),
      rhof_key_(model.Database().StorageKey("fluid density")),
      rhor_key_(model.Database().StorageKey("dry rock density")),
      rhob_key_(model.Database().StorageKey("bulk density")),
      phi_key_(model.Database().StorageKey("porosity")),
      k_key_(model.Database().StorageKey("permeability")),
      mu_key_(model.Database().StorageKey("fluid viscosity")),
      K_key_(model.Database().StorageKey("conductivity")),
      gf_key_(model.Database().StorageKey("gravity force")),
      gt_key_(model.Database().StorageKey("gravity term")),
      gravity_(PLAIN,PLAIN,PLAIN,0.,-1.,0.),
      acc_gravity_(std::fabs(acc_gravity)) // avoiding any sign confusions
{
    // TODO: add checks for types and placements of other variables

    if ( gf_key_.place != ELEMENT || gf_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "OverburdenPressureVisitor (constructor):",
                              "'gravity force' must be a VECTOR element variable." );

    if ( (gt_key_.place != ELEMENT_INTEGRATION_POINT and gt_key_.place != ELEMENT) || gt_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "OverburdenPressureVisitor (constructor):",
                               "'gravity term' must be a VECTOR variable placed on the element." );
}


template<size_t dim>
OverburdenPressureVisitor<dim>::~OverburdenPressureVisitor()
{
}




/**
    Loops over the element integration points and evaluates shear and tensile failure potential.
    Where failure occurred the variable 'failure' is set to 1.
    A distinction is made between tensile 'failure01' and shear failure 'failure' variables.
    
    @attention SKM 28/9/2014 - added case where stress is placed on the element
*/
template<size_t dim>
void OverburdenPressureVisitor<dim>::Visit( Element<dim>* e )
{
   // element properties
   // ------------------
   ScalarVariable  mu;
   e->PropertyValueAtBaryCenter( mu_key_, mu );
   // hydraulic conductivity computation
   e->Store( K_key_, makeScalar( PLAIN, e->Read(k_key_) / mu() ) );
   ScalarVariable  erho_fluid;
   e->PropertyValueAtBaryCenter( rhof_key_, erho_fluid );
  
   // gravity force (negative as it acts against Y-axis)
   const double dry_rhor = dryDensityFromBulkDensity( e->Read(rhob_key_), erho_fluid(), e->Read(phi_key_) );
   e->Store( rhor_key_, makeScalar( PLAIN, dry_rhor ) );
   gravity_(1) = dry_rhor * -acc_gravity_;
   e->Store( gf_key_, gravity_ );
  
   // gravity term
   gravity_(1) = -acc_gravity_ * erho_fluid() * e->Read(K_key_);
   e->Store( gt_key_, gravity_ );
    
} // end Visit(Element)

//template class OverburdenPressureVisitor<2U>;
template class OverburdenPressureVisitor<3U>;

} // end csmp
