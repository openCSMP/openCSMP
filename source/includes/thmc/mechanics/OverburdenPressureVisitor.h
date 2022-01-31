//
//  OverburdenPressureVisitor.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef OVERBURDEN_PRESSURE_VISITOR_H
#define OVERBURDEN_PRESSURE_VISITOR_H

#include "Visitor.h"

namespace csmp {

template<size_t> class VectorVariable;
template<size_t> class Element;
template<size_t> class Model;

/** 
     For a given Sv (overburden stress at the model top),
     this visitor computes:
     
     - the (hydraulic) 'conductivity' from the element 'permeability' and (nodal) 'fluid viscosity'
     
     - (nodal) 'gravity force' scalar variable as the product of g (specified by user) and the
       'dry rock density' specified as element property.
       
     - the (element) 'gravity term' vector variable used in for the integration
       of the fluid density via the corresponding PDE operator.
     
     @todo SKM: this visitor should include the influence of pore pressure.
 
*/
template<size_t dim>
class OverburdenPressureVisitor : public Visitor<dim> {
  public:
    /// for the computation of nodal 'gravity force'
    OverburdenPressureVisitor( Model<dim>&, double acc_gravity ); ///< local constant for reservoir

    virtual ~OverburdenPressureVisitor();
  
    /// this method computes all the element properties 'K' and 'gravity term', etc.
    virtual void Visit( Element<dim>* );

    /// nothing needs to be done at the level of the model
    virtual void Visit( Model<dim>* ) {}

  private:
    double dryDensityFromBulkDensity( double rho_bulk, double rho_fluid, double porosity );

  private:
    const csmp::Index    rhof_key_, rhor_key_,     ///< (nodal) fluid  and dry rock densities (scalars)
                         rhob_key_, phi_key_,      ///< bulk- (fluid+rock) density and porosity
                         k_key_, mu_key_, K_key_,  ///< (element) permeability, (nodal) 'fluid viscosity' and hydraulic 'conductivity' (scalar)
                         gf_key_, gt_key_;         ///< (nodal) gravity force and (elemental) 'gravity term' (vector)
    const double       acc_gravity_;
    VectorVariable<dim>  gravity_;
};


/// calculation ignoring the weight of air.
template<size_t dim>
inline double OverburdenPressureVisitor<dim>::dryDensityFromBulkDensity( double rho_bulk, double rho_fluid, double porosity )
 {
    return rho_bulk - porosity * rho_fluid;
 }


} // end csmp

#endif /* defined(OVERBURDEN_PRESSURE_VISITOR_H) */


