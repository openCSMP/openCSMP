//
//  FlatBoundaryStressVisitor.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Updated by Hossein Agheshlui on 09 April 16.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef FLAT_BOUNDARY_STRESS_VISITOR_H
#define FLAT_BOUNDARY_STRESS_VISITOR_H

#include "Visitor.h"
#include "VectorVariable.h"
#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Boundary;
template<size_t> class Model;

/**
     /// This visitor is to apply forces on flat boundary faces. 
     /// It accesses a flat face and computes the nodal forces based on the normal and shear stresses defined in the config file.

*/


// see DOxygen doc further below
template<size_t dim>
class FlatBoundaryStressVisitor : public Visitor<dim> {
  public:
    /// by default, any entries into the RHS force vectors will get overwritten
    explicit FlatBoundaryStressVisitor( const Model<dim>&, bool overwrite_force_vector=true );

    virtual ~FlatBoundaryStressVisitor();

    /// to zero out force vector, prior to accumulation of forces.
    virtual void Visit( Boundary<dim>* );
  
    /// computes nodal forces weighted by node
    virtual void Visit( Face<dim>* );

  private:
    const csmp::Index Sn_key_, Ss_key_; ///< Principal stresses, the vertical and two horizontal ones, and the azimuth of the Sh_max_. Anderson Stress convention.
    const csmp::Index F_key_;           ///< (nodal) 'force' (N)
    bool overwrite_previous_forces_;
    VectorVariable<dim>   vc_, vs_;
    VectorVariable<dim>   nrml_;
};

/** 
     Workaround for applying stress boundary conditions on a Model
     in the form of nodal forces.
     
     Visits the model Boundaries looking for the variables 
     'normal stress' and 'shear stress' and computes 
     nodal forces from these.
     
     Use this, for instance to compute and apply overburden stresses.
     
     @attention by default all forces that were previously assigned are 
     zeroed first.
     
     @attention compressive (inward pointing) stresses are positive.
     
     @author SKM 29/9/2014
*/


} // end csmp

#endif /* defined(BOUNDARY_STRESS_VISITOR_H) */
