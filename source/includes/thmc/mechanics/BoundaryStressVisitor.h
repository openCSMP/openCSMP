//
//  BoundaryStressVisitor.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef BOUNDARY_STRESS_VISITOR_H
#define BOUNDARY_STRESS_VISITOR_H

#include "Visitor.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Boundary;
template<size_t> class Model;
template<size_t> class VectorVariable;

// see DOxygen doc further below
template<size_t dim>
class BoundaryStressVisitor : public Visitor<dim> {
  public:
    /// by default, any entries into the RHS force vectors will get overwritten
    explicit BoundaryStressVisitor( const Model<dim>&, bool overwrite_force_vector=true );

    virtual ~BoundaryStressVisitor();

    /// to zero out force vector, prior to accumulation of forces.
    virtual void Visit( Boundary<dim>* );
  
    /// computes nodal forces weighted by node
    virtual void Visit( Face<dim>* );
  
  private:
    /// normals at the boundaries are not correct; this is a workaround
    void RetrieveNormalWithCorrectDirection( Face<dim>*,  VectorVariable<dim>& nrml );

  private:
    const csmp::Index Sn_key_, Ss_key_; ///< Face normal (scalar) - and shear (vector) stresses acting on the face.
    const csmp::Index F_key_;           ///< (nodal) 'force' (Pa)
    bool overwrite_previous_forces_;
    VectorVariable<dim>   vc_;
    VectorVariable<dim>   nrml_;
//    std::vector<double64> nrml_;
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
