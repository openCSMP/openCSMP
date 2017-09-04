//
//  BoundaryStressVisitor.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Updated by Hossein Agheshlui on 09 April 16.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef BOUNDARY_STRESS_VISITOR_H
#define BOUNDARY_STRESS_VISITOR_H

#include "Visitor.h"
#include "TensorVariable.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Boundary;
template<size_t> class Model;
template<size_t> class VectorVariable;

/**
     @brief BoundaryStressVisitor accesses a Boundary object that must be planar (flat),
     computing nodal forces based on the normal and shear stresses applied to that Boundary
     element. If the boundary is inclined or vertical and body forces are enabled,
     the visitor takes the vertical variation in the boundary normal stress into account.

     @note This visitor allows for a workaround for applying stress boundary conditions
     if one does not want to use the Neumann condition functionality of the PDE_Integrator.
 
     Input are the Face variables 'normal stress' and 'shear stress.' The nodal
     forces are computed from these.
     
     Use this, for instance to compute and apply overburden stresses.
     
     @attention by default all forces that were previously assigned are 
     zeroed first.
     
     @attention compressive (inward pointing) stresses are positive.
     
     @attention this visitor assumes that the far-field stress regime is 
     Andersonian, i.e. one of the principal stresses is oriented perpendicular
     to the earth surface.
     
     @author SKM 29/9/2014
*/
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
	  const csmp::Index Sv_key_, ///< key to variable vertical stress (a scalar)
                Ss_Hmax_key_,  ///< maximum horizontal stress
                Ss_Hmin_key_;  ///< minimum horizontal stress
	  const csmp::Index F_key_;  ///< (nodal) 'force' (N)
	  VectorVariable<dim>  Svv_, SHv_, Shv_, vc_;
	  VectorVariable<dim>  nrml_;
	  TensorVariable<dim>  stressTensor_;
	  bool overwrite_previous_forces_;
};

} // end csmp

#endif /* defined(BOUNDARY_STRESS_VISITOR_H) */
