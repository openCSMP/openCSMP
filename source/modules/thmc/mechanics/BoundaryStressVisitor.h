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

template<uint32_t> class Element;
template<uint32_t> class Boundary;
template<uint32_t> class Model;
template<uint32_t> class VectorVariable;

enum class BOUNDARY_STRESS { SINGLE_VALUED, DEPTH_DEPENDENT /* TODO: TRIANGULAR? */ };

/**
     @brief BoundaryStressVisitor accesses a Boundary object that must be a surface,
     computing nodal forces based on the normal and shear stresses applied to that Boundary
     element. If the boundary is inclined or vertical and body forces are enabled,
     the visitor takes the vertical variation in the boundary normal stress into account.

     @attention this visitor assumes that the driving stresses are Andersonian, i.e.
     that one of the principal stress axes is aligned with the vertical (Y)-axis of the model.

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
template<uint32_t dim>
class BoundaryStressVisitor : public Visitor<dim> {
  public:
    /// by default, any entries into the RHS force vectors will get overwritten
    BoundaryStressVisitor( const Model<dim>&,
                           const char* SV_variable, ///< should vary with depth
                           bool overwrite_force_vector=true );

    /// computations using absolute stress values on the boundary
    BoundaryStressVisitor( const Model<dim>&,
                           BOUNDARY_STRESS=BOUNDARY_STRESS::SINGLE_VALUED,
                           bool overwrite_force_vector=true );

    virtual ~BoundaryStressVisitor();

    /// to zero out force vector, prior to accumulation of forces.
    virtual void Visit( Boundary<dim>* );
  
    /// computes nodal forces weighted by node
    virtual void Visit( Face<dim>* );
  
  private:
    /// normals at the boundaries are not correct; this is a workaround
    void RetrieveNormalWithCorrectDirection( Face<dim>*,  VectorVariable<dim>& nrml );
  
    /// no vertical variation
    void ApplyConstantStressBoundaryConditions( Face<dim>* );
  
    /// variation with depth due to the action of body force = gravity
    void ApplyDepthDependentBoundaryConditions( Face<dim>* );
  

  private:
	  const csmp::Index    Sv_key_,               ///< vertical (overburden) stress (a scalar)
                         SH_max_key_,           ///< maximum horizontal stress expressed in terms of SV
                         Sh_min_key_,           ///< minimum horizontal stress Sh_min/SV
                         bstress_key_;          ///< stress acting on the boundary face (one value per boundary)
	  const csmp::Index    F_key_;                ///< (nodal) 'force' (N)
	  VectorVariable<dim>  vc_, force_;           ///< forces read from the nodes
	  VectorVariable<dim>  nrml_, bstress_;
	  TensorVariable<dim>  stressTensor_;
    ScalarVariable       SHmax_, Shmin_;        ///< stress ratios read from the boundary
    BOUNDARY_STRESS      stress_computation_;
	  bool                 overwrite_previous_forces_;
};

} // end csmp

#endif /* defined(BOUNDARY_STRESS_VISITOR_H) */
