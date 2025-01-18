//
//  AndersonianBoundaryStressVisitor.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Updated by Hossein Agheshlui on 09 April 16.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef ANDERSONIAN_BOUNDARY_STRESS_VISITOR_H
#define ANDERSONIAN_BOUNDARY_STRESS_VISITOR_H

#include "Visitor.h"
#include "TensorVariable.h"
#include "StressRotate.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Boundary;
template<uint32_t> class Model;
template<uint32_t> class VectorVariable;
class StressRegime;

/**
     @brief  Apply stresses to the boundaries of a three-dimensional geomechanic model.
 
     @author SKM 29/9/2014

     The stresses are entered via a StressRegime object, but the variable 'overburden pressure'
     storing the vertical stress Sv must be defined and initialised.
     
     The AndersonianBoundaryStressVisitor accesses all Boundary objects,
     and resolves the boundary normal stress and shear stresses using the
     normals of the boundary Face objects which must be outward pointing.
     
     In the corresponding calculations body forces must be enabled. 
 
     @note If a boundary is inclined or vertical, this is fine. In this case,
     the visitor takes the vertical variation in the boundary normal stress into account.
 
     Input variables are the stress state supplied by the StressRegime object.
     
     The visitor also computes the Face variables 'normal stress' and 'shear stress'
     on the visited boundaries. The nodal forces that it will apply
     are computed from these.
 
     @attention by default all forces that were previously assigned are 
     zeroed first.
     
     @attention Confining stresses are always compressive and result inward pointing forces.
     
     @attention this visitor assumes that the far-field stress regime is 
     Andersonian, i.e. one of the principal stresses is oriented perpendicular
     to the earth surface, i.e. that one of the principal stress axes is aligned
     with the vertical (Y)-axis of the CSMP model.
     
     @note This visitor allows for a workaround for applying stress boundary conditions
     if one does not want to use the Neumann condition functionality of the PDE_Integrator.
*/
class AndersonianBoundaryStressVisitor : public Visitor<3U> {
  public:
    /// by default, any entries into the RHS force vectors will get overwritten
    AndersonianBoundaryStressVisitor( const Model<3U>&,
                                      const StressRegime&,
                                      const char* SV_variable="overburden pressure", ///< should vary with depth
                                      bool overwrite_force_vector=false );

    virtual ~AndersonianBoundaryStressVisitor();

    /// to zero out force vector if overwrite_force_vector=true, prior to accumulation of forces
    virtual void Visit( Boundary<3U>* );
  
    /// computes nodal forces ensuing from stress and weighted by number of nodes of the boundary face
    virtual void Visit( Face<3U>* );
  
  private:
    const StressRegime&  stress_regime_;        ///< input Andersonian stress state
    StressRotate         Andersonian_stress_;   ///< the full tensor constructed for the stress regimes
  
	  const csmp::Index    Sv_key_,               ///< vertical (overburden) stress (a scalar)
                         Sn_key_,               ///< normal stress acting on boundary faces
                         Ss_key_;               ///< shear stress acting on boundary faces
	  const csmp::Index    F_key_;                ///< (nodal) 'force' (N)
	  VectorVariable<3U>   vc_, force_;           ///< forces read from the nodes
	  VectorVariable<3U>   nrml_;                 ///< outward-pointing unit normal to face
	  TensorVariable<3U>   far_field_stress_;     ///< as computed from stress regime
    ScalarVariable       Sv_, Sn_, Ss_;         ///< scalar vertical, normal and shear stresses
	  bool                 overwrite_previous_forces_;
};

} // end csmp

#endif /* ANDERSONIAN_BOUNDARY_STRESS_VISITOR_H */
