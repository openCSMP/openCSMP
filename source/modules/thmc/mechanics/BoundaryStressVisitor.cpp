// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  BoundaryStressVisitor.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Updated by Hossein Agheshlui on 09 April 16.
//

#include "BoundaryStressVisitor.h"
#include "Model.h"
#include "NodeManifold.h"
#include "Exception.h"
#include "Boundary.h"
#include "VectorVariable.h"
#include "variableOperations.h"
#include "Face.h"

using namespace std;

namespace csmp {

/** 
    For the computation of the nodal 'gravity force'
    
    Application level is model so that all boundaries are visited
    followed by a visitation of their faces.
    
    @param model computation domain.
    
    @param SV_variable denotes overburden = vertical stress 
    calculated by integration of the sediment column dry density 
    up to the earth' surface.
*/
template<uint32_t dim>
BoundaryStressVisitor<dim>::BoundaryStressVisitor( const Model<dim>& model,
                                                   const char* SV_variable,
                                                   bool overwrite_force_vector )
    : Visitor<dim>( BOUNDARY, FACE ),
      Sv_key_(model.Database().StorageKey(SV_variable)),
	    SH_max_key_(model.Database().StorageKey("max horizontal stress")),
	    Sh_min_key_(model.Database().StorageKey("min horizontal stress")),
      F_key_(model.Database().StorageKey("force")),
      stress_computation_(BOUNDARY_STRESS::SINGLE_VALUED),
      overwrite_previous_forces_(overwrite_force_vector)
{
    if ( Sv_key_.place != NODE || Sv_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                               SV_variable, "must be a SCALAR variable placed on the nodes." );

    if ( SH_max_key_.place != FACE || SH_max_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "the 'max horizontal stress' SHmax must be a SCALAR variable placed on the FACE." );

    if ( Sh_min_key_.place != FACE || Sh_min_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "the 'min horizontal stress' Shmin must be a SCALAR variable placed on the FACE." );

    if ( F_key_.place != NODE || F_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'force' must be a VECTOR variable placed on the nodes." );
}



template<uint32_t dim>
BoundaryStressVisitor<dim>::BoundaryStressVisitor( const Model<dim>& model,
                                                   BOUNDARY_STRESS stress_computation,
                                                   bool overwrite_force_vector )
    : Visitor<dim>( BOUNDARY, FACE ),
      bstress_key_(model.Database().StorageKey("boundary stress")),
      F_key_(model.Database().StorageKey("force")),
      stress_computation_(stress_computation),
      overwrite_previous_forces_(overwrite_force_vector)
{
    assert( stress_computation_ == BOUNDARY_STRESS::SINGLE_VALUED );
  
    if ( bstress_key_.place != BOUNDARY || bstress_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'boundary stress' must be a VECTOR variable placed on the boundary." );

    if ( F_key_.place != NODE || F_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'force' must be a VECTOR variable placed on the nodes." );
}






/**
    Collects the stress value placed on the boudary and 
    (optionally) zeros out the initial force values assigned to 
    this boundary.
*/
template<uint32_t dim>
void BoundaryStressVisitor<dim>::Visit( Boundary<dim>* b )
{
   if ( overwrite_previous_forces_ ) {
         std::cout <<"\nBoundaryStressVisitor<dim>::Visit(boundary): resetting previous 'force' conditions.\n";
         VectorVariable<dim> zero;
         zero = 0.;
         b->InputPropertyValue( "force", zero, COMPLETE );
      }
  
   if ( stress_computation_ == BOUNDARY_STRESS::SINGLE_VALUED )
     b->Read( bstress_key_, bstress_ );
  
   if ( stress_computation_ == BOUNDARY_STRESS::DEPTH_DEPENDENT ) {
       b->Read(SH_max_key_, SHmax_ );
       b->Read(Sh_min_key_, Shmin_ );
    }
}




/**
    Interpolates the applied
    element-centric forces to the nodes.
 
    An approach is used where the nodal contribution is
    weighted on the basis of the distance of the 
    barycenter to the node.
    
    @attention this is an integration
    (for uniform nodal forces, see option,
    we get non-uniform stresses)
 
*/
template<uint32_t dim>
void BoundaryStressVisitor<dim>::Visit( Face<dim>* f )
{
   bool Neumann_condition(false);

   // in the case of body forces
   if ( stress_computation_ == BOUNDARY_STRESS::DEPTH_DEPENDENT ) {
         for ( auto i{0U}; i<dim; ++i )
           if ( SHmax_.Flag() == NEUMANN or Shmin_.Flag() == NEUMANN ) {
                 Neumann_condition = true;
                 break;
             }
         if ( !Neumann_condition ) return;
     
         ApplyDepthDependentBoundaryConditions( f );
         return;
     }
  
   // when body forces are not turned on
   if ( stress_computation_ == BOUNDARY_STRESS::SINGLE_VALUED ) {
       for ( auto i{0U}; i<dim; ++i )
         if ( bstress_.Flag(i) == NEUMANN ) {
               Neumann_condition = true;
               break;
           }
       if ( !Neumann_condition ) return;

       ApplyConstantStressBoundaryConditions( f );
       return;
    }
  
   throw csmp::Exception( ERROR, "BoundaryStressVisitor<dim>::Visit:",
                            "boundary stress condition not implemented yet." );
} // end Visit(face)




/** 
    no vertical variation
*/
template<uint32_t dim>
void BoundaryStressVisitor<dim>::ApplyConstantStressBoundaryConditions( Face<dim>* f )
 {
   bool Neumann_condition(false);
   for ( auto i{0U}; i<dim; ++i )
     if ( bstress_.Flag(i) == NEUMANN ) {
           Neumann_condition = true;
           break;
       }
   if ( !Neumann_condition ) return;
   
   // going from specific- to area-integrated stress values to nodal forces
   const size_t  face_nodes(f->Nodes());
   force_ = bstress_ * (f->Area() / static_cast<double>(face_nodes));
  
   // adding the normal stress forces to the nodal forces of the nodes of the face
   for ( auto i{0U}; i<face_nodes; i++ ) {
        // the status of the variable is not touched
        f->N(i)->Read( F_key_, vc_ );
        vc_ += force_;
        f->N(i)->Store( F_key_, vc_ );
     }
}
 
 
  
/**
    variation with depth due to the action of body force = gravity
*/
template<uint32_t dim>
void BoundaryStressVisitor<dim>::ApplyDepthDependentBoundaryConditions( Face<dim>* f )
 {
   ScalarVariable Sv;
   f->Parent(INSIDE)->PropertyValueAtBaryCenter(Sv_key_, Sv );
   double SH = SHmax_() * Sv();
   double Sh = Shmin_() * Sv();
   assert( Sv() > 0. );
   assert( SH > 0. );
   assert( Sh > 0. );
   //double  Ss_magnitude = hypot(SH,Sh);
   //assert( Ss_magnitude );

   // if the values are negligibly small nothing needs to be done
   if ( fabs(Sv()) < numeric_limits<double>::epsilon() ) return;

   // creating the Andersonian stress tensor. SH always needs to be be greater than or equal to Sh.
   // TODO: check for rotations when one of the principal stresses is not aligned with a
   //       coordinate axis
   constexpr bool three_dimensional( dim == 3 );
   if ( three_dimensional ) {
        if ( Sv() >= SH )
          stressTensor_ = makeTensor(ANY,ANY,ANY, Sv(), 0., 0., 0., SH, 0., 0., 0., Sh);
        else if ( Sv() < SH && Sv() >= Sh )
          stressTensor_ = makeTensor(ANY,ANY,ANY, SH, 0., 0., 0., Sv(), 0., 0., 0., Sh);
        else if ( Sv() < Sh )
          stressTensor_ = makeTensor(ANY,ANY,ANY, SH, 0., 0., 0., Sh, 0., 0., 0., Sv());
        else throw csmp::Exception( ERROR, "BoundaryStressVisitor<dim>::Visit:",
                                   "the input values of Andersonian stresses are not correct.");
     }

   // getting the unit normal
   f->UnitNormal( nrml_ );

   // stress vector on the plane of the face (any orientation is OK, but normal must be outward pointing)
   VectorVariable<dim> faceStressVector = stressTensor_ * nrml_;
   
   // 2. ASSIGNMENT OF THE SHEAR STRESS COMPONENT
   // -------------------------------------------
   // going from specific- to area-integrated stress values to nodal forces
   const size_t  face_nodes(f->Nodes());
   faceStressVector = faceStressVector * (f->Area() / static_cast<double>(face_nodes));
  
   // adding the normal stress forces to the nodal forces of the nodes of the face
   for ( auto i{0U}; i<face_nodes; i++ ) {
        // the status of the variable is not touched
        f->N(i)->Read( F_key_, vc_ );
        vc_ += faceStressVector;
        f->N(i)->Store( F_key_, vc_ );
     }

 } // end (depth dependent computation)











/**  when normals at a boundary are not correct; this is a workaround.
  
     Idea: normal must be pointing from face-parent1 to face-parent2.
     This is established by comparisons of barycentre differences.
     
     @attention assumption is made that when a face is at the model 
     boundary only the first parent is assigned.

*/
template<uint32_t dim>
void BoundaryStressVisitor<dim>::RetrieveNormalWithCorrectDirection( Face<dim>* f, VectorVariable<dim>& nrml )
 {
    f->UnitNormal( nrml );
    const Point<dim> face_ctr(f->BaryCenter());
   
    // if the face is located at an external model boundary
    if ( f->Parent(OUTSIDE) == nullptr )
      {
          assert( f->Parent(INSIDE) != nullptr );
          const Point<dim> parent_ctr = f->Parent(INSIDE)->BaryCenter();
          // normalizing a copy of the unit normal with the distance between the 2 barycenters
          const double face_parent_ctr_distance(parent_ctr.DistanceTo(face_ctr));
          Point<dim> pnrml;
          for ( auto i{0U}; i<dim; i++ ) pnrml[i] = nrml[i];
          pnrml /= face_parent_ctr_distance;
          const double face_parent_nrml_tip_distance(parent_ctr.DistanceTo(face_ctr + pnrml));
          // now normal is added to face center, if this gets a point closer to parent center
          // the normal is inward pointing and needs correction
          if ( face_parent_nrml_tip_distance < face_parent_ctr_distance ) nrml *= -1.;
          return;
      }
   
   throw csmp::Exception( ERROR, "BoundaryStressVisitor<dim>::RetrieveNormalWithCorrectDirection:",
                         "case of boundary inside the model not handled yet." );

 } // end



//template class BoundaryStressVisitor<2U>;
template class BoundaryStressVisitor<3U>;

} // end csmp
