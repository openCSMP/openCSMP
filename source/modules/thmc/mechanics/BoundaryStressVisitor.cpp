//
//  BoundaryStressVisitor.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include "BoundaryStressVisitor.h"
#include "Model.h"
#include "Exception.h"
#include "Boundary.h"
#include "VectorVariable.h"
#include "Face.h"

namespace csmp {

/** 
    For the computation of the nodal 'gravity force'
    
    Application level is model so that all boundaries are visited
    followed by a visitation of their faces.
*/
template<size_t dim>
BoundaryStressVisitor<dim>::BoundaryStressVisitor( const Model<dim>& model,
                                                   bool overwrite_force_vector )
    : Visitor<dim>( BOUNDARY, FACE ),
      Sn_key_(model.Database().StorageKey("normal stress")),
      Ss_key_(model.Database().StorageKey("shear stress")),
      F_key_(model.Database().StorageKey("force")),
      overwrite_previous_forces_(overwrite_force_vector)
//      nrml_(dim)
{
    if ( Sn_key_.place != FACE || Sn_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'normal stress' must be a SCALAR face variable." );

    if ( Ss_key_.place != FACE || Ss_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'shear stress' must be a VECTOR face variable." );

    if ( F_key_.place != NODE || F_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'force' must be a VECTOR variable placed on the nodes." );
}




template<size_t dim>
BoundaryStressVisitor<dim>::~BoundaryStressVisitor()
{
}



/**
    To zero out initial values assigned to the boundary.
*/
template<size_t dim>
void BoundaryStressVisitor<dim>::Visit( Boundary<dim>* b )
{
   if ( overwrite_previous_forces_ ) {
         std::cout <<"\nBoundaryStressVisitor<dim>::Visit(boundary): resetting previous 'force' conditions.\n";
         VectorVariable<dim> zero;
         zero = 0.;
         b->InputPropertyValue( "force", zero, COMPLETE );
      }
}




/**
    This method interpolates the applied 
    element-centric forces to the nodes.
 
    An approach is used whre the contribution is
    weighted on the basis of the distance of the 
    barycenter to the node.
    
    @attention this is an integration
    (for uniform nodal forces, see option,
    we get non-uniform stresses)
 
*/
template<size_t dim>
void BoundaryStressVisitor<dim>::Visit( Face<dim>* f )
{
   // 1. NORMAL STRESS COMPONENT
   // --------------------------
//   double64  Sn = dynamic_cast<Element<dim>*>(f)->Read( Sn_key_ );
   double64  Sn = f->Read( Sn_key_ );
   // the shear stress in the plane of the face
   f->Read( Ss_key_, vc_ );
   const double64  Ss_magnitude = vc_.Length();

   // if the values are negligibly small nothing needs to be done
   if ( fabs(Sn) < std::numeric_limits<double64>::epsilon()
        and fabs(Ss_magnitude) < std::numeric_limits<double64>::epsilon() )
     return;

   // getting the unit normal and verifying that it is outward pointing
   // TODO: normal directions have to be fixed so that this is guaranteed
   // f->UnitNormal( nrml_ );
   RetrieveNormalWithCorrectDirection( f, nrml_ );

   // project stress on unit normal
   // negative value because nrml is outward pointing
   nrml_ *= -Sn;

   // going from specific- to areal stress on the face
   const bool interpolate_Sn_across_nodes(false);

   // smooth force-field on element boundary
   if ( interpolate_Sn_across_nodes ) {
        // adding the normal stress as force contribution to the nodes of the face
        for ( size_t i=0U; i<f->Nodes(); i++ ) {
             const double64 weighting(f->N(i)->Parents());
             f->N(i)->Read( F_key_, vc_ );
             vc_ += (nrml_ / weighting);
             f->N(i)->Store( F_key_, nrml_ );
          }
     }
   // standard FE accumulation of stresses across a boundary
   else {
        const double64  face_nodes(f->Nodes());
        Sn *= (f->Area() / face_nodes);
        // adding the normal stress as force contribution to the nodes of the face
        for ( size_t i=0U; i<f->Nodes(); i++ ) {
             f->N(i)->Read( F_key_, vc_ );
             vc_ += nrml_;
             f->N(i)->Store( F_key_, nrml_ );
          }
     }
  
/* NOT IMPLEMENTED YET

   // 2. ASSIGNMENT OF THE SHEAR STRESS COMPONENT
   // -------------------------------------------
   // renormalizing the unit normal
   //double64 sum(0.);
   //for ( size_t i=0U; i<dim; i++ ) sum += nrml_[i] * nrml_[i];
   //const double64 length = sqrt(sum);
   //for ( size_t i=0U; i<dim; i++ ) nrml_[i] /= length;
   nrml_.EuclideanNormalize();
  
   // computing and assigning the shear stresses
   dynamic_cast<Element<dim>*>(f)->Read( Ss_key_, vc_ ); // the shear stress in the plane of the face
   nrml_ = vc_.ProjectOnto( nrml_ );
  
   for ( size_t i=0U; i<f->Nodes(); i++ ) {
        f->N(i)->Read( F_key_, vc_ );
        vc_ += nrml_;
     }
   dynamic_cast<Element<dim>*>(f)->Store( F_key_, vc_ );
*/
} // end Visit(Element)






/*
    EARLIER VERSION THAT INTEGRATES ACROSS FACES
    
//
    This method attempts an integration of the forces at the nodes.
 
    Adds normal and shear stresses onto the nodes.
    
    @attention this means that previous values must make sense.
    
    The stresses are integrated over the boundary surface
    related to the face and the apportioned equally, 
    subdividing contributions by the number of its nodes.
    
    @attention this automatically takes care correctly of the 
    midside nodes.
 
//
template<size_t dim>
void BoundaryStressVisitor<dim>::Visit( Face<dim>* f )
{
   // NORMAL STRESS COMPONENT
   // -----------------------
//   double64  Sn = dynamic_cast<Element<dim>*>(f)->Read( Sn_key_ );
   double64  Sn = f->Read( Sn_key_ );
   // the shear stress in the plane of the face
   f->Read( Ss_key_, vc_ );
   const double64  Ss_magnitude = vc_.Length();

   // if the values are negligibly small nothing needs to be done
   if ( fabs(Sn) < std::numeric_limits<double64>::epsilon()
        and fabs(Ss_magnitude) < std::numeric_limits<double64>::epsilon() )
     return;
  
   // going from specific- to areal stress on the face
   const double64  face_nodes(f->Nodes());
   Sn *= (f->Area() / face_nodes);
  
   // getting the unit normal and verifying that it is outward pointing
   // TODO: normal directions have to be fixed so that this is guaranteed
   // f->UnitNormal( nrml_ );
   RetrieveNormalWithCorrectDirection( f, nrml_ );
  
   // project stress on unit normal; (-) because nrml is outward pointing
   nrml_ *= -Sn;

// vector version:   for ( size_t i=0U; i<dim; i++ ) nrml_[i] *= (Sn / face_nodes);
  
   // adding the computed normal stress to the forces already acting on the nodes
   for ( size_t i=0U; i<f->Nodes(); i++ ) {
         f->N(i)->Read( F_key_, vc_ );
         vc_ += nrml_;
         f->N(i)->Store( F_key_, nrml_ );
     }

   // ASSIGNMENT OF THE SHEAR STRESS COMPONENT
   // ----------------------------------------
   // renormalizing the unit normal
   //double64 sum(0.);
   //for ( size_t i=0U; i<dim; i++ ) sum += nrml_[i] * nrml_[i];
   //const double64 length = sqrt(sum);
   //for ( size_t i=0U; i<dim; i++ ) nrml_[i] /= length;
   nrml_.EuclideanNormalize();
  
   // computing and assigning the shear stresses
   dynamic_cast<Element<dim>*>(f)->Read( Ss_key_, vc_ ); // the shear stress in the plane of the face
   nrml_ = vc_.ProjectOnto( nrml_ );
  
   for ( size_t i=0U; i<f->Nodes(); i++ ) {
        f->N(i)->Read( F_key_, vc_ );
        vc_ += nrml_;
     }
   dynamic_cast<Element<dim>*>(f)->Store( F_key_, vc_ );

} // end Visit(Element)

*/



/**  when normals at a boundary are not correct; this is a workaround.
  
     Idea: normal must be pointing from face-parent1 to face-parent2.
     This is established by comparisons of barycentre differences.
     
     @attention assumption is made that when a face is at the model 
     boundary only the first parent is assigned.

*/
template<size_t dim>
void BoundaryStressVisitor<dim>::RetrieveNormalWithCorrectDirection( Face<dim>* f,  VectorVariable<dim>& nrml )
 {
    f->UnitNormal( nrml );
    const Point<dim> face_ctr(f->BaryCenter());
   
    // establishing whether the face is at a model boundary
    if ( f->Parent(INSIDE) == nullptr or f->Parent(OUTSIDE) == nullptr ) {
         if ( f->Parent(INSIDE) != nullptr ) {
              const Point<dim> parent_ctr = f->Parent(INSIDE)->BaryCenter();
              // normalizing a copy of the unit normal with the distance between the 2 barycenters
              const double64 face_parent_ctr_distance(parent_ctr.DistanceTo(face_ctr));
              Point<dim> pnrml;
              for ( size_t i=0U; i<dim; i++ ) pnrml[i] = nrml[i];
              pnrml /= face_parent_ctr_distance;
              const double64 face_parent_nrml_tip_distance(parent_ctr.DistanceTo(face_ctr + pnrml));
              // now normal is added to face center, if this gets a point closer to parent center
              // the normal is inward pointing and needs correction
              if ( face_parent_nrml_tip_distance < face_parent_ctr_distance ) nrml *= -1.;
              return;
           }
      }
   
   throw csmp::Exception( ERROR, "BoundaryStressVisitor<dim>::RetrieveNormalWithCorrectDirection:",
                         "case of boundary inside the model not handled yet." );

 } // end





//template class BoundaryStressVisitor<2U>;
//template class BoundaryStressVisitor<3U>;
template class BoundaryStressVisitor<3U>;

} // end csmp
