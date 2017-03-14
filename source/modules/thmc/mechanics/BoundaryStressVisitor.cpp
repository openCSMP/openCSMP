//
//  BoundaryStressVisitor.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Updated by Hossein Agheshlui on 09 April 16.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include "BoundaryStressVisitor.h"
#include "Model.h"
#include "Exception.h"
#include "Boundary.h"
#include "VectorVariable.h"
#include "variableOperations.h"
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
      Sv_key_(model.Database().StorageKey("vertical stress")),
	    Ss_Hmax_key_(model.Database().StorageKey("max horizontal stress")),
	    Ss_Hmin_key_(model.Database().StorageKey("min horizontal stress")),
      F_key_(model.Database().StorageKey("force")),
      overwrite_previous_forces_(overwrite_force_vector)
{
    if ( Sv_key_.place != FACE || Sv_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "BoundaryStressVisitor (constructor):",
                              "'normal stress' must be a SCALAR face variable." );

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
   double64  Sv = f->Parent(INSIDE)->Read( Sv_key_ );
   double64  SH = f->Parent(INSIDE)->Read(Ss_Hmax_key_);
   double64  Sh = f->Parent(INSIDE)->Read(Ss_Hmin_key_);
   double64  Ss_magnitude = hypot(SH,Sh);

   // if the values are negligibly small nothing needs to be done
   if ( fabs(Sv) < std::numeric_limits<double64>::epsilon()
        and fabs(Ss_magnitude) < std::numeric_limits<double64>::epsilon() )
     return;

   // creating the Andersonian stress tensor. SH always needs to be be greater than or equal to Sh.
   constexpr bool three_dimensional( dim == 3 );
   if ( three_dimensional ) {
        if ( Sv >= SH )
          stressTensor_ = makeTensor(ANY,ANY,ANY, Sv, 0., 0., 0., SH, 0., 0., 0., Sh);
        else if ( Sv < SH && Sv >= Sh )
          stressTensor_ = makeTensor(ANY,ANY,ANY, SH, 0., 0., 0., Sv, 0., 0., 0., Sh);
        else if ( Sv < Sh )
          stressTensor_ = makeTensor(ANY,ANY,ANY, SH, 0., 0., 0., Sh, 0., 0., 0., Sv);
        else throw csmp::Exception(ERROR, "BoundaryStressVisitor<dim>::Visit:",
                                  "the input values of andersonian stresses are not correct.");
     }

   // getting the unit normal and verifying that it is outward pointing
   // TODO: +*/normal directions have to be fixed so that this is guaranteed
   RetrieveNormalWithCorrectDirection( f, nrml_ );  //HA. TODO: Check to make sure that it is outward pointing.

   // Stress vector on an arbitrary plane
   VectorVariable<dim> T = stressTensor_ * nrml_;
   
   // project stress on unit normal
   // negative value because nrml is outward pointing
   nrml_ = nrml_ * -dotProduct(T,nrml_);
   VectorVariable<dim> Ts = T - nrml_;
   
   // going from specific- to areal stress on the face
   const size_t  face_nodes(f->Nodes());  //HA. This is added to be used throughout the method instead of reading it at few different places.
   nrml_ = nrml_ * (f->Area() / static_cast<double64> (face_nodes));
  
   // standard FE accumulation of stresses across a boundary
   // adding the normal stress as force contribution to the nodes of the face
   for ( size_t i=0U; i<face_nodes; i++ ) {
        // the status of the variable is not touched
        f->N(i)->Read( F_key_, vc_ );
        vc_ += nrml_;
        f->N(i)->Store( F_key_, vc_ );
     }
  
   // 2. ASSIGNMENT OF THE SHEAR STRESS COMPONENT
   // -------------------------------------------
  
   // computing and assigning the shear stresses
   // Shear stress vector is divided by the number of nodes and added to the exisiting force vector at the nodes. 
   Ts *= (f->Area() / static_cast<double64> (face_nodes));
  
   for (size_t i = 0U; i < face_nodes; i++) {
      f->N(i)->Read(F_key_, vc_);
      vc_ += Ts;
      f->N(i)->Store(F_key_, vc_);
   }

} // end Visit(face)








/**  when normals at a boundary are not correct; this is a workaround.
  
     Idea: normal must be pointing from face-parent1 to face-parent2.
     This is established by comparisons of barycentre differences.
     
     @attention assumption is made that when a face is at the model 
     boundary only the first parent is assigned.

*/
template<size_t dim>
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
   
   throw csmp::Exception( ERROR, "BoundaryStressVisitor<dim>::RetrieveNormalWithCorrectDirection:",
                         "case of boundary inside the model not handled yet." );

 } // end



//template class BoundaryStressVisitor<2U>;
template class BoundaryStressVisitor<3U>;

} // end csmp
