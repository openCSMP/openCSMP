//
//  AndersonianBoundaryStressVisitor.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include "AndersonianBoundaryStressVisitor.h"
#include "Model.h"
#include "Exception.h"
#include "Boundary.h"
#include "VectorVariable.h"
#include "variableOperations.h"
#include "Face.h"
#include "StressRegime.h"
#include "mechanics.h"

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
    
    TODO: the assumption is that the vertical stress is the largest principal stress, relax this
*/
AndersonianBoundaryStressVisitor::AndersonianBoundaryStressVisitor( const Model<3U>& model,
                                                                    const StressRegime& stress_regime,
                                                                    const char* SV_variable, ///< should vary with depth
                                                                    bool overwrite_force_vector )
    : Visitor<3U>( BOUNDARY, FACE ),
      stress_regime_(stress_regime),
      Andersonian_stress_( Point<3U>(0.,-1.,0.), // orientation vector for sigma1
                           Point<3U>(stress_regime.MinimumHorizontalStressUnitVector()[0],
                                     stress_regime.MinimumHorizontalStressUnitVector()[1],
                                     stress_regime.MinimumHorizontalStressUnitVector()[2]),
                           stress_regime.VerticalStressMagnitude(),
                           stress_regime.MaximumHorizontalStressMagnitude(),
                           stress_regime.MinimumHorizontalStressMagnitude() ),
      Sv_key_(model.Database().StorageKey(SV_variable)),
	    Sn_key_(model.Database().StorageKey("normal stress")),
	    Ss_key_(model.Database().StorageKey("shear stress")),
      F_key_(model.Database().StorageKey("force")),
      overwrite_previous_forces_(overwrite_force_vector)
{
    if ( Sv_key_.place != NODE || Sv_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "AndersonianBoundaryStressVisitor (constructor):",
                               SV_variable, "must be a SCALAR variable placed on the node." );

    if ( Sn_key_.place != FACE || Sn_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "AndersonianBoundaryStressVisitor (constructor):",
                              "the 'normal stress' must be a SCALAR variable placed on the FACE." );

    if ( Ss_key_.place != FACE || Ss_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "AndersonianBoundaryStressVisitor (constructor):",
                              "the 'shear stress' must be a SCALAR variable placed on the FACE." );

    if ( F_key_.place != NODE || F_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "AndersonianBoundaryStressVisitor (constructor):",
                              "'force' must be a VECTOR variable placed on the nodes." );
  
    // rotating the stress tensor around the y axis
    Andersonian_stress_.Rotate( 'y', stress_regime_.MaximumHorizontalStressTrend() );
    // initialising the far-field stress tensor
    Andersonian_stress_.CartesianStressTensor( far_field_stress_ );
}





AndersonianBoundaryStressVisitor::~AndersonianBoundaryStressVisitor()
{
}




/**
    Collects the stress value placed on the boudary and 
    (optionally) zeros out the initial force values assigned to 
    this boundary.
*/
void AndersonianBoundaryStressVisitor::Visit( Boundary<3U>* b )
{
   if ( overwrite_previous_forces_ ) {
         std::cout <<"\nAndersonianBoundaryStressVisitor<dim>::Visit(boundary): resetting previous 'force' conditions.\n";
         VectorVariable<3U> zero;
         zero = 0.;
         b->InputPropertyValue( "force", zero, COMPLETE );
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
void AndersonianBoundaryStressVisitor::Visit( Face<3U>* f )
{
   const double nodes(static_cast<double>(f->Nodes()));
  
   // 1. projecting the far-field stress onto the Face normal
   // -------------------------------------------------------
   f->UnitNormal( nrml_ );
   // stress vector on the plane of the face (any orientation is OK, but normal must be outward pointing)
   VectorVariable<3U> faceStressVector = far_field_stress_ * nrml_;
   // going from specific- to area-integrated stress values to nodal forces
   const size_t  face_nodes(f->Nodes());
   faceStressVector = faceStressVector * (f->Area() / static_cast<double>(face_nodes));
  
   // TODO: make sure that the stress is inward pointing relative to the boundary
   if ( dotProduct(faceStressVector,nrml_) >= 0. )
     faceStressVector *= -1.;
  
   // adding the stress forces augmented by the overburden stress to the nodal forces of the nodes of the face
   for ( auto i{0}; i<face_nodes; i++ ) {
        // overburden stress
        f->N(i)->Read( Sv_key_, Sv_ );
        Sv_ /= nodes;
        // the vertical stress acts opposite of the Y axis
        faceStressVector(1) += -Sv_();
        // the status of the variable is not touched
        f->N(i)->Read( F_key_, vc_ );
        vc_ += faceStressVector;
        f->N(i)->Store( F_key_, vc_ );
     }

   // 2. computing the face normal and shear stresses for display
   // -----------------------------------------------------------
   Point<3U> nrml(nrml_[0],nrml_[1],nrml_[2]);
   normalAndShearStressOnPlane( far_field_stress_, nrml, Sn_(), Ss_() );
  
   f->Store( Sn_key_, makeScalar( f->Status(Sn_key_), Sn_() ) );
   f->Store( Ss_key_, makeScalar( f->Status(Ss_key_), Ss_() ) );

} // end Visit(face)


} // end csmp
