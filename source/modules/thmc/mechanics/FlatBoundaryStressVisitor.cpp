// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  FlatBoundaryStressVisitor.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 9/29/14.
//  Updated by Hossein Agheshlui on 09 April 16.
//

#include <limits>

#include "CSMP_definitions.h"
#include "FlatBoundaryStressVisitor.h"
#include "Element.h"
#include "Model.h"
#include "Exception.h"
#include "Boundary.h"
#include "Face.h"
#include "InterFace.h"
#include "VectorVariable.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/** 
    For the computation of the nodal forces produced by boundary stresses

    Application level is model so that all boundaries are visited
    followed by a visitation of their faces.
*/
template<uint32_t dim>
FlatBoundaryStressVisitor<dim>::FlatBoundaryStressVisitor( const Model<dim>& model,
                                                           bool overwrite_force_vector )
    : Visitor<dim>( BOUNDARY, FACE ),
      Sn_key_(model.Database().StorageKey("normal stress")),
      Ss_key_(model.Database().StorageKey("shear stress")),
      F_key_(model.Database().StorageKey("force")),
      overwrite_previous_forces_(overwrite_force_vector)

{
    if ( Sn_key_.place != FACE || Sn_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "FlatBoundaryStressVisitor (constructor):",
                              "'normal stress' must be a SCALAR face variable." );

    if ( Ss_key_.place != FACE || Ss_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "FlatBoundaryStressVisitor (constructor):",
                              "'shear stress' must be a VECTOR face variable." );

    if ( F_key_.place != NODE || F_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "FlatBoundaryStressVisitor (constructor):",
                              "'force' must be a VECTOR variable placed on the nodes." );
     
}



/**
    To zero out initial values assigned to the boundary.
*/
template<uint32_t dim>
void FlatBoundaryStressVisitor<dim>::Visit( Boundary<dim>* b )
{
   if ( overwrite_previous_forces_ ) {
         std::cout <<"\nFlatBoundaryStressVisitor<dim>::Visit(boundary): resetting previous 'force' conditions.\n";
         VectorVariable<dim> zero;
         zero = 0.;
         b->InputPropertyValue( "force", zero, COMPLETE );
      }
}


/**
    This method interpolates the applied 
    element-centric forces to the nodes.
 
    @attention this is an integration
    (for uniform nodal forces, see option,
    we get non-uniform stresses)
 
*/
template<uint32_t dim>
void FlatBoundaryStressVisitor<dim>::Visit( Face<dim>* f )
{
   // 1. NORMAL STRESS COMPONENT
   // --------------------------

   const size_t elemType(f->InnerParent()->FE()->Interpolation()); //1 linear, 2 Quadratic
   if (elemType!=1 && elemType!=2)  {
      throw csmp::Exception(ERROR, "FlatBoundaryStressVisitor::Visit",
         "element type is not linear or quadratic.");
   }

   double Sn = f->Read( Sn_key_ ); // multiplied by -1. to make the compressive stresses positive
   // the shear stress in the plane of the face
   f->Read(Ss_key_, vs_);
   const double Ss_magnitude = vs_.Length();

   // stops the code if boundaries without boundary stresses (defined in the config file) are visited. 
   // This would accumulate NaN values to the rhs vector. 
   ErrorHandler&  csmp_error(ErrorHandler::Instance());
   if (isnan(Sn) || isnan(Ss_magnitude)) {
      csmp_error.Note(WARNING, "csmp::FlatBoundaryStressVisitor::Visit:", "Boundaries without boundary stresses are visited, returning NaN values.");
      terminate();
   }
   
   // if the values are negligibly small nothing needs to be done
   if (fabs(Sn) < std::numeric_limits<double>::epsilon()
      and fabs(Ss_magnitude) < std::numeric_limits<double>::epsilon())
         return;
   
   // using normal from the face. 
   f->UnitNormal(nrml_);
   nrml_ *= Sn;
   uint32_t  nodeNumbers(f->Nodes());
   uint32_t  BeginNodeNumber(0U);
   
   if (f->InnerParent()->Interpolation() == 2)  {
      nodeNumbers = f->FE()->MidSideNodes();
      BeginNodeNumber = f->FE()->CornerNodes();
   }
  
   vs_ *= (f->Area() / static_cast<double>(nodeNumbers));
   nrml_ *= (f->Area() / static_cast<double>(nodeNumbers));

   // adding the normal stress as force contribution to the nodes of the face
   // for linear elements, nodal forces are applied at corner nodes by a factor of A/nodeNumbers
   // for quadratic elements, nodal forces are applied at mid-side nodes by a factor of A/midsideNodeNumbers
   for (auto i = BeginNodeNumber; i < BeginNodeNumber+nodeNumbers; i++) {
      f->N(i)->Read(F_key_, vc_);
      vc_ += (nrml_ + vs_);
      f->N(i)->Store(F_key_, vc_);                                         // HA. Adds both shear and normal components to the nodal force vectors. 
      for (uint32_t j{0U}; j < dim; j++) f->N(i)->Status(F_key_, j, DIRICH);  // HA. Applies the forces as dirichlet boundary conditions. 
   }

   
} // end Visit(Element)


//template class FlatBoundaryStressVisitor<2U>;
template class FlatBoundaryStressVisitor<3U>;

} // end csmp
