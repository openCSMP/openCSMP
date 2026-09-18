// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "InterFaceFractureVisitor.h"
#include "Model.h"
#include "InterFace.h"
#include "Exception.h"


using namespace std;

namespace csmp {

template<uint32_t dim>
InterFaceFractureVisitor<dim>::InterFaceFractureVisitor( const PropertyDatabase<dim>& pref,
                                                         const char* displacement,
                                                         const char* aperture,
                                                         INTERFACE_SIDE aperture_side)
  : Visitor<dim>(SPLIT_BOUNDARY, INTER_FACE),
    aperture_key_(pref.StorageKey(aperture)),
    displacement_key_(pref.StorageKey(displacement)),
    aperture_side_(aperture_side),
    calculate_conductivity_(false)
{
  if (aperture_key_.place != NODE or displacement_key_.place != NODE)
    throw csmp::Exception(ERROR, "InterFaceApertureVisitor (constructor)",
                          "Apertures or Displacments must be defined on the nodes");
  if (aperture_key_.type != SCALAR or displacement_key_.type != VECTOR)
    throw csmp::Exception(ERROR, "InterFaceApertureVisitor (constructor)",
                          "aperture must be scalar, and displacement must be vector");
}

template<uint32_t dim>
InterFaceFractureVisitor<dim>::InterFaceFractureVisitor( const PropertyDatabase<dim>& pref,
                                                         const char* displacement,
                                                         const char* aperture,
                                                         INTERFACE_SIDE aperture_side,
                                                         const char* conductivity,
                                                         const char* viscosity,
                                                         double min_aperture)
  : Visitor<dim>(SPLIT_BOUNDARY, INTER_FACE),
    aperture_key_(pref.StorageKey(aperture)),
    displacement_key_(pref.StorageKey(displacement)),
    aperture_side_(aperture_side),
    conductivity_key_(pref.StorageKey(conductivity)),
    viscosity_key_(pref.StorageKey(viscosity)),
    min_aperture_(min_aperture),
    calculate_conductivity_(true)
{
  if (aperture_key_.place != NODE or displacement_key_.place != NODE)
    throw csmp::Exception(ERROR, "InterFaceApertureVisitor (constructor)",
                          "Apertures or Displacments must be defined on the nodes");
  if (aperture_key_.type != SCALAR or displacement_key_.type != VECTOR)
    throw csmp::Exception(ERROR, "InterFaceApertureVisitor (constructor)",
                          "aperture must be scalar, and displacement must be vector");

  if (conductivity_key_.place != NODE or viscosity_key_.place != NODE)
    throw csmp::Exception(ERROR, "InterFaceApertureVisitor (constructor)",
                          "viscosity or conductivity must be defined on the nodes");

  if (conductivity_key_.type != SCALAR or viscosity_key_.type != SCALAR)
    throw csmp::Exception(ERROR, "InterFaceApertureVisitor (constructor)",
                          "conductivity must be scalar, and viscosity must be scalar");


}


template<uint32_t dim>
InterFaceFractureVisitor<dim>::~InterFaceFractureVisitor()
{}


/**
    Takes the displacement defined on the nodes of the INSIDE and OUTSIDE
    faces of the InterFace object, dots them with the Interface unit normal
    (which by convention points inside to outside) to obtain the aperture,
    and stores the aperture on the middle Nodes of the InterFace object.

*/
template<uint32_t dim>
void InterFaceFractureVisitor<dim>::Visit( InterFace<dim>* f )
{

  if (aperture_side_ == MIDDLE && !f->HasInterveningElement())
    throw csmp::Exception( ERROR, "InterFaceFractureVisitor::Visit( interface )",
                           "Aperture property is designated to be stored on the MIDDLE, but MIDDLE Element does not exist in InterFace object.");

  if (calculate_conductivity_ && !f->HasInterveningElement())
    throw csmp::Exception( ERROR, "InterFaceFractureVisitor::Visit( interface )",
                           "Conductivity property is expected to be stored on the MIDDLE, but MIDDLE Element does not exist in InterFace object.");

  std::vector<VectorVariable<dim>> disp_plus, disp_minus;
  f->MatchingNodePropertyVector( displacement_key_, disp_plus, OUTSIDE);
  f->MatchingNodePropertyVector( displacement_key_, disp_minus, INSIDE);
  VectorVariable<dim> nrml(f->UnitNormal());

  if ( !calculate_conductivity_ ){    //Just calculate aperture

    for ( uint32_t i = 0u; i < f->FE()->Nodes() ; ++i){
      //Calculate aperture based on displacement
      double aperture         = nrml.DotProduct(disp_plus[i] - disp_minus[i]);
      //get node for storage
      Node<dim>* n = f->MatchingN(i,aperture_side_);
      //store new aperture value
      n->Store(aperture_key_, ScalarVariable(PLAIN, aperture));
    }

  } else {      //Calculate both aperture and conductivity

      ScalarVariable viscosity;
      for ( uint32_t i = 0u; i < f->FE()->Nodes() ; ++i){
        //Calculate aperture based on displacement
        double aperture         = nrml.DotProduct(disp_plus[i] - disp_minus[i]);
        //Store aperture on side
        Node<dim>* n = f->MatchingN(i,aperture_side_);
        n->Store(aperture_key_, ScalarVariable(PLAIN, aperture));    //store new aperture value

        //Read viscosity
        Node<dim>* n_middle = f->MatchingN(i,MIDDLE);
        n_middle->Read(viscosity_key_, viscosity);
        if (aperture < min_aperture_)
          aperture=min_aperture_;                 //setting to minimum aperture only for conductivity calculation
        //store conductivity
        n_middle->Store( conductivity_key_, ScalarVariable(PLAIN, aperture*aperture*aperture/12.0/viscosity())); //CUBIC LAW

      }
    }

}//end of Visit


template class InterFaceFractureVisitor<1U>;
template class InterFaceFractureVisitor<2U>;
template class InterFaceFractureVisitor<3U>;


} // csmp

