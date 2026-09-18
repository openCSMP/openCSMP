// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ControlVolumeElement.h"
#include "ElementFace.h"
#include "Element.h"
#include "Node.h"
#include "Exception.h"

namespace csmp {


// initializing the static elementID_Counter at file scope
template<uint32_t dim>
size_t ControlVolumeElement<dim>::controlVolumeElementID_Counter_(0U);
  
  
template<uint32_t dim>
ControlVolumeElement<dim>::ControlVolumeElement()
: active_(false)
{
}
    


template<uint32_t dim>    
ControlVolumeElement<dim>::~ControlVolumeElement()
{
}



template<uint32_t dim>
ControlVolumeElement<dim>::ControlVolumeElement(Element<dim>* element, size_t no_faces, size_t local_property_storage_size)
: active_(false)
{
    ConstructControlVolumeElement(element, no_faces, local_property_storage_size);
}



template<uint32_t dim>
ControlVolumeElement<dim>::ControlVolumeElement(Node<dim>* node, size_t no_faces, size_t local_property_storage_size)
: active_(false)
{
    ConstructControlVolumeElement(node, no_faces, local_property_storage_size);
}



template<uint32_t dim>
void ControlVolumeElement<dim>::ConstructControlVolumeElement(Element<dim>* element, size_t no_faces, size_t local_property_storage_size)
{

    if (element->FE()->IsSurface())
        elementDimension_ = 2;
    else if (element->FE()->IsLine())
        elementDimension_ = 1;
    else
        throw csmp::Exception( FATAL_ERROR, "ControlVolumeElement<dim>::ConstructControlVolumeElement",
                      "Unrecognized element type!" );
    
    controlVolumeElementID_ = controlVolumeElementID_Counter_++;
    elementPtr_ = element;
    nodePtr_ = 0;
    
    volume_ = element->Volume();
    
    controlVolumeElementFaces_.resize(no_faces);
    controlVolumeElementFacesNormalDirection_.resize(no_faces);
    for (size_t i = 0U; i < controlVolumeElementFaces_.size(); i++)
    {
        controlVolumeElementFaces_[i] = 0;
        controlVolumeElementFacesNormalDirection_[i] = 0; // formerly: numeric_limits<double>::quiet_NaN();
    }


    localPropertyStorage_.resize(local_property_storage_size);
    for (size_t i = 0U; i < localPropertyStorage_.size(); i++)
        localPropertyStorage_[i] = 0; // formerly: std::numeric_limits<double>::quiet_NaN();
}



template<uint32_t dim>
void ControlVolumeElement<dim>::ConstructControlVolumeElement(Node<dim>* node, size_t no_faces, size_t local_property_storage_size)
{
    
    controlVolumeElementID_ = controlVolumeElementID_Counter_++;
    elementPtr_ = 0;
    nodePtr_ = node;
    elementDimension_ = 0;    
    volume_ = 0.;
    
    controlVolumeElementFaces_.resize(no_faces);
    controlVolumeElementFacesNormalDirection_.resize(no_faces);
    for (size_t i = 0U; i < controlVolumeElementFaces_.size(); i++)
    {
        controlVolumeElementFaces_[i] = 0;
        controlVolumeElementFacesNormalDirection_[i] = std::numeric_limits<double>::quiet_NaN();
    }

    localPropertyStorage_.resize(local_property_storage_size);
    for (size_t i = 0U; i < localPropertyStorage_.size(); i++)
        localPropertyStorage_[i] = std::numeric_limits<double>::quiet_NaN();
}



template class ControlVolumeElement<2U>;


} // end csmp
