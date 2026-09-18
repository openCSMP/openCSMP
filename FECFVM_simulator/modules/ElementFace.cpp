// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ElementFace.h"
#include "ControlVolumeElement.h"
#include "Element.h"
#include "Node.h"
#include "Exception.h"


namespace csmp {


// initializing the static elementID_Counter at file scope
template<uint32_t dim>
size_t ElementFace<dim>::elementFaceID_Counter_(0U);



template<uint32_t dim>
ElementFace<dim>::ElementFace()
: active_(false)
{
}



template<uint32_t dim>
ElementFace<dim>::~ElementFace()
{
}
    


template<uint32_t dim>
ElementFace<dim>::ElementFace(ControlVolumeElement<dim>& CVE1, ControlVolumeElement<dim>& CVE2, size_t local_property_storage_size, BOX_BOUNDARY place)
: active_(false)
{
    ConstructFace(CVE1, CVE2, local_property_storage_size, place);
}



template<uint32_t dim>
void ElementFace<dim>::ConstructFace(ControlVolumeElement<dim>& CVE1, ControlVolumeElement<dim>& CVE2, size_t local_property_storage_size, BOX_BOUNDARY place)
{

    if ((CVE1.Dimension() == 0U && CVE2.Dimension() == 0U) && (place == NOT))
        throw csmp::Exception( FATAL_ERROR, "ElementFace<dim>::ConstructFace",
              "Interior zero dimentional Control Volume Element can not have a zero dimensional neighbor!" );

    if ((CVE1.Dimension() == 0U && CVE2.Dimension() == 2U) || (CVE1.Dimension() == 2U && CVE2.Dimension() == 0U))
        throw csmp::Exception( FATAL_ERROR, "ElementFace<dim>::ConstructFace",
              "Zero dimentional Control Volume Element can not have a two dimensional neighbor!" );
                 
    elementFaceID_ = elementFaceID_Counter_++;
    
    localPropertyStorage_.resize(local_property_storage_size);
    for (size_t i = 0; i < localPropertyStorage_.size(); i++)
        localPropertyStorage_[i] = std::strtod("NAN",NULL); 
        
    if (CVE1.ID() > CVE2.ID())
    {
        insideControlVolumeElement_ = &CVE1;
        outsideControlVolumeElement_ = &CVE2;
    }
    else
    {
        insideControlVolumeElement_ = &CVE2;
        outsideControlVolumeElement_ = &CVE1;
    }

    place_ = place;
    
    Point<dim> right_node_coordinate;
    Point<dim> left_node_coordinate;
    
    ControlVolumeElement<dim>* temp_inside(CVE1.Dimension() >= CVE2.Dimension() ? &CVE1 : &CVE2);
    ControlVolumeElement<dim>* temp_outside(CVE1.Dimension() >= CVE2.Dimension() ? &CVE2 : &CVE1);
    double normal_multiplier(temp_inside->ID() == insideControlVolumeElement_->ID() ? 1. : -1.);  
    
    if (place == NOT)
    {
        
        if (temp_inside->Dimension() == 2U)
        {
            std::vector<size_t> shared_nodes;
            
            for (size_t i = 0U; i < temp_inside->Nodes(); i++)
            {
                for (size_t j = 0U; j < temp_outside->Nodes(); j++)
                {
                    if (temp_inside->N(i)->Coordinate() == temp_outside->N(j)->Coordinate())
                    {
                        shared_nodes.push_back(i);
                        break;
                    }
                        
                }
                
            }

            if (shared_nodes[1U] - shared_nodes[0U] == 1U)
            {
                right_node_coordinate = temp_inside->N(shared_nodes[0U])->Coordinate();
                left_node_coordinate = temp_inside->N(shared_nodes[1U])->Coordinate();
            }
            else
            {
                right_node_coordinate = temp_inside->N(shared_nodes[1U])->Coordinate();
                left_node_coordinate = temp_inside->N(shared_nodes[0U])->Coordinate();
            }            

            elementFaceArea_ = right_node_coordinate.DistanceTo(left_node_coordinate);
            
            faceUnitNormal_[0U] = left_node_coordinate[1U] - right_node_coordinate[1U];
            faceUnitNormal_[1U] = right_node_coordinate[0U] - left_node_coordinate[0U];
            faceUnitNormal_.NormalizeLengthTo(1.);                
            faceUnitNormal_ *= normal_multiplier;
        
        } // end if inside element is a 2D element
        // if inside element is a 1D element
        else
        {
            // if outside element is a 1D element
            if (temp_outside->Dimension() == 1U)
            {
                Point<dim> inside_node;
                Point<dim> outside_node;
                Point<dim> shared_node;
                
                for (size_t i = 0U; i < temp_inside->Nodes(); i++)
                {
                    for (size_t j = 0U; j < temp_outside->Nodes(); j++)
                    {
                        if (temp_inside->N(i)->Coordinate() == temp_outside->N(j)->Coordinate())
                        {
                            shared_node = temp_inside->N(i)->Coordinate();
                            break;
                        }
                        
                    }
                    
                }
                
                inside_node = (temp_inside->N(0)->Coordinate() == shared_node ? temp_inside->N(1)->Coordinate() : 
                                                                                temp_inside->N(0)->Coordinate());
                outside_node = (temp_outside->N(0)->Coordinate() == shared_node ? temp_outside->N(1)->Coordinate() : 
                                                                                  temp_outside->N(0)->Coordinate());
                
                Point<dim> temp;
                
                // face normal is average of direction of two lines that meet
                temp = shared_node - inside_node;
                temp.NormalizeLengthTo(0.5);
                
                faceUnitNormal_ = temp;
                
                temp = outside_node - shared_node;
                temp.NormalizeLengthTo(0.5);
                
                faceUnitNormal_ += temp;
                faceUnitNormal_ *= normal_multiplier;
                
                // this should be modified via the interface of the class since for lower dimentional elements
                // the thickness should be determined as an input parameter
                elementFaceArea_ = std::strtod("NAN",NULL);
                
            }
            // if outside element is a 0D element
            else
            {
                if (temp_outside->N(0U)->Coordinate() == temp_inside->N(0U)->Coordinate())
                {
                    faceUnitNormal_ = temp_inside->N(0U)->Coordinate() - temp_inside->N(1U)->Coordinate();
                }
                else
                {
                    faceUnitNormal_ = temp_inside->N(1U)->Coordinate() - temp_inside->N(0U)->Coordinate();
                }
                
                faceUnitNormal_.NormalizeLengthTo(1.);                
                faceUnitNormal_ *= normal_multiplier;
                
                // this should be modified via the interface of the class since for lower dimentional elements
                // the thickness should be determined as an input parameter
                elementFaceArea_ = std::strtod("NAN",NULL);
            
            }
            
        } // end if inside element is a 1D element


    } // end if it is an inside element (palce = NOT)
    else
    {
        if (temp_inside->Dimension() == 2U)
        {
            std::vector<size_t> shared_nodes;
            
            for (size_t i = 0U; i < temp_inside->Nodes(); i++)
            {
                if (isLEFT(place) && isLEFT(temp_inside->N(i)->AtBoundary()))
                {
                    shared_nodes.push_back(i);
                    continue;
                }

                if (isRIGHT(place) && isRIGHT(temp_inside->N(i)->AtBoundary()))
                {
                    shared_nodes.push_back(i);
                    continue;
                }

                if (isTOP(place) && isTOP(temp_inside->N(i)->AtBoundary()))
                {
                    shared_nodes.push_back(i);
                    continue;
                }

                if (isBOTTOM(place) && isBOTTOM(temp_inside->N(i)->AtBoundary()))
                {
                    shared_nodes.push_back(i);
                    continue;
                }

            }

            if (shared_nodes[1U] - shared_nodes[0U] == 1U)
            {
                right_node_coordinate = temp_inside->N(shared_nodes[0U])->Coordinate();
                left_node_coordinate = temp_inside->N(shared_nodes[1U])->Coordinate();
            }
            else
            {
                right_node_coordinate = temp_inside->N(shared_nodes[1U])->Coordinate();
                left_node_coordinate = temp_inside->N(shared_nodes[0U])->Coordinate();
            }
            
            elementFaceArea_ = right_node_coordinate.DistanceTo(left_node_coordinate);
            
            faceUnitNormal_[0U] = left_node_coordinate[1U] - right_node_coordinate[1U];
            faceUnitNormal_[1U] = right_node_coordinate[0U] - left_node_coordinate[0U];
            faceUnitNormal_.NormalizeLengthTo(1.);                
                        
        } // end if it is a 2D element
        
        // if boundary element is a line element
        if (temp_inside->Dimension() == 1U)
        {
            // if line element is a boundary element
            if ((temp_inside->N(0U)->AtBoundary() != NOT) && (temp_inside->N(1U)->AtBoundary() != NOT))
            {
                if (isLEFT(place))
                {
                    if (temp_inside->N(0U)->Coordinate()[1U] > temp_inside->N(1U)->Coordinate()[1U])
                    {
                        right_node_coordinate = temp_inside->N(0U)->Coordinate();
                        left_node_coordinate = temp_inside->N(1U)->Coordinate();
                    }
                    else
                    {
                        right_node_coordinate = temp_inside->N(1U)->Coordinate();
                        left_node_coordinate = temp_inside->N(0U)->Coordinate();
                    }
                    
                }
                if (isRIGHT(place))
                {
                    if (temp_inside->N(0U)->Coordinate()[1U] < temp_inside->N(1U)->Coordinate()[1U])
                    {
                        right_node_coordinate = temp_inside->N(0U)->Coordinate();
                        left_node_coordinate = temp_inside->N(1U)->Coordinate();
                    }
                    else
                    {
                        right_node_coordinate = temp_inside->N(1U)->Coordinate();
                        left_node_coordinate = temp_inside->N(0U)->Coordinate();
                    }
                    
                }
                if (isTOP(place))
                {
                    if (temp_inside->N(0U)->Coordinate()[0U] > temp_inside->N(1U)->Coordinate()[0U])
                    {
                        right_node_coordinate = temp_inside->N(0U)->Coordinate();
                        left_node_coordinate = temp_inside->N(1U)->Coordinate();
                    }
                    else
                    {
                        right_node_coordinate = temp_inside->N(1U)->Coordinate();
                        left_node_coordinate = temp_inside->N(0U)->Coordinate();
                    }
                    
                }
                if (isBOTTOM(place))
                {
                    if (temp_inside->N(0U)->Coordinate()[0U] < temp_inside->N(1U)->Coordinate()[0U])
                    {
                        right_node_coordinate = temp_inside->N(0U)->Coordinate();
                        left_node_coordinate = temp_inside->N(1U)->Coordinate();
                    }
                    else
                    {
                        right_node_coordinate = temp_inside->N(1U)->Coordinate();
                        left_node_coordinate = temp_inside->N(0U)->Coordinate();
                    }
                    
                }
                
                elementFaceArea_ = right_node_coordinate.DistanceTo(left_node_coordinate);
                
                faceUnitNormal_[0U] = left_node_coordinate[1U] - right_node_coordinate[1U];
                faceUnitNormal_[1U] = right_node_coordinate[0U] - left_node_coordinate[0U];
                faceUnitNormal_.NormalizeLengthTo(1.);                
                
            } // end if line element is a boundary element
            // if line element has a node at a boundary
            else
            {
                Point<dim> interior_node(temp_inside->N(0U)->AtBoundary() == NOT ? temp_inside->N(0U)->Coordinate() : 
                                                                                   temp_inside->N(1U)->Coordinate());
                Point<dim> boundary_node(temp_inside->N(1U)->AtBoundary() == NOT ? temp_inside->N(0U)->Coordinate() : 
                                                                                   temp_inside->N(1U)->Coordinate());
                faceUnitNormal_ = boundary_node - interior_node;
                faceUnitNormal_.NormalizeLengthTo(1.);
                                
                // this should be modified via the interface of the class since for lower dimentional elements
                // the thickness should be determined as an input parameter
                elementFaceArea_ = std::strtod("NAN",NULL);
            
            } // end if line element has a node at a boundary
            
        } // end if boundary element is a line element

        // if boundary element is a point element
        if (temp_inside->Dimension() == 0U)
        {
            if (isLEFT(place)) faceUnitNormal_ = Point<dim>(-1., 0.);
            if (isRIGHT(place)) faceUnitNormal_ = Point<dim>(1., 0.);
            if (isTOP(place)) faceUnitNormal_ = Point<dim>(0., 1.);
            if (isBOTTOM(place)) faceUnitNormal_ = Point<dim>(0., -1.);
            
            // this should be modified via the interface of the class since for lower dimentional elements
            // the thickness should be determined as an input parameter
            elementFaceArea_ = std::strtod("NAN",NULL);
            
        } // end if boundary element is a point element
        
    } // end if it is a boundary element (palce != NOT)
    
}



template class ElementFace<2U>;


} // end csmp
