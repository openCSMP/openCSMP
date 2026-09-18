// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVE_Sector.h"
#include "Element.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "ControlVolumeElement.h"
#include "ElementFace.h"

namespace csmp {

template<uint32_t dim>
CVE_Sector<dim>::CVE_Sector(ControlVolumeElement<dim>& cvelm, Node<dim>& nd, double lower_dimensional_width)
 : CVE_ptr_(&cvelm)
{
    GenerateSector(cvelm, nd, lower_dimensional_width);
}




template<uint32_t dim>
CVE_Sector<dim>::~CVE_Sector()
{
}



template<uint32_t dim>
void CVE_Sector<dim>::GenerateSector(ControlVolumeElement<dim>& cvelm, Node<dim>& nd, double lower_dimensional_width)
{
    
    size_t equivalent_node;
    for (equivalent_node = 0U; equivalent_node < cvelm.Nodes(); equivalent_node++)
        if (nd.Coordinate() == cvelm.N(equivalent_node)->Coordinate()) break;
   
    if (equivalent_node >= cvelm.Nodes())
        throw csmp::Exception( FATAL_ERROR, "CVE_Sector::GenerateSector",
                               "The provided node does not belong to the control volume element" );
    switch (cvelm.Dimension())
    {
        case 0U:
        {
            faceArea_.resize(cvelm.Faces());
            CVE_SectorFaceConnection_.resize(cvelm.Faces());
            patchFaceNumber_.resize(cvelm.Faces());

            volume_ = cvelm.Volume();
            
            for (uint32_t f = 0U; f < cvelm.Faces(); f++)
            {
                    CVE_SectorFaceConnection_[f] = (unsigned char) f;
                    faceArea_[f] = cvelm.Face(f)->Area();
            }
            
            totalFacetsAreaVector_ = 0.;
            
            break;
        }
        
        case 1U:
        {

            if (nd.AtBoundary() == NOT)
            {
                uint32_t face_counter(0U);
                uint32_t face_num  = std::numeric_limits<uint32_t>::max();
                        
                // determining if the line sector has a 0D or 1D neighbor at the node nd
                for (uint32_t n = 0U; n < cvelm.Neighbors(); n++)
                {
                    if (((cvelm.Neighbor(n)->Dimension() == 0U) && (cvelm.Neighbor(n)->N(0U) == &nd)) ||
                        ((cvelm.Neighbor(n) != &cvelm) && (cvelm.Neighbor(n)->Dimension() == 1U) && (cvelm.Neighbor(n)->N(0U) == &nd)) ||
                        ((cvelm.Neighbor(n) != &cvelm) && (cvelm.Neighbor(n)->Dimension() == 1U) && (cvelm.Neighbor(n)->N(1U) == &nd)))
                    {
                        face_counter++;
                        face_num = n;
                        break;
                    }
                    
                }
                
                faceArea_.resize(2U + face_counter);
                CVE_SectorFaceConnection_.resize(2U + face_counter);
                patchFaceNumber_.resize(2U + face_counter);
                
                volume_ = 0.5 * cvelm.Volume();

                Point<dim> outside_node(&nd == cvelm.N(0U) ? 
                                        cvelm.N(1U)->Coordinate() : cvelm.N(0U)->Coordinate());
                totalFacetsAreaVector_ = outside_node - nd.Coordinate();
                totalFacetsAreaVector_.NormalizeLengthTo(lower_dimensional_width);

                uint32_t face_num_counter(0U);
                
                for (uint32_t n = 0U; n < cvelm.Neighbors(); n++)
                {
                    if (cvelm.Neighbor(n)->Dimension() == 2U)
                    {
                        CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) n;
                        faceArea_[face_num_counter] = 0.5 * cvelm.Face(n)->Area();
                        face_num_counter++;
                    }
                    
                }

                // if sector has a face with a 0D or 1D cve
                if (face_counter == 1U)
                {
                        CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) face_num;
                        faceArea_[face_num_counter] = cvelm.Face(face_num)->Area();            
                }
                
            }
            else
            {
                uint32_t point_face_num = std::numeric_limits<uint32_t>::max();
                
                // determining the boundary or point face associated with line element
                for (uint32_t n = 0U; n < cvelm.Neighbors(); n++)
                {
                    if (((cvelm.Neighbor(n)->Dimension() == 0U) && (cvelm.Neighbor(n)->N(0U) == &nd)) ||
                        (cvelm.Neighbor(n) == &cvelm))
                    {
                        point_face_num = n;
                        break;
                    }
                    
                }

                faceArea_.resize(3U);
                CVE_SectorFaceConnection_.resize(3U);
                patchFaceNumber_.resize(3U);
                
                volume_ = 0.5 * cvelm.Volume();

                Point<dim> outside_node(&nd == cvelm.N(0U) ? 
                                        cvelm.N(1U)->Coordinate() : cvelm.N(0U)->Coordinate());
                totalFacetsAreaVector_ = outside_node - nd.Coordinate();
                totalFacetsAreaVector_.NormalizeLengthTo(lower_dimensional_width);

                uint32_t face_num_counter{0U};
                
                for (uint32_t n = 0U; n < cvelm.Neighbors(); n++)
                {
                    if (cvelm.Neighbor(n)->Dimension() == 2U)
                    {
                        CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) n;
                        faceArea_[face_num_counter] = 0.5 * cvelm.Face(n)->Area();
                        face_num_counter++;
                    }
                    
                }

                // inserting the last face
                CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) point_face_num;
                faceArea_[face_num_counter] = cvelm.Face(point_face_num)->Area();            

            }

            break;
        }
        
        case 2U:
        {
            faceArea_.resize(2U);
            CVE_SectorFaceConnection_.resize(2U);
            patchFaceNumber_.resize(2U);
            
            Point<dim> elm_barycenter(cvelm.E()->BaryCenter());
            Point<dim> right_node_coordinate;
            Point<dim> left_node_coordinate;
            size_t right_node_num;
            size_t left_node_num;
            
            if (equivalent_node == 0U)
            {
                right_node_coordinate = cvelm.N(1U)->Coordinate();
                left_node_coordinate = cvelm.N(cvelm.Nodes() - 1U)->Coordinate();
                right_node_num = 1U;
                left_node_num = cvelm.Nodes() - 1U;
            }
            else if (equivalent_node == cvelm.Nodes() - 1U)
            {
                right_node_coordinate = cvelm.N(0U)->Coordinate();
                left_node_coordinate = cvelm.N(cvelm.Nodes() - 2U)->Coordinate();
                right_node_num = 0U;
                left_node_num = cvelm.Nodes() - 2U;
            }
            else
            {
                right_node_coordinate = cvelm.N(equivalent_node + 1U)->Coordinate();
                left_node_coordinate = cvelm.N(equivalent_node - 1U)->Coordinate();
                right_node_num = equivalent_node + 1U;
                left_node_num = equivalent_node - 1U;
            }

            volume_ = std::abs(crossProduct(elm_barycenter - nd.Coordinate(), 
                      (right_node_coordinate - nd.Coordinate()) / 2.)[0U]);
            volume_ += std::abs(crossProduct(elm_barycenter - nd.Coordinate(), 
                       (left_node_coordinate - nd.Coordinate()) / 2.)[0U]);
            volume_ *= 0.5;
            
            // calculating total facet area vector
            Point<dim> right_mid_point((right_node_coordinate + nd.Coordinate()) / 2.);
            Point<dim> left_mid_point((left_node_coordinate + nd.Coordinate()) / 2.);
            
            Point<dim> temp;
            temp[0U] = elm_barycenter[1U] - right_mid_point[1U];
            temp[1U] = -(elm_barycenter[0U] - right_mid_point[0U]);
            temp.NormalizeLengthTo(elm_barycenter.DistanceTo(right_mid_point));
            totalFacetsAreaVector_ = temp;
            
            temp[0U] = left_mid_point[1U] - elm_barycenter[1U];
            temp[1U] = -(left_mid_point[0U] - elm_barycenter[0U]);
            temp.NormalizeLengthTo(elm_barycenter.DistanceTo(left_mid_point));
            totalFacetsAreaVector_ += temp;

            // finding the faces for the sector
            size_t face_num_counter(0U);
            
            for (size_t f = 0U; f < cvelm.Faces(); f++)
            {
                if (cvelm.Face(f)->Placement() == NOT)
                {
                    ControlVolumeElement<dim>& neighbor(cvelm.Face(f)->InsideCVE() == &cvelm ?
                                                        *(cvelm.Face(f)->OutsideCVE()) : 
                                                        *(cvelm.Face(f)->InsideCVE()));
                    for (size_t n = 0U; n < neighbor.Nodes(); n++)
                    {
                        if (neighbor.N(n)->Coordinate() == nd.Coordinate())
                        {
                            for (size_t i = 0U; i < neighbor.Nodes(); i++)
                            {
                                if ((neighbor.N(i)->Coordinate() == right_node_coordinate) ||
                                    (neighbor.N(i)->Coordinate() == left_node_coordinate))
                                {
                                    CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) f;
                                    faceArea_[face_num_counter] = 0.5 * cvelm.Face(f)->Area();
                                    face_num_counter++;
                                    break;
                                }
                                
                            }
                            
                            break;
                        
                        }
                        
                    }
                    
                }
                else
                {
                    if (isLEFT(cvelm.Face(f)->Placement()) && isLEFT(nd.AtBoundary()) &&
                        (isLEFT(cvelm.N(right_node_num)->AtBoundary()) ||
                         isLEFT(cvelm.N(left_node_num)->AtBoundary())))
                        {
                            CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) f;
                            faceArea_[face_num_counter] = 0.5 * cvelm.Face(f)->Area();
                            face_num_counter++;
                            continue;
                        }
                        
                    if (isRIGHT(cvelm.Face(f)->Placement()) && isRIGHT(nd.AtBoundary()) &&
                        (isRIGHT(cvelm.N(right_node_num)->AtBoundary()) ||
                         isRIGHT(cvelm.N(left_node_num)->AtBoundary())))
                        {
                            CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) f;
                            faceArea_[face_num_counter] = 0.5 * cvelm.Face(f)->Area();
                            face_num_counter++;
                            continue;
                        }
                        
                    if (isTOP(cvelm.Face(f)->Placement()) && isTOP(nd.AtBoundary()) &&
                        (isTOP(cvelm.N(right_node_num)->AtBoundary()) ||
                         isTOP(cvelm.N(left_node_num)->AtBoundary())))
                        {
                            CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) f;
                            faceArea_[face_num_counter] = 0.5 * cvelm.Face(f)->Area();
                            face_num_counter++;
                            continue;
                        }
                        
                    if (isBOTTOM(cvelm.Face(f)->Placement()) && isBOTTOM(nd.AtBoundary()) &&
                        (isBOTTOM(cvelm.N(right_node_num)->AtBoundary()) ||
                         isBOTTOM(cvelm.N(left_node_num)->AtBoundary())))
                        {
                            CVE_SectorFaceConnection_[face_num_counter] = (unsigned char) f;
                            faceArea_[face_num_counter] = 0.5 * cvelm.Face(f)->Area();
                            face_num_counter++;
                            continue;
                        }

                }
                
                
            }
            
            break;
        }
            
        default:
            throw csmp::Exception( FATAL_ERROR, "CVE_Sector::GenerateSector",
                               "Unsupported control volume element dimension" );
            break;
    }

}




template class CVE_Sector<2U>;

} // end namespace csmp
