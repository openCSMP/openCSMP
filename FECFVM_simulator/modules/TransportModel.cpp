#include "TransportModel.h"
#include "Model.h"
#include "Node.h"
#include "CVE_Sector.h"
#include "ErrorHandler.h"
#include "Exception.h"

namespace csmp {


template<uint32_t dim>
TransportModel<dim>::TransportModel(Model<dim>& model, const char* region_name, const char* vol_mod,
                                    size_t element_local_storage_size, size_t boundary_element_local_storage_size,
                                    size_t face_local_storage_size, size_t boundary_face_local_storage_size)
 : modelRef_(model),
   regionName_(region_name),
   CVE_Collection_(0U),
   faceCollection_(0U)
{

    std::cout << "\nTransportModel<" << dim << ">::Constructor: Building the transport model.\n";
    
    CreateControlVolumeElementsAndFaces(regionName_.c_str(), element_local_storage_size, boundary_element_local_storage_size,
                                        face_local_storage_size, boundary_face_local_storage_size);
    CorrectVolumes(vol_mod);
    CreatePatches(vol_mod);
    AssignParentPatchesToCVEs();
    
    TransportModelReport();
    
}




template<uint32_t dim>
TransportModel<dim>::~TransportModel()
{
    std::cout << "\nDestruct Transport Model...\n";
}





template<uint32_t dim>
void TransportModel<dim>::CreateControlVolumeElementsAndFaces(const char* region_name, 
                                                              size_t element_local_storage_size, 
                                                              size_t boundary_element_local_storage_size,
                                                              size_t face_local_storage_size, 
                                                              size_t boundary_face_local_storage_size)
{
    std::cout << "\nTransportModel<" << dim << ">::CreateControlVolumeElementsAndFaces()\n";
    
    Region<2U>& mref = modelRef_.Region(region_name);
    
    std::set<Element<dim>*> interior_2D_elements;
    std::set<Element<dim>*> boundary_2D_elements;
    std::set<Element<dim>*> interior_1D_elements;
    std::set<Element<dim>*> boundary_1D_elements;
    std::set<Node<dim>*> interior_0D_elements;
    std::set<Node<dim>*> boundary_0D_elements;

    //  node min ID      node max ID            elm_ptr
    std::map<size_t, std::map<size_t, std::set<Element<dim>*> > > shared_faces_map;
    
    for ( auto eit = mref.CellsBegin(); eit != mref.CellsEnd(); eit++ )
    {
        if ( (*eit)->IsSurface() )
        {
            size_t left_boundary_nodes(0U);
            size_t right_boundary_nodes(0U);
            size_t top_boundary_nodes(0U);
            size_t bottom_boundary_nodes(0U);
            
            for (uint32_t n = 0U; n < (*eit)->Nodes(); n++)
            {
                if (isLEFT((*eit)->N(n)->AtBoundary())) left_boundary_nodes++;
                if (isRIGHT((*eit)->N(n)->AtBoundary())) right_boundary_nodes++;
                if (isTOP((*eit)->N(n)->AtBoundary())) top_boundary_nodes++;
                if (isBOTTOM((*eit)->N(n)->AtBoundary())) bottom_boundary_nodes++;
                
                uint32_t node_1(n);
                uint32_t node_2(n == 0U ? (*eit)->Nodes() - 1U : n - 1U);
                
                if ((*eit)->N(node_1)->Idx() < (*eit)->N(node_2)->Idx())
                    shared_faces_map[(*eit)->N(node_1)->Idx()][(*eit)->N(node_2)->Idx()].insert((*eit));
                else
                    shared_faces_map[(*eit)->N(node_2)->Idx()][(*eit)->N(node_1)->Idx()].insert((*eit));
                    
            }
            
            if ((left_boundary_nodes == 2U) || (right_boundary_nodes == 2U) ||
                (top_boundary_nodes == 2U) || (bottom_boundary_nodes == 2U)) boundary_2D_elements.insert((*eit));
            else interior_2D_elements.insert((*eit));
                
            continue;
            
        }
        
        if ((*eit)->IsLine())
        {

            if ((*eit)->N(0U)->Idx() < (*eit)->N(1U)->Idx())
                shared_faces_map[(*eit)->N(0U)->Idx()][(*eit)->N(1U)->Idx()].insert((*eit));
            else
                shared_faces_map[(*eit)->N(1U)->Idx()][(*eit)->N(0U)->Idx()].insert((*eit));

            if (((*eit)->N(0U)->AtBoundary() != NOT) || ((*eit)->N(1U)->AtBoundary() != NOT))
                boundary_1D_elements.insert((*eit));
            else 
                interior_1D_elements.insert((*eit));
        }
        
    } // end looping over elements
    
    // removing the 2D boundary elements from the set and add them to the interior set
    // where there is a line element at the boundary face and line element is a boundary element

    // stores number of faces between two nodes and parent elements
    std::map<size_t, std::map<size_t, std::set<Element<dim>*> > > boundary_element_faces;
    // stores number of boundary faces for each 2D element
    std::multiset<Element<dim>*> boundary_element_multiset;
    // stores 2D boundary elements
    std::set<Element<dim>*> boundary_element_set;
    // stores duplicated 2D boundary elements which should be inside element
    std::multiset<Element<dim>*> inside_element;
    
    // inserting boundary faces for 2D elements
    for ( auto eit = boundary_2D_elements.begin(); eit != boundary_2D_elements.end(); eit++)
    {
        for (uint32_t n = 0U; n < (*eit)->Nodes(); n++)
        {
            uint32_t node_1(n);
            uint32_t node_2(n == 0U ? (*eit)->Nodes() - 1U : n - 1U);
            
            if ((isLEFT((*eit)->N(node_1)->AtBoundary()) && isLEFT((*eit)->N(node_2)->AtBoundary())) ||
                (isRIGHT((*eit)->N(node_1)->AtBoundary()) && isRIGHT((*eit)->N(node_2)->AtBoundary())) ||
                (isTOP((*eit)->N(node_1)->AtBoundary()) && isTOP((*eit)->N(node_2)->AtBoundary())) ||
                (isBOTTOM((*eit)->N(node_1)->AtBoundary()) && isBOTTOM((*eit)->N(node_2)->AtBoundary())))
            {
                boundary_element_multiset.insert(*eit);
                boundary_element_set.insert(*eit);
                
                if ((*eit)->N(node_1)->Idx() < (*eit)->N(node_2)->Idx())
                    boundary_element_faces[(*eit)->N(node_1)->Idx()][(*eit)->N(node_2)->Idx()].insert(*eit);
                else
                    boundary_element_faces[(*eit)->N(node_2)->Idx()][(*eit)->N(node_1)->Idx()].insert(*eit);
            }
        }
    }
    
    // inserting 1D boundary elements which are aligned with the boundary
    for (typename std::set<Element<dim>*>::iterator eit = boundary_1D_elements.begin(); eit != boundary_1D_elements.end(); eit++)
    {
        if ((isLEFT((*eit)->N(0U)->AtBoundary()) && isLEFT((*eit)->N(1U)->AtBoundary())) ||
            (isRIGHT((*eit)->N(0U)->AtBoundary()) && isRIGHT((*eit)->N(1U)->AtBoundary())) ||
            (isTOP((*eit)->N(0U)->AtBoundary()) && isTOP((*eit)->N(1U)->AtBoundary())) ||
            (isBOTTOM((*eit)->N(0U)->AtBoundary()) && isBOTTOM((*eit)->N(1U)->AtBoundary())))
        {
            if ((*eit)->N(0U)->Idx() < (*eit)->N(1U)->Idx())
                boundary_element_faces[(*eit)->N(0U)->Idx()][(*eit)->N(1U)->Idx()].insert(*eit);
            else
                boundary_element_faces[(*eit)->N(1U)->Idx()][(*eit)->N(0U)->Idx()].insert(*eit);
        }
    }
    
    // detecting duplicated boundary elements
    for (typename std::map<size_t, std::map<size_t, std::set<Element<dim>*> > >::iterator min_n_id = boundary_element_faces.begin();
                  min_n_id != boundary_element_faces.end(); min_n_id++)
    {
        for (typename std::map<size_t, std::set<Element<dim>*> >::iterator max_n_id = min_n_id->second.begin();
                      max_n_id != min_n_id->second.end(); max_n_id++)
        {
            std::vector<Element<dim>*> e_ptr_2D;
            std::vector<Element<dim>*> e_ptr_1D;
            
            for (typename std::set<Element<dim>*>::iterator deit = max_n_id->second.begin();
                          deit != max_n_id->second.end(); deit++)
            {
                if ((*deit)->IsSurface())
                    e_ptr_2D.push_back((*deit));
                else
                    e_ptr_1D.push_back((*deit));
            }
            
            // if there is only one 2D element
            if (e_ptr_2D.size() == 1U && e_ptr_1D.size() == 0U)
            {
                continue;
            }
            // if there are one 2D elements and one 1D element
            else if (e_ptr_2D.size() == 1U && e_ptr_1D.size() == 1U)
            {
                inside_element.insert(e_ptr_2D[0U]);
                continue;
            }
            // too many faces (maybe because of duplicated elements)
            else
            {
                std::cout << "\nTransportModel::CreateControlVolumeElementsAndFaces()\n";
                std::cout << "Too many faces were found while detecting overlapping line and surface elements at the boundary.\n";
                std::cout << "Number of 2D Element: " << e_ptr_2D.size() << "\n";
                std::cout << "Number of 1D Element: " << e_ptr_1D.size() << "\n\n";
                for (size_t i = 0U; i < e_ptr_2D.size(); i++)
                {
                    std::cout << "2D Element ID " << e_ptr_2D[i]->Idx() << " Nodes\n";
                    for (uint32_t n = 0U; n < e_ptr_2D[i]->Nodes(); n++)
                    {
                        std::cout << "\tID " << e_ptr_2D[i]->N(n)->Idx() << 
                                     "\tx = " << e_ptr_2D[i]->N(n)->Coordinate()[0U] <<
                                     "\ty = " << e_ptr_2D[i]->N(n)->Coordinate()[1U] << "\n";
                    }

                }
                for (size_t i = 0U; i < e_ptr_1D.size(); i++)
                {
                    std::cout << "1D Element ID " << e_ptr_1D[i]->Idx() << " Nodes\n";
                    for (uint32_t n = 0U; n < e_ptr_1D[i]->Nodes(); n++)
                    {
                        std::cout << "\tID " << e_ptr_1D[i]->N(n)->Idx() << 
                                     "\tx = " << e_ptr_1D[i]->N(n)->Coordinate()[0U] <<
                                     "\ty = " << e_ptr_1D[i]->N(n)->Coordinate()[1U] << "\n";
                    }

                }

                std::cout << "\n";

                throw csmp::Exception( FATAL_ERROR, "TransportModel::CreateControlVolumeElementsAndFaces()",
                                      "Too many faces were found while detecting overlapping line and surface elements at the boundary." );
            }
        }
    }
  
    // removing duplicated boundary element from boundary set and inserting in interior set
    for (typename std::set<Element<dim>*>::iterator eit = boundary_element_set.begin(); 
                  eit != boundary_element_set.end(); eit++)
    {
        if (inside_element.count(*eit) == boundary_element_multiset.count(*eit))
          {
              interior_2D_elements.insert(*eit);
              boundary_2D_elements.erase(*eit);
          }
    }
    
    
    std::map<Node<dim>*, std::set<Element<dim>*> > node_element_connector;
    
    for (typename std::set<Element<dim>*>::iterator eit = interior_1D_elements.begin(); eit != interior_1D_elements.end(); eit++)
    {
        node_element_connector[(*eit)->N(0U)].insert((*eit));
        node_element_connector[(*eit)->N(1U)].insert((*eit));
    }
    
    for (typename std::set<Element<dim>*>::iterator eit = boundary_1D_elements.begin(); eit != boundary_1D_elements.end(); eit++)
    {
        node_element_connector[(*eit)->N(0U)].insert((*eit));
        node_element_connector[(*eit)->N(1U)].insert((*eit));
    }
    

    // 0D cves only exists inside the model. at the boundaries there is no 0D cve
    for (typename std::map<Node<dim>*, std::set<Element<dim>*> >::iterator nec = node_element_connector.begin();
                  nec != node_element_connector.end(); nec++)
    {
        if (nec->first->AtBoundary() == NOT)
        {
            if (nec->second.size() > 2U) interior_0D_elements.insert(nec->first);
        }
    }
    
        
    std::map<Element<dim>*, size_t> line_elements_faces;
    std::map<Node<dim>*, size_t> point_elements_faces;
    
    for (typename std::map<Node<dim>*, std::set<Element<dim>*> >::iterator neit = node_element_connector.begin();
         neit != node_element_connector.end(); neit++)
    {
        if (neit->first->AtBoundary() == NOT)
        {
            if (neit->second.size() > 1U)
            {
                for (typename std::set<Element<dim>*>::iterator eit = neit->second.begin(); eit != neit->second.end(); eit++)
                {
                    line_elements_faces[*eit]++;
                }
                
                if (neit->second.size() > 2U)
                    point_elements_faces[neit->first] = neit->second.size();
            }
        }
        else
        {
            for (typename std::set<Element<dim>*>::iterator eit = neit->second.begin(); eit != neit->second.end(); eit++)
            {
                line_elements_faces[*eit]++;
            }
            // there is no point cve at the boundary
        }
    }

    // assigning the interior and boundary element number
    interiorSurfaceCVE_ = interior_2D_elements.size();
    boundarySurfaceCVE_ = boundary_2D_elements.size();
    interiorLineCVE_ = interior_1D_elements.size();
    boundaryLineCVE_ = boundary_1D_elements.size();
    interiorPointCVE_ = interior_0D_elements.size();
    boundaryPointCVE_ = boundary_0D_elements.size();
    
    std::map<Element<dim>*, size_t> element_CVE_map;
    std::map<Node<dim>*, size_t> node_CVE_map;
    size_t counter(0U);

    // creating 2D interior CVEs
    for (typename std::set<Element<dim>*>::iterator eit = interior_2D_elements.begin(); eit != interior_2D_elements.end(); eit++)
    {
        CVE_Collection_.push_back(csmp::ControlVolumeElement<dim> ((*eit), (*eit)->Neighbors(), element_local_storage_size));
        element_CVE_map[*eit] = counter++;
    }

    // creating 1D interior CVEs
    for (typename std::set<Element<dim>*>::iterator eit = interior_1D_elements.begin(); eit != interior_1D_elements.end(); eit++)
    {
        const size_t faces(line_elements_faces[*eit] + 2U);
        CVE_Collection_.push_back(csmp::ControlVolumeElement<dim> ((*eit), faces, element_local_storage_size));
        element_CVE_map[*eit] = counter++;
    }

    // creating 0D interior CVEs
    for (typename std::set<Node<dim>*>::iterator nit = interior_0D_elements.begin(); nit != interior_0D_elements.end(); nit++)
    {
        const size_t faces(point_elements_faces[*nit]);
        // we create one extra storage for 0D cve, since there is no variable storage for that and we store saturation in it
        CVE_Collection_.push_back(csmp::ControlVolumeElement<dim> ((*nit), faces, element_local_storage_size + 1U));
        node_CVE_map[*nit] = counter++;
    }

    // creating 2D boundary CVEs
    for (typename std::set<Element<dim>*>::iterator eit = boundary_2D_elements.begin(); eit != boundary_2D_elements.end(); eit++)
    {
        CVE_Collection_.push_back(csmp::ControlVolumeElement<dim> ((*eit), (*eit)->Neighbors(), boundary_element_local_storage_size));
        element_CVE_map[*eit] = counter++;
    }

    // creating 1D boundary CVEs
    for (typename std::set<Element<dim>*>::iterator eit = boundary_1D_elements.begin(); eit != boundary_1D_elements.end(); eit++)
    {
        size_t faces(line_elements_faces[*eit] + 2U);
        CVE_Collection_.push_back(csmp::ControlVolumeElement<dim> ((*eit), faces, boundary_element_local_storage_size));
        element_CVE_map[*eit] = counter++;
    }

    // there are no 0D boundary CVEs

    
    std::map<size_t, Node<dim>*> id_node_map;
    for ( auto nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++)
    {
        id_node_map[(*nit)->Idx()] = *nit;
    }
    
    // creating interior line faces between 2D-2D and 1D-2D cves
    for (typename std::map<size_t, std::map<size_t, std::set<Element<dim>*> > >::iterator min_n_it = shared_faces_map.begin();
         min_n_it != shared_faces_map.end(); min_n_it++)
    {
        for (typename std::map<size_t, std::set<Element<dim>*> >::iterator max_n_it = min_n_it->second.begin();
             max_n_it != min_n_it->second.end(); max_n_it++)
        {
            size_t left_boundary_nodes(0U);
            size_t right_boundary_nodes(0U);
            size_t top_boundary_nodes(0U);
            size_t bottom_boundary_nodes(0U);
            
            if (isLEFT(id_node_map[min_n_it->first]->AtBoundary())) left_boundary_nodes++;
            if (isRIGHT(id_node_map[min_n_it->first]->AtBoundary())) right_boundary_nodes++;
            if (isTOP(id_node_map[min_n_it->first]->AtBoundary())) top_boundary_nodes++;
            if (isBOTTOM(id_node_map[min_n_it->first]->AtBoundary())) bottom_boundary_nodes++;

            if (isLEFT(id_node_map[max_n_it->first]->AtBoundary())) left_boundary_nodes++;
            if (isRIGHT(id_node_map[max_n_it->first]->AtBoundary())) right_boundary_nodes++;
            if (isTOP(id_node_map[max_n_it->first]->AtBoundary())) top_boundary_nodes++;
            if (isBOTTOM(id_node_map[max_n_it->first]->AtBoundary())) bottom_boundary_nodes++;

            if ((left_boundary_nodes != 2U) && (right_boundary_nodes != 2U) &&
                (top_boundary_nodes != 2U) && (bottom_boundary_nodes != 2U))
            {
                std::vector<Element<dim>*> e_ptr_2D;
                std::vector<Element<dim>*> e_ptr_1D;
                
                for ( auto eit = max_n_it->second.begin(); eit != max_n_it->second.end(); eit++)
                {
                    if ((*eit)->IsSurface())
                        e_ptr_2D.push_back((*eit));
                    else
                        e_ptr_1D.push_back((*eit));
                }
                
                // if the face is between two 2D elements
                if (e_ptr_2D.size() == 2U && e_ptr_1D.size() == 0U)
                {
                    ControlVolumeElement<dim>& CVE1(CVE_Collection_[element_CVE_map[e_ptr_2D[0U]]]);
                    ControlVolumeElement<dim>& CVE2(CVE_Collection_[element_CVE_map[e_ptr_2D[1U]]]);
                    faceCollection_.push_back(csmp::ElementFace<dim> (CVE1, CVE2, face_local_storage_size, NOT));
                    
                    continue;
                }
                // if the faces are between two 2D elements and one 1D element
                else if (e_ptr_2D.size() == 2U && e_ptr_1D.size() == 1U)
                {
                    ControlVolumeElement<dim>& SCVE1(CVE_Collection_[element_CVE_map[e_ptr_2D[0U]]]);
                    ControlVolumeElement<dim>& SCVE2(CVE_Collection_[element_CVE_map[e_ptr_2D[1U]]]);
                    ControlVolumeElement<dim>& LCVE(CVE_Collection_[element_CVE_map[e_ptr_1D[0U]]]);
                    faceCollection_.push_back(csmp::ElementFace<dim> (SCVE1, LCVE, face_local_storage_size, NOT));
                    faceCollection_.push_back(csmp::ElementFace<dim> (SCVE2, LCVE, face_local_storage_size, NOT));
                    
                    continue;
                }
                // too many faces (maybe because of duplicated elements)
                else
                {
                    std::cout << "\nTransportModel::CreateControlVolumeElementsAndFaces()\n";
                    std::cout << "Too many faces were found while creating faces inside the transport model.\n";
                    std::cout << "Number of 2D Element: " << e_ptr_2D.size() << "\n";
                    std::cout << "Number of 1D Element: " << e_ptr_1D.size() << "\n\n";
                    for (size_t i = 0U; i < e_ptr_2D.size(); i++)
                    {
                        std::cout << "2D Element ID " << e_ptr_2D[i]->Idx() << " Nodes\n";
                        for (uint32_t n = 0U; n < e_ptr_2D[i]->Nodes(); n++)
                        {
                            std::cout << "\tID " << e_ptr_2D[i]->N(n)->Idx() << 
                                "\tx = " << e_ptr_2D[i]->N(n)->Coordinate()[0U] <<
                                "\ty = " << e_ptr_2D[i]->N(n)->Coordinate()[1U] << "\n";
                        }
                    }
                    for (size_t i = 0U; i < e_ptr_1D.size(); i++)
                    {
                        std::cout << "1D Element ID " << e_ptr_1D[i]->Idx() << " Nodes\n";
                        for (uint32_t n = 0U; n < e_ptr_1D[i]->Nodes(); n++)
                        {
                            std::cout << "\tID " << e_ptr_1D[i]->N(n)->Idx() << 
                                "\tx = " << e_ptr_1D[i]->N(n)->Coordinate()[0U] <<
                                "\ty = " << e_ptr_1D[i]->N(n)->Coordinate()[1U] << "\n";
                        }
                    }
                    std::cout << "\n";

                    throw csmp::Exception( FATAL_ERROR, "TransportModel::CreateControlVolumeElementsAndFaces()",
                                          "Too many faces were found while creating faces inside the transport model" );
                }
                
            }
            // if both nodes are on the same boundary
            else
            {
                // we check if there is a line element which in this case is the boundary element
                std::vector<Element<dim>*> e_ptr_2D;
                std::vector<Element<dim>*> e_ptr_1D;
                
                for ( auto eit = max_n_it->second.begin(); eit != max_n_it->second.end(); eit++)
                {
                    if ((*eit)->IsSurface())
                        e_ptr_2D.push_back((*eit));
                    else
                        e_ptr_1D.push_back((*eit));
                }
                
                // if only one 2D element exists at the boundary this face is boundary face and we insert it later
                if (e_ptr_2D.size() == 1U && e_ptr_1D.size() == 0U)
                {
                    continue;
                }
                // if there are a line and 2D element at the boundary, 2D element is a interior element and we create the face for it
                else if (e_ptr_2D.size() == 1U && e_ptr_1D.size() == 1U)
                {
                    ControlVolumeElement<dim>& SCVE(CVE_Collection_[element_CVE_map[e_ptr_2D[0U]]]);
                    ControlVolumeElement<dim>& LCVE(CVE_Collection_[element_CVE_map[e_ptr_1D[0U]]]);
                    faceCollection_.push_back(csmp::ElementFace<dim> (SCVE, LCVE, face_local_storage_size, NOT));
                    continue;
                }
                // too many faces (maybe because of duplicated elements)
                else
                {
                    std::cout << "\nTransportModel::CreateControlVolumeElementsAndFaces()\n";
                    std::cout << "Too many faces were found while creating faces inside the transport model at the boundary.\n";
                    std::cout << "Number of 2D Element: " << e_ptr_2D.size() << "\n";
                    std::cout << "Number of 1D Element: " << e_ptr_1D.size() << "\n\n";
                    for (size_t i = 0U; i < e_ptr_2D.size(); i++)
                    {
                        std::cout << "2D Element ID " << e_ptr_2D[i]->Idx() << " Nodes\n";
                        for (uint32_t n = 0U; n < e_ptr_2D[i]->Nodes(); n++)
                        {
                            std::cout << "\tID " << e_ptr_2D[i]->N(n)->Idx() << 
                                "\tx = " << e_ptr_2D[i]->N(n)->Coordinate()[0U] <<
                                "\ty = " << e_ptr_2D[i]->N(n)->Coordinate()[1U] << "\n";
                        }
                    }
                    for (size_t i = 0U; i < e_ptr_1D.size(); i++)
                    {
                        std::cout << "1D Element ID " << e_ptr_1D[i]->Idx() << " Nodes\n";
                        for (uint32_t n = 0U; n < e_ptr_1D[i]->Nodes(); n++)
                        {
                            std::cout << "\tID " << e_ptr_1D[i]->N(n)->Idx() << 
                                "\tx = " << e_ptr_1D[i]->N(n)->Coordinate()[0U] <<
                                "\ty = " << e_ptr_1D[i]->N(n)->Coordinate()[1U] << "\n";
                        }
                    }

                    std::cout << "\n";

                    throw csmp::Exception( FATAL_ERROR, "TransportModel::CreateControlVolumeElementsAndFaces()",
                                          "Too many faces were found while creating faces inside the transport model at the boundary" );
                }                
            
            }
            
        }
        
    } // end of creating interior line faces between 2D-2D and 1D-2D cves

    
    // creating interior line faces between 1D-1D cves (this is also categorized as a line face, and it exists only when two 
    // line elements meet. This is only and only if 2 line elements meet. For the case that more than two line elements come 
    // together, a point cve exists at the meting point and we categorize these faces as point faces
    for (typename std::map<Node<dim>*, std::set<Element<dim>*> >::iterator nec = node_element_connector.begin();
                  nec != node_element_connector.end(); nec++)
    {
        if ((nec->first->AtBoundary() == NOT) && (nec->second.size() == 2U))
        {
            typename std::set<Element<dim>*>::iterator eit = nec->second.begin();
            
            ControlVolumeElement<dim>& LCVE1(CVE_Collection_[element_CVE_map[*eit]]);
            eit++;
            ControlVolumeElement<dim>& LCVE2(CVE_Collection_[element_CVE_map[*eit]]);
            
            faceCollection_.push_back(csmp::ElementFace<dim> (LCVE1, LCVE2, face_local_storage_size, NOT));
        }
    }
    
    interiorLineFaces_ = faceCollection_.size();

    // creating interior point faces
    for (typename std::map<Node<dim>*, std::set<Element<dim>*> >::iterator nec = node_element_connector.begin();
                  nec != node_element_connector.end(); nec++)
    {
        if ((nec->first->AtBoundary() == NOT) && (nec->second.size() > 2U))
        {
            ControlVolumeElement<dim>& PCVE(CVE_Collection_[node_CVE_map[nec->first]]);
            for (typename std::set<Element<dim>*>::iterator eit = nec->second.begin(); eit != nec->second.end(); eit++)
            {
                ControlVolumeElement<dim>& LCVE(CVE_Collection_[element_CVE_map[*eit]]);
                faceCollection_.push_back(csmp::ElementFace<dim> (LCVE, PCVE, face_local_storage_size, NOT));
            }
        }
    }
    
    interiorPointFaces_ = faceCollection_.size() - interiorLineFaces_;
    
    
    // creating boundary line faces
    for (typename std::map<size_t, std::map<size_t, std::set<Element<dim>*> > >::iterator min_n_it = shared_faces_map.begin();
         min_n_it != shared_faces_map.end(); min_n_it++)
    {
        for (typename std::map<size_t, std::set<Element<dim>*> >::iterator max_n_it = min_n_it->second.begin();
             max_n_it != min_n_it->second.end(); max_n_it++)
        {
            size_t left_boundary_nodes(0U);
            size_t right_boundary_nodes(0U);
            size_t top_boundary_nodes(0U);
            size_t bottom_boundary_nodes(0U);
            
            if (isLEFT(id_node_map[min_n_it->first]->AtBoundary())) left_boundary_nodes++;
            if (isRIGHT(id_node_map[min_n_it->first]->AtBoundary())) right_boundary_nodes++;
            if (isTOP(id_node_map[min_n_it->first]->AtBoundary())) top_boundary_nodes++;
            if (isBOTTOM(id_node_map[min_n_it->first]->AtBoundary())) bottom_boundary_nodes++;

            if (isLEFT(id_node_map[max_n_it->first]->AtBoundary())) left_boundary_nodes++;
            if (isRIGHT(id_node_map[max_n_it->first]->AtBoundary())) right_boundary_nodes++;
            if (isTOP(id_node_map[max_n_it->first]->AtBoundary())) top_boundary_nodes++;
            if (isBOTTOM(id_node_map[max_n_it->first]->AtBoundary())) bottom_boundary_nodes++;

            if ((left_boundary_nodes == 2U) || (right_boundary_nodes == 2U) ||
                (top_boundary_nodes == 2U) || (bottom_boundary_nodes == 2U))
            {
                BOX_BOUNDARY place{ NOT };
                
                if (isLEFT(id_node_map[min_n_it->first]->AtBoundary()) && isLEFT(id_node_map[max_n_it->first]->AtBoundary())) 
                    place = LEFT;
                if (isRIGHT(id_node_map[min_n_it->first]->AtBoundary()) && isRIGHT(id_node_map[max_n_it->first]->AtBoundary())) 
                    place = RIGHT;
                if (isTOP(id_node_map[min_n_it->first]->AtBoundary()) && isTOP(id_node_map[max_n_it->first]->AtBoundary())) 
                    place = TOP;
                if (isBOTTOM(id_node_map[min_n_it->first]->AtBoundary()) && isBOTTOM(id_node_map[max_n_it->first]->AtBoundary())) 
                    place = BOTTOM;
           
                std::vector<Element<dim>*> e_ptr_2D;
                std::vector<Element<dim>*> e_ptr_1D;
                
                for ( auto eit = max_n_it->second.begin(); eit != max_n_it->second.end(); eit++)
                {
                    if ((*eit)->IsSurface())
                        e_ptr_2D.push_back((*eit));
                    else
                        e_ptr_1D.push_back((*eit));
                }
                // if there is only one 2D element
                if (e_ptr_2D.size() == 1U && e_ptr_1D.size() == 0U)
                {
                    ControlVolumeElement<dim>& SCVE(CVE_Collection_[element_CVE_map[e_ptr_2D[0U]]]);
                    faceCollection_.push_back(csmp::ElementFace<dim> (SCVE, SCVE, boundary_face_local_storage_size, place));
                    
                    continue;
                
                }
                // if there are a line and surface element at the boundary the line element is a boundary element
                else if (e_ptr_2D.size() == 1U && e_ptr_1D.size() == 1U)
                {
                    ControlVolumeElement<dim>& LCVE(CVE_Collection_[element_CVE_map[e_ptr_1D[0U]]]);
                    faceCollection_.push_back(csmp::ElementFace<dim> (LCVE, LCVE, boundary_face_local_storage_size, place));
                    continue;
                }
                // too many faces (maybe because of duplicated elements)
                else
                {
                    std::cout << "\nTransportModel::CreateControlVolumeElementsAndFaces()\n";
                    std::cout << "Too many faces were found while creating faces at the boundary of transport model.\n";
                    std::cout << "Number of 2D Element: " << e_ptr_2D.size() << "\n";
                    std::cout << "Number of 1D Element: " << e_ptr_1D.size() << "\n\n";
                    for (size_t i = 0U; i < e_ptr_2D.size(); i++)
                    {
                        std::cout << "2D Element ID " << e_ptr_2D[i]->Idx() << " Nodes\n";
                        for (uint32_t n = 0U; n < e_ptr_2D[i]->Nodes(); n++)
                        {
                            std::cout << "\tID " << e_ptr_2D[i]->N(n)->Idx() << 
                                "\tx = " << e_ptr_2D[i]->N(n)->Coordinate()[0U] <<
                                "\ty = " << e_ptr_2D[i]->N(n)->Coordinate()[1U] << "\n";
                        }
                    }
                    for (size_t i = 0U; i < e_ptr_1D.size(); i++)
                    {
                        std::cout << "1D Element ID " << e_ptr_1D[i]->Idx() << " Nodes\n";
                        for (uint32_t n = 0U; n < e_ptr_1D[i]->Nodes(); n++)
                        {
                            std::cout << "\tID " << e_ptr_1D[i]->N(n)->Idx() << 
                                "\tx = " << e_ptr_1D[i]->N(n)->Coordinate()[0U] <<
                                "\ty = " << e_ptr_1D[i]->N(n)->Coordinate()[1U] << "\n";
                        }
                    }

                    std::cout << "\n";

                    throw csmp::Exception( FATAL_ERROR, "TransportModel::CreateControlVolumeElementsAndFaces()",
                                          "Too many faces were found while creating faces at the boundary of transport model" );
                }
            }
        }
    }

    // creating line faces for the line elements at the boundaries. line elements which are connected to 
    // the boundary by one node are boundary cves and they have a boundary face
    for (typename std::map<Node<dim>*, std::set<Element<dim>*> >::iterator nec = node_element_connector.begin();
                  nec != node_element_connector.end(); nec++)
    {
        if (nec->first->AtBoundary() != NOT)
        {
            for (typename std::set<Element<dim>*>::iterator eit = nec->second.begin(); eit != nec->second.end(); eit++)
            {
                ControlVolumeElement<dim>& LCVE(CVE_Collection_[element_CVE_map[*(eit)]]);
                faceCollection_.push_back(csmp::ElementFace<dim> (LCVE, LCVE, boundary_face_local_storage_size, nec->first->AtBoundary()));
            }
        }
    }
    
    boundaryLineFaces_ = faceCollection_.size() - interiorLineFaces_ - interiorPointFaces_;
    
    // there are no boundary point faces
    boundaryPointFaces_ = 0U;
    
    // at this stage faces know their parent elements, but elements don't know their faces
    // we make an inventory of element faces and assign the faces to the elements
    std::map<ControlVolumeElement<dim>*, std::set<ElementFace<dim>* > > CVE_face_map;
    
    for (typename std::vector<ElementFace<dim> >::iterator fit = faceCollection_.begin(); fit != faceCollection_.end(); fit++)
    {
        CVE_face_map[fit->InsideCVE()].insert(&(*fit));
        CVE_face_map[fit->OutsideCVE()].insert(&(*fit));
    }
    
    // checking and assigning the faces to the CVEs
    for (typename std::map<ControlVolumeElement<dim>*, std::set<ElementFace<dim>* > >::iterator cveit = CVE_face_map.begin();
         cveit != CVE_face_map.end(); cveit++)
    {
        assert (cveit->first->Faces() == cveit->second.size());
        
        std::vector<ElementFace<dim>* > f_ptr_l;
        std::vector<ElementFace<dim>* > f_ptr_p;
        
        for (typename std::set<ElementFace<dim>* >::iterator fit = cveit->second.begin(); fit != cveit->second.end(); fit++)
        {
            if (((*fit)->InsideCVE()->Dimension() == 0U) || ((*fit)->OutsideCVE()->Dimension() == 0U))
                f_ptr_p.push_back(*fit);
            else
                f_ptr_l.push_back(*fit);
        }
        
        std::vector<char> direction;
        for (size_t i = 0U; i < f_ptr_l.size(); i++)
        {
            if (f_ptr_l[i]->InsideCVE() ==  cveit->first) 
                direction.push_back(1);
            else
                direction.push_back(-1);
        }
        for (size_t i = 0U; i < f_ptr_p.size(); i++)
        {
            if (f_ptr_p[i]->InsideCVE() ==  cveit->first) 
                direction.push_back(1);
            else
                direction.push_back(-1);
        }
                
        size_t fcounter(0U);
        // putting line faces first
        for (size_t i = 0U; i < f_ptr_l.size(); i++, fcounter++)
          {
              cveit->first->Face(f_ptr_l[i], fcounter, direction[fcounter]);
          }
        // putting point faces at the end
        for (size_t i = 0U; i < f_ptr_p.size(); i++, fcounter++)
          {
              cveit->first->Face(f_ptr_p[i], fcounter, direction[fcounter]);
          }        
    }

} // end CreateControlVolumeElementsAndFaces








template<uint32_t dim>
void TransportModel<dim>::CreatePatches(const char* vol_mod)
{
    std::cout << "\nTransportModel<" << dim << ">::CreatePatches()\n";
    
    Index vol_mod_key(modelRef_.Database().StorageKey(vol_mod));
    
    interiorPatches_ = 0U;
    boundaryPatches_ = 0U;
    
    // making an inventory of a node's CVEs parent
    std::map<Node<dim>*, std::set<ControlVolumeElement<dim>* > > node_cve_map;
    
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = CVEsBegin(); cveit != CVEsEnd(); cveit++)
    {
        for (size_t n = 0U; n < cveit->Nodes(); n++)
        {
            node_cve_map[cveit->N(n)].insert(&(*cveit));
        }
        
    }
    
    patchCollection_.reserve(node_cve_map.size());
    
    for (typename std::map<Node<dim>*, std::set<ControlVolumeElement<dim>* > >::iterator ncveit = node_cve_map.begin(); 
         ncveit != node_cve_map.end(); ncveit++)
    {
        // first creating interior patches
        if (ncveit->first->AtBoundary() == NOT)
        {
            std::set<ControlVolumeElement<dim>* > cves;
            for (typename std::set<ControlVolumeElement<dim>* >::iterator cveit = ncveit->second.begin();
                 cveit != ncveit->second.end(); cveit++)
            {
                cves.insert(*cveit);
            }

            std::vector<ControlVolumeElement<dim>*> cves_vec; // stores CVEs in an order
            typename std::set<ControlVolumeElement<dim>* >::iterator cvesit;//(cves.begin());
            ControlVolumeElement<dim>* line_cve_ptr(0);
            ControlVolumeElement<dim>* point_cve_ptr(0);
            
            // finding a 1D dimensional CVE if exists
            for (cvesit = cves.begin(); cvesit != cves.end(); cvesit++)
            {
                if ((*cvesit)->Dimension() == 1U) break;
            }
            //while (((*cvesit)->Dimension() != 1U) && (cvesit != cves.end())) cvesit++;
            
            if (cvesit == cves.end())
            {
                cves_vec.push_back(*cves.begin());
                cves.erase(*cves.begin());
            }
            else
            {
                line_cve_ptr = *cvesit;
                cves.erase(*cvesit);
                
                // finding the 0D CVE if exists
                for ( cvesit = cves.begin(); cvesit != cves.end(); cvesit++)
                {
                    if ((*cvesit)->Dimension() == 0U)
                    {
                        point_cve_ptr = *cvesit;
                        cves.erase(*cvesit);
                        break;
                    }
                    
                }
                
                // finding the 2D neighbor of 1D CVE and inserting it in cves_vec as the first 2D cve
                for (size_t n = 0U; n < line_cve_ptr->Neighbors(); n++)
                {
                    if (line_cve_ptr->Neighbor(n)->Dimension() == 2U)
                    {
                        cves_vec.push_back(line_cve_ptr->Neighbor(n));
                        cves.erase(line_cve_ptr->Neighbor(n));
                        break;
                    }
                    
                }
                
            }

            size_t counter(0U);
            while (cves.size() != 0U)
            {
                for (size_t n = 0U; n < cves_vec[counter]->Neighbors(); n++)
                {
                    if (cves.find(cves_vec[counter]->Neighbor(n)) != cves.end())
                    {
                        cves_vec.push_back(cves_vec[counter]->Neighbor(n));
                        cves.erase(cves_vec[counter]->Neighbor(n));
                        counter++;
                        break;
                    }                
                }
                
            } // end while (cves is empty)


            // inserting line CVE and point CVE at the end if they exists
            if (line_cve_ptr != 0) cves_vec.push_back(line_cve_ptr);
            if (point_cve_ptr != 0) cves_vec.push_back(point_cve_ptr);
            
            // creating the patch and reserving storage for the sectors
            patchCollection_.push_back(CVE_Patch<dim> (*(ncveit->first)));
            patchCollection_[interiorPatches_].ReserveSectorsStorage(cves_vec.size());
            
            // creating sectors for the patch
            for (size_t s = 0U; s < cves_vec.size(); s++)
            {
                double width(cves_vec[s]->Dimension() == 1U ? cves_vec[s]->E()->Read(vol_mod_key) : 1.);
                patchCollection_[interiorPatches_].CreateAndAddSector(*cves_vec[s], width);
            }
            
            std::map<ElementFace<dim>*, size_t> face_number_map;
            size_t number_of_faces(0U);
            std::vector<char> cyclic_face_direction;
            
            // assigning numbers to the faces
            for (typename std::vector<CVE_Sector<dim> >::iterator sit = patchCollection_[interiorPatches_].SectorsBegin(); 
                 sit != patchCollection_[interiorPatches_].SectorsEnd(); sit++)
            {
                for (size_t f = 0U; f < sit->Faces(); f++)
                {
                    if (face_number_map.find(sit->Face(f)) == face_number_map.end())
                    {
                        sit->PatchFaceNumber(f, number_of_faces);
                        face_number_map[sit->Face(f)] = number_of_faces;
                        number_of_faces++;
                        
                        // assigning cyclic face number to the face in the patch
                        // the faces direction is in the same direction of sectors visitation in the 
                        // sectorCollection_ vector in the CVE_Patch class (if sectors are stored clockwise cyclic
                        // face direction is pointing in clockwise direction for the right face of sectors)
                        
                        // for the first sector we treat it seperately
                        if (sit == patchCollection_[interiorPatches_].SectorsBegin())
                        {
                            if (sit->Face(f)->InsideCVE() == sit->CVE())
                            {
                                if (sit->Face(f)->OutsideCVE() == patchCollection_[interiorPatches_].SectorAccessor(1U)->CVE())
                                {
                                    cyclic_face_direction.push_back(1);
                                }
                                else
                                {
                                    cyclic_face_direction.push_back(-1);
                                }
                                
                            }
                            else
                            {
                                if (sit->Face(f)->InsideCVE() == patchCollection_[interiorPatches_].SectorAccessor(1U)->CVE())
                                {
                                    cyclic_face_direction.push_back(-1);
                                }
                                else
                                {
                                    cyclic_face_direction.push_back(1);
                                }
                            
                            }
                            
                        }
                        else
                        {
                            if (sit->Face(f)->InsideCVE() == sit->CVE())
                            {
                                cyclic_face_direction.push_back(1);
                            }
                            else
                            {
                                cyclic_face_direction.push_back(-1);
                            }
                        
                        }
                        
                    }
                    else
                    {
                        sit->PatchFaceNumber(f, face_number_map[sit->Face(f)]);
                    }
                    
                }
                
            }
            
            patchCollection_[interiorPatches_].ReserveFacesStorage(number_of_faces);
            
            // creating a vector with the same sequence of sectors face number to insert in the same sequence 
            // in the patch face pointers container
            std::vector<ElementFace<dim>*> patch_faces(number_of_faces);
            for (typename std::map<ElementFace<dim>*, size_t>::iterator fit = face_number_map.begin(); fit != face_number_map.end(); fit++)
            {
                patch_faces[fit->second] = fit->first;
            }
            
            // inserting faces pointer in the patch face pointers container
            for (typename std::vector<ElementFace<dim>*>::iterator fit = patch_faces.begin(); fit != patch_faces.end(); fit++)
            {
                patchCollection_[interiorPatches_].PushbackFace(*fit);
            }

            // assigning cyclic face direction to the faces in the patch
            for (size_t f = 0U; f < number_of_faces; f++)
            {
                patchCollection_[interiorPatches_].FaceCyclicDirection(f, cyclic_face_direction[f]);
            }

            
            patchCollection_[interiorPatches_].InteriorFaces(number_of_faces);
            patchCollection_[interiorPatches_].BoundaryFaces(0);
            
            interiorPatches_++;
            
        } // end if it is not boundary node
        
    } // end for loop for node_cve_map (creating interior patches)

    
    // for loop for creating boundary patches
    for (typename std::map<Node<dim>*, std::set<ControlVolumeElement<dim>* > >::iterator ncveit = node_cve_map.begin(); 
         ncveit != node_cve_map.end(); ncveit++)
    {
        // boundary patches
        if (ncveit->first->AtBoundary() != NOT)
        {
            std::set<ControlVolumeElement<dim>* > cves;
            for (typename std::set<ControlVolumeElement<dim>* >::iterator cveit = ncveit->second.begin();
                 cveit != ncveit->second.end(); cveit++)
            {
                cves.insert(*cveit);
            }
            
            std::vector<ControlVolumeElement<dim>*> cves_vec; // stores CVEs in an order

            ControlVolumeElement<dim>* point_cve_ptr(0);
            
            // finding the 0D CVE if exists
            for (typename std::set<ControlVolumeElement<dim>* >::iterator cvesit = cves.begin(); 
                 cvesit != cves.end(); cvesit++)
            {
                if ((*cvesit)->Dimension() == 0U)
                {
                    point_cve_ptr = *cvesit;
                    cves.erase(*cvesit);
                    break;
                }
                
            }
            
            // finding a boundary CVE and inserting in the cves_vec as the first cve
            size_t left_boundary(isLEFT(ncveit->first->AtBoundary()));
            size_t right_boundary(isRIGHT(ncveit->first->AtBoundary()));
            size_t top_boundary(isTOP(ncveit->first->AtBoundary()));
            size_t bottom_boundary(isBOTTOM(ncveit->first->AtBoundary()));
            
            for (typename std::set<ControlVolumeElement<dim>* >::iterator cvesit = cves.begin(); cvesit != cves.end(); cvesit++)
            {
                // lower dimensional element can not be a boundary element
                if ((*cvesit)->Dimension() == 2U)
                {
                    for (size_t f = 0U; f < (*cvesit)->Faces(); f++)
                    {
                        if (isLEFT((*cvesit)->Face(f)->Placement())) left_boundary++;
                        if (isRIGHT((*cvesit)->Face(f)->Placement())) right_boundary++;
                        if (isTOP((*cvesit)->Face(f)->Placement())) top_boundary++;
                        if (isBOTTOM((*cvesit)->Face(f)->Placement())) bottom_boundary++;
                        
                        if ((left_boundary == 2U) || (right_boundary == 2U) || (top_boundary == 2U) || (bottom_boundary == 2U))
                        {
                            cves_vec.push_back(*cvesit);
                            cves.erase(*cvesit);
                            break;
                        }
                        
                    }
                
                }
                
                if (cves_vec.size() == 1U) break;
                
            }
            
            if (cves_vec.size() == 0U)
            {
                throw csmp::Exception( FATAL_ERROR, "TransportModel::CreatePatches()",
                                      "Could not find boundary CVE. Check the model maybe you have a lower dimensional element at the boundary!" );
            }
            
            
            // inserting the rest of sectors in cves_vec
            size_t counter(0U);
            while (cves.size() != 0U)
            {
                for (size_t n = 0U; n < cves_vec[counter]->Neighbors(); n++)
                {
                    if (cves.find(cves_vec[counter]->Neighbor(n)) != cves.end())
                    {
                        cves_vec.push_back(cves_vec[counter]->Neighbor(n));
                        cves.erase(cves_vec[counter]->Neighbor(n));
                        counter++;
                        break;
                    }                
                }
                
            } // end while (cves is empty)

            // inserting point CVE at the end if it exists
            if (point_cve_ptr != 0) cves_vec.push_back(point_cve_ptr);
            
            // creating the patch and reserving storage for the sectors
            patchCollection_.push_back(CVE_Patch<dim> (*(ncveit->first)));
            patchCollection_[interiorPatches_ + boundaryPatches_].ReserveSectorsStorage(cves_vec.size());
            
            // creating sectors for the patch
            for (size_t s = 0U; s < cves_vec.size(); s++)
            {
                double width(cves_vec[s]->Dimension() == 1U ? cves_vec[s]->E()->Read(vol_mod_key) : 1.);
                patchCollection_[interiorPatches_ + boundaryPatches_].CreateAndAddSector(*cves_vec[s], width);
            }
            
            std::map<ElementFace<dim>*, size_t> face_number_map;
            size_t number_of_faces(0U);
            std::vector<char> cyclic_face_direction;
           
            // assigning numbers to the interior faces
            for (typename std::vector<CVE_Sector<dim> >::iterator sit = patchCollection_[interiorPatches_ + boundaryPatches_].SectorsBegin(); 
                 sit != patchCollection_[interiorPatches_ + boundaryPatches_].SectorsEnd(); sit++)
            {
                for (size_t f = 0U; f < sit->Faces(); f++)
                {
                    if (sit->Face(f)->Placement() == NOT)
                    {
                        if (face_number_map.find(sit->Face(f)) == face_number_map.end())
                        {
                            sit->PatchFaceNumber(f, number_of_faces);
                            face_number_map[sit->Face(f)] = number_of_faces;
                            number_of_faces++;

                            // assigning cyclic face number to the face in the patch
                            // the faces direction is in the same direction of sectors visitation in the 
                            // sectorCollection_ vector in the CVE_Patch class (if sectors are stored clockwise cyclic
                            // face direction is pointing in clockwise direction for the right face of sectors)                            
                            if (sit->Face(f)->InsideCVE() == sit->CVE())
                            {
                                cyclic_face_direction.push_back(1);
                            }
                            else
                            {
                                cyclic_face_direction.push_back(-1);
                            }

                        }
                        else
                        {
                            sit->PatchFaceNumber(f, face_number_map[sit->Face(f)]);
                        }
                    
                    }
                    
                }
                
            }
            
            size_t interior_faces(number_of_faces);
            patchCollection_[interiorPatches_ + boundaryPatches_].InteriorFaces(interior_faces);
            
            // assigning numbers to the boundary faces
            for (typename std::vector<CVE_Sector<dim> >::iterator sit = patchCollection_[interiorPatches_ + boundaryPatches_].SectorsBegin(); 
                 sit != patchCollection_[interiorPatches_ + boundaryPatches_].SectorsEnd(); sit++)
            {
                for (size_t f = 0U; f < sit->Faces(); f++)
                {
                    if (sit->Face(f)->Placement() != NOT)
                    {
                        if (face_number_map.find(sit->Face(f)) == face_number_map.end())
                        {
                            sit->PatchFaceNumber(f, number_of_faces);
                            face_number_map[sit->Face(f)] = number_of_faces;
                            number_of_faces++;
                            
                            // assigning cyclic face number to the face in the patch
                            
                            // for the first boundary face which belongs to the first sector in the sectorCollection_ vector
                            // the cyclic face direction normal belongs to the imaginary element in the boundary
                            if (number_of_faces == interior_faces + 1U)
                            {
                                cyclic_face_direction.push_back(-1);
                            }
                            else
                            {
                                cyclic_face_direction.push_back(1);
                            }
                            
                        }
                        else
                        {
                            sit->PatchFaceNumber(f, face_number_map[sit->Face(f)]);
                        }
                    
                    }
                    
                }
                
            }
            
            patchCollection_[interiorPatches_ + boundaryPatches_].BoundaryFaces(number_of_faces - interior_faces);

            patchCollection_[interiorPatches_ + boundaryPatches_].ReserveFacesStorage(number_of_faces);

            // creating a vector with the same sequence of sectors face number to insert in the same sequence 
            // in the patch face pointers container
            std::vector<ElementFace<dim>*> patch_faces(number_of_faces);
            for (typename std::map<ElementFace<dim>*, size_t>::iterator fit = face_number_map.begin(); fit != face_number_map.end(); fit++)
            {
                patch_faces[fit->second] = fit->first;
            }
            
            // inserting faces pointer in the patch face pointers container
            for (typename std::vector<ElementFace<dim>*>::iterator fit = patch_faces.begin(); fit != patch_faces.end(); fit++)
            {
                patchCollection_[interiorPatches_ + boundaryPatches_].PushbackFace(*fit);
            }

            // assigning cyclic face direction to the faces in the patch
            for (size_t f = 0U; f < number_of_faces; f++)
            {
                patchCollection_[interiorPatches_ + boundaryPatches_].FaceCyclicDirection(f, cyclic_face_direction[f]);
            }
            
            boundaryPatches_++;
        
        }
        
    } // end for loop for node_cve_map (creating boundary patches)
    
} // end CreateControlVolumeElementsAndFaces







template<uint32_t dim>
void TransportModel<dim>::CorrectVolumes(const char* vol_mod)
{
    std::cout << "\nTransportModel<" << dim << ">::CorrectVolumes()\n";
    
    Index vol_mod_key(modelRef_.Database().StorageKey(vol_mod));
    
    double min_volume(std::numeric_limits<double>::max()); // this is used for assigning a volume to 0D cves
    
    // modifying interior surface and line CVE's volume
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = SurfaceCVEsBegin(); 
         cveit != PointCVEsBegin(); cveit++)
    {
        const double new_volume(cveit->E()->Read(vol_mod_key) * cveit->E()->Volume());
        cveit->Volume(new_volume);
        
        if (new_volume < min_volume) min_volume = new_volume;
    }
    // modifying boundary surface and line CVE's volume
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = PerimeterSurfaceCVEsBegin(); 
         cveit != PerimeterPointCVEsBegin(); cveit++)
    {
        const double new_volume(cveit->E()->Read(vol_mod_key) * cveit->E()->Volume());
        cveit->Volume(new_volume);
        
        if (new_volume < min_volume) min_volume = new_volume;
    }

    // looping over all faces and calculate face area for lower dimensional elements
    for (typename std::vector<ElementFace<dim> >::iterator fit = FacesBegin(); fit != FacesEnd(); fit++)
    {
        if (fit->InsideCVE()->Dimension() == 1U  && fit->OutsideCVE()->Dimension() == 1U)
        {
            double area(0.5 * (fit->InsideCVE()->E()->Read(vol_mod_key) + fit->OutsideCVE()->E()->Read(vol_mod_key)));
            fit->Area(area);
        }
        
        if (fit->InsideCVE()->Dimension() == 1U  && fit->OutsideCVE()->Dimension() == 0U)
        {
            fit->Area(fit->InsideCVE()->E()->Read(vol_mod_key));
        }
        
        if (fit->InsideCVE()->Dimension() == 0U  && fit->OutsideCVE()->Dimension() == 1U)
        {
            fit->Area(fit->OutsideCVE()->E()->Read(vol_mod_key));
        }
        
    }

    // modifying boundary point CVE boundary face area
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = PerimeterPointCVEsBegin(); 
         cveit != PointCVEsEnd(); cveit++)
    {
        double area(0.);
        size_t counter(0U);
        
        for (size_t f = 0U; f < cveit->Faces(); f++)
        {
            if (cveit->Face(f)->Placement() == NOT) 
            {
                area += cveit->Face(f)->Area();
                counter++;
            }
        }
        
        area /= (double) counter;

        for (size_t f = 0U; f < cveit->Faces(); f++)
        {
            if (cveit->Face(f)->Placement() != NOT) 
            {
                cveit->Face(f)->Area(area);
            }
        }
    }
    
    // assigning the lowest volume in model to the 0D cves
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = PointCVEsBegin(); 
         cveit != PerimeterSurfaceCVEsBegin(); cveit++)
    {
        cveit->Volume(min_volume);
    }

    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = PerimeterPointCVEsBegin(); 
         cveit != PointCVEsEnd(); cveit++)
    {
        cveit->Volume(min_volume);
    }
        
} // end CorrectVolumes





template<uint32_t dim>
void TransportModel<dim>::AssignParentPatchesToCVEs()
{
    std::cout << "\nTransportModel<" << dim << ">::AssignParentPatchesToCVEs()\n";

    typename std::map<ControlVolumeElement<dim>*, typename std::set<CVE_Patch<dim>*> > cve_patch_map;

    for (typename std::vector<CVE_Patch<dim> >::iterator pit = PatchesBegin(); pit != PatchesEnd(); ++pit)
    {
        for (typename std::vector<CVE_Sector<dim> >::iterator sit = pit->SectorsBegin(); sit != pit->SectorsEnd(); ++sit)
        {
            cve_patch_map[sit->CVE()].insert(&(*pit));
        }

    }

    // assigning the patches to the cves
    for (typename std::map<ControlVolumeElement<dim>*, typename std::set<CVE_Patch<dim>*> >::iterator
         cve_p_it = cve_patch_map.begin(); cve_p_it != cve_patch_map.end(); ++cve_p_it)
    {
        cve_p_it->first->ReservePatchStorage(cve_p_it->second.size());

        for (typename std::set<CVE_Patch<dim>*>::iterator pit = cve_p_it->second.begin(); 
             pit != cve_p_it->second.end(); ++pit)
        {
            cve_p_it->first->PushbackPatch(*pit);
        }

    }

} // end AssignParentPatchesToCVEs




template<uint32_t dim>
void TransportModel<dim>::TransportModelReport()
{
    std::cout << "\n*************************************\n";
    std::cout << "       Transport Model Report\n";
    std::cout << "*************************************\n\n";
    std::cout << "Interior 2D CVEs:       " << interiorSurfaceCVE_ << "\n";
    std::cout << "Interior 1D CVEs:       " << interiorLineCVE_ << "\n";
    std::cout << "Interior 0D CVEs:       " << interiorPointCVE_ << "\n";
    std::cout << "Boundary 2D CVEs:       " << boundarySurfaceCVE_ << "\n";
    std::cout << "Boundary 1D CVEs:       " << boundaryLineCVE_ << "\n";
    std::cout << "Boundary 0D CVEs:       " << boundaryPointCVE_ << "\n";
    std::cout << "-------------------------------------\n";
    std::cout << "Total CVEs:             " << CVEs() << "\n";
    std::cout << "-------------------------------------\n\n";

    std::cout << "Interior line faces:    " << interiorLineFaces_ << "\n";
    std::cout << "Interior point faces:   " << interiorPointFaces_ << "\n";
    std::cout << "Boundary line faces:    " << boundaryLineFaces_ << "\n";
    std::cout << "Boundary point faces:   " << boundaryPointFaces_ << "\n";
    std::cout << "-------------------------------------\n";
    std::cout << "Total faces:            " << Faces() << "\n";
    std::cout << "-------------------------------------\n\n";

    std::cout << "Interior patches:       " << interiorPatches_ << "\n";
    std::cout << "Boundary patches:       " << boundaryPatches_ << "\n";
    std::cout << "-------------------------------------\n";
    std::cout << "Total patches:          " << Patches() << "\n";
    std::cout << "-------------------------------------\n";
    
}




template class TransportModel<2U>;


} // end csmp
