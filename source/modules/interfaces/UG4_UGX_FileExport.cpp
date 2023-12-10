//
//  UG4_UGX_FileExport.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 3/12/2023.
//

#include "UG4_UGX_FileExport.h"
#include "XML_Document.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyDatabase.h"
#include "ErrorHandler.h"


using namespace std;

namespace csmp {

/**
      Helper function that spits out some reasonable colours for the visualisation of subsets inside of PROMESH
*/
tuple<float,float,float,float,std::string>  goodColour( int id )
{
   switch (id) { //       R    G    B    A(opaque)   colour name
        case  0: return { 1.0, 1.0, 1.0, 1.0, "white"};
        
        case  1: return { 0.0, 0.2, 1.0, 1.0, "blue"};
        case  2: return { 0.0, 0.8, 0.0, 1.0, "green"};

        case  3: return { 0.0, 1.0, 1.0, 1.0, "magenta"};
        case  4: return { 0.5, 0.0, 0.0, 1.0, "purple"};
        case  5: return { 0.5, 0.5, 0.0, 1.0, "khaki"};
        case  6: return { 0.6, 0.4, 0.2, 1.0, "brown"};
        case  7: return { 1.0, 0.0, 0.0, 1.0, "red"};
        case  8: return { 0.8, 0.5, 0.0, 1.0, "orange"};
        case  9: return { 1.0, 1.0, 0.0, 1.0, "yellow"};
        
        case 10: return { 0.0, 0.0, 0.0, 1.0, "black"};
     }
   throw std::invalid_argument("goodColours(id)");
}


template<uint32_t dim>
UG4_UGX_FileExport<dim>::UG4_UGX_FileExport( const Model<dim>& model )
 : trias_(0u), quads_(0u), tets_(0u), hexes_(0u), prisms_(0u), pyras_(0u)
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // if the model has P refinement or the order of the shape-function is greater than 1 (curvilinear elements) it cannot be output to UG4
    set<uint32_t> shape_function_order = model.Mesh().OrderOfShapeFunctions();
    if ( shape_function_order.size() > 1 || (*shape_function_order.begin()) != 1U ) {
         csmp_err.Note( ERROR, "UG4_UGX_FileExport (custom constructor)",
                       "UG4 - ugx format only supports output of straight-sided elements (shape-function order = 1; nothing was output");
         return;
      }
    
    // number elements globally in consecutive order
    model.Region("Model").UpdateMemberIndexes();
    // NumberElementsForUGX( model );

    cout<<"\nUG4_UGX_FileExport(ctor): collecting node numbers of cell edges into sorted 'edges_' vector...";
    size_t n_edges = CollectEdges( model, edges_ );
    assert( n_edges > model.Mesh().Elements() );
    cout <<" "<< n_edges <<" unique edges found"<<endl;
    
    cout<<"\nUG4_UGX_FileExport(ctor): collecting node numbers of cell faces into sorted 'faces_' vector...";
    size_t n_faces = CollectFaces( model, faces_ );
    assert( n_faces > 0U );
    cout <<" "<< n_faces <<" unique element faces found"<<endl;

    if constexpr( dim == 3U ) {
        cout<<"\nUG4_UGX_FileExport(ctor): collecting node numbers of cells into sorted 'volumes_' vector...";
        size_t n_volumes = CollectVolumes( model, volumes_ );
        assert( n_volumes > 0U && n_volumes <= model.Mesh().Elements() );
        cout <<" "<< n_volumes <<" unique element IDs found"<<endl;
      }
 
 } // end constructor
  
  
 
 
/** Naive version where the XML tags are handcrafted
 
  filename without extension which will be appended
  
  @attention method assumes that the nodes and elements have already been numbered from 0..n-1 following UG's conventions
  
*/
template<uint32_t dim>
bool UG4_UGX_FileExport<dim>::Write_UGX_FileASCII( const Model<dim>& model, const std::string& file_name ) const
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    ofstream ofs( (file_name + ".ugx").c_str() );
    
    // 1. File header
    // --------------
    ofs <<"<?xml version=\"1.0\" encoding=\"utf-8\"?>"<< endl;
    //    using the model name as grid name
    ofs <<"<grid name=\""<< model.Name() <<"\">"<< endl;
    
    // 2. Vertex coordinates (either groups of 2(2D) or 3(3D)
    // ------------------------------------------------------
    ofs <<"\t<vertices "<<"coords=\""<< dim <<"\">";
    const Region<dim>& model_domain(model.Region("Model"));
    bool print_whitespace{false};
    // remembering ofs state before setting it to scientific
    ios_base::fmtflags stream_state( ofs.flags() );
    // print in scientific format honpuring given precision
    ofs << scientific;
    for ( const auto& nit : model_domain.NodeVector() )
      for ( uint32_t i{0U}; i<dim; ++i ) {
           if ( print_whitespace ) ofs <<" ";
           else print_whitespace = true;
           ofs << (*nit)[i];
        }
    ofs.flags( stream_state );
    ofs <<"</vertices>"<< endl;
    
    // 3. Edges (ignoring line elements)
    // ---------------------------------
    ofs <<"\t<edges>";
    // printing the edge nodes (one pair per edge)
    print_whitespace = false;
    for ( const auto& it : edges_ ) {
          if ( print_whitespace ) ofs <<" ";
          else print_whitespace = true;
          ofs << it.first <<" "<< it.second;
      }
    ofs <<"</edges>"<< endl;
    
    // 4. Faces (triangles first and then quadrilaterals)
    // ---------------------------------------------------
    ofs <<"\t<triangles>";
    // printing the nodes
    if ( HasTriangles() ) {
        print_whitespace = false;
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsSurface() && isTriangular( it->FE_Type() ) ) {
               if ( print_whitespace ) ofs <<" ";
               else print_whitespace = true;
               ofs <<      it->N(0)->Idx();
               ofs <<" "<< it->N(1)->Idx();
               ofs <<" "<< it->N(2)->Idx();
            }
        ofs <<"</triangles>"<< endl;
      }
    // quads
    if ( HasQuadrilaterals() ) {
        print_whitespace = false;
        ofs <<"\t<quadrilaterals>";
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) ) {
               if ( print_whitespace ) ofs <<" ";
               else print_whitespace = true;
               ofs <<      it->N(0)->Idx();
               ofs <<" "<< it->N(1)->Idx();
               ofs <<" "<< it->N(2)->Idx();
               ofs <<" "<< it->N(3)->Idx();
            }
        ofs <<"</quadrilaterals>"<< endl;
      }
   
    // 5. Volumes (tet, hex, prism, pyra, octa)
    // ----------------------------------------
    if constexpr( dim == 3U ) {
        if ( HasTetrahedra() ) {
            print_whitespace = false;
            ofs <<"\t<tetrahedrons>";
            // printing the nodes (midside nodes last)
            for ( const auto& it : model_domain.CellVector() )
              if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
                   if ( print_whitespace ) ofs <<" ";
                   else print_whitespace = true;
                   ofs << it->N(0)->Idx();
                   for ( uint32_t i{1U}; i<it->Nodes(); ++i )
                     ofs <<" "<< it->N(i)->Idx();
                }
             ofs <<"</tetrahedrons>"<< endl;
          }
        // hexa
        if ( HasHexahedra() ) {
            print_whitespace = false;
            ofs <<"\t<hexahedrons>";
            for ( const auto& it : model_domain.CellVector() )
              if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
                   if ( print_whitespace ) ofs <<" ";
                   else print_whitespace = true;
                   ofs << it->N(0)->Idx();
                   for ( uint32_t i{1U}; i<it->Nodes(); ++i )
                     ofs <<" "<< it->N(i)->Idx();
                }
            ofs <<"</hexahedrons>"<< endl;
          }
        // prisms
        if ( HasPrisms() ) {
            print_whitespace = false;
            ofs <<"\t<prisms>";
            for ( const auto& it : model_domain.CellVector() )
              if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
                   if ( print_whitespace ) ofs <<" ";
                   else print_whitespace = true;
                   ofs << it->N(0)->Idx();
                   for ( uint32_t i{1U}; i<it->Nodes(); ++i )
                     ofs <<" "<< it->N(i)->Idx();
                }
            ofs <<"</prisms>"<< endl;
          }
        // pyramids
        if ( HasPyramids() ) {
            print_whitespace = false;
            ofs <<"\t<pyramids>";
            for ( const auto& it : model_domain.CellVector() )
              if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
                   if ( print_whitespace ) ofs <<" ";
                   else print_whitespace = true;
                   ofs << it->N(0)->Idx();
                   for ( uint32_t i{1U}; i<it->Nodes(); ++i )
                     ofs <<" "<< it->N(i)->Idx();
                }
            ofs <<"</pyramids>"<< endl;
          }
      }
    
    // 6. Property value attachments (cell properties only)
    // ----------------------------------------------------
    // (passOn="1" means that variable will be mapped to finer grids)
    // (global="1" means that the variable name will be known across UG)
    // the variable values must be in exactly in the same order as the elements inside of the 'faces' and 'volumes' arrays
    set<std::string>    output_props;
    VectorVariable<dim> vc;
    TensorVariable<dim> ts;
    size_t n_variables = model.Database().ListVariables( ELEMENT, output_props );
    assert( n_variables >= 1 );
    
    // 6.1 face attachments (first triangles, then quadrilaterals)
    // -----------------------------------------------------------
    vector<size_t>  cell_idx_in_attachment_order;
    CellIdxInFaceAttachmentsOrder( model, cell_idx_in_attachment_order );
    
    // where there is no property value, 'NaN' (not a number) will be written to file
    const string no_data_value{"NaN"};
    
    for ( const auto& pit : output_props ) {
         const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
         // skipping variable types that are not supported by UG
         if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
           continue;
         ofs <<"\t<face_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key) <<"\" ";
         ofs <<"passOn=\"1\" global=\"1\">";

         print_whitespace = false;
         for ( const auto& it : cell_idx_in_attachment_order ) {
              // there no cell index was recorded for the face the no data values is writen
              switch(var_key.type) {
                  case SCALAR: if ( !print_whitespace ) ofs << scientific << model_domain.E(it)->Read(var_key);
                               else ofs <<" "<< scientific << model_domain.E(it)->Read(var_key);
                               print_whitespace = true;
                    break;
                  case VECTOR:
                       model_domain.E(it)->Read(var_key,vc);
                       for ( uint32_t i{0u}; i<dim; ++i ) {
                             if ( print_whitespace ) ofs <<" ";
                             else print_whitespace = true;
                             ofs << scientific << vc[i];
                         }
                    break;
                  case TENSOR:
                       model_domain.E(it)->Read(var_key,ts);
                       for ( uint32_t i{0u}; i<dim; ++i )
                         for ( uint32_t j{0u}; j<dim; ++j ) {
                               if ( print_whitespace ) ofs <<" ";
                               else print_whitespace = true;
                               ofs << scientific << ts(i,j);
                           }
                    break;
                  default:
                    csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::Write_UGX_FileASCII", pit,
                                  "type of variable cannot be represented in UG4");
                }
            }

         ofs <<" </face_attachment>"<< endl;
         
      } // end for output properties


    // 6.2 volume attachments
    // ----------------------
    if constexpr ( dim == 3U )
      {
         CellIdxInVolumeAttachmentsOrder( model, cell_idx_in_attachment_order );

         for ( const auto& pit : output_props ) {
              const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
              // skipping variable types that are not supported by UG
              if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
                continue;
              ofs <<"\t<volume_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key);
              ofs <<"passOn=\"1\" global=\"1\">";

              print_whitespace = false;
              for ( const auto& it : cell_idx_in_attachment_order ) {
                    switch(var_key.type) {
                        case SCALAR: if ( !print_whitespace ) ofs << scientific << model_domain.E(it)->Read(var_key);
                                     else ofs <<" "<< scientific << model_domain.E(it)->Read(var_key);
                                     print_whitespace = true;
                          break;
                        case VECTOR:
                             model_domain.E(it)->Read(var_key,vc);
                             for ( uint32_t i{0u}; i<dim; ++i ) {
                                   if ( print_whitespace ) ofs <<" ";
                                   else print_whitespace = true;
                                   ofs << scientific << vc[i];
                               }
                          break;
                        case TENSOR:
                             model_domain.E(it)->Read(var_key,ts);
                             for ( uint32_t i{0u}; i<dim; ++i )
                               for ( uint32_t j{0u}; j<dim; ++j ) {
                                     if ( print_whitespace ) ofs <<" ";
                                     else print_whitespace = true;
                                     ofs << scientific << ts(i,j);
                                 }
                        break;
                            default:
                              csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::Write_UGX_FileASCII", pit,
                                            "type of variable cannot be represented in UG4");
                          }
                     }
              ofs <<" </volume_attachment>"<< endl;
           }
      }
    
    
    // 7. Subset (Region) Handler and subsets
    // --------------------------------------
    // (NB: subset do not have to include edges, faces, and volumes, but can only exist of vertices)
    // (NB: lower-dimensional element regions are communicated to UG via face-ony subsets)
    // (NB: each egde or node must only be present in one subset!!!)
    int colour{0}; // ten choices that will be used alternatingly, calling goodColour( int id )
    ofs <<"\t<subset_handler name=\"defSH\">"<< endl; // defSH means default subset handler (only one is needed)
    // model geometry-defining region by region
    for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit )
      {
         // resetting the colour that the subset will be visualised with in PROMESH
         if ( colour == 10 ) colour = 9;
         const size_t domain_index{ (*rit).second.DomainIndex() };
         const auto   rgba = goodColour( colour );
         ofs <<"\t\t<subset name=\""<< (*rit).first <<"\" color=\"";
         ofs << get<0>(rgba) <<" "<< get<1>(rgba) <<" "<< get<2>(rgba) <<" "<< get<3>(rgba);
         ofs <<"\" state=\""<< domain_index <<"\">"<< endl;
         // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
         ofs <<"\t\t\t<vertices>";
         print_whitespace = false;
         for ( const auto& nit : (*rit).second.NodeVector() ) {
              if ( print_whitespace ) ofs <<" ";
              else print_whitespace = true;
              ofs << nit->Idx();
           }
         ofs <<"</vertices>"<< endl;
         
        vector<size_t> ugx_id_numbers;

        // edges
         if ( with_region_edge_output_ ) {
             print_whitespace = false;
             // compiling the edges that form part of the region
             size_t n_edges = CollectEdgesInRegion( (*rit).second, ugx_id_numbers );
             if ( n_edges > 0 ) {
                 ofs <<"\t\t\t<edges>";
                 // printing only the IDs of global edges if they belong to this unique region
                 for ( const auto& it : ugx_id_numbers ) {
                      if ( print_whitespace ) ofs <<" ";
                      else print_whitespace = true;
                      ofs << it;
                   }
                 ofs <<"</edges>"<< endl;
               }
           }
         
         // faces
         print_whitespace = false;
         if ( with_region_face_output_ ) {
             size_t n_faces = CollectFacesInRegion( (*rit).second, ugx_id_numbers );
             if ( n_faces > 0 ) {
                 ofs <<"\t\t\t<faces>";
                 // printing only the IDs of global edges if they belong to this unique region
                 for ( const auto& it : ugx_id_numbers ) {
                      if ( print_whitespace ) ofs <<" ";
                      else print_whitespace = true;
                      ofs << it;
                   }
                 ofs <<"</faces>"<< endl;
               }
           }
         
         // volumes (3d only)
         if constexpr( dim == 3U ) {
             print_whitespace = false;
             if ( with_region_volume_output_ ) {
                 size_t n_volumes = CollectVolumesInRegion( (*rit).second, ugx_id_numbers );
                 if ( n_volumes > 0 ) {
                     ofs <<"\t\t\t<volumes>";
                     // printing only the IDs of global edges if they belong to this unique region
                     for ( const auto& it : ugx_id_numbers ) {
                          if ( print_whitespace ) ofs <<" ";
                          else print_whitespace = true;
                          ofs << it;
                       }
                     ofs <<" </volumes>"<< endl;
                   }
               }
           }
         ofs <<"\t\t</subset>"<< endl;
         colour++;
     }
      
     // closing the subset section
    ofs <<"\t</subset_handler>"<< endl;

    // considering the boundaries as subsets as well (this could be a separate subset handler)
    // TODO: loop over the outer model boundaries
    
    
    // closing the grid specification section
    ofs <<"</grid>"<< endl;
 
    cout <<"\n"<<"UG4_UGX_FileExport::Write_UGX_FileASCII: output file '"<< file_name <<".ugx' written successfully."<< endl;
    return true;
 
 } // end Write_UGX_FileASCII

 
  
  
  
  
  
/**
    Loops over model, making a vector of  the unique edges (node-iD pairs) in the model (corner nodes only).
   
    The elimination of duplicates and the sorting of the edges vector guarantees that the edges are always in the same order.
    
    @param model is used for data collection but not modified.
    @param edges uses  global node IDs for sorting.
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectEdges( const Model<dim>& model, vector<pair<size_t,size_t>>& edges )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // making sure that the vector is empty and resizing it
    if ( !edges.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectEdges", "supplied edge vector not empty; deleting previous content" );
         edges.clear();
         edges.reserve( model.Mesh().Elements() * 4U ); // approximate guess of the number of edges
      }
    
    // going over the unique regions in alphabetical order, tagging the edges found within them with region ids
    const Region<dim>& model_domain(model.Region("Model"));
    vector<uint32_t> snids;
    // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
    for ( const auto& it : model_domain.CellVector() ) {
      for ( uint32_t segm_id{0U}; segm_id<it->FE()->Segments(); ++segm_id ) {
             it->FE()->NodesOfSegment( segm_id, snids );
             assert( snids.size() == 2U );
             // getting the global node IDs
             pair<size_t,size_t> global_node_numbers{ it->N(snids[0])->Idx(), it->N(snids[1])->Idx() };
             // sorting the global node IDS in ascending order
             if ( global_node_numbers.first > global_node_numbers.second )
               swap(global_node_numbers.first,global_node_numbers.second);
             // inserting the edge node ids into the edges vector
             edges.emplace_back( global_node_numbers );
          }
      }
 
    // sorting and eliminating duplicate edges from the vector
    sort( edges.begin(), edges.end() );
    edges.erase( unique( edges.begin(), edges.end() ), edges.end() );

    return edges.size();
    
 } // end CollectEdges
 
 

 
 
 
/**
     @attention order of nodes in 'faces' does not matter, only the global ID numbers are required by UG4
     @attention 'faces' related to CSMP boundaries must not be collected extra at this stage because they are already included since they match up with the element faces
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectFaces( const Model<dim>& model, vector<set<size_t>>& faces )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // making sure that the vector is empty and resizing it
    if ( !faces.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectEdges", "supplied 'faces' vector not empty; deleting previous content" );
         faces.clear();
         faces.reserve( model.Mesh().Elements() * 4U ); // approximate guess of the number of edges
      }
    
    const Region<dim>& model_domain(model.Region("Model"));

    // in 2D triangular and quadrilateral elements make up the faces
    // -------------------------------------------------------------
    if constexpr( dim == 2U ) {
        // looping over the elements of the model, collecting their global node IDs into unique sets that will be stored in 'faces'
        for ( const auto& it : model_domain.CellVector() )
          // triangular elements
          if ( it->IsSurface() && isTriangular( it->FE_Type() ) )
            {
               vector<uint32_t>  fnids;
               it->FE()->CornerNodes( fnids );
               // getting the global node IDs
               set<size_t> global_node_numbers;
               for ( const auto& id : fnids ) global_node_numbers.insert( it->N(id)->Idx() );
               faces.push_back( global_node_numbers );
               trias_++;
            }
        // quadrilaterals
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) )
            {
               vector<uint32_t>  fnids;
               it->FE()->CornerNodes( fnids );
               set<size_t> global_node_numbers;
               for ( const auto& id : fnids ) global_node_numbers.insert( it->N(id)->Idx() );
               faces.push_back( global_node_numbers );
               quads_++;
            }
      }
    // three-dimensional models (now the faces of volume elements make up the faces)
    // -----------------------------------------------------------------------------
    else {
        // looping over the elements of the model, collecting their global node IDs into unique sets that will be stored in 'faces'
        // tetrahedra first
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    // getting the global node IDs
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // inserting the global face node ids into the faces vector
                    faces.emplace_back( global_node_numbers );
                 }
               tets_++;
            }
        // hexahedra second
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    faces.emplace_back( global_node_numbers );
                 }
               hexes_++;
            }
        // prism third
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    faces.emplace_back( global_node_numbers );
                 }
               prisms_++;
            }
        // pyramids fourth
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    faces.emplace_back( global_node_numbers );
                 }
               pyras_++;
            }
         // no octahedra in CSMP
      }

    // sorting and eliminating duplicate faces from the vector
    sort( faces.begin(), faces.end() );
    faces.erase( unique( faces.begin(), faces.end() ), faces.end() );

    return faces.size();
    
 } // end CollectFaces
 
 
 

 
 
 
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectVolumes( const Model<dim>& model, vector<set<size_t>>& volumes )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // making sure that this is a three-dimensional model
    if constexpr ( dim != 3U ) {
         csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::CollectVolumes", "this method is only for three-dimensional models" );
         return 0U;
      }
    // making sure that the vector is empty and resizing it
    if ( !volumes.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectVolumes", "supplied element vector not empty; deleting previous content" );
         volumes.clear();
         volumes.reserve( model.Mesh().Elements() * 4U ); // approximate guess of the number of edges
      }

    const Region<dim>& model_domain(model.Region("Model"));

    // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
    vector<uint32_t> cnids;
    // tetrahedra first
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           volumes_.emplace_back( global_node_numbers );
        }
    // hexahedra
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           volumes_.emplace_back( global_node_numbers );
        }
    // prisms
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           volumes_.emplace_back( global_node_numbers );
        }
    // pyramids
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           volumes_.emplace_back( global_node_numbers );
        }

    return volumes.size();
 
 } // end CollectVolumes










/**
      Method determines the position of the edges in the current region within the global 'edges' vector.
      
      @param region_edges is where the results are stored.
      @return number of unique edges in the region
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectEdgesInRegion( const Region<dim>& domain, std::vector<size_t>& region_edges ) const
 {
    if ( !region_edges.empty() ) region_edges.clear();
 
    vector<uint32_t> snids;
    for ( const auto& it : domain.CellVector() ) {
      for ( uint32_t segm_id{0U}; segm_id<it->Segments(); ++segm_id ) {
             it->FE()->NodesOfSegment( segm_id, snids );
             assert( snids.size() == 2U );
             // getting the global node IDs
             pair<size_t,size_t> global_node_numbers{ it->N(snids[0])->Idx(), it->N(snids[1])->Idx() };
             // sorting the global node IDS in ascending order
             if ( global_node_numbers.first > global_node_numbers.second )
               swap(global_node_numbers.first,global_node_numbers.second);
             // finding the edge in global edges_ vector and measuring its distance from the beginning of it = ID
             const auto edge_it = lower_bound( edges_.begin(), edges_.end(), global_node_numbers );
             if ( edge_it != edges_.end() )
               // NB: ptrdiff_t might be better in region_edges than size_t
               region_edges.push_back( distance( edges_.begin(), edge_it ) );
          }
      }
 
    // sorting and eliminating duplicate edges from the vector
    sort( region_edges.begin(), region_edges.end() );
    region_edges.erase( unique( region_edges.begin(), region_edges.end() ), region_edges.end() );
    region_edges.shrink_to_fit();

    return region_edges.size();
    
 } // end CollectEdgesInRegion




template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectFacesInRegion( const Region<dim>& domain, std::vector<size_t>& region_faces ) const
 {
    // if the model is one-dimensional
    if constexpr ( dim == 1U ) return 0U;
    
    if ( !region_faces.empty() ) region_faces.clear();
    
    // looping over the elements of the model, collecting their global node IDs into unique sets that will be stored in 'faces'
    if constexpr ( dim == 2U ) {
        vector<uint32_t> cnids;
        for ( const auto& it : domain.CellVector() )
          if ( it->IsSurface() ) {
               // getting the global node IDs of the triangle or quadrilateral
               it->FE()->CornerNodes(cnids);
               set<size_t> global_node_numbers;
               for ( const auto& id : cnids ) global_node_numbers.insert( it->N(id)->Idx() );
               // looking for the position of the face in the global 'faces_' vector
               const auto face_it = lower_bound( faces_.begin(), faces_.end(), global_node_numbers );
               if ( face_it != faces_.end() )
               region_faces.push_back( distance( faces_.begin(), face_it ) );
            }
      }

    if constexpr ( dim == 3U ) {
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() ) {
            for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id )
              {
                 // getting the global node IDs
                 set<size_t> global_node_numbers;
                 for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                 // looking for the position of the face in the global 'faces_' vector
                 const auto face_it = lower_bound( faces_.begin(), faces_.end(), global_node_numbers );
                 if ( face_it != faces_.end() )
                 region_faces.push_back( distance( faces_.begin(), face_it ) );
              }
           }
      }

    // sorting and eliminating duplicate edges from the vector
    sort( region_faces.begin(), region_faces.end() );
    region_faces.erase( unique( region_faces.begin(), region_faces.end() ), region_faces.end() );
    region_faces.shrink_to_fit();

    return region_faces.size();

 } // end CollectFacesInRegion






template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectVolumesInRegion( const Region<dim>& domain, std::vector<size_t>& region_volumes ) const
 {
    // volumes exist only in 3D models
    if constexpr( dim != 3U ) return 0U;
    
    if ( !region_volumes.empty() ) region_volumes.clear();

    // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
    vector<uint32_t> cnids;
    for ( const auto& it : domain.CellVector() ) {
             it->FE()->CornerNodes(cnids);
             // getting the global node IDs
             set<size_t> global_node_numbers;
             for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
             // searching for the element IDx in the global volumes vector
             const auto elmt_it = lower_bound( volumes_.begin(), volumes_.end(), global_node_numbers );
             if ( elmt_it != volumes_.end() )
               region_volumes.push_back( distance(volumes_.begin(),elmt_it) );
      }

    // sorting to allow that duplicate edges can be eliminated from the vector
    sort( region_volumes.begin(), region_volumes.end() );
    region_volumes.erase( unique( region_volumes.begin(), region_volumes.end() ), region_volumes.end() );
    region_volumes.shrink_to_fit();

    return region_volumes.size();
     
 } // end CollectVolumesInRegion








template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectEdgesInBoundary( const Boundary<dim>& domain, std::vector<size_t>& boundary_edges ) const
 {
    if ( !boundary_edges.empty() ) boundary_edges.clear();
 
    vector<uint32_t> snids;
    for ( const auto& it : domain.CellVector() ) {
      for ( uint32_t segm_id{0U}; segm_id<it->FE()->Segments(); ++segm_id ) {
             it->FE()->NodesOfSegment( segm_id, snids );
             assert( snids.size() == 2U );
             // getting the global node IDs
             pair<size_t,size_t> global_node_numbers{ it->N(snids[0])->Idx(), it->N(snids[1])->Idx() };
             // sorting the global node IDS in ascending order
             if ( global_node_numbers.first > global_node_numbers.second )
               swap(global_node_numbers.first,global_node_numbers.second);
             // finding the edge in global edges_ vector and measuring its distance from the beginning of it = ID
             const auto edge_it = lower_bound( edges_.begin(), edges_.end(), global_node_numbers );
             if ( edge_it != edges_.end() )
               // NB: ptrdiff_t might be better in region_edges than size_t
               boundary_edges.push_back( distance( edges_.begin(), edge_it ) );
          }
      }
 
    // sorting and eliminating duplicate edges from the vector
    sort( boundary_edges.begin(), boundary_edges.end() );
    boundary_edges.erase( unique( boundary_edges.begin(), boundary_edges.end() ), boundary_edges.end() );
    boundary_edges.shrink_to_fit();

    return boundary_edges.size();
    
 } // end CollectEdgesInBoundary




template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectFacesInBoundary( const Boundary<dim>& domain, std::vector<size_t>& boundary_faces ) const
 {
    if ( dim != 3U )
      throw csmp::Exception( ERROR, "UG4_UGX_FileExport<dim>::CollectFacesInBoundary",
                            "Boundary faces only exist in 3D models; only call method for three-dimensional models" );
 
    if ( !boundary_faces.empty() ) boundary_faces.clear();
    
    // looping over the elements of the model, collecting their global node IDs into unique sets that will be stored in 'faces'
    for ( const auto& it : domain.CellVector() ) {
        for ( uint32_t face_id{0U}; face_id<it->FE()->Faces(); ++face_id )
          {
             // getting the global node IDs
             set<size_t> global_node_numbers;
             for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
             // looking for the position of the face in the global 'faces_' vector
             const auto face_it = lower_bound( faces_.begin(), faces_.end(), global_node_numbers );
             if ( face_it != faces_.end() )
             boundary_faces.push_back( distance( faces_.begin(), face_it ) );
          }
     }

    // sorting and eliminating duplicate edges from the vector
    sort( boundary_faces.begin(), boundary_faces.end() );
    boundary_faces.erase( unique( boundary_faces.begin(), boundary_faces.end() ), boundary_faces.end() );
    boundary_faces.shrink_to_fit();

    return boundary_faces.size();

 } // end CollectFacesInBoundary





/**
    Retrieves which element Idx()  the global face number corresponds to and reports this result in the output vector.
    @note the element idx() must match the UGX order counting the triangles first and them the quadrilaters
    
    @param cell_idx_in_attachment_order via cell idx of surface elements (first triangles, then quadrilaterals) in the order in which they should appear in the attachment vector
    
    @attention attachments are present only for faces that originated from Element or Face objects!
*/
template<uint32_t dim>
void UG4_UGX_FileExport<dim>::CellIdxInFaceAttachmentsOrder( const Model<dim>& model, vector<size_t>& cell_idx_in_attachment_order ) const
 {
    const Region<dim>& model_domain = model.Region("Model");
    
    // the cell indices are initialised to UINT_MAX so that uninitialised indices can be detected immediately
    if ( !cell_idx_in_attachment_order.empty() ) cell_idx_in_attachment_order.clear();
    cell_idx_in_attachment_order.resize( model.Mesh().Elements(), numeric_limits<size_t>::max() );
    
    // looping over the elements of the model, sorting their global node IDs into search keys needed to find their position in 'faces'
    size_t n_faces{0U};
    
    // triangles first
    for ( const auto& it : model_domain.CellVector() )
       if ( it->IsSurface() && isTriangular( it->FE_Type() ) ) {
            // getting the elements global corner node IDs
            vector<uint32_t> nids;
            it->FE()->CornerNodes(nids);
            set<size_t> global_node_numbers;
            for ( const auto& id : nids ) global_node_numbers.insert( it->N(id)->Idx() );
              // looking for the position of the face in the global 'faces_' vector
              const auto face_it = lower_bound( faces_.begin(), faces_.end(), global_node_numbers );
              if ( face_it != faces_.end() ) {
#ifdef DEBUG
                   auto idx = distance(faces_.begin(),face_it);
                   if ( idx >= cell_idx_in_attachment_order.size() ) {
                        cerr <<"\n"<<"UG4_UGX_FileExport<dim>::CellIdxInFaceAttachmentsOrder: triangle idx=";
                        cerr << idx <<" vs "<< cell_idx_in_attachment_order.size() << endl;
                     }
                   if ( idx < 0 ) {
                        cerr <<"\n"<<"UG4_UGX_FileExport<dim>::CellIdxInFaceAttachmentsOrder: negative cell index ";
                        cerr << idx << endl;
                     }
#endif
                   cell_idx_in_attachment_order.at( distance(faces_.begin(),face_it) ) = it->Idx();
                   n_faces++;
                }
           }
           
    // quadrilaterals
    for ( const auto& it : model_domain.CellVector() )
       if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) ) {
          // getting the elements global corner node IDs
          vector<uint32_t> nids;
          it->FE()->CornerNodes(nids);
          set<size_t> global_node_numbers;
          for ( const auto& id : nids ) global_node_numbers.insert( it->N(id)->Idx() );
            // looking for the position of the face in the global 'faces_' vector
            const auto face_it = lower_bound( faces_.begin(), faces_.end(), global_node_numbers );
            if ( face_it != faces_.end() ) {
#ifdef DEBUG
                   auto idx = distance(faces_.begin(),face_it);
                   if ( idx >= cell_idx_in_attachment_order.size() ) {
                        cerr <<"\n"<<"UG4_UGX_FileExport<dim>::CellIdxInFaceAttachmentsOrder: quadrilateral idx=";
                        cerr << idx <<" vs "<< cell_idx_in_attachment_order.size() << endl;
                     }
                   if ( idx < 0 ) {
                        cerr <<"\n"<<"UG4_UGX_FileExport<dim>::CellIdxInFaceAttachmentsOrder: negative cell index ";
                        cerr << idx << endl;
                     }
#endif
                 cell_idx_in_attachment_order.at( distance(faces_.begin(),face_it) ) = it->Idx();
                 n_faces++;
              }
         }
         
                    
    cell_idx_in_attachment_order.resize(n_faces);
    cell_idx_in_attachment_order.shrink_to_fit();

    // checking the vector for gaps (it should be contiguous up n_faces-1)
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance());
    if ( *max_element(cell_idx_in_attachment_order.begin(),cell_idx_in_attachment_order.end()) == numeric_limits<size_t>::max() ) {
         if ( cell_idx_in_attachment_order.size() < 100 ) {
              cerr <<"\n\nvector of cell indices of triangle and quadrilateral faces:\n";
              for ( const auto& it : cell_idx_in_attachment_order ) cerr <<" "<< it;
              cerr << endl;
           }
         csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::CellIdxInFaceAttachmentsOrder",
                       "gap detected in face vector which should be contiguous");
      }
        
 } // end CellIdxInFaceAttachmentsOrder





    /// writes the cell idx of volume elements (tets, hexa, prims, pyra) in the order in which they should appear in the attachment vector
template<uint32_t dim>
void UG4_UGX_FileExport<dim>::CellIdxInVolumeAttachmentsOrder( const Model<dim>& model, vector<size_t>& cell_idx_in_attachment_order ) const
 {
    const Region<dim>& model_domain = model.Region("Model");

    // the cell indices are initialised to UINT_MAX so that uninitialised indices can be detected immediately
    if ( !cell_idx_in_attachment_order.empty() ) cell_idx_in_attachment_order.clear();
    cell_idx_in_attachment_order.resize( model.Mesh().Elements(), numeric_limits<size_t>::max() );
    
    // looping over the elements of the model, sorting their global node IDs into search keys needed to find their position in 'faces'
    size_t n_volumes{0U};
    
    // tetrahedra first
    for ( const auto& it : model_domain.CellVector() )
       if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
            // getting the elements global corner node IDs
            vector<uint32_t> nids;
            it->FE()->CornerNodes(nids);
            set<size_t> global_node_numbers;
            for ( const auto& id : nids ) global_node_numbers.insert( it->N(id)->Idx() );
              // looking for the position of the face in the global 'faces_' vector
              const auto face_it = lower_bound( volumes_.begin(), volumes_.end(), global_node_numbers );
              if ( face_it != volumes_.end() ) {
                   cell_idx_in_attachment_order[ distance(volumes_.begin(),face_it)-1 ] = it->Idx();
                   n_volumes++;
                }
           }
           
    // hexahedra
    for ( const auto& it : model_domain.CellVector() )
       if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
          // getting the elements global corner node IDs
          vector<uint32_t> nids;
          it->FE()->CornerNodes(nids);
          set<size_t> global_node_numbers;
          for ( const auto& id : nids ) global_node_numbers.insert( it->N(id)->Idx() );
            // looking for the position of the face in the global 'faces_' vector
            const auto face_it = lower_bound( volumes_.begin(), volumes_.end(), global_node_numbers );
            if ( face_it != volumes_.end() ) {
                 cell_idx_in_attachment_order[ distance(volumes_.begin(),face_it)-1 ] = it->Idx();
                 n_volumes++;
              }
         }

    // prism elements
    for ( const auto& it : model_domain.CellVector() )
       if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
          // getting the elements global corner node IDs
          vector<uint32_t> nids;
          it->FE()->CornerNodes(nids);
          set<size_t> global_node_numbers;
          for ( const auto& id : nids ) global_node_numbers.insert( it->N(id)->Idx() );
            // looking for the position of the face in the global 'faces_' vector
            const auto face_it = lower_bound( volumes_.begin(), volumes_.end(), global_node_numbers );
            if ( face_it != volumes_.end() ) {
                 cell_idx_in_attachment_order[ distance(volumes_.begin(),face_it)-1 ] = it->Idx();
                 n_volumes++;
              }
         }

    // pyramid elements
    for ( const auto& it : model_domain.CellVector() )
       if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
          // getting the elements global corner node IDs
          vector<uint32_t> nids;
          it->FE()->CornerNodes(nids);
          set<size_t> global_node_numbers;
          for ( const auto& id : nids ) global_node_numbers.insert( it->N(id)->Idx() );
            // looking for the position of the face in the global 'faces_' vector
            const auto face_it = lower_bound( volumes_.begin(), volumes_.end(), global_node_numbers );
            if ( face_it != volumes_.end() ) {
                 cell_idx_in_attachment_order[ distance(volumes_.begin(),face_it)-1 ] = it->Idx();
                 n_volumes++;
              }
         }

    // CSMP has no octahedra
                    
    cell_idx_in_attachment_order.resize(n_volumes);
    cell_idx_in_attachment_order.shrink_to_fit();

    // checking the vector for gaps (it should be contiguous up n_faces-1)
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance());
    if ( *max_element(cell_idx_in_attachment_order.begin(),cell_idx_in_attachment_order.end()) == numeric_limits<size_t>::max() ) {
         if ( cell_idx_in_attachment_order.size() < 100 ) {
              cerr <<"\n\nvector of cell indices of volumetric element types:\n";
              for ( const auto& it : cell_idx_in_attachment_order ) cerr <<" "<< it;
              cerr << endl;
           }
         csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::CellIdxInVolumeAttachmentsOrder",
                       "gap detected in face vector which should be contiguous");
      }

 } // end CellIdxInVolumeAttachmentsOrder




/** Numbers elements by  type, followig the convention adopted by the UGX format:
 
     edges, triangles, quadrilaterals, tetrahedra, hex, prisms, pyramids, octahedra
     
This is inferred from the code in the UGX reader listing the order of elements and attachments:

if(strcmp(name, "vertices") == 0)
bSuccess = create_vertices(vertices, grid, curNode, aaPos);
else if(strcmp(name, "edges") == 0)
bSuccess = create_edges(edges, grid, curNode, vertices);

else if(strcmp(name, "triangles") == 0)
bSuccess = create_triangles(faces, grid, curNode, vertices);
else if(strcmp(name, "quadrilaterals") == 0)
bSuccess = create_quadrilaterals(faces, grid, curNode, vertices);

else if(strcmp(name, "tetrahedrons") == 0)
bSuccess = create_tetrahedrons(volumes, grid, curNode, vertices);
else if(strcmp(name, "hexahedrons") == 0)
bSuccess = create_hexahedrons(volumes, grid, curNode, vertices);
else if(strcmp(name, "prisms") == 0)
bSuccess = create_prisms(volumes, grid, curNode, vertices);
else if(strcmp(name, "pyramids") == 0)
bSuccess = create_pyramids(volumes, grid, curNode, vertices);
else if(strcmp(name, "octahedrons") == 0)
bSuccess = create_octahedrons(volumes, grid, curNode, vertices);

else if(strcmp(name, "vertex_attachment") == 0)
bSuccess = read_attachment<Vertex>(grid, curNode);
else if(strcmp(name, "edge_attachment") == 0)
bSuccess = read_attachment<Edge>(grid, curNode);
else if(strcmp(name, "face_attachment") == 0)
bSuccess = read_attachment<Face>(grid, curNode);
else if(strcmp(name, "volume_attachment") == 0)
bSuccess = read_attachment<Volume>(grid, curNode);
 */
template<uint32_t dim>
void UG4_UGX_FileExport<dim>::NumberElementsForUGX( const Model<dim>& model ) const
 {
    // 2D elements first: triangles, quadrilaterals
    const Region<dim>& model_domain = model.Region("Model");
    
    size_t n_face{0U};
    
    // ignoring line elements
    
    // triangles first
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsSurface() && isTriangular( it->FE_Type() ) )
        it->Idx( n_face++ );
        
    // quadrilaterals
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) )
        it->Idx( n_face++ );
        
    if constexpr (dim == 3U ) {
    
    size_t n_volume{0U};
        
    // tetra
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) )
        it->Idx( n_volume++ );

    // hexa
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isHexahedral( it->FE_Type() ) )
        it->Idx( n_volume++ );

    // prisms
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) )
        it->Idx( n_volume++ );

    // pyramids
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) )
        it->Idx( n_volume++ );
        
    // CSMP has no octahedra
    
    } // end 3D
 
 } // end NumberElementsForUGX


/*
template<uint32_t dim>
bool UG4_UGX_FileExport<dim>::HasTriangles( const Model<dim>& model ) const
 {
    const Region<dim>& model_domain = model.Region("Model");
    // quadrilaterals
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsSurface() && isTriangle( it->FE_Type() ) )
        return true;
        
    return false;
 }



template<uint32_t dim>
bool UG4_UGX_FileExport<dim>::HasQuadrilaterals( const Model<dim>& model ) const
 {
    const Region<dim>& model_domain = model.Region("Model");
    // quadrilaterals
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) )
        return true;
        
    return false;
 }
 
 
 
template<uint32_t dim>
bool UG4_UGX_FileExport<dim>::HasHexahedra( const Model<dim>& model ) const
 {
    if constexpr( dim != 3U ) return false;
    const Region<dim>& model_domain = model.Region("Model");
    // quadrilaterals
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) )
        return true;
        
    return false;
 }



template<uint32_t dim>
bool UG4_UGX_FileExport<dim>::HasPrisms( const Model<dim>& model ) const
 {
    if constexpr( dim != 3U ) return false;
    const Region<dim>& model_domain = model.Region("Model");
    // prisms
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) )
        return true;
        
    return false;
 }



template<uint32_t dim>
bool UG4_UGX_FileExport<dim>::HasPyramids( const Model<dim>& model ) const
 {
    if constexpr( dim != 3U ) return false;
    const Region<dim>& model_domain = model.Region("Model");
    // pyramids
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) )
        return true;
        
    return false;
 }
*/


    // translation of CSMP variable types to UG variable types
template<uint32_t dim>
string UG4_UGX_FileExport<dim>::UG_VariableType( const csmp::Index& index ) const
 {
    switch(index.type) {
         case SCALAR:
           return "double";
         case VECTOR:
           if constexpr ( dim == 2U )
             return "vector2";
           else
             return "vector3";
         case TENSOR:
           if constexpr ( dim == 2U )
             return "matrix22";
           else
             return "matrix33";
         case ARRAY:
         case FLAGGEDARRAY:
           return string("vector") + to_string(index.dataDepth);
         default:
           return "UNDEFINED";
      }
    return "UNDEFINED";
      
 } // end


template class UG4_UGX_FileExport<2U>;
template class UG4_UGX_FileExport<3U>;

} // end csmp
