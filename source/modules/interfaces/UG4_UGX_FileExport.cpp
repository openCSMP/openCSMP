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




// custom constructor
template<uint32_t dim>
UG4_UGX_FileExport<dim>::UG4_UGX_FileExport( const Model<dim>& model )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // if the model has P refinement or the order of the shape-function is greater than 1 (curvilinear elements) it cannot be output to UG4
    set<uint32_t> shape_function_order = model.Mesh().OrderOfShapeFunctions();
    if ( shape_function_order.size() > 1 || (*shape_function_order.begin()) != 1U ) {
         csmp_err.Note( ERROR, "UG4_UGX_FileExport (custom constructor)",
                       "UG4 - ugx format only supports output of straight-sided elements (shape-function order = 1; nothing was output");
         return;
      }
    
    // number nodes and elements globally in consecutive order
    model.Region("Model").UpdateMemberIndexes();

    cout<<"\nUG4_UGX_FileExport(ctor): collecting node numbers of cell edges into sorted 'edges_' vector...";
    size_t n_edges = CollectEdges( model );
    assert( n_edges > model.Mesh().Elements() );
    cout <<"  done: "<< n_edges <<" unique edges found."<<endl;
    
    cout<<"\nUG4_UGX_FileExport(ctor): collecting node numbers of cell faces into sorted 'faces_' vector...";
    size_t n_faces = CollectFaces( model );
    assert( n_faces > 0U );
    cout <<"  done: "<< n_faces <<" unique element faces found."<<endl;

    if constexpr( dim == 3U ) {
        cout<<"\nUG4_UGX_FileExport(ctor): collecting node numbers of cells into sorted 'volumes_' vector...";
        size_t n_volumes = CollectVolumes( model );
        assert( n_volumes > 0U && n_volumes <= model.Mesh().Elements() );
        cout <<"  done: "<< n_volumes <<" unique element IDs found."<<endl;
      }
 
 } // end constructor
  
  
  
  /**
      For the allocatio of no-data values (which will be set to negative float max)
  */
  static bool equivalentEntityInCSMP( size_t idx ) {
       if ( idx == numeric_limits<size_t>::max() ) return false;
       return true;
    }
 
 
 
 /**
          Writes the global node numbers of the elements, faces, and interfaces in counter-clockwise order.
          
          @attention function only generates output for elements, faces or interfaces that are present in CSMP
 */
 template<uint32_t dim,template<uint32_t> class CELL>
 static void writeNodeIndices( ofstream& ofs, const ModelSubDomain<dim,CELL>& subdomain, size_t cell_idx, bool& print_whitespace ) {
     if ( equivalentEntityInCSMP(cell_idx) )
     for ( uint32_t i{0U}; i<subdomain.E(cell_idx)->Nodes(); ++i ) {
         if ( print_whitespace ) ofs <<" ";
         else print_whitespace = true;
         ofs << subdomain.E(cell_idx)->N(i)->Idx();
      }
  }

 
 
 /**
            Writes corresponding number of no-data values to file
 */
static void writeNoDataValue( ofstream& ofs, const csmp::Index& var_key, uint32_t dim, bool omit_1st_whitespace ) {
      // where there is no property value, the no-data value will be written to file
      constexpr double no_data_value{ -numeric_limits<float>::max() };
      assert( ofs.is_open() );
      assert( var_key.type != ARRAY && var_key.type != FLAGGEDARRAY );
      const auto prec = ofs.precision();
      ofs.precision(1);
      switch ( var_key.type ) {
          case SCALAR: if ( omit_1st_whitespace ) ofs << no_data_value;
                       else ofs <<" "<< no_data_value;
            break;
          case VECTOR: for ( uint32_t i{0U}; i<dim; ++i ) {
                            if ( omit_1st_whitespace && i == 0U )
                              ofs << no_data_value;
                            else ofs <<" "<< no_data_value;
                         }
            break;
          case TENSOR: for ( uint32_t i{0U}; i<dim; ++i )
                         for ( uint32_t j{0U}; j<dim; ++j ) {
                              if ( omit_1st_whitespace && i == 0U && j == 0U )
                                ofs << no_data_value;
                              else ofs <<" "<< no_data_value;
                           }
            break;
          default:
            cerr <<"\nwriteNoDataValue: variable type could not be parsed\n";
        }
      ofs.precision(prec);
   }
   
   
/**
    For cell-based variables (Element, Face, InterFace), writes its value, considering SCALAR, VECTOR and TENSOR variables and printing no_data values as appropriate
    
    @param ofs  output textfile stream UTF-8 on linux
    @param var_key variable accessor in CSMP
    @param subdomain Region, Boundary, or SplitBoundary which contains the variable values of interest stored on Element, Face or InterFace, respectively.
    @param cell_idx continous numbering of elements, faces and interfaces in CSMP
    @param print_whitespace to avoid printing blanks at the  beginning of attachment arrays
    
      @attention method has to switch print_whitespace on
*/
template<uint32_t dim,template<uint32_t> class CELL>
static void writeCellVariableValue( ofstream& ofs, const csmp::Index& var_key,
                                    const ModelSubDomain<dim,CELL>& subdomain, size_t cell_idx, bool& print_whitespace )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );

    if ( equivalentEntityInCSMP(cell_idx) ) {
        // there no cell index was recorded for the face the no data values is writen
        switch(var_key.type) {
            case SCALAR: {
                   if ( !print_whitespace ) ofs << scientific << subdomain.E(cell_idx)->Read(var_key);
                   else ofs <<" "<< scientific << subdomain.E(cell_idx)->Read(var_key);
                   print_whitespace = true;
                }
              break;
            case VECTOR: {
                   VectorVariable<dim> vc;
                   subdomain.E(cell_idx)->Read(var_key,vc);
                   for ( uint32_t i{0u}; i<dim; ++i ) {
                         if ( print_whitespace ) ofs <<" ";
                         else print_whitespace = true;
                         ofs << scientific << vc[i];
                     }
                }
              break;
            case TENSOR: {
                   TensorVariable<dim> ts;
                   subdomain.E(cell_idx)->Read(var_key,ts);
                   for ( uint32_t i{0u}; i<dim; ++i )
                     for ( uint32_t j{0u}; j<dim; ++j ) {
                           if ( print_whitespace ) ofs <<" ";
                           else print_whitespace = true;
                           ofs << scientific << ts(i,j);
                       }
                }
              break;
            default: {
                 cerr <<"\nvariable of type "<< parseType( var_key.type ) << endl;
                 csmp_err.Note( ERROR, "writeCellVariableValue", "type of variable has no equivalent representation in UG4");
              }
          }
      }
    else writeNoDataValue( ofs, var_key, dim, print_whitespace );
         
} // end writeVariableValue
  
//template void writeCellVariableValue<2,Element>( ostream&, const csmp::Index&, ModelSubDomain<2,Element>&, size_t, bool );
  
  
 
 
 
 /**
   Orders the unique regions in the model such that the highest-dimensional ones come first.
   
   @note The UGX format requires that subsets are listed in an order of decreasing dimensionality!
   
   @attention mixed-dimensional regions are skipped because UG4 does not support them.
 */
 template<uint32_t dim>
 static void uniqueRegionsOrderedByDimensionality( const Model<dim>& model, list<string>& unique_regions_ordered ) {
      csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
      
      // key value pairs of highest-dim / region name
      multimap<int,string, greater<int> >  ordered_regions;
      
      // examining the unique regions
      for ( auto rit=model.UniqueRegionsBegin(); rit!=model.UniqueRegionsEnd(); ++rit )
        {
           pair<CELL_SHAPE,bool> cell_characteristics = (*rit).second.SingleCellShapeDomain();
           // of this region consists of only a single element type (line, volume, surface) it is elible for output to UG
           if ( cell_characteristics.second ) {
                ordered_regions.insert( make_pair( cell_characteristics.first, (*rit).first ) );
             }
           else csmp_err.Note( ERROR, "uniqueRegionsOrderedByDimensionality", (*rit).first,
                              "unique csmp::Region is multi-dimensional and is therefore not output as a subset to UG" );
        }
      
      // output
      if ( !unique_regions_ordered.empty() ) unique_regions_ordered.clear();
      if ( ordered_regions.empty() )
        csmp_err.Note( ERROR, "uniqueRegionsOrderedByDimensionality", "no unique regions found that are eligible for output to UG" );
      for ( const auto& it : ordered_regions )
        unique_regions_ordered.push_back( it.second );
        
   } // end uniqueRegionsOrderedByDimensionality
 
 
 
 
/**
    Writes UGX file in text format (UTF-8), in which all CSMP elements are assigned property values in the form of face_attachments and volume_attachments.
    Naive version where the XML tags are handcrafted.
    
    @param model CSMP model
    @param file_name filename without extension. .UGX  will be appended.
   
    @attention method assumes that the nodes and elements have already been numbered from 0..n-1 following UG's conventions
  
    @attention no_data_value =  -numeric_limits<float>::max()  will be written for each edge, face, and volume for which there is no equivalent in CSMP model
    
    @section Export Sequence (using the keywords that were defined by S. Reiter for the UGX format
    
    vertices
    edges
    faces (triangles, quads)
    volumes (tetra, hexa, prism, pyra)
    
    vertex_attachments (Dirichlet conditions etc.)
    @todo  edge_attachments to implement the line element properties from CSM
    face_attachments (surfaces elements)
    volume_attachments (volume elements)
    
    subsetblock - csmp::Region named objects
    @todo   subsetblock - csmp::Boundary named objects
    @todo   setsetblock - csmp::SplitBoundary objects
    
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
    // (the numbering is implicit to the order printed here and must be followed by the property attachments)
    ofs <<"\t<triangles>";
    // printing the nodes
    if ( HasTriangles() ) {
        print_whitespace = false;
        for ( const auto& it : tria_faces_ ) {
//             for ( const auto nid : it.first ) {
//                  if ( print_whitespace ) ofs <<" ";
//                  else print_whitespace = true;
//                  ofs << model_domain.N(nid)->Idx();
//               }
             // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
             writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
          }
        ofs <<"</triangles>"<< endl;
      }
    // quads
    if ( HasQuadrilaterals() ) {
        ofs <<"\t<quadrilaterals>";
        print_whitespace = false;
        for ( const auto& it : quad_faces_ ) {
//             for ( const auto nid : it.first ) {
//                  if ( print_whitespace ) ofs <<" ";
//                  else print_whitespace = true;
//                  ofs << model_domain.N(nid)->Idx();
//               }
             // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
             writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
          }
        ofs <<"</quadrilaterals>"<< endl;
      }
   
    // 5. Volumes (tet, hex, prism, pyra, octa)
    // ----------------------------------------
    if constexpr( dim == 3U ) {
        if ( HasTetrahedra() ) {
            print_whitespace = false;
            ofs <<"\t<tetrahedrons>";
            for ( const auto& it : tetra_volumes_ ) {
                 // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
                 writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
              }
             ofs <<"</tetrahedrons>"<< endl;
          }
        // hexa
        if ( HasHexahedra() ) {
            print_whitespace = false;
            ofs <<"\t<hexahedrons>";
            for ( const auto& it : hexa_volumes_ ) {
                 // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
                 writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
              }
            ofs <<"</hexahedrons>"<< endl;
          }
        // prisms
        if ( HasPrisms() ) {
            print_whitespace = false;
            ofs <<"\t<prisms>";
            for ( const auto& it : prism_volumes_ ) {
                 // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
                 writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
              }
            ofs <<"</prisms>"<< endl;
          }
        // pyramids
        if ( HasPyramids() ) {
            print_whitespace = false;
            ofs <<"\t<pyramids>";
            for ( const auto& it : pyra_volumes_ ) {
                 // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
                 writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
              }
            ofs <<"</pyramids>"<< endl;
          }
      }
      
    // TODO: any 'vertex_attachment' properties should go here
    // TODO: any 'edge_attachment' properties should go here
    
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
    for ( const auto& pit : output_props ) {
         const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
         // skipping variable types that are not supported by UG
         if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
           continue;
         ofs <<"\t<face_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key) <<"\" ";
         ofs <<"passOn=\"1\" global=\"1\">";

         // 0. TODO: line element (edge) attachments

         // 1. triangle element/face attachments
         // ------------------------------------
         print_whitespace = false;
         for ( const auto& idx : tria_faces_ )
           writeCellVariableValue( ofs, var_key, model_domain, idx.second, print_whitespace );

         // 2. quadrilateral element/face attachments
         // -----------------------------------------
         for ( const auto& idx : quad_faces_ )
           writeCellVariableValue( ofs, var_key, model_domain, idx.second, print_whitespace );

         ofs <<"</face_attachment>"<< endl;
         
      } // end for output properties


    // 6.2 volume attachments
    // ----------------------
    if constexpr ( dim == 3U )
      {
         for ( const auto& pit : output_props ) {
              const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
              // skipping variable types that are not supported by UG
              if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
                continue;
              ofs <<"\t<volume_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key);
              ofs <<"\" passOn=\"1\" global=\"1\">";

              // 1. property attachments to tetrahedra
              print_whitespace = false;
              for ( const auto& idx : tetra_volumes_ )
                writeCellVariableValue( ofs, var_key, model_domain, idx.second, print_whitespace );

               // 2. property attachments to hexahedra
              for ( const auto& idx : hexa_volumes_ )
                writeCellVariableValue( ofs, var_key, model_domain, idx.second, print_whitespace );

               // 3. property attachments to prism elements
              for ( const auto& idx : prism_volumes_ )
                writeCellVariableValue( ofs, var_key, model_domain, idx.second, print_whitespace );

               // 4. property attachments to pyramid elements
              for ( const auto& idx : pyra_volumes_ )
                writeCellVariableValue( ofs, var_key, model_domain, idx.second, print_whitespace );

             ofs <<"</volume_attachment>"<< endl;
              
           } // output properties
      } // 3D
    
    
    // 7. Subset Handlers and subsets
    // =============================================================================================
    // RULE 1: the subsets must be listed in an order of decreasing dimensionality: volumes, surfaces, edges
    // RULE 2: nodes, edges, and faces that may be part of multiple subsets must only be listed in one of them
    // (NB: TODO: check what consequences ignoring Rule 2 may have later in UG
    // (NB: subset do not have to include edges, faces, and volumes, but can only exist of vertices)
    // (NB: lower-dimensional element regions are communicated to UG via face-ony subsets)
    // (NB: each egde or node must only be present in one subset!!!)
    int colour{0}; // ten choices that will be used alternatingly, calling goodColour( int id )
    
    // CSMP MODEL REGIONS
    // ------------------
    list<string>  unique_regions_ordered;
    uniqueRegionsOrderedByDimensionality( model, unique_regions_ordered );
    if ( !unique_regions_ordered.empty() )
      {
        ofs <<"\t<subset_handler name=\"Regions\">"<< endl; // defSH means default subset handler (only one is needed)
        // model geometry-defining region by region
        for ( const auto& rit : unique_regions_ordered )
          {
             const Region<dim>& subdomain = model.Region(rit);
             // resetting the colour that the subset will be visualised with in PROMESH
             if ( colour == 10 ) colour = 9;
             const size_t domain_index{ subdomain.DomainIndex() };
             const auto   rgba = goodColour( colour );
             ofs <<"\t\t<subset name=\""<< rit <<"\" color=\"";
             ofs << get<0>(rgba) <<" "<< get<1>(rgba) <<" "<< get<2>(rgba) <<" "<< get<3>(rgba);
             ofs <<"\" state=\""<< domain_index <<"\">"<< endl;
             // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
             ofs <<"\t\t\t<vertices>";
             print_whitespace = false;
             for ( const auto& nit : subdomain.NodeVector() ) {
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
                 size_t n_edges = CollectEdgesInRegion( subdomain, ugx_id_numbers );
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
             
             // faces (triangles first, then quads, all numbered continuously)
             print_whitespace = false;
             if ( with_region_face_output_ ) {
                 size_t n_faces = CollectFacesInRegion( subdomain, ugx_id_numbers );
                 if ( n_faces > 0 ) {
                     ofs <<"\t\t\t<faces>";
                     // printing the IDs of all  faces that belong to this unique region
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
                     size_t n_volumes = CollectVolumesInRegion( subdomain, ugx_id_numbers );
                     if ( n_volumes > 0 ) {
                         ofs <<"\t\t\t<volumes>";
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
        
      } // end: loop over eligible regions
      
    // TODO: Missing BOUNDARIES (Face collections)
    // TODO: Missing SPLITBOUNDARIES ( InterFace collections)


    // considering the boundaries as subsets as well (this could be a separate subset handler)
    // TODO: loop over the outer model boundaries
    
    // END
    // --------------------------------------
    // closing the grid specification section
    ofs <<"</grid>"<< endl;
 
    cout <<"\n"<<"UG4_UGX_FileExport::Write_UGX_FileASCII: output file '"<< file_name <<".ugx' written successfully."<< endl;
    return true;
 
 } // end Write_UGX_FileASCII

 
  
  
  
  
  
/**
    Loops over model, making a vector of  the unique edges (node-iD pairs) in the model (corner nodes only).
   
    The elimination of duplicates and the sorting of the edges vector guarantees that the edges are always in the same order.
    
    @param model is used for data collection but not modified.
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectEdges( const Model<dim>& model )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // making sure that the vector is empty and resizing it
    if ( !edges_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectEdges", "supplied edge vector not empty; deleting previous content" );
         edges_.clear();
         edges_.reserve( model.Mesh().Elements() * 4U ); // approximate guess of the number of edges
      }
    
    // going over the unique regions in alphabetical order, tagging the edges found within them with region ids
    const Region<dim>& model_domain(model.Region("Model"));
    vector<uint32_t> snids;
    // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
    for ( const auto& it : model_domain.CellVector() ) {
      for ( uint32_t segm_id{0U}; segm_id<it->Segments(); ++segm_id ) {
             it->FE()->NodesOfSegment( segm_id, snids );
             assert( snids.size() == 2U );
             // getting the global node IDs
             pair<size_t,size_t> global_node_numbers{ it->N(snids[0])->Idx(), it->N(snids[1])->Idx() };
             // sorting the global node IDS in ascending order
             if ( global_node_numbers.first > global_node_numbers.second )
               swap(global_node_numbers.first,global_node_numbers.second);
             // inserting the edge node ids into the edges vector
             edges_.push_back( global_node_numbers );
          }
      }
 
    // sorting and eliminating duplicate edges from the edges_ vector
    // --------------------------------------------------------------
    // (using lambda function to search only the first pair in the pair-value pair)
//    sort( edges_.begin(), edges_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
    sort( edges_.begin(), edges_.end() );
//    edges_.erase( unique( edges_.begin(), edges_.end(), [](auto &left, auto &right) { return left.first == right.first; } ), edges_.end() );
    edges_.erase( unique( edges_.begin(), edges_.end() ), edges_.end() );

    return edges_.size();
    
 } // end CollectEdges
 
 

 
 
/**
     @attention order of nodes in 'faces' does not matter, only the global ID numbers are required by UG4
     @attention 'faces' related to CSMP boundaries must not be collected extra at this stage because they are already included since they match up with the element faces
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectFaces( const Model<dim>& model )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // making sure that the vector is empty and resizing it
    if ( !tria_faces_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectFaces", "supplied 'triangular faces' vector not empty; deleting previous content" );
         tria_faces_.clear();
         tria_faces_.reserve( model.Mesh().Elements() ); // approximate guess of the number of edges
      }
    if ( !quad_faces_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectFaces", "supplied 'quadrilateral faces' vector not empty; deleting previous content" );
         quad_faces_.clear();
         quad_faces_.reserve( model.Mesh().Elements() ); // approximate guess of the number of edges
      }
    
    const Region<dim>& model_domain(model.Region("Model"));

    // in 2D, triangular and quadrilateral elements make up the faces
    // --------------------------------------------------------------
    for ( const auto& it : model_domain.CellVector() )
        // triangular elements
        if ( it->IsSurface() ) {
             if ( isTriangular( it->FE_Type() ) )
              {
                 vector<uint32_t>  fnids;
                 it->FE()->CornerNodes( fnids );
                 // getting the global node IDs
                 set<size_t> global_node_numbers;
                 for ( const auto& id : fnids ) global_node_numbers.insert( it->N(id)->Idx() );
                 tria_faces_.push_back( make_pair( global_node_numbers, it->Idx() ) );
              }
            // quadrilaterals
            else if ( isQuadrilateral( it->FE_Type() ) )
              {
                 vector<uint32_t>  fnids;
                 it->FE()->CornerNodes( fnids );
                 set<size_t> global_node_numbers;
                 for ( const auto& id : fnids ) global_node_numbers.insert( it->N(id)->Idx() );
                 quad_faces_.push_back( make_pair( global_node_numbers, it->Idx() ) );
              }
            else
              csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectFaces", "surface face type could not be identified" );
         }

    // three-dimensional models (now the faces of volume elements make up the faces)
    // -----------------------------------------------------------------------------
    // (at this point, the faces defined in the CSMP model are already collected)
    if constexpr ( dim == 3U ) {
        // csmp::Boundary: consist of triangular and quadrilateral faces in 3D models
        // --------------------------------------------------------------------------
        for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit ) {
              // triangular elements
              if ( isTriangular( (*fit).FE_Type() ) )
                {
                   vector<uint32_t>  fnids;
                   (*fit).FE()->CornerNodes( fnids );
                   // getting the global node IDs
                   set<size_t> global_node_numbers;
                   for ( const auto& id : fnids ) global_node_numbers.insert( (*fit).N(id)->Idx() );
                   tria_faces_.push_back( make_pair( global_node_numbers, (*fit).Idx() ) );
                }
              // quadrilaterals
              else if ( isQuadrilateral( (*fit).FE_Type() ) )
                {
                   vector<uint32_t>  fnids;
                   (*fit).FE()->CornerNodes( fnids );
                   set<size_t> global_node_numbers;
                   for ( const auto& id : fnids ) global_node_numbers.insert( (*fit).N(id)->Idx() );
                   quad_faces_.push_back( make_pair( global_node_numbers, (*fit).Idx() ) );
                }
            }
        // TODO: include split boundaries here
    
        // creating faces between volumetric elements
        // ------------------------------------------
        // (these have no equivalent in CSMP and their IDs are therefore set to ~UINTMAX)
        // tetrahedra first
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    // getting the global node IDs
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // inserting the global face node ids into the faces vector
                    tria_faces_.push_back( make_pair( global_node_numbers, numeric_limits<size_t>::max() ) );
                 }
            }
        // hexahedra second
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    quad_faces_.push_back( make_pair( global_node_numbers, numeric_limits<size_t>::max() ) );
                 }
            }
        // prism third
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // dealing with quad and tria faces
                    if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) )
                      tria_faces_.push_back( make_pair( global_node_numbers, numeric_limits<size_t>::max() ) );
                    else
                      quad_faces_.push_back( make_pair( global_node_numbers, numeric_limits<size_t>::max() ) );
                 }
            }
        // pyramids fourth
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) )
                      tria_faces_.push_back( make_pair( global_node_numbers, numeric_limits<size_t>::max() ) );
                    else
                      quad_faces_.push_back( make_pair( global_node_numbers, numeric_limits<size_t>::max() ) );
                 }
            }
         // no octahedra in CSMP
      } // end 3d

    // sorting and eliminating duplicate faces from the vectors
    // (lambda function to search only the first pair in the pair-value pair)
    sort( tria_faces_.begin(), tria_faces_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
    // eliminating duplicate faces
    tria_faces_.erase( unique( tria_faces_.begin(), tria_faces_.end(), [](auto &left, auto &right) { return left.first == right.first; } ), tria_faces_.end() );
    //                                                                                                                 ^^^
    sort( quad_faces_.begin(), quad_faces_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
    quad_faces_.erase( unique( quad_faces_.begin(), quad_faces_.end(), [](auto &left, auto &right) { return left.first == right.first; } ), quad_faces_.end() );

    return tria_faces_.size() + quad_faces_.size();
    
 } // end CollectFaces
 
 

 
 
 
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectVolumes( const Model<dim>& model )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // making sure that this is a three-dimensional model
    if constexpr ( dim != 3U ) {
         csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::CollectVolumes", "this method is only for three-dimensional models" );
         return 0U;
      }
    // making sure that the vector is empty and resizing it
    if ( !tetra_volumes_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectVolumes", "supplied tetrahedral element vector not empty; deleting previous content" );
         tetra_volumes_.clear();
         tetra_volumes_.reserve( model.Mesh().Elements() ); // approximate guess of the number of edges
      }
    if ( !hexa_volumes_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectVolumes", "supplied hexahedral element vector not empty; deleting previous content" );
         hexa_volumes_.clear();
         hexa_volumes_.reserve( model.Mesh().Elements() ); // approximate guess of the number of edges
      }
    if ( !prism_volumes_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectVolumes", "supplied prism element vector not empty; deleting previous content" );
         prism_volumes_.clear();
         prism_volumes_.reserve( model.Mesh().Elements() ); // approximate guess of the number of edges
      }
    if ( !pyra_volumes_.empty() ) {
         csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectVolumes", "supplied pyramid element vector not empty; deleting previous content" );
         pyra_volumes_.clear();
         pyra_volumes_.reserve( model.Mesh().Elements() ); // approximate guess of the number of edges
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
           tetra_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    // hexahedra
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           hexa_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    // prisms
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           prism_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    // pyramids
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
           it->FE()->CornerNodes(cnids);
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           pyra_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }

    // TODO: sorting to match order in volumes_?

    return tetra_volumes_.size() + hexa_volumes_.size() + prism_volumes_.size() + pyra_volumes_.size();
 
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
size_t UG4_UGX_FileExport<dim>::CollectFacesInRegion( const Region<dim>& domain, vector<size_t>& region_faces ) const
 {
    // if the model is one-dimensional
    if constexpr ( dim == 1U ) return 0U;
    
    if ( !region_faces.empty() ) region_faces.clear();
    
    // TRIANGULAR FACES first
    // ----------------------
    // looping over the elements of the model, collecting their global node IDs into unique sets that will be stored in 'faces'
    // dealing with all the triangular faces first
    vector<uint32_t> cnids;
    for ( const auto& it : domain.CellVector() )
      if ( it->IsSurface() && isTriangular( it->FE_Type() ) )
        {
            // for the elements in the region
            it->FE()->CornerNodes(cnids);
            set<size_t> global_node_numbers;
            for ( const auto& id : cnids ) global_node_numbers.insert( it->N(id)->Idx() );
            // finding their position in the global 'faces_' vector
            const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(), make_pair( global_node_numbers, it->Idx() ),
                                              [](auto &left, auto &right) { return left.first < right.first; }  );
            if ( face_it != tria_faces_.end() )
              region_faces.push_back( distance( tria_faces_.begin(), face_it ) );
         }

    // adding the triangular faces that exist between the volumetric elements in 3D
    if constexpr ( dim == 3U ) {
        // tetrahedra
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    // getting the global node IDs
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // finding the global number of the face in the faces_ vector
                    const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                                      make_pair( global_node_numbers, numeric_limits<size_t>::max() ),
                                                      [](auto &left, auto &right) { return left.first < right.first; }  );
                    // inserting the intervening face into the region_faces vector
                    if ( face_it != tria_faces_.end() )
                       region_faces.push_back( distance( tria_faces_.begin(), face_it ) );
                 }
            }
        // prisms
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // dealing with quad and tria faces
                    if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) ) {
                        // finding the global number of the face in the faces_ vector
                        const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                                          make_pair( global_node_numbers, numeric_limits<size_t>::max() ),
                                                          [](auto &left, auto &right) { return left.first < right.first; }  );
                        // inserting the intervening face into the region_faces vector
                        if ( face_it != tria_faces_.end() )
                          region_faces.push_back( distance( tria_faces_.begin(), face_it ) );
                      }
                 }
            }
        // pyramids fourth
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) ) {
                        // finding the global number of the face in the faces_ vector
                        const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                                          make_pair( global_node_numbers, numeric_limits<size_t>::max() ),
                                                          [](auto &left, auto &right) { return left.first < right.first; }  );
                        // inserting the intervening face into the region_faces vector
                        if ( face_it != tria_faces_.end() )
                            region_faces.push_back( distance( tria_faces_.begin(), face_it ) );
                     }
                 }
            }
         // no octahedra in CSMP
      } // 3D
      
      
    // QUADRILATERAL FACES
    // -------------------
    // order is important because tria_faces_ needs to be established before quadrilaterals are accumulated
    for ( const auto& it : domain.CellVector() )
      if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) ) {
          it->FE()->CornerNodes(cnids);
          set<size_t> global_node_numbers;
          for ( const auto& id : cnids ) global_node_numbers.insert( it->N(id)->Idx() );
          const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(), make_pair( global_node_numbers, it->Idx() ),
                                            [](auto &left, auto &right) { return left.first < right.first; }  );
          if ( face_it != quad_faces_.end() )
          region_faces.push_back( distance( quad_faces_.begin(), face_it ) + tria_faces_.size() );
       }

    // adding the quadrilateral faces that exist between the volumetric elements in 3D
    if constexpr ( dim == 3U ) {
        // hexahedra
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // finding the global number of the face in the faces_ vector
                    const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                                      make_pair( global_node_numbers, numeric_limits<size_t>::max() ),
                                                      [](auto &left, auto &right) { return left.first < right.first; }  );
                    // inserting the intervening face into the region_faces vector
                    if ( face_it != quad_faces_.end() )
                       region_faces.push_back( distance( quad_faces_.begin(), face_it ) + tria_faces_.size() );
                 }
            }
        // prisms
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    // dealing with quad and tria faces
                    if ( isQuadrilateral(it->FE()->ElementTypeOfFace(face_id)) ) {
                        // finding the global number of the face in the faces_ vector
                        const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                                          make_pair( global_node_numbers, numeric_limits<size_t>::max() ),
                                                          [](auto &left, auto &right) { return left.first < right.first; }  );
                        // inserting the intervening face into the region_faces vector
                        if ( face_it != quad_faces_.end() )
                          region_faces.push_back( distance( quad_faces_.begin(), face_it ) + tria_faces_.size() );
                      }
                 }
            }
        // pyramids fourth
        for ( const auto& it : domain.CellVector() )
          if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> global_node_numbers;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
                    if ( isQuadrilateral(it->FE()->ElementTypeOfFace(face_id)) ) {
                        // finding the global number of the face in the faces_ vector
                        const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                                          make_pair( global_node_numbers, numeric_limits<size_t>::max() ),
                                                          [](auto &left, auto &right) { return left.first < right.first; }  );
                        // inserting the intervening face into the region_faces vector
                        if ( face_it != quad_faces_.end() )
                          region_faces.push_back( distance( quad_faces_.begin(), face_it ) + tria_faces_.size() );
                      }
                 }
            }
         // no octahedra in CSMP
      } // 3D


    // removing duplicates that may arise due to the intervening faces
    if constexpr ( dim == 3U ) {
        sort( region_faces.begin(), region_faces.end() );
        region_faces.erase( unique( region_faces.begin(), region_faces.end() ), region_faces.end() );
      }
    region_faces.shrink_to_fit();

    return region_faces.size();

 } // end CollectFacesInRegion





/**
       Collecting the volumetric elements in the region in the order tetra, hexa, prism, pyra.
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectVolumesInRegion( const Region<dim>& domain, std::vector<size_t>& region_volumes ) const
 {
    // volumes exist only in 3D models
    if constexpr( dim != 3U ) return 0U;
    
    if ( !region_volumes.empty() ) region_volumes.clear();

    // looping over the element edges, collecting their global node IDs into unique sets that will be stored in 'edges'
    vector<uint32_t> cnids;
    // tetrahedral elements
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) )
        {
           it->FE()->CornerNodes(cnids);
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( tetra_volumes_.begin(), tetra_volumes_.end(), make_pair( global_node_numbers, it->Idx() ) );
           if ( elmt_it != tetra_volumes_.end() )
             region_volumes.push_back( distance(tetra_volumes_.begin(),elmt_it) );
        }
    // hexahedral elements
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isHexahedral( it->FE_Type() ) )
        {
           it->FE()->CornerNodes(cnids);
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( hexa_volumes_.begin(), hexa_volumes_.end(), make_pair( global_node_numbers, it->Idx() ) );
           if ( elmt_it != hexa_volumes_.end() )
             region_volumes.push_back( distance(hexa_volumes_.begin(),elmt_it) + region_volumes.size() );
        }
    // prism elements
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) )
        {
           it->FE()->CornerNodes(cnids);
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( prism_volumes_.begin(), prism_volumes_.end(), make_pair( global_node_numbers, it->Idx() ) );
           if ( elmt_it != prism_volumes_.end() )
             region_volumes.push_back( distance(prism_volumes_.begin(),elmt_it) + region_volumes.size() );
        }
    // pyramids
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) )
        {
           it->FE()->CornerNodes(cnids);
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( const auto& nit : cnids ) global_node_numbers.insert( it->N(nit)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( pyra_volumes_.begin(), pyra_volumes_.end(), make_pair( global_node_numbers, it->Idx() ) );
           if ( elmt_it != pyra_volumes_.end() )
             region_volumes.push_back( distance(pyra_volumes_.begin(),elmt_it) + region_volumes.size() );
        }

    // no duplicates contained here and no sorting needed because order does not matter

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
//             const auto edge_it = lower_bound( edges_.begin(), edges_.end(), make_pair( global_node_numbers, it->Idx() ) );
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
    
    // boundary consists of triangular and quadrilateral faces these are unique
    // triangular faces first
    for ( const auto& it : domain.CellVector() ) {
        for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id )
          {
             // getting the global node IDs
             set<size_t> global_node_numbers;
             for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
             // distinguishing triangular from quadrilateral faces
             if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) ) {
                  const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(), make_pair( global_node_numbers, it->Idx() ) );
                  if ( face_it != tria_faces_.end() )
                  boundary_faces.push_back( distance( tria_faces_.begin(), face_it ) );
               }
           }
      }
    // quadrilaterals
    for ( const auto& it : domain.CellVector() ) {
        for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id )
          {
             // getting the global node IDs
             set<size_t> global_node_numbers;
             for ( const auto& nit : it->CornerNodesOfFace(face_id) ) global_node_numbers.insert( nit->Idx() );
             // distinguishing triangular from quadrilateral faces
             if ( isQuadrilateral(it->FE()->ElementTypeOfFace(face_id)) ) {
                  const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(), make_pair( global_node_numbers, it->Idx() ) );
                  if ( face_it != quad_faces_.end() )
                  boundary_faces.push_back( distance( quad_faces_.begin(), face_it ) + boundary_faces.size() );
               }
           }
      }

    boundary_faces.shrink_to_fit();

    return boundary_faces.size();

 } // end CollectFacesInBoundary










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
