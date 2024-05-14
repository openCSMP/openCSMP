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
        
    // number nodes and elements and faces globally in consecutive order
    const Region<dim>& model_domain = model.Region("Model");
    model_domain.UpdateMemberIndexes();
    size_t counter{ model_domain.Cells() };
    // number faces and interfaces consecutively in the order in which they are stored in the MeshManager
    for ( auto it=model.Mesh().FacesBegin(); it!=model.Mesh().FacesEnd(); ++it ) { (*it).Idx( counter++ ); }
    for ( auto it=model.Mesh().InterFacesBegin(); it!=model.Mesh().InterFacesEnd(); ++it ) { (*it).Idx( counter++ ); }
 
 } // end constructor
  
  
  
/**
      Finds corresponding CSMP line element or returns UNSPECIFIED
*/
template<uint32_t dim>
long UG4_UGX_FileExport<dim>::EquivalentEdgeInCSMP( size_t ug_idx ) const {
     const auto map_it = csmp_edges_inverted_.find( ug_idx );
     if ( map_it == csmp_edges_inverted_.end() ) return UNSPECIFIED;
     return (*map_it).second;
  }

/**
      Finds corresponding CSMP line element or returns UNSPECIFIED
*/
template<uint32_t dim>
long UG4_UGX_FileExport<dim>::EquivalentFaceInCSMP( size_t ug_idx ) const {
     const auto map_it = csmp_faces_inverted_.find( ug_idx );
     if ( map_it == csmp_faces_inverted_.end() ) return UNSPECIFIED;
     return (*map_it).second;
  }


 
 /**
     Writes the CSMP global node indices (0..nodes-1) in their correct cell order to the output stream.
     The indices are read from the supplied array nids.
 */
 template<size_t n_nodes>
 static void writeNodeIndices( ofstream& ofs, const array<size_t,n_nodes>& nids, bool& print_whitespace )
  {
     for ( const auto& nid : nids ) {
           if ( print_whitespace ) ofs <<" ";
           else print_whitespace = true;
           ofs << nid;
        }
       
  } // end writeNodeIndices (array version)


/**
  Taking the indexing from the model domain the node indices are written to file.
*/
template<uint32_t dim>
static void writeNodeIndices( ofstream& ofs, const Region<dim>& domain, size_t vol_elmt_idx, bool& print_whitespace )
  {
     assert( vol_elmt_idx < domain.Cells() );
     assert( domain.E(vol_elmt_idx)->IsVolume() );
     
     for ( uint32_t i{0U}; i<domain.E(vol_elmt_idx)->Nodes(); ++i ) {
           if ( print_whitespace ) ofs <<" ";
           else print_whitespace = true;
           ofs << domain.E(vol_elmt_idx)->N(i)->Idx();
        }
       
  } // end writeNodeIndices ("Model" region version)


 
 
 /**
      Where there is no property value, the no-data value will be written to file.
 */
static void writeNoDataValue( ofstream& ofs, const csmp::Index& var_key, uint32_t dim, bool omit_1st_whitespace ) {
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
   


template<uint32_t dim>
void UG4_UGX_FileExport<dim>::WriteVertexVariableValue( ofstream& ofs, const csmp::Index& var_key,
                                                        const Region<dim>& domain, size_t cell_idx,
                                                        bool& print_whitespace ) const
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );

    if ( var_key.place != NODE ) {
        // there no cell index was recorded for the face the no data values is writen
        switch(var_key.type) {
            case SCALAR: {
                   double sc;
                   switch(var_key.place) {
                        case NODE:
                             assert( cell_idx < domain.Nodes() );
                             sc = domain.N(cell_idx)->Read(var_key);
                          break;
                        default:
                             sc = -numeric_limits<float>::max(); // no-data value
                     }
                   if ( !print_whitespace ) ofs << scientific << sc;
                   else ofs <<" "<< scientific << sc;
                   print_whitespace = true;
                }
              break;
            case VECTOR: {
                   VectorVariable<dim> vc;
                   switch(var_key.place) {
                        case NODE:
                             assert( cell_idx < domain.Nodes() );
                             domain.N(cell_idx)->Read(var_key,vc);
                          break;
                        default:
                             vc = -numeric_limits<float>::max(); // no-data value
                     }
                   for ( uint32_t i{0u}; i<dim; ++i ) {
                         if ( print_whitespace ) ofs <<" ";
                         else print_whitespace = true;
                         ofs << scientific << vc[i];
                     }
                }
              break;
            case TENSOR: {
                   TensorVariable<dim> ts;
                   switch(var_key.place) {
                        case NODE:
                             assert( cell_idx < domain.Nodes() );
                             domain.N(cell_idx)->Read(var_key,ts);
                          break;
                        default:
                             ts = -numeric_limits<float>::max(); // no-data value
                     }
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
                 csmp_err.Note( ERROR, "WriteFaceVariableValue", "type of variable has no equivalent representation in UG4");
              }
          }
      }
    else writeNoDataValue( ofs, var_key, dim, print_whitespace );
    
 } // end WriteVertexVariableValue





template<uint32_t dim>
void UG4_UGX_FileExport<dim>::WriteEdgeVariableValue( ofstream& ofs, const csmp::Index& var_key,
                                                      const Region<dim>& domain, size_t cell_idx,
                                                      bool& print_whitespace ) const
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    if ( !domain.E(cell_idx)->IsLine() ) {
         csmp_err.Note( ERROR, "UG4_UGX_FileExport<dim>::WriteEdgeVariableValue",
                        "Attempt to apply method to an non-line element; Nothing was done");
         return;
      }
    
     if ( var_key.place == ELEMENT ) {
        // there no cell index was recorded for the face the no data values is writen
        switch(var_key.type) {
            case SCALAR: {
                   double sc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < domain.Cells() );
                             sc = domain.E(cell_idx)->Read(var_key);
                          break;
                        default:
                             sc = -numeric_limits<float>::max(); // no-data value
                     }
                   if ( !print_whitespace ) ofs << scientific << sc;
                   else ofs <<" "<< scientific << sc;
                   print_whitespace = true;
                }
              break;
            case VECTOR: {
                   VectorVariable<dim> vc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < domain.Cells() );
                             domain.E(cell_idx)->Read(var_key,vc);
                          break;
                        default:
                             vc = -numeric_limits<float>::max(); // no-data value
                     }
                   for ( uint32_t i{0u}; i<dim; ++i ) {
                         if ( print_whitespace ) ofs <<" ";
                         else print_whitespace = true;
                         ofs << scientific << vc[i];
                     }
                }
              break;
            case TENSOR: {
                   TensorVariable<dim> ts;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < domain.Cells() );
                             domain.E(cell_idx)->Read(var_key,ts);
                          break;
                        default:
                             ts = -numeric_limits<float>::max(); // no-data value
                     }
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
                 csmp_err.Note( ERROR, "WriteEdgeVariableValue", "type of variable has no equivalent representation in UG4");
              }
          }
      }
    else writeNoDataValue( ofs, var_key, dim, print_whitespace );
         
} // end WriteEdgeVariableValue





/**
    For cell-based variables (Element, Face, InterFace), writes its value, considering SCALAR, VECTOR and TENSOR variables and printing no_data values as appropriate.
    The variable placement in the CSMP Index is the one which will be considered.
    
    @param ofs  output textfile stream UTF-8 on linux
    @param var_key variable accessor in CSMP
    @param mesh MeshManager which contains the variable values of interest stored on Element, Face or InterFace, respectively.
    @param ug_cell_idx continous numbering of faces in CSMP
    @param print_whitespace to avoid printing blanks at the  beginning of attachment arrays
    
      @attention method has to switch print_whitespace on
*/
template<uint32_t dim>
void UG4_UGX_FileExport<dim>::WriteFaceVariableValue( ofstream& ofs, const csmp::Index& var_key,
                                                      const MeshManager<dim>& mesh, size_t ug_cell_idx,
                                                      bool& print_whitespace ) const
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    const long cell_idx = EquivalentFaceInCSMP( ug_cell_idx );

    if ( cell_idx != UNSPECIFIED ) {
        // there no cell index was recorded for the face the no data values is writen
        switch(var_key.type) {
            case SCALAR: {
                   double sc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             sc = (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key);
                          break;
                        case FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces() );
                             sc = (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Read(var_key);
                          break;
                        case INTER_FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() );
                             sc = (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Read(var_key);
                          break;
                        default:
                             sc = -numeric_limits<float>::max(); // no-data value
                     }
                   if ( !print_whitespace ) ofs << scientific << sc;
                   else ofs <<" "<< scientific << sc;
                   print_whitespace = true;
                }
              break;
            case VECTOR: {
                   VectorVariable<dim> vc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key,vc);
                          break;
                        case FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces() );
                             (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Read(var_key,vc);
                          break;
                        case INTER_FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() );
                             (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Read(var_key,vc);
                          break;
                        default:
                             vc = -numeric_limits<float>::max(); // no-data value
                     }
                   for ( uint32_t i{0u}; i<dim; ++i ) {
                         if ( print_whitespace ) ofs <<" ";
                         else print_whitespace = true;
                         ofs << scientific << vc[i];
                     }
                }
              break;
            case TENSOR: {
                   TensorVariable<dim> ts;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key,ts);
                          break;
                        case FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces() );
                             (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Read(var_key,ts);
                          break;
                        case INTER_FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() );
                             (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Read(var_key,ts);
                          break;
                        default:
                             ts = -numeric_limits<float>::max(); // no-data value
                     }
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
                 csmp_err.Note( ERROR, "WriteFaceVariableValue", "type of variable has no equivalent representation in UG4");
              }
          }
      }
    else writeNoDataValue( ofs, var_key, dim, print_whitespace );
         
} // end WriteFaceVariableValue
  



/**
    Same as previous method, but for volumetric 3D elements only.
    @attention method assumes that there is an exact match between CSMP volume cell index and UG volume index.
*/
template<uint32_t dim>
void UG4_UGX_FileExport<dim>::WriteVolumeVariableValue( ofstream& ofs, const csmp::Index& var_key,
                                                        const MeshManager<dim>& mesh, size_t ug_cell_idx,
                                                        bool& print_whitespace ) const
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    const long cell_idx = ug_cell_idx;

    if ( cell_idx != UNSPECIFIED && var_key.place == ELEMENT ) {
        // there no cell index was recorded for the face the no data values is writen
        switch(var_key.type) {
            case SCALAR: {
                   double sc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             sc = (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key);
                          break;
                        default:
                             sc = -numeric_limits<float>::max(); // no-data value
                     }
                   if ( !print_whitespace ) ofs << scientific << sc;
                   else ofs <<" "<< scientific << sc;
                   print_whitespace = true;
                }
              break;
            case VECTOR: {
                   VectorVariable<dim> vc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key,vc);
                          break;
                        default:
                             vc = -numeric_limits<float>::max(); // no-data value
                     }
                   for ( uint32_t i{0u}; i<dim; ++i ) {
                         if ( print_whitespace ) ofs <<" ";
                         else print_whitespace = true;
                         ofs << scientific << vc[i];
                     }
                }
              break;
            case TENSOR: {
                   TensorVariable<dim> ts;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key,ts);
                          break;
                        default:
                             ts = -numeric_limits<float>::max(); // no-data value
                     }
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
                 csmp_err.Note( ERROR, "WriteVolumeVariableValue", "type of variable has no equivalent representation in UG4");
              }
          }
      }
    else writeNoDataValue( ofs, var_key, dim, print_whitespace );
         
} // end WriteVolumeVariableValue




  
 
 
 
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
 




 template<uint32_t dim>
 static void boundariesOrderedByDimensionality( const Model<dim>& model, list<string>& boundaries_ordered_by_dim ) {
      csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
      
      // key value pairs of highest-dim / region name
      multimap<int,string, greater<int> >  ordered_boundaries;
      
      // examining the unique regions
      for ( auto rit=model.BoundariesBegin(); rit!=model.BoundariesEnd(); ++rit )
        {
           pair<CELL_SHAPE,bool> cell_characteristics = (*rit).second.SingleCellShapeDomain();
           // of this region consists of only a single element type (line, volume, surface) it is elible for output to UG
           if ( cell_characteristics.second ) {
                ordered_boundaries.insert( make_pair( cell_characteristics.first, (*rit).first ) );
             }
           else csmp_err.Note( ERROR, "boundariesOrderedByDimensionality", (*rit).first,
                              "unique csmp::Boundary is multi-dimensional and is therefore not output as a subset to UG" );
        }
      
      // output
      if ( !boundaries_ordered_by_dim.empty() ) boundaries_ordered_by_dim.clear();
      if ( ordered_boundaries.empty() )
        csmp_err.Note( WARNING, "boundariesOrderedByDimensionality", "no unique boundaries found that are eligible for output to UG" );
      for ( const auto& it : ordered_boundaries )
        boundaries_ordered_by_dim.push_back( it.second );
        
   } // end boundariesOrderedByDimensionality
   
   


 
 
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
bool UG4_UGX_FileExport<dim>::Write_UGX_FileASCII( const Model<dim>& model, const std::string& file_name )
 {
    //csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    ofstream ofs( (file_name + ".ugx").c_str() );
    
    // 0. Creating the datastructures required by UG4
    // ----------------------------------------------
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
    // (consecutively numbered across the different element types)
    // (the numbering is implicit to the order printed here and must be followed by the property attachments as well)
    // 4.1 triangles that were already present in CSMP model
    // printing the nodes
    if ( HasTriangles() ) {
        ofs <<"\t<triangles>";
        print_whitespace = false;
        for ( const auto& it : tria_faces_ )
          // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
          writeNodeIndices( ofs, it.second, print_whitespace );
        ofs <<"</triangles>"<< endl;
      }

    // 4.2 Quadrilateral faces already present in CSMP model
    if ( HasQuadrilaterals() ) {
        ofs <<"\t<quadrilaterals>";
        print_whitespace = false;
        for ( const auto& it : quad_faces_ ) {
             // print node numbers in csmp order, but only for elements, faces, interfaces present in csmp mesh
             writeNodeIndices( ofs, it.second, print_whitespace );
          }
        ofs <<"</quadrilaterals>"<< endl;
      }


    // 5. Volumes
    // ----------
    // (consecutively numbered across the different element types)
    if constexpr ( dim == 3U ) {
        // tetrahedra
        if ( !tetra_volumes_.empty() ) {
             ofs <<"\t<tetrahedrons>";
             print_whitespace = false;
             for ( const auto& it : tetra_volumes_ )
               writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
             ofs <<"</tetrahedrons>"<< endl;
          }
        // hexahedra
        if ( !hexa_volumes_.empty() ) {
             ofs <<"\t<hexahedrons>";
             print_whitespace = false;
             for ( const auto& it : hexa_volumes_ )
               writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
             ofs <<"</hexahedrons>"<< endl;
          }
        // prisms
        if ( !prism_volumes_.empty() ) {
             ofs <<"\t<prisms>";
             print_whitespace = false;
             for ( const auto& it : prism_volumes_ )
               writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
             ofs <<"</prisms>"<< endl;
          }
        // pyramids
        if ( !pyra_volumes_.empty() ) {
             ofs <<"\t<pyramids>";
             print_whitespace = false;
             for ( const auto& it : pyra_volumes_ )
               writeNodeIndices( ofs, model_domain, it.second, print_whitespace );
             ofs <<"</pyramids>"<< endl;
          }
        // octahedra
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
    
    // 6.0 vertex attachments (assuming csmp nodes are equivalent to UG vertices)
    // --------------------------------------------------------------------------
    for ( const auto& pit : output_props ) {
         const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
         // skipping variable types that are not supported by UG
         if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
           continue;
         if ( var_key.place != NODE )
           continue;
         // vertex = node properties
         // ------------------------
             ofs <<"\t<vertex_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key) <<"\" ";
             ofs <<"passOn=\"1\" global=\"1\">";
             print_whitespace = false;
             for ( size_t n{0U}; n<model.Mesh().Nodes(); ++n )
               WriteVertexVariableValue( ofs, var_key, model_domain, n++, print_whitespace );
             ofs <<"</vertex_attachment>"<< endl;
      }
      
    // 6.1 edge attachments: CSMP line elements with element properties
    // ----------------------------------------------------------------
    // (no-data values are written everywhere else)
    for ( const auto& pit : output_props ) {
         const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
         // skipping variable types that are not supported by UG
         if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
           continue;
         if ( var_key.place != ELEMENT )
           continue;
         ofs <<"\t<edge_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key) <<"\" ";
         ofs <<"passOn=\"1\" global=\"1\">";
         print_whitespace = false;
         for ( size_t i{0U}; i<edges_.size(); i++ ) {
              const auto matching_edge = csmp_edges_inverted_.find(i);
              if ( matching_edge != csmp_edges_inverted_.end() )
                WriteEdgeVariableValue( ofs, var_key, model_domain, (*matching_edge).second, print_whitespace );
              else {
                   writeNoDataValue( ofs, var_key, dim, !print_whitespace );
                   print_whitespace = true;
                }
           }
         ofs <<"</edge_attachment>"<< endl;
     }
    
    // 6.2 face attachments (first triangles, then quadrilaterals)
    // -----------------------------------------------------------
    for ( const auto& pit : output_props ) {
         const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
         // skipping variable types that are not supported by UG
         if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
           continue;
         if constexpr ( dim == 3U )
           // In 3D, only variables placed on the Face will be considered
           if ( var_key.place != FACE )
             continue;

         ofs <<"\t<face_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key) <<"\" ";
         ofs <<"passOn=\"1\" global=\"1\">";

         // 6.1 triangle element/face attachments
         // -------------------------------------
         print_whitespace = false;
         for ( size_t n_ug_face{0U}; n_ug_face<tria_faces_.size(); ++n_ug_face )
           WriteFaceVariableValue( ofs, var_key, model.Mesh(), n_ug_face, print_whitespace );

         // 6.2 quadrilateral element/face attachments
         // ------------------------------------------
         for ( size_t n_ug_face{0U}; n_ug_face<quad_faces_.size(); ++n_ug_face )
           WriteFaceVariableValue( ofs, var_key, model.Mesh(), n_ug_face + tria_faces_.size(), print_whitespace );

         ofs <<"</face_attachment>"<< endl;

      } // end for output properties


    // 6.3 volume attachments
    // ----------------------
    // (the assumption here is that the CSMP model contains exactly the same volumetric elements as the UG model)
    if constexpr ( dim == 3U )
      {
         for ( const auto& pit : output_props ) {
              const csmp::Index var_key = model.Database().StorageKey(pit.c_str());
              // skipping variable types that are not supported by UG
              if ( var_key.type == ARRAY || var_key.type == FLAGGEDARRAY )
                continue;
              // only variables placed on the Element will be considered for now
              if ( var_key.place != ELEMENT )
                continue;
              ofs <<"\t<volume_attachment name=\""<< pit <<"\" type=\""<< UG_VariableType(var_key);
              ofs <<"\" passOn=\"1\" global=\"1\">";

              // 1. property attachments to tetrahedra
              print_whitespace = false;
              for ( const auto& idx : tetra_volumes_ )
                WriteVolumeVariableValue( ofs, var_key, model.Mesh(), idx.second, print_whitespace );

               // 2. property attachments to hexahedra
              for ( const auto& idx : hexa_volumes_ )
                WriteVolumeVariableValue( ofs, var_key, model.Mesh(), idx.second, print_whitespace );

               // 3. property attachments to prism elements
              for ( const auto& idx : prism_volumes_ )
                WriteVolumeVariableValue( ofs, var_key, model.Mesh(), idx.second, print_whitespace );

               // 4. property attachments to pyramid elements
              for ( const auto& idx : pyra_volumes_ )
                WriteVolumeVariableValue( ofs, var_key, model.Mesh(), idx.second, print_whitespace );

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
             const auto domain_index{ subdomain.DomainIndex() };
             const auto rgba = goodColour( colour );
             ofs <<"\t\t<subset name=\""<< rit <<"\" color=\"";
             ofs << get<0>(rgba) <<" "<< get<1>(rgba) <<" "<< get<2>(rgba) <<" "<< get<3>(rgba);
             ofs <<"\" state=\""<< domain_index <<"\">"<< endl;
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
                 n_edges = CollectEdgesInRegion( subdomain, ugx_id_numbers );
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
                 n_faces = CollectFacesInRegion( subdomain, ugx_id_numbers );
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
                         ofs <<"</volumes>"<< endl;
                       }
                   }
               }
             ofs <<"\t\t</subset>"<< endl;
             colour++;
         }
          
         // closing the subset section
        ofs <<"\t</subset_handler>"<< endl;
        
      } // end: loop over eligible regions
      
      
    // CSMP MODEL BOUNDARIES
    // ---------------------
    list<string>  boundaries_ordered_by_dimensionality;
    boundariesOrderedByDimensionality( model, boundaries_ordered_by_dimensionality );
    if ( !boundaries_ordered_by_dimensionality.empty() )
      {
        ofs <<"\t<subset_handler name=\"Boundaries\">"<< endl;
        // model geometry-defining region by region
        for ( const auto& bit : boundaries_ordered_by_dimensionality )
          {
             const Boundary<dim>& subdomain = model.Boundary(bit);
             // resetting the colour that the subset will be visualised with in PROMESH
             if ( colour == 10 ) colour = 9;
             const auto domain_index{ subdomain.DomainIndex() };
             const auto rgba = goodColour( colour );
             ofs <<"\t\t<subset name=\""<< bit <<"\" color=\"";
             ofs << get<0>(rgba) <<" "<< get<1>(rgba) <<" "<< get<2>(rgba) <<" "<< get<3>(rgba);
             ofs <<"\" state=\""<< domain_index <<"\">"<< endl;
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
             if ( with_boundary_edge_output_ ) {
                 print_whitespace = false;
                 // compiling the edges that form part of the region
                 size_t n_bedges = CollectEdgesInBoundary( subdomain, ugx_id_numbers );
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
             if constexpr ( dim == 3U ) {
                 if ( with_boundary_face_output_ ) {
                     size_t n_bfaces = CollectFacesInBoundary( subdomain, ugx_id_numbers );
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
               }
  
             ofs <<"\t\t</subset>"<< endl;
             colour++;
         }
          
         // closing the subset section
        ofs <<"\t</subset_handler>"<< endl;
        
      } // end: loop over eligible boundaries
      
      
    // TODO: Missing SPLITBOUNDARIES ( InterFace collections)

    // END
    // --------------------------------------
    // closing the grid specification section
    ofs <<"</grid>"<< endl;
 
    cout <<"\n"<<"UG4_UGX_FileExport::Write_UGX_FileASCII: output file '"<< file_name <<".ugx' written successfully."<< endl;
    return true;
 
 } // end Write_UGX_FileASCII

 
  
  
  
  
  
/**
    Loops over region "Model", making a vector of  unique edges (node-iD pairs) storing their corner nodes only.
   
    Elimination of duplicates and the sorting of the edges vector guarantees that the edges are always in the same order.
    The mapping between the UG4 edge IDs and CSMP line elements co-inciding with UG4 edges is created after the
    sorting and elimination of duplicates.
    
    @param model is used for data collection but not modified
    
    @todo deal with potential Edge objects in 3D models
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
    
    const Region<dim>& model_domain(model.Region("Model"));
    vector<uint32_t> snids;
    // looping over all equidimensional element edges in the Model,
    // collecting their global node IDs into unique pairs that will be stored in 'edges' vector as search keys
    if constexpr ( dim == 2U ) {
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsSurface() ) {
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
      }
    else if constexpr ( dim == 3U ) {
        for ( const auto& it : model_domain.CellVector() )
          if ( it->Volume() ) {
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
      }
      
    // sorting and eliminating duplicate edges from the edges_ vector
    // --------------------------------------------------------------
    sort( edges_.begin(), edges_.end() );
    edges_.erase( unique( edges_.begin(), edges_.end() ), edges_.end() );
    
    // singling out the edges which are also contained in the CSMP model
    // -----------------------------------------------------------------
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsLine() ) {
           // create search key
           pair<size_t,size_t> global_node_numbers{ it->N(0U)->Idx(), it->N(1U)->Idx() };
           // sort global node IDS in ascending order
           if ( global_node_numbers.first > global_node_numbers.second )
             swap(global_node_numbers.first,global_node_numbers.second);
             
           // search edge which must be contained in edge vector; recording its offset from origin in map
           const auto edge_it = lower_bound( edges_.begin(), edges_.end(), global_node_numbers );
           assert( edge_it != edges_.end() );
           // CSMP element idx, edge # in UG
           csmp_edges_.insert( make_pair( it->Idx(), distance(edges_.begin(), edge_it) ) );
       }
       
    // creating inverse map of csmp edges
    // ----------------------------------
    for ( const auto& it : csmp_edges_ )
      // edge_id, csmp_idx
      csmp_edges_inverted_.insert( make_pair( it.second, it.first ) );

    // finding those edges that have a representation in the CSMP model
    // ----------------------------------------------------------------
    if ( edges_in_csmp_.empty() ) edges_in_csmp_.clear();
    edges_in_csmp_.resize( edges_.size(), false );
    // initialising the vector
    for ( const auto& it : csmp_edges_ )
      edges_in_csmp_[ it.second ] = true;

    return edges_.size();
    
 } // end CollectEdges
 

// FOR EDGES with element ID stored within them
// (using lambda function to search only the first pair in the pair-value pair)
//    sort( edges_.begin(), edges_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
//    edges_.erase( unique( edges_.begin(), edges_.end(), [](auto &left, auto &right) { return left.first == right.first; } ), edges_.end() );

 
 
 
 
 
 
 /**
     This method collects all potential faces in the model into a face array: triangles first, quadrilaterals next.
     For UG faces that correspond to surface elements in CSMP model, their ID is recorded in a map for later retrieval.

     @attention Order of nodes in 'faces' does not matter, only the global ID numbers are required by UG4
     @attention 'faces' related to CSMP boundaries must not be collected extra at this stage because they are already included since they match up with the element faces
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectFaces( const Model<dim>& model )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );
    
    // emptying and resizing vectors
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
    // (NB: these faces will be unique and already have a number in the CSMP mesh.
    //      therefore they can be printed with correct orientation and associated properties)
    if constexpr ( dim == 2U )
      {
        for ( const auto& it : model_domain.CellVector() )
            if ( it->IsSurface() ) {
                 // triangular elements
                 if ( isTriangular( it->FE_Type() ) )
                  {
                     // getting the global node IDs
                     set<size_t> search_key;
                     array<size_t,3U> global_node_numbers;
                     for ( uint32_t i{0U}; i<3U; ++i ) {
                           global_node_numbers[i] = it->N(i)->Idx();
                           search_key.insert( global_node_numbers[i] );
                       }
                     tria_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                  }
                // quadrilaterals
                else if ( isQuadrilateral( it->FE_Type() ) )
                  {
                     set<size_t> search_key;
                     array<size_t,4U> global_node_numbers;
                     for ( uint32_t i{0U}; i<4U; ++i ) {
                           global_node_numbers[i] = it->N(i)->Idx();
                           search_key.insert( global_node_numbers[i] );
                       }
                     quad_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                  }
                else
                  csmp_err.Note( WARNING, "UG4_UGX_FileExport<dim>::CollectFaces", "surface face type could not be identified" );
             }
      } // end 2D
  
  
    // in 3D, faces are generated from the faces of the volumetric elements in the mesh
    // --------------------------------------------------------------------------------
    if constexpr ( dim == 3U )
      {
        // ----------------------
        // triangular faces
        // ----------------------
        // tetrahedra
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    // getting the global node IDs to form a search key
                    set<size_t> search_key;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) search_key.insert( nit->Idx() );
                    // getting the node IDs in their correct order
                    array<size_t,3U> global_node_numbers;
                    uint32_t n_node{ 0U };
                    for ( const auto& nid : it->FE()->CornerNodesOfFace(face_id) )
                      global_node_numbers[n_node++] = it->N(nid)->Idx();
                    // inserting the face into the face vector
                    tria_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                 }
            }
        // prisms
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> search_key;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) search_key.insert( nit->Idx() );
                    // dealing with triangular faces
                    if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) ) {
                         // getting the node IDs in their correct order
                         array<size_t,3U> global_node_numbers;
                         uint32_t n_node{ 0U };
                         for ( const auto& nid : it->FE()->CornerNodesOfFace(face_id) )
                           global_node_numbers[n_node++] = it->N(nid)->Idx();
                         tria_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                      }
                 }
            }
        // pyramids fourth
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> search_key;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) search_key.insert( nit->Idx() );
                    if ( isTriangular(it->FE()->ElementTypeOfFace(face_id)) ) {
                         array<size_t,3U> global_node_numbers;
                         uint32_t n_node{ 0U };
                         for ( const auto& nid : it->FE()->CornerNodesOfFace(face_id) )
                           global_node_numbers[n_node++] = it->N(nid)->Idx();
                         tria_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                      }
                 }
            }
         // no octahedra in CSMP

        // ------------------------
        // quadrilateral faces
        // ------------------------
        // hexahedra
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> search_key;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) search_key.insert( nit->Idx() );
                    array<size_t,4U> global_node_numbers;
                    uint32_t n_node{ 0U };
                    for ( const auto& nid : it->FE()->CornerNodesOfFace(face_id) )
                      global_node_numbers[n_node++] = it->N(nid)->Idx();
                    quad_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                 }
            }
        // prisms
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> search_key;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) search_key.insert( nit->Idx() );
                    // dealing with quad and tria faces
                    if ( isQuadrilateral(it->FE()->ElementTypeOfFace(face_id)) ) {
                         array<size_t,4U> global_node_numbers;
                         uint32_t n_node{ 0U };
                         for ( const auto& nid : it->FE()->CornerNodesOfFace(face_id) )
                           global_node_numbers[n_node++] = it->N(nid)->Idx();
                         quad_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                      }
                 }
            }
        // pyramids fourth
        for ( const auto& it : model_domain.CellVector() )
          if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
               for ( uint32_t face_id{0U}; face_id<it->Faces(); ++face_id ) {
                    set<size_t> search_key;
                    for ( const auto& nit : it->CornerNodesOfFace(face_id) ) search_key.insert( nit->Idx() );
                    if ( isQuadrilateral(it->FE()->ElementTypeOfFace(face_id)) ) {
                         array<size_t,4U> global_node_numbers;
                         uint32_t n_node{ 0U };
                         for ( const auto& nid : it->FE()->CornerNodesOfFace(face_id) )
                           global_node_numbers[n_node++] = it->N(nid)->Idx();
                         quad_faces_.push_back( make_pair( search_key, global_node_numbers ) );
                      }
                 }
            }

      } // end 3D

    // sorting and eliminating duplicate faces from the faces vectors
    // --------------------------------------------------------------
    // (lambda function to search only the first pair in the pair-value pair)
    sort( tria_faces_.begin(), tria_faces_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
    // eliminating duplicate faces
    tria_faces_.erase( unique( tria_faces_.begin(), tria_faces_.end(), [](auto &left, auto &right) { return left.first == right.first; } ), tria_faces_.end() );
    sort( quad_faces_.begin(), quad_faces_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
    quad_faces_.erase( unique( quad_faces_.begin(), quad_faces_.end(), [](auto &left, auto &right) { return left.first == right.first; } ), quad_faces_.end() );


    // finding those faces that have a representation in the CSMP model
    // ----------------------------------------------------------------
    for ( auto it=model.Mesh().ElementsBegin(); it!=model.Mesh().ElementsEnd(); ++it )
      if ( (*it).IsSurface() ) {
           // create search key from sorted global node Idx
           set<size_t> global_node_numbers;
           // searching for face if it is triangular
           if ( isTriangular((*it).FE_Type()) ) {
                for ( uint32_t i{0U}; i<3U; ++i ) global_node_numbers.insert( (*it).N(i)->Idx() );
                // search face within tria_faces_ vector; recording its offset from origin in map
                const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                                  make_pair( global_node_numbers, array<size_t,3>{} ),
                                                  [](auto &left, auto &right) { return left.first < right.first; } );
                assert( face_it != tria_faces_.end() );
                // CSMP element idx, face # in UG
                csmp_faces_.insert( make_pair( (*it).Idx(), distance(tria_faces_.begin(), face_it) ) );
             }
           else { // quadrilateral
                for ( uint32_t i{0U}; i<4U; ++i ) global_node_numbers.insert( (*it).N(i)->Idx() );
                const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                                  make_pair( global_node_numbers, array<size_t,4>{} ),
                                                  [](auto &left, auto &right) { return left.first < right.first; } );
                assert( face_it != quad_faces_.end() );
                // CSMP element idx, face # in UG
                csmp_faces_.insert( make_pair( (*it).Idx(), distance(quad_faces_.begin(), face_it) + tria_faces_.size() ) );
             }
       }
       
    if constexpr( dim == 3U ) {
        for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit )
          if ( (*fit).IsSurface() ) {
               // create search key from sorted global node Idx
               set<size_t> global_node_numbers;
               // searching for face if it is triangular
               if ( isTriangular((*fit).FE_Type()) ) {
                    for ( uint32_t i{0U}; i<3U; ++i ) global_node_numbers.insert( (*fit).N(i)->Idx() );
                    // search face within tria_faces_ vector; if found its offset from origin is stored in map
                    const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                                      make_pair( global_node_numbers, array<size_t,3>{} ),
                                                      [](auto &left, auto &right) { return left.first < right.first; } );
                    assert( face_it != tria_faces_.end() );
                    // CSMP element idx, face # in UG
                    csmp_faces_.insert( make_pair( (*fit).Idx(), distance(tria_faces_.begin(), face_it) ) );
                 }
               else { // quadrilateral
                    for ( uint32_t i{0U}; i<4U; ++i ) global_node_numbers.insert( (*fit).N(i)->Idx() );
                    const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                                      make_pair( global_node_numbers, array<size_t,4>{} ),
                                                      [](auto &left, auto &right) { return left.first < right.first; } );
                    assert( face_it != quad_faces_.end() );
                    // CSMP element idx, face # in UG
                    csmp_faces_.insert( make_pair( (*fit).Idx(), distance(quad_faces_.begin(), face_it) + tria_faces_.size() ) );
                 }
           }
       } // end 3D

    // creating inverse map of csmp edges
    // ----------------------------------
    for ( const auto& it : csmp_faces_ )
      // face_id, csmp_idx
      csmp_faces_inverted_.insert( make_pair( it.second, it.first ) );

    // creating vector (0..faces-1) that indicates whether a face has a csmp equivalent
    // --------------------------------------------------------------------------------
    if ( faces_in_csmp_.empty() ) faces_in_csmp_.clear();
    faces_in_csmp_.resize( tria_faces_.size() + quad_faces_.size(), false );
    // initialising the vector
    for ( const auto& it : csmp_faces_ )
      faces_in_csmp_[ it.second ] = true;

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
    // assumption: the first set of nodes in the element are the corner nodes
    // tetrahedra first
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) ) {
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<4U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           tetra_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    sort( tetra_volumes_.begin(), tetra_volumes_.end(), [](auto &left, auto &right) { return left.first < right.first; } );
    
    // hexahedra
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isHexahedral( it->FE_Type() ) ) {
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<8U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           hexa_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    sort( hexa_volumes_.begin(), hexa_volumes_.end(), [](auto &left, auto &right) { return left.first < right.first; } );

    // prisms
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) ) {
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<6U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           prism_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    sort( prism_volumes_.begin(), prism_volumes_.end(), [](auto &left, auto &right) { return left.first < right.first; } );

    // pyramids
    for ( const auto& it : model_domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) ) {
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<5U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           pyra_volumes_.push_back( make_pair( global_node_numbers, it->Idx() ) );
        }
    sort( pyra_volumes_.begin(), pyra_volumes_.end(), [](auto &left, auto &right) { return left.first < right.first; } );

    // TODO: sorting to match order in volumes_?

    return tetra_volumes_.size() + hexa_volumes_.size() + prism_volumes_.size() + pyra_volumes_.size();
 
 } // end CollectVolumes









/**
      Method determines the position of the edges in the current region within the global 'edges' vector.
      
      @param region_edges is where the results are stored.
      @return number of unique edges in the region
*/
template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectEdgesInRegion( const Region<dim>& domain, vector<size_t>& region_edges ) const
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
             assert( edge_it != edges_.end() );
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
    if constexpr ( dim == 2U ) {
        vector<uint32_t> cnids;
        for ( const auto& it : domain.CellVector() )
          if ( it->IsSurface() && isTriangular( it->FE_Type() ) )
            {
                // for the elements in the region
                set<size_t> search_key;
                for ( uint32_t i{0U}; i<3U; ++i ) search_key.insert( it->N(i)->Idx() );
                // finding their position in the global 'faces_' vector
                const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                                  make_pair( search_key, array<size_t,3>{} ),
                                                  [](auto &left, auto &right) { return left.first < right.first; }  );
                if ( face_it != tria_faces_.end() )
                  region_faces.push_back( distance( tria_faces_.begin(), face_it ) );
             }
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
        // hexahedra (no triangular faces here)
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
    if constexpr( dim == 2U ) {
        for ( const auto& it : domain.CellVector() )
          if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) ) {
              set<size_t> search_key;
              for ( uint32_t i{0U}; i<4U; ++i ) search_key.insert( it->N(i)->Idx() );
              const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                                make_pair( search_key, array<size_t,4>{} ),
                                                [](auto &left, auto &right) { return left.first < right.first; }  );
              if ( face_it != quad_faces_.end() )
                 region_faces.push_back( distance( quad_faces_.begin(), face_it ) + tria_faces_.size() );
           }
       } // end 2D

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
        region_faces.shrink_to_fit();
      }

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
    // tetrahedral elements
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isTetrahedral( it->FE_Type() ) )
        {
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<4U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( tetra_volumes_.begin(), tetra_volumes_.end(), make_pair( global_node_numbers, it->Idx() ),
                                             [](auto &left, auto &right) { return left.first < right.first; } );
           if ( elmt_it != tetra_volumes_.end() )
             region_volumes.push_back( distance(tetra_volumes_.begin(),elmt_it) );
        }
    // hexahedral elements
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isHexahedral( it->FE_Type() ) )
        {
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<8U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( hexa_volumes_.begin(), hexa_volumes_.end(), make_pair( global_node_numbers, it->Idx() ),
                                             [](auto &left, auto &right) { return left.first < right.first; } );
           if ( elmt_it != hexa_volumes_.end() )
             region_volumes.push_back( distance(hexa_volumes_.begin(),elmt_it) + tetra_volumes_.size() );
        }
    // prism elements
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isPrism( it->FE_Type() ) )
        {
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<6U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( prism_volumes_.begin(), prism_volumes_.end(), make_pair( global_node_numbers, it->Idx() ),
                                             [](auto &left, auto &right) { return left.first < right.first; } );
           if ( elmt_it != prism_volumes_.end() )
             region_volumes.push_back( distance(prism_volumes_.begin(),elmt_it) + tetra_volumes_.size() +  hexa_volumes_.size() );
        }
    // pyramids
    for ( const auto& it : domain.CellVector() )
      if ( it->IsVolume() && isPyramid( it->FE_Type() ) )
        {
           // getting the global node IDs
           set<size_t> global_node_numbers;
           for ( uint32_t i{0U}; i<5U; ++i ) global_node_numbers.insert( it->N(i)->Idx() );
           // searching for the element IDx in the global volumes vector
           const auto elmt_it = lower_bound( pyra_volumes_.begin(), pyra_volumes_.end(), make_pair( global_node_numbers, it->Idx() ),
                                             [](auto &left, auto &right) { return left.first < right.first; } );
           if ( elmt_it != pyra_volumes_.end() )
             region_volumes.push_back( distance(pyra_volumes_.begin(),elmt_it) + tetra_volumes_.size() +  hexa_volumes_.size() + prism_volumes_.size() );
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


/* TO SEARCH FOR SPECIFIC EDGES
    if constexpr( dim == 2U ) {
        for ( const auto& it : domain.CellVector() )
          if ( it->IsLine() ) {
               // making a search key for the face from the global node IDs
               pair<size_t,size_t> search_key{ it->N(0U)->Idx(), it->N(1U)->Idx() };
               if ( search_key.first > search_key.second ) swap( search_key.first, search_key.second );
               // searching the face in the faces_ vector
               const auto edge_it = lower_bound( edges_.begin(), edges_.end(), search_key );
               if ( edge_it != edges_.end() )
               boundary_edges.push_back( distance( edges_.begin(), edge_it ) );
            }
      }
*/









template<uint32_t dim>
size_t UG4_UGX_FileExport<dim>::CollectFacesInBoundary( const Boundary<dim>& domain, std::vector<size_t>& boundary_faces ) const
 {
    if ( dim != 3U )
      throw csmp::Exception( ERROR, "UG4_UGX_FileExport<dim>::CollectFacesInBoundary",
                            "Boundary faces only exist in 3D models; only call method for three-dimensional models" );
 
    if ( !boundary_faces.empty() ) boundary_faces.clear();
    
    // boundary consists of triangular and quadrilateral faces these are unique
    // triangular faces first
    for ( const auto& it : domain.CellVector() )
      if ( it->IsSurface() && isTriangular( it->FE_Type() ) ) {
           // making a search key for the face from the global node IDs
           set<size_t> search_key;
           for ( uint32_t i{0U}; i<3U; ++i ) search_key.insert( it->N(i)->Idx() );
           // searching the face in the faces_ vector
           const auto face_it = lower_bound( tria_faces_.begin(), tria_faces_.end(),
                                             make_pair( search_key, array<size_t,3>{} ),
                                             [](auto &left, auto &right) { return left.first < right.first; } );
           if ( face_it != tria_faces_.end() )
           boundary_faces.push_back( distance( tria_faces_.begin(), face_it ) );
        }
      
    // quadrilaterals
    for ( const auto& it : domain.CellVector() )
      if ( it->IsSurface() && isQuadrilateral( it->FE_Type() ) ) {
           set<size_t> search_key;
           for ( uint32_t i{0U}; i<4U; ++i ) search_key.insert( it->N(i)->Idx() );
           // searching the face in the faces_ vector
           const auto face_it = lower_bound( quad_faces_.begin(), quad_faces_.end(),
                                             make_pair( search_key, array<size_t,4>{} ),
                                             [](auto &left, auto &right) { return left.first < right.first; } );
           if ( face_it != quad_faces_.end() )
           boundary_faces.push_back( distance( quad_faces_.begin(), face_it ) + tria_faces_.size() );
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



// --------------------------------------------------------------------------
//
//     NOT CURRENTLY USED
//
// --------------------------------------------------------------------------

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
          
          Element<dim>* eptr = &(*next(mesh.ElementsBegin(),idx));
          
          NOT USED 
 */
 template<uint32_t dim>
 static void writeNodeIndices( ofstream& ofs, const MeshManager<dim>& mesh, size_t cell_idx, bool& print_whitespace ) {
     if ( cell_idx < mesh.Elements() )
       for ( uint32_t i{0U}; i<(*next(mesh.ElementsBegin(),cell_idx)).Nodes(); ++i ) {
           if ( print_whitespace ) ofs <<" ";
           else print_whitespace = true;
           ofs << (*next(mesh.ElementsBegin(),cell_idx)).N(i)->Idx();
        }
     else if ( cell_idx < mesh.Elements()+mesh.Faces() )
       for ( uint32_t i{0U}; i<(*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Nodes(); ++i ) {
           if ( print_whitespace ) ofs <<" ";
           else print_whitespace = true;
           ofs << (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).N(i)->Idx();
        }
     else if ( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() )
       for ( uint32_t i{0U}; i<(*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Nodes(); ++i ) {
           if ( print_whitespace ) ofs <<" ";
           else print_whitespace = true;
           ofs << (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).N(i)->Idx();
        }
     else {
       if ( !equivalentEntityInCSMP(cell_idx) )
         ErrorHandler::Instance().Note( ERROR, "writeNodeIndices", "cell index greater than size_t::max");
       else
         ErrorHandler::Instance().Note( ERROR, "writeNodeIndices", "cell index out of range");
       }
       
  } // end writeNodeIndices



/**
    For cell-based variables (Element, Face, InterFace), writes its value, considering SCALAR, VECTOR and TENSOR variables and printing no_data values as appropriate
    
    @param ofs  output textfile stream UTF-8 on linux
    @param var_key variable accessor in CSMP
    @param mesh MeshManager which contains the variable values of interest stored on Element, Face or InterFace, respectively.
    @param cell_idx continous numbering of elements, faces and interfaces in CSMP
    @param print_whitespace to avoid printing blanks at the  beginning of attachment arrays
    
      @attention method has to switch print_whitespace on
*/
template<uint32_t dim>
static void writeCellVariableValue( ofstream& ofs, const csmp::Index& var_key,
                                    const MeshManager<dim>& mesh, size_t cell_idx, bool& print_whitespace )
 {
    csmp::ErrorHandler& csmp_err( ErrorHandler::Instance() );

    if ( equivalentEntityInCSMP(cell_idx) ) {
        // there no cell index was recorded for the face the no data values is writen
        switch(var_key.type) {
            case SCALAR: {
                   double sc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             sc = (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key);
                          break;
                        case FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces() );
                             sc = (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Read(var_key);
                          break;
                        case INTER_FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() );
                             sc = (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Read(var_key);
                          break;
                        default:
                             sc = -numeric_limits<float>::max(); // no-data value
                     }
                   if ( !print_whitespace ) ofs << scientific << sc;
                   else ofs <<" "<< scientific << sc;
                   print_whitespace = true;
                }
              break;
            case VECTOR: {
                   VectorVariable<dim> vc;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key,vc);
                          break;
                        case FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces() );
                             (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Read(var_key,vc);
                          break;
                        case INTER_FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() );
                             (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Read(var_key,vc);
                          break;
                        default:
                             vc = -numeric_limits<float>::max(); // no-data value
                     }
                   for ( uint32_t i{0u}; i<dim; ++i ) {
                         if ( print_whitespace ) ofs <<" ";
                         else print_whitespace = true;
                         ofs << scientific << vc[i];
                     }
                }
              break;
            case TENSOR: {
                   TensorVariable<dim> ts;
                   switch(var_key.place) {
                        case ELEMENT:
                             assert( cell_idx < mesh.Elements() );
                             (*next(mesh.ElementsBegin(),cell_idx)).Read(var_key,ts);
                          break;
                        case FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces() );
                             (*next(mesh.FacesBegin(),cell_idx-mesh.Elements())).Read(var_key,ts);
                          break;
                        case INTER_FACE:
                             assert( cell_idx < mesh.Elements()+mesh.Faces()+mesh.InterFaces() );
                             (*next(mesh.InterFacesBegin(),cell_idx-mesh.Elements()-mesh.Faces())).Read(var_key,ts);
                          break;
                        default:
                             ts = -numeric_limits<float>::max(); // no-data value
                     }
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


template class UG4_UGX_FileExport<2U>;
template class UG4_UGX_FileExport<3U>;

} // end csmp
