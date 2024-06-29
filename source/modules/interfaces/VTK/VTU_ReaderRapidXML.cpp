//
//  XML_VTU_ReaderRapidXML.cpp
//
//  Created by Stephan Matthai on 26/1/2024.
//

#include "VTU_ReaderRapidXML.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <type_traits>
#
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

#include "CSMP_global_enumerations.h"
#include "VTK_Type.h"
#include "VSet.h"
#include "CSMP_ElementSpecifications.h"
#include "ModelTopology.h" // to form regions from 'region id' code in VTU file
#include "ErrorHandler.h"

// HELPER FUNCTIONS

/// for printing the different data stored in the VTU file
static void iterateSiblings(rapidxml::xml_node<>* node) {
    // Iterate through all sibling nodes, including the current node (parent)
    for (rapidxml::xml_node<>* sibling = node; sibling; sibling = sibling->next_sibling()) {
        // Process each sibling node here
        std::cout <<"\n"<<"XML-node name: " << sibling->name() << std::endl;
        // If you need to process attributes of each sibling, you can iterate through them similarly
        for (rapidxml::xml_attribute<>* attribute = sibling->first_attribute(); attribute; attribute = attribute->next_attribute()) {
            std::cout << "\tAttribute name: "<< attribute->name() <<", Value: '"<< attribute->value() <<"'"<< std::endl;
        }
    }
}

/// searching for a sibling dataset stored in the VTU file
static bool findSiblingNamed( rapidxml::xml_node<>* node, const char* name ) {
    std::string variable_name = name;
    // looping over xml DataArray nodes
    for (rapidxml::xml_node<>* sibling = node; sibling; sibling = sibling->next_sibling()) {
         // looping over their attributes searching the name of the dataset
         for (rapidxml::xml_attribute<>* attribute = sibling->first_attribute(); attribute; attribute = attribute->next_attribute())
           if ( variable_name == attribute->value() ) return true;
      }

   return false;
}

static rapidxml::xml_node<>* getSiblingNode( rapidxml::xml_node<>* node, const char* name ) {
    std::string variable_name = name;
    for (rapidxml::xml_node<>* sibling = node; sibling; sibling = sibling->next_sibling()) {
         for (rapidxml::xml_attribute<>* attribute = sibling->first_attribute(); attribute; attribute = attribute->next_attribute())
           if ( variable_name == attribute->value() ) return sibling;
      }

   return nullptr;
}

// CORE FUNCTION
using namespace std;

namespace csmp {

/**
   using VTK_TYPE = int;
*/
template<uint32_t dim>
int read_VTU_File( const char* fname, VSet<dim>& vset, ModelTopology& topology, bool verbose )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    try {
        // 1. Open and read the VTU-XML file into string-stream buffer
        // -----------------------------------------------------------
        // input files successfully tested:
        // string file_name{"blunt30deg_region_id.vtu"};
        // string file_name{"Greenshank_prism_mini-input_data_flow_domain.vtu"}; // OK
        // string file_name{"GFV3x3x03_PS2_PillarGridding_ProportionalLayering_quarter_cc_rt.vtu"};
        string file_name{ fname };
        file_name +=".vtu";
        ifstream    file(file_name.c_str());
        if (!file.is_open()) {
            csmp_error.Note( ERROR, "read_VTU_File:", file_name, "Unable to open XML file with the extension '.vtu'.");
            return 1;
        }

        // Read entire XML file into a string
        cout <<"\n"<<"readVTU_File: reading data from VTU file '"<< file_name <<"'"<<" in ASCII format."<< endl;
        string xmlContent((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
        file.close();
        
        // Parse the XML content using RapidXML
        rapidxml::xml_document<> doc;
        doc.parse<0>(&xmlContent[0]);

        // Access the root node called 'VTKFile'
        rapidxml::xml_node<>* vtkFileNode = doc.first_node("VTKFile");
        if (!vtkFileNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "Root node 'VTKFile' not found in file.");
            return 2;
        }

        // Access 'UnstructuredGrid' node
        rapidxml::xml_node<>* unstructuredGridNode = vtkFileNode->first_node("UnstructuredGrid");
        if (!unstructuredGridNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'UnstructuredGrid' node (record) not found.");
            return 3;
        }

        // 2. Access 'Piece' with vertex and cell data
        // -------------------------------------------------------------
        rapidxml::xml_node<>* pieceNode = unstructuredGridNode->first_node("Piece");
        if (!pieceNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'Piece' node (record) not found.");
            return 4;
        }
        // getting the number of node points and cells = finite elements
        const size_t n_vertices = atol( pieceNode->first_attribute("NumberOfPoints")->value() );
        const size_t n_cells    = atol( pieceNode->first_attribute("NumberOfCells")->value() );

        // 2.1 Accessing node 'Points' coordinates
        // ---------------------------------------
        rapidxml::xml_node<>* pointsNode = pieceNode->first_node("Points");
        if (!pointsNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'Points' node (record) not found.");
            return 5;
        }

        // Access DataArray with the node Points
        rapidxml::xml_node<>* dataArrayPointsNode = pointsNode->first_node("DataArray");
        if (!dataArrayPointsNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'DataArray' node (record) not found for point data.");
            return 6;
        }

        // Extracting the data from DataArray for the (vertex) Points
        deque<double> px, py, pz;
        double        coordinateValue[3];
        istringstream ss(dataArrayPointsNode->value());
        size_t counter{0u};
        while (ss >> coordinateValue[counter++] ) {
            if ( counter == 3U ) {
                 px.push_back( coordinateValue[0] );
                 py.push_back( coordinateValue[1] );
                 pz.push_back( coordinateValue[2] );  // Assuming 3 components with x, y, z
                 // resetting counter for next point
                 counter = 0U;
              }
        }
        if ( px.size() != n_vertices || py.size() != n_vertices || pz.size() != n_vertices ) {
             cout <<"\n\t"<<"n_vertices= "<< n_vertices <<" vs. "<< px.size() <<" vertices read."<< endl;
             csmp_error.Note( ERROR, "readVTU_File", "number of node points read does not match the number in the node attribute specification; trimming dequeues.");
             px.resize( n_vertices );
             py.resize( n_vertices );
             pz.resize( n_vertices );
          }

        // Print tge extracted vertex coordinates
        if ( verbose ) {
            cout <<"\n"<<"Extracted Point Coordinates:" << endl;
            for ( size_t i{0u}; i<px.size(); ++i ) {
                cout << "(" << px[i] << ", " << py[i] << ", " << pz[i] << ")" << endl;
            }
            cout << endl;
        }
        
        // 2.1 Accessing the 'Cells' / finite-element specifications
        // ---------------------------------------------------------
        // Access Cells node
        rapidxml::xml_node<>* cellsNode = pieceNode->first_node("Cells");
        if (!cellsNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'Cells' node not found.");
            return 7;
        }

        // Access DataArray node for connectivity
        rapidxml::xml_node<>* dataArrayConnectivityNode = cellsNode->first_node("DataArray");
        if (!dataArrayConnectivityNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'DataArray' node (record) not found for 'connectivity' data.");
            return 8;
        }

        // 2.1.1 Extracting data from DataArray 1 for connectivity
        // -------------------------------------------------------
        vector<vector<uint32_t>> connectivityValues;
        deque<int64_t> plist_data;
        {
          // Accessing attributes
          const char* typeAttrConnectivity = dataArrayConnectivityNode->first_attribute("type")->value();
          const char* nameAttrConnectivity = dataArrayConnectivityNode->first_attribute("Name")->value();
          const char* numComponentsAttrConnectivity = ( dataArrayConnectivityNode->first_attribute("NumberOfComponents") ) ?
                                                        dataArrayConnectivityNode->first_attribute("NumberOfComponents")->value() : nullptr;
          const char* formatAttrConnectivity = dataArrayConnectivityNode->first_attribute("format")->value();

          if ( verbose ) {
              cout << "Attributes for DataArray (connectivity='plist'):" << endl;
              cout << "type: " << typeAttrConnectivity << endl;
              cout << "Name: " << nameAttrConnectivity << endl;
              if ( numComponentsAttrConnectivity ) cout << "NumberOfComponents: " << numComponentsAttrConnectivity << endl;
              cout << "format: " << formatAttrConnectivity << endl;
            }

          // reading the nodes-per-element record to create 'plist'
          istringstream ssConnectivity(dataArrayConnectivityNode->value());
          int64_t valueNodeID;
          while (ssConnectivity >> valueNodeID) {
              plist_data.push_back(valueNodeID);
          }

          // Print extracted data
          if ( verbose ) {
              cout <<"\n"<<"Extracted Connectivity Values:" << endl;
              for (const auto& p : plist_data)
                cout << p << " ";
              cout << endl;
          }
        }

        // 2.1.2 Extracting data from DataArray 2 for offsets
        // --------------------------------------------------
        // Access DataArray node for offsets (beginnings of nodes of the next elements)
        rapidxml::xml_node<>* dataArrayOffsetsNode = cellsNode->first_node("DataArray")->next_sibling("DataArray");
        if (!dataArrayOffsetsNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'DataArray' node (record) not found for 'offsets' data.");
            return 9;
        }

        // Accessing attributes for offsets
        vector<uint32_t> offsetsValues;
        {
          const char* typeAttrOffsets = dataArrayOffsetsNode->first_attribute("type")->value();
          const char* nameAttrOffsets = dataArrayOffsetsNode->first_attribute("Name")->value(); // 'offsets'
          const char* numComponentsAttrOffsets = ( dataArrayOffsetsNode->first_attribute("NumberOfComponents") ) ?
                                                   dataArrayOffsetsNode->first_attribute("NumberOfComponents")->value() : nullptr;
          const char* formatAttrOffsets = dataArrayOffsetsNode->first_attribute("format")->value();

          if ( verbose ) {
              cout <<"\n"<<"Attributes for DataArray (offsets):" << endl;
              cout <<"\t"<< "type: " << typeAttrOffsets << endl;
              cout <<"\t"<< "Name: " << nameAttrOffsets << endl;
              if ( numComponentsAttrOffsets ) cout <<"\t"<< "NumberOfComponents: " << numComponentsAttrOffsets << endl;
              cout <<"\t"<< "format: " << formatAttrOffsets << endl;
            }

          // Extracting offsets data
          istringstream ssOffsets(dataArrayOffsetsNode->value());
          uint32_t valueOffsets;
          while (ssOffsets >> valueOffsets) {
              offsetsValues.push_back(valueOffsets);
          }

          if ( offsetsValues.size() != n_cells ) {
               cout <<"\n\t"<<"n_cells = "<< n_cells <<" vs. "<< offsetsValues.size() <<" cell 'offset' data read."<< endl;
               csmp_error.Note( ERROR, "readVTU_File", "number of offset records does not match the number of cells in the model");
            }
        }
        
        // Print extracted data for offsets
        if ( verbose ) {
            cout <<"\n"<<"Extracted Offsets Values:" << endl;
            for (const auto& value : offsetsValues) {
                cout << value << " ";
            }
            cout << endl;
        }

        // creating the 'npes' (number of nodes-per-element) and 'plist' (nodes-per-element) records from the plist_data and offset data
        deque<vector<size_t> >  plist;
        deque<uint32_t>         n_nodes_per_element;
        size_t                  prev_offset{0u};
        counter = 0u;
        for ( const auto& offset : offsetsValues ) {
             // nodes per element record
             n_nodes_per_element.push_back( static_cast<uint32_t>(offset - prev_offset) );
             // plist record
             vector<size_t> plist_entry;
             plist_entry.reserve(8u); // optimal for lin.-hex.
             for ( size_t i{prev_offset}; i<offset; ++i )
               plist_entry.push_back( plist_data[counter++] );
             plist.push_back( plist_entry );
             prev_offset = offset;
          }
        assert( counter == plist_data.size() );
        assert( plist.size() == n_cells );


        // 2.1.3 Extracting data from DataArray 3 for cell types
        // -----------------------------------------------------
        // Access DataArray node for offsets
        rapidxml::xml_node<>* dataArrayCellDataNode = cellsNode->first_node("DataArray")->next_sibling("DataArray")->next_sibling("DataArray");
        if (!dataArrayCellDataNode) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'DataArray' node (record) not found for cell 'types' data.");
            return 10;
        }

        // Accessing attributes
        deque<int8_t> valueCellTypes;
        {
          const char* typeAttrCellData = dataArrayCellDataNode->first_attribute("type")->value();
          const char* nameAttrCellData = dataArrayCellDataNode->first_attribute("Name")->value();
          const char* numComponentsAttrCellData = ( dataArrayCellDataNode->first_attribute("NumberOfComponents") ) ?
                                                   dataArrayCellDataNode->first_attribute("NumberOfComponents")->value() : nullptr;
          const char* formatAttrCellData = dataArrayCellDataNode->first_attribute("format")->value();

          if ( verbose ) {
               cout <<"\n"<<"Attributes for DataArray (types):" << endl;
               cout <<"\t"<<"type: " << typeAttrCellData << endl;
               cout <<"\t"<<"Name: " << nameAttrCellData << endl;
               if ( numComponentsAttrCellData ) cout <<"\t"<<"NumberOfComponents: " << numComponentsAttrCellData << endl;
               cout <<"\t"<<"format: " << formatAttrCellData << endl;
             }

          istringstream ssCellData(dataArrayCellDataNode->value());

          int valueCellType;
          while (ssCellData >> valueCellType) {
              valueCellTypes.push_back( csmpTypeFromVTK_TYPE( static_cast<VTK_TYPE>(valueCellType) ) );
          }

          // Print extracted finite-element type data
          if ( verbose ) {
               cout <<"\n"<<"Extracted Cell types:" << endl;
               for (const auto& value : valueCellTypes)
                 cout << parseFiniteElementType(value) << " ";
               cout << endl;
            }

          if ( valueCellTypes.size() != n_cells ) {
               cout <<"\n\t"<<"n_cells = "<< n_cells <<" vs. "<< valueCellTypes.size() <<" cells read."<< endl;
               csmp_error.Note( ERROR, "readVTU_File", "number of finite element type specifications does not match the number of finite elements in the model");
            }
        }
        plist_data.clear();
        offsetsValues.clear();


        // 3. Setting up the VSET to store the collected information
        // ---------------------------------------------------------
        cout <<"\n"<<"readVTU_File: setting up a csmp::VSet with the polygonal mesh data from the VTU file."<< endl;
        // getting number of neighbor elements per element from CSMP
        deque<uint32_t> n_neighbors_per_element; // no pfverts values, they will be created by CSMP elmt_specs;
        for ( const auto& cit : valueCellTypes )
          n_neighbors_per_element.push_back( CSMP_ElementSpecifications::NeighborsPerElementOfType( cit ) );
        
        // Configuring the VSET
        vset.Resize( valueCellTypes, n_nodes_per_element, n_neighbors_per_element, n_vertices, 0U, 0U );
        n_nodes_per_element.clear();
        n_neighbors_per_element.clear();

        // node coordinates from above
        vset.AddXYZ( px, py, pz );
        px.clear();
        py.clear();
        pz.clear();
        
        // finite-element / cell types
        vset.AddElementTypes( valueCellTypes.begin(), valueCellTypes.end() );
        valueCellTypes.clear();
        
        // nodes per element record
        vset.AddPlist( plist.begin(), plist.end() );
        plist.clear();
        
        // creating element neighbor information
        vset.EstablishElementConnectivity3D();
        const uint32_t MODEL_DIMENSION = vset.SpatialDimension();
        if ( MODEL_DIMENSION != dim )
          csmp_error.Note( INFO, "readVTU_File", to_string(MODEL_DIMENSION), "= spatial dimension of the model in the file is different from that expected by the reader.");


        // 4. creating regions stored in the model topology class
        // ------------------------------------------------------
        cout <<"\n\n"<<"readVTU_File: attempting to build a csmp::ModelTopology object, provided that 'region id' or 'attribute' information is contained in the VTU file."<< endl;

        // Access Cells node
        rapidxml::xml_node<>* cellData = pieceNode->first_node("CellData");
        if (!cellData) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'CellArray' node (record) not found.");
            return 11;
        }
        // Access DataArray node for connectivity
        rapidxml::xml_node<>* cellDataArray = cellData->first_node("DataArray");
        if (!cellDataArray) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'DataArray' node (record) not found for 'CellData'.");
            return 12;
        }
     
        // search for 'region id' or 'attribute' data inside of CellData to create CSMP regions from
        rapidxml::xml_node<>* topologyDataArray{ nullptr };

        if ( !findSiblingNamed( cellDataArray, "attribute" ) )
          csmp_error.Note( INFO, "readVTU_File", "the VTU file does not contain GEOS 'attribute' data to create ModelTopology (regions) from.");
        else topologyDataArray = getSiblingNode( cellDataArray, "attribute" );

        if ( !findSiblingNamed( cellDataArray, "region id" ) )
          csmp_error.Note( INFO, "readVTU_File", "the VTU file does not contain 'region id' to create ModelTopology (regions) from");
        else topologyDataArray = getSiblingNode( cellDataArray, "region id" );
 
        if ( !topologyDataArray )
          csmp_error.Note( WARNING, "readVTU_File", "there are no data in the VTU file to create csmp::Region objects from; step will be skipped and there will only be one region called 'Model'.");

         // Accessing attributes
        if ( topologyDataArray ) {
              cout <<"\n"<<"VTU file contains 'attribute' or 'region id' data that will be used to initialise csmp::ModelTopology..." << endl;
              const char* typeAttrCellPropertyData = topologyDataArray->first_attribute("type")->value();
              const char* nameAttrCellPropertyData = topologyDataArray->first_attribute("Name")->value();
              const char* numComponentsAttrCellPropertyData = ( topologyDataArray->first_attribute("NumberOfComponents") ) ?
                                                               topologyDataArray->first_attribute("NumberOfComponents")->value() : nullptr;
              const char* formatAttrCellPropertyData = topologyDataArray->first_attribute("format")->value();

              if ( verbose ) {
                   cout <<"\n"<<"Attributes for DataArray (property data):" << endl;
                   cout <<"\t"<<"type: " << typeAttrCellPropertyData << endl;
                   cout <<"\t"<<"Name: " << nameAttrCellPropertyData << endl;
                   if ( numComponentsAttrCellPropertyData ) cout <<"\t"<<"NumberOfComponents: " << numComponentsAttrCellPropertyData << endl;
                   cout <<"\t"<<"format: " << formatAttrCellPropertyData << endl;
                }

              istringstream ssCellPropertyData(topologyDataArray->value());

              // storing which elements belong to which unique regions
              long valueCellRegionProperty;
              map<long,set<size_t>> model_regions;
              counter = 0u;
              while (ssCellPropertyData >> valueCellRegionProperty) {
                  //                                              'region id' or 'attribute'           'cell id'
                  auto region_it = model_regions.insert( make_pair(valueCellRegionProperty, set<size_t>{counter} ) );
                  // if a region with this ID already exists, the cell is added to it
                  if ( region_it.second != true )
                    (*region_it.first).second.insert(counter);
                  counter++;
               }
              if ( model_regions.size() == 1U )
                cout <<"\n"<<"Identified 1 unique regions by its 'region id'"<< endl;
              else
                cout <<"\n"<<"Identified "<< model_regions.size() <<" unique regions by their '"<< nameAttrCellPropertyData <<"'"<< endl;

              // adding the regions to the model topology
              if ( !model_regions.empty() ) {
                    vector<size_t> region_cells;
                    for ( const auto& rit : model_regions ) {
                         region_cells.assign( rit.second.begin(), rit.second.end() );
                         // identifying the element type of the region assuming that they all consist of a single type of cell type (VTK_TYPE)
                         const int8_t fe_type = ( vset.HybridElementTypeMesh() == true ) ? vset.ElementType( region_cells[0u] ) : vset.ElementType(0u);
                         topology.AddDomain( (string("region") + to_string(rit.first)).c_str(), set<string>{parseFiniteElementType(fe_type)}, region_cells );
                         region_cells.clear();
                     }
                 }
              assert( topology.Cells() == vset.Cells() );
           } // end model topology
           

        // 5. Adding the rocktypes, material data and regions to VSet
        // ----------------------------------------------------------
        // Iterate through all attributes of the parent node, count and read them except for the 'region id' data (attribute in GEOS) that we already collected
        if ( verbose ) {
             cout <<"\n\n\n"<<"readVTU_File: CELL data stored in VTU file:";
             iterateSiblings( cellDataArray );
          }
        
                // Access PointData node
        rapidxml::xml_node<>* pointData = pieceNode->first_node("PointData");
        if (!pointData) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'PointData' node (record) not found.");
            return 13;
        }
        // Access point DataArray objects
        rapidxml::xml_node<>* pointDataArray = pointData->first_node("DataArray");
        if (!pointDataArray) {
            csmp_error.Note( ERROR, "read_VTU_File:", "'DataArray' node (record) not found for 'PointData'.");
            return 14;
        }
        if ( verbose ) {
             cout <<"\n\n\n"<<"readVTU_File: POINT data stored in VTU file:";
             iterateSiblings( pointDataArray );
          }

        // CELL PROPERTIES: reading property data placed on the cell
        cout <<"\n\n"<<"readVTU_File: reading VTK CELL-based data from the VTU file..." << endl;
        for (rapidxml::xml_node<>* sibling = cellDataArray; sibling; sibling = sibling->next_sibling()) {
            // Diagnostic attributes
            // - name of property
            const char* nameAttrCellPropertyData = sibling->first_attribute("Name")->value(); // property value
            // skipping property 'region id' since it was already read
            if ( string(nameAttrCellPropertyData) == "attribute" || string(nameAttrCellPropertyData) == "region id" ) continue;
            // - variable type (scalar, vector, tensor, array..
            const char* numComponentsAttrCellPropertyData = ( sibling->first_attribute("NumberOfComponents") ) ?
                                                             sibling->first_attribute("NumberOfComponents")->value() : nullptr;
            // reading the cell property data
            istringstream ssCellPropertyData(sibling->value());
            const int numberOfComponents = ( numComponentsAttrCellPropertyData ) ? atoi(numComponentsAttrCellPropertyData) : 1;
            const csmp::VARIABLE_TYPE var_type = inferTypeFromNumberOfComponents<dim>( numberOfComponents );
            double valueCellProperty;
            counter = 0U;
            
            // creating the appropriate datastructure to store them
            cout <<"\n\t\t"<<"reading CELL property: '"<< nameAttrCellPropertyData <<"' and adding it to the VSet."<< endl;
            switch( var_type ) {
                 case SCALAR: {
                        PropertyData elmt_prop( ELEMENT, SCALAR, dim );
                        elmt_prop.Reserve( vset.Cells() );
                        while (ssCellPropertyData >> valueCellProperty)
                          pushBack( elmt_prop, makeScalar( ANY, valueCellProperty ) );
                        vset.AddData( nameAttrCellPropertyData, elmt_prop );
                        // turning rocktype identifiers into PMTRL entries
                        if ( string("rocktype") == nameAttrCellPropertyData ) {
                             vector<int32_t> pmtrl;
                             pmtrl.reserve( n_cells );
                             for ( size_t i{0u}; i<n_cells; ++i )
                               // back from double to integer again (not ideal!)
                               pmtrl.push_back( static_cast<int32_t>(elmt_prop.Value(i)) );
                             vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
                             cout <<"\n\t\t"<<"found '"<< nameAttrCellPropertyData <<"' information and converted it into 'pmtrl' record in VSet."<< endl;
                          }
                     }
                   break;
                 case VECTOR: { // TODO: not tested yet
                        PropertyData elmt_prop( ELEMENT, VECTOR, dim );
                        elmt_prop.Reserve( vset.Cells() );
                        vector<double> value( dim );
                        while ( ssCellPropertyData >> value[counter++] ) {
                             if ( counter == dim ) {
                                  VectorVariable<dim> vc( value );
                                  pushBack( elmt_prop, vc );
                                  counter = 0U;
                               }
                          }
                        vset.AddData( nameAttrCellPropertyData, elmt_prop );
                     }
                   break;
                 case TENSOR: { // TODO: not tested yet
                        PropertyData elmt_prop( ELEMENT, TENSOR, dim );
                        elmt_prop.Reserve( vset.Cells() );
                        vector<double>      value( dim * dim );
                        TensorVariable<dim> ts;
                        while ( ssCellPropertyData >> value[counter++] ) {
                             if ( counter == dim * dim ) {
                                  uint32_t counter2{0u};
                                  for ( uint32_t i{0u}; i<dim; ++i )
                                    for ( uint32_t j{0u}; j<dim; ++j )
                                      ts(i,j) = value[counter2++];
                                  pushBack( elmt_prop, ts );
                                  counter = 0U;
                               }
                          }
                        vset.AddData( nameAttrCellPropertyData, elmt_prop );
                     }
                   break;
                 case ARRAY: {
                        PropertyData elmt_prop( ELEMENT, ARRAY, dim, numberOfComponents );
                        elmt_prop.Reserve( vset.Cells() );
                        vector<double> value( numberOfComponents );
                        ArrayVariable array( numberOfComponents );
                        while ( ssCellPropertyData >> array( static_cast<uint32_t>(counter++) ) ) {
                             if ( counter == numberOfComponents ) {
                                  pushBack( elmt_prop, array );
                                  counter = 0U;
                               }
                          }
                        vset.AddData( nameAttrCellPropertyData, elmt_prop );
                     }
                   break;
                 default:
                   cerr <<"\nReading cell property data: flagged-array or other variables not handled yet." << endl;
              }
        } // end cell properties



        // POINT PROPERTIES: reading property data placed on the cell
        cout <<"\n\n"<<"readVTU_File: reading VTK POINT-based data from the VTU file..." << endl;
        for (rapidxml::xml_node<>* sibling = pointDataArray; sibling; sibling = sibling->next_sibling()) {
            // Diagnostic attributes
            // - name of property
            const char* nameAttrPointPropertyData = sibling->first_attribute("Name")->value(); // property value
            // skipping property 'region id' since it was already read
            if ( string(nameAttrPointPropertyData) == "attribute" || string(nameAttrPointPropertyData) == "region id"  ) continue;
            // - variable type (scalar, vector, tensor, array.. (1=scalar, dim=vector, dim*dim=tensor..)
            const char* numComponentsAttrPointPropertyData = ( sibling->first_attribute("NumberOfComponents") ) ?
                                                              sibling->first_attribute("NumberOfComponents")->value() : nullptr;
            // reading the cell property data
            istringstream ssPointPropertyData(sibling->value());
            const int numberOfComponents = ( numComponentsAttrPointPropertyData ) ? atoi(numComponentsAttrPointPropertyData) : 1;
            const csmp::VARIABLE_TYPE var_type = inferTypeFromNumberOfComponents<dim>( numberOfComponents );
            double valuePointProperty;
            counter = 0U;
            
            // creating the appropriate datastructure to store them
            cout <<"\n\t\t"<<"reading POINT property: '"<< nameAttrPointPropertyData <<"' and adding it to the VSet."<< endl;
            switch( var_type ) {
                 case SCALAR: {
                        PropertyData node_prop( NODE, SCALAR, dim );
                        node_prop.Reserve( vset.Vertices() );
                        while (ssPointPropertyData >> valuePointProperty)
                          pushBack( node_prop, makeScalar( ANY, valuePointProperty ) );
                        vset.AddData( nameAttrPointPropertyData, node_prop );
                     }
                   break;
                 case VECTOR: { // TODO: not tested yet
                        PropertyData node_prop( NODE, VECTOR, dim );
                        node_prop.Reserve( vset.Vertices() );
                        vector<double> value( dim );
                        while ( ssPointPropertyData >> value[counter++] ) {
                             if ( counter == dim ) {
                                  VectorVariable<dim> vc( value );
                                  pushBack( node_prop, vc );
                                  counter = 0U;
                               }
                          }
                        vset.AddData( nameAttrPointPropertyData, node_prop );
                     }
                   break;
                 case TENSOR: { // TODO: not tested yet
                        PropertyData node_prop( NODE, TENSOR, dim );
                        node_prop.Reserve( vset.Vertices() );
                        vector<double>      value( dim * dim );
                        TensorVariable<dim> ts;
                        while ( ssPointPropertyData >> value[counter++] ) {
                             if ( counter == dim * dim ) {
                                  uint32_t counter2{0u};
                                  for ( uint32_t i{0u}; i<dim; ++i )
                                    for ( uint32_t j{0u}; j<dim; ++j )
                                      ts(i,j) = value[counter2++];
                                  pushBack( node_prop, ts );
                                  counter = 0U;
                               }
                          }
                        vset.AddData( nameAttrPointPropertyData, node_prop );
                     }
                   break;
                 case ARRAY: { // TODO: not tested yet
                        PropertyData node_prop( NODE, ARRAY, dim, numberOfComponents );
                        node_prop.Reserve( vset.Vertices() );
                        vector<double> value( numberOfComponents );
                        ArrayVariable array( numberOfComponents );
                        while ( ssPointPropertyData >> array( static_cast<uint32_t>(counter++) ) ) {
                             if ( counter == numberOfComponents ) {
                                   pushBack( node_prop, array );
                                  counter = 0U;
                               }
                          }
                        vset.AddData( nameAttrPointPropertyData, node_prop );
                     }
                   break;
                 default:
                   cerr <<"\nReading point properties (node data): flagged-array or other variables not handled yet." << endl;
              }
        } // end reading point properties

    cout <<"\n\n"<<"read_VTU_File: finished reading VTU file and extracted data into VSet successfully." << endl;


    } catch (const exception& e) {
        cerr << "\n\n"<<"read_VTU_File: Exception caught: " << e.what() << endl;
        return 15;
    }

    return 0;
    
} // end readVTU_File


template int read_VTU_File( const char* fname, VSet<2U>&, ModelTopology&,  bool );
template int read_VTU_File( const char* fname, VSet<3U>&, ModelTopology&,  bool );


} // end csmp
