#include <sstream>
#include "VTK_Interface.h"
#include "Element.h"
#include "Region.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_highLevelUtilities.h"
#include "OS_Utilities.h"

// revision history:
// SKM 1/7/2005 - replaced inefficient element/face list by vectors
//              - added constraint point output

using namespace std;

namespace csmp {

template<size_t dim>
VTK_Interface<dim>::VTK_Interface( const std::string& problemTitle,
                                   bool use_propblem_title_as_output_folder_name )
 : last_visualized_(static_cast<PLACEMENT>(UNSPECIFIED)), 
   node_output_of_element_data_(false),
   problemTitle_ ( problemTitle ),
   toFolder_     ( use_propblem_title_as_output_folder_name ),
   toSubFolder_  ( false )
 {

 }

template<size_t dim>
VTK_Interface<dim>::VTK_Interface( const std::string& problemTitle,
                                   const std::string& subFolderName,
                                   bool use_propblem_title_as_output_folder_name )
 : last_visualized_(static_cast<PLACEMENT>(UNSPECIFIED)),
   node_output_of_element_data_(false),
   problemTitle_ ( problemTitle ),
   subFolderName_( subFolderName ),
   toFolder_     ( use_propblem_title_as_output_folder_name ),
   toSubFolder_  ( true )
 {

 }

template<size_t dim>
VTK_Interface<dim>::~VTK_Interface()
 {

 }

template<size_t dim>
bool VTK_Interface<dim>::NodeOutputOfElementData() const 
 {
    return node_output_of_element_data_;
 }


template<size_t dim>
void VTK_Interface<dim>::NodeOutputOfElementData( bool yes_or_no )
 {
    node_output_of_element_data_ = yes_or_no;
 }

template<size_t dim>
std::string VTK_Interface<dim>::OutputFileAndSubFolderName( const std::string& fileName )
{
    if( toSubFolder_ || toFolder_ ){
        std::string output_file_name;
        output_file_name  = currentDirectorySymbol();
        if( toFolder_ ){
            output_file_name += problemTitle_;
            createDirectoryIfDoesntExist( output_file_name );
            output_file_name += directorySymbol();
            if( toSubFolder_ ){
                output_file_name += subFolderName_;
                createDirectoryIfDoesntExist( output_file_name );
                output_file_name += directorySymbol();
            }
        }else{
            output_file_name += subFolderName_;
            createDirectoryIfDoesntExist( output_file_name );
            output_file_name += directorySymbol();
        }
        if( fileName.find( output_file_name ) == std::string::npos ){
           output_file_name += fileName;
           return output_file_name;
        }
    }
    return fileName;
}


template<size_t dim>
template<class T>
void VTK_Interface<dim>::OutputNodeDataToVTK( const Model<dim>&  sg, 
                                              const std::string& file_name,
                                              T timestep )
  {
      OutputNodeDataToVTK( sg, "Model", file_name, timestep );
  }

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,int);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,int);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,int);

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,long);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,long);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,long);

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,size_t);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,size_t);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,size_t);

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,double);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,double);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,double);





/** output of variables placed on the nodes to VTK
*/
template<size_t dim>
template<class T>
void VTK_Interface<dim>::OutputNodeDataToVTK( const Model<dim>&  sg,
                                              const std::string& region,
                                              const std::string& initial_file_name,
                                              T timestep )
  {
     std::string file_name( this->OutputFileAndSubFolderName( initial_file_name ) );

     const Region<dim>&  gref(sg.Region(region));
     gref.UpdateMemberIndexes();

     // 0. creating list of consecutive element ID 0..n-1
     // -------------------------------------------------
     vector<size_t>  elmt_ids; elmt_ids.reserve( gref.Elements() );
     for ( size_t n=0U; n<elmt_ids.capacity(); n++ ) elmt_ids.push_back(n); 
       
     //   reading node properties in alphabetical order
     // -----------------------------------------------
     string       variable;
     set<string>  node_props;

     sg.Database().ListProperties( NODE, node_props );
     csmp::Index  prop_key(sg.Database().StorageKey((*node_props.begin()).c_str())); 
     

     // 1. getting new node mapping and updating storage if geometry has changed
     // ------------------------------------------------------------------------
     if ( last_visualized_ != NODE ) {
          NodeBasedTopology( gref, elmt_ids, plist, node_mapping );
          TransformPlist( gref, plist, transformed_plist );
          last_visualized_ = NODE;
       }


     // 2. opening data output file in ascii format
     // -------------------------------------------
     string region_name(region);
     replaceWhiteSpaceBy( region_name, '_' );
     char  outfile[NAME_STRING];
     if ( region_name != "Model" ) {
          strcpy( outfile, region_name.c_str() );
          strcat( outfile, "_" );
          strcat( outfile, region_name.c_str() );
       }
     else strcpy( outfile, file_name.c_str() );
     strcat( outfile, to_string(timestep).c_str() );
     strcat( outfile, ".vtk" );
       
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           throw csmp::Exception( CSMP_ERROR, "VTK_Interface<dim>::OutputDataToVTK", 
                                  "Output file could not be opened");
           return;
       }  
       
     // 2. writing the file header
     // --------------------------
     ofs.precision(15); // double precision has at list 15 significant digits
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): Complete set of NODE properties at timestep: "<< timestep << endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3. writing node coordinates & getting the first dataset
     // -------------------------------------------------------
     variable = (*node_props.begin());
     replaceWhiteSpaceBy( variable, '_' );
     NodeData( gref, prop_key, node_mapping, pxyz_data );
     
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << pxyz_data.size() <<" double"<< endl;
     for ( map<size_t,vector<double> >::const_iterator
           nit=pxyz_data.begin(); nit!=pxyz_data.end(); nit++ )
       {
          for ( size_t i=0U; i<dim; i++ ) ofs << (*nit).second[i] <<" ";
          if ( dim == 2U ) ofs << 0. <<"  ";
          ofs << endl;
       }
     ofs << endl;  
       
     // 4. getting total number of connections + numbers giving connections per element
     // -------------------------------------------------------------------------------
     size_t cell_list_size(0);                            
     for ( deque<vector<size_t> >::const_iterator 
           eit=transformed_plist.begin(); eit!=transformed_plist.end(); eit++ )
       cell_list_size += (*eit).size() + 1U;
     
     // 5. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< transformed_plist.size() <<" "<< cell_list_size << endl;
     for ( deque<vector<size_t> >::const_iterator
           it=transformed_plist.begin(); it!=transformed_plist.end(); it++ )
       {
          ofs << (*it).size() <<" ";
          for ( vector<size_t>::const_iterator 
                i=(*it).begin(); i!=(*it).end(); i++ ) ofs << *i <<" ";
          ofs << endl;
       }   
     ofs << endl;
     
     // 6. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< geometric_primitives_VTK.size() << endl;
     for ( typename deque<VTK_TYPE>::const_iterator
           vit=geometric_primitives_VTK.begin(); vit!=geometric_primitives_VTK.end(); vit++ ) ofs << (*vit) << endl;    
     ofs << endl;

     // 7. writing POINT_DATA point-type data values
     // --------------------------------------------
     size_t  line_break, offset(3); // offset for case where x,y,z are stored
     bool    first_iteration(true);

     // property after property
     for ( set<string>::const_iterator
           npit=node_props.begin(); npit!=node_props.end(); npit++ )
       {
          cout <<"\n\tOutputting property: '"<< (*npit) <<"' to VTK file..."<< endl;
          line_break = 1;
          // getting the property data, but only after first set was written
          if ( first_iteration ) {
               ofs <<"POINT_DATA "<< pxyz_data.size() << endl;
               ofs.setf( ios::scientific );
               first_iteration = false;
            }
          else {
               csmp::Index prop_key1(sg.Database().StorageKey((*npit).c_str()));
               variable = (*npit).c_str();
               replaceWhiteSpaceBy( variable, '_');
               // now only get data without coordinates
               RetrieveData( gref, prop_key1, node_mapping, pxyz_data );
               offset = 0; 
            }
       
          // writing the property data
          switch( prop_key.type )
            {
               case SCALAR:
                    ofs <<"SCALARS "<< variable <<" double"<< endl;
                    ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
                    for ( map<size_t,vector<double> >::const_iterator
                          nit=pxyz_data.begin(); nit!=pxyz_data.end(); nit++, line_break++ )
                      {
                         ofs << (*nit).second[offset] <<" ";
                         if ( line_break == 4 )
                           {
                              ofs << endl;
                              line_break = 0;
                           }
                      }
                    ofs << endl;
                 break;

               case VECTOR:
                    ofs <<"VECTORS "<< variable <<" double"<< endl;
                    for ( map<size_t,vector<double> >::const_iterator
                          nit=pxyz_data.begin(); nit!=pxyz_data.end(); nit++ )
                      {
                         // variables have always 3 components since view screen is 3D
                         for ( size_t j=offset; j<(*nit).second.size(); j++ ) 
                           ofs << (*nit).second[j] <<" ";
                         ofs << endl;
                      }
                break;
                
               case TENSOR:
                    ofs <<"TENSORS "<< variable <<" double"<< endl;
                    for ( map<size_t,vector<double> >::const_iterator
                          nit=pxyz_data.begin(); nit!=pxyz_data.end(); nit++ )
                      {
                         // always 3x3 components in the tensors
                         for ( size_t j=offset; j<(*nit).second.size(); j++ ) 
                           ofs << (*nit).second[j] <<" ";
                         ofs << endl;
                      }
                 break;
               default:
                 throw csmp::Exception( CSMP_ERROR, "VTK_Interface::OutputNodeDataToVTK:",
                                       "treatment of Array and FlaggedArray variables not handled yet.");
            }
            
       } // end for properties
       
     ofs.close();
     cout <<"\nVTK_Interface<"<<  dim;
     cout <<">::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,const std::string&,int);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,const std::string&,int);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,const std::string&,int);

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,const std::string&,long);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,const std::string&,long);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,const std::string&,long);

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,const std::string&,size_t);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,const std::string&,size_t);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,const std::string&,size_t);

template void VTK_Interface<1U>::OutputNodeDataToVTK(const Model<1U>&,const std::string&,const std::string&,double);
template void VTK_Interface<2U>::OutputNodeDataToVTK(const Model<2U>&,const std::string&,const std::string&,double);
template void VTK_Interface<3U>::OutputNodeDataToVTK(const Model<3U>&,const std::string&,const std::string&,double);










/** 

OutputDataToVTK writes the values of the specified distributed Model
variable to a VTK textfile of unstructured data.  

Note that in order to display element variables in VTK, flat shading must 
be set. Otherwise the variable values are again smeared out across 
neighboring nodes.   

@section arguments Input Arguments 

The first method argument is a reference to the Model object from 
which the data shall be retrieved. The second and third string arguments 
specify the name of the output file and the physical variable which shall
be output. The fourth argument specifies the simulation timestep which 
the dataset represents. This timestep is appended to the output file name.
The fifth method argument specifies whether OutputDataToVTK() should 
update the internally stored mesh topology or whether it can use again one 
which it created at an earlier output step. The default value here tells
the method to update the topology.   

@section implementation Implementation

From a 'plist', a unique list of nodes is created which are renumbered
from 0...n-1.  

The original 'plist' is reassembled by braking more complex elements down
into topological primitives which can be displayed in VTK. The new plist
then also contains the new node IDs which the VTK_Interface uses 
internally.   */
template<size_t dim>
template<class T>
void VTK_Interface<dim>::OutputDataToVTK( const Model<dim>&  sg,
                                          const std::string& file_name,
                                          const std::string& var_name,
                                          T    timestep,
                                          bool update_topology )
  {
     OutputDataToVTK( sg, "Model", file_name, var_name, timestep, update_topology );

  } // end OutputDataToVTK

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,unsigned int,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,unsigned int,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,unsigned int,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,int,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,int,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,int,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,long,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,long,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,long,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,size_t,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,size_t,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,size_t,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,double,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,double,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,double,bool);







/** To output subregions of a model represented by Region objects to a VTK
file in text format.  

@section arguments Input Arguments 

The first method argument is a reference to the Model object 
representing the model from which data shall be output to VTK. The following
string arguments give the names of (1) the region object inside the Model
to which the interface is applied, (2) the name of outputfile, 
and (3) the name of the output variable. 
The fifth method argument further supplies the current timestep of the
computation and the sixth argument tells the VTK_Interface whether
the topology has to be updated or whether it can be re-used. The default
value leads to updates.  

@section implementation Implementation

1. Makes a first map that contains a uniquely numbered list of the node coordinates
abd the old node IDs

2. Create a second map that contains new node IDs that correspond to the alreadzy
existing node coordinates

3. Build a third map for the Plist that contains the new node IDs of the new subtriangles

(4 subtriangles for the quadratic triangle , 6 subtriangles for the barycentric
quadratic triangle) in a counter clockwise order 

@section application Application

To visualize CSMP models using the Visualization Toolkit.
*/
template<size_t dim>
template<class T>
void VTK_Interface<dim>::OutputDataToVTK( const Model<dim>&  sg, 
                                          const std::string& group_name,
                                          const std::string& initial_file_name,
                                          const std::string& var_name,
                                          T     timestep,
                                          bool  update_topology )
  {
     std::string file_name( this->OutputFileAndSubFolderName( initial_file_name ) );

     csmp::Index  prop_key = sg.Database().StorageKey(var_name.c_str());

     if ( prop_key.place == REGION or prop_key.place == BOUNDARY )
       throw Exception( CSMP_ERROR, "VTK_Interface<dim>::OutputDataToVTK(region):",
                      "Thus far variables placed on REGION or BOUNDARY cannot be visualised with VTK (this could however be done with the VTK primitive POLYGONAL)");

     const Region<dim>&  gref(sg.Region(group_name));
     gref.UpdateMemberIndexes();  // renumber nodes and elements
     
     string variable(var_name);
     replaceWhiteSpaceBy( variable, '_' );

     // 0. Retrieving list of ID's of the Elements which make up the 
     //    target Region object.
     // ------------------------------------------------------------
     // finding the group in the group list
     vector<size_t>  elmt_ids;
     gref.MemberElementIndexes( elmt_ids );

     if ( elmt_ids.empty() )     
       throw csmp::Exception( CSMP_ERROR, "VTK_Interface<dim>::OutputDataToVTK(region)",
                              group_name, "region is empty");

     // 1. getting new node mapping and updating storage if geometry has changed
     // ------------------------------------------------------------------------
     if ( update_topology || last_visualized_ != prop_key.place ) 
       {
          // SKM fix for element based vector and tensor data
          if ( prop_key.place == ELEMENT and prop_key.type != SCALAR ) {
                geometric_primitives_VTK.clear();
                geometric_primitives_VTK.resize( plist.size() );
                fill( geometric_primitives_VTK.begin(), geometric_primitives_VTK.end(), VTK_VERTEX );
                transformed_plist.clear();
                // inserting just the element id into the plist
                size_t eidx(0U);
                for ( std::map<size_t,std::vector<size_t> >::const_iterator
                      eit=plist.begin(); eit!=plist.end(); eit++ )
                  // SKM FIX (only single point needs to be stored) transformed_plist.push_back((*eit).second);
                  transformed_plist.push_back(vector<size_t>(1,eidx++));
            }
          else if ( prop_key.place != ELEMENT_INTEGRATION_POINT ) {
               NodeBasedTopology( gref, elmt_ids, plist, node_mapping );
               TransformPlist( gref, plist, transformed_plist );
            }
          else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
               IntegrationPointBasedTopology( gref, elmt_ids, plist, node_mapping );
               // fill transform plist - AP Aug2006
               transformed_plist.clear();
               for ( std::map<size_t,std::vector<size_t> >::const_iterator
                     eit=plist.begin(); eit!=plist.end(); eit++ ) 
                  transformed_plist.push_back((*eit).second);
           }
          else
            throw csmp::Exception( CSMP_ERROR, "VTK_Interface<dim>::OutputDataToVTK(region)",
                                  "The placement of the variable was not recognized");
       }
     last_visualized_ = prop_key.place;

     // not the topology but only the data have to be updated
     if ( prop_key.place == NODE ) 
       NodeData( gref, prop_key, node_mapping, pxyz_data );
     else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) 
       IntegrationPointData( gref, prop_key, pxyz_data );
     else if ( prop_key.place == ELEMENT  or  prop_key.place == FACE  or prop_key.place == INTER_FACE  ) {
          if ( prop_key.type != SCALAR ) {
                PointBasedTopology( gref, elmt_ids, plist, node_mapping );
                ElementPointData( gref, prop_key, node_mapping, pxyz_data );
            }
          else NodeCoordinates( gref, node_mapping, pxyz_data );
       }

     // 2. write data into a stringstream
     // -------------------------------------------
     stringstream  ofs;
       
     // 3.1 writing the file header
     // ---------------------------
     ofs.precision(15); // double precision
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name;
     ofs <<", model domain: "<< group_name <<", timestep: "<< timestep << endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3.2 writing coordinates
     // ---------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << pxyz_data.size() <<" double"<< endl;
     
     typename map<size_t,vector<double> >::const_iterator  nit;
     const typename map<size_t,vector<double> >::const_iterator  nit_end(pxyz_data.end());
     for ( nit=pxyz_data.begin(); nit!=nit_end; nit++ )
       {
          for ( size_t i=0; i<dim; i++ ) ofs << (*nit).second[i] <<" ";
          if ( dim == 2 ) ofs << 0.0 <<"  ";
          ofs << endl;
       }
     ofs << endl;  
       
     // 3.3 getting total number of connections + numbers giving connections per element
     // -------------------------------------------------------------------------------
     size_t cell_list_size(0);                            
     deque<vector<size_t> >::const_iterator  eit;
     const deque<vector<size_t> >::const_iterator  eit_end(transformed_plist.end());
     for ( eit=transformed_plist.begin(); eit!=eit_end; eit++ )
       cell_list_size += (*eit).size() + 1;
     
     // 3.4 writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< transformed_plist.size() <<" "<< cell_list_size << endl;
     for ( eit=transformed_plist.begin(); eit!=eit_end; eit++ )
       {
          ofs << (*eit).size() <<" ";
          for ( vector<size_t>::const_iterator
                it=(*eit).begin(); it!=(*eit).end(); it++ ) ofs << *it <<" ";
          ofs << endl;
       }   
     ofs << endl;
     
     // 3.5 writing CELL_TYPES
     // ----------------------
     ofs <<"CELL_TYPES "<< geometric_primitives_VTK.size() << endl;
     deque<VTK_TYPE>::const_iterator  vit;
     const deque<VTK_TYPE>::const_iterator  vit_end(geometric_primitives_VTK.end());
     for ( vit=geometric_primitives_VTK.begin(); vit!=vit_end; vit++ ) 
       ofs << (*vit) << endl;    
     ofs << endl;

     // 3.6 writing POINT_DATA or CELL_DATA point-type data values
     // ----------------------------------------------------------
     int   line_break(1);
     size_t  offset(3U);
     if ( prop_key.type == SCALAR  and  (prop_key.place == ELEMENT  or prop_key.place == FACE  or prop_key.place == INTER_FACE ) ) {
          // resizing cell data record
          CellData( gref, prop_key, plist, pxyz_data );
//          assert( pxyz_data.size() == transformed_plist.size() );
          ofs <<"CELL_DATA "<< transformed_plist.size() << endl;
          offset = 0U;
       }
     else ofs <<"POINT_DATA "<< pxyz_data.size() << endl;
     ofs.setf( ios::scientific );
     switch( prop_key.type )
       {
          case SCALAR:
               ofs <<"SCALARS "<< variable <<" double"<< endl;
               ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
               for ( nit=pxyz_data.begin(); nit!=nit_end; nit++, line_break++ )
                 {
                    ofs << (*nit).second[offset] <<" ";
                    if ( line_break == 4 ) {
                         ofs << endl;
                         line_break = 0;
                      }
                 }
            break;

          case VECTOR:
               ofs <<"VECTORS "<< variable <<" double"<< endl;
               for ( nit=pxyz_data.begin(); nit!=nit_end; nit++ )
                 {
                    // variables have always 3 components since view screen is 3D
                    for ( size_t j=offset; j<(*nit).second.size(); j++ ) 
                      ofs << (*nit).second[j] <<" ";
                    ofs << endl;
                 }
           break;
          case TENSOR:
               ofs <<"TENSORS "<< variable <<" double"<< endl;
               for ( nit=pxyz_data.begin(); nit!=nit_end; nit++ )
                 {
                    // always 3x3 components in the tensors
                    for ( size_t j=offset; j<(*nit).second.size(); j++ ) 
                      ofs << (*nit).second[j] <<" ";
                    ofs << endl;
                 }
             break;
           default:
             throw csmp::Exception( CSMP_ERROR, "VTK_Interface<dim>::OutputDataToVTK(region):",
                                   "treatment of Array and FlaggedArray variables not handled yet.");
       }
     ofs << endl;
     
     // 2. opening data output file in ascii format
     // -------------------------------------------
     ostringstream  outfile;
     // if we are dealing with the Model as a whole, no name prefix is used
     if ( strncmp( group_name.c_str(), "Model", NAME_STRING ) == 0 ) outfile << file_name << timestep << ".vtk";
     else outfile << group_name << "_" << file_name << timestep << ".vtk";
     string file(outfile.str());
     for(string::iterator i = file.begin(); i < file.end(); i++ )
      if( *i == ' ' ) *i = '_'; 
     
     //write it into a file
     ofstream ofs_file;
     ofs_file.open( file.c_str(), ios::out|ios::trunc );
     if ( !ofs ) {
           throw csmp::Exception( CSMP_ERROR, 
                          "VTK_Interface<dim>::OutputDataToVTK(region)",
                          "Output file could not be opened");
           return;
       }  
     
     ofs_file << ofs.str();
     
     ofs_file.close();
     
     cout <<"\nVTK_Interface<"<<  dim;
     cout <<">::OutputDataToVTK(region): file '"<< outfile.str() <<"' written successfully."<< endl;

  } // end OutputDataToVTK (Region version)


template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,const std::string&,int,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,const std::string&,int,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,const std::string&,int,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,const std::string&,long,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,const std::string&,long,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,const std::string&,long,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,const std::string&,size_t,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,const std::string&,size_t,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,const std::string&,size_t,bool);

template void VTK_Interface<1U>::OutputDataToVTK(const Model<1U>&,const std::string&,const std::string&,const std::string&,double,bool);
template void VTK_Interface<2U>::OutputDataToVTK(const Model<2U>&,const std::string&,const std::string&,const std::string&,double,bool);
template void VTK_Interface<3U>::OutputDataToVTK(const Model<3U>&,const std::string&,const std::string&,const std::string&,double,bool);



/**
    deduces which vtkCellType (vtkCellType.h, vtk.org) is best suited
    to visualise the supplied CSMP finite element type.
    
    @attention method returns VTK_POLYGON if the correct element type cannot be found.
    
    @note the barycentric element types have a free floating node in the center. 
    in VTK these elements are called biquadratic but not all of them are there. Missing ones
    can be dealt with by adding an extra point to the normal primitive. Alternatively,
    one could display such elements using VTK_TRIANGLE_STRIP, i.e., a node by node
    linear tesselation.
    
    @note there are more VTK types that parse here, for instance the 20-noded vtkBiQuadraticQuadraticHexahedron; 
    @todo check list after major VTK updates to track changes
    
    @note to date, there are no cubic elements in VTK
    
    @author SKM 2/2/2016
*/
VTK_TYPE parseElementType( CSMP_FEM_TYPE etype )
 {
    if ( etype == LINEAR_BAR ) return VTK_LINE;
    if ( etype == QUADRATIC_BAR ) return VTK_QUADRATIC_EDGE;
//    if ( etype == CUBIC_BAR ) return VTK_POLYGON;
    if ( etype == LINEAR_TRIANGLE ) return VTK_TRIANGLE;
    if ( etype == LINEAR_TRIANGLE3D ) return VTK_TRIANGLE;
//    if ( etype == BARYCENTRIC_LINEAR_TRIANGLE ) return VTK_POLYGON;
    if ( etype == QUADRATIC_TRIANGLE ) return VTK_QUADRATIC_TRIANGLE;
//    if ( etype == BARYCENTRIC_QUADRATIC_TRIANGLE ) return VTK_POLYGON;
//    if ( etype == CUBIC_TRIANGLE ) return VTK_POLYGON;
    if ( etype == LINEAR_TETRAHEDRON ) return VTK_TETRA;
    if ( etype == QUADRATIC_TETRAHEDRON ) return VTK_QUADRATIC_TETRA;
//    if ( etype == BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_BAR ) return VTK_LINE;
    if ( etype == ISOPARAMETRIC_QUADRATIC_BAR ) return VTK_QUADRATIC_EDGE;
//    if ( etype == ISOPARAMETRIC_CUBIC_BAR ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_TRIANGLE ) return VTK_TRIANGLE;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_TRIANGLE ) return VTK_BIQUADRATIC_TRIANGLE;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TRIANGLE ) return VTK_QUADRATIC_TRIANGLE;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE ) return VTK_BIQUADRATIC_TRIANGLE;
//    if ( etype == ISOPARAMETRIC_CUBIC_TRIANGLE ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_TETRAHEDRON ) return VTK_TETRA;
    if ( etype == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON ) return VTK_QUADRATIC_TETRA;
//    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON ) return VTK_POLYGON;
//    if ( etype == ISOPARAMETRIC_CUBIC_TETRAHEDRON ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_PYRAMID ) return VTK_PYRAMID;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID13 ) return VTK_QUADRATIC_PYRAMID;
//    if ( etype == ISOPARAMETRIC_QUADRATIC_PYRAMID14 ) return VTK_POLYGON;
//    if ( etype == ISOPARAMETRIC_CUBIC_PYRAMID ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_PRISM ) return VTK_WEDGE;
    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM15 ) return VTK_QUADRATIC_WEDGE;
//    if ( etype == ISOPARAMETRIC_QUADRATIC_PRISM18 ) return VTK_POLYGON;
//    if ( etype == ISOPARAMETRIC_CUBIC_PRISM ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) return VTK_QUAD;
//    if ( etype == ISOPARAMETRIC_BARYCENTRIC_LINEAR_QUADRILATERAL ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_QUADRILATERAL ) return VTK_BIQUADRATIC_QUAD;
    if ( etype == ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9 ) return VTK_BIQUADRATIC_QUAD;
//    if ( etype == ISOPARAMETRIC_CUBIC_QUADRILATERAL ) return VTK_POLYGON;
    if ( etype == ISOPARAMETRIC_LINEAR_HEXAHEDRON ) return VTK_HEXAHEDRON;
    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20 ) return VTK_QUADRATIC_HEXAHEDRON;
//    if ( etype == ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27 ) return VTK_POLYGON;
//    if ( etype == ISOPARAMETRIC_CUBIC_HEXAHEDRON ) return VTK_POLYGON;
    if ( etype == ZERO_DIMENSIONAL_FACE ) return VTK_VERTEX;
    if ( etype == POINT_ELEMENT ) return VTK_VERTEX;
    if ( etype == POLYGONAL_ELEMENT ) return VTK_POLYGON;
    if ( etype == POLYHEDRAL_ELEMENT ) return VTK_POLYGON;
    if ( etype == EXPERIMENTAL_ELEMENT ) return VTK_POLYGON;

    return VTK_POLYGON;

 } // end parseElementType





/** The global element ID's are retained as keys of the 'plist' such that they
can be used at a later stage to query the Model.  
*/
template<size_t dim>
void VTK_Interface<dim>::NodeBasedTopology( const Region<dim>& sgref,
                                            const vector<size_t>& elmt_ids,
                                            map<size_t,vector<size_t> >& plist,
                                            map<size_t,size_t>& node_nums )
 {
    size_t             nodes(0U);
    vector<size_t>     pentry(3);
    
    plist.clear();
    node_nums.clear();

    // 1. creating unique node number list, and plist by looping over the selected elements
    // ------------------------------------------------------------------------------------
    vector<size_t>  empty_vec;       
    for ( typename vector<size_t>::const_iterator
          lit=elmt_ids.begin(); lit!=elmt_ids.end(); lit++ )
      {
         // getting node ids and renumbering them 0...n-1
         for ( size_t i=0U; i<sgref.E( (*lit) )->Nodes(); i++ ) {
              pair<typename map<size_t,size_t>::iterator,bool>
                node_it = node_nums.insert( make_pair( sgref.E( (*lit) )->N(i)->Idx(), nodes ) );
              if ( node_it.second == true ) nodes++;
           }
         
         // building the plist, minimizing the search by always using the smallest 
         // size of the map possible
         pentry.resize( sgref.E( (*lit) )->Nodes() );
         for ( size_t i=0U; i<sgref.E( (*lit) )->Nodes(); i++ ) {
              typename map<size_t,size_t>::const_iterator
                nit = node_nums.find( sgref.E( (*lit) )->N(i)->Idx() );
              pentry[i] = (*nit).second;
           }
         
         // inserting new element id and empty vector<double> into the plist
         pair<typename map<size_t,vector<size_t> >::iterator,bool>
           pit = plist.insert( make_pair( (*lit), empty_vec ) );
         if ( pit.second == true ) {
              // initializing plist vector
              (*pit.first).second.reserve( sgref.E( (*lit) )->Nodes() );
              for ( size_t i=0U; i<sgref.E( (*lit) )->Nodes(); i++ ) 
                (*pit.first).second.push_back( pentry[i] );
           }
      }   
     
 } // end NodeBasedTopology











/**

Since VTK cannot associate properties with cells, it must be tricked into
treating each cell separately. This is done (unfortunately) by introducing
redundant nodes such that each triangle / tetrahedron has its own,
unique set of nodes. The element properties are then assigned to all
the nodes of each element. To maximize efficiency in VTK flat shading
should be used to display this kind of element data.  

*/
template<size_t dim>
void VTK_Interface<dim>::ElementBasedTopology( const Region<dim>& sgref,
                                               const vector<size_t>& elmt_ids,
                                               map<size_t,vector<size_t> >& plist,
                                               map<size_t,size_t>& node_nums )
 {
    size_t                        nodes(0U);
    vector<size_t>                pentry(3);
    pair<size_t,vector<size_t> >  alpha_help;
    vector<size_t>                alpha_vec;
    
    plist.clear();
    node_nums.clear();

    // 1. creating unique node number list, and plist by looping over the selected elements
    // ------------------------------------------------------------------------------------
    for ( typename vector<size_t>::const_iterator 
          lit=elmt_ids.begin(); lit!=elmt_ids.end(); lit++ )
      {
         // getting node ids and renumbering them including duplicates 0...elmts * npe's
         pentry.resize( sgref.E( (*lit) )->Nodes() );
         for ( size_t i=0U; i<sgref.E( (*lit) )->Nodes(); i++ ) {
              node_nums[ nodes ] = sgref.E( (*lit) )->N(i)->Idx();
              pentry[i]          = nodes++;
           }
         // inserting new element id and empty vector<double> into the plist
         alpha_help.first  = (*lit);
         alpha_help.second = alpha_vec;
         pair<typename map<size_t,vector<size_t> >::iterator,bool> pit = plist.insert( alpha_help );
         assert( pit.second == true );
         
         // initializing plist vector
         (*pit.first).second.reserve( sgref.E( (*lit) )->Nodes() );
         for ( size_t i=0U; i<sgref.E( (*lit) )->Nodes(); i++ ) 
           (*pit.first).second.push_back( pentry[i] );
      }   

  } // end ElementBasedTopology
                                






/**
 
IntegrationPoints are visualized as point data. Such points are introduced
as nodes with locations that correspond to those of the integration points
of the corresponding finite elements. Thus, each element has its own 
unique set of constraint points. 

The constraint point properties are assigned later to each of these.  

@section implementation Implementation

The constraint points are treated like nodes which are not shared 
among the the elements.  

@section application Application

Output of variables that are stored at IntegrationPoints. These might
be strain, stress or similar computational variables. 

@todo SKM: sector and facet integration point placements have to be added
 
*/
template<size_t dim>
void VTK_Interface<dim>::IntegrationPointBasedTopology( const Region<dim>& sgref,
                                                       const vector<size_t>& elmt_ids,
                                                       map<size_t,vector<size_t> >&  clist,
                                                       map<size_t,size_t>&  cpoint_nums )
 {
    vector<size_t>                pentry(3);
    pair<size_t,vector<size_t> >  alpha_help;
    vector<size_t>                alpha_vec;
    
    // the clist will be used to hold the constraint point ids per element
    // -------------------------------------------------------------------
    clist.clear();
    cpoint_nums.clear();

    // 1. creating unique constraint point number list, and entering the cpoints into the 
    //    plist by looping over the selected elements
    // ------------------------------------------------------------------------------------
    size_t  cpoints(0U);
    for ( typename vector<size_t>::const_iterator
          lit=elmt_ids.begin(); lit!=elmt_ids.end(); lit++ )
      {
         // getting constraint point ids and data 
         pentry.resize( sgref.E( (*lit) )->IntegrationPoints() );
         for ( size_t i=0U; i<sgref.E( (*lit) )->IntegrationPoints(); i++ ) 
           {
              // associate the unique number of the constraint-point with the parent element 
              // for later retrieval in IntegrationPointData method
              cpoint_nums.insert( make_pair( cpoints, (*lit) ) ); 
              pentry[i] = cpoints++;
           }
         // inserting new element idx and empty vector<double> into the clist
         alpha_help.first  = (*lit);
         alpha_help.second = alpha_vec;
         pair<typename map<size_t,vector<size_t> >::iterator,bool> 
         cit = clist.insert( alpha_help );
         assert( cit.second == true );
         
         // initializing clist vector
         (*cit.first).second.reserve( sgref.E( (*lit) )->IntegrationPoints() );
         for ( size_t i=0U; i<sgref.E( (*lit) )->IntegrationPoints(); i++ ) 
           (*cit.first).second.push_back( pentry[i] );
      }  
      
    // 2. assigning the VTK type for the clist entries
    // -----------------------------------------------
    geometric_primitives_VTK.erase( geometric_primitives_VTK.begin(), geometric_primitives_VTK.end() );
    geometric_primitives_VTK.resize( clist.size() );
    fill( geometric_primitives_VTK.begin(), geometric_primitives_VTK.end(), VTK_POLY_VERTEX );

  } // end IntegrationPointBasedTopology









/**
    Rebuilds 'plist' and 'node_nums' creating a one-to-one mapping between elements and data points.
*/
template<size_t dim>
void VTK_Interface<dim>::PointBasedTopology( const Region<dim>& sgref,
                                             const vector<size_t>&  elmt_ids,
                                             map<size_t,vector<size_t> >&  plist,
                                             map<size_t,size_t>&  node_nums )
 {
    size_t          nodes(0U);
    vector<size_t>  pentry(1);
    
    plist.clear();
    node_nums.clear();
    
    // 1. creating unique point number list, and plist by looping over the selected elements
    // ------------------------------------------------------------------------------------
    for ( typename vector<size_t>::const_iterator
          lit=elmt_ids.begin(); lit!=elmt_ids.end(); lit++ )
      {
         // getting node ids and renumbering them including duplicates 0...elmts * npe's
         node_nums[ nodes ] = sgref.E( (*lit) )->Idx();
         pentry[0] = nodes++;
         // inserting new element id and single-element vector<double> into plist
         pair<typename map<size_t,vector<size_t> >::iterator,bool>
           pit = plist.insert( make_pair( (*lit), pentry ) );
         assert( pit.second == true );
      }   

  } // end PointBasedTopology






// only the coordinates of the nodes - no data
template<size_t dim>
void VTK_Interface<dim>::NodeCoordinates( const Region<dim>& sgref,
                                          const map<size_t,size_t>& node_nums,
                                          map<size_t,vector<double> >& pxyz_data )
 {
    assert( !node_nums.empty() );
 
    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );

    // extracting the node coordinates and property values from the Model
    // results are put into the first three elements of 'pxz_data'
    // --------------------------------------------------------------------------
    vector<double>  empty_vec;
    for ( typename map<size_t,size_t>::const_iterator
          nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
      {
         // inserting new element into the cordinate map using new node coordinate number
         pair<typename map<size_t,vector<double> >::iterator,bool>
         dit = pxyz_data.insert( make_pair( (*nit).second, empty_vec ) );
         assert( dit.second == true );
         
         // node coordinates (reserves storage for three coordinates and a scalar data value)
         (*dit.first).second.reserve(3U);
         (*dit.first).second.push_back( sgref.N( (*nit).first )->x() );
         if ( dim >= 2 ) (*dit.first).second.push_back( sgref.N( (*nit).first )->y() );
         else            (*dit.first).second.push_back( 0. );
         if ( dim == 3 ) (*dit.first).second.push_back( sgref.N( (*nit).first )->z() );
         else            (*dit.first).second.push_back( 0. );
      }
   
   assert( pxyz_data.size() == node_nums.size() );
     
 } // end NodeCoordinates







/** Method retrieves node coordinates and nodal property values from the
Model, for the property specified by the csmp::Index.  

@section arguments Input Arguments 

The first argument is a constant reference to the Model object,
the second specifies the property which will be output into the fourth
method argument. The third argument gives the mapping between Model
node numbers and the node numbers used by the VTK_Interface. Thus
the supplied map contains pairs of the old and the corresponding
new node numbers.  

A map which contains vectors holding the node coordinates (first three
entries and the property values. Depending whether the property is of
scalar, vector, or tensor type, the data segment of the vector<double> contains
1, 3, ord 9 entries respectively.  

*/
template<size_t dim>
void VTK_Interface<dim>::NodeData( const Region<dim>& sgref,
                                   const csmp::Index&  prop_key,
                                   const map<size_t,size_t>& node_nums,
                                   map<size_t,vector<double> >& pxyz_data )
 {
    if ( prop_key.place != NODE )
      throw csmp::Exception( CSMP_FATAL_ERROR, "VTK_Interface<dim>::NodeData", 
                                   "Method only applies to node properties" );
    ScalarVariable      sc;
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;
    
    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );

    // extracting the node coordinates and property values from the Model
    // results are put into the first three elements of 'pxz_data'
    // --------------------------------------------------------------------------
    vector<double>  empty_vec;
    for ( typename map<size_t,size_t>::const_iterator
          nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
      {
         // inserting new element into the cordinate map using new node coordinate number
         pair<typename map<size_t,vector<double> >::iterator,bool>
           dit = pxyz_data.insert( make_pair( (*nit).second, empty_vec ) );
         assert( dit.second == true );
         
         // node coordinates (reserves storage for three coordinates and a scalar data value)
         (*dit.first).second.reserve(4);
         (*dit.first).second.push_back( sgref.N( (*nit).first )->x() );
         (*dit.first).second.push_back( sgref.N( (*nit).first )->y() );
         if ( dim == 3 ) (*dit.first).second.push_back( sgref.N( (*nit).first )->z() );
         else            (*dit.first).second.push_back( 0. );
         
         // data values
         if ( prop_key.type == SCALAR )
           {
              sgref.N( (*nit).first )->Read( prop_key, sc );
              (*dit.first).second.push_back( sc() );
           }
         else if ( prop_key.type == VECTOR )
           {
              // storage for 3 coordinate values and 3 data values (in 2D 3rd place = 0.0) 
              (*dit.first).second.reserve(6);
              sgref.N( (*nit).first )->Read( prop_key, vc );
              // variables have always 3 components since view screen is 3D
              for ( size_t j=0U; j<2U; j++ ) (*dit.first).second.push_back( vc[j] );
              if ( dim == 3U ) (*dit.first).second.push_back( vc[2] );
              else             (*dit.first).second.push_back( 0. );
           }
         else if ( prop_key.type == TENSOR )
           {
             // storage for 3 coordinate values and 9 data values (in 2D 3rd row, column = 0.0) 
             (*dit.first).second.reserve(12);
              sgref.N( (*nit).first )->Read( prop_key, ts );
              // variables have always 3 components since view screen is 3D
              if ( dim == 3U )
              for ( size_t k=0U; k<3U; k++ )
                for ( size_t l=0U; l<3U; l++ ) (*dit.first).second.push_back( ts(k,l) );
              else {
                   (*dit.first).second.push_back( ts(0,0) );
                   (*dit.first).second.push_back( ts(0,1) );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( ts(1,0) );
                   (*dit.first).second.push_back( ts(1,1) );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( 0. );
                }
           }
      }
     
 } // end NodeData







/*  adds the variable values at the constraint points to the 'pxyz_data' array

clist contains the coordinates of the constraint points

*/
template<size_t dim>
void VTK_Interface<dim>::IntegrationPointData( const Region<dim>& gref,
                                              const csmp::Index&  prop_key,
                                              map<size_t,vector<double> >& pxyz_data )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( prop_key.place != ELEMENT_INTEGRATION_POINT )
      csmp_error.notice( CSMP_FATAL_ERROR, "VTK_Interface<dim>::IntegrationPointData", 
                                      "Method only applies to constraint point properties" );
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;
    
    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );

    // extracting constraint point coordinates and property values from Model
    // results are put into the first three elements of 'pxz_data'
    // --------------------------------------------------------------------------
    vector<double>  dentry;
    size_t            cpoints(0U);
    for ( typename vector<Element<dim>*>::const_iterator
          it=gref.ElementsBegin(); it!=gref.ElementsEnd(); it++ )
      {
         // getting the coordinates of the constraint point
         for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++, cpoints++ ) 
           {
              // inserting new point into the cordinate map using new constraint-point coordinate number
              pair<typename map<size_t,vector<double> >::iterator,bool>
              dit = pxyz_data.insert( make_pair( cpoints, dentry ) );
              assert( dit.second == true );
              // making space for missing coordinates and at least one scalar variable value
              (*dit.first).second.reserve(4U);
              Point<dim>  cp((**it).IntegrationPoint(i)); 
              (*dit.first).second.push_back( cp[0] );
              (*dit.first).second.push_back( cp[1] );
              (*dit.first).second.push_back( ((dim==3U) ? cp[2] : 0.) );
              
              // adding the variable data
              if ( prop_key.type == SCALAR ) (*dit.first).second.push_back( (*it)->Read( i, prop_key ) );
              else if ( prop_key.type == VECTOR ) {
                   // storage for 3 coordinate values and 3 data values (in 2D 3rd place = 0.0) 
                   (*dit.first).second.reserve(6U);
                   (*it)->Read( i, prop_key, vc );
                   // variables have always 3 components since view screen is 3D
                   for ( size_t j=0U; j<2U; j++ ) (*dit.first).second.push_back( vc[j] );
                   if ( dim == 3U ) (*dit.first).second.push_back( vc[2U] );
                   else             (*dit.first).second.push_back( 0. );
                 }
              else if ( prop_key.type == TENSOR ) {
                   // storage for 3 coordinate values and 9 data values (in 2D 3rd row, column = 0.0) 
                   (*dit.first).second.reserve(12U);
                   (*it)->Read( i, prop_key, ts );
                   // variables have always 3 components since view screen is 3D
                   if ( dim == 3U )
                    for ( size_t k=0U; k<3U; k++ )
                      for ( size_t l=0U; l<3U; l++ ) (*dit.first).second.push_back( ts(k,l) );
                   else {
                         (*dit.first).second.push_back( ts(0,0) );
                         (*dit.first).second.push_back( ts(0,1) );
                         (*dit.first).second.push_back( 0.0 );
                         (*dit.first).second.push_back( ts(1,0) );
                         (*dit.first).second.push_back( ts(1,1) );
                         (*dit.first).second.push_back( 0.0 );
                         (*dit.first).second.push_back( 0.0 );
                         (*dit.first).second.push_back( 0.0 );
                         (*dit.first).second.push_back( 0.0 );
                      }
                 } // tensor
              dentry.clear();
           } // cpoints
         
      } // cpoint map
     
 } // end IntegrationPointData








// just data without coordinates
template<size_t dim>
void VTK_Interface<dim>::RetrieveData( const Region<dim>& sgref,
                                      const csmp::Index&     prop_key,
                                      const map<size_t,size_t>& obj_nums,
                                      map<size_t,vector<double> >& sgdata )
 {
    if ( obj_nums.empty() ) {
         throw csmp::Exception( CSMP_ERROR, "VTK_Interface<dim>::RetrieveData",
                         "No output objects specified. Nothing was done." );
         return;
      }
    
    ScalarVariable      sc;
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;
    
    sgdata.erase( sgdata.begin(), sgdata.end() );

    for ( typename map<size_t,size_t>::const_iterator
          nit=obj_nums.begin(); nit!=obj_nums.end(); nit++ )
      {
         // inserting new element into the data map using new node number
         pair<typename map<size_t,vector<double> >::iterator,bool>
           dit = sgdata.insert( make_pair((*nit).second, vector<double>() ) );
         assert( dit.second == true );
         
         // node coordinates (reserves storage for a scalar data value)
         (*dit.first).second.reserve(1);
         
         // data values
         if ( prop_key.type == SCALAR )
           {
              sgref.N( (*nit).first )->Read( prop_key, sc );
              (*dit.first).second.push_back( sc() );
           }
         else if ( prop_key.type == VECTOR )
           {
             // storage for 3 data values (in 2D 3rd place = 0.0) 
             (*dit.first).second.reserve(3);
             sgref.N( (*nit).first )->Read( prop_key, vc );
              // variables have always 3 components since view screen is 3D
              if ( dim == 3U )
                for ( size_t j=0; j<3U; j++ ) (*dit.first).second.push_back( vc[j] );
              else {
                   (*dit.first).second.push_back( vc[0] );
                   (*dit.first).second.push_back( vc[1] );
                   (*dit.first).second.push_back( 0.0 );
                }
           }
         else if ( prop_key.type == TENSOR )
           {
             // storage for 9 data values (in 2D 3rd row, column = 0.0) 
             (*dit.first).second.reserve(9);
              sgref.N( (*nit).first )->Read( prop_key, ts );
              // variables have always 3 components since view screen is 3D
              if ( dim == 3U )
              for ( size_t k=0; k<3U; k++ )
                for ( size_t l=0; l<3U; l++ ) (*dit.first).second.push_back( ts(k,l) );
              else {
                   (*dit.first).second.push_back( ts(0,0) );
                   (*dit.first).second.push_back( ts(0,1) );
                   (*dit.first).second.push_back( 0.0 );
                   (*dit.first).second.push_back( ts(1,0) );
                   (*dit.first).second.push_back( ts(1,1) );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( 0. );
                   (*dit.first).second.push_back( 0. );
                }
           }
      }

 } // end RetrieveData








template<size_t dim>
void VTK_Interface<dim>::ElementData( const Region<dim>& sgref,
                                      const csmp::Index&     prop_key,
                                      const map<size_t,size_t>& node_nums,
                                      map<size_t,vector<double> >& pxyz_data )
 {
    if ( prop_key.place != ELEMENT and prop_key.place != REGION )
      throw csmp::Exception( CSMP_FATAL_ERROR, "VTK_Interface<dim>::ElementData", 
                                   "Method only applies to element or region properties" );
    ScalarVariable      sc;
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;
    
    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );
    
    // 1. extracting the node coordinates and property values from the Model
    //    results are put into the first three elements of 'pxz_data'
    // --------------------------------------------------------------------------
    vector<double>  empty_vec;
    for ( typename map<size_t,size_t >::const_iterator
          nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
      {
         // inserting new element into the cordinate map using new node coordinate number
         pair<typename map<size_t,vector<double> >::iterator,bool>
           dit = pxyz_data.insert( make_pair( (*nit).first, empty_vec ) );
         assert( dit.second == true );
         
         // node coordinates (reserves storage for three coordinates and a scalar data value)
         (*dit.first).second.reserve(4);
         (*dit.first).second.push_back( sgref.N( (*nit).first )->x() );
         (*dit.first).second.push_back( sgref.N( (*nit).first )->y() );
         if ( dim == 3U ) (*dit.first).second.push_back( sgref.N( (*nit).first )->z() );
         else             (*dit.first).second.push_back( 0. );
      }
      
    // 2. Looping through the plist, associating the element properties with the nodal pxyz_data     
    // data values
    // -----------------------------------------------------------------------------------------
    for ( typename map<size_t,vector<size_t> >::const_iterator
          it=plist.begin(); it!=plist.end(); it++ )
      {
         if ( prop_key.type == SCALAR )
           {
              sgref.E( (*it).first )->Read( prop_key, sc );
              // looping over the element's node 
              for ( typename vector<size_t>::const_iterator
                    pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   typename map<size_t,vector<double> >::iterator
                     dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.push_back( sc() );
                }
           }
         else if ( prop_key.type == VECTOR )
           {
              sgref.E( (*it).first )->Read( prop_key, vc );
              for ( typename vector<size_t>::const_iterator
                    pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   typename map<size_t,vector<double> >::iterator
                     dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.reserve(6U);
                   if ( dim == 3U )
                     for ( size_t j=0; j<3U; j++ ) (*dit2).second.push_back( vc[j] );
                   else {
                        (*dit2).second.push_back( vc[0] );
                        (*dit2).second.push_back( vc[1] );
                        (*dit2).second.push_back( 0. );
                     }
                }
             
           }
         else if ( prop_key.type == TENSOR )
           {
              sgref.E( (*it).first )->Read( prop_key, ts );
              for ( typename vector<size_t>::const_iterator
                    pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   typename map<size_t,vector<double> >::iterator
                     dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.reserve(12U);
                   if ( dim == 3U )
                     for ( size_t k=0; k<3U; k++ )
                       for ( size_t l=0; l<3U; l++ ) (*dit2).second.push_back( ts(k,l) );
                   else {
                      (*dit2).second.push_back( ts(0,0) );
                      (*dit2).second.push_back( ts(0,1) );
                      (*dit2).second.push_back( 0. );
                      (*dit2).second.push_back( ts(1,0) );
                      (*dit2).second.push_back( ts(1,1) );
                      (*dit2).second.push_back( 0. );
                      (*dit2).second.push_back( 0. );
                      (*dit2).second.push_back( 0. );
                      (*dit2).second.push_back( 0. );
                   }
                }
           }
      }
     
 } // end ElementData






/**
    Finds the element barycenters, storing these points in 'pxyz_data' as the
    locations associated with the 'node_nums' where the property values will be stored.
*/
template<size_t dim>
void VTK_Interface<dim>::ElementPointData( const Region<dim>& sgref,
                                           const csmp::Index&     prop_key,
                                           const map<size_t,size_t>& node_nums,
                                           map<size_t,vector<double> >& pxyz_data )
 {
    if ( prop_key.place != ELEMENT and prop_key.place != REGION )
      throw csmp::Exception( CSMP_FATAL_ERROR, "VTK_Interface<dim>::ElementPointData", 
                                   "Method only applies to element or region properties" );
    ScalarVariable      sc;
    VectorVariable<dim>  vc;
    TensorVariable<dim>  ts;
    
    pxyz_data.erase( pxyz_data.begin(), pxyz_data.end() );
    
    // 1. extracting the element-center coordinates and property values from the Model
    //    results are put into the first three elements of 'pxz_data'
    // -------------------------------------------------------------------------- 
    vector<double> empty_vec;
    for ( typename map<size_t,size_t>::const_iterator
          nit=node_nums.begin(); nit!=node_nums.end(); nit++ )
      {
         // inserting new element into the cordinate map using new node coordinate number
         pair<typename map<size_t,vector<double> >::iterator,bool>
           dit = pxyz_data.insert( make_pair( (*nit).first, empty_vec ) );
         assert( dit.second == true );
         
         // node coordinates (reserves storage for three coordinates and a scalar data value)
         (*dit.first).second.reserve(4);
         Point<dim>  bc(sgref.E( (*nit).second )->BaryCenter());
         (*dit.first).second.push_back( bc[0] );
         (*dit.first).second.push_back( bc[1] );
         if ( dim == 3 ) (*dit.first).second.push_back( bc[2] );
         else            (*dit.first).second.push_back( static_cast<double>(0.) );
      }
      
    // 2. Looping through the plist, associating the element properties with the nodal pxyz_data     
    // data values
    // -----------------------------------------------------------------------------------------
    for ( typename map<size_t,vector<size_t> >::const_iterator
          it=plist.begin(); it!=plist.end(); it++ )
      {
         if ( prop_key.type == SCALAR )
           {
              sgref.E( (*it).first )->Read( prop_key, sc );
              // looping over the element's node 
              for ( typename vector<size_t>::const_iterator
                    pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   typename map<size_t,vector<double> >::iterator dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.push_back( sc() );
                }
           }
         else if ( prop_key.type == VECTOR )
           {
              sgref.E( (*it).first )->Read( prop_key, vc );
              for ( typename vector<size_t>::const_iterator
                    pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   typename map<size_t,vector<double> >::iterator dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.reserve(6);
                   if ( dim == 3U )
                     for ( size_t j=0; j<3U; j++ ) (*dit2).second.push_back( vc[j] );
                   else {
                        (*dit2).second.push_back( vc[0] );
                        (*dit2).second.push_back( vc[1] );
                        (*dit2).second.push_back( static_cast<double>(0.) );
                     }
                }
             
           }
         else if ( prop_key.type == TENSOR )
           {
              sgref.E( (*it).first )->Read( prop_key, ts );
              for ( typename vector<size_t>::const_iterator
                    pit=(*it).second.begin(); pit!=(*it).second.end(); pit++ )
                {
                   typename map<size_t,vector<double> >::iterator dit2 = pxyz_data.find( (*pit) );
                   (*dit2).second.reserve(12);
                   if ( dim == 3U )
                     for ( size_t k=0U; k<3U; k++ )
                       for ( size_t l=0U; l<3U; l++ ) (*dit2).second.push_back( ts(k,l) );
                   else {
                      (*dit2).second.push_back( ts(0,0) );
                      (*dit2).second.push_back( ts(0,1) );
                      (*dit2).second.push_back( static_cast<double>(0.) );
                      (*dit2).second.push_back( ts(1,0) );
                      (*dit2).second.push_back( ts(1,1) );
                      (*dit2).second.push_back( static_cast<double>(0.) );
                      (*dit2).second.push_back( static_cast<double>(0.) );
                      (*dit2).second.push_back( static_cast<double>(0.) );
                      (*dit2).second.push_back( static_cast<double>(0.) );
                   }
                }
           }
      }
     
 } // end ElementPointData





/**
 
Depending on the element type, and the spatial dimensions of the model,
the VTK_Interface breaks the elements into simple triangles, quadrilaterals, 
tetrahedra, or hexahedra or pyramids for the the purpose of the visualization. 

@attention see "The Visualization Toolkit", >=4th ed., Cell types, p.129  

@test SKM 26/8/14, updated for nonlinear VTK>4.2 element types

*/
template<size_t dim>
void VTK_Interface<dim>::TransformPlist( const Region<dim>& sgref,
                                         const map<size_t,vector<size_t> >&  plist,
                                         deque<vector<size_t> >&  tdeque )
 {
    //                    triangle   tetrahedron  quadrilateral  hexahedron
    vector<size_t>     pentry(3), tentry(4),   qentry(4),     hentry(8);
 
    tdeque.erase( tdeque.begin(), tdeque.end() );
    geometric_primitives_VTK.erase( geometric_primitives_VTK.begin(), geometric_primitives_VTK.end() );
   
    for ( typename map<size_t,vector<size_t> >::const_iterator
          it=plist.begin(); it!=plist.end(); it++ )
      {
         switch ( sgref.E( (*it).first )->FE_Type() )
           {
              case LINEAR_BAR: // segment (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_LINE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_LINEAR_BAR: // segment (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_LINE);
                   tdeque.push_back( (*it).second );
                break;

              case ISOPARAMETRIC_QUADRATIC_BAR: // double segment
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_EDGE);
                   tdeque.push_back( (*it).second );
                break;
                
              case LINEAR_TRIANGLE: // triangle (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case LINEAR_TRIANGLE3D: // triangle (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_LINEAR_TRIANGLE: // triangle (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case LINEAR_TETRAHEDRON: // tetrahedron (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tdeque.push_back( (*it).second );
                break;
              
              case ISOPARAMETRIC_LINEAR_TETRAHEDRON: // tetrahedron (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tdeque.push_back( (*it).second );
                break;
              
              case QUADRATIC_TRIANGLE:
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_QUADRATIC_TRIANGLE: // quadratic triangle
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_TRIANGLE); // (=22)
                   tdeque.push_back( (*it).second );
                break;

              // there is not matching VTK element type, hence conversion into 6 triangles
              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE: // quadratic barycentric triangle (bubble node)
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_BIQUADRATIC_TRIANGLE); // (=32)
                   tdeque.push_back( (*it).second );
                break;

              case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: // quadratic tetrahedron -> 13 tetrahedra
                   // ANSYS 10-noded tetrahedron
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_TETRA); // (=24)
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON: // barycentric quadratic tetrahedron -> 13 tetrahedra
                   // basal 3 outer tetrahedra
                   // tetrahedron 1
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[0]; tentry[1]=(*it).second[4]; tentry[2]=(*it).second[6]; tentry[3]=(*it).second[7];
                   tdeque.push_back( tentry );
                   // tetrahedron 2
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[1]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[4]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                   // tetrahedron 3
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[2]; tentry[1]=(*it).second[6]; tentry[2]=(*it).second[5]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // top 3 tetrahedra
                   // tetrahedron 4
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[7]; tentry[1]=(*it).second[8]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // tetrahedron 5
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[8]; tentry[1]=(*it).second[9]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // tetrahedron 6
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[7]; tentry[1]=(*it).second[10]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // intermediate 3 side tetrahedra
                   // tetrahedron 7
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[7]; tentry[2]=(*it).second[8]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // tetrahedron 8
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[5]; tentry[1]=(*it).second[8]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // tetrahedron 9
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[6]; tentry[1]=(*it).second[7]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // oblique 3 tetrahedra in the lower part of the parent tetrahedron
                   // tetrahedron 10
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[6]; tentry[2]=(*it).second[7]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // tetrahedron 11
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                   // tetrahedron 12
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[6]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // lower plane basal tetrahedron
                   // tetrahedron 13
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[6]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                break;

              case ISOPARAMETRIC_LINEAR_QUADRILATERAL: // (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   tdeque.push_back( (*it).second );
                break;

              case ISOPARAMETRIC_LINEAR_PYRAMID: // (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_PYRAMID);
                   tdeque.push_back( (*it).second );
                break;

              case ISOPARAMETRIC_LINEAR_PRISM: // Wedge (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_WEDGE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_LINEAR_HEXAHEDRON: // linear hexahedron 
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   tdeque.push_back( (*it).second );
                break;
               
              case ISOPARAMETRIC_QUADRATIC_PRISM15:
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_WEDGE);
                   tdeque.push_back( (*it).second );
                break;
               
              case ISOPARAMETRIC_QUADRATIC_PYRAMID13:
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_PYRAMID);
                   tdeque.push_back( (*it).second );
                break;

              case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20: // quadratic hexahedron
                   // like the 20-noded ANSYS hexahedron
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_HEXAHEDRON); // (=25)
                   tdeque.push_back( (*it).second );
                 break;
               
              case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: // quadratic hexahedron
                   //Linear hex 1 => 0,8,20,11,   12,22,26,24
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[0]; hentry[1] = (*it).second[8]; hentry[2] = (*it).second[20];
                   hentry[3] = (*it).second[11]; hentry[4] = (*it).second[12]; hentry[5] = (*it).second[21];
                   hentry[6] = (*it).second[26]; hentry[7] = (*it).second[24];
                   tdeque.push_back( hentry );
                   
                    //Linear hex 2 => 8,1,9,20,   21,13,22,26
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[8]; hentry[1] = (*it).second[1]; hentry[2] = (*it).second[9];
                   hentry[3] = (*it).second[20]; hentry[4] = (*it).second[21]; hentry[5] = (*it).second[13];
                   hentry[6] = (*it).second[22]; hentry[7] = (*it).second[26];
                   tdeque.push_back( hentry );
                
                   //Linear hex 3 =>  20,9,2,10  26,22,14,23 
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[20]; hentry[1] = (*it).second[9]; hentry[2] = (*it).second[2];
                   hentry[3] = (*it).second[10]; hentry[4] = (*it).second[26]; hentry[5] = (*it).second[22];
                   hentry[6] = (*it).second[14]; hentry[7] = (*it).second[23];
                   tdeque.push_back( hentry );
                
                   //Linear hex 4 =>  11,20,10,3  24,26,23,15
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[11]; hentry[1] = (*it).second[20]; hentry[2] = (*it).second[10];
                   hentry[3] = (*it).second[3]; hentry[4] = (*it).second[24]; hentry[5] = (*it).second[26];
                   hentry[6] = (*it).second[23]; hentry[7] = (*it).second[15];
                   tdeque.push_back( hentry );
                   
                   // Second layer of hexes
                   //Linear hex 5 =>  12,21,26,24  4,16,25,19
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[12]; hentry[1] = (*it).second[21]; hentry[2] = (*it).second[26];
                   hentry[3] = (*it).second[24]; hentry[4] = (*it).second[4]; hentry[5] = (*it).second[16];
                   hentry[6] = (*it).second[25]; hentry[7] = (*it).second[19];
                   tdeque.push_back( hentry );
                   
                   //Linear hex 6 => 21,13,22,26,  16,5,17,25,
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[21]; hentry[1] = (*it).second[13]; hentry[2] = (*it).second[22];
                   hentry[3] = (*it).second[26]; hentry[4] = (*it).second[16]; hentry[5] = (*it).second[5];
                   hentry[6] = (*it).second[17]; hentry[7] = (*it).second[25];
                   tdeque.push_back( hentry );
                   
                   //Linear hex 7 => 26,22,14,23,  25,17,6,18
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[26]; hentry[1] = (*it).second[22]; hentry[2] = (*it).second[14];
                   hentry[3] = (*it).second[23]; hentry[4] = (*it).second[25]; hentry[5] = (*it).second[17];
                   hentry[6] = (*it).second[6]; hentry[7] = (*it).second[18];
                   tdeque.push_back( hentry );
                   
                   //Linear hex 8 => 24,26,23,15  19,25,18,7
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[24]; hentry[1] = (*it).second[26]; hentry[2] = (*it).second[23];
                   hentry[3] = (*it).second[15]; hentry[4] = (*it).second[19]; hentry[5] = (*it).second[25];
                   hentry[6] = (*it).second[18]; hentry[7] = (*it).second[7];
                   tdeque.push_back( hentry );
                break;

              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: // 8-noded quadratic quadrilateral
                   geometric_primitives_VTK.push_back(VTK_QUADRATIC_QUAD); // (=23)
                   tdeque.push_back( (*it).second  );
                break;
                
              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9: // quadratic quadrilateral
                   geometric_primitives_VTK.push_back(VTK_BIQUADRATIC_QUAD);
                   tdeque.push_back( (*it).second  );
                break;

              default:
                   cerr <<"\nElement Idx: "<< sgref.E( (*it).first )->Idx();
                   cerr <<" with "<< sgref.E( (*it).first )->Nodes() <<" nodes."<< endl;
                   throw csmp::Exception( CSMP_FATAL_ERROR, "VTK_Interface<dim>::TransformPlist", 
                                         "Unable to interpret how this element shall be broken in subelements");
           }
      }
 
 } // end TransformPlist





/** OLD STUFF (CONVERSION OF QUADRATIC ELEMENTS TO MULTIPLE LINEAR ELEMENTS

              case QUADRATIC_TRIANGLE: // quadratic triangle -> 4 triangles
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[3]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                   // triangle 2
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[4];
                   tdeque.push_back( pentry );
                   // triangle 3
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[5]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 4
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                break;

              case ISOPARAMETRIC_QUADRATIC_TRIANGLE: // quadratic triangle -> 4 triangles
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[3]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                   // triangle 2
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[4];
                   tdeque.push_back( pentry );
                   // triangle 3
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[5]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 4
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                break;

              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE: // quadratic barycentric triangle -> 6 triangles
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[3]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // triangle 2
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // triangle 3
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[6]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[4];
                   tdeque.push_back( pentry );
                   // triangle 4
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[6]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 5
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[5]; pentry[1] = (*it).second[6]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 6
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[6]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                break;

              case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: // quadratic tetrahedron -> 13 tetrahedra
                   // 4 corner tetrahedra
                   // tetrahedron 1
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[1]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[4]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                   // tetrahedron 2
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[2]; tentry[1]=(*it).second[6]; tentry[2]=(*it).second[5]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // tetrahedron 3
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[0]; tentry[1]=(*it).second[4]; tentry[2]=(*it).second[6]; tentry[3]=(*it).second[7];
                   tdeque.push_back( tentry );
                   // tetrahedron 4
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[7]; tentry[1]=(*it).second[8]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // 4 tetrahedra which make up central octahedron
                   // tetrahedron 5
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                   // tetrahedron 6
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[6]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // tetrahedron 7
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[6]; tentry[2]=(*it).second[7]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // tetrahedron 8
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[9]; tentry[2]=(*it).second[7]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                break;
                
              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON: // barycentric quadratic tetrahedron -> 13 tetrahedra
                   // basal 3 outer tetrahedra
                   // tetrahedron 1
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[0]; tentry[1]=(*it).second[4]; tentry[2]=(*it).second[6]; tentry[3]=(*it).second[7];
                   tdeque.push_back( tentry );
                   // tetrahedron 2
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[1]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[4]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                   // tetrahedron 3
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[2]; tentry[1]=(*it).second[6]; tentry[2]=(*it).second[5]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // top 3 tetrahedra
                   // tetrahedron 4
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[7]; tentry[1]=(*it).second[8]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // tetrahedron 5
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[8]; tentry[1]=(*it).second[9]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // tetrahedron 6
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[7]; tentry[1]=(*it).second[10]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[3];
                   tdeque.push_back( tentry );
                   // intermediate 3 side tetrahedra
                   // tetrahedron 7
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[7]; tentry[2]=(*it).second[8]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // tetrahedron 8
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[5]; tentry[1]=(*it).second[8]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // tetrahedron 9
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[6]; tentry[1]=(*it).second[7]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[9];
                   tdeque.push_back( tentry );
                   // oblique 3 tetrahedra in the lower part of the parent tetrahedron
                   // tetrahedron 10
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[6]; tentry[2]=(*it).second[7]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // tetrahedron 11
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[10]; tentry[3]=(*it).second[8];
                   tdeque.push_back( tentry );
                   // tetrahedron 12
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[6]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[9]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                   // lower plane basal tetrahedron
                   // tetrahedron 13
                   geometric_primitives_VTK.push_back(VTK_TETRA);
                   tentry[0]=(*it).second[4]; tentry[1]=(*it).second[5]; tentry[2]=(*it).second[6]; tentry[3]=(*it).second[10];
                   tdeque.push_back( tentry );
                break;

              case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: // quadratic hexahedron 
                   //Linear hex 1 => 0,8,20,11,   12,22,26,24
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[0]; hentry[1] = (*it).second[8]; hentry[2] = (*it).second[20];
                   hentry[3] = (*it).second[11]; hentry[4] = (*it).second[12]; hentry[5] = (*it).second[21];
                   hentry[6] = (*it).second[26]; hentry[7] = (*it).second[24];
                   tdeque.push_back( hentry );
                   
                    //Linear hex 2 => 8,1,9,20,   21,13,22,26
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[8]; hentry[1] = (*it).second[1]; hentry[2] = (*it).second[9];
                   hentry[3] = (*it).second[20]; hentry[4] = (*it).second[21]; hentry[5] = (*it).second[13];
                   hentry[6] = (*it).second[22]; hentry[7] = (*it).second[26];
                   tdeque.push_back( hentry );
                
                   //Linear hex 3 =>  20,9,2,10  26,22,14,23 
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[20]; hentry[1] = (*it).second[9]; hentry[2] = (*it).second[2];
                   hentry[3] = (*it).second[10]; hentry[4] = (*it).second[26]; hentry[5] = (*it).second[22];
                   hentry[6] = (*it).second[14]; hentry[7] = (*it).second[23];
                   tdeque.push_back( hentry );
                
                   //Linear hex 4 =>  11,20,10,3  24,26,23,15
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[11]; hentry[1] = (*it).second[20]; hentry[2] = (*it).second[10];
                   hentry[3] = (*it).second[3]; hentry[4] = (*it).second[24]; hentry[5] = (*it).second[26];
                   hentry[6] = (*it).second[23]; hentry[7] = (*it).second[15];
                   tdeque.push_back( hentry );
                   
                   // Second layer of hexes
                   //Linear hex 5 =>  12,21,26,24  4,16,25,19
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[12]; hentry[1] = (*it).second[21]; hentry[2] = (*it).second[26];
                   hentry[3] = (*it).second[24]; hentry[4] = (*it).second[4]; hentry[5] = (*it).second[16];
                   hentry[6] = (*it).second[25]; hentry[7] = (*it).second[19];
                   tdeque.push_back( hentry );
                   
                   //Linear hex 6 => 21,13,22,26,  16,5,17,25,
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[21]; hentry[1] = (*it).second[13]; hentry[2] = (*it).second[22];
                   hentry[3] = (*it).second[26]; hentry[4] = (*it).second[16]; hentry[5] = (*it).second[5];
                   hentry[6] = (*it).second[17]; hentry[7] = (*it).second[25];
                   tdeque.push_back( hentry );
                   
                   //Linear hex 7 => 26,22,14,23,  25,17,6,18
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[26]; hentry[1] = (*it).second[22]; hentry[2] = (*it).second[14];
                   hentry[3] = (*it).second[23]; hentry[4] = (*it).second[25]; hentry[5] = (*it).second[17];
                   hentry[6] = (*it).second[6]; hentry[7] = (*it).second[18];
                   tdeque.push_back( hentry );
                   
                   //Linear hex 8 => 24,26,23,15  19,25,18,7
                   geometric_primitives_VTK.push_back(VTK_HEXAHEDRON);
                   hentry[0] = (*it).second[24]; hentry[1] = (*it).second[26]; hentry[2] = (*it).second[23];
                   hentry[3] = (*it).second[15]; hentry[4] = (*it).second[19]; hentry[5] = (*it).second[25];
                   hentry[6] = (*it).second[18]; hentry[7] = (*it).second[7];
                   tdeque.push_back( hentry );
                break;

              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: // quadratic quadrilateral is split into 4 linear quads
           		// quad 1
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[0]; qentry[1] = (*it).second[4]; qentry[2] = (*it).second[8]; qentry[3] = (*it).second[7];
                   tdeque.push_back( qentry );
                   // quad 2
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[4]; qentry[1] = (*it).second[1]; qentry[2] = (*it).second[5]; qentry[3] = (*it).second[8];
                   tdeque.push_back( qentry );
                   // quad 3
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[8]; qentry[1] = (*it).second[5]; qentry[2] = (*it).second[2]; qentry[3] = (*it).second[6];
                   tdeque.push_back( qentry );
                   // quad 4
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[7]; qentry[1] = (*it).second[8]; qentry[2] = (*it).second[6]; qentry[3] = (*it).second[3];
                   tdeque.push_back( qentry );
                break;
                
              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9: // quadratic quadrilateral -> into 6 quadrilaterals
                   // quadrilateral 1
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[0]; qentry[1] = (*it).second[4]; qentry[2] = (*it).second[8]; qentry[3] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // quadrilateral 2
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[4]; qentry[1] = (*it).second[1]; qentry[2] = (*it).second[5]; qentry[3] = (*it).second[8];
                   tdeque.push_back( pentry );
                   // quadrilateral 3
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   pentry[0] = (*it).second[8]; qentry[1] = (*it).second[5]; qentry[2] = (*it).second[2]; qentry[3] = (*it).second[7];
                   tdeque.push_back( pentry );
                   // quadrilateral 4
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[6]; qentry[1] = (*it).second[8]; qentry[2] = (*it).second[7]; qentry[3] = (*it).second[3];
                   tdeque.push_back( pentry );
                break;
 
 END OLD STUFF*/





/*
template<size_t dim>
void VTK_Interface<dim>::TransformFacePlist( const Model<dim>&                 sg,
                                                const map<size_t,vector<size_t> >&  plist,
                                                deque<vector<size_t> >&        tdeque )

 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
    typename map<size_t,vector<size_t> >::const_iterator  it;
    vector<size_t>  pentry(3), lentry(2), qentry(4);
    bool               first_call(true);
 
    tdeque.erase( tdeque.begin(), tdeque.end() );
    geometric_primitives_VTK.erase( geometric_primitives_VTK.begin(), geometric_primitives_VTK.end() );

    for ( it=plist.begin(); it!=plist.end(); it++ )
      {
         switch ( sg.Mesh().F( (*it).first-1 )->FE_Type() )
           {
              case LINEAR_BAR: // segment (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_LINE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_LINEAR_BAR: // segment (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_LINE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_QUADRATIC_BAR: // double segment
                   geometric_primitives_VTK.push_back(VTK_LINE);
                   lentry[0] = (*it).second[0]; lentry[1] = (*it).second[2];  
                   tdeque.push_back( lentry );
                   geometric_primitives_VTK.push_back(VTK_LINE);
                   lentry[0] = (*it).second[2]; lentry[1] = (*it).second[1];  
                   tdeque.push_back( lentry );
                break;
                
              case ISOPARAMETRIC_LINEAR_TRIANGLE: // triangle (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case LINEAR_TRIANGLE: // triangle (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case LINEAR_TRIANGLE3D: // triangle (vectors are just copied over)
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   tdeque.push_back( (*it).second );
                break;
                
              case ISOPARAMETRIC_LINEAR_QUADRILATERAL: // quadrilateral
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   tdeque.push_back( (*it).second );
                break;
              
              case ISOPARAMETRIC_QUADRATIC_TRIANGLE: // quadratic triangle -> 4 triangles
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[3]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                   // triangle 2
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[4];
                   tdeque.push_back( pentry );
                   // triangle 3
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[5]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 4
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                break;

              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE: // quadratic barycentric triangle -> 6 triangles
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[3]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // triangle 2
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[3]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // triangle 3
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[6]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[4];
                   tdeque.push_back( pentry );
                   // triangle 4
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[6]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 5
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[5]; pentry[1] = (*it).second[6]; pentry[2] = (*it).second[2];
                   tdeque.push_back( pentry );
                   // triangle 6
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[6]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                break;
                
              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: // quadratic quadrilateral -> broken into 6 triangles
                   // triangle 1
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[0]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[7];
                   tdeque.push_back( pentry );
                   // triangle 2
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[4]; pentry[1] = (*it).second[1]; pentry[2] = (*it).second[5];
                   tdeque.push_back( pentry );
                   // triangle 3
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[5]; pentry[1] = (*it).second[2]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // triangle 4
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[7]; pentry[1] = (*it).second[6]; pentry[2] = (*it).second[3];
                   tdeque.push_back( pentry );
                   // triangle 5
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[7]; pentry[1] = (*it).second[4]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // triangle 6
                   geometric_primitives_VTK.push_back(VTK_TRIANGLE);
                   pentry[0] = (*it).second[4]; pentry[1] = (*it).second[5]; pentry[2] = (*it).second[6];
                   tdeque.push_back( pentry );
                break;
              
              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9: // quadratic quadrilateral -> into 6 quadrilaterals
                   // quadrilateral 1
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[0]; qentry[1] = (*it).second[4]; qentry[2] = (*it).second[8]; qentry[3] = (*it).second[6];
                   tdeque.push_back( pentry );
                   // quadrilateral 2
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[4]; qentry[1] = (*it).second[1]; qentry[2] = (*it).second[5]; qentry[3] = (*it).second[8];
                   tdeque.push_back( pentry );
                   // quadrilateral 3
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   pentry[0] = (*it).second[8]; qentry[1] = (*it).second[5]; qentry[2] = (*it).second[2]; qentry[3] = (*it).second[7];
                   tdeque.push_back( pentry );
                   // quadrilateral 4
                   geometric_primitives_VTK.push_back(VTK_QUAD);
                   qentry[0] = (*it).second[6]; qentry[1] = (*it).second[8]; qentry[2] = (*it).second[7]; qentry[3] = (*it).second[3];
                   tdeque.push_back( pentry );
                break;

              default:
                   cout <<"\nFace ID: "<< sg.Mesh().F( (*it).first-1 )->ID();
                   cout <<" with "<< sgref.E( (*it).first-1 )->Nodes() <<" nodes."<< endl;
                   throw csmp::Exception( CSMP_FATAL_ERROR, "VTK_Interface<dim>::TransformFacePlist", 
                                  "Unable to interpret how this face shall be broken in subelements");
           }
      }
 
 } // end TransformFacePlist
*/





/**
   the pxyz_data record is resized to the number of cells, no coordinates are stored
   only the cell centered variable values are recorded
   
   SKM fixed 1/9/2014 to include new VTK quadratic element types
*/ 
template<size_t dim>
void VTK_Interface<dim>::CellData( const Region<dim>& sgref,
                                   const csmp::Index& prop_key,
                                   const map<size_t,vector<size_t> >& plist,
                                   map<size_t,vector<double> >& cell_data )
 {
    if ( prop_key.place != ELEMENT  and  prop_key.place != INTER_FACE and prop_key.place != REGION )
      throw csmp::Exception( CSMP_FATAL_ERROR, "VTK_Interface<dim>::CellData", 
                                   "Method only applies to face, element or region properties" );
    
    cell_data.erase( cell_data.begin(), cell_data.end() );
      
    // 1. Looping over the plist, recording the desired data values
    // -----------------------------------------------------------------------------------------
    size_t   n(1U);
//?    size_t   n(0U);
    vector<double>  empty_vec;

    for ( typename map<size_t,vector<size_t> >::const_iterator 
          it=plist.begin(); it!=plist.end(); it++ )
      {
         // calculating into how many cells the element will be split up since the cell-data entry
         // must be multiplied by a corresponding number
         size_t  fragments(1U); 
         switch ( sgref.E( (*it).first )->FE_Type() )
           {
//              case ISOPARAMETRIC_QUADRATIC_BAR: fragments = 2U; // double segment
//                break;
//              case ISOPARAMETRIC_QUADRATIC_TRIANGLE: fragments = 4U; // quadratic triangle -> 4 triangles
//                break;
              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TRIANGLE: fragments = 6U; // quadratic barycentric triangle -> 6 triangles
                break;
//              case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON: fragments = 13U; // quadratic tetrahedron -> 13 tetrahedra
//                break;
              case ISOPARAMETRIC_BARYCENTRIC_QUADRATIC_TETRAHEDRON: fragments = 13U; // barycentric quadratic tetrahedron -> 13 tetrahedra
                break;
              case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON27: fragments = 8U; // quadratic hexahedron -> 8 hexahedra
                break;
//              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL: fragments = 4U; // quadratic quadrilateral is split into 4 linear quads
//                break;
              case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL9: fragments = 4U; // quadratic quadrilateral -> into 4 quadrilaterals
                break;
              default:
                // default=1: works for LINEAR_BAR, ISOPARAMETRIC_LINEAR_BAR, LINEAR_TRIANGLE, LINEAR_TRIANGLE3D
                // ISOPARAMETRIC_LINEAR_TRIANGLE, LINEAR_TETRAHEDRON, ISOPARAMETRIC_LINEAR_TETRAHEDRON,
                // ISOPARAMETRIC_LINEAR_QUADRILATERAL, ISOPARAMETRIC_LINEAR_PYRAMID, -WEDGE, -HEXAHEDRON 
                fragments = 1U;
           }
       
         // creating new records in 'cell_data'
         for ( size_t i=0U; i<fragments; i++ ) {
              pair<typename map<size_t,vector<double> >::iterator,bool> dit =
                                                   cell_data.insert( make_pair( n++, empty_vec ) );
              assert( dit.second == true );
         
              // initializing the data members of this record
              if ( prop_key.type == SCALAR ) {
                   (*dit.first).second.reserve(1U);
                   (*dit.first).second.push_back( sgref.E( (*it).first )->Read( prop_key ) );
                }
	            else if ( prop_key.type == VECTOR ) {
                  (*dit.first).second.resize( 3U, 0. ); // initialization of empty vector<double> to zero
	              VectorVariable<dim>  vc;
	              sgref.E( (*it).first )->Read( prop_key, vc );
                  for ( size_t j=0U; j<dim; j++ ) (*dit.first).second[j] = vc[j];
	              }
	            else if ( prop_key.type == TENSOR ) {
	              TensorVariable<dim>  ts;
	              sgref.E( (*it).first )->Read( prop_key, ts );
	              (*dit.first).second.reserve(9U);
	              if ( dim == 3U ) {
                       for ( size_t k=0U; k<dim; k++ )
	                     for ( size_t l=0U; l<dim; l++ ) (*dit.first).second.push_back( ts(k,l) );
	                }
	              else { // dim == 2 (since interface does not work for 1D elements)
	                   (*dit.first).second.push_back( ts(0,0) );
	                   (*dit.first).second.push_back( ts(0,1) );
                       (*dit.first).second.push_back( static_cast<double>(0.) );
	                   (*dit.first).second.push_back( ts(1,0) );
	                   (*dit.first).second.push_back( ts(1,1) );
                       (*dit.first).second.push_back( static_cast<double>(0.) );
                       (*dit.first).second.push_back( static_cast<double>(0.) );
                       (*dit.first).second.push_back( static_cast<double>(0.) );
                       (*dit.first).second.push_back( static_cast<double>(0.) );
	                }
	           }
	       }
      }
     
 } // end CellData




/**
    Outputs all subregion of model, with file_name preface and var-name attached
*/
template<size_t dim>
template<class T>
void VTK_Interface<dim>::OutputRegionByRegionToVTK( const Model<dim>& model,
                                                    const std::string& initial_file_name,
                                                    const std::string& var_name,
                                                    T timestep,
                                                    bool update_topology )
 {
    std::string file_name( this->OutputFileAndSubFolderName( initial_file_name ) );

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     
    if ( !model.Database().IsDefined(var_name.c_str()) ) {
         csmp_error.notice( CSMP_ERROR,
                            "VTK_Interface<dim>::OutputRegionByRegionToVTK",
                            var_name.c_str(),
                            "output property is not defined. Nothing was done.");
         return;
      }
    csmp::Index prop_key = model.Database().StorageKey(var_name.c_str());
    
    // unique regions first (model excluded)
    for ( typename std::map<std::string,csmp::Region<dim> >::const_iterator
          it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); it++ )
      OutputDataToVTK( model, (*it).first.c_str(), file_name, var_name, timestep, update_topology );

    // non-unique regions, if any
    for ( typename std::map<std::string,csmp::Region<dim> >::const_iterator
          it=model.RegionsBegin(); it!=model.RegionsEnd(); it++ )
      OutputDataToVTK( model, (*it).first.c_str(), file_name, var_name, timestep, update_topology );
 
 } // end OutputRegionByRegionToVTK

template void VTK_Interface<1U>::OutputRegionByRegionToVTK(const Model<1U>&,const std::string&,const std::string&,int,bool);
template void VTK_Interface<2U>::OutputRegionByRegionToVTK(const Model<2U>&,const std::string&,const std::string&,int,bool);
template void VTK_Interface<3U>::OutputRegionByRegionToVTK(const Model<3U>&,const std::string&,const std::string&,int,bool);

template void VTK_Interface<1U>::OutputRegionByRegionToVTK(const Model<1U>&,const std::string&,const std::string&,long,bool);
template void VTK_Interface<2U>::OutputRegionByRegionToVTK(const Model<2U>&,const std::string&,const std::string&,long,bool);
template void VTK_Interface<3U>::OutputRegionByRegionToVTK(const Model<3U>&,const std::string&,const std::string&,long,bool);

template void VTK_Interface<1U>::OutputRegionByRegionToVTK(const Model<1U>&,const std::string&,const std::string&,size_t,bool);
template void VTK_Interface<2U>::OutputRegionByRegionToVTK(const Model<2U>&,const std::string&,const std::string&,size_t,bool);
template void VTK_Interface<3U>::OutputRegionByRegionToVTK(const Model<3U>&,const std::string&,const std::string&,size_t,bool);

template void VTK_Interface<1U>::OutputRegionByRegionToVTK(const Model<1U>&,const std::string&,const std::string&,double,bool);
template void VTK_Interface<2U>::OutputRegionByRegionToVTK(const Model<2U>&,const std::string&,const std::string&,double,bool);
template void VTK_Interface<3U>::OutputRegionByRegionToVTK(const Model<3U>&,const std::string&,const std::string&,double,bool);


// METHODS FOR INDIVIDUAL ELEMENTS


/**
    For the user-defined discretised variable identified by key, method prints 
    the values stored on the given element to file.
*/
template<size_t dim>
void outputNodeDataToVTK( const Element<dim>& e, const csmp::Index& key,
                          const char* file_name, const char* variable_name )
  {
     if ( key.place != NODE ) {
          ErrorHandler::Instance().notice( CSMP_ERROR, "outputNodeDataToVTK:", "output variable must be placed on the node; nothing was done." );
          return;
       }

    string var_name(variable_name);
    replaceWhiteSpaceBy( var_name, '_' );

     string  outfile(file_name);
     // adding the element number
     outfile += to_string( e.Idx() );
     // adding the element type
     outfile += "_";
     outfile += parseFiniteElementType( e.FE()->ElementType() );
     outfile +=".vtk";

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs ) {
          ErrorHandler::Instance().notice( CSMP_ERROR, "outputNodeDataToVTK:", "output file could not be opened; nothing was done." );
          return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD;
     e.NodeCoordinateMatrix( COORD );
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     // NOTE: in VTK coordinates always are stored in single precision
     ofs <<"POINTS " << e.Nodes() <<" float"<< endl;
     for ( size_t i=0; i<e.Nodes(); i++ )
       {
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          if ( dim == 2 ) ofs << 0.;
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< e.Nodes()+1 << endl;
     ofs << e.Nodes() <<" ";
     for ( size_t i=0; i<e.Nodes(); i++ ) ofs << i <<" ";
     ofs << endl;

     // 4. writing CELL_TYPES - for the
     // ---------------------
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << parseElementType( e.FE_Type() ) << endl;
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     // Unfortunately the data can only be output as nodal variables
     ofs <<"POINT_DATA "<< e.Nodes() << endl;
     ofs.setf( ios::scientific );

     if ( key.type == SCALAR )
       {
           ofs <<"SCALARS "<< var_name <<" double"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA 1 x nodes
           for ( size_t i=0; i<e.Nodes(); i++ ) ofs << e.N(i)->Read(key) <<" ";
           ofs << endl;
       }
     else if ( key.type == VECTOR )
       {
          ofs <<"VECTORS "<< var_name <<" double"<< endl;
          vector<VectorVariable<dim> > data;
          e.NodePropertyVector( key, data );
          // matrix DATA is vec-dim x nodes
          for ( size_t i=0; i<e.Nodes(); i++ ) {
               for ( size_t j=0; j<dim; j++ ) ofs << data[i](j) <<" ";
               // variables have always 3 components since view screen is 3D
               if ( dim == 2 ) ofs << 0.;
               ofs << endl;
            }
       }
     else if ( key.type == TENSOR )
       {
          ofs <<"TENSORS "<< var_name <<" double"<< endl;
          vector<TensorVariable<dim> > data;
          e.NodePropertyVector( key, data );
          // matrix DATA is vec-dim x nodes
          for ( size_t i=0; i<e.Nodes(); i++ ) {
               for ( size_t j=0; j<dim; j++ ) {
                    for ( size_t k=0; k<dim; ++k ) {
                          ofs << data[i](j,k) <<" ";
                          if ( dim == 2 ) ofs << 0. <<" ";
                      }
                    if ( dim == 2 ) ofs << 0. <<" ";
                 }
               // variables have always 3 components since view screen is 3D
               if ( dim == 2 ) ofs << 0.;
               ofs << endl;
            }
       }
     else ErrorHandler::Instance().notice( CSMP_ERROR, "outputNodeDataToVTK:", parseType(key.type),
                                          "node property placement cannot be output; nothing was done." );
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricLinearTetrahedron::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

 } // end outputNodeDataToVTK

template void outputNodeDataToVTK( const Element<2>&, const csmp::Index&, const char*, const char* );
template void outputNodeDataToVTK( const Element<3>&, const csmp::Index&, const char*, const char* );



/**
 
Outputs supplied nodal values of a variable into an output file which
can be visualized with the visualization toolkit VTK. The file is in
text format. The ID of the parent element as well as the extension ".vtk" 
will be appended to the output file name.  

@section arguments Input Arguments 

The name of the data file and the name of the variable as it is later 
used in VTK must be supplied as output arguments. The third argument is
a data matrix which will contain n-colums=nodes and n-rows = dimensions
of data (n-rows=1 = scalar data, n-rows=3 = vector<double> data etc.) of
nodal variable data.  

@section application Application

Single elements are output to VTK in order to test the 
quality of interpolation and the computed properties directly.  

@section messages Messages 

The method will indicate if there is a problem in opening the output
file.  
*/
template<size_t dim>
void outputNodeDataToVTK( const Element<dim>& e,
                          const std::string& file_name, const std::string& var_name,
                          DenseMatrix<DM_MIN>& DATA ) 
  {
     if ( DATA.Rows() > dim )
       throw csmp::Exception( CSMP_FATAL_ERROR, "Element<dim>::OutputNodeDataToVTK", 
                      "Node DATA matrix can have at most model-dimensions rows" );

     if ( DATA.Cols() != e.Nodes() )
       throw csmp::Exception( CSMP_FATAL_ERROR, "Element<dim>::OutputNodeDataToVTK", 
                      "Node DATA matrix must have as many columns as the element has nodes" );
     
     e.CoordinateMatrix();
     e.FE()->OutputNodeDataToVTK( file_name, var_name, DATA );
  } 



template class VTK_Interface<2U>;
template class VTK_Interface<3U>;




/**
    Loops over the faces of the element, collecting the (outward-pointing) face normals
    for display in VTK.
    Each normal is scaled by the sqrt of the element volume and placed 
    with its origin on the barycenter of the corresponding face.
    
    method calls UnitNormalToFace(size_t face) to compute the unit normals.
    
    @author SKM 5/7/14
*/
template<size_t dim>
void outputFaceNormalsToVTK( const Element<dim>& e, const char* file )
 {
    vector<double64>    nrml;
    size_t              counter(0);
    map<size_t,pair<Point<dim>,Point<dim> > >  normals;

    // 0. generating face normals, scaling and storing them
    // -----------------------------------------------------------
    vector<size_t> fnids;
    for ( size_t i=0; i<e.Faces(); i++ ) {
         // computing and storing normal to face
         Point<dim> unormal = e.UnitNormalToFace( i );
         // scaling the normals (by empirical factor)
         unormal *= (sqrt(e.Volume()) / 10.);
         // finding root points for the normals = barycenters of faces
         Point<dim> bctr = e.FaceBaryCenter( i );
         // storing starting pointsd and normals in map
         normals[ counter++ ] = make_pair( bctr, unormal );
      }  
   
     // 1. writing the file header
     // --------------------------
     string  outfile(file);
     // adding the element number
     outfile += to_string( e.Idx() );
     // adding the element type
     outfile += "_";
     outfile += parseFiniteElementType( e.FE()->ElementType() );
     outfile +=".vtk";
     ofstream  ofs(outfile.c_str());
     ofs <<"# vtk DataFile Version 4.2"<< endl;
     ofs <<"normals to faces of element: '";
     ofs << "normals" <<"'"<< endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3. writing point coordinates in order
     // --------------------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << e.Faces() <<" float"<< endl;
     for ( typename map<size_t,pair<Point<dim>,Point<dim> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) {
          assert( dim != 1 );
          if ( dim == 2 ) ofs << (*it).second.first[0] <<" "<<  (*it).second.first[1] <<" "<<  0. << endl;
          else ofs << (*it).second.first[0] <<" "<<  (*it).second.first[1] <<" "<<  (*it).second.first[2] << endl;
       }
     ofs << endl;  
       
     // 4. writing CELLS (cell-size and member face (points))
     // -----------------------------------------------------
     // (each normal=Edge has two points plus a specifier)
     ofs <<"CELLS "<< e.Faces() <<" "<< (e.Faces() * 2) << endl;
     counter=0;
     for ( typename map<size_t,pair<Point<dim>,Point<dim> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) {
          ofs <<"1 "<< counter++ << endl;
       }
     ofs << endl;
     
     // 5. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< normals.size() << endl;
     // cell type VTK_VERTEX (=1)
     for ( typename map<size_t,pair<Point<dim>,Point<dim> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) ofs << 1 << endl;
     ofs << endl;

     // 6. writing POINT_DATA point-type data values
     // ----------------------------------------------------
     ofs <<"POINT_DATA "<< normals.size() << endl;
     ofs <<"VECTORS "<< "normals" <<" double"<< endl;
     for ( typename map<size_t,pair<Point<dim>,Point<dim> > >::const_iterator
           it=normals.begin(); it!=normals.end(); it++ ) {
          assert( dim != 1 );
          if ( dim == 2 ) ofs << (*it).second.second[0] <<" "<<  (*it).second.second[1] <<" "<<  0. << endl;
          else ofs << (*it).second.second[0] <<" "<<  (*it).second.second[1] <<" "<<  (*it).second.second[2] << endl;
       }
     ofs << endl;
     cout <<"\nVTK file '"<< outfile <<"' written successfully."<< endl;
      
 } // end outputFaceNormalsToVTK


template void outputFaceNormalsToVTK( const Element<2>&, const char* );
template void outputFaceNormalsToVTK( const Element<3>&, const char* );






/**
    The integration points are coloured by number and are then output to VTK file.
    Point 1=0, point2=1...n-1.
*/
template<size_t dim>
void outputIntegrationPointsToVTK( const Element<dim>& element, const char* file )
 {
     // 0. generating values, first integration point=0, second=1...
     // ------------------------------------------------------------
    vector<pair<Point<dim>,double64> >  ipoint_data(element.IntegrationPoints()); // locations in physical space
    // convention: integration points get values equivalent to their number 0..n-1
    double64 display_value(0.);
    for ( size_t i=0U; i<element.IntegrationPoints(); i++ ) {
         Point<dim> xyz = element.IntegrationPoint(i);
         ipoint_data[i] = make_pair( xyz, display_value );
         display_value += 1.;
      }
   
     // 1. writing the file header
     // --------------------------
     string  outfile(file);
     // adding the element number
     outfile += to_string( element.Idx() );
     // adding the element type
     outfile += "_";
     outfile += parseFiniteElementType( element.FE()->ElementType() );
     outfile +=".vtk";
     ofstream  ofs(outfile.c_str());
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"integration points colored by number: '";
     ofs << "point_number" <<"'"<< endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3. writing point coordinates in order
     // --------------------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << ipoint_data.size() <<" float"<< endl;
     for ( typename vector<pair<Point<dim>,double64> >::const_iterator
           it=ipoint_data.begin(); it!=ipoint_data.end(); it++ ) {
          assert( dim != 1 );
          if ( dim == 2 ) ofs << (*it).first[0] <<" "<< (*it).first[1] <<" "<< 0. << endl;
          else ofs << (*it).first[0] <<" "<< (*it).first[1] <<" "<< (*it).first[2] << endl;
       }
     ofs << endl;  
       
     // 4. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     // (each record has a cell number and the corresponding point specifier)
     ofs <<"CELLS "<< ipoint_data.size() <<" "<< (ipoint_data.size() * 2) << endl;
     for ( size_t counter(0); counter<ipoint_data.size(); counter++ ) {
          ofs << 1 <<" "<< counter << endl;
       }
     ofs << endl;
     
     // 5. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< ipoint_data.size() << endl;
     // cell type VTK_POINT (=1)
     for ( size_t counter(0); counter<ipoint_data.size(); counter++ ) ofs << 1 << endl;
     ofs << endl;

     // 6. writing POINT_DATA point-type data values
     // ----------------------------------------------------
     ofs <<"POINT_DATA "<< ipoint_data.size() << endl;
     ofs <<"SCALARS "<< "point_number" <<" double"<< endl;
     ofs <<"LOOKUP_TABLE default" << endl;
     for ( typename vector<pair<Point<dim>,double64> >::const_iterator
           it=ipoint_data.begin(); it!=ipoint_data.end(); it++ ) {
          ofs << (*it).second <<" "<< endl;
       }
     ofs << endl;
     cout <<"\nVTK file '"<< outfile <<"' written successfully."<< endl;
      
 } // end outputIntegrationPointsToVTK


template void outputIntegrationPointsToVTK( const Element<2>&, const char* );
template void outputIntegrationPointsToVTK( const Element<3>&, const char* );




/**
     For the element by element output of integration point values.
     
     File name is appended with element number and element type.
*/
template<size_t dim>
void outputIntegrationPointDataToVTK( const Element<dim>& element, const csmp::Index& key,
                                      const char* variable_name, const char* file )
 {
    string var_name(variable_name);
    replaceWhiteSpaceBy( var_name, '_' );
  
     // 0. generating values, first integration point=0, second=1...
     // ------------------------------------------------------------
    vector<pair<Point<dim>,double64> >  ipoint_data(element.IntegrationPoints()); // locations in physical space
    // convention: integration points get values equivalent to their number 0..n-1
    double64 display_value(0.);
    for ( size_t i=0U; i<element.IntegrationPoints(); i++ ) {
         Point<dim> xyz = element.IntegrationPoint(i);
         ipoint_data[i] = make_pair( xyz, display_value );
         display_value += 1.;
      }
   
     // 1. writing the file header
     // --------------------------
     string  outfile(file);
     // adding the element number
     outfile += to_string( element.Idx() );
     // adding the element type
     outfile += "_";
     outfile += parseFiniteElementType( element.FE()->ElementType() );
     outfile +=".vtk";
     ofstream  ofs(outfile.c_str());
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"integration points colored by number: '";
     ofs << var_name <<"'"<< endl;
     ofs <<"ASCII"<< endl << endl;
     
     // 3. writing point coordinates in order
     // --------------------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << ipoint_data.size() <<" float"<< endl;
     for ( typename vector<pair<Point<dim>,double64> >::const_iterator
           it=ipoint_data.begin(); it!=ipoint_data.end(); it++ ) {
          assert( dim != 1 );
          if ( dim == 2 ) ofs << (*it).first[0] <<" "<< (*it).first[1] <<" "<< 0. << endl;
          else ofs << (*it).first[0] <<" "<< (*it).first[1] <<" "<< (*it).first[2] << endl;
       }
     ofs << endl;  
       
     // 4. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     // (each record has a cell number and the corresponding point specifier)
     ofs <<"CELLS "<< ipoint_data.size() <<" "<< (ipoint_data.size() * 2) << endl;
     for ( size_t counter(0); counter<ipoint_data.size(); counter++ ) {
          ofs << 1 <<" "<< counter << endl;
       }
     ofs << endl;
     
     // 5. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< ipoint_data.size() << endl;
     // cell type VTK_POINT (=1)
     for ( size_t counter(0); counter<ipoint_data.size(); counter++ ) ofs << 1 << endl;
     ofs << endl;

     // 6. writing POINT_DATA point-type data values
     // ----------------------------------------------------
     ofs <<"POINT_DATA "<< element.Nodes() << endl;
     ofs.setf( ios::scientific );

     if ( key.type == SCALAR )
       {
           ofs <<"SCALARS "<< var_name <<" double"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA 1 x nodes
           for ( size_t i=0; i<element.IntegrationPoints(); i++ ) ofs << element.Read(i,key) <<" ";
           ofs << endl;
       }
     else if ( key.type == VECTOR )
       {
          ofs <<"VECTORS "<< var_name <<" double"<< endl;
          vector<VectorVariable<dim> > data;
          element.IntegrationPointPropertyVector( key, data );
          // matrix DATA is vec-dim x nodes
          for ( size_t i=0; i<element.Nodes(); i++ ) {
               for ( size_t j=0; j<dim; j++ ) ofs << data[i](j) <<" ";
               // variables have always 3 components since view screen is 3D
               if ( dim == 2 ) ofs << 0.;
               ofs << endl;
            }
       }
     else if ( key.type == TENSOR )
       {
          ofs <<"TENSORS "<< var_name <<" double"<< endl;
          vector<TensorVariable<dim> > data;
          element.IntegrationPointPropertyVector( key, data );
          // matrix DATA is vec-dim x nodes
          for ( size_t i=0; i<element.Nodes(); i++ ) {
               for ( size_t j=0; j<dim; j++ ) {
                    for ( size_t k=0; k<dim; ++k ) {
                          ofs << data[i](j,k) <<" ";
                          if ( dim == 2 ) ofs << 0. <<" ";
                      }
                    if ( dim == 2 ) ofs << 0. <<" ";
                 }
               // variables have always 3 components since view screen is 3D
               if ( dim == 2 ) ofs << 0.;
               ofs << endl;
            }
       }
     else ErrorHandler::Instance().notice( CSMP_ERROR, "outputIntegrationPointDataToVTK:", parseType(key.type),
                                          "integration-point property placement cannot be output; nothing was done." );
     ofs << endl;
     cout <<"\nVTK file '"<< outfile <<"' written successfully."<< endl;
      
 } // end outputIntegrationPointDataToVTK


template void outputIntegrationPointDataToVTK( const Element<2>&, const csmp::Index&, const char*, const char* );
template void outputIntegrationPointDataToVTK( const Element<3>&, const csmp::Index&, const char*, const char* );


} // end namespace csmp





