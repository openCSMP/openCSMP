//
//  SKUA_Interface.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/20/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#include "SKUA_Interface.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Point.h"
#include "CSMP_highLevelUtilities.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {


/**
    Writes file with the element barycentre coordinates, preceded by region-name(string), region-id(number), element-id(number).
    since SKUA doesn't import string type, a corresponding region-map and element-map will be generated
    to map their corresponding region-name and element-key from the region-id and the element-id.

    The file format is

    header line
    axes titles: region-name(string), region-id(number 0..n-1), element-id(number 0..n-1), x, y, z
    comma-seperated data...
    EOF
 
    @note to function correctly, the method needs the region id and the element invariable "element number"
    that must not change from run to run like CSMP's Idx() values do.
 
    @note ideally the region names should be echo'ed back from SKUA, but because this is a string and a unique integer
    code is written out and read in backwards as an alternative.
 
    @author SKM
*/
void SKUA_Interface::OutputElementNumbersAndBaryCentresRegionByRegion( const Model<3U>& model )
{
  ofstream ofs( string( model.Name() ) + "-element_barycentres.txt" );
  ofs << model.Name() << "-element_barycentres.txt  textfile with the element barcyentres of all elements in the mesh\n";
  ofs << "region-name,region-id,element number,x,y,z\n";
  printRangeOfVariable( model, "element number", true ); // == model.Region("Model").Cells() );
  size_t region_id = 0;
  
  const csmp::Index key_enr = model.Database().StorageKey("element number");
  
  for ( auto rit = model.UniqueRegionsBegin(); rit != model.UniqueRegionsEnd(); ++rit ) {
      for ( auto it = (*rit).second.CellsBegin(); it != (*rit).second.CellsEnd(); ++it )
        {
           // region name              id (0..n-1)       element number
           ofs << (*rit).first <<","<< region_id <<","<< static_cast<long>((*it)->Read( key_enr)) <<",";
           // point coordinates ordered and flipped to reflect SKUA's UTM lefthandrule coordinate system
           const Point<3U> bctr = (*it)->BaryCenter();
           //                   SKUA x             SKUA -y          SKUA -z
           ofs << scientific << bctr[0] << "," << -bctr[2] <<","<< -bctr[1] << endl;
        }
      region_id++;
   }
  
  ofs.close();
  cout << "\noutputElementBaryCentresRegionByRegion: " << model.Name() << "-element_barycentres.txt written successfully.\n";

} // end outputElementNumbersAndBaryCentresRegionByRegion






/**
     Using the supplied region names and "element number" values,
     maps the data from the input file onto the current model.
     The input file is expected to be named after the current model with the extension "-element_barycentres.txt".
     Its expected format is:
 
     header line
     axes titles: region-name(string), region-id(number 0..n-1), element-id(number 0..n-1), - comma-separated list of property names
     comma-seperated data...
     EOF

     @attention properties are mapped region by region so that values do not spill across region boundaries; not found non-NON-DATA values are reported
 
     @attention property values are mapped only if they are within range; out-of range values are reported, but ignored.
 
     @attention if a region does not exist, an error is reported
 
     @attention property values will only be written to elements whose property flag is not DIRICH or fixed
 
     @todo TODO: expand to handle tensor properties as well
 
     @author SKM
*/
template<uint32_t dim>
bool SKUA_Interface::ImportElementPropertyValuesFromSKUA( Model<dim>& model, const std::string& data_file )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    string  bc_file_name( model.Name() );
    bc_file_name += "-element_barycentre_properties.txt";
    // overwrite the default file name if a name argument was supplied to method
    if ( !data_file.empty() ) bc_file_name = data_file;

    ifstream  bc_ifs( bc_file_name );
    if ( !bc_ifs.is_open() ) {
      csmp_error.Note( ERROR, "importElementPropertyValuesFromSKUA:", bc_file_name, "could not be opened." );
      return false;
    }

    // 0. reading the contents of the file into maps of deques with the range-checked element data
    // -------------------------------------------------------------------------------------------
    const long LMAX( 1024 );
    char text_line[LMAX];
    string delimiter = ",";
    bc_ifs.getline( text_line, LMAX ); // reading the header line
    cout <<"\nimportElementPropertyValuesFromSKUA: reading text file: '"<< bc_file_name <<"'\n";
    cout <<"\n\tfile header: "<< text_line <<"\n";
    bc_ifs.getline( text_line, LMAX );
    cout <<"\tcontained variables: "<< text_line <<"\n";
   
    // 0.1 parsing the property list, checking whether the properties are available and what there ranges are
    // ------------------------------------------------------------------------------------------------------
    vector<string> properties = splitString( text_line, ',' );
    // expectations - DO NOT READ COORDINATE VALUES
    const size_t coords_to_read(0U);
    const size_t region_name_idx(0U), elmt_idx(2U), first_prop(3U+coords_to_read), items_per_line(properties.size());
    assert( items_per_line > first_prop );
    // checking that the properties do exist
    bool all_properties_known(true);
    for ( size_t i=first_prop; i<items_per_line; i++ )
       if ( !model.Database().IsDefined( properties[i].c_str() ) ) {
            csmp_error.Note( WARNING, "importElementPropertyValuesFromSKUA:", properties[i], "input property unknown to current model." );
            all_properties_known = false;
         }
    if ( !all_properties_known ) {
         csmp_error.Note( ERROR, "importElementPropertyValuesFromSKUA:", "unknown input properties; terminating mapping." );
         return false;
      }
    // getting bounds for the property values to be read
    vector<pair<double,double> > valid_ranges(items_per_line-first_prop,make_pair(-1.0e30,1.0e30));
    double min_val, max_val;
    for ( size_t i=first_prop; i<items_per_line; i++ ) {
         model.Database().RangeOf( properties[i].c_str(), min_val, max_val );
         valid_ranges[i-first_prop] = make_pair( min_val, max_val ); // ATTENTION - different range
      }
    // containers
    //  elmt-id, values of properties in order given in property string and boolean to tell whether value is to be mapped
    pair<size_t, vector<pair<bool, double> > > elmt_prop_values; //(UINT_MAX, items_per_line - first_prop);
    elmt_prop_values.second.resize( items_per_line - first_prop );

    // region by region values
    // region name, set of property values
    map<string,deque<pair<size_t,vector<pair<bool,double> > > > > region_data;
    // keeping statistics of range checks that failed and which properties were affected
    // property name, number of failures
    map<string,size_t> range_check_failures;
   
    // start to read properties
    bc_ifs.getline( text_line, LMAX );
    while ( !bc_ifs.eof() )
      {
        // tokenizing the text line
        vector<string> data_tokens = splitString( text_line, ',' );
        string region_name = data_tokens[region_name_idx];
        assert( !region_name.empty() );
        const long elmt_num = atol( data_tokens[elmt_idx].c_str() );
        assert( elmt_num >= 0U );
        elmt_prop_values.first = elmt_num;

        // reading the property values and performing range checks on the way
        for ( size_t i=first_prop; i<items_per_line; ++i ) {
              const double prop_value = stod( data_tokens[i] );
              // checking the range, including no-data values
              assert( !isnan(prop_value) );
              if ( prop_value < valid_ranges[i-first_prop].first || prop_value > valid_ranges[i-first_prop].second ||
                   static_cast<long>(prop_value) == -9999 || static_cast<long>(prop_value) == -99999  )
                {
                   if ( static_cast<long>(prop_value) != -9999 && static_cast<long>(prop_value) != -99999 ) {
                        cerr <<"\nOut-of-range value of '"<< properties[i] <<"' = "<< std::scientific << prop_value;
                        cerr <<" (region "<< region_name <<", element "<< elmt_num <<"), will be ignored.\n";
                        pair<map<string,size_t>::iterator,bool> it=range_check_failures.insert( make_pair(properties[i],1) );
                        if ( it.second == false ) (*it.first).second++; // incrementing the failure count
                     }                                             // false means that value will not be mapped
                   elmt_prop_values.second[i-first_prop] = make_pair( false, prop_value );
                }
              else elmt_prop_values.second[i-first_prop] = make_pair( true, prop_value );
          }
        
        // inserting the record into a new or existing map entry
        // region name, elements inside with their associated data
        // map<string,deque<pair<size_t,vector<pair<bool,double> > > > >
        map<string,deque<pair<size_t,vector<pair<bool,double> > > > >::iterator it = region_data.find(region_name);
        // if this is first element data set in the region data map
        if ( it == region_data.end() ) {
             deque<pair<size_t,vector<pair<bool,double> > > > new_data_set;
             new_data_set.push_back( elmt_prop_values );
             region_data.insert( make_pair(region_name,new_data_set) );
          }
        else (*it).second.push_back( elmt_prop_values );

        // reading next line
        bc_ifs.getline( text_line, LMAX );
      }
    bc_ifs.close();
   
   
    // 1. assigning the properties to the model, region-by-region, ignoring:
    //    - values that fall outside of the regions
    //    - NO_DATA values
    //    - out of range values
    // --------------------------------------------------------------------
    // 1.1 getting some property keys
    vector<csmp::Index> prop_keys;
    prop_keys.reserve( valid_ranges.size() );
    for ( size_t i=first_prop; i<items_per_line; i++ )
      prop_keys.push_back( model.Database().StorageKey(properties[i].c_str()) );
   
    // 1.2 establishing a mapping between current elements and id numbers
    const csmp::Index e_key = model.Database().StorageKey("element number");
   
    for ( map<string,deque<pair<size_t,vector<pair<bool,double> > > > >::iterator
          it=region_data.begin(); it!=region_data.end(); ++it )
      {
          // ignoring regions that do not exist in model
          if ( !model.ContainsRegion( (*it).first ) ) {
               csmp_error.Note( WARNING, "importElementPropertyValuesFromSKUA:", (*it).first,
                                 "region enlisted in file does not exist in current model; related data are ignored." );
               continue;
            }
          else if ( (*it).second.empty() )
            csmp_error.Note( WARNING, "importElementPropertyValuesFromSKUA:", (*it).first,
                              "region is enlisted in file, but without any assigneable values; nothing was done." );
          else {
               // finding the elements corresponding to element numbers in current model region
               Region<dim>& subdomain = model.Region((*it).first);
               map<size_t,Element<dim>*>  elmt_correspondance_map;
                for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
                  elmt_correspondance_map.insert( make_pair( (*eit)->Read(e_key), (*eit) ) );
            
               // assign properties element by element if these are valid
               for ( deque<pair<size_t,vector<pair<bool,double> > > >::iterator
                     et=(*it).second.begin(); et!=(*it).second.end(); ++et ) {
                     // searching for 'element number' in map
                     typename map<size_t,Element<dim>*>::iterator elmt_it = elmt_correspondance_map.find( (*et).first );
                     // if the element is within this region
                     if ( elmt_it != elmt_correspondance_map.end() ) {
                          Element<dim>* eptr = (*elmt_it).second;
                          // the property values are assigned to it
                          for ( size_t i{0U}; i<(*et).second.size(); ++i ) {
                               VARIABLE_FLAG flag = eptr->Status(prop_keys[i]);
                               if ( (*et).second[i].first == true && flag != DIRICH )
                                 eptr->Store( prop_keys[i], makeScalar( flag, (*et).second[i].second ) );
                            }
                       }
                 }
            }
      
      } // end loop over the regions
   
   // reporting
   if ( !range_check_failures.empty() ) {
        csmp_error.Note( WARNING, "importElementPropertyValuesFromSKUA:",
                              "range checks failed for several variable values." );
        for ( auto it=range_check_failures.begin(); it!=range_check_failures.end(); ++it )
          cerr <<"\n\t"<< (*it).second <<" range check failures occured for variable '"<< (*it).first <<"'";
        cerr << endl;
        return false;
     }
   
   cout <<"\nimportElementPropertyValuesFromSKUA: completed successfully.\n";
   
   return true;
   
} // end ImportElementPropertyValuesFromSKUA

template bool SKUA_Interface::ImportElementPropertyValuesFromSKUA( Model<2U>&, const std::string& );
template bool SKUA_Interface::ImportElementPropertyValuesFromSKUA( Model<3U>&, const std::string& );









/**  VariableToPointCloud

   Outputs any element variable from CSMP to SKUA point cloud (ASCII table) format.
   Variables discretized on volume, surface or line elements will be output 
   as point data placed on the barycenter of respective elements.  
   
   @attention methods outputs coordinates in SKUA format where y=-z  and  z=-y
   
*/   
void SKUA_Interface::VariableToPointCloud( const Model<3U>& model,
                                           const set<string>& regions_of_interest,
                                           const char* filename, const char* var ) const
 {
    const csmp::Index  var_key(model.Database().StorageKey(var));
    if ( var_key.place != ELEMENT )
      throw csmp::Exception( ERROR, "SKUA_Interface::VariableToPointCloud", "This method works only for 'Element' variables" );
   
    string var_name(var);
    replaceWhiteSpaceBy( var_name, '-' );
   
    // output file
    ofstream  ofs(filename);
    ofs <<"CSMP - SKUA_Interface::VariableToPointCloud:  output file '"<< filename <<"'\n";
   
    // setting to single precision (as SKUA only takes 7 decimal figures)
    ofs.setf(ios::scientific);
    const long prec = ofs.precision(7);
   
    // column titles
    ofs <<"\nmodel-region\tx\ty\tz\t";
    if ( var_key.type == SCALAR ) ofs << var_name;
    else if ( var_key.type == VECTOR ) ofs << var_name <<"[0]\t"<< var_name <<"[1]\t"<< var_name <<"[2]\n";
    else if ( var_key.type == TENSOR ) {
         int counter(0);
         for ( size_t i{0U}; i<3U; i++ )
           for ( size_t j{0U}; j<3U; j++ ) ofs << var_name <<"["<< counter++ <<"]\t";
      }
    else { // ARRAY variable
         for ( size_t i{0U}; i<var_key.index; i++ ) ofs << var_name <<"["<< i <<"]\t";
      }
    ofs <<"\n";
   
    // for all model regions
    for (  set<string>::const_iterator rt=regions_of_interest.begin(); rt!=regions_of_interest.end(); rt++ )
      {
         assert( model.ContainsRegion((*rt).c_str()) );
         const Region<3U>&  gref=model.Region((*rt).c_str());
         for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
           {
              ofs << (*rt) <<"\t";
              Point<3U> xyz((*it)->BaryCenter());
              // SKUA coordinates (x=x, y=-z, z=-y)
              ofs << xyz[0] <<"\t"<< -xyz[2] <<"\t"<< -xyz[1] <<"\t";
              if ( var_key.type == SCALAR ) ofs << (*it)->Read( var_key ) <<"\n";
              else if ( var_key.type == VECTOR ) {
                   VectorVariable<3U> vc;
                   (*it)->Read( var_key, vc );
                   ofs << vc[0] <<"\t"<< vc[1] <<"\t"<< vc[2] <<"\n";
                }
              else if ( var_key.type == TENSOR ) {
                   TensorVariable<3U> ts;
                   (*it)->Read( var_key, ts );
                   for ( size_t i{0U}; i<3U; i++ )
                     for ( size_t j{0U}; j<3U; j++ ) ofs << ts(i,j) <<"\t";
                   ofs <<"\n";
                }
              else { // ARRAY variable
                   ArrayVariable  ary;
                   (*it)->Read( var_key, ary );
                   for ( size_t i{0U}; i<ary.Size(); i++ ) ofs << ary[i] <<"\t";
                   ofs <<"\n";
                }
           }
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
   
    cout <<"\nSKUA_Interface::VariableToPointCloud: variable '";
    cout << var <<"' written successfully to file "<< filename << endl;
   
 } // end VariableToPointCloud







/// like previous method, but for multiple output variables discretized on the elements
void SKUA_Interface::VariablesToPointCloud( const Model<3U>& model,
                                            const set<string>& regions_of_interest,
                                            const char* filename, const set<string>& element_vars ) const
 {
    // processing variable names
    map<string,csmp::Index>  var_keys;
    for ( set<string>::const_iterator i=element_vars.begin(); i!=element_vars.end(); i++ )
      {
         const csmp::Index  var_key(model.Database().StorageKey((*i).c_str()));
         if ( var_key.place != ELEMENT )
            throw csmp::Exception( ERROR, "SKUA_Interface::VariablesToPointCloud", "This method works only for 'Element' variables" );
         if ( var_key.type != SCALAR )
            throw csmp::Exception( ERROR, "SKUA_Interface::VariablesToPointCloud", "This method works only for SCALAR variables" );
         string var_name((*i).c_str());
         replaceWhiteSpaceBy( var_name, '_' );
         var_keys.insert( make_pair(var_name,var_key) );
      }
   
    // output file
    ofstream  ofs(filename);
    ofs <<"CSMP - SKUA_Interface::VariablesToPointCloud:  output file '"<< filename <<"'\n";
   
    // setting to single precision (as SKUA only takes 7 decimal figures)
    ofs.setf(ios::scientific);
    const long prec = ofs.precision(7);
   
    // column titles
    ofs <<"\nmodel-region\tx\ty\tz\t";
    for ( map<string,csmp::Index>::const_iterator vit=var_keys.begin(); vit!=var_keys.end(); vit++ )
      ofs << (*vit).first <<"\t";
    ofs <<"\n";
   
    // for all model regions
    for ( set<string>::const_iterator rt=regions_of_interest.begin(); rt!=regions_of_interest.end(); rt++ )
      {
         assert( model.ContainsRegion((*rt).c_str()) );
         const Region<3U>&  gref=model.Region((*rt).c_str());
         for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
           {
              ofs << (*rt) <<"\t";
              Point<3U> xyz((*it)->BaryCenter());
              // SKUA coordinates (x=x, y=-z, z=-y)
              ofs << xyz[0] <<"\t"<< -xyz[2] <<"\t"<< -xyz[1] <<"\t";
              // writing the values of all the enlisted scalar element variables
              for ( map<string,csmp::Index>::const_iterator vit=var_keys.begin(); vit!=var_keys.end(); vit++ )
                ofs << (*it)->Read( (*vit).second ) <<"\t";
              ofs <<"\n";
           }
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
   
    cout <<"\nSKUA_Interface::VariablesToPointCloud: variables: ";
    for ( set<string>::const_iterator i=element_vars.begin(); i!=element_vars.end(); i++ ) cout <<"'"<< (*i) <<"'  ";
    cout <<"written successfully to file "<< filename << endl;
   
 } // end VariablesToPointCloud








/// creates point clouds around surface of interest using the 'thickness' attribute
// TODO: use UnitNormal() of surface elements to construct normals
void SKUA_Interface::SurfaceArrayVariableToPointCloud( const Model<3U>& model,
                                                       const set<string>& regions_of_interest,
                                                       const char* filename, const char* var ) const
 {
    const csmp::Index  var_key(model.Database().StorageKey(var));
    if ( var_key.place != ELEMENT )
      throw csmp::Exception( ERROR, "SKUA_Interface::SurfaceArrayVariableToPointCloud", "This method works only for 'Element' variables" );
    assert( var_key.type == ARRAY );
    string var_name(var);
    replaceWhiteSpaceBy( var_name, '-' );

    const csmp::Index  fth_key(model.Database().StorageKey("thickness"));
    assert( fth_key.type == SCALAR );
    Point<3U>      nrml;
    ArrayVariable  ary(6);

    // output file
    ofstream  ofs(filename);
    ofs <<"CSMP - SKUA_Interface::SurfaceArrayVariableToPointCloud:  output file '"<< filename <<"'\n";
   
    // setting to single precision (as SKUA only takes 7 decimal figures)
    ofs.setf(ios::scientific);
    const long prec = ofs.precision(7);
   
    // column titles
    ofs <<"\nmodel-region\tx\ty\tz\t"<< var_name <<"\n";

    // for all model regions
    for (  set<string>::const_iterator rt=regions_of_interest.begin(); rt!=regions_of_interest.end(); rt++ )
      {
         assert( model.ContainsRegion((*rt).c_str()) );
         const Region<3U>&  gref=model.Region((*rt).c_str());
         for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); ++it )
           {
              // each array variable entry is output as a singe line 
              assert( (*it)->FE()->IsSurface() );
             
              // getting the normal along which the output points will be created
              if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_QUADRILATERAL )
                nrml = normalAtFacetCenter( (*it)->N(0)->Coordinate(), (*it)->N(1)->Coordinate(),
                                            (*it)->N(2)->Coordinate(), (*it)->N(3)->Coordinate() );
              else
                if ( (*it)->FE_Type() == ISOPARAMETRIC_LINEAR_TRIANGLE )
                  nrml = normalOfTriangle( (*it)->N(0)->Coordinate(), (*it)->N(1)->Coordinate(), (*it)->N(2)->Coordinate() );
              else
              throw Exception( ERROR, "SKUA_Interface::SurfaceArrayVariableToPointCloud",
                              "attempt to compute normal on volume element rather than fault surface element");

              // reading the ARRAY variable and the thickness attribute
              (*it)->Read( var_key, ary );
              // writing out the first value representing center of surface
              Point<3U> xyz = (*it)->BaryCenter();
              // SKUA coordinates (x=x, y=-z, z=-y)
              ofs << (*rt) <<"\t" << xyz[0] <<"\t"<< -xyz[2] <<"\t"<< -xyz[1] <<"\t";

              // doing all subsequent points, assuming that values are symmetrically distributed around surface
              // and using the unit normal
              const double dx = (*it)->Read( fth_key ) / static_cast<double>(ary.Size()*2);
              for ( size_t i=1U; i<ary.Size(); i++ )
                {
                   Point<3U> out_pt =  (dx * i) * nrml;
                   Point<3U> ins_pt = (-dx * i) * nrml;
                   Point<3U> out = xyz + out_pt;
                   Point<3U> inp = xyz + ins_pt;
                   ofs << (*rt) <<"\t" << out[0] <<"\t"<< out[1] <<"\t"<< out[2] <<"\t"<< ary[i] <<"\n";
                   ofs << (*rt) <<"\t" << inp[0] <<"\t"<< inp[1] <<"\t"<< inp[2] <<"\t"<< ary[i] <<"\n";
                }
           }
      }

    ofs.unsetf( ios::scientific );
    ofs.precision(prec);
    ofs.close();
   
    cout <<"\nSKUA_Interface::SurfaceArrayVariableToPointCloud: variable '";
    cout << var <<"' written successfully to file "<< filename << endl;

 } // end SurfaceArrayVariableToPointCloud




// ======================================================================================================================

// PRIVATE AUXILIARY METHODS

// ======================================================================================================================


/**
    Enlists the unique 'element number' values of elements for which NO_DATA values (-9999, -99999)
    were reported in the input file in the target region.
 
    If the target region cannot be found of if there are no NO-DATA values in that region,
    the method returns false. Here, no-data values are expected for all the properties enlisted
    in the file header, e.g., rocktype, porosity and permeability.
 
    @attention the input data file is expected to have the extension:
    '"-element_barycentres_properties.txt"
 
    @author SKM
    @date 22/2/2019
*/
bool SKUA_Interface::Detect_NO_DATA_ElementsInDatasetFromSKUA( const string& input_txt_file,
                                                               const string& target_region,
                                                               std::set<size_t>& no_data_elmt_numbers )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    string datafile = input_txt_file + "-element_barycentres_properties.txt";

    ifstream  bc_ifs( datafile );
    if ( !bc_ifs.is_open() ) {
         csmp_error.Note( ERROR, "detect_NO_DATA_ElementsInDatasetFromSKUA:", datafile, "could not be opened." );
         return false;
      }
    if ( !no_data_elmt_numbers.empty() ) no_data_elmt_numbers.clear();


    // 0. scanning the contents of the file for the right region names and attributes
    // -------------------------------------------------------------------------------------------
    const long LMAX( 1024 );
    char text_line[LMAX];
    string delimiter = ",";
    bc_ifs.getline( text_line, LMAX ); // reading the header line
    cout <<"\ndetect_NO_DATA_ElementsInDatasetFromSKUA: reading text file: '"<< datafile <<"'\n";
    cout <<"\n\tfile header: "<< text_line <<"\n";
    bc_ifs.getline( text_line, LMAX );
    cout <<"\tcontained variables: "<< text_line <<"\n";
   
    // 0.1 parsing the property list, checking whether the properties are available and what there ranges are
    // ------------------------------------------------------------------------------------------------------
    vector<string> properties = splitString( text_line, ',' );
    // expectations
    const size_t coords_to_read(3U);
    const size_t region_name_idx(0U), elmt_idx(2U), first_prop(3U+coords_to_read), items_per_line(properties.size());
    assert( items_per_line > first_prop );
   
    // start to read properties
    bc_ifs.getline( text_line, LMAX );
    while ( !bc_ifs.eof() )
      {
        // tokenizing the text line
        vector<string> data_tokens = splitString( text_line, ',' );
        string region_name = data_tokens[region_name_idx];
        assert( !region_name.empty() );
        
        // skipping data from regions which are not of interest
        if ( region_name != target_region ) {
          bc_ifs.getline( text_line, LMAX );
          continue;
        }
        
        // else reading element numbers and storing them if all the property values are set to NO_DATA
        const long elmt_num = atol( data_tokens[elmt_idx].c_str() );
        assert( elmt_num >= 0U );

        // reading the property values and performing range checks on the way
        int no_data_count = 0.;
        for ( size_t i=first_prop; i<items_per_line; ++i ) {
              const double prop_value = stod( data_tokens[i] );
              // checking the value range, including no-data values
              assert( !isnan(prop_value) );
              if ( static_cast<long>(prop_value) == -9999 || static_cast<long>(prop_value) == -99999  )
                no_data_count++;
          }
        if ( no_data_count == (items_per_line - first_prop) )
          // for any NO_DATA element record insert element number into the element number set
          no_data_elmt_numbers.insert( static_cast<uint32_t>(elmt_num) );

        // reading next line
        bc_ifs.getline( text_line, LMAX );
      }
    bc_ifs.close();
   
   // reporting
   if ( no_data_elmt_numbers.empty() ) {
        csmp_error.Note( INFO, "detect_NO_DATA_ElementsInDatasetFromSKUA:", target_region,
                          "no empty element property records found for target region in datafile." );
        return false;
     }

   cout <<"\ndetect_NO_DATA_ElementsInDatasetFromSKUA: found "<< no_data_elmt_numbers.size();
   cout <<" empty property data records for region '"<< target_region <<"' in datafile. Function completed successfully.\n";
   
   return true;
   
} // end detect_NO_DATA_ElementsInDatasetFromSKUA






/**
removes the elements for which NO_DATA values (-9999, -99999) in the target region of the model.
After removing the elements, the unique regions, boundaries, and split boundaries are updated accordingly.

@author JC revised by SKM
@date 26/2/2019
@date 22/5/2021

*/
void SKUA_Interface::Erase_NO_DATA_ElementsFromModel( Model<3U>& model, const string& target_region, std::set<size_t>& no_data_elmt_numbers )
{
  csmp::Region<3U>& region = model.Region( target_region );
  
  vector<size_t>              element_ids( no_data_elmt_numbers.begin(), no_data_elmt_numbers.end() );
  vector<csmp::Element<3U>*>  ptrs_to_removed_elements;
  size_t n_removed_elmts    = region.RemoveByNumber( element_ids, ptrs_to_removed_elements );
  
  // reporting
  std::cout << "\nremove_NO_DATA_ElementsInModel: " << region.Name() << " (removed elements: " << n_removed_elmts << ")";
  std::cout << "\t" << region.Cells() << " elements remaining in '" << region.Name() << "'.\n";

  // updating the model, dependent on wether the removed elements were located only in a single unique region or across regions
  // if the region is unique ony that region needs to be modified
  if ( model.IsUnique(target_region) ) {
       // finding the target elements
       vector<Element<3U>*> elmt_ptrs;
       model.Mesh().DeleteCellsAndRepairConnnectivity( ptrs_to_removed_elements.begin(), ptrs_to_removed_elements.end() );
       return;
    }

  // if the region was non-unique, i.e., overlapping other regions, all regions the overlapped regions need to be rebuild
  // updating regions
  model.Mesh().UpdateConnectivity();
  
} // end Remove_NO_DATA_ElementsInModel





/**
    Reads file with the mapping from integer codes to region names and creates corresponding model regions.
    Elements with rocktype identifiers that are not contained in the regions file are ignore, but the
    unrecognised codes are reported to screen and file.
    The file has the following format:

    @code
    headline (printed to screen)
    identifier value  tab  region name (can contain white space)
    @endcode

    @attention if the rocktype file info file is omitted, the function will search for
    a file called:  "model_name-rocktype-identifiers.txt"

    @note only the first 2 tokens in each line of the file are read, allowing the user to add comments afterwards

    @author SKM
    @date 22/8/2018
*/
void  SKUA_Interface::ConvertRockTypesIntoRegions( Model<3U>& model, const string& rocktype_info_file )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. reading the rocktype identifier file (facies codes and corresponding facies names)
  // -------------------------------------------------------------------------------------
  string file_name = (rocktype_info_file == "\0")
    ? string( model.Name() ) + "-rocktype-identifiers.txt" : rocktype_info_file;

  std::ifstream ifs( file_name.c_str() );

  if ( !ifs.is_open() )
    csmp_error.Note( FATAL_ERROR,
                       "convertRockTypesIntoRegions:", file_name,
                       "ASCII rocktype identifier file could not be opened." );

  // 1.1 reading and discarding file header
  char         text_line[256];
  char*        token( 0 );
  const char* const delims = " ,\t,:,\n,\r";

  // 1.2 reading the names of the properties
  ifs.getline( text_line, 500 );
  cout << "\nconvertRockTypesIntoRegions: reading file '" << file_name << "' with header: ";
  cout << "\n\t" << text_line << endl;

  // 1.3 reading the rocktype identifiers from file
  map<int32_t, string>  rocktype_identifiers;

  while ( !ifs.eof() )
    {
      ifs.getline( text_line, 256 );
      if ( strlen( text_line ) == 0 ) break;
      // only the first 2 tokens are used alllowing the user to add comments afterwards
      // rocktype
      token = strtok( text_line, delims );
      int32_t rocktype = (token != NULL) ? atoi( token ) : UNSPECIFIED;
      // facies name / association
      token = strtok( NULL, delims );
      string rocktype_name = (token != NULL) ? token : "UNSPECIFIED";

      rocktype_identifiers.insert( make_pair( rocktype, rocktype_name ) );
    }
  ifs.close();


  // 2. creating unique model regions for each rocktype
  // -------------------------------------------------------------------------------------
  const csmp::Index  rrt_key( model.Database().StorageKey( "rocktype" ) );
  const bool         unique( true );

  cout << "\nconvertRockTypesIntoRegions: generating model regions from 'rocktype' integer codes using data from file: " << file_name;
  for ( auto it = rocktype_identifiers.begin(); it != rocktype_identifiers.end(); ++it ) {
    cout << "\n\t" << it->second.c_str();
    model.FormRegionFrom( it->second.c_str(), "rocktype", it->first, it->first, unique );
  }
  cout << "\n\n";

  set<int32_t>   unknown_identifiers;
  Region<3U>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.CellsBegin(); it != model_domain.CellsEnd(); ++it )
    {
      const int32_t rocktype = static_cast<int32_t>((*it)->Read( rrt_key ));
      // if the rocktype can be identified, we store the element id for the later creation of a region
      if ( rocktype_identifiers.find( rocktype ) == rocktype_identifiers.end() )
        unknown_identifiers.insert( rocktype );
    }
  if ( !unknown_identifiers.empty() ) {
    for ( auto p : unknown_identifiers ) cerr << p << " ";
    csmp_error.Note( WARNING, "convertRockTypesIntoRegions:",
                       "there were elements with unrecognized rocktype identifiers; they were ignored." );
  }

} // end convertRockTypesIntoRegions





} // end csmp
