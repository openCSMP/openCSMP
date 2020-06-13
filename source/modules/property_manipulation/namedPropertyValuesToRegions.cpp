//
//  namedPropertyValuesToRegions.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 9/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "namedPropertyValuesToRegions.h"
#include "Model.h"
#include "Region.h"

using namespace std;

namespace csmp {

/**
    Reads file with the mapping from integer codes to region names and creates corresponding model regions.
    Elements with property identifiers that are not contained in the regions file are ignore, but the
    unrecognised codes are reported to screen and file.
    The file has the following format:

    @code
    headline (printed to screen)
    identifier value  tab  region name (can contain white space)
    @endcode

    @attention if the rocktype file info file is omitted, the function will search for
    a file called:  "model_name-property-identifiers.txt"

    @note only the first 2 tokens in each line of the file are read, allowing the user to add comments afterwards

    @author SKM
    @date 22/8/2018
*/
template<size_t dim>
void  namedPropertyValuesToRegions( Model<dim>& model, const string& prop_name, const string& region_identifier_file )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. reading the rocktype identifier file (facies codes and corresponding facies names)
  // -------------------------------------------------------------------------------------
  string file_name = (region_identifier_file == "\0")
    ? string( model.Name() ) + "-region_identifiers.txt" : region_identifier_file;

  std::ifstream ifs( file_name.c_str() );

  if ( !ifs.is_open() )
    csmp_error.notice( FATAL_ERROR,
                       "namedPropertyValuesToRegions:", file_name,
                       "ASCII property identifier file could not be opened." );

  // 1.1 reading and discarding file header
  char         text_line[256];
  char*        token( 0 );
  const char* const delims = ",\t,:,\n,\r";

  // 1.2 reading the names of the properties
  ifs.getline( text_line, INFO_STRING );
  cout << "\nnamedPropertyValuesToRegions: reading file '" << file_name << "' with header: ";
  cout << "\n\t" << text_line << endl;

  // 1.3 reading the rocktype identifiers from file
  map<int32,string>  prop_value_region_name_mapping;

  while ( !ifs.eof() )
  {
    ifs.getline( text_line, INFO_STRING );
    if ( strlen( text_line ) == 0 ) break;
    // only the first 2 tokens are used alllowing the user to add comments afterwards
    // rocktype
    token = strtok( text_line, delims );
    int32 rocktype = (token != NULL) ? atoi( token ) : UNSPECIFIED;
    // facies name / association
    token = strtok( NULL, delims );
    string rocktype_name = (token != NULL) ? to_string( token ) : "UNSPECIFIED";

    prop_value_region_name_mapping.insert( make_pair( rocktype, rocktype_name ) );
  }
  ifs.close();


  // 2. creating unique model regions for each rocktype
  // -------------------------------------------------------------------------------------
  const csmp::Index  prop_key( model.Database().StorageKey( prop_name.c_str() ) );
  const bool         unique( true );

  cout << "\nnamedPropertyValuesToRegions: generating model regions from 'property' integer codes using data from file: " << file_name;
  for ( auto it = prop_value_region_name_mapping.begin(); it != prop_value_region_name_mapping.end(); ++it ) {
    cout << "\n\t" << it->second.c_str();
    model.FormRegionFrom( it->second.c_str(), prop_name.c_str(), it->first, it->first, unique );
  }
  cout << "\n\n";

  set<int32>   region_identifiers_without_name;
  Region<dim>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it )
  {
    const double64 value = (*it)->Read( prop_key );
    // checking that the property value can indeed be converted into an integer in a meaningful range
    // using the modulus operator % to determine whether the number has a decimal fraction
    if ( fmod( value, 1. ) != 0 ) {
         cerr <<"\n\tproperty value: "<< value << endl;
         csmp::Exception( ERROR, "namedPropertyValuesToRegions", "region identifier contains decimal places and can therefore not be converted to integer.");
      }
    const int32 region_identifier = static_cast<int32>(value);
    
    // if the rocktype can be identified, we store the element id for the later creation of a region
    if ( prop_value_region_name_mapping.find( region_identifier ) == prop_value_region_name_mapping.end() )
      region_identifiers_without_name.insert( region_identifier );
  }
  if ( !region_identifiers_without_name.empty() ) {
    for ( auto p : region_identifiers_without_name ) cerr << p << " ";
    csmp_error.notice( WARNING, "namedPropertyValuesToRegions:",
                       "there were elements with unrecognized property identifiers; they were ignored." );
  }

} // end namedPropertyValuesToRegions

template void namedPropertyValuesToRegions( Model<1U>&, const string&, const string& );
template void namedPropertyValuesToRegions( Model<2U>&, const string&, const string& );
template void namedPropertyValuesToRegions( Model<3U>&, const string&, const string& );




/// to remove NO_DATA values which were converted to NAN.
template<size_t dim>
void replaceElement_NAN_ValuesWith( Model<dim>& model, const char* element_var, double64 replacement_val )
 {
    const csmp::Index key(model.Database().StorageKey(element_var));
    if ( key.type!= SCALAR || key.place != ELEMENT )
      csmp::Exception( ERROR, "replaceElement_NAN_ValuesWith", "region identifier property must be a scalar placed on the element");
    
    Region<dim>& model_domain(model.Region("Model"));
    
    for ( typename vector<Element<dim>*>::iterator 
          it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        if ( isnan( (*it)->Read(key) ) )
          (*it)->Store( key, makeScalar(ANY,replacement_val) );
      }
    
 } // end replaceElement_NAN_ValuesWith

template void replaceElement_NAN_ValuesWith( Model<1U>&, const char*, double64 );
template void replaceElement_NAN_ValuesWith( Model<2U>&, const char*, double64 );
template void replaceElement_NAN_ValuesWith( Model<3U>&, const char*, double64 );

} // end csmp
