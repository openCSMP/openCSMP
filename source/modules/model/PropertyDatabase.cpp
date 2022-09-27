#include "PropertyDatabase.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "Standard_IO_Handler.h"
#include "BE_TokenIterator.h"
#include "IndexTracker.h"
#include "binaryReadWrite.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include <sstream>

using namespace std;

namespace csmp {


#define PHYSVARS_BFILE "CSMP_variables.dat"
#define MEMBERS 10 // maximum entries in one line of property description


/** Initializes variables database from text variables file
 */
template<uint32_t dim>
PropertyDatabase<dim>::PropertyDatabase( const char* variablesFileName )
    : vectorFlags(dim),
      tensorFlags(dim),
      physvarsFile((variablesFileName==NULL) ? "EmptyVariablesFile" : variablesFileName),
      variableCount_(),
      propList_(),
      indexTracker_(),
      verbose_(true)
 {
    Initialize(variablesFileName);
 }


/**
      Reader of binary variables file which also permits restriction to a subset of variables.
      
      @attention if the variablesFileName is empty, the variables file name is initialised to 'EmptyVariablesFile.'
*/
template<uint32_t dim>
PropertyDatabase<dim>::PropertyDatabase( const char* variablesFileName, const set<string>& subset_variables )
    : vectorFlags(dim),
      tensorFlags(dim),
      physvarsFile((variablesFileName==NULL) ? "EmptyVariablesFile" : variablesFileName),
      variableCount_(),
      propList_(),
      indexTracker_(),
      verbose_(true)
 {
    if ( !BinaryIn( variablesFileName, subset_variables ) )
      throw csmp::Exception( ERROR, "PropertyDatabase", "Not able to load from binary file" );
 }



/// default constructor: all variable counters are set to zero
template<uint32_t dim>
PropertyDatabase<dim>::PropertyDatabase()
 : vectorFlags(dim),
   tensorFlags(dim),
   physvarsFile("Empty_Database"),
   variableCount_(),
   propList_(),
   indexTracker_(),
   verbose_(true)
 {
    InitializeCount();
 }




/// @attention We do not copy the IndexTracker here!
template<uint32_t dim>
PropertyDatabase<dim>::PropertyDatabase( const PropertyDatabase<dim>& p )
  : vectorFlags(dim),
    tensorFlags(dim),
    physvarsFile(p.physvarsFile), 
    variableCount_( p.variableCount_),
    propList_(p.propList_),
    indexTracker_(),
    verbose_(true)
  {
    UpdateParametersAndDatabase();
  }













///  Initialization with variable specifications read from '*-variables.txt' ascii file.
template<uint32_t dim>
void PropertyDatabase<dim>::Initialize( const char* variables_file )
  {
    if (variables_file) { // this was the pre-existing logic(JEM - July 11-2014)
        physvarsFile = variables_file; // not sure why it was so, but I left it this way

        InitializeCount();

        cout <<"\nPropertyDatabase::Initialize: Initializing database from textfile: ";
        cout << physvarsFile << endl;

        TextToBinaryFile( physvarsFile.c_str() );

        FlushToScreen();
    }
  }



template<uint32_t dim>
void PropertyDatabase<dim>::InitializeVariableTypeCount( std::map<VARIABLE_TYPE,uint32_t>& typeCount )
  {
    set<VARIABLE_TYPE> types;
    variableTypeSet(types);
    for( set<VARIABLE_TYPE>::const_iterator it( types.begin() ); it != types.end(); ++it )
      typeCount[*it] = 0;
  }


/// Initializes private container to all available placements
template<uint32_t dim>
void PropertyDatabase<dim>::InitializeCount()
  {
    set<PLACEMENT> places;
    variablePlacementSet(places);
    for( set<PLACEMENT>::const_iterator it( places.begin() ); it != places.end(); ++it )
      InitializeVariableTypeCount( variableCount_[*it] );
  }



template<uint32_t dim>
PropertyDatabase<dim>& PropertyDatabase<dim>::operator=( const PropertyDatabase<dim>& p )
 {
    if ( &p == this ) return *this;

    physvarsFile = p.physvarsFile;
    variableCount_ = p.variableCount_;
    propList_  = p.propList_;
    verbose_ = p.verbose_;

    UpdateParametersAndDatabase();

    return *this;
 }


template<uint32_t dim>
PropertyDatabase<dim>::~PropertyDatabase()
  {
    indexTracker_.DetachFromAll();
  }




// property iterators
template<uint32_t dim>
std::map<std::string,Parameter>::const_iterator  PropertyDatabase<dim>::Begin() const
 {
    return propList_.begin();
 } 
 
 
template<uint32_t dim>
std::map<std::string,Parameter>::const_iterator  PropertyDatabase<dim>::End() const
 {
    return propList_.end();
 } 


template<uint32_t dim>
std::map<PLACEMENT,std::map<VARIABLE_TYPE,uint32_t> >::const_iterator PropertyDatabase<dim>::VariableCountEnd() const
  {
    return variableCount_.end();
  }


template<uint32_t dim>
std::map<PLACEMENT,std::map<VARIABLE_TYPE,uint32_t> >::const_iterator PropertyDatabase<dim>::VariableCountBegin() const
  {
    return variableCount_.begin();
  }

/**
 
Tells the user whether a variable exists in the property database. 

@param s The name or csmp::Index of the variable which shall be searched for in the
database. 

@return The boolean variable 'true' or 'false'.
*/
template<uint32_t dim>
bool PropertyDatabase<dim>::IsDefined( const char* s ) const
 {
    if ( propList_.find(std::string(s)) != propList_.end() ) return true;
     return false;
 }


/**
 
Returns the placement (Node, IntegrationPoint or Element) of the target
physical variable. 

@param s The name of the physical variable whose placement shall be determined.

@return The variable placement.

@section messages Messages

If the variable is not defined, an error will be reported. 
*/
template<uint32_t dim>
PLACEMENT PropertyDatabase<dim>::Placement( const char* s ) const 
  {
     auto iter(propList_.find(std::string(s)));
     
     if ( iter != propList_.end() ) 
       return (*iter).second.key.place;
     else 
       throw csmp::Exception( ERROR, "PropertyDatabase<dim>::Placement", "Unable to identify property", s );

     // shouldn't get here
     return MODEL;        
  } // end WhereIs   



/**
 
Returns the csmp::Index.index of the target physical variable. 

@param s The name of the physical variable whose index shall be determined.

@return The variable index which defines the number of the corresponding property
vector inside of the MemoryManager.  If undefined, numeric_limits<uint32_t>::max() is returned.

@section messages Messages

If the variable is not defined, an error will be reported. 
*/
template<uint32_t dim>
uint32_t PropertyDatabase<dim>::Index( const char* s ) const
  {
     auto iter(propList_.find(std::string(s)));
     
     if ( iter != propList_.end() ) 
       return (*iter).second.key.index;
     else 
       std::cerr <<"\nPropertyDatabase::Index: Unable to identify index of: " << s << std::endl;

     return numeric_limits<uint32_t>::max(); // ULONG_MAX
     
  } // end Index 



/**
 
Reports whether the target physical variable is of scalar, vector or
tensor type. 

@param s The name of the physical variable whose type shall be established.

@return The type of the physical variable.


@section messages Messages

If the variable is not defined, an error will be reported. 
*/
template<uint32_t dim>
VARIABLE_TYPE  PropertyDatabase<dim>::Type( const char* s ) const 
  {
     auto iter(propList_.find(std::string(s)));
     
     if ( iter != propList_.end() ) 
       return (*iter).second.key.type;
     else 
       std::cout <<"\n\nPropertyDatabase::Type: Unable to identify type of variable: '" << s <<"'\n";

     return static_cast<VARIABLE_TYPE>(-1);  
  } // end Type 

/**
 
Returns the number of components of a variable.

Scalar: 1
Vector: 2
Tensor: 3

Array:          variable (as defined in the variables file)
FlaggedArray:   variable (as defined in the variables file)


@param s The name of the physical variable whose type shall be established.

@return The number of components of a variable.

@section messages Messages

If the variable is not defined, an error will be reported. 
*/
template<uint32_t dim>
uint32_t  PropertyDatabase<dim>::Components( const char* s ) const
  {
     auto iter(propList_.find(std::string(s)));
     
     if ( iter != propList_.end() ) 
       return (*iter).second.key.dataDepth;
     else 
       std::cout <<"\n\nPropertyDatabase::Type: Unable to identify component count of variable: '" << s <<"'\n";

     return -1;  
  } // end Type 




/**
 
Returns the Index type for the target physical variable. This index is 
an efficient access key for the physical variable andcan be used in
repeated operations on the property storage, such as Read() / Write()
operations of variable values. 

@param s The name of the variable whose Index shall be retrieved.

@return csmp::Index initialised with variable specifications

A Index class object specifying the type, placement and property-array 
index of the target variable. If the variable cannot be found in the 
database the returned Index will be initialized using its default 
constructor. 

@section application Application

StorageKey() is used inside Interrelation, Algorithm, Visitor and other
class objects to derive Index keys for computations. 

@section messages Messages

If the variable does not exist in the property database, an error will be
reported. 
*/
template<uint32_t dim>
csmp::Index  PropertyDatabase<dim>::StorageKey( const char* s ) const
 {
    auto iter(propList_.find(std::string(s)));
    if ( iter == propList_.end() ) {
         std::string message("Variable '");
         message += s;
         message +="' is undefined";
         throw csmp::Exception( ERROR, "PropertyDatabase<dim>::StorageKey:", message );
      }
    if ( iter == propList_.end() ) 
      return csmp::Index();

    return (*iter).second.key;
 }

template<uint32_t dim>
csmp::Parameter  PropertyDatabase<dim>::Parameter( const char* s ) const
{
  auto iter(propList_.find(std::string(s)));
  if ( iter == propList_.end() ) {
    std::string message("Variable '");
    message += s;
    message +="' is undefined";
    throw csmp::Exception( ERROR, "PropertyDatabase<dim>::Parameter:", message );
  }
  if ( iter == propList_.end() ) 
    return csmp::Parameter();

  return (*iter).second;
}



/// reports how the variable is normally used (as input by user vs. computed)
template<uint32_t dim>
const char*  PropertyDatabase<dim>::Usage( const char* s ) const
 {
    auto iter(propList_.find(std::string(s)));
    if ( iter == propList_.end() ) {
         std::string message("Variable '");
         message += s;
         message +="' is undefined";
         throw csmp::Exception( ERROR, "PropertyDatabase<dim>::Usage:", message );
      }
    if ( iter == propList_.end() ) return "undefined";

    return (*iter).second.usage.c_str();
 }





/**
 Binary out.

 There is a special implementation here since we don't want to store the IndexTracker association with the Parameter Index objects.
 We write to binary following the below outline:
 -# Number of parameters (uint32_t)
 -# Each parameter
  - Name (char)
  - csmp::Parameter

 @author  P. Lang
 @date  10/16/2012

 @tparam  dim Spatial dimension
 @param [in,out]  fp  If non-null, the binary file pointed.

 @return  true if it succeeds, false if it fails.
 */
template<uint32_t dim>
bool PropertyDatabase<dim>::BinaryOut( fstream& fp ) const
  {
    size_t parameterCount( propList_.size() );
    fp.write( (char*) &parameterCount, sizeof(size_t) );

    for( auto& prop : propList_ )
      {
        binaryFileWrite( fp, prop.first.c_str() );
        prop.second.Out(fp);
      }

      return true; /// @todo (1-C) Meaningless return statement
  }


/// Writes to binary file, appends '_variables.dat' to fileName if necessary
template<uint32_t dim>
bool PropertyDatabase<dim>::BinaryOut( const char* fileName ) const
  {
  string outputFileName(fileName);

  fstream fp(outputFileName.c_str(), ios::out | ios::binary );
  if (!fp.is_open())
   return false;

  if( !BinaryOut(fp) )
    return false;

  fp.close();
  return true;
  }


/**
 Binary input. See BinaryOut

 @author  P. Lang
 @date  10/16/2012

 @tparam  dim Type of the dim.
 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.
 */
template<uint32_t dim>
bool PropertyDatabase<dim>::BinaryIn( fstream& fp, const set<string>& subset_variables )
{
  size_t parameterCount(0);
  fp.read( (char*) &parameterCount, sizeof(size_t) );

  char buf[255];
  for( auto i{0U}; i < parameterCount; ++i )
  {
    binaryFileRead( fp, buf );
    string parameterName(buf);
    csmp::Parameter parameter;
    parameter.In( fp );

    if ( !subset_variables.empty() ) {
        if ( subset_variables.find( parameterName ) != subset_variables.end() )
          propList_[parameterName] = parameter;
      }
    else propList_[parameterName] = parameter;
  }

  return true; /// @todo (1-C) Meaningless return statement
}


/// Reads from binary file, appends '_variables.dat' if necessary, and can read only a subset of variables from the variables as an option
template<uint32_t dim>
bool PropertyDatabase<dim>::BinaryIn( const char* fileName, const set<string>& subset_variables )
  {
    InitializeCount();

    string inputFileName(fileName);

    fstream fp(inputFileName.c_str(), ios::in | ios::binary);
    if (!fp.is_open())
      return false;

    if( !BinaryIn(fp, subset_variables ) )
      return false;

    fp.close();
    UpdateParametersAndDatabase();
    FlushToScreen();
    return true;
  }





template<uint32_t dim>
const char*  PropertyDatabase<dim>::VariablesFile() const
 {
   return physvarsFile.c_str();
  }

/** Tells the user whether a variable exists in the property database.

@param idx The csmp::Index of the variable which shall be searched for in the database. 

@return The boolean variable 'true' or 'false'.

*/
template<uint32_t dim>
bool PropertyDatabase<dim>::IsDefined( const csmp::Index& idx ) const
 {
     for ( auto& prop : propList_ )
       if ( prop.second.key == idx ) 
            return true;

     return false;
 }



/// Prints all parameter records in the database to std output.  
template<uint32_t dim>
void PropertyDatabase<dim>::FlushToScreen() const
 {
     map<string,csmp::Parameter>::const_iterator  iter;

     cout <<"\nPropertyDatabase::FlushToScreen: Current properties in alphabetical order: "<< endl;
     for ( auto& prop : propList_ )
       cout << prop.second;
   
     cout <<"\n\n";
     cout.flush();
    
 } // FlushToScreen



/** Prints the parameter record of the target variable to stdout.

@param propname name  of variable to be output

@section implementation Implementation

The method uses the 'Out()' interface of the Parameter class 
to print the variable record(s). 
*/
template<uint32_t dim>
void PropertyDatabase<dim>::FlushToScreen( const char* propname ) const
 {
     auto iter=propList_.find(string(propname));

     if ( iter == propList_.end() ) {
          cout <<"\nPropertyDatabase::FlushToScreen: Property '"<< propname <<"' is undefined."<< endl;
          return;
       }
       
     (*iter).second.Out();
     cout.flush();
    
 } // FlushToScreen


/**
 
Counts the total number of variables which reside in the property
database and uses the numbers to update its internal counters. 

@section application Application

Private method of the property database. 

@section messages Messages

If the property database is empty, a warning message will be reported. 
*/
template<uint32_t dim>
void PropertyDatabase<dim>::CountVariables()
  {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     
     if ( (propList_.empty()) ) {
          csmp_error.Note( WARNING, "PropertyDatabase<dim>::CountVariables", 
                                      "the property list is empty");
          return;
       }

     InitializeCount();
     
     for ( auto& prop : propList_ )
         ++( variableCount_[ prop.second.key.place ] [ prop.second.key.type ] );
  } // end CountVariables


/**
 
Once all the physical variables which describe a CSMP model have been 
successfully read from the input file, AssignVariableIndices() 
calculates CSMP indices for these which will be used by the 
MemoryManager as it builds the property storage. 

@section implementation Implementation

The PropertyDatabase stores the variable records in alphabetical order
inside of a map. The method iterates through this map assigning array
indices to each variable. These are later used as csmp::Index.index 
variables. 

@section application Application

Private method which is called when the property database is initialized
from file. 

@section messages Messages

If the property database is empty, a warning message will be reported. 
*/
template<uint32_t dim>
void PropertyDatabase<dim>::AssignVariableIndices() 
  {
     InitializeCount();
     
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     
     if ( (propList_.empty()) ) {
          csmp_error.Note( WARNING, "PropertyDatabase<dim>::AssignVariableIndices", 
                                      "the property list is empty");
          return;
       }
          
     for ( auto& prop : propList_ )
       prop.second.key.index = (variableCount_[prop.second.key.place][prop.second.key.type])++;
       
 } // AssignVariableIndices


/** Prints the Parameter records of all current physical variables to stdout. 

@section implementation Implementation
Method calls the Out() interface of the Parameter class. 

@section application Application
To get full descriptions of the characteristics of all model variables. 
*/
template<uint32_t dim>
void PropertyDatabase<dim>::ListVariables() const 
  {
     for ( auto& prop : propList_ )
       prop.second.Out();

  } // ListVariables


/**
    Prints an alphabetically ordered list of all current physical variables to 'stdout'.
    
    @test last modified: SKM 25/10/14: fixed segmentation fault caused by printing empty variable.
*/
template<uint32_t dim>
void PropertyDatabase<dim>::ListVariableNames() const
  {
     cout <<"\nPropertyDatabase::ListVariableNames: Variables in current database: "<< endl;
     for ( auto& prop : propList_ )
       // if there is a variable that has not been initialized
       if ( !prop.first.empty() )
         cout << prop.first << endl;

  } // CountVariables


/**

Reads a textfile defining the characteristics of each physical variable 
which shall be used in a computation. This 'text' file may be created 
by hand or exported as "tab-delimited" from an Excel spreadsheet variable 
database file. The successfully read input is translated into the 
default PHYSVARS_BFILE file.  

@section application Application

Private method which is used by the constructor of the property database. 

@section messages Messages

The method will report an error and return without reading the file if it
is unable to open the variable database input text file. 

@attention Variable enum int equivalents hardcoded here
 */
template<uint32_t dim>
void PropertyDatabase<dim>::TextToBinaryFile( const char* property_database_textfile, bool echo_to_screen )
  {
     ifstream  ifs( property_database_textfile );
     
     if ( !ifs.is_open() )
       throw csmp::Exception( ERROR, "PropertyDatabase<dim>::TextToBinaryFile",
                              property_database_textfile, "could not be opened; nothing was done" );

     if ( !propList_.empty() ) {
          cout <<"\nPropertyDatabase::TextToBinaryFile: property list was not empty; zapping it."<< endl;
          propList_.erase( propList_.begin(), propList_.end() );
       }  
       
     istreambuf_iterator<char>  ifsBegin(ifs), ifsEnd;
     Delimiters  delimiters("\t<>\n\r"); // check - this may be different from token to token
     TokenIterator<istreambuf_iterator<char>,Delimiters>  propertyIter( ifsBegin, ifsEnd, delimiters ),
                                                          propertiesEnd;   
     // lets read the first two lines
     bool  place_specified(false), usage_specified(false), explanation_specified(false), reference_specified(false); 
     
     while ( propertyIter!=propertiesEnd ) {
          if ( place_specified ) {
    	         if      ( *propertyIter == "usage" )       usage_specified       = true;
    	         else if ( *propertyIter == "explanation" ) explanation_specified = true;
    	         else if ( *propertyIter == "reference" )   reference_specified   = true;
    	         else break;
	          }
          if ( *propertyIter == "place" ) place_specified = true;
          propertyIter++;
       }

     // read data line by line 
     csmp::Parameter   new_param;
     while ( propertyIter!=propertiesEnd ) {
          // name	
          new_param.name = *propertyIter++;

          // variable notation
          new_param.notation = *propertyIter++;

          // unit
          new_param.unit = *propertyIter++;

          // variable type (scalar=1,vector=2,tensor=3, array=4..max<uint32_t>, flagged array=-1...-max<uint32_t>) and corresponding components/depth
          /// Roman, 2013: in addition to the old style the new way of reading type's is added
          /// The type's can be represented by name and size specification
          /// For types SCALAR, VECTOR, TENSOR there is no need to specify the size
          /// For ARRAY and FLAGGEDARRAY size should be specified after the "space" following the name of the variable,
          /// e.g. ARRAY 1000, FLAGGEDARRAY 2, but SCALAR, VECTOR, TENSOR
          /// Note: name of the variable is not case sensitive, i.e. one can put FlaggedArray or FLAGGEDARRAY or flaggedarray
          int varTypeQualifier( atoi((*propertyIter).c_str() ) );
          if( varTypeQualifier!=0)
              EstablishVariableTypeDependentProperties( varTypeQualifier, new_param.key );
          else
              EstablishVariableTypeDependentProperties( *propertyIter, new_param.key );

          propertyIter++;
          // physically realistic minimum and maximum values
          double vmin( atof( (*propertyIter++).c_str() ) ), 
                   vmax( atof( (*propertyIter++).c_str() ) );
          //cout<<"property name: "<< new_param.name<<endl;
          assert( vmin <= vmax );
          assert( vmin >= -1e50 and vmax <= 1e50 );
          new_param.min = vmin;
          new_param.max = vmax;

          // variable placement  
          new_param.key.place = parsePlacement( (*propertyIter++).c_str() );
          EstablishPlacementDependentProperties(  new_param.key.place, new_param.key );

          // optional specifications (usage, explanation, and reference)
          if ( usage_specified )       new_param.usage = *propertyIter++;
          if ( explanation_specified ) new_param.explanation = *propertyIter++;
          if ( reference_specified )   new_param.reference = *propertyIter++;
	      
            // storing new parameter
            propList_.emplace( new_param.name, new_param );

	        // testing
  	      if ( echo_to_screen ) new_param.Out();
       }
     
     ifs.close();
     
     // initialize indices and create binary database
     UpdateParametersAndDatabase(); 
     
  } // end TextToBinaryFile


/**

Interactive method allowing to a add a property at runtime, provided 
that the user knows which storage-vector index he can assign (this can 
be established by first calling ScalarProperties() or an equivalent method
to establish how many properties of the desired type there already
are. The index of the new property should use the next higher integer such
that the property can be allocated by the MemoryManager at the end
of the current property storage. 

@warning Any runtime change in variables invalidates existing Index objects! They need to be updated.

@section arguments Input Arguments

The user will be prompted to enter the Parameter record specifications. 

A csmp::Index for the new property. 

@section implementation Implementation

The method will first prompt the user to define the characteristics of 
the new variable. Then CountVariables() is called to update the variable
number lists inside of the PropertyDatabase. Finally, the updated variable
list will be written to file. 

@section application Application

To define a new variable at runtime. 

@section messages Messages

The method will report a warning and print the existing variable record
if the new variable already exists in the property database. 

@return The Index reflecting the new state ov LocalVariables and IntegrationPointVariables, including added property.

@todo (3-D) These cout should become either Exception or errors using ErrorHandler
@todo (3-D) Refactor to use overload internally
*/
template<uint32_t dim>
csmp::Index  PropertyDatabase<dim>::AddProperty() 
 {    
    string     new_prop;
    csmp::Parameter  added_prop;

    cout<<"\nEnter property you would like to add: ";
    cin >> new_prop;
    
    auto iter(propList_.find( new_prop ));

    if ( iter != propList_.end() ) {
         cout <<"\nWARNING, PropertyDatabase<dim>::AddProperty property '"<< new_prop <<"' already exists" << endl;
         (*iter).second.Out();
         return StorageKey( new_prop.c_str() );
      }
    else  
      {
         added_prop.DefineFromStdin();
         propList_[ added_prop.name ] = added_prop;
         
         VARIABLE_TYPE vType = (*propList_.find(new_prop)).second.key.type;

         /// Roman,2013: Added explicit way of reading the size of variable provided by user input in function added_prop.DefineFromStdin()
         uint32_t size = (*propList_.find(new_prop)).second.key.dataDepth;
         //EstablishVariableTypeDependentProperties( static_cast<int>(vType), added_prop.key );
         EstablishVariableTypeDependentProperties( static_cast<int>(vType), size, added_prop.key );
         EstablishPlacementDependentProperties( added_prop.key.place, added_prop.key );
      }
      
      // update indices and binary database
      UpdateParametersAndDatabase();

      // update existing index references
      UpdateIndexReferences();
      
      return StorageKey( new_prop.c_str() );          
 } // end AddProperty


/**

Adds a new physical variable at runtime to the property database. The
new variable is temporary in the sense that it will not be added to the
binary variable database file.

@warning Any runtime change in variables invalidates existing Index objects! The need to be refreshed.

@section arguments Input Arguments

The method takes four input arguments, the name of the new variable,
the unit descriptor, the variable type and its placement in the mesh.
This overload only accecpts SCALAR, VECTOR, TENSOR variables.
To add an ARRAY, use the corresponding overload.

@section arguments Output Arguments

A csmp::Index for the new property.

@section implementation Implementation

The method uses the interface of the Parameter class to define the
new record for the physical variable.

@section application Application

The method is used in the constructor of the PropertyHandle.

@section messages Messages

The method will report a warning and return without adding a new variable,
if another variable with the same name already exists in the property
database.

@return The Index reflecting the new state ov LocalVariables and IntegrationPointVariables, including added property.
*/
template<uint32_t dim>
csmp::Index PropertyDatabase<dim>::AddProperty( const char *property_name, const char* notation, const char *unit,
                                                VARIABLE_TYPE vtype, PLACEMENT vplace, uint32_t vsize,
                                                double vmin, double vmax, string usage )
{
    return AddProperty( property_name, notation, unit, VariableCount( vplace, vtype) , vtype, vplace, vsize, vmin, vmax, usage );
}


/**
 
Adds a new physical variable at runtime to the property database. The
new variable is temporary in the sense that it will not be added to the
binary variable database file. 

@warning Any runtime change in variables invalidates existing Index objects! They need to be refreshed.

@section arguments Input Arguments

The method takes five input arguments, the name of the new variable,
the unit descriptor, the csmp::Index.index integer, the variable type and 
its placement in the mesh. The csmp::Index.index integer must be established
by a call to ScalarProperties() or its vector or tensor variable equivalents. 
These methods will return the number of scalar, vector, or tensor properties,
respectively, and the user can now use this integer number as
index for the new variable which is defined.  

@section arguments Output Arguments

A csmp::Index for the new property. 

@section implementation Implementation

The method uses the interface of the Parameter class to define the 
new record for the physical variable. By default, this sets the depth/components
to 1(SCALAR), 2(VECTOR), 3(TENSOR), 4(ARRAY), 5(FLAGGEDARRAY).

@section application Application

The method is used in the constructor of the PropertyHandle. 

@section messages Messages

The method will report a warning and return without adding a new variable,
if another variable with the same name already exists in the property
database. 

@return The Index reflecting the new state ov LocalVariables and IntegrationPointVariables, including added property.

*/
template<uint32_t dim>
csmp::Index  PropertyDatabase<dim>::AddProperty( const char* name, const char* notation, const char* unit, uint32_t index,
                                                 VARIABLE_TYPE vtype, PLACEMENT place, uint32_t vsize,
                                                 double vmin, double vmax, string usage )
 {
    auto iter(propList_.find(string(name)));

    // checking for uniqueness of name
    if ( iter != propList_.end() ) {
         ErrorHandler::Instance().Note( WARNING, "PropertyDatabase<dim>::AddProperty:",
                                        name, "property already exists. Nothing was done.");
         (*iter).second.Out();
         return StorageKey(name);
      }

    // checking for uniqueness of notation
    {
      set<string>  variable_abbreviations;
      for ( const auto& it : propList_ ) variable_abbreviations.insert( it.second.notation );
      if ( variable_abbreviations.size() < propList_.size() )
           ErrorHandler::Instance().Note( WARNING, "PropertyDatabase<dim>::AddProperty:",
                                         "notation of stored variables in not unique");
                                         
      if ( variable_abbreviations.find(notation) != variable_abbreviations.end() ) {
           ErrorHandler::Instance().Note( WARNING, "PropertyDatabase<dim>::AddProperty:",
                                          notation, "notation for new variable is not unique; aborting operation");
           (*iter).second.Out();
           return csmp::Index();
        }
    }
   csmp::Parameter  added_prop;
   assert( strlen(name) > 0 );
   added_prop.name      = name;
   assert( strlen(notation) > 0 );
   added_prop.notation  = notation;
   added_prop.key.place = place;
   added_prop.key.type  = vtype;
   added_prop.key.index = index;
   added_prop.unit      = unit;
   added_prop.min       = vmin;
   added_prop.max       = vmax;
   added_prop.usage =usage;
   added_prop.reference ="not specified";
   added_prop.explanation ="new property defined at runtime";
   if (this->Verbose()) cout <<"\nINFO, PropertyDatabase<dim>::AddProperty adding new property: '"<< name <<"'\n";

   /// Roman,2013: Added explicit way of reading the size of variable
   EstablishVariableTypeDependentProperties( static_cast<int>(vtype), vsize, added_prop.key );
   EstablishPlacementDependentProperties(  added_prop.key.place, added_prop.key );

   propList_[ added_prop.name ] = added_prop;

   // update indices and binary database
   UpdateParametersAndDatabase();

   // update existing index references
   UpdateIndexReferences();

  // returning the up-to-date storage key for the newly created variable
  return StorageKey(name);

} // end AddProperty


/** Deletes the target property from the property database.

@warning Existing Index objects may be invalidated after invoking this function

@param s The name of the physical variable whose database entry shall be deleted.

@section messages Messages

A warning will be issued, if the target variable does not exist in the
property database.

*/
template<uint32_t dim>
void PropertyDatabase<dim>::DeleteProperty( const char* s ) 
 {
    string pName(s);
    auto iter(propList_.find(pName));

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( iter == propList_.end() ) { 
         csmp_error.Note( WARNING, "PropertyDatabase<dim>::DeleteProperty", const_cast<char*>(s), 
                                     "target property was not defined");
         return;
      }
    else {
         cout <<"\nINFO, PropertyDatabase<dim>::DeleteProperty: deleting '"<< s <<"'"<< endl;
         // memorizing where the variable was stored
         csmp::Index key=StorageKey(s);
      
         // deleting the property
         propList_.erase(iter);

         // update indices and binary database
         DetachIndices(pName);
         UpdateParametersAndDatabase();  

         // update existing index references
         UpdateIndexReferences();
      }                
 } // end DeleteProperty




/** Creates a list of Index class objects corresponding to the current
set of physical variables stored in the property database.

The first argument is a reference to an STL list of Index class
objects. If this list already contains data, it is erased before the
new indices are entered. 

*/
template<uint32_t dim>
void PropertyDatabase<dim>::ListKeys( list<csmp::Index>& keys ) const
 {
    keys.clear();
    for ( auto& prop : propList_ )
      keys.push_back( StorageKey( prop.second.name.c_str() ) );
 }


/** Creates a map of pairs of variable names and corresponding Index class objects. .
The variable names will be alphabetically ordered.

The first argument is a reference to the map which shall be filled with
variable records. If this map already contains any data, it will be 
erased before the new physical variable data are entered. 

*/
template<uint32_t dim>
void PropertyDatabase<dim>::ListVariables( map<string,csmp::Index>& props ) const
 {
    props.erase( props.begin(), props.end() );

    for ( auto& prop : propList_ )
      props[ prop.second.name ] = StorageKey( prop.second.name.c_str() );
 }

/// enlists properties with a specific placement
template<uint32_t dim>
uint32_t PropertyDatabase<dim>::ListVariables( PLACEMENT pl, map<string,csmp::Index>& props ) const
 {
    props.clear();

    for ( auto& prop : propList_ )
      {
        csmp::Index key = StorageKey( prop.second.name.c_str() );
        if( key.place == pl )
          props[ prop.second.name ] = key;
      }
    return static_cast<uint32_t>(props.size());
 }

/// reports properties with the target placement as a set
template<uint32_t dim>
uint32_t PropertyDatabase<dim>::ListVariables( PLACEMENT place, set<string>& props ) const
 {
    props.clear();

    for ( auto& prop : propList_ )
      if ( prop.second.key.place == place )
        props.insert( prop.second.name );

    return static_cast<uint32_t>(props.size());
 }

/// reports properties with a specific placement and type as a set
template<uint32_t dim>
uint32_t PropertyDatabase<dim>::ListVariables( PLACEMENT place, VARIABLE_TYPE vtype, set<string>& props ) const
  {
  props.clear();

  for ( auto& prop : propList_ )
    if ( prop.second.key.place == place && prop.second.key.type == vtype )
      props.insert( prop.second.name );
 
      return static_cast<uint32_t>(props.size());
 }


/**

Returns the unit record string of the Parameter record which describes 
the target variable. 

@param s The name of the target physical variable.

@return The unit string record in the Parameter object which stores the
specifications of the physical variable. If the variable record
cannot be found, a NULL pointer is returned. 

@section messages Messages

An error is reported if the variable cannot be found in the property
database. 
*/
template<uint32_t dim>
const char*  PropertyDatabase<dim>::Unit( const char* s ) const
  {
     auto iter = propList_.find(string(s));
     
     if ( iter != propList_.end() ) return (*iter).second.unit.c_str();
     else
     cout <<"\n\nPropertyDatabase::Unit: Unable to identify unit entry of: " << s << endl;

     return "undefined";  
     
  } // end Unit 


/**
 
Reports conversion factors which allow the user to convert variable values
from the SI to the CGS system or vice versa. 

@section arguments Input Arguments

The method expects three string arguments, the name of the system from
which the variable is converted (options are SI or CGS), the name of the 
system into which the variable is converted (SI or CGS), and the unit
which is converted. The SI unit options are the following: 
 
"oC","m","MPa","kg m-3","yr","m Ma-1","W m-2","m2","J m-1 K-1 s-1","kg m-2 s-1"
 
and the available CGS (and other) units which can be converted are:
 
"cm", "cm2", "cm-2", "m","feet","dyne cm-2 s-1","g cm-1 s-1", "g cm-2 s-1", 
"g cm-3", "cal", "cal g-1 oC-1", "cal cm-1 s-1 oC-1","cal cm-2 s-1","oC", 
"bars", "cm s-1".

@return  A floating point conversion factor.

@section application Application

To get conversion factors which allow to output variables in non SI
units for comparison with simulation results from non CSP simulation
tools. 

@section messages Messages

The method will report an error if it cannot identify one of its string
arguments. 

*/
template<uint32_t dim>
double  PropertyDatabase<dim>::UnitConversionFactor( const char* sys_in, 
                                                  const char* sys_out, 
                                                  const char* unit ) const
 {
   if ( (strcmp( "SI", sys_in )) == 0 ) 
      {
         if ( (strcmp( "CGS", sys_out )) == 0 )
           {
             if      ( (strcmp( "oC",         unit )) == 0 ) return 1.0;
             else if ( (strcmp( "m",          unit )) == 0 ) return 1.0e-3;
             else if ( (strcmp( "MPa",        unit )) == 0 ) return 1.0e+1;
             else if ( (strcmp( "kg m-3",     unit )) == 0 ) return 1.0e+3;
             else if ( (strcmp( "yr",         unit )) == 0 ) return 30758400.0;
             else if ( (strcmp( "m Ma-1",     unit )) == 0 ) return 3.2511444e-14;
             else if ( (strcmp( "W m-2",      unit )) == 0 ) return 0.239006e-4;
             else if ( (strcmp( "m2",         unit )) == 0 ) return 1.0e+4;
             else if ( (strcmp( "J m-1 K-1 s-1", unit )) == 0 ) return 0.239006e-2;
             else if ( (strcmp( "kg m-2 s-1", unit )) == 0 ) return 1.0e-1;
             else 
              {
                cout <<"\nPropertyDatabase::convert_unit: Failure to convert " << unit;
                cout <<" from SI to CGS format.";
              }
           }
      } // end if SI

    else if ( (strcmp( "CGS", sys_in )) == 0 ) 
      {
         if ( (strcmp( "SI", sys_out )) == 0 )
           {
             // conversion factors from:
             if      ( (strcmp( "cm",             unit )) == 0 ) return 1.0e+2;
             else if ( (strcmp( "cm2",            unit )) == 0 ) return 1.0e-4;
             else if ( (strcmp( "cm-2",           unit )) == 0 ) return 1.0e+4;
             else if ( (strcmp( "m",              unit )) == 0 ) return 1.0e+3;
             else if ( (strcmp( "feet",           unit )) == 0 ) return 0.3048;
             else if ( (strcmp( "dyne cm-2 s-1",  unit )) == 0 ) return 1.0e-1;
             else if ( (strcmp( "g cm-1 s-1",     unit )) == 0 ) return 1.0e-1;
             else if ( (strcmp( "g cm-2 s-1",     unit )) == 0 ) return 1.0e+1;
             else if ( (strcmp( "g cm-3",         unit )) == 0 ) return 1.0e+3;
             else if ( (strcmp( "cal",            unit )) == 0 ) return 4.184;
             else if ( (strcmp( "cal g-1 oC-1",   unit )) == 0 ) return 4.184e-3;
             else if ( (strcmp( "cal cm-1 s-1 oC-1",unit))== 0 ) return 4.184e-2;
             else if ( (strcmp( "cal cm-2 s-1",   unit )) == 0 ) return 4.184e-4;
             else if ( (strcmp( "oC",             unit )) == 0 ) return 1.0;
             else if ( (strcmp( "bars",           unit )) == 0 ) return 1.0e-1;
             else if ( (strcmp( "cm s-1",         unit )) == 0 ) return 1.0e-2;
             else 
              {
                cerr <<"\nPropertyDatabase::convert_unit: Failure to convert " << unit;
                cerr <<" from CGS to SI" << endl;
                return std::numeric_limits<double>::signaling_NaN();
              }
           }
      } // end if SI
 
   return std::numeric_limits<double>::signaling_NaN();

 } // end ConvertUnit


/**
 
Outputs the physically meaningful range of the target variable. This range
is equivalent to the range which was specified in the variable database
file. 

@param s The name of the physical variable whose range shall be established.

@param mn The second and third method arguments will store the mininum and maximum
permitted values of the target variable. 

@section messages Messages

If the target variable has not been defined in the property database,
an error will be reported. 
*/
template<uint32_t dim>
void PropertyDatabase<dim>::RangeOf( const char* s, double& mn, double& mx ) const
 {
     auto iter = propList_.find(string(s));
     if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
     else cout <<"\nPropertyDatabase::RangeOf: Unable to identify range of: " << s << endl;

 } //  end RangeOf



template<uint32_t dim>
double PropertyDatabase<dim>::LowerLimitOf( const char* property ) const
 {
     auto iter = propList_.find(string(property));
     if ( iter != propList_.end() ) return (*iter).second.MinValue();
     else cout <<"\nPropertyDatabase::LowerLimitOf: Unable to identify minimum value of: " << property << endl;
     return numeric_limits<double>::quiet_NaN();

 } //  end 


template<uint32_t dim>
double PropertyDatabase<dim>::UpperLimitOf( const char* property ) const
 {
     auto iter = propList_.find(string(property));
     if ( iter != propList_.end() ) return (*iter).second.MaxValue();
     else cout <<"\nPropertyDatabase::LowerLimitOf: Unable to identify maximum value of: " << property << endl;
     return numeric_limits<double>::quiet_NaN();

 } //  end 



/**
 
Tests whether the supplied value of the target variable is within the 
physical meaningful range which was specified for the variable in the 
database file.  

@param s The first method argument give the name of the target variable.

@param var The second method argument specifies the actual physical-variable value
which will be overwritten by the closest range bound if the target 
variable is out of range. 

@section implementation Implementation

Checks range of variable *s. A warning is issued if the variable is out 
of range and the variable is set to the closest range bound. 

@section messages Messages

If the target variable has not been defined in the property database,
an error will be reported. 
*/
template<uint32_t dim>
void PropertyDatabase<dim>::CheckRange( const char* s, double& var ) const
 {
    double mn, mx;
    auto iter = propList_.find(string(s));
    ErrorHandler& error_handler ( ErrorHandler::Instance() );
    if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
    else throw csmp::Exception( FATAL_ERROR, "PropertyDatabase<dim>::CheckRange", "property could not be identified");
    std::stringstream ss;
    string msg1,msg2,errmsg;
     if ( var < mn )
         {
            //cout <<"\n Error: PropertyDatabase<dim>::CheckRange:value too low"<< endl;
            //cout <<"\t User defined '"<< s <<"' minimum: "<< mn <<", versus: "<< var << endl;
            //var = mn;
            ss << var << ' ' << mn;
            ss >> msg1 >> msg2;
            errmsg=" variable : "+string(s)+" user defined : "+msg1+" while minimum was established at: "+msg2;
            error_handler.Note(FATAL_ERROR,"PropertyDatabase<dim>::CheckRange"," value is below the database minimum.",errmsg.c_str());
         }
      else if ( var > mx )
         {
            //cout <<"\n Error: PropertyDatabase<dim>::CheckRange: value too high"<< endl;
            //cout <<"\t User defined'"<< s <<"' maximum: "<< mx <<", versus: "<< var << endl;
            //var = mx;
            ss<<var<<' '<< mx;
            ss>>msg1>>msg2;
            errmsg=" variable : "+string(s)+" user defined : "+msg1+" while maximum was established at: "+msg2;
            error_handler.Note(FATAL_ERROR,"PropertyDatabase<dim>::CheckRange"," value is below the database minimum.",errmsg.c_str());
         }
        
 } // end CheckRange



/**
   returns true if tha value of the supplied variable is within the range stored in the database, else false
*/
template<uint32_t dim>
bool PropertyDatabase<dim>::CheckRange( const char* s, const ScalarVariable& var ) const
 {
    double mn, mx;
    auto iter = propList_.find(string(s));
    if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
    else throw csmp::Exception( FATAL_ERROR, "PropertyDatabase<dim>::CheckRange(scalar)", "property could not be identified");
    if ( isnan(var()) ) return false;
    if ( var() < mn ) return false;
    if ( var() > mx ) return false;

    return true;
   
 } // end CheckRange (scalar)


template<uint32_t dim>
bool PropertyDatabase<dim>::CheckRange( const char* s, const VectorVariable<dim>& var ) const
 {
    double mn, mx;
    auto iter = propList_.find(string(s));
    if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
    else throw csmp::Exception( FATAL_ERROR, "PropertyDatabase<dim>::CheckRange(vector)", "property could not be identified");
    if ( isnan(var.Length()) ) return false;
    if ( -var.Length() < mn ) return false;
    if ( var.Length() > mx ) return false;

    return true;
   
 } // end CheckRange (vector)



template<uint32_t dim>
bool PropertyDatabase<dim>::CheckRange( const char* s, const TensorVariable<dim>& var ) const
 {
    double mn, mx;
    auto iter = propList_.find(string(s));
    if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
    else throw csmp::Exception( FATAL_ERROR, "PropertyDatabase<dim>::CheckRange(tensor)", "property could not be identified");
    const double det = var.Determinant();
    if ( isnan(det) ) return false;
    if ( det < mn ) return false;
    else if ( det > mx ) return false;

    return true;
   
 } // end CheckRange (tensor)




template<uint32_t dim>
bool PropertyDatabase<dim>::CheckRange( const char* s, const ArrayVariable& var ) const
 {
    double mn, mx;
    auto iter = propList_.find(string(s));
    //ErrorHandler& error_handler ( ErrorHandler::Instance() );
    if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
    else throw csmp::Exception( FATAL_ERROR, "PropertyDatabase<dim>::CheckRange(array)", "property could not be identified");
    double ary_min(var[0]), ary_max(var[0]);
    for ( uint32_t i=1U; i<var.Size(); ++i ) {
        if ( isnan(var[i]) ) return false;
        ary_min = min( ary_min, var[i] );
        ary_max = max( ary_max, var[i] );
      }
    if ( ary_min < mn ) return false;
    else if ( ary_max > mx ) return false;

    return true;
   
 } // end CheckRange (array)




template<uint32_t dim>
bool PropertyDatabase<dim>::CheckRange( const char* s, const FlaggedArrayVariable& var ) const
 {
    double mn, mx;
    auto iter = propList_.find(string(s));
    //ErrorHandler& error_handler ( ErrorHandler::Instance() );
    if ( iter != propList_.end() ) (*iter).second.Range( mn, mx );
    else throw csmp::Exception( FATAL_ERROR, "PropertyDatabase<dim>::CheckRange(flagged array)", "property could not be identified");
    double ary_min(var[0]), ary_max(var[0]);
    for ( uint32_t i=1U; i<var.Size(); ++i ) {
        if ( isnan(var[i]) ) return false;
        ary_min = min( ary_min, var[i] );
        ary_max = max( ary_max, var[i] );
      }
    if ( ary_min < mn ) return false;
    else if ( ary_max > mx ) return false;

    return true;
   
 } // end CheckRange (array)




/**
     Sets of the range of the target variable at runtime in case relevant information becomes available.
*/
template<uint32_t dim>
void PropertyDatabase<dim>::SetRangeOf( const char* property_name, double vmin, double vmax )
 {
    auto iter = propList_.find(string(property_name));
    ErrorHandler& error_handler ( ErrorHandler::Instance() );
    if ( iter != propList_.end() )
      (*iter).second.Range( vmin, vmax );
    else
      error_handler.Note( ERROR, "PropertyDatabase<dim>::SetRangeOf:", property_name, "property could not be identified.");
 }




/** Finds the name of the physical variable on the basis of its Index.

@param idx The Index of the target physical variable whose name shall be
retrieved. 

@return The name of the physical variable returned in C-string format.
If the variable has not been defined in the property database, the 
output string will be 'undefined'.
*/
template<uint32_t dim>
const char* PropertyDatabase<dim>::Name( const csmp::Index& idx ) const
 {
    for ( auto& prop : propList_ )
      if ( prop.second.key == idx ) return prop.second.name.c_str();

    return "undefined";
   
 } 





template<uint32_t dim>
void PropertyDatabase<dim>::Out() const
 {
     cout <<"\nPropertyDatabase::Out: variables file '"<< physvarsFile <<"'"<< endl;
     cout <<"total number of stored properties: "<< VariableCount() << endl;
     if ( VariableCount() > 0 ) {
          cout <<"detailed variable counts:\n";
          for( map<PLACEMENT,map<VARIABLE_TYPE,uint32_t> >::const_iterator it( VariableCountBegin() ); it != VariableCountEnd(); ++it )
           for( map<VARIABLE_TYPE,uint32_t>::const_iterator iit( it->second.begin() ); iit != it->second.end(); ++iit )
             if ( VariableCount( it->first, iit->first ) > 0 )
               cout <<"\t"<< parseType(iit->first) << " variables on " << parsePlacement(it->first) << ": "  << VariableCount( it->first, iit->first ) << endl;

          cout <<"\nDetailed information on current properties in alphabetical order: "<< endl;
          for ( auto& prop : propList_ )
            cout << prop.first << ": \n" << prop.second;
          cout << endl;
          cout.flush();
       }
     
 }  // end Out


/// writes a csmp variables text file. NB: parameter is fully specified fileName(incl extension)
template<uint32_t dim>
bool PropertyDatabase<dim>::WriteVariablesFile( const char* fileName ) const
  {
    ofstream variablesFile( fileName );
    if ( variablesFile.is_open() )
      {
        // header
        variablesFile << "name\tMParameter\tunit\tindex\tmin.\tmax.\tplace\n";
        // properties
        for( auto& prop : propList_ )
          {
            const csmp::Index id( StorageKey( prop.second.name.data() ) );
            assert( prop.second.name == prop.first );
            string placementString( parsePlacement( id.place ) ); 
            
            variablesFile << prop.second.name << "\t" << prop.second.notation << "\t" << prop.second.unit << "\t"; 
              if( id.type != ARRAY && id.type != FLAGGEDARRAY )
                variablesFile << id.type;
              else
                variablesFile << id.dataDepth;
              variablesFile << "\t" << prop.second.min << "\t" << prop.second.max << "\t" << placementString;
            variablesFile << endl;
          } // all parameters
      } // if file

    variablesFile.close();
    return true;
  }




/**
   Writes templatized INDEX variable definitions of the current variables  into header file and instantiates their keys.
   If the variable set name is empty, the variables are inserted into the global namespace.
*/
template<uint32_t dim>
void PropertyDatabase<dim>::WriteVariableSetToHeaderFile( const char* header_file, const char* variable_set_name ) const
 {
    ErrorHandler& csmp_error{ ErrorHandler::Instance() };
 
    ofstream ofs( header_file );
    if ( !ofs.is_open() ) {
         csmp_error.Note( ERROR, "PropertyDatabase<dim>::WriteVariableSetToHeaderFile",
                          header_file, "input file could not be opened");
         return;
      }

    // 1. writing the file header and time when this file was created
    // --------------------------------------------------------------
    time_t now_time = chrono::system_clock::to_time_t(chrono::system_clock::now());

    if ( !variable_set_name )
      ofs <<"\\"<< header_file <<": written by PropertyDatabase<dim>::WriteVariableSetToHeaderFile: "<< ctime(&now_time) <<"\n";
    else {
         ofs <<"\\"<< header_file <<": variable set: "<< variable_set_name;
         ofs <<"written by PropertyDatabase<dim>::WriteVariableSetToHeaderFile: " << ctime(&now_time) <<"\n";
      }
    ofs << endl;
      
    // 2. Instantiating the variables, either globally or locally, creating names from their notation which must be unique
    // -------------------------------------------------------------------------------------------------------------------
    if ( variable_set_name ) {
         ofs <<"struct "<< variable_set_name <<" {";
      }
    // testing whether notation is unique
    size_t      n_variables{ propList_.size() };
    set<string> unique_variable_identifiers;
    for ( const auto& it : propList_ )
      unique_variable_identifiers.insert( it.second.notation );
        
    if ( unique_variable_identifiers.size() < n_variables ) {
         csmp_error.Note( ERROR, "PropertyDatabase<dim>::WriteVariableSetToHeaderFile",
                          header_file, "contained variable notation is non unique; keys could not be created");
         return;
      }

    // writing the variable keys
    for ( const auto& it : propList_ ) {
         ofs <<"csmp::INDEX<"<< parseType(it.second.key.type) <<","<< parsePlacement(it.second.key.place);
         ofs <<"> key_"<< it.second.notation <<"; \\ '"<< it.first << endl;
      }


    // 3. Providing a function or method to register the variables with the PropertyDatabase
    // -------------------------------------------------------------------------------------
    // function
    if ( !variable_set_name ) {
         ofs <<"\n"<<"template<uint32_t dim>"<< endl;
         ofs <<"void registerGlobalVariableINDEX_KeysWithPropertyDatabase( PropertyDatabase<dim>& db ) {"<< endl;
         ofs <<"double vmin, vmax;"<< endl;
         ofs <<"for ( const auto& it : propList_ ) {"<< endl;
         ofs <<"     it.second.Range( vmin, vmax );"<< endl;
         ofs <<"     db.AddProperty( it.first.c_str(), /* it.second.notation */ it.second.unit,"<< endl;
         ofs <<"                     it.second.key.type, it.second.key.place, it.second.key.dataDepth, vmin, vmax );"<< endl;
         ofs <<"  }"<< endl << endl;
         ofs <<"}"<< endl;
      }
    // member constructor of the variable set that registers it with the PropertyDatabase
    ofs <<"\n"<<"template<uint32_t dim>"<< endl;
    ofs <<"explicit "<< variable_set_name <<"( const PropertyDatabase<dim>& db )\n :"<< endl;
    // instantiating the variable keys
    for ( const auto& it : propList_ ) {
         ofs <<" key_"<< it.second.notation;
         ofs <<"( INDEX<"<< parseType(it.second.key.type) <<","<< parsePlacement(it.second.key.place) <<">( ";
         ofs <<"db.StorageKey\""<< it.first <<"\") ))"<< endl;
      }
   // finishing the in-class definition of the constructor
   ofs <<"  {"<< endl <<"  }"<< endl;

    // end of struct definition
    if ( variable_set_name ) {
         ofs <<"\n };"<< endl;
      }

 } // end WriteVariableSetToHeaderFile







/// Clears arrayLengths and inserts depth of ARRAY variables for place in order of index
template<uint32_t dim>
void PropertyDatabase<dim>::ArrayLengths( PLACEMENT place, std::vector<uint32_t>& arrayLengths ) const
  {
    arrayLengths.clear();
    map<uint32_t,uint32_t> offsetIndexMap;
    for ( auto& prop : propList_ )
      if ( prop.second.key.type  == ARRAY && prop.second.key.place == place  )
        offsetIndexMap[prop.second.key.index] = prop.second.key.dataDepth;
    for( auto& offs : offsetIndexMap )
      arrayLengths.push_back( offs.second );  
  }


template<uint32_t dim>
uint32_t PropertyDatabase<dim>::ArrayLengthTotal( PLACEMENT place ) const
  {
    vector<uint32_t> lengths;
    ArrayLengths( place, lengths );
    uint32_t totalLength(0);
    for( uint32_t i(0); i < lengths.size(); ++i )
      totalLength += lengths[i];
    return totalLength;
  }

/// Clears arrayLengths and inserts depth of FLAGGEDARRAY variables for place in order of index
template<uint32_t dim>
void PropertyDatabase<dim>::FlaggedArrayLengths( PLACEMENT place, std::vector<uint32_t>& arrayLengths ) const
  {
    arrayLengths.clear();
    map<uint32_t,uint32_t> offsetIndexMap;
    for ( auto& prop : propList_ )
      if ( prop.second.key.type  == FLAGGEDARRAY && prop.second.key.place == place  )
        offsetIndexMap[prop.second.key.index] = prop.second.key.dataDepth;
    for( auto& offs : offsetIndexMap )
      arrayLengths.push_back( offs.second );
  }


template<uint32_t dim>
uint32_t PropertyDatabase<dim>::FlaggedArrayLengthTotal( PLACEMENT place ) const
  {
    vector<uint32_t> lengths;
    FlaggedArrayLengths( place, lengths );
    uint32_t totalLength(0);
    for( uint32_t i(0); i < lengths.size(); ++i )
      totalLength += lengths[i];
    return totalLength;
  }


template<uint32_t dim>
void PropertyDatabase<dim>::UpdateParametersAndDatabase()
  {
    // updating the numbering of variables
    CountVariables();

    // assigning the variable indices to the properties in the list
    AssignVariableIndices();

    // writing local and integration point variables
    EstablishIndexLocalAndIntegrationPointVariables();

    // assigning flag & data offsets for arrays
    EstablishIndexOffsets();

    // attached indices to IndexTracker and vice versa
    AttachIndices();
  }


/// Creates IndexTracker <--> Index connection and updates offsets of all Index references
template<uint32_t dim>
void PropertyDatabase<dim>::AttachIndices()
  {
    // attach
    for( auto& prop : propList_ )
      {
        indexTracker_.Attach( &(prop.second.key), prop.second.name );
        prop.second.key.Attach(&indexTracker_);
      }
  }


/// Updates all register Index objects
template<uint32_t dim>
void PropertyDatabase<dim>::UpdateIndexReferences()
  {
    uint32_t i(0);
    // update
    for ( auto it(indexTracker_.IndicesBegin()); it != indexTracker_.IndicesEnd(); ++it, ++i )
      if ( IsDefined( it->second.c_str() ) ) {
           // (it->first)->UpdateData( propList_[it->second].key ); /// @todo (2-F) Buggy, only safe by design
           //            string,Parameter pairs
           auto key_it = propList_.find( (*it).second );
           if ( key_it == propList_.end() ) {
                cerr <<"\n"<< (*it).second;
                throw csmp::Exception( ERROR, "PropertyDatabase<dim>::UpdateIndexReferences:", "could not find variable in database when searching by name." );
             }
          else if ( (*it).first != nullptr )
            (*it).first->UpdateData( (*key_it).second.key );
        }
      else
        throw csmp::Exception( ERROR, "PropertyDatabase<dim>::UpdateIndexReferences:", "IndexTracker has Index objects that are linked to undefined properties." );
    
    if (this->Verbose()) cout << "\nPropertyDatabase<dim>::UpdateIndexReferences: updated " << i << " Index objects using IndexTracker.\n";
  }


/// Removes IndexTracker <--> Index connection
template<uint32_t dim>
void PropertyDatabase<dim>::DetachIndices( string parameterName )
  {
    indexTracker_.Detach(parameterName);  
  }


/// Sets local and integration point variables in indices for all parameters
template<uint32_t dim>
void PropertyDatabase<dim>::EstablishIndexLocalAndIntegrationPointVariables()
  {
    for( auto& prop : propList_ )
      {
        prop.second.key.localVariables = LocalVariablesAt( prop.second.key.place );
        prop.second.key.integrationPointVariables = IntegrationPointVariablesAt( prop.second.key.place );
      }
  }


template<uint32_t dim>
void PropertyDatabase<dim>::EstablishScalarOffsets( PLACEMENT whithin )
  {
    ///  Scalars: data offset = flag offset = index (as many flags as components)
    for( auto& prop : propList_ )
      if( prop.second.key.place == whithin && prop.second.key.type == SCALAR )
        prop.second.key.dataOffset = prop.second.key.flagOffset = prop.second.key.index;
  }


template<uint32_t dim>
void PropertyDatabase<dim>::EstablishVectorOffsets( PLACEMENT whithin )
  {
    const uint32_t scalars( VariableCount( whithin, SCALAR ) );

    ///  Vectors: data offset = flag offset = scalars at this placement plus index times dim (as many flags as components)
    for( auto& prop : propList_ )
      if( prop.second.key.place == whithin && prop.second.key.type == VECTOR )
        prop.second.key.dataOffset = prop.second.key.flagOffset = (scalars + prop.second.key.index*dim);
  }


template<uint32_t dim>
void PropertyDatabase<dim>::EstablishTensorOffsets( PLACEMENT whithin )
  {
    const uint32_t scalars( VariableCount( whithin, SCALAR ) );
    const uint32_t vectors( VariableCount( whithin, VECTOR ) );
    const uint32_t offsetToTensors( scalars + vectors*dim );

    ///  Tensors: data offset != flag offset (as many flags as dim)
    for( auto& prop : propList_ )
      if( prop.second.key.place == whithin && prop.second.key.type == TENSOR )
        {
          prop.second.key.dataOffset = offsetToTensors + prop.second.key.index*dim*dim;
          prop.second.key.flagOffset = offsetToTensors + prop.second.key.index*dim;
        }
  }


template<uint32_t dim>
void PropertyDatabase<dim>::EstablishArrayOffsets( PLACEMENT whithin )
  {
    const uint32_t scalars( VariableCount( whithin, SCALAR ) );
    const uint32_t vectors( VariableCount( whithin, VECTOR ) );
    const uint32_t tensors( VariableCount( whithin, TENSOR ) );
    const uint32_t offsetToArrayData( scalars + vectors*dim + tensors*dim*dim );
    const uint32_t offsetToArrayFlag( scalars + vectors*vectorFlags + tensors*tensorFlags );

    // we need to sort first acc to indexes
    map<uint32_t,csmp::Parameter*> sortedArrays;
    for( auto& prop : propList_ )
      if( prop.second.key.place == whithin && prop.second.key.type == ARRAY )
        sortedArrays[prop.second.key.index] = &(prop.second);

    // then we calculate the relative offset depending on the order and size of all arrays
    uint32_t currentDataOffset(0);
    for( auto& sarr : sortedArrays )
      {
        sarr.second->key.dataOffset = currentDataOffset;
        currentDataOffset += sarr.second->key.dataDepth;
        sarr.second->key.flagOffset = sarr.second->key.index;
      }

    // we adjust for the total offset resulting from scalars, vectors and tensors ahead
    for( auto& sarr : sortedArrays )
      {
        sarr.second->key.dataOffset += offsetToArrayData;
        sarr.second->key.flagOffset += offsetToArrayFlag;
      }
  }

template<uint32_t dim>
void PropertyDatabase<dim>::EstablishFlaggedArrayOffsets( PLACEMENT whithin )
  {
    const uint32_t scalars( VariableCount( whithin, SCALAR ) );
    const uint32_t vectors( VariableCount( whithin, VECTOR ) );
    const uint32_t tensors( VariableCount( whithin, TENSOR ) );
    const uint32_t arrays ( VariableCount( whithin, ARRAY ) );
    const uint32_t arraysDepth( ArrayLengthTotal  ( whithin ) );

    const uint32_t offsetToFlaggedArrayData( scalars + vectors*dim + tensors*dim*dim + arraysDepth);
    const uint32_t offsetToFlaggedArrayFlag( scalars + vectors*vectorFlags + tensors*tensorFlags + arrays );

    // we need to sort first acc to indexes
    map<uint32_t,csmp::Parameter*> sortedArrays;
    for( auto& prop : propList_ )
      if( prop.second.key.place == whithin && prop.second.key.type == FLAGGEDARRAY )
        sortedArrays[prop.second.key.index] = &(prop.second);

    // then we calculate the relative offset depending on the order and size of all arrays
    uint32_t currentOffset(0);
    for( auto& sarr : sortedArrays )
      {
        sarr.second->key.dataOffset = currentOffset;
        sarr.second->key.flagOffset = currentOffset;
        currentOffset += sarr.second->key.dataDepth;
      }

    // we adjust for the total offset resulting from scalars, vectors and tensors ahead
    for( auto& sarr : sortedArrays )
      {
        sarr.second->key.dataOffset += offsetToFlaggedArrayData;
        sarr.second->key.flagOffset += offsetToFlaggedArrayFlag;
      }
  }

template<uint32_t dim>
void PropertyDatabase<dim>::EstablishIndexOffsets()
  {
    for( auto it( VariableCountBegin() ); it != VariableCountEnd(); ++it )
      {
        EstablishScalarOffsets(it->first);   
        EstablishVectorOffsets(it->first);   
        EstablishTensorOffsets(it->first);   
        EstablishArrayOffsets(it->first);
        EstablishFlaggedArrayOffsets(it->first);
      } 
  }


/// Establishes integration point factors used in index arithmetic within LocalVariableStorage
template<uint32_t dim>
void PropertyDatabase<dim>::EstablishPlacementDependentProperties( PLACEMENT placement, csmp::Index& key ) const
  {
    if( placement == ELEMENT_INTEGRATION_POINT ||
        placement == FACE_INTEGRATION_POINT ||
        placement == INTER_FACE_INTEGRATION_POINT  )
      {
        key.offsetFactorSimplex = 0;
        key.offsetFactorSector  = 0; 
        key.ipFactorSimplex     = 1;
        key.ipFactorSector      = 0;
        key.ipFactorFacet       = 0;
      }
    else if( placement == SECTOR_INTEGRATION_POINT || 
             placement == FACE_SECTOR_INTEGRATION_POINT ||
             placement == INTER_FACE_SECTOR_INTEGRATION_POINT )
      {
        key.offsetFactorSimplex = 1;
        key.offsetFactorSector  = 0;
        key.ipFactorSimplex     = 0;
        key.ipFactorSector      = 1;
        key.ipFactorFacet       = 0;
      }
    else if( placement == FACET_INTEGRATION_POINT ||
             placement == FACE_FACET_INTEGRATION_POINT ||
             placement == INTER_FACE_FACET_INTEGRATION_POINT )
      {
        key.offsetFactorSimplex = 1;
        key.offsetFactorSector  = 1;
        key.ipFactorSimplex     = 0;
        key.ipFactorSector      = 0;
        key.ipFactorFacet       = 1;
      }
    else
      {
        key.offsetFactorSimplex = 0;
        key.offsetFactorSector  = 0; 
        key.ipFactorSimplex     = 0;
        key.ipFactorSector      = 0;
        key.ipFactorFacet       = 0;
      }
  }


/// Roman, 2013: Added the case of negative type for FlaggedArrayVariable
template<uint32_t dim>
void PropertyDatabase<dim>::EstablishVariableTypeDependentProperties( int vtype, csmp::Index& key ) const
  {
    if( vtype < 0 )
      {
        key.type = FLAGGEDARRAY;
        key.dataDepth = std::abs(vtype);
        key.flagDepth = key.dataDepth;
      }
    else if( vtype > 3 )
      {
        key.type = ARRAY;
        key.dataDepth = vtype;
        key.flagDepth = 1;
      }    
    else if( static_cast<VARIABLE_TYPE>(vtype) == SCALAR )
      {          
        key.type = SCALAR;
        key.dataDepth = 1; 
        key.flagDepth = 1;
      }
    else if( static_cast<VARIABLE_TYPE>(vtype) == VECTOR )
      {                 
        key.type = VECTOR;
        key.dataDepth = dim;   
        key.flagDepth = vectorFlags;  
      }
    else if( static_cast<VARIABLE_TYPE>(vtype) == TENSOR )
      {           
        key.type = TENSOR;
        key.dataDepth = dim*dim;   
        key.flagDepth = tensorFlags;  
      }
    else
      throw csmp::Exception( ERROR, "PropertyDatabase<dim>::EstablishVariableTypeDependentProperties(int,Index):",
                            "Variable type not supported/identified" );
  }





/// Roman, 2013: Added explicit way of establishing type and size of the variable dependent properties
template<uint32_t dim>
void PropertyDatabase<dim>::EstablishVariableTypeDependentProperties( int vtype, uint32_t vsize, csmp::Index& key ) const
  {
    if( static_cast<VARIABLE_TYPE>(vtype) == SCALAR )
      {
        key.type = SCALAR;
        key.dataDepth = 1;
        key.flagDepth = 1;
      }
    else if( static_cast<VARIABLE_TYPE>(vtype) == VECTOR )
      {
        key.type = VECTOR;
        key.dataDepth = dim;
        key.flagDepth = vectorFlags;
      }
    else if( static_cast<VARIABLE_TYPE>(vtype) == TENSOR )
      {
        key.type = TENSOR;
        key.dataDepth = dim*dim;
        key.flagDepth = tensorFlags;
      }
    else if( static_cast<VARIABLE_TYPE>(vtype) == ARRAY )
      {
        key.type = ARRAY;
        key.dataDepth = vsize;
        key.flagDepth = 1;
      }
    else if( static_cast<VARIABLE_TYPE>(vtype) == FLAGGEDARRAY )
      {
        key.type = FLAGGEDARRAY;
        key.dataDepth = vsize;
        key.flagDepth = key.dataDepth;
      }
    else
      throw csmp::Exception( ERROR, "PropertyDatabase<dim>::EstablishVariableTypeDependentProperties(int,size_t,Index):",
                            "Variable type not supported/identified" );
  }


/// Roman, 2013: Added new way of reading the type of the variable base on the name and size specification's
/// Available types:    SCALAR, VECTOR, TENSOR, ( the size defined automatically based on the dimension )
///                     ARRAY Size,             ( size of Array Variable with single flag specified by user, one should put space after the name ARRAY)
///                     FLAGGEDARRAY Size       ( size of Flagged Array Variable with multiple flag's specified by user, one should put space after the name FLAGGEDARRAY)
template<uint32_t dim>
void PropertyDatabase<dim>::EstablishVariableTypeDependentProperties( std::string vtype, csmp::Index& key ) const
  {
    std::string   vtype_name;
    unsigned long vtype_position;
    const char*   vtype_delims =" 1234567890";

    // Transform the name of the type to uppercase representation
    std::transform(vtype.begin(),vtype.end(),vtype.begin(),::toupper);

    vtype_position = vtype.find_first_of(vtype_delims);
    vtype_name     = vtype.substr(0,vtype_position);

    if( vtype_name      == "SCALAR" )
      {
        key.type = SCALAR;
        key.dataDepth = 1;
        key.flagDepth = 1;
      }
    else if( vtype_name == "VECTOR" )
      {
        key.type = VECTOR;
        key.dataDepth = dim;
        key.flagDepth = vectorFlags;
      }
    else if( vtype_name == "TENSOR" )
      {
        key.type = TENSOR;
        key.dataDepth = dim*dim;
        key.flagDepth = tensorFlags;
      }
    else if( vtype_name == "ARRAY" )
      {
        key.type = ARRAY;
        key.dataDepth = atoi( vtype.substr(vtype_position).c_str() );
        key.flagDepth = 1;
      }
    else if( vtype_name == "FLAGGEDARRAY" )
      {
        key.type = FLAGGEDARRAY;
        key.dataDepth = atoi( vtype.substr(vtype_position).c_str() );;
        key.flagDepth = key.dataDepth;
      }
    else
      throw csmp::Exception( ERROR, "PropertyDatabase<dim>::EstablishVariableTypeDependentProperties(string,Index):",
                            "Variable type not supported/identified", vtype.c_str() );
  }



template<uint32_t dim>
LocalVariables PropertyDatabase<dim>::LocalVariablesAt( PLACEMENT within ) const
  {
    if( within == ELEMENT_INTEGRATION_POINT ||
        within == FACET_INTEGRATION_POINT ||
        within == SECTOR_INTEGRATION_POINT  )
      return LocalVariablesAt(ELEMENT);
    else if( within == FACE_INTEGRATION_POINT || 
             within == FACE_SECTOR_INTEGRATION_POINT ||
             within == FACE_FACET_INTEGRATION_POINT )
      return LocalVariablesAt(FACE);
    else if( within == INTER_FACE_INTEGRATION_POINT ||
             within == INTER_FACE_SECTOR_INTEGRATION_POINT ||
             within == INTER_FACE_FACET_INTEGRATION_POINT )
      return LocalVariablesAt(INTER_FACE);

    const uint32_t scalars( VariableCount( within, SCALAR ) );
    const uint32_t vectors( VariableCount( within, VECTOR ) );
    const uint32_t tensors( VariableCount( within, TENSOR) );
    const uint32_t arrayDepth( ArrayLengthTotal(within) );
    const uint32_t arrays( VariableCount( within, ARRAY ) );
    const uint32_t flaggedArrayDepth  ( FlaggedArrayLengthTotal(within) );
    const uint32_t flaggedArrays      ( VariableCount( within, FLAGGEDARRAY ) );
    const uint32_t dataDepth( scalars + vectors*dim + tensors*dim*dim + arrayDepth + flaggedArrayDepth );
    const uint32_t flagDepth( scalars + vectors*vectorFlags + tensors*tensorFlags + arrays + flaggedArrayDepth);
    return LocalVariables( scalars, vectors, tensors, arrays, arrayDepth, flaggedArrays, flaggedArrayDepth, dataDepth, flagDepth );
  }


template<uint32_t dim>
IntegrationPointVariables PropertyDatabase<dim>::IntegrationPointVariablesAt( PLACEMENT within ) const
  {
    /// only for Element / Face / InterFace and their IntegrationPointVariables
    if( within == ELEMENT || within == ELEMENT_INTEGRATION_POINT || within == SECTOR_INTEGRATION_POINT || within == FACET_INTEGRATION_POINT )
      return ElementIntegrationPointVariables();
    else if( within == FACE || within == FACE_INTEGRATION_POINT || within == FACE_SECTOR_INTEGRATION_POINT || within == FACE_FACET_INTEGRATION_POINT )
      return FaceIntegrationPointVariables();
    else if( within == INTER_FACE || within == INTER_FACE_INTEGRATION_POINT || within == INTER_FACE_SECTOR_INTEGRATION_POINT || within == INTER_FACE_FACET_INTEGRATION_POINT )
      return InterFaceIntegrationPointVariables();
    /// for all others zero
    return IntegrationPointVariables();
  }


template<uint32_t dim>
IntegrationPointVariables PropertyDatabase<dim>::ElementIntegrationPointVariables() const
  {
    uint32_t scalarsSI(0), vectorsSI(0), tensorsSI(0), arrayCountSI(0), arrayLengthSI(0),flaggedArrayCountSI(0),flaggedArrayLengthSI(0);
    VariableCount( ELEMENT_INTEGRATION_POINT, scalarsSI, vectorsSI, tensorsSI, arrayCountSI, arrayLengthSI, flaggedArrayCountSI, flaggedArrayLengthSI );
    const uint32_t totalDataDepthSI( scalarsSI + vectorsSI*dim + tensorsSI*dim*dim + arrayLengthSI + flaggedArrayLengthSI );
    const uint32_t totalFlagDepthSI( scalarsSI + vectorsSI*vectorFlags + tensorsSI*tensorFlags + arrayCountSI + flaggedArrayLengthSI );
    uint32_t scalarsSE(0), vectorsSE(0), tensorsSE(0), arrayCountSE(0), arrayLengthSE(0),flaggedArrayCountSE(0),flaggedArrayLengthSE(0);
    VariableCount( SECTOR_INTEGRATION_POINT, scalarsSE, vectorsSE, tensorsSE, arrayCountSE, arrayLengthSE, flaggedArrayCountSE, flaggedArrayLengthSE );
    const uint32_t totalDataDepthSE( scalarsSE + vectorsSE*dim + tensorsSE*dim*dim + arrayLengthSE + flaggedArrayLengthSE);
    const uint32_t totalFlagDepthSE( scalarsSE + vectorsSE*vectorFlags + tensorsSE*tensorFlags + arrayCountSE + flaggedArrayLengthSE );
    uint32_t scalarsFA(0), vectorsFA(0), tensorsFA(0), arrayCountFA(0), arrayLengthFA(0),flaggedArrayCountFA(0),flaggedArrayLengthFA(0);
    VariableCount( FACET_INTEGRATION_POINT, scalarsFA, vectorsFA, tensorsFA, arrayCountFA, arrayLengthFA, flaggedArrayCountFA, flaggedArrayLengthFA );
    const uint32_t totalDataDepthFA( scalarsFA + vectorsFA*dim + tensorsFA*dim*dim + arrayLengthFA + flaggedArrayLengthFA );
    const uint32_t totalFlagDepthFA( scalarsFA + vectorsFA*vectorFlags + tensorsFA*tensorFlags + arrayCountFA + flaggedArrayLengthFA );
    
    return IntegrationPointVariables( scalarsSI, vectorsSI, tensorsSI, arrayCountSI, arrayLengthSI, flaggedArrayCountSI, flaggedArrayLengthSI,
                                      totalDataDepthSI, totalFlagDepthSI,
                                      scalarsSE, vectorsSE, tensorsSE, arrayCountSE, arrayLengthSE, flaggedArrayCountSE, flaggedArrayLengthSE,
                                      totalDataDepthSE, totalFlagDepthSE,
                                      scalarsFA, vectorsFA, tensorsFA, arrayCountFA, arrayLengthFA, flaggedArrayCountFA, flaggedArrayLengthFA,
                                      totalDataDepthFA, totalFlagDepthFA );
  }


template<uint32_t dim>
IntegrationPointVariables PropertyDatabase<dim>::FaceIntegrationPointVariables() const
  {
    uint32_t scalarsSI(0), vectorsSI(0), tensorsSI(0), arrayCountSI(0), arrayLengthSI(0),flaggedArrayCountSI(0),flaggedArrayLengthSI(0);
    VariableCount( FACE_INTEGRATION_POINT, scalarsSI, vectorsSI, tensorsSI, arrayCountSI, arrayLengthSI, flaggedArrayCountSI, flaggedArrayLengthSI );
    const uint32_t totalDataDepthSI( scalarsSI + vectorsSI*dim + tensorsSI*dim*dim + arrayLengthSI + flaggedArrayLengthSI );
    const uint32_t totalFlagDepthSI( scalarsSI + vectorsSI*vectorFlags + tensorsSI*tensorFlags + arrayCountSI + flaggedArrayLengthSI );
    uint32_t scalarsSE(0), vectorsSE(0), tensorsSE(0), arrayCountSE(0), arrayLengthSE(0),flaggedArrayCountSE(0),flaggedArrayLengthSE(0);
    VariableCount( FACE_SECTOR_INTEGRATION_POINT, scalarsSE, vectorsSE, tensorsSE, arrayCountSE, arrayLengthSE, flaggedArrayCountSE, flaggedArrayLengthSE );
    const uint32_t totalDataDepthSE( scalarsSE + vectorsSE*dim + tensorsSE*dim*dim + arrayLengthSE + flaggedArrayLengthSE);
    const uint32_t totalFlagDepthSE( scalarsSE + vectorsSE*vectorFlags + tensorsSE*tensorFlags + arrayCountSE + flaggedArrayLengthSE);
    uint32_t scalarsFA(0), vectorsFA(0), tensorsFA(0), arrayCountFA(0), arrayLengthFA(0),flaggedArrayCountFA(0),flaggedArrayLengthFA(0);
    VariableCount( FACE_FACET_INTEGRATION_POINT, scalarsFA, vectorsFA, tensorsFA, arrayCountFA, arrayLengthFA, flaggedArrayCountFA, flaggedArrayLengthFA );
    const uint32_t totalDataDepthFA( scalarsFA + vectorsFA*dim + tensorsFA*dim*dim + arrayLengthFA + flaggedArrayLengthFA );
    const uint32_t totalFlagDepthFA( scalarsFA + vectorsFA*vectorFlags + tensorsFA*tensorFlags + arrayCountFA + flaggedArrayLengthFA );

    return IntegrationPointVariables( scalarsSI, vectorsSI, tensorsSI, arrayCountSI, arrayLengthSI, flaggedArrayCountSI, flaggedArrayLengthSI,
                                      totalDataDepthSI, totalFlagDepthSI,
                                      scalarsSE, vectorsSE, tensorsSE, arrayCountSE, arrayLengthSE, flaggedArrayCountSE, flaggedArrayLengthSE,
                                      totalDataDepthSE, totalFlagDepthSE,
                                      scalarsFA, vectorsFA, tensorsFA, arrayCountFA, arrayLengthFA, flaggedArrayCountFA, flaggedArrayLengthFA,
                                      totalDataDepthFA, totalFlagDepthFA );
  }

template<uint32_t dim>
IntegrationPointVariables PropertyDatabase<dim>::InterFaceIntegrationPointVariables() const
  {
    uint32_t scalarsSI(0), vectorsSI(0), tensorsSI(0), arrayCountSI(0), arrayLengthSI(0),flaggedArrayCountSI(0),flaggedArrayLengthSI(0);
    VariableCount( INTER_FACE_INTEGRATION_POINT, scalarsSI, vectorsSI, tensorsSI, arrayCountSI, arrayLengthSI, flaggedArrayCountSI, flaggedArrayLengthSI );
    const uint32_t totalDataDepthSI( scalarsSI + vectorsSI*dim + tensorsSI*dim*dim + arrayLengthSI +flaggedArrayLengthSI );
    const uint32_t totalFlagDepthSI( scalarsSI + vectorsSI*vectorFlags + tensorsSI*tensorFlags + arrayCountSI + flaggedArrayLengthSI );
    uint32_t scalarsSE(0), vectorsSE(0), tensorsSE(0), arrayCountSE(0), arrayLengthSE(0),flaggedArrayCountSE(0),flaggedArrayLengthSE(0);
    VariableCount( INTER_FACE_SECTOR_INTEGRATION_POINT, scalarsSE, vectorsSE, tensorsSE, arrayCountSE, arrayLengthSE, flaggedArrayCountSE, flaggedArrayLengthSE );
    const uint32_t totalDataDepthSE( scalarsSE + vectorsSE*dim + tensorsSE*dim*dim + arrayLengthSE + flaggedArrayLengthSE);
    const uint32_t totalFlagDepthSE( scalarsSE + vectorsSE*vectorFlags + tensorsSE*tensorFlags + arrayCountSE + flaggedArrayLengthSE );
    uint32_t scalarsFA(0), vectorsFA(0), tensorsFA(0), arrayCountFA(0), arrayLengthFA(0),flaggedArrayCountFA(0),flaggedArrayLengthFA(0);
    VariableCount( INTER_FACE_FACET_INTEGRATION_POINT, scalarsFA, vectorsFA, tensorsFA, arrayCountFA, arrayLengthFA, flaggedArrayCountFA, flaggedArrayLengthFA );
    const uint32_t totalDataDepthFA( scalarsFA + vectorsFA*dim + tensorsFA*dim*dim + arrayLengthFA + flaggedArrayLengthFA );
    const uint32_t totalFlagDepthFA( scalarsFA + vectorsFA*vectorFlags + tensorsFA*tensorFlags + arrayCountFA + flaggedArrayLengthFA );

    return IntegrationPointVariables( scalarsSI, vectorsSI, tensorsSI, arrayCountSI, arrayLengthSI, flaggedArrayCountSI, flaggedArrayLengthSI,
                                      totalDataDepthSI, totalFlagDepthSI,
                                      scalarsSE, vectorsSE, tensorsSE, arrayCountSE, arrayLengthSE, flaggedArrayCountSE, flaggedArrayLengthSE,
                                      totalDataDepthSE, totalFlagDepthSE,
                                      scalarsFA, vectorsFA, tensorsFA, arrayCountFA, arrayLengthFA, flaggedArrayCountFA, flaggedArrayLengthFA,
                                      totalDataDepthFA, totalFlagDepthFA );
  }


template<uint32_t dim>
void PropertyDatabase<dim>::VariableCount( PLACEMENT within, uint32_t& scalars, uint32_t& vectors, uint32_t& tensors,
                                           uint32_t& arrayCount, uint32_t& arrayLength,
                                           uint32_t& flaggedArrayCount, uint32_t& flaggedArrayLength ) const
  {
    scalars             = VariableCount     ( within, SCALAR );
    vectors             = VariableCount     ( within, VECTOR );
    tensors             = VariableCount     ( within, TENSOR );
    arrayCount          = VariableCount     ( within, ARRAY  );
    arrayLength         = ArrayLengthTotal  ( within );
    flaggedArrayCount   = VariableCount     ( within, FLAGGEDARRAY  );
    flaggedArrayLength  = FlaggedArrayLengthTotal   ( within );
  }



template<uint32_t dim>
uint32_t PropertyDatabase<dim>::VariableCount() const
  {
    uint32_t variableCount(0);
    for( auto it( variableCount_.begin() );  it != variableCount_.end(); ++it )
      for( auto iit( it->second.begin() );  iit != it->second.end(); ++iit )
        variableCount += iit->second;
        
    return variableCount;
  }



template<uint32_t dim>
uint32_t PropertyDatabase<dim>::VariableCount( PLACEMENT variablePlacement, VARIABLE_TYPE variableType ) const
  {
    uint32_t variableCount(0);
    for( auto it( variableCount_.begin() );  it != variableCount_.end(); ++it )
      if( it->first == variablePlacement )
        for( auto iit( it->second.begin() ); iit != it->second.end(); ++iit )
          if( iit->first == variableType )
            variableCount += iit->second;
            
    return variableCount;
  }


template<uint32_t dim>
uint32_t PropertyDatabase<dim>::VariableCount( VARIABLE_TYPE variableType ) const
  {
    uint32_t variableCount(0);
    for( auto it( variableCount_.begin() ); it != variableCount_.end(); ++it )
         for( auto iit( it->second.begin() ); iit != it->second.end(); ++iit )
           if( iit->first == variableType )
             variableCount += iit->second;
             
    return variableCount;
  }


template<uint32_t dim>
uint32_t PropertyDatabase<dim>::VariableCount( PLACEMENT variablePlacement ) const
  {
    uint32_t variableCount(0);
    for( auto it( variableCount_.begin() ); it != variableCount_.end(); ++it )
      if( it->first == variablePlacement )
        for( auto iit( it->second.begin() ); iit != it->second.end(); ++iit )
          variableCount += iit->second;
          
    return variableCount;
  }


/// Clears propertyPlacements and inserts all PLACEMENT for all propertyNames
template<uint32_t dim>
void PropertyDatabase<dim>::ListVariables( const set<string>& propertyNames, set<PLACEMENT>& propertyPlacments ) const
{
  propertyPlacments.clear();
  for( set<string>::const_iterator it( propertyNames.begin() ); it != propertyNames.end(); ++it )
    {
      csmp::Index propKey( StorageKey( it->c_str() ) );
      propertyPlacments.insert( propKey.place );
    }
}



template class PropertyDatabase<1U>;
template class PropertyDatabase<2U>;
template class PropertyDatabase<3U>;


// non-member functions

/// counts variables in PropertyDatabase that are placed on finite volumes
template<uint32_t dim>
uint32_t finiteVolumeVariables( const PropertyDatabase<dim>& database )
 {
    uint32_t FV_variables = database.VariableCount( SECTOR_INTEGRATION_POINT );
    FV_variables    += database.VariableCount( FACET_INTEGRATION_POINT );
                      
    FV_variables    += database.VariableCount( FACE_SECTOR_INTEGRATION_POINT );
    FV_variables    += database.VariableCount( FACE_FACET_INTEGRATION_POINT );

    FV_variables    += database.VariableCount( INTER_FACE_SECTOR_INTEGRATION_POINT );
    FV_variables    += database.VariableCount( INTER_FACE_FACET_INTEGRATION_POINT );

    return FV_variables;
 }

template uint32_t finiteVolumeVariables( const PropertyDatabase<3>& );
template uint32_t finiteVolumeVariables( const PropertyDatabase<2>& );
template uint32_t finiteVolumeVariables( const PropertyDatabase<1>& );


} // end namespace csmp


