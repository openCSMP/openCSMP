/*! \file PropertyDatabase.h */

#ifndef CSMP_PROPERTY_DATABASE_H
#define CSMP_PROPERTY_DATABASE_H

#include <iostream>
#include <map>
#include <string>
#include <list>
#include <vector>
#include "Parameter.h"
#include "Index.h"
#include "IndexTracker.h"

namespace csmp {

                                        
/**
 
@brief Specifications of the variables used by the current model, i.e.
name, notation, SI-unit, placement, type etc.

@author S. K. Matthaei
@author P. Lang
@date 1994-2012

@section motivation Motivation
 
To manage access to physical variables in the numerical simulation.
Depending on their purpose such variables have a different type
such as scalar, vector, tensor or array.
In addition, they have a different placement in the mesh:
model, region, element node etc.
Thus, housekeeping is required to keep track of such different variables.
Firstly, one must be able to search for them by name. 
Secondly, their values have to be checked against a physically
meaningful range. Thirdly, one needs to know how a variable was 
defined (as a scalar, vector or tensor) and at which locations in a
discretization it shall be computed (nodes, integration points or
elements etc.). To accomplish this
there must be a unique access key for each variable 
in order to retrueve it in a
computationally efficient way. 
 
 
@section design Design Intent

The design intent of the PropertyDatabase is to have 
a single object in the CSMP code which manages the physical 
variables of a simulation and can be queried for their type,
placement, and property index. These three characteristics
combine into the structure csmp::Index
which serves the role of an access keys to distributed
variable values. 

The records in the PropertyDatabase are Parameter class instances
which contain the variable name, its csmp::Index, and a physically 
meaningful range for each variable. This range can be queried to test 
whether implausible results were computed. 

Finally, the database must ascertain uniqueness of the variable 
'string' identifiers and allow to retrieve these in a reverse approach,
starting from the inside of a computation in which only the efficient 
access keys are known. This capability is needed for error reporting. 

The creation or deletion of variables at runtime is a further task
of the PropertyDatabase.
 
 
@section structure Structure

The PropertyDatabase represents a binary tree of fixed-size
physical-variable records which are accessed via name strings. 
Access of variable specifications benefits from the tree structure, but
the reverse search for a variable name using the csmp::Index as a key
has a linear time complexity. Due to the typically small number of 
physical variables in a model, this is not important. 
 
 
@section participants Participants
 
The PropertyDatabase uses Parameter class objects as variable records.
These are accessed by the variable names which are stored in std::string
objects. 
 
 
@section collaborations Collaborations

The PropertyDatabase is a data member of the Model object
and collaborates with Index and Parameter records.

@attention The PropertyDatabase does not store any variable values!
These are stored in the LocalVariableStorage associated
with the basic tree objects: Node, Element, etc.

 
@section consequences Consequences

As is described in the CSMP User's Guide, properties of a physical 
variable and csmp::Index structures can be obtained from the 
PropertyDatabase. As additional functionality, the PropertyDatabase
can be used for unit conversions and to add or delete physical
variables at runtime.
 
 
@section implementation Implementation 

The PropertyDatabase contains an STL map of variable names serving
as keys to Parameter objects which hold the variable specifications. 
The variable names are stored in the map in alphabetical order and
the map is initialized from a text file when the PropertyDatabase is
constructed. After all variable records have been added to the map, they are
counted and csmp::Index.index integers are issued. These indices are used 
to access the right LocalVariableStorage object. 
 
 
@section examples Application Examples

An example of a variable entry in the text input file which is read by
the constructor to initialize the PropertyDatabase, is shown below: 
 
@code
name  parameter  unit  index  min.  max.  place  usage  explanation  reference
									
fluid pressure PF Pa 1 0.00E+00	1.00E+10 node	absolute value	fluid pressure at node	alpha version doc
@endcode
 
In this example, the string 'fluid pressure' defines the variable name which is
used in the top level interfaces of the Model, the Algorithm, Interrelation
or Visitor objects. The second variable 'PF' is an abbreviation which relates this
variable name to the notation which the modeler uses in his/her equations.
'index' defines the variable type: '1' for a scalar, '2' for a vector,
and '3' for a tensor variable. The following two floating point values in
scientific notation, define a physically meaningful range for the variable, 
and the third parameter denotes the variable placement in the finite-element 
mesh. 
 
The second examples illustrates how to get a csmp::Index from a variable name: 
 
@code
// get a reference to the property database from the Model
const PropertyDatabase<dim>&  p_ref = moab.Database();
// get the access key for the variable "fluid pressure"
csmp::Index property_key = p_ref.StorageKey("fluid pressure");
@endcode
 
This code retrieves the csmp::Index from the variable database. Typically,
one retrieves the property key only once if one needs to access 
variable values. Once the key is established it can be used instead
of the property name for fast property access. 

@section arrayVariables Array Variables

Worth mentioning is the variable file syntax for ArrayVariable handling. Array variables
are created by assigning a number greater than 3 where the variable type is specified. The number
indicates the size of the array, which as of now (9/3/2012) may not be changed at runtime.
@code
element array				EA	  X	  22	0.	1.	element
@endcode
here we instantiated an array on the element of size 22. So to recap the variable type options:
  1 .. ScalarVariable
  2 .. VectorVariable
  3 .. TensorVariable
  4 to max<size_t> ArrayVariable of that size

@todo (2-C) DocMe (Update to P. Lang 2012)
@todo (3-D) Refactor!! maps

*/
template<size_t dim>
class PropertyDatabase  {
 public:
   PropertyDatabase();
   explicit PropertyDatabase( const char* variablesFileName, bool isBinary = false );
   ~PropertyDatabase();
  
   PropertyDatabase<dim>& operator=( const PropertyDatabase<dim>& p );
   void Initialize( const char* variables_file );
   void DeepCopy( const PropertyDatabase<dim>& );

   bool BinaryOut( const char* fileName ) const;
   bool BinaryOut( FILE* fp ) const;
   bool BinaryIn( const char* fileName );
   bool BinaryIn( FILE* fp );

   const char*    VariablesFile() const;
   csmp::Index    StorageKey( const char* property_name ) const;
   csmp::Parameter  Parameter( const char* property_name ) const;
   PLACEMENT      Placement( const char* property_name ) const; 
   size_t         Index( const char* property_name ) const;
   VARIABLE_TYPE  Type( const char* property_name ) const;
   size_t         Components( const char* property_name ) const;
   const char*    Name( const csmp::Index& idx ) const;
   const char*    Usage( const char* property_name ) const;
   bool           IsDefined( const char* property_name ) const;
   bool           IsDefined( const csmp::Index& idx ) const;
   void           RangeOf( const char* property_name, double64& min, double64& max ) const;
   void           CheckRange( const char* property_name, double64& var ) const;
   const char*    Unit( const char* property_name ) const;
  
   LocalVariables             LocalVariablesAt( PLACEMENT within ) const;
   IntegrationPointVariables  IntegrationPointVariablesAt( PLACEMENT within ) const;
  
   void                       VariableCountAt( PLACEMENT within, size_t& scalars,
                                               size_t& vectors, size_t& tensors,
                                               size_t& arrayCount, size_t& arrayLength ) const;
   
   double64       UnitConversionFactor( const char* current_system, 
                                        const char* desired_system, 
                                        const char* unit ) const;

   /// adding a property at runtime from the command line (NB: you should call method of model to do this)
   csmp::Index    AddProperty();    

   csmp::Index    AddProperty(const char* property_name, const char* unit, size_t last_max_index,
                               VARIABLE_TYPE, PLACEMENT, size_t vsize = 1,
                               double64 vmin=-1.0e+30 , double64 vmax=1.0e+30 , std::string usage="???");

   csmp::Index    AddProperty( const char* property_name, const char* unit,
                               VARIABLE_TYPE, PLACEMENT, size_t vsize = 1,
                               double64 vmin=-1.0e+30 , double64 vmax=1.0e+30, std::string usage="???" );
                                   
   void           DeleteProperty( const char* property_name ); 
   
   std::map<std::string,csmp::Parameter>::const_iterator  Begin() const;
   std::map<std::string,csmp::Parameter>::const_iterator  End() const;

   std::map<PLACEMENT,std::map<VARIABLE_TYPE,size_t> >::const_iterator  VariableCountBegin() const; 
   std::map<PLACEMENT,std::map<VARIABLE_TYPE,size_t> >::const_iterator  VariableCountEnd() const; 

   void   FlushToScreen() const;
   void   FlushToScreen( const char* propname ) const;
   void   ListVariables() const;
   void   ListVariableNames() const;
   void   ListKeys( std::list<csmp::Index>& keys ) const;
   void   ListProperties( std::map<std::string,csmp::Index>& props ) const;
   void   ListProperties( const std::set<std::string>&, std::set<PLACEMENT>& props ) const;
   size_t ListProperties( PLACEMENT place, std::map<std::string,csmp::Index>& props ) const;
   size_t ListProperties( PLACEMENT place, std::set<std::string>& props ) const;
   size_t ListProperties( PLACEMENT, VARIABLE_TYPE, std::set<std::string>& props ) const;

   size_t VariableCount( PLACEMENT )                  const;
   size_t VariableCount( VARIABLE_TYPE )              const;
   size_t VariableCount( PLACEMENT, VARIABLE_TYPE )   const;
   size_t VariableCount()                             const;

   size_t ArrayLengthTotal( PLACEMENT place ) const;
   void   ArrayLengths( PLACEMENT place, std::vector<size_t>& arrayLengths ) const;

   size_t FlaggedArrayLengthTotal( PLACEMENT place ) const;
   void   FlaggedArrayLengths( PLACEMENT place, std::vector<size_t>& arrayLengths ) const;
   
   bool   WriteVariablesFile( const char* fileName ) const;

   void   Out() const { Out(std::cout); }
   void   Out(std::ostream& out) const;

   void   Verbose(bool verbose) { this->verbose_=verbose; }
   bool   Verbose() { return this->verbose_; }

 private:
   const size_t vectorFlags, tensorFlags;
   bool verbose_;
   std::string  physvarsFile;
   std::map<PLACEMENT,std::map<VARIABLE_TYPE,size_t> > variableCount_; ///< all placements and variable types in here
   std::map<std::string,csmp::Parameter>  propList_; ///< all parameters (and with those the indices)
   IndexTracker indexTracker_; ///< used to keep track of all index references and update indices after runtime changes

   PropertyDatabase( const PropertyDatabase<dim>& );
   void   InitializeCount();
   void   InitializeVariableTypeCount( std::map<VARIABLE_TYPE,size_t>& );
   void   UpdateParametersAndDatabase();
   void   CountVariables(); 
   void   AttachIndices();
   void   UpdateIndexReferences();
   void   DetachIndices( std::string parameterName );
   void   AssignVariableIndices();
   void   UpdateIndicesAfterDelete();
   void   TextToBinaryFile( const char* text_database_file, bool echo_to_screen=false );
   void   EstablishIndexLocalAndIntegrationPointVariables();
   void   EstablishIndexOffsets();
   void   EstablishScalarOffsets        ( PLACEMENT whithin );
   void   EstablishVectorOffsets        ( PLACEMENT whithin );
   void   EstablishTensorOffsets        ( PLACEMENT whithin );
   void   EstablishArrayOffsets         ( PLACEMENT whithin );
   void   EstablishFlaggedArrayOffsets  ( PLACEMENT whithin );
   void   EstablishVariableTypeDependentProperties( int vtype, csmp::Index& key ) const;
   void   EstablishVariableTypeDependentProperties( int vtype, size_t size, csmp::Index& key ) const;
   void   EstablishVariableTypeDependentProperties( std::string type, csmp::Index& key ) const;
   void   EstablishPlacementDependentProperties( PLACEMENT placement, csmp::Index& key ) const;
  
   void   VariableCount( PLACEMENT within, size_t& scalars, size_t& vectors, size_t& tensors,
                         size_t& arrayCount, size_t& arrayLength, size_t& flaggedArrayCount,
                         size_t& flaggedArrayLength ) const;
  
   IntegrationPointVariables  ElementIntegrationPointVariables() const;
   IntegrationPointVariables  FaceIntegrationPointVariables() const;
   IntegrationPointVariables  InterFaceIntegrationPointVariables() const;
};

<<<<<<< HEAD
=======



// property iterators
template<size_t dim>
inline std::map<std::string,Parameter>::const_iterator  PropertyDatabase<dim>::Begin() const
 {
    return propList_.begin();
 } 
 
 
template<size_t dim>
inline std::map<std::string,Parameter>::const_iterator  PropertyDatabase<dim>::End() const
 {
    return propList_.end();
 } 


template<size_t dim>
std::map<PLACEMENT,std::map<VARIABLE_TYPE,size_t> >::const_iterator PropertyDatabase<dim>::VariableCountEnd() const
  {
    return variableCount_.end();
  }


template<size_t dim>
std::map<PLACEMENT,std::map<VARIABLE_TYPE,size_t> >::const_iterator PropertyDatabase<dim>::VariableCountBegin() const
  {
    return variableCount_.begin();
  }

/**
 
Tells the user whether a variable exists in the property database. 

@param s The name or csmp::Index of the variable which shall be searched for in the
database. 

@return The boolean variable 'true' or 'false'.
*/
template<size_t dim>
inline bool PropertyDatabase<dim>::IsDefined( const char* s ) const
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
template<size_t dim>
inline PLACEMENT PropertyDatabase<dim>::Placement( const char* s ) const 
  {
     std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
     
     if ( iter != propList_.end() ) 
       return (*iter).second.key.place;
     else 
       throw csmp::Exception( CSMP_ERROR, "PropertyDatabase<dim>::Placement", "Unable to identify property", s );

     // shouldn't get here
     return MODEL;        
  } // end WhereIs   



/**
 
Returns the csmp::Index.index of the target physical variable. 

@param s The name of the physical variable whose index shall be determined.

@return The variable index which defines the number of the corresponding property
vector inside of the MemoryManager. 

@section messages Messages

If the variable is not defined, an error will be reported. 
*/
template<size_t dim>
inline size_t PropertyDatabase<dim>::Index( const char* s ) const 
  {
     std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
     
     if ( iter != propList_.end() ) 
       return (*iter).second.key.index;
     else 
       std::cout <<"\nPropertyDatabase::Index: Unable to identify index of: " << s << std::endl;

     return ULONG_MAX;  
  } // end Index 



/**
 
Reports whether the target physical variable is of scalar, vector or
tensor type. 

@param s The name of the physical variable whose type shall be established.

@return The type of the physical variable.


@section messages Messages

If the variable is not defined, an error will be reported. 
*/
template<size_t dim>
inline VARIABLE_TYPE  PropertyDatabase<dim>::Type( const char* s ) const 
  {
     std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
     
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
template<size_t dim>
inline size_t  PropertyDatabase<dim>::Components( const char* s ) const 
  {
     std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
     
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
template<size_t dim>
inline csmp::Index  PropertyDatabase<dim>::StorageKey( const char* s ) const
 {
    std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
    if ( iter == propList_.end() ) {
         std::string message("Variable '");
         message += s;
         message +="' is undefined";
         throw csmp::Exception( CSMP_ERROR, "PropertyDatabase<dim>::StorageKey:", message );
      }
    if ( iter == propList_.end() ) 
      return csmp::Index();

    return (*iter).second.key;
 }

template<size_t dim>
inline csmp::Parameter  PropertyDatabase<dim>::Parameter( const char* s ) const
{
  std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
  if ( iter == propList_.end() ) {
    std::string message("Variable '");
    message += s;
    message +="' is undefined";
    throw csmp::Exception( CSMP_ERROR, "PropertyDatabase<dim>::Parameter:", message );
  }
  if ( iter == propList_.end() ) 
    return csmp::Parameter();

  return (*iter).second;
}



/// reports how the variable is normally used (as input by user vs. computed)
template<size_t dim>
inline const char*  PropertyDatabase<dim>::Usage( const char* s ) const
 {
    std::map<std::string,csmp::Parameter>::const_iterator  iter(propList_.find(std::string(s)));
    if ( iter == propList_.end() ) {
         std::string message("Variable '");
         message += s;
         message +="' is undefined";
         throw csmp::Exception( CSMP_ERROR, "PropertyDatabase<dim>::Usage:", message );
      }
    if ( iter == propList_.end() ) return "undefined";

    return (*iter).second.usage.c_str();
 }




>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
} // csmp

#endif


