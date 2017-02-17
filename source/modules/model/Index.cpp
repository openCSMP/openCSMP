// ------------------------------------------------------------------
// Index.h
// Created: Mon May 11 15:28:50 MDT 1998
// (c) Copyright 1998 by Dr. Stephan Matthai
// ------------------------------------------------------------------
#include "Index.h"
#include <iomanip>
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

using namespace std;

namespace csmp {

/** comparator less

Provides required functionality for use in an STL container.

@attention Do not use this as an example for how to code
since it deduces information at runtime which should already be
known at compile time.
@param i a csmp::Index
@return true if smaller than argument

*/
bool Index::operator<( const csmp::Index& i ) const
 {
    const size_t lhs( (type+1) * (place+1) * (index+1U) * dataDepth * flagDepth  
                      * (offsetFactorSimplex+1) * (offsetFactorSector+1) *(flagOffset+1) );
   
    const size_t rhs( (i.type+1) * (i.place+1) * (i.index+1U) * i.dataDepth * i.flagDepth 
                      * (i.offsetFactorSimplex+1) * (i.offsetFactorSector+1) * (i.flagOffset+1));
   
    return (lhs < rhs);
 }


/**
@defgroup CSMPglobalFunctions Global CSMP Functions
Functions available globally in the csmp Namespace
 */


/**
@addtogroup CSMPglobalFunctions
@{
*/


/** Checks Variable for Type

@return the value of the corresponding csmp enumeration VARIABLE_TYPE.
This method uses RunTimeTypeIdentification with corresponding
implications for speed.

@todo achieve this with a compile time computation.

*/
template<typename Var>
VARIABLE_TYPE variableType( const Var& )
 {
    if ( typeid(Var) == typeid(ScalarVariable) )        return SCALAR;

    if ( typeid(Var) == typeid(VectorVariable<1U>) )    return VECTOR;
    if ( typeid(Var) == typeid(VectorVariable<2U>) )    return VECTOR;
    if ( typeid(Var) == typeid(VectorVariable<3U>) )    return VECTOR;

    if ( typeid(Var) == typeid(TensorVariable<1U>) )    return TENSOR;
    if ( typeid(Var) == typeid(TensorVariable<2U>) )    return TENSOR;
    if ( typeid(Var) == typeid(TensorVariable<3U>) )    return TENSOR;

    if ( typeid(Var) == typeid(ArrayVariable) )         return ARRAY;
    if ( typeid(Var) == typeid(FlaggedArrayVariable) )  return FLAGGEDARRAY;

    return static_cast<VARIABLE_TYPE>(UNSPECIFIED);
 }


/**
 @}
 */
// all these explicit instantiations
template VARIABLE_TYPE variableType( const ScalarVariable& );
template VARIABLE_TYPE variableType( const VectorVariable<1U>& );
template VARIABLE_TYPE variableType( const VectorVariable<2U>& );
template VARIABLE_TYPE variableType( const VectorVariable<3U>& );
template VARIABLE_TYPE variableType( const TensorVariable<1U>& );
template VARIABLE_TYPE variableType( const TensorVariable<2U>& );
template VARIABLE_TYPE variableType( const TensorVariable<3U>& );
template VARIABLE_TYPE variableType( const ArrayVariable& );
template VARIABLE_TYPE variableType( const FlaggedArrayVariable& );




/** Default constructor

  Assigns an Index of scalar type on the model.
  The index number itself will be set to the maximum
  and hence interpreted as not defined.
  Inline implementation.

*/
Index::Index() 
  : type(SCALAR), 
    place(UNDEFINED), 
    index(ULONG_MAX),
    dataDepth(0),
    flagDepth(0),
    dataOffset(0),
    flagOffset(0),
    offsetFactorSimplex(0),
    offsetFactorSector(0), 
    ipFactorSimplex(0),
    ipFactorSector(0),
    ipFactorFacet(0),
    localVariables(),
    integrationPointVariables(),
    indexTracker(nullptr)
  {}


/** Custom constructor

Neither sets local nor integration point variables and has no option for offsets/depth ( i.e. = 0 )

Assigns an Index with defined type and placement.
Inline implementation.
@param ty the variable type
@param pl the placement
@param idx the index

*/
Index::Index( VARIABLE_TYPE ty, PLACEMENT pl, size_t idx ) 
  : type(ty), place(pl),
    index(idx), dataDepth(0), flagDepth(0), dataOffset(0), 
    flagOffset(0), offsetFactorSimplex(0), offsetFactorSector(0), 
    ipFactorSimplex(0), ipFactorSector(0), ipFactorFacet(0),
    localVariables(), integrationPointVariables(), indexTracker(nullptr)
  {}


/** Custom constructor

This should be goto ctor, even though members may be initialized individually through struct
Assigns an Index with defined type and placement.
Inline implementation.

@param ty the variable type
@param pl the placement
@param idx the index
@param datadepth component count, size of variable (scalar 1, vector 2, tensor 3, array #components)
@param dataoffset Offset in of data in variable storage container (implicit dependence on LocalVariableStorage)
@param flagoffset Offset in of flag data in variable storage container (implicit dependence on LocalVariableStorage)
@param lvs LocalVariables of placement pl
@param ivs IntegrationPointVariables of elements
*/
Index::Index( VARIABLE_TYPE ty, PLACEMENT pl, size_t idx, 
                      size_t datadepth, size_t flagdepth, size_t dataoffset, size_t flagoffset,
                      const LocalVariables& lvs, const IntegrationPointVariables& ivs, 
                      size_t offsetfactorsimplex, size_t offsetfactorsector, size_t ipfactorsimplex, 
                      size_t ipfactorsector, size_t ipfactorfacet ) 
  : type(ty), place(pl),
    index(idx), dataDepth(datadepth), flagDepth(flagdepth), dataOffset(dataoffset), 
    flagOffset(flagoffset), offsetFactorSimplex(offsetfactorsimplex), offsetFactorSector(offsetfactorsector), 
    ipFactorSimplex(ipfactorsimplex), ipFactorSector(ipfactorsector), ipFactorFacet(ipfactorfacet),
    localVariables(lvs), integrationPointVariables(ivs), indexTracker(nullptr)
  {}

/** Copy constructor

Creates an index based on the provided parameter
@param idx an existing csmp::Index object
*/
Index::Index( const csmp::Index& idx ) 
  : type(idx.type), place(idx.place), index(idx.index),
    dataDepth(idx.dataDepth), flagDepth(idx.flagDepth), dataOffset(idx.dataOffset), flagOffset(idx.flagOffset),
    offsetFactorSimplex(idx.offsetFactorSimplex), offsetFactorSector(idx.offsetFactorSector), 
    ipFactorSimplex(idx.ipFactorSimplex), ipFactorSector(idx.ipFactorSector), ipFactorFacet(idx.ipFactorFacet),
    localVariables(idx.localVariables), integrationPointVariables(idx.integrationPointVariables),
    indexTracker(nullptr)
  {
  // attach itself to IndexTracker and vice versa
  if(idx.indexTracker)
    {
      Detach();
      idx.indexTracker->Attach( this, &idx ); 
      Attach( idx.indexTracker );
    }
  }


/// Assignment operator
Index&  Index::operator=( const csmp::Index& idx ) {
  if ( this != &idx ) {
    type                      = idx.type;
    place                     = idx.place;
    index                     = idx.index;
    dataDepth                 = idx.dataDepth;
    flagDepth                 = idx.flagDepth;
    dataOffset                = idx.dataOffset;
    flagOffset                = idx.flagOffset;
    offsetFactorSimplex       = idx.offsetFactorSimplex;
    offsetFactorSector        = idx.offsetFactorSector;
    ipFactorSimplex           = idx.ipFactorSimplex;
    ipFactorSector            = idx.ipFactorSector;
    ipFactorFacet             = idx.ipFactorFacet;
    localVariables            = idx.localVariables;
    integrationPointVariables = idx.integrationPointVariables;
    /// attach itself to IndexTracker and vice versa.
    if( idx.indexTracker )
      { 
        Detach();
        idx.indexTracker->Attach( this, &idx );
        Attach( idx.indexTracker );
      }
    }        

  return *this;
  }

/**
 @fn  void Index::UpdateData( const csmp::Index& idx )

 @brief Updates to the data described by idx without updating indexTracker

 @author  P. Lang
 @date  9/26/2012

 @param idx The index to be updated from
 */
void Index::UpdateData( const csmp::Index& idx )
  {
    if ( this != &idx ) {
        type                      = idx.type;
        place                     = idx.place;
        index                     = idx.index;
        dataDepth                 = idx.dataDepth;
        flagDepth                 = idx.flagDepth;
        dataOffset                = idx.dataOffset;
        flagOffset                = idx.flagOffset;
        offsetFactorSimplex       = idx.offsetFactorSimplex;
        offsetFactorSector        = idx.offsetFactorSector;
        ipFactorSimplex           = idx.ipFactorSimplex;
        ipFactorSector            = idx.ipFactorSector;
        ipFactorFacet             = idx.ipFactorFacet;
        localVariables            = idx.localVariables;
        integrationPointVariables = idx.integrationPointVariables;
        // MUST NOT BE UPSET: indexTracker = nullptr;
      }
  }




/** Boolean equals operator

@param i an csmp::Index object
@return boolean result of the comparison
@return true if equal
*/
bool Index::operator==( const csmp::Index& i ) const 
{
  return (  type==i.type && place==i.place && index==i.index && dataDepth==i.dataDepth && flagDepth==i.flagDepth && dataOffset==i.dataOffset 
            && flagOffset==i.flagOffset && offsetFactorSimplex==i.offsetFactorSimplex && offsetFactorSector==i.offsetFactorSector );
}



/** Boolean equals not operator

@return boolean result of the comparison
@return true if not equal
*/
bool Index::operator!=( const csmp::Index& i ) const 
{
  return ( type!=i.type || place!=i.place || index!=i.index || dataDepth!=i.dataDepth || flagDepth!=i.flagDepth || dataOffset!=i.dataOffset 
           || flagOffset!=i.flagOffset || offsetFactorSimplex!=i.offsetFactorSimplex || offsetFactorSector!=i.offsetFactorSector );
}


/** Boolean Function that checks if object has been defined

@return true if an index has been assigned
*/
bool Index::IsDefined() const { return (index != ULONG_MAX);  }


/// Registers an IndexTracker (does not detach from current!)
void Index::Attach( IndexTracker* newTracker ) { indexTracker = newTracker; }


/// DTor: Detaches intself from IndexTracker
Index::~Index()
  {
    if ( indexTracker != nullptr )
      indexTracker->Detach(this);
  }


/// Detaches itself from IndexTracker (i.e. called upon destruction)
void Index::Detach()
  {
    if ( indexTracker != nullptr ) indexTracker->Detach(this);
    indexTracker = nullptr;
  }




/** CSMP output feature

Out() returns all the internal data that Index needs to access 
variables. These include the type and placement of a variable in the 
finite element mesh and the associated data storage index for the CSMP
'MemoryManager'. 

*/
void Index::Out() const
{
  cout <<"\nIndex::Out: ";
  string str(parseType(type));
  cout <<"\n\tVARIABLE_TYPE: "<< str;
  str = parsePlacement(place);
  cout <<"\n\tPLACEMENT:        "<< str;
  cout <<"\n\tindex:            "<< index;
  cout <<"\n\t\tdata depth:       "<< dataDepth;
  cout <<"\n\t\tflag depth:       "<< flagDepth;
  cout <<"\n\t\tdata offset:      "<< dataOffset;
  cout <<"\n\t\tflag offset:      "<< flagOffset;
  cout <<"\n\t\toffset factor si: "<< offsetFactorSimplex;
  cout <<"\n\t\toffset factor se: "<< offsetFactorSector;
  cout <<"\n\t\tip factor FE:     "<< ipFactorSimplex <<" (number of integration points for finite element).";
  cout <<"\n\t\tip factor sector: "<< ipFactorSector <<" (number of integration points per finite volume sector).";
  cout <<"\n\t\tip factor facet:  "<< ipFactorFacet<<" (number of integration points per finite volume facet).";
  cout << endl;
}

bool Index::Out( FILE* fp ) const
{
fwrite( (void*) &type, sizeof(VARIABLE_TYPE), 1, fp );
fwrite( (void*) &place, sizeof(PLACEMENT), 1, fp );
fwrite( (void*) &index, sizeof(size_t), 1, fp );
fwrite( (void*) &dataDepth, sizeof(size_t), 1, fp );
fwrite( (void*) &flagDepth, sizeof(size_t), 1, fp );
fwrite( (void*) &dataOffset, sizeof(size_t), 1, fp );
fwrite( (void*) &flagOffset, sizeof(size_t), 1, fp );
fwrite( (void*) &offsetFactorSimplex, sizeof(size_t), 1, fp );
fwrite( (void*) &offsetFactorSector, sizeof(size_t), 1, fp );
fwrite( (void*) &ipFactorSimplex, sizeof(size_t), 1, fp );
fwrite( (void*) &ipFactorSector, sizeof(size_t), 1, fp );
fwrite( (void*) &ipFactorFacet, sizeof(size_t), 1, fp );
fwrite( (void*) &localVariables, sizeof(LocalVariables), 1, fp );
fwrite( (void*) &integrationPointVariables, sizeof(IntegrationPointVariables), 1, fp );
// we don't store pointers

return true; /// @todo (1-C) Meaningless return statement
}

bool Index::In( FILE* fp )
{
fread( (void*) &type, sizeof(VARIABLE_TYPE), 1, fp );
fread( (void*) &place, sizeof(PLACEMENT), 1, fp );
fread( (void*) &index, sizeof(size_t), 1, fp );
fread( (void*) &dataDepth, sizeof(size_t), 1, fp );
fread( (void*) &flagDepth, sizeof(size_t), 1, fp );
fread( (void*) &dataOffset, sizeof(size_t), 1, fp );
fread( (void*) &flagOffset, sizeof(size_t), 1, fp );
fread( (void*) &offsetFactorSimplex, sizeof(size_t), 1, fp );
fread( (void*) &offsetFactorSector, sizeof(size_t), 1, fp );
fread( (void*) &ipFactorSimplex, sizeof(size_t), 1, fp );
fread( (void*) &ipFactorSector, sizeof(size_t), 1, fp );
fread( (void*) &ipFactorFacet, sizeof(size_t), 1, fp );
fread( (void*) &localVariables, sizeof(LocalVariables), 1, fp );
fread( (void*) &integrationPointVariables, sizeof(IntegrationPointVariables), 1, fp );
indexTracker = nullptr; // we don't store pointers

return true; /// @todo (1-C) Meaningless return statement
}



/**
@defgroup CSMPglobalOperators Global CSMP Operators
Operators available globally in the csmp Namespace
 */


/**
@addtogroup CSMPglobalOperators
@{
*/

/** Custom global ostream operator

@attention Does not output IntegrationPointVariables / LocalVariables

Writes human-readable variable characteristics to given output stream
(use Out() to get all)
 
@param stream The output stream
@param o The csmp::Index object to output
@return Formatted output
*/
ostream& operator<<( ostream& stream, const csmp::Index& o )
 {
    string str(parseType(o.type));
    stream <<" VARIABLE_TYPE: "<< str;
    str  = parsePlacement(o.place);
    stream <<", PLACEMENT: "<< setw(16);
    stream << resetiosflags( ios::adjustfield );
    stream << setiosflags( ios::left ) << str;

#ifdef VARIABLE_STORAGE_DEBUG
    stream <<" index: "<< o.index <<" ";
    stream <<" data depth: "<< o.dataDepth <<" ";
    stream <<" flag depth: "<< o.flagDepth <<" ";
    stream <<" data offset: "<< o.dataOffset <<" ";
    stream <<" flag offset: "<< o.flagOffset <<" ";
    stream <<" offset factor si: "<< o.offsetFactorSimplex <<" ";
    stream <<" offset factor se: "<< o.offsetFactorSector <<" ";
    stream <<" ip factor si: "<< o.ipFactorSimplex <<" ";
    stream <<" ip factor se: "<< o.ipFactorSector <<" ";
    stream <<" ip factor fa: "<< o.ipFactorFacet <<" ";
#endif

    stream << resetiosflags( ios::adjustfield );
    return stream;
 }
 
 
 
 
/** Custom global istream operator

@attention Does not instantiate IntegrationPointVariables / LocalVariables

Provides formatted input to a given stream
@param stream The input stream
@param o The csmp::Index object to intput
@return Formatted intput
*/
istream& operator>>( istream& stream, csmp::Index& o )
 {
    string input;
    
    // type 
    stream >> input;
    o.type = parseType( input.c_str() );
    input.erase( input.begin(), input.end() );

    // placement 
    stream >> input;
    o.place = parsePlacement( input.c_str() );
    input.erase( input.begin(), input.end() );

    // index 
    stream >> o.index;

    // depth
    stream >> o.dataDepth;
    stream >> o.flagDepth;

    // offset
    stream >> o.dataOffset;
    stream >> o.flagOffset;
    stream >> o.offsetFactorSimplex;
    stream >> o.offsetFactorSector;
    stream >> o.ipFactorSimplex;
    stream >> o.ipFactorSector;
    stream >> o.ipFactorSimplex;

    return stream;     
 } // end 

/**
@}
*/


} // end namespace csmp
