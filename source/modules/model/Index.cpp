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
bool Index::operator<( const csmp::Index& i ) const // for stl
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
template<typename Var>  VARIABLE_TYPE variableType( const Var& )
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
    cout <<"\nVARIABLE_TYPE: "<< str;
    str = parsePlacement(place);
    cout <<"\nPLACEMENT:        "<< str;
    cout <<"\nindex:            "<< index << endl;
    cout <<"\ndata depth:       "<< dataDepth << endl;
    cout <<"\nflag depth:       "<< flagDepth << endl;
    cout <<"\ndata offset:      "<< dataOffset << endl;
    cout <<"\nflag offset:      "<< flagOffset << endl;
    cout <<"\noffset factor si: "<< offsetFactorSimplex << endl;
    cout <<"\noffset factor se: "<< offsetFactorSector << endl;
    cout <<"\nip factor si:     "<< ipFactorSimplex << endl;
    cout <<"\nip factor se:     "<< ipFactorSector << endl;
    cout <<"\nip factor fa:     "<< ipFactorFacet << endl;
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
  indexTracker = NULL; // we don't store pointers

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
