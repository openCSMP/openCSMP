//
//  PropertyData_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 7/04/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include <limits>
#include "PropertyData_Test.hpp"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

using namespace std;

namespace csmp {

void PropertyData_Test::run()
 {
    // TESTING THE PUBLIC INTERFACE OF CLASS
   
    // SCALAR VARIABLES
    // ----------------
    // constructor for all possible csmp variable types
    PropertyData dataset1( ELEMENT, SCALAR, 3U ), dataset3( ELEMENT, SCALAR, 3U );
   
    // stick in 3 scalars
    ScalarVariable sc1(PLAIN,2.), sc2(DIRICH,5.), sc3;
    dataset1.Reserve( 5, 5 );
    pushBack( dataset1, sc1 );
    pushBack( dataset1, sc2 );
    pushBack( dataset1, sc3 );
    // a copy
    PropertyData dataset2( dataset1 );
    // comparison
    dataset1.Out();
    dataset2.Out();
    // are they the same? (operator==) - should be fd
    _test( (dataset1 == dataset2) == false );
    // and a third one (asignment operator)
    dataset3 = dataset2;
    _test( dataset3.Size() == 3 );
    // should be empty
    dataset3.Clear();
    _test( dataset3.Size() == 0 );
   
    // resizing
    dataset1.Resize( 5, 5 );
    cout <<"\n\nPropertyData_Test::run: dataset1 after resizing to 5:";
    dataset1.Out();
    _test( dataset1.Size() == 5 );
 
    /// loosing entry 2 and 4, 5
    map<size_t,size_t> o_n_elmt_ids;
    o_n_elmt_ids.insert( make_pair(0,0) ); // keep first value: PLAIN:2
    o_n_elmt_ids.insert( make_pair(2,1) ); // keep third value: ANY:NaN
    dataset1.ReduceTo( o_n_elmt_ids );
    dataset1.Out();
    _test( dataset1.Size() == 2 );
 
    // accessors / mutators for values (operator[] cannot be overloaded for this)
    /// scalars
    _test( dataset1.Flag(0) == PLAIN );
    _test( dataset1.Flag(1) == ANY );
/*
    /// vectors and array variables
    VARIABLE_FLAG& Flag( size_t nth_value, size_t ith_dim );
    VARIABLE_FLAG  Flag( size_t nth_value, size_t ith_dim ) const;
    /// tensors
    VARIABLE_FLAG& Flag( size_t nth_value, size_t ith_row, size_t jth_col );
    VARIABLE_FLAG  Flag( size_t nth_value, size_t ith_row, size_t jth_col ) const;
 
    /// scalars
    double64&      Value( size_t nth_value );
    double64       Value( size_t nth_value ) const;
    /// vectors and array variables
    double64&      Value( size_t nth_value, size_t ith_dim );
    double64       Value( size_t nth_value, size_t ith_dim ) const;
    /// tensors
    double64&      Value( size_t nth_value, size_t ith_row, size_t jth_col );
    double64       Value( size_t nth_value, size_t ith_row, size_t jth_col ) const;
*/
    /// checks range; all values included
    pushBack( dataset1, sc1 );
    pushBack( dataset1, sc2 );
    // pushBack( dataset1, sc3 ); // NaN
    double64 tmin, tmax;
    dataset1.MinMaxOf( tmin, tmax );
    _equal( tmin, 2., 1e-10 );
    _equal( tmax, 5., 1e-10 );
  
    /// scales the current variable range to the new one; all values are scaled
    dataset1.ScaleRangeTo( 1., 20. );
    dataset1.MinMaxOf( tmin, tmax );
    _equal( tmin,  1., 1e-10 );
    _equal( tmax, 20., 1e-10 );

    /// add supplied value to all entries (look at conventions for vectors and tensors in the doc of these classes)
    dataset1.OffsetRangeBy( 15. );
    dataset1.MinMaxOf( tmin, tmax );
    _equal( tmin, 16., 1e-10 );
    _equal( tmax, 35., 1e-10 );
  
    cout <<"\n\nPropertyData_Test::run: scaled and offset dataset1:";
    dataset1.Out();
    cout <<"\n\nPropertyData_Test::run: sqrt of values ibn dataset1:";
    dataset1.TransformValues( std::sqrt );
    // eliminating the nan value
    dataset1.Value(1) = -1.;
    dataset1.Out();
   
    /// writing stored flag and data values to file
    cout <<"\n\nPropertyData_Test::run: dataset1 written to file and read back to memory:";
    PropertyData datasetN( dataset1 ); // backup copy

    std::FILE* out_fp = fopen( "PropertyData_Test", "wb" );
    dataset1.OutBinary( out_fp );
    fclose( out_fp );
  
    /// reading stored flag and data vaues from file
    std::FILE* in_fp = fopen( "PropertyData_Test", "rb" );
    datasetN = inBinaryPropertyData( in_fp );
    fclose( in_fp );
    datasetN.Out();
    _test( datasetN == dataset1 );
   
 } // end run


} // end csmp
