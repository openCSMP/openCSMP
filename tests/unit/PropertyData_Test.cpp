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

/// TESTING THE PUBLIC INTERFACE OF THE CLASS
void PropertyData_Test::run()
 {
    // ---------------------------------------------------------------------------
    // SCALAR VARIABLES
    // ---------------------------------------------------------------------------
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
    // are they the same? (operator==) - should be false because the dataset contains NAN values
    _test( !(dataset1 == dataset2) );
    // and a third one (asignment operator)
    dataset3 = dataset2;
    _test( dataset3.Size() == 3 );
    // should be empty
    dataset3.Clear();
    _test( dataset3.Size() == 0 );
    // resizing
    dataset1.Resize( 5, 5 );
    _info("PropertyData_Test::run: dataset1 after resizing to 5:");
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


/*  TESTING

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
    /// range check; all values included
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
  
    _info("PropertyData_Test::run: scaled and offset dataset1:");
    dataset1.Out();
    _info("PropertyData_Test::run: sqrt of values in dataset1:");
    dataset1.TransformValues( std::sqrt );
    // eliminating the nan value
    dataset1.Value(1) = -1.;
    dataset1.Out();
   
    /// writing stored flag and data values to file
    _info("PropertyData_Test::run: dataset1 written to file and read back to memory:");
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


    // ---------------------------------------------------------------------------
    // VECTOR VARIABLES
    // ---------------------------------------------------------------------------
    _info("PropertyData_Test: testing for VectorVariable:");
    // constructor for all possible csmp variable types
    const size_t DIM3(3);
    PropertyData dataset4( ELEMENT_INTEGRATION_POINT, VECTOR, DIM3 );
    // stick in 3 vectors
    VectorVariable<DIM3>  vc1(PLAIN,1.), vc2(DIRICH,2.), vc3(ANY,3.);
    dataset4.Reserve( 3 );
    pushBack( dataset4, vc1 );
    pushBack( dataset4, vc2 );
    pushBack( dataset4, vc3 );
    dataset4.Out();
    // a copy
    PropertyData dataset5( dataset4 );
    // comparison
    dataset5.Out();
    // are they the same? (operator==) - should be the same
    _test( (dataset4 == dataset5) == true );
    // reading element vc2
    VectorVariable<DIM3> vc_test;
    read( dataset5, 1, vc_test );
    vc_test.Out();
    _test( (vc_test == vc2) == true );
    // writing vc1 into position 3
    store( dataset5, 2, vc1 );
    dataset5.Out();
    _test( (dataset4 == dataset5) == false );
   

    // ---------------------------------------------------------------------------
    // TENSOR VARIABLES
    // ---------------------------------------------------------------------------
    _info("PropertyData_Test: testing for TensorVariable:");
    // constructor for all possible csmp variable types
    PropertyData dataset6( ELEMENT_INTEGRATION_POINT, TENSOR, DIM3 );
    // stick in 3 diagonal tensors
    TensorVariable<DIM3>  ts1(PLAIN,1.), ts2(DIRICH,2.), ts3(FIELD_DATA,3.);
    ts3(0,1) = ts3(1,0) =  1.;
    ts3(1,2) = ts3(2,1) = -3.;
    ts3(2,0) = ts3(0,2) =  1.;
    dataset6.Reserve( 3 );
    pushBack( dataset6, ts1 );
    pushBack( dataset6, ts2 );
    pushBack( dataset6, ts3 );
    dataset6.Out();
    // a copy
    PropertyData dataset7( dataset6 );
    // comparison
    dataset7.Out();
    // are they the same? (operator==) - should be the same
    _test( (dataset6 == dataset7) == true );

    // reading element ts3
    TensorVariable<DIM3> ts_test;
    read( dataset7, 2, ts_test );
    ts_test.Out();
    _test( (ts_test == ts3) == true );
    // (over) writing ts3 into position 3
    store( dataset7, 2, ts_test );
    // should not change the dataset
    _test( (dataset7 == dataset6) == true );


    // ---------------------------------------------------------------------------
    // ARRAY VARIABLES
    // ---------------------------------------------------------------------------
    _info("PropertyData_Test: testing for ArrayVariable:");
    // constructor for all possible csmp variable types
    const size_t array_length(4);
    PropertyData dataset8( NODE, ARRAY, DIM3, array_length );
    // stick in 2 arrays of 4
    ArrayVariable  ary1( array_length, 0., ANY ), ary2( 4, 1., PLAIN);
    dataset8.Reserve( 2 );
    pushBack( dataset8, ary1 );
    pushBack( dataset8, ary2 );
    dataset8.Out();
    // a copy
    PropertyData dataset9( dataset8 );
    // comparison
    dataset9.Out();
    // are they the same? (operator==) - should be the same
    _test( (dataset8 == dataset9) == true );
    // changing last element of ary1 to 4
    ary2(3) = 4.;
    // comparing with that element in the array storage
    dataset9.Value( 1, 3 ) = 4.;
    _equal( ary2(3), dataset9.Value(1,3), numeric_limits<double64>::epsilon() );
   
 } // end run


} // end csmp
