//
//  PropertyData_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 7/04/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "PropertyData_Test.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

using namespace std;

namespace csmp {

/// TESTING THE PUBLIC INTERFACE OF THE CLASS
void PropertyData_Test::run()
{
    // =========================================================================
    // SCALAR VARIABLES  (original tests preserved)
    // =========================================================================
    _info("PropertyData_Test: testing for ScalarVariable:");

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
    if ( verbose_ ) dataset1.Out();
    if ( verbose_ ) dataset2.Out();
    // are they the same? (operator==) - should be false because the dataset contains NAN values
    _test( !(dataset1 == dataset2) );
    // and a third one (assignment operator)
    dataset3 = dataset2;
    _test( dataset3.Size() == 3 );
    // should be empty
    dataset3.Clear();
    _test( dataset3.Size() == 0 );
    // resizing
    dataset1.Resize( 5, 5 );
    _info("PropertyData_Test::run: dataset1 after resizing to 5:");
    if ( verbose_ ) dataset1.Out();
    _test( dataset1.Size() == 5 );

    /// losing entry 2 and 4, 5
    map<size_t,size_t> o_n_elmt_ids;
    o_n_elmt_ids.insert( make_pair(0,0) ); // keep first value: PLAIN:2
    o_n_elmt_ids.insert( make_pair(2,1) ); // keep third value: ANY:NaN
    dataset1.ReduceTo( o_n_elmt_ids );
    if ( verbose_ ) dataset1.Out();
    _test( dataset1.Size() == 2 );

    // accessors / mutators for values (operator[] cannot be overloaded for this)
    /// scalars
    _test( dataset1.Flag(0) == PLAIN );
    _test( dataset1.Flag(1) == ANY );

    /// range check; all values included
    pushBack( dataset1, sc1 );
    pushBack( dataset1, sc2 );
    // pushBack( dataset1, sc3 ); // NaN
    double tmin, tmax;
    dataset1.MinMaxOf( tmin, tmax );
    _equal( tmin, 2., 1e-10 );
    _equal( tmax, 5., 1e-10 );

    /// scales the current variable range to the new one; all values are scaled
    dataset1.ScaleRangeTo( 1., 20. );
    dataset1.MinMaxOf( tmin, tmax );
    _equal( tmin,  1., 1e-10 );
    _equal( tmax, 20., 1e-10 );

    /// add supplied value to all entries
    dataset1.OffsetRangeBy( 15. );
    dataset1.MinMaxOf( tmin, tmax );
    _equal( tmin, 16., 1e-10 );
    _equal( tmax, 35., 1e-10 );

    _info("PropertyData_Test::run: scaled and offset dataset1:");
    if ( verbose_ ) dataset1.Out();
    _info("PropertyData_Test::run: sqrt of values in dataset1:");
    dataset1.TransformValues( std::sqrt );
    // eliminating the nan value
    dataset1.Value(1) = -1.;
    if ( verbose_ ) dataset1.Out();

    /// writing stored flag and data values to file
    _info("PropertyData_Test::run: dataset1 written to file and read back to memory:");
    PropertyData datasetN( dataset1 ); // backup copy

    std::fstream out_fp( "PropertyData_Test", ios::out | ios::binary );
    dataset1.OutBinary( out_fp );
    out_fp.close();

    /// reading stored flag and data values from file
    std::fstream in_fp("PropertyData_Test", ios::in | ios::binary);
    datasetN = inBinaryPropertyData( in_fp );
    in_fp.close();
    if ( verbose_ ) datasetN.Out();
    _test( datasetN == dataset1 );


    // =========================================================================
    // SCALAR — additional coverage
    // =========================================================================

    // --- Observers ---
    {
        PropertyData d( NODE, SCALAR, 2U );
        _test( d.Placement()  == NODE   );
        _test( d.Type()       == SCALAR );
        _test( d.Dim()        == 2U     );
        _test( d.Size()       == 0U     );
        _test( d.Components() == 1U     );
    }

    // --- flagOffset / valueOffset free functions ---
    _test( flagOffset(  SCALAR,      3U, 0U ) == 1U  );
    _test( flagOffset(  ARRAY,       3U, 0U ) == 1U  );
    _test( flagOffset(  VECTOR,      3U, 0U ) == 3U  );
    _test( flagOffset(  TENSOR,      2U, 0U ) == 2U  );
    _test( flagOffset(  FLAGGEDARRAY,3U, 5U ) == 5U  );
    _test( valueOffset( SCALAR,      3U, 0U ) == 1U  );
    _test( valueOffset( VECTOR,      3U, 0U ) == 3U  );
    _test( valueOffset( TENSOR,      3U, 0U ) == 9U  );
    _test( valueOffset( ARRAY,       3U, 7U ) == 7U  );
    _test( valueOffset( FLAGGEDARRAY,3U, 7U ) == 7U  );

    // --- Reserve(n_objects) then PushBack raw ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        d.Reserve( 3U );
        d.PushBack( PLAIN );  d.PushBack( 1.0 );
        d.PushBack( DIRICH ); d.PushBack( 2.0 );
        d.PushBack( ANY );    d.PushBack( 3.0 );
        _test( d.Size() == 3U );
        _test( d.Flag(0) == PLAIN  );
        _test( d.Flag(1) == DIRICH );
        _test( d.Flag(2) == ANY    );
        _equal( d.Value(0), 1.0, 1e-15 );
        _equal( d.Value(1), 2.0, 1e-15 );
        _equal( d.Value(2), 3.0, 1e-15 );
    }

    // --- Resize(n_objects): new entries initialised to ANY / NaN ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        d.Resize( 3U );
        _test( d.Size() == 3U );
        _test( d.Flag(0) == ANY );
        _test( std::isnan( d.Value(0) ) );
    }

    // --- Resize(n_flags, n_values) raw form ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        d.Resize( 2U, 2U );
        _test( d.Size() == 2U );
        _test( d.Flag(0) == ANY );
        _test( std::isnan( d.Value(0) ) );
    }

    // --- Flag and Value mutators ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        d.Resize( 2U );
        d.Flag(0)  = PLAIN;
        d.Value(0) = 42.0;
        _test(  d.Flag(0)  == PLAIN );
        _equal( d.Value(0), 42.0, 1e-15 );
    }

    // --- PushBackFrom ---
    {
        PropertyData src( ELEMENT, SCALAR, 3U );
        ScalarVariable s1(PLAIN, 7.0), s2(DIRICH, 8.0);
        pushBack( src, s1 );
        pushBack( src, s2 );

        PropertyData dst( ELEMENT, SCALAR, 3U );
        dst.PushBackFrom( src, 1U ); // copy second object
        _test(  dst.Flag(0)  == DIRICH );
        _equal( dst.Value(0), 8.0, 1e-15 );
    }

    // --- store / read round-trip for ScalarVariable ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        d.Resize( 3U );
        ScalarVariable s_in(PLAIN, 99.0);
        store( d, 1U, s_in );
        ScalarVariable s_out;
        read( d, 1U, s_out );
        _test(  s_out.Flag() == PLAIN );
        _equal( s_out(),       99.0, 1e-15 );
    }

    // --- ReduceTo with empty map clears container ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        ScalarVariable s(PLAIN, 1.0);
        pushBack( d, s );
        map<size_t,size_t> empty_map;
        d.ReduceTo( empty_map );
        _test( d.Size() == 0U );
    }

    // --- OutBinary returns true on a valid stream ---
    {
        PropertyData d( ELEMENT, SCALAR, 3U );
        ScalarVariable s(PLAIN, 3.14);
        pushBack( d, s );
        std::fstream fp( "PropertyData_scalar_rw", ios::out | ios::binary );
        _test( d.OutBinary( fp ) == true );
        fp.close();

        std::fstream fp2( "PropertyData_scalar_rw", ios::in | ios::binary );
        PropertyData d2 = inBinaryPropertyData( fp2 );
        fp2.close();
        _test( d2 == d );
    }


    // =========================================================================
    // VECTOR VARIABLES  (original tests preserved)
    // =========================================================================
    _info("PropertyData_Test: testing for VectorVariable:");

    const size_t DIM3(3);
    PropertyData dataset4( ELEMENT_INTEGRATION_POINT, VECTOR, DIM3 );
    VectorVariable<DIM3>  vc1(PLAIN,1.), vc2(DIRICH,2.), vc3(ANY,3.);
    dataset4.Reserve( 3 );
    pushBack( dataset4, vc1 );
    pushBack( dataset4, vc2 );
    pushBack( dataset4, vc3 );
    if ( verbose_ ) dataset4.Out();
    PropertyData dataset5( dataset4 );
    if ( verbose_ ) dataset5.Out();
    _test( (dataset4 == dataset5) == true );
    VectorVariable<DIM3> vc_test;
    read( dataset5, 1, vc_test );
    if ( verbose_ ) vc_test.Out();
    _test( (vc_test == vc2) == true );
    store( dataset5, 2, vc1 );
    if ( verbose_ ) dataset5.Out();
    _test( (dataset4 == dataset5) == false );


    // =========================================================================
    // VECTOR — additional coverage
    // =========================================================================

    // --- Observers ---
    {
        PropertyData d( NODE, VECTOR, 2U );
        _test( d.Type()       == VECTOR );
        _test( d.Dim()        == 2U     );
        _test( d.Components() == 2U     );
    }

    // --- Flag and Value accessors (two-argument forms) ---
    {
        PropertyData d( ELEMENT, VECTOR, DIM3 );
        d.Resize( 2U );
        d.Flag(  0U, 1U ) = DIRICH;
        d.Value( 0U, 2U ) = 7.5;
        _test(  d.Flag(  0U, 1U ) == DIRICH );
        _equal( d.Value( 0U, 2U ), 7.5, 1e-15 );
    }

    // --- store / read round-trip ---
    {
        PropertyData d( ELEMENT, VECTOR, DIM3 );
        d.Resize( 3U );
        VectorVariable<DIM3> v_in(PLAIN, 4.0);
        store( d, 1U, v_in );
        VectorVariable<DIM3> v_out;
        read( d, 1U, v_out );
        _test( v_out == v_in );
    }

    // --- PushBackFrom ---
    {
        PropertyData src( ELEMENT, VECTOR, DIM3 );
        pushBack( src, vc1 );
        pushBack( src, vc2 );
        PropertyData dst( ELEMENT, VECTOR, DIM3 );
        dst.PushBackFrom( src, 1U );
        VectorVariable<DIM3> v_out;
        read( dst, 0U, v_out );
        _test( v_out == vc2 );
    }

    // --- ReduceTo ---
    {
        PropertyData d( ELEMENT, VECTOR, DIM3 );
        pushBack( d, vc1 );
        pushBack( d, vc2 );
        pushBack( d, vc3 );
        map<size_t,size_t> m;
        m.insert( make_pair(0U, 0U) );
        m.insert( make_pair(2U, 1U) );
        d.ReduceTo( m );
        _test( d.Components() * 2U == d.Size() );
        VectorVariable<DIM3> v0, v1;
        read( d, 0U, v0 );
        read( d, 1U, v1 );
        _test( v0 == vc1 );
        _test( v1 == vc3 );
    }

    // --- MinMaxOf, OffsetRangeBy, ScaleRangeTo, TransformValues ---
    {
        PropertyData d( ELEMENT, VECTOR, DIM3 );
        pushBack( d, vc1 ); // all components 1.0
        pushBack( d, vc2 ); // all components 2.0
        pushBack( d, vc3 ); // all components 3.0
        double vmin, vmax;
        d.MinMaxOf( vmin, vmax );
        _equal( vmin, 1.0, 1e-10 );
        _equal( vmax, 3.0, 1e-10 );
        d.OffsetRangeBy( 10.0 );
        d.MinMaxOf( vmin, vmax );
        _equal( vmin, 11.0, 1e-10 );
        _equal( vmax, 13.0, 1e-10 );
        d.ScaleRangeTo( 0.0, 1.0 );
        d.MinMaxOf( vmin, vmax );
        _equal( vmin, 0.0, 1e-10 );
        _equal( vmax, 1.0, 1e-10 );
        d.TransformValues( std::sqrt );
        d.MinMaxOf( vmin, vmax );
        _equal( vmin, 0.0, 1e-10 );
        _equal( vmax, 1.0, 1e-10 );
    }

    // --- Binary round-trip ---
    {
        PropertyData d( ELEMENT, VECTOR, DIM3 );
        pushBack( d, vc1 );
        pushBack( d, vc2 );
        pushBack( d, vc3 );
        std::fstream fp( "PropertyData_vector_rw", ios::out | ios::binary );
        _test( d.OutBinary( fp ) == true );
        fp.close();
        std::fstream fp2( "PropertyData_vector_rw", ios::in | ios::binary );
        PropertyData d2 = inBinaryPropertyData( fp2 );
        fp2.close();
        _test( d2 == d );
    }


    // =========================================================================
    // TENSOR VARIABLES  (original tests preserved)
    // =========================================================================
    _info("PropertyData_Test: testing for TensorVariable:");

    PropertyData dataset6( ELEMENT_INTEGRATION_POINT, TENSOR, DIM3 );
    TensorVariable<DIM3>  ts1(PLAIN,1.), ts2(DIRICH,2.), ts3(FIELD_DATA,3.);
    ts3(0,1) = ts3(1,0) =  1.;
    ts3(1,2) = ts3(2,1) = -3.;
    ts3(2,0) = ts3(0,2) =  1.;
    dataset6.Reserve( 3 );
    pushBack( dataset6, ts1 );
    pushBack( dataset6, ts2 );
    pushBack( dataset6, ts3 );
    if ( verbose_ ) dataset6.Out();
    PropertyData dataset7( dataset6 );
    if ( verbose_ ) dataset7.Out();
    _test( (dataset6 == dataset7) == true );
    TensorVariable<DIM3> ts_test;
    read( dataset7, 2, ts_test );
    if ( verbose_ ) ts_test.Out();
    _test( (ts_test == ts3) == true );
    store( dataset7, 2, ts_test );
    _test( (dataset7 == dataset6) == true );


    // =========================================================================
    // TENSOR — additional coverage
    // =========================================================================

    // --- Observers ---
    {
        PropertyData d( NODE, TENSOR, 2U );
        _test( d.Type()       == TENSOR );
        _test( d.Dim()        == 2U     );
        _test( d.Components() == 4U     ); // 2*2
    }

    // --- Flag (diagonal) and Value (full) accessors ---
    {
        PropertyData d( ELEMENT, TENSOR, DIM3 );
        d.Resize( 1U );
        d.Flag(  0U, 1U, 1U ) = DIRICH;
        d.Value( 0U, 0U, 1U ) = 3.14;
        _test(  d.Flag(  0U, 1U, 1U ) == DIRICH );
        _equal( d.Value( 0U, 0U, 1U ), 3.14, 1e-15 );
    }

    // --- store / read round-trip ---
    {
        PropertyData d( ELEMENT, TENSOR, DIM3 );
        d.Resize( 2U );
        store( d, 0U, ts1 );
        store( d, 1U, ts3 );
        TensorVariable<DIM3> t0, t1;
        read( d, 0U, t0 );
        read( d, 1U, t1 );
        _test( t0 == ts1 );
        _test( t1 == ts3 );
    }

    // --- PushBackFrom ---
    {
        PropertyData src( ELEMENT, TENSOR, DIM3 );
        pushBack( src, ts1 );
        pushBack( src, ts2 );
        PropertyData dst( ELEMENT, TENSOR, DIM3 );
        dst.PushBackFrom( src, 1U );
        TensorVariable<DIM3> t_out;
        read( dst, 0U, t_out );
        _test( t_out == ts2 );
    }

    // --- ReduceTo ---
    {
        PropertyData d( ELEMENT, TENSOR, DIM3 );
        pushBack( d, ts1 );
        pushBack( d, ts2 );
        pushBack( d, ts3 );
        map<size_t,size_t> m;
        m.insert( make_pair(1U, 0U) );
        m.insert( make_pair(2U, 1U) );
        d.ReduceTo( m );
        TensorVariable<DIM3> t0, t1;
        read( d, 0U, t0 );
        read( d, 1U, t1 );
        _test( t0 == ts2 );
        _test( t1 == ts3 );
    }

    // --- Binary round-trip ---
    {
        PropertyData d( ELEMENT, TENSOR, DIM3 );
        pushBack( d, ts1 );
        pushBack( d, ts2 );
        pushBack( d, ts3 );
        std::fstream fp( "PropertyData_tensor_rw", ios::out | ios::binary );
        _test( d.OutBinary( fp ) == true );
        fp.close();
        std::fstream fp2( "PropertyData_tensor_rw", ios::in | ios::binary );
        PropertyData d2 = inBinaryPropertyData( fp2 );
        fp2.close();
        _test( d2 == d );
    }


    // =========================================================================
    // ARRAY VARIABLES  (original tests preserved)
    // =========================================================================
    _info("PropertyData_Test: testing for ArrayVariable:");

    const size_t array_length(4);
    PropertyData dataset8( NODE, ARRAY, DIM3, array_length );
    ArrayVariable  ary1( array_length, 0., ANY ), ary2( 4, 1., PLAIN);
    dataset8.Reserve( 2 );
    pushBack( dataset8, ary1 );
    pushBack( dataset8, ary2 );
    if ( verbose_ ) dataset8.Out();
    PropertyData dataset9( dataset8 );
    if ( verbose_ ) dataset9.Out();
    _test( (dataset8 == dataset9) == true );
    ary2(3) = 4.;
    dataset9.Value( 1, 3 ) = 4.;
    _equal( ary2(3), dataset9.Value(1,3), numeric_limits<double>::epsilon() );


    // =========================================================================
    // ARRAY — additional coverage
    // =========================================================================

    // --- Observers ---
    {
        PropertyData d( NODE, ARRAY, DIM3, array_length );
        _test( d.Type()       == ARRAY       );
        _test( d.Components() == array_length );
        _test( d.Dim()        == DIM3         );
    }

    // --- Flag accessor (single flag per array object) ---
    {
        PropertyData d( NODE, ARRAY, DIM3, array_length );
        d.Resize( 2U );
        d.Flag( 0U, 0U ) = PLAIN;
        _test( d.Flag( 0U, 0U ) == PLAIN );
    }

    // --- store / read round-trip ---
    {
        PropertyData d( NODE, ARRAY, DIM3, array_length );
        d.Resize( 2U );
        ArrayVariable a_in( array_length, 5.0, PLAIN );
        store( d, 0U, a_in );
        ArrayVariable a_out;
        read( d, 0U, a_out );
        _test( a_out.Flag() == PLAIN );
        for ( uint32_t i = 0U; i < array_length; ++i )
            _equal( a_out[i], 5.0, 1e-15 );
    }

    // --- PushBackFrom ---
    {
        PropertyData src( NODE, ARRAY, DIM3, array_length );
        pushBack( src, ary1 );
        pushBack( src, ary2 );
        PropertyData dst( NODE, ARRAY, DIM3, array_length );
        dst.PushBackFrom( src, 0U );
        ArrayVariable a_out;
        read( dst, 0U, a_out );
        _test( a_out.Flag() == ary1.Flag() );
        for ( uint32_t i = 0U; i < array_length; ++i )
            _equal( a_out[i], ary1[i], 1e-15 );
    }

    // --- ReduceTo ---
    {
        PropertyData d( NODE, ARRAY, DIM3, array_length );
        ArrayVariable a0( array_length, 1.0, PLAIN  );
        ArrayVariable a1( array_length, 2.0, DIRICH );
        ArrayVariable a2( array_length, 3.0, ANY    );
        pushBack( d, a0 );
        pushBack( d, a1 );
        pushBack( d, a2 );
        map<size_t,size_t> m;
        m.insert( make_pair(0U, 0U) );
        m.insert( make_pair(2U, 1U) );
        d.ReduceTo( m );
        ArrayVariable r0, r1;
        read( d, 0U, r0 );
        read( d, 1U, r1 );
        _equal( r0[0], 1.0, 1e-15 );
        _equal( r1[0], 3.0, 1e-15 );
    }

    // --- Binary round-trip ---
    {
        PropertyData d( NODE, ARRAY, DIM3, array_length );
        pushBack( d, ary1 );
        pushBack( d, ary2 );
        std::fstream fp( "PropertyData_array_rw", ios::out | ios::binary );
        _test( d.OutBinary( fp ) == true );
        fp.close();
        std::fstream fp2( "PropertyData_array_rw", ios::in | ios::binary );
        PropertyData d2 = inBinaryPropertyData( fp2 );
        fp2.close();
        _test( d2 == d );
    }


    // =========================================================================
    // FLAGGED ARRAY VARIABLES
    // =========================================================================
    _info("PropertyData_Test: testing for FlaggedArrayVariable:");

    const uint32_t fa_length(3);

    // --- Observers ---
    {
        PropertyData d( NODE, FLAGGEDARRAY, DIM3, fa_length );
        _test( d.Type()       == FLAGGEDARRAY );
        _test( d.Components() == fa_length    );
        _test( d.Dim()        == DIM3         );
    }

    // --- pushBack / Flag / Value accessors ---
    {
        PropertyData d( NODE, FLAGGEDARRAY, DIM3, fa_length );
        FlaggedArrayVariable fa1( fa_length );
        fa1(0) = 1.0;  fa1.Flag(0) = PLAIN;
        fa1(1) = 2.0;  fa1.Flag(1) = DIRICH;
        fa1(2) = 3.0;  fa1.Flag(2) = ANY;
        pushBack( d, fa1 );
        _test(  d.Flag(  0U, 0U ) == PLAIN  );
        _test(  d.Flag(  0U, 1U ) == DIRICH );
        _test(  d.Flag(  0U, 2U ) == ANY    );
        _equal( d.Value( 0U, 0U ), 1.0, 1e-15 );
        _equal( d.Value( 0U, 1U ), 2.0, 1e-15 );
        _equal( d.Value( 0U, 2U ), 3.0, 1e-15 );
    }

    // --- store / read round-trip ---
    {
        PropertyData d( NODE, FLAGGEDARRAY, DIM3, fa_length );
        d.Resize( 2U );
        FlaggedArrayVariable fa_in( fa_length );
        fa_in(0) = 10.0; fa_in.Flag(0) = PLAIN;
        fa_in(1) = 20.0; fa_in.Flag(1) = DIRICH;
        fa_in(2) = 30.0; fa_in.Flag(2) = ANY;
        store( d, 1U, fa_in );
        FlaggedArrayVariable fa_out;
        read( d, 1U, fa_out );
        _test(  fa_out.Flag(0) == PLAIN  );
        _test(  fa_out.Flag(1) == DIRICH );
        _test(  fa_out.Flag(2) == ANY    );
        _equal( fa_out[0], 10.0, 1e-15 );
        _equal( fa_out[1], 20.0, 1e-15 );
        _equal( fa_out[2], 30.0, 1e-15 );
    }

    // --- PushBackFrom ---
    {
        PropertyData src( NODE, FLAGGEDARRAY, DIM3, fa_length );
        FlaggedArrayVariable fa_a( fa_length ), fa_b( fa_length );
        for ( uint32_t i = 0U; i < fa_length; ++i ) {
            fa_a(i) = static_cast<double>(i);     fa_a.Flag(i) = PLAIN;
            fa_b(i) = static_cast<double>(i) * 2; fa_b.Flag(i) = DIRICH;
        }
        pushBack( src, fa_a );
        pushBack( src, fa_b );
        PropertyData dst( NODE, FLAGGEDARRAY, DIM3, fa_length );
        dst.PushBackFrom( src, 1U );
        FlaggedArrayVariable fa_out;
        read( dst, 0U, fa_out );
        for ( uint32_t i = 0U; i < fa_length; ++i ) {
            _test(  fa_out.Flag(i) == DIRICH );
            _equal( fa_out[i], fa_b[i], 1e-15 );
        }
    }

    // --- ReduceTo ---
    {
        PropertyData d( NODE, FLAGGEDARRAY, DIM3, fa_length );
        FlaggedArrayVariable fa0( fa_length ), fa1( fa_length ), fa2( fa_length );
        for ( uint32_t i = 0U; i < fa_length; ++i ) {
            fa0(i) = 1.0; fa0.Flag(i) = PLAIN;
            fa1(i) = 2.0; fa1.Flag(i) = DIRICH;
            fa2(i) = 3.0; fa2.Flag(i) = ANY;
        }
        pushBack( d, fa0 );
        pushBack( d, fa1 );
        pushBack( d, fa2 );
        map<size_t,size_t> m;
        m.insert( make_pair(0U, 0U) );
        m.insert( make_pair(2U, 1U) );
        d.ReduceTo( m );
        FlaggedArrayVariable r0, r1;
        read( d, 0U, r0 );
        read( d, 1U, r1 );
        for ( uint32_t i = 0U; i < fa_length; ++i ) {
            _test(  r0.Flag(i) == PLAIN );
            _equal( r0[i], 1.0, 1e-15 );
            _test(  r1.Flag(i) == ANY   );
            _equal( r1[i], 3.0, 1e-15 );
        }
    }

    // --- operator== distinguishes different flag patterns ---
    {
        PropertyData d1( NODE, FLAGGEDARRAY, DIM3, fa_length );
        PropertyData d2( NODE, FLAGGEDARRAY, DIM3, fa_length );
        FlaggedArrayVariable fa( fa_length );
        fa(0) = 1.0; fa.Flag(0) = PLAIN;
        fa(1) = 2.0; fa.Flag(1) = PLAIN;
        fa(2) = 3.0; fa.Flag(2) = PLAIN;
        pushBack( d1, fa );
        fa.Flag(2) = DIRICH;
        pushBack( d2, fa );
        _test( !(d1 == d2) );
    }

    // --- Binary round-trip ---
    {
        PropertyData d( NODE, FLAGGEDARRAY, DIM3, fa_length );
        FlaggedArrayVariable fa( fa_length );
        for ( uint32_t i = 0U; i < fa_length; ++i ) {
            fa(i) = static_cast<double>(i + 1);
            fa.Flag(i) = PLAIN;
        }
        pushBack( d, fa );
        std::fstream fp( "PropertyData_flaggedarray_rw", ios::out | ios::binary );
        _test( d.OutBinary( fp ) == true );
        fp.close();
        std::fstream fp2( "PropertyData_flaggedarray_rw", ios::in | ios::binary );
        PropertyData d2 = inBinaryPropertyData( fp2 );
        fp2.close();
        _test( d2 == d );
    }

} // end run

} // namespace csmp
