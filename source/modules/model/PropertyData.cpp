// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PropertyData.h"

#include "VectorVariable.h"
#include "TensorVariable.h"
#include "binaryReadWrite.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

// ============================================================================
// Constructor
// ============================================================================

PropertyData::PropertyData( PLACEMENT     place,
                             VARIABLE_TYPE type,
                             uint32_t      spatial_dimension,
                             uint32_t      array_length )
    : place_      ( place ),
      type_       ( type ),
      dim_        ( spatial_dimension ),
      flag_stride_( flagOffset( type, spatial_dimension, array_length ) ),
      data_stride_( valueOffset( type, spatial_dimension, array_length ) )
{}


// ============================================================================
// PushBackFrom
// ============================================================================

void PropertyData::PushBackFrom( const PropertyData& prop, size_t nth_value )
{
    assert( place_       == prop.place_ );
    assert( type_        == prop.type_ );
    assert( dim_         == prop.dim_ );
    assert( flag_stride_ == prop.flag_stride_ );
    assert( data_stride_ == prop.data_stride_ );

    for ( uint32_t i = 0U; i < flag_stride_; ++i )
        flags_.push_back( prop.flags_[ nth_value * flag_stride_ + i ] );

    for ( uint32_t i = 0U; i < data_stride_; ++i )
        data_.push_back( prop.data_[ nth_value * data_stride_ + i ] );
}


// ============================================================================
// ReduceTo
// ============================================================================

void PropertyData::ReduceTo( const map<size_t,size_t>& o_n_elmt_ids )
{
    if ( o_n_elmt_ids.empty() ) {
        cerr << "\nPropertyData::ReduceTo: supplied map is empty"
                " — container will be cleared.\n";
        Clear();
        return;
    }

    vector<VARIABLE_FLAG> new_flags( o_n_elmt_ids.size() * flag_stride_ );
    vector<double>        new_data ( o_n_elmt_ids.size() * data_stride_ );

    if ( type_ == SCALAR || type_ == ARRAY ) {
        // SCALAR: flag_stride_ == data_stride_ == 1.
        // ARRAY:  flag_stride_ == 1, data_stride_ == array_length.
        for ( const auto& [old_idx, new_idx] : o_n_elmt_ids ) {
            new_flags[ new_idx ] = flags_[ old_idx ];
            for ( uint32_t i = 0U; i < data_stride_; ++i )
                new_data[ new_idx * data_stride_ + i ] =
                    data_[ old_idx * data_stride_ + i ];
        }
    }
    else if ( type_ == VECTOR ) {
        // flag_stride_ == data_stride_ == dim_
        for ( const auto& [old_idx, new_idx] : o_n_elmt_ids ) {
            for ( uint32_t i = 0U; i < data_stride_; ++i ) {
                new_flags[ new_idx * flag_stride_ + i ] =
                    flags_[ old_idx * flag_stride_ + i ];
                new_data [ new_idx * data_stride_ + i ] =
                    data_ [ old_idx * data_stride_ + i ];
            }
        }
    }
    else if ( type_ == FLAGGEDARRAY ) {
        // flag_stride_ == data_stride_ == array_length
        for ( const auto& [old_idx, new_idx] : o_n_elmt_ids ) {
            for ( uint32_t i = 0U; i < data_stride_; ++i ) {
                new_flags[ new_idx * flag_stride_ + i ] =
                    flags_[ old_idx * flag_stride_ + i ];
                new_data [ new_idx * data_stride_ + i ] =
                    data_ [ old_idx * data_stride_ + i ];
            }
        }
    }
    else if ( type_ == TENSOR ) {
        // flag_stride_ == dim_, data_stride_ == dim_ * dim_
        for ( const auto& [old_idx, new_idx] : o_n_elmt_ids ) {
            for ( uint32_t i = 0U; i < flag_stride_; ++i )
                new_flags[ new_idx * flag_stride_ + i ] =
                    flags_[ old_idx * flag_stride_ + i ];
            for ( uint32_t i = 0U; i < data_stride_; ++i )
                new_data [ new_idx * data_stride_ + i ] =
                    data_ [ old_idx * data_stride_ + i ];
        }
    }

    flags_ = std::move( new_flags );
    data_  = std::move( new_data );

    // Release excess capacity.
    vector<VARIABLE_FLAG>( flags_ ).swap( flags_ );
    vector<double>( data_ ).swap( data_ );
}


// ============================================================================
// ScaleRangeTo
// ============================================================================

void PropertyData::ScaleRangeTo( double tmin, double tmax )
{
    double old_min, old_max;
    MinMaxOf( old_min, old_max );

    // Already within tolerance — nothing to do.
    if ( fabs( fabs(tmax) - fabs(old_max) ) < numeric_limits<double>::epsilon() &&
         fabs( fabs(tmin) - fabs(old_min) ) < numeric_limits<double>::epsilon() )
        return;

    const double old_range = old_max - old_min;
    const double new_range = tmax - tmin;

    transform( data_.begin(), data_.end(), data_.begin(),
               [tmin, old_min, old_range, new_range]( double val ) noexcept {
                   return tmin + ( (val - old_min) / old_range ) * new_range;
               } );

    // Epsilon-based post-check — avoids fragile exact equality after
    // floating-point arithmetic.
    double new_min, new_max;
    MinMaxOf( new_min, new_max );
    const double tol = fabs( new_range ) * numeric_limits<double>::epsilon() * 10.0;
    if ( fabs( new_min - tmin ) > tol || fabs( new_max - tmax ) > tol ) {
        cout << "\nPropertyData::ScaleRangeTo: scaling failed.\n";
        Out();
        throw range_error( "PropertyData::ScaleRangeTo" );
    }
}


// ============================================================================
// OutBinary
// ============================================================================

bool PropertyData::OutBinary( fstream& fp ) const
{
    const int8_t  var_placement   = static_cast<int8_t> ( place_       );
    const int8_t  var_type        = static_cast<int8_t> ( type_        );
    const int32_t var_dim         = static_cast<int32_t>( dim_         );
    const int32_t var_flag_stride = static_cast<int32_t>( flag_stride_ );
    const int32_t var_data_stride = static_cast<int32_t>( data_stride_ );

    fp.write( reinterpret_cast<const char*>(&var_placement),   sizeof(int8_t)  );
    fp.write( reinterpret_cast<const char*>(&var_type),        sizeof(int8_t)  );
    fp.write( reinterpret_cast<const char*>(&var_dim),         sizeof(int32_t) );
    fp.write( reinterpret_cast<const char*>(&var_flag_stride), sizeof(int32_t) );
    fp.write( reinterpret_cast<const char*>(&var_data_stride), sizeof(int32_t) );

    // Convert VARIABLE_FLAG -> int8_t for compact binary storage.
    vector<int8_t> flags;
    flags.reserve( flags_.size() );
    transform( flags_.begin(), flags_.end(), back_inserter( flags ),
               []( VARIABLE_FLAG f ) noexcept -> int8_t {
                   return static_cast<int8_t>( f );
               } );

    const bool ok_flags = binaryFileWrite( fp, flags );
    const bool ok_data  = binaryFileWrite( fp, data_ );

    return ok_flags && ok_data;
}


// ============================================================================
// inBinaryPropertyData
// ============================================================================

PropertyData inBinaryPropertyData( fstream& fp )
{
    int8_t  var_placement  ( UNSPECIFIED );
    int8_t  var_type       ( UNSPECIFIED );
    int32_t var_dim        ( 0 );
    int32_t var_flag_stride( 0 );
    int32_t var_data_stride( 0 );

    fp.read( reinterpret_cast<char*>(&var_placement),   sizeof(int8_t)  );
    fp.read( reinterpret_cast<char*>(&var_type),        sizeof(int8_t)  );
    fp.read( reinterpret_cast<char*>(&var_dim),         sizeof(int32_t) );
    fp.read( reinterpret_cast<char*>(&var_flag_stride), sizeof(int32_t) );
    fp.read( reinterpret_cast<char*>(&var_data_stride), sizeof(int32_t) );

    const VARIABLE_TYPE vtype = static_cast<VARIABLE_TYPE>( var_type );

    // Array length is encoded as the data stride for array types.
    const uint32_t array_length =
        ( vtype == ARRAY || vtype == FLAGGEDARRAY )
        ? static_cast<uint32_t>( var_data_stride )
        : 0U;

    PropertyData data( static_cast<PLACEMENT>( var_placement ),
                       vtype,
                       static_cast<uint32_t>( var_dim ),
                       array_length );

    // Read raw flag bytes and value doubles from file.
    vector<int8_t> raw_flags;
    binaryFileRead( fp, raw_flags );

    vector<double> values;
    binaryFileRead( fp, values );

    // Convert int8_t back to VARIABLE_FLAG.
    vector<VARIABLE_FLAG> flags;
    flags.reserve( raw_flags.size() );
    transform( raw_flags.begin(), raw_flags.end(), back_inserter( flags ),
               []( int8_t f ) noexcept -> VARIABLE_FLAG {
                   return static_cast<VARIABLE_FLAG>( f );
               } );

    data.Reserve( flags.size(), values.size() );
    for ( auto f : flags  ) data.PushBack( f );
    for ( auto v : values ) data.PushBack( v );

    return data;
}


// ============================================================================
// operator==
// ============================================================================

bool PropertyData::operator==( const PropertyData& d ) const noexcept
{
    return place_       == d.place_
        && type_        == d.type_
        && dim_         == d.dim_
        && flag_stride_ == d.flag_stride_
        && data_stride_ == d.data_stride_
        && flags_       == d.flags_
        && data_        == d.data_;
}


// ============================================================================
// Out
// ============================================================================

void PropertyData::Out() const
{
    cout << "\nPropertyData::Out: storage for "
         << data_.size() / data_stride_ << " "
         << string( parseType( type_ ) ) << " objects.\n"
         << "Data placement:    " << string( parsePlacement( place_ ) ) << "\n"
         << "Spatial dimension: " << dim_         << "\n"
         << "Flag stride:       " << flag_stride_ << "\n"
         << "Data stride:       " << data_stride_ << "\n";

    cout << "\nFlag values:\n";
    for ( auto f : flags_ )
        cout << parseStatus( f ) << " ";
    cout << "\n";

    cout << "\nVariable component values:\n";
    for ( auto v : data_ )
        cout << v << " ";
    cout << "\n";
}

} // namespace csmp

