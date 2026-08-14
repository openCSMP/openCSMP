#ifndef CSMP_PROPERTY_DATA_H
#define CSMP_PROPERTY_DATA_H

#include "CSMP_definitions.h"
#include "FlaggedArrayVariable.h"
#include "TensorVariable.h"

namespace csmp {

class PropertyData;

/// Append a whole csmp variable to a PropertyData container.
void pushBack( PropertyData&, const ScalarVariable& );
void pushBack( PropertyData&, const ArrayVariable& );
void pushBack( PropertyData&, const FlaggedArrayVariable& );
template<uint32_t dim> void pushBack( PropertyData&, const VectorVariable<dim>& );
template<uint32_t dim> void pushBack( PropertyData&, const TensorVariable<dim>& );

/// Store a csmp variable at a given position: store( data, i, makeScalar(ANY, 0.) );
template<typename csmp_var_type>
void store( PropertyData&, size_t position, const csmp_var_type& );
template<uint32_t dim>
void store( PropertyData&, size_t position, const VectorVariable<dim>& );
template<uint32_t dim>
void store( PropertyData&, size_t position, const TensorVariable<dim>& );

/// Read a csmp variable from a given position: VectorVariable<3U> a; read( data, i, a );
template<typename csmp_var_type>
void read( const PropertyData&, size_t position, csmp_var_type& );
template<uint32_t dim>
void read( const PropertyData&, size_t position, VectorVariable<dim>& );
template<uint32_t dim>
void read( const PropertyData&, size_t position, TensorVariable<dim>& );


/// Calculate the per-object flag count for a given variable type and geometry.
inline uint32_t flagOffset( VARIABLE_TYPE type,
                             uint32_t      spatial_dimension,
                             uint32_t      array_length ) noexcept;

/// Calculate the per-object value count for a given variable type and geometry.
inline uint32_t valueOffset( VARIABLE_TYPE type,
                              uint32_t      spatial_dimension,
                              uint32_t      array_length ) noexcept;


/**
  @brief Flexible-size container for CSMP data types, providing range-check
         and conversion functions.

  @author S.K. Matthai
  @date   8/4/2016

  @section motivation Motivation

  Efficient data transfer object with a compact binary interface for storing
  variables to disk.

  @section design Design Intent

  A more disk-storage-efficient alternative to FEM_Data in VSet, designed to
  streamline variable exchange and storage to disk and to allow variables to be
  read selectively from file.

  @section examples Application Examples

  Create a container for 1000 arrays, initialise it and write it to disk:

  @code
  PropertyData data( NODE, ARRAY, DIM3, 21 );
  data.Reserve( 1000 );
  ArrayVariable array( 21 );
  // ... fill array ...
  for ( auto& a : arrays )  pushBack( data, a );
  data.OutBinary( fp );
  @endcode

  @attention Use the free functions read() and store() to retrieve or set
  entire csmp variables; there are no member-level Read() or Store() accessors.
*/
class PropertyData {
  public:

    /// Construct a typed container for the given placement, variable type,
    /// spatial dimension, and (for array types) array length.
    PropertyData( PLACEMENT, VARIABLE_TYPE, uint32_t dim,
                  uint32_t array_length = 0U );

    /// Equality comparison.
    bool operator==( const PropertyData& ) const noexcept;

    // -------------------------------------------------------------------
    // Observers
    // -------------------------------------------------------------------

    PLACEMENT     Placement()  const noexcept;
    VARIABLE_TYPE Type()       const noexcept;
    uint32_t      Dim()        const noexcept;
    size_t        Size()       const noexcept;
    /// Number of value components stored per variable (equals array length
    /// for array types).
    uint32_t      Components() const noexcept;

    // -------------------------------------------------------------------
    // Container sizing
    // -------------------------------------------------------------------

    /// Clear all stored flags and values.
    void Clear() noexcept;

    /// Reserve storage for n_objects of the type this container was
    /// constructed for.
    void Reserve( size_t n_objects );

    /// Reserve raw storage without regard to variable-type layout.
    void Reserve( size_t flag_capacity, size_t value_capacity );

    /// Resize to hold n_objects; new entries are initialised to ANY / NaN.
    void Resize( size_t n_objects );

    /// Resize raw storage; new entries are initialised to ANY / NaN.
    void Resize( size_t n_flags, size_t n_values );

    /// Retain only those object indices present in the map
    /// (old index -> new index).
    void ReduceTo( const std::map<size_t,size_t>& o_n_elmt_ids );

    // -------------------------------------------------------------------
    // Raw inserters
    // Use the pushBack() free functions to insert whole csmp variables.
    // -------------------------------------------------------------------

    /// Raw flag inserter; caller must maintain correct stride.
    void PushBack( VARIABLE_FLAG f ) noexcept;

    /// Raw value inserter; caller must maintain correct stride.
    void PushBack( double v ) noexcept;

    /// Copy the nth object from another compatible PropertyData into this one.
    void PushBackFrom( const PropertyData&, size_t nth_value );

    // -------------------------------------------------------------------
    // Flag accessors / mutators
    // -------------------------------------------------------------------

    VARIABLE_FLAG& Flag( size_t nth_value ) noexcept;
    VARIABLE_FLAG  Flag( size_t nth_value ) const noexcept;

    VARIABLE_FLAG& Flag( size_t nth_value, uint32_t ith_dim ) noexcept;
    VARIABLE_FLAG  Flag( size_t nth_value, uint32_t ith_dim ) const noexcept;

    VARIABLE_FLAG& Flag( size_t nth_value,
                         uint32_t ith_row, uint32_t jth_col ) noexcept;
    VARIABLE_FLAG  Flag( size_t nth_value,
                         uint32_t ith_row, uint32_t jth_col ) const noexcept;

    // -------------------------------------------------------------------
    // Value accessors / mutators
    // -------------------------------------------------------------------

    double& Value( size_t nth_value ) noexcept;
    double  Value( size_t nth_value ) const noexcept;

    double& Value( size_t nth_value, uint32_t ith_dim ) noexcept;
    double  Value( size_t nth_value, uint32_t ith_dim ) const noexcept;

    double& Value( size_t nth_value,
                   uint32_t ith_row, uint32_t jth_col ) noexcept;
    double  Value( size_t nth_value,
                   uint32_t ith_row, uint32_t jth_col ) const noexcept;

    // -------------------------------------------------------------------
    // Range utilities
    // -------------------------------------------------------------------

    /// Compute the minimum and maximum of all stored values.
    void MinMaxOf( double& tmin, double& tmax ) const noexcept;

    /// Scale all stored values so that their range maps to [tmin, tmax].
    void ScaleRangeTo( double tmin, double tmax );

    /// Add a constant offset to every stored value.
    void OffsetRangeBy( double offset ) noexcept;

    /// Apply a unary math function (e.g. std::sqrt) to every stored value.
    void TransformValues( double (*f)(double) ) noexcept;

    // -------------------------------------------------------------------
    // I/O
    // -------------------------------------------------------------------

    /// Write all flag and value data to a binary file stream.
    bool OutBinary( std::fstream& fp ) const;

    /// Diagnostic dump to stdout.
    void Out() const;

  private:
    PLACEMENT                  place_       = ELEMENT; ///< placement of variable
    VARIABLE_TYPE              type_        = SCALAR;  ///< any of scalar..flagged array
    uint32_t                   dim_         = 3U;      ///< spatial dimension
    uint32_t                   flag_stride_ = 1U;      ///< per-object flag count
    uint32_t                   data_stride_ = 1U;      ///< per-object value count
    std::vector<VARIABLE_FLAG> flags_;                 ///< variable flags
    std::vector<double>        data_;                  ///< variable values
};


// ============================================================================
// Inline definitions — flagOffset / valueOffset
// ============================================================================

inline uint32_t flagOffset( VARIABLE_TYPE type,
                             uint32_t      spatial_dimension,
                             uint32_t      array_length ) noexcept
{
    if ( type == SCALAR || type == ARRAY  ) return 1U;
    if ( type == VECTOR || type == TENSOR ) return spatial_dimension;
    // FLAGGEDARRAY
    return array_length;
}

inline uint32_t valueOffset( VARIABLE_TYPE type,
                              uint32_t      spatial_dimension,
                              uint32_t      array_length ) noexcept
{
    if ( type == SCALAR ) return 1U;
    if ( type == VECTOR ) return spatial_dimension;
    if ( type == TENSOR ) return spatial_dimension * spatial_dimension;
    // ARRAY and FLAGGEDARRAY
    return array_length;
}


// ============================================================================
// Inline definitions — observers
// ============================================================================

inline PLACEMENT     PropertyData::Placement()  const noexcept { return place_; }
inline VARIABLE_TYPE PropertyData::Type()       const noexcept { return type_;  }
inline uint32_t      PropertyData::Dim()        const noexcept { return dim_;   }
inline size_t        PropertyData::Size()       const noexcept { return data_.size(); }
inline uint32_t      PropertyData::Components() const noexcept { return data_stride_; }


// ============================================================================
// Inline definitions — container sizing
// ============================================================================

inline void PropertyData::Clear() noexcept
{
    flags_.clear();
    data_.clear();
}

inline void PropertyData::Reserve( size_t n_objects )
{
    flags_.reserve( n_objects * flag_stride_ );
    data_.reserve(  n_objects * data_stride_ );
}

inline void PropertyData::Reserve( size_t flag_capacity, size_t value_capacity )
{
    flags_.reserve( flag_capacity );
    data_.reserve(  value_capacity );
}

inline void PropertyData::Resize( size_t n_objects )
{
    flags_.resize( n_objects * flag_stride_, ANY );
    data_.resize(  n_objects * data_stride_,
                   std::numeric_limits<double>::quiet_NaN() );
}

inline void PropertyData::Resize( size_t n_flags, size_t n_values )
{
    flags_.resize( n_flags,  ANY );
    data_.resize(  n_values, std::numeric_limits<double>::quiet_NaN() );
}


// ============================================================================
// Inline definitions — raw inserters
// ============================================================================

inline void PropertyData::PushBack( VARIABLE_FLAG f ) noexcept
{
    flags_.push_back( f );
}

inline void PropertyData::PushBack( double v ) noexcept
{
    data_.push_back( v );
}


// ============================================================================
// Inline definitions — flag accessors / mutators
// ============================================================================

inline VARIABLE_FLAG& PropertyData::Flag( size_t nth_value ) noexcept
{
    assert( type_ == SCALAR );
    assert( !flags_.empty() );
    assert( nth_value < flags_.size() );
    return flags_[ nth_value ];
}

inline VARIABLE_FLAG PropertyData::Flag( size_t nth_value ) const noexcept
{
    assert( type_ == SCALAR );
    assert( !flags_.empty() );
    assert( nth_value < flags_.size() );
    return flags_[ nth_value ];
}

inline VARIABLE_FLAG& PropertyData::Flag( size_t nth_value, uint32_t ith_dim ) noexcept
{
    assert( type_ == VECTOR || type_ == ARRAY || type_ == FLAGGEDARRAY );
    assert( !flags_.empty() );
    assert( nth_value * flag_stride_ + ith_dim < flags_.size() );
    assert( ( type_ == VECTOR       && ith_dim < dim_         ) ||
            ( type_ == ARRAY        && ith_dim == 0U          ) ||
            ( type_ == FLAGGEDARRAY && ith_dim < flag_stride_ ) );
    return flags_[ nth_value * flag_stride_ + ith_dim ];
}

inline VARIABLE_FLAG PropertyData::Flag( size_t nth_value, uint32_t ith_dim ) const noexcept
{
    assert( type_ == VECTOR || type_ == ARRAY || type_ == FLAGGEDARRAY );
    assert( !flags_.empty() );
    assert( nth_value * flag_stride_ + ith_dim < flags_.size() );
    assert( ( type_ == VECTOR       && ith_dim < dim_         ) ||
            ( type_ == ARRAY        && ith_dim == 0U          ) ||
            ( type_ == FLAGGEDARRAY && ith_dim < flag_stride_ ) );
    return flags_[ nth_value * flag_stride_ + ith_dim ];
}

inline VARIABLE_FLAG& PropertyData::Flag( size_t   nth_value,
                                           uint32_t ith_row,
                                           uint32_t jth_col ) noexcept
{
    assert( type_ == TENSOR );
    assert( !flags_.empty() );
    assert( nth_value * flag_stride_ + jth_col < flags_.size() );
    assert( ith_row < dim_ );
    assert( jth_col < dim_ );
    assert( jth_col == ith_row ); // only diagonal elements carry flags
    return flags_[ nth_value * flag_stride_ + jth_col ];
}

inline VARIABLE_FLAG PropertyData::Flag( size_t   nth_value,
                                          uint32_t ith_row,
                                          uint32_t jth_col ) const noexcept
{
    assert( type_ == TENSOR );
    assert( !flags_.empty() );
    assert( nth_value * flag_stride_ + jth_col < flags_.size() );
    assert( ith_row < dim_ );
    assert( jth_col < dim_ );
    assert( jth_col == ith_row ); // only diagonal elements carry flags
    return flags_[ nth_value * flag_stride_ + jth_col ];
}


// ============================================================================
// Inline definitions — value accessors / mutators
// ============================================================================

inline double& PropertyData::Value( size_t nth_value ) noexcept
{
    assert( type_ == SCALAR );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    return data_[ nth_value ];
}

inline double PropertyData::Value( size_t nth_value ) const noexcept
{
    assert( type_ == SCALAR );
    assert( !data_.empty() );
    assert( nth_value < data_.size() );
    return data_[ nth_value ];
}

inline double& PropertyData::Value( size_t nth_value, uint32_t ith_dim ) noexcept
{
    assert( type_ == VECTOR || type_ == ARRAY || type_ == FLAGGEDARRAY );
    assert( !data_.empty() );
    assert( nth_value * data_stride_ + ith_dim < data_.size() );
    assert( ( type_ == VECTOR       && ith_dim < dim_         ) ||
            ( type_ == ARRAY        && ith_dim < data_stride_ ) ||
            ( type_ == FLAGGEDARRAY && ith_dim < data_stride_ ) );
    return data_[ nth_value * data_stride_ + ith_dim ];
}

inline double PropertyData::Value( size_t nth_value, uint32_t ith_dim ) const noexcept
{
    assert( type_ == VECTOR || type_ == ARRAY || type_ == FLAGGEDARRAY );
    assert( !data_.empty() );
    assert( nth_value * data_stride_ + ith_dim < data_.size() );
    assert( ( type_ == VECTOR       && ith_dim < dim_         ) ||
            ( type_ == ARRAY        && ith_dim < data_stride_ ) ||
            ( type_ == FLAGGEDARRAY && ith_dim < data_stride_ ) );
    return data_[ nth_value * data_stride_ + ith_dim ];
}

inline double& PropertyData::Value( size_t   nth_value,
                                     uint32_t ith_row,
                                     uint32_t jth_col ) noexcept
{
    assert( type_ == TENSOR );
    assert( !data_.empty() );
    assert( nth_value * data_stride_ + ith_row * dim_ + jth_col < data_.size() );
    assert( ith_row < dim_ );
    assert( jth_col < dim_ );
    return data_[ nth_value * data_stride_ + ith_row * dim_ + jth_col ];
}

inline double PropertyData::Value( size_t   nth_value,
                                    uint32_t ith_row,
                                    uint32_t jth_col ) const noexcept
{
    assert( type_ == TENSOR );
    assert( !data_.empty() );
    assert( nth_value * data_stride_ + ith_row * dim_ + jth_col < data_.size() );
    assert( ith_row < dim_ );
    assert( jth_col < dim_ );
    return data_[ nth_value * data_stride_ + ith_row * dim_ + jth_col ];
}


// ============================================================================
// Inline definitions — MinMaxOf, OffsetRangeBy, TransformValues
// ============================================================================

inline void PropertyData::MinMaxOf( double& tmin, double& tmax ) const noexcept
{
    tmin = *std::min_element( data_.begin(), data_.end() );
    tmax = *std::max_element( data_.begin(), data_.end() );
}

inline void PropertyData::OffsetRangeBy( double offset ) noexcept
{
    std::transform( data_.begin(), data_.end(), data_.begin(),
                    [offset]( double v ) noexcept { return v + offset; } );
}

inline void PropertyData::TransformValues( double (*f)(double) ) noexcept
{
    std::transform( data_.begin(), data_.end(), data_.begin(), f );
}


// ============================================================================
// Inline definitions — pushBack free functions
// ============================================================================

inline void pushBack( PropertyData& data, const ScalarVariable& sc )
{
    assert( data.Type() == SCALAR );
    data.PushBack( sc.Flag() );
    data.PushBack( sc() );
}

inline void pushBack( PropertyData& data, const ArrayVariable& ar )
{
    assert( data.Type() == ARRAY );
    data.PushBack( ar.Flag() );
    for ( uint32_t i = 0U; i < ar.Size(); ++i )
        data.PushBack( ar[i] );
}

inline void pushBack( PropertyData& data, const FlaggedArrayVariable& fa )
{
    assert( data.Type() == FLAGGEDARRAY );
    // Note: PushBack( VARIABLE_FLAG ) and PushBack( double ) append to
    // separate internal vectors (flags_ and data_ respectively), so the
    // interleaved call order here is intentional and correct.
    for ( uint32_t i = 0U; i < fa.Size(); ++i ) {
        data.PushBack( fa.Flag(i) );
        data.PushBack( fa[i] );
    }
}

template<uint32_t dim>
inline void pushBack( PropertyData& data, const VectorVariable<dim>& vc )
{
    assert( data.Type() == VECTOR );
    assert( data.Dim()  == dim );
    for ( uint32_t i = 0U; i < dim; ++i ) {
        data.PushBack( vc.Flag(i) );
        data.PushBack( vc[i] );
    }
}

template<uint32_t dim>
inline void pushBack( PropertyData& data, const TensorVariable<dim>& ts )
{
    assert( data.Type() == TENSOR );
    assert( data.Dim()  == dim );
    for ( uint32_t i = 0U; i < dim; ++i ) {
        data.PushBack( ts.Flag(i) );
        for ( uint32_t j = 0U; j < dim; ++j )
            data.PushBack( ts(i,j) );
    }
}


// ============================================================================
// Inline definitions — store free functions
// ============================================================================

template<>
inline void store( PropertyData& data, size_t position, const ScalarVariable& sc )
{
    assert( data.Type() == SCALAR );
    data.Flag(  position ) = sc.Flag();
    data.Value( position ) = sc();
}

template<>
inline void store( PropertyData& data, size_t position, const ArrayVariable& ary )
{
    assert( data.Type() == ARRAY );
    data.Flag( position, 0U ) = ary.Flag();
    for ( uint32_t i = 0U; i < ary.Size(); ++i )
        data.Value( position, i ) = ary[i];
}

template<>
inline void store( PropertyData& data, size_t position, const FlaggedArrayVariable& ary )
{
    assert( data.Type() == FLAGGEDARRAY );
    for ( uint32_t i = 0U; i < ary.Size(); ++i ) {
        data.Flag(  position, i ) = ary.Flag(i);
        data.Value( position, i ) = ary[i];
    }
}

template<uint32_t dim>
inline void store( PropertyData& data, size_t position, const VectorVariable<dim>& vc )
{
    assert( data.Type() == VECTOR );
    assert( data.Dim()  == dim );
    for ( uint32_t i = 0U; i < dim; ++i ) {
        data.Flag(  position, i ) = vc.Flag(i);
        data.Value( position, i ) = vc[i];
    }
}

template<uint32_t dim>
inline void store( PropertyData& data, size_t position, const TensorVariable<dim>& ts )
{
    assert( data.Type() == TENSOR );
    assert( data.Dim()  == dim );
    for ( uint32_t i = 0U; i < dim; ++i ) {
        data.Flag( position, i, i ) = ts.Flag(i);
        for ( uint32_t j = 0U; j < dim; ++j )
            data.Value( position, i, j ) = ts(i,j);
    }
}


// ============================================================================
// Inline definitions — read free functions
// ============================================================================

template<>
inline void read( const PropertyData& data, size_t position, ScalarVariable& sc )
{
    assert( data.Type() == SCALAR );
    sc.Flag() = data.Flag(  position );
    sc()      = data.Value( position );
}

template<>
inline void read( const PropertyData& data, size_t position, ArrayVariable& ary )
{
    assert( data.Type() == ARRAY );
    ary.Flag() = data.Flag( position, 0U );
    const uint32_t array_size = data.Components();
    ary.Resize( array_size );
    for ( uint32_t i = 0U; i < array_size; ++i )
        ary(i) = data.Value( position, i );
}

template<>
inline void read( const PropertyData& data, size_t position, FlaggedArrayVariable& ary )
{
    assert( data.Type() == FLAGGEDARRAY );
    const uint32_t array_size = data.Components();
    ary.Resize( array_size );
    for ( uint32_t i = 0U; i < array_size; ++i ) {
        ary.Flag(i) = data.Flag(  position, i );
        ary(i)      = data.Value( position, i );
    }
}

template<uint32_t dim>
inline void read( const PropertyData& data, size_t position, VectorVariable<dim>& vc )
{
    assert( data.Type() == VECTOR );
    assert( data.Dim()  == dim );
    for ( uint32_t i = 0U; i < dim; ++i ) {
        vc.Flag(i) = data.Flag(  position, i );
        vc(i)      = data.Value( position, i );
    }
}

template<uint32_t dim>
inline void read( const PropertyData& data, size_t position, TensorVariable<dim>& ts )
{
    assert( data.Type() == TENSOR );
    assert( data.Dim()  == dim );
    for ( uint32_t i = 0U; i < dim; ++i ) {
        ts.Flag(i) = data.Flag( position, i, i );
        for ( uint32_t j = 0U; j < dim; ++j )
            ts(i,j) = data.Value( position, i, j );
    }
}


// ============================================================================
// Free function — read PropertyData record from binary file
// ============================================================================

/// Read a PropertyData record written by PropertyData::OutBinary() and return
/// an initialised container.
PropertyData inBinaryPropertyData( std::fstream& fp );

} // namespace csmp

#endif // CSMP_PROPERTY_DATA_H

