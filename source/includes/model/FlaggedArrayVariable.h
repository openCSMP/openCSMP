// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FLAGGED_ARRAY_VARIABLE_H
#define FLAGGED_ARRAY_VARIABLE_H

#include "CSMP_definitions.h"
#include "ScalarVariable.h"
#include "ArrayVariable.h"

namespace csmp {

class ScalarVariable;
class ArrayVariable;
template<uint32_t> class PropertyDatabase;

/**
@brief Storage for an arbitrary number of variables each with an
       independent flag; array size is determined at runtime.

@author Roman Manasipov
@author Philipp Lang
@date 2013
@author (refactored) 2024

@section motivation Motivation

Extends ArrayVariable by giving each element its own VARIABLE_FLAG.
This models chemical systems where some species concentrations are
fixed (e.g. DIRICH) and others are free (ANY).

@section usage Usage

@code
FlaggedArrayVariable fav;                   // 0 entries
FlaggedArrayVariable fav( 4, 10., DIRICH ); // 4 entries = 10, all DIRICH
FlaggedArrayVariable fav( "element array", model.Database() );
@endcode

@section flags Flag semantics

Each element has its own VARIABLE_FLAG. Operations in PropertyHandle3
that respect OutputCondition() only modify elements whose flag matches
the output condition. Use ApplyToFlaggedElements() for element-wise
flag-aware operations.

@section comparison Comparison operators

All six comparison operators (==, !=, <, <=, >, >=) compare values
element-by-element. Two arrays are equal if all corresponding values
are equal (flags are not compared). Ordering operators return true
only if the relation holds for ALL elements simultaneously.

@note The binary I/O methods Out(fstream&) and In(fstream&) are for
internal CSMP binary serialisation only. Do not use them directly.
*/
class FlaggedArrayVariable
{
public:
    static constexpr VARIABLE_TYPE VariableType = FLAGGEDARRAY;

    typedef std::vector<double> FlaggedArrayContainer;

    // -----------------------------------------------------------------------
    // Construction / destruction
    // -----------------------------------------------------------------------

    /** Default constructor: zero elements. */
    FlaggedArrayVariable() noexcept;

    /** Copy constructor. */
    FlaggedArrayVariable( const FlaggedArrayVariable& ) noexcept;

    /** Move constructor. */
    FlaggedArrayVariable( FlaggedArrayVariable&& other ) noexcept;

    /**
    Constructs with size and default values read from the
    PropertyDatabase entry for @p arrayPropertyName.
    */
    template<uint32_t dim>
    FlaggedArrayVariable( const char*                  arrayPropertyName,
                          const PropertyDatabase<dim>& pd,
                          double                       defaultValue =
                              std::numeric_limits<double>::quiet_NaN(),
                          VARIABLE_FLAG                flag = ANY );

    /**
    Constructs with @p arraySize elements, all set to @p defaultValue
    with flag @p flag.
    */
    FlaggedArrayVariable( uint32_t      arraySize,
                          double        defaultValue = 0.,
                          VARIABLE_FLAG flag = ANY ) noexcept;

    /**
    Constructs with size derived from @p arrayKey.dataDepth, all
    elements set to @p defaultValue with flag @p flag.
    */
    FlaggedArrayVariable( const Index&  arrayKey,
                          double        defaultValue = 0.,
                          VARIABLE_FLAG flag = ANY ) noexcept;

    // -----------------------------------------------------------------------
    // Assignment
    // -----------------------------------------------------------------------

    FlaggedArrayVariable& operator=( const FlaggedArrayVariable& ) noexcept;
    FlaggedArrayVariable& operator=( FlaggedArrayVariable&& other ) noexcept;

    /**
    Copies values from @p av into this variable's data, ignoring flags.
    Sizes must match.
    */
    void CopyValuesOnly( ArrayVariable& av ) noexcept;

    /** Sets all data elements to @p val. Flags are not modified. */
    FlaggedArrayVariable& operator=( double val ) noexcept;

    /** Sets all data elements to the scalar value. Flags are not modified. */
    FlaggedArrayVariable& operator=( const ScalarVariable& ) noexcept;

    // -----------------------------------------------------------------------
    // Arithmetic — returning new object (flags copied from *this)
    // -----------------------------------------------------------------------

    FlaggedArrayVariable operator+( double ) const noexcept;
    FlaggedArrayVariable operator-( double ) const noexcept;
    FlaggedArrayVariable operator*( double ) const noexcept;
    FlaggedArrayVariable operator/( double ) const noexcept;

    FlaggedArrayVariable operator+( const FlaggedArrayVariable& ) const noexcept;
    FlaggedArrayVariable operator-( const FlaggedArrayVariable& ) const noexcept;
    FlaggedArrayVariable operator*( const FlaggedArrayVariable& ) const noexcept;
    FlaggedArrayVariable operator/( const FlaggedArrayVariable& ) const noexcept;

    /**
    Returns a new FlaggedArrayVariable with each element raised to
    @p exponent. Flags are copied from *this.
    */
    FlaggedArrayVariable Pow( double exponent ) const noexcept;

    // -----------------------------------------------------------------------
    // Compound assignment — modifies *this in place (flags not modified)
    // -----------------------------------------------------------------------

    FlaggedArrayVariable& operator+=( double ) noexcept;
    FlaggedArrayVariable& operator-=( double ) noexcept;
    FlaggedArrayVariable& operator*=( double ) noexcept;
    FlaggedArrayVariable& operator/=( double ) noexcept;

    FlaggedArrayVariable& operator+=( const ScalarVariable& ) noexcept;
    FlaggedArrayVariable& operator-=( const ScalarVariable& ) noexcept;
    FlaggedArrayVariable& operator*=( const ScalarVariable& ) noexcept;
    FlaggedArrayVariable& operator/=( const ScalarVariable& ) noexcept;

    FlaggedArrayVariable& operator+=( const FlaggedArrayVariable& ) noexcept;
    FlaggedArrayVariable& operator-=( const FlaggedArrayVariable& ) noexcept;
    FlaggedArrayVariable& operator*=( const FlaggedArrayVariable& ) noexcept;
    FlaggedArrayVariable& operator/=( const FlaggedArrayVariable& ) noexcept;

    // -----------------------------------------------------------------------
    // Comparison (value-only; flags are not compared)
    // -----------------------------------------------------------------------

    /** True if all corresponding values are equal. */
    bool operator==( const FlaggedArrayVariable& ) const noexcept;
    bool operator!=( const FlaggedArrayVariable& ) const noexcept;

    /**
    Lexicographic ordering on values (flags ignored).
    Satisfies strict weak ordering for use in std::set / std::map.
    */
    bool operator< ( const FlaggedArrayVariable& ) const noexcept;
    bool operator<=( const FlaggedArrayVariable& ) const noexcept;
    bool operator> ( const FlaggedArrayVariable& ) const noexcept;
    bool operator>=( const FlaggedArrayVariable& ) const noexcept;

    // -----------------------------------------------------------------------
    // Element access
    // -----------------------------------------------------------------------

    double         operator[]( size_t i ) const noexcept;
    double&        operator()( size_t i )       noexcept;

    void           Component( uint32_t i, double val ) noexcept;
    double         Component( uint32_t i )       const noexcept;

    // -----------------------------------------------------------------------
    // Size and flags
    // -----------------------------------------------------------------------

    uint32_t       Size()   const noexcept;
    void           Resize( uint32_t newSize,
                           double   newValue =
                               std::numeric_limits<double>::quiet_NaN() ) noexcept;

    VARIABLE_FLAG  Flag( size_t i ) const noexcept;
    VARIABLE_FLAG& Flag( size_t i )       noexcept;
    void           Flag( size_t i, VARIABLE_FLAG flag ) noexcept;

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /**
    Returns true if all element values lie in [@p min, @p max].
    Flags are not considered.
    */
    bool IsWithinRange( double min, double max ) const noexcept;

    /** Returns true if any element value is NaN. */
    bool Has_NaN_Values() const noexcept;

    /**
    Returns the minimum and maximum element values into @p min and
    @p max. Handles empty arrays by returning NaN for both.
    */
    void MinMax( double& min, double& max ) const noexcept;

    /** Sorts data values in ascending order. Flags follow their values. */
    void Sort() noexcept;

    /**
    Returns the smallest value strictly greater than @p fromValue.
    Returns numeric_limits<double>::max() if no such value exists.
    O(n) linear scan — no allocation, no sort.
    */
    double NextLargestEntry( double fromValue ) const noexcept;

    /**
    Returns true if any element value is strictly greater than
    @p fromValue. O(n) worst case, O(1) best case.
    */
    bool HasLargerEntry( double fromValue ) const noexcept;

    // -----------------------------------------------------------------------
    // Standard range interface
    // -----------------------------------------------------------------------

    FlaggedArrayContainer::iterator       begin()        noexcept;
    FlaggedArrayContainer::iterator       end()          noexcept;
    FlaggedArrayContainer::const_iterator begin()  const noexcept;
    FlaggedArrayContainer::const_iterator end()    const noexcept;
    FlaggedArrayContainer::const_iterator cbegin() const noexcept;
    FlaggedArrayContainer::const_iterator cend()   const noexcept;

    /** Legacy uppercase iterators — prefer begin()/end(). */
    FlaggedArrayContainer::const_iterator Begin() const noexcept;
    FlaggedArrayContainer::const_iterator End()   const noexcept;

    // -----------------------------------------------------------------------
    // Output
    // -----------------------------------------------------------------------

    void Out( long digits = 3 ) const;
    bool Out( const char* filename, size_t precision = 9 ) const;

    /** Internal CSMP binary I/O — do not use directly. */
    bool Out( std::fstream& fp ) const;
    bool In ( std::fstream& fp );

private:
    std::vector<VARIABLE_FLAG>  flags_;
    FlaggedArrayContainer       data_;
};

std::ostream& operator<<( std::ostream& stream,
                           const FlaggedArrayVariable& o );


// ============================================================================
//  Inline definitions
// ============================================================================

inline FlaggedArrayVariable::FlaggedArrayVariable() noexcept
    : flags_(), data_()
{}

inline FlaggedArrayVariable::FlaggedArrayVariable(
    const FlaggedArrayVariable& other ) noexcept
    : flags_( other.flags_ ), data_( other.data_ )
{}

inline FlaggedArrayVariable::FlaggedArrayVariable(
    FlaggedArrayVariable&& other ) noexcept
    : flags_( std::move( other.flags_ ) ),
      data_(  std::move( other.data_  ) )
{}

inline FlaggedArrayVariable::FlaggedArrayVariable(
    uint32_t arraySize, double defaultValue, VARIABLE_FLAG flag ) noexcept
    : flags_( arraySize, flag ),
      data_(  arraySize, defaultValue )
{}

inline FlaggedArrayVariable::FlaggedArrayVariable(
    const Index& arrayKey, double defaultValue, VARIABLE_FLAG flag ) noexcept
    : flags_( static_cast<size_t>( arrayKey.dataDepth ), flag ),
      data_(  static_cast<size_t>( arrayKey.dataDepth ), defaultValue )
{}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator=(
    const FlaggedArrayVariable& other ) noexcept
{
    if ( this != &other )
    {
        data_  = other.data_;
        flags_ = other.flags_;
    }
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator=(
    FlaggedArrayVariable&& other ) noexcept
{
    if ( this != &other )
    {
        flags_ = std::move( other.flags_ );
        data_  = std::move( other.data_  );
    }
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator=(
    double val ) noexcept
{
    for ( auto& v : data_ ) v = val;
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator=(
    const ScalarVariable& sc ) noexcept
{
    for ( auto& v : data_ ) v = sc();
    return *this;
}

inline double FlaggedArrayVariable::operator[]( size_t i ) const noexcept
{
    assert( i < Size() );
    return data_[i];
}

inline double& FlaggedArrayVariable::operator()( size_t i ) noexcept
{
    assert( i < Size() );
    return data_[i];
}

inline void FlaggedArrayVariable::Component( uint32_t i, double val ) noexcept
{
    assert( i < Size() );
    data_[i] = val;
}

inline double FlaggedArrayVariable::Component( uint32_t i ) const noexcept
{
    assert( i < Size() );
    return data_[i];
}

inline uint32_t FlaggedArrayVariable::Size() const noexcept
{
    return static_cast<uint32_t>( data_.size() );
}

inline void FlaggedArrayVariable::Resize(
    uint32_t newSize, double newValue ) noexcept
{
    data_.resize ( newSize, newValue );
    flags_.resize( newSize, ANY );
}

inline VARIABLE_FLAG FlaggedArrayVariable::Flag( size_t i ) const noexcept
{
    return flags_[i];
}

inline VARIABLE_FLAG& FlaggedArrayVariable::Flag( size_t i ) noexcept
{
    return flags_[i];
}

inline void FlaggedArrayVariable::Flag( size_t i, VARIABLE_FLAG flag ) noexcept
{
    flags_[i] = flag;
}

inline bool FlaggedArrayVariable::IsWithinRange(
    double min, double max ) const noexcept
{
    if ( min >= max ) return false;
    for ( const double v : data_ )
        if ( v < min || v > max ) return false;
    return true;
}

inline bool FlaggedArrayVariable::Has_NaN_Values() const noexcept
{
    return std::any_of( data_.cbegin(), data_.cend(),
                        []( double v ) { return std::isnan(v); } );
}

inline void FlaggedArrayVariable::MinMax(
    double& min, double& max ) const noexcept
{
    if ( data_.empty() )
    {
        min = std::numeric_limits<double>::quiet_NaN();
        max = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const auto [minIt, maxIt] =
        std::minmax_element( data_.cbegin(), data_.cend() );
    min = *minIt;
    max = *maxIt;
}

inline double FlaggedArrayVariable::NextLargestEntry(
    double fromValue ) const noexcept
{
    double result = std::numeric_limits<double>::max();
    for ( const double v : data_ )
        if ( v > fromValue && v < result )
            result = v;
    return result;
}

inline bool FlaggedArrayVariable::HasLargerEntry(
    double fromValue ) const noexcept
{
    return std::any_of( data_.cbegin(), data_.cend(),
                        [fromValue]( double v ) { return v > fromValue; } );
}

inline void FlaggedArrayVariable::Sort() noexcept
{
    // Sort data and flags together by data value.
    // Build index array, sort by data, then reorder both.
    const uint32_t n = Size();
    std::vector<uint32_t> idx( n );
    std::iota( idx.begin(), idx.end(), 0U );
    std::sort( idx.begin(), idx.end(),
               [this]( uint32_t a, uint32_t b )
               { return data_[a] < data_[b]; } );

    FlaggedArrayContainer       sorted_data ( n );
    std::vector<VARIABLE_FLAG>  sorted_flags( n );
    for ( uint32_t i = 0; i < n; ++i )
    {
        sorted_data [i] = data_ [idx[i]];
        sorted_flags[i] = flags_[idx[i]];
    }
    data_  = std::move( sorted_data  );
    flags_ = std::move( sorted_flags );
}

inline FlaggedArrayVariable::FlaggedArrayContainer::iterator
FlaggedArrayVariable::begin() noexcept
{ return data_.begin(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::iterator
FlaggedArrayVariable::end() noexcept
{ return data_.end(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::const_iterator
FlaggedArrayVariable::begin() const noexcept
{ return data_.begin(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::const_iterator
FlaggedArrayVariable::end() const noexcept
{ return data_.end(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::const_iterator
FlaggedArrayVariable::cbegin() const noexcept
{ return data_.cbegin(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::const_iterator
FlaggedArrayVariable::cend() const noexcept
{ return data_.cend(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::const_iterator
FlaggedArrayVariable::Begin() const noexcept
{ return data_.begin(); }

inline FlaggedArrayVariable::FlaggedArrayContainer::const_iterator
FlaggedArrayVariable::End() const noexcept
{ return data_.end(); }

// --- Arithmetic returning new object ---

inline FlaggedArrayVariable FlaggedArrayVariable::operator+(
    double val ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( auto& v : result.data_ ) v += val;
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator-(
    double val ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( auto& v : result.data_ ) v -= val;
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator*(
    double val ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( auto& v : result.data_ ) v *= val;
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator/(
    double val ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( auto& v : result.data_ ) v /= val;
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::Pow(
    double exponent ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( uint32_t i = 0; i < Size(); ++i )
        result.data_[i] = std::pow( data_[i], exponent );
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator+(
    const FlaggedArrayVariable& other ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( uint32_t i = 0; i < Size(); ++i )
        result.data_[i] += other.data_[i];
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator-(
    const FlaggedArrayVariable& other ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( uint32_t i = 0; i < Size(); ++i )
        result.data_[i] -= other.data_[i];
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator*(
    const FlaggedArrayVariable& other ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( uint32_t i = 0; i < Size(); ++i )
        result.data_[i] *= other.data_[i];
    return result;
}

inline FlaggedArrayVariable FlaggedArrayVariable::operator/(
    const FlaggedArrayVariable& other ) const noexcept
{
    FlaggedArrayVariable result( *this );
    for ( uint32_t i = 0; i < Size(); ++i )
        result.data_[i] /= other.data_[i];
    return result;
}

// --- Compound assignment ---

inline FlaggedArrayVariable& FlaggedArrayVariable::operator+=(
    double val ) noexcept
{
    for ( auto& v : data_ ) v += val;
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator-=(
    double val ) noexcept
{
    for ( auto& v : data_ ) v -= val;
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator*=(
    double val ) noexcept
{
    for ( auto& v : data_ ) v *= val;
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator/=(
    double val ) noexcept
{
    for ( auto& v : data_ ) v /= val;
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator+=(
    const ScalarVariable& sc ) noexcept
{
    for ( auto& v : data_ ) v += sc();
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator-=(
    const ScalarVariable& sc ) noexcept
{
    for ( auto& v : data_ ) v -= sc();
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator*=(
    const ScalarVariable& sc ) noexcept
{
    for ( auto& v : data_ ) v *= sc();
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator/=(
    const ScalarVariable& sc ) noexcept
{
    for ( auto& v : data_ ) v /= sc();
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator+=(
    const FlaggedArrayVariable& other ) noexcept
{
    for ( uint32_t i = 0; i < Size(); ++i ) data_[i] += other.data_[i];
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator-=(
    const FlaggedArrayVariable& other ) noexcept
{
    for ( uint32_t i = 0; i < Size(); ++i ) data_[i] -= other.data_[i];
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator*=(
    const FlaggedArrayVariable& other ) noexcept
{
    for ( uint32_t i = 0; i < Size(); ++i ) data_[i] *= other.data_[i];
    return *this;
}

inline FlaggedArrayVariable& FlaggedArrayVariable::operator/=(
    const FlaggedArrayVariable& other ) noexcept
{
    for ( uint32_t i = 0; i < Size(); ++i ) data_[i] /= other.data_[i];
    return *this;
}

// --- Comparison ---

inline bool FlaggedArrayVariable::operator==(
    const FlaggedArrayVariable& other ) const noexcept
{
    return data_ == other.data_;
}

inline bool FlaggedArrayVariable::operator!=(
    const FlaggedArrayVariable& other ) const noexcept
{
    return data_ != other.data_;
}

inline bool FlaggedArrayVariable::operator<(
    const FlaggedArrayVariable& other ) const noexcept
{
    return data_ < other.data_;
}

inline bool FlaggedArrayVariable::operator<=(
    const FlaggedArrayVariable& other ) const noexcept
{
    return data_ <= other.data_;
}

inline bool FlaggedArrayVariable::operator>(
    const FlaggedArrayVariable& other ) const noexcept
{
    return data_ > other.data_;
}

inline bool FlaggedArrayVariable::operator>=(
    const FlaggedArrayVariable& other ) const noexcept
{
    return data_ >= other.data_;
}

} // namespace csmp

#endif // FLAGGED_ARRAY_VARIABLE_H
