// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_VECTOR_VARIABLE_H
#define CSMP_VECTOR_VARIABLE_H

#include "Point.h"
#include "ScalarVariable.h"
#include "VectorVariable1.h"
#include "VectorVariable2.h"

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace csmp {

/**
@brief 3D full specialisation of the VectorVariable class template.

@author S.K. Matthaei
@author Stephen G. Roberts
@author S. Geiger
@date 2001
@author (refactored) 2024

@section motivation Motivation

A plain STL vector cannot be used to implement this datatype since
each vector entry must carry a VARIABLE_FLAG for finite element
treatment (e.g. Dirichlet boundary conditions). All loops are
unrolled for performance.

@section operators Special operators

  operator&  — dot product (scalar product)
  operator%  — cross product
  Pow(e)     — raises each element to the power e
               (replaces the old operator^ which was ambiguous
               due to C++ operator precedence)

@note operator^ has been removed. Use Pow(double) instead.
*/
template<>
class VectorVariable<3U>
{
public:
    static constexpr VARIABLE_TYPE VariableType = VECTOR;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /** Default constructor: all elements NaN, all flags ANY. */
    VectorVariable() noexcept;

    /** Sets all elements to @p val with flag @p f. */
    VectorVariable( VARIABLE_FLAG f, double val ) noexcept;

    /** Full initialisation with independent flags and values. */
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3,
                    double val1, double val2, double val3 ) noexcept;

    /** Initialises from an STL vector. All flags set to ANY. */
    explicit VectorVariable( const std::vector<double>& v ) noexcept;

    /** Initialises from a Point. All flags set to ANY. */
    explicit VectorVariable( const csmp::Point<3U>& p ) noexcept;

    // -----------------------------------------------------------------------
    // Element access
    // -----------------------------------------------------------------------

    /** Read/write access to element @p i. */
    double&       operator()( uint32_t i )       noexcept;

    /** Read-only access to element @p i. */
    const double& operator()( uint32_t i ) const noexcept;

    /** Read-only access to element @p i (no flag returned). */
    double        operator[]( uint32_t i ) const noexcept;

    /** Alternative mutator. */
    void   Component( uint32_t i, double val ) noexcept;

    /** Alternative accessor. */
    [[nodiscard]] double Component( uint32_t i ) const noexcept;

    /** Number of elements (always 3). */
    static constexpr uint32_t Size() noexcept { return 3U; }

    // -----------------------------------------------------------------------
    // Assignment
    // -----------------------------------------------------------------------

    /** Sets all elements to @p val. Flags are not modified. */
    VectorVariable& operator=( double val )                noexcept;

    /** Sets all elements to the scalar value and all flags to the scalar's flag. */
    VectorVariable& operator=( const ScalarVariable& sc )  noexcept;

    /** Sets all elements from the point. Flags are not modified. */
    VectorVariable& operator=( const csmp::Point<3U>& pt ) noexcept;

    // -----------------------------------------------------------------------
    // Arithmetic — returning new vector (flags from *this)
    // -----------------------------------------------------------------------

    VectorVariable operator+( const VectorVariable& ) const noexcept;
    VectorVariable operator-( const VectorVariable& ) const noexcept;

    /** Element-by-element multiplication. */
    VectorVariable operator*( const VectorVariable& ) const noexcept;

    /** Element-by-element division. */
    VectorVariable operator/( const VectorVariable& ) const noexcept;

    /**
    Raises each element to the power @p exponent.
    Replaces the old operator^ which had ambiguous precedence.
    */
    VectorVariable Pow( double exponent ) const noexcept;

    // -----------------------------------------------------------------------
    // Compound assignment
    // -----------------------------------------------------------------------

    VectorVariable& operator+=( double val )            noexcept;
    VectorVariable& operator-=( double val )            noexcept;
    VectorVariable& operator*=( double val )            noexcept;
    VectorVariable& operator/=( double val )            noexcept;

    VectorVariable& operator+=( const ScalarVariable& ) noexcept;
    VectorVariable& operator-=( const ScalarVariable& ) noexcept;
    VectorVariable& operator*=( const ScalarVariable& ) noexcept;
    VectorVariable& operator/=( const ScalarVariable& ) noexcept;

    VectorVariable& operator+=( const VectorVariable& ) noexcept;
    VectorVariable& operator-=( const VectorVariable& ) noexcept;

    /** Element-by-element multiplication in place. */
    VectorVariable& operator*=( const VectorVariable& ) noexcept;

    /** Element-by-element division in place. */
    VectorVariable& operator/=( const VectorVariable& ) noexcept;

    // -----------------------------------------------------------------------
    // Comparison
    // -----------------------------------------------------------------------

    /** Element-by-element equality of values and flags. */
    bool operator==( const VectorVariable& ) const noexcept;

    /** Element-by-element inequality. */
    bool operator!=( const VectorVariable& ) const noexcept;

    /**
    Ordering by pointer address — satisfies strict weak ordering for
    use in STL associative containers. Does not compare values.
    */
    bool operator<( const VectorVariable& ) const noexcept;

    // -----------------------------------------------------------------------
    // Vector operations
    // -----------------------------------------------------------------------

    /** Dot product (scalar product). */
    double operator&( const VectorVariable& ) const noexcept;

    /** Cross product. */
    VectorVariable operator%( const VectorVariable& ) const noexcept;

    /** Euclidean length. */
    double Length() const noexcept;

    /**
    Returns the angle in degrees between this vector and @p v.
    Returns 90 if the denominator is zero.
    */
    double AngleTo( const VectorVariable& v ) const noexcept;

    /** Returns a csmp::Point initialised with the vector values. */
    Point<3U> P() const noexcept;

    /** Dot product with a Point. */
    double DotProduct( const csmp::Point<3U>& ) const noexcept;

    /** Dot product with a VectorVariable. */
    double DotProduct( const VectorVariable& ) const noexcept;

    /** Cross product with a Point. */
    VectorVariable CrossProduct( const csmp::Point<3U>& ) const noexcept;

    /** Cross product with a VectorVariable. */
    VectorVariable CrossProduct( const VectorVariable& ) const noexcept;

    /** Projects this vector onto @p v (given as std::vector). */
    VectorVariable ProjectOnto( const std::vector<double>& v ) const noexcept;

    /** Projects this vector onto @p v. */
    VectorVariable ProjectOnto( const VectorVariable& v ) const noexcept;

    // -----------------------------------------------------------------------
    // Modification
    // -----------------------------------------------------------------------

    /**
    Returns a new vector with elements in reversed order and flags
    swapped accordingly.
    */
    VectorVariable Flip() noexcept;

    /** Multiplies all elements by -1. */
    void Invert() noexcept;

    /** Normalises the vector to unit length. No-op if length is zero. */
    void EuclideanNormalize() noexcept;

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /**
    Returns true if the vector length lies within [@p vmin, @p vmax].
    */
    bool IsWithinRange( double vmin, double vmax ) const noexcept;

    /** Returns true if any element is NaN. */
    bool Has_NaN_Values() const noexcept;

    /** Returns the flag of element @p i. Const version. */
    VARIABLE_FLAG  Flag( uint32_t i = 0 ) const noexcept;

    /** Returns a reference to the flag of element @p i. */
    VARIABLE_FLAG& Flag( uint32_t i = 0 )       noexcept;

    // -----------------------------------------------------------------------
    // I/O
    // -----------------------------------------------------------------------

    /** Prompts the user to enter vector values from the command line. */
    void In();

    /** Prints the vector to stdout. */
    void Out() const noexcept;

    /** Reads the vector from a binary file stream. */
    bool In(  std::fstream& fp );

    /** Writes the vector to a binary file stream. */
    bool Out( std::fstream& fp ) const;

    friend class TensorVariable<3U>;

private:
    std::array<VARIABLE_FLAG, 3U> flag;
    std::array<double, 3U>        data;
};


// ============================================================================
//  Factory functions
// ============================================================================

VectorVariable<1U> makeVector( VARIABLE_FLAG,
                               double ) noexcept;

VectorVariable<2U> makeVector( VARIABLE_FLAG, VARIABLE_FLAG,
                               double, double ) noexcept;

VectorVariable<3U> makeVector( VARIABLE_FLAG, VARIABLE_FLAG, VARIABLE_FLAG,
                               double, double, double ) noexcept;

VectorVariable<3U> makeVector( const std::array<VARIABLE_FLAG,3U>&,
                               const std::array<double,3U>& ) noexcept;

VectorVariable<3U> makeVector( const std::vector<VARIABLE_FLAG>&,
                               const std::vector<double>& ) noexcept;


// ============================================================================
//  Stream output
// ============================================================================

template<uint32_t dim>
std::ostream& operator<<( std::ostream&, const VectorVariable<dim>& );


// ============================================================================
//  Binary operators — VectorVariable op ScalarVariable
// ============================================================================

VectorVariable<1U> operator+( const VectorVariable<1U>&, const ScalarVariable& ) noexcept;
VectorVariable<1U> operator-( const VectorVariable<1U>&, const ScalarVariable& ) noexcept;
VectorVariable<1U> operator*( const VectorVariable<1U>&, const ScalarVariable& ) noexcept;
VectorVariable<1U> operator/( const VectorVariable<1U>&, const ScalarVariable& ) noexcept;

VectorVariable<2U> operator+( const VectorVariable<2U>&, const ScalarVariable& ) noexcept;
VectorVariable<2U> operator-( const VectorVariable<2U>&, const ScalarVariable& ) noexcept;
VectorVariable<2U> operator*( const VectorVariable<2U>&, const ScalarVariable& ) noexcept;
VectorVariable<2U> operator/( const VectorVariable<2U>&, const ScalarVariable& ) noexcept;

VectorVariable<3U> operator+( const VectorVariable<3U>&, const ScalarVariable& ) noexcept;
VectorVariable<3U> operator-( const VectorVariable<3U>&, const ScalarVariable& ) noexcept;
VectorVariable<3U> operator*( const VectorVariable<3U>&, const ScalarVariable& ) noexcept;
VectorVariable<3U> operator/( const VectorVariable<3U>&, const ScalarVariable& ) noexcept;


// ============================================================================
//  Binary operators — VectorVariable op double
// ============================================================================

VectorVariable<1U> operator+( const VectorVariable<1U>&, double ) noexcept;
VectorVariable<1U> operator-( const VectorVariable<1U>&, double ) noexcept;
VectorVariable<1U> operator*( const VectorVariable<1U>&, double ) noexcept;
VectorVariable<1U> operator/( const VectorVariable<1U>&, double ) noexcept;

VectorVariable<2U> operator+( const VectorVariable<2U>&, double ) noexcept;
VectorVariable<2U> operator-( const VectorVariable<2U>&, double ) noexcept;
VectorVariable<2U> operator*( const VectorVariable<2U>&, double ) noexcept;
VectorVariable<2U> operator/( const VectorVariable<2U>&, double ) noexcept;

VectorVariable<3U> operator+( const VectorVariable<3U>&, double ) noexcept;
VectorVariable<3U> operator-( const VectorVariable<3U>&, double ) noexcept;
VectorVariable<3U> operator*( const VectorVariable<3U>&, double ) noexcept;
VectorVariable<3U> operator/( const VectorVariable<3U>&, double ) noexcept;


// ============================================================================
//  Binary operators — Point op VectorVariable
// ============================================================================

Point<1U> operator+( const Point<1U>&, const VectorVariable<1U>& ) noexcept;
Point<1U> operator-( const Point<1U>&, const VectorVariable<1U>& ) noexcept;
Point<1U> operator*( const Point<1U>&, const VectorVariable<1U>& ) noexcept;
Point<1U> operator/( const Point<1U>&, const VectorVariable<1U>& ) noexcept;

Point<2U> operator+( const Point<2U>&, const VectorVariable<2U>& ) noexcept;
Point<2U> operator-( const Point<2U>&, const VectorVariable<2U>& ) noexcept;
Point<2U> operator*( const Point<2U>&, const VectorVariable<2U>& ) noexcept;
Point<2U> operator/( const Point<2U>&, const VectorVariable<2U>& ) noexcept;

Point<3U> operator+( const Point<3U>&, const VectorVariable<3U>& ) noexcept;
Point<3U> operator-( const Point<3U>&, const VectorVariable<3U>& ) noexcept;
Point<3U> operator*( const Point<3U>&, const VectorVariable<3U>& ) noexcept;
Point<3U> operator/( const Point<3U>&, const VectorVariable<3U>& ) noexcept;


// ============================================================================
//  Inline definitions — VectorVariable<3U> members
// ============================================================================

inline VectorVariable<3U>::VectorVariable() noexcept
    : flag{ { ANY, ANY, ANY } },
      data{ { std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN() } }
{}

inline VectorVariable<3U>::VectorVariable(
    VARIABLE_FLAG f, double val ) noexcept
    : flag{ { f, f, f } },
      data{ { val, val, val } }
{}

inline VectorVariable<3U>::VectorVariable(
    VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3,
    double val1, double val2, double val3 ) noexcept
    : flag{ { f1, f2, f3 } },
      data{ { val1, val2, val3 } }
{}

inline VectorVariable<3U>::VectorVariable(
    const std::vector<double>& v ) noexcept
    : flag{ { ANY, ANY, ANY } },
      data{ { v[0], v[1], v[2] } }
{}

inline VectorVariable<3U>::VectorVariable(
    const csmp::Point<3U>& p ) noexcept
    : flag{ { ANY, ANY, ANY } },
      data{ { p[0], p[1], p[2] } }
{}

inline double& VectorVariable<3U>::operator()( uint32_t i ) noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nVectorVariable<3U>::operator(): violation i=" << i << "\n";
        return data[0];
    }
#endif
    return data[i];
}

inline const double& VectorVariable<3U>::operator()( uint32_t i ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nVectorVariable<3U>::operator() const: violation i=" << i << "\n";
        return data[0];
    }
#endif
    return data[i];
}

inline double VectorVariable<3U>::operator[]( uint32_t i ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nVectorVariable<3U>::operator[]: violation i=" << i << "\n";
        return data[0];
    }
#endif
    return data[i];
}

inline void VectorVariable<3U>::Component( uint32_t i, double val ) noexcept
{
    assert( i < 3U );
    data[i] = val;
}

inline double VectorVariable<3U>::Component( uint32_t i ) const noexcept
{
    assert( i < 3U );
    return data[i];
}

inline VARIABLE_FLAG& VectorVariable<3U>::Flag( uint32_t i ) noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nVectorVariable<3U>::Flag(): violation i=" << i << "\n";
        return flag[0];
    }
#endif
    return flag[i];
}

inline VARIABLE_FLAG VectorVariable<3U>::Flag( uint32_t i ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nVectorVariable<3U>::Flag() const: violation i=" << i << "\n";
        return flag[0];
    }
#endif
    return flag[i];
}

inline VectorVariable<3U>& VectorVariable<3U>::operator=( double val ) noexcept
{
    data[0] = data[1] = data[2] = val;
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator=(
    const ScalarVariable& sc ) noexcept
{
    flag[0] = flag[1] = flag[2] = sc.Flag();
    data[0] = data[1] = data[2] = sc();
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator=(
    const csmp::Point<3U>& pt ) noexcept
{
    data[0] = pt[0]; data[1] = pt[1]; data[2] = pt[2];
    return *this;
}

inline VectorVariable<3U> VectorVariable<3U>::operator+(
    const VectorVariable<3U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1], flag[2],
        data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2] );
}

inline VectorVariable<3U> VectorVariable<3U>::operator-(
    const VectorVariable<3U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1], flag[2],
        data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2] );
}

inline VectorVariable<3U> VectorVariable<3U>::operator*(
    const VectorVariable<3U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1], flag[2],
        data[0]*v.data[0], data[1]*v.data[1], data[2]*v.data[2] );
}

inline VectorVariable<3U> VectorVariable<3U>::operator/(
    const VectorVariable<3U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1], flag[2],
        data[0]/v.data[0], data[1]/v.data[1], data[2]/v.data[2] );
}

inline VectorVariable<3U> VectorVariable<3U>::Pow(
    double exponent ) const noexcept
{
    return VectorVariable( flag[0], flag[1], flag[2],
        std::pow( data[0], exponent ),
        std::pow( data[1], exponent ),
        std::pow( data[2], exponent ) );
}

inline VectorVariable<3U>& VectorVariable<3U>::operator+=(
    double val ) noexcept
{
    data[0]+=val; data[1]+=val; data[2]+=val;
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator-=(
    double val ) noexcept
{
    data[0]-=val; data[1]-=val; data[2]-=val;
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator*=(
    double val ) noexcept
{
    data[0]*=val; data[1]*=val; data[2]*=val;
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator/=(
    double val ) noexcept
{
    data[0]/=val; data[1]/=val; data[2]/=val;
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator+=(
    const ScalarVariable& sc ) noexcept
{
    data[0]+=sc(); data[1]+=sc(); data[2]+=sc();
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator-=(
    const ScalarVariable& sc ) noexcept
{
    data[0]-=sc(); data[1]-=sc(); data[2]-=sc();
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator*=(
    const ScalarVariable& sc ) noexcept
{
    data[0]*=sc(); data[1]*=sc(); data[2]*=sc();
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator/=(
    const ScalarVariable& sc ) noexcept
{
    data[0]/=sc(); data[1]/=sc(); data[2]/=sc();
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator+=(
    const VectorVariable<3U>& v ) noexcept
{
    data[0]+=v.data[0]; data[1]+=v.data[1]; data[2]+=v.data[2];
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator-=(
    const VectorVariable<3U>& v ) noexcept
{
    data[0]-=v.data[0]; data[1]-=v.data[1]; data[2]-=v.data[2];
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator*=(
    const VectorVariable<3U>& v ) noexcept
{
    data[0]*=v.data[0]; data[1]*=v.data[1]; data[2]*=v.data[2];
    return *this;
}

inline VectorVariable<3U>& VectorVariable<3U>::operator/=(
    const VectorVariable<3U>& v ) noexcept
{
    data[0]/=v.data[0]; data[1]/=v.data[1]; data[2]/=v.data[2];
    return *this;
}

inline bool VectorVariable<3U>::operator==(
    const VectorVariable<3U>& v ) const noexcept
{
    return ( flag == v.flag && data == v.data );
}

inline bool VectorVariable<3U>::operator!=(
    const VectorVariable<3U>& v ) const noexcept
{
    return ( flag != v.flag || data != v.data );
}

inline bool VectorVariable<3U>::operator<(
    const VectorVariable<3U>& v ) const noexcept
{
    return ( this < &v );
}

inline double VectorVariable<3U>::operator&(
    const VectorVariable<3U>& v ) const noexcept
{
    return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2];
}

inline VectorVariable<3U> VectorVariable<3U>::operator%(
    const VectorVariable<3U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1], flag[2],
        data[1]*v.data[2] - v.data[1]*data[2],
       -(data[0]*v.data[2] - v.data[0]*data[2]),
        data[0]*v.data[1] - v.data[0]*data[1] );
}

inline double VectorVariable<3U>::Length() const noexcept
{
    return std::hypot( data[0], data[1], data[2] );
}

inline Point<3U> VectorVariable<3U>::P() const noexcept
{
    return csmp::Point<3U>( data[0], data[1], data[2] );
}

inline bool VectorVariable<3U>::IsWithinRange(
    double vmin, double vmax ) const noexcept
{
    if ( data[0] < vmin || data[0] > vmax ) return false;
    if ( data[1] < vmin || data[1] > vmax ) return false;
    if ( data[2] < vmin || data[2] > vmax ) return false;
    return true;
}

inline bool VectorVariable<3U>::Has_NaN_Values() const noexcept
{
    return std::isnan( data[0] ) ||
           std::isnan( data[1] ) ||
           std::isnan( data[2] );
}

inline double VectorVariable<3U>::DotProduct(
    const csmp::Point<3U>& p ) const noexcept
{
    return data[0]*p[0] + data[1]*p[1] + data[2]*p[2];
}

inline double VectorVariable<3U>::DotProduct(
    const VectorVariable<3U>& v ) const noexcept
{
    return data[0]*v[0] + data[1]*v[1] + data[2]*v[2];
}

inline VectorVariable<3U> VectorVariable<3U>::CrossProduct(
    const csmp::Point<3U>& p ) const noexcept
{
    return VectorVariable<3U>( flag[0], flag[1], flag[2],
        data[1]*p[2] - data[2]*p[1],
        data[2]*p[0] - data[0]*p[2],
        data[0]*p[1] - data[1]*p[0] );
}

inline VectorVariable<3U> VectorVariable<3U>::CrossProduct(
    const VectorVariable<3U>& v ) const noexcept
{
    return VectorVariable<3U>( flag[0], flag[1], flag[2],
        data[1]*v[2] - data[2]*v[1],
        data[2]*v[0] - data[0]*v[2],
        data[0]*v[1] - data[1]*v[0] );
}

inline VectorVariable<3U> VectorVariable<3U>::ProjectOnto(
    const std::vector<double>& v ) const noexcept
{
    const double ratio = ( data[0]*v[0] + data[1]*v[1] + data[2]*v[2] )
                       / ( v[0]*v[0]   + v[1]*v[1]   + v[2]*v[2]   );
    return VectorVariable<3U>( flag[0], flag[1], flag[2],
        v[0]*ratio, v[1]*ratio, v[2]*ratio );
}

inline VectorVariable<3U> VectorVariable<3U>::ProjectOnto(
    const VectorVariable<3U>& v ) const noexcept
{
    const double ratio = ( data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2] )
                       / ( v.data[0]*v.data[0] + v.data[1]*v.data[1] + v.data[2]*v.data[2] );
    return VectorVariable<3U>( flag[0], flag[1], flag[2],
        v.data[0]*ratio, v.data[1]*ratio, v.data[2]*ratio );
}

inline void VectorVariable<3U>::Invert() noexcept
{
    data[0] *= -1.0;
    data[1] *= -1.0;
    data[2] *= -1.0;
}

inline void VectorVariable<3U>::EuclideanNormalize() noexcept
{
    const double norm = std::sqrt( data[0]*data[0]
                                 + data[1]*data[1]
                                 + data[2]*data[2] );
    if ( norm == 0.0 ) return;
    data[0] /= norm;
    data[1] /= norm;
    data[2] /= norm;
}

// ============================================================================
//  Inline definitions — factory functions
// ============================================================================

inline VectorVariable<1U> makeVector(
    VARIABLE_FLAG fx, double vx ) noexcept
{
    return VectorVariable<1U>( fx, vx );
}

inline VectorVariable<2U> makeVector(
    VARIABLE_FLAG fx, VARIABLE_FLAG fy,
    double vx, double vy ) noexcept
{
    return VectorVariable<2U>( fx, fy, vx, vy );
}

inline VectorVariable<3U> makeVector(
    VARIABLE_FLAG fx, VARIABLE_FLAG fy, VARIABLE_FLAG fz,
    double vx, double vy, double vz ) noexcept
{
    return VectorVariable<3U>( fx, fy, fz, vx, vy, vz );
}

inline VectorVariable<3U> makeVector(
    const std::array<VARIABLE_FLAG,3U>& flags,
    const std::array<double,3U>&        vals ) noexcept
{
    return VectorVariable<3U>( flags[0], flags[1], flags[2],
                               vals[0],  vals[1],  vals[2] );
}

// ============================================================================
//  Inline definitions — binary operators with double
// ============================================================================

inline VectorVariable<3U> operator+(
    const VectorVariable<3U>& a, double b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)+b, a(1)+b, a(2)+b );
}

inline VectorVariable<3U> operator-(
    const VectorVariable<3U>& a, double b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)-b, a(1)-b, a(2)-b );
}

inline VectorVariable<3U> operator*(
    const VectorVariable<3U>& a, double b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)*b, a(1)*b, a(2)*b );
}

inline VectorVariable<3U> operator/(
    const VectorVariable<3U>& a, double b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)/b, a(1)/b, a(2)/b );
}

// ============================================================================
//  Inline definitions — binary operators with ScalarVariable
// ============================================================================

inline VectorVariable<3U> operator+(
    const VectorVariable<3U>& a, const ScalarVariable& b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)+b(), a(1)+b(), a(2)+b() );
}

inline VectorVariable<3U> operator-(
    const VectorVariable<3U>& a, const ScalarVariable& b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)-b(), a(1)-b(), a(2)-b() );
}

inline VectorVariable<3U> operator*(
    const VectorVariable<3U>& a, const ScalarVariable& b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)*b(), a(1)*b(), a(2)*b() );
}

inline VectorVariable<3U> operator/(
    const VectorVariable<3U>& a, const ScalarVariable& b ) noexcept
{
    return VectorVariable<3U>( a.Flag(0), a.Flag(1), a.Flag(2),
        a(0)/b(), a(1)/b(), a(2)/b() );
}

// ============================================================================
//  Inline definitions — Point op VectorVariable (3D)
// ============================================================================

inline Point<3U> operator+(
    const Point<3U>& p, const VectorVariable<3U>& vc ) noexcept
{
    return Point<3U>( p[0]+vc[0], p[1]+vc[1], p[2]+vc[2] );
}

inline Point<3U> operator-(
    const Point<3U>& p, const VectorVariable<3U>& vc ) noexcept
{
    return Point<3U>( p[0]-vc[0], p[1]-vc[1], p[2]-vc[2] );
}

inline Point<3U> operator*(
    const Point<3U>& p, const VectorVariable<3U>& vc ) noexcept
{
    return Point<3U>( p[0]*vc[0], p[1]*vc[1], p[2]*vc[2] );
}

inline Point<3U> operator/(
    const Point<3U>& p, const VectorVariable<3U>& vc ) noexcept
{
    return Point<3U>( p[0]/vc[0], p[1]/vc[1], p[2]/vc[2] );
}

} // namespace csmp

#endif // CSMP_VECTOR_VARIABLE_H

