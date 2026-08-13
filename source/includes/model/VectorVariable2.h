#ifndef CSMP_VECTOR_VARIABLE2_H
#define CSMP_VECTOR_VARIABLE2_H

#include "CSMP_definitions.h"
#include "Point.h"
#include "ScalarVariable.h"

namespace csmp {

/**
@brief 2D full specialisation of the VectorVariable class template.

@author S.K. Matthaei
@author S. Geiger
@author Stephen G. Roberts
@date 2001
@author (refactored) 2024

@note operator^ has been removed. Use Pow(double) instead.
*/
template<>
class VectorVariable<2U>
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
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2,
                    double val1, double val2 ) noexcept;

    /** Initialises from an STL vector. All flags set to ANY. */
    explicit VectorVariable( const std::vector<double>& v ) noexcept;

    /** Initialises from a Point. All flags set to ANY. */
    explicit VectorVariable( const csmp::Point<2U>& p ) noexcept;

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

    /** Number of elements (always 2). */
    static constexpr uint32_t Size() noexcept { return 2U; }

    // -----------------------------------------------------------------------
    // Assignment
    // -----------------------------------------------------------------------

    /** Sets all elements to @p val. Flags are not modified. */
    VectorVariable& operator=( double val )                noexcept;

    /** Sets all elements from the point. Flags are not modified. */
    VectorVariable& operator=( const csmp::Point<2U>& p ) noexcept;

    /** Sets all elements to the scalar value and all flags to the scalar's flag. */
    VectorVariable& operator=( const ScalarVariable& sc ) noexcept;

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
    // Queries
    // -----------------------------------------------------------------------

    /** Returns the flag of element @p i. Const version. */
    VARIABLE_FLAG  Flag( uint32_t i = 0 ) const noexcept;

    /** Returns a reference to the flag of element @p i. */
    VARIABLE_FLAG& Flag( uint32_t i = 0 )       noexcept;

    /** Euclidean length. */
    double Length() const noexcept;

    /**
    Returns the angle in degrees between this vector and @p v.
    Returns 90 if the denominator is zero.
    */
    double AngleTo( const VectorVariable& v ) const noexcept;

    /** Returns a csmp::Point initialised with the vector values. */
    Point<2U> P() const noexcept;

    /**
    Returns true if both element values lie within [@p vmin, @p vmax].
    */
    bool IsWithinRange( double vmin, double vmax ) const noexcept;

    /** Returns true if any element is NaN. */
    bool Has_NaN_Values() const noexcept;

    // -----------------------------------------------------------------------
    // Vector operations
    // -----------------------------------------------------------------------

    /** Dot product with a Point. */
    double DotProduct( const csmp::Point<2U>& p ) const noexcept;

    /** Dot product with a VectorVariable. */
    double DotProduct( const VectorVariable& v ) const noexcept;

    /**
    2D cross product: returns a vector with x=0 and
    y = this.x * p.y - this.y * p.x.
    */
    VectorVariable CrossProduct( const csmp::Point<2U>& p ) const noexcept;

    /** 2D cross product with a VectorVariable. */
    VectorVariable CrossProduct( const VectorVariable& v ) const noexcept;

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

    friend class TensorVariable<2U>;

private:
    std::array<VARIABLE_FLAG, 2U> flag;
    std::array<double, 2U>        data;
};


// ============================================================================
//  Inline definitions
// ============================================================================

inline VectorVariable<2U>::VectorVariable() noexcept
    : flag{ { ANY, ANY } },
      data{ { std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN() } }
{}

inline VectorVariable<2U>::VectorVariable(
    VARIABLE_FLAG f, double val ) noexcept
    : flag{ { f, f } },
      data{ { val, val } }
{}

inline VectorVariable<2U>::VectorVariable(
    VARIABLE_FLAG f1, VARIABLE_FLAG f2,
    double val1, double val2 ) noexcept
    : flag{ { f1, f2 } },
      data{ { val1, val2 } }
{}

inline VectorVariable<2U>::VectorVariable(
    const std::vector<double>& v ) noexcept
    : flag{ { ANY, ANY } },
      data{ { v[0], v[1] } }
{}

inline VectorVariable<2U>::VectorVariable(
    const csmp::Point<2U>& p ) noexcept
    : flag{ { ANY, ANY } },
      data{ { p[0], p[1] } }
{}

inline double& VectorVariable<2U>::operator()( uint32_t i ) noexcept
{
    return data[i < 1U ? 0U : 1U];
}

inline const double& VectorVariable<2U>::operator()( uint32_t i ) const noexcept
{
    return data[i < 1U ? 0U : 1U];
}

inline double VectorVariable<2U>::operator[]( uint32_t i ) const noexcept
{
    return data[i < 1U ? 0U : 1U];
}

inline void VectorVariable<2U>::Component( uint32_t i, double val ) noexcept
{
    data[i < 1U ? 0U : 1U] = val;
}

inline double VectorVariable<2U>::Component( uint32_t i ) const noexcept
{
    return data[i < 1U ? 0U : 1U];
}

inline VARIABLE_FLAG& VectorVariable<2U>::Flag( uint32_t i ) noexcept
{
    return flag[i < 1U ? 0U : 1U];
}

inline VARIABLE_FLAG VectorVariable<2U>::Flag( uint32_t i ) const noexcept
{
    return flag[i < 1U ? 0U : 1U];
}

inline VectorVariable<2U>& VectorVariable<2U>::operator=(
    double val ) noexcept
{
    data[0] = data[1] = val;
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator=(
    const csmp::Point<2U>& p ) noexcept
{
    data[0] = p[0]; data[1] = p[1];
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator=(
    const ScalarVariable& sc ) noexcept
{
    flag[0] = flag[1] = sc.Flag();
    data[0] = data[1] = sc();
    return *this;
}

inline VectorVariable<2U> VectorVariable<2U>::operator+(
    const VectorVariable<2U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1],
        data[0]+v.data[0], data[1]+v.data[1] );
}

inline VectorVariable<2U> VectorVariable<2U>::operator-(
    const VectorVariable<2U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1],
        data[0]-v.data[0], data[1]-v.data[1] );
}

inline VectorVariable<2U> VectorVariable<2U>::operator*(
    const VectorVariable<2U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1],
        data[0]*v.data[0], data[1]*v.data[1] );
}

inline VectorVariable<2U> VectorVariable<2U>::operator/(
    const VectorVariable<2U>& v ) const noexcept
{
    return VectorVariable( flag[0], flag[1],
        data[0]/v.data[0], data[1]/v.data[1] );
}

inline VectorVariable<2U> VectorVariable<2U>::Pow(
    double exponent ) const noexcept
{
    return VectorVariable( flag[0], flag[1],
        std::pow( data[0], exponent ),
        std::pow( data[1], exponent ) );
}

inline VectorVariable<2U>& VectorVariable<2U>::operator+=(
    double val ) noexcept
{
    data[0]+=val; data[1]+=val;
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator-=(
    double val ) noexcept
{
    data[0]-=val; data[1]-=val;
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator*=(
    double val ) noexcept
{
    data[0]*=val; data[1]*=val;
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator/=(
    double val ) noexcept
{
    data[0]/=val; data[1]/=val;
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator+=(
    const ScalarVariable& sc ) noexcept
{
    data[0]+=sc(); data[1]+=sc();
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator-=(
    const ScalarVariable& sc ) noexcept
{
    data[0]-=sc(); data[1]-=sc();
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator*=(
    const ScalarVariable& sc ) noexcept
{
    data[0]*=sc(); data[1]*=sc();
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator/=(
    const ScalarVariable& sc ) noexcept
{
    data[0]/=sc(); data[1]/=sc();
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator+=(
    const VectorVariable<2U>& v ) noexcept
{
    data[0]+=v.data[0]; data[1]+=v.data[1];
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator-=(
    const VectorVariable<2U>& v ) noexcept
{
    data[0]-=v.data[0]; data[1]-=v.data[1];
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator*=(
    const VectorVariable<2U>& v ) noexcept
{
    data[0]*=v.data[0]; data[1]*=v.data[1];
    return *this;
}

inline VectorVariable<2U>& VectorVariable<2U>::operator/=(
    const VectorVariable<2U>& v ) noexcept
{
    data[0]/=v.data[0]; data[1]/=v.data[1];
    return *this;
}

inline bool VectorVariable<2U>::operator==(
    const VectorVariable<2U>& v ) const noexcept
{
    return ( flag == v.flag && data == v.data );
}

inline bool VectorVariable<2U>::operator!=(
    const VectorVariable<2U>& v ) const noexcept
{
    return ( flag != v.flag || data != v.data );
}

inline bool VectorVariable<2U>::operator<(
    const VectorVariable<2U>& v ) const noexcept
{
    return ( this < &v );
}

inline double VectorVariable<2U>::Length() const noexcept
{
    return std::hypot( data[0], data[1] );
}

inline Point<2U> VectorVariable<2U>::P() const noexcept
{
    return csmp::Point<2U>( data[0], data[1] );
}

inline bool VectorVariable<2U>::IsWithinRange(
    double vmin, double vmax ) const noexcept
{
    if ( data[0] < vmin || data[0] > vmax ) return false;
    if ( data[1] < vmin || data[1] > vmax ) return false;
    return true;
}

inline bool VectorVariable<2U>::Has_NaN_Values() const noexcept
{
    return std::isnan( data[0] ) || std::isnan( data[1] );
}

inline double VectorVariable<2U>::DotProduct(
    const csmp::Point<2U>& p ) const noexcept
{
    return data[0]*p[0] + data[1]*p[1];
}

inline double VectorVariable<2U>::DotProduct(
    const VectorVariable<2U>& v ) const noexcept
{
    return data[0]*v[0] + data[1]*v[1];
}

inline VectorVariable<2U> VectorVariable<2U>::CrossProduct(
    const csmp::Point<2U>& p ) const noexcept
{
    return VectorVariable<2U>( flag[0], flag[1],
        0.0,
        data[0]*p[1] - data[1]*p[0] );
}

inline VectorVariable<2U> VectorVariable<2U>::CrossProduct(
    const VectorVariable<2U>& v ) const noexcept
{
    return VectorVariable<2U>( flag[0], flag[1],
        0.0,
        data[0]*v[1] - data[1]*v[0] );
}

inline VectorVariable<2U> VectorVariable<2U>::ProjectOnto(
    const std::vector<double>& v ) const noexcept
{
    const double ratio = ( data[0]*v[0] + data[1]*v[1] )
                       / ( v[0]*v[0]   + v[1]*v[1]   );
    return VectorVariable<2U>( flag[0], flag[1],
        v[0]*ratio, v[1]*ratio );
}

inline VectorVariable<2U> VectorVariable<2U>::ProjectOnto(
    const VectorVariable<2U>& v ) const noexcept
{
    const double ratio = ( data[0]*v.data[0] + data[1]*v.data[1] )
                       / ( v.data[0]*v.data[0] + v.data[1]*v.data[1] );
    return VectorVariable<2U>( flag[0], flag[1],
        v.data[0]*ratio, v.data[1]*ratio );
}

inline void VectorVariable<2U>::Invert() noexcept
{
    data[0] *= -1.0;
    data[1] *= -1.0;
}

inline void VectorVariable<2U>::EuclideanNormalize() noexcept
{
    const double norm = std::hypot( data[0], data[1] );
    if ( norm == 0.0 ) return;
    data[0] /= norm;
    data[1] /= norm;
}

} // namespace csmp

#endif // CSMP_VECTOR_VARIABLE2_H

