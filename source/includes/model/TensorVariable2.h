#ifndef CSMP_TENSOR_VARIABLE2_H
#define CSMP_TENSOR_VARIABLE2_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "VectorVariable2.h"

namespace csmp {

/**
@brief 2D full specialisation of the TensorVariable class template.

@author S.K. Matthaei
@author S. Geiger
@author Stephen G. Roberts
@date 2001
@author (refactored) 2024

@section squaring Squaring operations

Three distinct squaring operations are provided to avoid ambiguity:

  HadamardSquared() — component-wise: T_ij^2
  MatrixSquared()   — matrix product: result_ij = sum_k T_ik * T_kj
  DoubleContraction() — scalar T:T = sum_ij T_ij^2

@section notes Notes

IsWithinRange checks all four elements (not just the diagonal),
consistent with the 2D case where off-diagonal elements are
physically meaningful.
*/
template<>
class TensorVariable<2U>
{
public:
    static constexpr VARIABLE_TYPE VariableType = TENSOR;

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /** Default constructor: all elements NaN, all flags ANY. */
    TensorVariable() noexcept;

    /**
    Constructs an isotropic diagonal tensor with all diagonal elements
    equal to @p val and all flags equal to @p f. Off-diagonal elements
    are zero.
    */
    TensorVariable( VARIABLE_FLAG f, double val ) noexcept;

    /**
    Full initialisation with a single flag for both diagonal elements.
    */
    TensorVariable( VARIABLE_FLAG f,
                    double v11, double v12,
                    double v21, double v22 ) noexcept;

    /**
    Full initialisation with independent flags for each diagonal element.
    */
    TensorVariable( const VARIABLE_FLAG& f11, const VARIABLE_FLAG& f22,
                    const double& v11, const double& v12,
                    const double& v21, const double& v22 ) noexcept;

    // -----------------------------------------------------------------------
    // Element access
    // -----------------------------------------------------------------------

    /** Read/write access to element (i,j). */
    double&       operator()( uint32_t i, uint32_t j )       noexcept;

    /** Read-only access to element (i,j). */
    const double& operator()( uint32_t i, uint32_t j ) const noexcept;

    /**
    Alternative mutator: accesses elements 0..3 sequentially row by row.
    Element k maps to row k/2, column k%2.
    */
    void   Component( uint32_t k, double val ) noexcept;

    /**
    Alternative accessor: accesses elements 0..3 sequentially row by row.
    */
    [[nodiscard]] double Component( uint32_t k ) const noexcept;

    /** Number of entries in the tensor (2×2 = 4). */
    static constexpr uint32_t Size() noexcept { return 4U; }

    // -----------------------------------------------------------------------
    // Arithmetic — returning new tensor (flags from *this)
    // -----------------------------------------------------------------------

    TensorVariable operator+( double val )             const noexcept;
    TensorVariable operator-( double val )             const noexcept;
    TensorVariable operator*( double val )             const noexcept;
    TensorVariable operator/( double val )             const noexcept;

    /** Element-by-element addition. */
    TensorVariable operator+( const TensorVariable& )  const noexcept;

    /** Element-by-element subtraction. */
    TensorVariable operator-( const TensorVariable& )  const noexcept;

    /**
    Matrix multiplication: result_ij = sum_k this_ik * ts_kj.
    Uses unrolled loops for performance.
    */
    TensorVariable operator*( const TensorVariable& )  const noexcept;

    /** Matrix-vector multiplication: result = this * vc (column vector). */
    VectorVariable<2U> operator*( const VectorVariable<2U>& ) const noexcept;

    /** Matrix-point multiplication: result = this * v. */
    Point<2U>          operator*( const Point<2U>& )           const noexcept;

    /** Element-by-element division. */
    TensorVariable operator/( const TensorVariable& )  const noexcept;

    // -----------------------------------------------------------------------
    // Compound assignment
    // -----------------------------------------------------------------------

    TensorVariable& operator+=( double val )           noexcept;
    TensorVariable& operator-=( double val )           noexcept;
    TensorVariable& operator*=( double val )           noexcept;
    TensorVariable& operator/=( double val )           noexcept;

    TensorVariable& operator+=( const ScalarVariable& ) noexcept;
    TensorVariable& operator-=( const ScalarVariable& ) noexcept;
    TensorVariable& operator*=( const ScalarVariable& ) noexcept;
    TensorVariable& operator/=( const ScalarVariable& ) noexcept;

    /** Element-by-element addition in place. */
    TensorVariable& operator+=( const TensorVariable& ) noexcept;

    /** Element-by-element subtraction in place. */
    TensorVariable& operator-=( const TensorVariable& ) noexcept;

    /** Element-by-element division in place. */
    TensorVariable& operator/=( const TensorVariable& ) noexcept;

    /**
    Matrix multiplication in place: this = this * ts.
    Creates a temporary tensor internally.
    Prefer operator*() when assigning to a new variable.
    */
    TensorVariable& operator*=( const TensorVariable& ) noexcept;

    // -----------------------------------------------------------------------
    // Assignment
    // -----------------------------------------------------------------------

    /** Sets all elements to @p val. Flags are not modified. */
    TensorVariable& operator=( double val )                    noexcept;

    /**
    Sets all elements to the scalar value and all flags to the
    scalar's flag.
    */
    TensorVariable& operator=( const ScalarVariable& )         noexcept;

    /**
    Assigns the vector to the diagonal elements of the tensor.
    Off-diagonal elements are set to zero. Flags are taken from
    the vector components.
    */
    TensorVariable& operator=( const VectorVariable<2U>& )     noexcept;

    // -----------------------------------------------------------------------
    // Comparison
    // -----------------------------------------------------------------------

    /** Element-by-element equality of values and flags. */
    bool operator==( const TensorVariable& ) const noexcept;

    /** Element-by-element inequality. */
    bool operator!=( const TensorVariable& ) const noexcept;

    /**
    Ordering by pointer address — satisfies strict weak ordering for
    use in STL associative containers. Does not compare values.
    */
    bool operator<( const TensorVariable& )  const noexcept;

    // -----------------------------------------------------------------------
    // Squaring operations (three distinct semantics)
    // -----------------------------------------------------------------------

    /**
    Component-wise (Hadamard) squaring in place: T_ij = T_ij^2.
    This is NOT the matrix product. Use MatrixSquared() for T*T.
    */
    TensorVariable& HadamardSquared() noexcept;

    /**
    Matrix product of this tensor with itself in place:
    result_ij = sum_k T_ik * T_kj.
    Use HadamardSquared() for component-wise squaring.
    */
    TensorVariable& MatrixSquared() noexcept;

    /**
    Double contraction T:T = sum_ij T_ij^2.
    Returns the Frobenius norm squared as a scalar.
    The tensor itself is not modified.
    */
    [[nodiscard]] double DoubleContraction() const noexcept;

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /** Returns the flag of diagonal element @p i (0 or 1). Const version. */
    VARIABLE_FLAG  Flag( uint32_t i = 0 ) const noexcept;

    /** Returns a reference to the flag of diagonal element @p i. */
    VARIABLE_FLAG& Flag( uint32_t i = 0 )       noexcept;

    /**
    Returns true if all four elements lie within [@p vmin, @p vmax].
    All elements are checked (unlike the 3D version which checks only
    the diagonal).
    */
    bool IsWithinRange( double vmin, double vmax ) const noexcept;

    /**
    Returns true if any diagonal element is NaN.
    Off-diagonal elements are not checked.
    */
    bool Has_NaN_Values() const noexcept;

    /** Returns the smallest element in the tensor. */
    double MinElement() const noexcept;

    /** Returns the largest element in the tensor. */
    double MaxElement() const noexcept;

    /** Returns the determinant of the tensor. */
    double Determinant() const noexcept;

    /** Returns the trace (sum of diagonal elements). */
    double Trace() const noexcept;

    // -----------------------------------------------------------------------
    // Matrix operations
    // -----------------------------------------------------------------------

    /** Returns the adjoint (classical adjugate) matrix. */
    TensorVariable Adjoint()    const noexcept;

    /**
    Returns the inverse of the tensor.
    Prints a warning and returns a NaN tensor if the determinant is zero.
    */
    TensorVariable Inverse()    const noexcept;

    /** Returns the transpose of the tensor. */
    TensorVariable Transposed() const noexcept;

    /** Converts this tensor to the identity matrix. */
    void Identity() noexcept;

    // -----------------------------------------------------------------------
    // Diagonal operations
    // -----------------------------------------------------------------------

    /** Assigns @p f_00 and @p f_11 to the diagonal elements. */
    void DiagonalValues( double f_00, double f_11 ) noexcept;

    /** Assigns the first two elements of @p v to the diagonal. */
    void DiagonalValues( const std::vector<double>& v ) noexcept;

    /** Assigns the vector components to the diagonal. Flags are also copied. */
    void DiagonalValues( const VectorVariable<2U>& v ) noexcept;

    /** Assigns the vector to row @p i of the tensor. */
    void AssignToRow(    uint32_t i, VectorVariable<2U>& vc ) noexcept;

    /** Assigns the vector to column @p j of the tensor. */
    void AssignToColumn( uint32_t j, VectorVariable<2U>& vc ) noexcept;

    /** Returns row @p iRow as a VectorVariable. */
    VectorVariable<2U> Row(    uint32_t iRow ) const noexcept;

    /** Returns column @p iCol as a VectorVariable. */
    VectorVariable<2U> Column( uint32_t iCol ) const noexcept;

    // -----------------------------------------------------------------------
    // Eigenvalue / eigenvector methods
    // -----------------------------------------------------------------------

    /**
    Computes eigenvalues of this 2×2 tensor via the quadratic formula.
    Returns eigenvalues in descending order.
    Returns false and sets NaN if the discriminant is negative.
    */
    bool EigenValues( VectorVariable<2U>& eigenvalues ) const;

    /**
    Overload returning eigenvalues into a std::vector<double>.
    Delegates to the VectorVariable overload.
    */
    bool EigenValues( std::vector<double>& eigenvalues ) const;

    /**
    Computes eigenvalues and eigenvectors of a symmetric 2×2 tensor.
    Off-diagonal elements are symmetrised by averaging.
    Eigenvectors are stored as columns of @p eigenVecs.
    Eigenvalues are returned in descending order.

    @param bNormalize  If true, eigenvectors are normalised to unit length.
    */
    bool EigenSymmetric( VectorVariable<2U>& eigenVals,
                         TensorVariable<2U>& eigenVecs,
                         bool                bNormalize = false ) const;

    /**
    Shorthand delegating to EigenSymmetric for interface consistency
    with the 3D specialisation.
    */
    bool Eigen( VectorVariable<2U>& eigenVals,
                TensorVariable<2U>& eigenVecs,
                bool                normalise ) const;

    /**
    Computes eigenvalues and eigenvectors for a weakly non-symmetric
    2×2 tensor by decomposing into symmetric and skew-symmetric parts.
    Returns true if ||W||_F / ||S||_F <= @p tolerance.
    Returns false and issues WARNING if asymmetry exceeds tolerance.

    @param tolerance  Maximum permitted asymmetry ratio. Default 1e-6.
    */
    bool EigenWeaklyNonSymmetric( VectorVariable<2U>& eigenVals,
                                  TensorVariable<2U>& eigenVecs,
                                  double              tolerance = 1.0e-6 ) const;

    // -----------------------------------------------------------------------
    // I/O
    // -----------------------------------------------------------------------

    /** Prompts the user to enter tensor values from the command line. */
    void In();

    /** Prints the tensor to stdout. */
    void Out() const noexcept;

    /** Reads the tensor from a binary file stream. */
    bool In(  std::fstream& fp );

    /** Writes the tensor to a binary file stream. */
    bool Out( std::fstream& fp ) const;

private:
    std::array<VARIABLE_FLAG, 2U>           flag;
    std::array<std::array<double, 2U>, 2U>  data;
};


// ============================================================================
//  Free functions
// ============================================================================

/**
Vector-matrix multiplication: result = vc^T * ts
(equivalent to ts^T * vc as a column vector).
*/
VectorVariable<2U> operator*( const VectorVariable<2U>& vc,
                               const TensorVariable<2U>& ts ) noexcept;

/** Point-matrix multiplication. */
Point<2U> operator*( const Point<2U>& vc,
                     const TensorVariable<2U>& ts ) noexcept;


// ============================================================================
//  Inline definitions
// ============================================================================

inline TensorVariable<2U>::TensorVariable() noexcept
    : flag{ { ANY, ANY } },
      data{ { { std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN() },
              { std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN() } } }
{}

inline TensorVariable<2U>::TensorVariable(
    VARIABLE_FLAG f, double val ) noexcept
    : flag{ { f, f } },
      data{ { { val, 0. }, { 0., val } } }
{}

inline TensorVariable<2U>::TensorVariable(
    VARIABLE_FLAG f,
    double v11, double v12,
    double v21, double v22 ) noexcept
    : flag{ { f, f } },
      data{ { { v11, v12 }, { v21, v22 } } }
{}

inline TensorVariable<2U>::TensorVariable(
    const VARIABLE_FLAG& f11, const VARIABLE_FLAG& f22,
    const double& v11, const double& v12,
    const double& v21, const double& v22 ) noexcept
    : flag{ { f11, f22 } },
      data{ { { v11, v12 }, { v21, v22 } } }
{}

inline double& TensorVariable<2U>::operator()(
    uint32_t i, uint32_t j ) noexcept
{
#ifndef NDEBUG
    if ( i >= 2U ) {
        std::cerr << "\nTensorVariable<2U>::operator(): row violation i="
                  << i << "\n";
        return data[0][0];
    }
    if ( j >= 2U ) {
        std::cerr << "\nTensorVariable<2U>::operator(): col violation j="
                  << j << "\n";
        return data[0][0];
    }
#endif
    return data[i][j];
}

inline const double& TensorVariable<2U>::operator()(
    uint32_t i, uint32_t j ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 2U ) {
        std::cerr << "\nTensorVariable<2U>::operator(): row violation i="
                  << i << "\n";
        return data[0][0];
    }
    if ( j >= 2U ) {
        std::cerr << "\nTensorVariable<2U>::operator(): col violation j="
                  << j << "\n";
        return data[0][0];
    }
#endif
    return data[i][j];
}

inline void TensorVariable<2U>::Component(
    uint32_t k, double val ) noexcept
{
    assert( k < Size() );
    data[k / 2][k % 2] = val;
}

inline double TensorVariable<2U>::Component( uint32_t k ) const noexcept
{
    assert( k < Size() );
    return data[k / 2][k % 2];
}

inline VARIABLE_FLAG TensorVariable<2U>::Flag( uint32_t i ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 2U ) {
        std::cerr << "\nTensorVariable<2U>::Flag(): violation i=" << i << "\n";
        return flag[0];
    }
#endif
    return flag[i];
}

inline VARIABLE_FLAG& TensorVariable<2U>::Flag( uint32_t i ) noexcept
{
#ifndef NDEBUG
    if ( i >= 2U ) {
        std::cerr << "\nTensorVariable<2U>::Flag(): violation i=" << i << "\n";
        return flag[0];
    }
#endif
    return flag[i];
}

inline bool TensorVariable<2U>::IsWithinRange(
    double vmin, double vmax ) const noexcept
{
    if ( data[0][0] < vmin || data[0][0] > vmax ) return false;
    if ( data[0][1] < vmin || data[0][1] > vmax ) return false;
    if ( data[1][0] < vmin || data[1][0] > vmax ) return false;
    if ( data[1][1] < vmin || data[1][1] > vmax ) return false;
    return true;
}

inline bool TensorVariable<2U>::Has_NaN_Values() const noexcept
{
    return std::isnan( data[0][0] ) || std::isnan( data[1][1] );
}

inline double TensorVariable<2U>::Determinant() const noexcept
{
    return data[0][0] * data[1][1] - data[0][1] * data[1][0];
}

inline double TensorVariable<2U>::Trace() const noexcept
{
    return data[0][0] + data[1][1];
}

inline bool TensorVariable<2U>::operator==(
    const TensorVariable<2U>& ts ) const noexcept
{
    return ( data == ts.data && flag == ts.flag );
}

inline bool TensorVariable<2U>::operator!=(
    const TensorVariable<2U>& ts ) const noexcept
{
    return !( *this == ts );
}

inline bool TensorVariable<2U>::operator<(
    const TensorVariable<2U>& t ) const noexcept
{
    return ( this < &t );
}

// -----------------------------------------------------------------------
//  Squaring operations
// -----------------------------------------------------------------------

inline TensorVariable<2U>& TensorVariable<2U>::HadamardSquared() noexcept
{
    data[0][0] *= data[0][0];
    data[0][1] *= data[0][1];
    data[1][0] *= data[1][0];
    data[1][1] *= data[1][1];
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::MatrixSquared() noexcept
{
    return *this *= *this;
}

inline double TensorVariable<2U>::DoubleContraction() const noexcept
{
    return data[0][0]*data[0][0] + data[0][1]*data[0][1]
         + data[1][0]*data[1][0] + data[1][1]*data[1][1];
}

// -----------------------------------------------------------------------
//  Arithmetic — returning new tensor
// -----------------------------------------------------------------------

inline TensorVariable<2U> TensorVariable<2U>::operator+(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]+val, data[0][1]+val,
        data[1][0]+val, data[1][1]+val );
}

inline TensorVariable<2U> TensorVariable<2U>::operator-(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]-val, data[0][1]-val,
        data[1][0]-val, data[1][1]-val );
}

inline TensorVariable<2U> TensorVariable<2U>::operator*(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]*val, data[0][1]*val,
        data[1][0]*val, data[1][1]*val );
}

inline TensorVariable<2U> TensorVariable<2U>::operator/(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]/val, data[0][1]/val,
        data[1][0]/val, data[1][1]/val );
}

inline TensorVariable<2U> TensorVariable<2U>::operator+(
    const TensorVariable<2U>& t ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]+t.data[0][0], data[0][1]+t.data[0][1],
        data[1][0]+t.data[1][0], data[1][1]+t.data[1][1] );
}

inline TensorVariable<2U> TensorVariable<2U>::operator-(
    const TensorVariable<2U>& t ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]-t.data[0][0], data[0][1]-t.data[0][1],
        data[1][0]-t.data[1][0], data[1][1]-t.data[1][1] );
}

inline TensorVariable<2U> TensorVariable<2U>::operator*(
    const TensorVariable<2U>& ts ) const noexcept
{
    TensorVariable<2U> temp;
    // Unrolled matrix multiplication
    temp.data[0][0] = data[0][0]*ts.data[0][0] + data[0][1]*ts.data[1][0];
    temp.data[0][1] = data[0][0]*ts.data[0][1] + data[0][1]*ts.data[1][1];
    temp.data[1][0] = data[1][0]*ts.data[0][0] + data[1][1]*ts.data[1][0];
    temp.data[1][1] = data[1][0]*ts.data[0][1] + data[1][1]*ts.data[1][1];
    temp.flag[0] = flag[0];
    temp.flag[1] = flag[1];
    return temp;
}

inline TensorVariable<2U> TensorVariable<2U>::operator/(
    const TensorVariable<2U>& ts ) const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0]/ts.data[0][0], data[0][1]/ts.data[0][1],
        data[1][0]/ts.data[1][0], data[1][1]/ts.data[1][1] );
}

inline VectorVariable<2U> TensorVariable<2U>::operator*(
    const VectorVariable<2U>& vc ) const noexcept
{
    return VectorVariable<2U>( vc.Flag(0), vc.Flag(1),
        data[0][0]*vc[0] + data[0][1]*vc[1],
        data[1][0]*vc[0] + data[1][1]*vc[1] );
}

inline Point<2U> TensorVariable<2U>::operator*(
    const Point<2U>& v ) const noexcept
{
    return Point<2U>(
        data[0][0]*v[0] + data[0][1]*v[1],
        data[1][0]*v[0] + data[1][1]*v[1] );
}

// -----------------------------------------------------------------------
//  Compound assignment
// -----------------------------------------------------------------------

inline TensorVariable<2U>& TensorVariable<2U>::operator+=(
    double val ) noexcept
{
    data[0][0]+=val; data[0][1]+=val;
    data[1][0]+=val; data[1][1]+=val;
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator-=(
    double val ) noexcept
{
    data[0][0]-=val; data[0][1]-=val;
    data[1][0]-=val; data[1][1]-=val;
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator*=(
    double val ) noexcept
{
    data[0][0]*=val; data[0][1]*=val;
    data[1][0]*=val; data[1][1]*=val;
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator/=(
    double val ) noexcept
{
    data[0][0]/=val; data[0][1]/=val;
    data[1][0]/=val; data[1][1]/=val;
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator+=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]+=sc(); data[0][1]+=sc();
    data[1][0]+=sc(); data[1][1]+=sc();
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator-=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]-=sc(); data[0][1]-=sc();
    data[1][0]-=sc(); data[1][1]-=sc();
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator*=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]*=sc(); data[0][1]*=sc();
    data[1][0]*=sc(); data[1][1]*=sc();
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator/=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]/=sc(); data[0][1]/=sc();
    data[1][0]/=sc(); data[1][1]/=sc();
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator+=(
    const TensorVariable<2U>& ts ) noexcept
{
    data[0][0]+=ts.data[0][0]; data[0][1]+=ts.data[0][1];
    data[1][0]+=ts.data[1][0]; data[1][1]+=ts.data[1][1];
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator-=(
    const TensorVariable<2U>& ts ) noexcept
{
    data[0][0]-=ts.data[0][0]; data[0][1]-=ts.data[0][1];
    data[1][0]-=ts.data[1][0]; data[1][1]-=ts.data[1][1];
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator/=(
    const TensorVariable<2U>& ts ) noexcept
{
    data[0][0]/=ts.data[0][0]; data[0][1]/=ts.data[0][1];
    data[1][0]/=ts.data[1][0]; data[1][1]/=ts.data[1][1];
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator*=(
    const TensorVariable<2U>& ts ) noexcept
{
    TensorVariable<2U> temp;
    // Unrolled matrix multiplication
    temp.data[0][0] = data[0][0]*ts.data[0][0] + data[0][1]*ts.data[1][0];
    temp.data[0][1] = data[0][0]*ts.data[0][1] + data[0][1]*ts.data[1][1];
    temp.data[1][0] = data[1][0]*ts.data[0][0] + data[1][1]*ts.data[1][0];
    temp.data[1][1] = data[1][0]*ts.data[0][1] + data[1][1]*ts.data[1][1];
    temp.flag[0] = flag[0];
    temp.flag[1] = flag[1];
    return *this = std::move( temp );
}

// -----------------------------------------------------------------------
//  Assignment
// -----------------------------------------------------------------------

inline TensorVariable<2U>& TensorVariable<2U>::operator=(
    double val ) noexcept
{
    data[0][0]=val; data[0][1]=val;
    data[1][0]=val; data[1][1]=val;
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator=(
    const ScalarVariable& sc ) noexcept
{
    flag[0] = flag[1] = sc.Flag();
    data[0][0]=sc(); data[0][1]=sc();
    data[1][0]=sc(); data[1][1]=sc();
    return *this;
}

inline TensorVariable<2U>& TensorVariable<2U>::operator=(
    const VectorVariable<2U>& vc ) noexcept
{
    flag[0] = vc.Flag(0);
    flag[1] = vc.Flag(1);
    data[0][0] = vc(0); data[0][1] = 0.0;
    data[1][0] = 0.0;   data[1][1] = vc(1);
    return *this;
}

// -----------------------------------------------------------------------
//  Matrix operations
// -----------------------------------------------------------------------

inline void TensorVariable<2U>::Identity() noexcept
{
    data[0][0] = 1.0; data[0][1] = 0.0;
    data[1][0] = 0.0; data[1][1] = 1.0;
}

inline TensorVariable<2U> TensorVariable<2U>::Transposed() const noexcept
{
    return TensorVariable( flag[0], flag[1],
        data[0][0], data[1][0],
        data[0][1], data[1][1] );
}

inline TensorVariable<2U> TensorVariable<2U>::Adjoint() const noexcept
{
    return TensorVariable( flag[0], flag[1],
         data[1][1], -data[1][0],
        -data[0][1],  data[0][0] );
}

inline TensorVariable<2U> TensorVariable<2U>::Inverse() const noexcept
{
    double det = Determinant();
    if ( det == 0.0 )
    {
        std::cerr << "\nTensorVariable<2U>::Inverse: Determinant = 0\n";
        return TensorVariable<2U>();
    }
    det = 1.0 / det;
    return TensorVariable( flag[0], flag[1],
         det *  data[1][1],
         det * -data[0][1],
         det * -data[1][0],
         det *  data[0][0] );
}

inline double TensorVariable<2U>::MinElement() const noexcept
{
    double me = data[0][0];
    if ( data[0][1] < me ) me = data[0][1];
    if ( data[1][0] < me ) me = data[1][0];
    if ( data[1][1] < me ) me = data[1][1];
    return me;
}

inline double TensorVariable<2U>::MaxElement() const noexcept
{
    double me = data[0][0];
    if ( data[0][1] > me ) me = data[0][1];
    if ( data[1][0] > me ) me = data[1][0];
    if ( data[1][1] > me ) me = data[1][1];
    return me;
}

// -----------------------------------------------------------------------
//  Diagonal operations
// -----------------------------------------------------------------------

inline void TensorVariable<2U>::DiagonalValues(
    double f_00, double f_11 ) noexcept
{
    data[0][0] = f_00;
    data[1][1] = f_11;
}

inline void TensorVariable<2U>::DiagonalValues(
    const std::vector<double>& v ) noexcept
{
    data[0][0] = v[0];
    data[1][1] = v[1];
}

inline void TensorVariable<2U>::DiagonalValues(
    const VectorVariable<2U>& v ) noexcept
{
    data[0][0] = v[0];
    data[1][1] = v[1];
    flag[0] = v.Flag(0);
    flag[1] = v.Flag(1);
}

inline void TensorVariable<2U>::AssignToRow(
    uint32_t iRow, VectorVariable<2U>& vc ) noexcept
{
    flag[ iRow == 0U ? 0U : 1U ] = vc.Flag( iRow == 0U ? 0U : 1U );
    data[iRow][0] = vc[0];
    data[iRow][1] = vc[1];
}

inline void TensorVariable<2U>::AssignToColumn(
    uint32_t iCol, VectorVariable<2U>& vc ) noexcept
{
    flag[ iCol == 0U ? 0U : 1U ] = vc.Flag( iCol == 0U ? 0U : 1U );
    data[0][iCol] = vc[0];
    data[1][iCol] = vc[1];
}

inline VectorVariable<2U> TensorVariable<2U>::Row(
    uint32_t iRow ) const noexcept
{
    return VectorVariable<2U>( flag[iRow], flag[iRow],
                               data[iRow][0], data[iRow][1] );
}

inline VectorVariable<2U> TensorVariable<2U>::Column(
    uint32_t iCol ) const noexcept
{
    return VectorVariable<2U>( flag[iCol], flag[iCol],
                               data[0][iCol], data[1][iCol] );
}

inline bool TensorVariable<2U>::Eigen(
    VectorVariable<2U>& eigenVals,
    TensorVariable<2U>& eigenVecs,
    bool                normalise ) const
{
    return EigenSymmetric( eigenVals, eigenVecs, normalise );
}

// -----------------------------------------------------------------------
//  Free function inline definitions
// -----------------------------------------------------------------------

inline VectorVariable<2U> operator*(
    const VectorVariable<2U>& vc,
    const TensorVariable<2U>& ts ) noexcept
{
    return VectorVariable<2U>( vc.Flag(0), vc.Flag(1),
        ts(0,0)*vc[0] + ts(1,0)*vc[1],
        ts(0,1)*vc[0] + ts(1,1)*vc[1] );
}

inline Point<2U> operator*(
    const Point<2U>& vc,
    const TensorVariable<2U>& ts ) noexcept
{
    return Point<2U>(
        ts(0,0)*vc[0] + ts(1,0)*vc[1],
        ts(0,1)*vc[0] + ts(1,1)*vc[1] );
}

} // namespace csmp

#endif // CSMP_TENSOR_VARIABLE2_H

