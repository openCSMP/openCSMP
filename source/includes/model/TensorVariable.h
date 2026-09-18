// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_TENSOR_VARIABLE_H
#define CSMP_TENSOR_VARIABLE_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable1.h"
#include "TensorVariable2.h"

namespace csmp {

/**
@brief Full 3D specialisation of the TensorVariable class template.

@author S.K. Matthai
@author Stephen G. Roberts
@date 2001
@author (refactored) 2024

@section motivation Motivation

Template specialisation implementing optimised 3×3 tensors in CSMP.
Each diagonal element is associated with a VARIABLE_FLAG which can
indicate to finite-element computations that the corresponding degree
of freedom is enabled or restricted.

@section squaring Squaring operations

Three distinct squaring operations are provided to avoid ambiguity:

  Squared()       — component-wise (Hadamard): T_ij^2
  MatrixSquared() — matrix product: result_ij = sum_k T_ik * T_kj
  DoubleContraction() — scalar T:T = sum_ij T_ij^2 (Frobenius norm squared)

@section implementation Implementation

Using the template specialisation mechanism of C++. Specialisations
exist for the 1D, 2D, and 3D cases. All performance-critical methods
are inlined with unrolled loops.
*/
template<>
class TensorVariable<3U>
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
    equal to @p val and all flags equal to @p flag. Off-diagonal
    elements are zero.
    */
    TensorVariable( VARIABLE_FLAG flag, double val ) noexcept;

    /**
    Full initialisation with a single flag for all diagonal elements.
    */
    TensorVariable( VARIABLE_FLAG f,
                    double v11, double v12, double v13,
                    double v21, double v22, double v23,
                    double v31, double v32, double v33 ) noexcept;

    /**
    Full initialisation with independent flags for each diagonal element.
    */
    TensorVariable( VARIABLE_FLAG f11, VARIABLE_FLAG f22, VARIABLE_FLAG f33,
                    double v11, double v12, double v13,
                    double v21, double v22, double v23,
                    double v31, double v32, double v33 ) noexcept;

    // -----------------------------------------------------------------------
    // Element access
    // -----------------------------------------------------------------------

    /** Read/write access to element (i,j). */
    double&       operator()( uint32_t i, uint32_t j )       noexcept;

    /** Read-only access to element (i,j). */
    const double& operator()( uint32_t i, uint32_t j ) const noexcept;

    /**
    Alternative mutator: accesses elements 0..8 sequentially row by row.
    Element k maps to row k/3, column k%3.
    */
    void   Component( uint32_t k, double val ) noexcept;

    /**
    Alternative accessor: accesses elements 0..8 sequentially row by row.
    */
    [[nodiscard]] double Component( uint32_t k ) const noexcept;

    /** Number of entries in the tensor (3×3 = 9). */
    static constexpr uint32_t Size() noexcept { return 9U; }

    /**
    Sets all elements to @p newValue. The first argument is unused
    (retained for interface compatibility with other variable types).
    */
    void Resize( uint32_t, double newValue =
                     std::numeric_limits<double>::quiet_NaN() ) noexcept;

    // -----------------------------------------------------------------------
    // Arithmetic — returning new tensor (flags from *this)
    // -----------------------------------------------------------------------

    TensorVariable operator+( double val )              const noexcept;
    TensorVariable operator-( double val )              const noexcept;
    TensorVariable operator*( double val )              const noexcept;
    TensorVariable operator/( double val )              const noexcept;

    /** Element-by-element addition. */
    TensorVariable operator+( const TensorVariable& )   const noexcept;

    /** Element-by-element subtraction. */
    TensorVariable operator-( const TensorVariable& )   const noexcept;

    /**
    Matrix multiplication: result_ij = sum_k this_ik * ts_kj.
    Uses unrolled loops for performance.
    */
    TensorVariable operator*( const TensorVariable& )   const noexcept;

    /** Matrix-vector multiplication: result = this * vc (column vector). */
    VectorVariable<3U> operator*( const VectorVariable<3U>& ) const noexcept;

    /** Matrix-point multiplication: result = this * v. */
    Point<3U>          operator*( const Point<3U>& )           const noexcept;

    /** Element-by-element division. */
    TensorVariable operator/( const TensorVariable& )   const noexcept;

    // -----------------------------------------------------------------------
    // Compound assignment
    // -----------------------------------------------------------------------

    TensorVariable& operator+=( double val )            noexcept;
    TensorVariable& operator-=( double val )            noexcept;
    TensorVariable& operator*=( double val )            noexcept;
    TensorVariable& operator/=( double val )            noexcept;

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
    TensorVariable& operator=( const VectorVariable<3U>& )     noexcept;

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
    Equivalent to *this = *this * *this.
    Use HadamardSquared() for component-wise squaring.
    */
    TensorVariable& MatrixSquared() noexcept;

    /**
    Double contraction T:T = sum_ij T_ij^2.
    Returns the Frobenius norm squared as a scalar.
    This is a reduction — the tensor itself is not modified.
    */
    [[nodiscard]] double DoubleContraction() const noexcept;

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /**
    Returns the flag of diagonal element @p i (0, 1, or 2).
    Const version.
    */
    VARIABLE_FLAG  Flag( uint32_t i = 0 ) const noexcept;

    /**
    Returns a reference to the flag of diagonal element @p i.
    Non-const version.
    */
    VARIABLE_FLAG& Flag( uint32_t i = 0 )       noexcept;

    /**
    Returns true if all diagonal elements lie within [@p vmin, @p vmax].
    Off-diagonal elements are not checked (assumes diagonal dominance).
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

    /** Assigns @p f_00, @p f_11, @p f_22 to the diagonal elements. */
    void DiagonalValues( double f_00, double f_11, double f_22 ) noexcept;

    /** Assigns the first three elements of @p v to the diagonal. */
    void DiagonalValues( const std::vector<double>& v ) noexcept;

    /** Assigns the vector components to the diagonal. Flags are also copied. */
    void DiagonalValues( const VectorVariable<3U>& v ) noexcept;

    /** Assigns the vector to row @p i of the tensor. */
    void AssignToRow(    uint32_t i, VectorVariable<3U>& v ) noexcept;

    /** Assigns the vector to column @p j of the tensor. */
    void AssignToColumn( uint32_t j, VectorVariable<3U>& v ) noexcept;

    /** Returns row @p iRow as a VectorVariable. */
    VectorVariable<3U> Row(    uint32_t iRow ) const;

    /** Returns column @p iCol as a VectorVariable. */
    VectorVariable<3U> Column( uint32_t iCol ) const;

    // -----------------------------------------------------------------------
    // Eigenvalue / eigenvector methods
    // -----------------------------------------------------------------------

    /**
    Computes eigenvalues for a symmetric positive definite tensor using
    Cardano's formula. Returns eigenvalues in descending order.
    Returns false if the discriminant is positive (complex roots —
    should not occur for a symmetric matrix).
    */
    bool EigenValuesPositiveDefiniteSymmetricMatrix(
             double& eigenValue0,
             double& eigenValue1,
             double& eigenValue2 ) const;

    /**
    Delegates to EigenValuesPositiveDefiniteSymmetricMatrix.
    Returns eigenvalues into @p vecEigenvalues.
    */
    bool EigenValues( VectorVariable<3U>& vecEigenvalues ) const;

    /**
    Delegates to EigenValuesPositiveDefiniteSymmetricMatrix.
    Returns eigenvalues into @p vecEigenvalues.
    */
    bool EigenValues( std::vector<double>& vecEigenvalues ) const;

    /**
    Computes eigenvalues and eigenvectors for a symmetric tensor using
    Householder tridiagonalisation + symmetric tridiagonal QL algorithm.
    Off-diagonal elements are symmetrised by averaging.
    Eigenvalues are returned in descending order.
    Eigenvectors are stored as columns of @p eigenVecs.
    */
    bool EigenSymmetric( VectorVariable<3U>& eigenVals,
                         TensorVariable<3U>& eigenVecs ) const;

    /**
    Computes eigenvalues and eigenvectors for a weakly non-symmetric
    tensor by decomposing into symmetric and skew-symmetric parts.
    Returns true if ||W||_F / ||S||_F <= @p tolerance.
    Returns false and issues WARNING if asymmetry exceeds tolerance.

    @param tolerance  Maximum permitted asymmetry ratio. Default 1e-6.
    */
    bool EigenWeaklyNonSymmetric( VectorVariable<3U>& eigenVals,
                                  TensorVariable<3U>& eigenVecs,
                                  double              tolerance = 1.0e-6 ) const;

    /**
    Computes eigenvalues and eigenvectors assuming the tensor is
    symmetric, using Cardano's formula for eigenvalues and
    LU-backsubstitution for eigenvectors.
    If @p bNormalize is true, eigenvectors are normalised to unit length.
    Returns false if the matrix rank is zero.
    */
    bool Eigen( VectorVariable<3U>& evals,
                TensorVariable<3U>& evecs,
                bool                bNormalize ) const;

    // -----------------------------------------------------------------------
    // I/O
    // -----------------------------------------------------------------------

    /** Prompts the user to enter tensor values from the command line. */
    void In();

    /** Prints the tensor to stdout. */
    void Out() const;

    /** Reads the tensor from a binary file stream. */
    bool In(  std::fstream& fp );

    /** Writes the tensor to a binary file stream. */
    bool Out( std::fstream& fp ) const;

private:
    std::array<VARIABLE_FLAG, 3U>           flag;
    std::array<std::array<double, 3U>, 3U>  data;
};


// ============================================================================
//  Free functions
// ============================================================================

/**
Vector-matrix multiplication: result = vc^T * ts
(equivalent to ts^T * vc as a column vector).
*/
VectorVariable<3U> operator*( const VectorVariable<3U>& vc,
                               const TensorVariable<3U>& ts ) noexcept;

/** Point-matrix multiplication. */
Point<3U> operator*( const Point<3U>& vc,
                     const TensorVariable<3U>& ts ) noexcept;

/** Stream output operator. */
template<uint32_t dim>
std::ostream& operator<<( std::ostream& stream,
                           const TensorVariable<dim>& o );

/** Factory function for 2D tensors. */
TensorVariable<2U> makeTensor( VARIABLE_FLAG, VARIABLE_FLAG,
                               double, double,
                               double, double ) noexcept;

/** Factory function for 3D tensors. */
TensorVariable<3U> makeTensor( VARIABLE_FLAG, VARIABLE_FLAG, VARIABLE_FLAG,
                               double, double, double,
                               double, double, double,
                               double, double, double ) noexcept;


// ============================================================================
//  Inline definitions
// ============================================================================

inline TensorVariable<3U>::TensorVariable() noexcept
    : flag{ { ANY, ANY, ANY } },
      data{ { { std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN() },
              { std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN() },
              { std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN() } } }
{}

inline TensorVariable<3U>::TensorVariable(
    VARIABLE_FLAG f, double val ) noexcept
    : flag{ { f, f, f } },
      data{ { { val, 0., 0. },
              { 0., val, 0. },
              { 0., 0., val } } }
{}

inline TensorVariable<3U>::TensorVariable(
    VARIABLE_FLAG f,
    double v11, double v12, double v13,
    double v21, double v22, double v23,
    double v31, double v32, double v33 ) noexcept
    : flag{ { f, f, f } },
      data{ { { v11, v12, v13 },
              { v21, v22, v23 },
              { v31, v32, v33 } } }
{}

inline TensorVariable<3U>::TensorVariable(
    VARIABLE_FLAG f11, VARIABLE_FLAG f22, VARIABLE_FLAG f33,
    double v11, double v12, double v13,
    double v21, double v22, double v23,
    double v31, double v32, double v33 ) noexcept
    : flag{ { f11, f22, f33 } },
      data{ { { v11, v12, v13 },
              { v21, v22, v23 },
              { v31, v32, v33 } } }
{}

inline double& TensorVariable<3U>::operator()(
    uint32_t i, uint32_t j ) noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nTensorVariable<3U>::operator(): row violation i="
                  << i << "\n";
        return data[0][0];
    }
    if ( j >= 3U ) {
        std::cerr << "\nTensorVariable<3U>::operator(): col violation j="
                  << j << "\n";
        return data[0][0];
    }
#endif
    return data[i][j];
}

inline const double& TensorVariable<3U>::operator()(
    uint32_t i, uint32_t j ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nTensorVariable<3U>::operator(): row violation i="
                  << i << "\n";
        return data[0][0];
    }
    if ( j >= 3U ) {
        std::cerr << "\nTensorVariable<3U>::operator(): col violation j="
                  << j << "\n";
        return data[0][0];
    }
#endif
    return data[i][j];
}

inline void TensorVariable<3U>::Component(
    uint32_t k, double val ) noexcept
{
    assert( k < Size() );
    data[k / 3][k % 3] = val;
}

inline double TensorVariable<3U>::Component( uint32_t k ) const noexcept
{
    assert( k < Size() );
    return data[k / 3][k % 3];
}

inline void TensorVariable<3U>::Resize(
    uint32_t, double newValue ) noexcept
{
    data[0][0] = data[0][1] = data[0][2] =
    data[1][0] = data[1][1] = data[1][2] =
    data[2][0] = data[2][1] = data[2][2] = newValue;
}

inline VARIABLE_FLAG TensorVariable<3U>::Flag( uint32_t i ) const noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nTensorVariable<3U>::Flag(): violation i=" << i << "\n";
        return flag[0];
    }
#endif
    return flag[i];
}

inline VARIABLE_FLAG& TensorVariable<3U>::Flag( uint32_t i ) noexcept
{
#ifndef NDEBUG
    if ( i >= 3U ) {
        std::cerr << "\nTensorVariable<3U>::Flag(): violation i=" << i << "\n";
        return flag[0];
    }
#endif
    return flag[i];
}

inline bool TensorVariable<3U>::IsWithinRange(
    double vmin, double vmax ) const noexcept
{
    if ( data[0][0] < vmin || data[0][0] > vmax ) return false;
    if ( data[1][1] < vmin || data[1][1] > vmax ) return false;
    if ( data[2][2] < vmin || data[2][2] > vmax ) return false;
    return true;
}

inline bool TensorVariable<3U>::Has_NaN_Values() const noexcept
{
    return std::isnan( data[0][0] ) ||
           std::isnan( data[1][1] ) ||
           std::isnan( data[2][2] );
}

inline double TensorVariable<3U>::Trace() const noexcept
{
    return data[0][0] + data[1][1] + data[2][2];
}

inline double TensorVariable<3U>::Determinant() const noexcept
{
    double det =  data[0][0] * ( data[1][1] * data[2][2]
                               - data[2][1] * data[1][2] );
    det        -= data[0][1] * ( data[1][0] * data[2][2]
                               - data[2][0] * data[1][2] );
    det        += data[0][2] * ( data[1][0] * data[2][1]
                               - data[2][0] * data[1][1] );
    return det;
}

inline bool TensorVariable<3U>::operator==(
    const TensorVariable<3U>& ts ) const noexcept
{
    return ( data == ts.data && flag == ts.flag );
}

inline bool TensorVariable<3U>::operator!=(
    const TensorVariable<3U>& ts ) const noexcept
{
    return !( *this == ts );
}

inline bool TensorVariable<3U>::operator<(
    const TensorVariable<3U>& t ) const noexcept
{
    return ( this < &t );
}

// -----------------------------------------------------------------------
//  Squaring operations
// -----------------------------------------------------------------------

inline TensorVariable<3U>& TensorVariable<3U>::HadamardSquared() noexcept
{
    data[0][0] *= data[0][0];
    data[0][1] *= data[0][1];
    data[0][2] *= data[0][2];
    data[1][0] *= data[1][0];
    data[1][1] *= data[1][1];
    data[1][2] *= data[1][2];
    data[2][0] *= data[2][0];
    data[2][1] *= data[2][1];
    data[2][2] *= data[2][2];
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::MatrixSquared() noexcept
{
    return *this *= *this;
}

inline double TensorVariable<3U>::DoubleContraction() const noexcept
{
    return data[0][0]*data[0][0] + data[0][1]*data[0][1] + data[0][2]*data[0][2]
         + data[1][0]*data[1][0] + data[1][1]*data[1][1] + data[1][2]*data[1][2]
         + data[2][0]*data[2][0] + data[2][1]*data[2][1] + data[2][2]*data[2][2];
}

// -----------------------------------------------------------------------
//  Arithmetic — returning new tensor
// -----------------------------------------------------------------------

inline TensorVariable<3U> TensorVariable<3U>::operator+(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]+val, data[0][1]+val, data[0][2]+val,
        data[1][0]+val, data[1][1]+val, data[1][2]+val,
        data[2][0]+val, data[2][1]+val, data[2][2]+val );
}

inline TensorVariable<3U> TensorVariable<3U>::operator-(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]-val, data[0][1]-val, data[0][2]-val,
        data[1][0]-val, data[1][1]-val, data[1][2]-val,
        data[2][0]-val, data[2][1]-val, data[2][2]-val );
}

inline TensorVariable<3U> TensorVariable<3U>::operator*(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]*val, data[0][1]*val, data[0][2]*val,
        data[1][0]*val, data[1][1]*val, data[1][2]*val,
        data[2][0]*val, data[2][1]*val, data[2][2]*val );
}

inline TensorVariable<3U> TensorVariable<3U>::operator/(
    double val ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]/val, data[0][1]/val, data[0][2]/val,
        data[1][0]/val, data[1][1]/val, data[1][2]/val,
        data[2][0]/val, data[2][1]/val, data[2][2]/val );
}

inline TensorVariable<3U> TensorVariable<3U>::operator+(
    const TensorVariable<3U>& t ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]+t.data[0][0], data[0][1]+t.data[0][1], data[0][2]+t.data[0][2],
        data[1][0]+t.data[1][0], data[1][1]+t.data[1][1], data[1][2]+t.data[1][2],
        data[2][0]+t.data[2][0], data[2][1]+t.data[2][1], data[2][2]+t.data[2][2] );
}

inline TensorVariable<3U> TensorVariable<3U>::operator-(
    const TensorVariable<3U>& t ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]-t.data[0][0], data[0][1]-t.data[0][1], data[0][2]-t.data[0][2],
        data[1][0]-t.data[1][0], data[1][1]-t.data[1][1], data[1][2]-t.data[1][2],
        data[2][0]-t.data[2][0], data[2][1]-t.data[2][1], data[2][2]-t.data[2][2] );
}

inline TensorVariable<3U> TensorVariable<3U>::operator/(
    const TensorVariable<3U>& ts ) const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0]/ts.data[0][0], data[0][1]/ts.data[0][1], data[0][2]/ts.data[0][2],
        data[1][0]/ts.data[1][0], data[1][1]/ts.data[1][1], data[1][2]/ts.data[1][2],
        data[2][0]/ts.data[2][0], data[2][1]/ts.data[2][1], data[2][2]/ts.data[2][2] );
}

inline TensorVariable<3U> TensorVariable<3U>::operator*(
    const TensorVariable<3U>& ts ) const noexcept
{
    TensorVariable<3U> temp;
    // Unrolled matrix multiplication — row 0
    temp.data[0][0] = data[0][0]*ts.data[0][0] + data[0][1]*ts.data[1][0] + data[0][2]*ts.data[2][0];
    temp.data[0][1] = data[0][0]*ts.data[0][1] + data[0][1]*ts.data[1][1] + data[0][2]*ts.data[2][1];
    temp.data[0][2] = data[0][0]*ts.data[0][2] + data[0][1]*ts.data[1][2] + data[0][2]*ts.data[2][2];
    // row 1
    temp.data[1][0] = data[1][0]*ts.data[0][0] + data[1][1]*ts.data[1][0] + data[1][2]*ts.data[2][0];
    temp.data[1][1] = data[1][0]*ts.data[0][1] + data[1][1]*ts.data[1][1] + data[1][2]*ts.data[2][1];
    temp.data[1][2] = data[1][0]*ts.data[0][2] + data[1][1]*ts.data[1][2] + data[1][2]*ts.data[2][2];
    // row 2
    temp.data[2][0] = data[2][0]*ts.data[0][0] + data[2][1]*ts.data[1][0] + data[2][2]*ts.data[2][0];
    temp.data[2][1] = data[2][0]*ts.data[0][1] + data[2][1]*ts.data[1][1] + data[2][2]*ts.data[2][1];
    temp.data[2][2] = data[2][0]*ts.data[0][2] + data[2][1]*ts.data[1][2] + data[2][2]*ts.data[2][2];
    temp.flag[0] = flag[0];
    temp.flag[1] = flag[1];
    temp.flag[2] = flag[2];
    return temp;
}

inline VectorVariable<3U> TensorVariable<3U>::operator*(
    const VectorVariable<3U>& vc ) const noexcept
{
    return VectorVariable<3U>( vc.Flag(0), vc.Flag(1), vc.Flag(2),
        data[0][0]*vc[0] + data[0][1]*vc[1] + data[0][2]*vc[2],
        data[1][0]*vc[0] + data[1][1]*vc[1] + data[1][2]*vc[2],
        data[2][0]*vc[0] + data[2][1]*vc[1] + data[2][2]*vc[2] );
}

inline Point<3U> TensorVariable<3U>::operator*(
    const Point<3U>& v ) const noexcept
{
    return Point<3U>(
        data[0][0]*v[0] + data[0][1]*v[1] + data[0][2]*v[2],
        data[1][0]*v[0] + data[1][1]*v[1] + data[1][2]*v[2],
        data[2][0]*v[0] + data[2][1]*v[1] + data[2][2]*v[2] );
}

// -----------------------------------------------------------------------
//  Compound assignment
// -----------------------------------------------------------------------

inline TensorVariable<3U>& TensorVariable<3U>::operator+=(
    double val ) noexcept
{
    data[0][0]+=val; data[0][1]+=val; data[0][2]+=val;
    data[1][0]+=val; data[1][1]+=val; data[1][2]+=val;
    data[2][0]+=val; data[2][1]+=val; data[2][2]+=val;
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator-=(
    double val ) noexcept
{
    data[0][0]-=val; data[0][1]-=val; data[0][2]-=val;
    data[1][0]-=val; data[1][1]-=val; data[1][2]-=val;
    data[2][0]-=val; data[2][1]-=val; data[2][2]-=val;
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator*=(
    double val ) noexcept
{
    data[0][0]*=val; data[0][1]*=val; data[0][2]*=val;
    data[1][0]*=val; data[1][1]*=val; data[1][2]*=val;
    data[2][0]*=val; data[2][1]*=val; data[2][2]*=val;
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator/=(
    double val ) noexcept
{
    data[0][0]/=val; data[0][1]/=val; data[0][2]/=val;
    data[1][0]/=val; data[1][1]/=val; data[1][2]/=val;
    data[2][0]/=val; data[2][1]/=val; data[2][2]/=val;
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator+=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]+=sc(); data[0][1]+=sc(); data[0][2]+=sc();
    data[1][0]+=sc(); data[1][1]+=sc(); data[1][2]+=sc();
    data[2][0]+=sc(); data[2][1]+=sc(); data[2][2]+=sc();
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator-=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]-=sc(); data[0][1]-=sc(); data[0][2]-=sc();
    data[1][0]-=sc(); data[1][1]-=sc(); data[1][2]-=sc();
    data[2][0]-=sc(); data[2][1]-=sc(); data[2][2]-=sc();
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator*=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]*=sc(); data[0][1]*=sc(); data[0][2]*=sc();
    data[1][0]*=sc(); data[1][1]*=sc(); data[1][2]*=sc();
    data[2][0]*=sc(); data[2][1]*=sc(); data[2][2]*=sc();
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator/=(
    const ScalarVariable& sc ) noexcept
{
    data[0][0]/=sc(); data[0][1]/=sc(); data[0][2]/=sc();
    data[1][0]/=sc(); data[1][1]/=sc(); data[1][2]/=sc();
    data[2][0]/=sc(); data[2][1]/=sc(); data[2][2]/=sc();
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator+=(
    const TensorVariable<3U>& ts ) noexcept
{
    data[0][0]+=ts.data[0][0]; data[0][1]+=ts.data[0][1]; data[0][2]+=ts.data[0][2];
    data[1][0]+=ts.data[1][0]; data[1][1]+=ts.data[1][1]; data[1][2]+=ts.data[1][2];
    data[2][0]+=ts.data[2][0]; data[2][1]+=ts.data[2][1]; data[2][2]+=ts.data[2][2];
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator-=(
    const TensorVariable<3U>& ts ) noexcept
{
    data[0][0]-=ts.data[0][0]; data[0][1]-=ts.data[0][1]; data[0][2]-=ts.data[0][2];
    data[1][0]-=ts.data[1][0]; data[1][1]-=ts.data[1][1]; data[1][2]-=ts.data[1][2];
    data[2][0]-=ts.data[2][0]; data[2][1]-=ts.data[2][1]; data[2][2]-=ts.data[2][2];
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator/=(
    const TensorVariable<3U>& ts ) noexcept
{
    data[0][0]/=ts.data[0][0]; data[0][1]/=ts.data[0][1]; data[0][2]/=ts.data[0][2];
    data[1][0]/=ts.data[1][0]; data[1][1]/=ts.data[1][1]; data[1][2]/=ts.data[1][2];
    data[2][0]/=ts.data[2][0]; data[2][1]/=ts.data[2][1]; data[2][2]/=ts.data[2][2];
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator*=(
    const TensorVariable<3U>& ts ) noexcept
{
    TensorVariable<3U> temp;
    // Unrolled matrix multiplication — row 0
    temp.data[0][0] = data[0][0]*ts.data[0][0] + data[0][1]*ts.data[1][0] + data[0][2]*ts.data[2][0];
    temp.data[0][1] = data[0][0]*ts.data[0][1] + data[0][1]*ts.data[1][1] + data[0][2]*ts.data[2][1];
    temp.data[0][2] = data[0][0]*ts.data[0][2] + data[0][1]*ts.data[1][2] + data[0][2]*ts.data[2][2];
    // row 1
    temp.data[1][0] = data[1][0]*ts.data[0][0] + data[1][1]*ts.data[1][0] + data[1][2]*ts.data[2][0];
    temp.data[1][1] = data[1][0]*ts.data[0][1] + data[1][1]*ts.data[1][1] + data[1][2]*ts.data[2][1];
    temp.data[1][2] = data[1][0]*ts.data[0][2] + data[1][1]*ts.data[1][2] + data[1][2]*ts.data[2][2];
    // row 2
    temp.data[2][0] = data[2][0]*ts.data[0][0] + data[2][1]*ts.data[1][0] + data[2][2]*ts.data[2][0];
    temp.data[2][1] = data[2][0]*ts.data[0][1] + data[2][1]*ts.data[1][1] + data[2][2]*ts.data[2][1];
    temp.data[2][2] = data[2][0]*ts.data[0][2] + data[2][1]*ts.data[1][2] + data[2][2]*ts.data[2][2];
    temp.flag[0] = flag[0];
    temp.flag[1] = flag[1];
    temp.flag[2] = flag[2];
    return *this = std::move( temp );
}

// -----------------------------------------------------------------------
//  Assignment
// -----------------------------------------------------------------------

inline TensorVariable<3U>& TensorVariable<3U>::operator=(
    double val ) noexcept
{
    data[0][0]=val; data[0][1]=val; data[0][2]=val;
    data[1][0]=val; data[1][1]=val; data[1][2]=val;
    data[2][0]=val; data[2][1]=val; data[2][2]=val;
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator=(
    const ScalarVariable& sc ) noexcept
{
    flag[0] = flag[1] = flag[2] = sc.Flag();
    data[0][0]=sc(); data[0][1]=sc(); data[0][2]=sc();
    data[1][0]=sc(); data[1][1]=sc(); data[1][2]=sc();
    data[2][0]=sc(); data[2][1]=sc(); data[2][2]=sc();
    return *this;
}

inline TensorVariable<3U>& TensorVariable<3U>::operator=(
    const VectorVariable<3U>& vc ) noexcept
{
    flag[0] = vc.Flag(0);
    flag[1] = vc.Flag(1);
    flag[2] = vc.Flag(2);
    data[0][0] = vc(0); data[0][1] = 0.0;   data[0][2] = 0.0;
    data[1][0] = 0.0;   data[1][1] = vc(1); data[1][2] = 0.0;
    data[2][0] = 0.0;   data[2][1] = 0.0;   data[2][2] = vc(2);
    return *this;
}

// -----------------------------------------------------------------------
//  Matrix operations
// -----------------------------------------------------------------------

inline void TensorVariable<3U>::Identity() noexcept
{
    data[0][0] = 1.0; data[0][1] = 0.0; data[0][2] = 0.0;
    data[1][0] = 0.0; data[1][1] = 1.0; data[1][2] = 0.0;
    data[2][0] = 0.0; data[2][1] = 0.0; data[2][2] = 1.0;
}

inline TensorVariable<3U> TensorVariable<3U>::Transposed() const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
        data[0][0], data[1][0], data[2][0],
        data[0][1], data[1][1], data[2][1],
        data[0][2], data[1][2], data[2][2] );
}

inline TensorVariable<3U> TensorVariable<3U>::Adjoint() const noexcept
{
    return TensorVariable( flag[0], flag[1], flag[2],
         data[1][1]*data[2][2] - data[1][2]*data[2][1],
        -data[1][0]*data[2][2] + data[1][2]*data[2][0],
         data[1][0]*data[2][1] - data[2][0]*data[1][1],
        -data[0][1]*data[2][2] + data[0][2]*data[2][1],
         data[0][0]*data[2][2] - data[0][2]*data[2][0],
        -data[0][0]*data[2][1] + data[0][1]*data[2][0],
         data[0][1]*data[1][2] - data[0][2]*data[1][1],
        -data[0][0]*data[1][2] + data[0][2]*data[1][0],
         data[0][0]*data[1][1] - data[0][1]*data[1][0] );
}

inline TensorVariable<3U> TensorVariable<3U>::Inverse() const noexcept
{
    double det = Determinant();
    if ( det == 0.0 )
    {
        std::cerr << "\nTensorVariable<3U>::Inverse: Determinant = 0\n";
        return TensorVariable<3U>();
    }
    det = 1.0 / det;
    return TensorVariable( flag[0], flag[1], flag[2],
        det * ( data[1][1]*data[2][2] - data[1][2]*data[2][1] ),
        det * (-data[0][1]*data[2][2] + data[0][2]*data[2][1] ),
        det * ( data[0][1]*data[1][2] - data[0][2]*data[1][1] ),
        det * (-data[1][0]*data[2][2] + data[1][2]*data[2][0] ),
        det * ( data[0][0]*data[2][2] - data[0][2]*data[2][0] ),
        det * (-data[0][0]*data[1][2] + data[0][2]*data[1][0] ),
        det * ( data[1][0]*data[2][1] - data[2][0]*data[1][1] ),
        det * (-data[0][0]*data[2][1] + data[0][1]*data[2][0] ),
        det * ( data[0][0]*data[1][1] - data[0][1]*data[1][0] ) );
}

inline double TensorVariable<3U>::MinElement() const noexcept
{
    double me = data[0][0];
    if ( data[0][1] < me ) me = data[0][1];
    if ( data[1][0] < me ) me = data[1][0];
    if ( data[1][1] < me ) me = data[1][1];
    if ( data[0][2] < me ) me = data[0][2];
    if ( data[1][2] < me ) me = data[1][2];
    if ( data[2][0] < me ) me = data[2][0];
    if ( data[2][1] < me ) me = data[2][1];
    if ( data[2][2] < me ) me = data[2][2];
    return me;
}

inline double TensorVariable<3U>::MaxElement() const noexcept
{
    double me = data[0][0];
    if ( data[0][1] > me ) me = data[0][1];
    if ( data[1][0] > me ) me = data[1][0];
    if ( data[1][1] > me ) me = data[1][1];
    if ( data[0][2] > me ) me = data[0][2];
    if ( data[1][2] > me ) me = data[1][2];
    if ( data[2][0] > me ) me = data[2][0];
    if ( data[2][1] > me ) me = data[2][1];
    if ( data[2][2] > me ) me = data[2][2];
    return me;
}

// -----------------------------------------------------------------------
//  Diagonal operations
// -----------------------------------------------------------------------

inline void TensorVariable<3U>::DiagonalValues(
    double f_00, double f_11, double f_22 ) noexcept
{
    data[0][0] = f_00;
    data[1][1] = f_11;
    data[2][2] = f_22;
}

inline void TensorVariable<3U>::DiagonalValues(
    const std::vector<double>& v ) noexcept
{
    data[0][0] = v[0];
    data[1][1] = v[1];
    data[2][2] = v[2];
}

inline void TensorVariable<3U>::DiagonalValues(
    const VectorVariable<3U>& v ) noexcept
{
    data[0][0] = v[0];
    data[1][1] = v[1];
    data[2][2] = v[2];
    flag[0] = v.Flag(0);
    flag[1] = v.Flag(1);
    flag[2] = v.Flag(2);
}

inline void TensorVariable<3U>::AssignToRow(
    uint32_t iRow, VectorVariable<3U>& vc ) noexcept
{
    if      ( iRow == 0U ) flag[0] = vc.Flag(0);
    else if ( iRow == 1U ) flag[1] = vc.Flag(1);
    else                   flag[2] = vc.Flag(2);
    data[iRow][0] = vc[0];
    data[iRow][1] = vc[1];
    data[iRow][2] = vc[2];
}

inline void TensorVariable<3U>::AssignToColumn(
    uint32_t iCol, VectorVariable<3U>& vc ) noexcept
{
    if      ( iCol == 0U ) flag[0] = vc.Flag(0);
    else if ( iCol == 1U ) flag[1] = vc.Flag(1);
    else                   flag[2] = vc.Flag(2);
    data[0][iCol] = vc[0];
    data[1][iCol] = vc[1];
    data[2][iCol] = vc[2];
}

inline VectorVariable<3U> TensorVariable<3U>::Row(
    uint32_t iRow ) const
{
    return VectorVariable<3U>( flag[iRow], flag[iRow], flag[iRow],
                               data[iRow][0], data[iRow][1], data[iRow][2] );
}

inline VectorVariable<3U> TensorVariable<3U>::Column(
    uint32_t iCol ) const
{
    return VectorVariable<3U>( flag[iCol], flag[iCol], flag[iCol],
                               data[0][iCol], data[1][iCol], data[2][iCol] );
}

inline bool TensorVariable<3U>::EigenValues(
    VectorVariable<3U>& ev ) const
{
    return EigenValuesPositiveDefiniteSymmetricMatrix( ev(0), ev(1), ev(2) );
}

inline bool TensorVariable<3U>::EigenValues(
    std::vector<double>& ev ) const
{
    return EigenValuesPositiveDefiniteSymmetricMatrix( ev[0], ev[1], ev[2] );
}

// -----------------------------------------------------------------------
//  Free function inline definitions
// -----------------------------------------------------------------------

inline VectorVariable<3U> operator*(
    const VectorVariable<3U>& vc,
    const TensorVariable<3U>& ts ) noexcept
{
    return VectorVariable<3U>( vc.Flag(0), vc.Flag(1), vc.Flag(2),
        ts(0,0)*vc[0] + ts(1,0)*vc[1] + ts(2,0)*vc[2],
        ts(0,1)*vc[0] + ts(1,1)*vc[1] + ts(2,1)*vc[2],
        ts(0,2)*vc[0] + ts(1,2)*vc[1] + ts(2,2)*vc[2] );
}

inline Point<3U> operator*(
    const Point<3U>& v,
    const TensorVariable<3U>& ts ) noexcept
{
    return Point<3U>(
        ts(0,0)*v[0] + ts(1,0)*v[1] + ts(2,0)*v[2],
        ts(0,1)*v[0] + ts(1,1)*v[1] + ts(2,1)*v[2],
        ts(0,2)*v[0] + ts(1,2)*v[1] + ts(2,2)*v[2] );
}

inline TensorVariable<3U> makeTensor(
    VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3,
    double v11, double v12, double v13,
    double v21, double v22, double v23,
    double v31, double v32, double v33 ) noexcept
{
    return TensorVariable<3U>( f1, f2, f3,
                               v11, v12, v13,
                               v21, v22, v23,
                               v31, v32, v33 );
}

} // namespace csmp

#endif // CSMP_TENSOR_VARIABLE_H

