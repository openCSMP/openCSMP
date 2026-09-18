// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_DENSE_MATRIX_H
#define CSMP_DENSE_MATRIX_H

#define DENSE_MATRIX_USED_TOGETHER_WITH_CSMP

#include "CSMP_definitions.h"
#include "compareFloats.h"

#ifdef DENSE_MATRIX_USED_TOGETHER_WITH_CSMP
#include "Point.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#endif

namespace csmp {

/*! \file DenseMatrix.h
@addtogroup CSMPglobalEnums
@{ */

/// choice of fixed sizes 'mn_max'
enum CSMP_DMAT_SIZE { DM1=1, DM2=2, DM3=3U, DM4=4U, DM6=6U, DM12=12U, DM_MIN=36U, DM_MAX=81U };

/** @} */

/**
    Fixed-size full-storage matrix template as a workhorse for finite element integral computations
    and accumulation of contributions to the global solution matrix.

    Involves: M = dense matrix, v = vector, _T = transposed,

    Interoperable with basic CSMP variable types:
      T = tensor matrix(dim,dim), default vector is a column vector (v = { c0,c2..cn-1 },
      vc = vector(dim) )

    @author S. K. Matthai
    @author S. Geiger
    @author S. G. Roberts
    @date 2001
*/
template<uint32_t mn_max>
class DenseMatrix {
  public:
    DenseMatrix();
    DenseMatrix( const std::initializer_list<std::initializer_list<double>>& );
    DenseMatrix( uint32_t m, uint32_t n );
    DenseMatrix( uint32_t m, uint32_t n, double val );
    // rule of zero — compiler-generated copy/move/destructor are correct

    uint32_t Rows() const noexcept;
    uint32_t Cols() const noexcept;
    void Resize( uint32_t m, uint32_t n );

    double&       operator()( uint32_t m, uint32_t n )       noexcept;
    const double& operator()( uint32_t m, uint32_t n ) const noexcept;

    // ----------------------------------------------------------------------
    // Cross-Capacity Matrix Arithmetic 
    // Templated on 'other_max' so DenseMatrix<36> can add DenseMatrix<12>
    // ----------------------------------------------------------------------
    template<uint32_t other_max>
    DenseMatrix& operator+=( const DenseMatrix<other_max>& rhs ) noexcept;
    
    template<uint32_t other_max>
    DenseMatrix& operator-=( const DenseMatrix<other_max>& rhs ) noexcept;
    
    template<uint32_t other_max>
    DenseMatrix& operator*=( const DenseMatrix<other_max>& rhs );
    
    // Standard vector/scalar operations
    DenseMatrix& operator*=( const std::vector<double>& );
    DenseMatrix& operator*=( const double* );
    DenseMatrix& operator*=( double ) noexcept;

    void Identity() noexcept;
    void AssignToDiagonal( double val ) noexcept;
    void AssignToDiagonalAndZeroOffDiagonal( uint32_t diag_elmts, double val );
    void RowCondenseTo( std::vector<double>& ) const;
    void Zero()                          noexcept;
    void ZeroRow( uint32_t row )         noexcept;
    void ZeroCol( uint32_t col )         noexcept;
    void Fill( double val )              noexcept;
    void FillRow( uint32_t row, double val ) noexcept;
    void FillCol( uint32_t col, double val ) noexcept;

    // ----------------------------------------------------------------------
    // Cross-Capacity Matrix Transpositions
    // ----------------------------------------------------------------------
    template<uint32_t other_max>
    void Transposed( DenseMatrix<other_max>& res ) const noexcept;
    
    template<uint32_t other_max>
    void TransposedProduct( DenseMatrix<other_max>& res ) const;
    
    template<uint32_t other_max1, uint32_t other_max2>
    void MultiplyWithTransposedOf( const DenseMatrix<other_max1>& A, DenseMatrix<other_max2>& res ) const;
    
    template<uint32_t other_max1, uint32_t other_max2>
    void MultiplyTransposedOfWith( const DenseMatrix<other_max1>& A, DenseMatrix<other_max2>& res ) const;

    double RowSum( uint32_t row ) const;
    double ColSum( uint32_t col ) const;
    double NormL1()               const noexcept;
    double NormL_Infinity()       const noexcept;

    void In();
    void Out( long digits=5L ) const;

    /// read-only view of the underlying fixed-size storage array
    const std::array<std::array<double,mn_max>,mn_max>& Data() const noexcept { return data; }

#ifdef DENSE_MATRIX_USED_TOGETHER_WITH_CSMP
    // ----------------------------------------------------------------------
    // Unified CSMP Integration using template parameter 'dim'
    // ----------------------------------------------------------------------
    template<uint32_t dim>
    requires ( mn_max >= dim )
    void AssignRow( uint32_t i, const Point<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    void AssignCol( uint32_t j, const Point<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    void AssignToDiagonal( const Point<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    DenseMatrix& operator*=( const Point<dim>& ) noexcept;


    void AssignToDiagonal( uint32_t diag_elmts, const ScalarVariable& );
    DenseMatrix& operator*=( const ScalarVariable& ) noexcept;


    template<uint32_t dim>
    requires ( mn_max >= dim )
    void AssignRow( uint32_t i, const VectorVariable<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    void AssignCol( uint32_t j, const VectorVariable<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    void AssignToDiagonal( const VectorVariable<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    DenseMatrix& operator*=( const VectorVariable<dim>& );

    // Standard ArrayVariables (un-dimensioned)
    void AssignRow( uint32_t i, const ArrayVariable& );
    void AssignCol( uint32_t j, const ArrayVariable& );
    void AssignToDiagonal( const ArrayVariable& );
    DenseMatrix& operator*=( const ArrayVariable& );

    void AssignRow( uint32_t i, const FlaggedArrayVariable& );
    void AssignCol( uint32_t j, const FlaggedArrayVariable& );
    void AssignToDiagonal( const FlaggedArrayVariable& );
    DenseMatrix& operator*=( const FlaggedArrayVariable& );


    template<uint32_t dim>
    requires ( mn_max >= dim )
    void ExportTo( TensorVariable<dim>& ) const;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    DenseMatrix& operator=( const TensorVariable<dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    DenseMatrix& operator*=( const TensorVariable<dim>& );


    template<uint32_t dim>
    requires ( mn_max >= dim )
    DenseMatrix& operator=( const std::array<std::array<double,dim>,dim>& ) noexcept;

    template<uint32_t dim>
    requires ( mn_max >= dim )
    DenseMatrix& operator*=( const std::array<std::array<double,dim>,dim>& );

#endif // DENSE_MATRIX_USED_TOGETHER_WITH_CSMP 

  private:
    // Allow DenseMatrix of different sizes to access each other's private data
    template<uint32_t> friend class DenseMatrix;

    std::array<std::array<double,mn_max>,mn_max> data;
    uint32_t rows, cols;

    bool CheckRange( uint32_t m, uint32_t n, const char* originator ) const noexcept;
    
    // CheckSizes also needs to become a template if it inspects the other matrix
    template<uint32_t other_max>
    bool CheckSizes( const DenseMatrix<other_max>& mat, const char* originator ) const noexcept;
};


// =============================================================================
// Associated binary operators — declarations
// =============================================================================

/// same-type comparators
template<uint32_t mn_max>
bool operator==( const DenseMatrix<mn_max>&, const DenseMatrix<mn_max>& ) noexcept;

template<uint32_t mn_max>
bool operator!=( const DenseMatrix<mn_max>&, const DenseMatrix<mn_max>& ) noexcept;

/// cross-type comparators
template<uint32_t mn_max_A, uint32_t mn_max_B>
bool operator==( const DenseMatrix<mn_max_A>&, const DenseMatrix<mn_max_B>& ) noexcept;

template<uint32_t mn_max_A, uint32_t mn_max_B>
bool operator!=( const DenseMatrix<mn_max_A>&, const DenseMatrix<mn_max_B>& ) noexcept;

// =============================================================================
// Inline performance-critical member functions — defined outside the class body
// =============================================================================

template<uint32_t mn_max>
inline uint32_t DenseMatrix<mn_max>::Rows() const noexcept { return rows; }

template<uint32_t mn_max>
inline uint32_t DenseMatrix<mn_max>::Cols() const noexcept { return cols; }

template<uint32_t mn_max>
template<uint32_t other_max>
inline bool DenseMatrix<mn_max>::CheckSizes( const DenseMatrix<other_max>& mat, const char* ) const noexcept
 {
    // called only inside #ifdef DEBUG guards — assert() gives zero overhead
    // in release builds and a clear diagnostic message in debug builds
    assert( rows == mat.rows
            && "DenseMatrix::CheckSizes: row mismatch — see originator for context" );
    assert( cols == mat.cols
            && "DenseMatrix::CheckSizes: col mismatch — see originator for context" );
    return true;
 }

template<uint32_t mn_max>
inline bool DenseMatrix<mn_max>::CheckRange( uint32_t m, uint32_t n, const char* ) const noexcept
 {
    assert( m < rows
            && "DenseMatrix::CheckRange: row index out of range — see originator for context" );
    assert( n < cols
            && "DenseMatrix::CheckRange: col index out of range — see originator for context" );
    return true;
 }

template<uint32_t mn_max>
inline double& DenseMatrix<mn_max>::operator()( uint32_t m, uint32_t n ) noexcept
 {
#ifdef DEBUG
    CheckRange( m, n, "DenseMatrix::operator()" );
#endif
    return data[m][n];
 }

template<uint32_t mn_max>
inline const double& DenseMatrix<mn_max>::operator()( uint32_t m, uint32_t n ) const noexcept
 {
#ifdef DEBUG
    CheckRange( m, n, "DenseMatrix::operator() const" );
#endif
    return data[m][n];
 }

template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::Zero() noexcept
 {
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = 0.0;
 }

template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::Fill( double val ) noexcept
 {
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = val;
 }

template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::Identity() noexcept
 {
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = ( i == j ) ? 1.0 : 0.0;
 }

template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::AssignToDiagonal( double val ) noexcept
 {
    for ( uint32_t i{0U}; i < rows; ++i )
      data[i][i] = val;
 }

template<uint32_t mn_max>
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( double val ) noexcept
 {
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] *= val;
    return *this;
 }

template<uint32_t mn_max>
template<uint32_t mn_max_B>
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator+=( const DenseMatrix<mn_max_B>& mat ) noexcept
 {
#ifdef DEBUG
    CheckSizes( mat, "DenseMatrix::operator+=" );
#endif
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] += mat(i, j);
    return *this;
 }

template<uint32_t mn_max>
template<uint32_t mn_max_B>
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator-=( const DenseMatrix<mn_max_B>& mat ) noexcept
 {
#ifdef DEBUG
    CheckSizes( mat, "DenseMatrix::operator-=" );
#endif
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] -= mat(i, j);
    return *this;
 }

template<uint32_t mn_max>
template<uint32_t mn_max_B>
inline void DenseMatrix<mn_max>::Transposed( DenseMatrix<mn_max_B>& M ) const noexcept
 {
    M.Resize( cols, rows );
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        M(j, i) = data[i][j];
 }
 

template<uint32_t mn_max>
template<uint32_t other_max>
inline void DenseMatrix<mn_max>::TransposedProduct( DenseMatrix<other_max>& res ) const
 {
    res.Resize( cols, cols );
    for ( uint32_t i{0U}; i < cols; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        {
          double sum{0.0};
          for ( uint32_t k{0U}; k < rows; ++k )
            sum += data[k][i] * data[k][j];
          res.data[j][i] = sum;
        }
 }
 
 
 
template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::ZeroRow( uint32_t row ) noexcept
 {
    if ( row < rows )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[row][j] = 0.0;
 }



template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::ZeroCol( uint32_t col ) noexcept
 {
    if ( col < cols )
      for ( uint32_t i{0U}; i < rows; ++i )
        data[i][col] = 0.0;
 }



template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::FillRow( uint32_t row, double val ) noexcept
 {
    if ( row < rows )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[row][j] = val;
 }



template<uint32_t mn_max>
inline void DenseMatrix<mn_max>::FillCol( uint32_t col, double val ) noexcept
 {
    if ( col < cols )
      for ( uint32_t i{0U}; i < rows; ++i )
        data[i][col] = val;
 }


// =============================================================================
// Cross-Capacity Matrix Arithmetic 
// =============================================================================

// ----------------------------------------------------------------------
// operator*=
// ----------------------------------------------------------------------
template<uint32_t mn_max>
template<uint32_t other_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const DenseMatrix<other_max>& rhs )
 {
    if ( cols != rhs.rows )
      {
        std::cerr << "\nDenseMatrix<" << mn_max << ">::operator*=: incompatible sizes" << std::endl;
        throw std::length_error("DenseMatrix::operator*=: incompatible sizes");
      }
    
    // Ensure the left-hand matrix has enough capacity to hold the new columns!
    if ( rhs.cols > mn_max )
      {
        std::cerr << "\nDenseMatrix<" << mn_max << ">::operator*=: capacity exceeded" << std::endl;
        throw std::length_error("DenseMatrix::operator*=: right-hand cols exceeds left-hand capacity");
      }
    
    DenseMatrix<mn_max> temp;
    temp.Resize( rows, rhs.cols );
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < rhs.cols; ++j )
        {
          double sum{0.0};
          for ( uint32_t k{0U}; k < cols; ++k )
            sum += data[i][k] * rhs.data[k][j];
          temp.data[i][j] = sum;
        }

    rows = temp.rows;
    cols = temp.cols;
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = temp.data[i][j];

    return *this;
 }
 


// =============================================================================
// std::array Integration with Template-Head Constraints
// =============================================================================

template<uint32_t mn_max>
template<uint32_t dim> requires ( mn_max >= dim )
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator=( const std::array<std::array<double, dim>, dim>& arr ) noexcept
 {
    rows = cols = dim;
    if constexpr ( dim == 1U ) {
        data[0][0] = arr[0][0];
    } else if constexpr ( dim == 2U ) {
        data[0][0] = arr[0][0];  data[0][1] = arr[0][1];
        data[1][0] = arr[1][0];  data[1][1] = arr[1][1];
    } else if constexpr ( dim == 3U ) {
        data[0][0] = arr[0][0];  data[0][1] = arr[0][1];  data[0][2] = arr[0][2];
        data[1][0] = arr[1][0];  data[1][1] = arr[1][1];  data[1][2] = arr[1][2];
        data[2][0] = arr[2][0];  data[2][1] = arr[2][1];  data[2][2] = arr[2][2];
    } else {
        for ( uint32_t i{0U}; i < dim; ++i )
          for ( uint32_t j{0U}; j < dim; ++j )
            data[i][j] = arr[i][j];
    }
    return *this;
 }

template<uint32_t mn_max>
template<uint32_t dim> requires ( mn_max >= dim )
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const std::array<std::array<double, dim>, dim>& arr )
 {
    if ( cols != dim )
      throw std::range_error("DenseMatrix::operator*=(array): cols != dim");

    if constexpr ( dim == 1U ) {
        if ( rows != cols )
          for ( uint32_t i{0U}; i < rows; ++i ) data[i][0] *= arr[0][0];
        else
          data[0][0] *= arr[0][0];
        cols = 1U;
        return *this;
    } else {
        if ( rows != cols )
          {
            DenseMatrix<mn_max> temp( rows, dim );
            for ( uint32_t i{0U}; i < rows; ++i )
              for ( uint32_t j{0U}; j < dim; ++j )
                {
                  temp.data[i][j] = 0.0;
                  for ( uint32_t k{0U}; k < dim; ++k )
                    temp.data[i][j] += data[i][k] * arr[k][j];
                }
            return *this = temp;
          }

        // Unrolled square cases with general loop fallback
        if constexpr ( dim == 2U ) {
            const double t00 = data[0][0]*arr[0][0] + data[0][1]*arr[1][0];
            const double t01 = data[0][0]*arr[0][1] + data[0][1]*arr[1][1];
            const double t10 = data[1][0]*arr[0][0] + data[1][1]*arr[1][0];
            const double t11 = data[1][0]*arr[0][1] + data[1][1]*arr[1][1];
            data[0][0] = t00;  data[0][1] = t01;
            data[1][0] = t10;  data[1][1] = t11;
        } else if constexpr ( dim == 3U ) {
            const double t00 = data[0][0]*arr[0][0] + data[0][1]*arr[1][0] + data[0][2]*arr[2][0];
            const double t01 = data[0][0]*arr[0][1] + data[0][1]*arr[1][1] + data[0][2]*arr[2][1];
            const double t02 = data[0][0]*arr[0][2] + data[0][1]*arr[1][2] + data[0][2]*arr[2][2];
            const double t10 = data[1][0]*arr[0][0] + data[1][1]*arr[1][0] + data[1][2]*arr[2][0];
            const double t11 = data[1][0]*arr[0][1] + data[1][1]*arr[1][1] + data[1][2]*arr[2][1];
            const double t12 = data[1][0]*arr[0][2] + data[1][1]*arr[1][2] + data[1][2]*arr[2][2];
            const double t20 = data[2][0]*arr[0][0] + data[2][1]*arr[1][0] + data[2][2]*arr[2][0];
            const double t21 = data[2][0]*arr[0][1] + data[2][1]*arr[1][1] + data[2][2]*arr[2][1];
            const double t22 = data[2][0]*arr[0][2] + data[2][1]*arr[1][2] + data[2][2]*arr[2][2];
            data[0][0] = t00;  data[0][1] = t01;  data[0][2] = t02;
            data[1][0] = t10;  data[1][1] = t11;  data[1][2] = t12;
            data[2][0] = t20;  data[2][1] = t21;  data[2][2] = t22;
        } else {
            DenseMatrix<mn_max> temp( rows, dim );
            for ( uint32_t i{0U}; i < rows; ++i )
              for ( uint32_t j{0U}; j < dim; ++j )
                {
                  temp.data[i][j] = 0.0;
                  for ( uint32_t k{0U}; k < dim; ++k )
                    temp.data[i][j] += data[i][k] * arr[k][j];
                }
            return *this = temp;
        }
        return *this;
    }
 }


// =============================================================================
// CSMP Integration Block - Manually Unrolled for Performance
// =============================================================================
#ifdef DENSE_MATRIX_USED_TOGETHER_WITH_CSMP

template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::AssignRow( uint32_t i, const Point<dim>& pt ) noexcept
 {
    cols = dim;
    if constexpr ( dim == 1U ) {
        data[i][0] = pt[0];
    } else if constexpr ( dim == 2U ) {
        data[i][0] = pt[0];  data[i][1] = pt[1];
    } else if constexpr ( dim == 3U ) {
        data[i][0] = pt[0];  data[i][1] = pt[1];  data[i][2] = pt[2];
    }
 }



template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::AssignCol( uint32_t j, const Point<dim>& pt ) noexcept
 {
    rows = dim;
    if constexpr ( dim == 1U ) {
        data[0][j] = pt[0];
    } else if constexpr ( dim == 2U ) {
        data[0][j] = pt[0];  data[1][j] = pt[1];
    } else if constexpr ( dim == 3U ) {
        data[0][j] = pt[0];  data[1][j] = pt[1];  data[2][j] = pt[2];
    }
 }



template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::AssignToDiagonal( const Point<dim>& pt ) noexcept
 {
    rows = cols = dim;
    if constexpr ( dim == 1U ) {
        data[0][0] = pt[0];
    } else if constexpr ( dim == 2U ) {
        data[0][0] = pt[0];  data[0][1] = 0.0;
        data[1][0] = 0.0;    data[1][1] = pt[1];
    } else if constexpr ( dim == 3U ) {
        data[0][0] = pt[0];  data[0][1] = 0.0;    data[0][2] = 0.0;
        data[1][0] = 0.0;    data[1][1] = pt[1];  data[1][2] = 0.0;
        data[2][0] = 0.0;    data[2][1] = 0.0;    data[2][2] = pt[2];
    }
 }




template<uint32_t mn_max>
template<uint32_t dim> requires ( mn_max >= dim )
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const Point<dim>& pt ) noexcept
{
    for ( uint32_t i{0U}; i < rows; ++i ) {
        double row_sum = 0.0;
        
        // 1. Calculate the dot product for the transformation
        if constexpr ( dim == 1U ) {
            row_sum = data[i][0] * pt[0];
        } else if constexpr ( dim == 2U ) {
            row_sum = (data[i][0] * pt[0]) + (data[i][1] * pt[1]);
        } else if constexpr ( dim == 3U ) {
            row_sum = (data[i][0] * pt[0]) + (data[i][1] * pt[1]) + (data[i][2] * pt[2]);
        }
        
        // 2. Store the resulting coordinate in the first column
        data[i][0] = row_sum;
        
        // 3. Zero out the remaining columns
        if constexpr ( dim > 1U ) data[i][1] = 0.0;
        if constexpr ( dim > 2U ) data[i][2] = 0.0;
    }
    cols = 1U; // <-- CRITICAL: keep state consistent!
    return *this;
}



template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::AssignRow( uint32_t i, const VectorVariable<dim>& vc ) noexcept
 {
    cols = dim;
    if constexpr ( dim == 1U ) {
        data[i][0] = vc[0];
    } else if constexpr ( dim == 2U ) {
        data[i][0] = vc[0];  data[i][1] = vc[1];
    } else if constexpr ( dim == 3U ) {
        data[i][0] = vc[0];  data[i][1] = vc[1];  data[i][2] = vc[2];
    }
 }



template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::AssignCol( uint32_t j, const VectorVariable<dim>& vc ) noexcept
 {
    rows = dim;
    if constexpr ( dim == 1U ) {
        data[0][j] = vc[0];
    } else if constexpr ( dim == 2U ) {
        data[0][j] = vc[0];  data[1][j] = vc[1];
    } else if constexpr ( dim == 3U ) {
        data[0][j] = vc[0];  data[1][j] = vc[1];  data[2][j] = vc[2];
    }
 }



template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::AssignToDiagonal( const VectorVariable<dim>& vc ) noexcept
 {
    rows = cols = dim;
    if constexpr ( dim == 1U ) {
        data[0][0] = vc[0];
    } else if constexpr ( dim == 2U ) {
        data[0][0] = vc[0];  data[0][1] = 0.0;
        data[1][0] = 0.0;    data[1][1] = vc[1];
    } else if constexpr ( dim == 3U ) {
        data[0][0] = vc[0];  data[0][1] = 0.0;    data[0][2] = 0.0;
        data[1][0] = 0.0;    data[1][1] = vc[1];  data[1][2] = 0.0;
        data[2][0] = 0.0;    data[2][1] = 0.0;    data[2][2] = vc[2];
    }
 }



template<uint32_t mn_max>
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const ScalarVariable& sc ) noexcept
 {
    // Note: Kept strictly on mn_max since a true math scalar has no spatial <dim>
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] *= sc();
    return *this;
 }



template<uint32_t mn_max>
template<uint32_t dim> requires ( mn_max >= dim )
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const VectorVariable<dim>& vc )
{
    for ( uint32_t i{0U}; i < rows; ++i ) {
        double row_sum = 0.0;
        
        // 1. Calculate the dot product for the current row
        if constexpr ( dim == 1U ) {
            row_sum = data[i][0] * vc[0];
        } else if constexpr ( dim == 2U ) {
            row_sum = (data[i][0] * vc[0]) + (data[i][1] * vc[1]);
        } else if constexpr ( dim == 3U ) {
            row_sum = (data[i][0] * vc[0]) + (data[i][1] * vc[1]) + (data[i][2] * vc[2]);
        }
        
        // 2. Store the resulting vector component in the first column
        data[i][0] = row_sum;
        
        // 3. (Optional but recommended) Zero out the remaining columns 
        // to prevent stale data from lingering in the matrix structure.
        if constexpr ( dim > 1U ) data[i][1] = 0.0;
        if constexpr ( dim > 2U ) data[i][2] = 0.0;
    }
    cols = 1U; // <-- CRITICAL: keep state consistent!
    return *this;
}



template<uint32_t mn_max>
template<uint32_t dim>
requires ( mn_max >= dim )
inline DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator=( const TensorVariable<dim>& ts ) noexcept
 {
    rows = cols = dim;
    if constexpr ( dim == 1U ) {
        data[0][0] = ts(0,0);
    } else if constexpr ( dim == 2U ) {
        data[0][0] = ts(0,0);  data[0][1] = ts(0,1);
        data[1][0] = ts(1,0);  data[1][1] = ts(1,1);
    } else if constexpr ( dim == 3U ) {
        data[0][0] = ts(0,0);  data[0][1] = ts(0,1);  data[0][2] = ts(0,2);
        data[1][0] = ts(1,0);  data[1][1] = ts(1,1);  data[1][2] = ts(1,2);
        data[2][0] = ts(2,0);  data[2][1] = ts(2,1);  data[2][2] = ts(2,2);
    }
    return *this;
 }
 
// =============================================================================
// ExportTo & Operator*= for TensorVariable
// =============================================================================

template<uint32_t mn_max>
template<uint32_t dim> requires ( mn_max >= dim )
inline void DenseMatrix<mn_max>::ExportTo( TensorVariable<dim>& ts ) const
 {
    if ( rows != dim || rows != cols )
      throw std::range_error("DenseMatrix::ExportTo: wrong matrix size");

    if constexpr ( dim == 1U ) {
        ts(0,0) = data[0][0];
    } else if constexpr ( dim == 2U ) {
        ts(0,0) = data[0][0];  ts(0,1) = data[0][1];
        ts(1,0) = data[1][0];  ts(1,1) = data[1][1];
    } else if constexpr ( dim == 3U ) {
        ts(0,0) = data[0][0];  ts(0,1) = data[0][1];  ts(0,2) = data[0][2];
        ts(1,0) = data[1][0];  ts(1,1) = data[1][1];  ts(1,2) = data[1][2];
        ts(2,0) = data[2][0];  ts(2,1) = data[2][1];  ts(2,2) = data[2][2];
    }
 }

template<uint32_t mn_max>
template<uint32_t dim> requires ( mn_max >= dim )
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const TensorVariable<dim>& ts )
 {
    if ( cols != dim )
      throw std::range_error("DenseMatrix::operator*=(TensorVariable): cols != dim");

    if constexpr ( dim == 1U ) {
        if ( rows != cols )
          for ( uint32_t i{0U}; i < rows; ++i ) data[i][0] *= ts(0,0);
        else
          data[0][0] *= ts(0,0);
        cols = 1U;
        return *this;
    } else {
        if ( rows != cols )
          {
            DenseMatrix<mn_max> temp( rows, dim );
            for ( uint32_t i{0U}; i < rows; ++i )
              for ( uint32_t j{0U}; j < dim; ++j )
                {
                  temp.data[i][j] = 0.0;
                  for ( uint32_t k{0U}; k < dim; ++k )
                    temp.data[i][j] += data[i][k] * ts(k,j);
                }
            return *this = temp;
          }

        // unrolled square case
        if constexpr ( dim == 2U ) {
            const double t00 = data[0][0]*ts(0,0) + data[0][1]*ts(1,0);
            const double t01 = data[0][0]*ts(0,1) + data[0][1]*ts(1,1);
            const double t10 = data[1][0]*ts(0,0) + data[1][1]*ts(1,0);
            const double t11 = data[1][0]*ts(0,1) + data[1][1]*ts(1,1);
            data[0][0] = t00;  data[0][1] = t01;
            data[1][0] = t10;  data[1][1] = t11;
        } else if constexpr ( dim == 3U ) {
            const double t00 = data[0][0]*ts(0,0) + data[0][1]*ts(1,0) + data[0][2]*ts(2,0);
            const double t01 = data[0][0]*ts(0,1) + data[0][1]*ts(1,1) + data[0][2]*ts(2,1);
            const double t02 = data[0][0]*ts(0,2) + data[0][1]*ts(1,2) + data[0][2]*ts(2,2);
            const double t10 = data[1][0]*ts(0,0) + data[1][1]*ts(1,0) + data[1][2]*ts(2,0);
            const double t11 = data[1][0]*ts(0,1) + data[1][1]*ts(1,1) + data[1][2]*ts(2,1);
            const double t12 = data[1][0]*ts(0,2) + data[1][1]*ts(1,2) + data[1][2]*ts(2,2);
            const double t20 = data[2][0]*ts(0,0) + data[2][1]*ts(1,0) + data[2][2]*ts(2,0);
            const double t21 = data[2][0]*ts(0,1) + data[2][1]*ts(1,1) + data[2][2]*ts(2,1);
            const double t22 = data[2][0]*ts(0,2) + data[2][1]*ts(1,2) + data[2][2]*ts(2,2);
            data[0][0] = t00;  data[0][1] = t01;  data[0][2] = t02;
            data[1][0] = t10;  data[1][1] = t11;  data[1][2] = t12;
            data[2][0] = t20;  data[2][1] = t21;  data[2][2] = t22;
        }
        return *this;
    }
 }


#endif // DENSE_MATRIX_USED_TOGETHER_WITH_CSMP


// =============================================================================
// Free operator implementations
// =============================================================================

template<uint32_t mn_max>
inline bool operator==( const DenseMatrix<mn_max>& MA, const DenseMatrix<mn_max>& MB ) noexcept
 {
    if ( &MA == &MB ) return true;
    if ( MA.Rows() != MB.Rows() || MA.Cols() != MB.Cols() ) return false;
    for ( uint32_t i{0U}; i < MA.Rows(); ++i )
      for ( uint32_t j{0U}; j < MA.Cols(); ++j )
        if ( !essentiallyEqual( MA(i,j), MB(i,j) ) ) return false;
    return true;
 }

template<uint32_t mn_max>
inline bool operator!=( const DenseMatrix<mn_max>& MA, const DenseMatrix<mn_max>& MB ) noexcept
 { return !( MA == MB ); }

template<uint32_t mn_max_A, uint32_t mn_max_B>
inline bool operator==( const DenseMatrix<mn_max_A>& MA, const DenseMatrix<mn_max_B>& MB ) noexcept
 {
    if ( MA.Rows() != MB.Rows() || MA.Cols() != MB.Cols() ) return false;
    for ( uint32_t i{0U}; i < MA.Rows(); ++i )
      for ( uint32_t j{0U}; j < MA.Cols(); ++j )
        if ( !essentiallyEqual( MA(i,j), MB(i,j) ) ) return false;
    return true;
 }

template<uint32_t mn_max_A, uint32_t mn_max_B>
inline bool operator!=( const DenseMatrix<mn_max_A>& MA, const DenseMatrix<mn_max_B>& MB ) noexcept
 { return !( MA == MB ); }


template<uint32_t mn_max>
inline DenseMatrix<mn_max> operator+( const DenseMatrix<mn_max>& a, const DenseMatrix<mn_max>& b )
 {
    DenseMatrix<mn_max> temp(a);
    temp += b;
    return temp;
 }


template<uint32_t mn_max>
inline DenseMatrix<mn_max> operator-( const DenseMatrix<mn_max>& a, const DenseMatrix<mn_max>& b )
 {
    DenseMatrix<mn_max> temp(a);
    temp -= b;
    return temp;
 }


template<uint32_t mn_max>
inline DenseMatrix<mn_max> operator*( const DenseMatrix<mn_max>& a, const DenseMatrix<mn_max>& b )
 {
#if defined(DEBUG) && defined(CSMP_DENSE_MATRIX_DEBUG)
    if ( a.Cols() != b.Rows() )
      {
        std::cerr <<"\nDenseMatrix<"<< mn_max <<"> operator*: incompatible sizes" << std::endl;
        a.Out(3L);
        b.Out(3L);
        throw std::length_error("DenseMatrix operator*: incompatible sizes");
      }
#endif
    DenseMatrix<mn_max> temp( a.Rows(), b.Cols() );
    for ( uint32_t i{0U}; i < a.Rows(); ++i )
      for ( uint32_t j{0U}; j < b.Cols(); ++j )
        {
          temp(i,j) = 0.0;
          for ( uint32_t k{0U}; k < b.Rows(); ++k )
            temp(i,j) += a(i,k) * b(k,j);
        }
    return temp;
 }


template<uint32_t mn_max_A, uint32_t mn_max_B>
inline DenseMatrix<(mn_max_A > mn_max_B ? mn_max_A : mn_max_B)>
  operator*( const DenseMatrix<mn_max_A>& A, const DenseMatrix<mn_max_B>& B )
 {
    if ( A.Cols() != B.Rows() )
      throw std::range_error("DenseMatrix operator*: A.cols != B.rows");

    constexpr uint32_t mn_max_C = ( mn_max_A > mn_max_B ? mn_max_A : mn_max_B );
    DenseMatrix<mn_max_C> C( A.Rows(), B.Cols() );
    C.Zero();
    for ( uint32_t i{0U}; i < A.Rows(); ++i )
      for ( uint32_t j{0U}; j < B.Cols(); ++j )
        for ( uint32_t k{0U}; k < A.Cols(); ++k )
          C(i,j) += A(i,k) * B(k,j);
    return C;
 }


template<uint32_t mn_max>
inline  std::vector<double> operator*( const DenseMatrix<mn_max>& mat, const std::vector<double>& vec )
 {
    assert( mat.Cols() == vec.size() );
    std::vector<double> temp( mat.Rows(), 0.0 );
    for ( uint32_t i{0U}; i < mat.Rows(); ++i )
      for ( uint32_t j{0U}; j < mat.Cols(); ++j )
        temp[i] += mat(i,j) * vec[j];
    return temp;
 }


template<uint32_t mn_max>
inline DenseMatrix<mn_max> operator*( const std::vector<double>& vec, const DenseMatrix<mn_max>& mat )
 {
    if ( vec.size() != mat.Rows() )
      {
        std::cerr <<"\noperator*(vector, DenseMatrix): incompatible sizes" << std::endl;
        throw std::length_error("operator*(vector, DenseMatrix): incompatible sizes");
      }
    DenseMatrix<mn_max> temp( vec.size(), mat.Cols() );
    for ( uint32_t i{0U}; i < vec.size(); ++i )
      for ( uint32_t j{0U}; j < mat.Cols(); ++j )
        {
          temp(i,j) = 0.0;
          for ( uint32_t k{0U}; k < mat.Rows(); ++k )
            temp(i,j) += vec[k] * mat(k,j);
        }
    return temp;
 }
 
} // namespace csmp

#endif // CSMP_DENSE_MATRIX_H

