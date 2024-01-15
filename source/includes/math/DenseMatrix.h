#ifndef CSMP_DENSE_MATRIX_H
#define CSMP_DENSE_MATRIX_H

#include <stdexcept>
#include <typeinfo>
#include <array>
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

///choice of fixed sizes 'mn_max'
enum CSMP_DMAT_SIZE { DM3=3U, DM4=4U, DM6=6U, DM12=12U, DM_MIN=36U, DM_MAX=81U };

/** @} */

/** 
    Fixed-size full-storage matrix template as a workhorse for finite element integral computations
    and accumulation of contributions to the global solution matrix.
 
    Involves: M = dense matrix, v = vector, _T = transposed,
 
    Interoperable with with basic CSMP variable types
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
    DenseMatrix( uint32_t m, uint32_t n );
    DenseMatrix( uint32_t m, uint32_t n, double val );
    DenseMatrix( const DenseMatrix& );
    DenseMatrix( DenseMatrix&& )= default;
    ~DenseMatrix();
    uint32_t Rows() const;
    uint32_t Cols() const;
    void     Resize( uint32_t m, uint32_t n );
    /// accessors M(i,j)
    double&       operator()( uint32_t m, uint32_t n );
    const double& operator()( uint32_t m, uint32_t n ) const;
    /// assignment
    DenseMatrix& operator=( const DenseMatrix& );
    DenseMatrix& operator=( DenseMatrix&& ) = default;
    DenseMatrix& operator=( double );
    DenseMatrix& operator+=( const DenseMatrix& );
    DenseMatrix& operator-=( const DenseMatrix& );
    DenseMatrix  operator+( const DenseMatrix& ) const;
    DenseMatrix  operator-( const DenseMatrix& ) const;
    /// matrix - matrix multiplication -> M(A.rows,B.cols) (NB: creates temporary M)
    DenseMatrix& operator*=( const DenseMatrix& );
    /// matrix vector multiplication
    DenseMatrix& operator*=( const std::vector<double>& );
    /// rhs is a C-array
    DenseMatrix& operator*=( const double* );
    DenseMatrix& operator*=( double );
    /// C = A B, matrix - matrix multiplication -> M(A.rows,B.cols) returns temporary matrix
    DenseMatrix  operator*( const DenseMatrix& ) const;
    /// creates identify matrix (Kronecker-Delta) with ones in the diagonal and zeros everywhere else
    void Identity();
    /// assigns 'value' to the i=j elements  o fmatrix
    void AssignToDiagonal( double val );
    /// resizes matrix to 'diag_elmts', assigns 'value' to the diagonal elements matrix diagonal while zeroing out all other values
    void AssignToDiagonalAndZeroOffDiagonal( uint32_t diag_elmts, double val );
    /// returns vec = Mat * unity vector
    void RowCondenseTo( std::vector<double>& ) const;
    void Zero();
    void ZeroRow( uint32_t row );
    void ZeroCol( uint32_t col );
    void Fill( double val );
    void FillRow( uint32_t row, double val );
    void FillCol( uint32_t col, double val );
    /// RES = M^T into its argument
    void Transposed( DenseMatrix& ) const;
    /// RES = M^T M
    void TransposedProduct( DenseMatrix& ) const;
    /// RES = A B^T
    void MultiplyWithTransposedOf( const DenseMatrix&, DenseMatrix& ) const;
    /// RES = A^T B
    void MultiplyTransposedOfWith( const DenseMatrix&, DenseMatrix& ) const;
    double   RowSum( uint32_t row ) const;
    double   ColSum( uint32_t col ) const;
    /// unscaled L1 matrix norm
    double   NormL1() const;
    /// L_infinity matrix norm
    double   NormL_Infinity() const;

    void In();
    void Out( long digits=5L ) const;

#ifdef DENSE_MATRIX_USED_TOGETHER_WITH_CSMP
    /// assignment of point coordinates to matrix rows or columns
    void AssignRow( uint32_t i, const Point<1U>& );
    void AssignRow( uint32_t i, const Point<2U>& );
    void AssignRow( uint32_t i, const Point<3U>& );
    void AssignCol( uint32_t j, const Point<1U>& );
    void AssignCol( uint32_t j, const Point<2U>& );
    void AssignCol( uint32_t j, const Point<3U>& );
    void AssignToDiagonal( const Point<1U>& );
    void AssignToDiagonal( const Point<2U>& );
    void AssignToDiagonal( const Point<3U>& );
    /// matrix point multiplication
    DenseMatrix& operator*=( const Point<1U>& );
    DenseMatrix& operator*=( const Point<2U>& );
    DenseMatrix& operator*=( const Point<3U>& );

    void AssignToDiagonal( uint32_t diag_elmts, const ScalarVariable& );
    DenseMatrix& operator*=( const ScalarVariable& );

    void AssignRow( uint32_t i, const VectorVariable<1U>& );
    void AssignRow( uint32_t i, const VectorVariable<2U>& );
    void AssignRow( uint32_t i, const VectorVariable<3U>& );
    void AssignCol( uint32_t j, const VectorVariable<1U>& );
    void AssignCol( uint32_t j, const VectorVariable<2U>& );
    void AssignCol( uint32_t j, const VectorVariable<3U>& );
    void AssignToDiagonal( const VectorVariable<1U>& );
    void AssignToDiagonal( const VectorVariable<2U>& );
    void AssignToDiagonal( const VectorVariable<3U>& );
    /// matrix vector multiplication
    DenseMatrix& operator*=( const VectorVariable<1U>& );
    DenseMatrix& operator*=( const VectorVariable<2U>& );
    DenseMatrix& operator*=( const VectorVariable<3U>& );

    void AssignRow( uint32_t i, const ArrayVariable& );
    void AssignCol( uint32_t j, const ArrayVariable& );
    void AssignToDiagonal( const ArrayVariable& );
    /// matrix array multiplication
    DenseMatrix& operator*=( const ArrayVariable& );

    void AssignRow( uint32_t i, const FlaggedArrayVariable& );
    void AssignCol( uint32_t j, const FlaggedArrayVariable& );
    void AssignToDiagonal( const FlaggedArrayVariable& );
    /// matrix FlaggedArray multiplication
    DenseMatrix& operator*=( const FlaggedArrayVariable& );

    void ExportTo( TensorVariable<1U>& ) const;
    void ExportTo( TensorVariable<2U>& ) const;
    void ExportTo( TensorVariable<3U>& ) const;
    DenseMatrix& operator=( const TensorVariable<1U>& );
    DenseMatrix& operator=( const TensorVariable<2U>& );
    DenseMatrix& operator=( const TensorVariable<3U>& );
    /// matrix tensor multiplication
    DenseMatrix& operator*=( const TensorVariable<1U>& );
    DenseMatrix& operator*=( const TensorVariable<2U>& );
    DenseMatrix& operator*=( const TensorVariable<3U>& );
#endif
  
  private:
   std::array<std::array<double,mn_max>, mn_max> data;
   uint32_t  rows, cols;
   
   bool  CheckRange( uint32_t m, uint32_t n, const char* originator ) const;
   bool  CheckSizes( const DenseMatrix& mat, const char* originator ) const;
};

// associated operators

/// comparitor of elements using compareFloats
template<uint32_t mn_max>
bool operator==( const DenseMatrix<mn_max>&, const DenseMatrix<mn_max>& );


template<uint32_t mn_max>
DenseMatrix<mn_max>  operator+( const DenseMatrix<mn_max>& a, 
                                const DenseMatrix<mn_max>& b );

template<uint32_t mn_max>
DenseMatrix<mn_max>  operator-( const DenseMatrix<mn_max>& a, 
                                const DenseMatrix<mn_max>& b );

template<uint32_t mn_max>
DenseMatrix<mn_max>  operator*( const DenseMatrix<mn_max>& a, 
                                const DenseMatrix<mn_max>& b );

template<uint32_t mn_max>
std::vector<double>  operator*( const DenseMatrix<mn_max>& mat,
                                  const std::vector<double>& );

/// v^T = (v^T M)^T
template<uint32_t mn_max>
DenseMatrix<mn_max>  operator*( const std::vector<double>& v,
                                const DenseMatrix<mn_max>& M );

} // csmp

#endif

