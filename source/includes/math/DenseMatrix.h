#ifndef CSMP_DENSE_MATRIX_H
#define CSMP_DENSE_MATRIX_H

#define USED_TOGETHER_WITH_CSMP

#include <stdexcept>
#include <typeinfo>
#include "CSMP_definitions.h"
#ifdef USED_TOGETHER_WITH_CSMP
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
template<size_t mn_max>
class DenseMatrix {
  public:
    DenseMatrix();
    DenseMatrix( size_t m, size_t n );
    DenseMatrix( size_t m, size_t n, double64 val );
    DenseMatrix( const DenseMatrix& );
    DenseMatrix( DenseMatrix&& );
    ~DenseMatrix();
    size_t Rows() const;
    size_t Cols() const;
    void Resize( size_t m, size_t n );
    /// accessors M(i,j)
    double64&       operator()( size_t m, size_t n );
    const double64& operator()( size_t m, size_t n ) const;
    /// assignment
    DenseMatrix& operator=( const DenseMatrix& );
    DenseMatrix& operator=( DenseMatrix&& );
    DenseMatrix& operator=( double64 );
    DenseMatrix& operator+=( const DenseMatrix& );
    DenseMatrix& operator-=( const DenseMatrix& );
    DenseMatrix  operator+( const DenseMatrix& ) const;
    DenseMatrix  operator-( const DenseMatrix& ) const;
    /// matrix - matrix multiplication -> M(A.rows,B.cols) (NB: creates temporary M)
    DenseMatrix& operator*=( const DenseMatrix& );
    /// matrix vector multiplication
    DenseMatrix& operator*=( const std::vector<double64>& );
    /// rhs is a C-array
    DenseMatrix& operator*=( const double64* );
    DenseMatrix& operator*=( double64 );
    /// C = A B, matrix - matrix multiplication -> M(A.rows,B.cols) returns temporary matrix
    DenseMatrix  operator*( const DenseMatrix& ) const;


    void Identity();
    void AssignToDiagonal( size_t diag_elmts, const double64& sc );
    /// returns vec = Mat * unity vector
    void RowCondenseTo( std::vector<double64>& ) const;
    void Zero();
    void ZeroRow( size_t row );
    void ZeroCol( size_t col );
    void Fill( double64 val );
    void FillRow( size_t row, double64 val );
    void FillCol( size_t col, double64 val );
    /// RES = M^T into its argument
    void Transposed( DenseMatrix& ) const;
    /// RES = M^T M
    void TransposedProduct( DenseMatrix& ) const;
    /// RES = A B^T
    void MultiplyWithTransposedOf( const DenseMatrix&, DenseMatrix& ) const;
    /// RES = A^T B
    void MultiplyTransposedOfWith( const DenseMatrix&, DenseMatrix& ) const;
    double64   RowSum( size_t row ) const;
    double64   ColSum( size_t col ) const;
    /// unscaled L1 matrix norm
    double64   NormL1() const;
    /// L_infinity matrix norm
    double64   NormL_Infinity() const;

    void In();
    void Out( long digits=5L ) const;

#ifdef USED_TOGETHER_WITH_CSMP
    /// assignment of point coordinates to matrix rows or columns
    void AssignRow( size_t i, const Point<1U>& );
    void AssignRow( size_t i, const Point<2U>& );
    void AssignRow( size_t i, const Point<3U>& );
    void AssignCol( size_t j, const Point<1U>& );
    void AssignCol( size_t j, const Point<2U>& );
    void AssignCol( size_t j, const Point<3U>& );
    void AssignToDiagonal( const Point<1U>& );
    void AssignToDiagonal( const Point<2U>& );
    void AssignToDiagonal( const Point<3U>& );
    /// matrix point multiplication
    DenseMatrix& operator*=( const Point<1U>& );
    DenseMatrix& operator*=( const Point<2U>& );
    DenseMatrix& operator*=( const Point<3U>& );

    void AssignToDiagonal( size_t diag_elmts, const ScalarVariable& );
    DenseMatrix& operator*=( const ScalarVariable& );

    void AssignRow( size_t i, const VectorVariable<1U>& );
    void AssignRow( size_t i, const VectorVariable<2U>& );
    void AssignRow( size_t i, const VectorVariable<3U>& );
    void AssignCol( size_t j, const VectorVariable<1U>& );
    void AssignCol( size_t j, const VectorVariable<2U>& );
    void AssignCol( size_t j, const VectorVariable<3U>& );
    void AssignToDiagonal( const VectorVariable<1U>& );
    void AssignToDiagonal( const VectorVariable<2U>& );
    void AssignToDiagonal( const VectorVariable<3U>& );
    /// matrix vector multiplication
    DenseMatrix& operator*=( const VectorVariable<1U>& );
    DenseMatrix& operator*=( const VectorVariable<2U>& );
    DenseMatrix& operator*=( const VectorVariable<3U>& );

    void AssignRow( size_t i, const ArrayVariable& );
    void AssignCol( size_t j, const ArrayVariable& );
    void AssignToDiagonal( const ArrayVariable& );
    /// matrix array multiplication
    DenseMatrix& operator*=( const ArrayVariable& );

    void AssignRow( size_t i, const FlaggedArrayVariable& );
    void AssignCol( size_t j, const FlaggedArrayVariable& );
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
   std::array<std::array<double64,mn_max>, mn_max> data;
   size_t    rows, cols;
   
   bool  CheckRange( size_t m, size_t n, const char* originator ) const;
   bool  CheckSizes( const DenseMatrix& mat, const char* originator ) const;
};

// associated operators

template<size_t mn_max>
DenseMatrix<mn_max>  operator+( const DenseMatrix<mn_max>& a, 
                                const DenseMatrix<mn_max>& b );

template<size_t mn_max>
DenseMatrix<mn_max>  operator-( const DenseMatrix<mn_max>& a, 
                                const DenseMatrix<mn_max>& b );

template<size_t mn_max>
DenseMatrix<mn_max>  operator*( const DenseMatrix<mn_max>& a, 
                                const DenseMatrix<mn_max>& b );

template<size_t mn_max>
std::vector<double64>  operator*( const DenseMatrix<mn_max>& mat,
                                  const std::vector<double64>& );

/// v^T = (v^T M)^T
template<size_t mn_max>
DenseMatrix<mn_max>  operator*( const std::vector<double64>& v,
                                const DenseMatrix<mn_max>& M );



/**
 
@class DenseMatrix  DenseMatrix "applied_math/DenseMatrix.h"
@author S.K. Matthaei
@author S. Geiger
@author Stephen G. Roberts
@date 2001

@section motivation Motivation

Class for operations on dense small matrices. Efficiency decreases if
the matrices have few non-zero entries.  

*/


/// constructor (i,j)
template<size_t mn_max>
inline DenseMatrix<mn_max>::DenseMatrix( size_t m, size_t n )
 : rows(m), cols(n)
 {
 }



/// operator (i,j)
template<size_t mn_max>
inline double64& DenseMatrix<mn_max>::operator()( size_t m, size_t n )
 {
    CheckRange( m, n, "DenseMatrix<mn_max>::operator()");
    return data[m][n];
 }



/// operator (i,j) const
template<size_t mn_max>
inline const double64& DenseMatrix<mn_max>::operator()( size_t m, size_t n ) const
 {
    CheckRange( m, n, "DenseMatrix<mn_max>::operator()");
    return data[m][n];
 }


 
// Rows()
template<size_t mn_max>
inline size_t DenseMatrix<mn_max>::Rows() const { return rows; }


// Cols()
template<size_t mn_max>
inline size_t DenseMatrix<mn_max>::Cols() const { return cols; }


/**  
    Resizes DenseMatrix without allocation of new memory.
    If the capacitiy is exceeded an exception is thrown,
    but only if the code is compiled in debug mode.
*/
template<size_t mn_max>
inline void DenseMatrix<mn_max>::Resize( size_t m, size_t n )
 {
#ifndef NDEBUG 
    if ( m > mn_max ) {
         std::cerr <<"\nDenseMatrix<"<< mn_max;
         std::cerr <<">::Resize: Requested m-rows exceed matric capacity (";
         std::cerr << m <<" versus "<< rows <<")."<< std::endl;
         throw std::length_error("DenseMatrix<mn_max>::Resize");
      }
    if ( n > mn_max ) {
         std::cerr <<"\nDenseMatrix<"<<  mn_max;
         std::cerr <<">::Resize: Requested n-columns exceed matric capacity (";
         std::cerr << n <<" versus "<< cols <<")."<< std::endl;
         throw std::length_error("DenseMatrix<mn_max>::Resize");
      }
#endif
    rows = m;
    cols = n;
 }


#ifndef NDEBUG 
/// indices checking but only in the debug version
template<size_t mn_max>
inline bool DenseMatrix<mn_max>::CheckRange( size_t m, size_t n, 
                                             const char* originator ) const
 {
    if ( m >= rows ) {
         std::cerr <<"\n"<< originator <<" row index violation, index="<< m;
         std::cerr <<" versus, row-max=" << rows << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckRange");
         return false;
      }
    if ( n >= cols ) {
         std::cerr <<"\n"<< originator <<" column index violation, index="<< n;
         std::cerr <<" versus, column-max=" << cols << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckRange");
         return false;
      }
    return true;
 }
#else
template<size_t mn_max>
inline bool DenseMatrix<mn_max>::CheckRange( size_t, size_t, 
                                             const char* ) const
 {
    return true;
 }
#endif




#ifndef NDEBUG 
template<size_t mn_max>
/// checks (in DEBUG mode) whether the sizes of the matrices on either side of the expression match
inline bool DenseMatrix<mn_max>::CheckSizes( const DenseMatrix& mat, 
                                                const char* originator ) const
#else
template<size_t mn_max>
inline bool DenseMatrix<mn_max>::CheckSizes( const DenseMatrix&, 
                                                const char* ) const
#endif
 {
#ifndef NDEBUG 
    if ( rows != mat.rows ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; rows1="<< rows;
         std::cerr <<" versus, rows2=" << mat.rows << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckSizes");
         return false;
      }
    if ( cols != mat.cols ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; columns1="<< cols;
         std::cerr <<" versus, columns2=" << mat.cols << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckSizes");
         return false;
      }
#endif
    return true;
 }


} // csmp

#endif

