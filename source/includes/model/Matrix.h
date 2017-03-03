#ifndef CSMP_MATRIX_H
#define CSMP_MATRIX_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class VectorVariable;
template<size_t> class TensorVariable;

/** 
    @brief Experimental implementation of class for a dense matrix
    greater than 50x50 matrix that can change its size.
    
    Supported operations (M = dense matrix, v = vector, _T = transposed,
    T = tensor matrix(dim,dim), default vector is a column vector (v = { c0,c2..cn-1 },
    vc = vector(dim) )
    
    @author S. K. Matthai
    @author S. Geiger
    @author S. G. Roberts
    @date 2001
    
    @todo SKM - implement non-member global methods for multiplication etc.
*/
class Matrix {
  public:
    Matrix();
    Matrix( size_t m, size_t n );
    Matrix( size_t m, size_t n, double64 val );
    Matrix( const Matrix& M );
    ~Matrix();
    size_t Rows() const;
    size_t Cols() const;
    void Resize( size_t m, size_t n );
    /// accessors M(i,j)
    double64&       operator()( size_t m, size_t n );
    const double64& operator()( size_t m, size_t n ) const;
    /// assignment
    Matrix& operator=( const Matrix& M );
    Matrix& operator=( double64 val );

    Matrix& operator=( const TensorVariable<2U>& M );
    Matrix& operator=( const TensorVariable<3U>& M );
    /// matrix tensor multiplication
    Matrix& operator*=( double64 val );
    Matrix& operator*=( const TensorVariable<2U>& M );
    Matrix& operator*=( const TensorVariable<3U>& M );
    /// matrix vector multiplication
    Matrix& operator*=( const VectorVariable<2U>& vc );
    Matrix& operator*=( const VectorVariable<3U>& vc );

    Matrix& operator+=( const Matrix& M );
    Matrix& operator-=( const Matrix& M );
    ///  matrix - matrix multiplication -> M(A.rows,B.cols) (operator creates temporary M)
    Matrix& operator*=( const Matrix& M );
    /// matrix vector multiplication
    Matrix& operator*=( const std::vector<double64>& v );
    /// matrix C-array vector multiplication
    Matrix& operator*=( const double64* v );
    /// adding of same-size matrices
    Matrix  operator+( const Matrix& M ) const;
    Matrix  operator-( const Matrix& M ) const;
    /// C = A B, matrix - matrix multiplication -> M(A.rows,B.cols)
    Matrix  operator*( const Matrix& M ) const;
    void Identity();
    
    void AssignToDiagonal( const VectorVariable<2U>& vc );
    void AssignToDiagonal( const VectorVariable<3U>& vc );
    void ExportTo( TensorVariable<2U>& ts ) const;
    void ExportTo( TensorVariable<3U>& ts ) const;
    /// returns vec = Mat * unity vector
    void RowCondenseTo( std::vector<double64>& vec ) const;

    void Zero();
    void ZeroRow( size_t row );
    void ZeroCol( size_t col );
    void Fill( double64 val );
    void FillRow( size_t row, double64 val );
    void FillCol( size_t col, double64 val );
    std::vector<double64> ReturnRow( size_t row ) const;
    std::vector<double64> ReturnCol( size_t col ) const;
    void AssignToRow( size_t row, std::vector<double64>& vec );
    void AssignToCol( size_t row, std::vector<double64>& vec );
    void Inversed( Matrix& RES ) const;
    /// RES = M^T into its argument
    void Transposed( Matrix& RES ) const;
    /// RES = M^T M
    void TransposedProduct( Matrix& RES ) const;
    /// RES = A B^T
    void MultiplyWithTransposedOf( const Matrix& B, Matrix& RES ) const;
    /// RES = A^T B
    void MultiplyTransposedOfWith( const Matrix& B, Matrix& RES ) const;
    double64   RowSum( size_t row ) const;
    double64   ColSum( size_t col ) const;
    /// unscaled matrix norms
    double64   NormL1() const;
    double64   NormL_Infinity() const;

    Matrix  Minor( const size_t row, const size_t col );

    void In();
    void Out(long digits=5L) const { Out(std::cout, digits); }
    void Out( std::ostream& os, long digits=5L ) const;

    /// LU decomposition and back-substitution
    void LUDecomposition(std::vector<size_t>& index, double64& d);
    void LUBackSubstitution(std::vector<size_t>& index, std::vector<double64>& b);
  
  private:
   size_t                               rows, cols;
   std::vector<std::vector<double64> >  data;
   
   bool  CheckRange( size_t m, size_t n, const char* originator ) const;
   bool  CheckSizes( const Matrix& mat, const char* originator ) const;
};


std::vector<double64>  operator*( const Matrix& mat, const std::vector<double64>& v );


/// v^T = (v^T M)^T
Matrix  operator*( const std::vector<double64>& v, const Matrix& M );

} // end csmp

#endif
