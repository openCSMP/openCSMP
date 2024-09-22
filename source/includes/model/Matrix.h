#ifndef CSMP_MATRIX_H
#define CSMP_MATRIX_H

#include "CSMP_definitions.h"
#include "DynamicArray2D.h"

namespace csmp {

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
    Matrix() {}
    Matrix( size_t m, size_t n );
    Matrix( size_t m, size_t n, double val );

    size_t Rows() const;
    size_t Cols() const;
    void Resize( size_t m, size_t n );
    
    /// accessors M(i,j)
    double&       operator()( size_t m, size_t n );
    const double& operator()( size_t m, size_t n ) const;
    /// assignment
    Matrix& operator=( double val );

    Matrix& operator+=( const Matrix& M );
    Matrix& operator-=( const Matrix& M );
    ///  matrix - matrix multiplication -> M(A.rows,B.cols) (operator creates temporary M)
    Matrix& operator*=( const Matrix& M );
    Matrix& operator*=( double val );
    /// matrix vector multiplication
    Matrix& operator*=( const std::vector<double>& v );
    /// adding of same-size matrices
    Matrix  operator+( const Matrix& M ) const;
    Matrix  operator-( const Matrix& M ) const;
    /// C = A B, matrix - matrix multiplication -> M(A.rows,B.cols)
    Matrix  operator*( const Matrix& M ) const;

    void Identity();
    
    /// returns vec = Mat * unity vector
    void RowCondenseTo( std::vector<double>& vec ) const;

    void Zero();
    void ZeroRow( size_t row );
    void ZeroCol( size_t col );
    void Fill( double val );
    void FillRow( size_t row, double val );
    void FillCol( size_t col, double val );
    std::vector<double> ReturnRow( size_t row ) const;
    std::vector<double> ReturnCol( size_t col ) const;
    void AssignToRow( size_t row, std::vector<double>& vec );
    void AssignToCol( size_t row, std::vector<double>& vec );

    /// RES = M^T into its argument
    void Transposed( Matrix& RES ) const;

    /// RES = M^T M
    void TransposedProduct( Matrix& RES ) const;

    /// RES = A B^T
    void MultiplyWithTransposedOf( const Matrix& B, Matrix& RES ) const;

    /// RES = A^T B
    void MultiplyTransposedOfWith( const Matrix& B, Matrix& RES ) const;
    double   RowSum( size_t row ) const;
    double   ColSum( size_t col ) const;
    /// unscaled matrix norms
    double   NormL1() const;
    double   NormL_Infinity() const;

    Matrix  Minor( const size_t row, const size_t col );

    void In();
    void Out( long digits=5L ) const;

    /// LU decomposition and back-substitution
    void LUDecomposition(std::vector<uint32_t>& index, double& d);
    void LUBackSubstitution(std::vector<uint32_t>& index, std::vector<double>& b);
  
  private:
   DynamicArray2D<double>  data;
   
   bool  CheckRange( size_t m, size_t n, const char* originator ) const;
   bool  CheckSizes( const Matrix& mat, const char* originator ) const;
};


std::vector<double>  operator*( const Matrix& mat, const std::vector<double>& v );


/// v^T = (v^T M)^T
Matrix  operator*( const std::vector<double>& v, const Matrix& M );

} // end csmp

#endif
