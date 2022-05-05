#ifndef CSMP_SPARSE_MATRIX_H
#define CSMP_SPARSE_MATRIX_H

#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cassert>

namespace csmp {

/**

@brief For the accumulation of element matrices in finite element
calculations.
 
@author S.K. Matthai
@date 2001

*/
class SparseMatrix {

  public:

    typedef std::vector<std::map<size_t,double> >::const_iterator    rowsConstIterator;
    typedef std::map<size_t,double>::const_iterator                  colsConstIterator;
    typedef std::vector<std::map<size_t,double> >::iterator          rowsIterator;
    typedef std::map<size_t,double>::iterator                        colsIterator;
    SparseMatrix(); 
    explicit SparseMatrix( size_t m_x_n );
    SparseMatrix( const SparseMatrix& sp );
    ~SparseMatrix();
    SparseMatrix& operator=( const SparseMatrix& sp );
    SparseMatrix& operator+=( SparseMatrix mat ); /// operator accumulates sparse matrices.  Initially created for OpenMP features.
    double            operator()( size_t, size_t ) const;
    double            At( size_t, size_t ) const;
    
    // this is not safe, nor efficient because we are dealing with a map! - use the row iterator RowBegin() and RowEnd() instead
//    const std::map<size_t,double>& Row( size_t i ) const { assert( i < data.size() ); return data[i]; }

    /// resets the rows=columns of the square matrix, retaining potential extra capacity of the vector used
    void                Resize( size_t n_x_m, bool preserve_allocated_memory=true );
  
    /// reports the non-zero elements currently stored in the sparse matrix
    size_t              Entries() const;
    size_t              Rows() const;
    size_t              Cols() const;
    rowsConstIterator   Begin() const;
    rowsConstIterator   End() const;
    colsConstIterator   RowBegin( size_t i) const;
    colsConstIterator   RowEnd( size_t i ) const;

    /// zeroing out rows in the context of parallel computations
    void      RemoveHalo( int32_t nrhalo );
    void      RemoveEntry( size_t, size_t );
    void      Erase();
    void      Zero();
    void      ZeroRow( size_t row );
    /// beware! - this is not an efficient operation
    void      ZeroColumn( size_t col );

    void      Add( size_t, size_t, double val );
    void      Assign( size_t, size_t, double val );
    void      MultiplyEntryWith( size_t i, size_t j, double val );
    void      MultiplyWith(const std::vector<double>& vec, std::vector<double>& res);

    /// returns j's of non-zero column entries in row
    void      ColumnIndices( size_t row, std::vector<uint32_t>& indices ) const;
    /// if suspected that a recent operation modified number of entries
    size_t    RecountEntries() const;

    double  InfinityNorm() const;

    bool      Symmetric() const;
    bool      ZeroesInDiagonal() const;
    bool      DiagonallyPositive() const;
    /// writes text matrix of zero's and one's to visualise sparsity pattern
    void      SparsityPattern( const char* txtfile ) const;
    
    template<class cspMat1,class cspMat2> 
    void      Assign( const cspMat1& idx, const cspMat2& data );
             
    /// output in the form of C-style arrays indexed 0...n-1
    void      OutCompressedRowFormat( int32_t*  ia, int32_t*  ja, double* a, bool reallocate=true ) const;
    void      OutCompressedRowFormat( long*  ia, long*  ja, double* a, bool reallocate=true ) const;
    /// output to Fortran arrays indexed 1...n
    void      OutCompressedRowFormat1_n( int32_t*  ia, int32_t*  ja, double* a, bool reallocate=true ) const;

    void      OutCompressedRowFormatParallel( std::vector<int32_t>&, std::vector<int32_t>&,
                                              std::vector<double>, int32_t, bool ) const;

    /// read (square) matrix from text file with zeros
    void      In( const char* file_name_without_extension );

    void      Out( long precis=5L ) const;
    void      Out( const char* file ) const;

	void OutForMatlab(const char * file) const;
    
  private:
    std::vector<std::map<size_t,double> >  data;
    size_t entries;

 };

} // end csmp

#endif 



