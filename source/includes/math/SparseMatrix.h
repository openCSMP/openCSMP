#ifndef CSMP_SPARSE_MATRIX_H
#define CSMP_SPARSE_MATRIX_H

#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include "CSMP_number_types.h"
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

    typedef std::vector<std::map<size_t,double64> >::const_iterator    rowsConstIterator;
    typedef std::map<size_t,double64>::const_iterator                  colsConstIterator;
    typedef std::vector<std::map<size_t,double64> >::iterator          rowsIterator;
    typedef std::map<size_t,double64>::iterator                        colsIterator;
    SparseMatrix(); 
    explicit SparseMatrix( size_t m_x_n );
    SparseMatrix( const SparseMatrix& sp );
    ~SparseMatrix();
    SparseMatrix& operator=( const SparseMatrix& sp );
    SparseMatrix& operator+=( SparseMatrix mat ); /// operator accumulates sparse matrices.  Initially created for OpenMP features.
    double64            operator()( size_t, size_t ) const;
    double64            At( size_t, size_t ) const;
    const std::map<size_t,double64>& Row( size_t i ) const { assert( i < data.size() ); return data[i]; }

    void                Resize( size_t n_x_m );
  
    /// reports the non-zero elements currently stored in the sparse matrix
    size_t              Entries() const;
    size_t              Rows() const;
    size_t              Cols() const;
    rowsConstIterator   Begin() const;
    rowsConstIterator   End() const;
    colsConstIterator   RowBegin( size_t i) const;
    colsConstIterator   RowEnd( size_t i ) const;

    /// zeroing out rows in the context of parallel computations
    void      RemoveHalo( int32 nrhalo );
    void      RemoveEntry( size_t, size_t );
    void      Erase();
    void      Zero();
    void      ZeroRow( size_t row );
    /// beware! - this is not an efficient operation
    void      ZeroColumn( size_t col );

    void      Add( size_t, size_t, double64 val );
    void      Assign( size_t, size_t, double64 val );
    void      MultiplyEntryWith( size_t i, size_t j, double64 val );
    void      MultiplyWith(const std::vector<double64>& vec, std::vector<double64>& res);

    /// returns j's of non-zero column entries in row
    void      ColumnIndices( size_t row, std::vector<size_t>& indices ) const;
    /// if suspected that a recent operation modified number of entries
    size_t    RecountEntries() const;

    double64  InfinityNorm() const;

    bool      Symmetric() const;
    bool      ZeroesInDiagonal() const;
    bool      DiagonallyPositive() const;
    /// writes text matrix of zero's and one's to visualise sparsity pattern
    void      SparsityPattern( const char* txtfile ) const;
    
    template<class cspMat1,class cspMat2> 
    void      Assign( const cspMat1& idx, const cspMat2& data );
             
    /// output in the form of C-style arrays indexed 0...n-1
    void      OutCompressedRowFormat( int32* ia, int32* ja, double64* a, bool reallocate=true ) const;
    void      OutCompressedRowFormat( long*  ia, long*  ja, double64* a, bool reallocate=true ) const;
    /// output to Fortran arrays indexed 1...n
    void      OutCompressedRowFormat1_n( int32* ia, int32* ja, double64* a, bool reallocate=true ) const;

    void      OutCompressedRowFormatParallel( std::vector<int32>&, std::vector<int32>&,
                                         std::vector<double64>, int32, bool ) const;

    /// read (square) matrix from text file with zeros
    void      In( const char* file_name_without_extension );

    void      Out( long precis=5L ) const;
    void      Out( const char* file ) const;

	void OutForMatlab(const char * file) const;
    
  private:
    std::vector<std::map<size_t,double64> >  data;
    size_t entries;

 };

} // end csmp

#endif 



