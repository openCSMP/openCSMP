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
#include <unordered_set>


/**
 
 @brief For the accumulation of element matrices in finite element
 calculations.
 
 @author A.J. Bromage
 @date 2017
 
 */
namespace csmp {

  struct SymmetricSparseMatrixEdge {
    size_t i, j;
    SymmetricSparseMatrixEdge( size_t i, size_t j ) : i(i), j(j) { }
  };
  
}
namespace std {
  template<>
  struct hash<csmp::SymmetricSparseMatrixEdge> {
    typedef csmp::SymmetricSparseMatrixEdge argument_type;
    typedef std::size_t result_type;
    
    inline size_t operator()( const csmp::SymmetricSparseMatrixEdge& x )
    {
      std::hash<size_t> h;
      return h(x.i) * 2000007949 + h(x.j);
    }
  };
}

namespace csmp {

class SymmetricSparseMatrix {
  private:
      struct Impl;
      std::unique_ptr<Impl> pimpl_;

  public:

    class Builder {
    public:
      Builder( size_t m_x_n );
      ~Builder();

      void AddElement( size_t i, size_t j );

    private:
      friend class SymmetricSparseMatrix;

      size_t m_x_n_;
      std::unordered_set< SymmetricSparseMatrixEdge > pattern_;
    };

    SymmetricSparseMatrix& operator=( const SymmetricSparseMatrix& sp ) = delete;
    SymmetricSparseMatrix( const SymmetricSparseMatrix& sp ) = delete;
    SymmetricSparseMatrix( SymmetricSparseMatrix&& sp ) = delete;

    explicit SymmetricSparseMatrix( const Builder& builder ); 
    ~SymmetricSparseMatrix();

    double64            operator()( size_t, size_t ) const;
    double64&           operator()( size_t, size_t );
    double64            At( size_t, size_t ) const;
    double64&           At( size_t, size_t );

    void Reset( Builder& builder );

    void      Zero();

    double64  InfinityNorm() const;

    bool      ZeroesInDiagonal() const;
    bool      DiagonallyPositive() const;

    /// output in the form of C-style arrays indexed 0...n-1
    void      OutCompressedRowFormat( int32* ia, int32* ja, double64* a, bool reallocate=true ) const;
    void      OutCompressedRowFormat( long*  ia, long*  ja, double64* a, bool reallocate=true ) const;
    /// output to Fortran arrays indexed 1...n
    void      OutCompressedRowFormat1_n( int32* ia, int32* ja, double64* a, bool reallocate=true ) const;

    void      OutCompressedRowFormatParallel( std::vector<int32>&, std::vector<int32>&,
                                         std::vector<double64>, int32, bool ) const;

    void      Out( long precis=5L ) const;
    void      Out( const char* file ) const;
 };

} // end csmp

#endif 



