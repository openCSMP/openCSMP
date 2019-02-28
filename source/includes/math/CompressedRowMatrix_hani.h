#ifndef COMPRESSED_ROW_MATRIX_HANI_H
#define COMPRESSED_ROW_MATRIX_HANI_H

#include "SparseMatrix.h"
#include "DenseMatrix.h"

namespace csmp {

class SparseMatrix;

struct Entry {
  Entry( size_t j, double64 v ) : ind( j ), val( v ) {}
  Entry() {}
  size_t ind;
  double64 val;
};



/**

@brief CSMP's implementation of a compressed-row storage container for sparse matrices.

@author S.K. Matthai
@author SS. Geiger
@author Coumou, D.
@author Stephen G. Roberts
@date 2006

This class implements the classic condensed row storage format for sparse matrices
that is common to many linear algebra packages.
This matrix is always square.

Currently used to convert csmp::SparseMatrix objects before they are passed to the SAMG solver.

*/

struct CompressedRowMatrix_hani {
  explicit CompressedRowMatrix_hani( csmp::SparseMatrix& );
  CompressedRowMatrix_hani();
  ~CompressedRowMatrix_hani();
  CompressedRowMatrix_hani( const CompressedRowMatrix_hani& );
  CompressedRowMatrix_hani( std::string& CRMfname );
  CompressedRowMatrix_hani& operator=( const CompressedRowMatrix_hani& );

  double64  operator()( uint32, uint32 ) const;

  void Initialize( const csmp::SparseMatrix& );
  void InitializePointBased( const SparseMatrix&, size_t nsys );

  size_t Rows() const { return (ia.size() - 1U); }
  size_t RowsJV() const { return JV_.size(); }
  size_t TotalExistingEntries() const { return ja.size(); }

  void Out() const;
  void Out( const std::string& outfile ) const;

  void Resize( size_t rows );
  void SetDiricheletRHS_CRM( std::vector<double64>& rhs, std::vector<Entry>& dirich );
  void mapToGlobal( std::vector<double64>& sol, std::vector<Entry>& dirich );
  void Add( std::vector<size_t>& rows, std::vector<size_t>& cols, DenseMatrix<DM_MIN>& vals, double64 factor );
  void SetIndicatorPBA( int32 nsys, std::vector<int32>& iu );
  void SetIndicatorUBA( int32 nsys, std::vector<int32>& iu );
  void SetPBA( int32 nsys, std::vector<double64>& u_, std::vector<double64>& f_, std::vector<double64>& x, std::vector<double64>& b );
  void AddToCRM( Entry& entry );

  std::vector<std::vector<Entry>> JV_;
  std::map<size_t, double64>      Dirich_;
  std::vector<int32>              mapDirich_;
  std::vector<int32>            ia, ///< ia(ilo) and the last row ends at position ia(ihi+1)-1 (see next).
    ja; ///< ja - pointer array pointing to the column indices. that is, for each matrix element a(j) with ia(ilo)<=j<=ia(ihi+1)-1, ja(j) contains the column index of that element. since, within each row, the diagonal element is stored first (see above), we always have ja(ia(i))=i.
  std::vector<double64>  a, trh_, tx_;  ///< array containing the rows of the matrix, one after the other, each row starting with its diagonal element. the first row starts at position
  int32 nnzi, dof_;
};

/// to print the vector of diagonal elements in the matrix
void print_crm( std::vector<std::pair<std::pair<uint32, uint32>, std::vector<bool> > >& );

} // end csmp

#endif
