#ifndef COMPRESSED_ROW_MATRIX_H
#define COMPRESSED_ROW_MATRIX_H

#include "SparseMatrix.h"

namespace csmp {

class SparseMatrix;

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
class CompressedRowMatrix {
public:
    explicit CompressedRowMatrix( csmp::SparseMatrix& );
    CompressedRowMatrix();
    ~CompressedRowMatrix();
    CompressedRowMatrix( const CompressedRowMatrix& );
    CompressedRowMatrix& operator=( const CompressedRowMatrix& );

    double64  operator()( uint32, uint32 ) const;

    void Initialize( const csmp::SparseMatrix& );
    void InitializePointBased( const SparseMatrix&, size_t nsys );

    size_t Rows() const { return (ia.size()-1U); }
    size_t Cols() const { return (ja.size()-1U); }
    size_t TotalExistingEntries() const { return ja.size(); }

    void Out() const;
    void Out( const std::string& outfile ) const;

    std::vector<int32>     ia, ///< ia(ilo) and the last row ends at position ia(ihi+1)-1 (see next).
                           ja; ///< ja - pointer array pointing to the column indices. that is, for each matrix element a(j) with ia(ilo)<=j<=ia(ihi+1)-1, ja(j) contains the column index of that element. since, within each row, the diagonal element is stored first (see above), we always have ja(ia(i))=i.
    std::vector<double64>  a;  ///< array containing the rows of the matrix, one after the other, each row starting with its diagonal element. the first row starts at position
};

/// to print the vector of diagonal elements in the matrix
void print( std::vector<std::pair<std::pair<uint32,uint32>, std::vector<bool> > >& );

} // end csmp

#endif
