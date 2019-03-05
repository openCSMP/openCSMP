#ifndef COMPRESSED_SPARSE_ROW_MATRIX_H
#define COMPRESSED_SPARSE_ROW_MATRIX_H

#include "CompressedRowMatrix.h"
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
class CompressedSparseRowMatrix : public CompressedRowMatrix
{
public:
    CompressedSparseRowMatrix();

    void      Resize( size_t n_x_m ) {};

    void      Add( size_t, size_t, double64 val ) {};    
    void      ZeroRow( size_t row ) {};
    
    void      Erase() {};
    void      Zero() {};
    
    double64  At( size_t row, size_t col) const { return (*this)(row, col);  };
};

} // end csmp

#endif
