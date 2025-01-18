#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include "CSMP_definitions.h"

namespace csmp {

class SparseMatrix;

/** 
    CSR matrix class created to emulate CompressedRowMatrix class
    when using OpenMP
     
    @author Julian Mindel
*/
struct CSRMatrix {
    explicit CSRMatrix( const csmp::SparseMatrix& );
    CSRMatrix();
    ~CSRMatrix();
    CSRMatrix( const CSRMatrix& );
    CSRMatrix& operator=( const CSRMatrix& );

    double  operator()( uint32_t, uint32_t ) const;

    void Initialize( const csmp::SparseMatrix & );
    void InitializePointBased( const SparseMatrix&, size_t nsys );

    size_t Rows() const { return (ia.size()-1U); }
    size_t TotalExistingEntries() const { return ja.size(); }

    void Out() const;

    /// refer to SAMG's documentation
    std::vector<int32_t>  ia,
                          ja;
    std::vector<double>   a;
};

} // end csmp

#endif // CSRMATRIX_H
