#ifndef COMPRESSED_ROW_MATRIX_H
#define COMPRESSED_ROW_MATRIX_H

#include "SparseMatrix.h"
#include "Index.h"
#include "Parameter.h"

namespace csmp {

class SparseMatrix;

/**

@brief CSMP's implementation of a compressed-row storage container for sparse matrices.

@author S.K. Matthai
@author SS. Geiger
@author Coumou, D.
@author Stephen G. Roberts
@date 2006
@author Qi Shao
@date 2022(modified)

This class implements the classic condensed row storage format for sparse matrices 
that is common to many linear algebra packages. 
This matrix is always square. 

Currently used to convert csmp::SparseMatrix objects before they are passed to the SAMG solver.
Now also directly used in PDE_Integrator (2022, Qi Shao)

*/
class CompressedRowMatrix {
  public:
    CompressedRowMatrix() {}
    ~CompressedRowMatrix() = default;

    /// this constructor constructs a compressed row matrix from the supplied sparse matrix by calling the Initialise function
    /// the constructed compressed row matrix is in SAMG format which has the diagonal term at the beginning of each row and fortran indexes.
    explicit CompressedRowMatrix( csmp::SparseMatrix&);

    /// custom move constructor
    CompressedRowMatrix(std::vector<int32_t>&& ia, std::vector<int32_t>&& ja, std::vector<double>&& a);

    /// copy constructor
    CompressedRowMatrix( const CompressedRowMatrix& );

    /// move constructor
    CompressedRowMatrix( CompressedRowMatrix&& ) noexcept ;

    CompressedRowMatrix& operator=( const CompressedRowMatrix& );
    CompressedRowMatrix& operator=( CompressedRowMatrix&& ) noexcept ;

    /// standard accessor of matrix elements (asserts i,j in debug mode)
    double  operator()( size_t, size_t ) const;

    /// range-checked accessor for matrix elements
    double  At( uint32_t, uint32_t ) const;

    /// use to assign matrix elements during accumulation; @attention matrix sparsity pattern must be established before
    void   Assign( uint32_t i, uint32_t j, double val );
    void   MultiplyEntryWith( uint32_t i, uint32_t j, double val );
    void   Add( uint32_t i, uint32_t j, double val );
    void   ZeroRow( size_t row );
    void   Zero();
    void   AddRowByAnotherRow(uint32_t i, uint32_t j);
    void   AssignRowByAnotherRow(uint32_t i, uint32_t j);

    /// Initialize function called by constructor to initialize a compressed row matrix from a supplied sparse matrix
    /// it also converts the matrix to SAMG format which has the diagonal term at the beginning of each row and fortran indexes.
    void Initialize( const csmp::SparseMatrix& );
    void InitializePointBased( const SparseMatrix&, size_t nsys );

    size_t Rows() const;
    size_t Cols() const;
    /// total number of available storage sites for non-zero matrix elements
    size_t NonZeroEntries() const;

    /// this function erases all elements in ia, ja and a
    void Erase();

    /// check whether matrix has been converted to SAMG format
    bool IsFormattedForSAMG() const;

    void Out() const;
    void Out( const std::string& outfile ) const;

    /// this function converts the supplies sparse matrix to compressed row matrix, but does not convert it to SAMG format
    void ConvertFromSparseMatrix(const csmp::SparseMatrix& );

    /// this function converts a compressed row matrix to SAMG format
    void ConvertToSAMGFormat();

  protected:
    /// ia stores the starting point/index for each row; its size equals num_rows + 1
    /// ja stores the column indices of all non-zero entries
    /// a stores the values of all non-zero entries
    std::vector<int32_t>  ia, ///< ia(ilo) and the last row ends at position ia(ihi+1)-1 (see next).
                          ja; ///< ja - pointer array pointing to the column indices. that is, for each matrix element a(j) with ia(ilo)<=j<=ia(ihi+1)-1, ja(j) contains the column index of that element. since, within each row, the diagonal element is stored first (see above), we always have ja(ia(i))=i.
    std::vector<double>   a;  ///< array containing the rows of the matrix, one after the other, each row starting with its diagonal element. the first row starts at position
    static const bool     verbose_ = true;

    template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
    friend void generateSparsityPatternEliminatingEssentialConditions( CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const COMPUTATION_DOMAIN<dim>&);

    friend class EigenSolver;
    friend class SAMG_Solver;
};

/// to print the vector of diagonal elements in the matrix
void print( std::vector<std::pair<std::pair<uint32_t,uint32_t>, std::vector<bool> > >& );


/// this standalone function generates sparsity pattern for a compressed row matrix from supplied test_operands_, DOF_indexes_ and computational domain, taking into account elimination of essential conditions
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void generateSparsityPatternEliminatingEssentialConditions( CompressedRowMatrix&, const std::map<Parameter,size_t>&, std::vector<size_t>&, const COMPUTATION_DOMAIN<dim>&);



} // end csmp

#endif
