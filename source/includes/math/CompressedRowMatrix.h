#ifndef COMPRESSED_ROW_MATRIX_H
#define COMPRESSED_ROW_MATRIX_H

#include "SparseMatrix.h"
#include "Index.h"
#include "Parameter.h"
#include "ModelSubDomain.h"

namespace csmp {

class SparseMatrix;

/**
    Non-zero elements in the matrix as identified by node-to-node connectivity,
    but also taking into account the degrees of freedom (DOF) in the system, i.e.
    that certain rows and columns do not contain entries because they were assigned
    essential (Dirichlet) conditions in the right-hand vector.

    'ia' stores the starting point/index for each row; its size equals num_rows + 1
    'ja' stores the column indices of all non-zero entries
    'a' stores the values of all non-zero entries
*/
struct SparsityPattern {
    std::vector<int32_t> ia;
    std::vector<int32_t> ja;
    std::vector<double>  a;
};

/**

@brief CSMP's implementation of a compressed-row storage container for sparse matrices (jagged matrix).

@author S.K. Matthai
@author SS. Geiger
@author Coumou, D.
@author Stephen G. Roberts
@date 2006
@author Qi Shao
@date 2022(modified)

This class implements the  condensed row storage format for sparse matrices
that is common to many linear algebra packages, and using the i=0..m-1,  n=0..n-1 indexing that is standard to C++.
This matrix is always square, and the rows start with the diagonal element (like in Trilinos, SAMG, or PETSc

@attention to use this CRM with the SAMG Solver which relies on 1..n indexing, call the method  ConvertToSAMGFormat().

@attention the conversion of the matrix into a point-based format (grouping unknowns by nodes), must be called before
the matrix is converted to SAMG format.

*/
class CompressedRowMatrix {
  public:
    /// to generate CRM without having to initialise it immediately
    CompressedRowMatrix() {}
 
    /// Construct via move semantics directly from a generated pattern
    CompressedRowMatrix( SparsityPattern&& pattern ) noexcept
        : ia(std::move(pattern.ia)),
          ja(std::move(pattern.ja)), 
          a(std::move(pattern.a)) 
    {}
    
    /// this constructor constructs a compressed row matrix from the supplied sparse matrix by calling the Initialise function
    /// the constructed compressed row matrix is in SAMG format which has the diagonal term at the beginning of each row and fortran indexes.
    CompressedRowMatrix( const csmp::SparseMatrix&, bool use_SAMG_format );

    /// constructor for generating sparsity pattern for testing
    CompressedRowMatrix( std::vector<int32_t>& ia, std::vector<int32_t>& ja, std::vector<double>& a );

    /// Initialize function called by constructor to initialize a compressed row matrix from a supplied sparse matrix
    void Initialize( const csmp::SparseMatrix& );
    
    /// Initialse matrix using Fortran indexing and putting diag-element first into each row
    void InitializeSAMG( const csmp::SparseMatrix& );
    
    /// Initialise for SAMG point-based approach where variables are grouped by node
    void InitializePointBasedSAMG( const SparseMatrix&, size_t nsys );
 
    /// Diagnostic method to verify sparsity pattern integrity; @return true if sparsity pattern is valid, false otherwise
    bool VerifySparsityPattern() const;
    bool VerifySparsityPatternSAMG() const;

    /// this function erases all elements in ia, ja and a
    void Erase() noexcept;


    // METHODS REQUIRING SPARSITY PATTERN

    size_t Rows() const noexcept { if (ia.empty()) return 0; return ia.size()-1U; }
    size_t Cols() const noexcept { if (ia.empty()) return 0; return ia.size()-1U; } // square matrix!
    
    /// total of stored non-zero entries in the compressed row format
    size_t Entries() const noexcept { return ja.size(); }
 
    /// standard accessor of matrix elements (asserts i,j in debug mode), i=0..m-1,  n=0..n-1
    double  operator()( size_t, size_t ) const noexcept;

    /// range-checked accessor for matrix elements, i=0..m-1,  n=0..n-1
    double  At( size_t, size_t ) const noexcept;

    /// overwrite potential entries in the matrix
    void   Assign( size_t i, size_t j, double val );

    /// add to entries that have already been assigned so that their value is non-NaN
    void   Add( size_t i, size_t j, double val );
    
    /// MainAccumulation method adds only to i,j matrix entries that at contained in the pre-allocated sparsity pattern
    bool   AddIf( size_t i, size_t j, double val ) noexcept;

    /// multiply entries by value
    void   MultiplyEntryWith( size_t i, size_t j, double val );
    
    /** matrix-vector multiplication: result = this * vector
        @param vector input vector of size Cols()
        @param result output vector of size Rows(); will be resized if necessary
        @attention vector must have size equal to Cols(), result will be resized to Rows()
    */
    void MultiplyWith( const std::vector<double>& vector, std::vector<double>& result ) const;

    void   ZeroRow( size_t row );
    void   Zero() noexcept;
    void   AddRowByAnotherRow(size_t i, size_t j);
    void   AssignRowByAnotherRow(size_t i, size_t j);
    
    /// places the diagonal element at the beginning of each row which can be used for efficiency gains
    void   ReorderDiagonalFirst();

    void Out() const;
    void Out( long digits ) const;
    void Out( const std::string& outfile ) const;
    
    // FOR SAMG

    /// check whether matrix has been converted to SAMG format
    bool IsFormattedForSAMG() const noexcept;

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

    friend class EigenSolver;
    friend class SAMG_Solver;
    friend class CompressedRowMatrix_Test;
    friend class PDE_Integrator_Test;
};


/// Constructs a CompressedRowMatrix from a dense matrix stored in row-major format.
CompressedRowMatrix makeCompressedRowMatrix( const std::vector<double>& dense, size_t n, double tol = 0.0 );
                                                 
/// to print the vector of diagonal elements in the matrix (see CSMP_mathUtilities)
void print( std::vector<std::pair<std::pair<size_t,size_t>, std::vector<bool> > >& );


/**
    Generate sparsity pattern for CRM from supplied test_operands_, DOF_indexes_ and computational domain, eliminating essential conditions.
    The matrix remains in the i=0..m-1,  n=0..n-1 indexing format.
    
    @test OK
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE>
SparsityPattern generateSparsityPattern( const std::map<Parameter, size_t>& test_operands,
                                         const std::vector<size_t>& DOF_indexes,
                                         const ModelSubDomain<dim, CELLTYPE>& comp_domain,
                                         double non_zero_value = 0.0 );


// INLINE FUNCTIONS

/**
    Check whether the compressed row matrix has been converted into SAMG format.
*/
inline bool CompressedRowMatrix::IsFormattedForSAMG() const noexcept {
  if (ia.empty()) return false;
  if( *(ia.end() - 1) > static_cast<int32_t>(ja.size()) ) return true;
  return false;
}


inline double CompressedRowMatrix::operator()(size_t i, size_t j) const noexcept
{
    // use proper dimension checks.
    assert(i < Rows());
    assert(j < Cols());
    return this->At(i, j);
}


/**
 * Provides member access for CRMs with different indexing schemes.
 *
 * Matrix access is range checked, reporting out of bound and missing element errors
 *
 * @param i row index
 * @param j column index
 * @return double matrix element
 */
inline double CompressedRowMatrix::At(size_t i, size_t j) const noexcept
{
    assert(i < Rows() && j < Cols());

    const size_t row_start = static_cast<size_t>(ia[i]);
    const size_t row_end   = static_cast<size_t>(ia[i + 1]);

    if (row_start == row_end) return 0.0;

    if (IsFormattedForSAMG()) {
        // SAMG is 1-based: ia values are 1-based, ja entries are 1-based.
        // The first element is still usually the diagonal.
        if (static_cast<size_t>(ja[row_start - 1]) - 1 == j) {
            return a[row_start - 1];
        }

        for (size_t index = row_start; index < row_end; ++index) {
            if (static_cast<size_t>(ja[index - 1]) - 1 == j) 
                return a[index - 1];
        }
    } 
    else {
        // 1. CSMP TRUE DIAGONAL-FIRST OPTIMIZATION (O(1))
        if (i == j) {
            if (static_cast<size_t>(ja[row_start]) == i) {
                return a[row_start];
            }
        }

        // 2. SEARCH OFF-DIAGONALS
        const int32_t target_j = static_cast<int32_t>(j);
        if (j < i) {
            // Forward scan (skip diagonal)
            for (size_t index = row_start + 1; index < row_end; ++index) {
                if (ja[index] == target_j) return a[index];
            }
        } else {
            // Backward scan (stop before diagonal)
            for (size_t index = row_end - 1; index > row_start; --index) {
                if (ja[index] == target_j) return a[index];
            }
        }
    }

    return 0.0;
}



/**
   When calling Assign on a row for the first time, the very first element assigned to that row will determine the "front" of the row's memory.
   To ensure the diagonal-first requirement is met for the test, one must assign the diagonal elements first for each row.
 */
inline void CompressedRowMatrix::Assign(size_t i, size_t j, double val)
{
    assert( IsFormattedForSAMG() == false );
    assert( i < Rows() );
    assert( j < Cols() );

    const size_t row_start = static_cast<size_t>(ia[i]);
    const size_t row_end   = static_cast<size_t>(ia[i+1]);

    if (row_start == row_end) {
        throw runtime_error("CompressedRowMatrix::Assign: Error: Row is empty.");
    }

    // 1. TRUE DIAGONAL-FIRST OPTIMIZATION (O(1) execution)
    if (i == j) {
        if (static_cast<size_t>(ja[row_start]) == i) {
            a[row_start] = val;
            return;
        }
    }

    // 2. SEARCH OFF-DIAGONALS
    // Since the rest of the row is sorted ascending, the directional split is efficient.
    const int32_t target_j = static_cast<int32_t>(j);

    if (j < i) {
        // Forward scan (skip the diagonal at row_start)
        for (size_t index = row_start + 1; index < row_end; ++index) {
            if (ja[index] == target_j) {
                a[index] = val; 
                return;
            }
        }
    } else {
        // Backward scan
        for (size_t index = row_end - 1; index > row_start; --index) {
            if (ja[index] == target_j) {
                a[index] = val; 
                return;
            }
        }
    }

    // 3. ERROR HANDLING
    cerr << "\nCompressedRowMatrix::Assign: Error: Cannot find target element i=" 
         << i << ", j=" << j << endl;
    cerr << "candidate col IDs in the row are:" << endl;
    for (size_t index = row_start; index < row_end; ++index) {
        cerr << ja[index] << ", ";
    }
    cerr << endl;
    
    throw runtime_error("CompressedRowMatrix::Assign: Error: Cannot find target element in the compressed row matrix.");
}



inline void CompressedRowMatrix::Add(size_t i, size_t j, double val)
{
    // 1. Quick exit for zero values
    if (val == 0.0) return;

    assert( IsFormattedForSAMG() == false );
    assert( i < Rows() );
    assert( j < Cols() );

    const size_t row_start = static_cast<size_t>(ia[i]);
    const size_t row_end   = static_cast<size_t>(ia[i+1]);

    if (row_start == row_end) {
        throw runtime_error("CompressedRowMatrix::Add: Error: Row is empty.");
    }

    // 2. TRUE DIAGONAL-FIRST OPTIMIZATION (O(1) execution)
    if (i == j) {
        if (static_cast<size_t>(ja[row_start]) == i) {
            a[row_start] += val;
            return;
        }
    }

    // 3. SEARCH OFF-DIAGONALS
    // Since the rest of the row is sorted ascending, the directional split is efficient.
    const int32_t target_j = static_cast<int32_t>(j);

    if (j < i) {
        // Forward scan (skip the diagonal at row_start)
        for (size_t index = row_start + 1; index < row_end; ++index) {
            if (ja[index] == target_j) {
                a[index] += val; 
                return;
            }
        }
    } else {
        // Backward scan (stop before row_start)
        for (size_t index = row_end - 1; index > row_start; --index) {
            if (ja[index] == target_j) {
                a[index] += val; 
                return;
            }
        }
    }

    // 4. ERROR HANDLING
    cerr << "\nCompressedRowMatrix::Add: Error: Cannot find target element i=" 
         << i << ", j=" << j << endl;
    cerr << "candidate col IDs in the row are:" << endl;
    for (size_t index = row_start; index < row_end; ++index) {
        cerr << ja[index] << ", ";
    }
    cerr << endl;

    throw runtime_error("CompressedRowMatrix::Add: Error: Cannot find target element in the compressed row matrix.");
}



inline bool CompressedRowMatrix::AddIf( size_t i, size_t j, double val ) noexcept
 {
    assert( IsFormattedForSAMG() == false );
    assert( i < Rows() );
    assert( j < Cols() );

    // Basic bounds check to prevent row access violations
    if (i >= this->Rows()) return false;

    // Search for column j in the pre-allocated row i
    const auto row_start = static_cast<size_t>(ia[i]);
    const auto row_end   = static_cast<size_t>(ia[i+1]);
    
    for ( auto k = row_start; k < row_end; ++k ) {
        if (ja[k] == static_cast<int32_t>(j)) {
            a[k] += val;
            return true; // Successfully added
        }
    }
    return false; // Column not in sparsity pattern; value ignored
}



/// this function only sets all values in vector a to zeros, but keeps the indexes ia and ja unchanged.
inline void CompressedRowMatrix::Zero() noexcept
{
  std::fill(a.begin(), a.end(), 0.);
}


} // end csmp

#endif
