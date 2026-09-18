// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CompressedRowMatrix.h"
#include "SparseMatrix.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "NodeManifold.h"

using namespace std;

namespace csmp {

/**
    Generates sparsity pattern based on mesh connectivity information (i.e., node-tonode connectivity), degrees of freedom per variable, and potential Dirichlet boundary conditions.
    Manifold nodes and the connectivity that these create are also taken into account.
    
    From this information,  the vectors  'ia' and 'ja' in the compressed row storage are initialised
    The values stored in 'a' are initialised to zeros
    Elimination of essential conditions is taken into account.
    
    @note By using this element-centric branch expansion, your sparsity pattern automatically scales its complexity depending on the types of elements meeting at the SplitBoundary:
    This ensures that regardless of whether a triangle meets a quadrilateral or a prism meets a hexahedron along an internal geological contact,
    the linear algebra engine allocates the exact matrix profile needed to integrate the trans-interface weak form without throwing an out-of-bounds error.
 */
template<uint32_t dim, template<uint32_t> class CELLTYPE>
SparsityPattern generateSparsityPattern( const map<Parameter, size_t>& test_operands,
                                         const vector<size_t>& DOF_indexes,
                                         const ModelSubDomain<dim, CELLTYPE>& domain,
                                         double initial_value )
{
    //------------------------------------------------------------------
    // count number of active equations
    //------------------------------------------------------------------

    size_t n_eq = 0;

    for (size_t v : DOF_indexes)
        if (v != NULL_IDX)
            n_eq = max(n_eq, v + 1);

     //------------------------------------------------------------------
    // pass 1 & 2 combined: Use an array of vectors to build ja directly
    //------------------------------------------------------------------
    // This is faster than a 2-pass over all elements if the matrix is highly sparse.
    // A vector of vectors is used because we don't know row lengths yet.
    // TODO: can we do better than first assembling this?
    vector<vector<int32_t>> dynamic_ja(n_eq);
    vector<size_t> eq;
    eq.reserve(256);

    // 1. Accumulate connections
    for (const auto* cell : domain.CellVector())
    {
        eq.clear();

        for (const auto& [parameter, offset] : test_operands)
        {
            cell->ActiveEquationIndices( parameter.key, DOF_indexes, offset, eq );
        }

        // Add connections. eq is already unique
        for (size_t row : eq)
        {
            for (size_t col : eq)
            {
                dynamic_ja[row].push_back(static_cast<int32_t>(col));
            }
        }
    }

    SparsityPattern sp;
    sp.ia.resize(n_eq + 1);
    sp.ia[0] = 0;

    // 2. Sort, deduplicate per row, and build ia
    size_t total_nnz = 0;
    for (size_t i = 0; i < n_eq; ++i)
    {
        auto& row_cols = dynamic_ja[i];
        
        // Sort and remove duplicate columns from adjacent elements
        sort(row_cols.begin(), row_cols.end());
        auto unique_end = unique(row_cols.begin(), row_cols.end());
        row_cols.erase(unique_end, row_cols.end());

        // Now we know the exact length of this row
        total_nnz += row_cols.size();
        sp.ia[i + 1] = static_cast<int32_t>(total_nnz);
    }

    // 3. Allocate final contiguous arrays
    sp.ja.reserve(total_nnz);
    sp.a.resize(total_nnz, initial_value);

    // 4. Flatten dynamic_ja into sp.ja and move diagonal to front
    for (size_t row = 0; row < n_eq; ++row)
    {
        auto& row_cols = dynamic_ja[row];
        
        // Find diagonal
        auto diag = lower_bound(row_cols.begin(), row_cols.end(), static_cast<int32_t>(row));
        
        // Rotate diagonal to front if it exists
        if (diag != row_cols.end() && *diag == static_cast<int32_t>(row))
        {
            rotate(row_cols.begin(), diag, diag + 1);
        }

        // Append to contiguous ja array
        sp.ja.insert(sp.ja.end(), row_cols.begin(), row_cols.end());
    }
    
  return sp;

} // end generateSparsityPattern1

template SparsityPattern generateSparsityPattern<1U,Element>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<1U,Element>&, double );
template SparsityPattern generateSparsityPattern<2U,Element>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<2U,Element>&, double );
template SparsityPattern generateSparsityPattern<3U,Element>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<3U,Element>&, double );
  
template SparsityPattern generateSparsityPattern<1U,Face>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<1U,Face>&, double );
template SparsityPattern generateSparsityPattern<2U,Face>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<2U,Face>&, double );
template SparsityPattern generateSparsityPattern<3U,Face>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<3U,Face>&, double );
  
template SparsityPattern generateSparsityPattern<1U,InterFace>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<1U,InterFace>&, double );
template SparsityPattern generateSparsityPattern<2U,InterFace>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<2U,InterFace>&, double );
template SparsityPattern generateSparsityPattern<3U,InterFace>( const map<Parameter,size_t>&, const vector<size_t>&, const ModelSubDomain<3U,InterFace>&, double );








CompressedRowMatrix::CompressedRowMatrix( vector<int32_t>& input_ia,
                                          vector<int32_t>& input_ja,
                                          vector<double>& input_a )
  : ia(input_ia),
    ja(input_ja),
    a(input_a)
{
  if(Rows()>=numeric_limits<int32_t>::max()) throw out_of_range("CompressedRowMatrix(ctor): matrix size too large for SAMG solver");
  if(verbose_) cout<<"\nCompressedRowMatrix: called custom constructor"<<endl;
}




/** Construct a CRM form a Sparsematrix
*/
CompressedRowMatrix::CompressedRowMatrix( const SparseMatrix& spmat, bool use_SAMG_format )
 {
    if ( use_SAMG_format ) InitializeSAMG( spmat );
    else Initialize( spmat );
 }




/**
 * @brief Constructs a CompressedRowMatrix from a dense matrix stored
 *        in row-major format.
 *
 * Uses the constructor:
 *   CompressedRowMatrix( vector<int32_t>& ia,
 *                        vector<int32_t>& ja,
 *                        vector<double>&  a  )
 *
 * The sparsity pattern includes ALL entries whose absolute value
 * exceeds the supplied tolerance, plus the diagonal (always included
 * to guarantee a valid sparsity pattern even for zero diagonal entries).
 *
 * @param dense   Row-major dense matrix, size n*n
 * @param n           Matrix dimension
 * @param tol       Drop tolerance for off-diagonal entries (default 0.0)
 * @return      Initialised CompressedRowMatrix
 */
CompressedRowMatrix makeCompressedRowMatrix( const vector<double>& dense, size_t n, double tol )
{
    if ( dense.size() != n * n )
        throw invalid_argument(
            "MakeCompressedRowMatrix: dense.size() != n*n" );

    vector<int32_t> ia, ja;
    vector<double>  a;

    // ia has n+1 entries: ia[0]=0, ia[i] = start of row i
    ia.reserve( n + 1 );
    ia.push_back( 0 );

    for ( size_t row = 0; row < n; ++row )
    {
        // Always insert diagonal first (CRS convention used by this class)
        // Then insert off-diagonal non-zeros in column order
        //
        // Pass 1: always insert diagonal first
        ja.push_back( static_cast<int32_t>(row) );
        a .push_back( dense[row * n + row] );

        // Pass 2: collect off-diagonal non-zeros into a temporary buffer
        std::vector<std::pair<int32_t, double>> off_diag;
        off_diag.reserve(n - 1);

        for ( size_t col = 0; col < n; ++col )
        {
            if ( col == row ) continue; // skip diagonal; already inserted

            const double val = dense[row * n + col];
            if ( std::abs(val) > tol )
            {
                off_diag.emplace_back( static_cast<int32_t>(col), val );
            }
        }

        // Pass 3: sort by column index ascending to satisfy At()/Assign() scan assumptions
        std::sort( off_diag.begin(), off_diag.end(),
                   []( const auto& lhs, const auto& rhs )
                   {
                       return lhs.first < rhs.first;
                   });

        // Pass 4: append sorted off-diagonal entries
        for ( const auto& [col_idx, val] : off_diag )
        {
            ja.push_back( col_idx );
            a .push_back( val );
        }

        ia.push_back( static_cast<int32_t>(ja.size()) );
    }

    return CompressedRowMatrix( ia, ja, a );
}




void CompressedRowMatrix::MultiplyEntryWith( size_t i, size_t j, double val )
{
  //Check if the compressed row matrix is converted into SAMG format.
  //This operation only applies to the matrix in its original form.
  assert( IsFormattedForSAMG() == false );
  assert( i < Rows() );
  assert( j < Cols() );

  if ( i >= Rows() ) {
    cerr <<"\nCompressedRowMatrix::MultiplyEntryWith("<< i <<","<< j <<","<< val <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::MultiplyEntryWith");
  }
  if ( i >= Cols() ) {
    cout <<"\nCompressedRowMatrix::MultiplyEntryWith("<< i <<","<< j <<","<< val <<"): ";
    cout <<"Column access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::MultiplyEntryWith");
  }


  if(j <= i) {
    for ( auto index = ia[i]; index < ia[i + 1]; index++) {
      if (ja[static_cast<size_t>(index)] == static_cast<int32_t>(j) ) {a[static_cast<size_t>(index)] *= val; return;}
    }
  } else {
    for ( auto index = ia[i+1]-1; index >= ia[i]; index--) {
      if (ja[static_cast<size_t>(index)] == static_cast<int32_t>(j) ) {a[static_cast<size_t>(index)] *= val; return;}
    }
  }

  cerr <<"\nCompressedRowMatrix::MultiplyEntryWith: Error: Cannot find target element in the compressed row matrix."<<endl;
  throw runtime_error("CompressedRowMatrix::MultiplyEntryWith: Error: Cannot find target element in the compressed row matrix.");

} // end MultiplyEntryWith






/**
    Matrix-vector multiplication using compressed row storage format.
    
    Computes: result = this * vector
    
    For each row i, the operation is:
    result[i] = sum over j of (a[i,j] * vector[j])
    
    In compressed row format:
    - ia[i] points to the start of row i in arrays ja and a
    - ia[i+1] points to the start of row i+1
    - ja[k] contains the column index for element a[k]
    - a[k] contains the matrix element value
    
    @param x input vector of size Cols()
    @param result output vector; will be resized to Rows()
    
    @attention vector size must equal Cols(), otherwise exception is thrown
    @attention result vector is resized to Rows() and all entries are zeroed
    
    @author S.K. Matthai
    @date 2024
    
    @section performance Performance Notes
    - Time complexity: O(nnz) where nnz is number of non-zero entries
    - Space complexity: O(Rows()) for result vector
    - Cache-friendly: sequential access to matrix arrays
    - Suitable for large sparse systems
*/
void CompressedRowMatrix::MultiplyWith( const vector<double>& x,
                                        vector<double>& result ) const
{
    // Ensure we are NOT in SAMG mode
    assert( IsFormattedForSAMG() == false );

    if ( x.size() != this->Cols() )
    {
        ostringstream oss;
        oss << "CompressedRowMatrix::MultiplyWith: vector size (" << x.size()
            << ") does not match matrix columns (" << this->Cols() << ")";
        throw csmp::Exception( ERROR, "CompressedRowMatrix::MultiplyWith", oss.str().c_str() );
    }

    const size_t num_rows = this->Rows();
    if ( num_rows == 0 || this->Cols() == 0 )
    {
        result.clear();
        return;
    }

    result.assign( num_rows, 0.0 );

    // MATRIX-VECTOR MULTIPLICATION
    for ( size_t i = 0; i < num_rows; ++i )
    {
        double row_sum = 0.0;
        
        // In Standard C++ mode:
        // ia[i] is the physical start index in vectors 'a' and 'ja'
        // ia[i+1] is the physical end index (exclusive)
        const size_t row_start = static_cast<size_t>(ia[i]);
        const size_t row_end   = static_cast<size_t>(ia[i + 1]);

        for ( size_t k = row_start; k < row_end; ++k )
        {
            // ja[k] is already 0-based column index
            const size_t col = static_cast<size_t>(ja[k]);
            
            // Debug check for valid column indexing
            assert( col < this->Cols() );
            
            // Standard dot product: A[i,col] * x[col]
            row_sum += a[k] * x[col];
        }
        
        result[i] = row_sum;
    }
}







void CompressedRowMatrix::ZeroRow( size_t row )
{
  assert( IsFormattedForSAMG() == false );
  assert( row < Rows() );
  assert( row >= 0);

  if ( row >= Rows() ) {
    cerr <<"\nCompressedRowMatrix::ZeroRow("<< row <<"): ";
    cerr <<"Row access index out of range."<< endl;
    throw range_error("CompressedRowMatrix::ZeroRow");
  }

  //Check if the compressed row matrix is converted into SAMG format.
  //This operation should only apply to the matrix in its original form.
  if(IsFormattedForSAMG()) {
    cerr <<"\nCompressedRowMatrix::ZeroRow: Error: this operation should perform before the matrix is turned into SAMG format."<<endl;
    throw runtime_error("CompressedRowMatrix::ZeroRow:: Error: this operation should perform before the matrix is turned into SAMG format.");
  }

  for ( auto index=ia[row]; index <ia[row+1]; index++ ) a[ static_cast<size_t>(index) ] = 0.;
}



//this function erase all elements in ia, ja and a
void CompressedRowMatrix::Erase() noexcept
{
  ia.clear();
  ja.clear();
  a.clear();
}



/// this function add elements in one row by corresponding elements in another row, but excluding the diagonal elements
void CompressedRowMatrix::AddRowByAnotherRow( size_t i, size_t j )
{
    assert( IsFormattedForSAMG() == false );
    
    // 1. Corrected Bounds Checking
    if (i >= Rows() || j >= Rows()) {
        cerr << "\nCompressedRowMatrix::AddRowByAnotherRow: Row index out of range (i=" 
             << i << ", j=" << j << ")." << endl;
        throw range_error("CompressedRowMatrix::AddRowByAnotherRow");
    }

    // 2. Physical representation of logical column indices i and j
    const int32_t physical_i = static_cast<int32_t>(i);
    const int32_t physical_j = static_cast<int32_t>(j);

    // 3. Row Iteration (using format-aware pointers)
    const size_t row_j_start = static_cast<size_t>(ia[j]);
    const size_t row_j_end   = static_cast<size_t>(ia[j + 1]);

    for (size_t index = row_j_start; index < row_j_end; ++index) {
        const auto   col_id_physical = ja[index];
        const double value = a[index];

        // 4. Skip Diagonal Elements
        // We skip the element if it's the diagonal of row j (col == j)
        // OR if it's the diagonal of the target row i (col == i)
        if (col_id_physical != physical_i && col_id_physical != physical_j) {
            
            // 5. Use the format-aware Add() method
            // Convert physical col_id back to logical 0-based for the Add call
            auto col_logical = static_cast<size_t>(col_id_physical);
            this->Add(i, col_logical, value);
        }
    }
}



/// this function assign elements in one row by corresponding elements in another row, but excluding the diagonal elements
void CompressedRowMatrix::AssignRowByAnotherRow(size_t i, size_t j)
{
    assert( IsFormattedForSAMG() == false );

    // 1. Corrected Bounds Checking (Fixed the 'j' check)
    if (i >= Rows() || j >= Rows()) {
        cerr << "\nCompressedRowMatrix::AssignRowByAnotherRow: Row index out of range (i=" 
             << i << ", j=" << j << ")." << endl;
        throw range_error("CompressedRowMatrix::AssignRowByAnotherRow");
    }

    // 2. Format Awareness
    const bool isSAMG = IsFormattedForSAMG();
    const int32_t offset = isSAMG ? 1 : 0;

    // Physical representation of logical column indices i and j
    const int32_t physical_i = static_cast<int32_t>(i) + offset;
    const int32_t physical_j = static_cast<int32_t>(j) + offset;

    // 3. Row Iteration
    // We navigate the 'j' row. Pointers in ia are adjusted by the offset.
    const size_t row_j_start = static_cast<size_t>(ia[j] - offset);
    const size_t row_j_end   = static_cast<size_t>(ia[j + 1] - offset);

    for (size_t index = row_j_start; index < row_j_end; ++index) {
        const int32_t col_id_physical = ja[index];
        const double value = a[index];

        // 4. Skip Diagonal Elements
        // Do not assign to the diagonal of row i, and do not copy the diagonal of row j.
        if (col_id_physical != physical_i && col_id_physical != physical_j) {
            
            // 5. Use format-aware Assign()
            // Convert physical col_id back to logical 0-based for the call
            size_t col_logical = static_cast<size_t>(col_id_physical - offset);
            this->Assign(i, col_logical, value);
        }
    }
}



/**
 
Initialises the public CompressedRowMatrix vectors ia, ja, a 
from the given csmp::SparseMatrix.  

@param A The fully accumulated global solution matrix.

@section implementation Implementation

If nnu denotes the number of rows (variables), the non-zero entries 
of the i-th row (1 ≤ i ≤ nnu) are stored in a(j) where

ia(i) ≤ j ≤ ia(i+1)-1.

In particular, according to the above-mentioned convention about the location of the diagonal element,
a(ia(i)) contains the diagonal entry of row i. Note that ia(1) = 1 and  ia(nnu+1) = nna+1 where nna denotes
the total number of matrix entries stored.

ia, ja define the order of the off-diagonal elements in the solution
matrix storage vector a (that holds them row by row). 

The pointer vector ja has to be defined so that ja(j) (1 ≤ j ≤ nna) equals the original matrix’ column index of
a(j), ie, a(j) corresponds to the variable u(ja(j)). In particular, since a(ia(i)) contains the diagonal entry of row i,
we have ja(ia(i))=i.

While, in the compressed row format, all diagonal entries must be 
stored, off-diagonal entries are only stored if they are non-zero.  

@section application Application

Transfer of the global solution matrix to conventional solvers.  

@section messages Messages

Reports if the solution matrix contains zero diagonal entries.

@attention converts  SparseMatrix ( 0..n-1 ) into CSM in SAMG format (1..n)

*/
void CompressedRowMatrix::InitializeSAMG( const SparseMatrix& A )
 {
      ia.resize( (A.Rows() + 1U) ); ia.shrink_to_fit();
      // ja is constructed with zero diagonal entries
      ja.resize( A.Entries(), 0 );  ja.shrink_to_fit();
      // 'a' stores the non-zero entries of the sparse matrix, row after row
      a.resize( ja.size() );        a.shrink_to_fit();

      // looping over all rows intializing ja and testing for diagonal entries which are zero
      // here n counts from 0 to j=nnu, i.e. all non-zero elements in the matrix
      int32_t n(0U);
      ia[0] = 0;

      for ( size_t i{0}; i < A.Rows(); i++ )
       {
          int32_t  diag(UNSPECIFIED);
          // looping over the non-zero elements row i
          for ( auto rit=A.RowBegin(i); rit!=A.RowEnd(i); rit++ ) {
               // copying A's entry row(i) into the compressed row storage vector 'a'
               a[ static_cast<size_t>(n) ]  = (*rit).second;
               // recording the corresponding column index in 'ja'
               // (NB: rit.first points to matrix column index from 0..rows-1)
               ja[ static_cast<size_t>(n) ] = static_cast<int32_t>((*rit).first);
               // if i=j, i.e., if this is a diagonal elemnt, its position is recorded by 'diag'
               // if the diagonal element is zero, however, it will not have been stored in 'ja'
               // so that this situation is never encountered and diag remains UNSPECIFIED
               if ( ja[ static_cast<size_t>(n) ] == static_cast<int32_t>(i) ) diag = n;
               n++;
	          }
          if ( diag == UNSPECIFIED ) {
               cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
               cout <<"\nSparseMatrix (rows=columns="<< A.Rows() <<") Zero entries (i=j): "<< endl;
               cout.setf(ios::scientific);
               long prec = cout.precision(15U);
               for ( size_t i2{0}; i2 < A.Rows(); i2++ )
                 if ( fabs(A(i2,i2)) < numeric_limits<double>::epsilon() )
                   cout <<"\n\t"<< i2 <<": "<< A(i2,i2);
               cout << endl;
               cout.unsetf( ios::scientific );
               cout.precision(prec);
               A.Out();
               throw underflow_error("CompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal.");
            }
        
          // setting matrix such that diagonal element is at the beginning of next row 
          ia[i+1U] = n;
          // inserting the diagonal elements at the beginning of each row
          const int32_t istart = ia[i];
          int32_t jatemp = ja[ static_cast<size_t>(istart) ];
          double  atemp  = a[ static_cast<size_t>(istart) ];
          int32_t dindex = diag;
          ja[ static_cast<size_t>(istart) ]   = ja[ static_cast<size_t>(dindex) ];
          a[ static_cast<size_t>(istart) ]    = a[ static_cast<size_t>(dindex) ];
          a[ static_cast<size_t>(dindex) ]    = atemp;
          ja[ static_cast<size_t>(dindex) ]   = jatemp;
       }

     // converting C++ array indices (0..n-1) into Fortran indices (1..n) 
     //for ( vector<int32_t>::iterator it=ia.begin(); it!=ia.end(); it++ ) (*it)++;
     for (auto& it : ia) it++;
     //for ( vector<int32_t>::iterator it=ja.begin(); it!=ja.end(); it++ ) (*it)++;
     for (auto& it : ja) it++;
    
} // end InitialiseSAMG



/**
  Initialisation of zero-based CRM.
  
   @attention this function is similar to the previous one, but does not convert compressed row matrix to the SAMG format
*/
/**
  Initialisation of zero-based CRM.
  The diagonal element is moved to the front of each row block.
  
  @attention this function is similar to the previous one, but does not convert compressed row matrix to the SAMG format
*/
void CompressedRowMatrix::Initialize(const SparseMatrix& A)
{
    const size_t num_rows = A.Rows();
    const size_t num_entries = A.Entries();

    ia.resize(num_rows + 1U);
    ia.shrink_to_fit();
    ja.resize(num_entries, 0);
    ja.shrink_to_fit();
    a.resize(num_entries);
    a.shrink_to_fit();

    uint32_t n = 0U;
    ia[0] = 0;

    for (size_t i = 0; i < num_rows; i++) 
    {
        int32_t diag_pos = UNSPECIFIED;
        const uint32_t row_start_pos = n;

        // 1. Copy entries from SparseMatrix into CRM vectors
        for (auto rit = A.RowBegin(i); rit != A.RowEnd(i); ++rit) {
            a[n] = rit->second;
            ja[n] = static_cast<int32_t>(rit->first);

            // Check if this is the diagonal element
            if (ja[n] == static_cast<int32_t>(i)) {
                diag_pos = static_cast<int32_t>(n);
            }
            n++;
        }

        // 2. Verify diagonal existence
        if (diag_pos == UNSPECIFIED) {
            cerr << "\nCompressedRowMatrix::Initialize: Error: Zero/Missing diagonal at row " << i << endl;
            // Detailed error reporting (limiting output for large matrices)
            cerr.setf(ios::scientific);
            long prec = cerr.precision(15U);
            for (size_t i2 = 0; i2 < num_rows; i2++) {
                if (fabs(A(i2, i2)) < numeric_limits<double>::epsilon()) {
                    cerr << "\tRow " << i2 << " diagonal: " << A(i2, i2) << "\n";
                }
            }
            cerr.unsetf(ios::scientific);
            cerr.precision(prec);
            throw underflow_error("CompressedRowMatrix::Initialize: Zero value(s) in matrix diagonal.");
        }

        // 3. Move the diagonal element to the front of the row
        // Standard C++ 0-based swap
        if (static_cast<uint32_t>(diag_pos) != row_start_pos) {
            swap(a[row_start_pos], a[static_cast<size_t>(diag_pos)]);
            swap(ja[row_start_pos], ja[static_cast<size_t>(diag_pos)]);
        }

        ia[i + 1U] = static_cast<int32_t>(n);
    }
    
} // end Initialize



/**
 
Initialises the public CompressedRowMatrix vectors ia, ja, a for given 
SparseMatrix in case the SAMG point-based approach.

*/
void CompressedRowMatrix::InitializePointBasedSAMG( const SparseMatrix& A, size_t system_size )
 {
      // resize internal storage
      ia.resize( (A.Rows() + 1) );
      ja.resize( A.Entries() );
      a.resize( ja.size() );
   
      int32_t  i, j, k, diag;
      bool     zero_diag_element(false);
	  
      vector<int32_t>  temp( ja.size() ); // auxilary vector

      const int32_t nnu = static_cast<int32_t>(A.Rows());
      const int32_t nsys = static_cast<int32_t>(system_size);
      
	    for ( k = 0U; k < nnu; k++)
	      temp[ static_cast<size_t>(k) ] = k%(nnu/nsys)*nsys+k/(nnu/nsys);
	  		
      for ( i=j=0, ia[0]=0; i < nnu; i++ )
       {
	        int32_t row = i%nsys*(nnu/nsys)+i/nsys; // amending the order rows will be written in a[]
          auto rit = A.RowBegin( static_cast<size_t>(row) );
          for ( diag=-1; rit!=A.RowEnd( static_cast<size_t>(row) ); rit++ )
            {
               a[ static_cast<size_t>(j) ]  = (*rit).second;
               ja[ static_cast<size_t>(j) ] = temp[(*rit).first];
               // rit.first points to matrix entries indexed from 0..rows-1
               if ( ja[ static_cast<size_t>(j) ] == temp[ static_cast<size_t>(row) ] ) diag = j;
               j++;
            }
          if ( diag == -1 ) zero_diag_element = true;
          
          ia[ static_cast<size_t>(i+1) ] = j;
			
		      // inserting the diagonal elements at the beginning of each row
          const auto istart = ia[ static_cast<size_t>(i) ];
          swap( ja[static_cast<size_t>(istart)], ja[static_cast<size_t>(diag)] );
          swap(  a[static_cast<size_t>(istart)],  a[static_cast<size_t>(diag)] );
       }

      if ( zero_diag_element ) {
          cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
          throw underflow_error("SparseM atrix::OutCompressedRowFormat");
       }
	   
     // converting C array indices (0..n-1) into Fortran indices (1..n) 
     for ( vector<int32_t>::iterator l=ia.begin(); l!=ia.end(); ++l )  (*l)++;
     for ( vector<int32_t>::iterator l=ja.begin(); l!=ja.end(); ++l )  (*l)++;

} // end InitialisePointBased (SAMG)



/**
    Reorganised sparse matrix so that each row vector starts with the diagonal matrix element.
    The method is efficient because of:
    
    In-Place Swaps: swap, which for double and int32_t is extremely fast (just a few register moves). No temporary vectors or reallocations are needed.
    
    Early Exit: The if (ja[row_start] == i) check ensures that rows already correctly formatted (which should be most of them if your generator is working) are processed in $O(1)$ time.
    
    Cache Locality: Because it iiterates through ja and a linearly within each row,
    CPU cache hits are likely.
    
    No Full Sort: A full sort on every row would be $O(N \cdot K \log K)$, where $K$ is entries per row. This linear scan is $O(N \cdot K)$.

     @note modifies the matrix in-place. It assumes the matrix is using 0-based indexing.
*/
void CompressedRowMatrix::ReorderDiagonalFirst()
{
    // Safety check: SAMG mode usually implies 1-based indexing, 
    // which would break this logic.
    assert(IsFormattedForSAMG() == false);

    const size_t num_rows = this->Rows();
    if (num_rows == 0 || ia.empty()) return;

    for (size_t i = 0; i < num_rows; ++i)
    {
        const size_t row_start = static_cast<size_t>(ia[i]);
        const size_t row_end   = static_cast<size_t>(ia[i + 1]);

        // If the diagonal is already at the start, skip to next row
        if (row_start < row_end && static_cast<size_t>(ja[row_start]) == i)
        {
            continue;
        }

        // Search for the diagonal element in the current row
        bool found = false;
        for (size_t k = row_start + 1; k < row_end; ++k)
        {
            if (static_cast<size_t>(ja[k]) == i)
            {
                // Found it! Swap it into the first position of the row
                // We use swap for both indices and values
                swap(ja[k], ja[row_start]);
                swap(a[k], a[row_start]);
                found = true;
                break;
            }
        }

        // Some FE formulations (like certain
        // boundary conditions) might result in a missing diagonal.
        // This logs a warning if the matrix is expected to be complete.
        if (!found) {
            // Optional: Handle missing diagonal (e.g., log warning or throw)
            cerr << "WARNING: CompressedRowMatrix::ReorderDiagonalFirst: Missing diagonal at row " << i << endl;
        }
    }
    
} // end ReorderDiagonalFirst




/// this function converts compressed row matrix to the SAMG format
void CompressedRowMatrix::ConvertToSAMGFormat()
{
  if(IsFormattedForSAMG()) {
    cout<<"Already in SAMG format, nothing was done"<<endl;
    return;
  }

  for ( size_t i{0lu}; i < ia.size()-1; i++ ) {
    int32_t  diag(UNSPECIFIED);
    for( int32_t n=ia[i]; n<ia[i+1]; n++ ){
      if ( ja[ static_cast<size_t>(n) ] == static_cast<int32_t>(i) ) {diag = n; break;}
    }
    if ( diag == UNSPECIFIED ) {
      cout <<"\nCompressedRowMatrix::ConvertToSAMGFormat: Error: no diagonal element can be found in row: "<<i<<endl;
      Out();
      throw runtime_error("CompressedRowMatrix::ConvertToSAMGFormat: Error: no diagonal element can be found.");
    }

    // setting matrix such that diagonal element is at the beginning of next row
    // inserting the diagonal elements at the beginning of each row
    const auto istart = ia[i];
    int32_t jatemp = ja[ static_cast<size_t>(istart) ];
    double  atemp  = a[ static_cast<size_t>(istart) ];
    size_t dindex = static_cast<size_t>(diag);
    ja[ static_cast<size_t>(istart) ]   = ja[ dindex ];
    a[ static_cast<size_t>(istart) ]     = a[ dindex ];
    a[ dindex ]    = atemp;
    ja[ dindex ]   = jatemp;
  }

  // converting C++ array indices (0..n-1) into Fortran indices (1..n)
  for (auto& it : ia) (it)++;
  for (auto& it : ja) (it)++;

} // end ConvertToSAMGFormat




/**
 * Verifies the sparsity pattern of the CompressedRowMatrix against standard 
 * C++ / CSMP conventions (0-based indexing).
 * 
 * Standard Conventions (0-based):
 * ------------------------------
 * - ia: Row pointer vector of size (nnu + 1).
 * - ja: Column index vector of size (nna).
 * - a:  Values vector of size (nna).
 * 
 * Logic:
 * 1. ia[0] must be 0.
 * 2. ia[nnu] must be equal to nna (the total number of non-zeros).
 * 3. For any row i (0 <= i < nnu):
 *    - Row entries are located at index j where: ia[i] <= j < ia[i+1].
 *    - The diagonal element must be the FIRST entry of the row: ja[ia[i]] == i.
 *    - All column indices ja[j] must be within range [0, nnu-1].
 * 
 * @return true if sparsity pattern is valid, false otherwise.
 */
bool CompressedRowMatrix::VerifySparsityPattern() const
{
    // 1. Basic Dimension Checks
    if (ia.empty() || ja.empty() || a.empty())
    {
        cerr << "ERROR: VerifySparsityPattern: One or more matrix containers are empty." << endl;
        return false;
    }
    
    if (ja.size() != a.size())
    {
        cerr << "ERROR: VerifySparsityPattern: ja.size() != a.size()." << endl;
        return false;
    }

    const size_t nnu = ia.size() - 1; // Number of rows
    const size_t nna = ja.size();     // Total non-zeros
    
    cout << "\n=== CSMP Sparsity Pattern Verification (0-based) ===" << endl;
    cout << "Rows (nnu): " << nnu << " | Non-zeros (nna): " << nna << endl;

    // 2. Check ia Boundaries (0-based convention)
    if (ia[0] != 0)
    {
        cerr << "ERROR: ia[0] = " << ia[0] << " but must be 0 for standard indexing." << endl;
        return false;
    }

    if (static_cast<size_t>(ia[nnu]) != nna)
    {
        cerr << "ERROR: ia[nnu] = " << ia[nnu] << " but must be nna = " << nna << endl;
        return false;
    }

    // 3. Check Monotonicity, Diagonals, and Column Ranges
    for (size_t i = 0; i < nnu; ++i)
    {
        const size_t row_start = static_cast<size_t>(ia[i]);
        const size_t row_end   = static_cast<size_t>(ia[i+1]);

        // A: Check ia Monotonicity
        if (row_start > row_end)
        {
            cerr << "ERROR: ia is not monotonic at row " << i 
                      << ": ia[" << i << "]=" << row_start 
                      << " > ia[" << i+1 << "]=" << row_end << endl;
            return false;
        }

        // B: Check Diagonal Convention: ja[ia[i]] == i
        // Every row must have at least one entry (the diagonal)
        if (row_start == row_end)
        {
             cerr << "ERROR: Row " << i << " is empty. SAMG/Point-based requires a diagonal." << endl;
             return false;
        }

        if (ja[row_start] != static_cast<int32_t>(i))
        {
            cerr << "ERROR: Diagonal not first in row " << i 
                      << ": ja[" << row_start << "] = " << ja[row_start] 
                      << " but must be " << i << endl;
            return false;
        }

        // C: Check all Column Indices in row
        for (size_t k = row_start; k < row_end; ++k)
        {
            const int32_t col = ja[k];
            if (col < 0 || col >= static_cast<int32_t>(nnu))
            {
                cerr << "ERROR: ja[" << k << "] = " << col 
                          << " in row " << i 
                          << " is out of valid column range [0, " << nnu - 1 << "]" << endl;
                return false;
            }
        }
    }

    // 4. Print Preview (First 5 rows)
    const size_t rows_to_print = min(size_t(5), nnu);
    cout << "\nStructure Preview (First " << rows_to_print << " rows):" << endl;
    for (size_t i = 0; i < rows_to_print; ++i)
    {
        cout << "Row " << i << ": ";
        for ( auto k = ia[i]; k < ia[i+1]; ++k )
        {
            cout << "(" << ja[static_cast<size_t>(k)] << "," << a[static_cast<size_t>(k)] << ") ";
        }
        cout << endl;
    }

    cout << "\n=== Sparsity Pattern Verification PASSED ===" << endl;
    return true;
}



/**
    Verifies the sparsity pattern of the CompressedRowMatrix against SAMG Fortran conventions.
    
    SAMG conventions (1-based indexing):
  
    If nnu denotes the number of rows (variables), the non‐zero entries of the i‐th row (1 i nnu) are
    stored in a(j) where
    
    ia(i) j ia(i+1)‐1.
    
    In particular, according to the above‐mentioned convention about the location of the diagonal
    element, a(ia(i)) contains the diagonal entry of row i. Note that ia(1)=1 and ia(nnu+1)=nna+1 where
    nna denotes the total number of matrix entries stored.
    The pointer vector ja has to be defined so that ja(j) (1 j nna) equals the original matrix' column
    index of a(j), ie, a(j) corresponds to the variable u(ja(j)). In particular, since a(ia(i)) contains the
    diagonal entry of row i, we have ja(ia(i))=i.
    
    Summarizing, for any 1 i nnu and ia(i) j ia(i+1)‐1, we have a(j) = A(i,ja(j)) and the i‐th equa‐
    tion of Au=f reads:
    
    a(j) u(ja(j)) = f(i)
    j1 j j2
    
    where u(i) and f(i) denote the i‐th component of u and f, respectively, and j1=ia(i), j2=ia(i+1)‐1.
    
    @attention the arrays are vector objects defined as:
    
    vector<int32_t>  ia,   ia(ilo) and the last row ends at position ia(ihi+1)-1 (see next).
                     ja;   ja - pointer array pointing to the column indices. that is, for each matrix element a(j) with ia(ilo)<=j<=ia(ihi+1)-1, ja(j) contains the column index of that element. since, within each row, the diagonal element is stored first (see above), we always have ja(ia(i))=i.
    vector<double>   a;   array containing the rows of the matrix, one after the other, each row starting with its diagonal element. the first row starts at position
    
    In summary,
    - ia(i) has to be defined for all 1 i nnu+1!
    - The order of rows has to be such that the i‐th row (equation) corresponds to the i‐th variable.
    - The order of entries within each row is arbitrary except that the diagonal entry has to be first.
    - For symmetric matrices A, all (non‐zero) entries need to be stored (not just a triangular part).
    - ia has nnu+1 entries, where nnu is the number of rows
    - ia(1) = 1 (first row starts at index 1)
    - ia(nnu+1) = nna+1 (one past the last entry)
    - Row i entries: a(j) where ia(i) <= j <= ia(i+1)-1
    - Diagonal element is first in each row: ja(ia(i)) = i
    - ja(j) contains the column index of a(j), 1 <= ja(j) <= nnu
    
    @return true if sparsity pattern is valid, false otherwise
*/
bool CompressedRowMatrix::VerifySparsityPatternSAMG() const
{
    // ========================================================================
    // BASIC CHECKS (Using size_t for vector sizes)
    // ========================================================================
    if (ia.empty() || ja.empty() || a.empty())
    {
        cout << "ERROR: CompressedRowMatrix::VerifySparsityPattern: One or more matrix containers (ia, ja, a) are empty." << endl;
        return false;
    }
    
    if (ja.size() != a.size())
    {
        cout << "ERROR: CompressedRowMatrix::VerifySparsityPattern: ja.size() (" << ja.size() << ") != a.size() (" << a.size() << ")" << endl;
        return false;
    }

    const size_t nnu = ia.size() - 1; // Number of rows
    const size_t nna = ja.size();    // Total non-zeros (nna)
    
    cout << "\n=== CompressedRowMatrix Sparsity Pattern Verification ===" << endl;
    cout << "nnu (rows): " << nnu << " | nna (entries): " << nna << endl;

    // ========================================================================
    // CHECK ia BOUNDARIES (SAMG 1-based convention)
    // ========================================================================
    if (ia[0] != 1)
    {
        cout << "ERROR: ia[0] = " << ia[0] << " but must be 1 (SAMG 1-based start)." << endl;
        return false;
    }

    // ia[nnu] is the last element (index nnu in 0-based is the (nnu+1)-th element)
    if (static_cast<size_t>(ia[nnu]) != nna + 1)
    {
        cout << "ERROR: ia(nnu+1) = " << ia[nnu] << " but must be nna+1 = " << nna + 1 << endl;
        return false;
    }

    // ========================================================================
    // CHECK ia MONOTONICITY & ja COLUMN INDICES
    // ========================================================================
    for (size_t i = 0; i < nnu; ++i)
    {
        // 1. Check Monotonicity
        if (ia[i] > ia[i+1])
        {
            cout << "ERROR: ia is not monotonic at row " << i + 1 
                      << ": ia(" << i + 1 << ")=" << ia[i] 
                      << " > ia(" << i + 2 << ")=" << ia[i + 1] << endl;
            return false;
        }

        // 2. Check Diagonal Convention: ja(ia(i)) = i
        // row_start_1based is the value in ia, e.g., 1. 
        // In C++, this is index row_start_1based - 1.
        const int32_t row_start_1based = ia[i];
        const size_t diag_idx_0based = static_cast<size_t>(row_start_1based) - 1;

        if (ja[diag_idx_0based] != static_cast<int32_t>(i + 1))
        {
            cout << "ERROR: Diagonal not first in row " << i + 1 
                      << ": ja[" << diag_idx_0based << "] = " << ja[diag_idx_0based] 
                      << " but must be " << i + 1 << endl;
            return false;
        }

        // 3. Check all Column Indices in row
        const int32_t row_end_1based = ia[i+1];
        for (int32_t k = row_start_1based; k < row_end_1based; ++k)
        {
            size_t k_0based = static_cast<size_t>(k) - 1;
            int32_t col = ja[k_0based];
            if (col < 1 || col > static_cast<int32_t>(nnu))
            {
                cout << "ERROR: ja(" << k << ") = " << col 
                          << " in row " << i + 1 
                          << " out of range [1, " << nnu << "]" << endl;
                return false;
            }
        }
    }

    // ========================================================================
    // PRINT PREVIEW (First 5 rows)
    // ========================================================================
    const size_t rows_to_print = min(size_t(5), nnu);
    cout << "\nStructure Check (First " << rows_to_print << " rows):" << endl;
    for (size_t i = 0; i < rows_to_print; ++i)
    {
        cout << "Row " << i + 1 << ": ";
        for (int32_t k = ia[i]; k < ia[i+1]; ++k)
        {
            size_t idx = static_cast<size_t>(k) - 1;
            cout << "(" << ja[idx] << "," << a[idx] << ") ";
        }
        cout << endl;
    }

    cout << "\n=== Sparsity Pattern Verification PASSED ===" << endl;
    return true;

} // end VerifySparsityPatternSAMG




/** Outputs matrix to screen.
*/
void CompressedRowMatrix::Out() const
{
    const bool SAMG_format = IsFormattedForSAMG();
    // Offset is only applied if SAMG format is active (1-based)
    const int32_t offset = SAMG_format ? 1 : 0;

    cout << flush << "\nCompressedRowMatrix::Out: " << endl;
    if (SAMG_format) 
        cout << "Format: SAMG (1-based indices, diagonal first)" << endl;
    else 
        cout << "Format: Standard C++ (0-based indices)" << endl;

    // ... (Vector prints remain the same) ...

    // Use ia.size() - 1 to ensure we iterate through all rows
    const size_t num_rows = ia.size() - 1;

    for (size_t i = 0; i < num_rows; ++i) {
        // Calculate physical 0-based indices into the ja/a vectors
        // We cast first, then subtract to avoid signed underflow issues
        size_t row_start = static_cast<size_t>(ia[i]) - static_cast<size_t>(offset);
        size_t row_end   = static_cast<size_t>(ia[i + 1]) - static_cast<size_t>(offset);

        cout << "Row " << i << ": ";
        
        // This loop now correctly captures every element from [row_start, row_end)
        for (size_t k = row_start; k < row_end; ++k) {
            cout << ja[k] << ":" << a[k] << "  ";
        }
        cout << "\n";
    }
    cout << endl << flush;
}


void CompressedRowMatrix::Out(long digits) const
{
    const bool SAMG_format = IsFormattedForSAMG();
    const int32_t offset = SAMG_format ? 1 : 0;

    cout << flush << "\nCompressedRowMatrix::Out (Precision: " << digits << "):" << endl;
    if (SAMG_format) 
        cout << "Format: SAMG (1-based indices, diagonal first)" << endl;
    else 
        cout << "Format: Standard C++ (0-based indices)" << endl;

    cout << "\nrow index vector 'ia' (size=" << ia.size() << "):\n";
    for (auto it : ia) cout << it << " ";
    
    cout << "\ncolumn index vector 'ja' (size=" << ja.size() << "):\n";
    for (auto it : ja) cout << it << " ";
    
    cout << "\nmatrix elements 'a' (size=" << a.size() << "):\n";

    // --- Format and Precision Setup ---
    long prec = static_cast<long>(cout.precision());
    if (digits > 0) {
        cout.precision(digits);
        cout.setf(ios::scientific);
    }

    // --- Row Printing Loop ---
    for (size_t i = 0; i < Rows(); ++i) {
        // Calculate physical C++ vector indices based on internal format
        size_t row_start = static_cast<size_t>(ia[i] - offset);
        size_t row_end   = static_cast<size_t>(ia[i + 1] - offset);

        cout << "Row " << i << ": ";
        for (size_t k = row_start; k < row_end; ++k) {
            // Print the column index as stored, then the value
            cout << ja[k] << ":" << a[k] << "  ";
        }
        cout << "\n";
    }

    // --- Restore Console State ---
    if (digits > 0) {
        cout.unsetf(ios::scientific);
        cout.precision(prec);
    }
    cout << endl << flush;
}


/** Outputs matrix to text file.
*/
void CompressedRowMatrix::Out(const string& outfile) const
{
    ofstream ofs(outfile);
    if (!ofs.is_open()) {
        cerr << "CompressedRowMatrix::Out: Error: Could not open file " << outfile << endl;
        return; 
    }

    const bool SAMG_format = IsFormattedForSAMG();
    const int32_t offset = SAMG_format ? 1 : 0;

    ofs << "CompressedRowMatrix::Out: " << outfile << endl;
    if (SAMG_format) 
        ofs << "Format: SAMG (1-based indices, diagonal first)" << endl;
    else 
        ofs << "Format: Standard C++ (0-based indices)" << endl;

    ofs << "\nrow index vector 'ia' (size=" << ia.size() << "):\n";
    for (auto it : ia) ofs << it << " ";
    
    ofs << "\ncolumn index vector 'ja' (size=" << ja.size() << "):\n";
    for (auto it : ja) ofs << it << " ";
    
    ofs << "\nmatrix elements 'a' (size=" << a.size() << "):\n";

    // Set high precision for file storage (standard for earth science data)
    const auto original_precision = ofs.precision();
    ofs.precision(15);
    ofs << scientific;

    for (size_t i = 0; i < Rows(); ++i) {
        // Calculate physical indices normalized to 0-based vector access
        size_t row_start = static_cast<size_t>(ia[i] - offset);
        size_t row_end   = static_cast<size_t>(ia[i + 1] - offset);

        for (size_t k = row_start; k < row_end; ++k) {
            // Write column index (as stored) and value
            ofs << ja[k] << ":" << a[k] << " ";
        }
        ofs << "\n";
    }

    // Restore stream state
    ofs.unsetf(ios::scientific);
    ofs.precision(original_precision);
    ofs << endl;
}

} // end csmp

