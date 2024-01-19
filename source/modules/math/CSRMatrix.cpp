#include "CSRMatrix.h"
#include "SparseMatrix.h"
#include "CSMP_global_enumerations.h"

using namespace std;

namespace csmp {

CSRMatrix::CSRMatrix()
{
}

CSRMatrix::CSRMatrix( const csmp::SparseMatrix& spmat )
{
    Initialize( spmat );
}


CSRMatrix::~CSRMatrix()
{
}

CSRMatrix::CSRMatrix( const CSRMatrix& crm )
{
    *this = crm;
}

CSRMatrix&  CSRMatrix::operator=( const CSRMatrix& crm )
{
    if ( &crm != this ) {
        ia = crm.ia;
        ja = crm.ja;
        a  = crm.a;
    }
    return *this;
}


double  CSRMatrix::operator()( uint32_t i, uint32_t j ) const
{
    //Since the off diagonal elements are stored after each diagonal element in the sparse matrix,
    //we must perform a small search until we match the requested column for given row i.

    for ( uint32_t  index=static_cast<uint32_t>(ia[i]); index <ia[i+1]; index++ )
        if ( ja[index-1]  ==  static_cast<int32_t>(j+1)  ) return a[ index - 1  ];
    return 0.;
}

void CSRMatrix::Initialize( const SparseMatrix& A )
 {
      ia.resize( (A.Rows() + 1U) ); vector<int32_t>( ia ).swap( ia );
      // ja is constructed with zero diagonal entries
      ja.resize( A.Entries(), 0 );  vector<int32_t>( ja ).swap( ja );
      // 'a' stores the non-zero entries of the sparse matrix, row after row
      a.resize( ja.size() );        vector<double>( a ).swap( a );

      // looping over all rows intializing ja and testing for diagonal entries which are zero
      // here n counts from 0 to j=nnu, i.e. all non-zero elements in the matrix
      uint32_t n(0U);
      ia[0] = 0;

#if defined(_OPENMP ) 
#pragma omp parallel for
      //first touch of ia.
      for ( int32_t i=0U; i < A.Rows(); i++ )
      {
          ia[i]       = 0;
      }
#endif

      // now create ia.
      for ( size_t i=0U; i < A.Rows(); i++ )
      {
          // looping over the non-zero elements row i
          for ( map<size_t,double>::const_iterator
                rit=A.RowBegin(i); rit!=A.RowEnd(i); rit++ ) {
              n++;
          }
          // setting matrix such that diagonal element is at the beginning of next row
          ia[i+1U]       = static_cast<int32_t>(n);
      }

#if defined(_OPENMP ) 
#pragma omp parallel for // algorithm has been done this way to complete idea of "first touch".
#endif
      // now first touch ja and a;
      for ( int32_t i=0U; i < A.Rows(); i++ )
      {
          // looping over the non-zero elements row i
          for ( size_t j = ia[i]; j < ia[i+1]; j++) {  //ia[i+1]-1 because indexes have not been shifted by 1 yet (happens below)
              ja[j] = 0;
              a[j]=0;
          }
      }

      n=0U;
      for ( size_t i=0U; i < A.Rows(); i++ )
       {
          int32_t  diag(UNSPECIFIED);
          // looping over the non-zero elements row i
          for ( map<size_t,double>::const_iterator
                rit=A.RowBegin(i); rit!=A.RowEnd(i); rit++ ) {
               // copying A's entry row(i) into the compressed row storage vector 'a'
               a[n]  = (*rit).second;
               // recording the corresponding column index in 'ja'
               // (NB: rit.first points to matrix column index from 0..rows-1)
               ja[n] = static_cast<int32_t>((*rit).first);
               // if i=j, i.e., if this is a diagonal elemnt, its position is recorded by 'diag'
               // if the diagonal element is zero, however, it will not have been stored in 'ja'
               // so that this situation is never encountered and diag remains UNSPECIFIED
               if ( ja[n] == static_cast<int32_t>(i) ) diag = static_cast<int32_t>(n);
               n++;
              }
          if ( diag == UNSPECIFIED ) {
               cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
               cout <<"\nSparseMatrix (rows=columns="<< A.Rows() <<") Zero entries (i=j): "<< endl;
               cout.setf(ios::scientific);
               long prec = cout.precision(15U);
               for ( size_t j=0U; j < A.Rows(); j++ )
                 if ( std::fabs(A(i,j)) < std::numeric_limits<double>::epsilon() )
                   cout <<"\n\t"<< j <<": "<< A(i,i);
               cout << endl;
               cout.unsetf( ios::scientific );
               cout.precision(prec);
               A.Out();
               throw underflow_error("CompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal.");
            }

          // setting matrix such that diagonal element is at the beginning of next row
          ia[i+1U]       = static_cast<int32_t>(n);
          // inserting the diagonal elements at the beginning of each row
          const uint32_t istart = static_cast<uint32_t>(ia[i]);
          int32_t jatemp   = ja[ istart ];
          double atemp = a[ istart ];
          int32_t dindex   = diag;
          ja[ istart ]   = ja[ dindex ];
          a[ istart]     = a[ dindex ];
          a[ dindex ]    = atemp;
          ja[ dindex ]   = jatemp;
       }

     // converting C++ array indices (0..n-1) into Fortran indices (1..n)
     for ( vector<int32_t>::iterator it=ia.begin(); it!=ia.end(); it++ ) (*it)++;
     for ( vector<int32_t>::iterator it=ja.begin(); it!=ja.end(); it++ ) (*it)++;

}  // end Initialize


/**
Initialises the public CSRMatrix vectors ia, ja, a for given
SparseMatrix in case the Point-based approach is selected.
*/
void CSRMatrix::InitializePointBased( const SparseMatrix& A, size_t nsys ) 
{
    // resize internal storage
    ia.resize( (A.Rows() + 1U) );
    ja.resize( A.Entries() );
    a.resize( ja.size() );
    size_t nnu_ = A.Rows();

    map<size_t,double>::const_iterator rit;
    long      i, j, k, row;
    int32_t     diag, istart, jatemp;
    double  atemp;
    bool      zero_diag_element(false);

    std::vector<int32_t>     temp; // auxilary vector
    temp.resize( nnu_ );
    for ( k = 0U; k < nnu_; k++)
        temp[k] = static_cast<int32_t>(k%(nnu_/nsys)*nsys+k/(nnu_/nsys));

    for ( i=j=0U, ia[0]=0; i < nnu_; i++ )
    {
        row = i%nsys*(nnu_/nsys)+i/nsys; // amending the order rows will be written in a[]
        for ( diag=-1, rit=A.RowBegin(row); rit!=A.RowEnd(row); rit++ ) {
            a[j]  = (*rit).second;
            ja[j] = temp[(*rit).first];
            // rit.first points to matrix entries indexed from 0..rows-1
            if ( ja[j] == temp[row] ) diag = static_cast<int32_t>(j);
            j++;
        }
        if ( diag == -1 ) zero_diag_element = true;

        ia[i+1]    = static_cast<int32_t>(j);

        // inserting the diagonal elements at the beginning of each row
        istart     = ia[i];
        jatemp     = ja[ static_cast<uint32_t>(istart) ];
        atemp      = a[ static_cast<uint32_t>(istart) ];
        ja[static_cast<uint32_t>(istart)] = ja[ static_cast<uint32_t>(diag) ];
        a[static_cast<uint32_t>(istart)]  = a[ static_cast<uint32_t>(diag) ];
        a[static_cast<uint32_t>(diag)]    = atemp;
        ja[static_cast<uint32_t>(diag)]   = jatemp;
    }

    if ( zero_diag_element ) {
        cout <<"\nCSRMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
        throw underflow_error("SparseM atrix::OutCompressedRowFormat");
    }

    // converting C array indices (0..n-1) into Fortran indices (1..n)
    for ( auto n=ia.begin(); n!=ia.end(); n++ )  (*n)++;
    for ( auto n=ja.begin(); n!=ja.end(); n++ )  (*n)++;

}  // end InitializePointBased


/** Outputs matrix to screen.
*/
void CSRMatrix::Out() const
{
    cout << flush <<"\nCSRMatrix::Out: "<< endl;
    cout <<"\nrow index vector 'ia' with size = "<<ia.size()<<"\n";
    for ( vector<int32_t>::const_iterator it=ia.begin(); it!=ia.end(); it++ )
        cout << *it <<" ";
    cout <<"\ncolumn index vector 'ja' with size = "<<ja.size()<<"\n";
    for ( vector<int32_t>::const_iterator it=ja.begin(); it!=ja.end(); it++ )
        cout << *it <<" ";
    cout <<"\nmatrix elements 'a' with size = "<<a.size()<<"\n";
    for ( size_t n=0U; n<ja.size(); n++ ) {
        cout << ja[n] <<":"<< a[n] <<" ";
        if ( n < ja.size()-1U and ja[n+1] < ja[n] ) cout << endl;
    }
    cout << endl;
    cout.flush();
}


} // end csmp

