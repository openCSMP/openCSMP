#include <cmath>
#include <limits>
#include <cassert>
#include "CompressedRowMatrix.h"
#include "SparseMatrix.h"
#include "Exception.h"


using namespace std;

namespace csmp {

/*! \file CSMP_mathUtilities.cpp */

/**
@addtogroup CSMPglobalFunctions
@{
*/

void print(  vector<pair<pair<uint32,uint32>,vector<bool> > >&  v )
 {
       cout <<"\nvector of off-diagonal elements:\n";
       int32 n(0);
       
         for ( vector<pair<pair<uint32,uint32>,vector<bool> > >::iterator it=v.begin(); it!=v.end(); it++ )
            {
                 cout <<"\nBlock "<< n++ <<" range: "<< (*it).first.first <<" - "<< (*it).first.second << endl;
                 cout <<"boolean vector (size="<< (*it).second.size() <<"):\n";
                 for ( vector<bool>::const_iterator i=(*it).second.begin(); i!=(*it).second.end(); i++ )
                   if ( *i ) cout <<" true  ";
                   else  cout <<"false ";
                   
                cout << endl << endl;
            }
  }
/**
@}
*/

CompressedRowMatrix::CompressedRowMatrix()
 {
 }



CompressedRowMatrix::CompressedRowMatrix( csmp::SparseMatrix& spmat )
 {
    Initialize( spmat );
 }


CompressedRowMatrix::~CompressedRowMatrix()
 {
 }

CompressedRowMatrix::CompressedRowMatrix( const CompressedRowMatrix& crm )
 {
      *this = crm;
 }

CompressedRowMatrix&  CompressedRowMatrix::operator=( const CompressedRowMatrix& crm )
  {
       if ( &crm != this ) {
             ia = crm.ia;
             ja = crm.ja;
             a  = crm.a;
         }
       return *this;
  }



/*
 Julian Mindel:  I proceeded to comment out the old code which contained the version of the () operator used before
 I have left it here below in the comment section for legacy purposes.

double64  CompressedRowMatrix::operator()( uint32 i, uint32 j ) const
{
  // if the diagonal element is requested
  if ( i == j ) return  a[ static_cast<uint32>(ia[i]) ];

  // now all row elements are stored to the right of the diagonal (by convention)
  if (  a.size() - static_cast<uint32>(ia[i]) == ia.size() or  ia[i+1] - ia[i] == static_cast<int32>(ia.size()) )
    return  a[ static_cast<uint32>(ia[i]) + j - 1U ];

  // i is the diagonal element of the matrix
  for ( uint32  index=static_cast<uint32>(ia[i]); index <= ia.size(); index++ )
    if ( ja[index]  ==  static_cast<int32>(j)  ) return a[ index  ];

  return 0.;
}
*/

double64  CompressedRowMatrix::operator()( uint32 i, uint32 j ) const
{
	assert( i < ia.size()-1U );
	assert( j < ia.size()-1U );

	//Since the off diagonal elements are stored after each diagonal element in the sparse matrix,
	//we must perform a small search until we match the requested column for given row i.
	for ( uint32  index=static_cast<uint32>(ia[i]); index <ia[i+1]; index++ )
    if ( ja[index-1]  ==  static_cast<int32>(j+1)  ) return a[ index - 1  ];

	return 0.;
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
*/

void CompressedRowMatrix::Initialize( const SparseMatrix& A ) 
 {
      ia.resize( (A.Rows() + 1U) ); vector<int32>( ia ).swap( ia );
      // ja is constructed with zero diagonal entries
      ja.resize( A.Entries(), 0 );  vector<int32>( ja ).swap( ja );
      // 'a' stores the non-zero entries of the sparse matrix, row after row
      a.resize( ja.size() );        vector<double64>( a ).swap( a );

      // looping over all rows intializing ja and testing for diagonal entries which are zero
      // here n counts from 0 to j=nnu, i.e. all non-zero elements in the matrix
      uint32 n(0U);
      ia[0] = 0;

      for ( size_t i=0U; i < A.Rows(); i++ )
       {
          int32  diag(UNSPECIFIED);
          // looping over the non-zero elements row i
          for ( map<size_t,double64>::const_iterator
                rit=A.RowBegin(i); rit!=A.RowEnd(i); rit++ ) {
               // copying A's entry row(i) into the compressed row storage vector 'a'
               a[n]  = (*rit).second;
               // recording the corresponding column index in 'ja'
               // (NB: rit.first points to matrix column index from 0..rows-1)
               ja[n] = static_cast<int32>((*rit).first);
               // if i=j, i.e., if this is a diagonal elemnt, its position is recorded by 'diag'
               // if the diagonal element is zero, however, it will not have been stored in 'ja'
               // so that this situation is never encountered and diag remains UNSPECIFIED
               if ( ja[n] == static_cast<int32>(i) ) diag = static_cast<int32>(n);
               n++;
	          }
          if ( diag == UNSPECIFIED ) {
               cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
               cout <<"\nSparseMatrix (rows=columns="<< A.Rows() <<") Zero entries (i=j): "<< endl;
               cout.setf(ios::scientific);
               long prec = cout.precision(15U);
               for ( size_t i=0U; i < A.Rows(); i++ )
                 if ( std::fabs(A(i,i)) < std::numeric_limits<double64>::epsilon() )
                   cout <<"\n\t"<< i <<": "<< A(i,i);
               cout << endl;
               cout.unsetf( ios::scientific );
               cout.precision(prec);
               A.Out();
               throw underflow_error("CompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal.");
            }
        
          // setting matrix such that diagonal element is at the beginning of next row 
          ia[i+1U]       = static_cast<int32>(n); 
          // inserting the diagonal elements at the beginning of each row
          const uint32 istart = static_cast<uint32>(ia[i]);
          int32 jatemp   = ja[ istart ];
          double64 atemp = a[ istart ];
          int32 dindex   = diag;
          ja[ istart ]   = ja[ dindex ];
          a[ istart]     = a[ dindex ];
          a[ dindex ]    = atemp;
          ja[ dindex ]   = jatemp;
       }

     // converting C++ array indices (0..n-1) into Fortran indices (1..n) 
     for ( vector<int32>::iterator it=ia.begin(); it!=ia.end(); it++ ) (*it)++;
     for ( vector<int32>::iterator it=ja.begin(); it!=ja.end(); it++ ) (*it)++;

}  // end Initialize






/**
 
Initialises the public CompressedRowMatrix vectors ia, ja, a for given 
SparseMatrix in case the Point-based approach is selected.
*/
void CompressedRowMatrix::InitializePointBased( const SparseMatrix& A, size_t nsys ) 
 {
      // resize internal storage
      ia.resize( (A.Rows() + 1U) );
      ja.resize( A.Entries() );
      a.resize( ja.size() );
   
      map<size_t,double64>::const_iterator rit;
      long      i, j, k, row;
      int32     diag;
      bool      zero_diag_element(false);
	  
      std::vector<int32>  temp( ja.size() ); // auxilary vector

      const size_t nnu_(A.Rows());
	    for ( k = 0U; k < nnu_; k++)
	      temp[k] = static_cast<int32>(k%(nnu_/nsys)*nsys+k/(nnu_/nsys));
	  		
      for ( i=j=0U, ia[0]=0; i < nnu_; i++ )
       {
	        row = i%nsys*(nnu_/nsys)+i/nsys; // amending the order rows will be written in a[]
          for ( diag=-1, rit=A.RowBegin(row); rit!=A.RowEnd(row); rit++ )
            {
               a[j]  = (*rit).second;
               ja[j] = temp[(*rit).first];
               // rit.first points to matrix entries indexed from 0..rows-1
               if ( ja[j] == temp[row] ) diag = static_cast<int32>(j);
               j++;
            }
          if ( diag == -1 ) zero_diag_element = true;
          
          ia[i+1] = static_cast<int32>(j);
			
		      // inserting the diagonal elements at the beginning of each row
          const int32 istart = ia[i];
          ja[static_cast<size_t>(istart)] = ja[ static_cast<size_t>(diag) ];
          a[static_cast<size_t>(istart)]  = a[ static_cast<size_t>(diag) ];
          a[static_cast<size_t>(diag)]    = a[ static_cast<size_t>(istart) ];
          ja[static_cast<size_t>(diag)]   = ja[ static_cast<size_t>(istart) ];
       }

      if ( zero_diag_element ) {
          cout <<"\nCompressedRowMatrix::Initialize: Error: Zero value(s) in matrix diagonal: ";
          throw underflow_error("SparseM atrix::OutCompressedRowFormat");
       }
	   
     // converting C array indices (0..n-1) into Fortran indices (1..n) 
     for ( auto i : ia ) i++;
     for ( auto j : ja ) j++;

}  // end InitializePointBased









/** Outputs matrix to screen.
*/
void CompressedRowMatrix::Out() const
 {
    cout << flush <<"\nCompressedRowMatrix::Out: "<< endl;
    cout <<"\nrow index vector 'ia' with size = "<<ia.size()<<"\n";
    for ( vector<int32>::const_iterator it=ia.begin(); it!=ia.end(); it++ ) 
      cout << *it <<" ";
    cout <<"\ncolumn index vector 'ja' with size = "<<ja.size()<<"\n";
    for ( vector<int32>::const_iterator it=ja.begin(); it!=ja.end(); it++ ) 
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

