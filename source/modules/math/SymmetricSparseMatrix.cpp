#include "SymmetricSparseMatrix.h"
#include <fstream>
#include <ctime>


// Detect what FPU operations are available.

// AVX2 gives 256 bit registers
#ifdef __AVX2__
#define CSMP_AVX2_AVAILABLE
#include <immintrin.h>
#endif

// SSE 4.1 gives pcmpeqq, i.e. _mm_cmpeq_epi64
#ifdef __SSE4_1__
#define CSMP_SSE4_1_AVAILABLE
#include <smmintrin.h>
#endif

// SSE2 is the minimum one can expect on x86-64
#ifdef __SSE2__
#define CSMP_SSE2_AVAILABLE
#include <emmintrin.h>
#endif


using namespace std;

namespace csmp {


SymmetricSparseMatrix::Builder::Builder( size_t m_x_n )
    : m_x_n_(m_x_n)
{
}
  
  SymmetricSparseMatrix::Builder::~Builder(  )
  {
  }
  
  bool operator==(const SymmetricSparseMatrixEdge& lhs, const SymmetricSparseMatrixEdge& rhs)
  {
    return lhs.i == rhs.i && lhs.j == rhs.j;
  }


void
SymmetricSparseMatrix::Builder::AddElement( size_t i, size_t j )
{
    pattern_.insert( SymmetricSparseMatrixEdge(i,j) );
    pattern_.insert( SymmetricSparseMatrixEdge(j,i) );
}


struct SymmetricSparseMatrix::Impl {

    struct RowDescriptor {
        size_t* index_begin_;
        size_t* index_end_;
        double64* elmt_begin_;
    };

    size_t n_;

    std::vector<RowDescriptor> rows_;
    std::vector<size_t> indices_;
    std::vector<double64> elements_;

    Impl( const Builder& builder );

    double64* At( size_t i, size_t j )
    {
        auto& row = rows_[i];
        for (size_t* c = row.index_begin_; c != row.index_end_; ++c) {
            if (*c == j) {
                return &row.elmt_begin_[c - row.index_begin_];
            }
        }
        return 0;
    }

    const double64* At( size_t i, size_t j ) const
    {
        auto& row = rows_[i];
        for (size_t* c = row.index_begin_; c != row.index_end_; ++c) {
            if (*c == j) {
                return &row.elmt_begin_[c - row.index_begin_];
            }
        }
        throw 0;
    }

};


#if 0
SymmetricSparseMatrix::SymmetricSparseMatrix( Builder& builder )
{
    Reset(builder);
}
#endif


SymmetricSparseMatrix::~SymmetricSparseMatrix()
{
}

double64
SymmetricSparseMatrix::operator()( size_t i, size_t j ) const
{
    auto elmt = pimpl_->At(i, j);
    assert( elmt );
    return *elmt;
}


double64&
SymmetricSparseMatrix::operator()( size_t i, size_t j )
{
    auto elmt = pimpl_->At(i, j);
    assert( elmt );
    return *elmt;
}


double64
SymmetricSparseMatrix::At( size_t i, size_t j ) const
{
    const double* elmt = pimpl_->At(i, j);
    assert( elmt );
    return *elmt;
}


double64&
SymmetricSparseMatrix::At( size_t i, size_t j )
{
    double64* elmt = pimpl_->At(i, j);
    assert( elmt );
    return *elmt;
}



#if 0
void
SymmetricSparseMatrix::Reset( Builder& builder )
{
}
#endif

void
SymmetricSparseMatrix::Zero()
{
    for (auto& e : pimpl_->elements_) {
        e = 0.0;
    }

}

#if 0
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
#endif





#if 0
/**

IN ORDER TO IDENTIFY EACH ELEMENT IN A, THE USER HAS TO 
PROVIDE TWO integer ARRAYS IA AND JA. IF NNU DENOTES THE TOTAL 
NUMBER OF UNKNOWNS, THE NON-ZERO ENTRIES OF ANY ROW I OF L
(1<=I<=NNU) ARE STORED IN A(J) WHERE THE RANGE OF J
IS GIVEN BY 

             IA(I) <= J <= IA(I+1)-1.

THUS, IA(I) POINTS TO THE POSITION OF THE DIAGONAL ENTRY OF
ROW I WITHIN THE VECTOR A. IN PARTICULAR,

             IA(1) = 1 ,  IA(NNU+1) = 1 + NNA

WHERE NNA DENOTES THE TOTAL NUMBER OF MATRIX ENTRIES STORED.
THE POINTER VECTOR JA HAS TO BE DEFINED SUCH THAT
ANY ENTRY A(J) CORRESPONDS TO THE UNKNOWN U(JA(J)), I.E.,
JA(J) identifies THE COLUMN INDEX OF A(J).
IN PARTICULAR, A(IA(I)) IS THE DIAGONAL ENTRY OF ROW I
AND CORRESPONDS TO THE UNKNOWN U(I): JA(IA(I))=I.

IN THIS TERMINOLOGY, THE I-TH EQUATION READS AS FOLLOWS
(FOR ANY I WITH  1.LE.I.LE.NNU):

          F(I) =        SUM      A(J) * U(JA(J))
          
                   J1 <= J <= J2

WHERE F(I) DENOTES THE I-TH COMPONENT OF THE RIGHT HAND
SIDE AND

              J1 = IA(I) ,  J2 = IA(I+1)-1. 

NOTES: THE ENTRY IA(NNU+1) HAS TO POINT TO THE FIRST FREE
ENTRY IN VECTORS A AND JA, RESPECTIVELY. OTHERWISE, 
AMG CANNOT KNOW THE LENGTH OF THE LAST MATRIX ROW.

THE INPUT VECTORS A, IA AND JA ARE CHANGED BY AMG1R5
SO, AFTER RETURN FROM AMG1R5, THE PACKAGE MUST NOT
BE CALLED A SECOND TIME WITHOUT HAVING NEWLY DEFINED
THE INPUT VECTORS AND USING ISWTCH=4. OTHERWISE, THE
SETUP PHASE WILL FAIL.
*/
void SymmetricSparseMatrix::OutCompressedRowFormat1_n( int32* ia, int32* ja, 
                                              double64* a, bool reallocate ) const
 {
   // memory allocation
   if ( reallocate ) {
        delete[] ia;
        delete[] ja;
        delete[] a;
        // ia indexing is from 1...rows+1
        ia = new int32[ Rows() + 2 ];
        // fortran indexing will be used
        ja = new int32[ Entries() + 1 ];
        a  = new double64[ Entries() + 1 ];
     }

   map<size_t,double64>::const_iterator rit;
   uint32         i;
   int32          j, diag, istart, jatemp;
   double64       atemp;
   bool           zero_diag_element(false);

   for ( i=0, ia[0]=j=0; i<Rows(); i++ )
    {
       for ( diag=-1, rit=data[i].begin(); rit!=data[i].end(); rit++ ) {
            a[j]  = (*rit).second;
            ja[j] = static_cast<int32>((*rit).first);
            if ( ja[j] == static_cast<int32>(i) ) diag = j;
            j++;
         }
       if ( diag == -1 ) zero_diag_element = true;
        
       // putting diagonal element at the beginning of each row 
       ia[i+1]    = j;
       // inserting the diagonal elements at the beginning of each row
       istart     = ia[i];
       jatemp     = ja[istart];
       atemp      = a[istart];
       ja[istart] = ja[diag];
       a[istart]  = a[diag];
       a[diag]    = atemp;
       ja[diag]   = jatemp;
    }

   if ( zero_diag_element ) {
       cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat_n: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SymmetricSparseMatrix::OutCompressedRowFormat");
    }

   // converting C arrays (0..n-1) to Fortran indexing (1..n) 
   for ( istart=static_cast<int32>(Rows()); istart>=0; istart-- ) 
     ia[istart+1]    = ia[istart]+1;
   for ( istart=static_cast<int32>(Entries())-1; istart>=0; istart-- ) {
        ja[istart+1] = ja[istart]+1;
        a[istart+1]  = a[istart];
     }

} // end OutCompressedRowFormat_n



void SymmetricSparseMatrix::OutCompressedRowFormat( int32* ia, int32* ja, 
                                           double64* a, bool reallocate ) const
 {
   cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Copying sparse matrix to AMG vectors..."<< endl;
 
   clock_t ticks = clock();

   // memory allocation
   if ( reallocate ) {
        delete[] ia;
        delete[] ja;
        delete[] a;
        // ia indexing is from 0...rows
        ia = new int32[ Rows() + 1 ];
        // fortran indexing will be used
        ja = new int32[ entries ];
        a  = new double64[ entries ];
     }

   size_t    i;
   int32     j, diag, istart, jatemp;
   double64  atemp;
   bool      zero_diag_element(false);

   for ( i=0U, ia[0]=j=0; i<Rows(); i++ )
    {
       map<size_t,double64>::const_iterator rit=data[i].begin();
       for ( diag=-1; rit!=data[i].end(); rit++ ) {
            a[j]  = (*rit).second;
            // rit.first points to matrix entries indexed from 0..rows-1
            ja[j] = static_cast<int32>((*rit).first);
            if ( ja[j] == static_cast<int32>(i) ) diag = j;
            j++;
         }
       if ( diag == -1 ) zero_diag_element = true;
        
       // setting matrix such that diagonal element is at the beginning of each row 
       ia[i+1]    = j; 
       // inserting the diagonal elements at the beginning of each row
       istart     = ia[i];
       jatemp     = ja[istart];
       atemp      = a[istart];
       ja[istart] = ja[diag];
       a[istart]  = a[diag];
       a[diag]    = atemp;
       ja[diag]   = jatemp;
    }

   if ( zero_diag_element ) {
       cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SymmetricSparseMatrix::OutCompressedRowFormat");
    }

   // converting C array indices (0..n-1) into Fortran indices (1..n) 
   for ( i=0; i<=Rows(); i++ )   ia[i]++;
   for ( i=0; i<Entries(); i++ ) ja[i]++;

   ticks = clock() - ticks;
  
   cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormat (C-format)



void SymmetricSparseMatrix::OutCompressedRowFormatParallel( std::vector<int32>& ia, std::vector<int32>& ja, 
                                                   std::vector<double64> a, int32 nrhalo, bool reallocate ) const
 {
   cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormatParallel: Copying sparse matrix to AMG vectors..."<< endl;
 
   clock_t ticks = clock();

   // memory allocation
   if ( reallocate ) {
        ia.clear();
        ja.clear();
        a.clear();
        ia.resize(Rows()+1-nrhalo); vector<int32>( ia ).swap( ia );
        ja.resize(entries);         vector<int32>( ja ).swap( ja );
        a.resize(entries);          vector<double64>( a ).swap( a);
     }

   size_t    i(0U), j(0U), istart;
   int32     diag, jatemp;
   double64  atemp;
   bool      zero_diag_element(false);

   ia[0] = 0;

   for ( vector<int32>::const_iterator it1=ia.begin(); it1!=ia.end(); it1++, i++ )
      {
         assert( i < data.size() );
         assert( !data[i].empty() );
         map<size_t,double64>::const_iterator rit =data[i].begin();
         for ( diag=-1; rit!=data[i].end(); rit++ )
           {
              a[j]  = (*rit).second;
              // rit.first points to matrix entries indexed from 0..rows-1
              ja[j] = static_cast<int32>((*rit).first);
              if ( ja[j] == static_cast<int32>(i) ) diag = static_cast<int32>(j);
              j++;
           }
      
        if ( diag == -1 ) zero_diag_element = true;
        
        // setting matrix such that diagonal element is at the beginning of each row 
        ia[i+1]    = static_cast<int32>(j); 
        // inserting the diagonal elements at the beginning of each row
        istart     = static_cast<uint32>(ia[i]);
        jatemp     = ja[istart];
        atemp      = a[istart];
        ja[istart] = ja[ static_cast<uint32>(diag) ];
        a[istart]  = a[ static_cast<uint32>(diag) ];
        a[ static_cast<uint32>(diag) ]  = atemp;
        ja[ static_cast<uint32>(diag) ] = jatemp;
     }

   if ( zero_diag_element ) {
       cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormatParallel: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SymmetricSparseMatrix::OutCompressedRowFormatParallel");
    }

   // converting C array indices (0..n-1) into Fortran indices (1..n) 
   for (std::vector<int32>::iterator it1=ia.begin(); it1!=ia.end(); it1++ ) (*it1)++;
   for (std::vector<int32>::iterator it1=ja.begin(); it1!=ja.end(); it1++ ) (*it1)++;

   ticks = clock() - ticks;
  
   cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormatParallel: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormatParallel (C-format)





/// This is the version which is used for the old AMG solver
void SymmetricSparseMatrix::OutCompressedRowFormat( long* ia, long* ja, 
                                               double64* a, bool reallocate ) const
 {
   cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Copying sparse matrix to AMG vectors..."<< endl;
 
   clock_t ticks = clock();

   // memory allocation for the Fortran-style arrays
   if ( reallocate ) {
        delete[] ia;
        delete[] ja;
        delete[] a;
        // ia indexing is from 0...rows+1
        ia = new long[ Rows() + 1 ];
        // fortran indexing will be used
        ja = new long[ entries ];
        a  = new double64[ entries ];
     }

   size_t    i;
   long      j, diag, istart, jatemp;
   double64  atemp;
   bool      zero_diag_element(false);

   // 1. NB: The AMG will automatically correct the offset from 0..n-1 to 1..n
   // but it will not augment the actual entries in the arrays. Thus, this
   // is done in 2
   for ( i=0U, ia[0]=j=0; i<Rows(); i++ )
    {
       map<size_t,double64>::const_iterator rit=data[i].begin();
       for ( diag=-1; rit!=data[i].end(); rit++ ) {
            a[j]  = (*rit).second;
            ja[j] = static_cast<long>((*rit).first);
            if ( ja[j] == static_cast<long>(i) ) diag = j;
            j++;
         }
       if ( diag == -1 ) {
            if ( !zero_diag_element ) 
              cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Zero diagonal element at i=j: ";
            cout << i <<" ";
            zero_diag_element = true;
         }
        
       // setting matrix such that diagonal element is at the beginning of each row 
       ia[i+1]    = j; 
       // inserting the diagonal elements at the beginning of each row
       istart     = ia[i];
       jatemp     = ja[istart];
       atemp      = a[istart];
       ja[istart] = ja[diag];
       a[istart]  = a[diag];
       a[diag]    = atemp;
       ja[diag]   = jatemp;
    }

   if ( zero_diag_element ) {
       cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SymmetricSparseMatrix::OutCompressedRowFormat");
    }

   // 2. Augmenting C array entries (0..n-1) to Fortran indexing (1..n) 
   for ( i=0; i<=Rows(); i++ )   ia[i]++;
   for ( i=0; i<Entries(); i++ ) ja[i]++;

   ticks = clock() - ticks;
  
   cout <<"\nSymmetricSparseMatrix::OutCompressedRowFormat: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormat (C-format)
 
 
void SymmetricSparseMatrix::In( const char* file_name_without_extension )
 {
      string in_file( file_name_without_extension );  in_file +=".txt";
      ifstream  ifs( in_file.c_str() );
      
      if ( !ifs.is_open() ) {
            cout <<"\n SymmetricSparseMatrix::In: failed to read '"<< in_file <<"'; nothing was done\n";
            return;
        }
        
      // read nrows = ncolums
      size_t    msize;
      double64  fdata;
      
      ifs >> msize;
      Resize( msize );
      
      for ( size_t i=0U; i<msize; i++ )
        for ( size_t j=0U; j<msize; j++ ) {
               ifs >> fdata;
               if ( fdata != 0. ) Add( i, j, fdata );
           }
 }


/**
    Prints matrix dimensions, # of non-zero elements and (non-zero)
    entries row by row.
    Each entry is prefaced by its column index.
*/
void SymmetricSparseMatrix::Out( long digits ) const
 {
    if ( data.empty() ) {
         cerr <<"\nSymmetricSparseMatrix::Out: Matrix is empty."<< endl;
         return;
      }

    long    prec(cout.precision(digits));
    size_t pcols(1);

    cout <<"\nSymmetricSparseMatrix::Out: entries (non-zero elements): "<< Entries();
    cout <<"\nrows: "<< data.size() <<", columns: "<< data.size() << endl;
    
    if ( digits != 0 ) cout.setf(ios::scientific);

    // for all rows
    for ( size_t i=0U; i<data.size(); i++ )
      // for all column entries
      for ( map<size_t,double64>::const_iterator
            ditc=data[i].begin(); ditc!=data[i].end(); ditc++ ) {
         // print the column index
         cout <<"("<< i <<","<< (*ditc).first;
         // prints values with extra spaces to achieve an alignment even if there are
         // negative elements
         if ( (*ditc).second > 0 ) cout <<"):  "<< (*ditc).second <<" ";
         else                      cout <<"): "<< (*ditc).second <<" ";
         // wraps the lines if there are more than 10 entries per line
         if ( pcols == 10 || pcols == data[i].size() ) {
              cout << endl;
              pcols = 0;
           }
         pcols++;
      }
    cout << endl << endl; 

    if ( digits != 0 ) {
         cout.unsetf( ios::scientific );
         cout.precision(prec);
      }
 }



void SymmetricSparseMatrix::Out( const char* file ) const
 {
    ofstream  ofs(file);

    if ( data.empty() ) {
         ofs <<"\nSymmetricSparseMatrix::Out: Matrix is empty."<< endl;
         return;
      }

    long         prec;
    const long   digits(3);
    size_t pcols(1);

    ofs <<"\nSymmetricSparseMatrix::Out: Entries: "<< Entries();
    ofs <<"\nrows: "<< data.size() <<", columns: "<< data.size() << endl;
    
    if ( digits != 0 ) { 
        ofs.setf(ios::scientific);
        prec = ofs.precision(digits);
     }

    for ( size_t i=0; i<data.size(); i++ )
      for ( map<size_t,double64>::const_iterator
            ditc=data[i].begin(); ditc!=data[i].end(); ditc++ ) {
         ofs <<"("<< i <<","<< (*ditc).first;
         ofs <<"): "<< (*ditc).second <<" ";
         if ( pcols == 10 || pcols == data[i].size() ) {
              ofs << endl;
              pcols = 0;
           }
         pcols++;
      }
    ofs << endl << endl; 

    if ( digits != 0 ) {
         ofs.unsetf( ios::scientific );
         ofs.precision(prec);
      }
 }
#endif


} // csmp

