#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Matrix.h"
#include <fstream>
#include <ctime>

using namespace std;

namespace csmp {

SparseMatrix::SparseMatrix()
 : entries(0)
 {
 }
 

SparseMatrix::SparseMatrix( size_t m_x_n )
 : entries(0), data(m_x_n)
 {
 }
 
 
SparseMatrix::SparseMatrix( const SparseMatrix& sp )
 {
    *this = sp;
 }
 
SparseMatrix& SparseMatrix::operator=( const SparseMatrix& sp )
 {
    if ( &sp != this ) {
         entries = sp.entries;
         data    = sp.data;
      }
    return *this;
 }

SparseMatrix::~SparseMatrix()
 {
 }


void SparseMatrix::Resize( size_t n_x_m, bool preserve_allocated_memory )
 {
    data.resize( n_x_m );
    if ( !preserve_allocated_memory )
      vector<map<size_t,double> >( data ).swap( data );
 }






size_t SparseMatrix::Entries() const 
 { 
    return entries; 
 }

/// iterator over the matrix rows
std::vector<std::map<size_t,double> >::const_iterator  SparseMatrix::Begin() const
 {
    return data.begin();
 }


/// iterator over the rows of the matrix
std::vector<std::map<size_t,double> >::const_iterator  SparseMatrix::End() const
 {
    return data.end();
 }


std::map<size_t,double>::const_iterator  SparseMatrix::RowBegin( size_t i ) const
 {
    assert( i < Rows() );
    return data[i].begin();
 }


std::map<size_t,double>::const_iterator  SparseMatrix::RowEnd( size_t i ) const
 {
    assert( i < Rows() );
    return data[i].end();
 }

size_t SparseMatrix::Rows() const
  {
     return data.size();
  }


size_t SparseMatrix::Cols() const
  {
     return data.size();
  }


double SparseMatrix::operator()( size_t i, size_t j ) const
 {
    if ( i >= data.size() ) {
         std::cerr <<"\nSparseMatrix::operator("<< i <<","<< j <<") const: ";
         std::cerr <<"Row access index out of range."<< std::endl;
         return std::numeric_limits<double>::signaling_NaN();
      }
    std::map<size_t,double>::const_iterator ditc = data[i].find(j);
    if ( ditc == data[i].end() ) return static_cast<double>(0.0);
    return (*ditc).second;
 }
  
    
double SparseMatrix::At( size_t i, size_t j ) const
 {
    return (*this)(i,j);
 }




void SparseMatrix::RemoveHalo( int nrhalo )
 {
   for ( size_t row=data.size()-nrhalo; row!=data.size(); row++ )
     ZeroRow(row);
 }

/// remove a single element from the matrix
void SparseMatrix::RemoveEntry( size_t i, size_t j )
 {
    std::map<size_t,double>::iterator dit = data[i].find(j);
    if ( dit != data[i].end() ) {
         data[i].erase( dit );
         entries--;
      }
 }

void SparseMatrix::Erase()
 {
    data.erase( data.begin(), data.end() );
    entries = 0;
 }

void SparseMatrix::Zero()
{
    Erase();
}

void SparseMatrix::ZeroRow( size_t row )
 {
    entries -= data[row].size();
    data[row].erase( data[row].begin(), data[row].end() );
 }


void SparseMatrix::ZeroColumn( size_t col )
 {
    for ( vector<map<size_t,double> >::iterator it=data.begin(); it!=data.end(); it++ ) {
        map<size_t,double>::iterator dit=(*it).find(col);
        if ( dit != (*it).end() ) {
             (*it).erase( dit );
             entries--;
          }
     }
 }

void SparseMatrix::Assign( size_t i, size_t j, double val )
 {
    if ( i >= data.size() ) {
         std::cout <<"\nSparseMatrix::Assign("<< i <<","<< j <<"): ";
         std::cout <<"Row access index out of range."<< std::endl;
         throw std::range_error("SparseMatrix::Assign");
      }
    // more costly than insert, operator[] overwrites or adds element
    const size_t n = data[i].size();
    data[i][j] = val;
    if ( data[i].size() > n ) entries++;

    // remove entry if val is zero
    if ( !( val > static_cast<double>(0.) || val < static_cast<double>(0.) ) ){
        data[i].erase( j );
        entries--;
    }
 }


/// more efficient additions to matrix than attainable with operator()
void SparseMatrix::Add( size_t i, size_t j, double val )
 {
    // zero elements are not stored
    if ( !( val > static_cast<double>(0.) || val < static_cast<double>(0.) ) ) return; 

    if ( i >= data.size() ) {
         std::cout <<"\nSparseMatrix::Add("<< i <<","<< j <<","<< val <<"): ";
         std::cout <<"Row access index out of range."<< std::endl;
         throw std::range_error("SparseMatrix::Add");
      }
 
    // if matrix element does not exist this statement inserts it
    pair<std::map<size_t,double>::iterator,bool>
      addit = data[i].insert( std::make_pair(j,val) );
    // if matrix element already exists, the value is added to it
    if ( !addit.second ) {
        (*addit.first).second += val;
         if ( !( (*addit.first).second > static_cast<double>(0.) ||
          (*addit.first).second < static_cast<double>(0.) ) ) {
             data[i].erase( (*addit.first).first );
             entries--;
          }
      }
    else entries++;
 }

/// Operator accumulates sparse matrices.  Initially created for OpenMP features.
/// All need to have the same number of row entries. Column entries may be differenct
/// @author Julian E. Mindel
SparseMatrix& SparseMatrix::operator+=( SparseMatrix mat )
  {
    if ( this->Rows() != mat.Rows() ) {
        std::cout <<" Added sparse matrices must have the same number of rows (even if they are empty)."<< std::endl;
        throw std::range_error("SparseMatrix::operator+=");
    }

    rowsIterator recipient_iterator=this->data.begin();
    for ( rowsConstIterator rit = mat.Begin(); rit!= mat.End();rit++){
        for (colsConstIterator cit = rit->begin(); cit!= rit->end();cit++){
            // if matrix element does not exist in the recipient this statement inserts it
            pair<std::map<size_t,double>::iterator,bool>
                    addit = recipient_iterator->insert( std::make_pair(cit->first,cit->second) );
            // if matrix element already exists, the value is added to it
            if ( !addit.second ) {
                (*addit.first).second += cit->second;
                if ( !( (*addit.first).second > static_cast<double>(0.) ||
                        (*addit.first).second < static_cast<double>(0.) ) ) {
                    recipient_iterator->erase( (*addit.first).first );
                    entries--;
                }
            }
            else entries++;
        }
        recipient_iterator++;
    }

    return *this;
  }

/// multiply an element of the matrix with val; if result is zero, the element is eliminated
void SparseMatrix::MultiplyEntryWith( size_t i, size_t j, double val )
 {
    std::map<size_t,double>::iterator dit = data[i].find(j);
    // if the element does not exist, i.e. is zero, the result would also be zero
    // so nothing needs to be done
    if ( dit == data[i].end() ) return;
    // removing element which becomes zero when multiplied with 'val'
    if ( std::fabs(val) <= std::numeric_limits<double>::epsilon() ) {
          data[i].erase( j );
          entries--;
         return;
      }
    (*dit).second *= val;
}


void SparseMatrix::MultiplyWith( const std::vector<double>& vec, std::vector<double>& res ) {
  assert(vec.size() == Cols());
  
  res.resize(Rows());
  vector<double>( res ).swap( res );
  fill(res.begin(), res.end(), 0.);
  
  size_t  i(0U);
  for ( vector<map<size_t,double> >::iterator it = data.begin(); it != data.end(); ++it, ++i )
    for ( map<size_t,double>::iterator dit = it->begin(); dit != it->end(); ++dit )
      res[i] += dit->second * vec[dit->first];
}



/// j's of cols with data
void SparseMatrix::ColumnIndices( size_t row, vector<uint32_t>& indices ) const
 {
    size_t i(0);

    indices.resize( data[row].size() );
    vector<uint32_t>( indices ).swap( indices );

    for ( map<size_t,double>::const_iterator
          ditc=data[row].begin(); ditc!=data[row].end(); ditc++ )
      indices[i++] = (*ditc).first;
 }

size_t SparseMatrix::RecountEntries() const
 {
    size_t current_entries(0U);

    for ( auto n=0U; n<data.size(); n++ )
      for ( map<size_t,double>::const_iterator
            rit=data[n].begin(); rit!=data[n].end(); rit++ )
        current_entries++;

    if ( current_entries != entries ) {
         cout <<"\nSparseMatrix::RecountEntries: Entries was not uptodate anymore (old vs. new): ";
         cout << entries <<" vs. "<< current_entries << endl;
      }

    return current_entries;

 } // end RecountEntries



double SparseMatrix::InfinityNorm() const
 {
    std::map<size_t,double>::const_iterator it;
    double sum(0.0), maxsum(0.0);

    // loop over all rows
    for ( size_t i{0U}; i<data.size(); i++ ) {
        sum = 0.;
        // loop over all columns and add up entries
        for ( it = data[i].begin(); it != data[i].end(); it++ ) {
            sum += std::abs((*it).second);
          }
        // take maximum rowsum
        if ( sum > maxsum ) maxsum = sum;
      }
    return maxsum;
 }


bool SparseMatrix::Symmetric() const
 {
    for ( size_t i{0U}; i<data.size(); i++ ) 
      {
          if ( data[i].empty() ) {
               cout <<"\nSparseMatrix::Symmetric: Matrix contains zero rows."<< endl;
               return false;
            }
          for ( map<size_t,double>::const_iterator ditc=data[i].begin(); ditc!=data[i].end(); ditc++ )
            // only for non-diagonal elements
            if ( i != (*ditc).first ) {
                  // checking whether there is a matrix element with exchanged indices
                  if ( data[ (*ditc).first ].find(i) != data[ (*ditc).first ].end() ) {
                       // if such an element exists a comparison of element values is made
                       if ( (*this)( (*ditc).first, i ) != (*ditc).second ) return false;
                    }
                  else return false;
              }
      }
    return true;
    
 } // end Symmetric

bool SparseMatrix::ZeroesInDiagonal() const
 {
    if ( entries < Rows() ) return true;
    
    map<size_t,double>::const_iterator  ditc;
    for ( size_t i{0U}; i<data.size(); i++ )
      if ( data[i].empty() ||
          (ditc=data[i].find(i)) == data[i].end() ||
          (*ditc).second == static_cast<double>(0.) ) {
           cout <<"\nSparseMatrix::ZeroesInDiagonal: Zero entry at ("<< i <<","<< i <<").";
           return true;
        }
    
    return false;
 }


/// if there are negative elements in the diagonal of the matrix
bool SparseMatrix::DiagonallyPositive() const
 {
    map<size_t,double>::const_iterator  ditc;
    for ( size_t i{0U}; i<data.size(); i++ ) {
         if ( (ditc=data[i].find(i)) == data[i].end() ) {
               cout <<"\nSparseMatrix::DiagonallyPositive: Warning: ";
               cout <<" This test can't be performed since there are zeroes in matrix diagonal."<< endl;
               return true;
            }
         if ( (*ditc).second < static_cast<double>(0.) ) return false;
      }
    
    return true;
 }

void SparseMatrix::SparsityPattern( const char* txtfile ) const
 {
    if ( data.empty() ) {
         cout <<"\nSparseMatrix::SparsityPattern: Matrix is empty."<< endl;
         return;
      }

    string  fname(txtfile);
    fname +=".txt";
    ofstream  ofs(fname.c_str());

     for ( size_t i{0U}; i<Rows(); i++ )
       {
          for ( size_t j{0U}; j<Cols(); j++ )
            if ( At(i,j) != 0. ) ofs << 1 <<" ";
            else                 ofs << 0 <<" ";
          ofs << endl;
       }
     ofs << endl;
     ofs.close();

    ofs <<"\nSparseMatrix::SparsityPattern: Matrix written successfully to: "<< fname << endl;

 } // end SparsityPattern





/// inserts (DenseMatrix) 'data' in sparse matrix at location identified by position matrix 'idx'
template<class cspMat1, class cspMat2>
void SparseMatrix::Assign( const cspMat1& idx, const cspMat2& d )
 {
    for ( size_t i{0U}; i<idx.Rows(); ++i )
      for ( size_t j{0U}; j<idx.Cols(); ++j )
        Add( static_cast<uint32_t>(idx(i,j)), static_cast<uint32_t>(idx(i,j)), d(i,j) );

 } // end Assign

//template
//void SparseMatrix::Assign<DenseMatrix<DM_MAX>,DenseMatrix<DM_MAX> >(
//const DenseMatrix<DM_MAX>&,
//const DenseMatrix<DM_MAX>& );

//template
//void SparseMatrix::Assign<DenseMatrix<DM_MIN>,DenseMatrix<DM_MIN> >(
//const DenseMatrix<DM_MIN>&,
//const DenseMatrix<DM_MIN>& );

template
void SparseMatrix::Assign<DenseMatrix<DM4>,Matrix >(
const DenseMatrix<DM4>&,
const Matrix& );


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
void SparseMatrix::OutCompressedRowFormat1_n( int32_t*  ia, int32_t*  ja, 
                                              double* a, bool reallocate ) const
 {
   // memory allocation
   if ( reallocate ) {
        delete[] ia;
        delete[] ja;
        delete[] a;
        // ia indexing is from 1...rows+1
        ia = new int32_t[ Rows() + 2 ];
        // fortran indexing will be used
        ja = new int32_t[ Entries() + 1 ];
        a  = new double[ Entries() + 1 ];
     }

   map<size_t,double>::const_iterator rit;
   uint32_t         i;
   int32_t          j, diag, istart, jatemp;
   double       atemp;
   bool           zero_diag_element(false);

   for ( i=0, ia[0]=j=0; i<Rows(); i++ )
    {
       for ( diag=-1, rit=data[i].begin(); rit!=data[i].end(); rit++ ) {
            a[j]  = (*rit).second;
            ja[j] = static_cast<int32_t>((*rit).first);
            if ( ja[j] == static_cast<int32_t>(i) ) diag = j;
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
       cout <<"\nSparseMatrix::OutCompressedRowFormat_n: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SparseMatrix::OutCompressedRowFormat");
    }

   // converting C arrays (0..n-1) to Fortran indexing (1..n) 
   for ( istart=static_cast<int32_t>(Rows()); istart>=0; istart-- ) 
     ia[istart+1]    = ia[istart]+1;
   for ( istart=static_cast<int32_t>(Entries())-1; istart>=0; istart-- ) {
        ja[istart+1] = ja[istart]+1;
        a[istart+1]  = a[istart];
     }

} // end OutCompressedRowFormat_n



void SparseMatrix::OutCompressedRowFormat( int32_t*  ia, int32_t*  ja, 
                                           double* a, bool reallocate ) const
 {
   cout <<"\nSparseMatrix::OutCompressedRowFormat: Copying sparse matrix to AMG vectors..."<< endl;
 
   clock_t ticks = clock();

   // memory allocation
   if ( reallocate ) {
        delete[] ia;
        delete[] ja;
        delete[] a;
        // ia indexing is from 0...rows
        ia = new int32_t[ Rows() + 1 ];
        // fortran indexing will be used
        ja = new int32_t[ entries ];
        a  = new double[ entries ];
     }

   size_t    i;
   int32_t     j, diag, istart, jatemp;
   double  atemp;
   bool      zero_diag_element(false);

   for ( i=0U, ia[0]=j=0; i<Rows(); i++ )
    {
       map<size_t,double>::const_iterator rit=data[i].begin();
       for ( diag=-1; rit!=data[i].end(); rit++ ) {
            a[j]  = (*rit).second;
            // rit.first points to matrix entries indexed from 0..rows-1
            ja[j] = static_cast<int32_t>((*rit).first);
            if ( ja[j] == static_cast<int32_t>(i) ) diag = j;
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
       cout <<"\nSparseMatrix::OutCompressedRowFormat: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SparseMatrix::OutCompressedRowFormat");
    }

   // converting C array indices (0..n-1) into Fortran indices (1..n) 
   for ( i=0; i<=Rows(); i++ )   ia[i]++;
   for ( i=0; i<Entries(); i++ ) ja[i]++;

   ticks = clock() - ticks;
  
   cout <<"\nSparseMatrix::OutCompressedRowFormat: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormat (C-format)



void SparseMatrix::OutCompressedRowFormatParallel( std::vector<int32_t>& ia, std::vector<int32_t>& ja, 
                                                   std::vector<double> a, int32_t nrhalo, bool reallocate ) const
 {
   cout <<"\nSparseMatrix::OutCompressedRowFormatParallel: Copying sparse matrix to AMG vectors..."<< endl;
 
   clock_t ticks = clock();

   // memory allocation
   if ( reallocate ) {
        ia.clear();
        ja.clear();
        a.clear();
        ia.resize(Rows()+1-nrhalo); vector<int32_t>( ia ).swap( ia );
        ja.resize(entries);         vector<int32_t>( ja ).swap( ja );
        a.resize(entries);          vector<double>( a ).swap( a);
     }

   size_t    i(0U), j(0U), istart;
   int32_t     diag, jatemp;
   double  atemp;
   bool      zero_diag_element(false);

   ia[0] = 0;

   for ( vector<int32_t>::const_iterator it1=ia.begin(); it1!=ia.end(); it1++, i++ )
      {
         assert( i < data.size() );
         assert( !data[i].empty() );
         map<size_t,double>::const_iterator rit =data[i].begin();
         for ( diag=-1; rit!=data[i].end(); rit++ )
           {
              a[j]  = (*rit).second;
              // rit.first points to matrix entries indexed from 0..rows-1
              ja[j] = static_cast<int32_t>((*rit).first);
              if ( ja[j] == static_cast<int32_t>(i) ) diag = static_cast<int32_t>(j);
              j++;
           }
      
        if ( diag == -1 ) zero_diag_element = true;
        
        // setting matrix such that diagonal element is at the beginning of each row 
        ia[i+1]    = static_cast<int32_t>(j); 
        // inserting the diagonal elements at the beginning of each row
        istart     = static_cast<uint32_t>(ia[i]);
        jatemp     = ja[istart];
        atemp      = a[istart];
        ja[istart] = ja[ static_cast<uint32_t>(diag) ];
        a[istart]  = a[ static_cast<uint32_t>(diag) ];
        a[ static_cast<uint32_t>(diag) ]  = atemp;
        ja[ static_cast<uint32_t>(diag) ] = jatemp;
     }

   if ( zero_diag_element ) {
       cout <<"\nSparseMatrix::OutCompressedRowFormatParallel: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SparseMatrix::OutCompressedRowFormatParallel");
    }

   // converting C array indices (0..n-1) into Fortran indices (1..n) 
   for (std::vector<int32_t>::iterator it1=ia.begin(); it1!=ia.end(); it1++ ) (*it1)++;
   for (std::vector<int32_t>::iterator it1=ja.begin(); it1!=ja.end(); it1++ ) (*it1)++;

   ticks = clock() - ticks;
  
   cout <<"\nSparseMatrix::OutCompressedRowFormatParallel: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormatParallel (C-format)





/// This is the version which is used for the old AMG solver
void SparseMatrix::OutCompressedRowFormat( long* ia, long* ja, 
                                               double* a, bool reallocate ) const
 {
   cout <<"\nSparseMatrix::OutCompressedRowFormat: Copying sparse matrix to AMG vectors..."<< endl;
 
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
        a  = new double[ entries ];
     }

   size_t    i;
   long      j, diag, istart, jatemp;
   double  atemp;
   bool      zero_diag_element(false);

   // 1. NB: The AMG will automatically correct the offset from 0..n-1 to 1..n
   // but it will not augment the actual entries in the arrays. Thus, this
   // is done in 2
   for ( i=0U, ia[0]=j=0; i<Rows(); i++ )
    {
       map<size_t,double>::const_iterator rit=data[i].begin();
       for ( diag=-1; rit!=data[i].end(); rit++ ) {
            a[j]  = (*rit).second;
            ja[j] = static_cast<long>((*rit).first);
            if ( ja[j] == static_cast<long>(i) ) diag = j;
            j++;
         }
       if ( diag == -1 ) {
            if ( !zero_diag_element ) 
              cout <<"\nSparseMatrix::OutCompressedRowFormat: Zero diagonal element at i=j: ";
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
       cout <<"\nSparseMatrix::OutCompressedRowFormat: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SparseMatrix::OutCompressedRowFormat");
    }

   // 2. Augmenting C array entries (0..n-1) to Fortran indexing (1..n) 
   for ( i=0; i<=Rows(); i++ )   ia[i]++;
   for ( i=0; i<Entries(); i++ ) ja[i]++;

   ticks = clock() - ticks;
  
   cout <<"\nSparseMatrix::OutCompressedRowFormat: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormat (C-format)






 
 
 
 
 
void SparseMatrix::In( const char* file_name_without_extension )
 {
      string in_file( file_name_without_extension );  in_file +=".txt";
      ifstream  ifs( in_file.c_str() );
      
      if ( !ifs.is_open() ) {
            cout <<"\n SparseMatrix::In: failed to read '"<< in_file <<"'; nothing was done\n";
            return;
        }
        
      // read nrows = ncolums
      size_t    msize;
      double  fdata;
      
      ifs >> msize;
      Resize( msize );
      
      for ( size_t i{0U}; i<msize; i++ )
        for ( size_t j{0U}; j<msize; j++ ) {
               ifs >> fdata;
               if ( fdata != 0. ) Add( i, j, fdata );
           }
 }


/**
    Prints matrix dimensions, # of non-zero elements and (non-zero)
    entries row by row.
    Each entry is prefaced by its column index.
*/
void SparseMatrix::Out( long digits ) const
 {
    if ( data.empty() ) {
         cerr <<"\nSparseMatrix::Out: Matrix is empty."<< endl;
         return;
      }

    long    prec(cout.precision(digits));
    size_t pcols(1);

    cout <<"\nSparseMatrix::Out: entries (non-zero elements): "<< Entries();
    cout <<"\nrows: "<< data.size() <<", columns: "<< data.size() << endl;
    
    if ( digits != 0 ) cout.setf(ios::scientific);

    // for all rows
    for ( size_t i{0U}; i<data.size(); i++ )
      // for all column entries
      for ( map<size_t,double>::const_iterator
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



void SparseMatrix::Out( const char* file ) const
 {
    ofstream  ofs(file);

    if ( data.empty() ) {
         ofs <<"\nSparseMatrix::Out: Matrix is empty."<< endl;
         return;
      }

    long         prec;
    const long   digits(3);
    size_t pcols(1);

    ofs <<"\nSparseMatrix::Out: Entries: "<< Entries();
    ofs <<"\nrows: "<< data.size() <<", columns: "<< data.size() << endl;
    
    if ( digits != 0 ) { 
        ofs.setf(ios::scientific);
        prec = ofs.precision(digits);
     }

    for ( size_t i{0U}; i<data.size(); i++ )
      for ( map<size_t,double>::const_iterator
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

void SparseMatrix::OutForMatlab(const char* file) const
{
	ofstream  ofs(file);

	long         prec;
	const long   digits(3);

	if (digits != 0) {
		ofs.setf(ios::scientific);
		prec = ofs.precision(digits);
	}

	for (auto i = 0; i<data.size(); i++)
		for (map<size_t, double>::const_iterator
			ditc = data[i].begin(); ditc != data[i].end(); ditc++) {
			ofs << i + 1 << " " << ((*ditc).first) + 1 << " " << (*ditc).second << "\n";
		}

	if (digits != 0) {
		ofs.unsetf(ios::scientific);
		ofs.precision(prec);
	}
}
} // csmp

