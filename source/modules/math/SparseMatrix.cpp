#include "SparseMatrix.h"
#include "DenseMatrix.h"
#include "Matrix.h"
#include "compareFloats.h"
#include <fstream>
#include <ctime>

using namespace std;

namespace csmp {

SparseMatrix::SparseMatrix()
 {
 }
 

SparseMatrix::SparseMatrix( size_t m_x_n )
 : data_(m_x_n)
 {
 }
 
 
SparseMatrix::SparseMatrix( const SparseMatrix& sp )
 : data_(sp.data_),
   entries_(sp.entries_)
 {
 }
 
 
SparseMatrix::SparseMatrix( SparseMatrix&& sp )
 : data_(move(sp.data_)),
   entries_(move(sp.entries_))
 {
 }
 


SparseMatrix& SparseMatrix::operator=( const SparseMatrix& sp )
 {
    if ( &sp != this ) {
         entries_ = sp.entries_;
         data_    = sp.data_;
      }
    return *this;
 }


SparseMatrix::~SparseMatrix()
 {
 }


void SparseMatrix::Resize( size_t n_x_m, bool preserve_allocated_memory )
 {
    data_.resize( n_x_m );
    if ( !preserve_allocated_memory )
      data_.shrink_to_fit();
 }



size_t SparseMatrix::Entries() const 
 { 
    return entries_;
 }


/// iterator over the matrix rows
vector<map<size_t,double> >::const_iterator  SparseMatrix::Begin() const
 {
    return data_.begin();
 }


/// iterator over the rows of the matrix
vector<map<size_t,double> >::const_iterator  SparseMatrix::End() const
 {
    return data_.end();
 }


map<size_t,double>::const_iterator  SparseMatrix::RowBegin( size_t i ) const
 {
    assert( i < Rows() );
    return data_[i].begin();
 }


map<size_t,double>::const_iterator  SparseMatrix::RowEnd( size_t i ) const
 {
    assert( i < Rows() );
    return data_[i].end();
 }

size_t SparseMatrix::Rows() const
  {
     return data_.size();
  }


size_t SparseMatrix::Cols() const
  {
     return data_.size();
  }


double SparseMatrix::operator()( size_t i, size_t j ) const
 {
    assert( i < data_.size() );
    assert( j < data_.size() );
    auto ditc = data_[i].find(j);
    if ( ditc == data_[i].end() ) return 0.;
    return (*ditc).second;
 }
  
    
double SparseMatrix::At( size_t i, size_t j ) const
 {
    if ( i >= data_.size() ) {
         cerr <<"\nSparseMatrix::At("<< i <<","<< j <<") const: ";
         cerr <<"ERROR, Row access index greater or equal to rows ("<< Rows() <<"). Returning: NaN"<< endl;
         return numeric_limits<double>::signaling_NaN();
      }
    if ( j >= data_.size() ) {
         cerr <<"\nSparseMatrix::At("<< i <<","<< j <<") const: ";
         cerr <<"ERROR, Column access index greater or equal to cols ("<< Cols() <<"). Returning: NaN"<< endl;
         return numeric_limits<double>::signaling_NaN();
      }
    if ( data_[i].find(j) == data_[i].end() ) {
         cerr <<"\nSparseMatrix::At("<< i <<","<< j <<") const: ";
         cerr <<"ERROR, Matrix entry does not exist. Returning: 0."<< endl;
        return 0.;
      }
    return (*this)(i,j);
 }




void SparseMatrix::RemoveHalo( int nrhalo )
 {
   for ( size_t row=data_.size()-nrhalo; row!=data_.size(); row++ )
     ZeroRow(row);
 }

/// remove a single element from the matrix
void SparseMatrix::RemoveEntry( size_t i, size_t j )
 {
    map<size_t,double>::iterator dit = data_[i].find(j);
    if ( dit != data_[i].end() ) {
         data_[i].erase( dit );
         entries_--;
      }
 }

void SparseMatrix::Erase()
 {
    data_.erase( data_.begin(), data_.end() );
    entries_ = 0;
 }


void SparseMatrix::Zero()
{
    Erase();
}


// deletes all entries from the row
void SparseMatrix::ZeroRow( size_t row )
 {
    entries_ -= data_[row].size();
    data_[row].erase( data_[row].begin(), data_[row].end() );
 }


void SparseMatrix::ZeroColumn( size_t col )
 {
    for ( auto& it : data_ ) {
        auto dit=it.find(col);
        if ( dit != it.end() ) {
             it.erase( dit );
             entries_--;
          }
     }
 }


void SparseMatrix::Assign( size_t i, size_t j, double val )
 {
    if ( !(val != 0.) ) return;
    
    if ( i >= data_.size() ) {
         cerr <<"\nSparseMatrix::Assign("<< i <<","<< j <<"): ";
         cerr <<"Row index out of range."<< endl;
         throw range_error("SparseMatrix::Assign");
      }
    if ( j >= data_.size() ) {
         cerr <<"\nSparseMatrix::Assign("<< i <<","<< j <<"): ";
         cerr <<"Column index out of range."<< endl;
         throw range_error("SparseMatrix::Assign");
      }
    
    auto it = data_[i].insert( make_pair( j, val ) );
    // if insertion was successful, there are now more non-zero elements
    if ( it.second ) entries_++;
    else  // value is overwritten
      (*it.first).second = val;
 }




/// more efficient additions to matrix than attainable with operator()
void SparseMatrix::Add( size_t i, size_t j, double val )
 {
    // zero elements are not stored
    if ( !(val != 0.) ) return;

    if ( i >= data_.size() ) {
         cerr <<"\nSparseMatrix::Add("<< i <<","<< j <<","<< val <<"): ";
         cerr <<"Row index out of range."<< endl;
         throw range_error("SparseMatrix::Add");
      }
    if ( j >= data_.size() ) {
         cerr <<"\nSparseMatrix::Add("<< i <<","<< j <<","<< val <<"): ";
         cerr <<"Column index out of range."<< endl;
         throw range_error("SparseMatrix::Add");
      }
 
    // if matrix element does not exist this statement inserts it
    auto addit = data_[i].insert( make_pair(j,val) );
    // if matrix element already exists, the value is added to it
    if ( !addit.second ) {
        (*addit.first).second += val;
         // if a < b = false and a > b = false then a == b
         if ( !((*addit.first).second > 0.) &&
              !((*addit.first).second < 0.) ) {
             data_[i].erase( (*addit.first).first );
             entries_--;
          }
      }
    else entries_++;
 }



/// Operator accumulates sparse matrices.  Initially created for OpenMP features.
/// All need to have the same number of row entries. Column entries may be differenct
/// @author Julian E. Mindel
SparseMatrix& SparseMatrix::operator+=( const SparseMatrix& mat )
  {
    if ( this->Rows() != mat.Rows() ) {
        cerr <<" Added sparse matrices must have the same number of rows (even if they are empty)."<< endl;
        throw range_error("SparseMatrix::operator+=");
    }

    rowsIterator recipient_iterator=data_.begin();
    for ( rowsConstIterator rit = mat.Begin(); rit!= mat.End();rit++){
        for (colsConstIterator cit = rit->begin(); cit!= rit->end();cit++){
            // if matrix element does not exist in the recipient this statement inserts it
            auto addit = recipient_iterator->insert( make_pair(cit->first,cit->second) );
            // if matrix element already exists, the value is added to it
            if ( !addit.second ) {
                (*addit.first).second += cit->second;
                if ( !((*addit.first).second > 0.) &&
                     !((*addit.first).second < 0.) ) {
                    recipient_iterator->erase( (*addit.first).first );
                    entries_--;
                }
            }
            else entries_++;
        }
        recipient_iterator++;
    }

    return *this;
  }



/// multiply an element of the matrix with val; if result is zero, the element is eliminated
void SparseMatrix::MultiplyEntryWith( size_t i, size_t j, double val )
 {
    if ( i >= data_.size() ) {
         cerr <<"\nSparseMatrix::MultiplyEntryWith("<< i <<","<< j <<","<< val <<"): ";
         cerr <<"Row index out of range."<< endl;
         throw range_error("SparseMatrix::MultiplyEntryWith");
      }
 
    map<size_t,double>::iterator dit = data_[i].find(j);
    // if the element does not exist, i.e. is zero, the result would also be zero
    // so nothing needs to be done
    if ( dit == data_[i].end() ) return;
    
    // multiplying entry with 'val'
    (*dit).second *= val;
    
    // removing entry if it has become zero by multiplication with 'val'
    if ( !((*dit).second > 0.) &&
         !((*dit).second < 0.) ) {
          data_[i].erase( j );
          entries_--;
      }
}


void SparseMatrix::MultiplyWith( const vector<double>& vec, vector<double>& res )
  {
    if ( vec.size() != Cols() ) {
         cerr <<"\nSparseMatrix::MultiplyWith: first argument vector must have the same size as matrix columns. ";
         cerr <<"Row index out of range."<< endl;
         throw range_error("SparseMatrix::MultiplyWith");
      }
  
  res.resize(Rows());
  vector<double>( res ).swap( res );
  fill(res.begin(), res.end(), 0. );
  
  size_t i{0U};
  for ( const auto& it : data_ ) {
      for ( const auto& dit : it )
        res[i] += dit.second * vec[ dit.first ];
      i++;
    }
}



/// j's of cols with data
void SparseMatrix::ColumnIndices( size_t row, vector<size_t>& indices ) const
 {
    if ( row >= data_.size() ) {
         cerr <<"\nSparseMatrix::ColumnIndices("<< row <<",..): ";
         cerr <<"Row index out of range."<< endl;
         throw range_error("SparseMatrix::ColumnIndices");
      }
      
    indices.resize( data_[row].size() );

    size_t i{0U};
    for ( const auto& ditc : data_[row] )
      indices[i++] = ditc.first;
 }



size_t SparseMatrix::RecountEntries() const
 {
    size_t current_entries{0U};

    for ( const auto& row : data_ )
      current_entries += row.size();

    if ( current_entries != entries_ ) {
         cout <<"\nSparseMatrix::RecountEntries: private 'entries_' is not uptodate anymore (old vs. recount): ";
         cout << entries_ <<" vs. "<< current_entries << endl;
      }

    return current_entries;

 } // end RecountEntries



double SparseMatrix::InfinityNorm() const
 {
    double maxsum{0.};

    // loop over all rows
    for ( const auto& vit : data_ ) {
         double sum{0.};
         // loop over all columns, adding up entries
         for ( const auto& it : vit ) sum += fabs(it.second);  // take maximum rowsum
         maxsum = max( sum, maxsum );
      }
    return maxsum;
 }



bool SparseMatrix::Symmetric() const
 {
    for ( size_t i{0U}; i<data_.size(); i++ )
      {
          if ( data_[i].empty() ) {
               cerr <<"\nSparseMatrix::Symmetric: Matrix contains zero rows."<< endl;
               return false;
            }
          for ( const auto& ditc : data_[i] )
            // only for non-diagonal elements
            if ( i != ditc.first ) {
                  // checking whether there is a matrix element with exchanged indices
                  if ( data_[ ditc.first ].find(i) != data_[ ditc.first ].end() ) {
                       // if such an element exists a comparison of element values is made
                       if ( (*this)( ditc.first, i ) != ditc.second ) return false;
                    }
                  else return false;
              }
      }
    return true;
    
 } // end Symmetric



bool SparseMatrix::ZeroesInDiagonal() const
 {
    if ( entries_ < Rows() ) return true;
    
    size_t i{0U};
    for ( const auto& it : data_ ) {
         if ( it.empty() ||
              it.find(i) == it.end() ||
              essentiallyEqual( (*it.find(i)).second, 0. ) ) {
                //cout <<"\nSparseMatrix::ZeroesInDiagonal: Zero entry at ("<< i <<","<< i <<").";
                return true;
            }
         i++;
      }
    return false;
 }


/// if there are negative elements in the diagonal of the matrix
bool SparseMatrix::DiagonallyPositive() const
 {
    for ( size_t i{0U}; i<data_.size(); i++ ) {
         auto diagonal_elmt = data_[i].find(i);
         if ( diagonal_elmt == data_[i].end() ) {
               cerr <<"\nSparseMatrix::DiagonallyPositive: Warning: ";
               cerr <<" This test can't be performed since there are zeroes in matrix diagonal."<< endl;
               return true;
            }
         if ( (*diagonal_elmt).second < 0. ) return false;
      }
    
    return true;
 }



void SparseMatrix::SparsityPattern( const char* txtfile ) const
 {
    if ( data_.empty() ) {
         cerr <<"\nSparseMatrix::SparsityPattern: Matrix is empty."<< endl;
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
    for ( uint32_t i{0U}; i<idx.Rows(); ++i )
      for ( uint32_t j{0U}; j<idx.Cols(); ++j )
        Add( idx(i,j), idx(i,j), d(i,j) );

 } // end Assign

template void SparseMatrix::Assign<DenseMatrix<DM4>,Matrix>( const DenseMatrix<DM4>&, const Matrix& );


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
   int32_t  j, diag, istart, jatemp;
   double   atemp;
   bool     zero_diag_element(false);

   size_t i{0U};
   for ( ia[0]=j=0; i<Rows(); i++ )
    {
       for ( diag=-1, rit=data_[i].begin(); rit!=data_[i].end(); rit++ ) {
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
       cerr <<"\nSparseMatrix::OutCompressedRowFormat_n: Error: Zero value(s) in matrix diagonal: ";
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
        ja = new int32_t[ entries_ ];
        a  = new double[ entries_ ];
     }

   int32_t j, diag, istart, jatemp;
   double  atemp;
   bool    zero_diag_element(false);

   size_t i{0U};
   for ( ia[0]=j=0; i<Rows(); i++ )
    {
       map<size_t,double>::const_iterator rit=data_[i].begin();
       for ( diag=-1; rit!=data_[i].end(); rit++ ) {
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
       cerr <<"\nSparseMatrix::OutCompressedRowFormat: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SparseMatrix::OutCompressedRowFormat");
    }

   // converting C array indices (0..n-1) into Fortran indices (1..n) 
   for ( i=0; i<=Rows(); i++ )   ia[i]++;
   for ( i=0; i<Entries(); i++ ) ja[i]++;

   ticks = clock() - ticks;
  
   cout <<"\nSparseMatrix::OutCompressedRowFormat: Copied matrix to vectors. ";
   cout <<"CPU ticks used for operation: "<< ticks << endl << endl;

} // end OutCompressedRowFormat (C-format)






void SparseMatrix::OutCompressedRowFormatParallel( vector<int32_t>& ia, vector<int32_t>& ja,
                                                   vector<double> a, int32_t nrhalo, bool reallocate ) const
 {
   cout <<"\nSparseMatrix::OutCompressedRowFormatParallel: Copying sparse matrix to AMG vectors..."<< endl;
 
   clock_t ticks = clock();

   // memory allocation
   if ( reallocate ) {
        ia.clear();
        ja.clear();
        a.clear();
        ia.resize(Rows()+1-nrhalo); vector<int32_t>( ia ).swap( ia );
        ja.resize(entries_);        vector<int32_t>( ja ).swap( ja );
        a.resize(entries_);         vector<double>( a ).swap( a);
     }

   size_t   i(0U), j(0U), istart;
   int32_t  diag, jatemp;
   double   atemp;
   bool     zero_diag_element(false);

   ia[0] = 0;

   for ( auto it1=ia.begin(); it1!=ia.end(); it1++, i++ )
      {
         assert( i < data_.size() );
         assert( !data_[i].empty() );
         map<size_t,double>::const_iterator rit =data_[i].begin();
         for ( diag=-1; rit!=data_[i].end(); rit++ )
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
       cerr <<"\nSparseMatrix::OutCompressedRowFormatParallel: Error: Zero value(s) in matrix diagonal: ";
       throw underflow_error("SparseMatrix::OutCompressedRowFormatParallel");
    }

   // converting C array indices (0..n-1) into Fortran indices (1..n) 
   for ( auto it1=ia.begin(); it1!=ia.end(); it1++ ) (*it1)++;
   for ( auto it1=ja.begin(); it1!=ja.end(); it1++ ) (*it1)++;

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
        ja = new long[ entries_ ];
        a  = new double[ entries_ ];
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
       map<size_t,double>::const_iterator rit=data_[i].begin();
       for ( diag=-1; rit!=data_[i].end(); rit++ ) {
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
       cerr <<"\nSparseMatrix::OutCompressedRowFormat: Error: Zero value(s) in matrix diagonal: ";
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
            cerr <<"\n SparseMatrix::In: failed to read '"<< in_file <<"'; nothing was done\n";
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
    if ( data_.empty() ) {
         cerr <<"\nSparseMatrix::Out: Matrix is empty."<< endl;
         return;
      }

    long    prec(cout.precision(digits));
    size_t pcols(1);

    cout <<"\nSparseMatrix::Out: entries (non-zero elements): "<< Entries();
    cout <<"\nrows: "<< data_.size() <<", columns: "<< data_.size() << endl;
    
    if ( digits != 0 ) cout.setf(ios::scientific);

    // for all rows
    for ( size_t i{0U}; i<data_.size(); i++ )
      // for all column entries
      for ( auto ditc=data_[i].begin(); ditc!=data_[i].end(); ditc++ ) {
         // print the column index
         cout <<"("<< i <<","<< (*ditc).first;
         // prints values with extra spaces to achieve an alignment even if there are
         // negative elements
         if ( (*ditc).second > 0 ) cout <<"):  "<< (*ditc).second <<" ";
         else                      cout <<"): "<< (*ditc).second <<" ";
         // wraps the lines if there are more than 10 entries per line
         if ( pcols == 10 || pcols == data_[i].size() ) {
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

    if ( data_.empty() ) {
         ofs <<"\nSparseMatrix::Out: Matrix is empty."<< endl;
         return;
      }

    long         prec;
    const long   digits(3);
    size_t pcols(1);

    ofs <<"\nSparseMatrix::Out: Entries: "<< Entries();
    ofs <<"\nrows: "<< data_.size() <<", columns: "<< data_.size() << endl;
    
    if ( digits != 0 ) { 
        ofs.setf(ios::scientific);
        prec = ofs.precision(digits);
     }

    for ( size_t i{0U}; i<data_.size(); i++ )
      for ( auto ditc=data_[i].begin(); ditc!=data_[i].end(); ditc++ ) {
         ofs <<"("<< i <<","<< (*ditc).first;
         ofs <<"): "<< (*ditc).second <<" ";
         if ( pcols == 10 || pcols == data_[i].size() ) {
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

    for ( size_t i{0U}; i<data_.size(); i++)
      for ( auto ditc = data_[i].begin(); ditc != data_[i].end(); ditc++)
        ofs << i + 1 << " " << ((*ditc).first) + 1 << " " << (*ditc).second << "\n";

    if (digits != 0) {
        ofs.unsetf(ios::scientific);
        ofs.precision(prec);
      }
  }

} // csmp

