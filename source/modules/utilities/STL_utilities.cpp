#include "STL_utilities.h"

using namespace std;

namespace csmp {


void printVector( const std::vector<uint32_t>& v )
 {
     cout <<"\nuint vector: size/capacity: "<< v.size() <<"/"<< v.capacity() <<"members:"<< endl;
     for ( size_t i{0U}; i<v.size(); i++ )
       cout << v[i] <<" ";
     cout << endl;
 }



void printVector( const std::vector<double>& v )
 {
     cout <<"\ndouble vector: size/capacity: "<< v.size() <<"/"<< v.capacity() <<"members:"<< endl;
     for ( size_t i{0U}; i<v.size(); i++ )
       cout << v[i] <<" ";
     cout << endl;
 }


/** Prints argument vector of vectors to screen for examination.

@section arguments Input Arguments 

The target vector.  
*/
void printVector( const vector<vector<double> >& array )
 {
    vector<vector<double> >::const_iterator    it;
    vector<double>::const_iterator             tit;
    size_t n;

    cout.setf(ios::scientific);
    long prec = cout.precision(2L);
    
    for ( n=0, it=array.begin(); it!=array.end(); it++, n++ ) {
         for ( tit=(*it).begin(); tit!=(*it).end(); tit++ ) cout << (*tit) <<" "; 
         if ( n < 3 ) cout <<"  ";
         else {
              cout << endl;
              n = 0;
           }
      }
      
    cout << endl;  
    cout.unsetf( ios::scientific );
    cout.precision(prec);

 } // end printVector





/**

Prints argument vector of vectors to screen for examination.

@section arguments Input Arguments 

The target vector.  
*/
void printVector( const char* headline, const vector<vector<double> >& array )
 {
    vector<vector<double> >::const_iterator    it;
    vector<double>::const_iterator             tit;
    size_t n;

    cout <<"\n"<< headline << endl;
    cout.setf(ios::scientific);
    long prec = cout.precision(2L);
    
    for ( n=0U, it=array.begin(); it!=array.end(); it++, n++ ) {
         for ( tit=(*it).begin(); tit!=(*it).end(); tit++ ) cout << (*tit) <<" "; 
         if ( n < 3U ) cout <<"  ";
         else {
              cout << endl;
              n = 0U;
           }
      }
      
    cout << endl;  
    cout.unsetf( ios::scientific );
    cout.precision(prec);

 } // end Print



/**
      reads vector<vector> from filestream where the elements of the vector are sequential
*/
template<typename T>
void readVectorOfVectors( ifstream& ifs, size_t total_items, size_t entries_per_vector, deque<vector<T> >& file_records )
 {
    file_records.clear();
    // not in deque: file_records.reserve( total_items / entries_per_vector );
    // reading plist
    size_t item=0U;
    while ( item < total_items )
      {
         // node ID's in file range 0...nodes-1
         vector<T> data;
         data.reserve( entries_per_vector );
         int32_t id;
         for ( size_t i{0U}; i<entries_per_vector; ++i ) {
              ifs >> id;
              assert( id >= 0 && id << total_items ); // assumption that there are not more nodes that elements*nodes_per_element
              data.push_back( id );
              item++;
           }
         file_records.emplace_back( data );
      }

    if ( item != total_items )
      throw std::range_error( "readVectorOfVectors: File record of vector<vector<typename>> was not correctly read." );

 } // end read inlined vector

template void readVectorOfVectors( ifstream&, size_t, size_t, deque<vector<uint32_t> >& );


} // end namespace csmp
