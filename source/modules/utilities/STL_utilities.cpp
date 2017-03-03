#include "STL_utilities.h"

using namespace std;

namespace csmp {


void printVector( const std::vector<size_t>& v )
 {
     cout <<"\nuint vector: size/capacity: "<< v.size() <<"/"<< v.capacity() <<"members:"<< endl;
     for ( size_t i=0U; i<v.size(); i++ )
       cout << v[i] <<" ";
     cout << endl;
 }



void printVector( const std::vector<double>& v )
 {
     cout <<"\ndouble vector: size/capacity: "<< v.size() <<"/"<< v.capacity() <<"members:"<< endl;
     for ( size_t i=0U; i<v.size(); i++ )
       cout << v[i] <<" ";
     cout << endl;
 }


/** Prints argument vector of vectors to screen for examination.

@section arguments Input Arguments 

The target vector.  
*/
void printVector( const vector<vector<double64> >& array )
 {
    vector<vector<double64> >::const_iterator    it;
    vector<double64>::const_iterator             tit;
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
void printVector( const char* headline, const vector<vector<double64> >& array )
 {
    vector<vector<double64> >::const_iterator    it;
    vector<double64>::const_iterator             tit;
    size_t n;

    cout <<"\n"<< headline << endl;
    cout.setf(ios::scientific);
    long prec = cout.precision(2L);
    
    for ( n=0, it=array.begin(); it!=array.end(); it++, n++ ) {
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

} // end namespace csmp
