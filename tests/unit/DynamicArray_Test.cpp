#include "DynamicArray_Test.h"
#include "vector"
#include <iomanip>
#include <iostream>
#include "DynamicArray2D.h"
#include "DynamicArray3D.h"


using namespace std;

namespace csmp {

void DynamicArray_Test::run()
 {
    Test_DynamicArray2D();
    Test_DynamicArray3D();
  
  } // end run


void DynamicArray_Test::Test_DynamicArray2D()
  {
    size_t rows{3u}, cols{2u};
    
    // default constructor
    DynamicArray2D<double>  a( rows, cols );
    // setting all values to 4.
    for ( auto& it : a ) it = 4.;
    if ( verbose_ ) {
         cout <<"\n"<<"Dynamic array 'a' initialised to '4.' using range loop:"<< endl;
         a.out();
      }

    // make a copy of 'a'
    DynamicArray2D<double>  b = a;
    // set all values of b to 7.
    for ( auto& rit : b ) rit = 7.;
    if ( verbose_ ) {
         cout <<"\n"<<"Dynamic array 'b'. Values changed to 7"<< endl;
         b.out();
      }
    // comparitor operator==
    _test( a != b );
    
    // C-array for reference (3(i) x 4(j))
     int c[3][4] = { { 1, 2, 3, 4 },
                     { 5, 6, 7, 8 },
                     { 9, 10, 11, 12 } };

     DynamicArray2D<int>  ca( 3, 4 );
     ca(0,0) =  1.; ca(0,1) =  2.; ca(0,2) =  3.; ca(0,3) =  4.;
     ca(1,0) =  5.; ca(1,1) =  6.; ca(1,2) =  7.; ca(1,3) =  8.;
     ca(2,0) =  9.; ca(2,1) = 10.; ca(2,2) = 11.; ca(2,3) = 12.;
    // dimensions of 'ca'
    _test( ca.rows()  == 3u );
    _test( ca.cols()  == 4u );
     
     // testing that the arrays are equal
     for (int i = 0; i < ca.rows(); i++)
       for (int j = 0; j < ca.cols(); j++)
         _test( ca(i,j) == c[i][j] );
     

     if ( verbose_ ) {
          cout <<"\n"<<"Dynamic array 'ca' mimicking C array:"<< endl;
          ca.out();
          cout <<"\n"<<"C-array for reference:\n";
          const int crows{3}, ccols{4};
          for (int i = 0; i < crows; i++) {
               printf("row %d:\n", i );
               for (int j = 0; j < ccols; j++) printf("%d ", c[i][j]);
               printf("\n");
             }
          printf("\n");
          
          cout <<"\n"<<"printing array 'ca' in reverse order:\n";
          int counter=0;
          for ( auto it=ca.rbegin(); it!=ca.rend(); ++it ) {
               cout <<" "<< (*it);
               counter++;
               if ( counter == ccols ) {
                    cout << endl;
                    counter = 0;
                 }
            }
         }
 
     // testing initialiser list
     DynamicArray2D<int> d{ { 11, 12, 13 },
                            { 21, 22, 23 },
                            { 31, 32, 33 } };  // initializer_list
    _test( d(1,1) == 22 );
    _test( d(2,0) == 31 );
    // 1d-iterators
    if ( verbose_ ) {
        cout <<"\n"<<"printing array (3x3):\n";
        d.out();
      }
      
    // provoking error
#ifdef DEBUG
    try {
        d.at(3,15);
     }
    catch( out_of_range& e ) {
        if ( verbose_ )
          cout <<"\n"<<"DynamicArray2D::at: '"<< e.what() <<"' threw correct exception"<< endl;
     }
#endif

    // resize operation array 'ca'
    rows  = ca.rows();
    cols  = ca.cols();
    ca.resize( rows+1, cols+1 );
    _test( ca.size()     == (rows+1) * (cols+1) );
    _test( ca.capacity() == (rows+1) * (cols+1) );
    if ( verbose_ ) {
         cout <<"\n"<<"array 'ca', increased by 1 row and 1 column\n";
         ca.out();
      }

    DynamicArray2D<double> e;
    _test( e.empty() );

    // data swapping (ca = 5 x 6)
    ca.swap( ca );
    _test( ca.rows() == 5 );
    _test( ca.cols() == 4 );
    swap( a, b );
    _test( a(0,0) == 7. );
    
 } // end Test_DynamicArray2D







void DynamicArray_Test::Test_DynamicArray3D()
  {
    size_t depth{4u}, rows{3u}, cols{2u};
    // default-insert rows*cols*pillars values
    DynamicArray3D<double>  a( depth, rows, cols );
    // setting all values to 4. ?
    for ( auto& it : a ) it = 4.;
    if ( verbose_ ) {
         cout <<"\n"<<"Dynamic array 'a' initialised to '4.' using range loop:"<< endl;
         a.out();
      }

    // copy initialized matrix rows*cols
    DynamicArray3D<double>  b( depth, rows, cols, 7. );
    if ( verbose_ ) {
         cout <<"\n"<<"Dynamic array 'b' constructed with 7'ths:"<< endl;
         b.out();
      }
  
    // comparitor operator==
    _test( a != b );
    
    // C-array for reference (3(i) x 3(j) x 2 (depth)
     int     c[2][3][4] = { { { 1, 2, 3, 4 },
                              { 5, 6, 7, 8 },
                              { 9, 10, 11, 12 } },
                            { { 13, 14, 15, 16 },
                              { 17, 18, 19, 20 },
                              { 21, 22, 23, 24 } } };

     DynamicArray3D<int>  ca( 2, 3, 4 );
     // slice 1
     ca(0,0,0) =  1; ca(0,0,1) =  2; ca(0,0,2) =  3; ca(0,0,3) =  4;
     ca(0,1,0) =  5; ca(0,1,1) =  6; ca(0,1,2) =  7; ca(0,1,3) =  8;
     ca(0,2,0) =  9; ca(0,2,1) = 10; ca(0,2,2) = 11; ca(0,2,3) = 12;
     // slice 2
     ca(1,0,0) = 13; ca(1,0,1) = 14; ca(1,0,2) = 15; ca(1,0,3) = 16;
     ca(1,1,0) = 17; ca(1,1,1) = 18; ca(1,1,2) = 19; ca(1,1,3) = 20;
     ca(1,2,0) = 21; ca(1,2,1) = 22; ca(1,2,2) = 23; ca(1,2,3) = 24;
     
     // testing
     for (int i = 0; i < ca.depth(); i++)
       for (int j = 0; j < ca.rows(); j++)
         for (int k = 0; k < ca.cols(); k++)
           _test( ca(i,j,k) == c[i][j][k] );
 
     if ( verbose_ ) {
         cout <<"\n"<<"Dynamic array 'ca' mimicking C array:"<< endl;
         ca.out();
         cout <<"\n"<<"C-array for reference:\n";
         const int cdepth{2}, crows{3}, ccols{4};
         for (int i = 0; i < cdepth; i++) {
            printf("Depth %d:\n", i );
            for (int j = 0; j < crows; j++) {
                for (int k = 0; k < ccols; k++) {
                    printf("%d ", c[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
          }
        cout <<"\n"<<"printing array 'ca' in reverse order:\n";
        int counter=0;
        for ( auto it=ca.rbegin(); it!=ca.rend(); ++it ) {
             cout <<" "<< (*it);
             counter++;
             if ( counter == 3 ) {
                  cout << endl;
                  counter = 0;
               }
          }
       }
 
    // dimensions of 'ca'
    _test( ca.depth() == 2U );
    _test( ca.rows()  == 3u );
    _test( ca.cols()  == 4u );
 
    /**
       matrix constructor taking vector-of-vectors-style initialiser list
       supporting initialisation by:
    */
    DynamicArray3D<int> d{ { { 111, 121, 131 },
                             { 211, 221, 231 },
                             { 311, 321, 331 }, },
                           { { 112, 122, 132 },
                             { 212, 222, 232 },
                             { 312, 322, 332 } } };  // initializer_list
    // testing
    _test( d.depth() == 2 );
    _test( d.rows() == 3 );
    _test( d.cols() == 3 );
    // individual elements
    _test( d(0,1,0) == 211 );
    _test( d(0,2,1) == 321 );
    _test( d(1,2,2) == 332 );

    if ( verbose_ ) {
        cout <<"\n"<<"printing array constructed from initialiser list (2x3x3):\n";
        d.out();
      }

    // provoking error
    d.at(3,15,2);

    // resize operation array 'ca'
    depth = ca.depth();
    rows  = ca.rows();
    cols  = ca.cols();
    ca.resize( depth, rows+1, cols+1 );
    _test( ca.size()     == (rows+1) * (cols+1) * (depth) );
    _test( ca.capacity() == (rows+1) * (cols+1) * (depth) );
    if ( verbose_ ) {
         cout <<"\n"<<"array 'ca', increased by 1 row and 1 column\n";
         ca.out();
      }
     
     // testing the layout of the array in memory
     {
       DynamicArray3D<int>  ia( 2, 3, 4 );
         // slice 1, row 1 and following
         ia(0,0,0) =  1; ia(0,0,1) =  1; ia(0,0,2) =  1; ia(0,0,3) =  1;
         ia(0,1,0) =  2; ia(0,1,1) =  2; ia(0,1,2) =  2; ia(0,1,3) =  2;
         ia(0,2,0) =  3; ia(0,2,1) =  3; ia(0,2,2) =  3; ia(0,2,3) =  3;
         // slice 2
         ia(1,0,0) =  4; ia(1,0,1) =  4; ia(1,0,2) =  4; ia(1,0,3) =  4;
         ia(1,1,0) =  5; ia(1,1,1) =  5; ia(1,1,2) =  5; ia(1,1,3) =  5;
         ia(1,2,0) =  6; ia(1,2,1) =  6; ia(1,2,2) =  6; ia(1,2,3) =  6;
         
         // correspoding vector
         if ( verbose_ ) {
             ia.out();
             cout <<"\ncontent of 'ia's storage vector: \n";
             int dep{0}, row{0}, col{0};
             for ( const auto& it : ia ) {
                  cout <<" "<< it;
                  if ( col++ == ia.cols()-1 ) {
                       cout << endl;
                       col = 0;
                       row++;
                    }
                  if ( row == ia.rows() ) {
                       cout << endl;
                       row = 0;
                       dep++;
                    }
                  if ( dep == ia.rows() * ia.cols() ) cout << endl;
               }
             cout << endl;
           }
        /*
            stores column after column, depth 0..1
        */
           
     } // end memory layout test

    DynamicArray3D<double> e;
    _test( e.empty() );

    // data swapping
    d.swap( d );
    swap( a, b ); // TODO: test
    
  
  } // end run

} // end csmp

