#include "LoopUnroller.h"

using namespace std;

namespace csmp {

LoopUnroller::LoopUnroller() {}

LoopUnroller::~LoopUnroller() {}

/// simple i loop, e.g. data[i] += obj.data[j]
void LoopUnroller::Unroll( const char* file,
                               int isz,
                               const char* var1, const char* var2,
                               const char* expression )
 {
    ofstream ofs(file);
    
    ofs <<"\nLoopUnroller::Unroll: Unrolled loop(i..."<< isz <<")";
    ofs <<", expression: "<< var1;
    ofs <<"[i] "<< expression <<" "<< var2 <<"[i]"<< endl;
    
    for ( int i=0; i<isz; i++ )
      {
         ofs <<"\n  "<< var1 <<"["<< i <<"] "<< expression;
         ofs <<" "<< var2 <<"["<< i <<"];";
      }
      
    cout <<"\nLoopUnroller::Unroll: '"<< file <<"' written successfully.\n";
 }



/// simple i,j loop, e.g. data[i][j] += obj.data[j][j]
void LoopUnroller::Unroll( const char* file,
                               int isz, int jsz,
                               const char* var1, const char* var2,
                               const char* expression )
 {
    ofstream ofs(file);
    
    ofs <<"\nLoopUnroller::Unroll: Unrolled loop(i..."<< isz <<",j..."<< jsz <<")";
    ofs <<", expression: "<< var1;
    ofs <<"[i][j] "<< expression <<" "<< var2 <<"[j][j]"<< endl;
    
    for ( int i=0; i<isz; i++ )
      {
 //        ofs <<"\n\/\/ row "<< i; to insert C++ comment (??/ trigraph for backslash)
         ofs <<"\n// row "<< i;
         for ( int j=0; j<jsz; j++ )
           {
              ofs <<"\n  "<< var1 <<"["<< i <<"]["<< j <<"] "<< expression;
              ofs <<" "<< var2 <<"["<< i <<"]["<< j <<"];";
           }
      }
      
    cout <<"\nLoopUnroller::Unroll: '"<< file <<"' written successfully.\n";
 }



/// version with inbuilt checks for the current size of the matrix
/// assumes that loop variables are called i,j
/// simple i,j loop, e.g. data[i][j] += obj.data[j][j]
void LoopUnroller::Unroll( const char* file,
                               const char* rows, const char* cols,
                               int isz, int jsz,
                               const char* var1, const char* var2,
                               const char* expression )
 {
    ofstream ofs(file);
    
    ofs <<"\nLoopUnroller::Unroll: Unrolled loop(i..."<< isz <<",j..."<< jsz <<")";
    ofs <<", expression: "<< var1;
    ofs <<"[i][j] "<< expression <<" "<< var2 <<"[j][j]"<< endl;
    
    for ( int i=0; i<isz; i++ )
      {
         ofs <<"\nLoopUnroller::Unroll: row "<< i;
         ofs <<"\nif ( i < "<< rows <<" ) {";
         for ( int j=0; j<jsz; j++ )
           {
              ofs <<"\n     if ( j < "<< cols <<" ) "<< var1;
              ofs <<"["<< i <<"]["<< j <<"] "<< expression;
              ofs <<" "<< var2 <<"["<< i <<"]["<< j <<"];";
           }
         ofs <<"\n  }";
      }
      
    cout <<"\nLoopUnroller::Unroll: '"<< file <<"' written successfully.\n";
 }


} // csmp
