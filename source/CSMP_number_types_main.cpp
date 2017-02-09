#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <cmath>


using namespace  std;

int main()
  {
     // 0. Checking the sizes of the inbuilt types
     // ------------------------------------------
     int char_size(sizeof(char)),
         uchar_size(sizeof(unsigned char)),
         short_size(sizeof(short)),
         ushort_size(sizeof(unsigned short)),
         long_size(sizeof(long)),
         ulong_size(sizeof(unsigned long)),
         int_size(sizeof(int)),
         uint_size(sizeof(unsigned int)),
         llong_size(sizeof(long long)),
         ullong_size(sizeof(unsigned long long)),
         ldouble_size(sizeof(long double)),
         float_size(sizeof(float)),
         double_size(sizeof(double));

     // 1. echoing the sizes of the different types to stdout
     // -----------------------------------------------------
     cout <<"\nmain: Sizes of inbuilt number types:  "<< endl;
     cout <<"\nSize of 'char' (bytes):               "<< char_size;
     cout <<"\nSize of 'unsigned char' (bytes):      "<< uchar_size;
     cout <<"\nSize of 'short' (bytes):              "<< short_size;
     cout <<"\nSize of 'unsigned short' (bytes):     "<< ushort_size;
     cout <<"\nSize of 'int' (bytes):                "<< int_size;
     cout <<"\nSize of 'unsigned int' (bytes):       "<< uint_size;
     cout <<"\nSize of 'long' (bytes):               "<< long_size;
     cout <<"\nSize of 'unsigned long' (bytes):      "<< ulong_size;
     cout <<"\nSize of 'float' (bytes):              "<< float_size;
     cout <<"\nSize of 'double' (bytes):             "<< double_size;
     cout <<"\nSize of 'long double' (bytes):        "<< ldouble_size;
     // special types which may or may not be there
     cout <<"\nSize of 'long long' (bytes):          "<< llong_size;
     cout <<"\nSize of 'unsigned long long' (bytes): "<< ullong_size;
     cout << endl;
     cout.flush();

     cout <<"\nSize of the typical STL container::size_type = size_t is: "<< sizeof(size_t);
     cout <<" versus vector<double>::size_type: "<< sizeof(vector<double>::size_type)  << endl;
     cout <<"\n\nUsing 'size_t' for loop integers."<< endl;
     assert( sizeof(size_t) == sizeof(vector<double>::size_type) );
     // if this assertion fails please contact: stephan.matthai@unimelb.edu.au

     // 2. writing a corresponding header file for the
     // CSMP port to this platform
     // --------------------------
     ofstream ofs("CSMP_number_types.h");
     ofs <<"#ifndef CSMP_NUMBER_TYPES_H"<< endl;
     ofs <<"#define CSMP_NUMBER_TYPES_H"<< endl << endl;

     ofs <<"/**"                                                << endl;
     ofs <<"@file CSMP_number_types.h"                          << endl;
     ofs <<"@author S.K. Matthai"                               << endl << endl;
     ofs <<"// ------------------------------------------------"<< endl << endl;
     ofs <<"// CSMP number types for current computing platform"<< endl << endl;
     ofs <<"// ------------------------------------------------"<< endl << endl;
     ofs <<"*/"                                                 << endl << endl;
     ofs <<"namespace csmp{ "<< endl << endl;

     ofs <<"/**"                                                << endl;
     ofs <<"@addtogroup CSMPglobalTypedefs"                     << endl;
     ofs <<"@{"                                                 << endl;
     ofs <<"*/\n"                                                 << endl;

     // -----------------------------------------------
     // char & short in CSP = 8 & 16 bits
     // -----------------------------------------------
     if ( char_size == 1 )
       {
          ofs <<"typedef char                 char8;"    << endl;
          ofs <<"typedef unsigned char        uchar8;"   << endl;
          ofs <<"typedef char                 int8;"     << endl;
          ofs <<"typedef unsigned char        uint8;"    << endl;
       }
     else
     cerr <<"\nmain: CSMP 'char8' etc. types are not available on this platform."<< endl;

      // ---------------------------------------------------------
     // short & unsigned short = 16 bits
     // --------------------------------
     // CSP integers should be the size of a machine word
     // on a 32 bit machine (2 words on a 64 bit machine)
     // size of loop variables: = int typically 4 bytes (32 bits)
     // ---------------------------------------------------------
     if ( short_size == 2 )
       {
          ofs <<"typedef short                int16;"   << endl;
          ofs <<"typedef unsigned short       uint16;"  << endl;
       }
     else if ( int_size == 2 )
       {
          ofs <<"typedef int                  int16;"   << endl;
          ofs <<"typedef unsigned int         uint16;"  << endl;
       }
     else
     cerr <<"\nmain: CSMP 'int16' and 'uint16' types are not available on this platform."<< endl;

     // ---------------------------------------------------------
     // int & unsigned int = 32 bits
     // ----------------------------
     // CSP integers should be the size of a machine word
     // on a 32 bit machine (2 words on a 64 bit machine)
     // size of loop variables: = int typically 4 bytes (32 bits)
     // ---------------------------------------------------------
     if ( int_size == 4 )
       {
          ofs <<"typedef int                  int32;"   << endl;
          ofs <<"typedef unsigned int         uint32;"  << endl;
       }
     else if ( short_size == 4 )
       {
          ofs <<"typedef short                int32;"   << endl;
          ofs <<"typedef unsigned short       uint32;"  << endl;
       }

     // --------------------------------------------------
     // long & unsigned long = 64 bits
     // --------------------------------------------------
     if ( int_size == 8 )
       {
          ofs <<"typedef int                  long64;"  << endl;
          ofs <<"typedef unsigned int         ulong64;" << endl;
       }
     else if ( long_size == 8 )
       {
          ofs <<"typedef long                 long64;"  << endl;
          ofs <<"typedef unsigned long        ulong64;" << endl;
       }
     else if ( llong_size == 8 )
       {
          ofs <<"typedef long long            long64;"  << endl;
          ofs <<"typedef unsigned long long   ulong64;" << endl;
       }
     else
     cerr <<"\nmain: CSMP 'long64' type is not available on this platform."<< endl;

     // --------------------------------------------------
     // long long & unsigned long long = 128 bits
     // --------------------------------------------------
     if ( long_size == 16 )
       {
          ofs <<"typedef long                 long128;"  << endl;
          ofs <<"typedef unsigned long        ulong128;" << endl;
       }
     else if ( llong_size == 16 ) // if long long is defined
       {
          ofs <<"typedef long long            long128;"  << endl;
          ofs <<"typedef unsigned long long   ulong128;" << endl;
       }
     else
     cerr <<"\nmain: CSMP 'long128' type is not available on this platform."<< endl;

     // --------------------------------------------------
     // float  = 32 bits
     // --------------------------------------------------
     if      ( float_size  == 4 )
       ofs <<"typedef float                float32;" << endl;
     else if ( double_size == 4 )
       ofs <<"typedef double               float32;" << endl;
     else
     cerr <<"\nmain: CSMP 'float32' type is not available on this platform."<< endl;

     // --------------------------------------------------
     // double = 64 bits
     // --------------------------------------------------
     if      ( double_size == 8 )
       ofs <<"typedef double               double64;"  << endl;
     else if ( float_size ==  8 )
       ofs <<"typedef float                double64;"  << endl;
     else
     cerr <<"\nmain: CSMP 'double64' type is not available on this platform."<< endl;

     // --------------------------------------------------
     // long double = 128 bits
     // --------------------------------------------------
     if ( ldouble_size == 16 )
       ofs <<"typedef long double          double128;"  << endl;
     else
     cerr <<"\nmain: CSMP 'double128' type is not available on this platform."<< endl;

     ofs <<"\n/**" << endl;
     ofs <<"@}" << endl;
     ofs << "*/" << endl;

     ofs <<"\n } // end namespace csmp"<< endl << endl;

     ofs <<"\n#endif"<< endl;

     cout <<"\nmain: CSMP header file '"<< "CSMP_number_types.h" <<"' written successfully."<< endl;

     return 0;

  } // end main
