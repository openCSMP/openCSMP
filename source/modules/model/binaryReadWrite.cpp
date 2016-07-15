#include "binaryReadWrite.h"

using namespace std;

namespace csmp {

bool skm_C_fwrite( FILE* fp, const char* str )
 {
     if ( fp == NULL ) {
          cout <<"\nskm_C_fwrite (const char*): ERROR: invalid file pointer."<< endl;
          return false;
       }
     // writing the size of the object
     size_t  characters = strlen(str);
     fwrite( (void*) &characters, sizeof(size_t), 1, fp );
     
     // writing the character string
     fwrite( (void*) str, sizeof(char) * characters, 1, fp );
    
     return true;
 }
 
 

bool skm_C_fread( FILE* fp, char str[] )
 {
     if ( fp == NULL ) {
          cout <<"\nskm_C_fread(char[]): ERROR: invalid file pointer."<< endl;
          return false;
       }
     // read size of the record and assert this 
     size_t  characters(0);
     char       buf[200];

     if ( !fread( (void*) &characters, sizeof(size_t), 1, fp ) )
       {
          cout <<"\nskm_C_fread(char[]): ERROR: could not read string length."<< endl;
          return false;
       }
     // reading the character string
     if ( characters != fread( (void*) buf, sizeof(char), characters, fp ) )
       {
          cout <<"\nskm_C_fread(char[]) ERROR: incorrect number of characters were read: ";
          cout <<"\nIndicated number: "<< characters <<", actual number read: "<< strlen(buf) << endl;
          return false;
       }
     // null terminate string and copy to 'str' argument
     buf[characters] = '\0';  
     strcpy( str, buf );  
    
     return true;
 } 
 
 
} // end namespace csmp
 
 
 
 
 
 
 
 
 
 
 
 
 
