#include "binaryReadWrite.h"
#include "Exception.h"

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
     fwrite( (void*) str, sizeof(char), characters, fp );
    
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
     char    buf[INFO_STRING];

     if ( !fread( (void*) &characters, sizeof(size_t), 1, fp ) )
       {
          cout <<"\nskm_C_fread(char[]): ERROR: could not read string length."<< endl;
          return false;
       }

     if (characters > INFO_STRING) {
       throw csmp::Exception(ERROR, "skm_C_fread", "Binary file appears to be corrupt");
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



  BinaryFileSectionRead::BinaryFileSectionRead(FILE* fp, const char* header)
      : fp_(fp)
  {
    char readhdr[CSMP_BINARY_FILE_HDR_SIZE];
    size_t hdrlen = strlen(header);
    assert(hdrlen <= CSMP_BINARY_FILE_HDR_SIZE);
    memset(hdr_, 0, sizeof(hdr_));
    memcpy(hdr_, header, std::min(hdrlen, (size_t)CSMP_BINARY_FILE_HDR_SIZE));
    fread(readhdr, sizeof(char), CSMP_BINARY_FILE_HDR_SIZE, fp);
    fread(&offset_, sizeof(offset_), 1, fp);
    assert(!memcmp(hdr_, readhdr, CSMP_BINARY_FILE_HDR_SIZE));
    sectoffset_ = ftell(fp);
  }

  BinaryFileSectionRead::~BinaryFileSectionRead()
  {
      long off = ftell(fp_);
      if (off != offset_) {
        throw csmp::Exception(FATAL_ERROR, "BinaryFileSectionRead", hdr_, "Binary file appears to be corrupt");
      }
  }

  BinaryFileSectionWrite::BinaryFileSectionWrite(FILE* fp, const char* header)
      : fp_(fp)
  {
    char hdr[8];
    size_t hdrlen = strlen(header);
    assert(hdrlen <= sizeof(hdr));
    memset(hdr, 0, sizeof(hdr));
    memcpy(hdr, header, std::min(hdrlen, sizeof(hdr)));
    fwrite(hdr, sizeof(char), sizeof(hdr) / sizeof(char), fp);
    offset_ = ftell(fp);
    fwrite(&offset_, sizeof(offset_), 1, fp);
  }

  BinaryFileSectionWrite::~BinaryFileSectionWrite()
  {
      long off = ftell(fp_);
      fseek(fp_, offset_, SEEK_SET);
      fwrite(&off, sizeof(off), 1, fp_);
      fseek(fp_, off, SEEK_SET);
  }
 
 
} // end namespace csmp
