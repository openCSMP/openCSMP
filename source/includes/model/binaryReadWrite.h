#ifndef CSMP_READ_WRITE_H
#define CSMP_READ_WRITE_H

#include "CSMP_definitions.h"

// STL vectors & deques

namespace csmp {

    /**
    @file binaryReadWrite.h
    */

#define CSMP_BINARY_FILE_HDR_SIZE  8

/// Helper class to read sections from a binary file
class BinaryFileSectionRead
{
public:
    BinaryFileSectionRead(FILE* fp, const char* header);

    ~BinaryFileSectionRead();

private:
    char hdr_[CSMP_BINARY_FILE_HDR_SIZE+1];
    FILE* fp_;

    // fseek()/ftell() returns long
    long offset_, sectoffset_;
};


/// Helper class to write sections to a binary file
class BinaryFileSectionWrite {
public:
    BinaryFileSectionWrite(FILE* fp, const char* header);

    ~BinaryFileSectionWrite();

private:
    FILE* fp_;
    long offset_;
};
    

    /**
    @addtogroup CSMPglobalFunctions
    */


template<class T>
bool skm_C_fwrite( std::FILE* fp, const std::vector<T>& stl_ctner );

template<class T>
bool skm_C_fread( std::FILE* fp, std::vector<T>& stl_ctner );

template<typename T>
bool skm_C_fwrite( std::FILE* fp, const std::deque<T>& stl_ctner );

template<typename T>
bool skm_C_fread( std::FILE* fp, std::deque<T>& stl_ctner );

template<class T>
bool skm_C_fwrite( std::FILE* fp, const std::deque<std::vector<T> >& stl_ctner );

template<class T>
bool skm_C_fread( std::FILE* fp, std::deque<std::vector<T> >& stl_ctner );

// maps

template<class M, class T>
bool skm_C_fwrite( std::FILE* fp, const std::map<M,T>& stl_ctner );
  
template<class M, class T>
bool skm_C_fwrite( std::FILE* fp, const std::unordered_map<M,T>& stl_ctner );

template<class M, class T>
bool skm_C_fread( std::FILE* fp, std::map<M,T>& stl_ctner );

// maps of vectors

template<class M, class T>
bool skm_C_fwrite( std::FILE* fp, const std::map<M,std::vector<T> >& stl_ctner );

template<class M, class T>
bool skm_C_fread( std::FILE* fp, std::map<M,std::vector<T> >& stl_ctner );


// character strings

bool skm_C_fwrite( std::FILE* fp, const char* str );

bool skm_C_fread( std::FILE* fp, char str[] );



/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL deque to a binary file. The data segment is preceded by an
'unsigned long' number which determines the number of objects which are 
written to file. 

@section arguments Input Arguments 

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL deque which will be written to file.  

@return returns a boolean indicating whether all data have been
written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite(). 

@section application Application

To efficiently write deque data to a binary file. 

@section messages Messages 

If the file pointer is invalid, method will quit, reporting an error. 
*/
template<class T>
bool skm_C_fwrite( std::FILE* fp, const std::vector<T>& stl_ctner )
 {
     if ( fp == NULL )
       {
          std::cerr <<"\nbool skm_C_fwrite: ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     typename std::vector<T>::const_iterator  it;
     size_t                     bytes    = sizeof(T);
     size_t                  elements = stl_ctner.size();
     
     // writing the size of the object
     fwrite( (void*) &elements, sizeof(size_t), 1, fp );
     
     // writing all elements
     for ( it=stl_ctner.begin(); it!=stl_ctner.end(); it++ )
       std::fwrite( (void*) &(*it), bytes, 1, fp );
    
     return true;
 }
 
/// deque version
template<typename T>
bool skm_C_fwrite( std::FILE* fp, const std::deque<T>& stl_ctner )
 {
     if ( fp == NULL )
       {
          std::cerr <<"\nbool skm_C_fwrite: ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     typename std::deque<T>::const_iterator  it;
     size_t                     bytes    = sizeof(T);
     size_t                  elements = stl_ctner.size();
     
     // writing the size of the object
     fwrite( (void*) &elements, sizeof(size_t), 1, fp );
     
     // writing all elements
     for ( it=stl_ctner.begin(); it!=stl_ctner.end(); it++ )
       std::fwrite( (void*) &(*it), bytes, 1, fp );
    
     return true;
 }
 

/**

Reads data from a binary file into an STL deque which is erased before
adding the data, if it already contained any data. How many 
records are read is specified by an 'unsigned long' number which precedes
the dataset. 

@section arguments Input Arguments 

A pointer to a binary file opened in read binary mode ("rb"), and a 
reference to an STL deque. 

@return returns the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the 
file pointer is invalid or less data are read, the function returns false. 

@section implementation Implementation

Uses the ANSI C function fread(). The deque which will hold the data
is erased if it is not empty, and storage is reserved for the new number 
of elements which is read from file. Then the elements are read and the
deque is returned.  

@section application Application

Efficiently read data into STL deques. 

@section messages Messages 

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this 
number does not match the number of records which were actually read. 
*/
template<class T>
bool skm_C_fread( std::FILE* fp, std::vector<T>& stl_ctner )
 {
     if ( !stl_ctner.empty() )
       stl_ctner.erase( stl_ctner.begin(), stl_ctner.end() );
       
     if ( fp == NULL )
       {
          std::cerr <<"\nbool skm_C_fread: ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     size_t     bytes    = sizeof(T), counter(0);
     size_t  elements(0), i;
     T          val;
     
     // read size of the record and check it 
     if ( !std::fread( (void*) &elements, sizeof(size_t), 1, fp ) ) {
          std::cerr <<"\nbool skm_C_fread: ERROR: could not read record length."<< std::endl;
          return false;
       }
     if ( elements > 0 ) {
          stl_ctner.reserve( elements );
          // writing all elements
          for ( i=0; i<elements; i++ ) {
               // counting the successfully read elements
               counter += std::fread( (void*) &val, bytes, 1, fp );
               stl_ctner.push_back( val );
            }
       }
     if ( counter != elements )
       {
          std::cerr <<"\nbool skm_C_fread: ERROR: incorrect number of records were read: ";
          std::cerr <<"\nIndicated number: "<< elements <<", actual number read: "<< counter << std::endl;
          return false;
       }
     return true;
 }


/// deque version
template<typename T>
bool skm_C_fread( std::FILE* fp, std::deque<T>& stl_ctner )
 {
     if ( !stl_ctner.empty() )
       stl_ctner.erase( stl_ctner.begin(), stl_ctner.end() );
       
     if ( fp == NULL ) {
          std::cerr <<"\nbool skm_C_fread: ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     size_t     bytes    = sizeof(T), counter(0);
     size_t  elements(0), i;
     T          val;
     
     // read size of the record and assert this 
     if ( !std::fread( (void*) &elements, sizeof(size_t), 1, fp ) )
       {
          std::cout <<"\nbool skm_C_fread: ERROR: could not read record length."<< std::endl;
          return false;
       }
     if ( elements > 0 )
       {
          // writing all elements
          for ( i=0; i<elements; i++ )
            {
               // counting the successfully read elements
               counter += std::fread( (void*) &val, bytes, 1, fp );
               stl_ctner.push_back( val );
            }
       }
     if ( counter != elements )
       {
          std::cerr <<"\nbool skm_C_fread: ERROR: incorrect number of records were read: ";
          std::cerr <<"\nIndicated number: "<< elements <<", actual number read: "<< counter << std::endl;
          return false;
       }
     return true;
 }



/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL deque<vector<T> > to a binary file. The data segment is preceded 
by an 'unsigned long' number which determines the number of vector objects 
which are written to file. 

@section arguments Input Arguments 

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL vector of vectors which will be written to 
file.  

@return returns a boolean indicating whether all vector data have 
been written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite(). 

@section application Application

To efficiently write vector data to a binary file. 

@section messages Messages 

If the file pointer is invalid, method will quit, reporting an error. 
*/
template<class T>
bool skm_C_fwrite( std::FILE* fp, const std::deque<std::vector<T> >& stl_ctner )
 {
     if ( fp == NULL ) {
          std::cerr <<"\nbool skm_C_fwrite(const deque<vector<T> >&): ";
          std::cerr <<"ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     
     // writing the number of vector objects
     size_t  elements(stl_ctner.size());
     std::fwrite( (void*) &elements, sizeof(size_t), 1, fp );
     
     // writing all elements
     for ( typename std::deque<std::vector<T> >::const_iterator
           it=stl_ctner.begin(); it!=stl_ctner.end(); it++ )
       skm_C_fwrite( fp, (*it) );
    
     return true;
 }
 
 

/**

Reads a vector of STL vector data from a binary file. The STL container
is erased before adding the data, if it already contained any data. How many 
records are read is specified by an 'unsigned long' number which precedes
the dataset. 

@section arguments Input Arguments 

A pointer to a binary file opened in read binary mode ("rb"), and a 
reference to an STL vector of vector. 

@return returns the boolean 'true' if the number of data records 
which precedes the dataset in the file has been read correctly. If the 
file pointer is invalid or less data are read, the function returns false. 

@section implementation Implementation

Uses the ANSI C function fread(). The vector<vector<T> > which will hold the 
data is erased if it is not empty. Storage is reserved for the new number 
of elements which is read from file. Then the elements are read and the
container is returned.  

@section application Application

Efficiently read data into STL vector of vectors. 

@section messages Messages 

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this 
number does not match the number of records which were actually read. 
*/
template<class T>
bool skm_C_fread( std::FILE* fp, std::deque<std::vector<T> >& stl_ctner )
 {
     if ( !stl_ctner.empty() )
       stl_ctner.erase( stl_ctner.begin(), stl_ctner.end() );
       
     if ( fp == NULL ) {
          std::cout <<"\nbool skm_C_fread(deque<vector<T> >&): ";
          std::cout <<"ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     size_t     counter(0);
     size_t  elements(0), i;
     std::vector<T>  val;
     
     // 1. reading number of vector records and assert this reading
     if ( !std::fread( (void*) &elements, sizeof(size_t), 1, fp ) ) {
          std::cerr <<"\nbool skm_C_fread(deque<vector<T> >&): ";
          std::cerr <<"ERROR: could not read record length."<< std::endl;
          return false;
       }
     if ( elements > 0U ) {
          //stl_ctner.reserve( elements );
          // 2. reading all the vector records
          for ( i=0; i<elements; i++ ) {
               // counting the successfully read elements
               if ( skm_C_fread( fp, val ) ) {
                    counter++;
                    stl_ctner.push_back( val );
                 }
            }
       }
     if ( counter != elements ) {
          std::cerr <<"\nbool skm_C_fread(deque<vector<T> >&): ";
          std::cerr <<"ERROR: incorrect number of records were read: ";
          std::cerr <<"\nIndicated number: "<< elements <<", actual number read: "<< counter << std::endl;
          return false;
       }
       
     return true;
 }
 


/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL map to a binary file. The data segment is preceded 
by an 'unsigned long' number which determines the number of objects 
which are written to file. Each record is preceded by the 
corresponding map key.  

@section arguments Input Arguments 

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL map which will be written to 
file.  

@return a boolean indicating whether all data have
been written correctly to the file. 

@section implementation Implementation

Uses the ANSI standard C function fwrite(). NOTE that the map key must
not contain any dynamically allocated data. Otherwise these data are
sliced of the binary record. 

@section application Application

To efficiently write map data to a binary file. 

@section messages Messages 

If the file pointer is invalid, method will quit, reporting an error. 
*/
template<class M, class T>
bool skm_C_fwrite( std::FILE* fp, const std::map<M,T>& stl_ctner )
 {
     if ( fp == NULL )
       {
          std::cerr <<"\nbool skm_C_fwrite(const map<M,T,less<M> >&): ";
          std::cerr <<"ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     typename std::map<M,T>::const_iterator  it;
     size_t      elements = stl_ctner.size();
     size_t         bytesM    = sizeof(M);
     size_t         bytesT    = sizeof(T);
     M              key;
     T              val;
     
     // writing the number of vector objects
     std::fwrite( (void*) &elements, sizeof(size_t), 1, fp );
     
     // writing all key-value pairs
     for ( it=stl_ctner.begin(); it!=stl_ctner.end(); it++ )
       {
          key = (*it).first;
          val = (*it).second;
          std::fwrite( (void*) &key, bytesM, 1, fp );
          std::fwrite( (void*) &val, bytesT, 1, fp );
       }
     return true;
 }
 
 
  

/**

Reads an STL map from a binary file. The map
is erased before adding the data if it already contains data. How many 
records are read is specified by an 'unsigned long' number which precedes
the dataset. 

@section arguments Input Arguments 

A pointer to a binary file opened in read binary mode ("rb"), and a 
reference to an STL map. 

@return skm_C_fread() returns the boolean 'true' if the number of data records 
which precedes the dataset in the file has been read correctly. If the 
file pointer is invalid or less data are read, the function returns false. 

@section implementation Implementation

Uses the ANSI C function fread(). The map which will hold the 
data is erased if it is not empty. Storage is reserved for the new number 
of elements which is read from file. Then the elements are read and the
container is returned.  

@section application Application

Efficiently read data into an STL map. 

@section messages Messages 

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this 
number does not match the number of records which were actually read. 
*/
template<class M, class T>
bool skm_C_fread( std::FILE* fp, std::map<M,T>& stl_ctner )
 {
     if ( !stl_ctner.empty() )
       stl_ctner.erase( stl_ctner.begin(), stl_ctner.end() );
       
     if ( fp == NULL )
       {
          std::cerr <<"\nbool skm_C_fread(map<M,T>&): ";
          std::cerr <<"ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     size_t         bytesM = sizeof(M),
                    bytesT = sizeof(T), 
                    counterM(0), counterT(0);
     size_t      elements(0), i;
     M              key;
     T              val;
     
     // 1. read number of record in the map and assert reading
     if ( !fread( (void*) &elements, sizeof(size_t), 1, fp ) ) {
          std::cout <<"\nbool skm_C_fread(map<M,T>&): ";
          std::cout <<"ERROR: could not read record length."<< std::endl;
          return false;
       }
     if ( elements > 0 ) {
          // 2. reading all map records
          for ( i=0; i<elements; i++ )
            {
               // reading key
               if ( std::fread( (void*) &key, bytesM, 1, fp ) ) counterM++;
               // reading value
               if ( std::fread( (void*) &val, bytesT, 1, fp ) ) counterT++;
               // storing value in map after key
               stl_ctner[ key ] = val;
            }
       }
     if ( counterM != elements || counterT != elements ) {
          std::cerr <<"\nbool skm_C_fread(vector<map<M,T>&): ";
          std::cerr <<"ERROR: incorrect number of map records were read: ";
          std::cerr <<"\nIndicated number: "<< elements <<", actual number read: "<< counterM << std::endl;
          return false;
       }
     return true;
 }
  
 
  /**
   
   Uses the C-style fwrite() function for binary file IO to write the contents
   of an STL unordered_map to a binary file. The data segment is preceded
   by an 'unsigned long' number which determines the number of objects
   which are written to file. Each record is preceded by the
   corresponding map key.
   
   @section arguments Input Arguments
   
   A file pointer of a binary file opened in write mode (e.g., "wb") and a
   constant reference to the STL map which will be written to
   file.
   
   @return a boolean indicating whether all data have
   been written correctly to the file.
   
   @section implementation Implementation
   
   Uses the ANSI standard C function fwrite(). NOTE that the map key must
   not contain any dynamically allocated data. Otherwise these data are
   sliced of the binary record.
   
   @section application Application
   
   To efficiently write map data to a binary file.
   
   @section messages Messages
   
   If the file pointer is invalid, method will quit, reporting an error.
   */
  template<class M, class T>
  bool skm_C_fwrite( std::FILE* fp, const std::unordered_map<M,T>& stl_ctner )
  {
    if ( fp == NULL )
    {
      std::cerr <<"\nbool skm_C_fwrite(const map<M,T,less<M> >&): ";
      std::cerr <<"ERROR: invalid file pointer."<< std::endl;
      return false;
    }
    typename std::unordered_map<M,T>::const_iterator  it;
    size_t      elements = stl_ctner.size();
    size_t         bytesM    = sizeof(M);
    size_t         bytesT    = sizeof(T);
    M              key;
    T              val;
    
    // writing the number of vector objects
    std::fwrite( (void*) &elements, sizeof(size_t), 1, fp );
    
    // writing all key-value pairs
    for ( it=stl_ctner.begin(); it!=stl_ctner.end(); it++ )
    {
      key = (*it).first;
      val = (*it).second;
      std::fwrite( (void*) &key, bytesM, 1, fp );
      std::fwrite( (void*) &val, bytesT, 1, fp );
    }
    return true;
  }
  
  
  
  
  /**
   
   Reads an STL unordered_map from a binary file. The map
   is erased before adding the data if it already contains data. How many
   records are read is specified by an 'unsigned long' number which precedes
   the dataset.
   
   @section arguments Input Arguments
   
   A pointer to a binary file opened in read binary mode ("rb"), and a
   reference to an STL map.
   
   @return skm_C_fread() returns the boolean 'true' if the number of data records
   which precedes the dataset in the file has been read correctly. If the
   file pointer is invalid or less data are read, the function returns false.
   
   @section implementation Implementation
   
   Uses the ANSI C function fread(). The map which will hold the
   data is erased if it is not empty. Storage is reserved for the new number
   of elements which is read from file. Then the elements are read and the
   container is returned.
   
   @section application Application
   
   Efficiently read data into an STL map.
   
   @section messages Messages
   
   The function will return false if (1) the file pointer is invalid, (2)
   the number of data records cannot be read correctly, and (3) if this
   number does not match the number of records which were actually read.
   */
  template<class M, class T>
  bool skm_C_fread( std::FILE* fp, std::unordered_map<M,T>& stl_ctner )
  {
    if ( !stl_ctner.empty() )
      stl_ctner.erase( stl_ctner.begin(), stl_ctner.end() );
    
    if ( fp == NULL )
    {
      std::cerr <<"\nbool skm_C_fread(map<M,T>&): ";
      std::cerr <<"ERROR: invalid file pointer."<< std::endl;
      return false;
    }
    size_t         bytesM = sizeof(M),
    bytesT = sizeof(T),
    counterM(0), counterT(0);
    size_t      elements(0), i;
    M              key;
    T              val;
    
    // 1. read number of record in the map and assert reading
    if ( !fread( (void*) &elements, sizeof(size_t), 1, fp ) ) {
      std::cout <<"\nbool skm_C_fread(map<M,T>&): ";
      std::cout <<"ERROR: could not read record length."<< std::endl;
      return false;
    }
    if ( elements > 0 ) {
      // 2. reading all map records
      for ( i=0; i<elements; i++ )
      {
        // reading key
        if ( std::fread( (void*) &key, bytesM, 1, fp ) ) counterM++;
        // reading value
        if ( std::fread( (void*) &val, bytesT, 1, fp ) ) counterT++;
        // storing value in map after key
        stl_ctner[ key ] = val;
      }
    }
    if ( counterM != elements || counterT != elements ) {
      std::cerr <<"\nbool skm_C_fread(vector<map<M,T>&): ";
      std::cerr <<"ERROR: incorrect number of map records were read: ";
      std::cerr <<"\nIndicated number: "<< elements <<", actual number read: "<< counterM << std::endl;
      return false;
    }
    return true;
  }
  
  


/**

Uses the C-style fwrite() function for binary file IO to write the contents
of an STL map of vector<T> to a binary file. The data segment is preceded 
by an 'size_t' number which determines the number of vector objects 
which are written to file. Each vector record is preceded by the 
corresponding map key.  

@section arguments Input Arguments 

A file pointer of a binary file opened in write mode (e.g., "wb") and a
constant reference to the STL map of vectors which will be written to 
file.  

@return a boolean indicating whether all map data have
been written correctly to the file.

@section implementation Implementation

Uses the ANSI standard C function fwrite(). NOTE that the map key must
not contain any dynamically allocated data. Otherwise these data are
sliced of the binary record. 

@section application Application

To efficiently write a map of vector data to a binary file. 

@section messages Messages 

If the file pointer is invalid, method will quit, reporting an error. 
*/
template<class M, class T>
bool skm_C_fwrite( std::FILE* fp, const std::map<M,std::vector<T> >& stl_ctner )
 {
     if ( fp == NULL ) {
          std::cerr <<"\nbool skm_C_fwrite(const map<M,vector<T>,less<M> >&): ";
          std::cerr <<"ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     typename std::map<M,std::vector<T> >::const_iterator  it;
     size_t      elements = stl_ctner.size();
     size_t         bytes    = sizeof(M);
     M              val;
     
     // writing the number of vector objects
     std::fwrite( (void*) &elements, sizeof(size_t), 1, fp );
     
     // writing all elements
     for ( it=stl_ctner.begin(); it!=stl_ctner.end(); it++ )
       {
          val = (*it).first;
          std::fwrite( (void*) &val, bytes, 1, fp );
          skm_C_fwrite( fp, (*it).second );
       }
     return true;
 }
 
 

/**

Reads a map of STL vectors and their keys from a binary file. The map
is erased before adding the data if it already contained any data. How many 
records are read is specified by an 'size_t' number which precedes
the dataset. 

@section arguments Input Arguments 

A pointer to a binary file opened in read binary mode ("rb"), and a 
reference to an STL map of vectors. 

@return the boolean 'true' if the number of data records
which precedes the dataset in the file has been read correctly. If the 
file pointer is invalid or less data are read, the function returns false. 

@section implementation Implementation

Uses the ANSI C function fread(). The map which will hold the 
data is erased if it is not empty. Storage is reserved for the new number 
of elements which is read from file. Then the elements are read and the
container is returned.  

@section application Application

Efficiently read data into an STL map of vectors. 

@section messages Messages 

The function will return false if (1) the file pointer is invalid, (2)
the number of data records cannot be read correctly, and (3) if this 
number does not match the number of records which were actually read. 
*/
template<class M, class T>
bool skm_C_fread( std::FILE* fp, std::map<M,std::vector<T> >& stl_ctner )
 {
     if ( !stl_ctner.empty() )
       stl_ctner.erase( stl_ctner.begin(), stl_ctner.end() );
       
     if ( fp == NULL ) {
          std::cerr <<"\nbool skm_C_fread(map<M,vector<T>,less<M> >&): ";
          std::cerr <<"ERROR: invalid file pointer."<< std::endl;
          return false;
       }
     size_t          bytes(sizeof(M)), counter(0);
     size_t       elements(0), i;
     std::vector<T>  val;
     M               key;
     
     // 1. read number of record in the map and assert reading
     if ( !std::fread( (void*) &elements, sizeof(size_t), 1, fp ) ) {
          std::cout <<"\nbool skm_C_fread(map<M,vector<T> >&): ";
          std::cout <<"ERROR: could not read record length."<< std::endl;
          return false;
       }
     if ( elements > 0 ) {
          // 2. reading all map records
          for ( i=0; i<elements; i++ )
            {
               // reading key
               fread( (void*) &key, bytes, 1, fp );
               // reading vector and counting the successfully read records
               if ( skm_C_fread( fp, val ) ) { 
                    counter++;
                    stl_ctner[ key ] = val;
                 }
            }
       }
     if ( counter != elements ) {
          std::cerr <<"\nbool skm_C_fread(vector<map<M,vector<T> >&): ";
          std::cerr <<"ERROR: incorrect number of map records were read: ";
          std::cerr <<"\nIndicated number: "<< elements <<", actual number read: "<< counter << std::endl;
          return false;
       }
     return true;
 }

/**
  @}
  */
  
 } // csmp


#endif







