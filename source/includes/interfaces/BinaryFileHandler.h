// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_BINARY_FILE_HANDLER_H
#define CSMP_BINARY_FILE_HANDLER_H

#include "CSMP_definitions.h"

namespace csmp {

/// writes objects and collections thereof to binary file
template<class S,class T>
class BinaryFileHandler {
    char  fname[NAME_STRING];
    std::FILE* fp;
  public:
    BinaryFileHandler();
    ~BinaryFileHandler();
     void SaveAs( const char* s, const std::vector<S>& obj_list );
     void SaveAs( const char* s, const std::map<S,T>& obj_list );
     bool Read( const char* s, std::vector<S>& obj_list );
     bool Read( const char* s, std::map<S,T>& obj_list );
};



template<class S,class T>
BinaryFileHandler<S,T>::BinaryFileHandler()
 {
 }


template<class S,class T>
BinaryFileHandler<S,T>::~BinaryFileHandler()
 {
 }


template<class S,class T>
void BinaryFileHandler<S,T>::SaveAs( const char* s, const std::vector<S>& obj_list )
 {
    unsigned long record_size;
    typename std::vector<S>::const_iterator iter;

    record_size = sizeof(S);

    // "wb" overwrite mode instead of "ab" for appending to file
    if ( (fp = std::fopen ( s,"wb")) == NULL ) 
     {
        std::cout <<"\nBinaryFileHandler::SaveAs: Object file '" << s;
        std::cout <<"' not accessible !";
        return;
     }
    /* Writing the contents of the class to binary datafile */
    for ( iter=obj_list.begin(); iter!=obj_list.end(); iter++ )
      {
        if ( (fwrite( &(*iter), record_size, 1, fp )) != 1 )
          std::cout <<"\nBinaryFileHandler::SaveAs: Unable to write data to datafile";
      }
    fclose( fp );
    
 } // end SaveAs(vector)




template<class S,class T>
void BinaryFileHandler<S,T>::SaveAs( const char* s, const std::map<S,T>& obj_list )
 {
    unsigned long record_size;
    typename std::map<S,T>::const_iterator iter;

    record_size = sizeof(T);

    // "wb" overwrite mode instead of "ab" for appending to file
    if ( (fp = std::fopen ( s,"wb")) == NULL ) 
     {
        std::cout <<"\nBinaryFileHandler::SaveAs: Object file '" << s;
        std::cout <<"' not accessible !";
        return;
     }
    /* Writing the contents of the class to binary datafile */
    for ( iter=obj_list.begin(); iter!=obj_list.end(); iter++ )
      {
        if ( (fwrite( &((*iter).second), record_size, 1, fp )) != 1 )
          std::cout <<"\nBinaryFileHandler::SaveAs: Unable to write data to datafile";
      }
    fclose( fp );
    
 } // end SaveAs(map)





/// Objects of type S are added in the vector<double> 'obj_list'
template<class S,class T>
bool BinaryFileHandler<S,T>::Read( const char* s, std::vector<S>& obj_list )
 {
    unsigned long record_size;
    S             object;

    record_size = sizeof(S);

    if ( (fp = std::fopen ( s,"rb")) == NULL ) 
     {
        std::cout <<"\nBinaryFileHandler::Read: Object file '" << s;
        std::cout <<"' not accessible !";
        return( false );
     }
     
    obj_list.erase( obj_list.begin(), obj_list.end() ); 
    /* Writing the contents of the class to binary datafile */
    while (!feof(fp))
       {
         if ( fread( &object, record_size, 1, fp ) == 1 )
           {
              obj_list.push_back( object );
           }
       }
    fclose( fp );
    return( true );
      
 } // end Read(vector)




/// Objects of type T are read and added in the map 'obj_list' with keys S
template<class S,class T>
bool BinaryFileHandler<S,T>::Read( const char* s, std::map<S,T>& obj_list )
 {
    unsigned long            record_size;
    typename std::map<S,T>::iterator  iter;
    T                        object;
    S                        key;

    record_size = sizeof(T);

    if ( (fp = std::fopen ( s,"rb")) == NULL ) 
     {
        std::cout <<"\nBinaryFileHandler::Read: Object file '" << s;
        std::cout <<"' not accessible !";
        return( false );
     }
    /* Writing the contents of the class to binary datafile */
    while (!feof(fp))
       {
         if ( fread( &object, record_size, 1, fp ) == 1 )
           {
              key           = object.name;
              obj_list[key] = object;
           }
       }
    fclose( fp );
    return( true );
      
 } // end Read(map)



} // end namespace csmp

#endif
