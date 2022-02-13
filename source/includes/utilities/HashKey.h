#ifndef CSMP_HASHKEY_H
#define CSMP_HASHKEY_H

#include "CSMP_definitions.h"

namespace csmp {

class HashKey {
    char  key_text[60];
    char  num[20];      // 20 = size of 8-byte size_t size_t as char string
    short hash_key_digits;
    std::set<uint32_t>            sorter;
    std::set<uint32_t>::iterator  sit;  
    size_t key;
  public:
    HashKey();
    ~HashKey() {};
    long Key( size_t a, size_t b );
    long Key( size_t a, size_t b, size_t c );
    long Key( size_t a, size_t b, size_t c, size_t d );
    
    void Key( size_t a, size_t b, char* s );
    void Key( size_t a, size_t b, size_t c, char* s );
    void Key( size_t a, size_t b, size_t c, size_t d, char* s );

    void Key( size_t a, size_t b, std::string& s );
    void Key( size_t a, size_t b, size_t c, std::string& s );
    void Key( size_t a, size_t b, size_t c, size_t d, std::string& s );
};



inline HashKey::HashKey()  
  : hash_key_digits(9)
 {
    if ( sizeof(size_t) ==  8 ) hash_key_digits = 19;
    if ( sizeof(size_t) == 12 ) hash_key_digits = 29;
 }


/** Key(a,b)

    Makes unique hashkey out of 2 unsigned integers. The integers
    are sorted in increasing order.

    @note use only unsigned integers here, because negative numbers cannot be concenated.
    */
inline long HashKey::Key( size_t a, size_t b )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    std::strcpy( key_text, "\0");
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         std::strcat( key_text, num );
      }
    if ( std::strlen(key_text) > static_cast<uint32_t>(hash_key_digits) )
      {
          std::cout <<"\nHashKey::Key: too many digits to convert string '";
          std::cout << key_text << "' to size_t."<< std::endl;
      }  
#ifdef SGI       
    return std::atoll( key_text );
#else    
    return std::atol( key_text );
#endif 
 }




inline long HashKey::Key( size_t a, size_t b, size_t c )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    sorter.insert(c);
    std::strcpy( key_text, "\0");
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         std::strcat( key_text, num );
      }
    if ( std::strlen(key_text) > static_cast<uint32_t>(hash_key_digits) )
      {
          std::cout <<"\nHashKey::Key: too many digits to convert string '";
          std::cout << key_text << "' to size_t."<< std::endl;
      } 
#ifdef SGI       
    return std::atoll( key_text );
#else    
    return std::atol( key_text );
#endif 
 }



inline long HashKey::Key( size_t a, size_t b, size_t c, size_t d )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    sorter.insert(c);
    sorter.insert(d);
    std::strcpy( key_text, "\0");
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         std::strcat( key_text, num );
      }
    if ( std::strlen(key_text) > static_cast<uint32_t>(hash_key_digits) )
      {
          std::cout <<"\nHashKey::Key: too many digits to convert string '";
          std::cout << key_text << "' to size_t."<< std::endl;
      }  
#ifdef SGI       
    return std::atoll( key_text );
#else    
    return std::atol( key_text );
#endif 
 }


// outputting strings


inline void HashKey::Key( size_t a, size_t b, char* s )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    std::strcpy( s, "\0");
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         std::strcat( s, num );
      }
 }


inline void HashKey::Key( size_t a, size_t b, size_t c, char* s )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    sorter.insert(c);
    std::strcpy( s, "\0");
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         std::strcat( s, num );
      }
 }



inline void HashKey::Key( size_t a, size_t b, size_t c, size_t d, char* s )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    sorter.insert(c);
    sorter.insert(d);
    std::strcpy( s, "\0");
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         std::strcat( s, num );
      }
 }


// string

inline void HashKey::Key( size_t a, size_t b, std::string& s )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    s = "\0";
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         s += num;
      }
 }


inline void HashKey::Key( size_t a, size_t b, size_t c, std::string& s )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    sorter.insert(c);
    s = "\0";
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         s += num;
      }
//    cout <<"\nHashKey::Key: "<< s.CharPointer() << endl;
 }



inline void HashKey::Key( size_t a, size_t b, size_t c, size_t d, std::string& s )
 {
    sorter.erase( sorter.begin(), sorter.end() );
    sorter.insert(a);
    sorter.insert(b);
    sorter.insert(c);
    sorter.insert(d);
    s = "\0";
    
    for ( sit=sorter.begin(); sit!=sorter.end(); sit++ )
      {
         std::sprintf( num, "%lu", *sit );
         s += num;
      }
 }

} // csp

#endif


