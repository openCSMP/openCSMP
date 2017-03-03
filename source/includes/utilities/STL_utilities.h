#ifndef CSMP_STL_UTILITIES_H
#define CSMP_STL_UTILITIES_H

#include "CSMP_definitions.h"

namespace csmp {

void printVector( const std::vector<size_t>& array );
void printVector( const std::vector<double>& array );
void printVector( const std::vector<std::vector<double64> >& array );
void printVector( const char* headline, const std::vector<std::vector<double64> >& array );

template<typename T1, typename T2>
void printMap( const std::map<T1,T2>& );



template<typename T1, typename T2>
inline void printMap( const std::map<T1,T2>& m )
 {
   for ( typename std::map<T1,T2>::const_iterator
         bit=m.begin(); bit!=m.end(); bit++ ) {
      std::cout << "Key: "<< (*bit).first<<" Value: " << (*bit).second<<std::endl;
   }

 } // end printMap

} // end namespace csmp


#endif
