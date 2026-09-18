// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_STL_UTILITIES_H
#define CSMP_STL_UTILITIES_H

#include "CSMP_definitions.h"

namespace csmp {

void printVector( const std::vector<uint32_t>& array );
void printVector( const std::vector<double>& array );
void printVector( const std::vector<std::vector<double> >& array );
void printVector( const char* headline, const std::vector<std::vector<double> >& array );

template<typename T1, typename T2>
void printMap( const std::map<T1,T2>& );


/// reads vector<vector> from filestream where the elements of the vector are sequential (more methods in readWriteBinary.h)
template<typename T>
void readVectorOfVectors( std::ifstream&, size_t total_items, size_t entries_per_vector, std::deque<std::vector<T> >& file_records );




template<typename T1, typename T2>
inline void printMap( const std::map<T1,T2>& m )
 {
   for ( typename std::map<T1,T2>::const_iterator
         bit=m.begin(); bit!=m.end(); bit++ ) {
      std::cout << "Key: "<< (*bit).first<<" Value: " << (*bit).second<<std::endl;
   }

 } // end printMap


/**
  sort() followed by unique(), erasing any non-unique items.
 */
template<typename Container>
inline void sortAndUnique(Container& container)
{
  std::sort(container.begin(), container.end());
  auto newend = std::unique(container.begin(), container.end());
  container.erase(newend, container.end());
}

/**
  Ensure that a container (typically a std::vector) has enough space,
  geometrically resizing if needed.
 */
template<typename Container>
inline void ensureCapacity(Container& container, size_t spaceNeeded = 1)
{
  size_t oldSize = container.size();
  size_t oldCapacity = container.capacity();
  size_t newCapacity = oldSize + spaceNeeded;
  if (newCapacity < oldCapacity) {
    newCapacity = std::max<size_t>(newCapacity, (oldCapacity * 4 + 2) / 3);
    newCapacity = std::max<size_t>(newCapacity, 4);
    container.reserve(newCapacity);
  }
}

} // end namespace csmp


#endif
