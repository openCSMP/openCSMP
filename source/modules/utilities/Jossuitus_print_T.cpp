// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <iostream>

/**
- prints optional C-string optcstr followed by
- all elements of the collection coll
- separated by spaces
 */
template <class T>
void PRINT_ELEMENTS (const T& coll, const char* optcstr="")
{
    std::cout << optcstr;
    for ( typename T::const_iterator 
          pos=coll.begin(); pos!=coll.end(); ++pos ) {
        std::cout << *pos << ' ';
    }
    std::cout << std::endl;
}
