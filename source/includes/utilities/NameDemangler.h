//
//  NameDemangler.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 9/10/19.
//  Copyright © 2019 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_NAME_DEMANGLER_H
#define CSMP_NAME_DEMANGLER_H

#include <string>
#include <typeinfo>

std::string demangle( const char* name );

template <class T>
std::string type(const T& t) {

    return demangle(typeid(t).name());
}

#endif /* CSMP_NAME_DEMANGLER_H */
