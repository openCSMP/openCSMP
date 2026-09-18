// Copyright © 2019 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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
