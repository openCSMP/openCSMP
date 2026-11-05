// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Experimental_Example.h
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//

#ifndef EXPERIMENTAL_EXAMPLE_H
#define EXPERIMENTAL_EXAMPLE_H

#include "Example.h"

namespace csmp {

template<uint32_t> class Model;


class  Experimental_Example final : public Example {
public:
  void Run() override;
  void Specifications() override;
};

} // csmp

#endif // EXPERIMENTAL_EXAMPLE_H
