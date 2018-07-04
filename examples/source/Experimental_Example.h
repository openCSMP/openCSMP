//
//  Experimental_Example.h
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#ifndef EXPERIMENTAL_EXAMPLE_H
#define EXPERIMENTAL_EXAMPLE_H

#include "Example.h"

namespace csmp {

class  Experimental_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};

} // csmp

#endif // EXPERIMENTAL_EXAMPLE_H
