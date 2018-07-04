//
//  convertColorToPermeability.h
//  CPP11_Tests
//
//  Created by Stephan Matthai on 2/10/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#ifndef CONVERT_COLOR_TO_PERMEABILITY_H
#define CONVERT_COLOR_TO_PERMEABILITY_H

#include "CSMP_number_types.h"

namespace csmp {

  class Matrix;

  /// converts 256-color bitmap values into permeability, see documentation for values
  void convertColorToPermeability( double64, Matrix& );
  
 } // end

#endif
