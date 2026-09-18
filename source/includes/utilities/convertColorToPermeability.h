// Copyright (c) 2014 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONVERT_COLOR_TO_PERMEABILITY_H
#define CONVERT_COLOR_TO_PERMEABILITY_H

#include <iostream>

namespace csmp {

  class Matrix;

  /// converts 256-color bitmap values into permeability, see documentation for values
  void convertColorToPermeability( double, Matrix& );
  
 } // end

#endif
