//
//  vectorOperations.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 8/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef VECTOR_VARIABLE_OPERATIONS_H
#define VECTOR_VARIABLE_OPERATIONS_H

#include "VectorVariable.h"

namespace csmp {

template<size_t dim>
double valueAverage( const VectorVariable<dim>& );

} // end csmp

#endif /* VECTOR_VARIABLE_OPERATIONS_H */
