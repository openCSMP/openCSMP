#ifndef CSMP_GENERIC_TRANSPORT_SCHEME_H
#define CSMP_GENERIC_TRANSPORT_SCHEME_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"

namespace csmp {

  // Forward declares
  template<size_t dim> class Model;
  template<size_t dim> class Region;
  template<size_t dim> class Element;
  template<size_t dim> class Node;

  template<size_t dim> class MatrixOperator;
  template<size_t dim> class VectorOperator;
  template<size_t dim> class LinearSystemAccumulator;

  /// What mode a system accumulation performs in
  enum ACCUMULATION_MODE {
    ADD_ACCUMULATE,
    MULTIPLY_ACCUMULATE,
    ADD_LATER
  };


} // end csmp


#endif /* CSMP_GENERIC_TRANSPORT_SCHEME_H */
