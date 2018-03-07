#ifndef CSMP_VARIABLE_TYPE_TRAITS_H
#define CSMP_VARIABLE_TYPE_TRAITS_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include <type_traits>


namespace csmp {
  template<size_t dim, VARIABLE_TYPE vt>
  struct VariableTypeTraits
  {
  };

  template<size_t dim>
  struct VariableTypeTraits<dim,SCALAR>
  {
    typedef double64 ReturnType;
    typedef ScalarVariable VariableType;
  };

  template<size_t dim>
  struct VariableTypeTraits<dim,VECTOR>
  {
    typedef Point<dim> ReturnType;
    typedef VectorVariable<dim> VariableType;
  };

  template<size_t dim>
  struct VariableTypeTraits<dim,TENSOR>
  {
    typedef TensorVariable<dim> VariableType;
  };

  template<size_t dim>
  struct VariableTypeTraits<dim,ARRAY>
  {
    typedef std::vector<double64> ReturnType;
    typedef ArrayVariable VariableType;
  };

  template<size_t dim>
  struct VariableTypeTraits<dim,FLAGGEDARRAY>
  {
    typedef FlaggedArrayVariable VariableType;
  };
}

#endif
