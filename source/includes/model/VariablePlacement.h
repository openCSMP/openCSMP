#ifndef CSMP_VARIABLE_PLACEMENT_H_INCLUDED
#define CSMP_VARIABLE_PLACEMENT_H_INCLUDED

#include "VariableTypeTraits.h"

namespace csmp {

  template<size_t dim> class Element;
  template<size_t dim> class Node;
  template<size_t dim, PLACEMENT pl> class FiniteElementPlacement;
  template<size_t dim, PLACEMENT pl> struct FiniteElementPlacementCollection;
  template<size_t dim, PLACEMENT pl> class FiniteVolumePlacement;
  template<size_t dim, PLACEMENT pl> struct FiniteVolumePlacementCollection;
  template<size_t dim> struct NeighbourNodeCollection;

  template<size_t dim>
  void calculateN(const Element<dim>& e, const Point<dim>& p, double* coeff);

  template<size_t dim>
  void calculateDN(const Element<dim>& e, const Point<dim>& p, std::vector<double>* coeff);

  template<size_t dim>
  Point<dim> directedAreaOfFacet(const Element<dim>& e, size_t iFacet);
}

#endif
