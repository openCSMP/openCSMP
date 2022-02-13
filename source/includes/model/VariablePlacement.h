#ifndef CSMP_VARIABLE_PLACEMENT_H_INCLUDED
#define CSMP_VARIABLE_PLACEMENT_H_INCLUDED

#include "VariableTypeTraits.h"

namespace csmp {

  template<uint32_t dim> class Element;
  template<uint32_t dim> class Node;
  template<uint32_t dim, PLACEMENT pl> class FiniteElementPlacement;
  template<uint32_t dim, PLACEMENT pl> struct FiniteElementPlacementCollection;
  template<uint32_t dim, PLACEMENT pl> class FiniteVolumePlacement;
  template<uint32_t dim, PLACEMENT pl> struct FiniteVolumePlacementCollection;
  template<uint32_t dim> struct NeighbourNodeCollection;

  template<uint32_t dim>
  void calculateN(const Element<dim>& e, const Point<dim>& p, double* coeff);

  template<uint32_t dim>
  void calculateDN(const Element<dim>& e, const Point<dim>& p, std::vector<double>* coeff);

  template<uint32_t dim>
  Point<dim> directedAreaOfFacet(const Element<dim>& e, uint32_t iFacet);
}

#endif
