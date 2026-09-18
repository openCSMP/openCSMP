// Copyright © 2021 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_MESH_PATCH_H
#define CSMP_MESH_PATCH_H

#include "FiniteElement.h"
#include "plf_colony.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Node;

/**
      Contiguous mesh patch that can, but does not have to be integrated into a Model.
      The idea is to temporarily store parts of a mesh during rebuilding operations;
      Size info is used also to keep number information uptodate.
          
         @author SKM
         @date 23/5/21
*/
template<uint32_t dim>
class MeshPatch {
  public:
    /// specify whether this is a volume, surface or line element patch
    explicit MeshPatch( CELL_SHAPE cs ) : cell_dimension_{cs} {}
    
/// creates an interconnected element patch from the shared faces of the supplied elements; replicates nodes upon request, and returns perimeter nodes and size of patch
    size_t BuildInterveningPatch( const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,
                                                              std::pair<Element<dim>*,uint32_t> > >&,
                                  bool copy_nodes,
                                  std::vector<Node<dim>*>  perimeter_node_ptrs );
                                  
    /// how many different finite element types the patch contains
    size_t CellTypes() const;
    /// whether this is a LINE, SURFACE or VOLUME mesh
    CELL_SHAPE CellGeometry() const;
    size_t Cells() const;
    size_t Nodes() const;

    typename plf::colony<Element<dim>>::iterator ElementsBegin();
    typename plf::colony<Element<dim>>::iterator ElementsEnd();

    typename plf::colony<Node<dim>>::iterator NodesBegin();
    typename plf::colony<Node<dim>>::iterator NodesEnd();
    
  private:
    std::map<CSMP_FEM_TYPE,FiniteElement*>  fe_ptrs_;        ///< finite-element pointers used in the FE policy
    plf::colony<Node<dim>>                  nodes_;          ///<  nodes (if they are pre-existing, a light copy without properties is made)
    plf::colony<Element<dim>>               elements_;       ///<  elements
    CELL_SHAPE                              cell_dimension_; ///<  types of elements in the patch
};

} // end csmp

#endif /* CSMP_MESH_PATCH_ATTRIBUTES_H */
