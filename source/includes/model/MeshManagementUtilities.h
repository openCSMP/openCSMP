//
//  MeshManagementUtilities.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 3/7/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_MESH_MANAGEMENT_UTILITIES_H
#define CSMP_MESH_MANAGEMENT_UTILITIES_H

#include "CSMP_definitions.h"
#include "MeshPatchAttributes.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class MeshManager;

/// counts and returns current indices of elements that may give rise to problems during the assignment of boundary conditions
template<size_t dim>
size_t detectElementsWithAllNodesOnBoundary( const MeshManager<dim>&, std::set<size_t>& );

// TODO: implement
template<size_t dim,template<size_t> class CELL>
size_t detectDisconnectedCells( const MeshManager<dim>&, std::set<size_t>& );

/// finds the connected (contiguous) mesh patches in the supplied range of cells storing them in map with names that reflect their dimensionality and cell numbers
template<size_t dim, template<size_t> class CELL>
size_t  findStandAloneMeshPatches( typename std::deque<CELL<dim>*>::const_iterator begin,
                                   typename std::deque<CELL<dim>*>::const_iterator end,
                                   std::map<std::string,std::deque<CELL<dim>*> >& );

/// finds pointers to all contiguous regions in a mesh (collection of mesh patches) returning pointers to them so that they can be explored; returns patches found
template<size_t dim, template<size_t> class CELL>
size_t  findPointersToStandAloneMeshPatches( typename std::deque<CELL<dim>*>::const_iterator begin,
                                             typename std::deque<CELL<dim>*>::const_iterator end,
                                             std::map<Element<dim>*,MeshPatchAttributes>& );

/// using a breadth-first mesh traversal, finds all the Element, Face, or InterFace objects that belong to this contiguous mesh patch
// TODO: this method might also find mesh patches with different dimensional elements; test this and change if necessary
template<size_t dim,template<size_t> class CELL>          // //
void findContiguousMeshPatch( CELL<dim>* const cell_pointer, std::set<CELL<dim>*>& contiguous_subset_of_cells );

/// breadth-first mesh traversal starting at a Node; returns number of discovered nodes
template<size_t dim>
size_t findContiguousMeshPatch( csmp::Node<dim>* const node_pointer, std::deque<Element<dim>*>& elements );

/// relying on the parent element information from its nodes, method finds higher-dim neighbors of each element face and connects itself with them and vice versa; returns # found
template<size_t dim>
size_t connectNeighborsUsingNodeParents( Element<dim>* const );

// TODO: not sure how to do this in a generic way
template<size_t dim>
void updateParentElementConnectivity( Node<dim>* const );


/// TODO: Using the parent elements of its nodes, finds its higher-dimensional neighbor on inside or outside
template<size_t dim>
Element<dim>* const findInnerHigherDimensionalNeighborFromNodes( Element<dim>* const, INTERFACE_SIDE );

/// Connects nodes to Face, finding them by matching the faces of the supplied higher dimensional elements
template<size_t dim>
void findNodesViaHigherDimensionalNeighbors( const Element<dim>* const inner_nbor,
                                             const Element<dim>* const outer_neighbor,
                                             Face<dim>* const );

/// Connects nodes to InterFace, finding them by matching the faces of the supplied higher dimensional elements; face IDs are set as well
template<size_t dim>
void findNodesViaHigherDimensionalNeighbors( const Element<dim>* const inner_nbor,
                                             const Element<dim>* const outer_neighbor,
                                             InterFace<dim>* const );


/// Surt's method to efficiently erase vector Element from a pointer vector.
template<size_t dim>
void eraseElementPointerFromVector( std::vector<csmp::Element<dim>*>&, const Element<dim>* );

} // end csmp

#endif /* CSMP_MESH_MANAGEMENT_UTILITIES_H */
