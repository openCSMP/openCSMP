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
#include "plf_colony.h"

namespace csmp {

template<uint32_t> class Point;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;

template<uint32_t> class Model;
template<uint32_t> class MeshManager;
template<uint32_t> class MeshPatch;
template<uint32_t> class Region;
template<uint32_t> class Boundary;
template<uint32_t> class SplitBoundary;


// MESH CONNECTIVITY

/// retrieves and returns the first contiguous element patch that can be reached by mesh traversal from the starting element
template<uint32_t dim, template<uint32_t> class CELL>
void floodFill( CELL<dim>* const eptr, std::set<CELL<dim>*>& output_contiguous_subset );

/// checks all elements of the surface region for whether their neighbor elements have normals that deviate less than 90o from their normals
template<uint32_t dim>
bool checkNeighborNormalsForConsistentOrientation( const Region<dim>& );

/// using a breadth-first mesh traversal, finds all the Element, Face, or InterFace objects that belong to this contiguous mesh patch
template<uint32_t dim,template<uint32_t> class CELL>
size_t findContiguousMeshPatch( CELL<dim>* const entry_cell, std::set<CELL<dim>*>& contiguous_subset_of_cells );

/// finds the connected (contiguous) mesh patches in the supplied range of cells storing them in map with names that reflect their dimensionality and cell numbers
template<uint32_t dim, template<uint32_t> class CELL>
size_t  findStandAloneMeshPatches( typename plf::colony<CELL<dim>>::iterator begin,
                                   typename plf::colony<CELL<dim>>::iterator end,
                                   std::map<std::string,std::vector<CELL<dim>*> >& );

/// finds pointers to all contiguous regions in a mesh (collection of mesh patches) returning pointers to them so that they can be explored; returns patches found
template<uint32_t dim, template<uint32_t> class CELL>
size_t  findPointersToStandAloneMeshPatches( typename std::vector<CELL<dim>*>::const_iterator begin,
                                             typename std::vector<CELL<dim>*>::const_iterator end,
                                             std::map<CELL<dim>*,MeshPatch<dim>>& );


// DIAGNOSTICS

/// calculates the number of model cells that fall into the cell category indicated by placement; returns total number of cells in the model
template<uint32_t dim>
size_t currentCellTypes( const MeshManager<dim>&, PLACEMENT, size_t& volume_cells, size_t& surface_cells, size_t& line_cells );

/// determines whether mesh in model is built from finite elements with a local coordinate system
template<uint32_t dim>
bool isoparametricElementMesh( const Model<dim>& );

/// counts and returns current indices of elements that may give rise to problems during the assignment of boundary conditions
template<uint32_t dim>
size_t detectElementsWithAllNodesOnBoundary( const MeshManager<dim>&, std::set<size_t>& );

/// finds cells that have the same nodes and reports their numbers; verbose reports the duplicate cells
template<uint32_t dim, template<uint32_t> class CELL>
size_t detectDuplicateCells( typename std::vector<CELL<dim>*>::const_iterator begin,
                             typename std::vector<CELL<dim>*>::const_iterator end,
                             bool verbose );

/// Computes parent element barycentre-to-node distances for range of nodes;  returns them into vector [e1,e2...e_n,e_sum] with a length of parent elements+1
template<uint32_t dim>
void distancesAndWeights( typename std::vector<Node<dim>*>::const_iterator nodes_begin,
                          typename std::vector<Node<dim>*>::const_iterator nodes_end,
                          std::vector<std::vector<double> >& distances_and_weight );

/// container of element pointers and local face ids of elements contacting each other across a split boundary
typedef std::pair<std::pair<Element<3U>*, uint32_t>, std::pair<Element<3U>*, uint32_t> > OppositeElements;

/// find all elements in a model that contact eachother across split interfaces and are node-matched
template<uint32_t dim>
bool findSplitInterfaceElements( const Region<dim>&,
                                 std::vector<std::pair<std::pair<Element<dim>*, uint32_t>,
                                                       std::pair<Element<dim>*, uint32_t> > >& opposite_elmts_and_face_ids );
                                                       
/// finds node-matched internal split model boundaries, a lower prop value distinguishes the inside; if both values are the same, no boundary is inserted; returns number of interfaces recorded
template<uint32_t dim>
size_t findSplitInterfaceElements( const Model<dim>& model, const std::string& property_to_distinguish_regions,
                                   std::vector<std::pair<std::pair<Element<dim>*, uint32_t>,
                                                         std::pair<Element<dim>*, uint32_t> > >& opposite_elmts_and_face_ids );

/// for supplied edge nodes, find their volumetric parent elements; if find segment ids is on, their local numbers are assigned to Idx of the parent elements
size_t parentElementsSharingMultipleEdgeNodes( const std::vector<Node<3U>*>&  edge_nodes,
                                               std::map<Element<3>*,std::vector<Node<3>*> >& segm_parents,
                                               bool find_segment_ids );

/// finds node by point coordinate; returns -1 if not found; @attention tolerance needs to account for single-precision of CAD tools
template<uint32_t dim>
long  findNode( const Model<dim>&, const Point<dim>& pxyz, double tolerance, bool verbose = false );

/// find node by its position as identified from its coordinates: tolerance should take into account single-precision of CAD tools
size_t  findNode( const Model<1U>&, double nx, double tolerance );
/// 2D version
size_t  findNode( const Model<2U>&, double nx, double ny, double tolerance );
/// 3D version
size_t  findNode( const Model<3U>&, double nx, double ny, double nz, double tolerance );

/// prints sorted global element node numbers in a compact way
template<uint32_t dim, template<uint32_t> class CELL>
void printNodes( const CELL<dim>& );


// MESHING UTILITIES


// MESHING UTILITIES FOR INDIVIDUAL ELEMENTS/FACES/INTERFACES

// TODO: implement
//template<uint32_t dim,template<uint32_t> class CELL>
//size_t detectDisconnectedCells( const MeshManager<dim>&, std::set<uint32_t>& );

/// loops over the valid neighbors of the cell and sets their neighbor pointers to point to this cell to nullptr
template<uint32_t dim, template<uint32_t> class CELL>
void detachNeighborsFrom( CELL<dim>* const cell_to_detach_neighbors_from );

/// returns angle (in degrees) between the normals of the two cells, which must be surfaces (only in 3D)
template<template<uint32_t> class CELL>
double angleBetweenSurfaceCells( const CELL<3>* const cell1, const CELL<3>* const cell2 );

/// Line elements can exist in all 3 spatial dimensions.
template<uint32_t dim, template<uint32_t> class CELL>
double angleBetweenLineCells( const CELL<dim>* const cell1, const CELL<dim>* const cell2 );

/// traverses mesh via node neighbors and collects nodes into argument set; @return number of discovered nodes; requires node to parent connectivity
template<uint32_t dim>
size_t findInterconnectedNodeCluster( Node<dim>* const, std::set<Node<dim>*>& contiguous_set_of_nodes );

/// relying on the parent element information from its nodes, method finds higher-dim neighbors of each element face and connects itself with them and vice versa; returns # found
template<uint32_t dim>
size_t connectNeighborsUsingNodeParents( Element<dim>* const );

// TODO: implement: not sure how to do this in a generic way
template<uint32_t dim>
void updateParentElementConnectivity( Node<dim>* const );

/// Using the parent elements of its nodes, finds its higher-dimensional neighbor on inside or outside
template<uint32_t dim>
Element<dim>* const findInnerHigherDimensionalNeighborFromNodes( Element<dim>* const, INTERFACE_SIDE );

/// Connects nodes to Face, finding them by matching the faces of the supplied higher dimensional elements
template<uint32_t dim>
void findNodesViaHigherDimensionalNeighbors( Element<dim>* const inner_nbor,
                                             Element<dim>* const outer_neighbor,
                                             Face<dim>* const );

/// Connects nodes to InterFace, finding them by matching the faces of the supplied higher dimensional elements; face IDs are set as well
template<uint32_t dim>
void findNodesViaHigherDimensionalNeighbors( Element<dim>* const inner_nbor,
                                             Element<dim>* const outer_neighbor,
                                             InterFace<dim>* const );


/// Surt's method to efficiently erase vector Element from a pointer vector.
template<uint32_t dim>
void eraseElementPointerFromVector( std::vector<csmp::Element<dim>*>&, const Element<dim>* );

/// by comparison of node locations, finds overlapping cells and reports them
bool findCollocatedCells(); // TODO: not implemented yet

/// assuming that the elements are adjacent, method finds their faces that are in contact with one another from their shared nodes (faster)
template<uint32_t dim>
std::pair<size_t,size_t> findAdjacentFacesFromNeighbors( Element<dim>* const eptr1, Element<dim>* const eptr2 );

/// assuming that the elements are adjacent, method finds their faces that are in contact with one another from their shared nodes (slower)
template<uint32_t dim>
std::pair<size_t,size_t> findAdjacentElementFaces( Element<dim>* const eptr1, Element<dim>* const eptr2 );

/// returns true if the elements contain each others barycentre
template<uint32_t dim>
bool interPenetrating( const Element<dim>* const, const Element<dim>* const );

/// Tests whether a tetrahedron is degenerate because all of its vertices lie within a single plane; tolerance in meters.
bool hasNonManifoldVertices( const csmp::Element<3U>* const tptr, double tolerance=1.0e-5 );

/// detects whether the point is contained in any of the elements of the region
csmp::Element<3u>* const pointInVolumeElement( Region<3U>& region, const Point<3U>& query );

/// detect degenerate elements by using the node coordinates to check whether some nodes have the same location
template<uint32_t dim>
size_t collocatedNodes( const Element<dim>* const );


// UTILITIES FOR TESTING ETC

///  captures a snapshot of the current cell connectivity for the range of cells; developed for testing
template<uint32_t dim, template<uint32_t> class CELL>
void backupNeighborConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                 typename std::vector<CELL<dim>*>::const_iterator last,
                                 std::vector<std::vector<CELL<dim>*> >& nbor_pointers );

/// tests whether all the expected cell functionality is there and operational
template<uint32_t dim, template<uint32_t> class CELL>
bool integrityCheck( typename plf::colony<CELL<dim>>::const_iterator first,
                     typename plf::colony<CELL<dim>>::const_iterator last );

template<uint32_t dim, template<uint32_t> class CELL>
bool integrityCheck( const plf::colony<CELL<dim> >&,
                     typename std::vector<CELL<dim>*>::const_iterator first,
                     typename std::vector<CELL<dim>*>::const_iterator last );
                     
/// check whether the right types of elements have been connected with one another; counts violations also considering elements without minimum number of neighbors
template<uint32_t dim>
size_t connectivityCheck( typename std::vector<Element<dim>*>::const_iterator first,
                          typename std::vector<Element<dim>*>::const_iterator last );

/// finds the min max corners of the bounding box for the supplied range of nodes
template<uint32_t dim>
std::pair<Point<dim>,Point<dim>>  boundingBox( typename std::vector<Node<dim>*>::const_iterator first,
                                               typename std::vector<Node<dim>*>::const_iterator last );

/// pretty prints line elements as a chain from beginning to end; returns number of elements printed
template<uint32_t dim>
size_t  printLineElementRegion( const Model<dim>&, const char* region_name, bool renumber_nodes );

/// prints the idx and coordinates of the supplied nodes in a format that can be pasted into a spreadsheet
template<uint32_t dim>
void printNodeCoordinates( typename std::vector<Node<dim>*>::const_iterator first,
                           typename std::vector<Node<dim>*>::const_iterator last );


} // end csmp

#endif /* CSMP_MESH_MANAGEMENT_UTILITIES_H */
