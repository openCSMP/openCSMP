//
//  MeshManagementUtilities.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 3/7/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "MeshManagementUtilities.h"
#include "CSMP_highLevelUtilities.h" // perhaps the mesh related stuff should be moved here
#include "MeshManager.h"
#include "Element.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "Box.h"

using namespace std;

namespace csmp {

/**
Counts elements the nodes of which are all located on the model boundary.

@attention such elements typically give rise to problems with the assignment of
boundary conditions and should be eliminated.

@attention method will work only if the elements are stored in elmt_connector deque.

@author SKM
*/
template<size_t dim>
size_t detectElementsWithAllNodesOnBoundary( MeshManager<dim>& mmgr, set<size_t>& belmts )
{
  belmts.clear();

  assert( mmgr.Elements() > 0 );
  if ( mmgr.Elements() == 0 ) return 0U;

  // traversal of the existing mesh nodes to find all its elements
  size_t boundary_only_elements( 0U );
  for ( typename set<Element<dim>*>::const_iterator
       it=mmgr.ElementsBegin(); it!=mmgr.ElementsEnd(); ++it ) {
      const size_t nodes((*it)->Nodes());
      size_t       counter(0U);
      for ( size_t i=0U; i<nodes; ++i )
        if ( (*it)->N(i)->AtBoundary() != NOT ) counter++;
      if ( counter == nodes ) {
            belmts.insert( (*it)->Idx() );
            boundary_only_elements++;
        }
   }
  	
  return boundary_only_elements;
}




/**

Breadth first traversal of mesh that can contain elements of any dimension.
Visit all elements without relying on how they are stored. However, connectivity must be established.
The Idx numbering of elements and nodes is not altered by this method.

@return number of nodes that were discovered.

@attention this method assumes that all the nodes are connected to elements.

@attention nodes must have been assigned their parent elements for this method to work.

@note method was formerly called AccumulateAll()

@author SKM 9/20/2008, CSMP Castasegna workshop, Switzerland.
*/
template<size_t dim>
size_t findContiguousMeshPatch( csmp::Node<dim>* const root_node, std::deque<Element<dim>*>& elements )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( root_node == nullptr )
    csmp_error.notice( FATAL_ERROR, "findContiguousMeshPatch(node)", "Root node pointer is dangling!" );

  if ( root_node->Parent( 0 ) == nullptr )
    csmp_error.notice( FATAL_ERROR, "findContiguousMeshPatch(node)",
                       "Root node must have been assigned parent elements; else this method cannot operate." );

  if ( !elements.empty() )
    csmp_error.notice( WARNING, "findContiguousMeshPatch(node)",
                       "supplied element set not empty; deleting all its content." );
  elements.clear();

  // 1. traversal of the existing mesh nodes to find all its elements
  set<csmp::Element<dim>*> explored_elements;
  set<csmp::Node<dim>*>    discovered_nodes;
  deque<csmp::Node<dim>*>  current_nodes;
  // starting at the root element
  discovered_nodes.insert( root_node );
  current_nodes.push_back( root_node );

  // MESH TRAVERSAL
  while ( !current_nodes.empty() ) {
    const csmp::Node<dim>*  n_ptr( *current_nodes.begin() );
    // for all parent elements of the current node
    for ( size_t i = 0U; i<n_ptr->Parents(); i++ ) {
      // for all the nodes of each parent element
      for ( size_t j = 0U; j<n_ptr->Parent( i )->Nodes(); j++ )
        // if this node is not the one from which we started
        if ( j != n_ptr->ParentNodeNumber( i ) ) {
          pair<typename set<csmp::Node<dim>*>::iterator, bool>
            new_node = discovered_nodes.insert( n_ptr->Parent( i )->N( j ) );
          if ( new_node.second ) current_nodes.push_back( n_ptr->Parent( i )->N( j ) );
        }
      // storing the explored element
      explored_elements.insert( n_ptr->Parent( i ) );
    }
    // removing the node from the discovered (but not yet explored) deque
    current_nodes.pop_front();
  }

  // 2. assigning and trimming excess storage from the element pointer vector
  elements.assign( explored_elements.begin(), explored_elements.end() );

  return discovered_nodes.size();

} // end findContiguousMeshPatch (node)

template size_t findContiguousMeshPatch( csmp::Node<1U>* const, std::deque<Element<1U>*>& );
template size_t findContiguousMeshPatch( csmp::Node<2U>* const, std::deque<Element<2U>*>& );
template size_t findContiguousMeshPatch( csmp::Node<3U>* const, std::deque<Element<3U>*>& );







/**  findContiguousMeshPatch

Attempts a floodfill on the supplied set of elements, this so identified
contiguous model region is returned.

The method depends on correct neighbor information.

@section arguments Input Arguments

Input element range iterator: In this region it proceeds to identify
a contiguous domain.

@param cells_contiguous_subset the method returns a set of pointers to
those elements forming the first contiguous domain that
it was able to reach from the supplied iterator.

@section application Application

To break regions into contiguous subdomains.

@author SKM
@date 21/5/2021

*/
template<size_t dim,template<size_t> class CELL>
void findContiguousMeshPatch( CELL<dim>* const eptr, set<CELL<dim>*>& cells_contiguous_subset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( eptr == nullptr ) {
         csmp_error.notice( ERROR, "findContiguousMeshPatch (set)",
                           "supplied element pointer is a null pointer; nothing was done.");
         return;
      }
    // identifying the neighbors of the first element to be looked at
    deque<CELL<dim>*>  neighbor_cells;
    const size_t  neighbors(eptr->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( eptr->Neighbor(i) != nullptr )
        neighbor_cells.push_back( eptr->Neighbor(i) );
 
    // performing the floodfill, starting with an empty set
    if ( !cells_contiguous_subset.empty() )
      cells_contiguous_subset.clear();
   
    // insert the first element into the new subset
    cells_contiguous_subset.insert( static_cast<CELL<dim>*>(eptr) );
      
    // element set for subsequent passes
    deque<CELL<dim>*>  new_neighbor_cells;
   
     while( !neighbor_cells.empty() )
       {
          // 1. loop over those neighbors that are not already part of the deque
          for ( typename deque<CELL<dim>*>::const_iterator
                nit=neighbor_cells.begin(); nit!=neighbor_cells.end(); ++nit )
            // if the element has not already been dealt with
            if ( cells_contiguous_subset.find( (*nit) ) == cells_contiguous_subset.end() ) {
                const size_t  neighbors((*nit)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( (*nit)->Neighbor(j) != nullptr )
                    new_neighbor_cells.push_back( (*nit)->Neighbor(j) );
                cells_contiguous_subset.insert( (*nit) );
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_cells = new_neighbor_cells;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_cells.clear();
      }
     
 } // end findContiguousMeshPatch


template void findContiguousMeshPatch( Element<1U>* const, set<Element<1U>*>& );
template void findContiguousMeshPatch( Element<2U>* const, set<Element<2U>*>& );
template void findContiguousMeshPatch( Element<3U>* const, set<Element<3U>*>& );

template void findContiguousMeshPatch( Face<1U>* const, set<Face<1U>*>& );
template void findContiguousMeshPatch( Face<2U>* const, set<Face<2U>*>& );
template void findContiguousMeshPatch( Face<3U>* const, set<Face<3U>*>& );

template void findContiguousMeshPatch( InterFace<1U>* const, set<InterFace<1U>*>& );
template void findContiguousMeshPatch( InterFace<2U>* const, set<InterFace<2U>*>& );
template void findContiguousMeshPatch( InterFace<3U>* const, set<InterFace<3U>*>& );






/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these a pointer is inserted into the argument deque (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
             @author SKM 14/8/21
*/
template<size_t dim, template<size_t> class CELL>
size_t  findStandAloneMeshPatches( typename deque<CELL<dim>*>::const_iterator begin,
                                   typename deque<CELL<dim>*>::const_iterator end,
                                   map<string,deque<CELL<dim>*> >& mesh_patches )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( begin == end ) {
             csmp_error.notice( WARNING, "findStandAloneMeshPatches:",
                                         "input CELL pointer range is empty." );
            return 0U;
        }
     
      // getting a copy of the element pointers of the mesh
      deque<CELL<dim>*>  cells( begin, end );
      sort( cells.begin(), cells.end() );

      // detecting via a flood-fill whether the group can be partitioned, else nothing is done
      set<CELL<dim>*>  cells_contiguous_subset;
      findContiguousMeshPatch( (*cells.begin()), cells_contiguous_subset );
      
      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( cells.size() <= cells_contiguous_subset.size() ) {
           std::cout <<"\nModel<" << dim << ">::findStandAloneMeshPatches: ";
           std::cout <<"mesh is already contiguous, nothing was done."<< std::endl;
           return 0U;
        }
      // else, we already have found one patch
      size_t  n_patches(1);
    
      // creating new contiguous group from the element subset
      while ( !cells.empty() )
        {
           // creating name for contiguous patch from finite-element type
           string patch_name( parseFiniteElementType( (*cells_contiguous_subset.begin())->FE_Type() ) );
           // appending the patch number and the number of elements in patch
           patch_name +="_patch";
           patch_name += to_string( n_patches );
           patch_name +="__";
           patch_name += to_string( cells_contiguous_subset.size() );
           patch_name +="cells";
           if ( n_patches == 1U ) {
                 std::cout <<"\nfindStandAloneMeshPatches: ";
                 std::cout <<"mesh is divided into disconnected patch(es):\n";
             }
           std::cout <<"\t\t\t'"<< patch_name <<"'";
           std::cout <<" ("<< cells_contiguous_subset.size() <<" elmts)"<< std::endl;

           // storing away the current contiguous subset
           if ( !cells_contiguous_subset.empty() )
             {
                pair<typename map<string,deque<CELL<dim>*> >::iterator,bool>
                  insertion = mesh_patches.insert( make_pair( patch_name,
                                                   move( deque<CELL<dim>*>( cells_contiguous_subset.begin(),
                                                                            cells_contiguous_subset.end() ) ) ) );
                if ( insertion.second == false ) {
                     csmp_error.notice( ERROR, "findStandAloneMeshPatches:", patch_name,
                                               "could not be inserted into patch map" );
                  }
                n_patches++;
              }
           else
             csmp_error.notice( ERROR, "findStandAloneMeshPatches:", patch_name,
                                       "patch contains no elements, nothing was done" );
              
           // deleting the cells that constitute the contiguous subset from the cell storage
           cells.erase( remove_if( cells.begin(), cells.end(),
                                   [&](auto x){ return binary_search( std::begin(cells_contiguous_subset),
                                                                      std::end(cells_contiguous_subset), x ); } ),
                                   cells.end() );

           // if there are no more cells to process the job is done, else the next patch is searched
           if ( !cells.empty() )
             findContiguousMeshPatch( (*cells.begin()), cells_contiguous_subset );
         }

      return n_patches;
    
   } // end findStandAloneMeshPatches

// 3D version
template size_t  findStandAloneMeshPatches( deque<Element<3U>*>::const_iterator,
                                            deque<Element<3U>*>::const_iterator,
                                            map<string,deque<Element<3U>*> >& );

template size_t  findStandAloneMeshPatches( deque<Face<3U>*>::const_iterator,
                                            deque<Face<3U>*>::const_iterator,
                                            map<string,deque<Face<3U>*> >& );

template size_t  findStandAloneMeshPatches( deque<InterFace<3U>*>::const_iterator,
                                            deque<InterFace<3U>*>::const_iterator,
                                            map<string,deque<InterFace<3U>*> >& );













/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these a pointer is inserted into the argument deque (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
*/
template<size_t dim, template<size_t> class CELL>
size_t  findPointersToStandAloneMeshPatches( typename deque<CELL<dim>*>::const_iterator begin,
                                             typename deque<CELL<dim>*>::const_iterator end,
                                             map<CELL<dim>*,MeshPatchAttributes>& root_pointers )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( begin == end ) {
             csmp_error.notice( WARNING, "findPointersToStandAloneMeshPatches:",
                                         "input CELL pointer range is empty." );
            return 0U;
        }
     
      // getting a copy of the element pointers of the mesh
      deque<CELL<dim>*>  cells( begin, end );
      sort( cells.begin(), cells.end() );

      // detecting via a flood-fill whether the group can be partitioned, else nothing is done
      set<CELL<dim>*>  cells_contiguous_subset;
      findContiguousMeshPatch( (*cells.begin()), cells_contiguous_subset );
      
      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( cells.size() <= cells_contiguous_subset.size() ) {
           std::cout <<"\nModel<" << dim << ">::findPointersToStandAloneMeshPatches: ";
           std::cout <<"mesh is already contiguous, nothing was done."<< std::endl;
           return 0U;
        }

      // else partitions can be created
      std::string  group_name("meshPatch");
      std::string  subgroup_name;
      char         num[128];
      size_t       n_subgroups(1);
    
      // creating new contiguous group from the element subset
      while ( !cells.empty() )
        {
           // creating name of contiguous subgroup
           sprintf( num, "%lu", n_subgroups );
           subgroup_name = group_name + num;
           if ( n_subgroups == 1U ) {
                 std::cout <<"\findPointersToStandAloneMeshPatches: ";
                 std::cout <<"mesh is divided into the subregion(s):\n";
             }
           std::cout <<"\t\t\t'"<< subgroup_name <<"'";
           std::cout <<" ("<< cells_contiguous_subset.size() <<" elmts)"<< std::endl;
         
           // subtracting the elements that constitute the new group from the remaining element list
           cells.erase( remove_if( cells.begin(), cells.end(),
                                   [&](auto x){ return binary_search( std::begin(cells_contiguous_subset),
                                                                      std::end(cells_contiguous_subset), x ); } ),
                                   cells.end() );

           // computing the next subset
           if ( cells.empty() ) break;
           else {
                MeshPatchAttributes attributes( cells_contiguous_subset.size(),
                                               parseFiniteElementDimension((*cells.begin())->FE_Type()) );
                
                pair<typename map<CELL<dim>*,MeshPatchAttributes>::iterator,bool>
                  insertion = root_pointers.insert( make_pair( (*cells.begin()), attributes ) );
                if ( insertion.second == false ) {
                     csmp_error.notice( WARNING, "findPointersToStandAloneMeshPatches:",
                                                 "mesh patch could not be inserted into root cell map. Does it already exist?" );
                  }
                findContiguousMeshPatch( (*cells.begin()), cells_contiguous_subset );
             }
           n_subgroups++;
        }

      return n_subgroups;
    
   } // end findPointersToStandAloneMeshPatches

// 3D version
template size_t  findPointersToStandAloneMeshPatches( deque<Element<3U>*>::const_iterator,
                                                      deque<Element<3U>*>::const_iterator,
                                                      map<Element<3U>*,MeshPatchAttributes>& );

template size_t  findPointersToStandAloneMeshPatches( deque<Face<3U>*>::const_iterator,
                                                      deque<Face<3U>*>::const_iterator,
                                                      map<Face<3U>*,MeshPatchAttributes>& );

template size_t  findPointersToStandAloneMeshPatches( deque<InterFace<3U>*>::const_iterator,
                                                      deque<InterFace<3U>*>::const_iterator,
                                                      map<InterFace<3U>*,MeshPatchAttributes>& );





/** relying on the parent element information from its nodes, method finds the neighbor elements for each Face (or boundary) and connects itself them
 
        @return the number of neighbors that were identified
 */
template<size_t dim>
size_t connectNeighborsUsingNodeParents( Element<dim>* const eptr )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( eptr == nullptr ) {
             csmp_error.notice( ERROR, "connectNeighborsUsingNodeParents:",
                               "invalid element pointer" );
            return 0U;
        }
        
     // making search keys from the nodes
     // TODO: the number of nodes might be reduced by ignoring boundary nodes ?!
     set<Node<dim>*>               node_keys( eptr->NodesBegin(), eptr->NodesEnd() );
     //  key             face
     map<set<Node<dim>*>,size_t>  face_keys;
     const size_t                 n_nbors(face_keys.size());
     vector<size_t>               fnids;
     for ( size_t i=0U; i<n_nbors; ++i ) {
          eptr->FE()->NodesOfFace( i, fnids );
          set<Node<dim>*> face_key;
          for ( auto j : fnids ) face_key.insert( eptr->N(j) );
          face_keys.insert( make_pair( face_key, i ) );
          // setting nbor pointers to null
          eptr->Assign( i, static_cast<Element<dim>*>(nullptr) );
       }
       
     // finding the neighbor elements opposite to the element's faces
     // -------------------------------------------------------------
     size_t nbors_found(0U);
     // making a subset of the parent elements that share a sufficient number of nodes to qualify
     set<Element<dim>*> potential_nbors;
     size_t min_face_nodes = (dim != 3) ? 2 : 3;
     if ( dim == 1 ) min_face_nodes = 1;
     // using the element's nodes to find the neighbors
     for ( auto nit=eptr->NodesBegin(); nit!=eptr->NodesEnd(); ++nit )
       for (size_t i=0U; i<(*nit)->Parents(); ++i )
         {
            size_t counter(0U);
            for ( auto it=(*nit)->Parent(i)->NodesBegin(); it!=(*nit)->Parent(i)->NodesEnd(); ++it )
              if ( node_keys.find(*it) != node_keys.end() ) counter++;
            // if the element shares an equal or greater number of nodes than face nodes it is a potential neighbor
            if ( counter >= min_face_nodes ) potential_nbors.insert( (*nit)->Parent(i) );
         }
     
     // searching the subset of elements
     set<Node<dim>*> nbor_face_key;
     for ( auto it : potential_nbors )
       {
          // loop over faces until matching face is found; else report
          const size_t n_faces(it->Neighbors());
          for ( size_t i=0U; i<n_faces; ++i ) {
              it->FE()->NodesOfFace( i, fnids );
              for ( auto j : fnids ) nbor_face_key.insert( it->N(j) );
              // searching & assigning neighbors found
              auto nbor_it(face_keys.find(nbor_face_key));
              if ( nbor_it != face_keys.end() ) {
                   // assigning the neighbor
                   eptr->Assign( (*nbor_it).second, it );
                   // assigning the neighbor's neighbor-element pointer to this new element
                   it->Assign( i, eptr );
                   nbors_found++;
                   // only of face of the neighbor may be connected
                   break;
                }
              nbor_face_key.clear();
           }
       }
       
    return nbors_found;
        
 } // end connectNeighborsUsingNodeParents
 
template size_t connectNeighborsUsingNodeParents( Element<1U>* );
template size_t connectNeighborsUsingNodeParents( Element<2U>* );
template size_t connectNeighborsUsingNodeParents( Element<3U>* );




/**
    Assigns nodes to the Face finding them from the nodes of the higher dimensional neighbors that share the Face.
       
    Uses unordered set of sets to find the interface between the higher dimensional elements.
    
    The nodes of the face are assigned directly
 */
template<size_t dim>
void findNodesViaHigherDimensionalNeighbors( const Element<dim>* const inner_nbor,
                                             const Element<dim>* const outer_nbor,
                                             Face<dim>* const face )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( face == nullptr )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(Face)", "pointer to target Face is not initialised");
   if ( inner_nbor == nullptr )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(Face)", "pointer to inner higher-dim Element not initialised");
   if ( outer_nbor == nullptr )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(Face)", "pointer to outer higher-dim Element  not initialised");
    
   // 1. Creating a map of the faces of the outer element
   vector<size_t> fnids;
   //      face key
   set<set<Node<dim>*> > outer_elmt_faces;

   const size_t n_outer_elmt_faces(outer_nbor->Faces());
   for ( size_t i=0U; i<n_outer_elmt_faces; ++i ) {
        outer_nbor->FE()->NodesOfFace( i, fnids );
        // creating and recording the search key and face number
        set<Node<dim>*> face_nodes;
        const size_t n_face_nodes(fnids.size());
        for ( size_t j=0U; j<n_face_nodes; j++ ) face_nodes.insert( inner_nbor->N( fnids[j] ) );
        outer_elmt_faces.insert( face_nodes );
     }
     
   // 2. Searching the faces of the inner element that matches this face
   const size_t n_inner_elmt_faces(inner_nbor->Faces());
   for ( size_t i=0U; i<n_inner_elmt_faces; ++i ) {
        inner_nbor->FE()->NodesOfFace( i, fnids );
        // creating and recording the search key
        set<Node<dim>*> face_nodes;
        const size_t n_face_nodes(fnids.size());
        for ( size_t j=0U; j<n_face_nodes; j++ ) face_nodes.insert( inner_nbor->N( fnids[j] ) );
        // performing the search
        auto search_it = outer_elmt_faces.find( face_nodes );
        // if a matching face is found
        if ( search_it != outer_elmt_faces.end() ) {
             // assigning the nodes which are in the right order in fnids
             // (remember that the nodes of the Face should match the order at the inner face)
             size_t k(0U);
             for ( size_t j=0U; j<n_face_nodes; j++ )
               face->Assign( k++, inner_nbor->N( fnids[j] ) );
             // ending the search because only one matching neighbor is expected
             break;
          }
     }

 } // end findNodesViaHigherDimensionalNeighbors

template void findNodesViaHigherDimensionalNeighbors( const Element<1>* const, const Element<1>* const, Face<1>* const );
template void findNodesViaHigherDimensionalNeighbors( const Element<2>* const, const Element<2>* const, Face<2>* const );
template void findNodesViaHigherDimensionalNeighbors( const Element<3>* const, const Element<3>* const, Face<3>* const );



/**
    Assigns nodes to the InterFace finding them by searching for collocated nodes in the higher dimensional neiighbor elements that share the Face.
    For each of the nodes their local integer code in in the inner and outer element are stored. With these a search key is created to find the corresponding
    element faces that are needed to recreate the node order.
       
    Uses unordered set of sets to find the interface between the higher dimensional elements.
    
    The nodes of the face are assigned directly
 */
template<size_t dim>
void findNodesViaHigherDimensionalNeighbors( const Element<dim>* const inner_nbor,
                                             const Element<dim>* const outer_nbor,
                                             InterFace<dim>* const interface )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( interface == nullptr )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "pointer to target InterFace is not initialised");
   if ( inner_nbor == nullptr )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "pointer to inner higher-dim Element not initialised");
   if ( outer_nbor == nullptr )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "pointer to outer higher-dim Element  not initialised");
    
   // 1. creating a search map from the nodes of the outer element
   //  key    inner local id, outer local node id
   map<Point<dim>,pair<size_t,size_t> > outer_elmt_nodes;
   // OUTER ELEMENT
   const size_t n_nodes_outer_elmt(outer_nbor->Nodes());
   for ( size_t i=0U; i<n_nodes_outer_elmt; ++i )
     outer_elmt_nodes.insert( make_pair( outer_nbor->N(i)->Coordinate(), make_pair(UINT_MAX,i) ) );
   
   // 2. searching for the shared nodes
   const size_t n_nodes_inner_elmt(inner_nbor->Nodes());
   for ( size_t i=0U; i<n_nodes_inner_elmt; ++i ) {
        auto search_it=outer_elmt_nodes.find( inner_nbor->N(i)->Coordinate() );
        // if the node is shared between the elements its local id is stored
        if ( search_it != outer_elmt_nodes.end() )
          (*search_it).second.first = i;
     }
   
   // 3. creating face_node ID search keys for the inner and outer elements
   set<size_t> inner_nodes, outer_nodes;
   for ( auto it : outer_elmt_nodes )
     if ( it.second.first != UINT_MAX ) {
          inner_nodes.insert( it.second.first );
          outer_nodes.insert( it.second.second );
       }
   
   // 4. searching the faces of the higher-dimensional for the node keys
   vector<size_t> fnids;
   // INNER ELEMENT
   const size_t n_inner_elmt_faces(inner_nbor->Faces());
   bool  inner_face_found(false);
   for ( size_t i=0U; i<n_inner_elmt_faces; ++i ) {
        inner_nbor->FE()->NodesOfFace( i, fnids );
        // creating and recording the search key
        set<size_t> face_nodes( fnids.begin(), fnids.end() );
        // if it matches the interface, the nodes are assigned and the loop is stopped
        if ( face_nodes == inner_nodes ) {
             inner_face_found = true;
             const size_t n_fnids(fnids.size());
             size_t k(0U);
             for ( size_t j=0U; j<n_fnids; ++j )
               interface->Assign( k++, inner_nbor->N( fnids[j] ), INSIDE );
             interface->ParentFaceID( INSIDE, i );
             break;
          }
     }
   if ( !inner_face_found )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "failed to find nodes of inner higher-dim neighbor element");
    
   // OUTER ELEMENT
   const size_t n_outer_elmt_faces(outer_nbor->Faces());
   bool  outer_face_found(false);
   for ( size_t i=0U; i<n_outer_elmt_faces; ++i ) {
        outer_nbor->FE()->NodesOfFace( i, fnids );
        // creating and recording the search key
        set<size_t> face_nodes( fnids.begin(), fnids.end() );
        // if it matches the interface, the nodes are assigned and the loop is stopped
        if ( face_nodes == outer_nodes ) {
             outer_face_found = true;
             const size_t n_fnids(fnids.size());
             size_t k(0U);
             for ( size_t j=0U; j<n_fnids; ++j )
               interface->Assign( k++, outer_nbor->N( fnids[j] ), OUTSIDE );
             interface->ParentFaceID( OUTSIDE, i );
             break;
          }
     }
   if ( !outer_face_found )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "failed to find nodes of outer higher-dim neighbor element");

 } // end findNodesViaHigherDimensionalNeighbors

template void findNodesViaHigherDimensionalNeighbors( const Element<1>* const, const Element<1>* const, InterFace<1>* const );
template void findNodesViaHigherDimensionalNeighbors( const Element<2>* const, const Element<2>* const, InterFace<2>* const );
template void findNodesViaHigherDimensionalNeighbors( const Element<3>* const, const Element<3>* const, InterFace<3>* const );





/**
     Surt's code to efficiently remove a single element from a sorted vector, without preserving sorted order
     (//inline void erase_v4(std::vector<int> &vec, int value)
     
     https://stackoverflow.com/questions/26719144/how-to-erase-a-value-efficiently-from-a-sorted-vector/26720032
 */
template<size_t dim>
void eraseElementPointerFromVector( vector<csmp::Element<dim>*>& vec, const Element<dim>* eptr )
 {
    // get the range in 2*log2(N), N=vec.size()
    auto bounds = std::equal_range(vec.begin(), vec.end(), eptr );

    // calculate the index of the first to be deleted O(1)
    auto last = vec.end() - std::distance(bounds.first, bounds.second);

    // swap the 2 ranges O(equals) , equal = std::distance(bounds.first, bounds.last)
    std::swap_ranges(bounds.first, bounds.second, last);

    // erase the victims O(equals)
    vec.erase(last, vec.end());
}

template void eraseElementPointerFromVector( vector<csmp::Element<1U>*>&, const Element<1U>* );
template void eraseElementPointerFromVector( vector<csmp::Element<2U>*>&, const Element<2U>* );
template void eraseElementPointerFromVector( vector<csmp::Element<3U>*>&, const Element<3U>* );


/* CLIPPING OF FACE NUMBERING FUNCTION FOR INTERFACE

    const size_t neighbors2x( f->Neighbors() * 2 );
    for ( size_t j = 0U; j<neighbors2x; ++j ) {
      InterFace<dim>* const ptr( f->Neighbor( j ) );
      if ( ptr != nullptr ) {
        // building search maps that we will use to find the shared interfaces
        // key=pointset   face iD
        map<set<Point<dim> >, pair<INTERFACE_SIDE, size_t> >   inner_elmt_faces, outer_elmt_faces;
        vector<size_t>  nids;
        // first element
        Element<dim>* e1 = f->InnerParent();
        for ( size_t face = 0U; face<e1->Faces(); ++face ) {
          e1->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j = 0U; j<nids.size(); ++j )
            face_key.insert( e1->N( nids[j] )->Coordinate() );
          outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
        }
        // second element
        Element<dim>* e2 = f->OuterParent();
        for ( size_t face = 0U; face<e2->Faces(); ++face ) {
          e2->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j = 0U; j<nids.size(); ++j )
            face_key.insert( e2->N( nids[j] )->Coordinate() );
          inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
        }

        // 2. finding the shared faces
        bool found( false );
        long64 inner_face_id( -1 ), outer_face_id( -1 );
        for ( auto inner_face : inner_elmt_faces ) {
          for ( auto outer_face : outer_elmt_faces ) {
            if ( inner_face.first == outer_face.first ) {
              inner_face_id = inner_face.second.second;
              outer_face_id = outer_face.second.second;
              found = true;
              break;
            }
          }
          if ( found ) break;
        }

        vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
      }
      else
        vset.Pfvert( eidx, j, REGION_BOUNDARY );
    }

*/ // FACE NUMBERING




} // end csmp
