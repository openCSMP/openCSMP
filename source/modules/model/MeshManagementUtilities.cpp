//
//  MeshManagementUtilities.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 3/7/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "MeshManagementUtilities.h"
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


/**

Attempts a floodfill on the supplied set of elements, this so identified
contiguous model region is returned.

The method depends on correct neighbor information.

@section arguments Input Arguments

Input element range iterator: In this region it proceeds to identify
a contiguous domain.

@param elements_contiguous_subset the method returns a set of pointers to
those elements forming the first contiguous domain that
it was able to reach from the supplied iterator.

@section application Application

To break regions into contiguous subdomains.

*/
template<size_t dim>
void floodFill( Element<dim>* const eptr, set<Element<dim>*>& elements_contiguous_subset )
 {
    assert( eptr != nullptr );
    // identifying the neighbors of the first element to be looked at
    deque<Element<dim>*>  neighbor_elements;
    const size_t  neighbors(eptr->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( eptr->Neighbor(i) != nullptr )
        neighbor_elements.push_back( eptr->Neighbor(i) );
 
    // performing the floodfill, starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( static_cast<Element<dim>*>(eptr) );
      
    // element set for subsequent passes
    deque<Element<dim>*>  new_neighbor_elements;
   
     while( !neighbor_elements.empty() )
       {
          // 1. loop over those neighbors that are not already part of the deque
          for ( typename deque<Element<dim>*>::const_iterator
                nit=neighbor_elements.begin(); nit!=neighbor_elements.end(); ++nit )
            // if the element has not already been dealt with
            if ( elements_contiguous_subset.find( (*nit) ) == elements_contiguous_subset.end() ) {
                const size_t  neighbors((*nit)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( (*nit)->Neighbor(j) != nullptr )
                    new_neighbor_elements.push_back( (*nit)->Neighbor(j) );
                elements_contiguous_subset.insert( (*nit) );
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_elements.clear();
      }
     
 } // end floodFill


template void floodFill( Element<1U>* const, set<Element<1U>*>& );
template void floodFill( Element<2U>* const, set<Element<2U>*>& );
template void floodFill( Element<3U>* const, set<Element<3U>*>& );






/**
    Retrieves and returns the element ids of the first contiguous element patch
    that can be reached by mesh traversal from the starting element.
    
    @attention the elements in the region must have a unique numbering scheme.
*/
template<size_t dim>
void floodFillViaIndexes( const Region<dim>& gref, size_t starting_idx,
                          std::set<size_t>& elements_contiguous_subset )
 {
    assert( starting_idx < gref.Elements() );
    assert( gref.E(starting_idx) != NULL );
    // identifying the neighbors of the first element to be looked at
    deque<size_t>  neighbor_elements;
    const size_t  neighbors(gref.E(starting_idx)->Neighbors());
    for ( size_t i=0U; i<neighbors; i++ )
      if ( gref.E(starting_idx)->Neighbor(i) != NULL )
        neighbor_elements.push_back( gref.E(starting_idx)->Neighbor(i)->Idx() );
 
    // performing the floodfill, starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( starting_idx );
      
    // element set for subsequent passes
    deque<size_t>  new_neighbor_elements;
   
     while( !neighbor_elements.empty() )
       {
          // 1. loop over those neighbors that are not already part of the list
          for ( typename deque<size_t>::const_iterator
                idx=neighbor_elements.begin(); idx!=neighbor_elements.end(); ++idx )
            // if the element has not already been dealt with
            if ( elements_contiguous_subset.find( (*idx) ) == elements_contiguous_subset.end() ) {
                const size_t  neighbors(gref.E(*idx)->Neighbors());
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( size_t j=0U; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( gref.E(*idx)->Neighbor(j) != NULL ) {
                        assert( gref.E(*idx)->Neighbor(j)->Idx() < gref.Elements() );
                        new_neighbor_elements.push_back( gref.E(*idx)->Neighbor(j)->Idx() );
                        elements_contiguous_subset.insert( gref.E(*idx)->Neighbor(j)->Idx() );
                    }
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
            
          // 3. emptying neighbor set for the next loop
          new_neighbor_elements.clear();
      }
     
 } // end floodFillViaIndexes


template void floodFillViaIndexes( const Region<1U>&, size_t, set<size_t>& );
template void floodFillViaIndexes( const Region<2U>&, size_t, set<size_t>&  );
template void floodFillViaIndexes( const Region<3U>&, size_t, set<size_t>&  );




/**
    Connects Element objects to their same-dimensional neighbors in as much as is possible.
    
    Where there are no neighbors the neighbor pointers will be nulled.
    
    @attention The assumption is made that all nodes in the model have a unique numbering.
    
    @author SKM 2012
    
        @todo deal with manifolds, disambiguating them on the basis of element orientation (only elements int the same plane or aligned elements should be neighbors)

        @TODO: no need to use set, use sort() then unique() on the result vector, searching will be much faster
*/
template<size_t dim>
void  establishNeighborConnectivity( vector<Element<dim>*>& simplexVector, bool unassign_neighbors_outside, bool verbose )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( Element )", "supplied element vector is empty; nothing was done" );
         return;
      }

    cout << "\nestablishNeighborConnectivity( Element ): Establishing CSMP FE neighbor connectivity...\n";
 
    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    if (verbose) cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >  volume_neighbor_keys,
                                                           surface_neighbor_keys,
                                                           line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<Element<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( size_t face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == nullptr ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( Element )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }
           // creating face key from idx's of face
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( size_t j=0U; j<fnids.size(); j++ )
               key.insert( (*it)->N( fnids[j] ) );
             
           // inserting newly generated keys into multimap
           if ( (*it)->FE()->IsVolumeElement() )
               volume_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else // for all line elements
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           key.clear();

           // unassign neighbors outside of the provided vector range
           if ( unassign_neighbors_outside )
           {
             // remove element from neighbor list of its neighbors
             const size_t neighbors( (*it)->Neighbors() );
             for ( size_t neighbor = 0; neighbor < neighbors; ++neighbor )
               if ( (*it)->Neighbor( neighbor ) != NULL ) {
                 const size_t neighbor_neighbors( (*it)->Neighbor( neighbor )->Neighbors() );
                 for ( size_t i = 0; i < neighbor_neighbors; i++ )
                   if ( (*it)->Neighbor( neighbor )->Neighbor( i ) != NULL )
                     if ( (*it)->Neighbor( neighbor )->Neighbor( i ) == (*it) ) {
                       (*it)->NeighborElementVector()[i] = NULL;
                     }
                 
                 (*it)->NeighborElementVector()[neighbor] = NULL;
               }
             (*it)->NeighborElementVector().clear();
           }
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)

    if (verbose) cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(nullptr);
        Element<dim>* e2Ptr(nullptr);

        if (verbose) cout << "\n\t\tline elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                                 it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // line elements
    
    // 2.2 surface elements
    // --------------------
    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(nullptr);
        Element<dim>* e2Ptr(nullptr);

        if (verbose) cout << "\n\t\tsurface elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                                 it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                          face pointer                   nbor face idx        neighbor pointer
                    ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                    ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( dim == 3U and !volume_neighbor_keys.empty() ) {

        Element<dim>* e1Ptr(nullptr);
        Element<dim>* e2Ptr(nullptr);

        if (verbose) cout << "\n\t\tvolume elements...\n";

        //                key             n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,Element<dim>*> >::iterator it1(volume_neighbor_keys.begin()),
                                                                                 it2(volume_neighbor_keys.begin());

        it2++;

        while ( it2 != volume_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ( (*it1).second.second != 0 and (*it2).second.second != 0 ) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                   //                          face pointer                   nbor face idx        neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == volume_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=3
    
 } // end establishNeighborConnectivity

// explicit instantiations
template void establishNeighborConnectivity<1U>( std::vector<csmp::Element<1U>*>&, bool, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::Element<2U>*>&, bool, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::Element<3U>*>&, bool, bool );




/** CONNECTIVITY BETWEEN INTERFACES
       
         Inside and outside must be considered.
*/
template<size_t dim>
void  establishNeighborConnectivity( std::vector<csmp::InterFace<dim>*>& simplexVector, bool unassign_neighbors_outside, bool verbose )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( simplexVector.empty() ) {
         csmp_error.notice( WARNING, "establishNeighborConnectivity( InterFace ):", "supplied element vector is empty; nothing was done." );
         return;
      }

    if (verbose) cout << "\nestablishNeighborConnectivity( InterFace ): Establishing CSMP FE neighbor connectivity...\n";

    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    if (verbose)  cout << "  Building element face list...\n";

    //       key             face number neighbor
    multimap<set<Node<dim>*>,pair<size_t,InterFace<dim>*> >  surface_neighbor_keys,
                                                             line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<InterFace<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( size_t face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == NULL ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }

           // creating face keys from idx's of interface for INSIDE & OUTSIDE
           (*it)->CurrentSide( INSIDE ); // just for node-vector
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( size_t j=0U; j<fnids.size(); j++ )
             key.insert( (*it)->N( fnids[j] ) );

           // inserting newly generated key into multimap
           if ( (*it)->FE()->IsSurfaceElement() )
               surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else if ( (*it)->FE()->IsLineElement() )
               line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
           else{
               csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )", "supplied element vector contains volumetric element! nothing was done" );
               return;
           }
           key.clear();

           // unassign neirghbors outside of the provided vector range
           if ( unassign_neighbors_outside )
           {
             // remove element from neighbor list of its neighbors
             const size_t neighbors( (*it)->Neighbors() );
             for ( size_t neighbor = 0; neighbor < neighbors; ++neighbor )
               if ( (*it)->Neighbor( neighbor ) != NULL ) {
                 const size_t neighbor_neighbors( (*it)->Neighbor( neighbor )->Neighbors() );
                 for ( size_t i = 0; i < neighbor_neighbors; i++ )
                   if ( (*it)->Neighbor( neighbor )->Neighbor( i ) != NULL )
                     if ( (*it)->Neighbor( neighbor )->Neighbor( i ) == (*it) ) {
                       (*it)->NeighborElementVector()[i] = NULL;
                     }

                 (*it)->NeighborElementVector()[neighbor] = NULL;
               }
             (*it)->NeighborElementVector().clear();
           }
        }


    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    if (verbose) cout << "  Building element neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() ) {

        InterFace<dim>* e1Ptr(NULL);
        InterFace<dim>* e2Ptr(NULL);

        if (verbose) cout << "\n\t\tline elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,InterFace<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                                   it2(line_neighbor_keys.begin());

        it2++;

        while ( it2 != line_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning eachothers faces
                    //                        face pointer, neighbor pointer, side-of interface
                   (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                   (*it2).second.second->Assign( (*it2).second.first, e1Ptr );

                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=1

    // 2.2 surface elements
    // --------------------
    if ( dim >= 2U and !surface_neighbor_keys.empty() ) {

        InterFace<dim>* e1Ptr(NULL);
        InterFace<dim>* e2Ptr(NULL);

        if (verbose) cout << "\n\t\tsurface elements...";

        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,InterFace<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                                   it2(surface_neighbor_keys.begin());

        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          {
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first and ((*it1).second.second != 0 and (*it2).second.second != 0) )
                {
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                    // assigning eachothers faces
                    //                        face pointer, neighbor pointer, side-of interface
                    (*it1).second.second->Assign( (*it1).second.first, e2Ptr );
                    (*it2).second.second->Assign( (*it2).second.first, e1Ptr );


                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          }
      } // dim=2


 } // end establishNeighborConnectivity

template void establishNeighborConnectivity<1U>( std::vector<csmp::InterFace<1U>*>&, bool, bool );
template void establishNeighborConnectivity<2U>( std::vector<csmp::InterFace<2U>*>&, bool, bool );
template void establishNeighborConnectivity<3U>( std::vector<csmp::InterFace<3U>*>&, bool, bool );








// MESH DIAGNOSTICS

/**
   Tests whether a tetrahedron is degenerate because all of its vertices lie within a single plane
*/
bool hasNonManifoldVertices( const csmp::Element<3>* const tptr, double64 tolerance )
 {
    assert( tptr != nullptr );
    assert( tptr->FE()->ElementType() == LINEAR_TETRAHEDRON ||
            tptr->FE()->ElementType() == ISOPARAMETRIC_LINEAR_TETRAHEDRON );
            
    // making a tensor of the 3 vectors from vertex 0 to the other vertices
    // using the coordinates of the tetrahedral element
    tptr->CoordinateMatrix();
    TensorVariable<3U> M;
    Point<3U> vec = tptr->N(1)->Coordinate() - tptr->N(0)->Coordinate();
    M(0,0) = vec[0]; M(0,1) = vec[1]; M(0,2) = vec[2];
    vec = tptr->N(2)->Coordinate() - tptr->N(0)->Coordinate();
    M(1,0) = vec[0]; M(1,1) = vec[1]; M(1,2) = vec[2];
    vec = tptr->N(3)->Coordinate() - tptr->N(0)->Coordinate();
    M(2,0) = vec[0]; M(2,1) = vec[1]; M(2,2) = vec[2];
    
    if ( M.Determinant() <= tolerance ) return false;
    return true;
    
 } // end hasNonManifoldVertices




/**
    Computes barycentre-to-node distances for all parent elements of node and returns them into the supplied vector
    distances_and_weight vector [e1,e2...e_n,e_sum] with a size of parent elements+1
*/
template<size_t dim>
void distanceWeights( typename vector<Node<dim>*>::const_iterator nodes_begin,
                      typename vector<Node<dim>*>::const_iterator nodes_end,
                      vector<vector<double64> >& distances_and_weights )
 {
     distances_and_weights.resize(distance(nodes_begin,nodes_end));
     size_t  node_index(0U);
   
     while( nodes_begin != nodes_end )
       {
          double64 weight(0.);
          const Point<dim> npt((*nodes_begin)->Coordinate());
          const size_t parents((*nodes_begin)->Parents());
          distances_and_weights[node_index].reserve(parents+1U);
         
          for ( size_t i=0U; i<parents; i++ ) {
               // recording the node-to-barycentre distances
               distances_and_weights[node_index].push_back( npt.DistanceTo( (*nodes_begin)->Parent(i)->BaryCenter() ) );
               weight += distances_and_weights[node_index][i];
            }
          // storing the weights (=sum of the distances) in the last element of the vector
          distances_and_weights[node_index][parents] = weight;
          nodes_begin++;
          node_index++;
       }
   
 } // end distanceWeights

template void distanceWeights<1>( vector<Node<1>*>::const_iterator, vector<Node<1>*>::const_iterator, vector<vector<double64> >& );
template void distanceWeights<2>( vector<Node<2>*>::const_iterator, vector<Node<2>*>::const_iterator, vector<vector<double64> >& );
template void distanceWeights<3>( vector<Node<3>*>::const_iterator, vector<Node<3>*>::const_iterator, vector<vector<double64> >& );




/**
    finds those elements in a model that contact eachother across split interfaces.
    Only those elements are discovered that are node-matched.
 
    @param subdomain (non-unique) region which contains the split boundary
 
    @return returns false if none of the perimeter elements are node matched
 
    @attention method only looks at highest dimensional elements in the model;
    thus, lower dimensional elements are ignored and the neighborhood relations
    on either side of them are returned.
 
    @attention the elements inside the model that are located along the splitboundary
    do not have neighbor pointers yet.
 
    @author SKM
    @date 14/01/2018
 
    @TODO do we need to remember which side of the interface we are on?
    @todo test method on Chloe's dataset
*/
template<size_t dim>
bool findSplitInterfaceElements( const Region<dim>& subdomain,
                                 set<pair<pair<Element<dim>*,size_t>,pair<Element<dim>*,size_t> > >& interface_elmt_pairs )
 {
    // multi-element container for all elements that are located on split boundaries
    //       face search key      element      face number
    multimap<set<Point<dim> >,pair<Element<dim>*,size_t> > element_face_keys;
   
    // 1. for all elements on the perimeter of the model subdomain,
    //    generate keys from their node coordinates that are then matched with one-another
    //    in order to connect these elements
    for ( size_t n=subdomain.InteriorElements(); n<subdomain.Elements(); ++n )
        {
           // for those element faces that define the perimeter surface
           for ( size_t i=0U; i<subdomain.PerimeterFaces(n); ++i )
             {
                // get the local node numbers of the perimeter face
                vector<size_t> fnids;
                subdomain.E(n)->FE()->NodesOfFace( subdomain.PerimeterFace( n, i ), fnids );
                // add the corresponding node points to a set that will form the element face key
                pair<set<Point<dim> >,size_t> face_key;
                for ( size_t j=0U; j<fnids.size(); ++j )
                  face_key.first.insert( subdomain.E(n)->N( fnids[j] )->Coordinate() );
                // remembering the face id
                face_key.second = subdomain.PerimeterFace( n, i );
                // storing the key in the correspondance search map
                //                              node-coordinate set  element pointer   local face id
                element_face_keys.insert( make_pair( face_key.first, make_pair( subdomain.E(n), face_key.second ) ) );
             }
        }
   
    // 2. searching map for matching interface elements
    //    - a match is obtained if the element pointers are different
    //    - if there is a match, the element pair is added to the OppositeElement container
    if ( !interface_elmt_pairs.empty() ) interface_elmt_pairs.clear();
    // searching
    for ( auto it=element_face_keys.begin(); it!=element_face_keys.end(); ++it ) {
         // multimap iterator pair containing the range of shared keys
         //  face search key (point set), element(Element*), face number (size_t)
         auto result = element_face_keys.equal_range( (*it).first );
         // if more than one value was found
         if ( distance( result.first, result.second ) > 1U ) {
           // if there are any manyfolds, then stop this process and return false.
           //if ( distance( result.first, result.second ) == 2U ) continue;
           // advancing the result range iterator as necessary to find an Element different from (*it).second.first
           do result.first++;
           while ( (*result.first).second.first == (*it).second.first );
           
           // building search maps that we will use to find the shared interfaces
           // key=pointset   face iD
           map<set<Point<dim> >, pair<INTERFACE_SIDE, size_t> >   inner_elmt_faces, outer_elmt_faces;
           vector<size_t>  nids;
           // first element
           Element<dim>* e1 = (*it).second.first;
           for ( size_t face = 0U; face<e1->Faces(); ++face ) {
             e1->FE()->NodesOfFace( face, nids );
             set<Point<dim> >  face_key;
             for ( size_t j = 0U; j<nids.size(); ++j )
               face_key.insert( e1->N( nids[j] )->Coordinate() );
             outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
           }
           // second element
           Element<dim>* e2 = (*result.first).second.first;
           for ( size_t face = 0U; face<e2->Faces(); ++face ) {
             e2->FE()->NodesOfFace( face, nids );
             set<Point<dim> >  face_key;
             for ( size_t j = 0U; j<nids.size(); ++j )
               face_key.insert( e2->N( nids[j] )->Coordinate() );
             inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
           }

           // 2. finding the shared faces
           bool found( false );
           long64 inner_face_id(-1), outer_face_id(-1);
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

           // we store the matching element pair
           // set<pair<pair<Element<dim>*,size_t>,pair<Element<dim>*,size_t> > >
           if ( found )
           interface_elmt_pairs.insert( make_pair(
                                        make_pair( (*it).second.first, (*it).second.second ),
                                        make_pair( (*result.first).second.first, (*result.first).second.second ) )
                                      );
          }
      }
  
    // 3. checking the results
    if (  interface_elmt_pairs.empty() ) return false;
 
    return true;
    
 } // end findSplitInterfaceElements

template bool findSplitInterfaceElements( const Region<3U>&, set<pair<pair<Element<3U>*,size_t>,pair<Element<3U>*,size_t> > >& );
template bool findSplitInterfaceElements( const Region<2U>&, set<pair<pair<Element<2U>*,size_t>,pair<Element<2U>*,size_t> > >& );
template bool findSplitInterfaceElements( const Region<1U>&, set<pair<pair<Element<1U>*,size_t>,pair<Element<1U>*,size_t> > >& );






/**
    Checks whether region of interest contains elements of the type of interest.
*/
template<size_t dim>
bool containsElementsOfType( const Region<dim>& gref, ELEMENT_DIMENSION dimension )
 {
    for ( typename vector<Element<dim>*>::const_iterator
          it=gref.ElementsBegin(); it!=gref.ElementsEnd(); it++ )
      if ( parseFiniteElementDimension( (*it)->FE_Type() ) == dimension )
        return true;
      
    return false;
      
 } // end containsElementsOfTtype

template bool containsElementsOfType<1U>( const Region<1>&, ELEMENT_DIMENSION );
template bool containsElementsOfType<2U>( const Region<2>&, ELEMENT_DIMENSION );
template bool containsElementsOfType<3U>( const Region<3>&, ELEMENT_DIMENSION );



/**
    loops over the surface cells of the model subdomain,
    checking whether any of the projections of the normals of the neighbor
    cells onto the normal of the current element are negative.
    
    @return if any of the projections is negative, method returns false.
    A negative projection will mean that the normals of the neighboring cells
    are at > to 90^o to the current cell; this should not be the case
    unless the surface has cusps with in it.
    
    @return if the region does not consist entirely of surface elements
    an error is reported.
    
    @attention even if the normals have a consistent orientation, this method
    may return false if the surface contains a cusp (>90^o kink).
    
    @test SKM 30/3/2016 for 3D only
 
*/
template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<3U>&  subdomain )
 {
    vector<double64> normal(3U), nbor_normal(3U);
   
    size_t non_surface_elements(0U);
    for ( auto it=subdomain.ElementsBegin(); it!=subdomain.ElementsEnd(); ++it )
      // this method only considers surface elements
      if ( (*it)->IsSurfaceElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( size_t i=0U; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double64 result(0.);
                  for ( size_t j=0U; j<3U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_surface_elements++;
   
    if ( non_surface_elements > 0U )
      ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation (3D):",
                                       subdomain.Name(), "region contained not only surface elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation


template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<2U>&  subdomain )
 {
    vector<double64> normal(2U), nbor_normal(2U);
   
    size_t non_line_elements(0U);
    for ( auto it=subdomain.ElementsBegin(); it!=subdomain.ElementsEnd(); ++it )
      // this method only considers line elements
      if ( (*it)->IsLineElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( size_t i=0U; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double64 result(0.);
                  for ( size_t j=0U; j<2U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_line_elements++;
   
    if ( non_line_elements > 0U )
      ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation (2D):",
                                       subdomain.Name(), "region contained not only line elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation





/**
    Stub: ID there are no boundaries so this should be a compile time assert
    
    @todo TODO: use static_assert<> here on the template argument
*/
template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<1U>&  subdomain )
 {
    ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation (1D):",
                                     subdomain.Name(), "one-dimensional models have no boundaries." );
    return false;
   
 } // end checkNeighborNormalsForConsistentOrientation








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
