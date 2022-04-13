//
//  MeshManagementUtilities.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 3/7/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "MeshManagementUtilities.h"
#include "MeshManager.h"
#include "Model.h"
#include "Region.h"
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
template<uint32_t dim>
size_t detectElementsWithAllNodesOnBoundary( const MeshManager<dim>& mmgr, set<size_t>& belmts )
{
  belmts.clear();

  assert( mmgr.Elements() > 0 );
  if ( mmgr.Elements() == 0 ) return 0U;

  // traversal of the existing mesh nodes to find all its elements
  size_t boundary_only_elements( 0U );
  for ( typename plf::colony<Element<dim>>::const_iterator
        it=mmgr.ElementsBegin(); it!=mmgr.ElementsEnd(); ++it ) {
      const auto nodes((*it).Nodes());
      uint32_t   counter{0};
      for ( auto i{0U}; i<nodes; ++i )
        if ( (*it).N(i)->AtBoundary() != NOT ) counter++;
      if ( counter == nodes ) {
            belmts.insert( (*it).Idx() );
            boundary_only_elements++;
        }
    }
  return boundary_only_elements;
}

template size_t detectElementsWithAllNodesOnBoundary( const MeshManager<2U>&, set<size_t>& );
template size_t detectElementsWithAllNodesOnBoundary( const MeshManager<3U>&, set<size_t>& );









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
template<uint32_t dim,template<uint32_t> class CELL>
size_t findContiguousMeshPatch( CELL<dim>* const eptr, set<CELL<dim>*>& cells_contiguous_subset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( eptr == nullptr ) {
         csmp_error.notice( ERROR, "findContiguousMeshPatch", "entry cell pointer is a nullptr; nothing was done.");
         return 0;
      }

    // clean the set if it is not empty
    if ( !cells_contiguous_subset.empty() ) cells_contiguous_subset.clear();
    // insert the first element into the new subset
    cells_contiguous_subset.insert( eptr );
    
    vector<CELL<dim>*>  neighbor_cells( eptr->NeighborsBegin(), eptr->NeighborsEnd() );
    assert( neighbor_cells.size() >= 1 );
    const uint32_t      max_cell_nbors{6}; // assuming that hexa has the most
      
    while( !neighbor_cells.empty() )
      {
          vector<CELL<dim>*>  new_neighbor_cells;
          // 1. loop over those neighbors that are not already part of the deque
          for ( auto& nit : neighbor_cells )
             // if the cell has not been encountered before
             if ( nit && cells_contiguous_subset.find(nit) == cells_contiguous_subset.end() ) {
                 const auto  n_neighbors{ nit->Neighbors() };
                 assert( n_neighbors <= max_cell_nbors );
                 new_neighbor_cells.reserve( n_neighbors );
                 for ( auto j{0U}; j<n_neighbors; ++j )
                   if ( nit->Neighbor(j) != nullptr )
                     new_neighbor_cells.push_back( nit->Neighbor(j) );
                 // marking the neighbor cell as discovered
                 cells_contiguous_subset.insert( nit );
              }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_cells = new_neighbor_cells;
          new_neighbor_cells.clear();
       }
  
    return cells_contiguous_subset.size();
     
 } // end findContiguousMeshPatch


template size_t findContiguousMeshPatch( Element<1U>* const, set<Element<1U>*>& );
template size_t findContiguousMeshPatch( Element<2U>* const, set<Element<2U>*>& );
template size_t findContiguousMeshPatch( Element<3U>* const, set<Element<3U>*>& );

template size_t findContiguousMeshPatch( Face<1U>* const, set<Face<1U>*>& );
template size_t findContiguousMeshPatch( Face<2U>* const, set<Face<2U>*>& );
template size_t findContiguousMeshPatch( Face<3U>* const, set<Face<3U>*>& );

template size_t findContiguousMeshPatch( InterFace<1U>* const, set<InterFace<1U>*>& );
template size_t findContiguousMeshPatch( InterFace<2U>* const, set<InterFace<2U>*>& );
template size_t findContiguousMeshPatch( InterFace<3U>* const, set<InterFace<3U>*>& );





/** loops over the valid neighbors of the cell and sets their neighbor pointers to point to this cell to nullptr
 */
template<uint32_t dim, template<uint32_t> class CELL>
void detachNeighborsFrom( CELL<dim>* const eptr )
 {
    // nulling the connections of neighbor neighbor elements to this element
    // (neighbor pointer to this element is nulled)
    for ( auto i{0U}; i<eptr->Neighbors(); ++i )
      if ( eptr->Neighbor(i) != nullptr )
        for ( auto j{0U}; j<eptr->Neighbor(i)->Neighbors(); ++j )
        if ( eptr->Neighbor(i)->Neighbor(j) == eptr )
          eptr->Neighbor(i)->Neighbor(j)->Unassign( eptr );
 }

template void detachNeighborsFrom( Element<1>* const );
template void detachNeighborsFrom( Element<2>* const );
template void detachNeighborsFrom( Element<3>* const );

template void detachNeighborsFrom( Face<1>* const );
template void detachNeighborsFrom( Face<2>* const );
template void detachNeighborsFrom( Face<3>* const );

template void detachNeighborsFrom( InterFace<1>* const );
template void detachNeighborsFrom( InterFace<2>* const );
template void detachNeighborsFrom( InterFace<3>* const );



/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these, a pointer is inserted into the argument vector (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
             @author SKM 14/8/21
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  findStandAloneMeshPatches( typename plf::colony<CELL<dim>>::iterator begin,
                                   typename plf::colony<CELL<dim>>::iterator end,
                                   map<string,vector<CELL<dim>*> >& mesh_patches )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( begin == end ) {
             csmp_error.notice( WARNING, "findStandAloneMeshPatches:",
                                         "input CELL pointer range is empty." );
            return 0U;
        }
     
      // getting a copy of the element pointers of the mesh (without writing an adaptor from colony iterator to cell pointer)
      vector<CELL<dim>*>   cells;
      cells.reserve( distance(begin,end) );
      while ( begin != end ) {
          cells.push_back( &(*begin) );
          ++begin;
        }

      // detecting via a flood-fill whether the group can be partitioned, else nothing is done
      set<CELL<dim>*>  cells_contiguous_subset;
      findContiguousMeshPatch<dim>( (*cells.begin()), cells_contiguous_subset );

      // if no contiguous cells could be found
      if ( cells_contiguous_subset.empty() ) {
           csmp_error.notice( WARNING, "findStandAloneMeshPatches:",
                                       "No contiguous cells found. Is the neighbor connectivity missing? - nothing was done" );
           return 0U;
        }

      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( cells.size() <= cells_contiguous_subset.size() ) {
           std::cout <<"\nfindStandAloneMeshPatches: ";
           std::cout <<"mesh is already contiguous, nothing was done."<< std::endl;
           
           // creating name for contiguous patch from finite-element type
           string patch_name( parseFiniteElementType( (*cells_contiguous_subset.begin())->FE_Type() ) );
           // appending the patch number and the number of elements in patch
           patch_name +="_patch";
           patch_name += to_string( 1 );
           patch_name +="_";
           patch_name += to_string( cells_contiguous_subset.size() );
           patch_name +="cells";
           mesh_patches.insert( make_pair( patch_name, move( vector<CELL<dim>*>{ cells_contiguous_subset.begin(),
                                                                                 cells_contiguous_subset.end() } ) ) );
           return 1U;
        }
      // else, we already have found one patch
      size_t  n_patches(1);
    
      // creating new contiguous group from the element subset
      while ( !cells.empty() )
        {
           // storing away the current contiguous subset
           if ( !cells_contiguous_subset.empty() )
             {
                // creating name for contiguous patch from finite-element type
                string patch_name( parseFiniteElementType( (*cells_contiguous_subset.begin())->FE_Type() ) );
                // appending the patch number and the number of elements in patch
                patch_name +="_patch";
                patch_name += to_string( n_patches );
                patch_name +="_";
                patch_name += to_string( cells_contiguous_subset.size() );
                patch_name +="cells";
                if ( n_patches == 1U ) {
                     cout <<"\nfindStandAloneMeshPatches: ";
                     cout <<"mesh is divided into disconnected patch(es):\n";
                  }
                cout <<"\t\t\t'"<< patch_name <<"'";
                cout <<" ("<< cells_contiguous_subset.size() <<" elmts)"<< std::endl;

                pair<typename map<string,vector<CELL<dim>*> >::iterator,bool>
                  insertion = mesh_patches.insert( make_pair( patch_name,
                                                   move( vector<CELL<dim>*>( cells_contiguous_subset.begin(),
                                                                             cells_contiguous_subset.end() ) ) ) );
                if ( insertion.second == false ) {
                     csmp_error.notice( ERROR, "findStandAloneMeshPatches:", patch_name,
                                               "could not be inserted into patch map" );
                  }
                n_patches++;
              }
           else csmp_error.notice( ERROR, "findStandAloneMeshPatches:", "patch contains no elements, nothing was done" );
              
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
template size_t  findStandAloneMeshPatches( plf::colony<Element<3U>>::iterator,
                                            plf::colony<Element<3U>>::iterator,
                                            map<string,vector<Element<3U>*> >& );

template size_t  findStandAloneMeshPatches( plf::colony<Face<3U>>::iterator,
                                            plf::colony<Face<3U>>::iterator,
                                            map<string,vector<Face<3U>*> >& );

template size_t  findStandAloneMeshPatches( plf::colony<InterFace<3U>>::iterator,
                                            plf::colony<InterFace<3U>>::iterator,
                                            map<string,vector<InterFace<3U>*> >& );

// 2D version
template size_t  findStandAloneMeshPatches( plf::colony<Element<2U>>::iterator,
                                            plf::colony<Element<2U>>::iterator,
                                            map<string,vector<Element<2U>*> >& );

template size_t  findStandAloneMeshPatches( plf::colony<Face<2U>>::iterator,
                                            plf::colony<Face<2U>>::iterator,
                                            map<string,vector<Face<2U>*> >& );

template size_t  findStandAloneMeshPatches( plf::colony<InterFace<2U>>::iterator,
                                            plf::colony<InterFace<2U>>::iterator,
                                            map<string,vector<InterFace<2U>*> >& );




/**
      tested:OK
*/
template<uint32_t dim>
size_t findInterconnectedNodeCluster( Node<dim>* const nptr, std::set<Node<dim>*>& contiguous_set_of_nodes )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    if ( nptr == nullptr ) {
         csmp_error.notice( ERROR, "findInterconnectedNodeCluster", "entry node is a nullpointer; nothing was done");
         return 0U;
      }
    
    if ( !contiguous_set_of_nodes.empty() ) contiguous_set_of_nodes.clear();

    vector<Node<dim>*>  node_neighbors{ nptr }, new_node_nbors;
    contiguous_set_of_nodes.insert( nptr );
    // cerr <<"\n\nfindInterconnectedNodeCluster: traversing nodes, starting at: "<< nptr->Idx() <<": ";

    while ( !node_neighbors.empty() )
      {
        for ( const auto& nit : node_neighbors )
          {
             const size_t n_node_nbors{ nit->Neighbors() };
             new_node_nbors.reserve( n_node_nbors );
             for ( auto i{0U}; i < n_node_nbors; ++i ) {
                 Node<dim>* nbor_ptr = nit->Neighbor(i);
                 assert( nbor_ptr != nullptr );
                 if ( contiguous_set_of_nodes.find( nbor_ptr ) == contiguous_set_of_nodes.end() ) {
                      new_node_nbors.push_back( nbor_ptr );
                      contiguous_set_of_nodes.insert( nbor_ptr );
                      // cerr <<" "<< nbor_ptr->Idx();
                   }
               }
           }
        node_neighbors = new_node_nbors;
        new_node_nbors.clear();
     }
    
  return contiguous_set_of_nodes.size();
      
} // end findInterconnectedNodeCluster


template size_t findInterconnectedNodeCluster( Node<3>* const nptr, set<Node<3>*>& contiguous_set_of_nodes );
template size_t findInterconnectedNodeCluster( Node<2>* const nptr, set<Node<2>*>& contiguous_set_of_nodes );
template size_t findInterconnectedNodeCluster( Node<1>* const nptr, set<Node<1>*>& contiguous_set_of_nodes );





/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these a pointer is inserted into the argument deque (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  findPointersToStandAloneMeshPatches( typename vector<CELL<dim>*>::const_iterator begin,
                                             typename vector<CELL<dim>*>::const_iterator end,
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
      findContiguousMeshPatch<dim>( (*cells.begin()), cells_contiguous_subset );
      
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
                findContiguousMeshPatch<dim>( (*cells.begin()), cells_contiguous_subset );
             }
           n_subgroups++;
        }

      return n_subgroups;
    
   } // end findPointersToStandAloneMeshPatches

// 3D version
template size_t  findPointersToStandAloneMeshPatches( vector<Element<3U>*>::const_iterator,
                                                      vector<Element<3U>*>::const_iterator,
                                                      map<Element<3U>*,MeshPatchAttributes>& );

template size_t  findPointersToStandAloneMeshPatches( vector<Face<3U>*>::const_iterator,
                                                      vector<Face<3U>*>::const_iterator,
                                                      map<Face<3U>*,MeshPatchAttributes>& );

template size_t  findPointersToStandAloneMeshPatches( vector<InterFace<3U>*>::const_iterator,
                                                      vector<InterFace<3U>*>::const_iterator,
                                                      map<InterFace<3U>*,MeshPatchAttributes>& );





/** relying on the parent element information from its nodes, method finds the neighbor elements for each Face (or boundary) and connects itself them
 
        @return the number of neighbors that were identified
 */
 // TODO: simplify by using the corner nodes only for the face matching
template<uint32_t dim>
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
     map<set<Node<dim>*>,uint32_t>  face_keys;
     const size_t                   n_nbors(face_keys.size());
     vector<uint32_t>               fnids;
     for ( auto i{0U}; i<n_nbors; ++i ) {
          eptr->FE()->NodesOfFace( i, fnids );
          set<Node<dim>*> face_key;
          for ( auto j : fnids ) face_key.insert( eptr->N(j) );
          face_keys.insert( make_pair( face_key, i ) );
          // setting nbor pointers to null
          eptr->Assign( i, static_cast<Element<dim>*>(nullptr) );
       }
       
     // finding the neighbor elements opposite to the element's faces
     // -------------------------------------------------------------
     uint32_t nbors_found(0U);
     // making a subset of the parent elements that share a sufficient number of nodes to qualify
     set<Element<dim>*> potential_nbors;
     uint32_t min_face_nodes = (dim != 3) ? 2 : 3;
     if ( dim == 1 ) min_face_nodes = 1;
     // using the element's nodes to find the neighbors
     for ( auto nit=eptr->NodesBegin(); nit!=eptr->NodesEnd(); ++nit )
       for (auto i{0U}; i<(*nit)->Parents(); ++i )
         {
            uint32_t counter(0U);
            for ( auto it=(*nit)->Parent(i)->NodesBegin(); it!=(*nit)->Parent(i)->NodesEnd(); ++it )
              if ( node_keys.find(*it) != node_keys.end() ) counter++;
            // if the element shares an equal or greater number of nodes than face nodes it is a potential neighbor
            if ( counter >= min_face_nodes ) potential_nbors.insert( (*nit)->Parent(i) );
         }
     
     // searching the subset of elements
     set<Node<dim>*> nbor_face_key;
     for ( auto& it : potential_nbors )
       {
          // loop over faces until matching face is found; else report
          const auto n_faces(it->Neighbors());
          for ( auto i{0U}; i<n_faces; ++i ) {
              it->FE()->NodesOfFace( i, fnids );
              for ( auto& j : fnids ) nbor_face_key.insert( it->N(j) );
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
       
    Uses  set of sets to find the interface between the higher dimensional elements.
    
    The nodes of the face are assigned directly
 */
template<uint32_t dim>
void findNodesViaHigherDimensionalNeighbors( Element<dim>* const inner_nbor,
                                             Element<dim>* const outer_nbor,
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
   //      face key
   set<set<Node<dim>*> > outer_elmt_faces;

   const auto n_outer_elmt_faces(outer_nbor->Faces());
   for ( auto i{0U}; i<n_outer_elmt_faces; ++i ) {
        // creating and recording the search key and face number
        outer_elmt_faces.insert( outer_nbor->CornerNodesOfFace(i) );
     }
     
   // 2. Searching the faces of the inner element that matches this face
   const auto n_inner_elmt_faces(inner_nbor->Faces());
   for ( auto i{0U}; i<n_inner_elmt_faces; ++i ) {
        // creating the search key and performing the search
        auto search_it = outer_elmt_faces.find( inner_nbor->CornerNodesOfFace(i) );
        // if a matching face is found
        if ( search_it != outer_elmt_faces.end() ) {
             // assigning the nodes which are in the right order in fnids
             // (remember that the nodes of the Face should match the order at the inner face)
             vector<uint32_t> fnids;
             inner_nbor->FE()->NodesOfFace( i, fnids );
             uint32_t k(0U);
             for ( auto& j : fnids )
               face->Assign( k++, inner_nbor->N(j) );
             // ending the search because only one matching neighbor is expected
             break;
          }
     }

 } // end findNodesViaHigherDimensionalNeighbors

template void findNodesViaHigherDimensionalNeighbors( Element<1>* const, Element<1>* const, Face<1>* const );
template void findNodesViaHigherDimensionalNeighbors( Element<2>* const, Element<2>* const, Face<2>* const );
template void findNodesViaHigherDimensionalNeighbors( Element<3>* const, Element<3>* const, Face<3>* const );



/**
    Assigns nodes to the InterFace finding them by searching for collocated nodes in the higher dimensional neiighbor elements that share the Face.
    For each of the nodes their local integer code in in the inner and outer element are stored. With these a search key is created to find the corresponding
    element faces that are needed to recreate the node order.
       
    Uses unordered set of sets to find the interface between the higher dimensional elements.
    
    The nodes of the face are assigned directly
 */
template<uint32_t dim>
void findNodesViaHigherDimensionalNeighbors( Element<dim>* const inner_nbor,
                                             Element<dim>* const outer_nbor,
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
   map<Point<dim>,pair<uint32_t,uint32_t> > outer_elmt_nodes;
   // OUTER ELEMENT
   const auto n_nodes_outer_elmt(outer_nbor->Nodes());
   for ( auto i{0U}; i<n_nodes_outer_elmt; ++i )
     outer_elmt_nodes.insert( make_pair( outer_nbor->N(i)->Coordinate(), make_pair(numeric_limits<uint32_t>::max(),i) ) );
   
   // 2. searching for the shared nodes
   const auto n_nodes_inner_elmt(inner_nbor->Nodes());
   for ( auto i{0U}; i<n_nodes_inner_elmt; ++i ) {
        auto search_it=outer_elmt_nodes.find( inner_nbor->N(i)->Coordinate() );
        // if the node is shared between the elements its local id is stored
        if ( search_it != outer_elmt_nodes.end() )
          (*search_it).second.first = i;
     }
   
   // 3. creating face_node ID search keys for the inner and outer elements
   set<uint32_t> inner_nodes, outer_nodes;
   for ( auto& it : outer_elmt_nodes )
     if ( it.second.first != numeric_limits<uint32_t>::max() ) {
          inner_nodes.insert( it.second.first );
          outer_nodes.insert( it.second.second );
       }
   
   // 4. searching the faces of the higher-dimensional for the node keys
   vector<uint32_t> fnids;
   // INNER ELEMENT
   const auto n_inner_elmt_faces(inner_nbor->Faces());
   bool  inner_face_found(false);
   for ( auto i{0U}; i<n_inner_elmt_faces; ++i ) {
        inner_nbor->FE()->NodesOfFace( i, fnids );
        // creating and recording the search key
        set<uint32_t> face_nodes( fnids.begin(), fnids.end() );
        // if it matches the interface, the nodes are assigned and the loop is stopped
        if ( face_nodes == inner_nodes ) {
             inner_face_found = true;
             const size_t n_fnids(fnids.size());
             uint32_t k(0U);
             for ( auto j{0U}; j<n_fnids; ++j )
               interface->Assign( k++, inner_nbor->N( fnids[j] ), INSIDE );
             interface->ParentFaceID( INSIDE, i );
             break;
          }
     }
   if ( !inner_face_found )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "failed to find nodes of inner higher-dim neighbor element");
    
   // OUTER ELEMENT
   const auto n_outer_elmt_faces(outer_nbor->Faces());
   bool  outer_face_found(false);
   for ( auto i{0U}; i<n_outer_elmt_faces; ++i ) {
        outer_nbor->FE()->NodesOfFace( i, fnids );
        // creating and recording the search key
        set<uint32_t> face_nodes( fnids.begin(), fnids.end() );
        // if it matches the interface, the nodes are assigned and the loop is stopped
        if ( face_nodes == outer_nodes ) {
             outer_face_found = true;
             const auto n_fnids(fnids.size());
             auto k(0U);
             for ( auto j{0U}; j<n_fnids; ++j )
               interface->Assign( k++, outer_nbor->N( fnids[j] ), OUTSIDE );
             interface->ParentFaceID( OUTSIDE, i );
             break;
          }
     }
   if ( !outer_face_found )
     csmp_error.notice( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "failed to find nodes of outer higher-dim neighbor element");

 } // end findNodesViaHigherDimensionalNeighbors

template void findNodesViaHigherDimensionalNeighbors( Element<1>* const, Element<1>* const, InterFace<1>* const );
template void findNodesViaHigherDimensionalNeighbors( Element<2>* const, Element<2>* const, InterFace<2>* const );
template void findNodesViaHigherDimensionalNeighbors( Element<3>* const, Element<3>* const, InterFace<3>* const );




/*
    Finds the face between 2 elements (if any) via the neighbor connectivity.
*/
template<uint32_t dim>
pair<size_t,size_t> findAdjacentFacesFromNeighbors( Element<dim>* const eptr1, Element<dim>* const eptr2 )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   if ( eptr1 == nullptr ) {
        csmp_error.notice( ERROR, "findAdjacentFacesFromNeighbors", "null pointer to first element.");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
   if ( eptr2 == nullptr ) {
        csmp_error.notice( ERROR, "findAdjacentFacesFromNeighbors", "null pointer to second element.");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
   if ( eptr1 == eptr2 ) {
        csmp_error.notice( ERROR, "findAdjacentFacesFromNeighbors", "the supplied pointers point to the same element!");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
    
    // 1. finding which face of element 1 is shared with element 2
    size_t face_elmt1 = UNSPECIFIED;
    const size_t n_faces1(eptr1->Faces());
    for ( auto i{0U}; i<n_faces1; ++i )
      if ( eptr1->Neighbor(i) == eptr2 ) {
           face_elmt1 = i;
           break;
        }

    // 2. finding which face of element 2 is shared with element 1
    size_t face_elmt2 = UNSPECIFIED;
    const size_t n_faces2(eptr2->Faces());
    for ( auto i{0U}; i<n_faces2; ++i )
      if ( eptr2->Neighbor(i) == eptr1 ) {
           face_elmt2 = i;
           break;
        }

    return make_pair( face_elmt1, face_elmt2 );
    
 } // end findAdjacentElementFaces
 
template pair<size_t,size_t> findAdjacentFacesFromNeighbors( Element<3>* const, Element<3>* const );
template pair<size_t,size_t> findAdjacentFacesFromNeighbors( Element<2>* const, Element<2>* const );
template pair<size_t,size_t> findAdjacentFacesFromNeighbors( Element<1>* const, Element<1>* const );
 




/**
     Finds the adjacent faces of the supplied elements via their shared nodes.
     
     @return pair of the local face ID numbers of element one and two.

     @attention if neighbor connectivty exists, use matching neighbor pointers which is much faster!
     @attention if no shared face can be found, function returns UNSPECIFIED.
*/
template<uint32_t dim>
pair<size_t,size_t> findAdjacentElementFaces( Element<dim>* const eptr1, Element<dim>* const eptr2 )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   if ( eptr1 == nullptr ) {
        csmp_error.notice( ERROR, "findAdjacentElementFaces", "null pointer to first element.");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
   if ( eptr2 == nullptr ) {
        csmp_error.notice( ERROR, "findAdjacentElementFaces", "null pointer to second element.");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
   if ( eptr1 == eptr2 ) {
        csmp_error.notice( ERROR, "findAdjacentElementFaces", "the supplied pointers point to the same element!");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
    
    // finding the shared face by their nodes
     // 1. creating unique keys from the nodes of the first elements faces not located at a boundary
    vector<set<Node<dim>*> > e1_face_keys;
    const size_t n_faces(eptr1->Faces());
    e1_face_keys.reserve(n_faces);
    for ( auto i{0U}; i<n_faces; ++i )
      if ( eptr1->Neighbor(i) != nullptr ) {
          e1_face_keys.emplace_back( eptr1->CornerNodesOfFace(i) );
       }
      else e1_face_keys.emplace_back( set<Node<dim>*>{} );
    
    // 2. creating face keys for the outer element trying match them
    //    with the faces of the inner one
    set<uint32_t>    face_key_n;
    const size_t n_faces2(eptr2->Faces());
    for ( auto i{0U}; i<n_faces2; ++i ) {
         // is this a matching face
         auto fit = find( e1_face_keys.begin(), e1_face_keys.end(), eptr2->CornerNodesOfFace(i) );
         if ( fit != e1_face_keys.end() ) {
              size_t face_elmt1 = distance(e1_face_keys.begin(),fit);
              size_t face_elmt2 = i;
              return make_pair( face_elmt1, face_elmt2 );
           }
      }

    return make_pair( UNSPECIFIED, UNSPECIFIED );
    
 } // end findAdjacentElementFaces
 
template pair<size_t,size_t> findAdjacentElementFaces( Element<3>* const, Element<3>* const );
template pair<size_t,size_t> findAdjacentElementFaces( Element<2>* const, Element<2>* const );
template pair<size_t,size_t> findAdjacentElementFaces( Element<1>* const, Element<1>* const );
 






/**
     Surt's code to efficiently remove a single element from a sorted vector, without preserving sorted order
     (//inline void erase_v4(std::vector<int> &vec, int value)
     
     https://stackoverflow.com/questions/26719144/how-to-erase-a-value-efficiently-from-a-sorted-vector/26720032
 */
template<uint32_t dim>
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
template<uint32_t dim, template<uint32_t> class CELL>
void floodFill( CELL<dim>* const eptr, set<CELL<dim>*>& elements_contiguous_subset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( eptr == nullptr ) {
         csmp_error.notice( ERROR, "floodFill", "root cell pointer is a nullptr; nothing was done.");
         return;
      }
    // identifying the neighbors of the first element to be looked at
    vector<CELL<dim>*>  neighbor_elements;
    const size_t  neighbors(eptr->Neighbors());
    neighbor_elements.reserve( neighbors );
    for ( auto i{0U}; i<neighbors; i++ )
      if ( eptr->Neighbor(i) != nullptr )
        neighbor_elements.push_back( eptr->Neighbor(i) );
 
    // starting with an empty set
    if ( !elements_contiguous_subset.empty() )
      elements_contiguous_subset.clear();
   
    // insert the first element into the new subset
    elements_contiguous_subset.insert( eptr );
      
     while( !neighbor_elements.empty() )
       {
          // 0. element set for subsequent passes
          vector<CELL<dim>*>  new_neighbor_elements;
   
          // 1. loop over those neighbors that are not already part of the deque
          for ( const auto& nit : neighbor_elements )
            // if the element has not already been dealt with
            if ( nit != nullptr && elements_contiguous_subset.find( nit ) == elements_contiguous_subset.end() ) {
                const auto  n_neighbors(nit->Neighbors());
                new_neighbor_elements.reserve( n_neighbors );
                // adding its neighbor ids to the element list to be processed next, if they haven't been dealt with already
                for ( auto j{0U}; j<neighbors; ++j )
                  // if there is a neighbor whose neighbors have not been traversed, it is input in the list
                  if ( nit->Neighbor(j) != nullptr )
                    new_neighbor_elements.push_back( nit->Neighbor(j) );
                elements_contiguous_subset.insert( nit );
             }
 
          // 2. obtain a new set of neighbors that has to be visited in the next iteration
          neighbor_elements = new_neighbor_elements;
      }
     
 } // end floodFill


template void floodFill( Element<1U>* const, set<Element<1U>*>& );
template void floodFill( Element<2U>* const, set<Element<2U>*>& );
template void floodFill( Element<3U>* const, set<Element<3U>*>& );

template void floodFill( Face<1U>* const, set<Face<1U>*>& );
template void floodFill( Face<2U>* const, set<Face<2U>*>& );
template void floodFill( Face<3U>* const, set<Face<3U>*>& );

template void floodFill( InterFace<1U>* const, set<InterFace<1U>*>& );
template void floodFill( InterFace<2U>* const, set<InterFace<2U>*>& );
template void floodFill( InterFace<3U>* const, set<InterFace<3U>*>& );







/**
    Connects Element objects to their same-dimensional neighbors in as much as is possible.
    
    Where there are no neighbors the neighbor pointers will be nulled.
    
    @attention The assumption is made that all nodes in the model have a unique numbering.
    
    @author SKM 2012
    
        @todo deal with manifolds, disambiguating them on the basis of element orientation (only elements int the same plane or aligned elements should be neighbors)

        @TODO: no need to use set, use sort() then unique() on the result vector, searching will be much faster
*/
template<uint32_t dim>
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
    multimap<set<Node<dim>*>,pair<uint32_t,Element<dim>*> >  volume_neighbor_keys,
                                                             surface_neighbor_keys,
                                                             line_neighbor_keys;
    vector<uint32_t>               fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<Element<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( auto face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == nullptr ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( Element )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }
           // creating face key from idx's of face
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( auto j{0U}; j<fnids.size(); j++ )
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
             const auto neighbors( (*it)->Neighbors() );
             for ( auto neighbor = 0; neighbor < neighbors; ++neighbor )
               if ( (*it)->Neighbor( neighbor ) != NULL ) {
                 const auto neighbor_neighbors( (*it)->Neighbor( neighbor )->Neighbors() );
                 for ( auto i = 0; i < neighbor_neighbors; i++ )
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

        auto it1(line_neighbor_keys.begin()), it2(line_neighbor_keys.begin());

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

        auto it1(surface_neighbor_keys.begin()), it2(surface_neighbor_keys.begin());

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

        auto it1(volume_neighbor_keys.begin()), it2(volume_neighbor_keys.begin());

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
template<uint32_t dim>
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
    multimap<set<Node<dim>*>,pair<uint32_t,InterFace<dim>*> >  surface_neighbor_keys,
                                                               line_neighbor_keys;
    vector<uint32_t>               fnids;
    typename std::set<Node<dim>*>  key;

    for ( typename vector<InterFace<dim>*>::const_iterator it = simplexVector.begin(); it!=simplexVector.end(); it++ )
      for ( auto face=0U; face<(*it)->Faces(); face++ )
        {
           if ( (*it) == NULL ) {
                csmp_error.notice( ERROR, "establishNeighborConnectivity( InterFace )",
                                  "supplied element contains NULL pointer to elements; nothing was done" );
                return;
             }

           // creating face keys from idx's of interface for INSIDE & OUTSIDE
           (*it)->CurrentSide( INSIDE ); // just for node-vector
           (*it)->FE()->NodesOfFace( face, fnids );
           for ( size_t j{0U}; j<fnids.size(); j++ )
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
             const auto neighbors( (*it)->Neighbors() );
             for ( auto neighbor = 0; neighbor < neighbors; ++neighbor )
               if ( (*it)->Neighbor( neighbor ) != NULL ) {
                 const auto neighbor_neighbors( (*it)->Neighbor( neighbor )->Neighbors() );
                 for ( auto i = 0; i < neighbor_neighbors; i++ )
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

        auto it1(line_neighbor_keys.begin()), it2(line_neighbor_keys.begin());

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

        auto it1(surface_neighbor_keys.begin()), it2(surface_neighbor_keys.begin());

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
bool hasNonManifoldVertices( const csmp::Element<3>* const tptr, double tolerance )
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
    For the supplied range of nodes, for each node,
    computes parent element barycentre-to-node distances for all parent elements and returns them into the supplied vector
    distances_and_weight vector [e1,e2...e_n,e_sum] with a size of parent elements+1
*/
template<uint32_t dim>
void distancesAndWeights( typename vector<Node<dim>*>::const_iterator nodes_begin,
                          typename vector<Node<dim>*>::const_iterator nodes_end,
                          vector<vector<double> >& distances_and_weights )
 {
     distances_and_weights.resize(distance(nodes_begin,nodes_end));
     size_t  node_index(0U);
   
     while( nodes_begin != nodes_end )
       {
          double weight(0.);
          const Point<dim> npt((*nodes_begin)->Coordinate());
          const size_t parents((*nodes_begin)->Parents());
          distances_and_weights[node_index].reserve(parents+1U);
         
          for ( auto i{0U}; i<parents; i++ ) {
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

template void distancesAndWeights<1>( vector<Node<1>*>::const_iterator, vector<Node<1>*>::const_iterator, vector<vector<double> >& );
template void distancesAndWeights<2>( vector<Node<2>*>::const_iterator, vector<Node<2>*>::const_iterator, vector<vector<double> >& );
template void distancesAndWeights<3>( vector<Node<3>*>::const_iterator, vector<Node<3>*>::const_iterator, vector<vector<double> >& );




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
template<uint32_t dim>
bool findSplitInterfaceElements( const Region<dim>& subdomain,
                                 set<pair<pair<Element<dim>*,size_t>,pair<Element<dim>*,size_t> > >& interface_elmt_pairs )
 {
    // multi-element container for all elements that are located on split boundaries
    //       face search key      element      face number
    multimap<set<Point<dim> >,pair<Element<dim>*,size_t> > element_face_keys;
   
    // 1. for all elements on the perimeter of the model subdomain,
    //    generate keys from their node coordinates that are then matched with one-another
    //    in order to connect these elements
    for ( auto n=subdomain.InteriorCells(); n<subdomain.Cells(); ++n )
        {
           // for those element faces that define the perimeter surface
           for ( auto i{0U}; i<subdomain.PerimeterFaces(n); ++i )
             {
                // get the local node numbers of the perimeter face
                vector<uint32_t> fnids;
                subdomain.E(n)->FE()->NodesOfFace( subdomain.PerimeterFace( n, i ), fnids );
                // add the corresponding node points to a set that will form the element face key
                pair<set<Point<dim> >,size_t> face_key;
                for ( size_t j{0U}; j<fnids.size(); ++j )
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
           // key =pointset                         face iD
           map<set<Point<dim> >,pair<INTERFACE_SIDE,uint32_t> >   inner_elmt_faces, outer_elmt_faces;
           vector<uint32_t>  nids;
           // first element
           Element<dim>* e1 = (*it).second.first;
           for ( auto face = 0U; face<e1->Faces(); ++face ) {
             e1->FE()->NodesOfFace( face, nids );
             set<Point<dim> >  face_key;
             for ( auto j{0U}; j<nids.size(); ++j )
               face_key.insert( e1->N( nids[j] )->Coordinate() );
             outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
           }
           // second element
           Element<dim>* e2 = (*result.first).second.first;
           for ( auto face = 0U; face<e2->Faces(); ++face ) {
             e2->FE()->NodesOfFace( face, nids );
             set<Point<dim> >  face_key;
             for ( size_t j{0U}; j<nids.size(); ++j )
               face_key.insert( e2->N( nids[j] )->Coordinate() );
             inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
           }

           // 2. finding the shared faces
           bool found( false );
           int64_t  inner_face_id(-1), outer_face_id(-1);
           for ( auto& inner_face : inner_elmt_faces ) {
             for ( auto& outer_face : outer_elmt_faces ) {
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
template<uint32_t dim>
bool containsElementsOfType( const Region<dim>& gref, CELL_SHAPE dimension )
 {
    for ( auto it=gref.CellsBegin(); it!=gref.CellsEnd(); it++ )
      if ( parseFiniteElementDimension( (*it)->FE_Type() ) == dimension )
        return true;
      
    return false;
      
 } // end containsElementsOfTtype

template bool containsElementsOfType<1U>( const Region<1>&, CELL_SHAPE );
template bool containsElementsOfType<2U>( const Region<2>&, CELL_SHAPE );
template bool containsElementsOfType<3U>( const Region<3>&, CELL_SHAPE );



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
    vector<double> normal(3U), nbor_normal(3U);
   
    size_t non_surface_elements(0U);
    for ( auto it=subdomain.CellsBegin(); it!=subdomain.CellsEnd(); ++it )
      // this method only considers surface elements
      if ( (*it)->IsSurfaceElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( auto i{0U}; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double result(0.);
                  for ( size_t j{0U}; j<3U; ++ j )
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
    vector<double> normal(2U), nbor_normal(2U);
   
    size_t non_line_elements(0U);
    for ( auto it=subdomain.CellsBegin(); it!=subdomain.CellsEnd(); ++it )
      // this method only considers line elements
      if ( (*it)->IsLineElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( auto i{0U}; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double result(0.);
                  for ( size_t j{0U}; j<2U; ++ j )
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



template<template<uint32_t> class CELL>
double angleBetweenSurfaceCells( const CELL<3>* const cell1, const CELL<3>* const cell2 )
 {
    assert( cell1 != nullptr );
    assert( cell2 != nullptr );
    assert( cell1->FE()->IsSurfaceElement() );
    assert( cell2->FE()->IsSurfaceElement() );
    
    Point<3> nrml1( cell1->UnitNormal() );
    Point<3> nrml2( cell2->UnitNormal() );
    
    // cos theta = dot-product over cross-product (length1 * length2)
    // already in degrees
    return angleBetweenEdges( make_pair( Point<3U>{0.,0.,0.}, nrml1 ),
                              make_pair( Point<3U>{0.,0.,0.}, nrml2 ) );
    
 } // end angleBetweenSurfaceCells

template double angleBetweenSurfaceCells<Element>( const Element<3>* const, const Element<3>* const );
template double angleBetweenSurfaceCells<Face>( const Face<3>* const, const Face<3>* const );
template double angleBetweenSurfaceCells<InterFace>( const InterFace<3>* const, const InterFace<3>* const );




/**
    Line elements can exist in all 3 spatial dimensions.
*/
template<uint32_t dim, template<uint32_t> class CELL>
double angleBetweenLineCells( const CELL<dim>* const cell1, const CELL<dim>* const cell2 )
 {
    assert( cell1 != nullptr );
    assert( cell2 != nullptr );
    assert( cell1->FE()->IsLineElement() );
    assert( cell2->FE()->IsLineElement() );
    
    // cos theta = dot-product over cross-product (length1 * length2)
    // already in degrees
    return angleBetweenEdges( make_pair( cell1->N(0)->Coordinate(), cell1->N(1)->Coordinate() ),
                              make_pair( cell2->N(0)->Coordinate(), cell2->N(1)->Coordinate() ) );
    
 } // end angleBetweenSurfaceCells

template double angleBetweenLineCells<3,Element>( const Element<3>* const, const Element<3>* const );
template double angleBetweenLineCells<3,Face>( const Face<3>* const, const Face<3>* const );
template double angleBetweenLineCells<3,InterFace>( const InterFace<3>* const, const InterFace<3>* const );

template double angleBetweenLineCells<2,Element>( const Element<2>* const, const Element<2>* const );
template double angleBetweenLineCells<2,Face>( const Face<2>* const, const Face<2>* const );
template double angleBetweenLineCells<2,InterFace>( const InterFace<2>* const, const InterFace<2>* const );






/**
      Using, neighbor and node to parent relationships, finds edges of volumetric parent elements which share the node.
      Returns the node pointers into the output vector, using the same order in which they appear in the parent element.
      
      @param edge_nodes sorted vector of edge nodes.
      @param segm_parents map of the parents found with the nodes that they match along their segments.
      @param find_segment_ids if this is requested the segment number that the edge corresponds to is reported via the elements IDX.
      @return the number of the parent elements that could be associated with the nodes.
      
      @attention if something went wrong, the Idx() of the involved nodes and volumetric elements will carry the values UINT_MAX.
      
      @attention if one of the nodes cannot be associated with a parent, this is reported.
      
      Uses node-to-parent connectivity.
      
      @attention does not identify which segment of the parent element matches the nodes.
      
      @todo check wether this method breaks for quadratic and cubic elements?
*/
size_t parentElementsSharingMultipleEdgeNodes( const vector<Node<3U>*>&  edge_nodes,
                                               map<Element<3>*,vector<Node<3>*> >& segm_parents,
                                               bool find_segment_ids  )
 {
     // debugging output for Paraview
     cerr <<"\nx,y,z,node_id,boundary_flag,";
     for ( const auto& nit : edge_nodes ) {
          cerr <<"\n"<< nit->x() <<","<< nit->y() <<","<< nit->z() <<","<< nit->Idx() <<","<< parseBoundary(nit->AtBoundary()) <<",";
       }
     cerr << endl;
     
     // setting edge node Idx to UINT_MAX so that processed edge nodes can be detected
     // without the use of another container
     for ( const auto& nit : edge_nodes ) nit->Idx( numeric_limits<size_t>::max() );
     
     // loop over the parent elements of the nodes
     size_t temp_idx{0};
     for ( const auto& nit : edge_nodes ) {
          node_processing:
          const auto n_parents{ nit->Parents() };
          for ( auto i{0U}; i<n_parents; ++i )
            if ( nit->Parent(i)->IsVolumeElement() )
              {
                 // checking which of the neighbor nodes of the parent element
                 // are also contained in the shared nodes vector
                 const auto n_node_nbors{ nit->Neighbors() };
                 for( auto j{0U}; j<n_node_nbors; ++j ) {
                   assert( nit->Neighbor(j) != nullptr );
                   if (  nit->Neighbor(j)->Idx() == numeric_limits<size_t>::max() && // if the node has not been encountered before
                         binary_search( edge_nodes.begin(), edge_nodes.end(), nit->Neighbor(j) ) )
                     {
                        // we have found a volumetric element whose edge contains 2 of the edge_nodes
                        segm_parents.insert( make_pair( nit->Parent(i), vector<Node<3>*>{ nit, nit->Neighbor(j) } ) );
                        // there can only be one shared edge per volume element
                        nit->Idx( temp_idx++ );
                        nit->Neighbor(j)->Idx( temp_idx++ );
                        goto node_processing;
                     }
                 }
              }
       }
       
     // if not all nodes could be associated with volumetric elements
     if ( temp_idx < edge_nodes.size() ) {
          cerr <<"\nparentElementsSharingMultipleEdgeNodes: only "<< temp_idx;
          cerr <<" of the "<< edge_nodes.size() <<" could be matched with segment nodes of volumetric elements.\n";
       }
       
     // finding the segments of the parent elements that the node-pairs match with
     // the segment Idx values are assigned to the elements
     if ( find_segment_ids ) {
         for ( auto& it : segm_parents ) {
              it.first->Idx( numeric_limits<size_t>::max() );
              const auto n_segments{ it.first->Segments() };
              for ( auto segm{0}; segm < n_segments; ++segm ) {
                   vector<uint32_t> snids;
                   it.first->FE()->NodesOfSegment( segm, snids );
                   // making set of segment node pointers to search for
                   set<Node<3>*> segm_nodes;
                   const size_t n_segm_nodes{ snids.size() };
                   for ( size_t j{0}; j<n_segm_nodes; ++j )
                     segm_nodes.insert( it.first->N( snids[j] ) );
                   // searching segm_nodes for the nodes previously associated with the element
                   bool all_nodes_found{true};
                   const size_t n_edge_nodes{ it.second.size() };
                   for ( size_t j{0}; j<n_edge_nodes; ++j )
                     if ( segm_nodes.find( it.second[j] ) == segm_nodes.end() ) {
                          all_nodes_found = false;
                          break;
                       }
                   if ( all_nodes_found ) {
                        // assigning the segment number i to IDX
                        it.first->Idx( segm );
                        // onto the next element
                        break;
                     }
                }
              // element loop gets here only if edge cannot be matched with segment
              cerr <<"\nparentElementsSharingMultipleEdgeNodes: Segment ID could not be found for Element: ";
              cerr << parseFiniteElementType( it.first->FE_Type() ) <<" with the nodes:\n\t";
              for ( auto j{0U}; j<it.first->Nodes(); ++j ) cerr <<" "<< it.first->N(j)->Idx();
              cerr <<"\n\tEdge nodes: "<< it.second[0]->Idx() <<" "<< it.second[1]->Idx();
           }
       }
       
     return segm_parents.size();
     
 } // end parentElementsSharingMultipleEdgeNodes







template<uint32_t dim, template<uint32_t> class CELL>
void backupNeighborConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                 typename std::vector<CELL<dim>*>::const_iterator last,
                                 std::vector<std::vector<CELL<dim>*> >& nbor_pointers )
 {
    nbor_pointers.clear();
    nbor_pointers.reserve( distance(first,last) );
    assert( nbor_pointers.capacity() > 1 );
    
    while ( first != last ) {
         const size_t n_nbors{ (*first)->Neighbors() };
         vector<CELL<dim>*>  nbors( n_nbors, nullptr );
         for ( auto i{0U}; i<n_nbors; ++i )
           if ( (*first)->Neighbor(i) != nullptr )
             nbors.push_back( (*first)->Neighbor(i) );
         nbor_pointers.emplace_back( nbors );
         first++;
      }
 
 } // end backupNeighborConnectivity
 
template void backupNeighborConnectivity( typename vector<Element<3>*>::const_iterator,
                                          typename vector<Element<3>*>::const_iterator,
                                          vector<vector<Element<3>*> >& );

template void backupNeighborConnectivity( typename vector<Face<3>*>::const_iterator,
                                          typename vector<Face<3>*>::const_iterator,
                                          vector<vector<Face<3>*> >& );

template void backupNeighborConnectivity( typename vector<Element<2>*>::const_iterator,
                                          typename vector<Element<2>*>::const_iterator,
                                          vector<vector<Element<2>*> >& );

template void backupNeighborConnectivity( typename vector<Face<2>*>::const_iterator,
                                          typename vector<Face<2>*>::const_iterator,
                                          vector<vector<Face<2>*> >& );
                                          
                                          



/**
   Checks validity of FE policy, nodes, neighbors, node parents, node neighbors.
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool integrityCheck( typename plf::colony<CELL<dim>>::const_iterator first,
                     typename plf::colony<CELL<dim>>::const_iterator last )
 {
    if ( first == last ) return false;
    
    const size_t              n_cells( distance(first,last) );
    vector<const Node<dim>*>  shared_nodes;
    size_t                    issues{0};
    shared_nodes.reserve( n_cells * dim );
    string                    celltype("Element");
    if constexpr ( is_same< CELL<dim>,Face<dim> >::value ) celltype = "Face";
    if constexpr ( is_same< CELL<dim>,InterFace<dim> >::value ) celltype = "InterFace";

    auto   copy_of_first{ first};
    size_t max_cell_idx{0};

    // 1. checking that all cells stored in the container are valid
    // ------------------------------------------------------------
    const string check1("\nAre all cells stored in the container valid?\n");
    bool first_call{true};
    
    while ( first != last ) {
         max_cell_idx = max( max_cell_idx, (*first).Idx() );
         // FE policy
         if ( (*first).FE() == nullptr ) {
              if ( first_call ) { cerr << check1; first_call=false; }
              cerr <<"\n"<< celltype << (*first).Idx() <<": FE pointer corrupt.";
              issues++;
           }
         first++;
      }
       
       
    // 2. checking that all connections between nodes and cells are valid
    // ------------------------------------------------------------------
    const string check2("\nintegrityCheck: Are all the nodes connected to the cell valid?\n");
    // global cell-id,local node-id
    multimap<size_t,uint32_t>  missing_nodes;
    first = copy_of_first;
    first_call = true;
    
    while ( first != last ) {
         // connected nodes
         for ( auto i{0U}; i<(*first).Nodes(); ++i ) {
               if ( (*first).N(i) == nullptr ) {
                    if ( first_call ) { cerr << check2; first_call=false; }
                    cerr <<"\n"<< celltype << (*first).Idx() <<": node: "<< i <<": node pointer corrupt.";
                    missing_nodes.insert( make_pair( (*first).Idx(), i ) );
                    issues++;
                 }
               else shared_nodes.push_back( (*first).N(i) );
           }
         first++;
      }
    // TODO: write the missing nodes to a file
    
    
   // 3. checking node parent connectivity after removing duplicate nodes
   // -------------------------------------------------------------------
   const string check4("\nintegrityCheck: Are all the parent elements of the nodes valid?\n");
   first_call = true;
   sort( shared_nodes.begin(), shared_nodes.end() );
   shared_nodes.erase( unique(shared_nodes.begin(), shared_nodes.end()), shared_nodes.end() );
   
   for ( const auto& nit : shared_nodes ) {
       if ( nit == nullptr ) cerr <<"\ndetected 'nullptr' node.";
       else
         for ( auto i{0U}; i<nit->Parents(); ++i )
            if ( nit->Parent(i) == nullptr ||
                 nit->Parent(i)->FE() == nullptr ) {
                 if ( first_call ) { cerr << check4; first_call=false; }
                 cerr <<"\nNode "<< nit->Idx() <<": parents vector contains nullptr.";
                 issues++;
              }
     }


   // 4. node to node connectivity is tested
   // --------------------------------------
   const string check5("\nintegrityCheck: Are all node-to-node connections valid?\n");
   first_call = true;
   for ( const auto& nit : shared_nodes ) {
       if ( nit == nullptr ) cerr <<"\ndetected 'nullptr' node.";
       else
         for ( auto i{0U}; i<nit->Neighbors(); ++i )
            if ( nit->Neighbor(i) == nullptr ) {
                 if ( first_call ) { cerr << check5; first_call=false; }
                 cerr <<"\nneighbor "<< i <<" of Node "<< nit->Idx() <<": is corrupt.";
                 issues++;
              }
     }
   
   
    // 5. checking that all cells have at least one neighbor
    // -----------------------------------------------------
    const string check3("\nintegrityCheck: Are there cells without any neighbors?\n");
    multimap<size_t,uint32_t>  missing_nbors;
    first = copy_of_first;
    first_call = true;
   
    while ( first != last ) {
           // there should be at least one cell neighbor
           size_t n_valid_nbors{0};
           for ( auto i{0U}; i<(*first).Neighbors(); ++i ) {
                 if ( (*first).Neighbor(i) != nullptr ) n_valid_nbors++;
                 else {
                      missing_nbors.insert( make_pair( (*first).Idx(), i ) );
                   }
             }
           if ( n_valid_nbors == 0 ) {
                if ( first_call ) { cerr << check3; first_call=false; }
                cerr <<"\n"<< celltype <<" "<< parseFiniteElementType((*first).FE_Type()) <<":"<< (*first).Idx() <<": has no neighbors.";
                issues++;
             }
         first++;
      }
   // TODO: check the missing neighbors against position of elements and write results to file


    // 6. checking that all non-null neighbors of the cells are valid
    // --------------------------------------------------------------
    cerr <<"\nintegrityCheck: Are all non-null nodes & cell neighbors of the cell valid?\n";
    first = copy_of_first;
   
    while ( first != last ) {
           // printing message before potentially catastrophic failure occurs
           cerr <<"\t"<< parseAbbreviated_FE_Type( (*first).FE_Type() ) <<":"<< (*first).Idx() <<" ("<< celltype <<"), barycenter: "<< (*first).BaryCenter() <<", node flags: ";
           for ( auto i{0U}; i<(*first).Nodes(); ++i )
             cerr <<" "<< parseBoundary((*first).N(i)->AtBoundary());
           cerr << endl;
           // valid cell neighbors should not be corrupt
           for ( auto i{0U}; i<(*first).Neighbors(); ++i ) {
               if ( (*first).Neighbor(i) != nullptr ) {
                    if ( !(*first).FE() ) cerr <<"\nelement "<< (*first).Idx() <<" has corrupt FE pointer.";
                    if ( (*first).Neighbor(i)->Idx() >= max_cell_idx )
                      cerr <<"\nis element "<< (*first).Idx() <<" neighbor idx="<< (*first).Neighbor(i)->Idx() <<" really this large?";
                 }
             }
         first++;
      }
      
   if ( issues > 0 ) return false;
   return true;
 
 } // end integrityCheck
 
template bool integrityCheck<3,Element>( typename plf::colony<Element<3>>::const_iterator, typename plf::colony<Element<3>>::const_iterator );
template bool integrityCheck<2,Element>( typename plf::colony<Element<2>>::const_iterator, typename plf::colony<Element<2>>::const_iterator );
template bool integrityCheck<1,Element>( typename plf::colony<Element<1>>::const_iterator, typename plf::colony<Element<1>>::const_iterator );

template bool integrityCheck<3,Face>( typename plf::colony<Face<3>>::const_iterator, typename plf::colony<Face<3>>::const_iterator );
template bool integrityCheck<2,Face>( typename plf::colony<Face<2>>::const_iterator, typename plf::colony<Face<2>>::const_iterator );
template bool integrityCheck<1,Face>( typename plf::colony<Face<1>>::const_iterator, typename plf::colony<Face<1>>::const_iterator );

template bool integrityCheck<3,InterFace>( typename plf::colony<InterFace<3>>::const_iterator, typename plf::colony<InterFace<3>>::const_iterator );
template bool integrityCheck<2,InterFace>( typename plf::colony<InterFace<2>>::const_iterator, typename plf::colony<InterFace<2>>::const_iterator );
template bool integrityCheck<1,InterFace>( typename plf::colony<InterFace<1>>::const_iterator, typename plf::colony<InterFace<1>>::const_iterator );
                                          
/* OLD METHOD

    while ( first != last ) {
         // FE policy
         if ( (*first).FE() == nullptr ) {
              cerr <<"\n"<< celltype << (*first).Idx() <<": FE pointer corrupt.";
              issues++;
           }
         else {
             // connected nodes
             for ( auto i{0U}; i<(*first).Nodes(); ++i )
               if ( (*first).N(i) == nullptr ) {
                    cerr <<"\n"<< celltype << (*first).Idx() <<": node: "<< i <<": node pointer corrupt.";
                    issues++;
                 }
               else shared_nodes.push_back( (*first).N(i) );
             // there should be at least one neighbor
             size_t n_valid_nbors{0};
             for ( auto i{0U}; i<(*first).Neighbors(); ++i )
               if ( (*first).Neighbor(i) != nullptr )
                 n_valid_nbors++;
             if ( n_valid_nbors == 0 ) {
                  cerr <<"\n"<< celltype <<" "<< parseFiniteElementType((*first).FE_Type()) <<":"<< (*first).Idx() <<": has no neighbors.";
                  issues++;
               }
             // valid neighbors should not be corrupt
             const long one_billion{1000000000};
             cerr <<"\n"<< parseAbbreviated_FE_Type( (*first).FE_Type() ) <<":"<< (*first).Idx() <<" "<< (*first).BaryCenter();
             for ( auto i{0U}; i<(*first).Nodes(); ++i )
               cerr <<" "<< parseBoundary((*first).N(i)->AtBoundary());
             for ( auto i{0U}; i<(*first).Neighbors(); ++i ) {
                 if ( (*first).Neighbor(i) != nullptr ) {
                      if ( !(*first).FE() ) cerr <<"\nelement "<< (*first).Idx() <<" has corrupt FE pointer.";
                      if ( (*first).Neighbor(i)->Idx() > one_billion )
                        cerr <<"\nis element "<< (*first).Idx() <<" neighbor idx="<< (*first).Neighbor(i)->Idx() <<" really this large?";
                   }
               }
           }
         first++;
      }
*/

template<uint32_t dim, template<uint32_t> class CELL>
bool integrityCheck( const plf::colony<CELL<dim> >& cells,
                     typename std::vector<CELL<dim>*>::const_iterator first,
                     typename std::vector<CELL<dim>*>::const_iterator last )
 {
    string  celltype("Element");
    if constexpr ( is_same< CELL<dim>,Face<dim> >::value ) celltype = "Face";
    if constexpr ( is_same< CELL<dim>,InterFace<dim> >::value ) celltype = "InterFace";
    size_t issues{0};

    // --------------------------------------------------------------
    cerr <<"\nintegrityCheck: Are all non-null nodes & cell neighbors of the cell valid?\n";
   
    while ( first != last ) {
           auto it = (*cells.get_iterator(*first));
           // printing message before potentially catastrophic failure occurs
           cerr <<"\t"<< parseAbbreviated_FE_Type( it.FE_Type() ) <<":"<< it.Idx() <<" ("<< celltype <<"), barycenter: "<< it.BaryCenter() <<", node flags: ";
           for ( auto i{0U}; i<it.Nodes(); ++i )
             cerr <<" "<< parseBoundary( it.N(i)->AtBoundary() );
           cerr << endl;
           // valid cell neighbors should not be corrupt
           for ( auto i{0U}; i<it.Neighbors(); ++i ) {
               if ( it.Neighbor(i) != nullptr ) {
                    if ( !it.FE() ) { cerr <<"\nelement "<< it.Idx() <<" has corrupt FE pointer."; issues++; }
                    if ( it.Neighbor(i)->Idx() >= 1e6 ) {
                         cerr <<"\nis element "<< it.Idx() <<" neighbor idx="<< it.Neighbor(i)->Idx() <<" really this large?";
                         issues++;
                      }
                 }
             }
         first++;
      }
      
   return !( issues > 0 );
    
 } // end

template bool integrityCheck<3,Element>( const plf::colony<Element<3>>&, vector<Element<3>*>::const_iterator, vector<Element<3>*>::const_iterator );
template bool integrityCheck<2,Element>( const plf::colony<Element<2>>&, vector<Element<2>*>::const_iterator, vector<Element<2>*>::const_iterator );
template bool integrityCheck<1,Element>( const plf::colony<Element<1>>&, vector<Element<1>*>::const_iterator, vector<Element<1>*>::const_iterator );


/**
      Checks whether volume elements are indeed connected with volume ones, surface ones with surface elements, and line elements with line ones.
      Also checks that elements have a minimum number of neighbors:  line,triangle, tetra->1, quadrilaterals,prism->2, hexa->3, pyramid->1
*/
template<uint32_t dim>
size_t connectivityCheck( typename std::vector<Element<dim>*>::const_iterator first,
                          typename std::vector<Element<dim>*>::const_iterator last )
  {
     size_t issues{0};
     while( first != last ) {
          const CSMP_FEM_TYPE etype = (*first)->FE_Type();
          int n_connected_neighbors{0};
          
          // 1. checking that elements have equivalent types as neighbors
          for ( auto i{0U}; i<(*first)->Neighbors(); ++i )
            if ( (*first)->Neighbor(i) ) {
                 if constexpr ( dim == 3 ) {
                      if ( (*first)->IsVolumeElement() && !(*first)->Neighbor(i)->IsVolumeElement() ) {
                           cerr <<"\nconnectivityCheck: volume Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           issues++;
                        }
                      if ( (*first)->IsSurfaceElement() && !(*first)->Neighbor(i)->IsSurfaceElement() ) {
                           cerr <<"\nconnectivityCheck: surface Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           issues++;
                        }
                      if ( (*first)->IsLineElement() && !(*first)->Neighbor(i)->IsLineElement() ) {
                           cerr <<"\nconnectivityCheck: line Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           issues++;
                        }
                   }
                 if constexpr ( dim == 2 ) {
                      if ( (*first)->IsSurfaceElement() && !(*first)->Neighbor(i)->IsSurfaceElement() ) {
                           cerr <<"\nconnectivityCheck: surface Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a line element!";
                           issues++;
                        }
                   }
                 // n such issues in 1-dimensional model
                 n_connected_neighbors++;
              }
              
          // 2. checking that the elements have the minimum plausible number of neighbors
          if ( n_connected_neighbors == 0 ) { // no neighbors
               cerr <<"\nconnectivityCheck: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
               cerr <<" has "<< n_connected_neighbors <<" neighbor(s).";
               issues++;
            }
          else if ( n_connected_neighbors == 1 ) {
                // minimum number of neighbors is 1
                if ( !isLineElement(etype) &&
                     !isTriangular(etype)  &&
                     !isTetrahedral(etype) ) {
                    cerr <<"\nconnectivityCheck: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
                    cerr <<" has only "<< n_connected_neighbors <<" neighbor(s).";
                    issues++;
                  }
            }
          
          if ( isQuadrilateral(etype) || isPrism(etype) ) {
                if ( n_connected_neighbors < 2 ) {
                    cerr <<"\nconnectivityCheck: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
                    cerr <<" has only "<< n_connected_neighbors <<" neighbor(s).";
                    issues++;
                  }
            }
          
          if ( isHexahedral(etype) || isPyramid(etype) ) {
                if ( n_connected_neighbors < 3 ) {
                    cerr <<"\nconnectivityCheck: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
                    cerr <<" has only "<< n_connected_neighbors <<" neighbor(s).";
                    issues++;
                  }
              }
          
          first++;
       }
       
     return issues;
       
  } // end connectivityCheck

template size_t connectivityCheck<3>( vector<Element<3>*>::const_iterator, vector<Element<3>*>::const_iterator );
template size_t connectivityCheck<2>( vector<Element<2>*>::const_iterator, vector<Element<2>*>::const_iterator );
template size_t connectivityCheck<1>( vector<Element<1>*>::const_iterator, vector<Element<1>*>::const_iterator );





                                          
template<uint32_t dim>
pair<Point<dim>,Point<dim>>  boundingBox( typename vector<Node<dim>*>::const_iterator first, typename vector<Node<dim>*>::const_iterator last )
 {
    if ( first == last ) return make_pair( Point<dim>{}, Point<dim>{} );
    
    Point<dim> pmin((*first)->Coordinate()), pmax((*first)->Coordinate());
    first++;
    
    while ( first != last ) {
         pmin = min( pmin, (*first)->Coordinate() );
         pmax = max( pmax, (*first)->Coordinate() );
         first++;
      }
    
    return make_pair( pmin, pmax );
    
 } // end boundingBox

template pair<Point<3>,Point<3>>  boundingBox( vector<Node<3>*>::const_iterator, vector<Node<3>*>::const_iterator );
template pair<Point<2>,Point<2>>  boundingBox( vector<Node<2>*>::const_iterator, vector<Node<2>*>::const_iterator );
template pair<Point<1>,Point<1>>  boundingBox( vector<Node<1>*>::const_iterator, vector<Node<1>*>::const_iterator );


/* CLIPPING OF FACE NUMBERING FUNCTION FOR INTERFACE

    const size_t neighbors2x( f->Neighbors() * 2 );
    for ( size_t j{0U}; j<neighbors2x; ++j ) {
      InterFace<dim>* const ptr( f->Neighbor( j ) );
      if ( ptr != nullptr ) {
        // building search maps that we will use to find the shared interfaces
        // key=pointset   face iD
        map<set<Point<dim> >, pair<INTERFACE_SIDE, size_t> >   inner_elmt_faces, outer_elmt_faces;
        vector<uint32_t>  nids;
        // first element
        Element<dim>* e1 = f->InnerParent();
        for ( size_t face = 0U; face<e1->Faces(); ++face ) {
          e1->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j{0U}; j<nids.size(); ++j )
            face_key.insert( e1->N( nids[j] )->Coordinate() );
          outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
        }
        // second element
        Element<dim>* e2 = f->OuterParent();
        for ( size_t face = 0U; face<e2->Faces(); ++face ) {
          e2->FE()->NodesOfFace( face, nids );
          set<Point<dim> >  face_key;
          for ( size_t j{0U}; j<nids.size(); ++j )
            face_key.insert( e2->N( nids[j] )->Coordinate() );
          inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
        }

        // 2. finding the shared faces
        bool found( false );
        int64_t  inner_face_id( -1 ), outer_face_id( -1 );
        for ( auto& inner_face : inner_elmt_faces ) {
          for ( auto& outer_face : outer_elmt_faces ) {
            if ( inner_face.first == outer_face.first ) {
              inner_face_id = inner_face.second.second;
              outer_face_id = outer_face.second.second;
              found = true;
              break;
            }
          }
          if ( found ) break;
        }

        vset.Pfvert( eidx, j, static_cast<int32_t>(ptr->Idx()) );
      }
      else
        vset.Pfvert( eidx, j, REGION_BOUNDARY );
    }

*/ // FACE NUMBERING





/**
   returns true if the elements contain each others barycentre,
      only considers equidimensional elements.
*/
template<uint32_t dim>
bool interPenetrating( const Element<dim>* const elmt1, const Element<dim>* const elmt2 )
 {
    if constexpr ( dim == 3 )
      if ( !elmt1->IsVolumeElement() || !elmt2->IsVolumeElement() ) {
           cerr <<"\ninterPenetrating<3>: only works for equidimensional (volumetric) elements.\n";
           return false;
        }

    if constexpr ( dim == 2 )
      if ( !elmt1->IsSurfaceElement() || !elmt2->IsSurfaceElement() ) {
           cerr <<"\ninterPenetrating<2>: only works for equidimensional (surface) elements.\n";
        }

    if constexpr ( dim == 1 ) {
         Point<dim> bctr = elmt1->BaryCenter();
         if ( bctr > elmt2->N(0)->Coordinate() &&
              bctr < elmt2->N(1)->Coordinate() ) return true;
      }

    // is the barycentre of element 1 contained in element 2 (all intpol functions positive?)
    Point<dim> bctr = elmt1->BaryCenter();
    elmt2->N_AtGlobalPoint( elmt2->FE()->NRST, bctr.Coordinates() );
    bool is_contained{true};
    for ( const auto& ni : elmt2->FE()->NRST )
      if ( ni < 0. ) {
          is_contained = false;
          break;
        }
    if ( is_contained ) return true;
    
    // is the barycentre of element 2 contained in element 1 (all intpol functions positive?)
    bctr = elmt2->BaryCenter();
    elmt1->N_AtGlobalPoint( elmt1->FE()->NRST, bctr.Coordinates() );
    is_contained = true;
    for ( const auto& ni : elmt2->FE()->NRST )
      if ( ni < 0. ) {
          is_contained = false;
          break;
        }
    if ( is_contained ) return true;

    return false;

 } // end interPenetrating

template bool interPenetrating( const Element<3>* const, const Element<3>* const );
template bool interPenetrating( const Element<2>* const, const Element<2>* const );
template bool interPenetrating( const Element<1>* const, const Element<1>* const );


/**
@brief Find element which contains a given point.

Note that this only searches volumetric elements. It is not recommended
that you use this function if you need to search for many points. It also
may not work if any elements are concave (possible in the case of hexahedra).

@author  A.J. Bromage
@date    11/04/2018

@param [in] region  region to search
@param [in] query   query point

@return  the element which contains the point, or NULL if no element does

*/
Element<3u>* const pointInVolumeElement( Region<3u>& region, const Point<3u>& query )
    {
      std::vector<uint32_t> fnids;
      fnids.reserve(4);

      const auto eend = region.CellsEnd();
      for (auto eit = region.CellsBegin(); eit != eend; ++eit) {

        // 1. Volume elements only
        
        if (!(*eit)->IsVolumeElement()) continue;

        // 2. Test against axis-aligned bounding box

        double minx = +std::numeric_limits<double>::max();
        double miny = +std::numeric_limits<double>::max();
        double minz = +std::numeric_limits<double>::max();
        double maxx = -std::numeric_limits<double>::max();
        double maxy = -std::numeric_limits<double>::max();
        double maxz = -std::numeric_limits<double>::max();
        const auto iNrNodes = (*eit)->Nodes();
        for ( auto iNode = 0; iNode < iNrNodes; ++iNode ) {
          auto n = (*eit)->N(iNode);
          minx = std::min(minx, n->x());
          maxx = std::max(maxx, n->x());
          miny = std::min(miny, n->y());
          maxy = std::max(maxy, n->y());
          minz = std::min(minz, n->z());
          maxz = std::max(maxz, n->z());
        }

        if (query[0] < minx || query[0] > maxx
            || query[1] < miny || query[1] > maxy
            || query[2] < minz || query[2] > maxz) {
          continue;
        }

        // 3. Test against all faces

        bool reject = false;
        auto fe = (*eit)->FE();
        const auto iNrFaces = (*eit)->Faces();
        for ( auto iFace = 0; iFace < iNrFaces && !reject; ++iFace ) {
          fe->NodesOfFace(iFace, fnids);
          const size_t iNrFacePts = fnids.size();
          for (size_t iFacePt = 0; iFacePt < iNrFacePts; iFacePt += 2) {
            auto p0 = (*eit)->N(fnids[(iFacePt+0) % iNrFacePts])->Coordinate();
            auto p1 = (*eit)->N(fnids[(iFacePt+1) % iNrFacePts])->Coordinate();
            auto p2 = (*eit)->N(fnids[(iFacePt+2) % iNrFacePts])->Coordinate();

            auto normal = crossProduct(p2-p0, p1-p0);
            normal.NormalizeLengthTo(1.0f);
            const double queryDotNormal = dotProduct(query, normal);
            const double p0DotNormal = dotProduct(p0, normal);

            if (queryDotNormal < p0DotNormal) {
              reject = true;
              break;
            }
          }
        }
        if (!reject) {
          return *eit;
        }
        fnids.clear();
      }

      return 0;
    
  } // end pointInVolumeElement




template<uint32_t dim>
size_t collocatedNodes( const Element<dim>* const eptr )
 {
    assert( eptr != nullptr );
    const size_t n_nodes{ eptr->Nodes() };
    set<Point<dim> > node_points;
    
    // set admits only unique node coordinates
    for ( auto i{0U}; i<n_nodes; ++i )
      node_points.insert( eptr->N(i)->Coordinate() );
 
    return n_nodes - node_points.size();
    
 } // end collocatedNodes

template size_t collocatedNodes( const Element<1>* const );
template size_t collocatedNodes( const Element<2>* const );
template size_t collocatedNodes( const Element<3>* const );





/// pretty prints line elements as a chain from beginning to end
template<uint32_t dim>
size_t printLineElementRegion( const Model<dim>& model, const char* region_name, bool renumber_nodes )
 {
    const Region<dim>& line_domain( model.Region(region_name) );
    if ( renumber_nodes ) line_domain.RenumberNodes();
    
    const Element<dim>* eptr1 = (*line_domain.PerimeterCellsBegin());
    const Element<dim>* eptr2 = (*prev(line_domain.CellsEnd(),1));
    
    const Element<dim>* previous_ptr{nullptr};
    size_t              traversed_elmts{0U};
    
    bool forward = ( eptr1->Neighbor(0) != nullptr ) ? true : false;
    
    cout <<"\n\nprintLineElementRegion: '"<< region_name <<"': printing chain of elements ";
    if ( forward ) {
         cout <<"in forward direction:"<< endl;
         while( eptr1 != nullptr ) {
              if ( !eptr1->IsLineElement() )
                throw csmp::Exception( ERROR, "printLineElementRegion", "current element is not a line element; aborting printing" );
              // Ideally (where the numbers are nodes and the labels are BOX_BOUNDARY flags)
              // we should get something like: TOP 1--0 0--11 11--12...56--56 BOTTOM
              if ( !eptr1->Neighbor(1) )
                cout <<"  "<< parseBoundary( eptr1->N(0)->AtBoundary() ) <<" "<< eptr1->N(0)->Idx();
              else cout << eptr1->N(0)->Idx();
              cout <<"->-";
              if ( !eptr1->Neighbor(0) )
                cout << eptr1->N(1)->Idx() <<" "<< parseBoundary( eptr1->N(1)->AtBoundary() ) <<" ";
              else cout << eptr1->N(1)->Idx();
              cout <<" ";
              // advancing along the polyline
              previous_ptr = eptr1;
              eptr1 = eptr1->Neighbor(0);
              // if progress stops
              if ( previous_ptr == eptr1 ) {
                   cout <<"\n| broken chain; ending after "<< traversed_elmts;
                   cout <<" elements vs. "<< line_domain.Cells() <<" in total."<< endl;
                   break;
                }
              traversed_elmts++;
           }
      }
    // backward
    else {
         cout <<"from back to front:"<< endl;
         while( eptr2 != nullptr ) {
              if ( !eptr2->IsLineElement() )
                throw csmp::Exception( ERROR, "printLineElementRegion", "current element is not a line element; aborting printing" );
              if ( !eptr2->Neighbor(0) )
                cout <<"  "<< parseBoundary( eptr2->N(1)->AtBoundary() ) <<" "<< eptr2->N(1)->Idx();
              else cout << eptr2->N(1)->Idx();
              cout <<"-<-";
              if ( !eptr2->Neighbor(1) )
                cout << eptr1->N(0)->Idx() <<" "<< parseBoundary( eptr2->N(0)->AtBoundary() ) <<" ";
              else cout << eptr2->N(0)->Idx();
              cout <<" ";
              // advancing along the polyline
              previous_ptr = eptr2;
              eptr2 = eptr2->Neighbor(1);
              // if progress stops
              if ( previous_ptr == eptr2 ) {
                   cout <<"\n| broken chain; ending after "<< traversed_elmts;
                   cout <<" elements vs. "<< line_domain.Cells() <<" in total."<< endl;
                   break;
                }
              traversed_elmts++;
           }
      }
      
   // checking
   if ( traversed_elmts < line_domain.Cells() ) {
        cout <<"\n"<<"printLineElementRegion('"<< region_name <<"'): only "<< traversed_elmts <<" of "<< line_domain.Cells();
        cout <<" total were discovered by neighbor to neighbor traversal. Broken connectivity?"<< endl;
     }
  
   return traversed_elmts;
    
 } // end printLineElementRegion

template size_t printLineElementRegion( const Model<3U>&, const char*, bool );
template size_t printLineElementRegion( const Model<2U>&, const char*, bool );




/**
    Prints coordinates of range of nodes nodes in a spread-sheet plottable format
*/
template<uint32_t dim>
void printNodeCoordinates( typename vector<Node<dim>*>::const_iterator first,
                           typename vector<Node<dim>*>::const_iterator last )
 {
     size_t n_nodes = distance(first,last);
     if ( n_nodes == 0 ) {
          cerr <<"\nprintNodeCoordinates: supplied iterator range is empty, nothing could be printed."<< endl;
          return;
       }
       
     cout <<"\n"<<"printNodeCoordinates: printing range of "<< n_nodes <<" nodes. Idx followed by coordinate values:"<< endl;
     cout <<"\n"<<"idx, x";
     if constexpr ( dim == 2U ) cout <<", y"<< endl;
     if constexpr ( dim == 3U ) cout <<", y, z"<< endl;
     
     while ( first != last ) {
          assert( (*first) != nullptr );
          cout << (*first)->Idx() <<", "<< (*first)->x();
          for ( auto i{1U}; i<dim; i++ ) cout <<", "<< (*(*first))[i];
          cout << endl;
          first++;
       }
     
 } // end printNodeCoordinates

template void printNodeCoordinates<3U>( typename vector<Node<3U>*>::const_iterator, typename vector<Node<3U>*>::const_iterator );
template void printNodeCoordinates<2U>( typename vector<Node<2U>*>::const_iterator, typename vector<Node<2U>*>::const_iterator );
template void printNodeCoordinates<1U>( typename vector<Node<1U>*>::const_iterator, typename vector<Node<1U>*>::const_iterator );

} // end csmp
