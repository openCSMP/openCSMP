//
//  MeshManagementUtilities.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 3/7/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include <unordered_set>
#include "meshManagementUtilities.h"
#include "MeshManager.h"
#include "MeshPatch.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "Box.h"
#include "plf_colony.h"

using namespace std;

namespace csmp {


/**
     @note Model is used for convenience here. MeshManager would be mode appropriate.
*/
template<uint32_t dim>
size_t currentCellTypes( const MeshManager<dim>& mesh, PLACEMENT etype, size_t& volume_cells, size_t& surface_cells, size_t& line_cells )
 {
    size_t n_cells_model{ mesh.Elements() + mesh.Faces() + mesh.Interfaces() };
    volume_cells = surface_cells = line_cells = 0U;
    
    // elements
    if ( etype == ELEMENT ) {
         if constexpr ( dim == 3U ) {
             for ( auto it=mesh.ElementsBegin(); it!=mesh.ElementsEnd(); ++it ) {
                  if      ( (*it).IsVolume() ) volume_cells++;
                  else if ( (*it).IsSurface() ) surface_cells++;
                  else { assert( (*it).IsLine() ); line_cells++; }
               }
           }
         else if constexpr ( dim == 2U ) {
             for ( auto it=mesh.ElementsBegin(); it!=mesh.ElementsEnd(); ++it ) {
                  if ( (*it).IsSurface() ) surface_cells++;
                  else { assert( (*it).IsLine() ); line_cells++; }
               }
           }
         else // 1D
           line_cells = mesh.Elements();
      }
    // faces
    else if ( etype == FACE ) {
         surface_cells = count_if( mesh.FacesBegin(), mesh.FacesEnd(), []( const Face<dim>& e ){ return e.IsSurface(); } );
         line_cells = mesh.Faces() - surface_cells;
      }
    // interfaces
    else if ( etype == INTER_FACE ) {
         surface_cells = count_if( mesh.InterfacesBegin(), mesh.InterfacesEnd(), []( const InterFace<dim>& e ){ return e.IsSurface(); } );
         line_cells = mesh.Interfaces() - surface_cells;
      }
    else cerr <<"\n\n"<< "currentCellTypes: invalid celltype: "<< parsePlacement(etype) << endl;
 
   return n_cells_model;
 
 } // end currentCellTypes

template size_t currentCellTypes( const MeshManager<3U>&, PLACEMENT, size_t&, size_t&, size_t& );
template size_t currentCellTypes( const MeshManager<2U>&, PLACEMENT, size_t&, size_t&, size_t& );
template size_t currentCellTypes( const MeshManager<1U>&, PLACEMENT, size_t&, size_t&, size_t& );





template<uint32_t dim>
bool isoparametricElementMesh( const Model<dim>& sg )
 {
    string  etype(parseFiniteElementType((*sg.Region("Model").CellsBegin())->FE()->ElementType()));

    if ( etype.find("ISOPARAMETRIC") != std::string::npos ) return true;
    return false;

 } // end isoparametricElementMesh
 
template bool isoparametricElementMesh<1U>( const Model<1U>& );
template bool isoparametricElementMesh<2U>( const Model<2U>& );
template bool isoparametricElementMesh<3U>( const Model<3U>& );





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






/// finds cells that have the same nodes and reports their numbers
template<uint32_t dim, template<uint32_t> class CELL>
size_t detectDuplicateCells( typename vector<CELL<dim>*>::const_iterator first,
                             typename vector<CELL<dim>*>::const_iterator last,
                             bool verbose )
 {
    map<set<Node<dim>*>,set<CELL<dim>*> > potential_duplicates;
    size_t  n_duplicates{0ul};
    size_t  cell_counter{0ul};
    bool    first_issue{ true };
    
    while( first != last ) {
         if ( (*first) == nullptr ) {
               if ( first_issue ) cout <<"\n"<<"detectDuplicateCells: cell "<< cell_counter <<"=nullptr";
               else cout <<" cell "<< cell_counter <<"=nullptr "<< endl;
               first_issue = false;
               continue;
           }
         // creating cell keys from their node pointers
         set<Node<dim>*> node_set;
         const auto n_nodes{ (*first)->Nodes() };
         if constexpr( is_same<CELL<dim>,InterFace<dim>>::value ) {
              for ( uint32_t i{0U}; i<n_nodes/2u; i++ ) {
                   if ( (*first)->N(i) == nullptr ) {
                        if ( first_issue ) cout <<"\n"<<"detectDuplicateCells: cell "<< cell_counter <<", node"<< i <<"=nullptr";
                        else cout <<" cell "<< cell_counter <<", node"<< i <<"=nullptr";
                        first_issue = false;
                        continue;
                     }
                   node_set.insert( (*first)->N(i) );
                }
           }
         else if constexpr( is_same<CELL<dim>,Element<dim>>::value ||
                            is_same<CELL<dim>,Face<dim>>::value )  {
              for ( uint32_t i{0U}; i<n_nodes; i++ ) {
                   if ( (*first)->N(i) == nullptr ) {
                        if ( first_issue ) cout <<"\n"<<"detectDuplicateCells: cell "<< cell_counter <<", node"<< i <<"=nullptr";
                        else cout <<" cell "<< cell_counter <<", node"<< i <<"=nullptr";
                        first_issue = false;
                        continue;
                     }
                   node_set.insert( (*first)->N(i) );
                }
           }
         // recording the cells
         auto it = potential_duplicates.insert( make_pair( node_set, set<CELL<dim>*>{(*first)} ) );
         // if there is a cell with the same nodes but a different pointer, it is recorded
         if ( it.second == false ) {
              auto cit =(*it.first).second.insert( (*first) );
              if ( verbose && cit.second == false ) {
                  cout <<"\n"<<"detectDuplicateCells: input range contains multiple copies of:";
                  (*first)->Out();
                }
              n_duplicates++;
           }
         first++;
         cell_counter++;
      }
      
    // printing the duplicate cells if any
    if ( n_duplicates  > 0U && verbose ) {
         cout <<"\n\n"<<"detectDuplicateCells: found "<< n_duplicates <<" cells sharing all nodes in input range:";
         for ( auto pd : potential_duplicates )
           if ( pd.second.size() > 1U )
             (*pd.second.begin())->Out();
      }
  
    return n_duplicates;
    
 } // end detectDuplicateCells

template size_t detectDuplicateCells<3,Element>( vector<Element<3>*>::const_iterator, vector<Element<3>*>::const_iterator, bool );
template size_t detectDuplicateCells<2,Element>( vector<Element<2>*>::const_iterator, vector<Element<2>*>::const_iterator, bool );
template size_t detectDuplicateCells<1,Element>( vector<Element<1>*>::const_iterator, vector<Element<1>*>::const_iterator, bool );

template size_t detectDuplicateCells<3,Face>( vector<Face<3>*>::const_iterator, vector<Face<3>*>::const_iterator, bool );
template size_t detectDuplicateCells<2,Face>( vector<Face<2>*>::const_iterator, vector<Face<2>*>::const_iterator, bool );
template size_t detectDuplicateCells<1,Face>( vector<Face<1>*>::const_iterator, vector<Face<1>*>::const_iterator, bool );

template size_t detectDuplicateCells<3,InterFace>( vector<InterFace<3>*>::const_iterator, vector<InterFace<3>*>::const_iterator, bool );
template size_t detectDuplicateCells<2,InterFace>( vector<InterFace<2>*>::const_iterator, vector<InterFace<2>*>::const_iterator, bool );
template size_t detectDuplicateCells<1,InterFace>( vector<InterFace<1>*>::const_iterator, vector<InterFace<1>*>::const_iterator, bool );





/// finds nodes that are not connected to any elements, faces or interfaces. If there are any, it returns their number and pointers to them into the argument set
template<uint32_t dim>
size_t detectOrphanNodes( MeshManager<dim>& mesh, set<Node<dim>*>& orphan_nodes )
 {
    set<Node<dim>*> node_ptrs;
    
    // loop over elements, faces and interfaces, returning storing pointers to their nodes in a vector
    for ( auto it=mesh.ElementsBegin(); it!=mesh.ElementsEnd(); ++it ) {
         const auto n_nodes{ (*it).Nodes() };
         for ( uint32_t i{0U}; i<n_nodes; i++ )
           node_ptrs.insert( (*it).N(i) );
      }
    if ( mesh.Faces() > 0 )
      for ( auto it=mesh.FacesBegin(); it!=mesh.FacesEnd(); ++it ) {
           const auto n_nodes{ (*it).Nodes() };
           for ( uint32_t i{0U}; i<n_nodes; i++ )
             node_ptrs.insert( (*it).N(i) );
        }
    if ( mesh.Interfaces() > 0 )
      for ( auto it=mesh.InterfacesBegin(); it!=mesh.InterfacesEnd(); ++it ) {
           const auto n_nodes{ (*it).Nodes() };
           for ( uint32_t i{0U}; i<n_nodes; i++ )
             node_ptrs.insert( (*it).N(i) );
        }
    
    if ( node_ptrs.size() < mesh.Nodes() ) {
         orphan_nodes.clear();
         // aua! - checking which of the nodes in the node vector have no connections to elements, faces or interfaces
         for ( auto nit=mesh.NodesBegin(); nit!=mesh.NodesEnd(); ++nit )
           if ( node_ptrs.find( &(*nit) ) == node_ptrs.end() )
             orphan_nodes.insert( &(*nit) );
    
         return orphan_nodes.size();
      }
    
    return 0U;
    
 } // end detectOrphanNodes

template size_t detectOrphanNodes( MeshManager<3U>&, set<Node<3U>*>& );
template size_t detectOrphanNodes( MeshManager<2U>&, set<Node<2U>*>& );
template size_t detectOrphanNodes( MeshManager<1U>&, set<Node<1U>*>& );





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
         csmp_error.Note( ERROR, "findContiguousMeshPatch", "entry cell pointer is a nullptr; nothing was done.");
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
                 for ( uint32_t j{0U}; j<n_neighbors; ++j )
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





/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these, a pointer is inserted into the argument vector (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
             @author SKM 14/8/21
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  findContiguousMeshPatches( typename plf::colony<CELL<dim>>::iterator begin,
                                   typename plf::colony<CELL<dim>>::iterator end,
                                   map<string,vector<CELL<dim>*> >& mesh_patches )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( begin == end ) {
             csmp_error.Note( WARNING, "findContiguousMeshPatches:",
                                       "input CELL pointer range is empty." );
            return 0U;
        }
     
      // getting a copy of the element pointers of the mesh (without writing an adaptor from colony iterator to cell pointer)
      vector<CELL<dim>*>   cells;
      cells.reserve( static_cast<size_t>(distance(begin,end)) );
      while ( begin != end ) {
          cells.push_back( &(*begin) );
          ++begin;
        }

      // detecting via a flood-fill whether the group can be partitioned, else nothing is done
      set<CELL<dim>*>  cells_contiguous_subset;
      findContiguousMeshPatch<dim>( (*cells.begin()), cells_contiguous_subset );

      // if no contiguous cells could be found
      if ( cells_contiguous_subset.empty() ) {
           csmp_error.Note( WARNING, "findContiguousMeshPatches:",
                                     "No contiguous cells found. Is the neighbor connectivity missing? - nothing was done" );
           return 0U;
        }

      // if the first flood-fill reached all elements of the region or more on the outside it is contiguous
      if ( cells.size() <= cells_contiguous_subset.size() ) {
           std::cout <<"\n"<<"findContiguousMeshPatches: ";
           std::cout <<"mesh is already contiguous, nothing was done."<< std::endl;
           
           // creating name for contiguous patch from finite-element type
           string patch_name( parseFiniteElementType( (*cells_contiguous_subset.begin())->FE_Type() ) );
           // appending the patch number and the number of elements in patch
           patch_name +="_patch";
           patch_name += to_string( 1 );
           patch_name +="_";
           patch_name += to_string( cells_contiguous_subset.size() );
           patch_name +="cells";
           mesh_patches.insert( make_pair( patch_name, std::move( vector<CELL<dim>*>{ cells_contiguous_subset.begin(),
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
                     cout <<"\n"<<"findContiguousMeshPatches: ";
                     cout <<"mesh is divided into disconnected patch(es):\n";
                  }
                cout <<"\t\t\t'"<< patch_name <<"'";
                cout <<" ("<< cells_contiguous_subset.size() <<" elmts)"<< std::endl;

                pair<typename map<string,vector<CELL<dim>*> >::iterator,bool>
                  insertion = mesh_patches.insert( make_pair( patch_name,
                                                   std::move( vector<CELL<dim>*>( cells_contiguous_subset.begin(),
                                                                             cells_contiguous_subset.end() ) ) ) );
                if ( insertion.second == false ) {
                     csmp_error.Note( ERROR, "findContiguousMeshPatches:", patch_name,
                                               "could not be inserted into patch map" );
                  }
                n_patches++;
              }
           else csmp_error.Note( ERROR, "findContiguousMeshPatches:", "patch contains no elements, nothing was done" );
              
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
    
   } // end findContiguousMeshPatches

// 3D version
template size_t  findContiguousMeshPatches( plf::colony<Element<3U>>::iterator,
                                            plf::colony<Element<3U>>::iterator,
                                            map<string,vector<Element<3U>*> >& );

template size_t  findContiguousMeshPatches( plf::colony<Face<3U>>::iterator,
                                            plf::colony<Face<3U>>::iterator,
                                            map<string,vector<Face<3U>*> >& );

template size_t  findContiguousMeshPatches( plf::colony<InterFace<3U>>::iterator,
                                            plf::colony<InterFace<3U>>::iterator,
                                            map<string,vector<InterFace<3U>*> >& );

// 2D version
template size_t  findContiguousMeshPatches( plf::colony<Element<2U>>::iterator,
                                            plf::colony<Element<2U>>::iterator,
                                            map<string,vector<Element<2U>*> >& );

template size_t  findContiguousMeshPatches( plf::colony<Face<2U>>::iterator,
                                            plf::colony<Face<2U>>::iterator,
                                            map<string,vector<Face<2U>*> >& );

template size_t  findContiguousMeshPatches( plf::colony<InterFace<2U>>::iterator,
                                            plf::colony<InterFace<2U>>::iterator,
                                            map<string,vector<InterFace<2U>*> >& );
                                            
// 1D version
template size_t  findContiguousMeshPatches( plf::colony<Element<1U>>::iterator,
                                            plf::colony<Element<1U>>::iterator,
                                            map<string,vector<Element<1U>*> >& );

template size_t  findContiguousMeshPatches( plf::colony<Face<1U>>::iterator,
                                            plf::colony<Face<1U>>::iterator,
                                            map<string,vector<Face<1U>*> >& );

template size_t  findContiguousMeshPatches( plf::colony<InterFace<1U>>::iterator,
                                            plf::colony<InterFace<1U>>::iterator,
                                            map<string,vector<InterFace<1U>*> >& );





/**
      For the user-defined starting Node, finds the nodes which are connected to it.
      
      @param nptr pointer from which the breadth first search for neighbor nodes is started
      @param contiguous_set_of_nodes the  nodes forming an interconnected cluster
      @return how many nodes are in the interconnected cluster that was found
*/
template<uint32_t dim>
size_t findInterconnectedNodeCluster( Node<dim>* const nptr, std::set<Node<dim>*>& contiguous_set_of_nodes )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    if ( nptr == nullptr ) {
         csmp_error.Note( ERROR, "findInterconnectedNodeCluster", "entry node is a nullpointer; nothing was done");
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
             const auto n_node_nbors{ nit->Neighbors() };
             new_node_nbors.reserve( n_node_nbors );
             for ( uint32_t i{0U}; i < n_node_nbors; ++i ) {
                 Node<dim>* nbor_ptr = nit->Neighbor(i);
                 assert( nbor_ptr != nullptr );
                 // if the neighbor-node is not part of the cluster yet, it is added
                 if ( contiguous_set_of_nodes.find( nbor_ptr ) == contiguous_set_of_nodes.end() ) {
                      contiguous_set_of_nodes.insert( nbor_ptr );
                      // remembering the new nodes for the next search
                      new_node_nbors.push_back( nbor_ptr );
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



template<uint32_t dim>
size_t countCornerNodes(typename plf::colony<Element<dim>>::const_iterator elmts_begin,
                        typename plf::colony<Element<dim>>::const_iterator elmts_end ){

  set<Node<dim>*> corner_nodes;

  for (typename plf::colony<Element<dim>>::const_iterator eit = elmts_begin; eit != elmts_end; eit++ ){
    uint32_t n_corner_nds = eit->FE()->CornerNodes();
    for ( uint32_t n{0U}; n<n_corner_nds; n++ ){
      corner_nodes.insert( eit->N(n) );
    }
  }

  return corner_nodes.size();

}

template size_t countCornerNodes<3>( typename plf::colony<Element<3>>::const_iterator , typename plf::colony<Element<3>>::const_iterator  );
template size_t countCornerNodes<2>( typename plf::colony<Element<2>>::const_iterator , typename plf::colony<Element<2>>::const_iterator  );
template size_t countCornerNodes<1>( typename plf::colony<Element<1>>::const_iterator , typename plf::colony<Element<1>>::const_iterator  );




/**
    Recursive Depth-First Search (DFS) traversal function of mesh based on its node connectivity.
    Node connectivity-based traversal of mesh.
*/
/*
template<uint32_t dim>
void depthFirstSearch( const Node<dim>* node, std::unordered_set<const Node<dim>*>& visited )
 {
    if ( !node || visited.count(node) ) {
        return; // If node is null or already visited, return
      }

    // Mark the current node as visited
    visited.insert(node);

    // Recur for all the neighbors of this node
    for ( auto nit=node->NeighborsBegin(); nit!=node->NeighborsEnd(); ++nit ) {
         depthFirstSearch( (*nit), visited );
      }
      
} // end depthFirstSearch

template void depthFirstSearch<3>( const Node<3>*, std::unordered_set<const Node<3>*>& );
template void depthFirstSearch<2>( const Node<2>*, std::unordered_set<const Node<2>*>& );
template void depthFirstSearch<1>( const Node<1>*, std::unordered_set<const Node<1>*>& );
*/

/**
      Relying on the node-to-node connectivity, discovers all nodes in the supplied mesh.
*/
/*
template<uint32_t dim>
size_t findInterconnectedNodes( const plf::colony<Node<dim>>& nodes ) {
    std::unordered_set<const Node<dim>*> visited; // To keep track of visited nodes
    if (!nodes.empty()) {
         // Start DFS from the first node in the vector
         depthFirstSearch<dim>( (&(*nodes.begin())), visited );
      }
   return nodes.size();
}

template size_t findInterconnectedNodes( const plf::colony<Node<3>>& );
template size_t findInterconnectedNodes( const plf::colony<Node<2>>& );
template size_t findInterconnectedNodes( const plf::colony<Node<1>>& );
*/


// Template function to traverse the tree
template <uint32_t dim>
void depthFirstFilteredSearch( const Node<dim>* node, std::unordered_set<const Node<dim>*>& visited,
                               std::set<std::pair<const Node<dim>*, const Node<dim>*>>& validEdges ) {
    // Mark the current node as visited
    visited.insert(node);

    // Process the node (e.g., print it)
    // process(node);

    // Recur for all the nodes adjacent to this node that are connected by valid edges
    //for (auto& neighbor : node->neighbors) {
    for ( auto neighbor=node->NeighborsBegin(); neighbor!=node->NeighborsEnd(); ++neighbor ) {
        auto edge        = make_pair( node, (*neighbor) );
        //auto reverseEdge = make_pair( (*neighbor), node ); // For undirected trees
        //if (validEdges.count(edge) || validEdges.count(reverseEdge) ) {
        if ( validEdges.count(edge) ) {
            if (visited.find(*neighbor) == visited.end()) {
                depthFirstFilteredSearch( (*neighbor), visited, validEdges );
            }
        }
    }
}

template void depthFirstFilteredSearch<3>( const Node<3>*, unordered_set<const Node<3>*>&, set<pair<const Node<3>*,const Node<3>*>>& );
template void depthFirstFilteredSearch<2>( const Node<2>*, unordered_set<const Node<2>*>&, set<pair<const Node<2>*,const Node<2>*>>& );
template void depthFirstFilteredSearch<1>( const Node<1>*, unordered_set<const Node<1>*>&, set<pair<const Node<1>*,const Node<1>*>>& );



/**
      Relying on the node-to-node connectivity, discovers all nodes in the supplied mesh.
*/
template<uint32_t dim>
size_t findInterconnectedNodes( const plf::colony<Node<dim>>& nodes,
                                set<pair<const Node<dim>*,const Node<dim>*>>& validEdges )
 {
    std::unordered_set<const Node<dim>*> visited; // To keep track of visited nodes
    if (!nodes.empty()) {
         // Start DFS from the first node in the vector
         depthFirstFilteredSearch<dim>( (&(*nodes.begin())), visited, validEdges );
      }
   return nodes.size();
 }

template size_t findInterconnectedNodes( const plf::colony<Node<3>>&, set<pair<const Node<3>*,const Node<3>*>>& );
template size_t findInterconnectedNodes( const plf::colony<Node<2>>&, set<pair<const Node<2>*,const Node<2>*>>& );
template size_t findInterconnectedNodes( const plf::colony<Node<1>>&, set<pair<const Node<1>*,const Node<1>*>>& );




/**
      Traversing mesh to find patches that cannot be reached by neighborhood traversal.
      For each of these a pointer is inserted into the argument deque (root pointer container).
       
      returns number of stand-alone mesh patches found.
      
*/
template<uint32_t dim, template<uint32_t> class CELL>
size_t  findPointersToStandAloneMeshPatches( typename vector<CELL<dim>*>::const_iterator begin,
                                             typename vector<CELL<dim>*>::const_iterator end,
                                             map<CELL<dim>*,MeshPatch<dim>>& root_pointers )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( begin == end ) {
             csmp_error.Note( WARNING, "findPointersToStandAloneMeshPatches:",
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
      size_t       n_subgroups(1);
    
      // creating new contiguous group from the element subset
      while ( !cells.empty() )
        {
           // creating name of contiguous subgroup
           subgroup_name = group_name + to_string( n_subgroups );
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
                MeshPatch<dim> attributes( parseFiniteElementDimension((*cells.begin())->FE_Type()) );
                
                pair<typename map<CELL<dim>*,MeshPatch<dim>>::iterator,bool>
                  insertion = root_pointers.insert( make_pair( (*cells.begin()), attributes ) );
                if ( insertion.second == false ) {
                     csmp_error.Note( WARNING, "findPointersToStandAloneMeshPatches:",
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
                                                      map<Element<3U>*,MeshPatch<3U>>& );

template size_t  findPointersToStandAloneMeshPatches( vector<Face<3U>*>::const_iterator,
                                                      vector<Face<3U>*>::const_iterator,
                                                      map<Face<3U>*,MeshPatch<3U>>& );

template size_t  findPointersToStandAloneMeshPatches( vector<InterFace<3U>*>::const_iterator,
                                                      vector<InterFace<3U>*>::const_iterator,
                                                      map<InterFace<3U>*,MeshPatch<3U>>& );









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
   if ( face == nullptr ) {
        csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(Face)", "pointer to target Face is not initialised");
        return;
     }
   if ( inner_nbor == nullptr ) {
        csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(Face)", "pointer to inner higher-dim Element not initialised");
        return;
     }
   if ( outer_nbor == nullptr ) {
        csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(Face)", "pointer to outer higher-dim Element  not initialised");
        return;
     }
    
   // 1. Creating a map of the faces of the outer element
   //      face key
   set<set<Node<dim>*> > outer_elmt_faces;

   const auto n_outer_elmt_faces(outer_nbor->Faces());
   for ( uint32_t i{0U}; i<n_outer_elmt_faces; ++i ) {
        // creating and recording the search key and face number
        outer_elmt_faces.insert( outer_nbor->CornerNodesOfFace(i) );
     }
     
   // 2. Searching the faces of the inner element that matches this face
   const auto n_inner_elmt_faces(inner_nbor->Faces());
   for ( uint32_t i{0U}; i<n_inner_elmt_faces; ++i ) {
        // creating the search key and performing the search
        auto search_it = outer_elmt_faces.find( inner_nbor->CornerNodesOfFace(i) );
        // if a matching face is found
        if ( search_it != outer_elmt_faces.end() ) {
             // assigning the nodes which are in the right order in fnids
             // (remember that the nodes of the Face should match the order at the inner face)
             uint32_t k(0U);
             for ( const auto& j : inner_nbor->FE()->NodesOfFace(i) )
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
   if ( interface == nullptr ) {
        csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "pointer to target InterFace is not initialised");
        return;
     }
   if ( inner_nbor == nullptr ) {
        csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "pointer to inner higher-dim Element not initialised");
        return;
     }
   if ( outer_nbor == nullptr ) {
        csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "pointer to outer higher-dim Element  not initialised");
        return;
     }
    
   // 1. creating a search map from the nodes of the outer element
   //  key       inner local id, outer local node id
   map<Point<dim>,pair<uint32_t,uint32_t> > outer_elmt_nodes;
   // OUTER ELEMENT
   const auto n_nodes_outer_elmt(outer_nbor->Nodes());
   for ( uint32_t i{0U}; i<n_nodes_outer_elmt; ++i )
     outer_elmt_nodes.insert( make_pair( outer_nbor->N(i)->Coordinate(), make_pair(numeric_limits<uint32_t>::max(),i) ) );
   
   // 2. searching for the shared nodes
   const auto n_nodes_inner_elmt(inner_nbor->Nodes());
   for ( uint32_t i{0U}; i<n_nodes_inner_elmt; ++i ) {
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
   // INNER ELEMENT
   const auto n_inner_elmt_faces(inner_nbor->Faces());
   bool  inner_face_found(false);
   for ( uint32_t i{0U}; i<n_inner_elmt_faces; ++i ) {
        auto fnids = inner_nbor->FE()->NodesOfFace(i);
        // creating and recording the search key
        set<uint32_t> face_nodes( fnids.begin(), fnids.end() );
        // if it matches the interface, the nodes are assigned and the loop is stopped
        if ( face_nodes == inner_nodes ) {
             inner_face_found = true;
             uint32_t k(0U);
             for ( const auto& j : fnids )
               interface->Assign( k++, inner_nbor->N(j), INSIDE );
             interface->ParentFaceID( INSIDE, i );
             break;
          }
     }
   if ( !inner_face_found )
     csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "failed to find nodes of inner higher-dim neighbor element");
    
   // OUTER ELEMENT
   const auto n_outer_elmt_faces(outer_nbor->Faces());
   bool  outer_face_found(false);
   for ( uint32_t i{0U}; i<n_outer_elmt_faces; ++i ) {
        auto fnids = outer_nbor->FE()->NodesOfFace( i );
        // creating and recording the search key
        set<uint32_t> face_nodes( fnids.begin(), fnids.end() );
        // if it matches the interface, the nodes are assigned and the loop is stopped
        if ( face_nodes == outer_nodes ) {
             outer_face_found = true;
             auto k(0U);
             for ( const auto& j : fnids )
               interface->Assign( k++, outer_nbor->N(j), OUTSIDE );
             interface->ParentFaceID( OUTSIDE, i );
             break;
          }
     }
   if ( !outer_face_found )
     csmp_error.Note( ERROR, "findNodesViaHigherDimensionalNeighbors(InterFace)", "failed to find nodes of outer higher-dim neighbor element");

 } // end findNodesViaHigherDimensionalNeighbors

template void findNodesViaHigherDimensionalNeighbors( Element<1>* const, Element<1>* const, InterFace<1>* const );
template void findNodesViaHigherDimensionalNeighbors( Element<2>* const, Element<2>* const, InterFace<2>* const );
template void findNodesViaHigherDimensionalNeighbors( Element<3>* const, Element<3>* const, InterFace<3>* const );




/**
    Finds the face between 2 elements (if any) via the neighbor connectivity.
*/
template<uint32_t dim>
pair<uint32_t,uint32_t> findAdjacentFacesFromNeighbors( Element<dim>* const eptr1, Element<dim>* const eptr2 )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   const uint32_t unspecified = numeric_limits<uint32_t>::max();
   
   if ( eptr1 == nullptr ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNeighbors", "null pointer to first element.");
        return make_pair( unspecified, unspecified );
     }
   if ( eptr2 == nullptr ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNeighbors", "null pointer to second element.");
        return make_pair( unspecified, unspecified );
     }
   if ( eptr1 == eptr2 ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNeighbors", "the supplied pointers point to the same element!");
        return make_pair( unspecified, unspecified );
     }
    
    // 1. finding which face of element 1 is shared with element 2
    uint32_t face_of_elmt1 = unspecified;
    const auto n_faces1(eptr1->Faces());
    for ( uint32_t i{0U}; i<n_faces1; ++i )
      if ( eptr1->Neighbor(i) == eptr2 ) {
           face_of_elmt1 = i;
           break;
        }

   if ( face_of_elmt1 == unspecified ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNeighbors", "could not find matching face of Element 1");
        return make_pair( unspecified, unspecified );
     }

    // 2. finding which face of element 2 is shared with element 1
    uint32_t face_of_elmt2 = unspecified;
    const auto n_faces2(eptr2->Faces());
    for ( uint32_t i{0U}; i<n_faces2; ++i )
      if ( eptr2->Neighbor(i) == eptr1 ) {
           face_of_elmt2 = i;
           break;
        }

   if ( face_of_elmt2 == unspecified ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNeighbors", "could not find matching face of Element 2");
        return make_pair( unspecified, unspecified );
     }

    return make_pair( face_of_elmt1, face_of_elmt2 );
    
 } // end findAdjacentElementFaces
 
template pair<uint32_t,uint32_t> findAdjacentFacesFromNeighbors( Element<3>* const, Element<3>* const );
template pair<uint32_t,uint32_t> findAdjacentFacesFromNeighbors( Element<2>* const, Element<2>* const );
template pair<uint32_t,uint32_t> findAdjacentFacesFromNeighbors( Element<1>* const, Element<1>* const );
 




/**
     Finds the adjacent faces of the supplied elements via their shared nodes.
     
     @return pair of the local face ID numbers of element one and two.

     @attention if neighbor connectivty exists, use matching neighbor pointers which is much faster!
     @attention if no shared face can be found, function returns UNSPECIFIED.
*/
template<uint32_t dim>
pair<uint32_t,uint32_t> findAdjacentFacesFromNodes( Element<dim>* const eptr1, Element<dim>* const eptr2 )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   if ( eptr1 == nullptr ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNodes", "null pointer to first element.");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
   if ( eptr2 == nullptr ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNodes", "null pointer to second element.");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
   if ( eptr1 == eptr2 ) {
        csmp_error.Note( ERROR, "findAdjacentFacesFromNodes", "the supplied pointers point to the same element!");
        return make_pair( UNSPECIFIED, UNSPECIFIED );
     }
    
    // finding the shared face by their nodes
     // 1. creating unique keys from the nodes of the first elements faces not located at a boundary
    vector<set<Node<dim>*> > e1_face_keys;
    const size_t n_faces(eptr1->Faces());
    e1_face_keys.reserve(n_faces);
    for ( uint32_t i{0U}; i<n_faces; ++i )
      if ( eptr1->Neighbor(i) != nullptr ) {
          e1_face_keys.emplace_back( eptr1->CornerNodesOfFace(i) );
       }
      else e1_face_keys.emplace_back( set<Node<dim>*>{} );
    
    // 2. creating face keys for the outer element trying match them
    //    with the faces of the inner one
    set<uint32_t>    face_key_n;
    const size_t n_faces2(eptr2->Faces());
    for ( uint32_t i{0U}; i<n_faces2; ++i ) {
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
 
template pair<uint32_t,uint32_t> findAdjacentFacesFromNodes( Element<3>* const, Element<3>* const );
template pair<uint32_t,uint32_t> findAdjacentFacesFromNodes( Element<2>* const, Element<2>* const );
template pair<uint32_t,uint32_t> findAdjacentFacesFromNodes( Element<1>* const, Element<1>* const );
 






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
         csmp_error.Note( ERROR, "floodFill", "root cell pointer is a nullptr; nothing was done.");
         return;
      }
    // identifying the neighbors of the first element to be looked at
    vector<CELL<dim>*>  neighbor_elements;
    const size_t  neighbors(eptr->Neighbors());
    neighbor_elements.reserve( neighbors );
    for ( uint32_t i{0U}; i<neighbors; i++ )
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
                for ( uint32_t j{0U}; j<n_neighbors; ++j )
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
    Finds elements in region 'subdomain' that contact eachother across split interfaces via matching face nodes.
    Only those elements are discovered that are node-matched.
    
    The neighborhood relations of the elements discovered on either side of the internal boundary are returned in the element pair vector.
 
    @param subdomain (non-unique) region which is anticipated to contain a node-matched split boundary
 
    @return returns false if none of the perimeter elements has nodes that match another perimeter face
 
    @attention method only looks at elements with the same spatial dimension as the  model, i.e. volumes in 3D, surfaces in 2D etc.
    lower-dimensional elements are ignored.
 
    @attention the elements inside the model that are located along the splitboundary
    do not have neighbor pointers yet.
 
    @author SKM
    @date 14/01/2018
    @date 25/4/2022 refactored
 
    @TODO do we need to remember which side of the interface we are on?
    @todo test method on Chloe's dataset
*/
template<uint32_t dim>
bool findSplitInterfaceElements( const Region<dim>& subdomain,
                                 vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >& interface_elmt_pairs )
 {
    // get coordinate ranges to determine the tolerance of the point matching that will be performed
    /*
       pair<Point<dim>,Point<dim>> box = boundingBox<dim>( subdomain.PerimeterNodesBegin(), subdomain.NodesEnd() );
       const double tolerance = fabs( box.first.DistanceTo( box.second ) ) * numeric_limits<double>::epsilon();
       NOTE: a specific tolerance cannot be used unless a specific comparitor is defined for set below
    */
    
    // multi-element container for all elements that are located on split boundaries
    //       face search key      element      face number
    multimap<set<Point<dim> >,pair<Element<dim>*,uint32_t> > element_face_keys;
   
    // 1. for all elements on the perimeter of the model subdomain,
    //    generate keys from their node coordinates that will later be matched with one-another
    //    in order to connect these elements
    // -------------------------------------
    for ( auto n=subdomain.InteriorCells(); n<subdomain.Cells(); ++n )
        {
           // for those element faces that define the perimeter surface
           for ( uint32_t i{0U}; i<subdomain.PerimeterFaces(n); ++i )
             {
                pair<set<Point<dim> >,uint32_t> face_key;
                const uint32_t face_id{ subdomain.PerimeterFace(n,i) };
                // get the local corner node numbers of the perimeter face
                for ( const auto& nit : subdomain.E(n)->CornerNodesOfFace( face_id ) )
                  // add the corresponding node points to a set that will form the element face key
                  face_key.first.insert( nit->Coordinate() );
                // remembering the face id
                face_key.second = face_id;
                // storing the key in the correspondance search map
                //                              node-coordinate set  element pointer   local face id
                element_face_keys.insert( make_pair( face_key.first, make_pair( subdomain.E(n), face_key.second ) ) );
             }
        }
   
    // 2. searching map for matching interface elements
    // ------------------------------------------------
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
           for ( uint32_t face{0U}; face<e1->Faces(); ++face ) {
               set<Point<dim> >  face_key;
               for ( const auto& nit : e1->CornerNodesOfFace( face ) )
                 face_key.insert( nit->Coordinate() );
               outer_elmt_faces.emplace( make_pair( face_key, make_pair( INSIDE, face ) ) );
             }
           // second element
           Element<dim>* e2 = (*result.first).second.first;
           for ( uint32_t face{0U}; face<e2->Faces(); ++face ) {
               set<Point<dim> >  face_key;
               for ( const auto& nit : e2->CornerNodesOfFace( face ) )
                 face_key.insert( nit->Coordinate() );
               inner_elmt_faces.emplace( make_pair( face_key, make_pair( OUTSIDE, face ) ) );
             }

           // 2. finding the shared faces
           bool found( false );
           int64_t  inner_face_id(-1), outer_face_id(-1);
           for ( const auto& inner_face : inner_elmt_faces ) {
             for ( const auto& outer_face : outer_elmt_faces ) {
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
           interface_elmt_pairs.push_back( make_pair(
                                           make_pair( (*it).second.first, (*it).second.second ),
                                           make_pair( (*result.first).second.first, (*result.first).second.second ) )
                                         );
          }
      }
  
    // 3. checking the results
    if (  interface_elmt_pairs.empty() ) return false;
 
    return true;
    
 } // end findSplitInterfaceElements

template bool findSplitInterfaceElements( const Region<3U>&, vector<pair<pair<Element<3U>*,uint32_t>,pair<Element<3U>*,uint32_t> > >& );
template bool findSplitInterfaceElements( const Region<2U>&, vector<pair<pair<Element<2U>*,uint32_t>,pair<Element<2U>*,uint32_t> > >& );
template bool findSplitInterfaceElements( const Region<1U>&, vector<pair<pair<Element<1U>*,uint32_t>,pair<Element<1U>*,uint32_t> > >& );





/**
        Finds  node-matched Element faces in the supplied model domain, if any.
        The pairs of elements in the model that share these faces are recorded.
        Using property values, the ordering of the elements into inside and outside of the boundary is determined.
        Elements on the inside of this domain are enlisted first in the output Element-face number pairs.
        If neither of the elements is found in the model subdomain, an error is reported.
        
        @param model a  CSMP contiguous or discontiguous model
        @param property_to_distinguish_regions scalar placed on element: the inside of any split boundary will be where the value is lower; where the values are the same no interface will be created
        @param interface_elmt_pairs interface construction data consisting of pairs of elements that share a face across the boundary
        @return number of interfaces created.
*/
template<uint32_t dim>
size_t findSplitInterfaceElements( const Model<dim>& model, const string& property_to_distinguish_regions,
                                   vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >& interface_elmt_pairs )
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    // get coordinate ranges to determine the tolerance of the point matching that will be performed
    /*
       pair<Point<dim>,Point<dim>> box = boundingBox<dim>( subdomain.PerimeterNodesBegin(), subdomain.NodesEnd() );
       const double tolerance = fabs( box.first.DistanceTo( box.second ) ) * numeric_limits<double>::epsilon();
       NOTE: a specific tolerance cannot be used unless a specific comparitor is defined for set below
    */
    const Region<dim>&  model_domain{ model.Region("Model") };
    
    if ( !model.Database().IsDefined( property_to_distinguish_regions.c_str() ) ) {
         csmp_error.Note( ERROR, "findSplitInterfaceElements", property_to_distinguish_regions, "is not defined" );
         return 0U;
      }
    const csmp::Index reg_key = model.Database().StorageKey( property_to_distinguish_regions.c_str() );
    if ( reg_key.place != ELEMENT && reg_key.type != SCALAR )
      csmp_error.Note( ERROR, "findSplitInterfaceElements", property_to_distinguish_regions,
                      "must be scalar placed on the element" );
    
    // multi-element container for all elements that are located on split boundaries
    //   face search key      element      face number
    map<set<Point<dim>>,set<pair<Element<dim>*,uint32_t>>>  element_face_keys;
   
    // 1. for the faces of all elements on the perimeter of the model subdomain,
    //    generate keys from their node coordinates that will later be matched with one-another
    //    in order to connect these elements
    // -------------------------------------
    for ( auto n{model_domain.InteriorCells()}; n<model_domain.Cells(); ++n )
        {
           // making sure that the perimeter elements considered have the same dimension as the model
           // (else, they do share their nodes with higher dimensional elements and are therefore considered as well)
           assert( model_domain.E(n)->IsEquidimensional() );
           // for those element faces that define the perimeter surface
           for ( uint32_t i{0U}; i<model_domain.PerimeterFaces(n); ++i )
             {
                // create a search key for the face from the coordinates of the corner nodes
                pair<set<Point<dim> >,uint32_t> face_key;
                for ( const auto& nit : model_domain.E(n)->CornerNodesOfFace( model_domain.PerimeterFace( n, i ) ) )
                  face_key.first.insert( nit->Coordinate() );
                // remembering the face id
                face_key.second = model_domain.PerimeterFace( n, i );
                // storing the key in the correspondance search map
                set<pair<Element<dim>*,uint32_t>>  elmt_face{ { model_domain.E(n), face_key.second } };
                //                                        node-coordinate set  element pointer   local face id
                auto it = element_face_keys.insert( make_pair( face_key.first, elmt_face ) );
                // if a matching face is detected, its parent element and face id are added to the map
                if ( it.second == false ) {
                     (*it.first).second.insert( make_pair(model_domain.E(n), face_key.second) );
                  }
             }
        }
   
    // 2. creating interface construction data
    // ------------------------------------------------
    // map<set<Point<dim>>,set<pair<Element<dim>*,uint32_t>>>  element_face_keys;
    for ( auto& fit : element_face_keys )
      // if there are matching faces
      if ( fit.second.size() > 1U )
        {
           // interface pairs are generated of the first and the subsequent elmt-face pairs in the map
           auto inside_it = fit.second.begin();
           for ( auto outside_it{next(inside_it,1)}; outside_it!=fit.second.end(); ++outside_it ) {
                // we store the matching element pair
                // vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >  interface_elmt_pairs;
                interface_elmt_pairs.push_back( make_pair( (*inside_it), (*outside_it) ) );
             }
        }
     if (  interface_elmt_pairs.empty() ) {
          ErrorHandler::Instance().Note( ERROR, "findSplitInterfaceElements", "No node-matched interfaces were detected.");
          return 0U;
       }
       
        
    // 3. Re-ordering the element pairs such that elements which belong to inside region are first in the output pair
    // --------------------------------------------------------------------------------------------------------------
    set<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > > pairs_to_eliminate;
    for ( auto& it : interface_elmt_pairs ) {
         const double value1 = it.first.first->Read( reg_key );
         const double value2 = it.second.first->Read( reg_key );
         // value based re-ordering
         if ( value1 > value2 ) swap( it.first, it.second );
         // pairs of elements with the same value will be eliminated
         if ( !(value1 < value2) && !(value1 > value2) )
           pairs_to_eliminate.insert( make_pair( it.first, it.second ) );
      }
  
    // removing the interface without property jump from output vector
    // ---------------------------------------------------------------
    if ( pairs_to_eliminate.size() > 0U ) {
         cerr <<"\n\n"<<"findSplitInterfaceElements: removing "<< pairs_to_eliminate.size() <<" element pairs without ";
         cerr << property_to_distinguish_regions <<" contrast.";
         // removal
         interface_elmt_pairs.erase( remove_if( interface_elmt_pairs.begin(),
                                                interface_elmt_pairs.end(),
                                                [&](auto x) {
                                                     return (pairs_to_eliminate.find(x) != pairs_to_eliminate.end());
                                                  } ) );
      }
 
    return interface_elmt_pairs.size();
    
 } // end findSplitInterfaceElements

template size_t findSplitInterfaceElements( const Model<3U>&, const string&, vector<pair<pair<Element<3U>*,uint32_t>,pair<Element<3U>*,uint32_t> > >& );
template size_t findSplitInterfaceElements( const Model<2U>&, const string&, vector<pair<pair<Element<2U>*,uint32_t>,pair<Element<2U>*,uint32_t> > >& );
template size_t findSplitInterfaceElements( const Model<1U>&, const string&, vector<pair<pair<Element<1U>*,uint32_t>,pair<Element<1U>*,uint32_t> > >& );






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
    Loops over the surface cells of the model subdomain,
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
    size_t non_surface_elements(0U); // NB testing must involve only interior cells
    for ( auto it=subdomain.CellsBegin(); it!=subdomain.PerimeterCellsBegin(); ++it )
      // this method only considers surface elements
      if ( (*it)->IsSurface() ) {
           auto normal = (*it)->UnitNormal();
           const uint32_t neighbors{ (*it)->Neighbors() };
           for ( uint32_t i{0U}; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  auto nbor_normal = (*it)->Neighbor(i)->UnitNormal();
                  // projection
                  double scalar_product = dotProduct( normal, nbor_normal );
                  if ( scalar_product < 0. )
                    return false;
               }
         }
       else non_surface_elements++;
   
    if ( non_surface_elements > 0U )
      ErrorHandler::Instance().Note( ERROR, "checkNeighborNormalsForConsistentOrientation (3D):",
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
      if ( (*it)->IsLine() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( uint32_t i{0U}; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double result(0.);
                  for ( uint32_t j{0U}; j<2U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_line_elements++;
   
    if ( non_line_elements > 0U )
      ErrorHandler::Instance().Note( ERROR, "checkNeighborNormalsForConsistentOrientation (2D):",
                                       subdomain.Name(), "region contained not only line elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation






/**
    Tries to determine whether surface mesh has consistent element normal orientations (none of the element numberings is flipped).
    @pre surface patch must be contiguous.
*/
static bool doNormalsInContiguousSurfacePatchPointToSameSide( const Region<3U>& subdomain )
 {
    subdomain.RenumberNodes();
    
    // performing somewhat expensive normal computation only once on each node
    vector<Point<3U>> vertex_normals( subdomain.Nodes() );
    for ( const auto& nit : subdomain.NodeVector() )
      vertex_normals[ nit->Idx() ] = nit->VertexNormal( subdomain.NodesBegin(), subdomain.NodesEnd() );
      
    // 1. comparing the node normals with one-another does not really answer the question, it only detects discontinuites in surface derivative
    // (getting average normal orientation and then comparing individual normals with it)
    Point<3U>  avg_normal(0);
    for ( const auto& nrml : vertex_normals ) avg_normal += nrml;
    avg_normal /= static_cast<double>( vertex_normals.size() );
    for ( const auto& nrml : vertex_normals )
      if ( dotProduct( avg_normal, nrml ) < 0.5 ) {
           cout <<"\n"<<"doNormalsInContiguousSurfacePatchPointToSameSide: not all vertex normals point in about the same direction ";
           cout <<"(this may also be caused by a discontinuity in the surface)."<< endl;
           break;
        }
      
    // 2. comparing the node normals that are consistent by construction with the normals of the elements
    unordered_set<Element<3U>*> elmts_to_be_flipped;
    for ( const auto& eit : subdomain.CellVector() ) {
         if ( !eit->IsSurface() ) throw csmp::Exception( ERROR, "doNormalsInContiguousSurfacePatchPointToSameSide",
                                                        "detected Element that is not a surface element");
         Point enrml = eit->UnitNormal();
         for ( auto nit = eit->NodesBegin(); nit != eit->NodesEnd(); ++nit )
           if ( dotProduct( enrml, vertex_normals[ (*nit)->Idx() ] ) < 0. ) {
                elmts_to_be_flipped.insert( eit );
                break;
             }
       }
    
    // 3. reporting out
    if ( !elmts_to_be_flipped.empty() ) {
         cout <<"\n"<<"doNormalsInContiguousSurfacePatchPointToSameSide: detected "<< elmts_to_be_flipped.size();
         cout <<" inconsistently numbered elements: "<< endl;
         return false;
      }
      
    return true;
      
 } // end doNormalsInContiguousSurfacePatchPointToSameSide





/**
    Stub: ID there are no boundaries so this should be a compile time assert
    
    @todo TODO: use static_assert<> here on the template argument
*/
template<>
bool checkNeighborNormalsForConsistentOrientation( const Region<1U>& subdomain )
 {
    ErrorHandler::Instance().Note( ERROR, "checkNeighborNormalsForConsistentOrientation (1D):",
                                     subdomain.Name(), "one-dimensional models have no boundaries." );
    return false;
   
 } // end checkNeighborNormalsForConsistentOrientation



template<template<uint32_t> class CELL>
double angleBetweenSurfaceCells( const CELL<3>* const cell1, const CELL<3>* const cell2 )
 {
    assert( cell1 != nullptr );
    assert( cell2 != nullptr );
    assert( cell1->FE()->IsSurface() );
    assert( cell2->FE()->IsSurface() );
    
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
    assert( cell1->FE()->IsLine() );
    assert( cell2->FE()->IsLine() );
    
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
          for ( uint32_t i{0U}; i<n_parents; ++i )
            if ( nit->Parent(i)->IsVolume() )
              {
                 // checking which of the neighbor nodes of the parent element
                 // are also contained in the shared nodes vector
                 const auto n_node_nbors{ nit->Neighbors() };
                 for( uint32_t j{0U}; j<n_node_nbors; ++j ) {
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
              for ( uint32_t segm{0}; segm < n_segments; ++segm ) {
                   vector<uint32_t> snids;
                   it.first->FE()->NodesOfSegment( segm, snids );
                   // making set of segment node pointers to search for
                   set<Node<3>*> segm_nodes;
                   const auto n_segm_nodes{ snids.size() };
                   for ( uint32_t j{0}; j<n_segm_nodes; ++j )
                     segm_nodes.insert( it.first->N( snids[j] ) );
                   // searching segm_nodes for the nodes previously associated with the element
                   bool all_nodes_found{true};
                   const auto n_edge_nodes{ it.second.size() };
                   for ( uint32_t j{0}; j<n_edge_nodes; ++j )
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





/**

Returns the ID (0..n-1) of any node which is located within the tolerance of the
target coordinates. If the node cannot be found it returns max() of index type.

@section arguments Input Arguments

The current model that shall be searched for the node, the node
coordinates, and the tolerance which shall be applied in comparing the
supplied coordinates with those of the actual node points.

@return The Idx (0..n-1) of the node of interest or UINTMAX
(if this node does not exist) will
be returned.

@section application Application

The method is used to retrieve point locations from the mesh in order to
identify points that cannot be grouped into individual families using
the ANSYS mesher.

*/
size_t  findNode( const Model<3U>& sg, double nx, double ny, double nz,
                  double tolerance )
 {
    const Region<3>&  sgroup(sg.Region("Model"));
 
    for ( auto it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ ) {
         if ( fabs(nx-(*it)->x()) <= tolerance and
              fabs(ny-(*it)->y()) <= tolerance and
              fabs(nz-(*it)->z()) <= tolerance )
           return (*it)->Idx();
      }
 
    stringstream  out("The targeted node with the coordinate (x,y,z): ");
    out << nx <<" "<< ny <<" "<< nz <<" could not be found; ";
    out <<" returning node index="<< std::numeric_limits<uint32_t>::max() << endl;
    throw csmp::Exception( WARNING, "findNode", out.str() );
    
    return std::numeric_limits<uint32_t>::max();
     
} // end find_node



size_t  findNode( const Model<2U>& sg, double nx, double ny,
                  double tolerance )
 {
    const Region<2>&  sgroup(sg.Region("Model"));
 
    for ( auto it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ ) {
           if ( fabs(nx-(*it)->x()) <= tolerance and
                fabs(ny-(*it)->y()) <= tolerance )
             return (*it)->Idx();
      }
 
    stringstream  out("The targeted node with the coordinates (x,y): ");
    out << nx <<" "<< ny <<" could not be found; ";
    out <<" returning node index="<< UINT_MAX << endl;
    throw csmp::Exception( WARNING, "findNode", out.str() );
    
    return std::numeric_limits<uint32_t>::max();
     
} // end find_node



size_t  findNode( const Model<1U>& sg, double nx, double tolerance )
 {
    const Region<1>&  sgroup(sg.Region("Model"));
 
    for ( auto it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ )
      if ( fabs(nx-(*it)->x()) <= tolerance ) return (*it)->Idx();
         
    stringstream  out("The targeted node with the coordinate (x): ");
    out << nx <<" could not be found; ";
    out <<" returning node index="<< UINT_MAX << endl;
    throw csmp::Exception( WARNING, "findNode", out.str() );
    
    return std::numeric_limits<uint32_t>::max();
     
} // end find_node


/**
    Generic version for 1-3 dimensions, using Point object to identify the node
    location.
    
    If a node is found its local Idx() number is returned
    (care has to be taken that this index is a valid number.
    
    If the node cannot be found, -1, is returned.
    
    @attention if verbose is on and the point cannot be found this is reported.
    
    @author SKM 22/9/2014.
*/
template<uint32_t dim>
long  findNode( const Model<dim>& sg, const Point<dim>& pxyz, double tolerance, bool verbose )
 {
    const Region<dim>&  sgroup(sg.Region("Model"));
 
    for ( auto it=sgroup.NodesBegin(); it!=sgroup.NodesEnd(); it++ )
      if ( pxyz.CoincidesWithWithinTolerance( (*it)->Coordinate(), tolerance ) )
        return (*it)->Idx();

    ErrorHandler& csmp_error(ErrorHandler::Instance());
   
    if ( verbose ) {
         stringstream  out("The targeted node with the coordinate (x): ");
         out << pxyz <<" could not be found; ";
         out <<" returning node index="<< -1 << endl;
         csmp_error.Note( WARNING, "findNode:", out.str() );
      }
    return -1;
     
} // end find_node

template long findNode( const Model<1U>&, const Point<1U>&, double, bool );
template long findNode( const Model<2U>&, const Point<2U>&, double, bool );
template long findNode( const Model<3U>&, const Point<3U>&, double, bool );



/// prints sorted global element node numbers in a compact way
template<uint32_t dim, template<uint32_t> class CELL>
void printNodes( const CELL<dim>& c )
 {
    set<size_t> nodes;
    for ( uint32_t i{0U}; i<c.Nodes(); i++ ) nodes.insert( c.N(i)->Idx() );
    cout <<" "<< c.Idx() <<": ";
    for ( auto& it : nodes ) cout << it <<",";
    cout <<" ";
 }

template void printNodes( const Element<1U>& );
template void printNodes( const Element<2U>& );
template void printNodes( const Element<3U>& );
template void printNodes( const Face<1U>& );
template void printNodes( const Face<2U>& );
template void printNodes( const Face<3U>& );
template void printNodes( const InterFace<1U>& );
template void printNodes( const InterFace<2U>& );
template void printNodes( const InterFace<3U>& );





/**
    by comparison of barycentre locations, finds overlapping cells and reports them; returns true if collocated cells were found.
    
    @param collocated_cells_to_eliminate pointers to duplicate cells so that these can be eliminated
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool findCollocatedCells( typename std::vector<CELL<dim>*>::const_iterator first,
                          typename std::vector<CELL<dim>*>::const_iterator last,
                          std::vector<CELL<dim>*>& collocated_cells_to_eliminate )
 {
     collocated_cells_to_eliminate.clear();
     
     map<Point<dim>,set<CELL<dim>*>> collocated_cell_map;
     
     while( first != last ) {
          Point<dim> bctr = (*first)->BaryCenter();
          auto insertion = collocated_cell_map.insert( make_pair( bctr, set{ (*first) } ) );
          // if there is already a cell with this barycentre, the pointer to it is added to this key
          if ( insertion.second == false ) (*insertion.first).second.insert( (*first) );
          first++;
       }
     
     // analysing results
     bool found_collocated_cells{ false };
     for ( const auto& it : collocated_cell_map )
       // if there are collacated cells
       if ( it.second.size() > 1 )
         {  // recording these extra cells
            for ( auto cit=next(it.second.begin(),1); cit!=it.second.end(); ++cit )
              collocated_cells_to_eliminate.push_back( (*cit) );
            // reporting them
             if ( !found_collocated_cells ) {
                  cout <<"\n"<<"findCollocatedCells: detected collocated cells:"<< endl;
                  found_collocated_cells= true;
               }
             cout << it.first <<" cell indices: ";
             for ( const auto& c : it.second )
               cout << c->Idx() <<" ";
         }
     
     // the duplicated cells are expected to be unique
     if ( !collocated_cells_to_eliminate.empty() )
       collocated_cells_to_eliminate.shrink_to_fit();
     
     if ( collocated_cells_to_eliminate.empty() ) return false;
     return true;
 }

template
bool findCollocatedCells( typename vector<Element<3U>*>::const_iterator,
                          typename vector<Element<3U>*>::const_iterator,
                          vector<Element<3U>*>& );




template<uint32_t dim>
bool areUnitNormalsToFacesAreOutwardPointing( const Element<dim>* eptr )
 {
    assert( eptr != nullptr );
    Point<dim> ebctr = eptr->BaryCenter();
    bool nrmls_are_outward_pointing{true};
    
    for ( uint32_t face{0}; face<eptr->Faces(); ++face ) {
         // computing unit vector from barycentre to face barycenter
         Point<dim> ctr_face_vec = eptr->FaceBaryCenter(face) - ebctr;
         ctr_face_vec.NormalizeLengthTo(1.);
         // computing unit normal to face
         vector<double> unrml;
         eptr->UnitNormalToFace( face, unrml );
         Point<dim> face_nrml(unrml);
         // if the unit_vec is pointing in opposite direction of normal, normal is inward pointing
         if ( dotProduct(ctr_face_vec,face_nrml) < 0. ) {
              cout <<"\n\n"<<"Element "<< eptr->Idx() <<": face "<< face <<": normal "<< face_nrml <<" is inward pointing."<< endl;
              nrmls_are_outward_pointing = false;
           }
     }
   return nrmls_are_outward_pointing;
    
 } // end

template bool areUnitNormalsToFacesAreOutwardPointing( const Element<3>* );
template bool areUnitNormalsToFacesAreOutwardPointing( const Element<2>* );




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
    size_t                    issues{0};
    string                    celltype("Element");
    if constexpr ( is_same< CELL<dim>,Face<dim> >::value ) celltype = "Face";
    if constexpr ( is_same< CELL<dim>,InterFace<dim> >::value ) celltype = "InterFace";

    auto   copy_of_first{ first };
    size_t max_cell_idx{0};

    // 1. checking that all cells stored in the container are valid
    // ------------------------------------------------------------
    const string check1("\nintegrityCheck: Are all cells stored in the container valid?");
    bool first_call{true};
    
    while ( first != last ) {
         max_cell_idx = max( max_cell_idx, (*first).Idx() );
         // FE policy
         if ( (*first).FE() == nullptr ) {
              if ( first_call ) { cerr << check1; first_call=false; }
              cerr <<"\n\t"<< celltype << (*first).Idx() <<": FE pointer corrupt.";
              issues++;
           }
         first++;
      }
    if ( !first_call ) cerr << endl;
       
    // 2. checking that all connections between nodes and cells are valid
    // ------------------------------------------------------------------
    const string check2("\nintegrityCheck: Are all the nodes connected to the cell valid?");
    // global cell-id,local node-id
    multimap<size_t,uint32_t>  missing_nodes;
    first = copy_of_first;
    first_call = true;

    vector<const Node<dim>*>  shared_nodes;
    shared_nodes.reserve( n_cells * dim );

    while ( first != last ) {
         // connected nodes
         for ( uint32_t i{0U}; i<(*first).FE()->Nodes(); ++i ) {
               if ( (*first).N(i) == nullptr ) {
                    if ( first_call ) { cerr << check2; first_call=false; }
                    cerr <<"\n\t"<< celltype << (*first).Idx() <<": node: "<< i <<": node pointer corrupt.";
                    missing_nodes.insert( make_pair( (*first).Idx(), i ) );
                    issues++;
                 }
               else shared_nodes.push_back( (*first).N(i) );
           }
         first++;
      }
    if ( !first_call ) cerr << endl;
    
    
   // 3. checking node parent connectivity after removing duplicate nodes
   // -------------------------------------------------------------------
   const string check3("\nintegrityCheck: Are all the parent elements of the nodes valid?");
   first_call = true;

   sort( shared_nodes.begin(), shared_nodes.end() );
   shared_nodes.erase( unique(shared_nodes.begin(), shared_nodes.end()), shared_nodes.end() );
   
   for ( const auto& nit : shared_nodes ) {
       if ( nit == nullptr ) cerr <<"\n\t"<<"detected 'nullptr' node.";
       else
         for ( uint32_t i{0U}; i<nit->Parents(); ++i ) {
              if ( nit->Parent(i) == nullptr ) {
                   if ( first_call ) { cerr << check3; first_call=false; }
                   cerr <<"\n\t\t"<<"Node "<< nit->Idx() <<": parent-element vector "<< i <<" contains nullptr.";
                   issues++;
                }
              if ( nit->Parent(i) && nit->Parent(i)->FE() == nullptr ) {
                   if ( first_call ) { cerr << check3; first_call=false; }
                   cerr <<"\n\t\t"<<"Node "<< nit->Idx() <<": parent-element vector "<< i <<" contains Element whose FiniteElement pointer = nullptr.";
                   issues++;
                }
           }
     }
    if ( !first_call ) cerr << endl;


   // 4. node to node connectivity is tested
   // --------------------------------------
   const string check4("\nintegrityCheck: Are all node-to-node connections valid?");
   first_call = true;

   for ( const auto& nit : shared_nodes ) {
       if ( nit == nullptr ) cerr <<"\ndetected 'nullptr' node.";
       else
         for ( uint32_t i{0U}; i<nit->Neighbors(); ++i )
            if ( nit->Neighbor(i) == nullptr ) {
                 if ( first_call ) { cerr << check4; first_call=false; }
                 cerr <<"\n\tneighbor "<< i <<" of Node "<< nit->Idx() <<": is corrupt.";
                 issues++;
              }
     }
    if ( !first_call ) cerr << endl;
   
   
    // 5. checking that all cells have at least one neighbor
    // -----------------------------------------------------
    // if not, this is not necessarily an issue, only if the cells are equidimensional
    const string check5("\nintegrityCheck: Are there cells without any neighbors?");
    multimap<size_t,uint32_t>  missing_nbors;
    first = copy_of_first;
    first_call = true;
   
    while ( first != last ) {
           // there should be at least one cell neighbor
           size_t n_valid_nbors{0};
           for ( uint32_t i{0U}; i<(*first).Neighbors(); ++i ) {
                 if ( (*first).Neighbor(i) != nullptr ) n_valid_nbors++;
                 else missing_nbors.insert( make_pair( (*first).Idx(), i ) );
             }
           if ( n_valid_nbors == 0 ) {
                if ( first_call ) { cerr << check5; first_call=false; }
                cerr <<"\n\t"<< celltype <<" "<< parseFiniteElementType((*first).FE_Type()) <<":"<< (*first).Idx() <<": has no neighbors.";
                if ( (*first).IsEquidimensional() ) issues++;
             }
         first++;
      }
    if ( !first_call ) cerr << endl;


    // 6. checking that all non-null neighbors of the cells are valid
    // --------------------------------------------------------------
    string check6("\nintegrityCheck: Are all non-null nodes and cell neighbors valid?");
    first = copy_of_first;
    first_call = true;
   
    while ( first != last ) {
           // printing message before potentially catastrophic failure occurs
           // valid cell neighbors should not be corrupt
           for ( uint32_t i{0U}; i<(*first).Neighbors(); ++i ) {
               if ( (*first).Neighbor(i) != nullptr ) {
                    // Neighbor element is corrupted so that it no longer has no storage for nodes or elements
                    if ( (*first).Neighbor(i)->Nodes() == 0 || (*first).Neighbor(i)->Neighbors() == 0 ) {
                         if ( first_call ) { cerr << check6; first_call=false; }
                         cerr <<"\n\t"<< celltype <<": "<< (*first).Idx() <<": "<< parseAbbreviated_FE_Type((*first).FE_Type()) <<": neighbor: "<< i <<" has unitialised Node or Neighbor Storage.";
                         issues++;
                      }
                    // FE pointer invalid or FE pointer pointing to FiniteElement base class
                    if ( (*first).Neighbor(i)->FE() == nullptr ) {
                        if ( first_call ) { cerr << check6; first_call=false; }
                        cerr <<"\n\t"<< celltype <<": "<< (*first).Idx() <<": "<< parseAbbreviated_FE_Type((*first).FE_Type()) <<": neighbor: "<< i <<": FE pointer is NULL.";
                        issues++;
                      }
                    // FE pointer invalid or FE pointer pointing to FiniteElement base class
                    if ( (*first).Neighbor(i)->FE() && (*first).Neighbor(i)->FE_Type() == UNKNOWN ) {
                        if ( first_call ) { cerr << check6; first_call=false; }
                        cerr <<"\n\t"<< celltype <<": "<< (*first).Idx() <<": "<< parseAbbreviated_FE_Type((*first).FE_Type()) <<": neighbor: "<< i <<": FE pointer points to FE base class.";
                        issues++;
                      }
                    // Neighbor ID higher than elements in the mesh
                    if ( (*first).Neighbor(i)->Idx() > max_cell_idx ) {
                         if ( first_call ) { cerr << check6; first_call=false; }
                         cerr <<"\n\t"<< celltype <<": "<< (*first).Idx() <<": "<< parseAbbreviated_FE_Type((*first).FE_Type()) <<": neighbor: "<< i <<" has ID out of range.";
                         issues++;
                      }
                 }
              }
           first++;
      }
    if ( !first_call ) cerr << endl;
      
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

 /* TESTING
           cerr <<"\t"<< parseAbbreviated_FE_Type( (*first).FE_Type() ) <<":"<< (*first).Idx() <<" ("<< celltype <<"), barycenter: "<< (*first).BaryCenter() <<", node flags: ";
           for ( uint32_t i{0U}; i<(*first).Nodes(); ++i )
             cerr <<" "<< parseBoundary((*first).N(i)->AtBoundary());
           cerr << endl;
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
           for ( uint32_t i{0U}; i<it.Nodes(); ++i )
             cerr <<" "<< parseBoundary( it.N(i)->AtBoundary() );
           cerr << endl;
           // valid cell neighbors should not be corrupt
           for ( uint32_t i{0U}; i<it.Neighbors(); ++i ) {
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







template<uint32_t dim, template<uint32_t> class CELL>
size_t duplicatesCheck( typename std::vector<CELL<dim>*>::const_iterator first,
                        typename std::vector<CELL<dim>*>::const_iterator last )
 {
    string  celltype("Element");
    if constexpr ( is_same< CELL<dim>,Face<dim> >::value ) celltype = "Face";
    if constexpr ( is_same< CELL<dim>,InterFace<dim> >::value ) celltype = "InterFace";

    long n_cells = distance( first, last );
    
    unordered_set<CELL<dim>*>  unique_cell_set( first, last );
    size_t n_duplicates = n_cells - unique_cell_set.size();
    
    if ( n_duplicates > 0 ) {
         // --------------------------------------------------------------
         cout <<"\n"<<"duplicatesCheck<"<< celltype <<"<"<< dim <<">*> Are there any duplicate pointers in supplied iterator sequence?";
         cout <<"\n\n\t"<<"found "<< n_duplicates <<" "<< celltype <<" pointers."<< endl;
      }
      
    return n_duplicates;
    
 } // end

template size_t duplicatesCheck<3,Element>( vector<Element<3>*>::const_iterator, vector<Element<3>*>::const_iterator );
template size_t duplicatesCheck<2,Element>( vector<Element<2>*>::const_iterator, vector<Element<2>*>::const_iterator );
template size_t duplicatesCheck<1,Element>( vector<Element<1>*>::const_iterator, vector<Element<1>*>::const_iterator );
template size_t duplicatesCheck<3,Face>( vector<Face<3>*>::const_iterator, vector<Face<3>*>::const_iterator );
template size_t duplicatesCheck<2,Face>( vector<Face<2>*>::const_iterator, vector<Face<2>*>::const_iterator );
template size_t duplicatesCheck<1,Face>( vector<Face<1>*>::const_iterator, vector<Face<1>*>::const_iterator );
template size_t duplicatesCheck<3,InterFace>( vector<InterFace<3>*>::const_iterator, vector<InterFace<3>*>::const_iterator );
template size_t duplicatesCheck<2,InterFace>( vector<InterFace<2>*>::const_iterator, vector<InterFace<2>*>::const_iterator );
template size_t duplicatesCheck<1,InterFace>( vector<InterFace<1>*>::const_iterator, vector<InterFace<1>*>::const_iterator );





template<uint32_t dim, template<uint32_t> class CELL>
size_t setNeighborsWithInvalidFE_PointersTo_nullptr( typename plf::colony<CELL<dim>>::iterator first,
                                                     typename plf::colony<CELL<dim>>::iterator last  )
 {
    size_t n_broken_nbors{0ul};
    while( first != last ) {
         for ( uint32_t i{0u}; i<(*first).Neighbors(); ++i )
           if ( (*first).Neighbor(i) && !(*first).FE() ) {
                (*first).UnassignNeighbor(i);
                n_broken_nbors++;
             }
         first++;
      }
    return n_broken_nbors;
 } // end

template size_t setNeighborsWithInvalidFE_PointersTo_nullptr<3,Element>( plf::colony<Element<3>>::iterator, plf::colony<Element<3>>::iterator );
template size_t setNeighborsWithInvalidFE_PointersTo_nullptr<2,Element>( plf::colony<Element<2>>::iterator, plf::colony<Element<2>>::iterator );
template size_t setNeighborsWithInvalidFE_PointersTo_nullptr<1,Element>( plf::colony<Element<1>>::iterator, plf::colony<Element<1>>::iterator );






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
          assert( etype != UNKNOWN );
          int n_connected_neighbors{0};
          
          // 1. checking that elements have equivalent types as neighbors
          for ( uint32_t i{0U}; i<(*first)->Neighbors(); ++i )
            if ( (*first)->Neighbor(i) ) {
                 if constexpr ( dim == 3U ) {
                      if ( (*first)->IsVolume() && !(*first)->Neighbor(i)->IsVolume() ) {
                           cerr <<"\nconnectivityCheck: ERROR: volume Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           issues++;
                        }
                      if ( (*first)->IsSurface() && !(*first)->Neighbor(i)->IsSurface() ) {
                           cerr <<"\nconnectivityCheck: ERROR: surface Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           cerr <<"\n\n"<<"first element:";
                           (*first)->Out();
                           cerr <<"\n\n"<<"wrongly connected neighbor element:";
                           (*first)->Neighbor(i)->Out();
                           issues++;
                        }
                      if ( (*first)->IsLine() && !(*first)->Neighbor(i)->IsLine() ) {
                           cerr <<"\nconnectivityCheck: line Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           issues++;
                        }
                   }
                 if constexpr ( dim == 2U ) {
                      if ( (*first)->IsLine() && !(*first)->Neighbor(i)->IsLine() ) {
                           cerr <<"\nconnectivityCheck: ERROR: line Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a ";
                           cerr << parseAbbreviated_FE_Type( (*first)->Neighbor(i)->FE_Type() );
                           issues++;
                        }
                      if ( (*first)->IsSurface() && !(*first)->Neighbor(i)->IsSurface() ) {
                           cerr <<"\nconnectivityCheck: ERROR: surface Element "<< (*first)->Idx() <<": neighbor("<< i <<") is a line element!";
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
                    cerr <<"\nconnectivityCheck: ERROR: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
                    cerr <<" has only "<< n_connected_neighbors <<" neighbor(s).";
                    issues++;
                  }
            }
          
          if ( isPrism(etype) ) {
                if ( n_connected_neighbors <= 2 ) {
                    cerr <<"\nconnectivityCheck: ERROR: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
                    cerr <<" has only "<< n_connected_neighbors <<" neighbor(s).";
                    issues++;
                  }
            }
          
          if ( isHexahedral(etype) || isPyramid(etype) ) {
                if ( n_connected_neighbors <= 3 ) {
                    cerr <<"\nconnectivityCheck: ERROR: Element "<< (*first)->Idx() <<": "<< parseAbbreviated_FE_Type(etype);
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




// array-based version
template<uint32_t dim>
pair<array<double,dim>,array<double,dim>>  boundingBox1( typename vector<Node<dim>*>::const_iterator first,
                                                         typename vector<Node<dim>*>::const_iterator last )
 {
    if ( first == last ) return make_pair( array<double,dim>{}, array<double,dim>{} );
    
    array<double,dim> pmin, pmax;
    if constexpr ( dim == 3U ) {
         pmin[0]=1.0e30;  pmin[1]=1.0e30;  pmin[2]=1.0e30;
         pmax[0]=-1.0e30; pmax[1]=-1.0e30; pmax[2]=-1.0e30;
      }
    else if constexpr ( dim == 2U ) {
         pmin[0]=1.0e30; pmin[1]=1.0e30;
         pmax[0]=-1.0e30; pmax[1]=-1.0e30;
      }
    else { pmin[0]=1.0e30; pmax[0]=-1.0e30; }
    
    while ( first != last ) {
         pmin = min( pmin, (*first)->Coordinate().CoordinateArray() );
         pmax = max( pmax, (*first)->Coordinate().CoordinateArray() );
         first++;
      }
    
    return make_pair( pmin, pmax );
    
 } // end boundingBox

template pair<array<double,3>,array<double,3>>  boundingBox1<3>( vector<Node<3>*>::const_iterator, vector<Node<3>*>::const_iterator );
template pair<array<double,2>,array<double,2>>  boundingBox1<2>( vector<Node<2>*>::const_iterator, vector<Node<2>*>::const_iterator );
template pair<array<double,1>,array<double,1>>  boundingBox1<1>( vector<Node<1>*>::const_iterator, vector<Node<1>*>::const_iterator );






/**
   returns true if the elements contain each others barycentre,
      only considers equidimensional elements.
*/
template<uint32_t dim>
bool interPenetrating( const Element<dim>* const elmt1, const Element<dim>* const elmt2 )
 {
    if constexpr ( dim == 3 )
      if ( !elmt1->IsVolume() || !elmt2->IsVolume() ) {
           cerr <<"\ninterPenetrating<3>: only works for equidimensional (volumetric) elements.\n";
           return false;
        }

    if constexpr ( dim == 2 )
      if ( !elmt1->IsSurface() || !elmt2->IsSurface() ) {
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
        
        if (!(*eit)->IsVolume()) continue;

        // 2. Test against axis-aligned bounding box

        double minx = +std::numeric_limits<double>::max();
        double miny = +std::numeric_limits<double>::max();
        double minz = +std::numeric_limits<double>::max();
        double maxx = -std::numeric_limits<double>::max();
        double maxy = -std::numeric_limits<double>::max();
        double maxz = -std::numeric_limits<double>::max();
        const auto iNrNodes = (*eit)->Nodes();
        for ( uint32_t iNode = 0; iNode < iNrNodes; ++iNode ) {
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
        for ( uint32_t iFace = 0u; iFace < iNrFaces && !reject; ++iFace ) {
          fnids = fe->NodesOfFace(iFace);
          const size_t iNrFacePts = fnids.size();
          for ( uint32_t iFacePt = 0u; iFacePt < iNrFacePts; iFacePt += 2) {
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
    for ( uint32_t i{0U}; i<n_nodes; ++i )
      node_points.insert( eptr->N(i)->Coordinate() );
 
    return n_nodes - node_points.size();
    
 } // end collocatedNodes

template size_t collocatedNodes( const Element<1>* const );
template size_t collocatedNodes( const Element<2>* const );
template size_t collocatedNodes( const Element<3>* const );





/**
    Pretty prints line elements as a chain from beginning to end.
    
     Ideally (where the numbers are nodes and the labels are BOX_BOUNDARY flags)
     you should get something like: TOP 1->-0 0--11 11->-12...56->-56 BOTTOM.
*/
template<uint32_t dim>
size_t printLineElementRegion( const Model<dim>& model, const char* region_name, bool renumber_nodes )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const Region<dim>& line_domain( model.Region(region_name) );
    if ( renumber_nodes ) line_domain.RenumberNodes();
    
    const Element<dim>* eptr1 = (*line_domain.PerimeterCellsBegin());
    const Element<dim>* eptr2 = (*prev(line_domain.CellsEnd(),1));
    
    const Element<dim>* previous_ptr{nullptr};
    size_t              traversed_elmts{0U}, total_line_elmts{ line_domain.Cells() };
    
    bool forward = ( eptr1->Neighbor(0) != nullptr ) ? true : false;
    
    cout <<"\n\nprintLineElementRegion: '"<< region_name <<"': printing chain of elements ";
    if ( forward ) {
         cout <<"in forward direction:"<< endl;
         while( eptr1 != nullptr ) {
              if ( !eptr1->IsLine() )
                throw csmp::Exception( ERROR, "printLineElementRegion", "current element is not a line element; aborting printing" );
              // Ideally (where the numbers are nodes and the labels are BOX_BOUNDARY flags)
              // we should get something like: TOP 1->-0 0--11 11->-12...56->-56 BOTTOM
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
              if ( traversed_elmts > total_line_elmts ) {
                   csmp_error.Note( ERROR, "printLineElementRegion", "traversed more elements than in regions; check output for issues with line element connectivity" );
                   break;
                }
           }
      }
    // backward
    else {
         cout <<"from back to front:"<< endl;
         while( eptr2 != nullptr ) {
              if ( !eptr2->IsLine() )
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
              if ( traversed_elmts > total_line_elmts ) {
                   csmp_error.Note( ERROR, "printLineElementRegion", "traversed more elements than in regions; check output for issues with line element connectivity" );
                   break;
                }
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
     long n_nodes = distance(first,last);
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
          for ( uint32_t i{1U}; i<dim; i++ ) cout <<", "<< (*(*first))[i];
          cout << endl;
          first++;
       }
     
 } // end printNodeCoordinates

template void printNodeCoordinates<3U>( typename vector<Node<3U>*>::const_iterator, typename vector<Node<3U>*>::const_iterator );
template void printNodeCoordinates<2U>( typename vector<Node<2U>*>::const_iterator, typename vector<Node<2U>*>::const_iterator );
template void printNodeCoordinates<1U>( typename vector<Node<1U>*>::const_iterator, typename vector<Node<1U>*>::const_iterator );




/**
      returns the maximum distance between the nodes in the provided range
*/
template<uint32_t dim>
double maximumNodeSpacing( typename vector<Node<dim>*>::const_iterator first,
                           typename vector<Node<dim>*>::const_iterator last )
 {
    double node_spacing{0.};
    while( first != last ) {
          for ( uint32_t i{0u}; i<(*first)->Neighbors(); ++i  )
            node_spacing = max( node_spacing, distance( (*first)->Coordinate(), (*first)->Neighbor(i)->Coordinate() ) );
          first++;
       }
    
    return node_spacing;
       
 } // end maximumNodeSpacing

template double maximumNodeSpacing<1U>( typename vector<Node<1U>*>::const_iterator,
                                        typename vector<Node<1U>*>::const_iterator );
template double maximumNodeSpacing<2U>( typename vector<Node<2U>*>::const_iterator,
                                        typename vector<Node<2U>*>::const_iterator );
template double maximumNodeSpacing<3U>( typename vector<Node<3U>*>::const_iterator,
                                        typename vector<Node<3U>*>::const_iterator );



template<uint32_t dim>
double minimumNodeSpacing( typename vector<Node<dim>*>::const_iterator first,
                           typename vector<Node<dim>*>::const_iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    assert( first != last );
    
    double node_spacing{ 1.0e30 };
    while( first != last ) {
          for ( uint32_t i{0u}; i<(*first)->Neighbors(); ++i  ) {
               node_spacing = min( node_spacing, distance( (*first)->Coordinate(), (*first)->Neighbor(i)->Coordinate() ) );
               if ( node_spacing <= numeric_limits<double>::epsilon() ) {
                    cerr <<"\n\t"<<"minimumNodeSpacing: found collocated nodes:\n";
                    cerr <<"\n\t\t"<< (*first)->Idx() <<": "<< (*first)->Coordinate() << endl;
                    cerr <<"\n\t\t"<< (*first)->Neighbor(i)->Idx() <<": "<< (*first)->Neighbor(i)->Coordinate() << endl;
                    csmp_error.Note( ERROR, "minimumNodeSpacing", "found collocated nodes" );
                    node_spacing = 1e30;
                 }
            }
          first++;
       }
    
    return node_spacing;
       
 } // end minimumNodeSpacing

template double minimumNodeSpacing<1U>( typename vector<Node<1U>*>::const_iterator,
                                        typename vector<Node<1U>*>::const_iterator );
template double minimumNodeSpacing<2U>( typename vector<Node<2U>*>::const_iterator,
                                        typename vector<Node<2U>*>::const_iterator );
template double minimumNodeSpacing<3U>( typename vector<Node<3U>*>::const_iterator,
                                        typename vector<Node<3U>*>::const_iterator );


// Arithmetic mean: node that the obtained value will be affected by how many times a particular spacing gets evaluated
template<uint32_t dim>
double averageNodeSpacing( typename vector<Node<dim>*>::const_iterator first,
                           typename vector<Node<dim>*>::const_iterator last )
 {
    //ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    assert( first != last );
    
    const auto n_nodes{ distance(first,last) };
    double     sum_of_node_spacings{ 0. };
    
    while( first != last ) {
          double     nbor_spacing_average{0.};
          const auto n_nbor_nodes{ (*first)->Neighbors() };
          for ( uint32_t i{0u}; i<n_nbor_nodes; ++i  )
               nbor_spacing_average += distance( (*first)->Coordinate(), (*first)->Neighbor(i)->Coordinate() );
          nbor_spacing_average /= static_cast<double>(n_nbor_nodes);
          sum_of_node_spacings += nbor_spacing_average;
          first++;
       }
    
    // finding the arithmetic mean
    return sum_of_node_spacings / static_cast<double>(n_nodes);
       
 } // end averageNodeSpacing

template double averageNodeSpacing<1U>( typename vector<Node<1U>*>::const_iterator,
                                        typename vector<Node<1U>*>::const_iterator );
template double averageNodeSpacing<2U>( typename vector<Node<2U>*>::const_iterator,
                                        typename vector<Node<2U>*>::const_iterator );
template double averageNodeSpacing<3U>( typename vector<Node<3U>*>::const_iterator,
                                        typename vector<Node<3U>*>::const_iterator );





/// checks whether point is contained in any of the elements in supplied region returning 'nullptr' or the element in which it is contained
template<uint32_t dim>
const Element<dim>* isContainedIn( const Region<dim>& subdomain, const array<double,dim>& search_point ) {
     vector<double> xyz( search_point.begin(), search_point.end() );
     vector<double> N;
     // testing containment by a linear search that diagnoses whether point is contained
     // by determining whether all element interpolation function values are between zero and one.
     for ( const auto& it : subdomain.CellVector() )
      if ( it->IsEquidimensional() )
       {
           it->N_AtGlobalPoint( N, xyz );
           size_t counter{0ul};
           for( const auto& nval : N ) {
                if ( nval < 0. || nval > 1. ) break;
                counter++;
             }
           // if all N values are within [0..1] the point has been found
           if ( counter == N.size() )
             return it;
       }
     return nullptr;
  }

template const Element<2U>* isContainedIn<2U>( const Region<2U>&, const array<double,2U>& );
template const Element<3U>* isContainedIn<3U>( const Region<3U>&, const array<double,3U>& );



/**
      Finding the neighbors nodes of each node.
      
      This method is equivalent to creating a sparsity pattern for an accumulation.
      However, only up to a single mid-side node per segment is handled.
      
      @test OK SKM 8/12/21
*/
template<uint32_t dim>
void nodeNeighbors( const Region<dim>& subdomain, vector<set<size_t>>& node_neighbors, bool verbose )
 {
    if ( !node_neighbors.empty() ) node_neighbors.clear();
    node_neighbors.resize( subdomain.Nodes() );
    
    // 1. get continuous indices to access the sets contained in the node_neighbor vectors
    subdomain.RenumberNodes();
 
    vector<uint32_t> segm_nodes;

    const auto elmtsEnd{ subdomain.CellsEnd() };
    for ( auto it=subdomain.CellsBegin(); it!=elmtsEnd; ++it ) {
         const uint32_t n_segments{ (*it)->Segments() };
         for ( uint32_t segm_id{0u}; segm_id < n_segments; ++segm_id ) {
              (*it)->FE()->NodesOfSegment( segm_id, segm_nodes );
              // replacing local with global node ids
              for ( auto& sit : segm_nodes ) sit = static_cast<uint32_t>((*it)->N(sit)->Idx());
              // corner nodes
              size_t segm_node1{ subdomain.N(*segm_nodes.begin())->Idx() };
              size_t segm_node2{ subdomain.N(*next(segm_nodes.begin(),1))->Idx() };
              node_neighbors[ segm_node1 ].insert( segm_node2 );
              node_neighbors[ segm_node2 ].insert( segm_node1 );
              // if there is a mid-side node
              if ( segm_nodes.size() == 3U ) {
                   size_t segm_node3{ subdomain.N(*next(segm_nodes.begin(),2))->Idx() };
                   node_neighbors[ segm_node1 ].insert( segm_node3 );
                   node_neighbors[ segm_node3 ].insert( segm_node1 );
                }
              // if there are two midside nodes
              if ( segm_nodes.size() > 3U )
                throw csmp::Exception( ERROR, "nodeNeighbors", "method only handles a single segment midside node" );
           }
      }
      
    // printing the node-neighbor vector for testing
    if ( verbose ) {
      cout <<"\n\nnodeNeighbors: connectivity created for "<< node_neighbors.size() <<" nodes:";
      size_t node{0};
      for ( auto nit : node_neighbors ) {
           cout <<"\n\t" << node <<": ";
           for (  auto i : nit )
             cout << subdomain.N(i)->Idx() <<" ";
           cout <<" ("<< parseBoundary( subdomain.N(node)->AtBoundary() ) <<")";
           node++;
        }
      }
      
 } // end nodeNeighbors

template void nodeNeighbors( const Region<3>&, vector<set<size_t>>&, bool );
template void nodeNeighbors( const Region<2>&, vector<set<size_t>>&, bool );




} // end csmp
