//
//  MeshIterator.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 13/8/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "MeshIterator.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
MeshIterator<dim,CELL>::MeshIterator( Node<dim>* root_node )
 : discovered_nodes_({root_node}),
   current_nodes_({root_node}) 
 {
    // 1. traversal of existing contiguous mesh patch to find all its elements
    while ( !current_nodes_.empty() ) {
        const csmp::Node<dim>*  n_ptr( *current_nodes_.begin() );
        // for all parent elements of the current node
        const uint32_t parent_cells(n_ptr->Parents());
        for ( auto i{0U}; i<parent_cells; ++i ) {
            // for all the nodes of each parent element
            const uint32_t parent_nodes(n_ptr->Parent(i)->Nodes());
            for ( uint32_t j{0U}; j<parent_nodes; ++j )
              // if this node is not the one from which we started
              if ( j != n_ptr->ParentNodeNumber(i) ) {
                  pair<typename set<csmp::Node<dim>*>::iterator, bool>
                    new_node = discovered_nodes_.insert( n_ptr->Parent(i)->N( j ) );
                  if ( new_node.second ) current_nodes_.push_back( n_ptr->Parent( i )->N( j ) );
                }
            // storing the explored element
            discovered_cells_.insert( n_ptr->Parent(i) );
          }
      // removing the node from the discovered (but not yet explored) deque
      current_nodes_.pop_front();
    }

} // end constructor( root_node )



template class MeshIterator<1U,Element>;
template class MeshIterator<2U,Element>;
template class MeshIterator<3U,Element>;

/*
template class MeshIterator<1U,Face>;
template class MeshIterator<2U,Face>;
template class MeshIterator<3U,Face>;

template class MeshIterator<1U,InterFace>;
template class MeshIterator<2U,InterFace>;
template class MeshIterator<3U,InterFace>;
*/

} // end csmp
