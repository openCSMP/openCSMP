//
//  MeshIterator.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 13/8/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_MESH_ITERATOR_H
#define CSMP_MESH_ITERATOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Node;

/**
    Interator to perform a breadth-first traversal of contiguous mesh patch during its construction
    and then contain cell and node sets that can be used to iterate over the mesh tree.
    
    @author SKM
    @date 15/3/2017
    
    @todo needs to be completed to become useful
    @todo use this breadth-first graph traversal rather than any other method to discover patch of CELLs.
*/
template<uint32_t dim, template<uint32_t> class CELL> 
class MeshIterator {
  public:
    /// discovers interconnected mesh patch initiaising the member cells
    explicit MeshIterator( Node<dim>* mesh_node );
    
    typename std::set<CELL<dim>*>::const_iterator CellsBegin() const { return discovered_cells_.begin(); }
    typename std::set<CELL<dim>*>::const_iterator CellsEnd() const { return discovered_cells_.end(); }

    typename std::set<CELL<dim>*>::iterator CellsBegin() { return discovered_cells_.begin(); }
    typename std::set<CELL<dim>*>::iterator CellsEnd() { return discovered_cells_.end(); }   
     
    typename std::set<csmp::Node<dim>*>::const_iterator NodesBegin() const { return discovered_nodes_.begin(); }
    typename std::set<csmp::Node<dim>*>::const_iterator NodesEnd() const { return discovered_nodes_.end(); }

    typename std::set<csmp::Node<dim>*>::iterator NodesBegin() { return discovered_nodes_.begin(); }
    typename std::set<csmp::Node<dim>*>::iterator NodesEnd() { return discovered_nodes_.end(); }   

  private:
    // use deque contain made unique to store the cell references
    std::set<CELL<dim>*>          discovered_cells_;
    std::set<csmp::Node<dim>*>    discovered_nodes_;
    std::deque<csmp::Node<dim>*>  current_nodes_;
};

} // end csmp

#endif /* CSMP_MESH_ITERATOR_H */
