#ifndef CSMP_NODE_MANIFOLD_MANAGER_H
#define CSMP_NODE_MANIFOLD_MANAGER_H

#include "NodeManifold.h"

namespace csmp {

template<size_t> class NodeManifold;
template<size_t> class Node;

/**
@brief Policy of Model class for the management of manifolds that consist of multiplicated nodes

@note this class is not copy constructable because only a single instance is needed.
@note just storage management for manifolds; no targeted access.
*/

template<size_t dim>
class NodeManifoldManager {
  public:
    NodeManifoldManager()                                        = default;
    NodeManifoldManager( const NodeManifoldManager& )            = delete;
    NodeManifoldManager& operator=( const NodeManifoldManager& ) = delete;
    
    /// deleting the dynamically allocated manifolds
    ~NodeManifoldManager();
    
    /// for sorting the manifolds by their pointers
    typedef typename std::deque< NodeManifold<dim>* >::iterator manifoldIterator;
    typedef typename std::deque< NodeManifold<dim>* >::const_iterator manifoldConstIterator;

    manifoldIterator       ManifoldsBegin();
    manifoldIterator       ManifoldsEnd();
    manifoldConstIterator  ManifoldsBegin() const;
    manifoldConstIterator  ManifoldsEnd() const;

    size_t                 Manifolds() const;
    
     // find - no need for this because each node has direct access to connected manifolds

    /// creates a node manifold to which nodes can be added
    NodeManifold<dim>* const NewManifold( Node<dim>* const inside, Node<dim>* const outside, ManifoldType );
    
    /// puts the nodes inside of the manifolds into the ascending order of values of the user specified  variable
    void SortManifoldsByVariableValue( std::string var_name, const csmp::Index& var_index );

    /// deletes = erases manifold from storage container, reordering / compacting as necessary
    bool   Delete( NodeManifold<dim>* );
    
    /// scans for manifolds with a single Node only and deletes them
    size_t DeleteSingleNodeManifolds();

    /// prints manifolds and other stats to screen
    void Out() const;
    
    /// writes manifolds to binary file; relies on unique node indexes
    void OutputNodeManifoldsToBinary( const char* file_name,
                                      const char* node_sorting_variable ) const;

    /// Reads manifolds from binary file, using indices to create pointer connections; the name of the sorting variable for the manifolds is returned
    std::string InputNodeManifoldsFromBinary( std::deque<Node<dim>*>& nodes, const char* file_name );
    
  protected:
     std::deque< NodeManifold<dim>* >  node_manifolds_;                 ///<  sorted deque of pointers to manifolds created with new and delete
     std::string                       current_sort_variable_ = "none"; ///< TODO: update the sort order if some modifications are made
}; 

} // csmp


#endif
