#ifndef CSMP_NODE_MANIFOLD_MANAGER_H
#define CSMP_NODE_MANIFOLD_MANAGER_H

#include "CSMP_definitions.h"
#include "NodeManifold.h"

namespace csmp {

template<size_t> class NodeManifold;
template<size_t> class Node;

/**
@brief Policy of Model class for the management of manifolds that consist of multiplicated nodes


@note this class is not copy constructable because only a single instance is needed.
*/

template<size_t dim>
class NodeManifoldManager {
  public:
    NodeManifoldManager();
    NodeManifoldManager( const NodeManifoldManager& ) = delete;
    NodeManifoldManager& operator=( const NodeManifoldManager& ) = delete;
    virtual ~NodeManifoldManager();

    typedef typename std::deque< NodeManifold<dim>* >::iterator manifoldIterator;
    typedef typename std::deque< NodeManifold<dim>* >::const_iterator manifoldConstIterator;

    manifoldIterator       NodeManifoldsBegin();
    manifoldIterator       NodeManifoldsEnd();
    manifoldConstIterator  NodeManifoldsBegin() const;
    manifoldConstIterator  NodeManifoldsEnd() const;
    //manifoldIterator       Manifold( const Manifold<dim>* );
    size_t                 NodeManifolds() const;    

    // Manifold creation, modification & removal
    bool InsertNodeManifold( NodeManifold<dim>* );
    bool RemoveNodeManifold( NodeManifold<dim>* );

    // TODO: is this really meant to be call by value? - (causing 2 copy constructions)
    bool AddNodeManifold( std::vector<Node<dim>*> );

    // TODO: is this really meant to be call by value? - (causing 2 copy constructions)
    bool AddNodeManifold( std::vector<Node<dim>*>, PARENT_GEOMETRIC_ENTITY );

    bool AddNodeManifold( typename std::vector<Node<dim>*>::iterator nodesBegin,
                          typename std::vector<Node<dim>*>::iterator nodesEnd);

    bool AddNodeManifold( typename std::vector<Node<dim>*>::iterator nodesBegin,
                          typename std::vector<Node<dim>*>::iterator nodesEnd,
                          PARENT_GEOMETRIC_ENTITY);

    bool InsertNodeToManifold( NodeManifold<dim>*, Node<dim>* );

    bool RemoveNodeFromManifold( NodeManifold<dim>*, Node<dim>* );
    
    // prints manifolds and other stats to screen 
    void NodeanifoldsOut() const;                                     


  protected:
     std::deque< NodeManifold<dim>* >  node_manifolds_; ///< storage of the manifolds  

}; 

} // csmp


#endif
