#ifndef CSMP_NODE_MANIFOLD_H
#define CSMP_NODE_MANIFOLD_H

#include "CSMP_definitions.h"
#include "CSMP_global_enumerations.h"
#include "plf_colony.h"

namespace csmp {

template<uint32_t dim> class Node;
class Index;

// TODO: shall we include specific InterFace parent connectivity to NodeManifolds?
// TODO: do we need to track intersections between SB and Boundaries or lower-dim Regions?
  
/**
    @brief: a class that contains pointers to topologically (and most often also physically) collocated nodes.
    SplitBoundaries (node-matching standalone mesh patches connected by interface elements typically consist of NodeManifold objects.
    
    @attention the nodes on the perimeter of a SplitBoundary that terminates in a volume are not duplicated! - therefore they are not manifolds.
*/
template<uint32_t dim>
class NodeManifold {
    public:
      NodeManifold() = default;

      using manifold              = std::vector<Node<dim>*>;
      using manifoldIterator      = typename std::vector<Node<dim>*>::iterator;
      using manifoldConstIterator = typename std::vector<Node<dim>*>::const_iterator;

      /// main constructor in NodeManifoldManager: sorts created vector by Node pointers in ascending order, so that it can be searched for nodes using std::binary_search
      NodeManifold( plf::colony<Node<dim> >& nodes,
                    const std::vector<size_t>& manifold_nodes,
                    ManifoldType );
      
      /// constructs new Manifold from two nodes on the inside and outside of it; @attention nodes must be assigned to this manifold once it has been constructed
      NodeManifold( Node<dim>& node1, Node<dim>& node2, ManifoldType );

      /// constructs manifold from vector pointer and qualifier pairs
      NodeManifold( const manifold&, ManifoldType );

      NodeManifold( const NodeManifold& );
      NodeManifold( NodeManifold&& );
      NodeManifold& operator=( const NodeManifold& );
      NodeManifold& operator=( NodeManifold&& );

      ~NodeManifold();
      
      /// adds a node to the manifold storing the interface side, it is on; @note  this might also have implications for Manifold geometry to be addressed later
      bool Add( Node<dim>* );

      /// removes node from the current manifold and sets its manifold pointer to zero because a Node can only belong to a single manifold
      bool Remove( const Node<dim>* const );

      /// asscending sort (scalar on Node or Element only) - default is sorted by pointer in sequence entered
      void SortByVariableValue( const Index& scalar_node_variable );

      /// number of entries
      uint32_t Branches() const;

      /// access to node
      Node<dim>* const N( size_t branch ) const;
      
      /// reports manifold classifier that indicates the topologic position of the manifold
      ManifoldType GeometricClassifier() const;
      void GeometricClassifier( ManifoldType );
      
      /// checks whether all nodes in the  manifold have the same location using operator< of point
      bool AreNodesCollocated() const;
      
      /// outputs manifold state to data structure used to initialise VData
      std::pair<std::vector<size_t>,ManifoldType> Data() const;

      /// prints out state of the manifold
      void Out() const;

    private:
      manifold      branches_;         ///< vector of Node - classifier pairs
      ManifoldType  parent_geometry_;  ///< classifier for the manifold as a whole
};


/// check whether through modification of the manifold, the original topology identifier is no longer valid
template<uint32_t dim>
ManifoldType consistencyCheck( const NodeManifold<dim>& );
  
} // end csmp

#endif

