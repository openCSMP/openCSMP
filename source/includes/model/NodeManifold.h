#ifndef CSMP_NODE_MANIFOLD_H
#define CSMP_NODE_MANIFOLD_H

#include "CSMP_definitions.h"
#include "CSMP_global_enumerations.h"
#include "plf_colony.h"

namespace csmp {

class Index;

template<uint32_t> class NodeManifold;

template<uint32_t> class Node;

template<uint32_t> class InterFace;

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
      NodeManifold( plf::colony<Node<dim> >& nodes, const std::vector<size_t>& manifold_nodes );
      
      /// constructs new Manifold from two nodes on the inside and outside of it; @attention nodes must be assigned to this manifold once it has been constructed
      NodeManifold( Node<dim>& node1, Node<dim>& node2 );

      /// constructs manifold from vector pointer and qualifier pairs
      NodeManifold( const manifold& );

      /// adds a node to the manifold storing the interface side, it is on; @note  this might also have implications for Manifold geometry to be addressed later
      bool Add( Node<dim>* );

      /// removes node from the current manifold and sets its manifold pointer to zero because a Node can only belong to a single manifold
      bool Remove( const Node<dim>* const );

      /// asscending sort (scalar on Node or Element only) - default is sorted by pointer in sequence entered
      void SortByVariableValue( const Index& scalar_node_variable );

      ///Assign ------------------------------------------------------
      void Assign( Node<dim>*,
                   std::set<std::pair<InterFace<dim>*,std::pair<uint32_t,INTERFACE_SIDE>>> interface_indexes );


      ///Access -------------------------------------------------------

      /// number of entries
      uint32_t Branches() const noexcept;

      /// access to node
      Node<dim>* const N( size_t branch ) const  noexcept;

      /// Gives size of node interface parent map - should correspond with number of nodes(branches) when configured correctly
      //   MeshManager::ReplaceElementsByInterFaces
      uint32_t NodeMapSize() const;
      
      ///number of interfaces connected to a node
      //   MeshManager::ReplaceElementsByInterFaces
      uint32_t InterFaces( Node<dim>* const ) const noexcept;

      ///Access of single interface of node
      //   MeshManager::ReplaceElementsByInterFaces
      InterFace<dim>* I( Node<dim>* const, uint32_t i );

      ///Pair with interface, and corresponding position of node in nodeconnector of interface
      // MeshManager::ReplaceElementsByInterFaces
      std::pair<InterFace<dim>*, std::pair<uint32_t,INTERFACE_SIDE> >  InterFaceIndex( Node<dim>* const n, uint32_t i );

      /// InterFace vector of node
      // REMOVE-not used
      std::vector< std::pair<InterFace<dim>*, std::pair<uint32_t,INTERFACE_SIDE>>> InterFaceIndexVector( Node<dim>* const n );

     ///Query -----------------------------------------------------------

      /// reports geometric role of the manifold inferred from the TOPOTYPEs of the nodes that care contained in the
      ManifoldType Classify() const noexcept;
      
      /// checks whether all nodes in the  manifold have the same location using operator< of point
      bool AreNodesCollocated() const noexcept;
      
      /// outputs manifold state to data structure used to initialise VData
      std::pair<std::vector<size_t>,ManifoldType> Data() const noexcept;

      /// prints out state of the manifold
      void Out() const;

    private:
      manifold branches_;  ///< vector of Node pointers

      // IDEA: Nodes know their parent elements; NodeManifolds know their parent InterFaces
      // TODO: remove this huge-overhead structure (if needs be, put it as an auxiliary vector inside of MeshManager::ReplaceElementsByInterFaces()
      std::map< Node<dim>*,std::vector<std::pair<InterFace<dim>*,std::pair<uint32_t,INTERFACE_SIDE> >> >   node_parent_interface_map_;     ///interfaces and index for each node on manifold
};


/// check whether through modification of the manifold, the original topology identifier is no longer valid
template<uint32_t dim>
ManifoldType consistencyCheck( const NodeManifold<dim>&, bool verbose );
  
} // end csmp

#endif

