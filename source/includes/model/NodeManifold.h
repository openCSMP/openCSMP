#ifndef CSMP_NODE_MANIFOLD_H
#define CSMP_NODE_MANIFOLD_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t dim> class Node;
class Index;

// TODO: distinguish perimeter nodes in SplitBoundaries from interior nodes, using criteria below
// TODO: do we need to track intersections between SB and Boundaries or lower-dim Regions?
/**
  @note This classification is for NodeManifolds (topologically collocated Nodes) only)!
  @note Not all of these classifiers apply to manifolds, like endpoint
  @note NodeManifolds exist only at SplitBoundary objects
  @note SplitBoundary objects exist only inside of models
  @note since SplitBoundaries are surfaces, this is the highest dimension
  @note Classification applies to all nodes within Manifold simultaneously (including intervening ones)
*/
enum class ManifoldType : int8_t {
                                    POINT2,       ///<  intersection or contact point between 2 line SplitBoundaries
                                    POINT3,       ///<  3 line SplitBoundary objects in contact
                                    POINTX,       ///<  multiple SB in contact
                                    END_POINT,    ///< where a line SplitBoundary terminates in a volumetric region
                                    EDGE,         ///<  outer edge of a SplitBoundary terminating within volume
                                    BPOINT1,      ///< intersection between line SB and model boundary
                                    BEDGE,        ///< edge of SB located on model boundary
                                    BPOINT2,             ///< intersection of boundary edge with model boundary
                                    BPOINTX,
                                    INTERSECTION,        ///<  when 2 SplitBoundary objects meet
                                    MULTI_INTERSECTION,  ///<  more than 2 SplitBoundary objects cross
                                    INTERSECTION_POINT,  ///< edge point due the contact between multiple SB objects
                                    BINTERSECTION_POINT, ///< SB's touch each other on model boundary
                                    INTERFACE,           ///<  on SplitBoundary
                                    NOT_CLASSIFIED       ///<  an isolated Point that has therefore no classification
                                };

/// converts classifiers to strings so that they can be printed
std::string parse( ManifoldType );
  
/**
* @brief: a class contains pointers to coincident nodes at split-boundary
*/
template<size_t dim>
class NodeManifold {
    public:
      NodeManifold() = default;

      typedef typename std::vector<std::pair<Node<dim>*,INTERFACE_SIDE> >                   manifold;
      typedef typename std::vector<std::pair<Node<dim>*,INTERFACE_SIDE> >::iterator         manifoldIterator;
      typedef typename std::vector<std::pair<Node<dim>*,INTERFACE_SIDE> >::const_iterator   manifoldConstIterator;

      /// sorts created vector by Node pointers in ascending order, so that it can be searched for nodes using std::binary_search
      NodeManifold( const std::vector<Node<dim>*>&, const std::vector<INTERFACE_SIDE>&, ManifoldType );

      NodeManifold( const manifold&, ManifoldType );

      NodeManifold( const NodeManifold& )            = default;
      NodeManifold( NodeManifold&& )                 = default;
      NodeManifold& operator=( const NodeManifold& ) = default;
      NodeManifold& operator=( NodeManifold&& )      = default;

      ~NodeManifold();
      
      /// asscending sort (scalar on Node or Element only) - default is sorted by pointer in sequence entered
      void SortByVariableValue( const Index& scalar_node_variable );

      /// number of entries
      size_t Branches() const;

      /// access to node
      Node<dim>* const N( size_t branch );
      
      /// where the node resides
      INTERFACE_SIDE InterFaceSide( size_t branch ) const;
      
      /// returns the  Node of an intervening element; else returns nullptr
      Node<dim>* const MIDDLE_Node() const;
      
      /// retrieve node(s) with specific flags (manifold should only contain one intervening Node (MIDDLE)
      std::vector<Node<dim>*>  NodesLocatedAt( INTERFACE_SIDE ) const;

      /// reports manifold classifier that indicates the topologic position of the manifold
      ManifoldType GeometricClassifier() const;
      void GeometricClassifier( ManifoldType );

      /// adds a node to the manifold storing the interface side, it is on; note that this might also have implications for the Manifold geometry which is therefore required
      bool Add( Node<dim>*, INTERFACE_SIDE, ManifoldType );

      /// removes node from the current manifold and sets its manifold pointer to zero because a Node can only belong to a single manifold
      bool Remove( Node<dim>* );

      /// print out information
      void Out() const;

    private:
      manifold      branches_;         ///< vector of Node - classifier pairs
      ManifoldType  parent_geometry_;  ///< classifier for the manifold as a whole
};


// TODO: needs rigorous testing; most certainly incomplete
/// check whether through modification of the manifold, the original topology identifier is no longer valid
template<size_t dim>
ManifoldType consistencyCheck( const NodeManifold<dim>* const );
  
/// determines the type of manifold on the basis of its node members and their BOX_BOUNDARY flags
//template<size_t dim>
//ManifoldType manifoldType( const NodeManifold<dim>* const );


} // end csmp

#endif

