#ifndef CSMP_NODE_MANIFOLD_H
#define CSMP_NODE_MANIFOLD_H

#include "CSMP_definitions.h"


namespace csmp
{

template<size_t dim> class Node;
class Index;

/**
* @brief: a class contains pointers to coincident nodes at split-boundary
*/

enum PARENT_GEOMETRIC_ENTITY { AT_POINT, AT_EDGE, AT_SURFACE, AT_EDGE_INTERSECTION, NOT_SPECIFIED };

template<size_t dim>
class NodeManifold {
    public:
      NodeManifold();

      NodeManifold( const std::vector<Node<dim>*>& );

      NodeManifold( const std::vector<Node<dim>*>&, PARENT_GEOMETRIC_ENTITY );

      NodeManifold( typename std::vector<Node<dim>*>::iterator nodesBegin,
                    typename std::vector<Node<dim>*>::iterator nodesEnd );

      NodeManifold( typename std::vector<Node<dim>*>::iterator nodesBegin,
                    typename std::vector<Node<dim>*>::iterator nodesEnd,
                    PARENT_GEOMETRIC_ENTITY );

      /// destructor
      ~NodeManifold();
      
      /// copy constructor
      NodeManifold( const NodeManifold& );
      
      /// asignment constructor
      NodeManifold( NodeManifold&& );
      
      /// assigment operator
      NodeManifold& operator=( const NodeManifold& );
      NodeManifold& operator=( NodeManifold&& );
      
      /// asscending sort
      void SortByVariableIndex(const Index&);

      /// number of nodes
      size_t Nodes() const;

      /// access to nodes
      Node<dim>* N(const size_t&);

      PARENT_GEOMETRIC_ENTITY ParentGeometry(){return parent_geometry_;};

      bool Add( Node<dim>* );

      bool Delete( Node<dim>* );

      //bool Updated(){return updated_;}
      //void Updated(bool flag) {updated_ = flag;}

      /// print out information
      void Out();

    private:
      std::vector<Node<dim>*> nodes_;
      PARENT_GEOMETRIC_ENTITY parent_geometry_;
      //bool updated_ = false;
};
  
} // end csmp

#endif

