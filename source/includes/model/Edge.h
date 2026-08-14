//
//  Edge.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 19/10/2022.
//

#ifndef CSMP_EDGE_H
#define CSMP_EDGE_H

#include "LocalVariableStorage.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Face;
template<uint32_t> class Visitor;
class FiniteElementManager;
template<uint32_t> class FiniteVolumeStencilManager;

/**
    Line "face" that knows its surface Face neighbors.
    Supports computations on the boundaries of 3D models.
    Variable placement is FACE.
    Can be a manifold.
    
    @todo design issue: edges should only exist only in 3D models because, only there, they can have Face neighbors.
    @todo remove templatisation, making this a 3D only entity
*/
template<uint32_t dim>
class Edge : public FiniteElementPolicy<dim,Edge>,
             public FiniteVolumePolicy<dim,Edge>,
             public LocalVariableStorage<dim,Edge> {
  public:
    /// custom constructor, finds face ids and assigns nodes automatically
    Edge( FiniteElement*,
          const FiniteVolumeStencil<dim>*,
          Face<dim>* const inner_parent,
          Face<dim>* const outer_parent, ///< if any
          const LocalVariables&,
          const IntegrationPointVariables& );
          
    ~Edge() {}
  
    /// Local variable storage interface
    PLACEMENT Placement() const { return FACE; }

    uint32_t  Nodes() const { return static_cast<uint32_t>(node_connector_.size()); }
    uint32_t  Neighbors() const { return 2U; }
	  uint32_t  ConnectedNeighbors() const;
    
    /// sides of Face object by analogy with Element
    uint32_t  Faces() const { return 2U; }

    /// only constant iterators are provided because the user is not supposed to change the node pr neighbor connectivity (done by MeshManager)
    typename std::vector<Node<dim>*>::const_iterator NodesBegin() const;
    typename std::vector<Node<dim>*>::const_iterator NodesEnd() const;

    /// to apply visitors whose application level is Boundary and target is Edge
    void Accept( Visitor<dim>& );

    /// access the nodes that are connected to the Edge
    Node<dim>* const N( uint32_t n_local ) const;
  
    /// access the neighbor faces of this face
    Edge<dim>* const Neighbor( uint32_t ) const;

    /// initialisation of  nodes x dim matrix with the  coordinates of the Edge's nodes
    void NodeCoordinateMatrix( DenseMatrix<DM_MIN>& ) const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void   Idx( size_t ) const;
    size_t Idx() const;

    /// higher-dimensional Face located on side opposite to where the unit normal points; will always be present
    Face<dim>* const InnerParent() const;
  
    /// higher-dimensional Face located on the side of the face to which the unit normal points; @attention does not exist on model boundary
    Face<dim>* const OuterParent() const;

    Face<dim>* const Parent( INTERFACE_SIDE ) const;

    /// returns which Face of the higher dimensional inner neighbor Face this Edge shares its nodes with
    uint32_t InnerParentFaceID() const;
    uint32_t OuterParentFaceID() const;

    // ------------------------------------------------------------------------
    //  Face geometry operations
    // ------------------------------------------------------------------------
  
    /// returns length of the Edge; method assumes same role as Volume() for the element
    double      Length() const;

    /// the centre of gravity of the element
    Point<dim>  BaryCenter() const;
    
    void Out() const;

  private:
    Face<dim>*               innerParent_;    ///< higher-dimensional neighbor in opposite direction of unit normal (always there)
    Face<dim>*               outerParent_;    ///< higher-dimensional neighbor element in direction of interface normal
    std::vector<Node<dim>*>  node_connector_; ///< pointers to the nodes of the face
    mutable size_t  idx_;
    uint32_t   inner_parent_face_id_ = std::numeric_limits<uint32_t>::max(); ///< face number of inside higher-dimensional parent Face
    uint32_t   outer_parent_face_id_ = std::numeric_limits<uint32_t>::max(); ///< face number of outside higher-dimensional parent Face
    Edge<dim>* nbor_edge0_=nullptr;
    Edge<dim>* nbor_edge1_=nullptr;  ///< the neighbors of the edge
};

} // end csmp

#endif /* CSMP_EDGE_H */
