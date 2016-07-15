#ifndef CSMP_FACE_H
#define CSMP_FACE_H

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementTraits.h"
#include "FiniteVolumeTraits.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Visitor;
//class FiniteElementManager;

/// surface (3D) or line (2D) element connector for use at material interfaces
// SKM FIX removed FaceRemeshingTraits<dim,Face>

template<size_t dim>
class Face : public FiniteElementTraits<dim,Face>,
             public FiniteVolumeTraits<dim,Face>,
             public LocalVariableStorage<dim,Face<dim> >
{
  public:
  
    // ------------------------------------------------------------------------
    // Functionality used in Face construction process (in that order)
    // ------------------------------------------------------------------------

    /// constructs model-interior face as an exact copy of the supplied lower-dimensional element; no neighbor faces yet
    Face( const Element<dim>& dim_minus1_element,
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs model-boundary face as n-th (external) boundary face of the supplied higher dimensional parent element; no neighbor faces yet
    Face( Element<dim>& dim_dimensional_inner_parent_element,
          csmp::FiniteElement* FE_type_of_boundary_face,
          size_t n_boundary_face,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs model-edge line-element face connected with two volumetric elements at model boundary sharing its nodes
    Face( csmp::FiniteElement* FE_type_of_boundary_face,
          Element<dim>* const inner_parent,
          Element<dim>* const outer_parent,
          const std::vector<Node<dim>*>&  edge_nodes,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// constructs face as n-th (external) boundary face of the supplied higher dimensional parent element; no neighbor faces yet
    // TODO: test whether ever needed
/*
    Face( const FiniteElementManager& finiteElementManager,
          Element<dim>& dim_dimensional_inner_parent_element,
          size_t& nth_boundary_face, ///< takes target boundary face as input and returns number of discovered boundary faces as output
          const LocalVariables&,
          const IntegrationPointVariables& );
*/
    /// prefered custom constructor creates face with together with variable storage
    Face( csmp::FiniteElement*,
          csmp::FiniteVolumeStencil<dim>*,
          const LocalVariables&,
          const IntegrationPointVariables& );

    /// for reconstruction of model from binary file; with storage but without connectivity
    Face( size_t index,
          csmp::FiniteElement*,
          const LocalVariables&,
          const IntegrationPointVariables& );

    Face( const Face& );
    Face( Face&& );

    ~Face();

    /// @attention because of the pointers, this assignment makes sense only in the rarest cases
    Face& operator=( const Face<dim>& );

    /// connects face to the supplied node
    void Assign( size_t node, Node<dim>* const );
  
    /// assign higher-dimensional neighbor elements to either side of face (outside is optional); needs nodes to be assigned first
    void Assign( Element<dim>* const innerElement, Element<dim>* const outerElement );
  
    /// tell face about its face neighbors
    void Assign( size_t nbor, Face<dim>* const );
  
    /// endow Face with FiniteVolume functionality
    void Assign( const FiniteVolumeStencil<dim>* const );
  
    // TODO: SKM: deprecate once new functionality is available
    /// assigns inner parent and node indices after finding if boundary face by matching the provided face nodes
//    void Assign( Element<dim>* const parent, const std::vector<Node<dim>*>& faceNodes );

    /// compares faces with one-another
    bool operator==( const Face<dim>& ) const;

    /// Local variable storage interface
    PLACEMENT Placement() const { return FACE; }
  
    /// box boundary flagging is not carried over to face

    // ------------------------------------------------------------------------
    // Member functions
    // ------------------------------------------------------------------------

    /// to apply visitors whose application level is Boundary and target is Face
    void Accept( csmp::Visitor<dim>& );

    /// helper method for remeshing purposes; @todo move to remeshing policy
    typename  std::vector<csmp::Face<dim>*>& NeighborElementVector();

    /// access the nodes that are connected to the Face
    csmp::Node<dim>*  N( size_t n_local ) const;
  
    /// access the meighbor faces of this face
    csmp::Face<dim>*  Neighbor( size_t ) const;

    /// access the finite element subclass that is associated with this face (triangle etc.)
    FiniteElement*    FE() const;
  
    /// access low-level finite volume functionality
    const FiniteVolumeStencil<dim>* FV_Stencil() const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void           Idx( size_t ) const;
    size_t         Idx() const;

    /// access the higher dimensional elements on either side of face; @attention returns nullptr if outside is not present
    Element<dim>*  Parent( INTERFACE_SIDE ) const;
  
    /// higher-dimensional element located on side opposite to where the unit normal points; will always be present
    Element<dim>*  InnerParent() const;
  
    /// higher-dimensional element located on the side of the face to which the unit normal points; @attention does not exist on model boundary
    Element<dim>*  OuterParent() const;
  
    /// returns which Face of the higher dimensional inner neighbor element this Face shares its nodes with
    size_t         InnerParentFaceNumber() const;
  
    /// returns the number of the desired node in the inner parent element of the Face
    size_t         ParentNodeNumber( size_t n_local ) const;

    // ------------------------------------------------------------------------
    // Geometry
    // ------------------------------------------------------------------------
  
    /// returns area of the face; method assumes same role as Volume() for the element
    double64       Area() const;
  
    /// not a face-normal vector, but the shortest path between the barycenters of face and element
    void           VectorToInnerElementBaryCenter( VectorVariable<dim>& ) const;
  
    //// to compute unit normal to Face, see FiniteElementTraits for following methods
    // Point<dim> UnitNormal() const;
    // void       UnitNormal( std::vector<double64>& nrml ) const;
    // void       UnitNormal( VectorVariable<dim>& nrml ) const;

    // ------------------------------------------------------------------------
    // Screen Output
    // ------------------------------------------------------------------------

    /// prints state of this object
    void  Out() const;

  private:
    Face();
  
    // TODO: SKM: deprecate this inefficient method
    void AssignFaceID( const std::vector<Node<dim>*>& faceNodes );
    //void CheckNodeOrderingAccordingToUnitNormalOrientation();
  
    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t                          idx_;
    csmp::FiniteElement*                    fptr_;
    const csmp::FiniteVolumeStencil<dim>*   fvptr_;
  
    std::vector<Node<dim>*>  node_connector_;  ///< pointers to the nodes of the face
    std::vector<Face<dim>*>  face_connector_;  ///< the (equidimensional) neighbors of the face
    // not references or constant pointers because these may need to change during remeshing
    Element<dim>*            innerParent_;     ///< higher-dimensional neighbor in opposite direction of unit normal (always there)
    Element<dim>*            outerParent_;     ///< (optional) higher-dimensional neighbor in direction of unit normal
};



/**
@class Face Face "main_library/Face.h"

@author P. Lang
@author S.K. Matthai
@author R. Manasipov
@date 2010,2014,2016

Faces are lower-dimensional element objects with a unit normal that point
from their inner to their outer higher-dimensional parent neighbors.

The inside of a Face is classified as the side opposite to the direction into which the
unit normal of the face points to.

Like Element objects, Face objects contain pointers to their equidimensional
neigbor elements, but - in addition - they also have pointers to their
higher-dimensional parent elements that lie on either side of them.

If a Face lies on a model boundary, there is no outside parent element.
Thus, when using face objects, one can only rely on the presence of the 
innerParent element.

@attention the indices of the nodes stored in parent_element_node_ids_.
It contains the local parent element node numbers that match the nodes of the face
as the inner parent element face pointers will always be defined. Since these nodes are
shared with the Face and the outer parent element, only one such container is needed.

@note By contrast with Elements, Face and InterFace objects are not registered as parents
elements within Node objects.

@note In Box shaped models, you can also rely on the BOX_BOUNDARY flag to determine whether 
a Face is located on an external model boundary: when flagged INTERNAL the face lies inside
the model.

@note The property placement FACE is unique and to recuperate ELEMENT properties when
working with a Face one needs to call on its higher-dimensional parents.

@section motivation Motivation

To represent internal boundaries / material interfaces in an mesh.
Faces permit to implement algorithms that couple different domains of a model together
efficiently.

Algorithms that require a special treatment of boundary elements can also be made more
efficiently and written more transparently by braking them up into interior and
boundary parts, rather than writing iffy code that checks each element
whether it is at a boundary or not.

*/


} // csmp

#endif



