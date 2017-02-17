#ifndef CSMP_ELEMENT_H
#define CSMP_ELEMENT_H

#include "LocalVariableStorage.h"
#include "FiniteElement.h"
#include "FiniteElementPolicy.h"
#include "FiniteVolumePolicy.h"
#include "ElementRemeshingTraits.h"

#include "Box.h"

namespace csmp {

struct  Index;
template<size_t> class Node;
template<size_t> class FiniteVolumeStencil;
template<size_t> class Visitor;

/**

@brief Object representation of a finite element Bridge pattern together with
FiniteElement class. 

@author S.K. Matthai
@author Stephen G. Roberts
@date 1999

@section motivation Motivation

The finite-element method uses areal/volumetric elements as smallest computational
units. Each element has 2 functions: (1) to store material properties and
nodes with associated properties; (2) To generate matrices which contribute to
a global system of algebraic equations which is solved to obtain the value
of the dependent spatially-distributed variable which typically is stored 
on the nodes.
 
 
@section design Design Intent

The Element class is a structural object that forms part of a Bridge design
pattern (pattern 151, Gamma et al. 1995). The bridge pattern is used to
allow that elements can be mass-contructed by a container class (an STL vector 
in the MeshManager), but that they have polymorphic properties at the same.
time. Also a mesh usually has many elements, but only a few element types.
This is achieved by "bridging" the Element objects to a few FiniteElement
object pointers stored in the global array globalFE_types[].

A class diagram of the Element-FiniteElement bridge pattern is shown in
the CSP3D User Guide in the developers section. 
 
 
@section Applicability

The Element class should be applicable in any type of finite-element 
computation. Higher-order elements will make use of the ability to store
retrieve integration point data.
 
 
@section structure Structure

The structure of the Element - FiniteElement bridge pattern is described
in Gamma et al. 1995 (p. 151). The element interface is broken down into
a part that echoes the interfaces of the FiniteElement base class from
which the desired element types are inherited thereby giving these inter-
faces polymorphic behaviour.  

The main operations which the element is used for are the retrieval of
property values, node coordinate matrices, and node ID number vectors which
are used in the assembly of the element equations. The element matrices are
square and theur size is equivalent to the product of the degrees of freedom
of each node and the number of nodes. 
 
 
@section partcipants Participants

The Element queries Nodes, IntegrationPoints, and the MemoryManager<dim> to
obtain the data which it needs in its computations. It cannot exist without
reference to the FiniteElement which implements half of its interface. 

The Element uses DenseMatrix<DM_MIN>, and STL vector objects to receive
and tranfer data.
 
 
@section collaborations Collaborations

Elements can collaborate with
Visitor objects through their Accept() interface. This combines the 
methods of the Visitor with the functionality of the Element. 

The main collaboration exists between Elements and Algorithms. The latter
query the Element for coordinate and property matrices and a node number
vector. The new Algorithm class queries Elements through  virtual methods
of the MathOperator base classes.
 
 
@section consequences Consequences

Elements are "finite-element-independent" mesh-building blocks. Because
variables are only accessed through the Elements but not stored there, the
mesh storage is also independent from the number of properties associated
with it, while encapsulation is maintained and the access hierarchy is
preserved.
 
 
@section implementation Implementation

The bridge between the Element and the FiniteElement is implemented as
a reference. It is initialized either through the default constructor
which accesses the globalFE_types[] array, or, through a constructor
which takes a reference to a finite element as argument.  

Nodes, IntegrationPoints, and Neighbors of an Element are 
accessed by pointers which are initialized through the ConnectTo...() 
interfaces. Methods with the same name interface the Element with the 
methods of the FiniteElement. The Element methods Read(), Store() and
Status() are used to access element properties. Alternatively, to 
property output into DenseMatrix<DM_MIN> and vector objects 
you can also get STL vectors
of scalar, vector, and tensor properties.
 
 
@section application Application Examples

The Visit(Element* e) method in "TransportVisitor.h" illustrates nicely
the functionality of the Element class. In this method element and node
properties and node property averages are obtained and written back to
the mesh.

@section IMPORTANT - FOR DEVELOPERS

Consider that users of the policies may also want to work on other 
incarnations of the Element, in particular ones with less policies.
This only works when the policies only depend on the functionality
of the element, but not on each other!

As a key difference between policies and traits, the latter have state.
Thus, the finite-element policy own the fe pointer and the fv policy ons the fv pointer
which is key to the design.

For the current design, however, this is only partially true because the 
FiniteVolumePolicy depends on finite element functionality.
Thus, one can have an element with FEM but without FVM, but not vice versa.

@todo (3) Write faster code to determine whether a Point is contained in a certain element (A)

*/
template<size_t dim>
class Element : public ElementRemeshingTraits<dim,Element>,    ///< TODO: @todo is this really needed? - else deprecate
                public FiniteElementPolicy<dim,Element>,
                public FiniteVolumePolicy<dim,Element>,        
                public LocalVariableStorage<dim,Element<dim> > ///< TODO: @todo fix template - template parameter
  {
  public:
    explicit Element( BOX_BOUNDARY bflag=IRREGULAR );
    
    /// constructor for testing element in isolation
    explicit Element( FiniteElement* );
    
    /// constructor for testing element in isolation
    Element( FiniteElement*,
             FiniteVolumeStencil<dim>* );
    
    Element( FiniteElement*, 
             FiniteVolumeStencil<dim>*, 
             const LocalVariables& element_props, 
             const IntegrationPointVariables& integration_point_props );
    
    /// for model reconstruction from CSMP native binary file
    Element( size_t index,
             FiniteElement*,
             const LocalVariables& element_props, 
             const IntegrationPointVariables& integration_point_props,
             BOX_BOUNDARY=NOT );
             
    Element( const Element& );
//    Element( Element&& );
    ~Element();

    Element&  operator=( const Element& );

    /// relation operators
    bool operator==( const Element<dim>& );

    /// Local variable storage interface
    PLACEMENT Placement() const { return ELEMENT; }

    void Accept( csmp::Visitor<dim>& );

    // ------------------------------------------------------------------------
    // Functionality of construction process
    // ------------------------------------------------------------------------

    /// (re)connect the element to its neighbors (during the model construction process or after remeshing)
    void Assign( size_t nbor, Element<dim>* const );
    
    /// (re)connect the element to its nodes (during the model construction process or after remeshing / split boundary creation)
    void Assign( size_t node, Node<dim>* const );

    // ------------------------------------------------------------------------
    // Member access
    // ------------------------------------------------------------------------

    /// number of nodes of this element
    size_t  Nodes() const     { return node_connector_.size(); };
    
    /// number of equidimensional neighbor elements of this element (not necessarily connected)
    size_t  Neighbors() const { return elmt_connector_.size(); };
    
    /// number of faces (side-surfaces) of the current element; for each element face, there can be a neighbor
    size_t  Faces() const { return elmt_connector_.size(); };

    typename std::vector<csmp::Node<dim>*>::iterator            NodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator            NodesEnd();
    typename std::vector<csmp::Element<dim>*>::iterator         NeighborsBegin();
    typename std::vector<csmp::Element<dim>*>::iterator         NeighborsEnd();

    typename std::vector<csmp::Node<dim>*>::const_iterator      NodesBegin()     const;
    typename std::vector<csmp::Node<dim>*>::const_iterator      NodesEnd()       const;
    typename std::vector<csmp::Element<dim>*>::const_iterator   NeighborsBegin() const;
    typename std::vector<csmp::Element<dim>*>::const_iterator   NeighborsEnd()   const;

    typename  std::vector<csmp::Node<dim>*>&                    NodeVector();
    typename  std::vector<csmp::Element<dim>*>&                 NeighborElementVector();

    /// accessor of the nodes of the current finite element
    csmp::Node<dim>*     N( size_t n_local ) const;
    
    /// accessor of the equidimensional neighbor elements of the current element (volume->volume, surface->surfaces element etc.)
    csmp::Element<dim>*  Neighbor( size_t )  const;

    /// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
    void         Idx( size_t ) const;
    size_t       Idx() const;
    
    /// is element located at an outside or internal model boundary; if it shares a face with a boundary, this is true
    void         AtBoundary( BOX_BOUNDARY b );
    BOX_BOUNDARY AtBoundary() const;

    // ------------------------------------------------------------------------
    // Functionality
    // ------------------------------------------------------------------------

    /// returns a vector of the property of interest discretized on the node
    template<class Var>
    void        NodePropertyVector( const csmp::Index&, std::vector<Var>& ) const;

    /// inputs node coordinates into supplied matrix
    void        NodeCoordinateMatrix( DenseMatrix<DM_MIN>& ) const;

    /// the centre of gravity of the elemt
    Point<dim>  BaryCenter() const;

    /// projects node points onto line returning max distance between them; vec direction can have any length
    double64    LengthInDirection( const VectorVariable<dim>& vecDirection ) const;
    

    // ------------------------------------------------------------------------
    // Screen Output
    // ------------------------------------------------------------------------

    void Out() const;

 private:

    // ------------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------------

    mutable size_t                 idx_;
    std::vector<Element<dim>*>     elmt_connector_; ///< neighbors
    std::vector<csmp::Node<dim>*>  node_connector_; ///< nodes
    BOX_BOUNDARY                   at_boundary_;
};


} // csmp

#endif



