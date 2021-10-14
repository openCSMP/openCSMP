#include "Element.h"
#include "Node.h"

#include "Exception.h"

#include "Visitor.h"

#include <cfloat>

using namespace std;

namespace csmp {


/**
Element stub used to model an element at the model boundary. The key
information here is the boundary flag.

Default constructor associates Element with instance of
fem_manager.DefaultElement() pointing to the default Element as specified in main.
Also, vectors of pointers to nodes, constraint points and neighbor elements
are constructed.

@section application Application

Elements are constructed inside the MeshManager. There are 2 methods
for VSets and VSets respectively. VSets (a class that will be
phased out in 1999) cannot hold different element types in one mesh. In
this case the Element default constructor will initialize the FiniteElement
pointer to zero.
*/
template<size_t dim>
Element<dim>::Element( BOX_BOUNDARY bflag )
  : idx_( UINT_MAX ),
    material_id_(UNSPECIFIED)
{
}



/**
Constructor associates Element with instance of specific finite element.
This initialized the bridge pattern. Also, vectors of pointers to nodes,
constraint points and neighbor elements are constructed.

@section input Input Arguments

Finite element, for example LinearTriangle. The element type determines how
many nodes and constraint points the Element will have and if affects the
element dimensionality which can be queried with Dim().

@section application Application

Elements are constructed inside the MeshManager. There are 2 methods
for VSets and VSets respectively. VSets (a class that will be
phased out in 1999) cannot hold different element types in one mesh. In
this case the Element default constructor will initialize the FiniteElement
reference to FiniteElementManager. In the case of VSets, specific
element types can be stored and this constructor is used to reference
them properly.
*/
template<size_t dim>
Element<dim>::Element( csmp::FiniteElement* f )
  : FiniteElementPolicy<dim, csmp::Element>( f ),
  idx_( UINT_MAX ),
  elmt_connector_( f->Neighbors(), nullptr ),
  node_connector_( f->Nodes(), nullptr ),
  material_id_(UNSPECIFIED)
{
  elmt_connector_.resize( f->Neighbors(), nullptr );
  node_connector_.resize( f->Nodes(), nullptr );
}



template<size_t dim>
Element<dim>::Element( csmp::FiniteElement* f,
                       const csmp::FiniteVolumeStencil<dim>* fvs )
  : FiniteElementPolicy<dim, csmp::Element>( f ),
  FiniteVolumePolicy<dim, ::csmp::Element>( fvs ),
  idx_( UINT_MAX ),
  elmt_connector_( f->Neighbors(), nullptr ),
  node_connector_( f->Nodes(), nullptr ),
  material_id_(UNSPECIFIED)
{
}



template<size_t dim>
Element<dim>::Element( csmp::FiniteElement* f,
                       const csmp::FiniteVolumeStencil<dim>* fvs,
                       const LocalVariables& ep,
                       const IntegrationPointVariables& cp )

  : FiniteElementPolicy<dim, csmp::Element>( f ),
    FiniteVolumePolicy<dim, ::csmp::Element>( fvs ),
    idx_( UINT_MAX ),
    elmt_connector_( f->Neighbors(), nullptr ),
    node_connector_( f->Nodes(), nullptr ),
    material_id_(UNSPECIFIED)
{
  if ( this->UsesLocalCoordinates() )
    this->ResizePropertyStorage( ep, cp );
  else
    this->ResizePropertyStorage( ep );
}



/**
For model reconstruction from binary file

*/
template<size_t dim>
Element<dim>::Element( size_t idx,
                       csmp::FiniteElement* f,
                       const FiniteVolumeStencil<dim>* s,
                       const LocalVariables& ep,
                       const IntegrationPointVariables& cp,
                       int32 material )

  : FiniteElementPolicy<dim, csmp::Element>( f ),
    FiniteVolumePolicy<dim, ::csmp::Element>( s ),
    idx_( idx ),
    elmt_connector_( f->Neighbors(), nullptr ),
    node_connector_( f->Nodes(), nullptr ),
    material_id_(material)
{
  if ( this->UsesLocalCoordinates() )
    this->ResizePropertyStorage( ep, cp );
  else
    this->ResizePropertyStorage( ep );
}




template<size_t dim>
Element<dim>::Element( const Element<dim>& el )
  : FiniteElementPolicy<dim, csmp::Element>( el.FE() ),
    FiniteVolumePolicy<dim, csmp::Element>( el.FV() ),
    idx_( el.idx_ ),
    elmt_connector_( el.elmt_connector_ ), // the pointers point to the same elements as for the original element
    node_connector_( el.node_connector_ ),  
    material_id_(el.material_id_)
{
  assert( !node_connector_.empty() /* detected unitialized element*/ );
  assert( !elmt_connector_.empty() /* detected unitialized element*/ );
  // variable storage: call of initialization function
  this->LVS( el.LVS() );
}



/// move constructor
template<size_t dim>
Element<dim>::Element( Element<dim>&& el )
  : FiniteElementPolicy<dim, csmp::Element>( move( el.FE() ) ),
    FiniteVolumePolicy<dim, csmp::Element>( move( el.FV() ) ),
    idx_( move( el.idx_ ) ),
    elmt_connector_( move( el.elmt_connector_ ) ),
    node_connector_( move( el.node_connector_ ) ),
    material_id_(el.material_id_)
{
  this->LVS( move( el.LVS() ) );
  el.AssignFiniteElementNullPtr();
  el.AssignFiniteVolumeNullPtr();
}




template<size_t dim>
Element<dim>::~Element()
 {
    // disconnecting the neighbor elements that are connected to this element
    for ( auto it : elmt_connector_ )
      if ( it != nullptr )
        for ( auto nit : it->elmt_connector_ )
          if ( nit == this ) {
               nit = nullptr;
               break;
            }
    // disconnecting the element from its nodes and neighbors
    for ( auto& it : elmt_connector_ ) it = nullptr;
    for ( auto& it : node_connector_ ) it = nullptr;
 }



template<size_t dim>
Element<dim>& Element<dim>::operator=( const Element<dim>& el )
{
  if ( &el != this ) {
    if ( el.FE() ) FiniteElementPolicy<dim, csmp::Element>::Assign( el.FE() );
    if ( el.FV() ) FiniteVolumePolicy<dim, csmp::Element>::AssignFiniteVolume( el.FV() );
    idx_ = el.idx_;
    elmt_connector_ = el.elmt_connector_;
    node_connector_ = el.node_connector_;
    material_id_    = el.material_id_;
    this->LVS( el.LVS() );
  }
  return *this;
}



/**
@note a temporary variable cannot be equivalent to lvalue!
*/
template<size_t dim>
Element<dim>& Element<dim>::operator=( Element<dim>&& el )
{
  // should never happen because a temporary variable cannot be an lvalue
  assert( &el != this );

  if ( el.FE() ) FiniteElementPolicy<dim, csmp::Element>::Assign( move( el.FE() ) );
  if ( el.FV() ) FiniteVolumePolicy<dim, csmp::Element>::AssignFiniteVolume( move( el.FV() ) );

  idx_ = move( el.idx_ );
  elmt_connector_ = move( el.elmt_connector_ );
  node_connector_ = move( el.node_connector_ );
  material_id_    = move(el.material_id_ );
  this->LVS( move( el.LVS() ) );

  el.AssignFiniteElementNullPtr();
  el.AssignFiniteVolumeNullPtr();

  return *this;
}



/// Roman, 2014
/// WARNING: this operator is used specifically in the process of creation of particular ??? Region ???
/// Therefore only important infromation for that process is taken into account in order to distinguish two Element's.
/// That must be FE_Type and attached Nodes
///   SKM revised 2021
template<size_t dim>
bool  Element<dim>::operator==( const Element<dim>& el )
{
  if ( &el != this )
    {
       if ( this->FE_Type() != el.FE_Type() ) return false;
       if ( node_connector_.size() != el.node_connector_.size() ) return false;
       // ignored: if ( idx_ == el.idx_ ) return true;
       set<Node<dim>*> nodes1( node_connector_.begin(), node_connector_.end() );
       set<Node<dim>*> nodes2( el.node_connector_.begin(), el.node_connector_.end() );
       vector<Node<dim>*> nodes_intersect;
       set_intersection( nodes1.begin(), nodes1.end(),
                         nodes2.begin(), nodes2.end(),
                         back_inserter( nodes_intersect ) );
        if ( nodes_intersect.size() == node_connector_.size() )
          return true;
    }
  return true;
}


/**
  less_than<> predicate for storage of Elements in STL container objects, including equal comparisons.
    @author SKM
    @date 6/9/2021
*/
template<size_t dim>
bool  Element<dim>::operator<( const Element<dim>& el )
{
  if ( &el != this )
    {
       Point<dim> barycentre(BaryCenter());
       Point<dim> barycentre_el(el.BaryCenter());
       return barycentre < barycentre_el;
    }
  return false;
}


// privatized to avoid use outside of MeshManager
template<size_t dim>
void* Element<dim>::operator new( size_t size )
  {
//      std::cout<< "\nElement<"<< dim <<">: called overloaded new operator.\n";
      //void * p = malloc(size); will also work fine
      return ::operator new(size);
  }
 

// privatized to avoid use outside of MeshManager
template<size_t dim>
void Element<dim>::operator delete( void* p )
  {
//     std::cout<< "\nElement<"<< dim <<">: called overloaded delete operator.\n";
     free(p);
     p = nullptr;
  }





/**
Visitor subclassed objects may be passed to the Element object to gain
access to its public interface. This then allows the visitor to query
the Elements methods for the data it needs for its own computations. Apart from just
visiting each Element sequentially, a Visitor can target the Elements
neigbor elements to which it simultaneously has access. This allows to
implement a mesh traversal technique like streamline routing.

@section implementation Implementation

The visitor is queried for its application Target. Depending on the result
it is applied to this Element or passed on to all of its Nodes,
IntegrationPoints or Neighbors.

@section application Application

The Accept method is used, for instance, by the TranportVisitor class.

@param vis A reference to a Visitor subclass.
*/
template<size_t dim>
void Element<dim>::Accept( csmp::Visitor<dim>& vis )
{
  if ( vis.ApplicationTarget() == ELEMENT ) {
    vis.Visit( this );
    return;
  }
  if ( vis.ApplicationTarget() == NODE ) {
    for ( typename vector<csmp::Node<dim>*>::iterator
          nit = node_connector_.begin(); nit != node_connector_.end(); nit++ ) (*nit)->Accept( vis );
      return;
  }
  throw logic_error( "Element<dim>::Accept: target of visitation unresolved." );

} // end Accept    


  /// number of nodes of this element
template<size_t dim>
size_t  Element<dim>::Nodes() const
{
  return node_connector_.size();
}

/// number of equidimensional neighbor elements of this element (not necessarily connected)
template<size_t dim>
size_t  Element<dim>::Neighbors() const
{
  return elmt_connector_.size();
}

/// number of equidimensional neighbor elements of this element (necessarily connected)
template<size_t dim>
size_t  Element<dim>::ConnectedNeighbors() const
{
  size_t nulls( 0 );
  for ( auto f : elmt_connector_ )
    if ( !f ) nulls++;
  return (elmt_connector_.size() - nulls);
}

/// number of faces (side-surfaces) of the current element; for each element face, there can be a neighbor
template<size_t dim>
size_t  Element<dim>::Faces() const
{
  return elmt_connector_.size();
}

// ITERATORS

/// iterator to the element nodes
template<size_t dim>
typename std::vector<csmp::Node<dim>*>::iterator  Element<dim>::NodesBegin()
{
  return node_connector_.begin();
}

template<size_t dim>
typename std::vector<csmp::Node<dim>*>::iterator  Element<dim>::NodesEnd()
{
  return node_connector_.end();
}

template<size_t dim>
typename std::vector<Element<dim>*>::iterator  Element<dim>::NeighborsBegin()
{
  return elmt_connector_.begin();
}

template<size_t dim>
typename std::vector<Element<dim>*>::iterator  Element<dim>::NeighborsEnd()
{
  return elmt_connector_.end();
}

template<size_t dim>
typename std::vector<csmp::Node<dim>*>::const_iterator  Element<dim>::NodesBegin() const
{
  return node_connector_.begin();
}

template<size_t dim>
typename std::vector<csmp::Node<dim>*>::const_iterator  Element<dim>::NodesEnd() const
{
  return node_connector_.end();
}

template<size_t dim>
typename std::vector<Element<dim>*>::const_iterator  Element<dim>::NeighborsBegin() const
{
  return elmt_connector_.begin();
}

template<size_t dim>
typename std::vector<Element<dim>*>::const_iterator  Element<dim>::NeighborsEnd() const
{
  return elmt_connector_.end();
}






// CONSTRUCTION PROCESS


template<size_t dim>
void Element<dim>::Assign( size_t i, Element<dim>* const e_ptr ) // neighbor elements
{
  assert( elmt_connector_.size() == this->Neighbors() );
  assert( i < this->Neighbors() );

  elmt_connector_[i] = e_ptr;
}


/**
    unassigns the neighbor element, setting the pointer in the 'elmt_connector' vector to null
*/
template<size_t dim>
void Element<dim>::Unassign( const Element<dim>* e_ptr ) 
  {
    for ( size_t i = 0U; i < elmt_connector_.size(); ++i ) {
        if ( e_ptr == nullptr || elmt_connector_[i] == nullptr )
          continue;
        if ( e_ptr == elmt_connector_[i] ) {
            elmt_connector_[i] = nullptr;
            break;
          }
      }
  }




template<size_t dim>
void Element<dim>::Assign( size_t i, csmp::Node<dim>* const nd_ptr )
{
  assert( node_connector_.size() == this->Nodes() );
  assert( i < this->Nodes() );
  assert( nd_ptr != nullptr );

  node_connector_[i] = nd_ptr;
}


template<size_t dim>
void Element<dim>::Unassign( csmp::Node<dim>* const nd_ptr )
  {
    for ( size_t i = 0U; i < node_connector_.size(); i++ ) {
        if ( nd_ptr == nullptr || node_connector_[i] == nullptr )
          continue;
        if ( (*nd_ptr) == (*node_connector_[i]) ) {
            node_connector_[i] = nullptr;
            break;
          }
      }
  }



template<size_t dim>
void  Element<dim>::Idx( size_t idx_to_assign ) const
{
  idx_ = idx_to_assign;
}


/// unique material identifier that matches number of parent unique region
template<size_t dim>
int32 Element<dim>::Material_ID() const
 {
    return material_id_;
 }
 
 
template<size_t dim>
void Element<dim>::Material_ID( int32 id )
 {
    material_id_ = id;
 }



// ACCESSORS

template<size_t dim>
size_t   Element<dim>::Idx() const
{
  return idx_;
}

template<size_t dim>
BOX_BOUNDARY  Element<dim>::AtBoundary( size_t boundary_face ) const
{
  return atBoundary( this, boundary_face );
}

template<size_t dim>
typename  std::vector<csmp::Node<dim>*>&  Element<dim>::NodeVector()
{
  return node_connector_;
}

template<size_t dim>
typename std::vector<csmp::Element<dim>*>&  Element<dim>::NeighborElementVector()
{
  return elmt_connector_;
}






/**
The suite of methods Node(), IntegrationPoint(),
Face() and Neighbor(),
enables access of these connected objects via
const references. A NULL pointer is returned by Neighbor() if the
Element lies at the model boundary and a Neighbor does not exist.
Also constraint points do not exist unless you defined these for the
specific FiniteElement which you are using.

@section input Input Arguments

An integer from 0...n-1, where n is the number of nodes per Element. Nodes
are numbered clockwise for regular-gridded meshes (since Y points downward)
and counter-clockwise for all other meshes.

@section implementation Implementation

A range check is performed.

@section application Application

Use only when the other methods of the Element fail to give you enough
information about the target object. For instance when you want to retrieve
a Node coordinate:

@code
double64 x = (*element.N(2))->x();
@endcode

@return return A pointer to the Target object, e.g. a Node.
*/
template<size_t dim>
csmp::Node<dim>*  Element<dim>::N( size_t n ) const
{
  assert( node_connector_.size() == this->Nodes() );
  assert( n < this->Nodes() );
  return node_connector_[n];
}


/**
Returns a pointer to neighbouring element i
(see CSMP's FEM type conventions to understand the neighbor
numbering scheme).

@attention always check whether the neighbor pointer is valid
before you are trying to use it.
*/
template<size_t dim>
csmp::Element<dim>*  Element<dim>::Neighbor( size_t n ) const
{
  assert( elmt_connector_.size() == this->Neighbors() );
  assert( n < this->Neighbors() );
  return elmt_connector_[n];
}


/**
Returns a matrix with 'nodes'-rows and 'coordinate-directions' columns.
This Meschach++ matrix defines the positions of the elements nodes for
the finite-element matrix assembly. Since the number of element nodes
may vary among different elements types, the number of rows in XY may
also vary from element to element.

@param XY A DenseMatrix<DM_MIN> class object (value type fT). This matrix is dynamically
resized if necessary but must have been constructed with a finite size
before passing it to CoordinateMatrix().

@return Element The node coordinates are returned into the supplied matrix.

@section application Application

Finite-element forms of differential equations require the global node
coordinates of the element to calculate the element constribution to the
global solution matrix. If the element uses local coordinates, the global
node coordinates will still be required to compute Jacobian (coordinate-
transformation) matrix.
*/
template<size_t dim>
void  Element<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
{
  const size_t n_nodes( Nodes() );
  XY.Resize( n_nodes, dim );
  for ( size_t i = 0U; i<n_nodes; ++i )
    XY.AssignRow( i, N( i )->Coordinate() );

} // end CoordinateMatrix







/**
BaryCenter() calculates the node coordinate average for the element. This
coordinate value is equivalent to the center of gravity of the Element
type.

@section input Input Arguments

The barycentre is returned into a CSMP vector variable a reference to which
is supplied as single argument.

@section application Application

The element barycentre could be used for instance to output element material
properties from the Model as point data. This is done if you use the
Model OutputToTextFile() methods.

@return The VARIABLE_FLAG of the returned vector variable will not be changed by
BaryCentre().

@note could be madfe more

@test O.K. SKM25/8/14 after refactoring loop

*/
template<size_t dim>
Point<dim>  Element<dim>::BaryCenter() const
{
  Point<dim>    pt( N( 0U )->Coordinate() );
  const size_t  n_nodes( Nodes() );
  for ( size_t i = 1U; i<n_nodes; ++i )
    pt += N( i )->Coordinate();

  return pt / static_cast<double64>(Nodes());
}



/**
Measures the length of the element (in physical space) in a certain direction.

@param vecDirection direction in which the element is to be measured

@return Length of the element.

@sectin implementation Implementation

The main idea is to measure the height of the oriented bounding box
(oriented by the direction of vecDirection) which surrounds the element.
This is done by projecting each node onto the direction vector and calculating
the difference between the largest and smallest magnitude. The distance between both
projections will be given by the absolute value of this difference.

@section application Application

The length of the element in a certain dimension is used to weigh the error of the element,
when evaluating the quality of a certain mesh.

@todo this method could be optimised if only the corner nodes would be used

*/
template<size_t dim>
double64  Element<dim>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
{
  double64 fMinTemp( static_cast<double64>(DBL_MAX) );
  double64 fMaxTemp( static_cast<double64>(-DBL_MAX) );

  // this normalisation is necessary because the vector variable
  // being any physical quantity may have any magnitude
  const double64 fMagnitudeOfDirection( vecDirection.Length() );
  // avoid division by zero
  assert( fMagnitudeOfDirection >= numeric_limits<double64>::epsilon() );

  const size_t n_nodes( Nodes() );
  for ( size_t i = 0; i<n_nodes; ++i ) {
    // fTemp is the projection of the vector (0,0,0)-node(i) on the vector direction
    double64 fTemp( vecDirection.DotProduct( N( i )->Coordinate() ) );
    fTemp /= fMagnitudeOfDirection;

    // update minimum value
    fMinTemp = std::min( fMinTemp, fTemp );
    // update maximum value
    fMaxTemp = std::max( fMaxTemp, fTemp );
  }

  //substract magnitudes
  return fMaxTemp - fMinTemp;
}




/**
returns property values at the nodes
*/
template<size_t dim>
template< class Var>
void  Element<dim>::NodePropertyVector( const csmp::Index& idx, std::vector<Var>& V ) const
{
  if ( idx.place != NODE ) {
    std::cerr << "\nElement<" << dim;
    std::cerr << ">::NodePropertyVector: Requested property ";
    std::cerr << "is not placed on the nodes; property Index: " << std::endl;
    idx.Out();
    return;
  }

  // resizing V if necessary
  const size_t  n_nodes( Nodes() );
  V.resize( n_nodes );

  for ( size_t i = 0U; i<n_nodes; i++ )
    N( i )->Read( idx, V[i] );
}

// scalar
template void  Element<1U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>& ) const;
template void  Element<2U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>& ) const;
template void  Element<3U>::NodePropertyVector( const csmp::Index&, std::vector<ScalarVariable>& ) const;
// vector
template void  Element<1U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<1U> >& ) const;
template void  Element<2U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<2U> >& ) const;
template void  Element<3U>::NodePropertyVector( const csmp::Index&, std::vector<VectorVariable<3U> >& ) const;
// tensor
template void  Element<1U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<1U> >& ) const;
template void  Element<2U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<2U> >& ) const;
template void  Element<3U>::NodePropertyVector( const csmp::Index&, std::vector<TensorVariable<3U> >& ) const;
// array
template void  Element<1U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>& ) const;
template void  Element<2U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>& ) const;
template void  Element<3U>::NodePropertyVector( const csmp::Index&, std::vector<ArrayVariable>& ) const;
// flagged array
template void  Element<1U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>& ) const;
template void  Element<2U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>& ) const;
template void  Element<3U>::NodePropertyVector( const csmp::Index&, std::vector<FlaggedArrayVariable>& ) const;


// OUTPUT


/// prints Element internal data and those of connected objects.
template<size_t dim>
void Element<dim>::Out() const
{
  string str( parseBoundary( atBoundary(this) ) );
  cout << "\n\nElement<" << dim << ">::Out: number: " << idx_;
  cout <<"\n\tMaterial ID: "<< material_id_;
  cout << " (" << parseFiniteElementType( this->FE_Type() ) << " = ";
  if ( this->IsLineElement() )    cout << "line element";
  else if ( this->IsSurfaceElement() ) cout << "surface element";
  else if ( this->IsVolumeElement() )  cout << "volume element";
  cout << "), boundary flag: " << str << "\n";

  cout << "\n\tconnected nodes (indices : boundary flags):  ";
  for ( size_t i = 0U; i<this->Nodes(); i++ ) {
    str = parseBoundary( N( i )->AtBoundary() );
    cout << N( i )->Idx() << ":" << str << "  ";
  }
  cout << endl;

  cout << "\n\tconnected neighbors (finite element types : boundary flags):\n";
  for ( size_t i = 0U; i<this->Neighbors(); i++ )
    if ( Neighbor( i ) != nullptr ) {
      cout << "\t\t" << Neighbor( i )->Idx() << ": ";
      cout << parseFiniteElementType( Neighbor( i )->FE_Type() ) << ": ";
      str = parseBoundary( Neighbor( i )->AtBoundary(i) );
      cout << str << endl;
    }
    else cout << "\t\tnone.\n";

    // barycentre
    Point<dim>  pt( this->BaryCenter() );
    if ( dim == 1U )
      cout << "\n\tbarycentre at (xyz): " << pt[0] << endl;
    else if ( dim == 2U )
      cout << "\n\tbarycentre at (xyz): " << pt[0] << ", " << pt[1] << endl;
    else
      cout << "\n\tbarycentre at (xyz): " << pt[0] << ", " << pt[1] << ", " << pt[2] << endl;

    // length, area, volue
    const double64 volume( this->Volume() );
    if ( this->IsLineElement() ) {
      if ( volume > 0. ) cout << "\n\tlength: " << volume << endl;
      else cerr << "\n\tlength: ERROR (negative value indicates numbering problem): " << volume << endl;
    }
    if ( this->IsSurfaceElement() ) {
      if ( volume > 0. ) cout << "\n\tarea: " << volume << endl;
      else cerr << "\n\tarea: ERROR (negative value indicates numbering problem): " << volume << endl;
    }
    else if ( this->IsVolumeElement() ) {
      if ( volume > 0. ) cout << "\n\tvolume: " << volume << endl;
      else cerr << "\n\tvolume: ERROR (negative value indicates numbering problem): " << volume << endl;
    }

    // inner radius
    if ( this->IsSurfaceElement() )
      cout << "\n\tradius of inscribed circle: " << this->InnerRadius() << endl;
    else if ( this->IsVolumeElement() )
      cout << "\n\tradius of inscribed sphere: " << this->InnerRadius() << endl;

    // aspect ratio
    if ( !this->IsLineElement() ) cout << "\n\taspect ratio (b-box):   " << this->AspectRatio() << endl;

} // end Out



template class Element<1U>;
template class Element<2U>;
template class Element<3U>;

} // end namespace csmp

