#include "Element.h"
#include "Node.h"

#include "Exception.h"

#include "Visitor.h"

#include <cfloat>

using namespace std;

namespace csmp {



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
template<uint32_t dim>
Element<dim>::Element( csmp::FiniteElement* f )
  : FiniteElementPolicy<dim, csmp::Element>( f ),
    idx_( numeric_limits<size_t>::max() ),
    material_id_(UNSPECIFIED),
    region_id_(UNSPECIFIED),
    elmt_connector_( f->Neighbors(), nullptr ),
    node_connector_( f->Nodes(), nullptr )
{
  assert( f != nullptr );
}



template<uint32_t dim>
Element<dim>::Element( csmp::FiniteElement* f,
                       const csmp::FiniteVolumeStencil<dim>* fvs )
  : FiniteElementPolicy<dim, csmp::Element>( f ),
    FiniteVolumePolicy<dim, ::csmp::Element>( fvs ),
    idx_( numeric_limits<size_t>::max() ),
    material_id_(UNSPECIFIED),
    region_id_(UNSPECIFIED),
    elmt_connector_( f->Neighbors(), nullptr ),
    node_connector_( f->Nodes(), nullptr )
{
   assert( f != nullptr );
}



template<uint32_t dim>
Element<dim>::Element( csmp::FiniteElement* f,
                       const csmp::FiniteVolumeStencil<dim>* fvs,
                       const LocalVariables& ep,
                       const IntegrationPointVariables& cp )

  : FiniteElementPolicy<dim, csmp::Element>( f ),
    FiniteVolumePolicy<dim, ::csmp::Element>( fvs ),
    idx_( numeric_limits<size_t>::max() ),
    material_id_(UNSPECIFIED),
    region_id_(UNSPECIFIED),
    elmt_connector_( f->Neighbors(), nullptr ),
    node_connector_( f->Nodes(), nullptr )
{
  assert( f != nullptr );
  if ( this->UsesLocalCoordinates() ) {
       this->ResizePropertyStorage( ep, cp );
    }
  else
    this->ResizePropertyStorage( ep );
}



/**
For model reconstruction from binary file

*/
template<uint32_t dim>
Element<dim>::Element( size_t idx,
                       csmp::FiniteElement* f,
                       const FiniteVolumeStencil<dim>* s,
                       const LocalVariables& ep,
                       const IntegrationPointVariables& cp,
                       int32_t material )

  : FiniteElementPolicy<dim, csmp::Element>( f ),
    FiniteVolumePolicy<dim, ::csmp::Element>( s ),
    idx_( idx ),
    material_id_(material),
    region_id_(UNSPECIFIED),
    elmt_connector_( f->Neighbors(), nullptr ),
    node_connector_( f->Nodes(), nullptr )
{
  assert( f != nullptr );
  if ( this->UsesLocalCoordinates() ) {
       this->ResizePropertyStorage( ep, cp );
    }
  else
    this->ResizePropertyStorage( ep );
}


template<uint32_t dim>
Element<dim>::Element( const Element<dim>& el )
  : FiniteElementPolicy<dim, csmp::Element>( el.FE() ),
    FiniteVolumePolicy<dim, csmp::Element>( el.FV() ),
    idx_( el.idx_ ),
    material_id_(el.material_id_),
    region_id_(el.region_id_),
    elmt_connector_( el.elmt_connector_ ),
    node_connector_( el.node_connector_ )
{
  // variable storage: call of initialization function
  this->LVS( el.LVS() );
}

/*
template<uint32_t dim>
Element<dim>::Element( Element<dim>&& el )
  : FiniteElementPolicy<dim, csmp::Element>( el.FE() ),
    FiniteVolumePolicy<dim, csmp::Element>( el.FV() ),
    // in-built types are copied
    idx_{ el.idx_ },
    material_id_{ el.material_id_},
    region_id_{ el.region_id_ },
    elmt_connector_{ el.elmt_connector_ },
    node_connector_{ el.node_connector_ }
{
   // nulling the pointers in the dying object
   el.AssignFiniteElementNullPtr();
   el.AssignFiniteVolumeNullPtr();
   
   this->LVS( std::move(el.LVS()) );
    
//  cout <<"\nElement(ctor): moved element: "<< Idx();
}
*/



template<uint32_t dim>
Element<dim>& Element<dim>::operator=( const Element<dim>& el )
{
  if ( &el != this ) {
    FiniteElementPolicy<dim, csmp::Element>::Assign( el.FE() );
    FiniteVolumePolicy<dim, csmp::Element>::AssignFiniteVolume( el.FV() );
    idx_ = el.idx_;
    elmt_connector_ = el.elmt_connector_;
    node_connector_ = el.node_connector_;
    material_id_    = el.material_id_;
    region_id_      = el.region_id_;
    this->LVS( el.LVS() );
  }
  return *this;
}


/*
template<uint32_t dim>
Element<dim>& Element<dim>::operator=( Element<dim>&& el )
{
  // should never happen because a temporary variable cannot be an lvalue
  assert( &el != this );

  FiniteElementPolicy<dim, csmp::Element>::Assign( el.FE() );
  FiniteVolumePolicy<dim, csmp::Element>::AssignFiniteVolume( el.FV() );

  idx_            = el.idx_;
  elmt_connector_ = el.elmt_connector_;
  node_connector_ = el.node_connector_;
  material_id_    = el.material_id_;
  region_id_      = el.region_id_;
  
  this->LVS( std::move( el.LVS() ) );

//  cerr <<"\nElement: move-assigned element: "<< Idx();

  return *this;
}
*/


template<uint32_t dim>
bool  Element<dim>::operator==( const Element<dim>& el ) const
{
  if ( &el != this )
    {
       if ( this->FE_Type() != el.FE_Type() ) return false;
       if ( this->FV()      != el.FV() ) return false;
       if ( node_connector_ != el.node_connector_ ) return false;
       if ( elmt_connector_ != el.elmt_connector_ ) return false;
       if ( material_id_    != el.material_id_ ) return false;
       if ( region_id_      != el.region_id_ ) return false;
       // ignored mutable index: if ( idx_ == el.idx_ ) return true;
    }
  return true;
}


/**
  less_than<> predicate for storage of Elements in STL container objects, including equal comparisons.
  
  @note Uses comparitor on barycentre point object to order the elements.
  
    @author SKM
    @date 6/9/2021
*/
template<uint32_t dim>
bool  Element<dim>::operator<( const Element<dim>& el ) const
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
template<uint32_t dim>
void* Element<dim>::operator new( size_t size )
  {
//      std::cout<< "\nElement<"<< dim <<">: called overloaded new operator.\n";
      //void * p = malloc(size); will also work fine
      return ::operator new(size);
  }
 

// privatized to avoid use outside of MeshManager
template<uint32_t dim>
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
template<uint32_t dim>
void Element<dim>::Accept( csmp::Visitor<dim>& vis )
{
  if ( vis.ApplicationTarget() == ELEMENT ) {
    vis.Visit( this );
    return;
  }
  if ( vis.ApplicationTarget() == NODE ) {
    for ( auto& nit : node_connector_ ) nit->Accept( vis );
      return;
  }
  throw logic_error( "Element<dim>::Accept: target of visitation unresolved." );

} // end Accept    


  /// number of nodes of this element
template<uint32_t dim>
uint32_t  Element<dim>::Nodes() const
{
  return static_cast<uint32_t>(node_connector_.size());
}

/// number of equidimensional neighbor elements of this element (not necessarily connected)
template<uint32_t dim>
uint32_t  Element<dim>::Neighbors() const
{
  return static_cast<uint32_t>(elmt_connector_.size());
}

/// number of equidimensional neighbor elements of this element (necessarily connected)
template<uint32_t dim>
uint32_t  Element<dim>::ConnectedNeighbors() const
{
  uint32_t nulls{0u};
  for ( auto& f : elmt_connector_ )
    if ( !f ) nulls++;
  return static_cast<uint32_t>(elmt_connector_.size() - nulls);
}

/// number of faces (side-surfaces) of the current element; for each element face, there can be a neighbor
template<uint32_t dim>
uint32_t  Element<dim>::Faces() const
{
  return static_cast<uint32_t>(elmt_connector_.size());
}

// ITERATORS

/// iterator to the element nodes

template<uint32_t dim>
typename std::vector<csmp::Node<dim>*>::const_iterator  Element<dim>::NodesBegin() const
{
  return node_connector_.begin();
}

template<uint32_t dim>
typename std::vector<csmp::Node<dim>*>::const_iterator  Element<dim>::NodesEnd() const
{
  return node_connector_.end();
}




/// assuming that the corner nodes are the first, the midside nodes the second, and the .. in the elements node set according to CSMP_FEM_conventions.pdf
template<uint32_t dim>
typename vector<csmp::Node<dim>*>::const_iterator Element<dim>::CornerNodesBegin() const
{
  return node_connector_.begin();
}

template<uint32_t dim>
typename vector<csmp::Node<dim>*>::const_iterator Element<dim>::CornerNodesEnd() const
{
  return next(node_connector_.begin(),this->FE()->CornerNodes());
}




template<uint32_t dim>
typename std::vector<Element<dim>*>::const_iterator  Element<dim>::NeighborsBegin() const
{
  return elmt_connector_.begin();
}

template<uint32_t dim>
typename std::vector<Element<dim>*>::const_iterator  Element<dim>::NeighborsEnd() const
{
  return elmt_connector_.end();
}






// CONSTRUCTION PROCESS


/**
     Assign equidimensional neighbor elements to element.
*/
template<uint32_t dim>
void Element<dim>::Assign( uint32_t i, Element<dim>* const e_ptr )
{
  assert( i < elmt_connector_.size() );
  // assert( e_ptr != nullptr ); null is a legitimate assignment if their is no neighbor
  elmt_connector_[i] = e_ptr;
}



/**
    Unassigns the neighbor element, setting the pointer in the 'elmt_connector' vector to nullptr.
    Neither the argument pointer nor the Element pointed to are modified.
    
        @return whether removal was successful

*/
template<uint32_t dim>
bool Element<dim>::Unassign( const Element<dim>* const e_ptr )
  {
    if ( e_ptr == nullptr ) return false;
    const auto n_nbors = static_cast<uint32_t>(elmt_connector_.size());
    for ( uint32_t i{0U}; i < n_nbors; ++i )
      if ( elmt_connector_[i] == e_ptr ) {
          elmt_connector_[i] = nullptr;
          return true; // a removal was made
        }
    return false;
  }

/**
      Fast version which assumes that the number of the neighbor is known
*/
template<uint32_t dim>
void Element<dim>::UnassignNeighbor( uint32_t nbor )
  {
     assert( nbor < elmt_connector_.size() );
     elmt_connector_[nbor] = nullptr;
  }


template<uint32_t dim>
void Element<dim>::Assign( uint32_t i, csmp::Node<dim>* const nd_ptr )
{
  assert( i < node_connector_.size() );
  assert( nd_ptr != nullptr );
  node_connector_[i] = nd_ptr;
}




template<uint32_t dim>
void Element<dim>::Unassign( const csmp::Node<dim>* const nd_ptr )
  {
    if ( nd_ptr == nullptr ) return;
    const auto n_nodes = static_cast<uint32_t>(node_connector_.size());
    for ( uint32_t i{0U}; i < n_nodes; i++ )
      if ( nd_ptr == node_connector_[i] ) {
          node_connector_[i] = nullptr;
          break;
        }
  }



template<uint32_t dim>
void  Element<dim>::Idx( size_t idx_to_assign ) const
{
  idx_ = idx_to_assign;
}


/// unique material identifier that matches number of parent unique region
template<uint32_t dim>
int32_t Element<dim>::Material_ID() const
 {
    return material_id_;
 }
 
 
template<uint32_t dim>
void Element<dim>::Material_ID( int32_t id )
 {
    material_id_ = id;
 }


template<uint32_t dim>
int32_t Element<dim>::Region_ID() const
 {
    return region_id_;
 }
 
 
template<uint32_t dim>
void Element<dim>::Region_ID( int32_t id )
 {
    region_id_ = id;
 }



// ACCESSORS

template<uint32_t dim>
size_t   Element<dim>::Idx() const
{
  return idx_;
}




template<uint32_t dim>
BOX_BOUNDARY  Element<dim>::AtBoundary( uint32_t boundary_face ) const
{
   return atBoundary( this, boundary_face );
}

template<uint32_t dim>
typename  std::vector<csmp::Node<dim>*>&  Element<dim>::NodeVector()
{
  return node_connector_;
}

template<uint32_t dim>
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
double x = (*element.N(2))->x();
@endcode

@return return A pointer to the Target object, e.g. a Node.
*/
/*
template<uint32_t dim>
const csmp::Node<dim>*  Element<dim>::N( uint32_t n ) const
{
  assert( n < node_connector_.size() );
  return node_connector_[n];
}
*/

template<uint32_t dim>
csmp::Node<dim>*  const Element<dim>::N( uint32_t n ) const
{
  assert( n < node_connector_.size() );
  return node_connector_[n];
}


/**
Returns a pointer to neighbouring element i
(see CSMP's FEM type conventions to understand the neighbor
numbering scheme).

@attention always check whether the neighbor pointer is valid
before you are trying to use it.
*/
template<uint32_t dim>
csmp::Element<dim>*  const Element<dim>::Neighbor( uint32_t n ) const
{
  assert( n < elmt_connector_.size() );
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

@note Element The node coordinates are returned into the supplied matrix.

@section application Application

Finite-element forms of differential equations require the global node
coordinates of the element to calculate the element constribution to the
global solution matrix. If the element uses local coordinates, the global
node coordinates will still be required to compute Jacobian (coordinate-
transformation) matrix.
*/
template<uint32_t dim>
void  Element<dim>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
{
  const uint32_t n_nodes{ Nodes() };
  XY.Resize( n_nodes, dim );
  for ( uint32_t i{0U}; i<n_nodes; ++i )
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
template<uint32_t dim>
Point<dim>  Element<dim>::BaryCenter() const
{
  Point<dim>  pt( N( 0U )->Coordinate() );
  const uint32_t  n_nodes{ Nodes() };
  for ( uint32_t i{1U}; i<n_nodes; ++i )
    pt += N( i )->Coordinate();

  return pt / static_cast<double>(Nodes());
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
template<uint32_t dim>
double  Element<dim>::LengthInDirection( const VectorVariable<dim>& vecDirection ) const
{
  double fMinTemp( static_cast<double>(DBL_MAX) );
  double fMaxTemp( static_cast<double>(-DBL_MAX) );

  // this normalisation is necessary because the vector variable
  // being any physical quantity may have any magnitude
  const double fMagnitudeOfDirection( vecDirection.Length() );
  // avoid division by zero
  assert( fMagnitudeOfDirection >= numeric_limits<double>::epsilon() );

  const auto n_nodes( Nodes() );
  for ( uint32_t i{0u}; i<n_nodes; ++i ) {
    // fTemp is the projection of the vector (0,0,0)-node(i) on the vector direction
    double fTemp( vecDirection.DotProduct( N( i )->Coordinate() ) );
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
template<uint32_t dim>
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
  const auto  n_nodes{ Nodes() };
  V.resize( n_nodes );

  for ( uint32_t i{0U}; i<n_nodes; i++ )
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
template<uint32_t dim>
void Element<dim>::Out() const
{
  cout << "\n\nElement<" << dim << ">::Out: number: " << idx_;
  cout <<"\n\tMaterial ID: "<< material_id_;
  cout << " (" << parseFiniteElementType( this->FE_Type() ) << " = ";
  if ( this->IsLine() )    cout << "line element";
  else if ( this->IsSurface() ) cout << "surface element";
  else if ( this->IsVolume() )  cout << "volume element";
  cout <<")\n";

  cout << "\n\tconnected nodes (indices : boundary flags):  ";
  string str("undefined");
  for ( uint32_t i{0U}; i<this->Nodes(); i++ ) {
    str = parseBoundary( N( i )->AtBoundary() );
    cout << N( i )->Idx() << ":" << str << "  ";
  }
  cout << endl;

  cout << "\n\tconnected neighbors (finite element types : boundary flags):\n";
  for ( uint32_t i{0U}; i<this->Neighbors(); i++ ) {
        if ( Neighbor( i ) != nullptr ) {
            cout << "\t\t"<<"elmt "<< Neighbor( i )->Idx() << ": ";
            cout << parseFiniteElementType( Neighbor( i )->FE_Type() ) << " ";
          }
        else cout << "\t\tnone ("<< parseBoundary( AtBoundary(i) ) <<").";
        cout << endl;
     }
     
    // barycentre
    Point<dim>  pt( this->BaryCenter() );
    if ( dim == 1U )
      cout << "\n\tbarycentre at (xyz): " << pt[0] << endl;
    else if ( dim == 2U )
      cout << "\n\tbarycentre at (xyz): " << pt[0] << ", " << pt[1] << endl;
    else
      cout << "\n\tbarycentre at (xyz): " << pt[0] << ", " << pt[1] << ", " << pt[2] << endl;

    // length, area, volue
    const double volume( this->Volume() );
    if ( this->IsLine() ) {
      if ( volume > 0. ) cout << "\n\tlength: " << volume << endl;
      else cerr << "\n\tlength: ERROR (negative value indicates numbering problem): " << volume << endl;
    }
    if ( this->IsSurface() ) {
      if ( volume > 0. ) cout << "\n\tarea: " << volume << endl;
      else cerr << "\n\tarea: ERROR (negative value indicates numbering problem): " << volume << endl;
    }
    else if ( this->IsVolume() ) {
      if ( volume > 0. ) cout << "\n\tvolume: " << volume << endl;
      else cerr << "\n\tvolume: ERROR (negative value indicates numbering problem): " << volume << endl;
    }

    // inner radius
    if ( this->IsSurface() )
      cout << "\n\tradius of inscribed circle: " << this->InnerRadius() << endl;
    else if ( this->IsVolume() )
      cout << "\n\tradius of inscribed sphere: " << this->InnerRadius() << endl;

    // aspect ratio
    if ( !this->IsLine() ) cout << "\n\taspect ratio (b-box):   " << this->AspectRatio() << endl;

} // end Out



template class Element<1U>;
template class Element<2U>;
template class Element<3U>;




// NON-MEMBER FUNCTIONS

/// returns the nodes that are shared between the two elements; if these are the node-set of a shared face the returned boolean is set to true
template<uint32_t dim>
pair<vector<Node<dim>*>,bool>  sharedNodes( const Element<dim>* const eptr1, const Element<dim>* const eptr2 )
 {
    if ( !eptr1 || !eptr2 ) return make_pair( vector<Node<dim>*>(), false );
 
    const bool describes_face = ( (eptr1->IsLine() && eptr2->IsLine()) ||
                                  (eptr1->IsSurface() && eptr2->IsSurface()) ||
                                  (eptr1->IsVolume() && eptr2->IsVolume()) ) ? true : false;

     // agrrrh! - we need sorted vectors, so we might as well get sorted node sets
     set<Node<dim>*> nodes1( eptr1->NodesBegin(), eptr1->NodesEnd() ), nodes2( eptr2->NodesBegin(), eptr2->NodesEnd() );
     // store the nodes shared between the first and the second element
     vector<Node<dim>*> shared_nodes;
     set_intersection( nodes1.begin(), nodes1.end(),
                       nodes2.begin(), nodes2.end(),
                       back_inserter(shared_nodes) );

     return make_pair( shared_nodes, describes_face );
 
 } // end sharedNodes

template pair<vector<Node<3>*>,bool>  sharedNodes( const Element<3>* const, const Element<3>* const );
template pair<vector<Node<2>*>,bool>  sharedNodes( const Element<2>* const, const Element<2>* const );
template pair<vector<Node<1>*>,bool>  sharedNodes( const Element<1>* const, const Element<1>* const );


} // end namespace csmp

