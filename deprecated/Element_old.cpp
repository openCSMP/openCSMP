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
Element<dim>::Element( BOX_BOUNDARY bflag  )
  : at_boundary_(bflag),
    idx_(UINT_MAX),
    fptr_(0),
    fvptr_(0)
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
  : at_boundary_(NOT),
    idx_(UINT_MAX),
    fptr_(f),
    fvptr_(0),
    elmt_connector_(f->Neighbors(),nullptr),
    node_connector_(f->Nodes(),nullptr)
 {
 }

template<size_t dim>
Element<dim>::Element( csmp::FiniteElement* f, 
                       csmp::FiniteVolumeStencil<dim>* fvs )	
  : at_boundary_(NOT),
    idx_(UINT_MAX),
    fptr_(f),
    fvptr_(fvs),
    elmt_connector_(f->Neighbors(),nullptr),
    node_connector_(f->Nodes(),nullptr)
 {
 }                 



template<size_t dim>
Element<dim>::Element( csmp::FiniteElement* f,
                       csmp::FiniteVolumeStencil<dim>* fvs,
                       const LocalVariables& ep,
                       const IntegrationPointVariables& cp )
 
  : at_boundary_(NOT),
    idx_(UINT_MAX),
    fptr_(f),
    fvptr_(fvs),
    elmt_connector_(f->Neighbors(),nullptr),
    node_connector_(f->Nodes(),nullptr)
 {
    // variable storage is resized here because the
    // finite element pointer must be initialised first
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
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
                       const LocalVariables& ep,
                       const IntegrationPointVariables& cp,
                       BOX_BOUNDARY boundary_flag )
 
  : idx_(idx),
    fptr_(f),
    fvptr_(nullptr),
    elmt_connector_(f->Neighbors(),nullptr),
    node_connector_(f->Nodes(),nullptr),
    at_boundary_(boundary_flag)
 {
    assert( fptr_ );
    if ( fptr_->UsesLocalCoordinates() )
        this->ResizePropertyStorage( ep, cp );
    else
        this->ResizePropertyStorage( ep );
 }




template<size_t dim>
Element<dim>::Element( const Element<dim>& el )  
 : at_boundary_   (el.at_boundary_   ),
   idx_           (el.idx_           ),
   fptr_          (el.fptr_          ),
   fvptr_         (el.fvptr_         ),
   elmt_connector_(el.elmt_connector_), // watch out where the pointers point to
   node_connector_(el.node_connector_)  // watch out where the pointers point to
 {
   assert( fptr_ != nullptr /* detected unitialized element*/ );
   assert( !node_connector_.empty() /* detected unitialized element*/ );
   assert( !elmt_connector_.empty() /* detected unitialized element*/ );
   // variable storage: call of initialization function
   this->LVS( el.LVS() );
 }



/// move constructor
template<size_t dim>
Element<dim>::Element( Element<dim>&& el )
 : at_boundary_   {el.at_boundary_   },
   idx_           {el.idx_           },
   fptr_          {el.fptr_          },
   fvptr_         {el.fvptr_         },
   elmt_connector_{el.elmt_connector_},
   node_connector_{el.node_connector_}
 {
   this->LVS( move(el.LVS()) );
   el.fptr_  = nullptr;
   el.fvptr_ = nullptr;
 }


template<size_t dim>
Element<dim>::~Element()
 {
//    if ( idx_ == UINT_MAX ) cerr <<"x ";
//    else cerr << idx_ <<" ";
 }



template<size_t dim>
Element<dim>& Element<dim>::operator=( const Element<dim>& el )
 {
    if ( &el != this ) {
        at_boundary_    = el.at_boundary_;
        idx_            = el.idx_;
        fptr_           = el.fptr_;
        fvptr_          = el.fvptr_;
        elmt_connector_ = el.elmt_connector_;
        node_connector_ = el.node_connector_;
        this->LVS( el.LVS() );
      }
    return *this;
 }



/// Roman, 2014
/// WARNING: this operator is used specifically in the process of creation of particular Region.
/// Therefore only important infromation for that process is taken into account in order to distinguish two Element's.
/// That must be FE_Type and attached Nodes
// TODO: check whether this is still needed after refactoring
template<size_t dim>
bool  Element<dim>::operator==( const Element<dim>& el )
 {
    if ( &el != this )
    {
        if( this->FE_Type() == el.FE_Type() )
        {
            if( node_connector_.size() == el.node_connector_.size() )
            {
                if( node_connector_.empty() )
                {
                    if( idx_ == el.idx_ )
                        return true;
                    return false;
                }

                std::set<Node<dim>*> nodes1(node_connector_.begin(), node_connector_.end());
                std::set<Node<dim>*> nodes2(el.node_connector_.begin(), el.node_connector_.end());
                std::vector<Node<dim>*> nodes_intersect;
                std::set_intersection( nodes1.begin(), nodes1.end(),
                                       nodes2.begin(), nodes2.end(),
                                       std::back_inserter(nodes_intersect) );
                if( nodes_intersect.size() == node_connector_.size() )
                    return true;

                return false;
            }
            return false;
        }
        return false;
    }
    return true;
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
         vis.Visit(this);
         return;
      }
    if ( vis.ApplicationTarget() == NODE ) {
         for ( typename vector<csmp::Node<dim>*>::iterator
               nit=node_connector_.begin(); nit!=node_connector_.end(); nit++ ) (*nit)->Accept( vis );
         return;
      }
    throw logic_error("Element<dim>::Accept: target of visitation unresolved.");

 } // end Accept    


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
void Element<dim>::Assign( const csmp::FiniteVolumeStencil<dim>* const stencil_ptr )
 {
    assert( stencil_ptr != NULL );
    fvptr_ = stencil_ptr;
 }

template<size_t dim>
void Element<dim>::Assign(FiniteElement * fem_ptr )
 {
    assert( fem_ptr != NULL );
    fptr_ = fem_ptr;
 }

template<size_t dim>
void Element<dim>::Assign( size_t i, Element<dim>* const e_ptr ) // neighbor elements
 {
    assert( fptr_ != NULL );
    assert( elmt_connector_.size() == this->Neighbors() );
    assert( i < this->Neighbors() );

    elmt_connector_[i] = e_ptr;
 }


template<size_t dim>
void Element<dim>::Assign( size_t i, csmp::Node<dim>* const nd_ptr )
 {
    assert( fptr_ != NULL );
    assert( node_connector_.size() == this->Nodes() );
    assert( i < this->Nodes() );
    assert( nd_ptr != NULL );

    node_connector_[i] = nd_ptr;
 }


template<size_t dim>
void  Element<dim>::Idx( size_t idx_to_assign ) const
 {
    idx_ = idx_to_assign;
 }

template<size_t dim>
void  Element<dim>::AtBoundary( BOX_BOUNDARY b )
 {
    at_boundary_ = b;
 }






// ACCESSORS

template<size_t dim>
size_t   Element<dim>::Idx() const
 {
    return idx_;
 }

template<size_t dim>
BOX_BOUNDARY  Element<dim>::AtBoundary() const
 {
    return at_boundary_;
 }

// methods supporting the bridge to the subclasses of FiniteElement
template<size_t dim>
FiniteElement*  Element<dim>::FE() const
 {
    return fptr_;
 }

template<size_t dim>
const FiniteVolumeStencil<dim>*  Element<dim>::FV_Stencil() const
 {
    return fvptr_;
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














// OUTPUT


/// prints Element internal data and those of connected objects.
template<size_t dim>
void Element<dim>::Out() const
 {
    string str(parseBoundary(at_boundary_));
    cout <<"\n\nElement<"<< dim <<">::Out: number: "<< idx_;
    cout <<" ("<< parseFiniteElementType(this->FE_Type()) <<" = ";
    if      ( this->IsLineElement() )    cout <<"line element";
    else if ( this->IsSurfaceElement() ) cout <<"surface element";
    else if ( this->IsVolumeElement() )  cout <<"volume element";
    cout <<"), boundary flag: "<< str <<"\n";

    cout <<"\n\tconnected nodes (indices : boundary flags):  ";
    for ( size_t i=0U; i<this->Nodes(); i++ ) {
         str = parseBoundary(N(i)->AtBoundary());
         cout << N(i)->Idx() <<":"<< str <<"  ";
      }
    cout << endl;
    
    cout <<"\n\tconnected neighbors (finite element types : boundary flags):\n";
    for ( size_t i=0U; i<this->Neighbors(); i++ )
      if ( Neighbor(i) != NULL ) {
           cout <<"\t\t"<< Neighbor(i)->Idx() <<": ";
           cout << parseFiniteElementType( Neighbor(i)->FE_Type()) <<": ";
           str = parseBoundary(Neighbor(i)->AtBoundary());
           cout << str << endl;
        }
      else cout <<"none.\n";
   
    // barycentre
    Point<dim>  pt(this->BaryCenter());
    if ( dim == 1U )
       cout <<"\n\tbarycentre at (xyz): "<< pt[0] << endl;
    else if ( dim == 2U )                        
       cout <<"\n\tbarycentre at (xyz): "<< pt[0] <<", "<< pt[1] << endl;
    else                          
       cout <<"\n\tbarycentre at (xyz): "<< pt[0] <<", "<< pt[1] <<", "<< pt[2] << endl;

    // length, area, volue
    const double64 volume(this->Volume());
    if (  this->IsLineElement() ) {
         if ( volume > 0. ) cout <<"\n\tlength: "<< volume << endl;
         else cerr <<"\n\tlength: ERROR (negative value indicates numbering problem): "<< volume << endl;
      }
    if (  this->IsSurfaceElement() ) {
         if ( volume > 0. ) cout <<"\n\tarea: "<< volume << endl;
         else cerr <<"\n\tarea: ERROR (negative value indicates numbering problem): "<< volume << endl;
      }
    else if ( this->IsVolumeElement() ) {
         if ( volume > 0. ) cout <<"\n\tvolume: "<< volume << endl;
         else cerr <<"\n\tvolume: ERROR (negative value indicates numbering problem): "<< volume << endl;
      }
   
    // inner radius
    if (  this->IsSurfaceElement() )
      cout <<"\n\tradius of inscribed circle: "<< this->InnerRadius() << endl;
    else if ( this->IsVolumeElement() )
      cout <<"\n\tradius of inscribed sphere: "<< this->InnerRadius() << endl;
   
    // aspect ratio
    if ( !this->IsLineElement() ) cout <<"\n\taspect ratio (b-box):   "<< this->AspectRatio() << endl;
   
 } // end Out





template class Element<1U>;
template class Element<2U>;
template class Element<3U>;

} // end namespace csmp






