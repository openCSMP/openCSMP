#include "ElementRemeshingTraits.h"
#include "Element.h"

namespace csmp {

template<size_t dim, template<size_t> class SIMPLEX>
ElementRemeshingTraits<dim,SIMPLEX>::ElementRemeshingTraits()
{

}


/** Detaches the element from its neighbors and from its nodes (removes itself as parent)
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ElementRemeshingTraits<dim,SIMPLEX>::Detach()
 {
    UnassignNodes();
    UnassignNeighbors();
 }

/// TODO: @todo: SKM: questionable logic -> remove
template<size_t dim, template<size_t> class SIMPLEX>
void  ElementRemeshingTraits<dim,SIMPLEX>::UnassignNodes(  ) // unassign nodes of current element
 {
    SIMPLEX<dim>* eptr( static_cast<SIMPLEX<dim>* >(this) );

    // remove itself from the parent containers of all its nodes
    for( size_t node = 0; node < eptr->Nodes(); ++node )
        eptr->N(node)->Unassign( eptr );

    // deleting node connectivity
    eptr->NodeVector().clear();
 }

/// TODO: @todo: SKM: questionable logic -> remove
template<size_t dim, template<size_t> class SIMPLEX>
void  ElementRemeshingTraits<dim,SIMPLEX>::UnassignNeighbors(  ) // unassign neighbors of current element
 {
    SIMPLEX<dim>* eptr( static_cast<SIMPLEX<dim>* >(this) );

    // remove element from neighbor list of its neighbors
    const size_t neighbors( eptr->Neighbors() );
    for( size_t neighbor = 0; neighbor < neighbors; ++neighbor )
      if( eptr->Neighbor( neighbor ) != NULL ){
          const size_t neighbor_neighbors( eptr->Neighbor( neighbor )->Neighbors() );
          for( size_t i = 0; i < neighbor_neighbors; i++ )
              if( eptr->Neighbor( neighbor )->Neighbor( i ) != NULL )
                  if( eptr->Neighbor( neighbor )->Neighbor( i ) == eptr )
                      eptr->Neighbor( neighbor )->UnassignNeighbor( i );
          eptr->UnassignNeighbor( neighbor );
      }
    eptr->NeighborElementVector().clear();
 }

/// TODO: @todo: SKM: questionable logic -> remove
template<size_t dim, template<size_t> class SIMPLEX>
void  ElementRemeshingTraits<dim,SIMPLEX>::UnassignNeighbor( size_t i ) // unassign neighbor of current element
 {
    SIMPLEX<dim>* eptr( static_cast<SIMPLEX<dim>* >(this) );

    assert( eptr->FE() != NULL );
    assert( eptr->NeighborElementVector().size() == eptr->Neighbors() );
    assert( i < eptr->Neighbors() );

    eptr->NeighborElementVector()[ i ] = NULL;
 }



/// Change from clockwise to counter clockwise and vice versa
template<size_t dim, template<size_t> class SIMPLEX>
void  ElementRemeshingTraits<dim,SIMPLEX>::LowDimRevertNodeNumbering()
  {
    SIMPLEX<dim>* eptr( static_cast<SIMPLEX<dim>* >(this) );

    Node<dim>* node_buf(NULL);
    const size_t nodes( eptr->Nodes() );
    const size_t upTo( nodes/2 ); // integer division to index entry before median
    for( size_t i(0); i < upTo; ++i )
    {
        node_buf = eptr->N(i);
        eptr->Assign( i, static_cast<csmp::Node<dim>*>( eptr->N( nodes-(i+1) ) ) );
        eptr->Assign( nodes-(i+1), static_cast<csmp::Node<dim>*>(node_buf) );
    }
    eptr->FE()->CurrentID( FiniteElement::InitialID() );
  }

template<size_t dim, template<size_t> class SIMPLEX>
void  ElementRemeshingTraits<dim,SIMPLEX>::LowDimRevertNeighborNumbering()
  {
    SIMPLEX<dim>* eptr( static_cast<SIMPLEX<dim>* >(this) );

    SIMPLEX<dim>* elmnt_buf(NULL);
    const size_t neighbors( eptr->Neighbors() );
    const size_t upTo( neighbors/2 ); // integer division to index entry before median
    for( size_t i(0); i < upTo; ++i )
    {
        elmnt_buf = eptr->Neighbor(i);
        eptr->Assign( i, static_cast<csmp::Element<dim>*>( eptr->Neighbor( neighbors-(i+1) ) ) );
        eptr->Assign( neighbors-(i+1), static_cast<csmp::Element<dim>*>(elmnt_buf) );
    }
    eptr->FE()->CurrentID( FiniteElement::InitialID() );
  }

template class ElementRemeshingTraits<1U,Element>;
template class ElementRemeshingTraits<2U,Element>;
template class ElementRemeshingTraits<3U,Element>;

  } // end csmp
