#include "InterFaceRemeshingTraits.h"
#include "InterFace.h"

namespace csmp {


template<size_t dim, template<size_t> class SIMPLEX>
InterFaceRemeshingTraits<dim,SIMPLEX>::InterFaceRemeshingTraits()
 {

 }

template<size_t dim, template<size_t> class SIMPLEX>
void  InterFaceRemeshingTraits<dim,SIMPLEX>::UnassignNeighbor( size_t i ) // unassign neighbour interface
 {
    SIMPLEX<dim>* fptr( static_cast<SIMPLEX<dim>* >(this) );

    assert( fptr->FE() != NULL );
    assert( fptr->NeighborElementVector().size() == fptr->Neighbors() );
    assert( i < fptr->Neighbors() );

    fptr->NeighborElementVector()[ i ] = NULL;
 }

template<size_t dim, template<size_t> class SIMPLEX>
void  InterFaceRemeshingTraits<dim,SIMPLEX>::UnassignNeighbors(  ) // unassign neighbour interfaces
 {
    SIMPLEX<dim>* fptr( static_cast<SIMPLEX<dim>* >(this) );

    // remove element from neighbor list of its neighbors
    for( size_t neighbor = 0; neighbor < fptr->Neighbors(); ++neighbor )
      if( fptr->Neighbor( neighbor ) != NULL )
          for( size_t i = 0; i < fptr->Neighbor( neighbor )->Neighbors(); i++ )
              if( fptr->Neighbor(neighbor)->Neighbor( i ) != NULL )
                  if( fptr->Neighbor(neighbor)->Neighbor( i ) == fptr )
                      fptr->Neighbor(neighbor)->UnassignNeighbor( i );
 }

/** Detaches the interface from its neighbors
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  InterFaceRemeshingTraits<dim,SIMPLEX>::Detach()
 {
    SIMPLEX<dim>* fptr( static_cast<SIMPLEX<dim>* >(this) );

    // remove element from neighbor list of its neighbors
    for( size_t neighbor = 0; neighbor < fptr->Neighbors(); ++neighbor )
      if( fptr->Neighbor( neighbor ) != NULL )
          for( size_t i = 0; i < fptr->Neighbor( neighbor )->Neighbors(); i++ )
              if( fptr->Neighbor(neighbor)->Neighbor( i ) != NULL )
                  if( fptr->Neighbor(neighbor)->Neighbor( i ) == fptr )
                      fptr->Neighbor(neighbor)->UnassignNeighbor( i );
 }

/// Exchanges parents (as a consequence orientation in space)
template<size_t dim, template<size_t> class SIMPLEX>
void  InterFaceRemeshingTraits<dim,SIMPLEX>::Flip()
 {
    SIMPLEX<dim>* fptr( static_cast<SIMPLEX<dim>* >(this) );
    fptr->Assign( fptr->InnerParent( ), fptr->OuterParent() );
 }

template class InterFaceRemeshingTraits<1U,InterFace>;
template class InterFaceRemeshingTraits<2U,InterFace>;
template class InterFaceRemeshingTraits<3U,InterFace>;

} // end csmp
