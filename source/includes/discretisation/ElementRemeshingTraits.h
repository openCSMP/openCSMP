#ifndef ELEMENT_REMESHING_TRAITS_H
#define ELEMENT_REMESHING_TRAITS_H

#include "Node.h"

namespace csmp {

/**
    Local operations on mesh
    @author Roman
    @date 8/8/2014
*/
template<size_t dim, template<size_t> class SIMPLEX>
class ElementRemeshingTraits
  {
  public:

    ElementRemeshingTraits();

    /// detaches element from its neighbors and from its nodes (removes itself as parent)
    void UnassignNodes( );
    void UnassignNeighbors( );
    void UnassignNeighbor( size_t );
    void Detach();

    /// revert node and neighbors in case if low dimensional element oriented not correctly
    void LowDimRevertNodeNumbering( );
    void LowDimRevertNeighborNumbering( );

  private:

    ElementRemeshingTraits( const SIMPLEX<dim>& );

  };
 
} // end csmp

#endif
