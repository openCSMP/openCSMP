#ifndef INTERFACE_REMESHING_TRAITS_H
#define INTERFACE_REMESHING_TRAITS_H

#include "Node.h"

namespace csmp {

/**

    Local operations on mesh
    @author Roman
    @date 8/8/2014
*/
template<size_t dim, template<size_t> class SIMPLEX>
class InterFaceRemeshingTraits
  {
  public:

    InterFaceRemeshingTraits();
    
    /// unassign neighbor interface
    void UnassignNeighbor( size_t nbor );
    void UnassignNeighbors(  );

    /// detaches interface from its neighbors
    void Detach();

    /// exchange parents
    void Flip();

  private:

    InterFaceRemeshingTraits( const SIMPLEX<dim>& );

  };
 
} // end csmp

#endif
