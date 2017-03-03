#ifndef CSMP_MESH_DIAGNOSTICS_H
#define CSMP_MESH_DIAGNOSTICS_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t dim> class Model;

/**
    @brief Common checks performed on a newly constructed mesh to detect issues.
    
    MeshDiagnostics determines node-numbering and other issues such as duplicate
    elements in a mesh.
 
    Reports element volume ranges.
    
    @author Stephan Matthai
    @date 1999
*/
template<size_t dim>
class MeshDiagnostics {
  public:
    void ElementVolumeRange( const Model<dim>&, double64& vmin, double64& vmax ) const;
    
    /// checks/reports negative element orientations, flat elements, zero volumes etc.
    bool ScrutinizeMesh( Model<3U>& ) const;

  private:
    /// recomputes a consistent node-numbering for all surfaces in the mesh
    void FixFiniteElementNeighborOrientationOfSurfaceMeshes( Model<dim>& ) const;
};


/// checks input mesh for duplicate elements
template<size_t dim>
bool detectDuplicateElement( const Model<dim>& );



} // csmp

#endif
