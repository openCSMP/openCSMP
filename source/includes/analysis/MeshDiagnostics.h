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
    void ElementVolumeRange( const Model<dim>&, double& vmin, double& vmax ) const;
  
    /// summary: diagnostics including checks/reports negative element orientations, flat elements, zero volumes etc.
    bool ScrutinizeMesh( Model<3U>& ) const;

    /// tests for negative element volumes
    bool DetectPotentiallyMisnumberedElements( const Model<dim>& ) const;
  
    /// detects whether different Dirichlet conditions have been assigned to a single finite element
    bool DetectConflictingDirichletConditions( const Model<dim>&, const char* variable_of_interest, VARIABLE_FLAG status=DIRICH ) const;
  
    /// elements where all nodes have a status constraint so that they do not participate in the computation
    bool DetectOverConstrainedElements( const Model<dim>& model, const char* variable_of_interest, VARIABLE_FLAG status=DIRICH ) const;
  
    /// checks whether the any value of a computed node variable (P,T,C) lies outside of the range of the values in its neighborhood
    bool DetectNonMonotonicity( const Model<dim>&, const char* variable_of_interest ) const;

  private:
    /// recomputes a consistent node-numbering for all surfaces in the mesh
    void FixFiniteElementNeighborOrientationOfSurfaceMeshes( Model<dim>& ) const;
};


/// checks input mesh for duplicate elements
template<size_t dim>
bool detectDuplicateElements( const Model<dim>& );



} // csmp

#endif
