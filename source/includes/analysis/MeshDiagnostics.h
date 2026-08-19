#ifndef CSMP_MESH_DIAGNOSTICS_H
#define CSMP_MESH_DIAGNOSTICS_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t dim> class Model;


/**
    @brief Structured result of a mesh scrutiny operation.

    Each counter reflects a distinct category of mesh defect.
    The overall @c problems flag is true if any defect was found.
*/
struct MeshDiagnosticsResult {
    bool     problems              = false; ///< true if any defect was detected
    size_t   coincident_node_pairs = 0U;    ///< node pairs within the same element that coincide
    size_t   negative_volumes      = 0U;    ///< elements with negative volume/area/length
    size_t   zero_volumes          = 0U;    ///< elements with zero or near-zero volume/area/length
    size_t   erratic_volumes       = 0U;    ///< elements with implausibly large volume/area/length
    size_t   erratic_coordinates   = 0U;    ///< nodes with coordinates outside the validity bounds
    size_t   large_aspect_ratios   = 0U;    ///< elements whose segment length ratio exceeds threshold
    double   sl_min                = 0.;    ///< shortest segment length found
    double   sl_max                = 0.;    ///< longest segment length found
    double   aspect_ratio_min      = 0.;    ///< minimum per-element segment length ratio
    double   aspect_ratio_max      = 0.;    ///< maximum per-element segment length ratio
};


/**
    @brief Mesh quality diagnostics for CSMP models of any spatial dimension.

    Checks node coordinates, segment lengths, element aspect ratios,
    and element volumes/areas/lengths. Problematic elements are collected
    into named regions and written to VTK for visualisation.
*/

/**
    @brief Common checks performed on a newly constructed mesh to detect issues.
    
    MeshDiagnostics determines node-numbering and other issues such as duplicate
    elements in a mesh.
 
    Reports element volume ranges.
    
    @author Stephan Matthai
    @date 1999
*/
template<uint32_t dim>
class MeshDiagnostics {
  public:
    /**
        @brief Scrutinises the mesh of the supplied model.

        Validity bounds for coordinates and volumes are derived automatically
        from the model bounding box. The aspect ratio warning threshold and
        the coincidence tolerance can be overridden.

        @param model              The model to scrutinise.
        @param aspect_ratio_warn  Per-element segment length ratio above which
                                  a warning is issued (default 10).
        @param coincidence_tol    Distance below which two nodes in the same
                                  element are considered coincident (default 1e-10).
        @return                   A @c MeshDiagnosticsResult summarising all findings.
    */
    MeshDiagnosticsResult ScrutinizeMesh( Model<dim>& model,
                                          double aspect_ratio_warn = 10.,
                                          double coincidence_tol   = 1.0e-10 ) const;
    
    /// reports volume (3D), area (2D) or length (1D) of the elements in the mesh
    void ElementVolumeRange( const Model<dim>&, double& vmin, double& vmax ) const;
  
    /// creates a range of quality variables, computes them and outputs them from the target region to VTU; creates diagnostic regions for this purpose
    bool ComputeQualityMetricsAndOutputToVTU( Model<dim>&, const std::string& region );

    /// tests for negative element volumes
    bool DetectPotentiallyMisnumberedElements( const Model<dim>& ) const;
  
    /// detects whether different Dirichlet conditions have been assigned to a single finite element
    bool DetectConflictingDirichletConditions( const Model<dim>&, const char* variable_of_interest, VARIABLE_FLAG status=DIRICH ) const;
  
    /// elements where all nodes have a status constraint so that they do not participate in the computation
    bool DetectOverConstrainedElements( const Model<dim>& model, const char* variable_of_interest, VARIABLE_FLAG status=DIRICH ) const;
  
    /// checks whether the any value of a computed node variable (P,T,C) lies outside of the range of the values in its neighborhood
    bool DetectNonMonotonicity( const Model<dim>&, const char* variable_of_interest ) const;

    /// recomputes a consistent node-numbering for all surfaces in the mesh
    void FixFiniteElementNeighborOrientationOfSurfaceMeshes( Model<dim>& ) const;
};


/// checks input mesh for duplicate elements
/// TODO: generalise to cells
template<uint32_t dim>
bool detectCollocatedElements( const Model<dim>& );



} // csmp

#endif
