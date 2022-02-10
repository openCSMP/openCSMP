#ifndef CSMP_FINITE_VOLUME_STENCIL_H
#define CSMP_FINITE_VOLUME_STENCIL_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"
#include "Point.h"

namespace csmp {

template<size_t> class FiniteVolumeStencilManager;


/**

@brief Finite-element specific finite volume stencil (partition into facets and sectors)
of the finite element for which the stencil was initialized.

@author S.K. Matthai
@author A.Paluszny
@date 2003

@section motivation Motivation
 
Class provides a whole set of functionality, related to construction of the
node-centered FV virtual cells for realisation of the CVPEM in CSP.
 
Tesarius: sector - sector part of the finite element (FE) , attributed to the FE node.
Number of sectors corresponds to the number of nodes. Sector is a topological
hexahedron and is defined in physical and parametric spaces of the FE.
Physical and parametric spaces of the FE are related by the mapping,
characterized by Jacobian transformation.
 
facet  - four-side face inside the FE, dividing  it to volumetric sectors.
 
Other set of methods here is providing accessors for the stencil dynamic data:
integration points in the parametric space of the finite element, integration
weights, normals and correction multiplies for unification of the parametric
paces of different elements.
 
@section applicability Applicability

Class is applicable to any FE model, containing isoparametric linear
elements. It can work on the mono- and poly-element type meshes, containing
tetrahedrons, hexahedrons, prisms  and pyramids.
 
@section participants Participants

Gets it data from IntegrationPointsAndWeights.
 

 
@section implementation Implementation

Class provides facet / sector integration methods for finite volumes
associated with specific finite elements. For tabulation of the points
and weights if uses class FV_IntegrationPointsAndWeights ().

*/
template<size_t dim>
class FiniteVolumeStencil {
  public:
    explicit FiniteVolumeStencil( const char* csp_finite_element_type );
    FiniteVolumeStencil( const FiniteVolumeStencil& );
    FiniteVolumeStencil& operator=( const FiniteVolumeStencil& );
    ~FiniteVolumeStencil();
    void      Initialize( const char* csp_finite_element_type );
    void      Resize( size_t isrf, size_t srfs_per_node, size_t ivol, size_t spts, size_t vpts );
  
    /// number of facets that delimited the FV sector on the insider of the parent finite element
    size_t    FacetsPerSector( size_t iSector ) const;
  
    /// n-th facet that delimits the FV sector in the inside of the finite element
    size_t    FacetSurroundingSector( size_t iSector, size_t n ) const;

    /// reports the type of a facet
    FV_FACET_TYPE FacetType( size_t iFacet ) const;
  
    /// choice of rst-integration point coordinate of the facet integration point
    double  FacetIntegrationPoint( size_t iFacet, size_t ip, size_t r_or_s_or_t ) const;
  
    /// rst coordinate of facet integration point in parametric space
    void      FacetIntegrationPoint( size_t iFacet, size_t ip, std::vector<double>& rst ) const;
  
    /// returns rst coordinate of facet integration point in parametric space
    const Point<dim>&  FacetIntegrationPoint( size_t iFacet, size_t ip ) const;
  
    /// integration weight of given integration point correlated with facet area and number of facet integration points
    double  FacetIntegrationWeight( size_t iFacet, size_t ip ) const;
  
    /// scale factor to match sector integration weights between adjacent finite elements that have a different volume in rst space
    double  FacetProjectionWeight( size_t iFacet, size_t ip ) const;
  
    /// nodes that sit on the opposite sides of the facet; outside is direction into which facet normal points
    void      FacetEdgeNodes( size_t iFacet, size_t& inside_node, size_t& outside_node ) const;
  
    /// opposite side of the outward pointing normal of the facet
    size_t    InsideNode( size_t iFacet ) const;

    /// same side as the outward pointing normal of the facet
    size_t    OutsideNode( size_t iFacet ) const;
  
    /// point on the finite-element edge that is touched by the facet
	  const Point<dim>& FacetEdgeMidPoint(size_t iFacet ) const;
  
    /// of finite element in parametric space
    const Point<dim>& Barycenter() const;
  
    double  SectorIntegrationPoint( size_t iSector, size_t ip, size_t r_or_s_or_t ) const;
    void      SectorIntegrationPoint( size_t iSector, size_t ip, std::vector<double>& rst ) const;
    const Point<dim>&  SectorIntegrationPoint( size_t iSector, size_t ip ) const;
    double  SectorIntegrationWeight( size_t iSector, size_t ip ) const;
    const Point<dim>&  FacetPoint( size_t iFacet, size_t iPoint ) const;
    size_t    FacetPoints( size_t iFacet ) const { return facet_points[iFacet].size(); }

    /// corner points of the (hexahedral) sector in parametric space
    const Point<dim>&  SectorPoint( size_t iSector, size_t iPoint ) const;
  
    /// pairs of sector edge points definining the intersection lines of the sector with the finite element faces
    std::pair<Point<dim>,Point<dim> >  SectorEdgePoints( size_t iSector, size_t iEdge ) const;
  
    /// pairs of the ids (0..n-1) of the sector edge point definining the intersection lines of the sector with the finite element faces
    const std::pair<size_t,size_t>&  SectorEdge( size_t iSector, size_t iEdge ) const;
  
    /// number of corner points of the finite volume sector
    size_t    SectorPoints( size_t iSector ) const;
  
    /// number of edges of the finite volume sector
    size_t    SectorEdges( size_t iSector ) const; 
  
    /// number of finite-volume facets in the parent finite element
    size_t    Facets() const;
  
    /// number of finite-volume sectors in the parent finite element
    size_t    Sectors() const;
  
    /// number of quadrature points of the finite volume facet i
    size_t    IntegrationPointsPerFacet( size_t iFacet=0U ) const;
  
    /// number of quadrature points of the finite volume sector i
    size_t    IntegrationPointsPerSector( size_t iSector=0U ) const;
  
    /// returns unit normal to facet in parametric space
    const Point<dim>&  UnitParametricNormalTo( size_t iFacet ) const;
  
    /// returns x,y or z component of unit normal in parametric space
    double           UnitParametricNormalComponent( size_t iFacet, size_t x_or_y_or_z ) const;

    std::pair<double,double>    FacetNormalTransformationNodeWeights( size_t iFacet, size_t iNode ) const;
  
    /// returns the parent element of the finite volume stencil
    const std::string& ParentElement() const { return parent_element_; }

    /// reports whether the parent element of the stencil is a line, surface or volume
    CELL_SHAPE  Geometry() const;

    /// returns the data (private members) stored in this finite volume stencil
    void  Out() const;
    
  private:
    FiniteVolumeStencil() { /* do not use this default constructor */ };
    // TODO: flatten all these multidimensional arrays for faster access; use specific type for this
    std::vector<std::pair<size_t,size_t> >  edges_of_element;             ///< = facets
    std::vector<std::vector<size_t> >       facets_surrounding_node;      ///< [node][facet]
    std::vector<std::vector<Point<dim> > >  facet_integration_points;     ///< [isrf][spts][dim]
    std::vector<std::vector<double> >       facet_integration_weights;    ///< [isrf][spts]
    std::vector<std::vector<double> >       facet_projection_weights;     ///< [isrf][spts]
    std::vector<std::vector<double> >       facet_normals;                ///< [isrf][dim] //[node*3][dim] -3d
    std::vector<Point<dim> >                facet_parametric_normals;     ///< [isrf][dim] //[node*3][dim] -3d
    std::vector<std::vector<std::pair<double,double>>> facet_normal_xforms; ///< [isrf][node]
    std::vector<std::vector<Point<dim> > >  sector_integration_points;    ///< [ivol][vpts][dim]
    std::vector<std::vector<double> >       sector_integration_weights;   ///< [ivol][vpts]
    std::vector<Point<dim> >                facet_edge_midpoints;
    Point<dim>					          	        barycenter;
    std::vector<std::vector<Point<dim> > >  facet_points;
    ///< for each sector, all the points delimiting the sector in parametric space
    std::vector<FV_FACET_TYPE>              facet_types;                  ///< [isrf]
    std::vector<std::vector<Point<dim> > >  sector_points_;         
    ///< for each sector, for each sector edge, the indices of the end points
    std::vector<std::vector<std::pair<size_t,size_t> > >  sector_edges_;   
    std::string                             parent_element_;              ///< name of parent finite element
    CELL_SHAPE                              space_dimension_;             ///< line, surface, or volumetric parent element
    
    friend class FiniteVolumeStencilManager<dim>;
};

} // end namespace csmp

#endif
