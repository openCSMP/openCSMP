#ifndef CSMP_FINITE_VOLUME_STENCIL_REFACTORED_H
#define CSMP_FINITE_VOLUME_STENCIL_REFACTORED_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"
#include "DynamicArray2D.h"
#include "DynamicArray3D.h"
#include "Point.h"

// VERSION WITH SINGLE INTEGRATION POINTS PER FACET AND SECTOR
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
#undef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
#endif

namespace csmp {

template<uint32_t> class FiniteVolumeStencilManager;
template<uint32_t> class FV_IntegrationPointsAndWeights;
template<uint32_t> struct TempVecs;

/**

@brief Finite-element specific finite volume stencil (partition into facets and sectors)
of the finite element for which the stencil was initialized.

@author S.K. Matthai
@author A.Paluszny
@date 2003

@section motivation Motivation
 
Class provides a whole set of functionality, related to construction of the
node-centered FV virtual cells for realisation of the CVPEM in CSP.
 
Tesarius: sector - sector part of the finite element (FE), attributed to the FE node.
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
template<uint32_t dim>
class FiniteVolumeStencil {
  public:
    explicit FiniteVolumeStencil( const char* csp_finite_element_type );
    // TODO: refactor and test: rule of zero should work
    FiniteVolumeStencil( const FiniteVolumeStencil& );
    FiniteVolumeStencil& operator=( const FiniteVolumeStencil& );
    ~FiniteVolumeStencil() = default;

    /// number of facets that delimited the FV sector on the insider of the parent finite element
    uint32_t  FacetsPerSector( uint32_t iSector ) const;
  
    /// n-th facet that delimits the FV sector in the inside of the finite element
    uint32_t  FacetSurroundingSector( uint32_t iSector, uint32_t n ) const;

    /// reports the type of a facet
    FV_FACET_TYPE FacetType( uint32_t iFacet ) const;
  
    /// choice of rst-integration point coordinate of the facet integration point
    double  FacetIntegrationPoint( uint32_t iFacet, uint32_t ip, uint32_t r_or_s_or_t ) const;
  
    /// rst coordinate of facet integration point in parametric space
    void    FacetIntegrationPoint( uint32_t iFacet, uint32_t ip, std::vector<double>& rst ) const;
  
    /// returns rst coordinate of facet integration point in parametric space
    const Point<dim>  FacetIntegrationPoint( uint32_t iFacet, uint32_t ip ) const;
  
    /// integration weight of given integration point correlated with facet area and number of facet integration points
    double  FacetIntegrationWeight( uint32_t iFacet, uint32_t ip ) const;
  
    /// scale factor to match sector integration weights between adjacent finite elements that have a different volume in rst space
    double  FacetProjectionWeight( uint32_t iFacet, uint32_t ip ) const;
  
    /// nodes that sit on the opposite sides of the facet; outside is direction into which facet normal points
    void    FacetEdgeNodes( uint32_t iFacet, uint32_t& inside_node, uint32_t& outside_node ) const;
  
    /// opposite side of the outward pointing normal of the facet
    uint32_t InsideNode( uint32_t iFacet ) const;

    /// same side as the outward pointing normal of the facet
    uint32_t OutsideNode( uint32_t iFacet ) const;
  
    /// point on the finite-element edge that is touched by the facet
	  const Point<dim>  FacetEdgeMidPoint( uint32_t iFacet ) const;
  
    /// of finite element in parametric space
    const Point<dim>  Barycenter() const;
  
    double             SectorIntegrationPoint( uint32_t iSector, uint32_t ip, uint32_t r_or_s_or_t ) const;
    void               SectorIntegrationPoint( uint32_t iSector, uint32_t ip, std::vector<double>& rst ) const;
    const Point<dim>   SectorIntegrationPoint( uint32_t iSector, uint32_t ip ) const;
    double             SectorIntegrationWeight( uint32_t iSector, uint32_t ip ) const;
    const Point<dim>   FacetPoint( uint32_t iFacet, uint32_t iPoint ) const;
    uint32_t           FacetPoints( uint32_t iFacet ) const { return static_cast<uint32_t>(facet_points_[iFacet].size()); }

    /// corner points of the (hexahedral) sector in parametric space
    const Point<dim>   SectorPoint( uint32_t iSector, uint32_t iPoint ) const;
  
    /// pairs of sector edge points definining the intersection lines of the sector with the finite element faces
    std::pair<Point<dim>,Point<dim> >  SectorEdgePoints( uint32_t iSector, uint32_t iEdge ) const;
  
    /// pairs of the ids (0..n-1) of the sector edge point definining the intersection lines of the sector with the finite element faces
    const std::pair<uint32_t,uint32_t>  SectorEdge( uint32_t iSector, uint32_t iEdge ) const;
  
    /// number of corner points of the finite volume sector
    uint32_t    SectorPoints( uint32_t iSector ) const;
  
    /// number of edges of the finite volume sector
    uint32_t    SectorEdges( uint32_t iSector ) const;
  
    /// number of finite-volume facets in the parent finite element
    uint32_t    Facets() const;
  
    /// number of finite-volume sectors in the parent finite element
    uint32_t    Sectors() const;
  
    /// number of quadrature points of the finite volume facet i
    uint32_t    IntegrationPointsPerFacet( uint32_t iFacet=0U ) const;
  
    /// number of quadrature points of the finite volume sector i
    uint32_t    IntegrationPointsPerSector( uint32_t iSector=0U ) const;
  
    /// returns unit normal to facet in parametric space
    const Point<dim>  UnitParametricNormalTo( uint32_t iFacet ) const;
  
    /// returns x,y or z component of unit normal in parametric space
    double  UnitParametricNormalComponent( uint32_t iFacet, uint32_t x_or_y_or_z ) const;

    std::pair<double,double>  FacetNormalTransformationNodeWeights( uint32_t iFacet, uint32_t iNode ) const;
  
    /// returns the parent element of the finite volume stencil
    const std::string  ParentElement() const { return parent_element_; }

    /// reports whether the parent element of the stencil is a line, surface or volume
    CELL_SHAPE  Geometry() const;

    /// returns the data (private members) stored in this finite volume stencil
    void  Out() const;
    
  private:
    void  Resize( uint32_t isrf, uint32_t srfs_per_node, uint32_t ivol, uint32_t spts, uint32_t vpts );
    void  ResizeTemporaryVectors( TempVecs<dim>& tvecs,
                                  uint32_t n_isrf,
                                  uint32_t srfs_per_node,
                                  uint32_t n_ivol,
                                  uint32_t n_spts,
                                  uint32_t n_vpts );

    void  Initialize( const char* csp_finite_element_type );
    void  InitializeDataMembers( TempVecs<dim>&, const FV_IntegrationPointsAndWeights<dim>& );

  private:
    FiniteVolumeStencil() = default; // only for FiniteVolumeStencilManager
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    DynamicArray3D<double>    facet_integration_points_;     ///< [isrf][spts][rst]
    DynamicArray3D<double>    sector_integration_points_;    ///< [ivol][vpts][rst]
    DynamicArray2D<double>    facet_integration_weights_;    ///< [isrf][spts]
    DynamicArray2D<double>    facet_projection_weights_;     ///< [isrf][spts]
    DynamicArray2D<double>    sector_integration_weights_;   ///< [ivol][vpts]
#else
    DynamicArray2D<double>    facet_integration_points_;     ///< [isrf][rst]
    DynamicArray2D<double>    sector_integration_points_;    ///< [ivol][rst]
    std::vector<double>       facet_integration_weights_;    ///< [isrf]
    std::vector<double>       facet_projection_weights_;     ///< [isrf]
    std::vector<double>       sector_integration_weights_;   ///< [ivol]
#endif
    DynamicArray2D<double>    facet_edge_midpoints_;         ///< [isrf][rst]
    DynamicArray2D<double>    facet_normals_;                ///< [isrf][dim]
    DynamicArray2D<double>    facet_parametric_normals_;     ///< [isrf][dim]

    DynamicArray2D<uint32_t>  edges_of_element_;             ///< [isrf][2=pair]  element segments = number of facets, each identified by the nodes it connects

    // jagged arrays (where rows have different number of columns) TODO: not refactored yet
    
// DynamicArray3D<double>    facet_normal_xforms_; ///< jagged array for PYRAMID(triangular facets) !- else [isrf][node][2=pair]
    std::vector<std::vector<std::pair<double,double>>> facet_normal_xforms_;   ///< [isrf][node]
// DynamicArray3D<double>    facet_points_;        ///< jagged array for PYRAMID(triangular facets)! - else [isrf][spts][rst]
    std::vector<std::vector<Point<dim> > >      facet_points_;
// DynamicArray3D<double>    sector_points_;       ///< jagged arrat for PYRAMID(triangular facets) !- else [ivol][spts][rst]
    std::vector<std::vector<Point<dim> > >      sector_points_;
//DynamicArray2D<uint32_t>  facets_surrounding_node_; ///< jagged array [ivol][isrf]
    std::vector<std::vector<uint32_t> >  facets_surrounding_node_;      ///< jagged array[node][facet - number is variable]

    ///< for each sector, for each sector edge, pair of local indices of the end nodes of the edge
//DynamicArray3D<uint32_t>  edges_of_sector_; ///< [sector][edges][node1,node2] @attention does not work because it is a JAGGED ARRAY !!!
    std::vector<std::vector<std::pair<uint32_t,uint32_t> > > edges_of_sectors_; ///< edges of FV sector on outside of the finite element

    ///< for each sector, all the points delimiting the sector in parametric space
    std::vector<FV_FACET_TYPE>  facet_types_;      ///< [isrf]

    Point<dim>					      barycenter_;
    std::string               parent_element_;     ///< name of parent finite element
    CSMP_FEM_TYPE             parent_element_type_;
    CELL_SHAPE                space_dimension_;    ///< line, surface, or volumetric parent element
    
    friend class FiniteVolumeStencilManager<dim>;
};






// INLINE FUNCTIONS

/**
    Irrespective of the dimension of the model,
    this method reports whether the stencil corresponds to a line, surface or volumetric
    finite element.
*/
template<uint32_t dim>
inline CELL_SHAPE  FiniteVolumeStencil<dim>::Geometry() const
 {
    return space_dimension_;
 }



/**

Returns the number of internal facets - division walls inside the FE.

@return Returns the number of the internal facets.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
*/
template<uint32_t dim>
inline uint32_t  FiniteVolumeStencil<dim>::Facets() const
 {
    return static_cast<uint32_t>(edges_of_element_.rows());
 }
 
 
 
/**

Returns number of facets, opposite to the given node - this corresponds to
the number of facets, adding flux to the given node.

@return the number of the facets in front of the node of given element type.

@section implementation Implementation

Accesses the private data of the class, retrieving the value.
*/
template<uint32_t dim>
inline uint32_t  FiniteVolumeStencil<dim>::FacetsPerSector( uint32_t iSector ) const
 {
    assert( iSector < facets_surrounding_node_.size() );
    return static_cast<uint32_t>(facets_surrounding_node_[iSector].size());
/*
    if constexpr ( dim == 1u ) return 1U;
    if constexpr ( dim == 2U ) {
         if ( space_dimension_ == LINE ) return 1U;
         return 2U;
      }
    
    // special cases of pyramid apex and lower dimensional elements
    if constexpr( dim == 3U ) {
         if ( space_dimension_ == VOLUME ) {
              if ( parent_element_type_ == ISOPARAMETRIC_LINEAR_PYRAMID && iSector == 4U ) {
 #ifdef PYRAMID_TRIANGULAR_FACETS
                   return 9U; // apex node 4
 #else
                   return 4U; // Special case of the Pyramid vertex 4 is surrounded by 4 facets
 #endif
                 }
              return dim;
           }
         if ( space_dimension_ == SURFACE ) return 2U;
         if ( space_dimension_ == LINE ) return 1U;
      }
    return dim;
*/
 }



/**

Returns number of internal FV sectors, given FE is divided to.

@return the number of sectors.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
*/
template<uint32_t dim>
inline uint32_t  FiniteVolumeStencil<dim>::Sectors() const
 {
    // since this is a vector of vectors of points per sector
//    return static_cast<uint32_t>(sector_points_.depth());
    return static_cast<uint32_t>(sector_points_.size());
 }



/**

Returns number of integration points per facet.

@return the number of integration points per facet.

@section implementation Implementation 

Accesses the private data of the class, retreiving the value.
For current implementation this value is 1.

@attention assumes that the number of integration points per facet is the same for all facets

*/
template<uint32_t dim>
inline uint32_t  FiniteVolumeStencil<dim>::IntegrationPointsPerFacet( uint32_t iFacet ) const
 {
    assert( iFacet < facet_integration_weights_.size() );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return static_cast<uint32_t>(facet_integration_weights_.cols());
#else
    return 1U;
#endif
 }



/**

Returns number of integration points per sectortric sector of given FE
type.

@return the number of integration points per volumetric sector.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
For current implementation this value is 1.
*/
template<uint32_t dim>
inline uint32_t  FiniteVolumeStencil<dim>::IntegrationPointsPerSector( uint32_t iSector ) const
 {
    assert( iSector < sector_integration_weights_.size() );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return static_cast<uint32_t>(sector_integration_weights_.cols());
#else
    return 1U;
#endif
  }




/**

Returns to argument rst the reference to STL vector of the parametric coordinates
of the given facet integration point.

@param iFacet index (No) of the internal facet isnside the FE
@param ip index (No) of the integration point on the facet

@return reference to STL vector of the parametric coordinates 
of the given facet integration point

@section implementation Implementation

Accesses the private data of the tabulated points, associated with the facet.

@section application  Application

For current implementation index ip is constrained to 0 only, i.e. the
FVPEM method is working with 1 facet integration point only.
*/
template<uint32_t dim>
inline const Point<dim>  FiniteVolumeStencil<dim>::FacetIntegrationPoint( uint32_t iFacet,
                                                                          uint32_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    if constexpr ( dim == 1u )
      return Point<dim>( facet_integration_points_(iFacet,ip,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_integration_points_(iFacet,ip,0u), facet_integration_points_(iFacet,ip,1u) );
    else if constexpr ( dim == 3u ) {
       facet_integration_points_.out();
       return Point<dim>( facet_integration_points_(iFacet,ip,0u),
                          facet_integration_points_(iFacet,ip,1u),
                          facet_integration_points_(iFacet,ip,2u) );
      }
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    if constexpr ( dim == 1u )
      return Point<dim>( facet_integration_points_(iFacet,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_integration_points_(iFacet,0u), facet_integration_points_(iFacet,1u) );
    else if constexpr ( dim == 3u ) {
       facet_integration_points_.out();
       return Point<dim>( facet_integration_points_(iFacet,0u),
                          facet_integration_points_(iFacet,1u),
                          facet_integration_points_(iFacet,2u) );
      }
#endif
}



template<uint32_t dim>
inline void FiniteVolumeStencil<dim>::FacetIntegrationPoint( uint32_t iFacet,
                                                             uint32_t ip,
                                                             std::vector<double>& rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    rst.resize(dim);

#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    if constexpr ( dim == 1u )
      rst[0u] = facet_integration_points_(iFacet,ip,0u);
    else if constexpr ( dim == 2u ) {
        rst[0u] = facet_integration_points_(iFacet,ip,0u);
        rst[1u] = facet_integration_points_(iFacet,ip,1u);
      }
    else if constexpr ( dim == 3u ) {
       facet_integration_points_.out();
       rst[0u] = facet_integration_points_(iFacet,ip,0u);
       rst[1u] = facet_integration_points_(iFacet,ip,1u);
       rst[2u] = facet_integration_points_(iFacet,ip,2u);
     }
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    if constexpr ( dim == 1u )
      rst[0u] = facet_integration_points_(iFacet,0u);
    else if constexpr ( dim == 2u ) {
        rst[0u] = facet_integration_points_(iFacet,0u);
        rst[1u] = facet_integration_points_(iFacet,1u);
      }
    else if constexpr ( dim == 3u ) {
       facet_integration_points_.out();
       rst[0u] = facet_integration_points_(iFacet,0u);
       rst[1u] = facet_integration_points_(iFacet,1u);
       rst[2u] = facet_integration_points_(iFacet,2u);
     }
#endif
 }



// SKM addon
template<uint32_t dim>
inline double FiniteVolumeStencil<dim>::FacetIntegrationPoint( uint32_t iFacet,
                                                               uint32_t ip,
                                                               uint32_t rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return facet_integration_points_(iFacet,ip,rst);
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    return facet_integration_points_(iFacet,rst);
#endif
 }



/**

Defines FE edge, associated with that facet of interest. For a given
facet method returns local indices of the FE start node and 
the end node. The end node is the one on the side of the facet into
which the outward pointing normal points.  

@param iFacet index of the internal dividing facet in the FE.
@param estart local index of the FE edge start  node
@param eend   local index of the FE edge end node

@section implementation Implementation

Accesses the private pair of nodes, associated with the facet.
*/
template<uint32_t dim>
inline void FiniteVolumeStencil<dim>::FacetEdgeNodes( uint32_t iFacet,
                                                      uint32_t& estart,
                                                      uint32_t& eend ) const
{
    assert( iFacet < Facets() );
    estart=edges_of_element_(iFacet,0U);
    eend=edges_of_element_(iFacet,1U);
}


template<uint32_t dim>
inline uint32_t FiniteVolumeStencil<dim>::InsideNode( uint32_t iFacet ) const
{
    assert( iFacet < Facets() );
    return edges_of_element_(iFacet,0u);
}


template<uint32_t dim>
inline uint32_t FiniteVolumeStencil<dim>::OutsideNode( uint32_t iFacet ) const
{
    assert( iFacet < Facets() );
    return edges_of_element_(iFacet,1U);
}



/**

Returns to argument rst the reference to STL vector of the parametric coordinates
of the given sector integration point, associated with given sector "volume".

@param iSector number of the internal sector isnside the FE
@param ip index (No) of the integration point on the facet

@return reference to STL vector of the parametric coordinates
of the given facet integration point

@section implementation Implementation

Accesses the private data of the tabulated points, associated with the sector.

@section application Application

For current implementation index ip is constrained to 0 only, i.e. the
FVPEM method is working with 1 facet integration point only.
*/
template<uint32_t dim>
inline const Point<dim>  FiniteVolumeStencil<dim>::SectorIntegrationPoint( uint32_t iSector,
                                                                           uint32_t ip ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    if constexpr ( dim == 1u )
      return Point<dim>( sector_integration_points_(iSector,ip,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( sector_integration_points_(iSector,ip,0u), sector_integration_points_(iSector,ip,1u) );
    else if constexpr ( dim == 3u )
       return Point<dim>( sector_integration_points_(iSector,ip,0u),
                          sector_integration_points_(iSector,ip,1u),
                          sector_integration_points_(iSector,ip,2u) );
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    if constexpr ( dim == 1u )
      return Point<dim>( sector_integration_points_(iSector,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( sector_integration_points_(iSector,0u), sector_integration_points_(iSector,1u) );
    else if constexpr ( dim == 3u )
       return Point<dim>( sector_integration_points_(iSector,0u),
                          sector_integration_points_(iSector,1u),
                          sector_integration_points_(iSector,2u) );
#endif
 }



template<uint32_t dim>
inline void FiniteVolumeStencil<dim>::SectorIntegrationPoint( uint32_t iSector,
                                                              uint32_t ip,
                                                              std::vector<double>& rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
    rst.resize(dim);
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    if constexpr ( dim == 1u )
      rst[0u] = sector_integration_points_(iSector,ip,0u);
    else if constexpr ( dim == 2u ) {
        rst[0u] = sector_integration_points_(iSector,ip,0u);
        rst[1u] = sector_integration_points_(iSector,ip,1u);
      }
    else if constexpr ( dim == 3u ) {
       rst[0u] = sector_integration_points_(iSector,ip,0u);
       rst[1u] = sector_integration_points_(iSector,ip,1u);
       rst[2u] = sector_integration_points_(iSector,ip,2u);
     }
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    if constexpr ( dim == 1u )
      rst[0u] = sector_integration_points_(iSector,0u);
    else if constexpr ( dim == 2u ) {
        rst[0u] = sector_integration_points_(iSector,0u);
        rst[1u] = sector_integration_points_(iSector,1u);
      }
    else if constexpr ( dim == 3u ) {
       rst[0u] = sector_integration_points_(iSector,0u);
       rst[1u] = sector_integration_points_(iSector,1u);
       rst[2u] = sector_integration_points_(iSector,2u);
     }
#endif
 }



template<uint32_t dim>
inline double FiniteVolumeStencil<dim>::SectorIntegrationPoint( uint32_t iSector,
                                                                uint32_t ip,
                                                                uint32_t rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return sector_integration_points_(iSector,ip,rst);
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    return sector_integration_points_(iSector,rst);
#endif
 }


/**

Returns for the given facet integration point integration weight,
associated with given facet "facet".

@param iFacet index (No) of the internal facet  isnside the FE
@param ip index (No) of the integration point on the facet

@return the value of the facet integration weight.

@section implementation Implementation

Accesses the private data of the tabulated points, associated with the facet.

@section application Application

For current implementation index ip is constrained to 0 only, i.e. the
FVPEM method is working with 1 facet integration point only.
*/    
template<uint32_t dim>
inline double FiniteVolumeStencil<dim>::FacetIntegrationWeight( uint32_t iFacet,
                                                                uint32_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return facet_integration_weights_(iFacet,ip);
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    return facet_integration_weights_[iFacet];
#endif
 }
 


/**

Returns for the given sector integration point projection weight,
associated with given facet "facet". Projection weight is an extra multiplier,
which is bringing to the same scale the size of the unit parametric
spaces of different element types.

@param iFacet index (No) of the internal facet  isnside the FE
@param ip index (No) of the integration point on the facet

@return the value of the facet projection weight.

@section implementation Implementation

Accesses the private data of the tabulated points, associated with the facet.

@section application Application

For current implementation index ip is constrained to 0 only, i.e. the
FVPEM method is working with 1 facet integration point only.
*/
 
template<uint32_t dim>
inline double FiniteVolumeStencil<dim>::FacetProjectionWeight( uint32_t iFacet, uint32_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return facet_projection_weights_(iFacet,ip);
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    return facet_projection_weights_[iFacet];
#endif
 }



/**

Returns for the given sector integration point integration weight,
associated with given sector "volume".

@param iSector corresponding to the finite element node with the same number
@param ip represents the n-th integration point on the facet

@return the value of the sector integration weight=volume in parametric space.

@section implementation Implementation

Accesses the private data of the tabulated points, associated with the sector.

@section application Application

For current implementation index ip is constrained to 0 only, i.e. the
FVPEM method is working with 1 facet integration point only.
*/
template<uint32_t dim>
inline double FiniteVolumeStencil<dim>::SectorIntegrationWeight( uint32_t iSector, uint32_t ip ) const
 {
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    return sector_integration_weights_(iSector,ip);
#else
    assert( ip == 0U ); // restrictive assumpting of this scheme for testing
    return sector_integration_weights_[iSector];
#endif
 }
 
 
 
  // overloaded method, returns specific index of a facet 
template<uint32_t dim>
inline uint32_t  FiniteVolumeStencil<dim>::FacetSurroundingSector( uint32_t iSector, uint32_t n ) const
 {
    assert( iSector < Sectors() );
    assert( n < FacetsPerSector(iSector) );
//    return facets_surrounding_node_(iSector,n);
    return facets_surrounding_node_[iSector][n];
 }

 


/**

Returns the normal vector in parametric space of the FE.
 
@param iFacet is the number 0..n-1 of the internal facet inside of the FE

@return to argument "nrml" reference to STL vector of the normal vector
 coordinates in parametric space of the FE. 

@section implementation Implementation
 
Accesses the private data of the class, retreiving tabulated vector.

@attention here parametric space has same dimensions as physical space, but third dimension may be zero
*/
template<uint32_t dim>
inline const Point<dim>  FiniteVolumeStencil<dim>::UnitParametricNormalTo( uint32_t iFacet ) const
{
    assert( iFacet < facet_types_.size() );
    if constexpr ( dim == 1u )
      return Point<dim>( facet_parametric_normals_(iFacet,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_parametric_normals_(iFacet,0u), facet_parametric_normals_(iFacet,1u) );
    else if constexpr ( dim == 3u )
      return Point<dim>( facet_parametric_normals_(iFacet,0u),
                         facet_parametric_normals_(iFacet,1u),
                         facet_parametric_normals_(iFacet,2u) );
}



template<uint32_t dim>
inline double FiniteVolumeStencil<dim>::UnitParametricNormalComponent( uint32_t iFacet,
                                                                uint32_t x_or_y_or_z ) const
{
   assert( iFacet < facet_types_.size() );
   assert( x_or_y_or_z < dim );
   return facet_parametric_normals_(iFacet,x_or_y_or_z);
}




/**

Returns the component of the transformation from nodes to normals in physical space.
 
@param iFacet is the number 0..n-1 of the internal facet inside of the FE
@param iNode is the number 0..n-1 of the FE node

@return a pair of weights, one for each tangent vector

@section implementation Implementation

Accesses the private data of the class, retreiving tabulated vector.

*/
template<uint32_t dim>
inline std::pair<double,double>  FiniteVolumeStencil<dim>::FacetNormalTransformationNodeWeights(
                                                               uint32_t iFacet, uint32_t iNode ) const
{
//    assert( iFacet < facet_normal_xforms_.depth() );
//    assert( iNode < sector_integration_weights_.size() );
//    return make_pair( facet_normal_xforms_(iFacet,iNode,0u), facet_normal_xforms_(iFacet,iNode,1u) );
    assert( iFacet < facet_normal_xforms_.size() );
    assert( iNode < facet_normal_xforms_[iFacet].size() );
    return facet_normal_xforms_[iFacet][iNode];
}



template<uint32_t dim>
inline const Point<dim>  FiniteVolumeStencil<dim>::FacetEdgeMidPoint( uint32_t iFacet ) const
 {
    assert( iFacet < facet_types_.size() );
    if constexpr ( dim == 1u )
      return Point<dim>( facet_edge_midpoints_(iFacet,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_edge_midpoints_(iFacet,0u), facet_edge_midpoints_(iFacet,1u) );
    else if constexpr ( dim == 3u )
      return Point<dim>( facet_edge_midpoints_(iFacet,0u),
                         facet_edge_midpoints_(iFacet,1u),
                         facet_edge_midpoints_(iFacet,2u) );
}


template<uint32_t dim>
inline const Point<dim>  FiniteVolumeStencil<dim>::Barycenter() const
{
	 return barycenter_;
}


template<uint32_t dim>
inline const Point<dim>  FiniteVolumeStencil<dim>::FacetPoint( uint32_t iFacet, uint32_t iPoint ) const
 {
    assert( iFacet < facet_points_.size() );
    assert( iPoint < facet_points_[iPoint].size() );
    return facet_points_[iFacet][iPoint];
    /*
    assert( iFacet < facet_points_.depth() );
    assert( iPoint < facet_points_.rows() );
    if constexpr ( dim == 1u )
      return Point<dim>( facet_points_(iFacet,iPoint,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_points_(iFacet,iPoint,0u), facet_points_(iFacet,iPoint,1u) );
    else if constexpr ( dim == 3u )
      return Point<dim>( facet_points_(iFacet,iPoint,0u),
                         facet_points_(iFacet,iPoint,1u),
                         facet_points_(iFacet,iPoint,2u) );
    */
 }



template<uint32_t dim>
inline FV_FACET_TYPE  FiniteVolumeStencil<dim>::FacetType( uint32_t iFacet ) const
{
   assert( iFacet < facet_types_.size() );
	 return facet_types_[iFacet];
}

    
    
template<uint32_t dim>
inline const csmp::Point<dim>  FiniteVolumeStencil<dim>::SectorPoint( uint32_t iSector, uint32_t iPoint ) const
 {
    assert( iSector < sector_points_.size() );
    assert( iPoint < sector_points_[iSector].size() );
    return sector_points_[iSector][iPoint];
    /*
    assert( iSector < sector_points_.depth() );
    assert( iPoint < sector_points_.rows() );
    if constexpr ( dim == 1u )
      return Point<dim>( sector_points_(iSector,iPoint,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( sector_points_(iSector,iPoint,0u), sector_points_(iSector,iPoint,1u) );
    else if constexpr ( dim == 3u )
      return Point<dim>( sector_points_(iSector,iPoint,0u),
                         sector_points_(iSector,iPoint,1u),
                         sector_points_(iSector,iPoint,2u) );
    */
 }
  
  
  
template<uint32_t dim>
inline std::pair<csmp::Point<dim>,csmp::Point<dim> >  FiniteVolumeStencil<dim>::SectorEdgePoints(
                                                               uint32_t iSector, uint32_t iEdge ) const
{
   assert( iSector < edges_of_sectors_.size() );
   assert( iEdge < edges_of_sectors_[iSector].size() );
   const uint32_t pt1(edges_of_sectors_[iSector][iEdge].first);
   const uint32_t pt2(edges_of_sectors_[iSector][iEdge].second);
   
	 return std::make_pair( sector_points_[iSector][pt1], sector_points_[iSector][pt2] );
   /*
    assert( iSector < edges_of_sectors_.size() );
    assert( iEdge < edges_of_sectors_[iSector].size() );
    const uint32_t pt1(edges_of_sectors_[iSector][iEdge].first);
    const uint32_t pt2(edges_of_sectors_[iSector][iEdge].second);
   
    if constexpr ( dim == 1u )
      return make_pair( Point<1U>(sector_points_(iSector,pt1,0u)), Point<1U>(sector_points_(iSector,pt2,0u)) );
    else if constexpr ( dim == 2u )
      return make_pair( Point<2U>(sector_points_(iSector,pt1,0u),sector_points_(iSector,pt1,1u)),
                        Point<2U>(sector_points_(iSector,pt2,0u),sector_points_(iSector,pt2,1u)) );
    else if constexpr ( dim == 3u )
      return make_pair( Point<3U>(sector_points_(iSector,pt1,0u),sector_points_(iSector,pt1,1u), sector_points_(iSector,pt1,2u)),
                        Point<3U>(sector_points_(iSector,pt2,0u),sector_points_(iSector,pt2,1u), sector_points_(iSector,pt2,2u)) );
  */
}



template<uint32_t dim>
inline const std::pair<uint32_t,uint32_t>  FiniteVolumeStencil<dim>::SectorEdge( uint32_t iSector,
                                                                                 uint32_t iEdge ) const
{
   assert( iSector < edges_of_sectors_.size() );
   assert( iEdge < edges_of_sectors_[iSector].size() );
	 return edges_of_sectors_[iSector][iEdge];
}
    
    
/// number of corner points of the finite volume sector
template<uint32_t dim>
inline uint32_t FiniteVolumeStencil<dim>::SectorPoints( uint32_t iSector ) const
{
   assert( iSector < sector_points_.size() );
   // assumes that all sectors have the same number of points
   return static_cast<uint32_t>(sector_points_[iSector].size());
}


template<uint32_t dim>
inline uint32_t FiniteVolumeStencil<dim>::SectorEdges( uint32_t iSector ) const
{
   assert( iSector < edges_of_sectors_.size() );
   return static_cast<uint32_t>(edges_of_sectors_[iSector].size());
}
    
 
} // end namespace csmp

#endif
