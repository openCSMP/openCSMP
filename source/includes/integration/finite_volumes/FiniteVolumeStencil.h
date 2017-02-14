#ifndef CSMP_FINITE_VOLUME_STENCIL_H
#define CSMP_FINITE_VOLUME_STENCIL_H

#include "CSMP_definitions.h"
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
    FiniteVolumeStencil( const FiniteVolumeStencil& fvs );
    FiniteVolumeStencil& operator=( const FiniteVolumeStencil& fvs );
    ~FiniteVolumeStencil();
    void      Initialize( const char* csp_finite_element_type );
    void      Resize( size_t isrf, size_t srfs_per_node, size_t ivol, size_t spts, size_t vpts );
  
    /// number of facets that delimited the FV sector on the insider of the parent finite element
    size_t    FacetsPerSector( size_t iSector ) const;
  
    /// n-th facet that delimits the FV sector in the inside of the finite element
    size_t    FacetSurroundingSector( size_t iSector, size_t n ) const;
  
    /// choice of rst-integration point coordinate of the facet integration point
    double64  FacetIntegrationPoint( size_t iFacet, size_t ip, size_t r_or_s_or_t ) const;
  
    /// rst coordinate of facet integration point in parametric space
    void      FacetIntegrationPoint( size_t iFacet, size_t ip, std::vector<double64>& rst ) const;
  
    /// returns rst coordinate of facet integration point in parametric space
    const Point<dim>&  FacetIntegrationPoint( size_t iFacet, size_t ip ) const;
  
    /// integration weight of given integration point correlated with facet area and number of facet integration points
    double64  FacetIntegrationWeight( size_t iFacet, size_t ip ) const;
  
    /// scale factor to match sector integration weights between adjacent finite elements that have a different volume in rst space
    double64  FacetProjectionWeight( size_t iFacet, size_t ip ) const;
  
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
  
    double64  SectorIntegrationPoint( size_t iSector, size_t ip, size_t r_or_s_or_t ) const;
    void      SectorIntegrationPoint( size_t iSector, size_t ip, std::vector<double64>& rst ) const;
    const Point<dim>&  SectorIntegrationPoint( size_t iSector, size_t ip ) const;
    double64  SectorIntegrationWeight( size_t iSector, size_t ip ) const;
    const Point<dim>&  FacetPoint( size_t iFacet, size_t iPoint ) const;
    size_t    FacetPoints( size_t iFacet ) const { return facet_points[iFacet].size(); }

    // sector intersection line with finite element faces
    const Point<dim>&  SectorPoint( size_t iSector, size_t iPoint ) const;              
    std::pair<Point<dim>,Point<dim> > SectorEdgePoints( size_t iSector, size_t iEdge ) const;
    const std::pair<size_t,size_t>&  SectorEdge( size_t iSector, size_t iEdge ) const;
    size_t    SectorPoints( size_t iSector ) const; 
    size_t    SectorEdges( size_t iSector ) const; 
    
    size_t    Facets() const;
    size_t    Sectors() const;
    size_t    IntegrationPointsPerFacet( size_t iFacet=0U ) const;
    size_t    IntegrationPointsPerSector( size_t iSector=0U ) const;
  
    /// returns unit normal to facet in parametric space
    const Point<dim>&  UnitParametricNormalTo( size_t iFacet ) const;
  
    /// returns x,y or z component of unit normal in parametric space
    double64           UnitParametricNormalComponent( size_t iFacet, size_t x_or_y_or_z ) const;
  
    /// returns the parent element of the finite volume stencil
    const std::string& ParentElement() const { return parent_element_; }
  
    /// returns the data (private members) stored in this finite volume stencil
    void Out() const { Out(std::cout); }
    void  Out(std::ostream& os) const;
    
  private:
    FiniteVolumeStencil() { /* do not use this default constructor */ };
    std::vector<std::pair<size_t,size_t> >  edges_of_element;             // = facets
    std::vector<std::vector<size_t> >       facets_surrounding_node;      // [node][facet]
    std::vector<std::vector<Point<dim> > >  facet_integration_points;     // [isrf][spts][dim]
    std::vector<std::vector<double64> >     facet_integration_weights;    // [isrf][spts]
    std::vector<std::vector<double64> >     facet_projection_weights;     // [isrf][spts]
    std::vector<std::vector<double64> >     facet_normals;                // [isrf][dim] //[node*3][dim] -3d
    std::vector<Point<dim> >                facet_parametric_normals;     // [isrf][dim] //[node*3][dim] -3d
    std::vector<std::vector<Point<dim> > >  sector_integration_points;    // [ivol][vpts][dim]
    std::vector<std::vector<double64> >     sector_integration_weights;   // [ivol][vpts]
    std::vector<Point<dim> >                facet_edge_midpoints;
    Point<dim>					          	        barycenter;
    std::vector<std::vector<Point<dim> > >  facet_points;
    // for each sector, all the points delimiting the sector in parametric space
    std::vector<std::vector<Point<dim> > >  sector_points_;         
    // for each sector, for each sector edge, the indices of the end points
    std::vector<std::vector<std::pair<size_t,size_t> > >  sector_edges_;   
    std::string                             parent_element_;
    
    friend class FiniteVolumeStencilManager<dim>;
};





/**

Returns the number of internal facets - division walls inside the FE.

@return Returns the number of the internal facets.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
*/
template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::Facets() const
 {
    return edges_of_element.size();
 }
 
 
 
/**

Returns number of facets, opposite to the given node - this corresponds to
the number of facets, adding flux to the given node.

@return the number of the facets in front of the node of given element type.

@section implementation Implementation

Accesses the private data of the class, retrieving the value.
*/
template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::FacetsPerSector( size_t iSector ) const
 {
    assert( iSector < facets_surrounding_node.size() );

    return facets_surrounding_node[iSector].size();
 }



/**

Returns number of internal FV sectors, given FE is divided to.

@return the number of sectors.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
*/
template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::Sectors() const
 {
    // since this is a vector of vectors of points per sector
    return sector_integration_points.size();
 }



/**

Returns number of integration points per facet.

@return the number of integration points per facet.

@section implementation Implementation 

Accesses the private data of the class, retreiving the value.
For current implementation this value is 1.

*/
template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::IntegrationPointsPerFacet( size_t iFacet ) const
 {
    assert( iFacet < facet_integration_weights.size() );

    // assuming that each facet has the same number of ip's
    return facet_integration_weights[iFacet].size();
 }



/**

Returns number of integration points per sectortric sector of given FE
type.

@return the number of integration points per volumetric sector.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
For current implementation this value is 1.
*/
template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::IntegrationPointsPerSector( size_t iSector ) const
 {
    assert( iSector < sector_integration_weights.size() );

    return sector_integration_weights[iSector].size();
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
template<size_t dim>
inline const Point<dim>&  FiniteVolumeStencil<dim>::FacetIntegrationPoint( size_t iFacet, 
                                                                           size_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    return facet_integration_points[iFacet][ip];
 }



template<size_t dim>
inline void FiniteVolumeStencil<dim>::FacetIntegrationPoint( size_t iFacet, 
                                                             size_t ip, 
                                                             std::vector<double64>& rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    rst = facet_integration_points[iFacet][ip].Coordinates();
 }



// SKM addon
template<size_t dim>
inline double64 FiniteVolumeStencil<dim>::FacetIntegrationPoint( size_t iFacet, 
                                                                 size_t ip, 
                                                                 size_t rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );

    return facet_integration_points[iFacet][ip][rst];
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
template<size_t dim>
inline void FiniteVolumeStencil<dim>::FacetEdgeNodes( size_t iFacet, 
                                                      size_t& estart, 
                                                      size_t& eend ) const
{
    assert( iFacet < Facets() );

    estart=edges_of_element[iFacet].first;
    eend= edges_of_element[iFacet].second;
}


template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::InsideNode( size_t iFacet ) const
{
    assert( iFacet < Facets() );

    return edges_of_element[iFacet].first;
}


template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::OutsideNode( size_t iFacet ) const
{
    assert( iFacet < Facets() );

    return edges_of_element[iFacet].second;
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
template<size_t dim>
inline const Point<dim>&  FiniteVolumeStencil<dim>::SectorIntegrationPoint( size_t iSector, 
                                                                            size_t ip ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    return sector_integration_points[iSector][ip];
 }


template<size_t dim>
inline void FiniteVolumeStencil<dim>::SectorIntegrationPoint( size_t iSector, 
                                                              size_t ip, 
                                                              std::vector<double64>& rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    rst = sector_integration_points[iSector][ip].Coordinates();
 }



template<size_t dim>
inline double64 FiniteVolumeStencil<dim>::SectorIntegrationPoint( size_t iSector, 
                                                                  size_t ip, 
                                                                  size_t rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    return sector_integration_points[iSector][ip][rst];
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
template<size_t dim>
inline double64 FiniteVolumeStencil<dim>::FacetIntegrationWeight( size_t iFacet, 
                                                                  size_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );

    return facet_integration_weights[iFacet][ip];
 }


/**

Returns for the given sector integration point projection weight,
associated with given facet "facet". Projection weight is an extra multiplier,
which is bringing to the same scale the size of the unit parametric
spaces of different element types.

@param facet index (No) of the internal facet  isnside the FE
@param ip index (No) of the integration point on the facet

@return the value of the facet projection weight.

@section implementation Implementation

Accesses the private data of the tabulated points, associated with the facet.

@section application Application

For current implementation index ip is constrained to 0 only, i.e. the
FVPEM method is working with 1 facet integration point only.
*/
 
template<size_t dim>
inline double64 FiniteVolumeStencil<dim>::FacetProjectionWeight( size_t iFacet, size_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    
    return facet_projection_weights[iFacet][ip];
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
template<size_t dim>
inline double64 FiniteVolumeStencil<dim>::SectorIntegrationWeight( size_t iSector, size_t ip ) const
 {
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    return sector_integration_weights[iSector][ip];
 }
 
 
 
  // overloaded method, returns specific index of a facet 
template<size_t dim>
inline size_t  FiniteVolumeStencil<dim>::FacetSurroundingSector( size_t iSector, size_t n ) const
 {
    assert( iSector < Sectors() );
    assert( n < FacetsPerSector(iSector) );

    return facets_surrounding_node[iSector][n];
 }

 
 

 





/**

Returns STL vector of the normal vector coordinates in parametric space of the FE.
 
@param iFacet is the number 0..n-1 of the internal facet inside of the FE

@return to argument "nrml" reference to STL vector of the normal vector
 coordinates in parametric space of the FE. 

@section implementation Implementation
 
Accesses the private data of the class, retreiving tabulated vector.
*/ 
template<size_t dim>
inline const Point<dim>& FiniteVolumeStencil<dim>::UnitParametricNormalTo( size_t iFacet ) const
{
    assert( iFacet < edges_of_element.size() );

    return facet_parametric_normals[iFacet];
}


template<size_t dim>
inline double64 FiniteVolumeStencil<dim>::UnitParametricNormalComponent( size_t iFacet, 
                                                         size_t x_or_y_or_z ) const
{
    assert( iFacet < edges_of_element.size() );

   return facet_parametric_normals[iFacet][x_or_y_or_z];
} 


template<size_t dim>
inline const Point<dim>& FiniteVolumeStencil<dim>::FacetEdgeMidPoint(size_t iFacet ) const
{
    assert( iFacet < edges_of_element.size() );

	return facet_edge_midpoints[iFacet];
}

template<size_t dim>
inline const Point<dim>& FiniteVolumeStencil<dim>::Barycenter() const
{
	 return barycenter;
}
   
   
template<size_t dim>
inline const Point<dim>&  FiniteVolumeStencil<dim>::FacetPoint( size_t iFacet, size_t iPoint ) const
{
	 return facet_points[iFacet][iPoint];
}

    
    
template<size_t dim>
inline const csmp::Point<dim>&  FiniteVolumeStencil<dim>::SectorPoint( size_t iSector, size_t iPoint ) const
{
	 return sector_points_[iSector][iPoint];
}
  
template<size_t dim>
inline std::pair<csmp::Point<dim>,csmp::Point<dim> >  FiniteVolumeStencil<dim>::SectorEdgePoints( size_t iSector, size_t iEdge ) const
{
   const size_t pt1(sector_edges_[iSector][iEdge].first);
   const size_t pt2(sector_edges_[iSector][iEdge].second);
   
	 return std::make_pair( sector_points_[iSector][pt1], sector_points_[iSector][pt2] );
}       

template<size_t dim>
inline const std::pair<size_t,size_t>&  FiniteVolumeStencil<dim>::SectorEdge( size_t iSector, size_t iEdge ) const
{
	 return sector_edges_[iSector][iEdge];
}       
    
template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::SectorPoints( size_t iSector ) const 
{ return sector_points_[iSector].size(); } 

template<size_t dim>
inline size_t FiniteVolumeStencil<dim>::SectorEdges( size_t iSector ) const 
{ return sector_edges_[iSector].size(); } 
    
    
    
} // end namespace csmp

#endif
