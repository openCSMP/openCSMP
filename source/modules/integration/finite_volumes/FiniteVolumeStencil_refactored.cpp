#include "FiniteVolumeStencil_refactored.h"
#include "FiniteElement.h"
#include "FV_IntegrationPointsAndWeights.h"

using namespace std;

namespace csmp { 
 

/**
 
 Parametrised constructor for the FiniteVolumeStencil1 class, setting up
 stencil for given CSMP type of the FE. 

@param csp_finite_element_type is the @enum CSMP_FEM_TYPE

 const char* csp_finite_element_type - pointer to the character array of the FE type.  

@section application Application 

 Standard as for any parametrised constructor. Method parameter should be CSP
 FE type.
*/
template<uint32_t dim>
FiniteVolumeStencil1<dim>::FiniteVolumeStencil1( const char* csp_finite_element_type )
 : parent_element_("not initialized"), parent_element_type_(UNKNOWN)
 {
    Initialize( csp_finite_element_type );
 }




template<uint32_t dim>
FiniteVolumeStencil1<dim>::~FiniteVolumeStencil1()
 {
 }


/**
 
 Parametrised constructor for the FiniteVolumeStencil1 classpoints.

 const FiniteVolumeStencil1<dim>& fvs - reference to the FiniteVolumeStencil1
 object.  

@section application Application 

 Standard as for any parametrised constructor. Method parametr should be previously
 initialised as well.
*/
template<uint32_t dim>
FiniteVolumeStencil1<dim>::FiniteVolumeStencil1( const FiniteVolumeStencil1<dim>& fvs )
 : edges_of_element_(fvs.edges_of_element_),
   facets_surrounding_node_(fvs.facets_surrounding_node_),
   facet_integration_points_(fvs.facet_integration_points_),
   facet_integration_weights_(fvs.facet_integration_weights_),
   facet_normals_(fvs.facet_normals_),
   facet_parametric_normals_(fvs.facet_parametric_normals_),
   facet_normal_xforms_(fvs.facet_normal_xforms_),
   sector_integration_points_(fvs.sector_integration_points_),
   sector_integration_weights_(fvs.sector_integration_weights_),
   parent_element_(fvs.parent_element_),
   parent_element_type_(fvs.parent_element_type_),
   facet_edge_midpoints_(fvs.facet_edge_midpoints_),
   barycenter_(fvs.barycenter_),
   facet_points_(fvs.facet_points_),
   facet_types_(fvs.facet_types_),
   sector_points_(fvs.sector_points_),
   sector_edges_(fvs.sector_edges_),
   space_dimension_(fvs.space_dimension_)
{
}
 
 
/**
 
 Operator "=" for FiniteVolumeStencil1 object.

 const FiniteVolumeStencil1<dim>& fvs - reference to the FiniteVolumeStencil1
 object.  

@section application Application 

 Standard as for any "equals" operator.
*/ 
template<uint32_t dim>
FiniteVolumeStencil1<dim>&  FiniteVolumeStencil1<dim>::operator=( const FiniteVolumeStencil1<dim>& fvs )
 {
     if ( &fvs != this ) {
          edges_of_element_           = fvs.edges_of_element_;
          facets_surrounding_node_    = fvs.facets_surrounding_node_;   // [node][facet]
          facet_integration_points_   = fvs.facet_integration_points_;  // [isrf][spts][dim]
          facet_integration_weights_  = fvs.facet_integration_weights_; // [isrf][spts]
          facet_normals_              = fvs.facet_normals_;             // [isrf][dim]
          facet_normal_xforms_        = fvs.facet_normal_xforms_;       // [isrf][node]
          facet_parametric_normals_   = fvs.facet_parametric_normals_;  // [isrf][dim]
          sector_integration_points_  = fvs.sector_integration_points_;   // [ivol][vpts][dim]
          sector_integration_weights_ = fvs.sector_integration_weights_;  // [ivol][vpts]
          parent_element_             = fvs.parent_element_;
          parent_element_type_        = fvs.parent_element_type_;
          facet_edge_midpoints_       = fvs.facet_edge_midpoints_;
     	    barycenter_                 = fvs.barycenter_;
    	    facet_points_  		          = fvs.facet_points_;
    	    facet_types_  		          = fvs.facet_types_;
    	    sector_points_              = fvs.sector_points_;
     	    sector_edges_               = fvs.sector_edges_;
          space_dimension_            = fvs.space_dimension_;
       }
     return *this;
 }


/**
 
@brief The method is setting up dynamic storage vectors for integration points and weights.

 uint32_t n_isrf          - number of internal facets inside FE;
 uint32_t srfs_per_node   - number of internal facets in front of the node;
 uint32_t n_ivol          - number of internal sectors inside FE;
 uint32_t n_spts          - number of integration points per facet;
 uint32_t n_vpts          - number of integration points per sector;
 
@section implementation Implementation 

 The method is setting up dynamic storage vectors for the tabulated points, 
 associated with the sectors.

@section application Application 

 Can be used at the stencil initialization procedure.

@section messages Messages 
*/
template<uint32_t dim>
void FiniteVolumeStencil1<dim>::Resize( uint32_t n_isrf,
                                        uint32_t srfs_per_node,
                                        uint32_t n_ivol,
                                        uint32_t n_spts,
                                        uint32_t n_vpts )
 {
    facet_integration_points_.resize( n_isrf, n_spts, dim );  // [isrf][spts][dim]
    facet_integration_weights_.resize( n_isrf, n_spts );      // [isrf][spts]
    facet_normals_.resize( n_isrf, dim );
    facet_parametric_normals_.resize(  n_isrf, dim );
    facet_normal_xforms_.resize( n_isrf, n_ivol, 2u );        // [isrf][node][2] (was pair)
    facet_projection_weights_.resize( n_isrf, n_spts );       // [isrf][spts]
    facets_surrounding_node_.resize( n_ivol, srfs_per_node ); // [node][facet], ivol = nodes
    facet_edge_midpoints_.resize( n_isrf, dim );
    facet_points_.resize( n_isrf, n_spts, dim );
    facet_types_.resize( n_isrf );
    edges_of_element_.resize( n_isrf, 2u ); // 2 edges per surface (formerly pair)

 } // end Resize



/**
    Irrespective of the dimension of the model,
    this method reports whether the stencil corresponds to a line, surface or volumetric
    finite element.
*/
template<uint32_t dim>
CELL_SHAPE  FiniteVolumeStencil1<dim>::Geometry() const
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
uint32_t  FiniteVolumeStencil1<dim>::Facets() const
 {
    return static_cast<uint32_t>(edges_of_element_.size());
 }
 
 
 
/**

Returns number of facets, opposite to the given node - this corresponds to
the number of facets, adding flux to the given node.

@return the number of the facets in front of the node of given element type.

@section implementation Implementation

Accesses the private data of the class, retrieving the value.
*/
template<uint32_t dim>
uint32_t  FiniteVolumeStencil1<dim>::FacetsPerSector( uint32_t iSector ) const
 {
    assert( iSector < facets_surrounding_node_.size() );
    
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
                   return 12U;
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
 }



/**

Returns number of internal FV sectors, given FE is divided to.

@return the number of sectors.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
*/
template<uint32_t dim>
uint32_t  FiniteVolumeStencil1<dim>::Sectors() const
 {
    // since this is a vector of vectors of points per sector
    return static_cast<uint32_t>(sector_integration_points_.size());
 }



/**

Returns number of integration points per facet.

@return the number of integration points per facet.

@section implementation Implementation 

Accesses the private data of the class, retreiving the value.
For current implementation this value is 1.

*/
template<uint32_t dim>
uint32_t  FiniteVolumeStencil1<dim>::IntegrationPointsPerFacet( uint32_t /* iFacet */ ) const
 {
    // assert( iFacet < facet_integration_weights_.size() );
    // TODO: check the assumption that each facet has the same number of ip's
    return static_cast<uint32_t>(facet_integration_weights_.cols());
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
uint32_t  FiniteVolumeStencil1<dim>::IntegrationPointsPerSector( uint32_t /* iSector */ ) const
 {
    // assert( iSector < sector_integration_weights_.size() );
    // TODO: check the assumption that each facet has the same number of ip's
    return static_cast<uint32_t>(sector_integration_weights_.cols());
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
const Point<dim>  FiniteVolumeStencil1<dim>::FacetIntegrationPoint( uint32_t iFacet,
                                                                    uint32_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    if constexpr ( dim == 1u )
      return Point<dim>( facet_integration_points_(iFacet,ip,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_integration_points_(iFacet,ip,0u), facet_integration_points_(iFacet,ip,1u) );
    else if constexpr ( dim == 3u )
       return Point<dim>( facet_integration_points_(iFacet,ip,0u),
                          facet_integration_points_(iFacet,ip,1u),
                          facet_integration_points_(iFacet,ip,2u) );
}



template<uint32_t dim>
void FiniteVolumeStencil1<dim>::FacetIntegrationPoint( uint32_t iFacet,
                                                       uint32_t ip,
                                                       vector<double>& rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    rst.resize(dim);

    if constexpr ( dim == 1u )
      rst[0u] = facet_integration_points_(iFacet,ip,0u);
    else if constexpr ( dim == 2u ) {
        rst[0u] = facet_integration_points_(iFacet,ip,0u);
        rst[1u] = facet_integration_points_(iFacet,ip,1u);
      }
    else if constexpr ( dim == 3u ) {
       rst[0u] = facet_integration_points_(iFacet,ip,0u);
       rst[1u] = facet_integration_points_(iFacet,ip,1u);
       rst[2u] = facet_integration_points_(iFacet,ip,2u);
     }
 }



// SKM addon
template<uint32_t dim>
double FiniteVolumeStencil1<dim>::FacetIntegrationPoint( uint32_t iFacet,
                                                         uint32_t ip,
                                                         uint32_t rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    return facet_integration_points_(iFacet,ip,rst);
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
void FiniteVolumeStencil1<dim>::FacetEdgeNodes( uint32_t iFacet,
                                               uint32_t& estart,
                                               uint32_t& eend ) const
{
    assert( iFacet < Facets() );

    estart=edges_of_element_(iFacet,0u);
    eend=edges_of_element_(iFacet,1U);
}


template<uint32_t dim>
uint32_t FiniteVolumeStencil1<dim>::InsideNode( uint32_t iFacet ) const
{
    assert( iFacet < Facets() );
    return edges_of_element_(iFacet,0u);
}


template<uint32_t dim>
uint32_t FiniteVolumeStencil1<dim>::OutsideNode( uint32_t iFacet ) const
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
const Point<dim>  FiniteVolumeStencil1<dim>::SectorIntegrationPoint( uint32_t iSector,
                                                                     uint32_t ip ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
    if constexpr ( dim == 1u )
      return Point<dim>( sector_integration_points_(iSector,ip,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( sector_integration_points_(iSector,ip,0u), sector_integration_points_(iSector,ip,1u) );
    else if constexpr ( dim == 3u )
       return Point<dim>( sector_integration_points_(iSector,ip,0u),
                          sector_integration_points_(iSector,ip,1u),
                          sector_integration_points_(iSector,ip,2u) );
 }


template<uint32_t dim>
void FiniteVolumeStencil1<dim>::SectorIntegrationPoint( uint32_t iSector,
                                                        uint32_t ip,
                                                        vector<double>& rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
    rst.resize(dim);
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
 }



template<uint32_t dim>
double FiniteVolumeStencil1<dim>::SectorIntegrationPoint( uint32_t iSector,
                                                          uint32_t ip,
                                                          uint32_t rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
    return sector_integration_points_(iSector,ip,rst);
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
double FiniteVolumeStencil1<dim>::FacetIntegrationWeight( uint32_t iFacet,
                                                          uint32_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    return facet_integration_weights_(iFacet,ip);
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
double FiniteVolumeStencil1<dim>::FacetProjectionWeight( uint32_t iFacet, uint32_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    return facet_projection_weights_(iFacet,ip);
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
double FiniteVolumeStencil1<dim>::SectorIntegrationWeight( uint32_t iSector, uint32_t ip ) const
 {
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );
    return sector_integration_weights_(iSector,ip);
 }
 
 
 
  // overloaded method, returns specific index of a facet 
template<uint32_t dim>
uint32_t  FiniteVolumeStencil1<dim>::FacetSurroundingSector( uint32_t iSector, uint32_t n ) const
 {
    assert( iSector < Sectors() );
    assert( n < FacetsPerSector(iSector) );
    return facets_surrounding_node_(iSector,n);
 }

 


/**

Returns the normal vector in parametric space of the FE.
 
@param iFacet is the number 0..n-1 of the internal facet inside of the FE

@return to argument "nrml" reference to STL vector of the normal vector
 coordinates in parametric space of the FE. 

@section implementation Implementation
 
Accesses the private data of the class, retreiving tabulated vector.
*/ 
template<uint32_t dim>
const Point<dim>  FiniteVolumeStencil1<dim>::UnitParametricNormalTo( uint32_t iFacet ) const
{
    assert( iFacet < edges_of_element_.size() );
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
double FiniteVolumeStencil1<dim>::UnitParametricNormalComponent( uint32_t iFacet,
                                                                 uint32_t x_or_y_or_z ) const
{
   assert( iFacet < edges_of_element_.size() );
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
pair<double,double>  FiniteVolumeStencil1<dim>::FacetNormalTransformationNodeWeights( uint32_t iFacet, uint32_t iNode ) const
{
    assert( iFacet < facet_normal_xforms_.size() );
    return make_pair( facet_normal_xforms_(iFacet,iNode,0u), facet_normal_xforms_(iFacet,iNode,1u) );
}


template<uint32_t dim>
const Point<dim>  FiniteVolumeStencil1<dim>::FacetEdgeMidPoint( uint32_t iFacet ) const
 {
    assert( iFacet < edges_of_element_.size() );
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
const Point<dim>  FiniteVolumeStencil1<dim>::Barycenter() const
{
	 return barycenter_;
}


template<uint32_t dim>
const Point<dim>  FiniteVolumeStencil1<dim>::FacetPoint( uint32_t iFacet, uint32_t iPoint ) const
 {
    assert( iFacet < facet_points_.rows() );
    assert( iPoint < facet_points_.cols() );
    if constexpr ( dim == 1u )
      return Point<dim>( facet_points_(iFacet,iPoint,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( facet_points_(iFacet,iPoint,0u), facet_points_(iFacet,iPoint,1u) );
    else if constexpr ( dim == 3u )
      return Point<dim>( facet_points_(iFacet,iPoint,0u),
                         facet_points_(iFacet,iPoint,1u),
                         facet_points_(iFacet,iPoint,2u) );
 }


template<uint32_t dim>
FV_FACET_TYPE  FiniteVolumeStencil1<dim>::FacetType( uint32_t iFacet ) const
{
   assert( iFacet < facet_types_.size() );
	 return facet_types_[iFacet];
}

    
    
template<uint32_t dim>
const csmp::Point<dim>  FiniteVolumeStencil1<dim>::SectorPoint( uint32_t iSector, uint32_t iPoint ) const
 {
    assert( iSector < sector_points_.rows() );
    assert( iPoint < sector_points_.cols() );
    if constexpr ( dim == 1u )
      return Point<dim>( sector_points_(iSector,iPoint,0u) );
    else if constexpr ( dim == 2u )
      return Point<dim>( sector_points_(iSector,iPoint,0u), sector_points_(iSector,iPoint,1u) );
    else if constexpr ( dim == 3u )
      return Point<dim>( sector_points_(iSector,iPoint,0u),
                         sector_points_(iSector,iPoint,1u),
                         sector_points_(iSector,iPoint,2u) );
 }
  
  
  
template<uint32_t dim>
std::pair<csmp::Point<dim>,csmp::Point<dim> >  FiniteVolumeStencil1<dim>::SectorEdgePoints( uint32_t iSector, uint32_t iEdge ) const
{
    assert( iSector < sector_edges_.rows() );
    assert( iEdge < sector_edges_.cols() );
    const uint32_t pt1(sector_edges_(iSector,iEdge,0u));
    const uint32_t pt2(sector_edges_(iSector,iEdge,1u));
   
    if constexpr ( dim == 1u )
      return make_pair( Point<1U>(sector_points_(iSector,pt1,0u)), Point<1U>(sector_points_(iSector,pt2,0u)) );
    else if constexpr ( dim == 2u )
      return make_pair( Point<2U>(sector_points_(iSector,pt1,0u),sector_points_(iSector,pt1,1u)),
                        Point<2U>(sector_points_(iSector,pt2,0u),sector_points_(iSector,pt2,1u)) );
    else if constexpr ( dim == 3u )
      return make_pair( Point<3U>(sector_points_(iSector,pt1,0u),sector_points_(iSector,pt1,1u), sector_points_(iSector,pt1,2u)),
                        Point<3U>(sector_points_(iSector,pt2,0u),sector_points_(iSector,pt2,1u), sector_points_(iSector,pt2,2u)) );
}



template<uint32_t dim>
const pair<uint32_t,uint32_t>  FiniteVolumeStencil1<dim>::SectorEdge( uint32_t iSector, uint32_t iEdge ) const
{
	 return make_pair( sector_edges_(iSector,iEdge,1u), sector_edges_(iSector,iEdge,1u) );
}
    
    
template<uint32_t dim>
uint32_t FiniteVolumeStencil1<dim>::SectorPoints( uint32_t iSector ) const
// TODO: check assumption that all sectors have the same number of points!
{ return static_cast<uint32_t>(sector_points_.rows()); }


template<uint32_t dim>
uint32_t FiniteVolumeStencil1<dim>::SectorEdges( uint32_t iSector ) const
// TODO: check assumption that all sectors have the same number of edges!
{ return static_cast<uint32_t>(sector_edges_.depth()*sector_edges_.rows()); }
    
  




/**

 This method is the main workhorse that sets up the stencil for the specific 
 finite element type.  

@param csp_finite_element_type - pointer to the character array of the FE type.

@section implementation Implementation 

 Resizes and defines the private data of the tabulated points, weights and normals, 
 associated with the specific type of the finite element type.

@section application Application 

 Should be called at the moment of stencil creation for specific type of the 
 finite element within CSP. Alternatively, could be used for re-definition of
 internal FV stencil data as well.  
 
 Method uses CSP object FV_IntegrationPointsAndWeights to generate tabulated
 integration points, weights, normals for any type of element. Currently, works
 only with linear isoparametric finite elements.  
 
 
@section messages Messages 
*/
template<uint32_t dim>
void FiniteVolumeStencil1<dim>::Initialize( const char* csp_finite_element_type )
 {
   // uses the types specified in the FiniteElement.h header file
   const string csp_fem_type(csp_finite_element_type);
   bool debug=false;

   uint32_t NumOfInternalFacets=6U;
   uint32_t NumOfInternalFacetsPerNode=3U;
   uint32_t NumOfInternalVolumes=4U;
   uint32_t NumOfIPperVolume=1U;
   uint32_t NumOfIPperFacet=1U;
   
   if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_BAR) )
     {
          NumOfInternalFacets=1U;
          NumOfInternalFacetsPerNode=1U;
          NumOfInternalVolumes=2U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );

          FV_IntegrationPointsAndWeights<dim> FV(ISOPARAMETRIC_LINEAR_BAR);
          InitializeDataMembers( FV );
  
          parent_element_      = "ISOPARAMETRIC_LINEAR_BAR";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_BAR;
          space_dimension_     = LINE;

          if(debug) {
             cout<<"\nISOPARAMETRIC_LINEAR_BAR:" <<endl;
             Out();         
           }
          return;
      }
      
   if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_TRIANGLE) ) {

          NumOfInternalFacets=3U;
          NumOfInternalFacetsPerNode=2U;
          NumOfInternalVolumes=3U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );

             
          FV_IntegrationPointsAndWeights<dim> FV(ISOPARAMETRIC_LINEAR_TRIANGLE);
          InitializeDataMembers( FV );
             
          parent_element_ = "ISOPARAMETRIC_LINEAR_TRIANGLE";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_TRIANGLE;
          space_dimension_ = SURFACE;
     
          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_TRIANGLE:" <<endl;
              Out();         
          }
           
          return;
      }
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_TETRAHEDRON) ) {

         NumOfInternalFacets=6U;
         NumOfInternalFacetsPerNode=3U;
         NumOfInternalVolumes=4U;

         Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );

         FV_IntegrationPointsAndWeights<dim> FVT(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
         InitializeDataMembers( FVT );
             
         parent_element_ = "ISOPARAMETRIC_LINEAR_TETRAHEDRON";
         parent_element_type_ = ISOPARAMETRIC_LINEAR_TETRAHEDRON;
         space_dimension_ = VOLUME;

             if(debug) {
                cout<<"\nISOPARAMETRIC_LINEAR_TETRAHEDRON: "<<endl;
                Out();         
             }
          return;
     }
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_QUADRILATERAL) ) {

         NumOfInternalFacets=4U;
         NumOfInternalFacetsPerNode=2U;
         NumOfInternalVolumes=4U; // == VolMultipliers

         Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVQ(ISOPARAMETRIC_LINEAR_QUADRILATERAL);
          InitializeDataMembers( FVQ );

          parent_element_      = "ISOPARAMETRIC_LINEAR_QUADRILATERAL";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_QUADRILATERAL;
          space_dimension_     = SURFACE;

          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_QUADRILATERAL:"<<endl;
              Out();         
           }
         return;
      }
      
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_HEXAHEDRON) ) {

          NumOfInternalFacets=12U;
          NumOfInternalFacetsPerNode=3U;
          NumOfInternalVolumes=8U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVH(ISOPARAMETRIC_LINEAR_HEXAHEDRON);
          InitializeDataMembers( FVH );

          parent_element_      = "ISOPARAMETRIC_LINEAR_HEXAHEDRON";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_HEXAHEDRON;
          space_dimension_     = VOLUME;

         if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_HEXAHEDRON:"<<endl;
              Out();         
          }
         return;
      }
      
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_PYRAMID) ) {

 #ifdef PYRAMID_TRIANGULAR_FACETS
          NumOfInternalFacets=12U;
          NumOfInternalFacetsPerNode=8U; // Special case of the Pyramid vertex require 8
          NumOfInternalVolumes=5U;       // == VolMultipliers
 #else
          NumOfInternalFacets=8U;
          NumOfInternalFacetsPerNode=4U; // Special case of the Pyramid vertex require 4
          NumOfInternalVolumes=5U;       // == VolMultipliers
 #endif
          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVPy(ISOPARAMETRIC_LINEAR_PYRAMID);
          InitializeDataMembers( FVPy );

          parent_element_ = "ISOPARAMETRIC_LINEAR_PYRAMID";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_PYRAMID;
          space_dimension_ = VOLUME;

         if ( debug ) {
              cout<<"\nISOPARAMETRIC_LINEAR_PYRAMID:"<<endl;
              Out();         
        }
       return;
       
    }
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_PRISM) ) {

          NumOfInternalFacets=9U;
          NumOfInternalFacetsPerNode=3U;
          NumOfInternalVolumes=6U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVP(ISOPARAMETRIC_LINEAR_PRISM);
          InitializeDataMembers( FVP );
		  
          parent_element_      = "ISOPARAMETRIC_LINEAR_PRISM";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_PRISM;
          space_dimension_     = VOLUME;

          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_PRISM:"<<endl;
              Out();         
           }
          return;
      }
    else {
         cerr <<"\nFiniteVolumeStencil1::Initialize: Element type not recognized: ";
         cerr << csp_finite_element_type << endl;
      }
      
 } // end Initialize





template<uint32_t dim>
void FiniteVolumeStencil1<dim>::InitializeDataMembers( const FV_IntegrationPointsAndWeights<dim>& FV )
 {

      // facet integration weights
      {
         vector<vector<double>>  facet_integration_weights;
         FV.FacetIntegrationWeights(facet_integration_weights);
         facet_integration_weights_.resize( facet_integration_weights.size(), facet_integration_weights[0].size() );
         uint32_t i{0u};
         for ( const auto& iit : facet_integration_weights ) {
             uint32_t j{0u};
             for ( const auto& jit : iit )
               facet_integration_weights_(i,++j) = jit;
             ++i;
          }
      }

      // TODO: missing: vector<std::vector<double> >  facet_projection_weights;
      // facet projection weights
      /*
      {
         vector<vector<double>>  facet_projection_weights;
         FV.FacetProjectionWeights(facet_projection_weights);
         facet_projection_weights_.resize( facet_projection_weights.size(), facet_projection_weights[0].size() );
         uint32_t i{0u};
         for ( const auto& iit : facet_projection_weights ) {
             uint32_t j{0u};
             for ( const auto& jit : iit )
               facet_projection_weights_(i,++j) = jit;
             ++i;
          }
      }
      */

      // facet normals
      {
         vector<Point<dim>>  facet_normals; // TODO: why do these get initialised here?
         FV.FacetNormals(facet_normals);
         facet_normals_.resize( facet_normals.size(), dim );
         uint32_t i{0u};
         for ( const auto& iit : facet_normals ) {
             for ( uint32_t j{0u}; j<dim; ++j )
               facet_normals_(i,j) = iit[j];
             ++i;
          }
      }

      // facet parametric normals
      {
         vector<Point<dim>>  facet_parametric_normals;
         FV.FacetNormals(facet_parametric_normals);
         facet_parametric_normals_.resize( facet_parametric_normals.size(), dim );
         uint32_t i{0u};
         for ( const auto& iit : facet_parametric_normals ) {
             for ( uint32_t j{0u}; j<dim; ++j )
               facet_parametric_normals_(i,j) = iit[j];
             ++i;
          }
      }

      // facet edge midpoints
      {
         vector<Point<dim>> facet_edge_midpoints;
         FV.EdgeMidpoints(facet_edge_midpoints);
         facet_edge_midpoints_.resize( facet_edge_midpoints.size(), dim );
         uint32_t i{0u};
         for ( const auto& iit : facet_edge_midpoints ) {
             for ( uint32_t j{0u}; j<dim; ++j )
               facet_edge_midpoints_(i,j) = iit[j];
             ++i;
          }
      }

      // sector integration weights
      {
         vector<vector<double> > sector_integration_weights;
         FV.SectorIntegrationWeights(sector_integration_weights);
         sector_integration_weights_.resize( sector_integration_weights.size(), sector_integration_weights[0].size() );
         uint32_t i{0u};
         for ( const auto& iit : sector_integration_weights ) {
             uint32_t j{0u};
             for ( const auto& jit : iit )
               sector_integration_weights_(i,++j) = jit;
             ++i;
          }
      }
      
    // THREE-DIMENSIONAL ARRAYS
     
     // facet transformations (parametric to physical)
     {
        vector<vector<pair<double,double>>>  facet_normal_xforms;
        FV.FacetNormalTransformations(facet_normal_xforms);
         facet_normal_xforms_.resize( facet_normal_xforms.size(),
                                      facet_normal_xforms[0].size(), 2u );
         uint32_t i{0u};
         for ( const auto& iit : facet_normal_xforms ) {
             uint32_t j{0u};
             for ( const auto& jit : iit ) {
                   facet_normal_xforms_(i,j,0u) = jit.first;
                   facet_normal_xforms_(i,j,1u) = jit.second;
                 ++j;
               }
             ++i;
          }
     }
      
     // facet points
     {
        vector<vector<Point<dim>>>  facet_points;
        FV.FacetPoints(facet_points);
        facet_points_.resize( facet_points.size(), facet_points[0].size(), dim );
        uint32_t i{0u};
        for ( const auto& iit : facet_points ) {
            uint32_t j{0u};
            for ( const auto& jit : iit ) {
                for ( uint32_t k{0u}; k<dim; ++k )
                  facet_points_(i,j,k) = jit[j];
                ++j;
              }
            ++i;
         }
      }

     // facet integration points
     {
        vector<vector<Point<dim>>>  facet_integration_points;
        FV.FacetIntegrationPoints(facet_integration_points);
        facet_integration_points_.resize( facet_integration_points.size(),
                                          facet_integration_points[0].size(), dim );
        uint32_t i{0u};
        for ( const auto& iit : facet_integration_points ) {
            uint32_t j{0u};
            for ( const auto& jit : iit ) {
                for ( uint32_t k{0u}; k<dim; ++k )
                  facet_integration_points_(i,j,k) = jit[j];
                ++j;
              }
            ++i;
         }
      }

     // sector points
     {
        vector<vector<Point<dim>>>  sector_points;
        FV.SectorPoints(sector_points);
        sector_points_.resize( sector_points.size(), sector_points[0].size(), dim );
        uint32_t i{0u};
        for ( const auto& iit : sector_points ) {
            uint32_t j{0u};
            for ( const auto& jit : iit ) {
                for ( uint32_t k{0u}; k<dim; ++k )
                  sector_points_(i,j,k) = jit[j];
                ++j;
              }
            ++i;
         }
      }

     // sector integration points
     {
        vector<vector<Point<dim>>>  sector_integration_points;
        FV.SectorIntegrationPoints(sector_integration_points);
        sector_integration_points_.resize( sector_integration_points.size(),
                                           sector_integration_points[0].size(), dim );
        uint32_t i{0u};
        for ( const auto& iit : sector_integration_points ) {
            uint32_t j{0u};
            for ( const auto& jit : iit ) {
                for ( uint32_t k{0u}; k<dim; ++k )
                  sector_integration_points_(i,j,k) = jit[j];
                ++j;
              }
            ++i;
         }
      }

      // edges of elements
      {
         vector<pair<uint32_t,uint32_t> >  edges_of_element;
         FV.EdgePairs(edges_of_element);
         edges_of_element_.resize( edges_of_element.size(), 2u );
         uint32_t i{0u};
         for ( const auto& iit : edges_of_element ) {
               edges_of_element_(i,0u) = iit.first;
               edges_of_element_(i,1u) = iit.second;
             ++i;
          }
      }

      // facets surrounding node
      {
         vector<vector<uint32_t> >  facets_surrounding_node;
         FV.FacetsSurroundingNode(facets_surrounding_node);
         facets_surrounding_node_.resize( facets_surrounding_node.size(), facets_surrounding_node[0].size() );
         uint32_t i{0u};
         for ( const auto& iit : facets_surrounding_node ) {
             uint32_t j{0u};
             for ( const auto& jit : iit )
               sector_integration_weights_(i,++j) = jit;
             ++i;
          }
      }

     // sector edges
     {
        vector<vector<pair<uint32_t,uint32_t> > >  sector_edges;
        FV.SectorEdgePairs(sector_edges);
        sector_edges_.resize( sector_edges.size(), sector_edges[0].size(), 2u );
        uint32_t i{0u};
        for ( const auto& iit : sector_edges ) {
            uint32_t j{0u};
            for ( const auto& jit : iit ) {
                  sector_edges_(i,j,0u) = jit.first;
                  sector_edges_(i,j,1u) = jit.second;
                ++j;
              }
            ++i;
         }
      }

    FV.Barycenter(barycenter_);
    FV.FacetTypes(facet_types_);
          
 } // InitializeDataMembers




// SKM modified version TODO: assumption is that the arrays do not have different row lenghts
template<uint32_t dim>
void FiniteVolumeStencil1<dim>::Out() const
{
   cout<<"\nFiniteVolumeStencil1<dim>::Out: \nVolume Weights ";
   cout <<"("<< sector_integration_weights_.rows() <<","<< sector_integration_weights_.cols() <<"):";
   for( uint32_t i{0U}; i<sector_integration_weights_.size(); i++)
      for( uint32_t j{0U}; j<sector_integration_weights_.cols(); j++)
        cout<<" "<< sector_integration_weights_(i,j);
   
   cout<<endl<<endl<<" Volume IPs "<<endl;
   cout <<"("<< sector_integration_points_.rows() <<","<< sector_integration_points_.cols() <<",rst):";
   for( uint32_t i{0U}; i<sector_integration_points_.rows(); i++) {
     for( uint32_t j{0U}; j<sector_integration_points_.cols(); j++) {
        for( uint32_t k{0U}; k<dim; k++) cout<<" "<< sector_integration_points_(i,j,k);
          cout<<endl;
        }
   }
                
   cout<<endl<<" Facet Weights ";
   cout <<"("<< facet_integration_weights_.rows() <<","<< facet_integration_weights_.cols() <<"):";
   for( uint32_t i{0U}; i<facet_integration_weights_.rows(); i++) {
     for( uint32_t j{0U}; j<facet_integration_weights_.cols(); j++) {
       cout<<" "<< facet_integration_weights_(i,j);
      }
   }
   cout<<endl;
   
   cout<<endl<<" Facet IPs: "<<endl;
   cout <<"("<< facet_integration_points_.rows() <<","<< facet_integration_points_.cols() <<",rst):";
   for( uint32_t i{0U}; i<facet_integration_points_.rows(); i++) {
     for( uint32_t j{0U}; j<facet_integration_points_.cols(); j++) {
       for( uint32_t k{0U}; k<dim; k++) cout<<" "<< facet_integration_points_(i,j,k);
       cout<<endl;
     }
   }

  cout<<endl<<" Facet parametric normals inside the element "<<endl;
  cout <<"("<< facet_parametric_normals_.rows() <<",rst):";
  for( uint32_t i{0U}; i<facet_parametric_normals_.rows(); i++) {
    for( uint32_t j{0U}; j<dim; j++) {
       cout<<" "<< facet_parametric_normals_(i,j);
     }
     cout<<endl;
   }
     
   if ( !facet_normals_.empty() && facet_normals_.cols() == 0U ) {
       cout<<endl<<" Facet physical normals inside the element (may not be initialised): "<<endl;
       cout <<"("<< facet_normals_.rows() <<",xyz):";
       for( uint32_t i{0U}; i<facet_normals_.rows(); i++) {
         for( uint32_t j{0U}; j<dim; j++) {
             cout<<" "<< facet_normals_(i,j);
           }
         cout<<endl;
       }
     }

    cout<<endl<<" Facet surrounding node (corresponding to sector) in the element "<<endl;
    cout <<"("<< Sectors() <<","<< dim <<"):";
    for( uint32_t i{0U}; i<Sectors(); i++) {
      for( uint32_t j{0U}; j<FacetsPerSector(i); j++) {
        cout<<" "<< FacetSurroundingSector(i,j);
      }
      cout<<endl;
    }
    
    cout<<endl<<" Edges pairs in the element "<<endl;
    cout <<"("<< edges_of_element_.rows() <<",2):";
    for( uint32_t i{0U}; i<edges_of_element_.rows(); i++) {
      cout<<" "<< edges_of_element_(i,0U) <<", "<< edges_of_element_(i,1U) <<" ";
    }
  cout << endl << endl; 
   
} // end Out


template class FiniteVolumeStencil1<1U>; 
template class FiniteVolumeStencil1<2U>; 
template class FiniteVolumeStencil1<3U>; 


} // end namespace csmp
