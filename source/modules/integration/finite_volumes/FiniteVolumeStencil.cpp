#include "FiniteVolumeStencil.h"
#include "FiniteElement.h"
#include "FV_IntegrationPointsAndWeights.h"

using namespace std;

namespace csmp { 
 

/**
 
 Parametrised constructor for the FiniteVolumeStencil class, setting up
 stencil for given CSMP type of the FE. 

@param csp_finite_element_type is the @enum CSMP_FEM_TYPE

 const char* csp_finite_element_type - pointer to the character array of the FE type.  

@section application Application 

 Standard as for any parametrised constructor. Method parameter should be CSP
 FE type.
*/
template<size_t dim>
FiniteVolumeStencil<dim>::FiniteVolumeStencil( const char* csp_finite_element_type )
 : parent_element_("not initialized")
 {
    Initialize( csp_finite_element_type );
 }




template<size_t dim>
FiniteVolumeStencil<dim>::~FiniteVolumeStencil()
 {
 }


/**
 
 Parametrised constructor for the FiniteVolumeStencil classpoints.

 const FiniteVolumeStencil<dim>& fvs - reference to the FiniteVolumeStencil
 object.  

@section application Application 

 Standard as for any parametrised constructor. Method parametr should be previously
 initialised as well.
*/
template<size_t dim>
FiniteVolumeStencil<dim>::FiniteVolumeStencil( const FiniteVolumeStencil<dim>& fvs )
 : edges_of_element(fvs.edges_of_element),
   facets_surrounding_node(fvs.facets_surrounding_node),
   facet_integration_points(fvs.facet_integration_points),
   facet_integration_weights(fvs.facet_integration_weights),
   facet_normals(fvs.facet_normals),
   facet_parametric_normals(fvs.facet_parametric_normals),
   facet_normal_xforms(fvs.facet_normal_xforms),
   sector_integration_points(fvs.sector_integration_points),
   sector_integration_weights(fvs.sector_integration_weights),
   parent_element_(fvs.parent_element_),
   facet_edge_midpoints(fvs.facet_edge_midpoints),
   barycenter(fvs.barycenter),
   facet_points(fvs.facet_points),
   facet_types(fvs.facet_types),
   sector_points_(fvs.sector_points_),
   sector_edges_(fvs.sector_edges_),
   space_dimension_(fvs.space_dimension_)
{
}
 
 
/**
 
 Operator "=" for FiniteVolumeStencil object.

 const FiniteVolumeStencil<dim>& fvs - reference to the FiniteVolumeStencil
 object.  

@section application Application 

 Standard as for any "equals" operator.
*/ 
template<size_t dim>
FiniteVolumeStencil<dim>&  FiniteVolumeStencil<dim>::operator=( const FiniteVolumeStencil<dim>& fvs )
 {
     if ( &fvs != this ) {
          edges_of_element           = fvs.edges_of_element;
          facets_surrounding_node    = fvs.facets_surrounding_node;   // [node][facet]
          facet_integration_points   = fvs.facet_integration_points;  // [isrf][spts][dim]
          facet_integration_weights  = fvs.facet_integration_weights; // [isrf][spts]
          facet_normals              = fvs.facet_normals;             // [isrf][dim]
          facet_normal_xforms        = fvs.facet_normal_xforms;       // [isrf][node]
          facet_parametric_normals   = fvs.facet_parametric_normals;  // [isrf][dim]
          sector_integration_points  = fvs.sector_integration_points;   // [ivol][vpts][dim]
          sector_integration_weights = fvs.sector_integration_weights;  // [ivol][vpts]
          parent_element_            = fvs.parent_element_;
          facet_edge_midpoints       = fvs.facet_edge_midpoints;
     	    barycenter                 = fvs.barycenter;
    	    facet_points  		         = fvs.facet_points;
    	    facet_types  		         = fvs.facet_types;
    	    sector_points_             = fvs.sector_points_;
     	    sector_edges_              = fvs.sector_edges_;
          space_dimension_           = fvs.space_dimension_;
       }
     return *this;
 }


/**
 
@brief The method is setting up dynamic storage vectors for integration points and weights.

 size_t n_isrf          - number of internal facets inside FE;
 size_t srfs_per_node   - number of internal facets in front of the node;
 size_t n_ivol          - number of internal sectors inside FE; 
 size_t n_spts          - number of integration points per facet;
 size_t n_vpts          - number of integration points per sector;
 
@section implementation Implementation 

 The method is setting up dynamic storage vectors for the tabulated points, 
 associated with the sectors.

@section application Application 

 Can be used at the stencil initialization procedure.

@section messages Messages 
*/
template<size_t dim>
void FiniteVolumeStencil<dim>::Resize( size_t n_isrf,
                                        size_t srfs_per_node,
                                        size_t n_ivol, 
                                        size_t n_spts, 
                                        size_t n_vpts
                                     )
 {
    facet_integration_points.resize( n_isrf );
    facet_integration_weights.resize( n_isrf );
    facet_normals.resize( n_isrf );
    facet_normal_xforms.resize( n_isrf );
    facet_parametric_normals.resize( n_isrf );
    facet_projection_weights.resize( n_isrf );    
    facets_surrounding_node.resize( n_ivol ); // ivol = nodes
    facet_edge_midpoints.resize( n_isrf );
    facet_points.resize( n_isrf );
    facet_types.resize( n_isrf );
    edges_of_element.resize(n_isrf);

    for ( size_t i=0; i<n_ivol; i++ ) 
      facets_surrounding_node[i].resize( srfs_per_node );
    
    for ( size_t i=0; i<n_isrf; i++ ) {
         facet_integration_points[i].resize( n_spts );
         facet_projection_weights[i].resize( n_spts );
         facet_integration_weights[i].resize( n_spts );
         facet_normal_xforms[i].resize( n_ivol ); // ivol = nodes
      } 
    
    sector_integration_points.resize( n_ivol );
    sector_integration_weights.resize( n_ivol );
    //par_sector_integration_points.resize( n_ivol );

    for ( size_t i=0; i<n_ivol; i++ ) {
         sector_integration_points[i].resize( n_vpts );
         //par_sector_integration_points[i].resize( n_vpts );
           //par_sector_integration_points[i][j].resize(pdim);
         sector_integration_weights[i].resize(n_vpts);
      }

 } // end Resize



/**
    Irrespective of the dimension of the model,
    this method reports whether the stencil corresponds to a line, surface or volumetric
    finite element.
*/
template<size_t dim>
ELEMENT_DIMENSION FiniteVolumeStencil<dim>::Geometry() const
 {
    return space_dimension_;
 }



/**

Returns the number of internal facets - division walls inside the FE.

@return Returns the number of the internal facets.

@section implementation Implementation

Accesses the private data of the class, retreiving the value.
*/
template<size_t dim>
size_t FiniteVolumeStencil<dim>::Facets() const
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
size_t FiniteVolumeStencil<dim>::FacetsPerSector( size_t iSector ) const
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
size_t FiniteVolumeStencil<dim>::Sectors() const
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
size_t FiniteVolumeStencil<dim>::IntegrationPointsPerFacet( size_t iFacet ) const
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
size_t FiniteVolumeStencil<dim>::IntegrationPointsPerSector( size_t iSector ) const
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
const Point<dim>&  FiniteVolumeStencil<dim>::FacetIntegrationPoint( size_t iFacet, 
                                                                           size_t ip ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    return facet_integration_points[iFacet][ip];
 }



template<size_t dim>
void FiniteVolumeStencil<dim>::FacetIntegrationPoint( size_t iFacet, 
                                                             size_t ip, 
                                                             std::vector<double64>& rst ) const
 {
    assert( iFacet < Facets() );
    assert( ip < IntegrationPointsPerFacet(iFacet) );
    rst = facet_integration_points[iFacet][ip].Coordinates();
 }



// SKM addon
template<size_t dim>
double64 FiniteVolumeStencil<dim>::FacetIntegrationPoint( size_t iFacet, 
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
void FiniteVolumeStencil<dim>::FacetEdgeNodes( size_t iFacet, 
                                                      size_t& estart, 
                                                      size_t& eend ) const
{
    assert( iFacet < Facets() );

    estart=edges_of_element[iFacet].first;
    eend= edges_of_element[iFacet].second;
}


template<size_t dim>
size_t FiniteVolumeStencil<dim>::InsideNode( size_t iFacet ) const
{
    assert( iFacet < Facets() );

    return edges_of_element[iFacet].first;
}


template<size_t dim>
size_t FiniteVolumeStencil<dim>::OutsideNode( size_t iFacet ) const
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
const Point<dim>&  FiniteVolumeStencil<dim>::SectorIntegrationPoint( size_t iSector, 
                                                                            size_t ip ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    return sector_integration_points[iSector][ip];
 }


template<size_t dim>
void FiniteVolumeStencil<dim>::SectorIntegrationPoint( size_t iSector, 
                                                              size_t ip, 
                                                              std::vector<double64>& rst ) const
{
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    rst = sector_integration_points[iSector][ip].Coordinates();
 }



template<size_t dim>
double64 FiniteVolumeStencil<dim>::SectorIntegrationPoint( size_t iSector, 
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
double64 FiniteVolumeStencil<dim>::FacetIntegrationWeight( size_t iFacet, 
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
double64 FiniteVolumeStencil<dim>::FacetProjectionWeight( size_t iFacet, size_t ip ) const
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
double64 FiniteVolumeStencil<dim>::SectorIntegrationWeight( size_t iSector, size_t ip ) const
 {
    assert( iSector < Sectors() );
    assert( ip < IntegrationPointsPerSector(iSector) );

    return sector_integration_weights[iSector][ip];
 }
 
 
 
  // overloaded method, returns specific index of a facet 
template<size_t dim>
size_t  FiniteVolumeStencil<dim>::FacetSurroundingSector( size_t iSector, size_t n ) const
 {
    assert( iSector < Sectors() );
    assert( n < FacetsPerSector(iSector) );

    return facets_surrounding_node[iSector][n];
 }

 
 

 





/**

Returns the normal vector in parametric space of the FE.
 
@param iFacet is the number 0..n-1 of the internal facet inside of the FE

@return to argument "nrml" reference to STL vector of the normal vector
 coordinates in parametric space of the FE. 

@section implementation Implementation
 
Accesses the private data of the class, retreiving tabulated vector.
*/ 
template<size_t dim>
const Point<dim>& FiniteVolumeStencil<dim>::UnitParametricNormalTo( size_t iFacet ) const
{
    assert( iFacet < edges_of_element.size() );

    return facet_parametric_normals[iFacet];
}


template<size_t dim>
double64 FiniteVolumeStencil<dim>::UnitParametricNormalComponent( size_t iFacet, 
                                                                         size_t x_or_y_or_z ) const
{
    assert( iFacet < edges_of_element.size() );

   return facet_parametric_normals[iFacet][x_or_y_or_z];
} 


/**

Returns the component of the transformation from nodes to normals in physical space.
 
@param iFacet is the number 0..n-1 of the internal facet inside of the FE

@param iNode is the number 0..n-1 of the FE node

@return a pair of weights, one for each tangent vector

@section implementation Implementation

Accesses the private data of the class, retreiving tabulated vector.

*/

template<size_t dim>
std::pair<double64,double64>
FiniteVolumeStencil<dim>::FacetNormalTransformationNodeWeights( size_t iFacet, size_t iNode ) const
{
    assert( iFacet < facet_normal_xforms.size() );
    assert( iNode < facet_normal_xforms[iFacet].size() );
    return facet_normal_xforms[iFacet][iNode];
}


template<size_t dim>
const Point<dim>& FiniteVolumeStencil<dim>::FacetEdgeMidPoint(size_t iFacet ) const
{
    assert( iFacet < edges_of_element.size() );

	return facet_edge_midpoints[iFacet];
}

template<size_t dim>
const Point<dim>& FiniteVolumeStencil<dim>::Barycenter() const
{
	 return barycenter;
}


template<size_t dim>
const Point<dim>&  FiniteVolumeStencil<dim>::FacetPoint( size_t iFacet, size_t iPoint ) const
{
	 return facet_points[iFacet][iPoint];
}


template<size_t dim>
FV_FACET_TYPE  FiniteVolumeStencil<dim>::FacetType( size_t iFacet ) const
{
     assert( iFacet < facet_types.size() );

	 return facet_types[iFacet];
}

    
    
template<size_t dim>
const csmp::Point<dim>&  FiniteVolumeStencil<dim>::SectorPoint( size_t iSector, size_t iPoint ) const
{
	 return sector_points_[iSector][iPoint];
}
  
template<size_t dim>
const std::pair<csmp::Point<dim>,csmp::Point<dim> >&  FiniteVolumeStencil<dim>::SectorEdgePoints( size_t iSector, size_t iEdge ) const
{
   const size_t pt1(sector_edges_[iSector][iEdge].first);
   const size_t pt2(sector_edges_[iSector][iEdge].second);
   
	 return std::move(std::make_pair( sector_points_[iSector][pt1], sector_points_[iSector][pt2] ));
}       

template<size_t dim>
const std::pair<size_t,size_t>&  FiniteVolumeStencil<dim>::SectorEdge( size_t iSector, size_t iEdge ) const
{
	 return sector_edges_[iSector][iEdge];
}       
    
template<size_t dim>
size_t FiniteVolumeStencil<dim>::SectorPoints( size_t iSector ) const 
{ return sector_points_[iSector].size(); } 

template<size_t dim>
size_t FiniteVolumeStencil<dim>::SectorEdges( size_t iSector ) const 
{ return sector_edges_[iSector].size(); } 
    
  




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
template<size_t dim>
void FiniteVolumeStencil<dim>::Initialize( const char* csp_finite_element_type )
 {
   // uses the types specified in the FiniteElement.h header file
   const string csp_fem_type(csp_finite_element_type);
   bool debug=false;

   size_t NumOfInternalFacets=6U;
   size_t NumOfInternalFacetsPerNode=3U;
   size_t NumOfInternalVolumes=4U;
   size_t NumOfIPperVolume=1U;
   size_t NumOfIPperFacet=1U;
   
   
   if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_BAR) ) {

          NumOfInternalFacets=1U;
          NumOfInternalFacetsPerNode=1U;
          NumOfInternalVolumes=2U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );

          FV_IntegrationPointsAndWeights<dim> FV(ISOPARAMETRIC_LINEAR_BAR);
             
          FV.SectorIntegrationWeights(sector_integration_weights); 

          FV.SectorIntegrationPoints(sector_integration_points);
             
          FV.FacetIntegrationWeights(facet_integration_weights); 
             
          FV.FacetIntegrationPoints(facet_integration_points);
             
          FV.FacetNormals(facet_parametric_normals);

          FV.FacetNormalTransformations(facet_normal_xforms);

          FV.FacetsSurroundingNode(facets_surrounding_node);
             
          FV.EdgePairs(edges_of_element);

		      FV.EdgeMidpoints(facet_edge_midpoints);
		  
		      FV.Barycenter(barycenter);
    
    	    FV.FacetPoints(facet_points);
    	      	  
    	    FV.FacetTypes(facet_types);
    	      	  
    	    FV.SectorPoints(sector_points_);
    	    
     	    FV.SectorEdgePairs(sector_edges_);
     
          space_dimension_ = LINE;

          if(debug) {
             cout<<"\nISOPARAMETRIC_LINEAR_BAR:" <<endl;
             Out();         
           }
           
          parent_element_ = "ISOPARAMETRIC_LINEAR_BAR";
          return;
      }
   if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_TRIANGLE) ) {

          NumOfInternalFacets=3U;
          NumOfInternalFacetsPerNode=2U;
          NumOfInternalVolumes=3U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );

             
          FV_IntegrationPointsAndWeights<dim> FV(ISOPARAMETRIC_LINEAR_TRIANGLE);
             
          FV.SectorIntegrationWeights(sector_integration_weights); 

          FV.SectorIntegrationPoints(sector_integration_points);
             
          FV.FacetIntegrationWeights(facet_integration_weights); 
             
          FV.FacetIntegrationPoints(facet_integration_points);
             
          FV.FacetNormals(facet_parametric_normals);

          FV.FacetNormalTransformations(facet_normal_xforms);

          FV.FacetsSurroundingNode(facets_surrounding_node);
             
          FV.EdgePairs(edges_of_element);

          FV.EdgeMidpoints(facet_edge_midpoints);
		  
          FV.Barycenter(barycenter);

          FV.FacetPoints(facet_points);
       
          FV.FacetTypes(facet_types);

    	    FV.SectorPoints(sector_points_);
    	    
     	    FV.SectorEdgePairs(sector_edges_); 
    	  
          space_dimension_ = SURFACE;
     
          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_TRIANGLE:" <<endl;
              Out();         
          }
           
          parent_element_ = "ISOPARAMETRIC_LINEAR_TRIANGLE";
          return;
      }
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_TETRAHEDRON) ) {

         NumOfInternalFacets=6U;
         NumOfInternalFacetsPerNode=3U;
         NumOfInternalVolumes=4U;

         Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );

         FV_IntegrationPointsAndWeights<dim> FVT(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
         
             FVT.SectorIntegrationWeights(sector_integration_weights); 

             FVT.SectorIntegrationPoints(sector_integration_points);
             
             FVT.FacetIntegrationWeights(facet_integration_weights); 
             
             FVT.FacetIntegrationPoints(facet_integration_points);
             
             FVT.FacetNormals(facet_parametric_normals);

             FVT.FacetNormalTransformations(facet_normal_xforms);

             FVT.FacetsSurroundingNode(facets_surrounding_node);
             
             FVT.EdgePairs(edges_of_element);
             
             FVT.ProjectionWeights(facet_projection_weights);
             
     		     FVT.EdgeMidpoints(facet_edge_midpoints);
		  
             FVT.Barycenter(barycenter);
             
             FVT.FacetPoints(facet_points);
             
             FVT.FacetTypes(facet_types);
             
    	       FVT.SectorPoints(sector_points_);
    	    
        	   FVT.SectorEdgePairs(sector_edges_); 
             
             space_dimension_ = VOLUME;

             if(debug) {
                cout<<"\nISOPARAMETRIC_LINEAR_TETRAHEDRON: "<<endl;
                Out();         
             }
       parent_element_ = "ISOPARAMETRIC_LINEAR_TETRAHEDRON";
       return;
     }
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_QUADRILATERAL) ) {

         NumOfInternalFacets=4U;
         NumOfInternalFacetsPerNode=2U;
         NumOfInternalVolumes=4U; // == VolMultipliers

         Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVQ(ISOPARAMETRIC_LINEAR_QUADRILATERAL);
          
          FVQ.SectorIntegrationWeights(sector_integration_weights); 

          FVQ.SectorIntegrationPoints(sector_integration_points);
         
          FVQ.FacetIntegrationWeights(facet_integration_weights); 
             
          FVQ.FacetIntegrationPoints(facet_integration_points);
             
          FVQ.FacetNormals(facet_parametric_normals);

          FVQ.FacetNormalTransformations(facet_normal_xforms);

          FVQ.FacetsSurroundingNode(facets_surrounding_node);
            
          FVQ.EdgePairs(edges_of_element);

		      FVQ.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVQ.Barycenter(barycenter);

          FVQ.FacetPoints(facet_points);

          FVQ.FacetTypes(facet_types);

    	    FVQ.SectorPoints(sector_points_);
    	    
     	    FVQ.SectorEdgePairs(sector_edges_); 

          space_dimension_ = SURFACE;

          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_QUADRILATERAL:"<<endl;
              Out();         
           }

         parent_element_ = "ISOPARAMETRIC_LINEAR_QUADRILATERAL";
         return;
        
      }
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_HEXAHEDRON) ) {

          NumOfInternalFacets=12U;
          NumOfInternalFacetsPerNode=3U;
          NumOfInternalVolumes=8U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVH(ISOPARAMETRIC_LINEAR_HEXAHEDRON);
         
          FVH.SectorIntegrationWeights(sector_integration_weights); 

          FVH.SectorIntegrationPoints(sector_integration_points);
             
          FVH.FacetIntegrationWeights(facet_integration_weights); 
             
          FVH.FacetIntegrationPoints(facet_integration_points);
             
          FVH.FacetNormals(facet_parametric_normals);

          FVH.FacetNormalTransformations(facet_normal_xforms);
             
          FVH.FacetsSurroundingNode(facets_surrounding_node);

          FVH.EdgePairs(edges_of_element);
          
          FVH.ProjectionWeights(facet_projection_weights);

		      FVH.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVH.Barycenter(barycenter);
          
          FVH.FacetPoints(facet_points);

          FVH.FacetTypes(facet_types);

    	    FVH.SectorPoints(sector_points_);
    	    
     	    FVH.SectorEdgePairs(sector_edges_); 

          space_dimension_ = VOLUME;

         if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_HEXAHEDRON:"<<endl;
              Out();         
          }

      parent_element_ = "ISOPARAMETRIC_LINEAR_HEXAHEDRON";
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
         
          FVPy.SectorIntegrationWeights(sector_integration_weights); 

          FVPy.SectorIntegrationPoints(sector_integration_points);
             
          FVPy.FacetIntegrationWeights(facet_integration_weights); 
             
          FVPy.FacetIntegrationPoints(facet_integration_points);
             
          FVPy.FacetNormals(facet_parametric_normals);

          FVPy.FacetNormalTransformations(facet_normal_xforms);
             
          FVPy.FacetsSurroundingNode(facets_surrounding_node);
          
          FVPy.EdgePairs(edges_of_element);
          
          FVPy.ProjectionWeights(facet_projection_weights);

		      FVPy.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVPy.Barycenter(barycenter);
          
          FVPy.FacetPoints(facet_points);

          FVPy.FacetTypes(facet_types);

    	    FVPy.SectorPoints(sector_points_);
    	    
     	    FVPy.SectorEdgePairs(sector_edges_); 

          space_dimension_ = VOLUME;

         if ( debug ) {
              cout<<"\nISOPARAMETRIC_LINEAR_PYRAMID:"<<endl;
              Out();         
        }

      parent_element_ = "ISOPARAMETRIC_LINEAR_PYRAMID";
      return;     
    }  
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_PRISM) ) {

          NumOfInternalFacets=9U;
          NumOfInternalFacetsPerNode=3U;
          NumOfInternalVolumes=6U;

          Resize(NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                 NumOfInternalVolumes, NumOfIPperVolume, NumOfIPperFacet );
                 
          FV_IntegrationPointsAndWeights<dim> FVP(ISOPARAMETRIC_LINEAR_PRISM);
         
          FVP.SectorIntegrationWeights(sector_integration_weights); 

          FVP.SectorIntegrationPoints(sector_integration_points);
             
          FVP.FacetIntegrationWeights(facet_integration_weights); 
             
          FVP.FacetIntegrationPoints(facet_integration_points);
             
          FVP.FacetNormals(facet_parametric_normals);

          FVP.FacetNormalTransformations(facet_normal_xforms);
             
          FVP.FacetsSurroundingNode(facets_surrounding_node);
          
          FVP.EdgePairs(edges_of_element);
          
          FVP.ProjectionWeights(facet_projection_weights);

		      FVP.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVP.Barycenter(barycenter);

          FVP.FacetPoints(facet_points);
		  
          FVP.FacetTypes(facet_types);
		  
          space_dimension_ = VOLUME;

         if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_PRISM:"<<endl;
              Out();         
          }

      parent_element_ = "ISOPARAMETRIC_LINEAR_PRISM";
      return; 
      }
    else {
         cerr <<"\nFiniteVolumeStencil::Initialize: Element type not recognized: ";
         cerr << csp_finite_element_type << endl;
      }
      
 } // end Initialize






// SKM modified version
template<size_t dim>
void FiniteVolumeStencil<dim>::Out() const
{
   cout<<"\nFiniteVolumeStencil<dim>::Out: \nVolume Weights: ";
   for(size_t i=0; i<sector_integration_weights.size(); i++)
      for(size_t j=0; j<sector_integration_weights[i].size(); j++) 
        cout<<" "<<sector_integration_weights[i][j];
   
   cout<<endl<<endl<<" Volume IPs: "<<endl;
   for(size_t i=0; i<sector_integration_points.size(); i++) {
     for(size_t j=0; j<sector_integration_points[i].size(); j++) {
        for(size_t k=0; k<dim; k++) cout<<" "<<sector_integration_points[i][j][k];
          cout<<endl;
        }
   }
                
   cout<<endl<<" Facet Weights: ";
   for(size_t i=0; i<facet_integration_weights.size(); i++) {
     for(size_t j=0; j<facet_integration_weights[i].size(); j++) {
       cout<<" "<<facet_integration_weights[i][j] ;
      }
   }
   cout<<endl;
   
   cout<<endl<<" Facet IPs: "<<endl;
   for(size_t i=0; i<facet_integration_points.size(); i++) {
     for(size_t j=0; j<facet_integration_points[i].size(); j++) {
       for(size_t k=0; k<dim; k++) cout<<" "<<facet_integration_points[i][j][k];
       cout<<endl;
     }
   }

  cout<<endl<<" Facet parametric normals inside the element: "<<endl;
   for(size_t i=0; i<facet_parametric_normals.size(); i++) {
     for(size_t j=0; j<dim; j++) {
       cout<<" "<<facet_parametric_normals[i][j] ;
     }
     cout<<endl;
   }
     
   cout<<endl<<" Facet physical normals inside the element: "<<endl;
   for(size_t i=0; i<facet_normals.size(); i++) {
     for(size_t j=0; j<dim; j++) {
       cout<<" "<<facet_normals[i][j] ;
     }
     cout<<endl;
   }

    cout<<endl<<" Facet surrounding node in the element: "<<endl;
    for(size_t i=0; i<Sectors(); i++) {
      for(size_t j=0; j<facets_surrounding_node[i].size(); j++) {
        cout<<" "<<facets_surrounding_node[i][j] ;
      }
      cout<<endl;
    }
    
    cout<<endl<<" Edges pairs in the element "<<endl;
    for(size_t i=0; i<edges_of_element.size(); i++) {
      cout<<" "<< edges_of_element[i].first <<", "<< edges_of_element[i].second <<" ";
    }
  cout << endl << endl; 
   
} // end Out


template class FiniteVolumeStencil<1U>; 
template class FiniteVolumeStencil<2U>; 
template class FiniteVolumeStencil<3U>; 


} // end namespace csmp
