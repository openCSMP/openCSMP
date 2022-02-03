#include "FV_IntegrationPointsAndWeights.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<size_t dim>
FV_IntegrationPointsAndWeights<dim>::FV_IntegrationPointsAndWeights()
 {
 } // end constructor



/**

Public constructor of the FV_IntegrationPointsAndWeights class. 
As default constructor of the class is declares private, this constructor
should be used in all cases.

@param typeOfElement type of Finite Element in CSMP
*/
template<size_t dim>
FV_IntegrationPointsAndWeights<dim>::FV_IntegrationPointsAndWeights( CSMP_FEM_TYPE typeOfElement )
 {
    ErrorHandler& error_handler ( ErrorHandler::Instance() );
    if (error_handler.Verbose()){
        cout <<"\nFV_IntegrationPointsAndWeights (constructor): Constructing stencil for ";
        cout << parseFiniteElementType(typeOfElement) <<endl;
    }
    m_typeOfSubdividedElement=typeOfElement;
       
    switch( typeOfElement ) {
         case ISOPARAMETRIC_LINEAR_BAR:
             CreateDataFor_ISOPARAMETRIC_LINEAR_BAR(); 
         break;
    
         case ISOPARAMETRIC_LINEAR_TRIANGLE:
             CreateDataFor_ISOPARAMETRIC_LINEAR_TRI();
         break;

         case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
             CreateDataFor_ISOPARAMETRIC_LINEAR_QUAD();
         break;
         
         case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
              CreateDataFor_ISOPARAMETRIC_LINEAR_TET();
         break;  
         
         case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
              CreateDataFor_ISOPARAMETRIC_LINEAR_HEX();
         break;
         
         case ISOPARAMETRIC_LINEAR_PRISM:
              CreateDataFor_ISOPARAMETRIC_LINEAR_PRISM();
         break;
         
         case ISOPARAMETRIC_LINEAR_PYRAMID:
              CreateDataFor_ISOPARAMETRIC_LINEAR_PYR();
         break;  
          
         default:
         cout <<"\nFV_IntegrationPointsAndWeights (constructor): FiniteElementType not recognized: ";
         cout << parseFiniteElementType(typeOfElement) << endl;
      }
      
 } // end constructor




template<size_t dim>
FV_IntegrationPointsAndWeights<dim>::~FV_IntegrationPointsAndWeights()
{
}





/**

 Method resizes private data arrays of the tabulated data of the  
 class  

@param n_isrf number of internal facets
@param srfs_per_node number of facets, confronting the node
@param n_ivol number of internal sectors
@param n_spts number of facet integration points
@param n_vpts number of volume integration points

@section application Application

The method is used in FiniteVolumeStencil class to resize the internal arrays 
of the FV_IntegrationPointsAndWeights class.
*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::Resize( size_t n_isrf,
                                                  size_t srfs_per_node, 
                                                  size_t n_ivol, 
                                                  size_t n_spts, 
                                                  size_t n_vpts ) 
 {
    m_facet_integration_points.resize( n_isrf );   vector<vector<Point<dim> > >( m_facet_integration_points ).swap( m_facet_integration_points );
    m_facet_integration_weights.resize( n_isrf );  vector<vector<double> >( m_facet_integration_weights ).swap( m_facet_integration_weights );

    m_projection_weights.resize( n_isrf );         vector<vector<double> >(m_projection_weights).swap(m_projection_weights);

    m_facet_normals.resize( n_isrf );              vector<Point<dim> >(m_facet_normals).swap(m_facet_normals);
    m_facets_surrounding_node.resize( n_ivol );    vector<vector<size_t> >(m_facets_surrounding_node).swap(m_facets_surrounding_node);
    
    for ( size_t i=0U; i<n_ivol; i++ ) {
         m_facets_surrounding_node[i].resize( srfs_per_node );
         vector<size_t>(m_facets_surrounding_node[i]).swap(m_facets_surrounding_node[i]);
      }
    
    for ( size_t i=0U; i<n_isrf; i++ ) {
         m_facet_integration_points[i].resize( n_spts );
         vector<Point<dim> >(m_facet_integration_points[i]).swap(m_facet_integration_points[i]);
         m_facet_integration_weights[i].resize(n_spts);
         vector<double>(m_facet_integration_weights[i]).swap(m_facet_integration_weights[i]);
         m_projection_weights[i].resize(n_spts);
         vector<double>(m_projection_weights[i]).swap(m_projection_weights[i]);
      } 
    
    m_volume_integration_points1.resize( n_ivol );     
    vector<vector<Point<dim> > >(m_volume_integration_points1).swap(m_volume_integration_points1);
    m_volume_integration_weights.resize( n_ivol );     
    vector<vector<double> >(m_volume_integration_weights).swap(m_volume_integration_weights);
    m_par_volume_integration_points.resize( n_ivol );  
    vector<vector<Point<dim> > >(m_par_volume_integration_points).swap(m_par_volume_integration_points);

    for ( size_t i=0U; i<n_ivol; i++ ) {
         m_volume_integration_points1[i].resize( n_vpts );
         vector<Point<dim> >(m_volume_integration_points1[i]).swap(m_volume_integration_points1[i]);
         m_par_volume_integration_points[i].resize( n_vpts );
         vector<Point<dim> >(m_par_volume_integration_points[i]).swap(m_par_volume_integration_points[i]);
         m_volume_integration_weights[i].resize(n_vpts);
         vector<double>(m_volume_integration_weights[i]).swap(m_volume_integration_weights[i]);
      } 
    
//    m_facet_edge_midpoints.resize( n_isrf );  vector<Point<dim> >(m_facet_edge_midpoints).swap(m_facet_edge_midpoints);
//    m_facet_points.resize( n_isrf );  vector<vector<Point<dim> > >(m_facet_points).swap(m_facet_points);

    // SKM_FIX for each facet there is one edge with 2 subedges
    m_facet_edge_midpoints.resize( n_isrf*2 );  vector<Point<dim> >(m_facet_edge_midpoints).swap(m_facet_edge_midpoints);

    m_facet_points.resize( n_isrf );  vector<vector<Point<dim> > >(m_facet_points).swap(m_facet_points);
    m_facet_types.resize( n_isrf );
    m_sector_points.resize( n_ivol ); 
    m_edge_of_sectors.resize( n_ivol ); 
     
} // end Resize





template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::SectorIntegrationPoints( 
                          std::vector<std::vector<Point<dim> > >&  volume_integration_points ) const
{
  volume_integration_points = m_volume_integration_points1;
}


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::SectorIntegrationWeights( 
                                       std::vector<std::vector<double> >&  volume_integration_weights ) const
{
  volume_integration_weights = m_volume_integration_weights;
}


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetIntegrationPoints( 
                         std::vector<std::vector<Point<dim> > >&  facet_integration_points ) const
{
   facet_integration_points = m_facet_integration_points;
}


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetIntegrationWeights( 
                                      std::vector<std::vector<double> >&  facet_integration_weights ) const
{
   facet_integration_weights = m_facet_integration_weights;
}


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::ProjectionWeights( 
                                                std::vector<std::vector<double> >&  projection_weights ) const
{
   projection_weights = m_projection_weights;
}

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetNormalTransformations( 
                                                   std::vector<std::vector<std::pair<double,double>>>&  facet_normal_xforms ) const
{
   facet_normal_xforms = m_facet_normal_xforms;
}

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetNormals( 
                                                   std::vector<Point<dim> >&  facet_normals ) const
{
   facet_normals = m_facet_normals;
}


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetsSurroundingNode( 
                                   std::vector<std::vector<size_t> >& facets_surrounding_node ) const
{
   facets_surrounding_node = m_facets_surrounding_node;
}


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::EdgePairs( 
                                    std::vector<std::pair<size_t,size_t> >&  edges_of_element ) const
{
   edges_of_element = m_edges_of_element;
}

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::EdgeMidpoints(std::vector<Point<dim> >& rst_facet_edge_midpoints)
{
   rst_facet_edge_midpoints = m_facet_edge_midpoints;
}
		  
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::Barycenter( Point<dim>& rst_barycenter )
{
	rst_barycenter = m_barycenter;
}

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetPoints( std::vector<std::vector<Point<dim> > >& rst_facet_points )
{
   rst_facet_points = m_facet_points;
}

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::FacetTypes( std::vector<FV_FACET_TYPE>& facet_types )
{
   facet_types = m_facet_types;
}

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::SectorPoints(
                          std::vector<std::vector<Point<dim> > >& rst_sector_points)
{
   rst_sector_points = m_sector_points;
}    


template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::SectorEdgePairs( 
               std::vector<std::vector<std::pair<size_t,size_t> > >&  edges_of_sector ) const
{
   edges_of_sector = m_edge_of_sectors;
}
    


/**

Generates data needed in the subdivision of the isoparametric linear
bar element into CV sectors.
*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_BAR()
{
  const size_t NumOfIPperVolume(1U);
  const size_t NumberOfInternalFacets(1U);
  const size_t NumberOfInternalFacetsPerNode(1U);
  const size_t NumberOfInternalVolumes(2U);
  const size_t NumberOfIntegrationPointsPerFacet(1U);

  Resize( NumberOfInternalFacets, NumberOfInternalFacetsPerNode,
          NumberOfInternalVolumes, NumberOfIntegrationPointsPerFacet, NumOfIPperVolume ); 
  
  m_volume_integration_weights[0][0]  = 1.;
  m_volume_integration_weights[1][0]  = 1.;
  
  m_facet_integration_weights[0][0]   = 1.;
  m_facet_integration_points[0][0][0] = 0.;
  m_facet_normals[0][0]               = 1.; 
   
  // Subdivided line volume integration points (x,y,z, but dependent on element dimension, only some will be non zero)
  const double vip[2U][3U] = {
                      {-0.5, 0., 0. }, 
                      { 0.5, 0., 0. } 
                   }; 
    
  // in parametric space, the bar element has just a single dimension r
  for( size_t i=0; i<NumberOfInternalVolumes; i++ )
    m_volume_integration_points1[i][0][0U]=vip[i][0U];

  // Facet surrounding node
  m_facets_surrounding_node[0][0] = 0U;
  m_facets_surrounding_node[1][0] = 0U;

  // Edges
  m_edges_of_element.resize(1U);
  vector<pair<size_t,size_t> >(m_edges_of_element).swap(m_edges_of_element);
  m_edges_of_element[0]=make_pair(0U,1U);

  m_facet_edge_midpoints[0][0] = -1./2.;
  m_facet_edge_midpoints[1][0] =  1./2.;

  m_barycenter[0] = 0.;

  m_facet_types[0] = POINT_FACET;

  m_facet_points[0].resize(1U); //nr of points in this facet
  vector<Point<dim> >(m_facet_points[0]).swap(m_facet_points[0]);
  m_facet_points[0][0] = m_barycenter; 

   m_facet_normal_xforms.resize(1u);
   m_facet_normal_xforms[0].resize(2u);
   m_facet_normal_xforms[0][0] = std::make_pair(-0.5, -0.5);
   m_facet_normal_xforms[0][1] = std::make_pair(0.5, 0.5);

  m_sector_points[0].resize(2U);
  m_sector_points[0][0] = -1.; 
  m_sector_points[0][1] = m_barycenter; 
  
  m_sector_points[1].resize(2U); 
  m_sector_points[1][0] = m_barycenter; 
  m_sector_points[1][1] = 1.; 

  m_edge_of_sectors.resize(2U);
  for ( size_t jSector=0U; jSector < 2U; jSector++ ) {
      m_edge_of_sectors[jSector].resize(1U);
      // in this case each sector only has one edge
      m_edge_of_sectors[jSector][0] = make_pair(0U,1U); 
   }
}





/**

Generates data for subdivision of isoparametric linear
triangle to CV sectors.

*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_TRI()
{
  // Total number of line segments in the triangle is 9, but internal - 3, by one IP 
  const size_t NumOfIPperVolume(1U);
  Resize(3,2,3,1,NumOfIPperVolume); 

  // Volume (area) of each sector 
  for(size_t i=0; i<3U; i++) //where 3 is the number of sectors
    for(size_t j=0;j<NumOfIPperVolume;j++)
       // Hard coded, check multiple points!
       m_volume_integration_weights[i][j]=1./6.;
     
  // facet integration points
  const double fip[3U][3U]=
     {
         { 5./12., 1./6.,  0. }, 
         { 5./12., 5./12., 0. }, // point adjacent to hypotenuse
         { 1./6.,  5./12., 0. }, 
     }; 

   // products of detJ's, weights & conversion factors
   const size_t iNumberOfFacets(3U);
  
   m_facet_integration_weights[0][0]=sqrt(5.)/6.;
   m_facet_integration_weights[1][0]=sqrt(2.)/6.;
   m_facet_integration_weights[2][0]=sqrt(5.)/6.;

   // fill in vectors of facet integration points
   for(size_t i=0;i<iNumberOfFacets;i++)
    for (size_t j=0;j<dim; j++)
      m_facet_integration_points[i][0][j]=fip[i][j];  
  
   // facet unit normals  
   const double cosAlpha1(static_cast<double>(2./sqrt(5.)));
   const double sinAlpha1(static_cast<double>(1./sqrt(5.)));
   const double cosAlpha2(static_cast<double>(sqrt(2.)/2.));
   const double sinAlpha2(cosAlpha2);

   m_facet_normals[0][0]= cosAlpha1; m_facet_normals[0][1]= sinAlpha1;
   m_facet_normals[1][0]=-cosAlpha2; m_facet_normals[1][1]= sinAlpha2;
   m_facet_normals[2][0]=-sinAlpha1; m_facet_normals[2][1]=-cosAlpha1; 
   
   // triangle sector integration points
   const double vip[3U][3U] =
    {
      { 7./36.,  7./36., 0. },
      { 11./18., 7./36., 0. },
      { 7./36., 11./18., 0. }
    }; 
    
    for(size_t i=0U;i<3U;i++)
      for(size_t j=0U;j<dim; j++) m_volume_integration_points1[i][0][j]=vip[i][j];
      
    // facets surrounding node
    m_facets_surrounding_node[0][0]=0U; m_facets_surrounding_node[0][1]=2U;     
    m_facets_surrounding_node[1][0]=1U; m_facets_surrounding_node[1][1]=0U;
    m_facets_surrounding_node[2][0]=2U; m_facets_surrounding_node[2][1]=1U;       
    
    // Edges
    m_edges_of_element.resize(3U);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,0U);
    
    //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 1./2.; m_facet_edge_midpoints[0][1] = 0.;
    m_facet_edge_midpoints[1][0] = 1./2.; m_facet_edge_midpoints[1][1] = 1./2.;
    m_facet_edge_midpoints[2][0] = 0.;    m_facet_edge_midpoints[2][1] = 1./2.;
    
    //m_barycenter
    m_barycenter[0] = 1./3.;
    m_barycenter[1] = 1./3.;
    
    // tested: 25/Oct/2017 AJB
    double facet_normal_transform_scale[3] = {
        4.0 * sqrt(5.0) / 3.0,
        8.0 / (3.0 * sqrt(2.0)),
        4.0 * sqrt(5.0) / 3.0
    };
    const double facet_normal_transforms[3][2][3] = {
        { { -1.0/6.0,-1.0/6.0,1.0/3.0 }, { -5.0/12.0,1.0/3.0,1.0/12.0 } },
        { { 1.0/3.0,-1.0/6.0,-1.0/6.0 }, { 1.0/12.0,-5.0/12.0,1.0/3.0 } },
        { { -1.0/6.0,1.0/3.0,-1.0/6.0 }, { 1.0/3.0,1.0/12.0,-5.0/12.0 } }
    };
    m_facet_normal_xforms.resize(3u);
    for ( size_t iFacet=0U; iFacet<3U; iFacet++ ) {
        m_facet_types[iFacet] = UNIT_LINEAR_FACET;

        m_facet_points[iFacet].resize(2U);
        vector<Point<dim> >(m_facet_points[iFacet]).swap(m_facet_points[iFacet]);
        m_facet_points[iFacet][0] = m_facet_edge_midpoints[iFacet];
        m_facet_points[iFacet][1] = m_barycenter;

        m_facet_normal_xforms[iFacet].resize(3u);
        for (size_t j = 0; j < 3; ++j) {
            m_facet_normal_xforms[iFacet][j]
                = std::make_pair(
                        facet_normal_transforms[iFacet][0][j],
                        facet_normal_transforms[iFacet][1][j] * facet_normal_transform_scale[iFacet]);
        }
      }

    
    //tested: 30/Oct/2007 Hamid
    Point<dim> p;
    // sector points numbered counter-clockwise from the outside looking in
    // and starting with the sector (node) point
    for ( size_t iSector=0U; iSector<3U; iSector++ ) { 
         // each sector has four points 
         m_sector_points[iSector].resize(4U); 
         // generating sector node coordinates assuming that triangle has a 90o
         // angle at first node         
         p[0] =  iSector % 2; 
         p[1] = (iSector/2)*((iSector+1) % 2);
         m_sector_points[iSector][0] = p;
         m_sector_points[iSector][1] = m_facet_edge_midpoints[iSector];        
         m_sector_points[iSector][2] = m_barycenter;                           
         m_sector_points[iSector][3] = m_facet_edge_midpoints[(2+iSector) % 3];
      }
    
    // for all sectors  
    m_edge_of_sectors.resize(3U);
    for ( size_t jSector = 0U; jSector < 3U; jSector++ ) {
         m_edge_of_sectors[jSector].resize(4U);
         for ( size_t jEdge = 0U; jEdge < 4U; jEdge++ )
           m_edge_of_sectors[jSector][jEdge] = make_pair(jEdge,(jEdge+1) % 4); 
      }
    
} // end CreateDataFor_ISOPARAMETRIC_LINEAR_TRIANGLE





/**

Generates data for subdivision of isoparametric linear
quadrilateral to the CV sectors.
*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_QUAD()
{
  const size_t NumOfInternalFacets(4U);
  const size_t NumOfInternalFacetsPerNode(2U);
  const size_t NumOfInternalVolumes(4U); // == VolMultipliers
  const size_t NumOfIPperVolume(1U);
  const size_t NumOfIPperFacet(1U);

  //Resize SKM data structure  int_surf, srf_per_node, int_vols, int_pts_per_surf, 
   Resize( NumOfInternalFacets,
           NumOfInternalFacetsPerNode,
           NumOfInternalVolumes,
           NumOfIPperFacet,
           NumOfIPperVolume );
 
  // Map volume multipliers to new vectors
  for(size_t i=0;i<NumOfInternalVolumes;i++)
    for(size_t j=0;j<NumOfIPperVolume;j++)
      m_volume_integration_weights[i][j]=1.;
     
  // facet integration points
  const double fip[4U][3U]=
     {
         { 0.0,-0.5, 0. },
         { 0.5, 0.0, 0. },
         { 0.0, 0.5, 0. },
         {-0.5, 0.0, 0. }
     }; 

   // products of detJ's, weights & conversion factors
   
   // Fill in vectors of face variables
   for(size_t i=0;i<NumOfInternalFacets; i++ ){
     m_facet_integration_weights[i][0]=1.0;
    	for(size_t j=0;j<dim; j++) {
           m_facet_integration_points[i][0][j]=fip[i][j];  
        }
    }
 
   // Subdivided Quad Volume integration points
   const double vip[4U][3U] =
    {
      {-0.5,-0.5, 0. },
      {0.5 ,-0.5, 0. },
      {0.5 , 0.5, 0. },
      {-0.5, 0.5, 0. }
    }; 
    
   for(size_t i=0U; i<NumOfInternalVolumes; i++)
      for(size_t j=0U; j<dim; j++ ) m_volume_integration_points1[i][0][j]=vip[i][j];

   m_facet_normals[0][0]= 1.0; m_facet_normals[0][1]= 0.0;
   m_facet_normals[1][0]= 0.0; m_facet_normals[1][1]= 1.0;
   m_facet_normals[2][0]=-1.0; m_facet_normals[2][1]= 0.0;
   m_facet_normals[3][0]= 0.0; m_facet_normals[3][1]=-1.0;
     
   // facets surrounding node
   m_facets_surrounding_node[0][0]=0U; m_facets_surrounding_node[0][1]=3U;     
   m_facets_surrounding_node[1][0]=1U; m_facets_surrounding_node[1][1]=0U;
   m_facets_surrounding_node[2][0]=2U; m_facets_surrounding_node[2][1]=1U;
   m_facets_surrounding_node[3][0]=3U; m_facets_surrounding_node[3][1]=2U;
   
    // Edges
    m_edges_of_element.resize(4U);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,3U);
    m_edges_of_element[3]=make_pair(3U,0U);
    
    //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 0.; m_facet_edge_midpoints[0][1] = -1.;
    m_facet_edge_midpoints[1][0] = 1.; m_facet_edge_midpoints[1][1] = 0.;
    m_facet_edge_midpoints[2][0] = 0.; m_facet_edge_midpoints[2][1] = 1.;
    m_facet_edge_midpoints[3][0] =-1.; m_facet_edge_midpoints[3][1] = 0.;
    
    //m_barycenter
    m_barycenter[0] = 0.;
    m_barycenter[1] = 0.;

    const double facet_normal_transforms[4][2][4] = {
        { { -1.0/4.0, -1.0/4.0, 1.0/4.0, 1.0/4.0 }, { -3.0/8.0, 3.0/8.0, 1.0/8.0, -1.0/8.0 } },
        { { 1.0/4.0, -1.0/4.0, -1.0/4.0, 1.0/4.0 }, { -1.0/8.0, -3.0/8.0, 3.0/8.0, 1.0/8.0 } },
        { { 1.0/4.0, 1.0/4.0, -1.0/4.0, -1.0/4.0 }, { 1.0/8.0, -1.0/8.0, -3.0/8.0, 3.0/8.0 } },
        { { -1.0/4.0, 1.0/4.0, 1.0/4.0, -1.0/4.0 }, { 3.0/8.0, 1.0/8.0, -1.0/8.0, -3.0/8.0 } }
    };
    m_facet_normal_xforms.resize(4u);

    for ( size_t iFacet = 0U; iFacet < 4; iFacet++ )
       {
         m_facet_types[iFacet] = UNIT_LINEAR_FACET;

         m_facet_normal_xforms[iFacet].resize(4u);
         for (size_t j = 0; j < 4; ++j) {
             m_facet_normal_xforms[iFacet][j]
                 = std::make_pair(
                        facet_normal_transforms[iFacet][0][j],
                        facet_normal_transforms[iFacet][1][j]);
         }

         m_facet_points[iFacet].resize(2U);
         vector<Point<dim> >(m_facet_points[iFacet]).swap(m_facet_points[iFacet]);
         m_facet_points[iFacet][0] = m_facet_edge_midpoints[iFacet];
         m_facet_points[iFacet][1] = m_barycenter;
       } 
    
    // tested: 08/Nov/2007 Hamid	
    Point<dim> p;
    for ( size_t iSector=0U; iSector<4U; iSector++ ) {  
         m_sector_points[iSector].resize(4U);               
         m_sector_points[iSector][1] = m_facet_edge_midpoints[iSector];         
         m_sector_points[iSector][2] = m_barycenter;                             
         m_sector_points[iSector][3] = m_facet_edge_midpoints[(3+iSector) % 4];  
      }
 
    p[0] = -1.; p[1] = -1.;
    m_sector_points[0][0] = p;
    p[0] =  1.; p[1] = -1.;
    m_sector_points[1][0] = p;
    p[0] =  1.; p[1] =  1.;
    m_sector_points[2][0] = p;
    p[0] = -1.; p[1] =  1.;
    m_sector_points[3][0] = p;

    m_edge_of_sectors.resize(4U);
    for ( size_t jSector = 0U; jSector < 4U; jSector++ ) {
        m_edge_of_sectors[jSector].resize(4U);
        for ( size_t jEdge = 0U; jEdge < 4U; jEdge++ )
          m_edge_of_sectors[jSector][jEdge] = make_pair(jEdge,(jEdge+1)%4); 
      }
    
} // end CreateDataFor_ISOPARAMETRIC_LINEAR_QUAD




/**

Generates stencil data for subdivision of unit isoparametric linear
tetrahedron to the CV sectors.
*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_TET()
 {
    const size_t NumOfInternalFacets(6U);
    const size_t NumOfInternalFacetsPerNode(3U);
    const size_t NumOfInternalVolumes(4U); // == VolMultipliers
  	const size_t NumOfIPperVolume(1U);
    const size_t NumOfIPperFacet(1U);

   Resize( NumOfInternalFacets,
            NumOfInternalFacetsPerNode,
            NumOfInternalVolumes,
            NumOfIPperFacet,
            NumOfIPperVolume );
 
    // Map volume multipliers to new vectors
    for(size_t i=0;i<NumOfInternalVolumes;i++)
      for(size_t j=0;j<NumOfIPperVolume;j++)
       // Hard coded, check multiple points!
       m_volume_integration_weights[i][j]=1./24.;
	
    const double fip[6U][3U] =
    {
     {13./36.,5./36.,5./36.},
     {13./36.,13./36.,5./36.},
     {5./36.,13./36.,5./36.},
     {5./36.,5./36.,13./36.},
     {13./36.,5./36.,13./36.},
     {5./36.,13./36.,13./36.}
    };
    
    // Mapping to the new vectors    		
    for(size_t i=0;i<NumOfInternalFacets;i++)
      for(size_t j=0;j<dim; j++)
    	  m_facet_integration_points[i][0][j]=fip[i][j];
    		  		
    // Subdivided Tet Volume integration points
    const double vip[4U][3U] =
    {
     {23./144.,23./144.,23./144.},
     {25./48.,23./144.,23./144.},
     {23./144.,25./48.,23./144.},
     {23./144.,23./144.,25./48.}
    };
    
    // Map volume IP to vectors 		
    for(size_t i=0;i<NumOfInternalVolumes;i++)
     for(size_t j=0;j<NumOfIPperVolume; j++)
      for(size_t k=0;k<dim; k++)
        m_volume_integration_points1[i][j][k]=vip[i][k];		
    
    // products of detJ's, weights & conversion factors
   m_facet_integration_weights[0][0]=m_facet_integration_weights[2][0]=m_facet_integration_weights[3][0]=sqrt(6.)/24;
   m_facet_integration_weights[4][0]=m_facet_integration_weights[5][0]=m_facet_integration_weights[1][0]=sqrt(2.)/24;
  
   //Projection weights of the transformation  tet->hex
   m_projection_weights[0][0]=0.0757606720336959499;
   m_projection_weights[1][0]=0.0537296378927313301;
   m_projection_weights[2][0]=0.0531111233171085023;
   m_projection_weights[3][0]=0.0401725120584219275;
   m_projection_weights[4][0]=0.0406761158970540965;
   m_projection_weights[5][0]=0.0240110680599455104;

    // Facet normals pointing outside the internal face, all together 6 normals!
    // Normal 0
    m_facet_normals[0][0]=sqrt(6.)/3.;m_facet_normals[0][1]=sqrt(6.)/6.;m_facet_normals[0][2]=sqrt(6.)/6.;
    // Normal 1
    m_facet_normals[1][0]=-sqrt(2.)/2.;m_facet_normals[1][1]=sqrt(2.)/2.;m_facet_normals[1][2]=0.;
    // Normal 2
    m_facet_normals[2][0]=-sqrt(6.)/6.;m_facet_normals[2][1]=-sqrt(6.)/3.;m_facet_normals[2][2]=-sqrt(6.)/6.;
    // Normal 3
    m_facet_normals[3][0]=sqrt(6.)/6.;m_facet_normals[3][1]=sqrt(6.)/6.;m_facet_normals[3][2]=sqrt(6.)/3.;
    // Normal 4
    m_facet_normals[4][0]=-sqrt(2.)/2.;m_facet_normals[4][1]=0.;m_facet_normals[4][2]=sqrt(2.)/2.;
    // Normal 5
    m_facet_normals[5][0]=0;m_facet_normals[5][1]=-sqrt(2.)/2;m_facet_normals[5][2]=sqrt(2.)/2;

    // facets surrounding node
    m_facets_surrounding_node[0][0]=0; m_facets_surrounding_node[0][1]=2; m_facets_surrounding_node[0][2]=3;    
    m_facets_surrounding_node[1][0]=1; m_facets_surrounding_node[1][1]=0; m_facets_surrounding_node[1][2]=4;
    m_facets_surrounding_node[2][0]=2; m_facets_surrounding_node[2][1]=1; m_facets_surrounding_node[2][2]=5;
    m_facets_surrounding_node[3][0]=5; m_facets_surrounding_node[3][1]=4; m_facets_surrounding_node[3][2]=3;
    
    // Edges
    m_edges_of_element.resize(6U);
    vector<pair<size_t,size_t> >(m_edges_of_element).swap(m_edges_of_element);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,0U);
    m_edges_of_element[3]=make_pair(0U,3U);
    m_edges_of_element[4]=make_pair(1U,3U);    
    m_edges_of_element[5]=make_pair(2U,3U);            
  
    //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 1./3.; m_facet_edge_midpoints[0][1] = 0.;    m_facet_edge_midpoints[0][2] = 1./3.; 
    m_facet_edge_midpoints[1][0] = 1./3.; m_facet_edge_midpoints[1][1] = 1./3.; m_facet_edge_midpoints[1][2] = 1./3.; 
    m_facet_edge_midpoints[2][0] = 0.;    m_facet_edge_midpoints[2][1] = 1./3.; m_facet_edge_midpoints[2][2] = 1./3.; 
    
    m_facet_edge_midpoints[3] = m_facet_edge_midpoints[2]; 
    m_facet_edge_midpoints[4] = m_facet_edge_midpoints[0];
    m_facet_edge_midpoints[5] = m_facet_edge_midpoints[1];
    
    //m_barycenter
    m_barycenter[0] = 1./4.;
    m_barycenter[1] = 1./4.;
	  m_barycenter[2] = 1./4.;
    
    Point<dim>  pt_c123, pt_c234, pt_c124, pt_c134, pt_c23, pt_c13, pt_c14, pt_c34, pt_c24, pt_c12;
    pt_c123[0]=1./3.;pt_c123[1]=1./3.;pt_c123[2]=0.;
    pt_c234[0]=1./3.;pt_c234[1]=1./3.;pt_c234[2]=1./3.;
    pt_c124[0]=1./3.;pt_c124[1]=0.;pt_c124[2]=1./3.;
    pt_c134[0]=0.;pt_c134[1]=1./3.;pt_c134[2]=1./3.;
    
    pt_c23[0]=1./2.;pt_c23[1]=1./2.;pt_c23[2]=0.;
    pt_c13[0]=0.;pt_c13[1]=1./2.;pt_c13[2]=0.;
    pt_c14[0]=0.;pt_c14[1]=0.;pt_c14[2]=1./2.;
    pt_c34[0]=0.;pt_c34[1]=1./2.;pt_c34[2]=1./2.;
    pt_c24[0]=1./2.;pt_c24[1]=0.;pt_c24[2]=1./2.;
    pt_c12[0]=1./2.;pt_c12[1]=0.;pt_c12[2]=0.;
   
    m_facet_types[0] = QUADRILATERAL_FACET;

    m_facet_points[0].resize(4U);
    vector<Point<dim> >(m_facet_points[0]).swap(m_facet_points[0]);
    m_facet_points[0][0] = pt_c12;
    m_facet_points[0][1] = pt_c123; 
    m_facet_points[0][2] = m_barycenter;
    m_facet_points[0][3] = pt_c124;
    
    m_facet_types[1] = QUADRILATERAL_FACET;

    m_facet_points[1].resize(4U);
    vector<Point<dim> >(m_facet_points[1]).swap(m_facet_points[1]);
    m_facet_points[1][0] = pt_c23;
    m_facet_points[1][1] = pt_c123; 
    m_facet_points[1][2] = m_barycenter;
    m_facet_points[1][3] = pt_c234;
    
    m_facet_types[2] = QUADRILATERAL_FACET;

    m_facet_points[2].resize(4U);
    vector<Point<dim> >(m_facet_points[2]).swap(m_facet_points[2]);
    m_facet_points[2][0] = pt_c13;
    m_facet_points[2][1] = pt_c123; 
    m_facet_points[2][2] = m_barycenter;
    m_facet_points[2][3] = pt_c134;
    
    m_facet_types[3] = QUADRILATERAL_FACET;

    m_facet_points[3].resize(4U);
    vector<Point<dim> >(m_facet_points[3]).swap(m_facet_points[3]);
    m_facet_points[3][0] = pt_c14;
    m_facet_points[3][1] = pt_c124; 
    m_facet_points[3][2] = m_barycenter;
    m_facet_points[3][3] = pt_c134;
 
    m_facet_types[4] = QUADRILATERAL_FACET;

    m_facet_points[4].resize(4U);
    vector<Point<dim> >(m_facet_points[4]).swap(m_facet_points[4]);
    m_facet_points[4][0] = pt_c124;
    m_facet_points[4][1] = pt_c24; 
    m_facet_points[4][2] = pt_c234;
    m_facet_points[4][3] = m_barycenter;
 
    m_facet_types[5] = QUADRILATERAL_FACET;

    m_facet_points[5].resize(4U);
    vector<Point<dim> >(m_facet_points[5]).swap(m_facet_points[5]);
    m_facet_points[5][0] = m_barycenter;
    m_facet_points[5][1] = pt_c234; 
    m_facet_points[5][2] = pt_c34;
    m_facet_points[5][3] = pt_c134;

    // tested: 18/Oct/2017 AJB
    const double facet_normal_transforms[6][2][4] = {
        { { -1.0/8.0,-1.0/8.0,-1.0/24.0,7.0/24.0 }, { -1.0/8.0,-1.0/8.0,7.0/24.0,-1.0/24.0 } },
        { { -1.0/24.0,-1.0/8.0,-1.0/8.0,7.0/24.0 }, { 7.0/24.0,-1.0/8.0,-1.0/8.0,-1.0/24.0 } },
        { { -1.0/8.0,-1.0/24.0,-1.0/8.0,7.0/24.0 }, { -1.0/8.0,7.0/24.0,-1.0/8.0,-1.0/24.0 } },
        { { -1.0/8.0,-1.0/24.0,7.0/24.0,-1.0/8.0 }, { -1.0/8.0,7.0/24.0,-1.0/24.0,-1.0/8.0 } },
        { { -1.0/24.0,-1.0/8.0,7.0/24.0,-1.0/8.0 }, { -7.0/24.0,1.0/8.0,1.0/24.0,1.0/8.0 } },
        { { 1.0/24.0,-7.0/24.0,1.0/8.0,1.0/8.0 }, { -7.0/24.0,1.0/24.0,1.0/8.0,1.0/8.0 } }
    };
    m_facet_normal_xforms.resize(6u);
    for (size_t f = 0; f < 6u; ++f) {
        m_facet_normal_xforms[f].resize(4u);
        for (size_t n = 0; n < 4u; ++n) {
            m_facet_normal_xforms[f][n] = std::make_pair(
                    facet_normal_transforms[f][0][n],
                    facet_normal_transforms[f][1][n]
                );
        }
    }

    // tested: 08/Nov/2007 Hamid	
    Point<dim> p;
    m_sector_points[0].resize(8U); 
    p[0]=0.;p[1]=0.;p[2]=0.;  
    m_sector_points[0][0] = p;
    m_sector_points[0][1] = pt_c12; 
    m_sector_points[0][2] = pt_c123;
    m_sector_points[0][3] = pt_c13; 
    m_sector_points[0][4] = pt_c14; 
    m_sector_points[0][5] = pt_c124;
    m_sector_points[0][6] = m_barycenter;
    m_sector_points[0][7] = pt_c134; 

    m_sector_points[1].resize(8U); 
    p[0]=1.;p[1]=0.;p[2]=0.;  
    m_sector_points[1][0] = p;
    m_sector_points[1][1] = pt_c23; 
    m_sector_points[1][2] = pt_c123;
    m_sector_points[1][3] = pt_c12; 
    m_sector_points[1][4] = pt_c24; 
    m_sector_points[1][5] = pt_c234;
    m_sector_points[1][6] = m_barycenter; 
    m_sector_points[1][7] = pt_c124; 
  
    m_sector_points[2].resize(8U); 
    p[0]=0.;p[1]=1.;p[2]=0.;  
    m_sector_points[2][0] = p;
    m_sector_points[2][1] = pt_c13;
    m_sector_points[2][2] = pt_c123;
    m_sector_points[2][3] = pt_c23; 
    m_sector_points[2][4] = pt_c34; 
    m_sector_points[2][5] = pt_c134;
    m_sector_points[2][6] = m_barycenter;
    m_sector_points[2][7] = pt_c234;
  
    m_sector_points[3].resize(8U); 
    p[0]=0.;p[1]=0.;p[2]=1.;  
    m_sector_points[3][0] = p;
    m_sector_points[3][1] = pt_c24; 
    m_sector_points[3][2] = pt_c234;
    m_sector_points[3][3] = pt_c34; 
    m_sector_points[3][4] = pt_c14; 
    m_sector_points[3][5] = pt_c124;
    m_sector_points[3][6] = m_barycenter;
    m_sector_points[3][7] = pt_c134;

    // for all sectors
    m_edge_of_sectors.resize(4U); 
    for ( size_t iSector=0U; iSector<4U; iSector++ ) 
      { 
         m_edge_of_sectors[iSector].resize(12U);
         for ( size_t iLine=0U; iLine<4U; iLine++ )
           m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine+1)%4); 
         for ( size_t iLine=4U; iLine<8U; iLine++ )
           m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine-3)%4+4);  
         for ( size_t iLine=8U; iLine<12U; iLine++ ) 
           m_edge_of_sectors[iSector][iLine]=make_pair(iLine-8,iLine-4);  
      }  
    
} // CreateDataFor_ISOPARAMETRIC_LINEAR_TET()





/**
 
Generates data for subdivision of unit isoparametric linear
hexahedron to the CV sectors.  
*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_HEX()
 {
    const size_t NumOfInternalFacets(12U);
    const size_t NumOfInternalFacetsPerNode(3U);
    const size_t NumOfInternalVolumes(8U); // == VolMultipliers
  	const size_t NumOfIPperVolume(1U);
    const size_t NumOfIPperFacet(1U);

   Resize( NumOfInternalFacets,
            NumOfInternalFacetsPerNode,
            NumOfInternalVolumes,
            NumOfIPperFacet,
            NumOfIPperVolume );

    // Map volume multipliers to new vectors
    for(size_t i=0;i<NumOfInternalVolumes;i++)
      for(size_t j=0;j<NumOfIPperVolume;j++)
       // Hard coded, check multiple points!
       m_volume_integration_weights[i][j]=1.0;
	
    // Faces 
    // products of detJ's, weights & conversion factors
    //Map to new vectors
    for( size_t i=0;i<NumOfInternalFacets;i++)
      for ( size_t j=0;j<NumOfIPperFacet; j++) {
           m_projection_weights[i][j]=1.;
           m_facet_integration_weights[i][j]=1.;
        }

	// face integration points (3D element)
	const double fip[12U][3U] =
    {
     // Lower
     {0., - 0.5,-0.5},
     {0.5,  0.0, -0.5},    
     {0.0,  0.5, -0.5},
     {-0.5, 0.,-0.5},
     //Mid
     {-0.5,-0.5,0.0},
     { 0.5,-0.5,0.0 },
     { 0.5, 0.5,0.0},
     {-0.5, 0.5,0.0},
     //Up
     {0.0,-0.5,0.5},
     {0.5,0.0,0.5},
     {0.0,0.5,0.5},
     {-0.5,0.0,0.5},
    };
	
     // Map to new vectors
    for(size_t i=0;i<NumOfInternalFacets;i++)
     for(size_t j=0;j<NumOfIPperFacet; j++)
      for(size_t k=0;k<dim; k++) 
         m_facet_integration_points[i][j][k]=fip[i*NumOfIPperFacet+j][k];
    		
    // Subdivided volume integration points for hexahedron
    const double vip[8U][3U]=
        {
            {-0.5,-0.5,-0.5},
            {0.5,-0.5,-0.5}, 
            {0.5,0.5,-0.5},
            {-0.5,0.5,-0.5}, 

            {-0.5,-0.5,0.5},
            {0.5,-0.5,0.5},
            {0.5,0.5,0.5},
            {-0.5,0.5,0.5} 
        }; 
  
    // m_volume_integration_points1;   // [ivol][vpts][dim]
    for( size_t i=0U; i<NumOfInternalVolumes; i++ )
      for( size_t j=0U; j<NumOfIPperVolume; j++ )
    	for( size_t k=0U; k<dim; k++)
    	  m_volume_integration_points1[i][j][k]=vip[i*NumOfIPperVolume+j][k];

    // Facet normals pointing outside the internal face, all together 6 normals!
    // Normal 0, +r facet IP[0] OK
    //start 0., - 0.5,-0.5
    m_facet_normals[0][0]=1.0; m_facet_normals[0][1]=0.0; m_facet_normals[0][2]=0.0; 
    // Normal 1, +s
    // start -0.5, 0., -0.5
    m_facet_normals[1][0]=0.0; m_facet_normals[1][1]=1.0; m_facet_normals[1][2]=0.0;
    // Normal 2 +t
    // start -0.5,-0.5, 0.0
    m_facet_normals[2][0]=-1.0; m_facet_normals[2][1]=0.0; m_facet_normals[2][2]=0.0;
     // Normal 3 +S
    // start 0.5, 0.,-0.5
    m_facet_normals[3][0]=0.0; m_facet_normals[3][1]=-1.0; m_facet_normals[3][2]=0.0;
   
    // Normal 4, +t
    // start 0.5,-0.5,0.
    m_facet_normals[4][0]=0.0; m_facet_normals[4][1]=0.0; m_facet_normals[4][2]=1.0;
    // Normal 5, -r
    // start 0.0,0.5,-0.5
    m_facet_normals[5][0]=0.0; m_facet_normals[5][1]=0.0; m_facet_normals[5][2]=1.0;
    // Normal 6, +t
    // start 0.5,0.5,0.0
    m_facet_normals[6][0]=0.0; m_facet_normals[6][1]=0.0; m_facet_normals[6][2]=1.0;
    // Normal 7, +t    // start -0.5,0.5,0.0
    m_facet_normals[7][0]=0.0; m_facet_normals[7][1]=0.0; m_facet_normals[7][2]=1.0;

    // Normal 8,+r start 0.0,-0.5,0.5
    m_facet_normals[8][0]=1.0; m_facet_normals[8][1]=0.0; m_facet_normals[8][2]=0.0;
    // Normal 9, +s  start -0.5,0.0,0.5
    m_facet_normals[9][0]=0.0; m_facet_normals[9][1]=1.0; m_facet_normals[9][2]=0.0;
    // Normal 10 +s start 0.5,0.0,0.5
    m_facet_normals[10][0]=-1.0; m_facet_normals[10][1]=0.0; m_facet_normals[10][2]=0.0;
    // Normal 11 -r   start 0.0,0.5,0.5
    m_facet_normals[11][0]=0.0;m_facet_normals[11][1]=-1.0; m_facet_normals[11][2]=0.0;

   
    // facets surrounding node
    m_facets_surrounding_node[0][0]=0; m_facets_surrounding_node[0][1]=3; m_facets_surrounding_node[0][2]=4;    
    m_facets_surrounding_node[1][0]=1; m_facets_surrounding_node[1][1]=0; m_facets_surrounding_node[1][2]=5;
    m_facets_surrounding_node[2][0]=2; m_facets_surrounding_node[2][1]=1; m_facets_surrounding_node[2][2]=6;
    m_facets_surrounding_node[3][0]=3; m_facets_surrounding_node[3][1]=2; m_facets_surrounding_node[3][2]=7;
    m_facets_surrounding_node[4][0]=11; m_facets_surrounding_node[4][1]=8;  m_facets_surrounding_node[4][2]=4;    
    m_facets_surrounding_node[5][0]=8;  m_facets_surrounding_node[5][1]=9;  m_facets_surrounding_node[5][2]=5;
    m_facets_surrounding_node[6][0]=9;  m_facets_surrounding_node[6][1]=10; m_facets_surrounding_node[6][2]=6;
    m_facets_surrounding_node[7][0]=10; m_facets_surrounding_node[7][1]=11; m_facets_surrounding_node[7][2]=7;
  

    //Edges   
    m_edges_of_element.resize(12U);
    vector<pair<size_t,size_t> >(m_edges_of_element).swap(m_edges_of_element);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,3U);
    m_edges_of_element[3]=make_pair(3U,0U);
 
    m_edges_of_element[4]=make_pair(0U,4U);    
    m_edges_of_element[5]=make_pair(1U,5U);
    m_edges_of_element[6]=make_pair(2U,6U);
    m_edges_of_element[7]=make_pair(3U,7U);            
 
    m_edges_of_element[8]=make_pair(4U,5U);    
    m_edges_of_element[9]=make_pair(5U,6U);
    m_edges_of_element[10]=make_pair(6U,7U);
    m_edges_of_element[11]=make_pair(7U,4U);            
 
    //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 0.; m_facet_edge_midpoints[0][1] = -1.; m_facet_edge_midpoints[0][2] = 0.; 
    m_facet_edge_midpoints[1][0] = 1.; m_facet_edge_midpoints[1][1] = 0.; m_facet_edge_midpoints[1][2] = 0.; 
    m_facet_edge_midpoints[2][0] = 0.; m_facet_edge_midpoints[2][1] = 1.; m_facet_edge_midpoints[2][2] = 0.; 
    m_facet_edge_midpoints[3][0] = -1.; m_facet_edge_midpoints[3][1] = 0.; m_facet_edge_midpoints[3][2] = 0.; 

    m_facet_edge_midpoints[4] = m_facet_edge_midpoints[3]; 
    m_facet_edge_midpoints[5] = m_facet_edge_midpoints[0]; 
    m_facet_edge_midpoints[6] = m_facet_edge_midpoints[1];
    m_facet_edge_midpoints[7] = m_facet_edge_midpoints[2];
 
    m_facet_edge_midpoints[8][0] = 0.; m_facet_edge_midpoints[8][1] = 0.; m_facet_edge_midpoints[8][2] = 1.; 
    m_facet_edge_midpoints[11] = m_facet_edge_midpoints[10] = m_facet_edge_midpoints[9] = m_facet_edge_midpoints[8]; 
    
    //m_barycenter
    m_barycenter[0] = 0.;
    m_barycenter[1] = 0.;
	  m_barycenter[2] = 0.;
	
	  Point<dim> pt_c1234, pt_c2367, pt_c3478, pt_c1458, pt_c5678, pt_c1256;
    Point<dim> pt_c12, pt_c23, pt_c34, pt_c14, pt_c15, pt_c26, pt_c37, pt_c48, pt_c56, pt_c67, pt_c78, pt_c58;
    
    pt_c1234[0]=0.;pt_c1234[1]=0.;pt_c1234[2]=-1.;
    pt_c2367[0]=1.;pt_c2367[1]=0.;pt_c2367[2]=0.;
    pt_c3478[0]=0.;pt_c3478[1]=1.;pt_c3478[2]=0.;
    pt_c1458[0]=-1.;pt_c1458[1]=0.;pt_c1458[2]=0.;
    pt_c5678[0]=0.;pt_c5678[1]=0.;pt_c5678[2]=1.;
    pt_c1256[0]=0.;pt_c1256[1]=-1.;pt_c1256[2]=0.;
    pt_c12[0]=0.;pt_c12[1]=-1.;pt_c12[2]=-1.;
    pt_c23[0]=1.;pt_c23[1]=0.;pt_c23[2]=-1.;
    pt_c34[0]=0.;pt_c34[1]=1.;pt_c34[2]=-1.;
    pt_c14[0]=-1.;pt_c14[1]=0.;pt_c14[2]=-1.;
    pt_c15[0]=-1.;pt_c15[1]=-1.;pt_c15[2]=0.;
    pt_c26[0]=1.;pt_c26[1]=-1.;pt_c26[2]=0.;
    pt_c37[0]=1.;pt_c37[1]=1.;pt_c37[2]=0.;
    pt_c48[0]=-1.;pt_c48[1]=1.;pt_c48[2]=0.;
    pt_c56[0]=0.;pt_c56[1]=-1.;pt_c56[2]=1.;
    pt_c67[0]=1.;pt_c67[1]=0.;pt_c67[2]=1.;
    pt_c78[0]=0.;pt_c78[1]=1.;pt_c78[2]=1.;
    pt_c58[0]=-1.;pt_c58[1]=0.;pt_c58[2]=1.;
    
    m_facet_types[0] = QUADRILATERAL_FACET;

    m_facet_points[0].resize(4);
    vector<Point<dim> >(m_facet_points[0]).swap(m_facet_points[0]);
    m_facet_points[0][0] = pt_c12;
    m_facet_points[0][1] = pt_c1234; 
    m_facet_points[0][2] = m_barycenter;
    m_facet_points[0][3] = pt_c1256;

    m_facet_types[1] = QUADRILATERAL_FACET;

    m_facet_points[1].resize(4);
    vector<Point<dim> >(m_facet_points[1]).swap(m_facet_points[1]);
    m_facet_points[1][0] = pt_c23;
    m_facet_points[1][1] = pt_c1234; 
    m_facet_points[1][2] = m_barycenter;
    m_facet_points[1][3] = pt_c2367;

    m_facet_types[2] = QUADRILATERAL_FACET;

    m_facet_points[2].resize(4);
    vector<Point<dim> >(m_facet_points[2]).swap(m_facet_points[2]);
    m_facet_points[2][0] = pt_c34;
    m_facet_points[2][1] = pt_c1234; 
    m_facet_points[2][2] = m_barycenter;
    m_facet_points[2][3] = pt_c3478;

    m_facet_types[3] = QUADRILATERAL_FACET;

    m_facet_points[3].resize(4);
    vector<Point<dim> >(m_facet_points[3]).swap(m_facet_points[3]);
    m_facet_points[3][0] = pt_c14;
    m_facet_points[3][1] = pt_c1234; 
    m_facet_points[3][2] = m_barycenter;
    m_facet_points[3][3] = pt_c1458;

    m_facet_types[4] = QUADRILATERAL_FACET;

    m_facet_points[4].resize(4);
    vector<Point<dim> >(m_facet_points[4]).swap(m_facet_points[4]);
    m_facet_points[4][0] = pt_c15;
    m_facet_points[4][1] = pt_c1256; 
    m_facet_points[4][2] = m_barycenter;
    m_facet_points[4][3] = pt_c1458;

    m_facet_types[5] = QUADRILATERAL_FACET;

    m_facet_points[5].resize(4);
    vector<Point<dim> >(m_facet_points[5]).swap(m_facet_points[5]);
    m_facet_points[5][0] = pt_c26;
    m_facet_points[5][1] = pt_c2367; 
    m_facet_points[5][2] = m_barycenter;
    m_facet_points[5][3] = pt_c1256;

    m_facet_types[6] = QUADRILATERAL_FACET;

    m_facet_points[6].resize(4);
    vector<Point<dim> >(m_facet_points[6]).swap(m_facet_points[6]);
    m_facet_points[6][0] = pt_c37;
    m_facet_points[6][1] = pt_c3478; 
    m_facet_points[6][2] = m_barycenter;
    m_facet_points[6][3] = pt_c2367;

    m_facet_types[7] = QUADRILATERAL_FACET;

    m_facet_points[7].resize(4);
    vector<Point<dim> >(m_facet_points[7]).swap(m_facet_points[7]);
    m_facet_points[7][0] = pt_c48;
    m_facet_points[7][1] = pt_c1458; 
    m_facet_points[7][2] = m_barycenter;
    m_facet_points[7][3] = pt_c3478;

    m_facet_types[8] = QUADRILATERAL_FACET;

    m_facet_points[8].resize(4);
    vector<Point<dim> >(m_facet_points[8]).swap(m_facet_points[8]);
    m_facet_points[8][0] = pt_c56;
    m_facet_points[8][1] = pt_c1256; 
    m_facet_points[8][2] = m_barycenter;
    m_facet_points[8][3] = pt_c5678;

    m_facet_types[9] = QUADRILATERAL_FACET;

    m_facet_points[9].resize(4);
    vector<Point<dim> >(m_facet_points[9]).swap(m_facet_points[9]);
    m_facet_points[9][0] = pt_c67;
    m_facet_points[9][1] = pt_c2367; 
    m_facet_points[9][2] = m_barycenter;
    m_facet_points[9][3] = pt_c5678;

    m_facet_types[10] = QUADRILATERAL_FACET;

    m_facet_points[10].resize(4);
    vector<Point<dim> >(m_facet_points[10]).swap(m_facet_points[10]);
    m_facet_points[10][0] = pt_c78;
    m_facet_points[10][1] = pt_c3478; 
    m_facet_points[10][2] = m_barycenter;
    m_facet_points[10][3] = pt_c5678;

    m_facet_types[11] = QUADRILATERAL_FACET;

    m_facet_points[11].resize(4);
    vector<Point<dim> >(m_facet_points[11]).swap(m_facet_points[11]);
    m_facet_points[11][0] = pt_c58;
    m_facet_points[11][1] = pt_c1458; 
    m_facet_points[11][2] = m_barycenter;
    m_facet_points[11][3] = pt_c5678;

    // tested: 18/Oct/2017 AJB
    const double facet_normal_transforms[12][2][8] = {
        { { -3.0/16.0, -3.0/16.0, -1.0/16.0, -1.0/16.0, 3.0/16.0, 3.0/16.0, 1.0/16.0, 1.0/16.0 },
          { -3.0/16.0, -3.0/16.0, 3.0/16.0, 3.0/16.0, -1.0/16.0, -1.0/16.0, 1.0/16.0, 1.0/16.0 } },
        { { -1.0/16.0, -3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0, 3.0/16.0, 1.0/16.0 },
          { 3.0/16.0, -3.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0, -1.0/16.0, 1.0/16.0 } },
        { { -1.0/16.0, -1.0/16.0, -3.0/16.0, -3.0/16.0, 1.0/16.0, 1.0/16.0, 3.0/16.0, 3.0/16.0 },
          { 3.0/16.0, 3.0/16.0, -3.0/16.0, -3.0/16.0, 1.0/16.0, 1.0/16.0, -1.0/16.0, -1.0/16.0 } },
        { { -3.0/16.0, -1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0, 1.0/16.0, 3.0/16.0 },
          { -3.0/16.0, 3.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0, 1.0/16.0, -1.0/16.0 } },
        { { -3.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0 },
          { -3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0 } },
        { { 3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0 },
          { -1.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0 } },
        { { 1.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0 },
          { 1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0 } },
        { { -1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0 },
          { 3.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0 } },
        { { -1.0/16.0, -1.0/16.0, 1.0/16.0, 1.0/16.0, -3.0/16.0, -3.0/16.0, 3.0/16.0, 3.0/16.0 },
          { 3.0/16.0, 3.0/16.0, 1.0/16.0, 1.0/16.0, -3.0/16.0, -3.0/16.0, -1.0/16.0, -1.0/16.0 } },
        { { 1.0/16.0, -1.0/16.0, -1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0, -3.0/16.0, 3.0/16.0 },
          { 1.0/16.0, 3.0/16.0, 3.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0, -3.0/16.0, -1.0/16.0 } },
        { { 1.0/16.0, 1.0/16.0, -1.0/16.0, -1.0/16.0, 3.0/16.0, 3.0/16.0, -3.0/16.0, -3.0/16.0 },
          { 1.0/16.0, 1.0/16.0, 3.0/16.0, 3.0/16.0, -1.0/16.0, -1.0/16.0, -3.0/16.0, -3.0/16.0 } },
        { { -1.0/16.0, 1.0/16.0, 1.0/16.0, -1.0/16.0, -3.0/16.0, 3.0/16.0, 3.0/16.0, -3.0/16.0 },
          { 3.0/16.0, 1.0/16.0, 1.0/16.0, 3.0/16.0, -3.0/16.0, -1.0/16.0, -1.0/16.0, -3.0/16.0 } }
    };
    m_facet_normal_xforms.resize(12u);
    for (size_t f = 0; f < 12u; ++f) {
        m_facet_normal_xforms[f].resize(8u);
        for (size_t n = 0; n < 8u; ++n) {
            m_facet_normal_xforms[f][n] = std::make_pair(
                    facet_normal_transforms[f][0][n],
                    facet_normal_transforms[f][1][n]
                );
        }
    }

    // tested: 07/Nov/2007 Hamid
    Point<dim> p;
    m_sector_points[0].resize(8U); 
    p[0]=-1.;p[1]=-1.;p[2]=-1.;  
    m_sector_points[0][0] = p;
    m_sector_points[0][1] = pt_c12;
    m_sector_points[0][2] = pt_c1234;
    m_sector_points[0][3] = pt_c14; 
    m_sector_points[0][4] = pt_c15; 
    m_sector_points[0][5] = pt_c1256;
    m_sector_points[0][6] = m_barycenter;
    m_sector_points[0][7] = pt_c1458; 
    
    m_sector_points[1].resize(8U); 
    p[0]=1.;p[1]=-1.;p[2]=-1.;  
    m_sector_points[1][0] = p;
    m_sector_points[1][1] = pt_c23;
    m_sector_points[1][2] = pt_c1234;
    m_sector_points[1][3] = pt_c12; 
    m_sector_points[1][4] = pt_c26; 
    m_sector_points[1][5] = pt_c2367;
    m_sector_points[1][6] = m_barycenter;
    m_sector_points[1][7] = pt_c1256;
  
    m_sector_points[2].resize(8U); 
    p[0]=1.;p[1]=1.;p[2]=-1.;  
    m_sector_points[2][0] = p;
    m_sector_points[2][1] = pt_c34; 
    m_sector_points[2][2] = pt_c1234;
    m_sector_points[2][3] = pt_c23; 
    m_sector_points[2][4] = pt_c37; 
    m_sector_points[2][5] = pt_c3478;
    m_sector_points[2][6] = m_barycenter; 
    m_sector_points[2][7] = pt_c2367; 
  
    m_sector_points[3].resize(8U); 
    p[0]=-1.;p[1]=1.;p[2]=-1.;  
    m_sector_points[3][0] = p;
    m_sector_points[3][1] = pt_c14; 
    m_sector_points[3][2] = pt_c1234;
    m_sector_points[3][3] = pt_c34; 
    m_sector_points[3][4] = pt_c48; 
    m_sector_points[3][5] = pt_c1458;
    m_sector_points[3][6] = m_barycenter;
    m_sector_points[3][7] = pt_c3478;
    
    m_sector_points[4].resize(8U); 
    p[0]=-1.;p[1]=-1.;p[2]=1.;  
    m_sector_points[4][0] = p;
    m_sector_points[4][1] = pt_c56; 
    m_sector_points[4][2] = pt_c5678;
    m_sector_points[4][3] = pt_c58; 
    m_sector_points[4][4] = pt_c15; 
    m_sector_points[4][5] = pt_c1256;
    m_sector_points[4][6] = m_barycenter;
    m_sector_points[4][7] = pt_c1458;
    
    m_sector_points[5].resize(8U); 
    p[0]=1.;p[1]=-1.;p[2]=1.;  
    m_sector_points[5][0] = p;
    m_sector_points[5][1] = pt_c67; 
    m_sector_points[5][2] = pt_c5678;
    m_sector_points[5][3] = pt_c56; 
    m_sector_points[5][4] = pt_c26; 
    m_sector_points[5][5] = pt_c2367;
    m_sector_points[5][6] = m_barycenter;
    m_sector_points[5][7] = pt_c1256;
  
    m_sector_points[6].resize(8U); 
    p[0]=1.;p[1]=1.;p[2]=1.;  
    m_sector_points[6][0] = p;
    m_sector_points[6][1] = pt_c78; 
    m_sector_points[6][2] = pt_c5678;
    m_sector_points[6][3] = pt_c67; 
    m_sector_points[6][4] = pt_c37; 
    m_sector_points[6][5] = pt_c3478;
    m_sector_points[6][6] = m_barycenter;
    m_sector_points[6][7] = pt_c2367;
  
    m_sector_points[7].resize(8U);
    p[0]=-1.;p[1]=1.;p[2]=1.;  
    m_sector_points[7][0] = p;
    m_sector_points[7][1] = pt_c58; 
    m_sector_points[7][2] = pt_c5678;
    m_sector_points[7][3] = pt_c78; 
    m_sector_points[7][4] = pt_c48; 
    m_sector_points[7][5] = pt_c1458;
    m_sector_points[7][6] = m_barycenter;
    m_sector_points[7][7] = pt_c3478; 

    m_edge_of_sectors.resize(8U);
    for ( size_t iSector=0U; iSector<8U; iSector++ ) 
      { 
         m_edge_of_sectors[iSector].resize(12U);
         for ( size_t iLine=0U; iLine<4U; iLine++ )
            m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine+1)%4); 
          for ( size_t iLine=4U; iLine<8U; iLine++ )
            m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine-3)%4+4);  
          for ( size_t iLine=8U; iLine<12U; iLine++ ) 
            m_edge_of_sectors[iSector][iLine]=make_pair(iLine-8,iLine-4);  
      }  

} // CreateDataFor_ISOPARAMETRIC_LINEAR_HEX()







/**
 
Generates data for subdivision of unit isoparametric linear
prism to the CV sectors.
*/
template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_PRISM()
 {
//	  const size_t NumberOfPhysicalDimensions(3U);
	  enum { NumberOfPhysicalDimensions = 3, 
           NumOfInternalFacets        = 9,
           NumOfInternalFacetsPerNode = 3,
           NumOfInternalVolumes       = 6, // == VolMultipliers
  	       NumOfIPperVolume           = 1,
           NumOfIPperFacet            = 1 };

    Resize( NumOfInternalFacets,
            NumOfInternalFacetsPerNode,
            NumOfInternalVolumes,
            NumOfIPperFacet,
            NumOfIPperVolume );

    // Map volume multipliers to new vectors
    for(size_t i=0;i<NumOfInternalVolumes;i++)
      for(size_t j=0;j<NumOfIPperVolume;j++)
       // Hard coded, check multiple points!
       m_volume_integration_weights[i][j]=1./6.;

    // products of detJ's, weights & conversion factors
    m_facet_integration_weights[0][0]=m_facet_integration_weights[2][0]=m_facet_integration_weights[6][0]=m_facet_integration_weights[8][0]=sqrt(5.)/6.;
    m_facet_integration_weights[1][0]=m_facet_integration_weights[7][0]=sqrt(2.)/6.;
    m_facet_integration_weights[3][0]=m_facet_integration_weights[4][0]=m_facet_integration_weights[5][0]=1./6.;
  
    // projection weights
    m_projection_weights[0][0]=0.20833323822827704;
    m_projection_weights[1][0]=0.145806006795383303;
    m_projection_weights[2][0]=0.145828940771960402;
    m_projection_weights[3][0]=0.197916676224492188;
    m_projection_weights[4][0]=0.197916157089035288;
    m_projection_weights[5][0]=0.104250300417627992;
    m_projection_weights[6][0]=0.208333238228277012;
    m_projection_weights[7][0]=0.145806006795383247;
    m_projection_weights[8][0]=0.145828940771960402;

	// face integration points
    const double fip[9U][NumberOfPhysicalDimensions] =
    {
      {5./12.,1./6.,-1./2.},
      {5./12.,5./12.,-1./2.},
      {1./6.,5./12.,-1./2.},
      {7./36.,7./36.,0.},
      {11./18.,7./36.,0.},
      {7./36.,11./18.,0.},
      {5./12.,1./6.,1./2.},
      {5./12.,5./12.,1./2.},
      {1./6.,5./12.,1./2.}
    }; //end of faces integration points initialisation list

    // Map to new vectors
    for(size_t i=0;i<NumOfInternalFacets;i++)
     for(size_t j=0;j<NumOfIPperFacet; j++)
      for(size_t k=0;k<NumberOfPhysicalDimensions; k++) 
         m_facet_integration_points[i][j][k]=fip[i*NumOfIPperFacet+j][k];
    		
    // Subdivided Tet Volume integration points
    const double vip[6U][NumberOfPhysicalDimensions]=
    {
     {5./24.,5./24.,-1./2.},
     {7./12.,5./24.,-1./2.}, 
     {5./24.,7./12.,-1./2.},
     {5./24.,5./24.,1./2.}, 
     {7./12.,5./24.,1./2.},
     {5./24.,7./12.,1./2.} 
    }; 
    
    // m_volume_integration_points1;   // [ivol][vpts][dim]
    for(size_t i=0;i<NumOfInternalVolumes;i++)
      for(size_t j=0;j<NumOfIPperVolume;j++)
    	for(size_t k=0;k<NumberOfPhysicalDimensions; k++)
    		m_volume_integration_points1[i][j][k]=vip[i*NumOfIPperVolume+j][k];
    		
    const double cosAlpha1(static_cast<double>(2./sqrt(5.0)));
    const double sinAlpha1(static_cast<double>(1./sqrt(5.0)));
   
    const double cosAlpha2(static_cast<double>(sqrt(2.0)/2.0));
    const double sinAlpha2(cosAlpha2);


   // Facet normals pointing outside the internal face, all together 6 normals!
    // Normal 0, +r facet IP[0] OK
    m_facet_normals[0][0]=cosAlpha1; m_facet_normals[0][1]=sinAlpha1; m_facet_normals[0][2]=0.; 
    // Normal 1 On the bisector
    m_facet_normals[1][0]=-cosAlpha2; m_facet_normals[1][1]=sinAlpha2; m_facet_normals[1][2]=0.;
    // Normal 2 -s Check orientation!
    m_facet_normals[2][0]=-sinAlpha1;  m_facet_normals[2][1]=-cosAlpha1;m_facet_normals[2][2]=0.;
    // Normal 3 +t, Up
    m_facet_normals[3][0]=0.0; m_facet_normals[3][1]=0.0; m_facet_normals[3][2]=1.;
    // Normal 4+t, Up
    m_facet_normals[4][0]=0.0; m_facet_normals[4][1]=0.0; m_facet_normals[4][2]=1.;
    // Normal 5+t, Up
    m_facet_normals[5][0]=0.0; m_facet_normals[5][1]=0.0; m_facet_normals[5][2]=1.;
    // Normal 6, +t  
    m_facet_normals[6]=m_facet_normals[0];
    // Normal 7,+r 
    m_facet_normals[7]=m_facet_normals[1];
    // Normal 8 -s Check orientation!
    m_facet_normals[8]=m_facet_normals[2];

    // facets surrounding node
    m_facets_surrounding_node[0][0]=0; m_facets_surrounding_node[0][1]=2; m_facets_surrounding_node[0][2]=3;    
    m_facets_surrounding_node[1][0]=1; m_facets_surrounding_node[1][1]=0; m_facets_surrounding_node[1][2]=4;
    m_facets_surrounding_node[2][0]=2; m_facets_surrounding_node[2][1]=1; m_facets_surrounding_node[2][2]=5;
    m_facets_surrounding_node[3][0]=8; m_facets_surrounding_node[3][1]=6; m_facets_surrounding_node[3][2]=3;
    m_facets_surrounding_node[4][0]=6; m_facets_surrounding_node[4][1]=7; m_facets_surrounding_node[4][2]=4;    
    m_facets_surrounding_node[5][0]=7; m_facets_surrounding_node[5][1]=8; m_facets_surrounding_node[5][2]=5;
    
    // Edges
    m_edges_of_element.resize(9U);
    vector<pair<size_t,size_t> >(m_edges_of_element).swap(m_edges_of_element);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,0U);

    m_edges_of_element[3]=make_pair(0U,3U);
    m_edges_of_element[4]=make_pair(1U,4U);    
    m_edges_of_element[5]=make_pair(2U,5U);
    
    m_edges_of_element[6]=make_pair(3U,4U);
    m_edges_of_element[7]=make_pair(4U,5U);            
    m_edges_of_element[8]=make_pair(5U,3U);    

    //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 1./2.; m_facet_edge_midpoints[0][1] =    0.; m_facet_edge_midpoints[0][2] = 0.; 
    m_facet_edge_midpoints[1][0] = 1./2.; m_facet_edge_midpoints[1][1] = 1./2.; m_facet_edge_midpoints[1][2] = 0.; 
    m_facet_edge_midpoints[2][0] = 0.;    m_facet_edge_midpoints[2][1] = 1./2.; m_facet_edge_midpoints[2][2] = 0.; 
    
   	m_facet_edge_midpoints[3] = m_facet_edge_midpoints[2]; 
    m_facet_edge_midpoints[4] = m_facet_edge_midpoints[0]; 
    m_facet_edge_midpoints[5] = m_facet_edge_midpoints[1]; 
    
    m_facet_edge_midpoints[6][0] = 1./3.; m_facet_edge_midpoints[6][1] =  1./3.; m_facet_edge_midpoints[6][2] = 1.; 
    m_facet_edge_midpoints[8] = m_facet_edge_midpoints[7] = m_facet_edge_midpoints[6]; 
    
    //m_barycenter 
    m_barycenter[0] = 1./3.;
    m_barycenter[1] = 1./3.;
	  m_barycenter[2] = 0.;

	  Point<dim> pt_c1245, pt_c2356, pt_c1346, pt_c123, pt_c456;
    Point<dim> pt_c12, pt_c13, pt_c23, pt_c25, pt_c14, pt_c45, pt_c56, pt_c46, pt_c36;
    
    pt_c1245[0]=1./2.;pt_c1245[1]=0.;pt_c1245[2]=0.;
    pt_c2356[0]=1./2.;pt_c2356[1]=1./2.;pt_c2356[2]=0.;
    pt_c1346[0]=0.;pt_c1346[1]=1./2.;pt_c1346[2]=0.;
    pt_c123[0]=1./3.;pt_c123[1]=1./3.;pt_c123[2]=-1.;
    pt_c456[0]=1./3.;pt_c456[1]=1./3.;pt_c456[2]=1.;
   
    pt_c12[0]=1./2.;pt_c12[1]=0.;pt_c12[2]=-1.;
    pt_c13[0]=0.;pt_c13[1]=1./2.;pt_c13[2]=-1.;
    pt_c23[0]=1./2.;pt_c23[1]=1./2.;pt_c23[2]=-1.;
    pt_c25[0]=1.;pt_c25[1]=0.;pt_c25[2]=0.;
    pt_c14[0]=0.;pt_c14[1]=0.;pt_c14[2]=0.;
    pt_c45[0]=1./2.;pt_c45[1]=0.;pt_c45[2]=1.;
    pt_c56[0]=1./2.;pt_c56[1]=1./2.;pt_c56[2]=1.;
    pt_c46[0]=0.;pt_c46[1]=1./2.;pt_c46[2]=1.;
    pt_c36[0]=0.;pt_c36[1]=1.;pt_c36[2]=0.;

    m_facet_types[0] = QUADRILATERAL_FACET;

	m_facet_points[0].resize(4);
    vector<Point<dim> >(m_facet_points[0]).swap(m_facet_points[0]);
    m_facet_points[0][0] = pt_c12;
    m_facet_points[0][1] = pt_c123; 
    m_facet_points[0][2] = m_barycenter;
    m_facet_points[0][3] = pt_c1245;

    m_facet_types[1] = QUADRILATERAL_FACET;

    m_facet_points[1].resize(4);
    vector<Point<dim> >(m_facet_points[1]).swap(m_facet_points[1]);
    m_facet_points[1][0] = pt_c123;
    m_facet_points[1][1] = m_barycenter; 
    m_facet_points[1][2] = pt_c2356;
    m_facet_points[1][3] = pt_c23;

    m_facet_types[2] = QUADRILATERAL_FACET;

	m_facet_points[2].resize(4);
    vector<Point<dim> >(m_facet_points[2]).swap(m_facet_points[2]);
    m_facet_points[2][0] = pt_c123;
    m_facet_points[2][1] = m_barycenter; 
    m_facet_points[2][2] = pt_c1346;
    m_facet_points[2][3] = pt_c13;

    m_facet_types[3] = QUADRILATERAL_FACET;

	m_facet_points[3].resize(4);
    vector<Point<dim> >(m_facet_points[3]).swap(m_facet_points[3]);
    m_facet_points[3][0] = pt_c14;
    m_facet_points[3][1] = pt_c1245; 
    m_facet_points[3][2] = m_barycenter;
    m_facet_points[3][3] = pt_c1346;

    m_facet_types[4] = QUADRILATERAL_FACET;

	m_facet_points[4].resize(4);
    vector<Point<dim> >(m_facet_points[4]).swap(m_facet_points[4]);
    m_facet_points[4][0] = pt_c25;
    m_facet_points[4][1] = pt_c2356; 
    m_facet_points[4][2] = m_barycenter;
    m_facet_points[4][3] = pt_c1245;

    m_facet_types[5] = QUADRILATERAL_FACET;

	m_facet_points[5].resize(4);
    vector<Point<dim> >(m_facet_points[5]).swap(m_facet_points[5]);
    m_facet_points[5][0] = pt_c36;
    m_facet_points[5][1] = pt_c1346; 
    m_facet_points[5][2] = m_barycenter;
    m_facet_points[5][3] = pt_c2356;

    m_facet_types[6] = QUADRILATERAL_FACET;

	m_facet_points[6].resize(4);
    vector<Point<dim> >(m_facet_points[6]).swap(m_facet_points[6]);
    m_facet_points[6][0] = pt_c45;
    m_facet_points[6][1] = pt_c1245; 
    m_facet_points[6][2] = m_barycenter;
    m_facet_points[6][3] = pt_c456;

    m_facet_types[7] = QUADRILATERAL_FACET;

    m_facet_points[7].resize(4);
    vector<Point<dim> >(m_facet_points[7]).swap(m_facet_points[7]);
    m_facet_points[7][0] = pt_c56;
    m_facet_points[7][1] = pt_c2356; 
    m_facet_points[7][2] = m_barycenter;
    m_facet_points[7][3] = pt_c456;

    m_facet_types[8] = QUADRILATERAL_FACET;

	m_facet_points[8].resize(4);
    vector<Point<dim> >(m_facet_points[8]).swap(m_facet_points[8]);
    m_facet_points[8][0] = pt_c46;
    m_facet_points[8][1] = pt_c1346; 
    m_facet_points[8][2] = m_barycenter;
    m_facet_points[8][3] = pt_c456;

    // tested: 18/Oct/2017 AJB
    const double facet_normal_transforms[9][2][6] = {
        { { -5.0/24.0, -5.0/24.0, -1.0/12.0, 5.0/24.0, 5.0/24.0, 1.0/12.0 },
          { -1.0/8.0, -1.0/8.0, 1.0/4.0, -1.0/24.0, -1.0/24.0, 1.0/12.0 } },
        { { -1.0/4.0, 1.0/8.0, 1.0/8.0, -1.0/12.0, 1.0/24.0, 1.0/24.0 },
          { -1.0/12.0, -5.0/24.0, -5.0/24.0, 1.0/12.0, 5.0/24.0, 5.0/24.0 } },
        { { 1.0/8.0, -1.0/4.0, 1.0/8.0, 1.0/24.0, -1.0/12.0, 1.0/24.0 },
          { -5.0/24.0, -1.0/12.0, -5.0/24.0, 5.0/24.0, 1.0/12.0, 5.0/24.0 } },
        { { -1.0/6.0, -1.0/24.0, 5.0/24.0, -1.0/6.0, -1.0/24.0, 5.0/24.0 },
          { -1.0/6.0, 5.0/24.0, -1.0/24.0, -1.0/6.0, 5.0/24.0, -1.0/24.0 } },
        { { 5.0/24.0, -1.0/6.0, -1.0/24.0, 5.0/24.0, -1.0/6.0, -1.0/24.0 },
          { -1.0/24.0, -1.0/6.0, 5.0/24.0, -1.0/24.0, -1.0/6.0, 5.0/24.0 } },
        { { -1.0/24.0, 5.0/24.0, -1.0/6.0, -1.0/24.0, 5.0/24.0, -1.0/6.0 },
          { 5.0/24.0, -1.0/24.0, -1.0/6.0, 5.0/24.0, -1.0/24.0, -1.0/6.0 } },
        { { -1.0/24.0, -1.0/24.0, 1.0/12.0, -1.0/8.0, -1.0/8.0, 1.0/4.0 },
          { 5.0/24.0, 5.0/24.0, 1.0/12.0, -5.0/24.0, -5.0/24.0, -1.0/12.0 } },
        { { 1.0/12.0, -1.0/24.0, -1.0/24.0, 1.0/4.0, -1.0/8.0, -1.0/8.0 },
          { 1.0/12.0, 5.0/24.0, 5.0/24.0, -1.0/12.0, -5.0/24.0, -5.0/24.0 } },
        { { -1.0/24.0, 1.0/12.0, -1.0/24.0, -1.0/8.0, 1.0/4.0, -1.0/8.0 },
          { 5.0/24.0, 1.0/12.0, 5.0/24.0, -5.0/24.0, -1.0/12.0, -5.0/24.0 } }
    };
    m_facet_normal_xforms.resize(9u);
    for (size_t f = 0; f < 9u; ++f) {
        m_facet_normal_xforms[f].resize(6u);
        for (size_t n = 0; n < 6u; ++n) {
            m_facet_normal_xforms[f][n] = std::make_pair(
                    facet_normal_transforms[f][0][n],
                    facet_normal_transforms[f][1][n]
                );
        }
    }

    // tested: 08/Nov/2007 Hamid	
    Point<dim> p;    
    m_sector_points[0].resize(8U); 
    p[0] = 0.; p[1] = 0.; p[2] = -1.;
    m_sector_points[0][0] = p;
    m_sector_points[0][1] = pt_c12; 
    m_sector_points[0][2] = pt_c123;
    m_sector_points[0][3] = pt_c13; 
    m_sector_points[0][4] = pt_c14; 
    m_sector_points[0][5] = pt_c1245;
    m_sector_points[0][6] = m_barycenter;
    m_sector_points[0][7] = pt_c1346;
    
    m_sector_points[1].resize(8U); 
    p[0] = 1.; p[1] = 0.; p[2] = -1.;
    m_sector_points[1][0] = p;
    m_sector_points[1][1] = pt_c23; 
    m_sector_points[1][2] = pt_c123;
    m_sector_points[1][3] = pt_c12; 
    m_sector_points[1][4] = pt_c25; 
    m_sector_points[1][5] = pt_c2356;
    m_sector_points[1][6] = m_barycenter;
    m_sector_points[1][7] = pt_c1245;
  
    m_sector_points[2].resize(8U); 
    p[0] = 0.; p[1] = 1.; p[2] = -1.;
    m_sector_points[2][0] = p;
    m_sector_points[2][1] = pt_c13; 
    m_sector_points[2][2] = pt_c123;
    m_sector_points[2][3] = pt_c23; 
    m_sector_points[2][4] = pt_c36; 
    m_sector_points[2][5] = pt_c1346;
    m_sector_points[2][6] = m_barycenter;
    m_sector_points[2][7] = pt_c2356;
    
    m_sector_points[3].resize(8U); 
    p[0] = 0.; p[1] = 0.; p[2] = 1.;
    m_sector_points[3][0] = p;
    m_sector_points[3][1] = pt_c45;
    m_sector_points[3][2] = pt_c456; 
    m_sector_points[3][3] = pt_c46; 
    m_sector_points[3][4] = pt_c14; 
    m_sector_points[3][5] = pt_c1245;
    m_sector_points[3][6] = m_barycenter;
    m_sector_points[3][7] = pt_c1346; 
    
    m_sector_points[4].resize(8U); 
    p[0] = 1.; p[1] = 0.; p[2] = 1.;
    m_sector_points[4][0] = p;
    m_sector_points[4][1] = pt_c56;
    m_sector_points[4][2] = pt_c456;
    m_sector_points[4][3] = pt_c45; 
    m_sector_points[4][4] = pt_c25; 
    m_sector_points[4][5] = pt_c2356;
    m_sector_points[4][6] = m_barycenter;
    m_sector_points[4][7] = pt_c1245;
  
    m_sector_points[5].resize(8U); 
    p[0] = 0.; p[1] = 1.; p[2] = 1.;
    m_sector_points[5][0] = p;
    m_sector_points[5][1] = pt_c46;
    m_sector_points[5][2] = pt_c456;
    m_sector_points[5][3] = pt_c56; 
    m_sector_points[5][4] = pt_c36; 
    m_sector_points[5][5] = pt_c1346;
    m_sector_points[5][6] = m_barycenter;
    m_sector_points[5][7] = pt_c2356;

    m_edge_of_sectors.resize(6U); 
    for ( size_t iSector=0U; iSector<6U; iSector++ ) 
      { 
        m_edge_of_sectors[iSector].resize(12U);
        for ( size_t iLine=0U; iLine<4U; iLine++ )
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine+1)%4);  //   [iSector][iLine] 
        for ( size_t iLine=4U; iLine<8U; iLine++ )  
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine-3)%4+4);  
        for ( size_t iLine=8U; iLine<12U; iLine++ ) 
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine-8,iLine-4);  
      }  
  
    
} // CreateDataFor_ISOPARAMETRIC_LINEAR_PRISM




/**
 
Generates data for subdivision of unit isoparametric linear
pyramid to the CV sectors.

NB: The pyramid is a special element type since it has two sector
integration points in the apex, but only a single one in each of the
other 4 sectors.
*/
#ifndef PYRAMID_TRIANGULAR_FACETS

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_PYR()
 {
    const size_t NumOfInternalFacets(8U);
    const size_t NumOfInternalFacetsPerNode(4U); //Special case of 4 facets?
    const size_t NumOfInternalVolumes(5U); // == VolMultipliers
  	const size_t NumOfIPperVolume(1U);
    const size_t NumOfIPperFacet(1U);

    Resize( NumOfInternalFacets,
            NumOfInternalFacetsPerNode,
            NumOfInternalVolumes,
            NumOfIPperFacet,
            NumOfIPperVolume );
            
    for ( size_t i=0U; i<NumOfInternalVolumes; i++ )	
      m_volume_integration_weights[i][0]=1./4.;//0.24537037;//
    
    m_volume_integration_weights[4][0]=1./3.;//0.351851852;//
    
    // Subdivided pyramid volume integration points
    const double vip[5U][3U]=
      {
         {-31./72.,-31./72.,23./144.},
         {31./72.,-31./72.,23./144.},
         {31./72.,31./72.,23./144.},
         {-31./72.,31./72.,23./144.},
         {0.,0.,25./48.}
      }; 
    
    // map to m_volume_integration_points1;   // [ivol][vpts][dim]
    for( size_t i=0U; i<NumOfInternalVolumes; i++ )
      for( size_t j=0U; j<NumOfIPperVolume; j++ )
    	for( size_t k=0U; k<3U; k++ )
    	  m_volume_integration_points1[i][j][k]=vip[i*NumOfIPperVolume+j][k];


	  //face integration points
    const double fip[8U][3U] =
   	 {
  	   {0.,-4./9., 5./36.},
         {4./9.,0.,5./36.},
         {0.,4./9.,5./36.},
         {-4./9.,0.,5./36.},
  	   {-7./24.,-7./24.,17./48.},
         {7./24.,-7./24.,17./48.},
         {7./24.,7./24.,17./48.},
         {-7./24.,7./24.,17./48.}
    }; //end of faces integration points initialisation list
    
     for(size_t i=0;i<NumOfInternalFacets;i++)
      for(size_t j=0;j<NumOfIPperFacet; j++)
       for(size_t k=0;k<3U; k++) 
          m_facet_integration_points[i][j][k]=fip[i*NumOfIPperFacet+j][k];    		
   
     // the facet integration weights is the area of the facet in parametric space
     m_facet_integration_weights[0][0]=m_facet_integration_weights[1][0]=m_facet_integration_weights[2][0]=m_facet_integration_weights[3][0]=0.25;
     m_facet_integration_weights[4][0]=m_facet_integration_weights[5][0]=m_facet_integration_weights[6][0]=m_facet_integration_weights[7][0]=sqrt(2.)/4.;
      
     // Projection multipliers from pyramid to unit cube
     m_projection_weights[0][0]=m_projection_weights[1][0]=m_projection_weights[2][0]=m_projection_weights[3][0]=0.3648003127564263;
     m_projection_weights[4][0]=m_projection_weights[5][0]=m_projection_weights[6][0]=m_projection_weights[7][0]=0.2085503283484722;
    
    // Facet normals pointing outside the internal face, all together 6 normals!
    //normal 0 
    m_facet_normals[0][0]=1.0; m_facet_normals[0][1]=0.0; m_facet_normals[0][2]=0.0; 
    //normal 1
    m_facet_normals[1][0]=0.0; m_facet_normals[1][1]=1.0; m_facet_normals[1][2]=0.0;
    //normal 2
    m_facet_normals[2][0]=-1.0; m_facet_normals[2][1]=0.0; m_facet_normals[2][2]=0.0;      
	  //normal 3
    m_facet_normals[3][0]=0.0; m_facet_normals[3][1]=-1.0; m_facet_normals[3][2]=0.0;      
    //normal 4
    m_facet_normals[4][0]=sqrt(8.)/12.;m_facet_normals[4][1]=sqrt(8.)/12.;m_facet_normals[4][2]=sqrt(8.)/3.;
    //normal 5
    m_facet_normals[5][0]=-sqrt(8.)/12.;m_facet_normals[5][1]=sqrt(8.)/12.;m_facet_normals[5][2]=sqrt(8.)/3.;
    //normal 6
    m_facet_normals[6][0]=-sqrt(8.)/12.;m_facet_normals[6][1]=-sqrt(8.)/12.;m_facet_normals[6][2]=sqrt(8.)/3.;
    //normal 7 
    m_facet_normals[7][0]=sqrt(8.)/12.;m_facet_normals[7][1]=-sqrt(8.)/12.;m_facet_normals[7][2]=sqrt(8.)/3.;
     
    // Facet surrounding node
    m_facets_surrounding_node[0].resize(3U);
    m_facets_surrounding_node[0][0]=0; m_facets_surrounding_node[0][1]=3; m_facets_surrounding_node[0][2]=4;    
    m_facets_surrounding_node[1].resize(3U);
    m_facets_surrounding_node[1][0]=1; m_facets_surrounding_node[1][1]=0; m_facets_surrounding_node[1][2]=5;
    m_facets_surrounding_node[2].resize(3U);
    m_facets_surrounding_node[2][0]=2; m_facets_surrounding_node[2][1]=1; m_facets_surrounding_node[2][2]=6;
	  m_facets_surrounding_node[3].resize(3U);
    m_facets_surrounding_node[3][0]=3; m_facets_surrounding_node[3][1]=2; m_facets_surrounding_node[3][2]=7;
	  m_facets_surrounding_node[4].resize(4U); // apex node
    m_facets_surrounding_node[4][0]=4; m_facets_surrounding_node[4][1]=7; m_facets_surrounding_node[4][2]=6;
    m_facets_surrounding_node[4][3]=5;

    //Edges   
    m_edges_of_element.resize(8U);
    vector<pair<size_t,size_t> >(m_edges_of_element).swap(m_edges_of_element);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,3U);
    m_edges_of_element[3]=make_pair(3U,0U);
 
    m_edges_of_element[4]=make_pair(0U,4U);    
    m_edges_of_element[5]=make_pair(1U,4U);
    m_edges_of_element[6]=make_pair(2U,4U);
    m_edges_of_element[7]=make_pair(3U,4U);            
 
    //Facet mid-edge points 
    //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 0.;     m_facet_edge_midpoints[0][1] = -2./3.; m_facet_edge_midpoints[0][2] = 1./3.; 
    m_facet_edge_midpoints[1][0] = 2./3.;  m_facet_edge_midpoints[1][1] = 0.;     m_facet_edge_midpoints[1][2] = 1./3.; 
    m_facet_edge_midpoints[2][0] = 0.;     m_facet_edge_midpoints[2][1] =  2./3.; m_facet_edge_midpoints[2][2] = 1./3.; 
    m_facet_edge_midpoints[3][0] = -2./3.; m_facet_edge_midpoints[3][1] = 0.;     m_facet_edge_midpoints[3][2] = 1./3.; 
    
    m_facet_edge_midpoints[4] = m_facet_edge_midpoints[3]; 
    m_facet_edge_midpoints[5] = m_facet_edge_midpoints[0]; 
    m_facet_edge_midpoints[6] = m_facet_edge_midpoints[1]; 
    m_facet_edge_midpoints[7] = m_facet_edge_midpoints[2]; 
    
    //m_barycenter 
    m_barycenter[0] = 0.;
    m_barycenter[1] = 0.;
	  m_barycenter[2] = 1./4.; 

    Point<dim> pt_c125, pt_c145, pt_c235, pt_c345, pt_c1234;
  	Point<dim> pt_c12, pt_c23, pt_c34, pt_c14, pt_c15, pt_c25, pt_c35, pt_c45;
  	
  	pt_c12[0]=0.;pt_c12[1]=-1.;pt_c12[2]=0.;
  	pt_c23[0]=1.;pt_c23[1]=0.;pt_c23[2]=0.;
  	pt_c34[0]=0.;pt_c34[1]=1.;pt_c34[2]=0.;
  	pt_c14[0]=-1.;pt_c14[1]=0.;pt_c14[2]=0.;
  	pt_c15[0]=-1./2.;pt_c15[1]=-1./2.;pt_c15[2]=1./2.;
  	pt_c25[0]=1./2.;pt_c25[1]=-1./2.;pt_c25[2]=1./2.;
  	pt_c35[0]=1./2.;pt_c35[1]=1./2.;pt_c35[2]=1./2.;
  	pt_c45[0]=-1./2.;pt_c45[1]=1./2.;pt_c45[2]=1./2.;
  	
  	pt_c125[0]=0.;pt_c125[1]=-2./3.;pt_c125[2]=1./3.;
  	pt_c145[0]=-2./3.;pt_c145[1]=0.;pt_c145[2]=1./3.;
  	pt_c235[0]=2./3.;pt_c235[1]=0.;pt_c235[2]=1./3.;
  	pt_c345[0]=0.;pt_c345[1]=2./3.;pt_c345[2]=1./3.;
  	pt_c1234[0]=0.;pt_c1234[1]=0.;pt_c1234[2]=0.;

    m_facet_types[0] = QUADRILATERAL_FACET;

  	m_facet_points[0].resize(4);
    vector<Point<dim> >(m_facet_points[0]).swap(m_facet_points[0]);
    m_facet_points[0][0] = pt_c12;
    m_facet_points[0][1] = pt_c1234; 
    m_facet_points[0][2] = m_barycenter;
    m_facet_points[0][3] = pt_c125;

    m_facet_types[1] = QUADRILATERAL_FACET;

	m_facet_points[1].resize(4);
    vector<Point<dim> >(m_facet_points[1]).swap(m_facet_points[1]);
    m_facet_points[1][0] = pt_c23;
    m_facet_points[1][1] = pt_c1234; 
    m_facet_points[1][2] = m_barycenter;
    m_facet_points[1][3] = pt_c235;    

    m_facet_types[2] = QUADRILATERAL_FACET;

  	m_facet_points[2].resize(4);
    vector<Point<dim> >(m_facet_points[2]).swap(m_facet_points[2]);
    m_facet_points[2][0] = pt_c34;
    m_facet_points[2][1] = pt_c1234; 
    m_facet_points[2][2] = m_barycenter;
    m_facet_points[2][3] = pt_c345;    

    m_facet_types[3] = QUADRILATERAL_FACET;

	m_facet_points[3].resize(4);
	vector<Point<dim> >(m_facet_points[3]).swap(m_facet_points[3]);
    m_facet_points[3][0] = pt_c14;
    m_facet_points[3][1] = pt_c1234; 
    m_facet_points[3][2] = m_barycenter;
    m_facet_points[3][3] = pt_c145;    

    m_facet_types[4] = QUADRILATERAL_FACET;

	m_facet_points[4].resize(4);
    vector<Point<dim> >(m_facet_points[4]).swap(m_facet_points[4]);
    m_facet_points[4][0] = pt_c15;
    m_facet_points[4][1] = pt_c125; 
    m_facet_points[4][2] = m_barycenter;
    m_facet_points[4][3] = pt_c145;    

    m_facet_types[5] = QUADRILATERAL_FACET;

	m_facet_points[5].resize(4);
    vector<Point<dim> >(m_facet_points[5]).swap(m_facet_points[5]);
    m_facet_points[5][0] = pt_c25;
    m_facet_points[5][1] = pt_c235; 
    m_facet_points[5][2] = m_barycenter;
    m_facet_points[5][3] = pt_c125;    

    m_facet_types[6] = QUADRILATERAL_FACET;

	m_facet_points[6].resize(4);
    vector<Point<dim> >(m_facet_points[6]).swap(m_facet_points[6]);
    m_facet_points[6][0] = pt_c35;
    m_facet_points[6][1] = pt_c345; 
    m_facet_points[6][2] = m_barycenter;
    m_facet_points[6][3] = pt_c235;    

    m_facet_types[7] = QUADRILATERAL_FACET;

	m_facet_points[7].resize(4);
    vector<Point<dim> >(m_facet_points[7]).swap(m_facet_points[7]);
    m_facet_points[7][0] = pt_c45;
    m_facet_points[7][1] = pt_c145; 
    m_facet_points[7][2] = m_barycenter;
    m_facet_points[7][3] = pt_c345;    

    // tested: 07/Nov/2007 Hamid	
    Point<dim> p;
    m_sector_points[0].resize(8U); 
    p[0]=-1.;p[1]=-1.;p[2]=0.;  
    m_sector_points[0][0] = p;
    m_sector_points[0][1] = pt_c12; 
    m_sector_points[0][2] = pt_c1234;
    m_sector_points[0][3] = pt_c14; 
    m_sector_points[0][4] = pt_c15; 
    m_sector_points[0][5] = pt_c125;
    m_sector_points[0][6] = m_barycenter;
    m_sector_points[0][7] = pt_c145;
    
    m_sector_points[1].resize(8U); 
    p[0]=1.;p[1]=-1.;p[2]=0.;  
    m_sector_points[1][0] = p;
    m_sector_points[1][1] = pt_c23; 
    m_sector_points[1][2] = pt_c1234;
    m_sector_points[1][3] = pt_c12; 
    m_sector_points[1][4] = pt_c25; 
    m_sector_points[1][5] = pt_c235 ;
    m_sector_points[1][6] = m_barycenter;
    m_sector_points[1][7] = pt_c125; 
  
    m_sector_points[2].resize(8U); 
    p[0]=1.;p[1]=1.;p[2]=0.;  
    m_sector_points[2][0] = p;
    m_sector_points[2][1] = pt_c34; 
    m_sector_points[2][2] = pt_c1234;
    m_sector_points[2][3] = pt_c23; 
    m_sector_points[2][4] = pt_c35; 
    m_sector_points[2][5] = pt_c345;
    m_sector_points[2][6] = m_barycenter;
    m_sector_points[2][7] = pt_c235; 
    
    m_sector_points[3].resize(8U); 
    p[0]=-1.;p[1]=1.;p[2]=0.;  
    m_sector_points[3][0] = p;
    m_sector_points[3][1] = pt_c14; 
    m_sector_points[3][2] = pt_c1234;
    m_sector_points[3][3] = pt_c34; 
    m_sector_points[3][4] = pt_c45; 
    m_sector_points[3][5] = pt_c145;
    m_sector_points[3][6] = m_barycenter;
    m_sector_points[3][7] = pt_c345;

    m_sector_points[4].resize(10U);
    p[0]=0.;p[1]=0.;p[2]=1.;  
    m_sector_points[4][0] = pt_c15; 
    m_sector_points[4][1] = pt_c125;
    m_sector_points[4][2] = pt_c25; 
    m_sector_points[4][3] = pt_c235;
    m_sector_points[4][4] = pt_c35; 
    m_sector_points[4][5] = pt_c345;
    m_sector_points[4][6] = pt_c45; 
    m_sector_points[4][7] = pt_c145;
    m_sector_points[4][8] = p;
    m_sector_points[4][9] = m_barycenter; 
    
    m_edge_of_sectors.resize(5U); 
    for ( size_t iSector=0U; iSector<4U; iSector++ ) 
      { 
        m_edge_of_sectors[iSector].resize(13U);
        for ( size_t iLine=0U; iLine<4U; iLine++ )
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine+1)%4);  //   [iSector][iLine] 
        for ( size_t iLine=4U; iLine<8U; iLine++ )
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine-3)%4+4);  
        for ( size_t iLine=8U; iLine<12U; iLine++ ) 
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine-8,iLine-4);  
        m_edge_of_sectors[iSector][12]=make_pair(4U,6U);  
      }  
    m_edge_of_sectors[4].resize(20U);
    for ( size_t iLine=0U; iLine<8U; iLine++ )
      m_edge_of_sectors[4][iLine]=make_pair(iLine,(iLine+1)%8);  //   [iSector][iLine] 
    size_t j=0U;
    for ( size_t iLine=8U; iLine<12U; iLine++,j++,j++ )
      m_edge_of_sectors[4][iLine]=make_pair(8U,j);  
    j=0U;
    for ( size_t iLine=12U; iLine<20U; iLine++ , j++) 
      m_edge_of_sectors[4][iLine]=make_pair(9U,j);     


} // CreateDataFor_ISOPARAMETRIC_LINEAR_PYRAMID()

#else

template<size_t dim>
void FV_IntegrationPointsAndWeights<dim>::CreateDataFor_ISOPARAMETRIC_LINEAR_PYR()
{
    const size_t NumOfInternalFacets(12U);
    const size_t NumOfInternalFacetsPerNode(8U); //Special case of 8 facets at the apex
    const size_t NumOfInternalVolumes(5U); // == VolMultipliers
  	const size_t NumOfIPperVolume(1U);
    const size_t NumOfIPperFacet(1U);

    Resize( NumOfInternalFacets,
            NumOfInternalFacetsPerNode,
            NumOfInternalVolumes,
            NumOfIPperFacet,
            NumOfIPperVolume );
            
    for ( size_t i=0U; i<NumOfInternalVolumes; i++ )	
      m_volume_integration_weights[i][0]=1./4.;
    
    m_volume_integration_weights[4][0]=1./3.;
    
    // Subdivided pyramid volume integration points
    const double vip[5U][3U]=
      {
         {-31./72.,-31./72.,23./144.},
         {31./72.,-31./72.,23./144.},
         {31./72.,31./72.,23./144.},
         {-31./72.,31./72.,23./144.},
         {0.,0.,25./48.}
      }; 
    
    // map to m_volume_integration_points1;   // [ivol][vpts][dim]
    for( size_t i=0U; i<NumOfInternalVolumes; i++ )
      for( size_t j=0U; j<NumOfIPperVolume; j++ )
    	for( size_t k=0U; k<3U; k++ )
    	  m_volume_integration_points1[i][j][k]=vip[i*NumOfIPperVolume+j][k];

	  //face integration points
    const double fip[12U][3U] =
   	 {
  	   {0.,-4./9., 5./36.},
       {4./9.,0.,5./36.},
       {0.,4./9.,5./36.},
       {-4./9.,0.,5./36.},
  	   {-7./18.,-1./6.,13./36.}, 
       {-1./6.,-7./18.,13./36.}, 
       { 1./6.,-7./18.,13./36.}, 
       { 7./18.,-1./6.,13./36.}, 
       { 7./18.,1./6., 13./36.}, 
       { 1./6.,7./18., 13./36.}, 
       {-1./6.,7./18., 13./36.}, 
       {-7./18.,1./6., 13./36.} 
    }; //end of faces integration points initialisation list
    
   for(size_t i=0;i<NumOfInternalFacets;i++)
    for(size_t j=0;j<NumOfIPperFacet; j++)
     for(size_t k=0;k<3U; k++) 
        m_facet_integration_points[i][j][k]=fip[i*NumOfIPperFacet+j][k];    		
 
   // the facet integration weights is the area of the facet in parametric space
   m_facet_integration_weights[0][0]=m_facet_integration_weights[1][0]=m_facet_integration_weights[2][0]=m_facet_integration_weights[3][0]=0.25;
   m_facet_integration_weights[4][0]=m_facet_integration_weights[5][0]=m_facet_integration_weights[6][0]=m_facet_integration_weights[7][0]=sqrt(74.)/48.;
   m_facet_integration_weights[8][0]=m_facet_integration_weights[9][0]=m_facet_integration_weights[10][0]=m_facet_integration_weights[11][0]=sqrt(74.)/48.;
    
   // Projection multipliers from pyramid to unit cube
   m_projection_weights[0][0]=m_projection_weights[1][0]=m_projection_weights[2][0]=m_projection_weights[3][0]=0.3648003127564263;
   m_projection_weights[4][0]=m_projection_weights[5][0]=m_projection_weights[6][0]=m_projection_weights[7][0]=0.2085503283484722;
    
    // Facet normals pointing outside the internal face, all together 12 normals!
    //normal 0 
    m_facet_normals[0][0]=1.0; m_facet_normals[0][1]=0.0; m_facet_normals[0][2]=0.0; 
    //normal 1
    m_facet_normals[1][0]=0.0; m_facet_normals[1][1]=1.0; m_facet_normals[1][2]=0.0;
    //normal 2
    m_facet_normals[2][0]=-1.0; m_facet_normals[2][1]=0.0; m_facet_normals[2][2]=0.0;      
	  //normal 3
    m_facet_normals[3][0]=0.0; m_facet_normals[3][1]=-1.0; m_facet_normals[3][2]=0.0;      
    //normal 4
    m_facet_normals[4][0]=sqrt(74.)/74.; m_facet_normals[4][1]=3*sqrt(74.)/74.; m_facet_normals[4][2]=4*sqrt(74.)/37.;
    //normal 5
    m_facet_normals[5][0]=3*sqrt(74.)/74.; m_facet_normals[5][1]=sqrt(74.)/74.; m_facet_normals[5][2]=4*sqrt(74.)/37.;
    //normal 6
    m_facet_normals[6][0]=-3*sqrt(74.)/74.; m_facet_normals[6][1]=sqrt(74.)/74.; m_facet_normals[6][2]=4*sqrt(74.)/37.;
    //normal 7 
    m_facet_normals[7][0]=-sqrt(74.)/74.; m_facet_normals[7][1]=3*sqrt(74.)/74.; m_facet_normals[7][2]=4*sqrt(74.)/37.;
    //normal 8
    m_facet_normals[8][0]=-sqrt(74.)/74.; m_facet_normals[8][1]=-3*sqrt(74.)/74.; m_facet_normals[8][2]=4*sqrt(74.)/37.;
    //normal 9
    m_facet_normals[9][0]=-3*sqrt(74.)/74.; m_facet_normals[9][1]=-sqrt(74.)/74.; m_facet_normals[9][2]=4*sqrt(74.)/37.;
    //normal 10
    m_facet_normals[10][0]=3*sqrt(74.)/74.; m_facet_normals[10][1]=-sqrt(74.)/74.; m_facet_normals[10][2]=4*sqrt(74.)/37.;
    //normal 11 
    m_facet_normals[11][0]=sqrt(74.)/74. ; m_facet_normals[11][1]=-3*sqrt(74.)/74.; m_facet_normals[11][2]=4*sqrt(74.)/37.;
    
     
    // Facet surrounding node
    m_facets_surrounding_node[0].resize(4U);
    m_facets_surrounding_node[0][0]=0; m_facets_surrounding_node[0][1]=3; m_facets_surrounding_node[0][2]=4; m_facets_surrounding_node[0][3]=5;
    m_facets_surrounding_node[1].resize(4U);
    m_facets_surrounding_node[1][0]=1; m_facets_surrounding_node[1][1]=0; m_facets_surrounding_node[1][2]=6; m_facets_surrounding_node[1][3]=7;
    m_facets_surrounding_node[2].resize(4U);
    m_facets_surrounding_node[2][0]=2; m_facets_surrounding_node[2][1]=1; m_facets_surrounding_node[2][2]=8; m_facets_surrounding_node[2][3]=9;
	  m_facets_surrounding_node[3].resize(4U);
    m_facets_surrounding_node[3][0]=3; m_facets_surrounding_node[3][1]=2; m_facets_surrounding_node[3][2]=10; m_facets_surrounding_node[3][3]=11;
	  m_facets_surrounding_node[4].resize(8U); // apex node
    m_facets_surrounding_node[4][0]=4; m_facets_surrounding_node[4][1]=5; m_facets_surrounding_node[4][2]=6; m_facets_surrounding_node[4][3]=7;
    m_facets_surrounding_node[4][4]=8; m_facets_surrounding_node[4][5]=9; m_facets_surrounding_node[4][6]=10;m_facets_surrounding_node[4][7]=11;

    //Edges   
    m_edges_of_element.resize(12U);
    vector<pair<size_t,size_t> >(m_edges_of_element).swap(m_edges_of_element);
    m_edges_of_element[0]=make_pair(0U,1U);
    m_edges_of_element[1]=make_pair(1U,2U);
    m_edges_of_element[2]=make_pair(2U,3U);
    m_edges_of_element[3]=make_pair(3U,0U);
 
    m_edges_of_element[4]=make_pair(0U,4U);    
    m_edges_of_element[5]=make_pair(0U,4U);
    m_edges_of_element[6]=make_pair(1U,4U);
    m_edges_of_element[7]=make_pair(1U,4U);            
    m_edges_of_element[8]=make_pair(2U,4U);    
    m_edges_of_element[9]=make_pair(2U,4U);
    m_edges_of_element[10]=make_pair(3U,4U);
    m_edges_of_element[11]=make_pair(3U,4U);            
 
     //Facet mid-edge points 
    m_facet_edge_midpoints[0][0] = 0.;     m_facet_edge_midpoints[0][1] = -2./3.; m_facet_edge_midpoints[0][2] = 1./3.; 
    m_facet_edge_midpoints[1][0] = 2./3.;  m_facet_edge_midpoints[1][1] = 0.;     m_facet_edge_midpoints[1][2] = 1./3.; 
    m_facet_edge_midpoints[2][0] = 0.;     m_facet_edge_midpoints[2][1] =  2./3.; m_facet_edge_midpoints[2][2] = 1./3.; 
    m_facet_edge_midpoints[3][0] = -2./3.; m_facet_edge_midpoints[3][1] = 0.;     m_facet_edge_midpoints[3][2] = 1./3.; 
    
    m_facet_edge_midpoints[4] = m_facet_edge_midpoints[3]; 
    m_facet_edge_midpoints[5] = m_facet_edge_midpoints[0]; 
    m_facet_edge_midpoints[6] = m_facet_edge_midpoints[1]; 
    m_facet_edge_midpoints[7] = m_facet_edge_midpoints[2]; 
    
    //m_barycenter 
    m_barycenter[0] = 0.;
    m_barycenter[1] = 0.;
	  m_barycenter[2] = 1./4.; 

    Point<dim> pt_c125, pt_c145, pt_c235, pt_c345, pt_c1234;
	  Point<dim> pt_c12, pt_c23, pt_c34, pt_c14, pt_c15, pt_c25, pt_c35, pt_c45;
	
	  pt_c12[0]=0.;pt_c12[1]=-1.;pt_c12[2]=0.;
	  pt_c23[0]=1.;pt_c23[1]=0.;pt_c23[2]=0.;
	  pt_c34[0]=0.;pt_c34[1]=1.;pt_c34[2]=0.;
	  pt_c14[0]=-1.;pt_c14[1]=0.;pt_c14[2]=0.;
	  
	  pt_c15[0]=-1./2.;pt_c15[1]=-1./2.;pt_c15[2]=1./2.;
	  pt_c25[0]=1./2.;pt_c25[1]=-1./2.;pt_c25[2]=1./2.;
	  pt_c35[0]=1./2.;pt_c35[1]=1./2.;pt_c35[2]=1./2.;
	  pt_c45[0]=-1./2.;pt_c45[1]=1./2.;pt_c45[2]=1./2.;
	
	  pt_c125[0]=0.;pt_c125[1]=-2./3.;pt_c125[2]=1./3.;
	  pt_c145[0]=-2./3.;pt_c145[1]=0.;pt_c145[2]=1./3.;
	  pt_c235[0]=2./3.;pt_c235[1]=0.;pt_c235[2]=1./3.;
	  pt_c345[0]=0.;pt_c345[1]=2./3.;pt_c345[2]=1./3.;
	  pt_c1234[0]=0.;pt_c1234[1]=0.;pt_c1234[2]=0.;

    m_facet_types[0] = QUADRILATERAL_FACET;

	m_facet_points[0].resize(4);
    vector<Point<dim> >(m_facet_points[0]).swap(m_facet_points[0]);
    m_facet_points[0][0] = pt_c12;
    m_facet_points[0][1] = pt_c1234; 
    m_facet_points[0][2] = m_barycenter;
    m_facet_points[0][3] = pt_c125;

    m_facet_types[1] = QUADRILATERAL_FACET;

	m_facet_points[1].resize(4);
    vector<Point<dim> >(m_facet_points[1]).swap(m_facet_points[1]);
    m_facet_points[1][0] = pt_c23;
    m_facet_points[1][1] = pt_c1234; 
    m_facet_points[1][2] = m_barycenter;
    m_facet_points[1][3] = pt_c235;    

    m_facet_types[2] = QUADRILATERAL_FACET;

	m_facet_points[2].resize(4);
    vector<Point<dim> >(m_facet_points[2]).swap(m_facet_points[2]);
    m_facet_points[2][0] = pt_c34;
    m_facet_points[2][1] = pt_c1234; 
    m_facet_points[2][2] = m_barycenter;
    m_facet_points[2][3] = pt_c345;    

    m_facet_types[3] = QUADRILATERAL_FACET;

	m_facet_points[3].resize(4);
    vector<Point<dim> >(m_facet_points[3]).swap(m_facet_points[3]);
    m_facet_points[3][0] = pt_c14;
    m_facet_points[3][1] = pt_c1234; 
    m_facet_points[3][2] = m_barycenter;
    m_facet_points[3][3] = pt_c145;    

    // triangular facets

    m_facet_types[4] = TRIANGULAR_FACET;

	m_facet_points[4].resize(3);
    vector<Point<dim> >(m_facet_points[4]).swap(m_facet_points[4]);
    m_facet_points[4][0] = pt_c145;
    m_facet_points[4][1] = pt_c15; 
    m_facet_points[4][2] = m_barycenter;

    m_facet_types[5] = TRIANGULAR_FACET;

	m_facet_points[5].resize(3);
    vector<Point<dim> >(m_facet_points[5]).swap(m_facet_points[5]);
    m_facet_points[5][0] = pt_c15;
    m_facet_points[5][1] = pt_c125; 
    m_facet_points[5][2] = m_barycenter;

    m_facet_types[6] = TRIANGULAR_FACET;

	m_facet_points[6].resize(3);
    vector<Point<dim> >(m_facet_points[6]).swap(m_facet_points[6]);
    m_facet_points[6][0] = pt_c125;
    m_facet_points[6][1] = pt_c25; 
    m_facet_points[6][2] = m_barycenter;

    m_facet_types[7] = TRIANGULAR_FACET;

	m_facet_points[7].resize(3);
    vector<Point<dim> >(m_facet_points[7]).swap(m_facet_points[7]);
    m_facet_points[7][0] = pt_c25;
    m_facet_points[7][1] = pt_c235; 
    m_facet_points[7][2] = m_barycenter;

    m_facet_types[8] = TRIANGULAR_FACET;

	m_facet_points[8].resize(3);
    vector<Point<dim> >(m_facet_points[8]).swap(m_facet_points[8]);
    m_facet_points[8][0] = pt_c235;
    m_facet_points[8][1] = pt_c35; 
    m_facet_points[8][2] = m_barycenter;

    m_facet_types[9] = TRIANGULAR_FACET;

	m_facet_points[9].resize(3);
    vector<Point<dim> >(m_facet_points[9]).swap(m_facet_points[9]);
    m_facet_points[9][0] = pt_c35;
    m_facet_points[9][1] = pt_c345; 
    m_facet_points[9][2] = m_barycenter;

    m_facet_types[10] = TRIANGULAR_FACET;

	m_facet_points[10].resize(3);
    vector<Point<dim> >(m_facet_points[10]).swap(m_facet_points[10]);
    m_facet_points[10][0] = pt_c345;
    m_facet_points[10][1] = pt_c45; 
    m_facet_points[10][2] = m_barycenter;

    m_facet_types[11] = TRIANGULAR_FACET;

	m_facet_points[11].resize(3);
    vector<Point<dim> >(m_facet_points[11]).swap(m_facet_points[11]);
    m_facet_points[11][0] = pt_c45;
    m_facet_points[11][1] = pt_c145; 
    m_facet_points[11][2] = m_barycenter;


    // tested: 18/Oct/2017 AJB
    const double sqrt2 = 1.41421356237309504880168872420969807856967187537694807317667973799;
    const double facet_normal_transforms[12][2][5] = {
        { { -11.0/96.0, -11.0/96.0, -1.0/32.0, -1.0/32.0, 7.0/24.0 },
            { -19.0/96.0, -19.0/96.0, 7.0/32.0, 7.0/32.0, -1.0/24.0 } },
        { { -1.0/32.0, -11.0/96.0, -11.0/96.0, -1.0/32.0, 7.0/24.0 },
            { 7.0/32.0, -19.0/96.0, -19.0/96.0, 7.0/32.0, -1.0/24.0 } },
        { { -1.0/32.0, -1.0/32.0, -11.0/96.0, -11.0/96.0, 7.0/24.0 },
            { 7.0/32.0, 7.0/32.0, -19.0/96.0, -19.0/96.0, -1.0/24.0 } },
        { { -11.0/96.0, -1.0/32.0, -1.0/32.0, -11.0/96.0, 7.0/24.0 },
            { -19.0/96.0, 7.0/32.0, 7.0/32.0, -19.0/96.0, -1.0/24.0 } },
        { { -7.0*sqrt2/96.0, 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -7.0*sqrt2/96.0, -sqrt2/24.0 },
            { sqrt2/12.0, 0.0, 0.0, -sqrt2/6.0, sqrt2/12.0 } },
        { { -5.0*sqrt2/32.0, 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -sqrt2/8.0 },
            { -sqrt2/12.0, sqrt2/6.0, 0.0, 0.0, -sqrt2/12.0 } },
        { { -7.0*sqrt2/96.0, -7.0*sqrt2/96.0, 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -sqrt2/24.0 },
            { -sqrt2/6.0, sqrt2/12.0, 0.0, 0.0, sqrt2/12.0 } },
        { { 3.0*sqrt2/32.0, -5.0*sqrt2/32.0, 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -sqrt2/8.0 },
            { 0.0, -sqrt2/12.0, sqrt2/6.0, 0.0, -sqrt2/12.0 } },
        { { 3.0*sqrt2/32.0, -7.0*sqrt2/96.0, -7.0*sqrt2/96.0, 3.0*sqrt2/32.0, -sqrt2/24.0 },
            { 0.0, -sqrt2/6.0, sqrt2/12.0, 0.0, sqrt2/12.0 } },
        { { 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -5.0*sqrt2/32.0, 3.0*sqrt2/32.0, -sqrt2/8.0 },
            { 0.0, 0.0, -sqrt2/12.0, sqrt2/6.0, -sqrt2/12.0 } },
        { { 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -7.0*sqrt2/96.0, -7.0*sqrt2/96.0, -sqrt2/24.0 },
            { 0.0, 0.0, -sqrt2/6.0, sqrt2/12.0, sqrt2/12.0 } },
        { { 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, 3.0*sqrt2/32.0, -5.0*sqrt2/32.0, -sqrt2/8.0 },
            { sqrt2/6.0, 0.0, 0.0, -sqrt2/12.0, -sqrt2/12.0 } }
    };
    m_facet_normal_xforms.resize(12u);
    for (size_t f = 0; f < 12u; ++f) {
        m_facet_normal_xforms[f].resize(5u);
        for (size_t n = 0; n < 5u; ++n) {
            m_facet_normal_xforms[f][n] = std::make_pair(
                    facet_normal_transforms[f][0][n],
                    facet_normal_transforms[f][1][n]
                );
        }
    }

    // tested: 07/Nov/2007 Hamid	
    Point<dim> p;
    m_sector_points[0].resize(8U); 
    p[0]=-1.;p[1]=-1.;p[2]=0.;  
    m_sector_points[0][0] = p;
    m_sector_points[0][1] = pt_c12; 
    m_sector_points[0][2] = pt_c1234;
    m_sector_points[0][3] = pt_c14; 
    m_sector_points[0][4] = pt_c15; 
    m_sector_points[0][5] = pt_c125;
    m_sector_points[0][6] = m_barycenter;
    m_sector_points[0][7] = pt_c145;
    
    m_sector_points[1].resize(8U); 
    p[0]=1.;p[1]=-1.;p[2]=0.;  
    m_sector_points[1][0] = p;
    m_sector_points[1][1] = pt_c23; 
    m_sector_points[1][2] = pt_c1234;
    m_sector_points[1][3] = pt_c12; 
    m_sector_points[1][4] = pt_c25; 
    m_sector_points[1][5] = pt_c235 ;
    m_sector_points[1][6] = m_barycenter;
    m_sector_points[1][7] = pt_c125; 
  
    m_sector_points[2].resize(8U); 
    p[0]=1.;p[1]=1.;p[2]=0.;  
    m_sector_points[2][0] = p;
    m_sector_points[2][1] = pt_c34; 
    m_sector_points[2][2] = pt_c1234;
    m_sector_points[2][3] = pt_c23; 
    m_sector_points[2][4] = pt_c35; 
    m_sector_points[2][5] = pt_c345;
    m_sector_points[2][6] = m_barycenter;
    m_sector_points[2][7] = pt_c235; 
    
    m_sector_points[3].resize(8U); 
    p[0]=-1.;p[1]=1.;p[2]=0.;  
    m_sector_points[3][0] = p;
    m_sector_points[3][1] = pt_c14; 
    m_sector_points[3][2] = pt_c1234;
    m_sector_points[3][3] = pt_c34; 
    m_sector_points[3][4] = pt_c45; 
    m_sector_points[3][5] = pt_c145;
    m_sector_points[3][6] = m_barycenter;
    m_sector_points[3][7] = pt_c345;

    m_sector_points[4].resize(10U);
    p[0]=0.;p[1]=0.;p[2]=1.;  
    m_sector_points[4][0] = pt_c15; 
    m_sector_points[4][1] = pt_c125;
    m_sector_points[4][2] = pt_c25; 
    m_sector_points[4][3] = pt_c235;
    m_sector_points[4][4] = pt_c35; 
    m_sector_points[4][5] = pt_c345;
    m_sector_points[4][6] = pt_c45; 
    m_sector_points[4][7] = pt_c145;
    m_sector_points[4][8] = p;
    m_sector_points[4][9] = m_barycenter; 
    
    m_edge_of_sectors.resize(5U); 
    for ( size_t iSector=0U; iSector<4U; iSector++ ) 
      { 
        m_edge_of_sectors[iSector].resize(13U);
        for ( size_t iLine=0U; iLine<4U; iLine++ )
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine+1)%4);  //   [iSector][iLine] 
        for ( size_t iLine=4U; iLine<8U; iLine++ )
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine,(iLine-3)%4+4);  
        for ( size_t iLine=8U; iLine<12U; iLine++ ) 
          m_edge_of_sectors[iSector][iLine]=make_pair(iLine-8,iLine-4);  
        m_edge_of_sectors[iSector][12]=make_pair(4U,6U);
      }  
    m_edge_of_sectors[4].resize(20U);
    for ( size_t iLine=0U; iLine<8U; iLine++ )
      m_edge_of_sectors[4][iLine]=make_pair(iLine,(iLine+1)%8);  //   [iSector][iLine] 
    size_t j=0U;
    for ( size_t iLine=8U; iLine<12U; iLine++ ,j++,j++ )
      m_edge_of_sectors[4][iLine]=make_pair(8U,j);  
    j=0U;
    for ( size_t iLine=12U; iLine<20U; iLine++,j++ ) 
      m_edge_of_sectors[4][iLine]=make_pair(9U,j);  

}

#endif
 
template class FV_IntegrationPointsAndWeights<1U>;
template class FV_IntegrationPointsAndWeights<2U>;
template class FV_IntegrationPointsAndWeights<3U>;


} // end namespace csp
