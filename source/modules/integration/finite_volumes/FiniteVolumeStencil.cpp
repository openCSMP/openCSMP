// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FiniteVolumeStencil.h"
#include "FiniteElement.h"
#include "FV_IntegrationPointsAndWeights.h"

using namespace std;

namespace csmp { 
 


/**
 
 Parametrised constructor for the FiniteVolumeStencil classpoints.

 const FiniteVolumeStencil<dim>& fvs - reference to the FiniteVolumeStencil
 object.  

@section application Application 

 Standard as for any parametrised constructor. Method parametr should be previously
 initialised as well.
*/
template<uint32_t dim>
FiniteVolumeStencil<dim>::FiniteVolumeStencil( const FiniteVolumeStencil<dim>& fvs )
 : facet_integration_weights_(fvs.facet_integration_weights_),
   facet_projection_weights_(fvs.facet_projection_weights_),
   facet_normals_(fvs.facet_normals_),
   facet_parametric_normals_(fvs.facet_parametric_normals_),
   facet_edge_midpoints_(fvs.facet_edge_midpoints_),
   sector_integration_weights_(fvs.sector_integration_weights_),
   facet_normal_xforms_(fvs.facet_normal_xforms_),
   facet_points_(fvs.facet_points_),
   facet_integration_points_(fvs.facet_integration_points_),
   sector_points_(fvs.sector_points_),
   sector_integration_points_(fvs.sector_integration_points_),
   facets_surrounding_node_(fvs.facets_surrounding_node_),
   edges_of_element_(fvs.edges_of_element_),
   edges_of_sectors_(fvs.edges_of_sectors_),
   facet_types_(fvs.facet_types_),
   barycenter_(fvs.barycenter_),
   parent_element_(fvs.parent_element_),
   parent_element_type_(fvs.parent_element_type_),
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
template<uint32_t dim>
FiniteVolumeStencil<dim>&  FiniteVolumeStencil<dim>::operator=( const FiniteVolumeStencil<dim>& fvs )
 {
     if ( &fvs != this ) {
          edges_of_element_           = fvs.edges_of_element_;
          facets_surrounding_node_    = fvs.facets_surrounding_node_;   // [node][facet]
          facet_integration_points_   = fvs.facet_integration_points_;  // [isrf][spts][dim]
          facet_integration_weights_  = fvs.facet_integration_weights_; // [isrf][spts]
          facet_projection_weights_   = fvs.facet_projection_weights_;
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
     	    edges_of_sectors_           = fvs.edges_of_sectors_;
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
void FiniteVolumeStencil<dim>::Resize( uint32_t n_isrf,
                                       uint32_t srfs_per_node,
                                       uint32_t n_ivol,
                                       uint32_t n_spts,
                                       uint32_t n_vpts )
 {
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
    sector_integration_points_.resize( n_ivol, n_vpts, dim );
    sector_integration_weights_.resize( n_ivol, n_vpts );      // [isrf][spts]
    sector_points_.resize( n_ivol );
    facet_integration_points_.resize( n_isrf, n_spts, dim );
    facet_integration_weights_.resize( n_isrf, n_spts );
    facet_normal_xforms_.resize( n_isrf );                    // [isrf][node][2] (was pair)
    facet_projection_weights_.resize( n_isrf, n_spts );       // [isrf][spts]
    facets_surrounding_node_.resize( n_ivol );                // [node][facet], ivol = nodes
    facet_points_.resize( n_isrf );
#else
    sector_integration_points_.resize( n_ivol, dim );
    sector_integration_weights_.resize( n_ivol );      // [isrf][spts]
    sector_points_.resize( n_ivol );
    facet_integration_points_.resize( n_isrf, dim );
    facet_integration_weights_.resize( n_isrf );
    facet_normal_xforms_.resize( n_isrf );        // [isrf][node][2] (was pair)
    facet_projection_weights_.resize( n_isrf );       // [isrf]
    facets_surrounding_node_.resize( n_ivol ); // [node][facet], ivol = nodes
    facet_points_.resize( n_isrf );
#endif
    facet_normals_.resize( n_isrf, dim );
    facet_parametric_normals_.resize(  n_isrf, dim );
    facet_edge_midpoints_.resize( n_isrf, dim );
    facet_types_.resize( n_isrf );
    edges_of_element_.resize( n_isrf, 2u ); // 2 edges per surface (formerly pair)

 } // end Resize



/**
 
 Parametrised constructor for the FiniteVolumeStencil class, setting up
 stencil for given CSMP type of the FE. 

@param csp_finite_element_type is the @enum CSMP_FEM_TYPE

 const char* csp_finite_element_type - pointer to the character array of the FE type.  

@section application Application 

 Standard as for any parametrised constructor. Method parameter should be CSP
 FE type.
*/
template<uint32_t dim>
FiniteVolumeStencil<dim>::FiniteVolumeStencil( const char* csp_finite_element_type )
 : parent_element_("not initialized"), parent_element_type_(UNKNOWN)
 {
    Initialize( csp_finite_element_type );
 }




/**
Temporary storage used to convert vectors of vectors to the dynamic array based storage
*/
template<uint32_t dim>
struct TempVecs {
    std::vector<std::vector<Point<dim> > >      facet_integration_points;     ///< [isrf][spts][dim]
    std::vector<std::vector<double> >           facet_integration_weights;    ///< [isrf][spts]
    std::vector<Point<dim> >                    facet_normals;                ///< [isrf][dim] //[node*3][dim] -3d
    std::vector<std::vector<std::pair<double,double>>> facet_normal_xforms; ///< [isrf][node]
    std::vector<Point<dim> >                    facet_parametric_normals;     ///< [isrf][dim] //[node*3][dim] -3d
    std::vector<std::vector<double> >           facet_projection_weights;     ///< [isrf][spts]
    std::vector<std::vector<uint32_t> >         facets_surrounding_node;      ///< [node][facet]
    std::vector<Point<dim> >                    facet_edge_midpoints;
    std::vector<std::vector<Point<dim> > >      facet_points;
    std::vector<FV_FACET_TYPE>                  facet_types;                  ///< [isrf]
    std::vector<std::pair<uint32_t,uint32_t> >  edges_of_element;             ///< = facets
    std::vector<std::vector<Point<dim> > >      sector_integration_points;    ///< [ivol][vpts][dim]
    std::vector<std::vector<double> >           sector_integration_weights;   ///< [ivol][vpts]
    std::vector<std::vector<Point<dim> > >      sector_points;
    std::vector<std::vector<std::pair<uint32_t,uint32_t> > >  sector_edges;
};

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
void FiniteVolumeStencil<dim>::ResizeTemporaryVectors( TempVecs<dim>& tvecs,
                                                       uint32_t n_isrf,
                                                       uint32_t srfs_per_node,
                                                       uint32_t n_ivol,
                                                       uint32_t n_spts,
                                                       uint32_t n_vpts )
 {
    assert( n_spts == 1U ); // change projection_weights, facet_integration_weights & volume_integration_weights into 2D arrays if
    assert( n_vpts == 1U ); // there are more than 1 integration points per facet and sector
    tvecs.facet_integration_points.resize( n_isrf );
    tvecs.facet_integration_weights.resize( n_isrf );
    tvecs.facet_normals.resize( n_isrf );
    tvecs.facet_normal_xforms.resize( n_isrf );
    tvecs.facet_parametric_normals.resize( n_isrf );
    tvecs.facet_projection_weights.resize( n_isrf );
    tvecs.facets_surrounding_node.resize( n_ivol ); // ivol = nodes
    tvecs.sector_edges.resize( n_ivol ); // ivol = nodes
    tvecs.facet_edge_midpoints.resize( n_isrf );
    tvecs.facet_points.resize( n_isrf );
    tvecs.facet_types.resize( n_isrf );
    tvecs.edges_of_element.resize(n_isrf);

    for ( uint32_t i{0U}; i<n_ivol; i++ ) {
         tvecs.facets_surrounding_node[i].resize( srfs_per_node );
         tvecs.sector_edges[i].resize( srfs_per_node ); // no further resize, since vec[vec[pairs]]]
      }
    
    for ( uint32_t i{0U}; i<n_isrf; i++ ) {
         tvecs.facet_integration_points[i].resize( n_spts );
         tvecs.facet_projection_weights[i].resize( n_spts );
         tvecs.facet_integration_weights[i].resize( n_spts );
         tvecs.facet_points[i].resize( n_spts );
         tvecs.facet_normal_xforms[i].resize( n_ivol ); // ivol = nodes
      }
    
    tvecs.sector_integration_points.resize( n_ivol );
    tvecs.sector_integration_weights.resize( n_ivol );
    //par_sector_integration_points.resize( n_ivol );

    for ( uint32_t i{0U}; i<n_ivol; i++ ) {
         tvecs.sector_integration_points[i].resize( n_vpts );
         //par_sector_integration_points[i].resize( n_vpts );
           //par_sector_integration_points[i][j].resize(pdim);
         tvecs.sector_integration_weights[i].resize(n_vpts);
      }
      
    tvecs.sector_points.resize(n_ivol);

 } // end ResizeTemporaryVectors





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
void FiniteVolumeStencil<dim>::Initialize( const char* csp_finite_element_type )
 {
   // uses the types specified in the FiniteElement.h header file
   const string csp_fem_type(csp_finite_element_type);
   bool debug=false;

   uint32_t NumOfInternalFacets=6U;
   uint32_t NumOfInternalFacetsPerNode=3U;
   uint32_t NumOfInternalVolumes=4U;
   uint32_t NumOfIPperVolume=1U;
   uint32_t NumOfIPperFacet=1U;
   
   TempVecs<dim> temp_vecs;
   
   if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_BAR) )
     {
          parent_element_      = "ISOPARAMETRIC_LINEAR_BAR";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_BAR;
          space_dimension_     = LINE;

          NumOfInternalFacets=1U;
          NumOfInternalFacetsPerNode=1U;
          NumOfInternalVolumes=2U;

          FV_IntegrationPointsAndWeights<dim> FV(ISOPARAMETRIC_LINEAR_BAR);

          // containers in FV_IntegrationPointsAndWeights
          ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode,
                                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );
          // local ones
          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          InitializeDataMembers( temp_vecs, FV );

          if (debug) {
               cout<<"\nISOPARAMETRIC_LINEAR_BAR:" <<endl;
               Out();
             }
          return;
      }
      
   if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_TRIANGLE) ) {

          parent_element_ = "ISOPARAMETRIC_LINEAR_TRIANGLE";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_TRIANGLE;
          space_dimension_ = SURFACE;
     
          NumOfInternalFacets=3U;
          NumOfInternalFacetsPerNode=2U;
          NumOfInternalVolumes=3U;

          FV_IntegrationPointsAndWeights<dim> FV(ISOPARAMETRIC_LINEAR_TRIANGLE);

          ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode,
                                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          InitializeDataMembers( temp_vecs, FV );
             
          if (debug) {
                cout<<"\nISOPARAMETRIC_LINEAR_TRIANGLE:" <<endl;
                Out();
            }
          return;
      }
      
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_TETRAHEDRON) ) {

         parent_element_ = "ISOPARAMETRIC_LINEAR_TETRAHEDRON";
         parent_element_type_ = ISOPARAMETRIC_LINEAR_TETRAHEDRON;
         space_dimension_ = VOLUME;

         NumOfInternalFacets=6U;
         NumOfInternalFacetsPerNode=3U;
         NumOfInternalVolumes=4U;

         FV_IntegrationPointsAndWeights<dim> FVT(ISOPARAMETRIC_LINEAR_TETRAHEDRON);

         ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode,
                                 NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

         InitializeDataMembers( temp_vecs, FVT );
             
         if (debug) {
                cout<<"\nISOPARAMETRIC_LINEAR_TETRAHEDRON: "<<endl;
                Out();         
             }
          return;
     }
     
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_QUADRILATERAL) ) {

          parent_element_      = "ISOPARAMETRIC_LINEAR_QUADRILATERAL";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_QUADRILATERAL;
          space_dimension_     = SURFACE;

          NumOfInternalFacets=4U;
          NumOfInternalFacetsPerNode=2U;
          NumOfInternalVolumes=4U; // == VolMultipliers

          FV_IntegrationPointsAndWeights<dim> FVQ(ISOPARAMETRIC_LINEAR_QUADRILATERAL);

          ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode,
                                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );
                 
          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          InitializeDataMembers( temp_vecs, FVQ );

          if (debug) {
                cout<<"\nISOPARAMETRIC_LINEAR_QUADRILATERAL:"<<endl;
                Out();
             }
          return;
      }
      
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_HEXAHEDRON) ) {

          parent_element_      = "ISOPARAMETRIC_LINEAR_HEXAHEDRON";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_HEXAHEDRON;
          space_dimension_     = VOLUME;

          NumOfInternalFacets=12U;
          NumOfInternalFacetsPerNode=3U;
          NumOfInternalVolumes=8U;

          FV_IntegrationPointsAndWeights<dim> FVH(ISOPARAMETRIC_LINEAR_HEXAHEDRON);

          ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode, 
                                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );
                 
          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          InitializeDataMembers( temp_vecs, FVH );

         if (debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_HEXAHEDRON:"<<endl;
              Out();         
          }
         return;
      }
      
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_PYRAMID) ) {

          parent_element_ = "ISOPARAMETRIC_LINEAR_PYRAMID";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_PYRAMID;
          space_dimension_ = VOLUME;

 #ifdef PYRAMID_TRIANGULAR_FACETS
          NumOfInternalFacets=12U;
          NumOfInternalFacetsPerNode=8U; // Special case of the Pyramid vertex require 8
          NumOfInternalVolumes=5U;       // == VolMultipliers
 #else
          NumOfInternalFacets=8U;
          NumOfInternalFacetsPerNode=4U; // Special case of the Pyramid vertex require 4
          NumOfInternalVolumes=5U;       // == VolMultipliers
 #endif

          FV_IntegrationPointsAndWeights<dim> FVP(ISOPARAMETRIC_LINEAR_PYRAMID);

          ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode,
                                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );
                 
          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          InitializeDataMembers( temp_vecs, FVP );

         if ( debug ) {
              cout<<"\nISOPARAMETRIC_LINEAR_PYRAMID:"<<endl;
              Out();         
          }
         return;
      }
    
    if ( csp_fem_type == parseFiniteElementType(ISOPARAMETRIC_LINEAR_PRISM) ) {

          parent_element_      = "ISOPARAMETRIC_LINEAR_PRISM";
          parent_element_type_ = ISOPARAMETRIC_LINEAR_PRISM;
          space_dimension_     = VOLUME;

          NumOfInternalFacets=9U;
          NumOfInternalFacetsPerNode=3U;
          NumOfInternalVolumes=6U;

          FV_IntegrationPointsAndWeights<dim> FVP(ISOPARAMETRIC_LINEAR_PRISM);

          ResizeTemporaryVectors( temp_vecs, NumOfInternalFacets, NumOfInternalFacetsPerNode,
                                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );
                 
          Resize( NumOfInternalFacets, NumOfInternalFacetsPerNode,
                  NumOfInternalVolumes, NumOfIPperFacet, NumOfIPperVolume );

          InitializeDataMembers( temp_vecs, FVP );
		  
          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_PRISM:"<<endl;
              Out();         
           }
          return;
      }
    else {
         cerr <<"\nFiniteVolumeStencil::Initialize: Element type not recognized: ";
         cerr << csp_finite_element_type << endl;
      }
      
 } // end Initialize





template<uint32_t dim>
void FiniteVolumeStencil<dim>::InitializeDataMembers( TempVecs<dim>& tvecs,
                                                      const FV_IntegrationPointsAndWeights<dim>& FV )
 {
    // TODO: why is fictitious dimension needed for line and surface elements?
    const uint32_t parametric_space_dimension = dim; // static_cast<int8_t>(space_dimension_);
    
    // facet integration weights
    {
       FV.FacetIntegrationWeights(tvecs.facet_integration_weights);
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
       facet_integration_weights_.resize( tvecs.facet_integration_weights.size(),
                                          tvecs.facet_integration_weights[0].size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_integration_weights ) {
            uint32_t j{0u};
            for ( const auto& jit : iit )
               facet_integration_weights_(i,j++) = jit;
            i++;
         }
#else
       facet_integration_weights_.resize( tvecs.facet_integration_weights.size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_integration_weights )
         for ( const auto& jit : iit )
           facet_integration_weights_[i++] = jit;
#endif
    }

    // facet projection weights
    {
       FV.ProjectionWeights(tvecs.facet_projection_weights);
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
       facet_projection_weights_.resize( tvecs.facet_projection_weights.size(),
                                         tvecs.facet_projection_weights[0].size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_projection_weights ) {
            uint32_t j{0u};
            for ( const auto& jit : iit )
               facet_projection_weights_(i,j++) = jit;
            i++;
         }
#else
       facet_projection_weights_.resize( tvecs.facet_projection_weights.size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_projection_weights )
         for ( const auto& jit : iit )
           facet_projection_weights_[i++] = jit;
#endif
    }

    // facet normals
    {
       FV.FacetNormals(tvecs.facet_normals);
       facet_normals_.resize( tvecs.facet_normals.size(), parametric_space_dimension );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_normals ) {
           for ( uint32_t j{0u}; j<parametric_space_dimension; ++j )
             facet_normals_.at(i,j) = iit[j];
           ++i;
        }
    }

    // facet parametric normals
    {
       FV.FacetNormals(tvecs.facet_parametric_normals);
       facet_parametric_normals_.resize( tvecs.facet_parametric_normals.size(), parametric_space_dimension );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_parametric_normals ) {
           for ( uint32_t j{0u}; j<parametric_space_dimension; ++j )
             facet_parametric_normals_.at(i,j) = iit[j];
           ++i;
        }
    }

    // facet edge midpoints
    {
       FV.EdgeMidpoints(tvecs.facet_edge_midpoints);
       facet_edge_midpoints_.resize( tvecs.facet_edge_midpoints.size(), parametric_space_dimension );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_edge_midpoints ) {
           for ( uint32_t j{0u}; j<parametric_space_dimension; ++j )
             facet_edge_midpoints_.at(i,j) = iit[j];
           ++i;
        }
    }

    // sector integration weights
    {
       FV.SectorIntegrationWeights(tvecs.sector_integration_weights);
#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
       sector_integration_weights_.resize( tvecs.sector_integration_weights.size(),
                                           tvecs.sector_integration_weights[0].size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.sector_integration_weights ) {
            uint32_t j{0u};
            for ( const auto& jit : iit )
               sector_integration_weights_(i,j++) = jit;
            i++;
         }
#else
       sector_integration_weights_.resize( tvecs.sector_integration_weights.size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.sector_integration_weights )
         for ( const auto& jit : iit )
           sector_integration_weights_[i++] = jit;
#endif
    }
    
    
  // THREE-DIMENSIONAL ARRAYS
   
   // facet transformations (parametric to physical)
   {
      FV.FacetNormalTransformations(tvecs.facet_normal_xforms);
      facet_normal_xforms_ = tvecs.facet_normal_xforms;
      /*
       facet_normal_xforms_.resize( tvecs.facet_normal_xforms.size(),
                                    tvecs.facet_normal_xforms[0].size(), 2u ); // pair
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facet_normal_xforms ) {
           uint32_t j{0u};
           for ( const auto& jit : iit ) {
                 facet_normal_xforms_.at(i,j,0u) = jit.first;
                 facet_normal_xforms_.at(i,j,1u) = jit.second;
               ++j;
             }
           ++i;
        }
      */
   }
    
   // facet points
   {
      FV.FacetPoints(tvecs.facet_points);
      facet_points_ = tvecs.facet_points;
      /*
      facet_points_.resize( tvecs.facet_points.size(), tvecs.facet_points[0].size(), parametric_space_dimension );
      uint32_t i{0u};
      for ( const auto& iit : tvecs.facet_points ) {
          uint32_t j{0u};
          for ( const auto& jit : iit ) {
              for ( uint32_t k{0u}; k<parametric_space_dimension; ++k )
                facet_points_.at(i,j,k) = jit[j];
              ++j;
            }
          ++i;
       }
      */
    }

   // sector points
   {
      FV.SectorPoints(tvecs.sector_points);
      sector_points_ = tvecs.sector_points;
      /*
      sector_points_.resize( tvecs.sector_points.size(), tvecs.sector_points[0].size(), parametric_space_dimension );
      uint32_t i{0u};
      for ( const auto& iit : tvecs.sector_points ) {
          uint32_t j{0u};
          for ( const auto& jit : iit ) {
              for ( uint32_t k{0u}; k<parametric_space_dimension; ++k )
                sector_points_.at(i,j,k) = jit[j];
              ++j;
            }
          ++i;
       }
      */
    }

#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
   // facet integration points
   {
      FV.FacetIntegrationPoints(tvecs.facet_integration_points);
      facet_integration_points_.resize( tvecs.facet_integration_points.size(),
                                        tvecs.facet_integration_points[0].size(), parametric_space_dimension );

      for ( uint32_t i{0u}; i<tvecs.facet_integration_points.size(); ++i )
        for ( uint32_t j{0u}; j<tvecs.facet_integration_points[i].size(); ++j )
          for ( uint32_t k{0u}; k<parametric_space_dimension; ++k )
            facet_integration_points_.at(i,j,k) = tvecs.facet_integration_points[i][j][k];
    }

   // sector integration points
   {
      FV.SectorIntegrationPoints(tvecs.sector_integration_points);
      sector_integration_points_.resize( tvecs.sector_integration_points.size(),
                                         tvecs.sector_integration_points[0].size(), parametric_space_dimension );

      for ( uint32_t i{0u}; i<tvecs.sector_integration_points.size(); ++i )
        for ( uint32_t j{0u}; j<tvecs.sector_integration_points[i].size(); ++j )
          for ( uint32_t k{0u}; k<parametric_space_dimension; ++k )
            sector_integration_points_.at(i,j,k) = tvecs.sector_integration_points[i][j][k];
    }
#else // single integration point
   {
      FV.FacetIntegrationPoints(tvecs.facet_integration_points);
      facet_integration_points_.resize( tvecs.facet_integration_points.size(), parametric_space_dimension );

      for ( uint32_t i{0u}; i<tvecs.facet_integration_points.size(); ++i )
        for ( uint32_t j{0u}; j<tvecs.facet_integration_points[i].size(); ++j )
          for ( uint32_t k{0u}; k<parametric_space_dimension; ++k )
            facet_integration_points_.at(i,k) = tvecs.facet_integration_points[i][j][k];
    }
   {
      FV.SectorIntegrationPoints(tvecs.sector_integration_points);
      sector_integration_points_.resize( tvecs.sector_integration_points.size(), parametric_space_dimension );

      for ( uint32_t i{0u}; i<tvecs.sector_integration_points.size(); ++i )
        for ( uint32_t j{0u}; j<tvecs.sector_integration_points[i].size(); ++j )
          for ( uint32_t k{0u}; k<parametric_space_dimension; ++k )
            sector_integration_points_.at(i,k) = tvecs.sector_integration_points[i][j][k];
    }
#endif

    // edges of elements
    {
       FV.EdgePairs(tvecs.edges_of_element);
       edges_of_element_.resize( tvecs.edges_of_element.size(), 2u ); // pair
       uint32_t i{0u};
       for ( const auto& iit : tvecs.edges_of_element ) {
             edges_of_element_.at(i,0u) = iit.first;
             edges_of_element_.at(i,1u) = iit.second;
           ++i;
        }
    }

    // facets surrounding node
    {
       // jagged array
       FV.FacetsSurroundingNode(tvecs.facets_surrounding_node);
       facets_surrounding_node_ = tvecs.facets_surrounding_node;
       /*
       facets_surrounding_node_.resize( tvecs.facets_surrounding_node.size(),
                                        tvecs.facets_surrounding_node[0].size() );
       uint32_t i{0u};
       for ( const auto& iit : tvecs.facets_surrounding_node ) {
           uint32_t j{0u};
           for ( const auto& jit : iit ) {
                facets_surrounding_node_.at(i,j) = jit;
                ++j;
             }
           ++i;
        }
       */
    }

   // sector edges
   {
      FV.SectorEdgePairs(tvecs.sector_edges);
      // jagged array!
      edges_of_sectors_ = tvecs.sector_edges;
      /*
      edges_of_sectors_.resize( tvecs.sector_edges.size() );
      for ( uint32_t i{0u}; i<tvecs.sector_edges[i].size(); ++i )
        edges_of_sectors_[i].resize(tvecs.sector_edges[i].size());
        
      uint32_t i{0u};
      for ( const auto& iit : tvecs.sector_edges ) {
          uint32_t j{0u};
          for ( const auto& jit : iit ) {
               edges_of_sectors_[i][j++] = jit;
            }
          ++i;
       }
      */
    }

    FV.Barycenter(barycenter_);
    FV.FacetTypes(facet_types_);
          
 } // InitializeDataMembers






// SKM modified version TODO: check assumption that n-columns in these arrays is the same for all rows
template<uint32_t dim>
void FiniteVolumeStencil<dim>::Out() const
{
   cout<<"\nFiniteVolumeStencil<dim>::Out: \nsector integration points"; // 2d array
   cout<<endl<<endl<<" Volume IPs "<<endl; // 3d array

#ifdef FV_STENCIL_WITH_MORE_THAN_1_INTEGRATION_POINT_PER_FACET_OR_SECTOR
   // sectors
   cout <<"("<< sector_integration_points_.depth() <<","<< sector_integration_points_.rows() <<", rst):";
   for( uint32_t i{0U}; i<sector_integration_points_.depth(); i++)
     for( uint32_t j{0U}; j<sector_integration_points_.rows(); j++) {
             for( uint32_t k{0U}; k<dim; k++) cout<<" "<< sector_integration_points_(i,j,k);
             cout<<endl;
          }
   cout<<"\n Volume weights (sectors)"; // 2d array
   cout <<"("<< sector_integration_weights_.rows() <<","<< sector_integration_weights_.cols() <<"):";
   for( uint32_t i{0U}; i<sector_integration_weights_.rows(); i++)
      for( uint32_t j{0U}; j<sector_integration_weights_.cols(); j++)
        cout<<" "<< sector_integration_weights_(i,j);
   
   //facets
   cout<<endl<<" Facet IPs: "<<endl; // 3d array
   cout <<"("<< facet_integration_points_.depth() <<","<< facet_integration_points_.rows() <<",rst):";
   for( uint32_t i{0U}; i<facet_integration_points_.depth(); i++)
     for( uint32_t j{0U}; j<facet_integration_points_.rows(); j++) {
          for( uint32_t k{0U}; k<dim; k++) cout<<" "<< facet_integration_points_(i,j,k);
          cout<<endl;
       }

   cout<<endl<<" Facet weights (facets) "; // 2d array
   cout <<"("<< facet_integration_weights_.rows() <<","<< facet_integration_weights_.cols() <<"):";
   for( uint32_t i{0U}; i<facet_integration_weights_.rows(); i++)
     for( uint32_t j{0U}; j<facet_integration_weights_.cols(); j++) {
         cout<<" "<< facet_integration_weights_(i,j);
       }
   cout<<endl;

   cout<<endl<<" Facet projection weights (facets) "; // 2d array
   cout <<"("<< facet_projection_weights_.rows() <<","<< facet_projection_weights_.cols() <<"):";
   for( uint32_t i{0U}; i<facet_projection_weights_.size(); i++)
     for( uint32_t j{0U}; j<facet_integration_weights_.cols(); j++) {
         cout<<" "<< facet_projection_weights_(i,j);
      }
   cout<<endl;

#else
   // sectors
   // ---------------------------
   // integration points 2D array
   cout <<"("<< sector_integration_points_.rows() <<","<< sector_integration_points_.cols() <<",rst):";
   for( uint32_t i{0U}; i<sector_integration_points_.rows(); i++ )
     for( uint32_t j{0U}; j<sector_integration_points_.cols(); j++ ) {
             cout<<" "<< sector_integration_points_(i,j);
             cout<<endl;
          }
   // 1D vector
   cout <<"("<< sector_integration_weights_.size() <<"):";
   for( uint32_t i{0U}; i<sector_integration_weights_.size(); i++)
     cout<<" "<< sector_integration_weights_[i];
   
   // facets
   // --------------------------
   cout<<endl<<" Facet IPs: "<<endl; // 2d array
   cout <<"("<< facet_integration_points_.rows() <<","<< facet_integration_points_.cols() <<"=rst):";
   for( uint32_t i{0U}; i<facet_integration_points_.rows(); i++) {
     for( uint32_t j{0U}; j<facet_integration_points_.cols(); j++) {
         cout<<" "<< facet_integration_points_(i,j);
         cout<<endl;
       }
     }

   // facet (integration) weights 1D vector
   cout<<endl<<" Facet weights (facets) "; // 1d array
   cout <<"("<< facet_integration_weights_.size() <<"):";
   for( uint32_t i{0U}; i<facet_integration_weights_.size(); i++) {
        cout<<" "<< facet_integration_weights_[i];
     }
   cout<<endl;

   cout<<endl<<" Facet projection weights (facets) "; // 1d vector
   cout <<"("<< facet_projection_weights_.size() <<"):";
   for( uint32_t i{0U}; i<facet_projection_weights_.size(); i++) {
        cout<<" "<< facet_projection_weights_[i];
     }
   cout<<endl;

#endif

   cout<<endl<<" Facet parametric normals inside the element (vector) "<<endl; // 2d array
   cout <<"("<< facet_parametric_normals_.rows() <<",rst):";
   for( uint32_t i{0U}; i<facet_parametric_normals_.rows(); i++) {
       for( uint32_t j{0U}; j<dim; j++) cout<<" "<< facet_parametric_normals_(i,j);
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

   cout<<endl<<" Facets surrounding node (corresponding to sector) in the element "<<endl; // 2d array
   cout <<"("<< Sectors() <<","<< dim <<"):";
   for( uint32_t i{0U}; i<Sectors(); i++) {
     for( uint32_t j{0U}; j<FacetsPerSector(i); j++) {
       cout<<" "<< FacetSurroundingSector(i,j);
     }
     cout<<endl;
   }
  
   cout<<endl<<" Edge pairs in the element "<<endl; // 2d array
   cout <<"("<< edges_of_element_.rows() <<",2):";
   for( uint32_t i{0U}; i<edges_of_element_.rows(); i++) {
     cout<<" "<< edges_of_element_(i,0U) <<", "<< edges_of_element_(i,1U) <<" ";
   }
   cout << endl << endl;
   
} // end Out


template class FiniteVolumeStencil<1U>; 
template class FiniteVolumeStencil<2U>; 
template class FiniteVolumeStencil<3U>; 


} // end namespace csmp
