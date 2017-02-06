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
 
 Parametrised constructor for the FiniteVolumeStencil class.

 const FiniteVolumeStencil<dim>& fvs - reference to the FiniteVolumeStencil
 object.  

@section application Application 

 Standard as for any parametrised constructor. Method parametr should be previously
 initialised as well.
*/
template<size_t dim>
FiniteVolumeStencil<dim>::FiniteVolumeStencil( const FiniteVolumeStencil<dim>& fvs )
 {
    *this = fvs;
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
          facet_parametric_normals   = fvs.facet_parametric_normals;  // [isrf][dim]
          sector_integration_points  = fvs.sector_integration_points;   // [ivol][vpts][dim]
          sector_integration_weights = fvs.sector_integration_weights;  // [ivol][vpts]
          parent_element_            = fvs.parent_element_;
          facet_edge_midpoints       = fvs.facet_edge_midpoints;
     	    barycenter                 = fvs.barycenter;
    	    facet_points  		         = fvs.facet_points;
    	    sector_points_             = fvs.sector_points_;
     	    sector_edges_              = fvs.sector_edges_;
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
    facet_parametric_normals.resize( n_isrf );
    facet_projection_weights.resize( n_isrf );    
    facets_surrounding_node.resize( n_ivol ); // ivol = nodes
    facet_edge_midpoints.resize( n_isrf );
    facet_points.resize( n_isrf );
    edges_of_element.resize(n_isrf);

    for ( size_t i=0; i<n_ivol; i++ ) 
      facets_surrounding_node[i].resize( srfs_per_node );
    
    for ( size_t i=0; i<n_isrf; i++ ) {
         facet_integration_points[i].resize( n_spts );
         facet_projection_weights[i].resize( n_spts );
         facet_integration_weights[i].resize( n_spts );
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
             
          FV.FacetsSurroundingNode(facets_surrounding_node);
             
          FV.EdgePairs(edges_of_element);

		      FV.EdgeMidpoints(facet_edge_midpoints);
		  
		      FV.Barycenter(barycenter);
    
    	    FV.FacetPoints(facet_points);
    	      	  
    	    FV.SectorPoints(sector_points_);
    	    
     	    FV.SectorEdgePairs(sector_edges_); 

          if(debug) {
             cout<<"\nISOPARAMETRIC_LINEAR_BAR:" <<endl;
             Out(cout);         
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
             
          FV.FacetsSurroundingNode(facets_surrounding_node);
             
          FV.EdgePairs(edges_of_element);

		      FV.EdgeMidpoints(facet_edge_midpoints);
		  
		      FV.Barycenter(barycenter);

		      FV.FacetPoints(facet_points);
    	  
    	    FV.SectorPoints(sector_points_);
    	    
     	    FV.SectorEdgePairs(sector_edges_); 
    	  
          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_TRIANGLE:" <<endl;
              Out(cout);         
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
             
             FVT.FacetsSurroundingNode(facets_surrounding_node);
             
             FVT.EdgePairs(edges_of_element);
             
             FVT.ProjectionWeights(facet_projection_weights);
             
     		     FVT.EdgeMidpoints(facet_edge_midpoints);
		  
             FVT.Barycenter(barycenter);
             
             FVT.FacetPoints(facet_points);
             
    	       FVT.SectorPoints(sector_points_);
    	    
        	   FVT.SectorEdgePairs(sector_edges_); 
             
             if(debug) {
                cout<<"\nISOPARAMETRIC_LINEAR_TETRAHEDRON: "<<endl;
                Out(cout);         
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
             
          FVQ.FacetsSurroundingNode(facets_surrounding_node);
            
          FVQ.EdgePairs(edges_of_element);

		      FVQ.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVQ.Barycenter(barycenter);

          FVQ.FacetPoints(facet_points);

    	    FVQ.SectorPoints(sector_points_);
    	    
     	    FVQ.SectorEdgePairs(sector_edges_); 

          if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_QUADRILATERAL:"<<endl;
              Out(cout);         
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
             
          FVH.FacetsSurroundingNode(facets_surrounding_node);

          FVH.EdgePairs(edges_of_element);
          
          FVH.ProjectionWeights(facet_projection_weights);

		      FVH.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVH.Barycenter(barycenter);
          
          FVH.FacetPoints(facet_points);

    	    FVH.SectorPoints(sector_points_);
    	    
     	    FVH.SectorEdgePairs(sector_edges_); 

         if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_HEXAHEDRON:"<<endl;
              Out(cout);         
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
             
          FVPy.FacetsSurroundingNode(facets_surrounding_node);
          
          FVPy.EdgePairs(edges_of_element);
          
          FVPy.ProjectionWeights(facet_projection_weights);

		      FVPy.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVPy.Barycenter(barycenter);
          
          FVPy.FacetPoints(facet_points);

    	    FVPy.SectorPoints(sector_points_);
    	    
     	    FVPy.SectorEdgePairs(sector_edges_); 

         if ( debug ) {
              cout<<"\nISOPARAMETRIC_LINEAR_PYRAMID:"<<endl;
              Out(cout);         
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
             
          FVP.FacetsSurroundingNode(facets_surrounding_node);
          
          FVP.EdgePairs(edges_of_element);
          
          FVP.ProjectionWeights(facet_projection_weights);

		      FVP.EdgeMidpoints(facet_edge_midpoints);
		  
		      FVP.Barycenter(barycenter);

          FVP.FacetPoints(facet_points);
		  
         if(debug) {
              cout<<"\nISOPARAMETRIC_LINEAR_PRISM:"<<endl;
              Out(cout);         
          }

      parent_element_ = "ISOPARAMETRIC_LINEAR_PRISM";
      return; 
      }
    else {
         cout <<"\nFiniteVolumeStencil::Initialize: Element type not recognized: ";
         cout << csp_finite_element_type << endl;
      }
      
 } // end Initialize






// SKM modified version
template<size_t dim>
void FiniteVolumeStencil<dim>::Out(std::ostream& os) const
{
   os<<"\nFiniteVolumeStencil<dim>::Out: \nVolume Weights: ";
   for(size_t i=0; i<sector_integration_weights.size(); i++)
      for(size_t j=0; j<sector_integration_weights[i].size(); j++) 
        os<<" "<<sector_integration_weights[i][j];
   
   os<<endl<<endl<<" Volume IPs: "<<endl;
   for(size_t i=0; i<sector_integration_points.size(); i++) {
     for(size_t j=0; j<sector_integration_points[i].size(); j++) {
        for(size_t k=0; k<dim; k++) os<<" "<<sector_integration_points[i][j][k];
          os<<endl;
        }
   }
                
   os<<endl<<" Facet Weights: ";
   for(size_t i=0; i<facet_integration_weights.size(); i++) {
     for(size_t j=0; j<facet_integration_weights[i].size(); j++) {
       os<<" "<<facet_integration_weights[i][j] ;
      }
   }
   os<<endl;
   
   os<<endl<<" Facet IPs: "<<endl;
   for(size_t i=0; i<facet_integration_points.size(); i++) {
     for(size_t j=0; j<facet_integration_points[i].size(); j++) {
       for(size_t k=0; k<dim; k++) os<<" "<<facet_integration_points[i][j][k];
       os<<endl;
     }
   }

  os<<endl<<" Facet parametric normals inside the element: "<<endl;
   for(size_t i=0; i<facet_parametric_normals.size(); i++) {
     for(size_t j=0; j<dim; j++) {
       os<<" "<<facet_parametric_normals[i][j] ;
     }
     os<<endl;
   }
     
   os<<endl<<" Facet physical normals inside the element: "<<endl;
   for(size_t i=0; i<facet_normals.size(); i++) {
     for(size_t j=0; j<dim; j++) {
       os<<" "<<facet_normals[i][j] ;
     }
     os<<endl;
   }

    os<<endl<<" Facet surrounding node in the element: "<<endl;
    for(size_t i=0; i<Sectors(); i++) {
      for(size_t j=0; j<facets_surrounding_node[i].size(); j++) {
        os<<" "<<facets_surrounding_node[i][j] ;
      }
      os<<endl;
    }
    
    os<<endl<<" Edges pairs in the element "<<endl;
    for(size_t i=0; i<edges_of_element.size(); i++) {
      os<<" "<< edges_of_element[i].first <<", "<< edges_of_element[i].second <<" ";
    }
  os << endl << endl; 
   
} // end Out


template class FiniteVolumeStencil<1U>; 
template class FiniteVolumeStencil<2U>; 
template class FiniteVolumeStencil<3U>; 


} // end namespace csmp
