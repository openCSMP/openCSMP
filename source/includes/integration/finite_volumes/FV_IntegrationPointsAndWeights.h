#ifndef FV_INTEGRATION_POINTS_AND_WEIGHTS_H
#define FV_INTEGRATION_POINTS_AND_WEIGHTS_H

#include "FiniteElement.h"
#include "Point.h"

namespace csmp {

/**

@brief Data class that contains all information required to initialize the
FiniteElementStencil objects containing the node-centered FV parametric space representations
matching the different isoparametric finite element types in CSMP.

@author S.K. Matthaei
@author S. Geiger
@author A.A.Mezentsev
@date 2003

Class fills in the C-type arrays of tabulated integration points and weights 
for further use in inplementation of the CVPEM method in CSMP. Characteristic 
feature of the class is consistent support of the mixed (poly-element) mesh.
 
@section functionality Functionality

Generates positions of integration points and corresponding weights for
numerical integration routines  

@section motivation Motivation

Part of a consistent framework for generation of the Finite Element/Finite
Volume method in CSP. Uses the concept of 
 
@section implementation Implementation

Using the template specialization mechanism of C++, the class is implemented
for different floating point types. Currently the class is supporting linear
iso-parametric family of Finite Elements with reduced polynomial interpolation
basis .
 
@section application Application Examples

No specific application examples are required, as the method is basically 
called once in the FiniteVolumeStencil class of the great code CSMP.

*/
template<size_t dim>
class FV_IntegrationPointsAndWeights {
  public:
    explicit FV_IntegrationPointsAndWeights( CSMP_FEM_TYPE type=ISOPARAMETRIC_LINEAR_TRIANGLE );
    ~FV_IntegrationPointsAndWeights();
    
    // TODO: here is a problem because physical dimension can be different from parametric dimension,
    // but dim is always greater than parametric dimension
    // maybe one should only be able to get these point by point?
    void FacetIntegrationPoints( std::vector<std::vector<Point<dim> > >&  facet_integration_points ) const;   
    void FacetIntegrationWeights( std::vector<std::vector<double64> >&  facet_integration_weights ) const;

    void SectorIntegrationPoints( std::vector<std::vector<Point<dim> > >&  volume_integration_points ) const;
    void SectorIntegrationWeights( std::vector<std::vector<double64> >&  volume_integration_weights ) const;
    
    void ProjectionWeights( std::vector<std::vector<double64> >& projection_weights ) const;
    
    void FacetNormals( std::vector<Point<dim> >& facet_normals ) const;
    
    void FacetsSurroundingNode( std::vector<std::vector<size_t> >&  facets_surrounding_node ) const;
    
    void EdgePairs( std::vector<std::pair<size_t,size_t> >& edges_of_element ) const;
    
    void EdgeMidpoints(std::vector<Point<dim> >& rst_facet_edge_midpoints );
		  
	  void Barycenter( Point<dim>& rst_barycenter );
    
    void FacetPoints( std::vector<std::vector<Point<dim> > > & rst_facet_points );
    
    void FacetTypes( std::vector<FV_FACET_TYPE> & facet_types );
    
    void SectorPoints( std::vector<std::vector<Point<dim> > >& rst_sector_points ); 

    void SectorEdgePairs( std::vector<std::vector<std::pair<size_t,size_t> > >& edges_of_sector ) const;

  private:
    FV_IntegrationPointsAndWeights();
    void Resize( size_t n_isrf, size_t srfs_per_node, size_t n_ivol, size_t n_spts, size_t n_vpts );

    CSMP_FEM_TYPE  m_typeOfSubdividedElement; ///< finite-element type to which the stencil belongs
  
    // TODO: replace by static arrays; make constant and / or use initialiser lists to set them up
    std::vector<std::pair<size_t,size_t> >   m_edges_of_element;              ///< for each edge, local node numbers on inside and outside of facet
    std::vector<std::vector<size_t> >        m_facets_surrounding_node;       ///< [node][facet] = sector-delimiting facets
    std::vector<std::vector<Point<dim> > >   m_facet_integration_points;      ///< [isrf][spts][dim] = 1 for lowest-order integration
    std::vector<std::vector<double64> >      m_facet_integration_weights;     ///< [isrf][spts]
    std::vector< Point<dim> >                m_facet_normals;                 ///< [isrf][dim]
    std::vector<std::vector<Point<dim> > >   m_volume_integration_points1;    ///< [ivol][vpts][dim]
    std::vector<std::vector<double64> >      m_volume_integration_weights;    ///< [ivol][vpts]
    std::vector<std::vector<double64> >      m_projection_weights;            ///< [isrf][spts]
    std::vector<std::vector<Point<dim> > >   m_par_volume_integration_points; ///< [ivol][vpts][dim]
    std::vector<Point<dim> >                 m_facet_edge_midpoints;          ///< points where the facet touches the finite element edges
    Point<dim>					                     m_barycenter;                    ///< center of gravity of finite element in parametric space
    std::vector< std::vector<Point<dim> > >  m_facet_points;                  ///< corner points of facet (counter-clockwise looking in),  starting with B.C.
    std::vector< FV_FACET_TYPE >             m_facet_types;                   ///< types of facet
    std::vector< std::vector<Point<dim> > >  m_sector_points;                 ///< corner points of finite volume sector corresponding to node
    std::vector<std::vector<std::pair<size_t,size_t> > > m_edge_of_sectors;   ///< edges of FV sector on outside of the finite element
    
    void CreateDataFor_ISOPARAMETRIC_LINEAR_TRI();
    void CreateDataFor_ISOPARAMETRIC_LINEAR_QUAD();
    void CreateDataFor_ISOPARAMETRIC_LINEAR_TET();
    void CreateDataFor_ISOPARAMETRIC_LINEAR_HEX();
    void CreateDataFor_ISOPARAMETRIC_LINEAR_PRISM();
    void CreateDataFor_ISOPARAMETRIC_LINEAR_PYR();
    void CreateDataFor_ISOPARAMETRIC_LINEAR_BAR();
};

} // end namespace csmp

#endif
