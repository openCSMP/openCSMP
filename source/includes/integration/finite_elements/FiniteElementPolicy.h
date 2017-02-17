#ifndef FINITE_ELEMENT_POLICY_H
#define FINITE_ELEMENT_POLICY_H

#include "FiniteElement.h"

#include "TriangularFacet.h"
#include "QuadrilateralFacet.h"

#include "CSMP_definitions.h"
#include "CSMP_mathUtilities.h"

namespace csmp {

/// finite element policy for class Element
template<size_t dim, template<size_t> class SIMPLEX>
class FiniteElementPolicy {
  public:
    FiniteElementPolicy( FiniteElement* = nullptr );
    FiniteElementPolicy( const FiniteElementPolicy& p ) : fptr_(p.fptr_) {}
  
    /// for deferred assignment or changing the element at runtime
    void Assign( FiniteElement* fe_ptr );

    /// the type is an enumeration that is used in the generation of finite elements
    CSMP_FEM_TYPE  FE_Type() const;
  
    /// you are allowed to switch the element at runtime and modifies the volatile data it stores
    FiniteElement* FE() const;
    
    bool       IsLineElement() const;
    bool       IsSurfaceElement() const;
    bool       IsVolumeElement() const;

    size_t     Interpolation() const;
    bool       UsesLocalCoordinates() const;

    /// number of element integration points for current quadrature scheme
    size_t     IntegrationPoints() const;
  
    /// returns the location of the integration point in global coordinates
    Point<dim> IntegrationPoint( size_t ip ) const;
    double64   WeightAtIntegrationPoint( size_t i ) const;

    void       N_AtGlobalPoint( std::vector<double64>& N, const std::vector<double64>& xyz ) const;
    void       N_AtBaryCenter( std::vector<double64>& N ) const;
    void       N_AtIntegrationPoint( size_t ipoint, std::vector<double64>& N ) const;

    /// first (constant) derivatives of linear interpolation functions of analytically integrated simplex
    void       dN( DenseMatrix<DM_MIN>& ) const;
  
    /// first derivatives of interpolation functions at the given node (returned is determinant of Jacobian matrix at i,j)
    double64   dN_AtNode( DenseMatrix<DM_MIN>&, size_t nd, size_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at the (returned is determinant of Jacobian matrix at i,j)
    double64   dN_AtBaryCenter( DenseMatrix<DM_MIN>&, size_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at Gauss point gp (returned is determinant of Jacobian matrix at i,j)
    double64   dN_AtIntegrationPoint( DenseMatrix<DM_MIN>&, size_t gp, size_t dof=1 ) const;

    /// returned is determinant of Jacobian matrix at given Gauss point
    double64   det_JINV_AtIntegrationPoint( size_t ipoint ) const;

    /// interpolation function products matrix for analytically integrated element
    void       IntegralNN( DenseMatrix<DM_MIN>& M ) const;
  
    /// initialises the nodes x dim matrix XY stored in the connected finite element class 
    void       CoordinateMatrix() const;

    /// integrates a scalar property over the element and returns this value
    double64   PropertyIntegral( const csmp::Index& scalar_property ) const;

    /// interpolates node property values to the point of interest (in global coordinates); @attention costly for isoparametric elements
    template<class Var>
    void       PropertyValueAt( const csmp::Index&,
                                const std::vector<double64>& xyz, Var& ) const;

    /// returns the value of any property interpolated to the element's center of gravity
    template<class Var>
    void       PropertyValueAtBaryCenter( const csmp::Index&, Var& ) const;

    /// returns the value of any property interpolated to the integration point of interest
    template<class Var>
    void       PropertyValueAtIntegrationPoint( const csmp::Index&, size_t integration_point, Var& ) const;

    /// returns the value of a scalar property interpolated to the integration point of interest
    double64   PropertyValueAtIntegrationPoint( const csmp::Index&, size_t integration_point ) const;

    /// returns property values at the integration points
    template<class Var>
    void       IntegrationPointPropertyVector( const csmp::Index&, std::vector<Var>& ) const;

    /// uses linear extrapolation of property values from the integration points to the nodes
    void       ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                           const std::vector<double64>&integration_point_var,
                                                           std::vector<double64>& node_variable ) const;

    /// reports volume, area*1, and length*1*1 products for volumetric, surfacic and line elements respectively
    double64   Volume() const;

    /// ratio of longest to shortest edge=segment
    double64   AspectRatio() const;
  
    /// the radius of an inscribed circle / sphere
    double64   InnerRadius() const;
  
    /// number of edges = segments in the finite element of interest
    size_t     Segments() const;

    /// retrieve position of the segment (element edge) mid-points
    Point<dim> SegmentMidPoint( size_t segm ) const;

    /// returns the length of the nth segment=edge of the finite element
    double64   SegmentLength( size_t segm ) const;
  
    /// returns a vector with the lengths of all finite element segments = edges
    void       SegmentLengths( std::vector<double64>& ) const;

    /// retrieve barycenter of element face
    Point<dim> FaceBaryCenter( size_t segm ) const;

    /// returns area of face i (to be scaled with thickness attribute if this is a lower-dimensional element)
    double64   FaceArea( size_t face ) const;

    /// returns unit normal to face of element i
    Point<dim> UnitNormalToFace( size_t face ) const;

    /// returns unit normal to face of element i into STL vector
    void       UnitNormalToFace( size_t face, std::vector<double64>& ) const;
  
    /// returns unit normal to face of element i into supplied csmp VectorVariable
    void       UnitNormalToFace( size_t face, VectorVariable<dim>& ) const;
  
    /// computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into point
    Point<dim> UnitNormal() const;

    /// fastest: computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into argument vector
    void       UnitNormal( std::vector<double64>& nrml ) const;

    /// computes unit normal to element if it is planar (3D) or linear (2D), vector components are returned into VectorVariable
    void       UnitNormal( VectorVariable<dim>& nrml ) const;

  private:
    FiniteElementPolicy( const SIMPLEX<dim>& );
    csmp::FiniteElement*  fptr_;
};

} // end csmp

#endif
