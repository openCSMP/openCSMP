#ifndef FINITE_ELEMENT_TRAITS_H
#define FINITE_ELEMENT_TRAITS_H

#include "TriangularFacet.h"
#include "QuadrilateralFacet.h"

#include "CSMP_definitions.h"
#include "CSMP_mathUtilities.h"
#include "CSMP_highLevelUtilities.h"

#include "Exception.h"

namespace csmp {

template<size_t> class Element;

/// finite element policy for class Element 
template<size_t dim, template<size_t> class SIMPLEX>
class FiniteElementTraits {

  public:

    FiniteElementTraits();

    CSMP_FEM_TYPE  FE_Type() const;
    bool       IsLineElement() const;
    bool       IsSurfaceElement() const;
    bool       IsVolumeElement() const;
    size_t     Interpolation() const;
    bool       UsesLocalCoordinates() const;

    size_t     Nodes() const;
    size_t     Neighbors() const;
    size_t     Segments() const;
    size_t     Faces() const;

    /// initializes node coordinates in connected finite element
    void       CoordinateMatrix() const;
  
    /// inputs node coordinates into supplied matrix
    void       CoordinateMatrix( DenseMatrix<DM_MIN>& XY )  const;

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
    /// returns property values at the nodes
    template<class Var>
    void       NodePropertyVector( const csmp::Index&, std::vector<Var>& ) const;



    // Finite Element Geometry( can potentially go to elemnt geometry traits)

    double64   Volume() const;
    double64   AspectRatio() const;
    double64   InnerRadius() const;

    /// projects element onto supplied line in space
    double64   LengthInDirection( const VectorVariable<dim>& vecDirection ) const;

    /// returns the coordinates of the element barycenter
    Point<dim> BaryCenter() const;

    /// retrieve barycenter of element face
    Point<dim> FaceBaryCenter( size_t segm ) const;
    /// returns area of face i (to be scaled with thickness attribute if this is a lower-dimensional element)
    double64   FaceArea( size_t face ) const;

    /// retrieve position of the segment (element edge) mid-points
    Point<dim> SegmentMidPoint( size_t segm ) const;
    /// returns the length of all element segments=edges
    double64   SegmentLength( size_t segm ) const;
    void       SegmentLengths( std::vector<double64>& ) const;

    /// returns unit normal to face of element i
    Point<dim> UnitNormalToFace( size_t face ) const;

    /// returns unit normal to face of element i into STL vector
    void       UnitNormalToFace( size_t face, std::vector<double64>& ) const;
  
    /// returns unit normal to face of element i into supplied csmp VectorVariable
    void       UnitNormalToFace( size_t face, VectorVariable<dim>& ) const;
  
    /// computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into point
    Point<dim> UnitNormal() const;

    /// computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into argument vector
    void       UnitNormal( std::vector<double64>& nrml ) const;

    /// computes unit normal to element if it is planar (3D) or linear (2D), vector components are returned into VectorVariable
    void       UnitNormal( VectorVariable<dim>& nrml ) const;

    /// TODO: @todo: SKM: questionable logic -> remove
    bool       isCorrectUnitNormalOrientation( Element<dim>* innerParent ) const;
    /// TODO: @todo: SKM: questionable logic -> remove
    bool       isCorrectUnitNormalOrientation( Element<dim>* innerParent, size_t face_id ) const;
    /// TODO: @todo: SKM: questionable logic -> remove
    bool       FromInside ( Element<dim>* parentElement, VectorVariable<dim> face_unit_normal ) const;
    /// TODO: @todo: SKM: questionable logic -> remove
    bool       FromOutside( Element<dim>* parentElement, VectorVariable<dim> face_unit_normal ) const;

  private:

    FiniteElementTraits( const SIMPLEX<dim>& );

};

} // end csmp

#endif
