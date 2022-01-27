#ifndef FINITE_ELEMENT_POLICY_H
#define FINITE_ELEMENT_POLICY_H

#include "FiniteElement.h"

namespace csmp {

template<size_t> class Node;

/// finite element policy for class Element
template<size_t dim, template<size_t> class CELL>
class FiniteElementPolicy {
  public:
    FiniteElementPolicy( FiniteElement* = nullptr );
    FiniteElementPolicy( const FiniteElementPolicy& p ) : fptr_(p.fptr_) {}
  
    /// for deferred assignment or changing the element at runtime
    void Assign( FiniteElement* fe_ptr );
    void AssignFiniteElementNullPtr() { fptr_ = nullptr; }
    
    /// returns search key to match element faces; pointers in order so that they can be searched
    std::set<Node<dim>*> CornerNodesOfFace( size_t face_id );

    /// returns set of nodes that are neighbors of the target node in this element
    std::set<Node<dim>*> CornerNodesConnectedTo( size_t node_id );

    /// the type is an enumeration that is used in the generation of finite elements
    CSMP_FEM_TYPE  FE_Type() const;
  
    /// you are allowed to switch the element at runtime and modifies the volatile data it stores
    FiniteElement* FE() const;
    
    /// true for volumes in 3D, surfaces in 2D, and line elements in 1D, else this is a lower dimensional element
    bool       IsEquidimensional() const;
    
    bool       IsLineElement() const;
    bool       IsSurfaceElement() const;
    bool       IsVolumeElement() const;

    size_t     Interpolation() const;
    bool       UsesLocalCoordinates() const;

    /// number of element integration points for current quadrature scheme
    size_t     IntegrationPoints() const;
  
    /// returns the location of the integration point in global coordinates
    Point<dim> IntegrationPoint( size_t ip ) const;
    double     WeightAtIntegrationPoint( size_t i ) const;

    /// mapping of local points from local to global coordinates
    Point<dim> RstToXYZ( Point<dim> rst ) const;

    /// interpolation function values at local points
    void    N_At( const Point<dim>& rst ) const;
    void    N_At( const Point<dim>& rst, std::vector<double>& N )  const;
  
    /// for elements with local coordinates, this method uses an iterative approach to find the correspoding location in parametric space
    void    N_AtGlobalPoint( std::vector<double>& N, const std::vector<double>& xyz ) const;
    
    void    N_AtBaryCenter( std::vector<double>& N ) const;
    void    N_AtIntegrationPoint( size_t ipoint, std::vector<double>& N ) const;

    /// first (constant) derivatives of linear interpolation functions of analytically integrated simplex
    void	  Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K) const;

	  void    dN( DenseMatrix<DM_MIN>& ) const;
  
    /// first derivatives of interpolation functions at the given node (returned is determinant of Jacobian matrix at i,j)
    double   dN_AtNode( DenseMatrix<DM_MIN>&, size_t nd, size_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at the (returned is determinant of Jacobian matrix at i,j)
    double   dN_AtBaryCenter( DenseMatrix<DM_MIN>&, size_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at Gauss point gp (returned is determinant of Jacobian matrix at i,j)
    double   dN_AtIntegrationPoint( DenseMatrix<DM_MIN>&, size_t gp, size_t dof=1 ) const;

    /// returned is determinant of Jacobian matrix at given Gauss point
    double   det_JINV_AtIntegrationPoint( size_t ipoint ) const;

    /// interpolation function products matrix for analytically integrated element
    void     IntegralNN( DenseMatrix<DM_MIN>& M ) const;
  
    /// initialises the nodes x dim matrix XY stored in the connected finite element class 
    void     CoordinateMatrix() const;

    /// integrates a scalar property over the element and returns this value
    double   PropertyIntegral( const csmp::Index& scalar_property ) const;

    /// interpolates node property values to the point of interest (in global coordinates); @attention costly for isoparametric elements
    template<class Var>
    void       PropertyValueAt( const csmp::Index&,
                                const std::vector<double>& xyz, Var& ) const;
    /// scalar version
    double   PropertyValueAt( const csmp::Index&, const std::vector<double>& xyz ) const;
  
    /// returns the value of any property interpolated to the element's center of gravity
    template<class Var>
    void       PropertyValueAtBaryCenter( const csmp::Index&, Var& ) const;

    /// scalar version
    double   PropertyValueAtBaryCenter( const csmp::Index& ) const;
  
    /// returns the value of any property interpolated to the integration point of interest
    template<class Var>
    void       PropertyValueAtIntegrationPoint( const csmp::Index&, size_t integration_point, Var& ) const;

    /// returns the value of a scalar property interpolated to the integration point of interest
    double   PropertyValueAtIntegrationPoint( const csmp::Index&, size_t integration_point ) const;

    /// returns property values at the integration points
    template<class Var>
    void       IntegrationPointPropertyVector( const csmp::Index&, std::vector<Var>& ) const;

    /// uses linear extrapolation of property values from the integration points to the nodes
    void       ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                           const std::vector<double>&integration_point_var,
                                                           std::vector<double>& node_variable ) const;

    /// reports volume, area*1, and length*1*1 products for volumetric, surfacic and line elements respectively
    double   Volume() const;

    /// ratio of longest to shortest edge=segment
    double   AspectRatio() const;
  
    /// the radius of an inscribed circle / sphere
    double   InnerRadius() const;
  
    /// number of edges = segments in the finite element of interest
    size_t     Segments() const;

    /// retrieve position of the segment (element edge) mid-points
    Point<dim> SegmentMidPoint( size_t segm ) const;

    /// returns the length of the nth segment=edge of the finite element
    double     SegmentLength( size_t segm ) const;
  
    /// returns a vector with the lengths of all finite element segments = edges
    void       SegmentLengths( std::vector<double>& ) const;

    /// retrieve barycenter of element face
    Point<dim> FaceBaryCenter( size_t face ) const;

    /// returns area of face i (to be scaled with thickness attribute if this is a lower-dimensional element)
    double     FaceArea( size_t face ) const;

    /// returns unit normal to face of element i
    Point<dim> UnitNormalToFace( size_t face ) const;

    /// returns unit normal to face of element i into STL vector; fastest version!
    void       UnitNormalToFace( size_t face, std::vector<double>& ) const;
  
    /// returns unit normal to face of element i into supplied csmp VectorVariable
    void       UnitNormalToFace( size_t face, VectorVariable<dim>& ) const;
  
    /// computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into point
    Point<dim> UnitNormal() const;

    /// fastest: computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into argument vector
    void       UnitNormal( std::vector<double>& nrml ) const;

    /// computes unit normal to element if it is planar (3D) or linear (2D), vector components are returned into VectorVariable
    void       UnitNormal( VectorVariable<dim>& nrml ) const;

  private:
    explicit FiniteElementPolicy( const CELL<dim>& );
    csmp::FiniteElement*  fptr_ = nullptr;
};

} // end csmp

#endif
