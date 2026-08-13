#ifndef FINITE_ELEMENT_POLICY_H
#define FINITE_ELEMENT_POLICY_H

#include <set>
#include "FiniteElement.h"
#include <Eigen/Dense>

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;

/// finite element policy for class Element
template<uint32_t dim, template<uint32_t> class CELL>
class FiniteElementPolicy {
  public:
    explicit FiniteElementPolicy( const CELL<dim>& ) noexcept;
    /// takes pointer reference to avoid copying
    FiniteElementPolicy( FiniteElement* ) noexcept;
    FiniteElementPolicy( const FiniteElementPolicy& p ) noexcept : fptr_{p.fptr_} {}
    FiniteElementPolicy& operator=( const FiniteElementPolicy& p ) noexcept { if ( this != &p ) fptr_ = p.fptr_; return *this; }
    /// move semantics
    FiniteElementPolicy( FiniteElementPolicy&& p ) noexcept : fptr_{p.fptr_} { p.fptr_ = nullptr; }
    FiniteElementPolicy& operator=( FiniteElementPolicy&& p ) noexcept { std::swap(fptr_,p.fptr_); p.fptr_=nullptr; return *this; }
  
    /// for deferred assignment or changing the element at runtime
    void Assign( FiniteElement* fe_ptr );
    void AssignFiniteElementNullPtr() noexcept { fptr_ = nullptr; }
    
    /// the type is an enumeration that is used in the generation of finite elements
    CSMP_FEM_TYPE  FE_Type() const noexcept;
  
    /// you are allowed to switch the element at runtime and modifies the volatile data it stores
    FiniteElement* FE() const noexcept;
    
    /// true for volumes in 3D, surfaces in 2D, and line elements in 1D, else this is a lower dimensional element
    bool       IsEquidimensional() const noexcept;
    
    bool       IsLine() const noexcept;
    bool       IsSurface() const noexcept;
    bool       IsVolume() const noexcept;

    uint32_t   Interpolation() const noexcept;
    bool       UsesLocalCoordinates() const noexcept;

    /// number of element integration points for current quadrature scheme
    uint32_t   IntegrationPoints() const noexcept;
  
    /// returns the location of the integration point in global coordinates
    Point<dim> IntegrationPoint( uint32_t ip ) const;
    double     WeightAtIntegrationPoint( uint32_t i ) const;

    /// mapping of local points from local to global coordinates
    Point<dim> RstToXYZ( Point<dim> rst ) const;

    /// interpolation function values at local points
    void    N_At( const Point<dim>& rst ) const;
    void    N_At( const Point<dim>& rst, std::vector<double>& N )  const;
  
    /// for elements with local coordinates, this method uses an iterative approach to find the correspoding location in parametric space
    void    N_AtGlobalPoint( std::vector<double>& N, const std::vector<double>& xyz ) const;
    
    void    N_AtBaryCenter( std::vector<double>& N ) const;
    void    N_AtIntegrationPoint( uint32_t ipoint, std::vector<double>& N ) const;

    /// first (constant) derivatives of linear interpolation functions of analytically integrated simplex
    void	  Integral_dNT_K_dN( DenseMatrix<DM_MIN>& M, DenseMatrix<DM3>& MTRL ) const;

	  void    dN( DenseMatrix<DM_MIN>& ) const;
  
    /// first derivatives of interpolation functions at the given node (returned is determinant of Jacobian matrix at i,j)
    double  dN_AtNode( DenseMatrix<DM_MIN>&, uint32_t nd, uint32_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at the (returned is determinant of Jacobian matrix at i,j)
    double  dN_AtBaryCenter( DenseMatrix<DM_MIN>&, uint32_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at Gauss point gp (returned is determinant of Jacobian matrix at i,j)
    double  dN_AtIntegrationPoint( DenseMatrix<DM_MIN>&, uint32_t gp, uint32_t dof=1 ) const;

    /// returned is determinant of Jacobian matrix at given Gauss point
    double  det_J_AtIntegrationPoint( uint32_t ipoint ) const;

    /// interpolation function products matrix for analytically integrated element
    void    IntegralNN( DenseMatrix<DM_MIN>& M ) const;
  
    /// initialises the nodes x dim matrix XY stored in the connected finite element class 
    void    CoordinateMatrix() const;

    /// integrates a scalar property over the element and returns this value
    double  PropertyIntegral( const csmp::Index& scalar_property ) const;

    /// interpolates node property values to the point of interest (in global coordinates); @attention costly for isoparametric elements
    template<class Var> requires CsmpVariable<dim, Var>
    void    PropertyValueAt( const csmp::Index&, const std::vector<double>& xyz, Var& ) const;
    /// scalar version
    double  PropertyValueAt( const csmp::Index&, const std::vector<double>& xyz ) const;
  
    /// returns the value of any property interpolated to the element's center of gravity
    template<class Var> requires CsmpVariable<dim, Var>
    void    PropertyValueAtBaryCenter( const csmp::Index&, Var& ) const;

    /// scalar version
    double  PropertyValueAtBaryCenter( const csmp::Index& ) const;
  
    /// returns the value of any property interpolated to the integration point of interest
    template<class Var> requires CsmpVariable<dim, Var>
    void    PropertyValueAtIntegrationPoint( const csmp::Index& node_prop, uint32_t integration_point, Var& ) const;

    /// returns the value of a scalar property interpolated to the integration point of interest
    double  PropertyValueAtIntegrationPoint( const csmp::Index& node_prop, uint32_t integration_point ) const;
    
    /// computes the gradient of a property at the barycenter of the element, face or interface (for scalar this is a vector for vector this is a matrix); returns determinant of DN
    Eigen::Matrix<double,dim,Eigen::Dynamic> PropertyGradientAtBaryCenter( const Index& prop_key, double& detJ ) const;
    
    /// computes the gradient of a property at the integration point of the element, face or interface (for scalar this is a vector for vector this is a matrix); returns determinant of DN
    Eigen::Matrix<double,dim,Eigen::Dynamic> PropertyGradientAtIntegrationPoint( const Index& prop_key, uint32_t ip, double& detJ ) const;

    /// returns property values at the integration points
    template<class Var> requires CsmpVariable<dim, Var>
    void    IntegrationPointPropertyVector( const csmp::Index&, std::vector<Var>& ) const;

    /// uses linear extrapolation of property values from the integration points to the nodes
    void    ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                        const std::vector<double>&integration_point_var,
                                                        std::vector<double>& node_variable ) const;

    /// reports volume, area*1, and length*1*1 products for volumetric, surfacic and line elements respectively
    double  Volume() const;

    /// ratio of longest to shortest edge=segment
    double  AspectRatio() const;
  
    /// the radius of an inscribed circle / sphere
    double  InnerRadius() const;
  
    /// number of edges = segments in the finite element of interest
    uint32_t Segments() const noexcept;

    /// retrieve position of the segment (element edge) mid-points
    Point<dim> SegmentMidPoint( uint32_t segm ) const;

    /// returns the length of the nth segment=edge of the finite element
    double   SegmentLength( uint32_t segm ) const;
  
    /// returns a vector with the lengths of all finite element segments = edges
    void     SegmentLengths( std::vector<double>& ) const;

    /// returns search key to match element faces; pointers in order so that they can be searched
    std::set<Node<dim>*> CornerNodesOfFace( uint32_t face_id ) const;

    /// returns set of nodes that are neighbors of the target node in this element
    std::set<Node<dim>*> CornerNodesConnectedTo( uint32_t node_id ) const;

    /// retrieve barycenter of element face
    Point<dim> FaceBaryCenter( uint32_t face ) const;

    /// returns unit normal to face of element i
    Point<dim> UnitNormalToFace( uint32_t face ) const;

    /// returns unit normal to face of element i into STL vector; fastest version!
    void       UnitNormalToFace( uint32_t face, std::vector<double>& ) const;
  
    /// returns unit normal to face of element i into supplied csmp VectorVariable
    void       UnitNormalToFace( uint32_t face, VectorVariable<dim>& ) const;
  
    /// computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into point
    Point<dim> UnitNormal() const;

    /// fastest: computes unit normal to element if it is planar (3D) or linear (2D); vector components are returned into argument vector
    void       UnitNormal( std::vector<double>& nrml ) const;

    /// computes unit normal to element if it is planar (3D) or linear (2D), vector components are returned into VectorVariable
    void       UnitNormal( VectorVariable<dim>& nrml ) const;
    
    /// outputs finite element and discretised variable to VTK file
    void       OutputPropertyToVTK( const csmp::Index&, const char* file_name, const char* var_name ) const;
    
    friend class Element<dim>;
    friend class Face<dim>;
    friend class InterFace<dim>;

  private:
    void AccumulateGradient( Eigen::Matrix<double,dim,Eigen::Dynamic>& gradVar,
                             const DenseMatrix<DM_MIN>& DERIV, const csmp::Index& prop_key, uint32_t ip ) const;
  
    csmp::FiniteElement*  fptr_ = nullptr;
};

/**
     Computing the area of the finite element Face #
 */
template<uint32_t dim>
double faceArea( const Element<dim>& elmt, uint32_t element_face );

} // end csmp

#endif
