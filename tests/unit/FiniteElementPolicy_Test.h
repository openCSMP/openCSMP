#ifndef CSMP_FINITE_ELEMENT_POLICY_TEST_H
#define CSMP_FINITE_ELEMENT_POLICY_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Model;

/// Tests the finite element functionality made available via FiniteElementPolicy
class FiniteElementPolicy_Test : public Test
{
  public:
    FiniteElementPolicy_Test();
    ~FiniteElementPolicy_Test();
    
    virtual void run();
    
  private:
    Model<3U>* model_ptr_ = nullptr;
    Element<3U>* e_ptr_   = nullptr;

    void InitialiseModel();
    void CreateLinearScalarNodePropertyVariation();

    // FiniteElementPolicy interfaces trivial or tested elsewhere
    // ----------------------------------------------------------
    //std::set<Node<dim>*> CornerNodesOfFace( uint32_t face_id ) const;
    //std::set<Node<dim>*> CornerNodesConnectedTo( uint32_t node_id ) const;
    //CSMP_FE_TYPE  FE_Type() const;
    //FiniteElement* FE() const;
    //bool       IsEquidimensional() const;
    //bool       IsLine() const;
    //bool       IsSurface() const;
    //bool       IsVolume() const;
    //uint32_t   Interpolation() const;
    //bool       UsesLocalCoordinates() const;
    //uint32_t   IntegrationPoints() const;
    //Point<dim> IntegrationPoint( uint32_t ip ) const;
    //double     WeightAtIntegrationPoint( uint32_t i ) const;

    /// mapping of local points from local to global coordinates
    bool Test_RstToXYZ();

    // NOT YET
    /*
    /// interpolation function values at local points
    void    N_At( const Point<dim>& rst ) const;
    void    N_At( const Point<dim>& rst, std::vector<double>& N )  const;
  
    /// for elements with local coordinates, this method uses an iterative approach to find the correspoding location in parametric space
    void    N_AtGlobalPoint( std::vector<double>& N, const std::vector<double>& xyz ) const;
    
    void    N_AtBaryCenter( std::vector<double>& N ) const;
    void    N_AtIntegrationPoint( uint32_t ipoint, std::vector<double>& N ) const;

    /// first (constant) derivatives of linear interpolation functions of analytically integrated simplex
    void	  Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K) const;

	  void    dN( DenseMatrix<DM_MIN>& ) const;
  
    /// first derivatives of interpolation functions at the given node (returned is determinant of Jacobian matrix at i,j)
    double  dN_AtNode( DenseMatrix<DM_MIN>&, uint32_t nd, uint32_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at the (returned is determinant of Jacobian matrix at i,j)
    double  dN_AtBaryCenter( DenseMatrix<DM_MIN>&, uint32_t dof=1 ) const;
  
    /// first derivatives of interpolation functions at Gauss point gp (returned is determinant of Jacobian matrix at i,j)
    double  dN_AtIntegrationPoint( DenseMatrix<DM_MIN>&, uint32_t gp, uint32_t dof=1 ) const;

    /// returned is determinant of Jacobian matrix at given Gauss point
    double  det_JINV_AtIntegrationPoint( uint32_t ipoint ) const;

    /// interpolation function products matrix for analytically integrated element
    void    IntegralNN( DenseMatrix<DM_MIN>& M ) const;
  
    /// initialises the nodes x dim matrix XY stored in the connected finite element class
    void    CoordinateMatrix() const;
*/
    /// integrates a scalar property over the element and returns this value
    bool Test_PropertyIntegral();

    /// interpolates node property values to the point of interest (in global coordinates); @attention costly for isoparametric elements
    bool Test_PropertyValueAt();
  
    /// returns the value of any property interpolated to the element's center of gravity
    bool Test_PropertyValueAtBaryCenter();
  
    /// returns the value of any property interpolated to the integration point of interest
    bool Test_PropertyValueAtIntegrationPoint();

    /// returns property values at the integration points
    bool Test_IntegrationPointPropertyVector();

    /// uses linear extrapolation of property values from the integration points to the nodes
    bool Test_ExtrapolateIntegrationPointVariableToNodes();

    /// retrieve barycenter of element face
    bool Test_FaceBaryCenter();

    //double  Volume() const;
    //double  AspectRatio() const;
    //double  InnerRadius() const;
    //uint32_t Segments() const;
    //Point<dim> SegmentMidPoint( uint32_t segm ) const;
    //double   SegmentLength( uint32_t segm ) const;
    //void     SegmentLengths( std::vector<double>& ) const;

    /// returns unit normal to face of element i
    //Point<dim> UnitNormalToFace( uint32_t face ) const;
    //void       UnitNormalToFace( uint32_t face, std::vector<double>& ) const;
    //void       UnitNormalToFace( uint32_t face, VectorVariable<dim>& ) const;
    //Point<dim> UnitNormal() const;
    //void       UnitNormal( std::vector<double>& nrml ) const;
    //void       UnitNormal( VectorVariable<dim>& nrml ) const;

    const bool verbose_ = true;
};


}

#endif // CSMP_FINITE_ELEMENT_POLICY_TEST_H
