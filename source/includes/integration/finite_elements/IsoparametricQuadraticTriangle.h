#ifndef ISOPARAMETRIC_QUADRATIC_TRIANGLE_H
#define ISOPARAMETRIC_QUADRATIC_TRIANGLE_H

#include "FiniteElement.h"

namespace csmp {

/**

@class IsoparametricQuadraticTriangle  IsoparametricQuadraticTriangle "finite_elements/IsoparametricQuadraticTriangle.h"
@date 2001
@author S.K. Matthaei
@author S. Geiger



     @section motivation Motivation


Quadratic triangular finite-element with local interpolation function
and Jacobian transformation capability. Element numbering:


           / s-axis (of local coordinate system)
       3 /
       o
      / \
 6  o     o 5
   /       \
  o -- o -- o -----> r-axis
 1     4     2


The segment nodes are represented by IntegrationPoint variables. The local
coordinates of the integration points are stored in the array objects 'ri'
and 'si'.


The coordinates of the local coordinate system are chosen such that the
nodes have the following local coordinates:

N1 = 0,0, N2 = 1,0, N3 = 0,1, N4 = 0.5,0, N5 = 0.5,0.5, N6 = 0, 0.5


@section implementation Implementation

Since the quadratic triangle uses local shape functions, the integration
procedure is as follows.

The weights of the interpolation functions always must add up to 1. Thus,
since there are 3 integration points, each weight is 0.33...

*/
class IsoparametricQuadraticTriangle : public FiniteElement {

public:

    explicit IsoparametricQuadraticTriangle( uint32_t dimensions );

    /// sets the integration points in the inside of the triangle near the nodes (default is midside node integration)
    void GaussPointsNearCorners();
    void GaussRadau4PointIntegration();
    void Dimensions( uint32_t dimensions );

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void      EdgeLengths( std::vector<double>& vec );
    virtual void      NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;
    
    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node ) const;
    
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void      CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual void      MidSideNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  MidSideNodes() const { return 3U; }
    virtual uint32_t  CornerNodes() const  { return 3U; }
    virtual std::vector<double> UnitNormal() const;
    virtual void UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    virtual void      N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void      N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N );
    virtual void      N_AtBaryCenter( std::vector<double>& N );
    virtual void      JacobianAtIntegrationPoint( uint32_t ip );
    virtual double    JacobianInverse();
    virtual double    JacobianDeterminant();

    virtual double    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void      dN( DenseMatrix<DM_MIN>& M );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, uint32_t gauss_point );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t node );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void      Nrs( double r, double s, std::vector<double>& nrst ) const;
    virtual void      Nrs( double r, double s, double* nrst ) const;
    virtual void      dNr( double r, double s, std::vector<double>& dNr ) const;
    virtual void      dNs( double r, double s, std::vector<double>& dNs ) const;

    virtual void      ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                  const std::vector<double>& IVAR,
                                                                  std::vector<double>& NVAR ) const;

    virtual double    WeightAtIntegrationPoint( uint32_t i ) const;
    virtual void      IntegrationPoint( uint32_t i, std::vector<double>& xyz ) const; // in global coords

    virtual void      OutputNodeDataToVTK( const char* file_name,
                                           const char* var_name,
                                           DenseMatrix<DM_MIN>& DATA ) const;

    virtual void      ReferenceCoordinates( DenseMatrix<DM_MIN> & matCoords ) const;

  private:

    std::vector<double>   rr;   ///< // (Gauss) integration point coordinates
    std::vector<double>   ss;   ///< // (Gauss) integration point coordinates
    std::vector<double>   W;    ///< // (Gauss) integration point weights
    DenseMatrix<DM_MIN>   NXY;
    DenseMatrix<DM_MIN>   DN;
    DenseMatrix<DM_MIN>   JMAT;
    DenseMatrix<DM_MIN>   BEE;
    std::vector<double>   LXY;
    std::vector<double>   RS, EFG;
    bool                  use2Dto3Djacobi;
    static const uint32_t   parametricDimensions = 2U;

    double  Jacobi( const std::vector<double>& rs, std::vector<double>& EFG, DenseMatrix<DM_MIN>& J );
    double  Jacobi( const std::vector<double>& rs );

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz);
    void PhysicalToParametric(std::vector<double>& rst,const std::vector<double>& xyz);
};


} // end namespace csmp

#endif




