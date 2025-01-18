#ifndef ISOPARAMETRIC_LINEAR_TRIANGLE_H
#define ISOPARAMETRIC_LINEAR_TRIANGLE_H

#include "FiniteElement.h"

namespace csmp {

/// 3-noded triangle with a choice of Gaussian quadrature points and construction options for 2D and 3D models
class IsoparametricLinearTriangle : public FiniteElement {
  public:
    explicit IsoparametricLinearTriangle( uint32_t dimensions=2,
                                          uint32_t ipoints=3 );

    // standard interfaces
    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void      EdgeLengths( std::vector<double>& vec );
    virtual void      NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;

    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
    virtual void      CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual void      MidSideNodes( std::vector<uint32_t>& ids ) const;
    virtual void      CounterClockwiseNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  MidSideNodes() const { return 0; }
    virtual uint32_t  CornerNodes() const  { return 3U; }
    virtual std::vector<double> UnitNormal() const;
    virtual void      UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t /* segment */ ) const { return ISOPARAMETRIC_LINEAR_BAR; };

    // shape functions
    virtual void      N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void      N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N );
    virtual void      N_AtBaryCenter( std::vector<double>& N );
    virtual void      JacobianAtIntegrationPoint( uint32_t ip );
    virtual double    JacobianInverse();
    virtual void      JacobianAt( const std::vector<double>& rst );
    virtual double    JacobianDeterminant();

    // partial derivatives of individual interpolation functions at corresponding nodes
    virtual double    dN( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void      dN( DenseMatrix<DM_MIN>& M );

    // partial derivatives of all interpolation functions at indicated gauss point 'xy'
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, uint32_t gauss_point );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t node );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    // interpolation functions and their derivatives in local coordinates r,s,t.
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

    void              IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS);

    virtual void      ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;


  private:
    // (Gauss) integration point coordinates, and weights
    std::vector<double>   rr_;
    std::vector<double>   ss_;
    std::vector<double>   W_;
    DenseMatrix<DM_MIN>   NXY_;
    DenseMatrix<DM_MIN>   JMAT_;
    std::vector<double>   RS_;
    bool                  use2Dto3Djacobi_;

    static const uint32_t parametricDimensions_ = 2;

    void Dimensions( uint32_t dimensions );

    double  Jacobi( const std::vector<double>& rs,
                    std::vector<double>& EFG,
                    DenseMatrix<DM_MIN>& J );

    double  Jacobi( const std::vector<double>& rs );

    //Local &  Global coordinates
    void ParametricToPhysical( std::vector<double> &rst, std::vector<double> &xyz );
    void PhysicalToParametric( std::vector<double>& rst, const std::vector<double>& xyz );

};


/**

@class  IsoparametricLinearTriangle IsoparametricLinearTriangle "finite_elements/IsoparametricLinearTriangle.h"

@section motivation Motivation


Linear Triangle finite-element with local interpolation function
and Jacobian transformation capability. Element numbering:


  | s-axis (of local coordinate system)
  |
  o3
  |
  |
  |
  o ---- o -----> r-axis
 1          2


@section implementation Implementation

Since the linear triangle uses local shape functions, the integration
procedure is as follows.


The weights of the interpolation functions always must add up to 1. Thus,
since there are 4 integration points, each weight is 0.25
*/

} // end namespace csp

#endif




