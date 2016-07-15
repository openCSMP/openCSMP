#ifndef ISOPARAMETRIC_LINEAR_TRIANGLE_H
#define ISOPARAMETRIC_LINEAR_TRIANGLE_H

#include "FiniteElement.h"

namespace csmp {

/// 3-noded triangle with a choice of Gaussian quadrature points and construction options for 2D and 3D models
class IsoparametricLinearTriangle : public FiniteElement {
  public:
    explicit IsoparametricLinearTriangle( size_t dimensions=2,
                                          size_t ipoints=3 );
    ~IsoparametricLinearTriangle();

    // standard interfaces
    virtual double64    Volume();
    virtual double64    AspectRatio();
    virtual double64    InnerRadius();
    virtual void        EdgeLengths( std::vector<double64>& vec );
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 0; }
    virtual size_t      CornerNodes() const  { return 3U; }
    virtual void        UnitNormal( std::vector<double64>& vc ) const;
    virtual void        UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_LINEAR_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                  std::vector<size_t>& fnids );
    // shape functions
    virtual void        N( std::vector<double64>& N, const std::vector<double64>& xyz );
    virtual void        N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
    virtual void        N_AtBaryCenter( std::vector<double64>& N );
    virtual void        JacobianAtIntegrationPoint( size_t ip );
    virtual double64    JacobianInverse();
    virtual void        JacobianAt( const std::vector<double64>& rst );
    virtual double64    JacobianDeterminant();

    // partial derivatives of individual interpolation functions at corresponding nodes
    virtual double64    dN( DenseMatrix<DM_MIN>& dn, const std::vector<double64>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& M );

    // partial derivatives of all interpolation functions at indicated gauss point 'xy'
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    // interpolation functions and their derivatives in local coordinates r,s,t.
    virtual void        Nrs( double64 r, double64 s, std::vector<double64>& nrst ) const;
    virtual void        dNr( double64 r, double64 s, std::vector<double64>& dNr ) const;
    virtual void        dNs( double64 r, double64 s, std::vector<double64>& dNs ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                               const std::vector<double64>& IVAR,
                                                               std::vector<double64>& NVAR ) const;

    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double64>& xyz ) const; // in global coords

    virtual void        OutputNodeDataToVTK( const char* file_name,
                                      const char* var_name,
                                      DenseMatrix<DM_MIN>& DATA ) const;

    void                IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS);

    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;


  private:
    // (Gauss) integration point coordinates, and weights
    std::vector<double64>   rr;
    std::vector<double64>   ss;
    std::vector<double64>   W;
    DenseMatrix<DM_MIN>     NXY;
    DenseMatrix<DM_MIN>     DN;
    DenseMatrix<DM_MIN>     JMAT;
    std::vector<double64>   LXY;
    std::vector<double64>   RS;
    bool use2Dto3Djacobi;

    static const size_t parametricDimehsions=2;

    void Dimensions( size_t dimensions );

    double64  Jacobi( const std::vector<double64>& rs, std::vector<double64>& EFG,
                       DenseMatrix<DM_MIN>& J );

    double64  Jacobi( const std::vector<double64>& rs );

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double64> &rst, std::vector<double64> &xyz);
    void PhysicalToParametric(std::vector<double64>& rst,const std::vector<double64>& xyz);

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




