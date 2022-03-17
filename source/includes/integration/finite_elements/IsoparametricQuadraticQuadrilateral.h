#ifndef ISOPARAMETRIC_QUADRATIC_QUADRILATERAL_H
#define ISOPARAMETRIC_QUADRATIC_QUADRILATERAL_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticQuadrilateral : public FiniteElement {

public:

    explicit IsoparametricQuadraticQuadrilateral( uint32_t dimensions=3U );
    ~IsoparametricQuadraticQuadrilateral();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void        EdgeLengths( std::vector<double>& vec );
    virtual void        NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;
    virtual void        NodesOfFace( uint32_t face_id, std::vector<uint32_t>& fnids ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
    virtual void        CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual void        MidSideNodes( std::vector<uint32_t>& ids ) const;
    virtual void        CounterClockwiseNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t      MidSideNodes() const { return 4U; }
    virtual uint32_t      CornerNodes() const  { return 4U; }
    virtual void        UnitNormal( std::vector<double>& vc ) const;
    // TODO: virtual void        UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    virtual CSMP_FEM_TYPE ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<uint32_t>& bnodes,
                                                    std::vector<uint32_t>& fnids );

    virtual void      N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void      N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N );
    virtual void      N_AtBaryCenter( std::vector<double>& N );
    virtual void      JacobianAtIntegrationPoint( uint32_t ip );
    virtual double    JacobianInverse();

    virtual double    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void      dN( DenseMatrix<DM_MIN>& M );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, uint32_t gauss_point );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t node );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrs( double r, double s, std::vector<double>& N ) const;
    virtual void        Nrs( double r, double s, double* N ) const;
    virtual void        dNr( double r, double s, std::vector<double>& dNr ) const;
    virtual void        dNs( double r, double s, std::vector<double>& dNs ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                    const std::vector<double>& IVAR,
                                                                    std::vector<double>& NVAR ) const;

    virtual double    WeightAtIntegrationPoint( uint32_t i ) const;
    virtual void        IntegrationPoint( uint32_t i, std::vector<double>& xyz ) const;

    virtual void        OutputNodeDataToVTK( const char* file_name,
                                             const char* var_name,
                                            DenseMatrix<DM_MIN>& DATA ) const;

    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

  private:

    // (Gauss) integration point coordinates, and weights
    std::vector<double>   rr;
    std::vector<double>   ss;
    std::vector<double>   W;
    DenseMatrix<DM_MIN>     NXY;
    DenseMatrix<DM_MIN>     DN;
    DenseMatrix<DM_MIN>     JMAT;
    DenseMatrix<DM_MIN>     BEE;

    std::vector<double>   RS, EFG;
    bool                    use2Dto3Djacobi;

    static const uint32_t parametricDimensions=2;

    void Dimensions( uint32_t dimensions );

    double  Jacobi( const std::vector<double>& rs, std::vector<double>& EFG,
                       DenseMatrix<DM_MIN>& J );

    double  Jacobi( const std::vector<double>& rs );

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz);
    void PhysicalToParametric(std::vector<double>& rst,const std::vector<double>& xyz);

};


/**

@class IsoparametricQuadraticQuadrilateral  IsoparametricQuadraticQuadrilateral "finite_elements/IsoparametricQuadraticQuadrilateral.h"
@date 2002
@author S.K. Matthaei
@author S. Geiger

@section motivation Motivation

Quadratic Quadrilateral finite-element with local interpolation function
and Jacobian transformation capability. Element numbering:


           / s-axis (of local coordinate system)
         /
 4o---7o---o3
  |        |
  8   9o   6
  |        |
  o---5o---o -----> r-axis
  1        2


The segment nodes are represented by IntegrationPoint variables. The local
coordinates of the integration points are stored in the array objects 'ri'
and 'si'.


@section implementation Implementation

Since the quadratic triangle uses local shape functions, the integration
procedure is as follows.

The weights of the interpolation functions always must add up to 1. Thus,
since there are 4 integration points, each weight is 1.0 ...

*/

}// csmp namespace

#endif

