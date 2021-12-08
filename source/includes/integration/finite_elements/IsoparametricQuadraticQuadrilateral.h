#ifndef ISOPARAMETRIC_QUADRATIC_QUADRILATERAL_H
#define ISOPARAMETRIC_QUADRATIC_QUADRILATERAL_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticQuadrilateral : public FiniteElement {

public:

    explicit IsoparametricQuadraticQuadrilateral( size_t dimensions=3U );
    ~IsoparametricQuadraticQuadrilateral();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void        EdgeLengths( std::vector<double>& vec );
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual std::vector<size_t>  CornerNodesOfFace( size_t face_id ) const;
    virtual std::vector<size_t>  NodesConnectedTo( size_t node_id ) const;
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 4U; }
    virtual size_t      CornerNodes() const  { return 4U; }
    virtual void        UnitNormal( std::vector<double>& vc ) const;
    // TODO: virtual void        UnitNormalToFace( size_t face, std::vector<double>& unrml ) const;

    virtual CSMP_FEM_TYPE ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                    std::vector<size_t>& fnids );

    virtual void      N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void      N_AtIntegrationPoint( size_t ip, std::vector<double>& N );
    virtual void      N_AtBaryCenter( std::vector<double>& N );
    virtual void      JacobianAtIntegrationPoint( size_t ip );
    virtual double    JacobianInverse();

    virtual double    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void      dN( DenseMatrix<DM_MIN>& M );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrs( double r, double s, std::vector<double>& N ) const;
    virtual void        Nrs( double r, double s, double* N ) const;
    virtual void        dNr( double r, double s, std::vector<double>& dNr ) const;
    virtual void        dNs( double r, double s, std::vector<double>& dNs ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double>& IVAR,
                                                                    std::vector<double>& NVAR ) const;

    virtual double    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double>& xyz ) const;

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

    static const size_t parametricDimensions=2;

    void Dimensions( size_t dimensions );

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

