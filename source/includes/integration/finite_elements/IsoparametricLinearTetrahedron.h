#ifndef ISOPARAMETRIC_LINEAR_TETRAHEDRON_H
#define ISOPARAMETRIC_LINEAR_TETRAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

/// 4-noded tetrahedron with a choice of Gaussian quadrature points
class IsoparametricLinearTetrahedron : public FiniteElement {
  public:
    explicit IsoparametricLinearTetrahedron( uint32_t integrationPoints = 4 /* 1 or 4*/ );
    ~IsoparametricLinearTetrahedron();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void      EdgeLengths( std::vector<double>& vec );
    virtual void      CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  CornerNodes() const  { return 4U; }
    virtual void      MidSideNodes( std::vector<uint32_t>& ids ) const;
    virtual void      CounterClockwiseNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  MidSideNodes() const { return 0; }

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t /* segment */ ) const { return ISOPARAMETRIC_LINEAR_BAR; };

    virtual void   NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;

    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
  
    virtual void   UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    // shape functions
    virtual void   N(std::vector<double>& N, const std::vector<double>& xyz );
    virtual void   N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N );
    virtual void   N_AtBaryCenter( std::vector<double>& N );
    virtual void   JacobianAtIntegrationPoint( uint32_t ip );
    virtual void   JacobianAt( const std::vector<double>& rst );

    // partial derivatives of individual interpolation functions at corresponding nodes
    virtual double dN( DenseMatrix<DM_MIN>& DN, const std::vector<double>& xyz  );
    virtual void 	 dN( DenseMatrix<DM_MIN>& DN4 );
    virtual double dN_AtNode( DenseMatrix<DM_MIN>& DN, uint32_t node );
    // partial derivatives of all interpolation functions at indicated gauss point 'xy'
    virtual double dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, uint32_t gauss_point );
    virtual double dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    // interpolation functions and their derivatives in local coordinates r,s,t.
    virtual void   Nrst( double r, double s, double t, std::vector<double>& nrst ) const;
    virtual void   Nrst( double r, double s, double t, double* nrst ) const;
    virtual void   dNr( double r, double s, double t, std::vector<double>& dNr ) const;
    virtual void   dNs( double r, double s, double t, std::vector<double>& dNs ) const;
    virtual void   dNt( double r, double s, double t, std::vector<double>& dNt ) const;

    virtual double WeightAtIntegrationPoint( uint32_t i ) const;
    virtual void   IntegrationPoint( uint32_t i, std::vector<double>& xyz ) const; // in global coords

    virtual void   ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                               const std::vector<double>& IVAR,
                                                               std::vector<double>& NVAR ) const;

    virtual void   OutputNodeDataToVTK( const char* file_name,
                                        const char* var_name,
                                        DenseMatrix<DM_MIN>& DATA ) const;

    void           N(std::vector<double>& N, uint32_t& iterations,	double& distance,const std::vector<double>& xyz);

    virtual void   ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

  private:

    std::vector<double>   W;
    DenseMatrix<DM_MIN>   DN, NXYZ, IP;
    double accDistance;
    uint32_t totIterations;
    uint32_t nonConvergenceOfProjections;
    uint32_t projectionCalledNTimes;

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz);
    void PhysicalToParametric(std::vector<double>& rst,const std::vector<double>& xyz);
    inline uint32_t n( uint32_t i, uint32_t a ) const;

};

/**

@class IsoparametricLinearTetrahedron  IsoparametricLinearTetrahedron "finite_elements/IsoparametricLinearTetrahedron.h"
@date 1998
@author S.K. Matthai
@author Stephen G. Roberts */


} // end namespace csmp

#endif


