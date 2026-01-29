#ifndef ISOPARAMETRIC_QUADRATIC_PYRAMID_H
#define ISOPARAMETRIC_QUADRATIC_PYRAMID_H

#include "IsoparametricLinearPyramid.h"

namespace csmp {

class IsoparametricQuadraticPyramid : public FiniteElement {

public:

    IsoparametricQuadraticPyramid();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void      EdgeLengths( std::vector<double>& vec );
    virtual void      CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  CornerNodes() const  { return 5U; }
    virtual void      MidSideNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  MidSideNodes() const { return 8U; }
    virtual uint32_t  NodesPerFace( uint32_t face ) const noexcept;

    virtual void      ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                  const std::vector<double>& IVAR,
                                                                  std::vector<double>&       NVAR ) const;

    virtual CSMP_FEM_TYPE ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE ElementTypeOfSegment( uint32_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void        NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;

    virtual std::vector<uint32_t> NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t> CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t> NodesConnectedTo( uint32_t node_id ) const;

    // TODO: virtual void        UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    virtual void        N(std::vector<double>& N, const std::vector<double>& xyz );
    virtual void        N_AtIntegrationPoint( uint32_t ip, std::vector<double>& N );
    virtual void        N_AtBaryCenter( std::vector<double>& N );
    virtual void        JacobianAtIntegrationPoint( uint32_t ip );

    virtual double    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void 	      dN( DenseMatrix<DM_MIN>& DN8 );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t node );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, uint32_t gauss_point );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double r, double s, double t, std::vector<double>& nrst ) const;
    virtual void        Nrst( double r, double s, double t, double* nrst ) const;
    virtual void        dNr( double r, double s, double t, std::vector<double>& dNr ) const;
    virtual void        dNs( double r, double s, double t, std::vector<double>& dNs ) const;
    virtual void        dNt( double r, double s, double t, std::vector<double>& dNt ) const;

    virtual double    WeightAtIntegrationPoint( uint32_t i ) const;
    virtual void        IntegrationPoint( uint32_t i, std::vector<double>& xyz ) const;
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

    void                OutputNodeDataToVTK( const char* file_name,
                                             const char* var_name,
                                             DenseMatrix<DM_MIN>& DATA ) const;

  private:

    std::vector<double>   W;
    DenseMatrix<DM_MIN>     DN, NXYZ, IP;

    void GenerateIntegrationPoints(DenseMatrix<DM_MIN> &Ip, std::vector<double>& We);

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz);
    void PhysicalToParametric(std::vector<double>& rst,const std::vector<double>& xyz);

};

/**

@class IsoparametricQuadraticPyramid  IsoparametricQuadraticPyramid "finite_elements/IsoparametricQuadraticPyramid.h"
@date 1998
@author S.K. Matthaei
@author Stephen G. Roberts */

} // namespace csmp

#endif



