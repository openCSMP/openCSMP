#ifndef ISOPARAMETRIC_LINEAR_HEXAHEDRON_H
#define ISOPARAMETRIC_LINEAR_HEXAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricLinearHexahedron : public FiniteElement {

public:

    /// options are 1, 4, or 8 integration points
    explicit IsoparametricLinearHexahedron( size_t integrationPoints = 8 );
    ~IsoparametricLinearHexahedron();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void        EdgeLengths( std::vector<double>& vec );
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const  { return 8; }
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 0U; }

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_LINEAR_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                       std::vector<size_t>& fnids );

    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;

    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;

    virtual std::vector<size_t>  CornerNodesOfFace( size_t face_id ) const;
  
    virtual void        UnitNormalToFace( size_t face, std::vector<double>& unrml ) const;

    virtual void        N(std::vector<double>& N, const std::vector<double>& xyz );
    virtual void        N_AtIntegrationPoint( size_t ip, std::vector<double>& N );
    virtual void        N_AtBaryCenter( std::vector<double>& N );
    virtual void        JacobianAtIntegrationPoint( size_t ip );
    virtual void        JacobianAt( const std::vector<double>& rst );

    virtual double    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& DN8 );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double r, double s, double t, std::vector<double>& nrst ) const;
    virtual void        Nrst( double r, double s, double t, double* nrst ) const;
    virtual void        dNr( double r, double s, double t, std::vector<double>& dNr ) const;
    virtual void        dNs( double r, double s, double t, std::vector<double>& dNs ) const;
    virtual void        dNt( double r, double s, double t, std::vector<double>& dNt ) const;

    virtual double    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double>& xyz ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double>& IVAR,
                                                                    std::vector<double>& NVAR ) const;
                                                                    
    /// output node coordinates from parametric space to physical space
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

    /// project integration points from parametric space to physical the IPs
    void                IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS );

    /// generic projection from cube space to physical space ...
    void                ReferenceCubeParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz );

    void                OutputElementToRhino( const char* file_name, bool NodesOn, bool IP_PointsVolOn, bool IP_PointsFaceOn );

    virtual void        OutputNodeDataToVTK( const char* file_name, const char* var_name,
                                             DenseMatrix<DM_MIN>& DATA ) const;

	virtual void	   Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K);

  private:

    std::vector<double>   W, Vx_, Vy_, Vz_, V_;
    DenseMatrix<DM_MIN>     DN, NXYZ, IP, B, BT;

    double VolumeOfTetra (
                        size_t verticeIndex1,
                        size_t verticeIndex2,
                        size_t verticeIndex3,
                        size_t verticeIndex4
                        );

    // Correct indexes of VolHexaCSP_Ind
    double VolumeOfHexa();

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz);
    void PhysicalToParametric(std::vector<double>& rst,const std::vector<double>& xyz);

};



/**

@class IsoparametricLinearHexahedron  IsoparametricLinearHexahedron "finite_elements/IsoparametricLinearHexahedron.h"
@date 1998
@author S.K. Matthaei
@author Stephen G. Roberts
*/


} // end namespace csp

#endif


