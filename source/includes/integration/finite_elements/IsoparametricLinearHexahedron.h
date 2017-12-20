#ifndef ISOPARAMETRIC_LINEAR_HEXAHEDRON_H
#define ISOPARAMETRIC_LINEAR_HEXAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricLinearHexahedron : public FiniteElement {

public:

    /// options are 1, 4, or 8 integration points
    explicit IsoparametricLinearHexahedron( size_t integrationPoints = 8 );
    ~IsoparametricLinearHexahedron();

    virtual double64    Volume();
    virtual double64    AspectRatio();
    virtual double64    InnerRadius();
    virtual void        EdgeLengths( std::vector<double64>& vec );
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

    virtual void        UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    virtual void        N(std::vector<double64>& N, const std::vector<double64>& xyz );
    virtual void        N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
    virtual void        N_AtBaryCenter( std::vector<double64>& N );
    virtual void        JacobianAtIntegrationPoint( size_t ip );
    virtual void        JacobianAt( const std::vector<double64>& rst );

    virtual double64    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double64>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& DN8 );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double64 r, double64 s, double64 t, std::vector<double64>& nrst ) const;
    virtual void        Nrst( double64 r, double64 s, double64 t, double64* nrst ) const;
    virtual void        dNr( double64 r, double64 s, double64 t, std::vector<double64>& dNr ) const;
    virtual void        dNs( double64 r, double64 s, double64 t, std::vector<double64>& dNs ) const;
    virtual void        dNt( double64 r, double64 s, double64 t, std::vector<double64>& dNt ) const;

    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double64>& IVAR,
                                                                    std::vector<double64>& NVAR ) const;
                                                                    
    /// output node coordinates from parametric space to physical space
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

    /// project integration points from parametric space to physical the IPs
    void                IntegrationPointsFromParToPhys(DenseMatrix<DM_MIN>&  IPPHYS );

    /// generic projection from cube space to physical space ...
    void                ReferenceCubeParametricToPhysical(std::vector<double64> &rst, std::vector<double64> &xyz );

    void                OutputElementToRhino( const char* file_name, bool NodesOn, bool IP_PointsVolOn, bool IP_PointsFaceOn );

    virtual void        OutputNodeDataToVTK( const char* file_name, const char* var_name,
                                             DenseMatrix<DM_MIN>& DATA ) const;

  private:

    std::vector<double64>   W;
    DenseMatrix<DM_MIN>     DN, NXYZ, IP;

    double64 VolumeOfTetra (
                        size_t verticeIndex1,
                        size_t verticeIndex2,
                        size_t verticeIndex3,
                        size_t verticeIndex4
                        );

    // Correct indexes of VolHexaCSP_Ind
    double64 VolumeOfHexa();

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double64> &rst, std::vector<double64> &xyz);
    void PhysicalToParametric(std::vector<double64>& rst,const std::vector<double64>& xyz);

};



/**

@class IsoparametricLinearHexahedron  IsoparametricLinearHexahedron "finite_elements/IsoparametricLinearHexahedron.h"
@date 1998
@author S.K. Matthaei
@author Stephen G. Roberts
*/


} // end namespace csp

#endif


