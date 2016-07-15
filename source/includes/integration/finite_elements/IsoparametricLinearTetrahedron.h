#ifndef ISOPARAMETRIC_LINEAR_TETRAHEDRON_H
#define ISOPARAMETRIC_LINEAR_TETRAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

/// 4-noded tetrahedron with a choice of Gaussian quadrature points
class IsoparametricLinearTetrahedron : public FiniteElement {
  public:
    explicit IsoparametricLinearTetrahedron( size_t integrationPoints = 4 /* 1 or 4*/ );
    ~IsoparametricLinearTetrahedron();

    virtual double64    Volume();
    virtual double64    AspectRatio();
    virtual double64    InnerRadius();
    virtual void        EdgeLengths( std::vector<double64>& vec );
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const  { return 4U; }
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 0; }

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_LINEAR_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                    std::vector<size_t>& fnids );

    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;

    virtual void        UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    // shape functions
    virtual void        N(std::vector<double64>& N, const std::vector<double64>& xyz );
    virtual void        N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
    virtual void        N_AtBaryCenter( std::vector<double64>& N );
    virtual void        JacobianAtIntegrationPoint( size_t ip );
    virtual void        JacobianAt( const std::vector<double64>& rst );

    // partial derivatives of individual interpolation functions at corresponding nodes
    virtual double64    dN( DenseMatrix<DM_MIN>& DN, const std::vector<double64>& xyz  );
    virtual void 	      dN( DenseMatrix<DM_MIN>& DN4 );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& DN, size_t node );
    // partial derivatives of all interpolation functions at indicated gauss point 'xy'
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    // interpolation functions and their derivatives in local coordinates r,s,t.
    virtual void        Nrst( double64 r, double64 s, double64 t, std::vector<double64>& nrst ) const;
    virtual void        dNr( double64 r, double64 s, double64 t, std::vector<double64>& dNr ) const;
    virtual void        dNs( double64 r, double64 s, double64 t, std::vector<double64>& dNs ) const;
    virtual void        dNt( double64 r, double64 s, double64 t, std::vector<double64>& dNt ) const;

    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double64>& xyz ) const; // in global coords

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                    const std::vector<double64>& IVAR,
                                    std::vector<double64>& NVAR ) const;

    virtual void        OutputNodeDataToVTK( const char* file_name,
                                    const char* var_name,
                                    DenseMatrix<DM_MIN>& DATA ) const;

    void                N(std::vector<double64>& N, size_t& iterations,	double64& distance,const std::vector<double64>& xyz);

    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

  private:

    std::vector<double64>   W;
    DenseMatrix<DM_MIN>     DN, NXYZ, IP;
    double64 accDistance;
    size_t totIterations;
    size_t nonConvergenceOfProjections;
    size_t projectionCalledNTimes;

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double64> &rst, std::vector<double64> &xyz);
    void PhysicalToParametric(std::vector<double64>& rst,const std::vector<double64>& xyz);
    inline size_t n( size_t i, size_t a ) const;

};

/**

@class IsoparametricLinearTetrahedron  IsoparametricLinearTetrahedron "finite_elements/IsoparametricLinearTetrahedron.h"
@date 1998
@author S.K. Matthai
@author Stephen G. Roberts */


} // end namespace csmp

#endif


