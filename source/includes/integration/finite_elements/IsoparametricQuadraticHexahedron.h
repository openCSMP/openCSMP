#ifndef ISOPARAMETRIC_QUADRATIC_HEXAHEDRON_H
#define ISOPARAMETRIC_QUADRATIC_HEXAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticHexahedron : public FiniteElement {

public:
    /// 4, 6 or 8 integration points are possible
    explicit IsoparametricQuadraticHexahedron( size_t integrationPoints = 8 );
    ~IsoparametricQuadraticHexahedron();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void        EdgeLengths( std::vector<double>& vec );
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const  { return 8; }
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 12; }
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
  
    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double>& IVAR,
                                                                    std::vector<double>& NVAR ) const;

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                  std::vector<size_t>& fnids );
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual std::vector<size_t>  CornerNodesOfFace( size_t face_id ) const;
  
    virtual   void      N(std::vector<double>& N, const std::vector<double>& xyz );
    virtual   void      N_AtIntegrationPoint( size_t ip, std::vector<double>& N );
    virtual   void      N_AtBaryCenter( std::vector<double>& N );
    virtual   void      JacobianAtIntegrationPoint( size_t ip );

    virtual double    dN( DenseMatrix<DM_MIN>& dn, const std::vector<double>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& DN8 );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double r, double s, double t, std::vector<double>& nrst ) const;
    virtual void        Nrst( double r, double s, double t, double* nrst ) const;
    virtual void        dNr( double r, double s, double t, std::vector<double>& dNr ) const;
    virtual void        dNs( double r, double s, double t, std::vector<double>& dNs ) const;
    virtual void        dNt( double r, double s, double t, std::vector<double>& dNt ) const;

    virtual   void      IntegrationPoint( size_t i, std::vector<double>& xyz ) const;
    virtual double    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

    void OutputNodeDataToVTK( const char* file_name,
                              const char* var_name,
                              DenseMatrix<DM_MIN>& DATA ) const;
  private:

    std::vector<double>   W;
    DenseMatrix<DM_MIN>     DN, NXYZ, IP, JETAL;

    void ExtrapolateNodalValToIntegrationPoints( size_t NumOfVariables,
                                                 std::vector<double>& Ni,
                                                 DenseMatrix<DM_MIN>& IpValues );

    /// transformations of coordinatea: local->global directly, global->local iteratively
    void ParametricToPhysical(std::vector<double>& rst, std::vector<double>& xyz );
    void PhysicalToParametric(std::vector<double>& rst, const std::vector<double>& xyz );
};

/**

@class IsoparametricQuadraticHexahedron  IsoparametricQuadraticHexahedron "finite_elements/IsoparametricQuadraticHexahedron.h"
@date 1998
@author S.K. Matthaei
@author Stephen G. Roberts */


}// csp namespace

#endif


