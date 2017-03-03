#ifndef ISOPARAMETRIC_QUADRATIC_PRISM_H
#define ISOPARAMETRIC_QUADRATIC_PRISM_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticPrism : public FiniteElement {

public:

    IsoparametricQuadraticPrism();
    ~IsoparametricQuadraticPrism();

    virtual double64    Volume();

    virtual double64    AspectRatio();
    virtual double64    InnerRadius();
    virtual void        EdgeLengths( std::vector<double64>& vec );
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const  { return 8; }
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 12; }
    virtual void      CounterClockwiseNodes( std::vector<size_t>& ids ) const;

    virtual void      ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                  const std::vector<double64>& IVAR,
                                                                  std::vector<double64>&       NVAR ) const;

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                  std::vector<size_t>& fnids );
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;

    // TODO: virtual void        UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    virtual   void      N(std::vector<double64>& N, const std::vector<double64>& xyz );
    virtual   void      N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
    virtual   void      N_AtBaryCenter( std::vector<double64>& N );
    virtual   void      JacobianAtIntegrationPoint( size_t ip );

    virtual double64    dN_At( DenseMatrix<DM_MIN>& dn, const std::vector<double64>& xyz  );
    virtual void 	      dN( DenseMatrix<DM_MIN>& DN8 );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double64 r, double64 s, double64 t, std::vector<double64>& nrst ) const;
    virtual void        dNr( double64 r, double64 s, double64 t, std::vector<double64>& dNr ) const;
    virtual void        dNs( double64 r, double64 s, double64 t, std::vector<double64>& dNs ) const;
    virtual void        dNt( double64 r, double64 s, double64 t, std::vector<double64>& dNt ) const;

    virtual void        IntegralNN( DenseMatrix<DM_MIN>& M );
    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

    void                OutputNodeDataToVTK( const char* file_name,
                                             const char* var_name,
                                             DenseMatrix<DM_MIN>& DATA ) const;
  private:

    std::vector<double64>   W;
    DenseMatrix<DM_MIN>     DN,NXYZ,IP,PROP;

    void ExtrapolateNodalValToIntegrationPoints(
                                            size_t NumOfVariables,
                                            std::vector<double64> &Ni,
                                            DenseMatrix<DM_MIN> &IpValues
                                            );

    void GenerateIntegrationPoints(DenseMatrix<DM_MIN> &Ip, std::vector<double64>& We);

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double64> &rst, std::vector<double64> &xyz);
    void PhysicalToParametric(std::vector<double64>& rst,const std::vector<double64>& xyz);

};

/**

@class IsoparametricQuadraticPrism  IsoparametricQuadraticPrism "finite_elements/IsoparametricQuadraticPrism.h"
@date 1998
@author S.K. Matthai
@author Stephen G. Roberts */

} // end namespace csmp

#endif


