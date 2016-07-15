#ifndef ISOPARAMETRIC_LINEAR_PRISM_H
#define ISOPARAMETRIC_LINEAR_PRISM_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricLinearPrism : public FiniteElement {

  public:

    explicit IsoparametricLinearPrism( size_t integrationPoints = 6 /* 1 or 6*/ );
    ~IsoparametricLinearPrism();

    virtual double64    Volume();
    virtual double64    AspectRatio();
    virtual double64    InnerRadius();
    virtual void        EdgeLengths( std::vector<double64>& vec );
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const  { return 6; }
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

    virtual void        N(std::vector<double64>& N, const std::vector<double64>& xyz );
    virtual void        N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
    virtual void        N_AtBaryCenter( std::vector<double64>& N );
    virtual void        JacobianAtIntegrationPoint( size_t ip );
    virtual void        JacobianAt( const std::vector<double64>& rst );

    virtual double64    dN( DenseMatrix<DM_MIN>& dn, const std::vector<double64>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& DN6 );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double64>& IVAR,
                                                                    std::vector<double64>& NVAR ) const;

    virtual void        OutputNodeDataToVTK( const char* file_name,
                                             const char* var_name,
                                             DenseMatrix<DM_MIN>& DATA ) const;

    /// interpolation functions in parametric space
    virtual void        Nrst(double64 r, double64 s, double64 t, std::vector<double64>& N ) const;

    virtual void        dNr (
                            double64 r,
                            double64 s,
                            double64 t,
                            std::vector<double64>& dNr ) const;

    virtual void        dNs(
                            double64 r,
                            double64 s,
                            double64 t,
                            std::vector<double64>& dNs ) const;

    virtual void        dNt (
                            double64 r,
                            double64 s,
                            double64 t,
                            std::vector<double64>& dNt ) const;

    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

  private:

    std::vector<double64>   W;
    DenseMatrix<DM_MIN>     DN,NXYZ,IP;

    double64 AreaOfBase (
                       std::vector<double64> &Vertice1XYZ,
                       std::vector<double64> &Vertice2XYZ,
                       std::vector<double64> &Vertice3XYZ
                      );

    double64 VolumeOfRegularPrism();

    void GenerateIntegrationPoints( DenseMatrix<DM_MIN>& ipoints, std::vector<double64>& weights );

    void ParametricToPhysical( std::vector<double64> &rst, std::vector<double64> &xyz );
    void PhysicalToParametric( std::vector<double64>& rst,const std::vector<double64>& xyz );
};

/**

@class IsoparametricLinearPrism  IsoparametricLinearPrism "finite_elements/IsoparametricLinearPrism.h"
@date 1998
@author S.K. Matthai
@author Stephen G. Roberts */


} // end namespace csmp

#endif


