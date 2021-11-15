#ifndef ISOPARAMETRIC_QUADRATIC_TETRAHEDRON_H
#define ISOPARAMETRIC_QUADRATIC_TETRAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticTetrahedron : public FiniteElement {
  public:
    IsoparametricQuadraticTetrahedron();
    ~IsoparametricQuadraticTetrahedron();

    virtual double    Volume();
    virtual double    AspectRatio();
    virtual double    InnerRadius();
    virtual void        EdgeLengths( std::vector<double>& vec );
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes, std::vector<size_t>& fnids );
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 6U; }
    virtual size_t      CornerNodes() const  { return 4U; }
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return ISOPARAMETRIC_QUADRATIC_BAR; };

    virtual void        UnitNormalAtFaceBarycenter( size_t face, std::vector<double>& nrml );
    //virtual void      UnitNormalToFace( size_t face, std::vector<double>& unrml ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double>& IVAR,
                                                                    std::vector<double>&       NVAR ) const;

    virtual   void      N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual   void      N_AtIntegrationPoint( size_t ip, std::vector<double>& N );
    virtual   void      N_AtBaryCenter( std::vector<double>& N );
    virtual   void      JacobianAtIntegrationPoint( size_t ip );
    virtual   void      JacobianAt(const std::vector<double>& rst);

    virtual double    dN_At( DenseMatrix<DM_MIN>& DN2, const std::vector<double>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& M );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double r, double s, double t, std::vector<double>& nrst ) const;
    virtual void        Nrst( double r, double s, double t, double* nrst ) const;
    virtual void        dNr( double r, double s, double t, std::vector<double>& dNr ) const;
    virtual void        dNs( double r, double s, double t, std::vector<double>& dNs ) const;
    virtual void        dNt( double r, double s, double t, std::vector<double>& dNt ) const;

    virtual double    WeightAtIntegrationPoint( size_t i ) const;
    virtual   void      IntegrationPoint( size_t i, std::vector<double>& xyz ) const;

    virtual void        OutputNodeDataToVTK( const char* file_name,
                                             const char* var_name,
                                             DenseMatrix<DM_MIN>& DATA ) const;
                                             
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

    /// underintegrate element by just using a single quadrature point located at the barycentre
    void QuadratureRules_1Point();

    /// integrate element by using a 4 quadrature points
    void QuadratureRules_4Points();
    void QuadratureRules_5Points();
    void QuadratureRules_11Points();
    void QuadratureRules_14Points();
    void QuadratureRules_24Points();
    void QuadratureRules_45Points();

  private:
    /// returns interpolation function values at barycenter of the element faces
    std::array<double, 10> N_AtFaceBarycenter(size_t i) const;

    /// returns coordinates of barycenters of the element faces
    void FaceBarycenterCoordinates( size_t face, std::vector<double>& barycenterCoord );
  
  private:
    // (Gauss) integration point coordinates, and weights
    DenseMatrix<DM_MIN>  NXYZ, DN, IP, BEE;
    std::vector<double> W;
    double RST[3], LXY[4];
    double accDistance;
    size_t totIterations;
    size_t nonConvergenceOfProjections;
    size_t projectionCalledNTimes;
    // 3x3 matrix operations
    double Determinant( DenseMatrix<DM_MIN>& M ) const;
    double InvertMatrix( DenseMatrix<DM_MIN>& M ) const;

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double> &rst, std::vector<double> &xyz);
    void PhysicalToParametric(std::vector<double>& rst,const std::vector<double>& xyz);
    size_t n( size_t i, size_t a ) const;
};


/**

@class IsoparametricQuadraticTetrahedron  IsoparametricQuadraticTetrahedron "finite_elements/IsoparametricQuadraticTetrahedron.h"
@date 2000
@author S.K. Matthai
@author S. Geiger

Inherited from the FiniteElement base class.


@section collaboration Collaboration

In an element-centered finite volume computation, the IsoparametricQuadraticTetrahedron
collaborates with the QuadraticTriangle element.


@section implementation Implementation

Quadratic tetrahedral finite-element with local interpolation function
and Jacobian transformation capability. The node numbering of the
tetrahedron with regard to the local coordinate frame is:

          ^ t
          |
          |

        4 o
          | \
          |   \
        8 o     o 10
          |      \
   9 o    |        \
        1 o -- o --- o 3  --> s
         /     7
       /
   5 o     o 6
   /
  /
o 2
r

Note that this numbering scheme is different from that used by ICEMCFD
where node 7 (0..n-1) is the midside node of the edge 0-3, node 8 of
edge 1-3, and node 9 of edge 2-3.
*/

///  loop mapping over the nodes of the element. @test AM o.k.
inline size_t IsoparametricQuadraticTetrahedron::n( size_t i, size_t a ) const {
     if ( i+a >= 4 ) return i+a-4;
     return i+a;
}


} // end namespace csmp

#endif




