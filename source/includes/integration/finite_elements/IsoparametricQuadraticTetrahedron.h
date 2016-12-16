#ifndef ISOPARAMETRIC_QUADRATIC_TETRAHEDRON_H
#define ISOPARAMETRIC_QUADRATIC_TETRAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticTetrahedron : public FiniteElement {
  public:
    IsoparametricQuadraticTetrahedron();
    ~IsoparametricQuadraticTetrahedron();

    virtual double64    Volume();
    virtual double64    AspectRatio();
    virtual double64    InnerRadius();
    virtual void        EdgeLengths( std::vector<double64>& vec );
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

    virtual void        UnitNormalAtFaceBarycenter( size_t face, std::vector<double64>& nrml );
    //virtual void      UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                    const std::vector<double64>& IVAR,
                                                                    std::vector<double64>&       NVAR ) const;

    virtual   void      N( std::vector<double64>& N, const std::vector<double64>& xyz );
    virtual   void      N_AtIntegrationPoint( size_t ip, std::vector<double64>& N );
    virtual   void      N_AtBaryCenter( std::vector<double64>& N );
    virtual   void      JacobianAtIntegrationPoint( size_t ip );
    virtual   void      JacobianAt(const std::vector<double64>& rst);

    virtual double64    dN_At( DenseMatrix<DM_MIN>& DN2, const std::vector<double64>& xyz  );
    virtual void        dN( DenseMatrix<DM_MIN>& M );
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t gauss_point );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        Nrst( double64 r, double64 s, double64 t, std::vector<double64>& nrst ) const;
    virtual void        dNr( double64 r, double64 s, double64 t, std::vector<double64>& dNr ) const;
    virtual void        dNs( double64 r, double64 s, double64 t, std::vector<double64>& dNs ) const;
    virtual void        dNt( double64 r, double64 s, double64 t, std::vector<double64>& dNt ) const;

    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual   void      IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;

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
    std::array<double64, 10> N_AtFaceBarycenter(size_t i) const;

    /// returns coordinates of barycenters of the element faces
    void FaceBarycenterCoordinates( size_t face, std::vector<double64>& barycenterCoord );
  
  private:
    // (Gauss) integration point coordinates, and weights
    DenseMatrix<DM_MIN>  NXYZ, DN, IP, BEE;
    std::vector<double64> W;
    double64 RST[3], LXY[4];
    double64 accDistance;
    size_t totIterations;
    size_t nonConvergenceOfProjections;
    size_t projectionCalledNTimes;
    // 3x3 matrix operations
    double64 Determinant( DenseMatrix<DM_MIN>& M ) const;
    double64 InvertMatrix( DenseMatrix<DM_MIN>& M ) const;

    //Local &  Global coordinates
    void ParametricToPhysical(std::vector<double64> &rst, std::vector<double64> &xyz);
    void PhysicalToParametric(std::vector<double64>& rst,const std::vector<double64>& xyz);
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




