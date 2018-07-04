#ifndef ISOPARAMETRIC_LINEAR_LINE_ELEMENT_H
#define ISOPARAMETRIC_LINEAR_LINE_ELEMENT_H

#include "FiniteElement.h"
#include "MJL_Edge.h"
#include "MJL_Edge3D.h"

namespace csmp {

/// line element with linear interpolation functions for 1D, 2D and 3D models and variable number of integration points
class IsoparametricLinearLineElement : public FiniteElement {
  public:
    explicit IsoparametricLinearLineElement( size_t dimensions = 2, size_t ips = 1  );
    ~IsoparametricLinearLineElement();

    virtual double64  Volume();
    virtual double64  InnerRadius();
  
    /// assuming an unscaled unit thickness, the aspect ratio of this element is equivalent to its length
    virtual double64  AspectRatio();
  
    virtual void      CornerNodes( std::vector<size_t>& ids ) const;
    virtual void      MidSideNodes( std::vector<size_t>& ids ) const;
    virtual void      NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual size_t    CornerNodes() const;
    virtual size_t    MidSideNodes() const;
    virtual void      NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual void      CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual void      EdgeLengths( std::vector<double64>& vec );
    virtual void      ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                  std::vector<size_t>& fnids );

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;

    virtual void      N( std::vector<double64>& FN, const std::vector<double64>& xy );
    virtual void      IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;
    virtual double64  WeightAtIntegrationPoint( size_t i ) const;
    virtual void      N_AtIntegrationPoint( size_t IP, std::vector<double64>& N );
    virtual  void     N_AtBaryCenter( std::vector<double64>& N );
    virtual void      JacobianAtIntegrationPoint( size_t IP );
    virtual void      JacobianAt( const std::vector<double64>& rst );
  
    /// inverts current JAC matrix and returns the determinant of J
    virtual double64  JacobianInverse();
    virtual double64  JacobianDeterminant();

    /// constant partial derivatives of interpolation functions on the element
    virtual void      Nr( double64 r, std::vector<double64>& nr ) const;
    virtual void      Nr( double64 r, double64* nr ) const;
    virtual void      dNr( double64 r, std::vector<double64>& dnr ) const;
    virtual void      dN( DenseMatrix<DM_MIN>& M );
    virtual double64  dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64  dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t IP );
    virtual double64  dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void      IntegralN( DenseMatrix<DM_MIN>& M );

    virtual void      UnitNormal( std::vector<double64>& vc ) const;
  
    /// 2 normals located on the nodes and aligned with the elements
    virtual void      UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    virtual void      ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                  const std::vector<double64>& IVAR,
                                                                  std::vector<double64>& NVAR ) const;
    virtual void      OutputNodeDataToVTK( const char* file_name,
                                           const char* var_name,
                                           DenseMatrix<DM_MIN>& DATA ) const;

    virtual void ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

  private:
    /// Jacobians and their determinants as scale factors: dx = J * dr
    double64         JacobianFor( const std::vector<double64>& DNR, size_t spatial_dimension ) const;
    double64         Jacobian1D( const std::vector<double64>& DNR ) const;
    double64         Jacobian2D( const std::vector<double64>& DNR ) const;
    double64         Jacobian3D( const std::vector<double64>& DNR ) const;

  private:
    std::vector<double64> IP, W;         ///< gauss-point coordinates and weights
    double64              NX[2];         ///< node coordinates in parametric space
    double64              current_detJ;  ///< current determinant of Jacobian matrix
};

/**

@class IsoparametricLinearLineElement  IsoparametricLinearLineElement "finite_elements/IsoparametricLinearLineElement.h"
@date 2001
@author S.K. Matthaei */

} // end namespace csmp

#endif




