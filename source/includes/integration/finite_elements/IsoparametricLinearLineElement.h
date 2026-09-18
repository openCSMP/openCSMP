// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ISOPARAMETRIC_LINEAR_LINE_ELEMENT_H
#define ISOPARAMETRIC_LINEAR_LINE_ELEMENT_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"

namespace csmp {

/// line element with linear interpolation functions for 1D, 2D and 3D models and variable number of integration points
class IsoparametricLinearLineElement final : public FiniteElement {
  public:
    explicit IsoparametricLinearLineElement( uint32_t dimensions = 2, uint32_t ips = 2  );

    virtual double  Volume();
    virtual double  InnerRadius();
  
    /// assuming an unscaled unit thickness, the aspect ratio of this element is equivalent to its length
    virtual double  AspectRatio();
  
    virtual void      CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual void      MidSideNodes( std::vector<uint32_t>& ids ) const;
    virtual void      NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;
    virtual uint32_t  CornerNodes() const;
    virtual uint32_t  MidSideNodes() const;

    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
    virtual void      CounterClockwiseNodes( std::vector<uint32_t>& ids ) const;
    virtual void      EdgeLengths( std::vector<double>& vec );

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const;

    virtual void      N( std::vector<double>& FN, const std::vector<double>& xy );
    virtual void      IntegrationPoint( uint32_t i, std::vector<double>& xyz ) const;
    virtual double    WeightAtIntegrationPoint( uint32_t i ) const;
    virtual void      N_AtIntegrationPoint( uint32_t IP, std::vector<double>& N );
    virtual  void     N_AtBaryCenter( std::vector<double>& N );
    virtual void      JacobianAtIntegrationPoint( uint32_t IP );
    virtual void      JacobianAt( const std::vector<double>& rst );
  
    /// inverts current JAC matrix and returns the determinant of J
    virtual double  JacobianInverse();
    virtual double  JacobianDeterminant();

    /// constant partial derivatives of interpolation functions on the element
    virtual void    Nr( double r, std::vector<double>& nr ) const;
    virtual void    Nr( double r, double* nr ) const;
    virtual void    dNr( double r, std::vector<double>& dnr ) const;
    virtual void    dN( DenseMatrix<DM_MIN>& M );
    virtual double  dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t node );
    virtual double  dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, uint32_t IP );
    virtual double  dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void    IntegralN( DenseMatrix<DM_MIN>& M );

    virtual std::vector<double>  UnitNormal() const;
  
    /// 2 normals located on the nodes and aligned with the elements
    virtual void    UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    virtual void    ExtrapolateIntegrationPointVariableToNodes( uint32_t nvars,
                                                                const std::vector<double>& IVAR,
                                                                std::vector<double>& NVAR ) const;
    virtual void    OutputNodeDataToVTK( const char* file_name,
                                         const char* var_name,
                                         DenseMatrix<DM_MIN>& DATA ) const;

    virtual void ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;

  private:
    /// Jacobians and their determinants as scale factors: dx = J * dr
    double         JacobianFor( const std::vector<double>& DNR, uint32_t spatial_dimension ) const;
    double         Jacobian1D( const std::vector<double>& DNR ) const;
    double         Jacobian2D( const std::vector<double>& DNR ) const;
    double         Jacobian3D( const std::vector<double>& DNR ) const;

  private:
    std::vector<double> IP, W;         ///< gauss-point coordinates and weights
    double              NX[2];         ///< node coordinates in parametric space
    double              current_detJ;  ///< current determinant of Jacobian matrix
};

/**

@class IsoparametricLinearLineElement  IsoparametricLinearLineElement "finite_elements/IsoparametricLinearLineElement.h"
@date 2001
@author S.K. Matthaei */

} // end namespace csmp

#endif




