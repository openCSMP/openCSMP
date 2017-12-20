#ifndef ISOPARAMETRIC_QUADRATIC_LINE_ELEMENT_H
#define ISOPARAMETRIC_QUADRATIC_LINE_ELEMENT_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticLineElement : public FiniteElement {

public:

    explicit IsoparametricQuadraticLineElement( size_t dimensions=2 );
    ~IsoparametricQuadraticLineElement();

    virtual double64    Volume();
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const { return 4U; }
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 4U; }
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual void        N_AtBaryCenter( std::vector<double64>& N );
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual void        EdgeLengths( std::vector<double64>& vec );
    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                  std::vector<size_t>& fnids );

    // TODO: virtual void        UnitNormalToFace( size_t face, std::vector<double64>& unrml ) const;

    virtual double64    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        N_AtIntegrationPoint( size_t IP, std::vector<double64>& N );
    virtual void        JacobianAtIntegrationPoint( size_t IP );
    virtual double64    JacobianInverse(); // returns determinant J for values of previous function
    
    virtual void        dN( DenseMatrix<DM_MIN>& M );
    virtual double64    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double64    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t IP );
    virtual double64    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        IntegralN( DenseMatrix<DM_MIN>& M );

    virtual void        UnitNormal( std::vector<double64>& vc ) const;
    
    virtual void        IntegrationPoint( size_t i, std::vector<double64>& xyz ) const;
    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                  const std::vector<double64>& IVAR, 
                                                                  std::vector<double64>& NVAR ) const; 
    virtual void        OutputNodeDataToVTK( const char* file_name,
                                           const char* var_name, 
                                           DenseMatrix<DM_MIN>& DATA ) const;

    virtual void        Nr( double64 r, std::vector<double64>& nr ) const;
    virtual void        Nr( double64 r, double64* nr ) const;
    virtual void        dNr( double64 r, std::vector<double64>& dnr ) const;
    
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;
	  
  private: 

    // Methods return determinants of Jacobians as scale factors: dx = J * dr
    double64            JacobianFor( const std::vector<double64>& DNR, size_t spatial_dimension ) const;
    double64            Jacobian1D( const std::vector<double64>& DNR ) const;
    double64            Jacobian2D( const std::vector<double64>& DNR ) const;
    double64            Jacobian3D( const std::vector<double64>& DNR ) const;

    std::vector<double64> IP, W; // gauss-point coordinates and weights
    double64            NX[3];
    mutable double64    DNXYZ[3];
    double64            current_detJ;
};

/**

@class IsoparametricQuadraticLineElement  IsoparametricQuadraticLineElement "finite_elements/IsoparametricQuadraticLineElement.h"
@date 2002
@author S.K. Matthaei */

} // end namespace csmp

#endif




