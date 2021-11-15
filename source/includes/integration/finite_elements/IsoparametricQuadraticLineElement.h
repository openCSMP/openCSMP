#ifndef ISOPARAMETRIC_QUADRATIC_LINE_ELEMENT_H
#define ISOPARAMETRIC_QUADRATIC_LINE_ELEMENT_H

#include "FiniteElement.h"

namespace csmp {

class IsoparametricQuadraticLineElement : public FiniteElement {

public:

    explicit IsoparametricQuadraticLineElement( size_t dimensions=2 );
    ~IsoparametricQuadraticLineElement();

    virtual double    Volume();
    virtual void        CornerNodes( std::vector<size_t>& ids ) const;
    virtual size_t      CornerNodes() const { return 4U; }
    virtual void        MidSideNodes( std::vector<size_t>& ids ) const;
    virtual size_t      MidSideNodes() const { return 4U; }
    virtual void        NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void        NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual void        N_AtBaryCenter( std::vector<double>& N );
    virtual void        CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual void        EdgeLengths( std::vector<double>& vec );
    virtual void        ConsecutiveNodesAtBoundary( const std::vector<size_t>& bnodes,
                                                  std::vector<size_t>& fnids );

    // TODO: virtual void        UnitNormalToFace( size_t face, std::vector<double>& unrml ) const;

    virtual double    WeightAtIntegrationPoint( size_t i ) const;
    virtual void        N_AtIntegrationPoint( size_t IP, std::vector<double>& N );
    virtual void        JacobianAtIntegrationPoint( size_t IP );
    virtual double    JacobianInverse(); // returns determinant J for values of previous function
    
    virtual void        dN( DenseMatrix<DM_MIN>& M );
    virtual double    dN_AtNode( DenseMatrix<DM_MIN>& M, size_t node );
    virtual double    dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& M, size_t IP );
    virtual double    dN_AtBarycenter( DenseMatrix<DM_MIN>& M );

    virtual void        IntegralN( DenseMatrix<DM_MIN>& M );

    virtual void        UnitNormal( std::vector<double>& vc ) const;
    
    virtual void        IntegrationPoint( size_t i, std::vector<double>& xyz ) const;
    virtual void        ExtrapolateIntegrationPointVariableToNodes( size_t nvars,
                                                                  const std::vector<double>& IVAR, 
                                                                  std::vector<double>& NVAR ) const; 
    virtual void        OutputNodeDataToVTK( const char* file_name,
                                           const char* var_name, 
                                           DenseMatrix<DM_MIN>& DATA ) const;

    virtual void        Nr( double r, std::vector<double>& nr ) const;
    virtual void        Nr( double r, double* nr ) const;
    virtual void        dNr( double r, std::vector<double>& dnr ) const;
    
    virtual void        ReferenceCoordinates(DenseMatrix<DM_MIN> & matCoords) const;
	  
  private: 

    // Methods return determinants of Jacobians as scale factors: dx = J * dr
    double            JacobianFor( const std::vector<double>& DNR, size_t spatial_dimension ) const;
    double            Jacobian1D( const std::vector<double>& DNR ) const;
    double            Jacobian2D( const std::vector<double>& DNR ) const;
    double            Jacobian3D( const std::vector<double>& DNR ) const;

    std::vector<double> IP, W; // gauss-point coordinates and weights
    double            NX[3];
    mutable double    DNXYZ[3];
    double            current_detJ;
};

/**

@class IsoparametricQuadraticLineElement  IsoparametricQuadraticLineElement "finite_elements/IsoparametricQuadraticLineElement.h"
@date 2002
@author S.K. Matthaei */

} // end namespace csmp

#endif




