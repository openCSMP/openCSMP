#ifndef LINEAR_TRIANGLE_H
#define LINEAR_TRIANGLE_H

#include "FiniteElement.h"

namespace csmp {

/// Straight-sided, analytically integrated triangular element for 2D calculations
class LinearTriangle : public FiniteElement {
  public:
    LinearTriangle();

    virtual double       Volume();
    virtual double       AspectRatio();
    virtual double       InnerRadius();
    virtual void           EdgeLengths( std::vector<double>& vec );
    virtual void           NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;

    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
    virtual uint32_t       CornerNodes() const  { return 3U; }
    virtual void           CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t /* segment */ ) const { return LINEAR_BAR; };

    /// outward-pointing normals
    virtual void           UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    /// element interpolation functions at point in global coordinates
    virtual void           N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void           N_AtBaryCenter( std::vector<double>& N );

    /// first derivatives of element interpolation functions
    virtual void           dN( DenseMatrix<DM_MIN>& );
    virtual   double     dN_At( DenseMatrix<DM_MIN>&, const std::vector<double>& xyz );
    virtual   double     dN_AtBarycenter( DenseMatrix<DM_MIN>& );

    /// integrals over interpolation functions products returned into DenseMatrix
    virtual void           IntegralNN( DenseMatrix<DM_MIN>& );
    
    virtual void           OutputNodeDataToVTK( const char* file_name,
                                                const char* var_name,
                                                DenseMatrix<DM_MIN>& DATA ) const;
  private:
    void TestFunctionCoefficients( DenseMatrix<DM_MIN>& M );

    double a[3U], b[3U], c[3U]; ///< test function coefficients
};

} // end namespace csmp

#endif




