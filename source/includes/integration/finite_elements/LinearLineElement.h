#ifndef LINEAR_LINE_ELEMENT_H
#define LINEAR_LINE_ELEMENT_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"

namespace csmp {

class LinearLineElement final : public FiniteElement {
  public:
    explicit LinearLineElement( uint32_t dimensions=2 );

    virtual double   Volume();
    virtual void     CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t CornerNodes() const { return 2U; }

    virtual void     NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;

    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;

    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t ) const { return LINEAR_BAR; }

    virtual void     N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void     dN( DenseMatrix<DM_MIN>& M );
    virtual double   dN_At( DenseMatrix<DM_MIN>& DN2, const std::vector<double>& xyz  );
    virtual double   dN_AtNode( DenseMatrix<DM_MIN>& M, uint32_t node );
    virtual double   dN_AtBarycenter( DenseMatrix<DM_MIN>& M );
	  virtual void     N_AtBaryCenter(std::vector<double>& N);

    virtual void     IntegralN( DenseMatrix<DM_MIN>& M );
    virtual void     IntegralNN( DenseMatrix<DM_MIN>& M );

    /// returns the unit normal computed as the cross-product between the line and the supplied vector vc
    virtual std::vector<double> UnitNormal() const;
    
    /// 2 normals located on the nodes and aligned with the elements
    virtual void     UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;

    virtual void     OutputNodeDataToVTK( const char* file_name,
                                          const char* var_name,
                                          DenseMatrix<DM_MIN>& DATA ) const;
};

/**

@class LinearLineElement  LinearLineElement "finite_elements/LinearLineElement.h"
@date 2002
@author S.K. Matthaei
@author S. Geiger */

} // end csmp


#endif





