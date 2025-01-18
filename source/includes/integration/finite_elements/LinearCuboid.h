#ifndef LINEAR_CUBOID_H
#define LINEAR_CUBOID_H

#include "FiniteElement.h"

namespace csmp {

/**
    @class LinearCuboid
    @author Hani Akbari (Jan. 2017)
    This class implements linear cuboid element, so each element has 8 nodes.
*/
class LinearCuboid : public FiniteElement {
  public:
    LinearCuboid();

    virtual double   Volume();
    virtual double   AspectRatio();
    virtual double   InnerRadius();
    virtual uint32_t   CornerNodes() const { return 8U; }
    virtual void	     EdgeLengths( std::vector<double>& );
    virtual void       CornerNodes(std::vector<uint32_t>& ids) const;
    virtual void       NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids) const;

    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
    virtual void       UnitNormalToFace( uint32_t face, std::vector<double>& unrml) const;
    virtual void	     OutputNodeDataToVTK(const char * file_name, const char * var_name, DenseMatrix<DM_MIN>& DATA) const;
    virtual void       N(std::vector<double>& M, const std::vector<double>& xyz);
    virtual void       N_AtGlobalPoint(std::vector<double>& M, const std::vector<double>& xyz);
    virtual void       N_AtBaryCenter(std::vector<double>& M);
    virtual void       IntegralNN(DenseMatrix<DM_MIN>& N);
    virtual void       IntegraldNdN(DenseMatrix<DM_MIN>& DN);
    virtual double     dN_At(DenseMatrix<DM_MIN>& B, const std::vector<double>& xyz);

    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t) const { return LINEAR_BAR; };
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t) const { return LINEAR_RECTANGLE; };
    // void TestElementIntegrals(DenseMatrix<DM_MIN>& XY);
  private:
    void dN_Partial_At(std::vector<double>& V, const std::vector<double>& xyz, uint32_t partial);
    void MidSideNodes(std::vector<uint32_t>& ids) const;
    void CenterOfFacePoints(DenseMatrix<DM12>& XF);
    void MidSegmentPoints(DenseMatrix<DM12>& XS);
    DenseMatrix<DM12> M1_, M2_, M3_; // Auxiliary Dense Matrix
    std::vector<double> V_; // Auxiliary vector
};

} // end csmp

#endif // LINEAR_CUBOID_H
