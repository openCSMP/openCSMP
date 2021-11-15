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
    virtual ~LinearCuboid();
    virtual double   Volume();
    virtual double   AspectRatio();
    virtual double   InnerRadius();
    virtual size_t     CornerNodes() const { return 8U; }
    virtual void	     EdgeLengths( std::vector<double>& );
    virtual void       CornerNodes(std::vector<size_t>& ids) const;
    virtual void       NodesOfSegment(size_t segm_id, std::vector<size_t>& snids) const;
    virtual void       NodesOfFace(size_t face_id, std::vector<size_t>& fnids) const;
    virtual std::vector<size_t>  CornerNodesOfFace( size_t face_id ) const;
    virtual void       UnitNormalToFace(size_t face, std::vector<double>& unrml) const;
    virtual void	     OutputNodeDataToVTK(const char * file_name, const char * var_name, DenseMatrix<DM_MIN>& DATA) const;
    virtual void       N(std::vector<double>& M, const std::vector<double>& xyz);
    virtual void       N_AtGlobalPoint(std::vector<double>& M, const std::vector<double>& xyz);
    virtual void       N_AtBaryCenter(std::vector<double>& M);
    virtual void       IntegralNN(DenseMatrix<DM_MIN>& N);
    virtual void       IntegraldNdN(DenseMatrix<DM_MIN>& DN);
    virtual double     dN_At(DenseMatrix<DM_MIN>& B, const std::vector<double>& xyz);
    virtual void       CounterClockwiseNodes(std::vector<size_t>& ids) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment(size_t) const { return LINEAR_BAR; };
    virtual CSMP_FEM_TYPE  ElementTypeOfFace(size_t) const { return LINEAR_RECTANGLE; };
    // void TestElementIntegrals(DenseMatrix<DM_MIN>& XY);
  private:
    void dN_Partial_At(std::vector<double>& V, const std::vector<double>& xyz, size_t partial);
    void MidSideNodes(std::vector<size_t>& ids) const;
    void CenterOfFacePoints(DenseMatrix<DM12>& XF);
    void MidSegmentPoints(DenseMatrix<DM12>& XS);
    DenseMatrix<DM12> M1_, M2_, M3_; // Auxiliary Dense Matrix
    std::vector<double> V_; // Auxiliary vector
};

} // end csmp

#endif // LINEAR_CUBOID_H
