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
	virtual double64   Volume();
	virtual double64   AspectRatio();
	virtual double64   InnerRadius();
	virtual size_t     CornerNodes() const { return 8U; }
	virtual void	     EdgeLengths( std::vector<double64>& );
	virtual void       CornerNodes(std::vector<size_t>& ids) const;
	virtual void       NodesOfSegment(size_t segm_id, std::vector<size_t>& snids) const;
	virtual void       NodesOfFace(size_t face_id, std::vector<size_t>& fnids) const;
	virtual void       UnitNormalToFace(size_t face, std::vector<double64>& unrml) const;
	virtual void	     OutputNodeDataToVTK(const char * file_name, const char * var_name, DenseMatrix<DM_MIN>& DATA) const;
	virtual void       N(std::vector<double64>& M, const std::vector<double64>& xyz);
	virtual void       N_AtGlobalPoint(std::vector<double64>& M, const std::vector<double64>& xyz);
	virtual void       N_AtBaryCenter(std::vector<double64>& M);
	virtual void       IntegralNN(DenseMatrix<DM_MIN>& N);
	virtual void       IntegraldNdN(DenseMatrix<DM_MIN>& DN);
	virtual double64   dN_At(DenseMatrix<DM_MIN>& B, const std::vector<double64>& xyz);
	virtual void       CounterClockwiseNodes(std::vector<size_t>& ids) const;
	virtual CSMP_FEM_TYPE  ElementTypeOfSegment(size_t) const { return LINEAR_BAR; };
	virtual CSMP_FEM_TYPE  ElementTypeOfFace(size_t) const { return LINEAR_RECTANGLE; };
	// void TestElementIntegrals(DenseMatrix<DM_MIN>& XY);
private:
	void dN_Partial_At(std::vector<double64>& V, const std::vector<double64>& xyz, size_t partial);
	void MidSideNodes(std::vector<size_t>& ids) const;
	void CenterOfFacePoints(DenseMatrix<DM12>& XF);
	void MidSegmentPoints(DenseMatrix<DM12>& XS);
	DenseMatrix<DM12> M1_, M2_, M3_; // Auxiliary Dense Matrix
	std::vector<double64> V_; // Auxiliary vector
};

} // end csmp

#endif // LINEAR_CUBOID_H
