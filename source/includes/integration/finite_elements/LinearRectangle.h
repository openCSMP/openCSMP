#ifndef LINEAR_RECTANGLE_H
#define LINEAR_RECTANGLE_H

#include "FiniteElement.h"

namespace csmp {

class LinearRectangle : public FiniteElement {
	public:
		LinearRectangle(size_t dim);
		virtual ~LinearRectangle();
		virtual void IntegralNN(DenseMatrix<DM_MIN>& M);
		virtual double Volume();
		virtual void   EdgeLengths(std::vector<double>& );
		virtual size_t CornerNodes() const { return 4U; }
		virtual size_t MidSideNodes() const { return 0; }
		virtual void CornerNodes(std::vector<size_t>& ids) const;
		virtual void NodesOfSegment(size_t segm_id, std::vector<size_t>& snids) const;
		virtual void NodesOfFace(size_t face_id, std::vector<size_t>& fnids) const;
    virtual std::vector<size_t>  CornerNodesOfFace( size_t face_id ) const;  
		virtual void Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K);
		virtual void N(std::vector<double>& N, const std::vector<double>& xyz);
		virtual void dN(DenseMatrix<DM_MIN>& DN);
		virtual void N_AtBaryCenter(std::vector<double>& N);
		virtual void UnitNormal(std::vector<double>& vc) const;
		virtual CSMP_FEM_TYPE  ElementTypeOfSegment(size_t) const { return LINEAR_BAR; };
		virtual CSMP_FEM_TYPE  ElementTypeOfFace(size_t) const { return LINEAR_BAR; };
	private:
		void Nrs( std::vector<double>& N, DenseMatrix<DM_MIN>& RS, const std::vector<double>& rs );
		DenseMatrix<DM_MIN>   M1_;
		std::vector<double> V_; 
    
}; // end class linear rectangle

} // end csmp

#endif // LINEAR_RECTANGLE_H

