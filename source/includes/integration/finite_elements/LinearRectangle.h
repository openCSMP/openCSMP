// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LINEAR_RECTANGLE_H
#define LINEAR_RECTANGLE_H

#include "FiniteElement.h"

namespace csmp {

class LinearRectangle final : public FiniteElement {
	public:
		explicit LinearRectangle( uint32_t dim);

		virtual void IntegralNN( DenseMatrix<DM_MIN>& MASS_MATRIX );
		virtual double Volume();
		virtual void   EdgeLengths(std::vector<double>& );
		virtual uint32_t CornerNodes() const noexcept { return 4U; }
		virtual uint32_t MidSideNodes() const noexcept { return 0; }
		virtual void CornerNodes(std::vector<uint32_t>& ids) const noexcept;
		virtual void NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids) const noexcept;
    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const noexcept;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const noexcept;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const noexcept;
		virtual void Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM3>& K) const noexcept;
		virtual void N(std::vector<double>& N, const std::vector<double>& xyz);
		virtual void dN(DenseMatrix<DM_MIN>& DN);
		virtual void N_AtBaryCenter(std::vector<double>& N);
		virtual std::vector<double> UnitNormal() const;
		virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t) const { return LINEAR_BAR; };
		virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t) const { return LINEAR_BAR; };

	private:
		void nrs( std::vector<double>& N, DenseMatrix<DM_MIN>& RS, const std::vector<double>& rs );
		DenseMatrix<DM_MIN> M1_;
		mutable std::vector<double> V_;
    
}; // end class linear rectangle

} // end csmp

#endif // LINEAR_RECTANGLE_H

