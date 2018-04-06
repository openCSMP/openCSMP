/**          __________ 
^s         3|         |2
|           |         |
|_____>r   0|_________|1
*/

#include "LinearRectangle.h"

using namespace std;

namespace csmp {

LinearRectangle::LinearRectangle(size_t dims) :FiniteElement(LINEAR_RECTANGLE, false, false, 1U)
	{
		dim = dims;       /**< spatial dimension of element */
		itp = 1;       /**< degree of interpolation */
		npf = 2;       /**< nodes per face */
		npe = 4;       /**< nodes per element  */
		spe = 4;       /**< segments per element  */
		fpe = 4;       /**< faces per element  */
		epe = 4;       /**< neighbors of element */
		nne = 4;       /**< typical number of elements that share each node */
		cne = 0;       /**< typical number of elements that share each integration point */
		gpe = 0;       /**< Gauss points per element for numerical integration */

		XY.Resize(npe, dim);
		UsesLocalCoordinates(false);
		Isoparametric(false);
		SurfaceElement();
		ElementType(LINEAR_RECTANGLE);
	}

	LinearRectangle::~LinearRectangle() {}

	void LinearRectangle::IntegralNN(DenseMatrix<DM_MIN>& M)
	{
		double64 vol = Volume();
		V_ = { vol / 9.,  -vol / 18.,  vol / 36., -vol / 18.,
			-vol / 18.,  vol / 9.,  -vol / 18.,  vol / 36.,
			vol / 36., -vol / 18.,  vol / 9.,  -vol / 18.,
			-vol / 18.,  vol / 36., -vol / 18.,  vol / 9. };
		size_t k(0);
		for (size_t i = 0; i < 4; ++i)
			for (size_t j = 0; j < 4; ++j)
				M(i, j) = V_[k++];
	}

	void LinearRectangle::N(std::vector<double64>& M, const std::vector<double64>& xyz)
	{
		if (dim == 2) { Nrs(M, XY, xyz); return; }
		size_t r = 0, s = 1; // Nodes on xy-plane
		if ((XY(0, 0) == XY(1, 0)) && (XY(1, 0) == XY(2, 0))) { r = 1; s = 2; } // nodes on yz-plane
		else if ((XY(0, 1) == XY(1, 1)) && (XY(1, 1) == XY(2, 1))) { r = 0; s = 2; } // nodes on xz-plane
		M1_.Resize(npe, 2);
		for (size_t i = 0; i < npe; ++i) { M1_(i, 0) = XY(i, r); M1_(i, 1) = XY(i, s); }
		Nrs(M, M1_, { xyz[r],xyz[s] });
	}

	// dN returns dN at barcy center of the element (Stephan's opinion)
	void LinearRectangle::dN(DenseMatrix<DM_MIN>& DN) 
	{
		size_t r = 0, s = 1; // Nodes on xy-plane
		if (dim == 3){
			if ((XY(0, 0) == XY(1, 0)) && (XY(1, 0) == XY(2, 0))) { r = 1; s = 2; } // nodes on yz-plane
			else if ((XY(0, 1) == XY(1, 1)) && (XY(1, 1) == XY(2, 1))) { r = 1; s = 2; } // nodes on xz-plane
		}

		//double64  rc = 0.5*(XY(2, r) + XY(0, r)), sc = 0.5*(XY(2, s) + XY(0, s)), // bary center
		double64  rc = XY(0, r) + (XY(2, r) - XY(0, r))/sqrt(3.), sc = XY(0, s) + (XY(2, s) - XY(0, s))/sqrt(3.), //
				area = (XY(1, r) - XY(0, r))*(XY(3, s) - XY(0, s));

		size_t indV = 0;
		V_ = { (sc - XY(2,s)) / area, -(sc - XY(3,s)) / area, (sc - XY(0,s)) / area, -(sc - XY(1,s)) / area,
			   (rc - XY(2,r)) / area, -(rc - XY(3,r)) / area, (rc - XY(0,r)) / area, -(rc - XY(1,r)) / area };

		DN.Resize(dim, npe);
		DN.Zero();
		for (auto i : { r , s })
			for (size_t j = 0; j < npe; ++j)
				DN(i, j) = V_[indV++];
	}

	void LinearRectangle::CornerNodes(std::vector<size_t>& ids) const { ids = { 0, 1, 2, 3 }; }

	void LinearRectangle::NodesOfSegment(size_t segm_id, std::vector<size_t>& snids) const

	{
		snids.resize(2);
		if (segm_id == 0) {
			snids[0] = 0;
			snids[1] = 1;
		}
		else if (segm_id == 1) {
			snids[0] = 1;
			snids[1] = 2;
		}
		else if (segm_id == 2) {
			snids[0] = 2;
			snids[1] = 3;
		}
		else if (segm_id == 3) {
			snids[0] = 3;
			snids[1] = 0;
		}
		else
			std::cout << "\nIsoparametricLinearQuadrilateral::NodesOfSegment: Erratic segment id requested: " << segm_id << std::endl;

	} // end NodesOfSegment

	void LinearRectangle::NodesOfFace(size_t face_id, std::vector<size_t>& fnids) const

	{
		fnids.resize(2);
		if (face_id == 0)
		{
			fnids[0] = 0;
			fnids[1] = 1;

		}
		else if (face_id == 1)
		{
			fnids[0] = 1;
			fnids[1] = 2;

		}
		else if (face_id == 2)
		{
			fnids[0] = 2;
			fnids[1] = 3;

		}
		else if (face_id == 3)
		{
			fnids[0] = 3;
			fnids[1] = 0;

		}
		else
			std::cout << "\nIsoparametricLinearQuadrilateral::NodesOfFace: Erratic input face ID: " << face_id << std::endl;
	}

	void LinearRectangle::Integral_dNT_K_dN(DenseMatrix<DM_MIN>& M, DenseMatrix<DM_MIN>& K)
	{
		size_t r = 0, s = 1; // Nodes on xy-plane
		if (dim == 3) {
			if ((XY(0, 0) == XY(1, 0)) && (XY(1, 0) == XY(2, 0))) { r = 1; s = 2; } // nodes on yz-plane
			else if ((XY(0, 1) == XY(1, 1)) && (XY(1, 1) == XY(2, 1))) { r = 1; s = 2; } // nodes on xz-plane
		}

		double64 KsLr2 = 0.5*K(s, s)*(XY(1, r) - XY(0, r))*(XY(1, r) - XY(0, r)), 
			     KrLs2 = 0.5*K(r, r)*(XY(3, s) - XY(0, s))*(XY(3, s) - XY(0, s)), 
			       vol = (XY(1, r) - XY(0, r))*(XY(3, s) - XY(0, s)), 
			         C = (KrLs2 + KsLr2)/vol/1.5;

		V_.resize(16);
		V_ = {  C,	             -C + KsLr2,	     (1 - 1.5*vol)*C,  -C + KrLs2,
			   -C + KsLr2,	      C,	            -C + KrLs2,	        (1 - 1.5*vol)*C,
			   (1 - 1.5*vol)*C,	 -C + KrLs2,	     C,	               -C + KsLr2,
			   -C + KrLs2,	      (1 - 1.5*vol)*C,	-C + KsLr2,		  	C};

		size_t vind = 0;
		M.Resize(4, 4);
		for (size_t i = 0; i < 4; ++i)
			for (size_t j = 0; j < 4; ++j)
				M(i, j) = V_[vind++];
	}

	void LinearRectangle::Nrs(std::vector<double64>& M, DenseMatrix<DM_MIN>& RS, const std::vector<double64>& rs)
	{
		double64 vol = (RS(1, 0) - RS(0, 0))*(RS(3, 1) - RS(0, 1));
		if (vol < 0) { cout << "\n negative vol Rectangle\n"; getchar(); }
		M.resize(4);
		M = { +vol*(rs[0] - RS(2, 0))*(rs[1] - RS(2, 1)),
			  -vol*(rs[0] - RS(3, 0))*(rs[1] - RS(3, 1)),
			  +vol*(rs[0] - RS(0, 0))*(rs[1] - RS(0, 1)),
			  -vol*(rs[0] - RS(1, 0))*(rs[1] - RS(1, 1)) };
	}
	
	void LinearRectangle::N_AtBaryCenter(std::vector<double64>& M) 
	{
		if (dim == 3) {
			N(M, { 0.5*(XY(2,0) + XY(0,0)), 0.5*(XY(2,1) + XY(0,1)), 0.5*(XY(2,2) + XY(0,2)) });
			return;
		}
		N(M, { 0.5*(XY(2,0) + XY(0,0)), 0.5*(XY(2,1) + XY(0,1)) });
	}

	double64 LinearRectangle::Volume() 
	{ 
		if (dim == 2) return (XY(1, 0) - XY(0, 0)) * (XY(3, 1) - XY(0, 1)); 

		// Vec1 = X1 - X0; Vec2 = X3 - X0 
		double64  X1 = XY(1, 0) - XY(0, 0), // X
			X2 = XY(3, 0) - XY(0, 0),
			Y1 = XY(1, 1) - XY(0, 1), // Y
			Y2 = XY(3, 1) - XY(0, 1),
			Z1 = XY(1, 2) - XY(0, 2), // Z
			Z2 = XY(3, 2) - XY(0, 2);

		return std::sqrt((X1*X1 + Y1*Y1 + Z1*Z1)*(X2*X2 + Y2*Y2 + Z2*Z2));
	}

	void LinearRectangle::EdgeLengths(std::vector<double64>& vec)
	{
		double64 L0 = sqrt((XY(0, 0) - XY(1, 0))*(XY(0, 0) - XY(1, 0)) 
			             + (XY(0, 1) - XY(1, 1))*(XY(0, 1) - XY(1, 1))
		                 + (XY(0, 2) - XY(1, 2))*(XY(0, 2) - XY(1, 2)));

		double64 L1 = sqrt((XY(0, 0) - XY(3, 0))*(XY(0, 0) - XY(3, 0))
					     + (XY(0, 1) - XY(3, 1))*(XY(0, 1) - XY(3, 1))
					     + (XY(0, 2) - XY(3, 2))*(XY(0, 2) - XY(3, 2)));

		vec = {L0 , L1, L0, L1};
	}

	void LinearRectangle::UnitNormal(std::vector<double64>& vc) const
	{
		if (dim == 2) {
			vc.resize(3); // points perpendicular to plane
			vc[0] = vc[1] = static_cast<double64>(0.0);
			vc[2] = static_cast<double64>(1.0);
			return;
		}
		// cross product of two vectors
		// Vec1 = X1 - X0; Vec2 = X3 - X0 
		double64  X1 = XY(1, 0) - XY(0, 0), // X
			X2 = XY(3, 0) - XY(0, 0),
			Y1 = XY(1, 1) - XY(0, 1), // Y
			Y2 = XY(3, 1) - XY(0, 1),
			Z1 = XY(1, 2) - XY(0, 2), // Z
			Z2 = XY(3, 2) - XY(0, 2);
		// Second : cross product (normal to quad)
		vc.resize(3);
		vc[0] = Z2*Y1 - Y2*Z1,
		vc[1] = Z1*X2 - X1*Z2,
		vc[2] = Y2*X1 - X2*Y1 ;
		// Third : normalization to unit length
		double64 length = sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);
		vc[0] /= length;
		vc[1] /= length;
		vc[2] /= length;
	}
} // end csmp
