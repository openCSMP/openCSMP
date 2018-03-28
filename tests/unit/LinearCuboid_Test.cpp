#include "LinearCuboid.h"
#include "LinearCuboid_Test.h"

using namespace std;

namespace csmp {
	LinearCuboid_Test::LinearCuboid_Test(bool verbose) :
		lcuboid_(new LinearCuboid),
		element_(new Element<3U>(lcuboid_)),
		verbose_(verbose),
		tolerance_factor_(1e-10)
	{
		n0.Idx(0);
		n1.Idx(1);
		n2.Idx(2);
		n3.Idx(3);
		n4.Idx(4);
		n5.Idx(5);
		n6.Idx(6);
		n7.Idx(7);

		element_->Idx(0);
		element_->Assign(0, &n0);
		element_->Assign(1, &n1);
		element_->Assign(2, &n2);
		element_->Assign(3, &n3);
		element_->Assign(4, &n4);
		element_->Assign(5, &n5);
		element_->Assign(6, &n6);
		element_->Assign(7, &n7);

		// 8 nodes each includes 3 coordinates :
		// Standard ordering of reference element of isoparametric linear hexahedron
		double tXY[8][3] = { -1,-1,1, 1,-1,1,  1,-1,-1, -1,-1,-1, -1,1,1, 1,1,1, 1,1,-1, -1,1,-1 };

		// Set Node coordinates
		for (size_t i = 0; i < 8; ++i)
			element_->N(i)->x(tXY[i][0]), element_->N(i)->y(tXY[i][1]), element_->N(i)->z(tXY[i][2]);

		// setting the id of the underlying finite element as it will be used to generate VTK output
		element_->FE()->CurrentID(1U);
		element_->CoordinateMatrix();

	}

	LinearCuboid_Test::~LinearCuboid_Test()
	{
	}

	void LinearCuboid_Test::run()
	{
#if 0
		DenseMatrix<DM_MIN>  XY(element_->Nodes(), 3U);
		element_->NodeCoordinateMatrix(XY);

		cout << "\n LinearCuboid_Test::run: Element coordinate matrix: \n" << endl;
		XY.Out();

		// IntegralNN over Element
		cout << "\n\nLinearCuboid_Test::run: Integral NN \n";
		DenseMatrix<DM_MIN> V;
		V.Resize(8, 8);
		element_->IntegralNN(V);
		cout << " Integral NN \n ";
		V.Out();

		cout << "\n\nLinearCuboid_Test::run: Integral dNdN \n";
		cout << "Unfortunately FiniteElementPolicy doesn't have IntegraldNdN!!\n";
		cout << "So we run it by a LinearCuboid object (not an Element object)\n";

		// Integral dNdN over Element		
		DenseMatrix<DM_MIN> W;
		W.Resize(8, 8);
		lcuboid_->IntegraldNdN(W);
		cout << "Integral DNDN \n";
		W.Out();

		cout << "\n Test Interpolation Function Values at Nodes : \n";
		TestInterpolationFunctionValues(*element_);
		cout << " End test. \n";

		cout << "\n Sum of Shape functions at barycenter : \n";
		TestSumShapesAtBaryCenter(*element_);
		cout << " End test. \n";
		
		//cerr << "\n\nCurrently LinearCuboid_Test ( Run ) terminates program. Hit Enter to exit!\n";
		//getchar();
		//exit(EXIT_SUCCESS);
#endif
	}



	void LinearCuboid_Test::TestInterpolationFunctionValues(const Element<3U>& e)
	{
		vector<double64> IPOL(element_->Nodes()), xyz(3U);
		for (size_t i = 0; i < e.Nodes(); i++) {
			Point<3U> pt = e.N(i)->Coordinate();
			xyz = pt.Coordinates();
			e.N_AtGlobalPoint(IPOL, xyz);
			// 1 at point
			_equal(IPOL[i], 1., tolerance_factor_);
			// zero everywhere else
			for (size_t j = 0; j < IPOL.size(); j++)
				if (j != i)
					_equal(IPOL[j], 0., tolerance_factor_);
			// sum = 1 (is given)
		}
	} // end

	void LinearCuboid_Test::TestSumShapesAtBaryCenter(const Element<3U>& e)
	{
		std::vector<double64> M;
		element_->N_AtBaryCenter(M);
		double sum = 0.;
		for (size_t i = 0; i < M.size(); ++i) sum += M[i];
		_equal(sum, 1., tolerance_factor_);
	}
}
