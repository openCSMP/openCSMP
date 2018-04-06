/**
Standard ordering:
Assume center of cuboid (barycenter) is located on (0,0,0).
Node 0 : x,y < 0 , z > 0
Node 1 : x,z > 0 , y < 0
Node 2 : x > 0 , y,z < 0 
Node 3 : x,y,z < 0 
Nodes 4,5,6,7 are above nodes 0,1,2,3 respectively. (change y > 0).

NodesOfSegment() and NodesOfFace() determine order of segments and faces respectively.
*/
#include "LinearCuboid.h"

using namespace std;

namespace csmp {

LinearCuboid::LinearCuboid() : FiniteElement(LINEAR_CUBOID,false,false,1U), V_(8)
	{
		dim = 3;       /**< spatial dimension of element */
		itp = 1;       /**< degree of interpolation */
		npf = 4;       /**< nodes per face */
		npe = 8;       /**< nodes per element  */
		spe = 12;       /**< segments per element  */
		fpe = 6;       /**< faces per element  */
		epe = 6;       /**< neighbors of element */
		nne = 4;       /**< typical number of elements that share each node */
		cne = 1;       /**< typical number of elements that share each integration point */
		gpe = 1;       /**< Gauss points per element for numerical integration */

		UsesLocalCoordinates(false);
		Isoparametric(false);
	}

	 LinearCuboid::~LinearCuboid() {}

	 double64 LinearCuboid::Volume()
	 {
		 return (XY(5, 0) - XY(3, 0))*(XY(5, 1) - XY(3, 1))*(XY(5, 2) - XY(3, 2));
	 }

	 double64 LinearCuboid::AspectRatio()
	 {
		 double64 dx = XY(5, 0) - XY(3, 0);
		 double64 dy = XY(5, 1) - XY(3, 1);
		 double64 dz = XY(5, 2) - XY(3, 2);
		 return max({dx,dy,dz}) / min({dx,dy,dz}) / 2.0;
	 }

	 /** In cuboid diameter of the biggest inside sphere is the minimum edge length */
	 double64 LinearCuboid::InnerRadius()
	 {
		 double64 dx = XY(5, 0) - XY(3, 0);
		 double64 dy = XY(5, 1) - XY(3, 1);
		 double64 dz = XY(5, 2) - XY(3, 2);
		 return min({ dx,dy,dz }) / 2.0;
	 }

	 /** To integrate multiplication of shape functions over element.
	     We use the values of shape functions at nodes, mid-sigments, face center and barycenter.
	 */
	void LinearCuboid::IntegralNN(DenseMatrix<DM_MIN>& V) {
		DenseMatrix<DM12>& X = M1_, VS = M2_, VF = M3_;

		MidSegmentPoints(X); // mide-sigments
		// VS : Values of Shapes at mid-sigments (X)
		VS.Resize(12, 8);
		for (auto i = 0; i < 12; ++i) {
			N(V_, { X(i,0),X(i,1),X(i,2) });
			for (auto j = 0; j < 8; ++j) VS(i, j) = V_[j];
		}

		CenterOfFacePoints(X); // center of each face 
		// VF : Values of Shapes at center of faces
		VF.Resize(6, 8);
		for (auto i = 0; i < 6; ++i) {
			N(V_, { X(i,0),X(i,1),X(i,2) });
			for (auto j = 0; j < 8; ++j) VF(i, j) = V_[j];
		}

		// M : Values of Shapes at center of element
		vector<double64>& M = V_;
		N(M, { 0.5*(XY(3,0) + XY(5, 0)), 0.5*(XY(3,1) + XY(5, 1)), 0.5*(XY(3,2) + XY(5, 2)) });
	
		//--- Computing integral N_i*N_j
		double64 vol = Volume();
		V.Resize(8, 8);
		for (auto i = 0; i < 8 ; ++i)
			for (auto j = i; j < 8; ++j) {
				double64 sum = 64*M[i]*M[j]; // for center
				// if (i == j) sum += 1; // for nodes ! We add a for loop instead of this if
				for (auto k = 0; k < 12; ++k) sum += 4 * VS(k, i)*VS(k, j); // for mid-segment
				for (auto k = 0; k < 6; ++k) sum += 16 * VF(k, i)*VF(k, j); // for center of face
				V(i, j) = sum*vol/216;
				V(j, i) = V(i,j);
			}
		for (auto i = 0; i < 8; ++i) V(i, i) += vol/216; // if(i==j)
	}

	void LinearCuboid::IntegraldNdN(DenseMatrix<DM_MIN>& V)
	{
		for (auto partial=0; partial < 3; ++partial)
		{
			DenseMatrix<DM12>& X = M1_, DS = M2_, DF = M3_;
		
			MidSegmentPoints(X); // mide-sigments.
			// DS = Values of Partial Shapes at mid-sigments
			DS.Resize(12, 8);
			for (auto i = 0; i < 12; ++i) {
				dN_Partial_At(V_, { X(i,0),X(i,1),X(i,2) }, partial);
				for (auto j = 0; j < 8; ++j) DS(i, j) = V_[j];
			}

			CenterOfFacePoints(X); // center of each face 
			// DF = Values of Partial Shapes at center of faces
			DF.Resize(6, 8);
			for (auto i = 0; i < 6; ++i) {
				dN_Partial_At(V_, { X(i,0), X(i,1), X(i,2) }, partial);
				for (auto j = 0; j < 8; ++j) DF(i, j) = V_[j];
			}

			// DN = Values of Partial Shapes at Nodes
			DenseMatrix<DM12>& DN = M1_;
			DN.Resize(8, 8);
			for (auto i = 0; i < 8; ++i) {
				dN_Partial_At(V_ , { XY(i,0), XY(i,1), XY(i,2) }, partial);
				for (auto j = 0; j < 8; ++j) DN(i, j) = V_[j];
			}		

			// M = Values of Partial Shapes at center of element
			vector<double64>& M = V_;
			dN_Partial_At(M, { 0.5*(XY(3,0) + XY(5, 0)), 0.5*(XY(3,1) + XY(5, 1)), 0.5*(XY(3,2) + XY(5, 2)) }, partial);

			//--- Computing integral grad N_i dot grad N_j
			V.Resize(8, 8);
			V.Zero();
			double64 vol = Volume();
			for (auto i = 0; i < 8; ++i)
				for (auto j = i; j < 8; ++j) {
					double64 sum = 64 * M[i] * M[j]; // for center
					for (auto k = 0; k < 8; ++k) sum += DN(k, i)*DN(k, j); // for nodes
					for (auto k = 0; k < 12; ++k) sum += 4 * DS(k, i)*DS(k, j); // for mid-segment
					for (auto k = 0; k < 6; ++k) sum += 16 * DF(k, i)*DF(k, j); // for center of face
					V(i, j) += sum*vol/216;
					V(j, i) = V(i,j);
				}
		}
	}

	/** 
	To calculate each shape function corresponding to a node, we use the symmetric point (with
	respect to barycenter) and care about sign of it. Dividing by volume() gives one at that node.
	*/
	void LinearCuboid::N(std::vector<double64>& N, const std::vector<double64>& xyz) 
	{
		double64 vol = Volume();
		for(auto i = 0; i < 8; ++i)
			V_[i] = (xyz[0] - XY(i, 0))*(xyz[1] - XY(i, 1))*(xyz[2] - XY(i, 2)) / vol;
		// 0->6 1->(-7) 2->4 3->(-5)
		N[0] =  V_[6]; N[6] = -V_[0];
		N[1] = -V_[7]; N[7] =  V_[1];
		N[2] =  V_[4]; N[4] = -V_[2];
		N[3] = -V_[5]; N[5] =  V_[3];
	}

	void LinearCuboid::N_AtGlobalPoint(std::vector<double64>& M, const std::vector<double64>& xyz)
	{
		return N(M, xyz);
	}

	void LinearCuboid::N_AtBaryCenter(std::vector<double64>& N)
	{
		// xyz Center of cuboid or center of gravity
		double64 xyz[] = { 0.5*(XY(3, 0) + XY(5, 0)), 0.5*(XY(3, 1) + XY(5, 1)), 0.5*(XY(3, 2) + XY(5, 2)) };
		double64 vol = Volume();
		for (auto i = 0; i<8; ++i)
			V_[i] = (xyz[0] - XY(i, 0))*(xyz[1] - XY(i, 1))*(xyz[2] - XY(i, 2)) / vol;
		N.resize(8);
		// 0->6 1->(-7) 2->4 3->(-5)
		N[0] = V_[6]; N[6] = -V_[0];
		N[1] = -V_[7]; N[7] = V_[1];
		N[2] = V_[4]; N[4] = -V_[2];
		N[3] = -V_[5]; N[5] = V_[3];
	}

	/** Outputs shape function derivative matrix of the form:
		dN0dx ... dN7dx
	DN = dN0dy ... dN7dy
		dN0dz ... dN7dz
	*/
	double64 LinearCuboid::dN_At(DenseMatrix<DM_MIN>& DN, const vector<double64>& xyz) 
	{
		// 0->6 1->(-7) 2->4 3->(-5) for N0 use node6 ....
		DN.Resize(3, 8); //DN_size = dim x npe
		const double64 vol = Volume(), sgn[] = {1.,-1.,-1.,1.,1.,-1.,-1.,1.};
		size_t ind[] = {6,7,4,5,2,3,0,1};
		for (auto k = 0; k < 8; ++k) {
			DN(0, k) = sgn[k]*(xyz[1] - XY(ind[k], 1))*(xyz[2] - XY(ind[k], 2)) / vol;
			DN(1, k) = sgn[k]*(xyz[0] - XY(ind[k], 0))*(xyz[2] - XY(ind[k], 2)) / vol;
			DN(2, k) = sgn[k]*(xyz[0] - XY(ind[k], 0))*(xyz[1] - XY(ind[k], 1)) / vol;
		}
		return 0.0; // To satisfy base function and compiler!
	}

	/** Standard ordering is un-counter clockwise */
	void LinearCuboid::CounterClockwiseNodes(std::vector<size_t>& ids) const
	{
		ids.resize(8);
		ids[0] = 0; ids[1] = 3; ids[2] = 2; ids[3] = 1;
		ids[4] = 4; ids[5] = 7; ids[6] = 6; ids[7] = 5;
	}

	void LinearCuboid::dN_Partial_At(vector<double64>& DN, const vector<double64>& xyz, size_t partial) 
	{
		const double64 vol = Volume(), sgn[] = { 1.,-1.,-1.,1.,1.,-1.,-1.,1. };
		const size_t ind[] = { 6,7,4,5,2,3,0,1 };
		switch (partial) {
			case 0:	for (auto k = 0; k < 8; ++k)
						DN[k] = sgn[k] * (xyz[1] - XY(ind[k], 1))*(xyz[2] - XY(ind[k], 2)) / vol;
			case 1: for (auto k = 0; k < 8; ++k)
						DN[k] = sgn[k] * (xyz[0] - XY(ind[k], 0))*(xyz[2] - XY(ind[k], 2)) / vol;
			case 2: for (auto k = 0; k < 8; ++k)
						DN[k] = sgn[k] * (xyz[0] - XY(ind[k], 0))*(xyz[1] - XY(ind[k], 1)) / vol;
		}	
	}

	void LinearCuboid::MidSideNodes(std::vector<size_t>& ids) const
	{
		cerr << "\nIsoparametricLinearCuboid::MidSideNodes WARNING: MidSideNodes not present " << endl;
		ids[0] = 0;
	}

	// mide-sigments by averaging of corresponding nodes
	void LinearCuboid::MidSegmentPoints(DenseMatrix<DM12>& XS)
	{
		XS.Resize(12, 3); // 12 npe , 3 coordinates
		for (auto i = 0; i < 3; ++i) {
			XS(0, i) = 0.5*(XY(0, i) + XY(1, i));
			XS(1, i) = 0.5*(XY(1, i) + XY(2, i));
			XS(2, i) = 0.5*(XY(2, i) + XY(3, i));
			XS(3, i) = 0.5*(XY(0, i) + XY(3, i));
			XS(4, i) = 0.5*(XY(0, i) + XY(4, i));
			XS(5, i) = 0.5*(XY(1, i) + XY(5, i));
			XS(6, i) = 0.5*(XY(2, i) + XY(6, i));
			XS(7, i) = 0.5*(XY(3, i) + XY(7, i));
			XS(8, i) = 0.5*(XY(4, i) + XY(5, i));
			XS(9, i) = 0.5*(XY(5, i) + XY(6, i));
			XS(10, i) = 0.5*(XY(6, i) + XY(7, i));
			XS(11, i) = 0.5*(XY(4, i) + XY(7, i));
		}
	}

	//center of each face by average of mid-sigments
	void LinearCuboid::CenterOfFacePoints(DenseMatrix<DM12>& XF)
	{
		XF.Resize(6, 3); // 6 number of faces, 3 coordinates
		for (auto i = 0; i < 3; ++i) {
			XF(0, i) = 0.5*(XY(0, i) + XY(2, i));
			XF(1, i) = 0.5*(XY(0, i) + XY(5, i));
			XF(2, i) = 0.5*(XY(1, i) + XY(6, i));
			XF(3, i) = 0.5*(XY(2, i) + XY(7, i));
			XF(4, i) = 0.5*(XY(0, i) + XY(7, i));
			XF(5, i) = 0.5*(XY(5, i) + XY(7, i));
		}
	}

	/** for cuboid and with standard ordering */
	void LinearCuboid::EdgeLengths(std::vector<double64>& v )
	{
		const double64 dx = XY(5, 0) - XY(3, 0);
		const double64 dy = XY(5, 1) - XY(3, 1);
		const double64 dz = XY(5, 2) - XY(3, 2);
		v = {dx,dz,dx,dz,dy,dy,dy,dy,dx,dz,dx,dz};
	}

	void LinearCuboid::CornerNodes(std::vector<size_t>& ids) const
	{
		ids.resize(8);
		ids[0] = 0; ids[1] = 1; ids[2] = 2; ids[3] = 3;
		ids[4] = 4; ids[5] = 5; ids[6] = 6; ids[7] = 7;
	}

	void LinearCuboid::NodesOfSegment(size_t segm_id, std::vector<size_t>& snids) const
	{
		snids.resize(2);
		switch (segm_id) {
		case 0: snids = {0,1}; break;
		case 1: snids = {1,2}; break;
		case 2: snids = {2,3}; break;
		case 3: snids = {0,3}; break;
		case 4: snids = {0,4}; break;
		case 5: snids = {1,5}; break;
		case 6: snids = {2,6}; break;
		case 7: snids = {3,7}; break;
		case 8: snids = {4,5}; break;
		case 9: snids = {5,6}; break;
		case 10: snids = {6,7}; break;
		case 11: snids = {4,7}; break;
		}
	}

	void LinearCuboid::NodesOfFace(size_t face_id, std::vector<size_t>& fnids) const
	{
		fnids.resize(4);
		switch (face_id) {
		case 0: fnids = {0,1,2,3}; break;
		case 1: fnids = {0,1,4,5}; break;
		case 2: fnids = {1,2,5,6}; break;
		case 3: fnids = {2,3,6,7}; break;
		case 4: fnids = {0,3,4,7}; break;
		case 5: fnids = {4,5,6,7}; break;
		}
	}

	void LinearCuboid::UnitNormalToFace(size_t face, std::vector<double64>& unrml) const
	{
		unrml.resize(3);
		switch (face) {
		case 0: unrml = {0, -1, 0}; break;
		case 1: unrml = {0,  0, 1}; break;
		case 2: unrml = {1,  0, 0}; break;
		case 3: unrml = {0,  0, -1}; break;
		case 4: unrml = {-1, 0, 0}; break;
		case 5: unrml = {0, 1, 0}; break;
		}
	}

	void LinearCuboid::OutputNodeDataToVTK(const char* file_name, const char* var_name,
		DenseMatrix<DM_MIN>& DATA) const
	{
		char  outfile[NAME_STRING], elmt[30];
		strcpy(outfile, file_name);
		sprintf(elmt, "%lu", CurrentID());
		strcat(outfile, elmt);
		strcat(outfile, ".vtk");

		// 0. opening data output file in ascii format
		ofstream ofs;
		ofs.open(outfile, ios::out | ios::trunc);
		if (!ofs)
		{
			cout << "\nLinearCuboid::OutputNodeDataToVTK ";
			cout << "Output file could not be opened." << endl;
			return;
		}

		// 1. writing the file header
		// --------------------------
		ofs << "# vtk DataFile Version 2.0" << endl;
		ofs << "Finite-element dataset (CSMP): variable: " << var_name << endl;
		ofs << "ASCII" << endl << endl;

		// 2. writing node coordinates
		// ---------------------------
		DenseMatrix<DM_MIN> COORD(XY);
		ofs << "DATASET UNSTRUCTURED_GRID" << endl;
		ofs << "POINTS " << npe << " float" << endl;
		for (size_t i = 0; i<npe; i++)
		{
			for (size_t j = 0; j<dim; j++) ofs << COORD(i, j) << " ";
			ofs << endl;
		}
		ofs << endl;

		// 3. writing CELLS (cell-size and member nodes (point)) -corect
		// Cells
		// -----------------------------------------------------
		ofs << "CELLS " << 1 << " " << 9 << endl;
		// (1+4)X8
		// the 4 corner hexahedra
		ofs << 8 << " 0 1 2 3 4 5 6 7" << endl;
		ofs << endl;

		// 4. writing CELL_TYPES - for the
		// ---------------------
		ofs << "CELL_TYPES " << 1 << endl;
		ofs << 12 << endl; // VTK_HEX
		ofs << endl;

		// 5. writing POINT_DATA point-type data values
		// --------------------------------------------
		// Unfortunately the data can only be output as nodal variables
		ofs << "POINT_DATA " << npe << endl;
		ofs.setf(ios::scientific);

		if (DATA.Rows() == 1)
		{
			ofs << "SCALARS " << var_name << " float" << endl;
			ofs << "LOOKUP_TABLE default" << endl; // table must always be created
												   // matrix DATA is 1x9
			for (size_t i = 0; i<DATA.Cols(); i++) ofs << DATA(0, i) << " ";
			ofs << endl;
		}
		else
		{
			ofs << "VECTORS " << var_name << " float" << endl;
			// variables have always 3 components since view screen is 3D
			// matrix DATA is vec-dim x 10
			for (size_t i = 0; i<DATA.Cols(); i++) {
				for (size_t j = 0; j<DATA.Rows(); j++) ofs << DATA(j, i) << "  ";
				ofs << endl;
			}
		}
		ofs << endl;
		ofs.close();
		cout << "\nLinearCuboid::OutputNodeDataToVTK: file '" << outfile << "' written successfully." << endl;

	} // end OutputNodeDataToVTK

	/*
	// This function is implemented in LinearCuboid_Test
	void LinearCuboid::TestElementIntegrals(DenseMatrix<DM_MIN>& tXY)
	{
		XY = tXY;
		cout << "Nodes Position\n";
		XY.Out();
		double64 dx = XY(5, 0) - XY(3, 0);
		double64 dy = XY(5, 1) - XY(3, 1);
		double64 dz = XY(5, 2) - XY(3, 2);
		double64 vol = dx*dy*dz;

		// IntegralNN over Element
		DenseMatrix<DM_MIN> V;
		V.Resize(8, 8);
		IntegralNN(V);
		cout << " Integral NN \n ";
		V.Out();


		// IntegralNN over Element
		DenseMatrix<DM_MIN> W;
		W.Resize(8, 8);
		IntegraldNdN(W);
		cout << " Integral DNDN \n";
		W.Out();
	}
	*/
}
