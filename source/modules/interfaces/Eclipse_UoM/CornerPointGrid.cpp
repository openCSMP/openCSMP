#include "CornerPointGrid_UoM.h"
#include "PropertyData.h"
#include "CSMP_highLevelUtilities.h"
#include "STL_utilities.h"
#include "Pillar_UoM.h"

#include "ErrorHandler.h"
#include "EclipseInterface_UoM.h"

#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPrism.h"

#define CSMP_SIZE_OF_ARRAY(a)  (sizeof(a) / sizeof(a[0]))

using namespace std;

namespace csmp {

	namespace eclipse {

		enum class FACE_TYPE {

			//       3_________2
			//    /|        /|
			//     0_|_______1 |
			//     | 7-------|-6
			//     |/        |/
			//     4_________5


			FULL_QUAD,
			SPLIT_02_OR_46,
			SPLIT_13_OR_57,
			SPLIT_X
		};

		struct CellGenerator {
			CornerPointGrid_UoM& grid;

			size_t elementID = 0;
			std::vector<Point<3u>> ordinaryNodes;
			std::deque<Point<3u>> extraNodes;

			std::map<size_t, vector<size_t>>  plist;
			std::vector<int32> fem_types;

			Pillar* p0;
			Pillar* p1;
			Pillar* p2;
			Pillar* p3;

			CellGenerator(CornerPointGrid_UoM& grid)
				: grid(grid)
			{
			}

			Point<3u> getGlobalNodeCoord(size_t node) const {
				if (node < ordinaryNodes.size()) {
					return ordinaryNodes[node];
				}
				else {
					return extraNodes[node - ordinaryNodes.size()];
				}
			}

			Point<3u> getNodeCoord(ColumnCell& cell, size_t vertex) const {
				switch (vertex)
				{
				case 0:
					return p0->GetPoint(cell.z[0][0]);
				case 1:
					return p1->GetPoint(cell.z[1][0]);
				case 2:
					return p2->GetPoint(cell.z[2][0]);
				case 3:
					return p3->GetPoint(cell.z[3][0]);
				case 4:
					return p0->GetPoint(cell.z[0][1]);
				case 5:
					return p1->GetPoint(cell.z[1][1]);
				case 6:
					return p2->GetPoint(cell.z[2][1]);
				case 7:
					return p3->GetPoint(cell.z[3][1]);
				}
				throw csmp::Exception(ERROR, "CornerPointGrid::CellGenerator::getNodeCoord",
					"Node id out of range");
			}

			size_t getNodeID(ColumnCell& cell, size_t vertex) const {
				switch (vertex)
				{
				case 0:
					return cell.z[0][0] + p0->FirstNodeNum();
				case 1:
					return cell.z[1][0] + p1->FirstNodeNum();
				case 2:
					return cell.z[2][0] + p2->FirstNodeNum();
				case 3:
					return cell.z[3][0] + p3->FirstNodeNum();
				case 4:
					return cell.z[0][1] + p0->FirstNodeNum();
				case 5:
					return cell.z[1][1] + p1->FirstNodeNum();
				case 6:
					return cell.z[2][1] + p2->FirstNodeNum();
				case 7:
					return cell.z[3][1] + p3->FirstNodeNum();
				}
				throw csmp::Exception(ERROR, "CornerPointGrid::CellGenerator::getNodeID",
					"Node id out of range");
			}


			void setIJ(size_t i, size_t j) {
				p0 = &grid(i + 0, j + 0);
				p1 = &grid(i + 1, j + 0);
				p2 = &grid(i + 1, j + 1);
				p3 = &grid(i + 0, j + 1);
			}

			size_t generateCellCentroid(ColumnCell& cell) {
				Point<3> p(0, 0, 0);
				p += p0->GetPoint(cell.z[0][0]);
				p += p0->GetPoint(cell.z[0][1]);
				p += p1->GetPoint(cell.z[1][0]);
				p += p1->GetPoint(cell.z[1][1]);
				p += p2->GetPoint(cell.z[2][0]);
				p += p2->GetPoint(cell.z[2][1]);
				p += p3->GetPoint(cell.z[3][0]);
				p += p3->GetPoint(cell.z[3][1]);
				extraNodes.push_back(p * 0.125);
				return ordinaryNodes.size() + extraNodes.size() - 1;
			}

			size_t vertexIDs[8];
			IsoparametricLinearHexahedron hexa;
			IsoparametricLinearPyramid pyra;
			IsoparametricLinearTetrahedron tetra;
			IsoparametricLinearPrism prism;

			bool ConstructHexahedron(ColumnCell& cell) {

				//std::cerr << "hexa: \n";				
				for (size_t i = 0; i < 8; ++i) {
					vertexIDs[i] = i;
					auto p = getNodeCoord(cell, vertexIDs[i]);					
					grid.ConvertFromReservoirToCSMPcoordinateSystem(p);					
					hexa.XYZ(i, 0, p[0]);
					hexa.XYZ(i, 1, p[1]);
					hexa.XYZ(i, 2, p[2]);
					//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
				}
				
				size_t iNrIps = hexa.IntegrationPoints();
				for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
					hexa.JacobianAtIntegrationPoint(iIp);
					double64 jacdet = hexa.JacobianDeterminant();
					//std::cerr << "jacdet = " << jacdet << '\n';
					if (jacdet <= 0) {
						return false;
					}
				}
				return true;
			}

			size_t EmitHexahedron(ColumnCell& cell) {
				auto ids = getGlobalIDList(cell, 8, vertexIDs);
				size_t elid = elementID++;
				plist.emplace(elid, ids);
				fem_types.push_back(ISOPARAMETRIC_LINEAR_HEXAHEDRON);
				return elid;
			}

			bool ConstructPyramidOnFace(ColumnCell& cell, size_t face0, size_t face1, size_t face2, size_t face3, size_t apex) {
				vertexIDs[0] = getNodeID(cell, face0);
				vertexIDs[1] = getNodeID(cell, face1);
				vertexIDs[2] = getNodeID(cell, face2);
				vertexIDs[3] = getNodeID(cell, face3);
				vertexIDs[4] = apex;

				//std::cerr << "pyramid: \n";
				for (size_t i = 0; i < 5; ++i) {
					auto p = getGlobalNodeCoord(vertexIDs[i]);
					//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
					grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
					pyra.XYZ(i, 0, p[0]);
					pyra.XYZ(i, 1, p[1]);
					pyra.XYZ(i, 2, p[2]);
					//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
				}

				size_t iNrIps = pyra.IntegrationPoints();
				for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
					pyra.JacobianAtIntegrationPoint(iIp);
					double64 jacdet = pyra.JacobianDeterminant();
					//std::cerr << "jacdet = " << jacdet << '\n';
					if (jacdet <= 0) {
						return false;
					}
				}
				return true;
			}

			size_t EmitPyramid(ColumnCell& cell) {
				std::vector<size_t> ids(&vertexIDs[0], &vertexIDs[5]);
				size_t elid = elementID++;
				plist.emplace(elid, ids);
				fem_types.push_back(ISOPARAMETRIC_LINEAR_PYRAMID);
				return elid;
			}


			bool ConstructTetrahedronOnFace(ColumnCell& cell, size_t face0, size_t face1, size_t face2, size_t apex) {
				vertexIDs[0] = getNodeID(cell, face0);
				vertexIDs[1] = getNodeID(cell, face1);
				vertexIDs[2] = getNodeID(cell, face2);
				vertexIDs[3] = apex;

				//std::cerr << "tetra: \n";
				for (size_t i = 0; i < 4; ++i) {
					auto p = getGlobalNodeCoord(vertexIDs[i]);
					//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
					grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
					tetra.XYZ(i, 0, p[0]);
					tetra.XYZ(i, 1, p[1]);
					tetra.XYZ(i, 2, p[2]);
					//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
				}

				size_t iNrIps = tetra.IntegrationPoints();
				for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
					tetra.JacobianAtIntegrationPoint(iIp);
					double64 jacdet = tetra.JacobianDeterminant();
					//std::cerr << "jacdet = " << jacdet << '\n';
					if (jacdet <= 0) {
						return false;
					}
				}
				return true;
			}

			size_t EmitTetrahedron(ColumnCell& cell) {
				std::vector<size_t> ids(&vertexIDs[0], &vertexIDs[4]);
				size_t elid = elementID++;
				plist.emplace(elid, ids);
				fem_types.push_back(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
				return elid;
			}


			bool ConstructPrism(ColumnCell& cell, size_t face0, size_t face1, size_t face2, size_t face3, size_t face4, size_t face5) {
				vertexIDs[0] = getNodeID(cell, face0);
				vertexIDs[1] = getNodeID(cell, face1);
				vertexIDs[2] = getNodeID(cell, face2);
				vertexIDs[3] = getNodeID(cell, face3);
				vertexIDs[4] = getNodeID(cell, face4);
				vertexIDs[5] = getNodeID(cell, face5);

				//std::cerr << "prism: \n";
				for (size_t i = 0; i < 6; ++i) {
					auto p = getGlobalNodeCoord(vertexIDs[i]);
					grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
					prism.XYZ(i, 0, p[0]);
					prism.XYZ(i, 1, p[1]);
					prism.XYZ(i, 2, p[2]);
					//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
				}

				size_t iNrIps = prism.IntegrationPoints();
				for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
					prism.JacobianAtIntegrationPoint(iIp);
					double64 jacdet = prism.JacobianDeterminant();
					//std::cerr << "jacdet = " << jacdet << '\n';
					if (jacdet <= 0) {
						return false;
					}
				}
				return true;
			}

			size_t EmitPrism(ColumnCell& cell) {
				std::vector<size_t> ids(&vertexIDs[0], &vertexIDs[6]);
				size_t elid = elementID++;
				plist.emplace(elid, ids);
				fem_types.push_back(ISOPARAMETRIC_LINEAR_PRISM);
				return elid;
			}


			// return a list of globalNodeID from local vertex index
			std::vector<size_t> getGlobalIDList(ColumnCell& cell, size_t size, const size_t* vertexIDs) {
				std::vector<size_t> nodeIDs;
				nodeIDs.reserve(size);
				for (int i = 0; i < size; ++i) {
					nodeIDs.push_back(getNodeID(cell, vertexIDs[i]));
				}
				return nodeIDs;
			}

			/// degenerates to one prism

			vector<size_t> degenerateToOnePrismAtEdge23(ColumnCell& cell) {    // 034 125
				static const size_t vertexIDs[6] = { 0, 3, 4, 1, 2, 5 };
				return getGlobalIDList(cell, 6, vertexIDs);
			}

			vector<size_t> degenerateToOnePrismAtEdge30(ColumnCell& cell) { //051 362
				static const size_t vertexIDs[6] = { 0, 5, 1, 3, 6, 2 };
				return getGlobalIDList(cell, 6, vertexIDs);
			}

			vector<vector<size_t>> degenerateToTwoTetrahedrasAtEdge02(ColumnCell& cell) {    // 1010
				vector<vector<size_t>> nodeLists;          // tetra 0125 0237
				nodeLists.reserve(2);

				static const size_t vertexIDs1[4] = { 0, 1, 2, 5 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

				static const size_t vertexIDs2[4] = { 0, 2, 3, 7 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

				return nodeLists;
			}

			vector<vector<size_t>> degenerateToTwoTetrahedraAtEdge13(ColumnCell& cell) {    // 0101
				vector<vector<size_t>> nodeLists;          // tetra 1304 1326
				nodeLists.reserve(2);

				static const size_t vertexIDs1[4] = { 1, 3, 0, 4 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

				static const size_t vertexIDs2[4] = { 1, 3, 2, 6 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

				return nodeLists;
			}

			vector<vector<size_t>> splitToFiveTetrahedrasAtTwoEdges0257(ColumnCell& cell) {
				vector<vector<size_t>> nodeLists;                // 0457 0125 0237 0257 2756
				nodeLists.reserve(5);

				static const size_t vertexIDs1[4] = { 0, 4, 5, 7 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

				static const size_t vertexIDs2[4] = { 0, 1, 2, 5 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

				static const size_t vertexIDs3[4] = { 0, 2, 3, 7 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs3));

				static const size_t vertexIDs4[4] = { 0, 2, 5, 7 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs4));

				static const size_t vertexIDs5[4] = { 2, 7, 5, 6 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs5));

				return nodeLists;
			}


			vector<vector<size_t>> splitToFiveTetrahedrasAtTwoEdges1347(ColumnCell& cell) {
				vector<vector<size_t>> nodeLists;                // 0134 1456 1236 1346 3467
				nodeLists.reserve(5);

				static const size_t vertexIDs1[4] = { 0, 1, 3, 4 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

				static const size_t vertexIDs2[4] = { 1, 4, 5, 6 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

				static const size_t vertexIDs3[4] = { 1, 2, 3, 6 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs3));

				static const size_t vertexIDs4[4] = { 1, 3, 4, 6 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs4));

				static const size_t vertexIDs5[4] = { 3, 4, 6, 7 };
				nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs5));

				return nodeLists;
			}



			FACE_TYPE getShapeOfTopFace(ColumnCell& cell) {
				// using reservor coord, smaller z -> higher

				switch (cell.classification) {
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_021_023: {
					return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_201_203: {
					return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_130_132: {
					return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_310_312: {
					return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_021_023: {
					return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_130_132: {
					return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_0: {
					if (cell.z[0][0] < cell.z[0][1]) {
						return FACE_TYPE::SPLIT_13_OR_57;
					}
					else return FACE_TYPE::FULL_QUAD;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_1: {
					if (cell.z[1][0] < cell.z[1][1]) {
						return FACE_TYPE::SPLIT_02_OR_46;
					}
					else return FACE_TYPE::FULL_QUAD;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_2: {
					if (cell.z[2][0] < cell.z[2][1]) {
						return FACE_TYPE::SPLIT_13_OR_57;
					}
					else return FACE_TYPE::FULL_QUAD;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3: {
					if (cell.z[3][0] < cell.z[3][1]) {
						return FACE_TYPE::SPLIT_02_OR_46;
					}
					else return FACE_TYPE::FULL_QUAD;
				}
				default:
					return FACE_TYPE::FULL_QUAD;
				}
			}

			FACE_TYPE getShapeOfBottomFace(ColumnCell& cell) {
				switch (cell.classification) {
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_021_023: {
					return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_201_203: {
					return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_130_132: {
					return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_310_312: {
					return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_021_023: {
					return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_130_132: {
					return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_0: {
					if (cell.z[0][0] < cell.z[0][1]) {
						return FACE_TYPE::FULL_QUAD;
					}
					else return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_1: {
					if (cell.z[1][0] < cell.z[1][1]) {
						return FACE_TYPE::FULL_QUAD;
					}
					else return FACE_TYPE::SPLIT_13_OR_57;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_2: {
					if (cell.z[2][0] < cell.z[2][1]) {
						return FACE_TYPE::FULL_QUAD;
					}
					else return FACE_TYPE::SPLIT_02_OR_46;
				}
				case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3: {
					if (cell.z[3][0] < cell.z[3][1]) {
						return FACE_TYPE::FULL_QUAD;
					}
					else return FACE_TYPE::SPLIT_13_OR_57;
				}
				default:
					return FACE_TYPE::FULL_QUAD;
				}
			}
		};

		// CORNER POINT GRID

		CornerPointGrid_UoM::CornerPointGrid_UoM()
			: NX_(0),
			NY_(0),
			NZ_(0),
			NX_x_NY_(0)
		{
		}

		CornerPointGrid_UoM::~CornerPointGrid_UoM()
		{
		}


		void CornerPointGrid_UoM::AssignDimensions(size_t nx, size_t ny, size_t nz)
		{
			NX_ = nx;
			NY_ = ny;
			NZ_ = nz;
		}


		void CornerPointGrid_UoM
			::InitializeGridSpecs()
		{
			std::cout << "InitializeGridSpecs\n";
			NX_x_NY_ = NX_ * NY_;
			elements_ = NX_x_NY_ * NZ_;
			if (cell_activity_.empty()) {
				cell_activity_.resize(elements_, 1);
				active_elements_ = elements_;
			}
			else {
				active_elements_ = 0;
				for (auto activity : cell_activity_) {
					if (activity) ++active_elements_;
				}
			}
		}

		/** JC: Why is it needed???
		x stays, y = -z and z=y.
		*/
		void CornerPointGrid_UoM
			::ConvertFromReservoirToCSMPcoordinateSystem(csmp::Point<3U>& pt) const
		{
			const double y(pt[1U]);
			pt[0U] = pt[0U];
			pt[1U] = -pt[2U];
			pt[2U] = y;
		};


		uint8_t&
			CornerPointGrid_UoM::CellActivity(size_t i, size_t j, size_t k)
		{
			assert(i < NX_ && j < NY_ && k < NZ_);
			return cell_activity_[i + j * NX_ + k * NX_x_NY_];
		}


		std::vector<uint8_t>&
			CornerPointGrid_UoM::GetCellActivity()
		{
			return cell_activity_;
		}



		/**
		Builds the pillars and columns.
		*/
		void CornerPointGrid_UoM::ConstructPillarsAndColumns(const std::vector<double64>& zcorn)
		{
			const size_t NXxNY = NX_ * NY_;
			size_t node_count = 0;

			{
				// 1. Build pillars
				std::vector<double64> zcoord;
				zcoord.reserve(4 * (NZ_ + 1));
				for (size_t i = 0; i <= NX_; ++i)
				{
					for (size_t j = 0; j <= NY_; ++j)
					{
						for (size_t k = 0; k < NZ_; ++k) {
							if (i < NX_ && j < NY_ && CellActivity(i, j, k)) {
								zcoord.push_back(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 0]); // t_nw
								zcoord.push_back(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 4]); // b_nw
							}

							if (i > 0 && j < NY_ && CellActivity(i - 1, j, k)) {
								zcoord.push_back(zcorn[((i - 1) + j * NX_ + k * NXxNY) * 8 + 1]); // t_ne
								zcoord.push_back(zcorn[((i - 1) + j * NX_ + k * NXxNY) * 8 + 5]); // b_ne
							}

							if (i < NX_ && j > 0 && CellActivity(i, j - 1, k)) {
								zcoord.push_back(zcorn[(i + (j - 1) * NX_ + k * NXxNY) * 8 + 2]); // t_sw
								zcoord.push_back(zcorn[(i + (j - 1) * NX_ + k * NXxNY) * 8 + 6]); // b_sw
							}

							if (i > 0 && j > 0 && CellActivity(i - 1, j - 1, k)) {
								zcoord.push_back(zcorn[((i - 1) + (j - 1) * NX_ + k * NXxNY) * 8 + 3]); // t_se
								zcoord.push_back(zcorn[((i - 1) + (j - 1) * NX_ + k * NXxNY) * 8 + 7]); // b_se
							}
						}

						sortAndUnique(zcoord);
						auto& p = (*this)(i, j);
						p.SetZCoords(zcoord);
						p.SetFirstNodeNum(node_count);
						node_count += zcoord.size();
						zcoord.clear();
					}
				}
			}
			ordinaryNodes_ = node_count;
			std::cerr << " " << node_count << " unique nodes detected\n";

			// 2. Build columns
			size_t cell_count = 0;
			size_t fully_degenerate_cells = 0;
			{
				std::vector<ColumnCell> cells;
				cells.reserve(NZ_);
				for (size_t i = 0; i < NX_; ++i)
				{
					for (size_t j = 0; j < NY_; ++j)
					{
						std::pair<size_t, size_t> index(i, j);
						Pillar& p0 = (*this)(i + 0, j + 0);
						Pillar& p1 = (*this)(i + 1, j + 0);
						Pillar& p2 = (*this)(i + 1, j + 1);
						Pillar& p3 = (*this)(i + 0, j + 1);

						for (size_t k = 0; k < NZ_; ++k)
						{
							if (!CellActivity(i, j, k)) {
								continue;
							}

							ColumnCell cell;
							cell.k = k;
							cell.z[0][0] = p0.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 0]); // t_nw
							cell.z[0][1] = p0.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 4]); // b_nw
							cell.z[1][0] = p1.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 1]); // t_ne
							cell.z[1][1] = p1.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 5]); // b_ne
							cell.z[2][0] = p2.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 3]); // t_se
							cell.z[2][1] = p2.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 7]); // b_se
							cell.z[3][0] = p3.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 2]); // t_sw
							cell.z[3][1] = p3.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 6]); // b_sw
							
						    // If all four corners are degenerate, the cell is fully degenerate.

							if (cell.z[0][0] == cell.z[0][1]
								&& cell.z[1][0] == cell.z[1][1]
								&& cell.z[2][0] == cell.z[2][1]
								&& cell.z[3][0] == cell.z[3][1]) {
								++fully_degenerate_cells;
								continue;
							}
							++cell_count;
							cells.push_back(cell);
						}

						columns_.emplace(index, cells);
						cells.clear();
					}
				}
			}
			std::cerr << " " << cell_count << " unique active cells detected\n";
			if (fully_degenerate_cells > 0) {
				std::cerr << " " << fully_degenerate_cells << " fully degenerate cells detected\n";
			}

			std::cerr << "Pillars and columns built\n";

		} // end ConstructPillarsAndColumns




		void
			CornerPointGrid_UoM::ConstructFiniteElementsFromColumns(VSet<3U>& vset)
		{
			vset.HybridElementTypeMesh(true);

			// 1. Classify the cells

			size_t skewCells = 0;
			for (auto& index_column : columns_) {
				auto& column = index_column.second;
				size_t i = index_column.first.first;
				size_t j = index_column.first.second;

				const size_t iNrCells = column.cells_.size();
				for (size_t iCell = 0; iCell < iNrCells; ++iCell) {
					auto& cell = column.cells_[iCell];

					Pillar& p0 = (*this)(i + 0, j + 0); //nw
					Pillar& p1 = (*this)(i + 1, j + 0); //ne
					Pillar& p2 = (*this)(i + 1, j + 1); //se
					Pillar& p3 = (*this)(i + 0, j + 1); //sw

					double64 z[4][2];
					z[0][0] = p0.GetZCoord(cell.z[0][0]);
					z[0][1] = p0.GetZCoord(cell.z[0][1]);
					z[1][0] = p1.GetZCoord(cell.z[1][0]);
					z[1][1] = p1.GetZCoord(cell.z[1][1]);
					z[2][0] = p2.GetZCoord(cell.z[2][0]);
					z[2][1] = p2.GetZCoord(cell.z[2][1]);
					z[3][0] = p3.GetZCoord(cell.z[3][0]);
					z[3][1] = p3.GetZCoord(cell.z[3][1]);

					//AB: Classify the degeneracy
					uint8_t classification = 0;
					for (size_t v = 0; v < 4; ++v) {
						if (cell.z[v][0] == cell.z[v][1]) {
							classification |= (1 << (3 - v));
						}
					}
					cell.classification = static_cast<ECLIPSE_CELL_CLASSIFICATION>(classification);

					//std::cout << "p0 t: "; p0.GetPoint(cell.z[0][0]).Out();
					//std::cout << "p0 b: "; p0.GetPoint(cell.z[0][1]).Out();
					//std::cout << "p1 t: "; p1.GetPoint(cell.z[1][0]).Out();
					//std::cout << "p1 b: "; p1.GetPoint(cell.z[1][1]).Out();
					//std::cout << "p2 t: "; p2.GetPoint(cell.z[2][0]).Out();
					//std::cout << "p2 b: "; p2.GetPoint(cell.z[2][1]).Out();
					//std::cout << "p3 t: "; p3.GetPoint(cell.z[3][0]).Out();
					//std::cout << "p3 b: "; p3.GetPoint(cell.z[3][1]).Out();
					//std::cout << "cell.classification: " << int(classification);
					//std::cout << "\n";
				}
			}


			// 2. Constructing Elements
			CellGenerator generator(*this);
			{
				generator.ordinaryNodes.reserve(ordinaryNodes_);
				for (size_t i = 0; i <= NX_; ++i) {
					for (size_t j = 0; j <= NY_; ++j) {
						Pillar& pillar = (*this)(i, j);
						for (size_t k = 0; k < pillar.GetNumPoints(); ++k) {
							generator.ordinaryNodes.push_back(pillar.GetPoint(k));
						}
					}
				}
			}

			size_t badHexahedra = 0;
			size_t badPyramids = 0;
			size_t badTetrahedra = 0;
			size_t badPrisms = 0;


			//Degeneration process starts!
			//Degeneration process starts!
			for (auto column = columns_.begin(); column != columns_.end(); column++) {
				Column& Col = column->second;
				for (size_t k = 0; k < Col.cells_.size(); k++) {
					ColumnCell&  cell = Col.cells_[k];
					pair<size_t, size_t> index = column->first;
					size_t i = index.first;
					size_t j = index.second;

					generator.setIJ(i, j);

					ECLIPSE_CELL_CLASSIFICATION cellType = cell.classification;

					/// referred to the classifications of the cell above and the cell beneath  
					ColumnCell* cellAbove = nullptr;
					ColumnCell* cellBeneath = nullptr;

					if (k > 0) {
						cellAbove = &Col.cells_[k - 1];
					}

					if (k < Col.cells_.size() - 1) {
						cellBeneath = &Col.cells_[k + 1];
					}

					assert(generator.elementID == generator.plist.size());
					assert(generator.elementID == generator.fem_types.size());

					switch (cellType) {
					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_HEXAHEDRON: {        // 0000
						bool faceAboveIsQuad = !cellAbove || generator.getShapeOfBottomFace(*cellAbove) == FACE_TYPE::FULL_QUAD;
						bool faceBeneathIsQuad = !cellBeneath || generator.getShapeOfTopFace(*cellBeneath) == FACE_TYPE::FULL_QUAD;

						if (faceAboveIsQuad) {
							if (faceBeneathIsQuad) {
								if (generator.ConstructHexahedron(cell)) {
									addElementToMap(i, j, k, generator.EmitHexahedron(cell));
								}
								else {
									++badHexahedra;
									continue;
								}
							}  // null - hex - null
							else {
								switch (generator.getShapeOfTopFace(*cellBeneath)) {
								case FACE_TYPE::SPLIT_02_OR_46:
								{
									size_t centroid = generator.generateCellCentroid(cell);

									// Top face
									if (generator.ConstructPyramidOnFace(cell, 0, 1, 2, 3, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
									}

									// Left face
									if (generator.ConstructPyramidOnFace(cell, 0, 4, 5, 1, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
									}

									// Right face
									if (generator.ConstructPyramidOnFace(cell, 3, 2, 6, 7, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
									}

									// Front face
									if (generator.ConstructPyramidOnFace(cell, 0, 3, 7, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
									}

									// Back face
									if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
									}

									// Bottom faces
									if (generator.ConstructTetrahedronOnFace(cell, 6, 5, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									if (generator.ConstructTetrahedronOnFace(cell, 6, 4, 7, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									break;
								}

								case FACE_TYPE::SPLIT_13_OR_57:
								{
									size_t centroid = generator.generateCellCentroid(cell);

									// Top face
									if (generator.ConstructPyramidOnFace(cell, 0, 1, 2, 3, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Left face
									if (generator.ConstructPyramidOnFace(cell, 0, 4, 5, 1, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Right face
									if (generator.ConstructPyramidOnFace(cell, 3, 2, 6, 7, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Front face
									if (generator.ConstructPyramidOnFace(cell, 0, 3, 7, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Back face
									if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Bottom faces
									if (generator.ConstructTetrahedronOnFace(cell, 7, 6, 5, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									if (generator.ConstructTetrahedronOnFace(cell, 7, 5, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									break;
								}

								default:
								{
									throw csmp::Exception(FATAL_ERROR,
										"CornerPointGrid::ConstructFiniteElementsFromColumns",
										"Hexahedron with unknown bottom face");
								}
								}
							}
						}
						else {
							if (faceBeneathIsQuad) {
								switch (generator.getShapeOfBottomFace(*cellAbove)) {
								case FACE_TYPE::SPLIT_02_OR_46:
								{
									size_t centroid = generator.generateCellCentroid(cell);

									// Left face
									if (generator.ConstructPyramidOnFace(cell, 0, 4, 5, 1, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Right face
									if (generator.ConstructPyramidOnFace(cell, 3, 2, 6, 7, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Front face
									if (generator.ConstructPyramidOnFace(cell, 0, 3, 7, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Back face
									if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Bottom face
									if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Top faces
									if (generator.ConstructTetrahedronOnFace(cell, 0, 1, 2, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									if (generator.ConstructTetrahedronOnFace(cell, 0, 2, 3, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									break;
								}

								case FACE_TYPE::SPLIT_13_OR_57:
								{
									size_t centroid = generator.generateCellCentroid(cell);

									// Left face
									if (generator.ConstructPyramidOnFace(cell, 0, 4, 5, 1, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Right face
									if (generator.ConstructPyramidOnFace(cell, 3, 2, 6, 7, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Front face
									if (generator.ConstructPyramidOnFace(cell, 0, 3, 7, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Back face
									if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Bottom face
									if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, centroid)) {
										addElementToMap(i, j, k, generator.EmitPyramid(cell));
									}
									else {
										++badPyramids;
										continue;
									}

									// Top faces
									if (generator.ConstructTetrahedronOnFace(cell, 1, 2, 3, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									if (generator.ConstructTetrahedronOnFace(cell, 1, 3, 0, centroid)) {
										addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
									}
									else {
										++badTetrahedra;
									}

									break;
								}

								default:;
								}
							}
							else {
								auto topShape = generator.getShapeOfBottomFace(*cellAbove);
								auto bottomShape = generator.getShapeOfTopFace(*cellBeneath);
								switch (topShape) {
								case FACE_TYPE::SPLIT_02_OR_46:
								{
									switch (bottomShape) {
									case FACE_TYPE::SPLIT_02_OR_46:
									{
										if (generator.ConstructPrism(cell, 6, 5, 4, 2, 1, 0)) {
											addElementToMap(i, j, k, generator.EmitPrism(cell));
										}
										else {
											++badPrisms;
										}

										if (generator.ConstructPrism(cell, 0, 2, 3, 4, 6, 7)) {
											addElementToMap(i, j, k, generator.EmitPrism(cell));
										}
										else {
											++badPrisms;
										}

										break;
									}

									case FACE_TYPE::SPLIT_13_OR_57:
									{
										size_t centroid = generator.generateCellCentroid(cell);

										if (generator.ConstructPyramidOnFace(cell, 0, 4, 5, 1, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
											continue;
										}

										// Right face
										if (generator.ConstructPyramidOnFace(cell, 3, 2, 6, 7, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
											continue;
										}

										// Front face
										if (generator.ConstructPyramidOnFace(cell, 0, 3, 7, 4, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
											continue;
										}

										// Back face
										if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
											continue;
										}

										// Top faces
										if (generator.ConstructTetrahedronOnFace(cell, 0, 1, 2, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										if (generator.ConstructTetrahedronOnFace(cell, 0, 2, 3, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										// Bottom faces
										if (generator.ConstructTetrahedronOnFace(cell, 7, 6, 5, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										if (generator.ConstructTetrahedronOnFace(cell, 7, 5, 4, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										break;
									}

									default:;
									}
									break;
								}

								case FACE_TYPE::SPLIT_13_OR_57:
								{
									switch (bottomShape) {
									case FACE_TYPE::SPLIT_02_OR_46:
									{
										size_t centroid = generator.generateCellCentroid(cell);

										// Left face
										if (generator.ConstructPyramidOnFace(cell, 0, 4, 5, 1, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
										}

										// Right face
										if (generator.ConstructPyramidOnFace(cell, 3, 2, 6, 7, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
										}

										// Front face
										if (generator.ConstructPyramidOnFace(cell, 0, 3, 7, 4, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
										}

										// Back face
										if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, centroid)) {
											addElementToMap(i, j, k, generator.EmitPyramid(cell));
										}
										else {
											++badPyramids;
										}

										// Top faces
										if (generator.ConstructTetrahedronOnFace(cell, 1, 2, 3, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										if (generator.ConstructTetrahedronOnFace(cell, 1, 3, 0, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										// Bottom faces
										if (generator.ConstructTetrahedronOnFace(cell, 6, 5, 4, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										if (generator.ConstructTetrahedronOnFace(cell, 6, 4, 7, centroid)) {
											addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
										}
										else {
											++badTetrahedra;
										}

										break;
									}

									case FACE_TYPE::SPLIT_13_OR_57:
									{
										if (generator.ConstructPrism(cell, 5, 4, 7, 1, 0, 3)) {
											addElementToMap(i, j, k, generator.EmitPrism(cell));
										}
										else {
											++badPrisms;
										}

										if (generator.ConstructPrism(cell, 7, 6, 5, 3, 2, 1)) {
											addElementToMap(i, j, k, generator.EmitPrism(cell));
										}
										else {
											++badPrisms;
										}
										break;
									}

									default:;
									}
									break;
								}
								default:;
								}
							}
						}
						break;
					}


					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_310_312: {      // 0001
						if (generator.ConstructPyramidOnFace(cell, 4, 5, 1, 0, generator.getNodeID(cell, 3))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						if (generator.ConstructPyramidOnFace(cell, 1, 5, 6, 2, generator.getNodeID(cell, 3))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_201_203: {      // 0010
						if (generator.ConstructPyramidOnFace(cell, 4, 5, 1, 0, generator.getNodeID(cell, 2))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						if (generator.ConstructPyramidOnFace(cell, 3, 7, 4, 0, generator.getNodeID(cell, 2))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_23: {          // 0011
						generator.plist.emplace(generator.elementID, generator.degenerateToOnePrismAtEdge23(cell));
						addElementToMap(i, j, k, generator.elementID);
						generator.fem_types.push_back(ISOPARAMETRIC_LINEAR_PRISM);
						++generator.elementID;
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_130_132: {      // 0100
						if (generator.ConstructPyramidOnFace(cell, 3, 7, 4, 0, generator.getNodeID(cell, 1))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						if (generator.ConstructPyramidOnFace(cell, 2, 6, 7, 3, generator.getNodeID(cell, 1))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_130_132: {    // 0101
						auto elementList = generator.degenerateToTwoTetrahedraAtEdge13(cell);
						generator.plist.emplace(generator.elementID, elementList[0]);
						addElementToMap(i, j, k, generator.elementID);
						generator.fem_types.push_back(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
						++generator.elementID;
						generator.plist.emplace(generator.elementID, elementList[1]);
						addElementToMap(i, j, k, generator.elementID);
						generator.fem_types.push_back(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
						++generator.elementID;
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_12: {          // 0111
						if (generator.ConstructPrism(cell, 0, 4, 1, 3, 7, 2)) {
							addElementToMap(i, j, k, generator.EmitPrism(cell));
						}
						else {
							++badPrisms;
						}
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_0: {          // 0111
																						 //if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 0))) {
																						 // addElementToMap(i, j, k, generator.EmitPyramid(cell));
																						 //}

						if (generator.ConstructTetrahedronOnFace(cell, 7, 5, 4, generator.getNodeID(cell, 0))) {
							addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
						}
						else {
							++badPyramids;
						}
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMIDS_021_023: {      // 1000
						if (generator.ConstructPyramidOnFace(cell, 5, 6, 2, 1, generator.getNodeID(cell, 0))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						if (generator.ConstructPyramidOnFace(cell, 7, 3, 2, 6, generator.getNodeID(cell, 0))) {
							addElementToMap(i, j, k, generator.EmitPyramid(cell));
						}
						else {
							++badPyramids;
						}

						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_03: {                    // 1001
						generator.plist.emplace(generator.elementID, generator.degenerateToOnePrismAtEdge30(cell));
						addElementToMap(i, j, k, generator.elementID);
						generator.fem_types.push_back(ISOPARAMETRIC_LINEAR_PRISM);
						++generator.elementID;
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_TETRAHEDRONS_021_023: {              // 1010
						auto elementList = generator.degenerateToTwoTetrahedrasAtEdge02(cell);
						generator.plist.emplace(generator.elementID, elementList[0]);
						addElementToMap(i, j, k, generator.elementID);
						generator.fem_types.push_back(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
						++generator.elementID;
						generator.plist.emplace(generator.elementID, elementList[1]);
						addElementToMap(i, j, k, generator.elementID);
						generator.fem_types.push_back(ISOPARAMETRIC_LINEAR_TETRAHEDRON);
						++generator.elementID;
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_1: {                    // 1011
																								   /*if (generator.ConstructHexahedron(cell, 7, 6, 5, 4, generator.getNodeID(cell, 1))) {
																								   addElementToMap(i, j, k, generator.EmitPyramid(cell));
																								   }*/
						if (generator.ConstructTetrahedronOnFace(cell, 7, 5, 4, generator.getNodeID(cell, 1))) {
							addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
						}
						else {
							++badPyramids;
						}
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PRISM_01: {                    // 1100
						if (generator.ConstructPrism(cell, 0, 3, 7, 1, 2, 6)) {
							addElementToMap(i, j, k, generator.EmitPrism(cell));
						}
						else {
							++badPrisms;
						}
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_2: {                    // 1101
																								   /*if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 2))) {
																								   addElementToMap(i, j, k, generator.EmitPyramid(cell));
																								   }*/
						if (generator.ConstructTetrahedronOnFace(cell, 7, 6, 5, generator.getNodeID(cell, 2))) {
							addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
						}
						else {
							++badPyramids;
						}
						break;
					}

					case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3: {                    // 1110
																								   /*if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 3))) {
																								   addElementToMap(i, j, k, generator.EmitPyramid(cell));
																								   }*/
						if (generator.ConstructTetrahedronOnFace(cell, 7, 6, 4, generator.getNodeID(cell, 3))) {
							addElementToMap(i, j, k, generator.EmitTetrahedron(cell));
						}
						else {
							++badPyramids;
						}
						break;
					}
					}
				}
			}

			std::cerr << "Bad hexahedra: " << badHexahedra << '\n';
			std::cerr << "Bad pyramids: " << badPyramids << '\n';
			std::cerr << "Bad tetrahedra: " << badTetrahedra << '\n';
			std::cerr << "Bad prisms: " << badPrisms << '\n';

			// 3. store node coordinates
			{
				deque<double64> x, y, z;
				for (size_t i = 0; i <= NX_; ++i) {
					for (size_t j = 0; j <= NY_; ++j) {
						Pillar& pillar = (*this)(i, j);
						for (size_t k = 0; k < pillar.GetNumPoints(); ++k) {
							csmp::Point<3u> point = pillar.GetPoint(k);
							// convert to CSMP coordinate
							ConvertFromReservoirToCSMPcoordinateSystem(point);
							x.push_back(point[0]);
							y.push_back(point[1]);
							z.push_back(point[2]);
						}
					}
				}
				for (auto p : generator.extraNodes) {
					ConvertFromReservoirToCSMPcoordinateSystem(p);
					x.push_back(p[0]);
					y.push_back(p[1]);
					z.push_back(p[2]);
				}
				vset.AddXYZ(x, y, z);
			}

			// 4. Set up the rest of the vset
			// JC: There is no neighbor information in the Eclipse data(*.grdecl). The information will be created later.
			vset.RemovePfverts();
			//vset.ResizePfverts(generator.fem_types.size());
			vset.ResizePlist(generator.plist.size());
			vset.AddPlist(generator.plist.begin(), generator.plist.end());
			vset.AddElementTypes(generator.fem_types.begin(), generator.fem_types.end());
		}

		/**
		MASTER METHOD for the creation of corner point grids
		from the original cell-centered grids read by the EclipseInterface.

		The PolygonGridManager is used to construct the pillars.

		@todo this is confused! - why should the corner point grid be responsible for building the CSMP model? - code should be in the EclipseModel

		@todo regions are not recognised properly.
		*/
		void CornerPointGrid_UoM::CreateModel(const std::string&     model_name,
			csmp::VSet<3U>&       vset,
			csmp::ModelTopology&   model_topology,
			const std::vector<double64>& zcorn,
			std::set<std::string>& regions,
			std::set<std::string>& faults,
			std::set<std::string>& wells,
			bool tetra_mesh, bool exclude_inactive_cells)
		{

			csmp::ErrorHandler& csmp_error(csmp::ErrorHandler::Instance());


			// temporary data
			deque<double64>              x, y, z;        // node x, y, z coordinates
			map<size_t, vector<size_t> >  plist;        // element list
			map<size_t, vector<long64> >  pfverts;        // neighbor list


														  // 1. Initialise grid specs
			InitializeGridSpecs();

			// 2. Construct pillars and columns
			ConstructPillarsAndColumns(zcorn);

			// 3. Construct FEs from columns and add to Vset
			ConstructFiniteElementsFromColumns(vset);

		}

		void CornerPointGrid_UoM::addElementToMap(size_t i, size_t j, size_t k, size_t elementID) {
			this->elementMap.emplace(ijk(i, j, k), elementID);
		};

		template<class VarType>
		void CornerPointGrid_UoM::WritePropertyToVSet(csmp::VSet<3U>&             vset,
			const std::vector<VarType>& prop_data,
			const std::string&          prop_name,
			const csmp::PLACEMENT&      prop_place) const
		{
			/// correcting data cell id's due to existance of embedded cells
			VarType var;
			var = 0.;
			size_t current_id = 0;
			std::vector<VarType> cell_data(vset.Elements(), var);
			for (auto& entry : elementMap) {
				auto coord = entry.first;
				size_t idx = coord.i + coord.j * NX_ + coord.k * NX_x_NY_;
				cell_data[entry.second] = prop_data[idx];
			}
			if (prop_place == csmp::NODE) {
				// Add nodal data to vset
				std::vector<VarType> nodal_data;
				// TODO: why node property?
				extrapolateElementToNodeProperty<3U, VarType>(vset, cell_data, nodal_data);
				//csmp::FEM_Data<VarType> property_values( prop_place, nodal_data );
				const size_t array_length = (var.Size() > 9U) ? var.Size() : 0U;
				assert(array_length == 1);
				PropertyData property_values(prop_place, VarType::VariableType, 3U, 0U);
				property_values.Reserve(nodal_data.size());
				for (const auto& it : nodal_data) pushBack(property_values, it);
				vset.AddData(prop_name.c_str(), property_values);
			}
			else {
				// Add cell data to vset
				// csmp::FEM_Data<VarType> property_values( prop_place, cell_data );
				PropertyData property_values(prop_place, VarType::VariableType, 3U, 0U);
				property_values.Reserve(cell_data.size());
				for (const auto& it : cell_data) pushBack(property_values, it);
				vset.AddData(prop_name.c_str(), property_values);
			}
		}

		template void CornerPointGrid_UoM::WritePropertyToVSet(csmp::VSet<3U>&, const std::vector<ScalarVariable>&, const std::string&, const csmp::PLACEMENT&) const;
		template void CornerPointGrid_UoM::WritePropertyToVSet(csmp::VSet<3U>&, const std::vector<VectorVariable<3U> >&, const std::string&, const csmp::PLACEMENT&) const;
		template void CornerPointGrid_UoM::WritePropertyToVSet(csmp::VSet<3U>&, const std::vector<TensorVariable<3U> >&, const std::string&, const csmp::PLACEMENT&) const;
		template void CornerPointGrid_UoM::WritePropertyToVSet(csmp::VSet<3U>&, const std::vector<ArrayVariable>&, const std::string&, const csmp::PLACEMENT&) const;
		template void CornerPointGrid_UoM::WritePropertyToVSet(csmp::VSet<3U>&, const std::vector<FlaggedArrayVariable>&, const std::string&, const csmp::PLACEMENT&) const;


		void CornerPointGrid_UoM::Resize(size_t i_pillar_max, size_t j_pillar_max)
		{
			pillars_.resize(i_pillar_max * j_pillar_max);
		}


		/** accessing the contained pillars
		*/
		Pillar & CornerPointGrid_UoM::operator()(size_t i, size_t j)
		{
			assert(i <= NX_);
			assert(j <= NY_);

			return pillars_[j * (NX_ + 1) + i];
		}

		// CELL CENTERED GRID

		CellCenteredGrid::CellCenteredGrid()
			:dx_(3U)
		{
		}

		CellCenteredGrid::~CellCenteredGrid()
		{
		}

		void CellCenteredGrid::Clear()
		{
			/// grid
			tops_.clear();
			dx_.clear();
			dy_.clear();
			dz_.clear();
		}



		size_t CellCenteredGrid::GetNumCells() const
		{
			return NX_*NY_*NZ_;
		}



		std::vector<csmp::ScalarVariable>& CellCenteredGrid::GetCellDepths()
		{
			return tops_;
		}



		std::vector<ScalarVariable>& CellCenteredGrid::GetCellSizes(size_t i)
		{
			if (i == 0)
				return dx_;
			if (i == 1)
				return dy_;
			return dz_;
		}



		void CellCenteredGrid::AssignDimensionX(size_t NX)
		{
			NX_ = NX;
		}



		void CellCenteredGrid::AssignDimensionY(size_t NY)
		{
			NY_ = NY;
		}



		void CellCenteredGrid::AssignDimensionZ(size_t NZ)
		{
			NZ_ = NZ;
		}

		// WELLS
		/// add well path based on symmetry assumption ( neighbouring cell defines the direction )
		/// by default well is assumed to be vertical
		void addWellPath(size_t NX, size_t NY, size_t NZ,
			const std::string& well_name,
			const std::vector<ijk>& cell_ids,
			std::map<std::string, EclipseWellPath>& well_path)
		{
			const size_t NX_x_NY(NX*NY);
			std::map<int64_t, CORNER_POINT_CELL_FACE_INDEX> face_map;
			face_map.insert(std::make_pair(-1, CORNER_POINT_CELL_FACE_Xminus));
			face_map.insert(std::make_pair(+1, CORNER_POINT_CELL_FACE_Xplus));
			face_map.insert(std::make_pair(-NX, CORNER_POINT_CELL_FACE_Yminus));
			face_map.insert(std::make_pair(+NX, CORNER_POINT_CELL_FACE_Yplus));
			face_map.insert(std::make_pair(-NX_x_NY, CORNER_POINT_CELL_FACE_Zminus));
			face_map.insert(std::make_pair(+NX_x_NY, CORNER_POINT_CELL_FACE_Zplus));

			EclipseWellPath  wpath;
			size_t num_cells(cell_ids.size());
			for (size_t cid = 0; cid<num_cells; ++cid)
			{
				wpath.path.emplace_back(cell_ids[cid], CORNER_POINT_CELL_FACE_Zminus, CORNER_POINT_CELL_FACE_Zplus);
			}
			if (wpath.path.size() > 1)
			{
				size_t nid(wpath.path.size() - num_cells + 1); /// neighbour is a next cell
				for (size_t cid = 0; cid <(num_cells - 1); ++cid, ++nid)
				{
					//              neighbor id          cell id
					int64_t face_id = ((int64_t)wpath.path[nid].cell.i - (int)cell_ids[cid].i)
						+ ((int64_t)wpath.path[nid].cell.j - (int64_t)cell_ids[cid].j) * (int64_t)NX
						+ ((int64_t)wpath.path[nid].cell.k - (int64_t)cell_ids[cid].k) * (int64_t)NX_x_NY;
					if (nid - 1 == 0) wpath.path[nid - 1].from = face_map[-face_id];
					wpath.path[nid - 1].to = face_map[face_id];
					wpath.path[nid].from = face_map[-face_id];
					wpath.path[nid].to = face_map[face_id];
				}
			}
			well_path.emplace(well_name, wpath);
		}

		/// add well path with explicitly specified faces
		void addWellPath(const std::string& well_name,
			const std::vector<size_t>& cell_ids,
			const std::vector<std::pair<size_t, size_t> >& face_ids,
			std::map<std::string, std::vector<std::pair<size_t, std::pair<size_t, size_t> > > >& well_path)
		{
			/// temp data
			std::pair<size_t, size_t> direction;
			std::pair<size_t, std::pair<size_t, size_t> > path;
			size_t cell_id;
			size_t num_cells(cell_ids.size());
			std::vector<std::pair<size_t, std::pair<size_t, size_t> > >   empty_path;
			well_path.insert(std::make_pair(well_name, empty_path));
			std::vector<std::pair<size_t, std::pair<size_t, size_t> > >& wpath(well_path[well_name]);
			for (size_t cid = 0; cid<num_cells; ++cid)
			{
				cell_id = cell_ids[cid];
				direction.first = face_ids[cid].first;
				direction.second = face_ids[cid].second;
				path.first = cell_id;
				path.second = direction;
				wpath.push_back(path);
			}
		}

		/// add well path with explicitly specified faces ( same for all cells )
		void addWellPath(const std::string& well_name,
			const std::vector<size_t>& cell_ids,
			std::pair<size_t, size_t> face_id,
			std::map<std::string, std::vector<std::pair<size_t, std::pair<size_t, size_t> > > >& well_path)
		{
			/// temp data
			std::pair<size_t, size_t> direction;
			std::pair<size_t, std::pair<size_t, size_t> > path;
			size_t cell_id;
			size_t num_cells(cell_ids.size());
			std::vector<std::pair<size_t, std::pair<size_t, size_t> > >   empty_path;
			well_path.insert(std::make_pair(well_name, empty_path));
			std::vector<std::pair<size_t, std::pair<size_t, size_t> > >& wpath(well_path[well_name]);
			for (size_t i = 0; i<num_cells; ++i)
			{
				cell_id = cell_ids[i];
				direction.first = face_id.first;
				direction.second = face_id.second;
				path.first = cell_id;
				path.second = direction;
				wpath.push_back(path);
			}
		}

		/// add well path with explicitly specified faces ( for single cell )
		void addWellPath(const std::string& well_name,
			size_t cell_id,
			std::pair<size_t, size_t> face_id,
			std::map<std::string, std::vector<std::pair<size_t, std::pair<size_t, size_t> > > >& well_path)
		{
			/// temp data
			std::pair<size_t, size_t> direction;
			std::pair<size_t, std::pair<size_t, size_t> > path;
			direction.first = face_id.first;
			direction.second = face_id.second;
			path.first = cell_id;
			path.second = direction;
			std::vector<std::pair<size_t, std::pair<size_t, size_t> > >   empty_path;
			well_path.insert(std::make_pair(well_name, empty_path));
			well_path[well_name].push_back(path);
		}

	} // eclipse

} // end namespace csmp
