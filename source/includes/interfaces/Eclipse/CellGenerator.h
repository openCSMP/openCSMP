#ifndef CELL_GENERATOR_UOM_H
#define CELL_GENERATOR_UOM_H

#include "CornerPointGrid.h"

#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearLineElement.h"

#include "ErrorHandler.h"


namespace csmp {

enum class FACE_TYPE : std::int8_t {

	//       0_________1
	//      /|        /|
	//     3_|_______2 |
	//     | 4-------|-5
	//     |/        |/
	//     7_________6

	FULL_QUAD,
	SPLIT_02_OR_46,
	SPLIT_13_OR_57,
	SPLIT_X
};

class CellGenerator {
public:
	CornerPointGrid& grid;

	size_t elementID = UINT_MAX;
	std::vector<Point<3u>> ordinaryNodes;
	std::deque<Point<3u>> extraNodes;

	std::map<size_t, std::vector<size_t>>  plist;
	std::vector<int8_t> fem_types;

	Pillar* p0;
	Pillar* p1;
	Pillar* p2;
	Pillar* p3;

	CellGenerator(CornerPointGrid& grid) : grid(grid) {}

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


	void setIJ(uint32_t i, uint32_t j) {
		p0 = &grid(i + 0u, j + 0u);
		p1 = &grid(i + 1u, j + 0u);
		p2 = &grid(i + 1u, j + 1u);
		p3 = &grid(i + 0u, j + 1u);
	}


  struct vec3 {
    double x, y, z;

    vec3( double x0, double y0, double z0){
      x = x0; y = y0; z = z0;
    }

    vec3 operator-( vec3 p ) const {
      return vec3{ x - p.x, y - p.y, z - p.z };
    }

    vec3 cross( vec3 p ) const {
      return vec3{
        y * p.z - p.y * z,
        z * p.x - p.z * x,
        x * p.y - p.x * y
      };
    }

    double dot( vec3 p ) const {
      return x * p.x + y * p.y + z * p.z;
    }

    double norm() const {
      return std::sqrt( x*x + y*y + z*z );
    }
  };

  struct face3 {
    std::vector<vec3> v;

    vec3 normal() const {
      assert( v.size() > 2 );
      vec3 dir1 = v[1] - v[0];
      vec3 dir2 = v[2] - v[0];
      vec3 n = dir1.cross( dir2 );
      double d = n.norm();
      return vec3{ n.x / d, n.y / d, n.z / d };
    }
  };


  bool is_inside_polyheadron( vec3 const& cp, std::vector<vec3> const& pts )
  {
    std::vector<face3> polyhedron{ // faces with 4 points
      face3{ { pts[7], pts[6], pts[2], pts[3] } }, // front
      face3{ { pts[4], pts[0], pts[1], pts[5] } }, // back
      face3{ { pts[7], pts[3], pts[0], pts[4] } }, // left
      face3{ { pts[6], pts[5], pts[1], pts[2] } }, // right
      face3{ { pts[3], pts[2], pts[1], pts[0] } }, // top
      face3{ { pts[7], pts[4], pts[5], pts[6] } }, // bottom
    };

    for ( face3 const& f : polyhedron ) {
      vec3 p2f = f.v[0] - cp;  // an arbitrary point on face
      double d = p2f.dot( f.normal() );
      d /= p2f.norm();

      double bound = -1e-15; // to exclude boundaries
      if ( d < bound )
        return false;
    }

    return true;
  }

	size_t generateCellCentroid(ColumnCell& cell) {    
		Point<3> pts[8], cp(0, 0, 0), tcp(0, 0, 0), bcp(0, 0, 0);
    
		pts[0] = p0->GetPoint(cell.z[0][0]);
    tcp += pts[0];
    pts[1] = p1->GetPoint( cell.z[1][0] );
    tcp += pts[1];
    pts[2] = p2->GetPoint( cell.z[2][0] );
    tcp += pts[2];
    pts[3] = p3->GetPoint( cell.z[3][0] );
    tcp += pts[3];
    pts[4] = p0->GetPoint(cell.z[0][1]); 
    bcp += pts[4];    
    pts[5] = p1->GetPoint(cell.z[1][1]);
    bcp += pts[5];    
    pts[6] = p2->GetPoint(cell.z[2][1]);
    bcp += pts[6];    
    pts[7] = p3->GetPoint(cell.z[3][1]);
    bcp += pts[7];    

    tcp = tcp * 0.25;
    bcp = bcp * 0.25;
    cp = (tcp + bcp) * 0.5; // in normal cases
    Point<3> rcp0 = (pts[0] + pts[6])*0.5; // in exceptial cases such as skewed cells
    Point<3> rcp1 = (pts[2] + pts[4])*0.5;
    Point<3> rcp2 = (pts[1] + pts[7])*0.5;
    Point<3> rcp3 = (pts[3] + pts[5])*0.5;

    std::vector<vec3> cp_vecs;
    cp_vecs.push_back( vec3( cp[0], cp[1], cp[2] ));
    cp_vecs.push_back( vec3( rcp0[0], rcp0[1], rcp0[2] ) );
    cp_vecs.push_back( vec3( rcp1[0], rcp1[1], rcp1[2] ) );
    cp_vecs.push_back( vec3( rcp2[0], rcp2[1], rcp2[2] ) );
    cp_vecs.push_back( vec3( rcp3[0], rcp3[1], rcp3[2] ) );

    std::vector<vec3> vecs; vecs.reserve( 8 );
    for(auto p: pts) vecs.push_back( vec3( p[0], p[1], p[2] ) );
    
    // find a valid centroid which is inside of its cell
    bool found (false);
    for ( auto v : cp_vecs ) {      
      if ( is_inside_polyheadron( v, vecs ) ){
        cp = Point<3>( v.x, v.y, v.z ); 
        found = true;
        break; 
      }
    }
    
    // the centroid is not inside of the cell, return zero, otherwise return its new node index
    if ( found ) { 
      size_t idx = ordinaryNodes.size() + extraNodes.size();
      extraNodes.push_back( cp );
      return idx;
    }
    else
      return found;
	}

	size_t  vertexIDs[8];
	IsoparametricLinearHexahedron hexa;
	IsoparametricLinearPyramid pyra;
	IsoparametricLinearTetrahedron tetra;
	IsoparametricLinearPrism prism;
	IsoparametricLinearLineElement line;

	bool ConstructLine(ColumnCell& cell) {
		vertexIDs[0] = 0;
		auto point_p1 = getNodeCoord(cell, vertexIDs[0]);
		grid.ConvertFromReservoirToCSMPcoordinateSystem(point_p1);

		vertexIDs[1] = 4;
    auto point_p2 = getNodeCoord(cell, vertexIDs[1]);
		grid.ConvertFromReservoirToCSMPcoordinateSystem(point_p2);

		line.XY.Resize(2U, 3U);
		line.XY.AssignRow(0U, point_p1);
		line.XY.AssignRow(0U, point_p2);
				
		return true;
	}

	size_t EmitLine(ColumnCell& cell) {
		auto ids = getGlobalIDList(cell, 2, vertexIDs);
		size_t elid = elementID++;
		plist.emplace(elid, ids);
		fem_types.push_back(ISOPARAMETRIC_LINEAR_BAR);
		return elid;
	}

	bool ConstructHexahedron(ColumnCell& cell) {
		//std::cerr << "hexa: \n";				
		for ( uint32_t i = 0; i < 8; ++i) {
			vertexIDs[i] = i;
			auto p = getNodeCoord(cell, vertexIDs[i]);
			grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
			hexa.XYZ(i, 0, p[0]);
			hexa.XYZ(i, 1, p[1]);
			hexa.XYZ(i, 2, p[2]);
			//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
		}

		uint32_t iNrIps = hexa.IntegrationPoints();
		for (uint32_t iIp = 0u; iIp < iNrIps; ++iIp) {
			hexa.JacobianAtIntegrationPoint(iIp);
			double jacdet = hexa.JacobianDeterminant();
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
		for (uint32_t i = 0; i < 5; ++i) {
			auto p = getGlobalNodeCoord(vertexIDs[i]);
			//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
			grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
			pyra.XYZ(i, 0, p[0]);
			pyra.XYZ(i, 1, p[1]);
			pyra.XYZ(i, 2, p[2]);
			//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
		}

    uint32_t iNrIps = pyra.IntegrationPoints();
		for (uint32_t iIp = 0; iIp < iNrIps; ++iIp) {
			pyra.JacobianAtIntegrationPoint(iIp);
			double jacdet = pyra.JacobianDeterminant();
			//std::cerr << "jacdet = " << jacdet << '\n';
			if (jacdet <= 0) {
				return false;
			}
		}
		return true;
	}

	size_t EmitPyramid(ColumnCell& ) {
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
		for (uint32_t i = 0; i < 4; ++i) {
			auto p = getGlobalNodeCoord(vertexIDs[i]);
			//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
			grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
			tetra.XYZ(i, 0, p[0]);
			tetra.XYZ(i, 1, p[1]);
			tetra.XYZ(i, 2, p[2]);
			//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
		}

    uint32_t iNrIps = tetra.IntegrationPoints();
		for (uint32_t iIp = 0; iIp < iNrIps; ++iIp) {
			tetra.JacobianAtIntegrationPoint(iIp);
			double jacdet = tetra.JacobianDeterminant();
			//std::cerr << "jacdet = " << jacdet << '\n';
			if (jacdet <= 0) {
				return false;
			}
		}
		return true;
	}

	size_t EmitTetrahedron(ColumnCell& ) {
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
		for ( uint32_t i = 0; i < 6; ++i) {
			auto p = getGlobalNodeCoord(vertexIDs[i]);
			grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
			prism.XYZ(i, 0, p[0]);
			prism.XYZ(i, 1, p[1]);
			prism.XYZ(i, 2, p[2]);
			//std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
		}

		auto iNrIps = prism.IntegrationPoints();
		for ( uint32_t iIp = 0; iIp < iNrIps; ++iIp) {
			prism.JacobianAtIntegrationPoint(iIp);
			double jacdet = prism.JacobianDeterminant();
			//std::cerr << "jacdet = " << jacdet << '\n';
			if (jacdet <= 0) {
				return false;
			}
		}
		return true;
	}

	size_t EmitPrism(ColumnCell& ) {
		std::vector<size_t> ids(&vertexIDs[0], &vertexIDs[6]);
		size_t elid = elementID++;
		plist.emplace(elid, ids);
		fem_types.push_back(ISOPARAMETRIC_LINEAR_PRISM);
		return elid;
	}


	// return a list of globalNodeID from local vertex index
	  std::vector<size_t> getGlobalIDList(ColumnCell& cell, size_t list_size, const size_t* vertexID ) {
        std::vector<size_t> nodeIDs;
        nodeIDs.reserve(list_size);
        for (size_t i = 0; i < list_size; ++i) {
          nodeIDs.push_back(getNodeID(cell, vertexID[i]));
        }
        return nodeIDs;
	}

	/// degenerates to one prism

	std::vector<size_t> degenerateToOnePrismAtEdge23(ColumnCell& cell) {    // 034 125
    const size_t  vertices[6] = { 0, 3, 4, 1, 2, 5 };
		return getGlobalIDList(cell, 6, vertices);
	}

	std::vector<size_t> degenerateToOnePrismAtEdge30(ColumnCell& cell) { //051 362
		const size_t  vertices[6] = { 0, 5, 1, 3, 6, 2 };
		return getGlobalIDList(cell, 6, vertices );
	}

	std::vector<std::vector<size_t>> degenerateToTwoTetrahedrasAtEdge02(ColumnCell& cell) {    // 1010
		std::vector<std::vector<size_t>> nodeLists;          // tetra 0125 0237
		nodeLists.reserve(2);

		static const size_t  vertexIDs1[4] = { 0, 1, 2, 5 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

		static const size_t  vertexIDs2[4] = { 0, 2, 3, 7 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

		return nodeLists;
	}

	std::vector<std::vector<size_t>> degenerateToTwoTetrahedraAtEdge13(ColumnCell& cell) {    // 0101
		std::vector<std::vector<size_t>> nodeLists;          // tetra 1304 1326
		nodeLists.reserve(2);

		static const size_t  vertexIDs1[4] = { 1, 3, 0, 4 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

		static const size_t  vertexIDs2[4] = { 1, 3, 2, 6 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

		return nodeLists;
	}

	std::vector<std::vector<size_t>> splitToFiveTetrahedrasAtTwoEdges0257(ColumnCell& cell) {
		std::vector<std::vector<size_t>> nodeLists;                // 0457 0125 0237 0257 2756
		nodeLists.reserve(5);

		static const size_t  vertexIDs1[4] = { 0, 4, 5, 7 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

		static const size_t  vertexIDs2[4] = { 0, 1, 2, 5 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

		static const size_t  vertexIDs3[4] = { 0, 2, 3, 7 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs3));

		static const size_t  vertexIDs4[4] = { 0, 2, 5, 7 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs4));

		static const size_t  vertexIDs5[4] = { 2, 7, 5, 6 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs5));

		return nodeLists;
	}

	std::vector<std::vector<size_t>> splitToFiveTetrahedrasAtTwoEdges1347(ColumnCell& cell) {
		std::vector<std::vector<size_t>> nodeLists;                // 0134 1456 1236 1346 3467
		nodeLists.reserve(5);

		static const size_t  vertexIDs1[4] = { 0, 1, 3, 4 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs1));

		static const size_t  vertexIDs2[4] = { 1, 4, 5, 6 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs2));

		static const size_t  vertexIDs3[4] = { 1, 2, 3, 6 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs3));

		static const size_t  vertexIDs4[4] = { 1, 3, 4, 6 };
		nodeLists.push_back(getGlobalIDList(cell, 4, vertexIDs4));

		static const size_t  vertexIDs5[4] = { 3, 4, 6, 7 };
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
      if ( cell.z[0][0] < cell.z[0][1] ) {
        return FACE_TYPE::SPLIT_13_OR_57;
        //return FACE_TYPE::SPLIT_02_OR_46;
      }
      else return FACE_TYPE::FULL_QUAD;
    }
    case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_1: {
      if ( cell.z[1][0] < cell.z[1][1] ) {
        return FACE_TYPE::SPLIT_02_OR_46;
        //return FACE_TYPE::SPLIT_13_OR_57;
      }
      else return FACE_TYPE::FULL_QUAD;
    }
    case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_2: {
      if ( cell.z[2][0] < cell.z[2][1] ) {
        return FACE_TYPE::SPLIT_13_OR_57;
        //return FACE_TYPE::SPLIT_02_OR_46;
      }
      else return FACE_TYPE::FULL_QUAD;
    }
    case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3: {
      if ( cell.z[3][0] < cell.z[3][1] ) {
        return FACE_TYPE::SPLIT_02_OR_46;
        //return FACE_TYPE::SPLIT_13_OR_57;
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

} // end namespace csmp

#endif
