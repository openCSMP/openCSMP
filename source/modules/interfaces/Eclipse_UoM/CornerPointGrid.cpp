#include "CornerPointGrid_UoM.h"
#include "PropertyData.h"
#include "CSMP_highLevelUtilities.h"
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
          p1 = &grid(i + 0, j + 1);
          p2 = &grid(i + 1, j + 1);
          p3 = &grid(i + 1, j + 0);
      }

      size_t generateCellCentroid(ColumnCell& cell) {
        Point<3> p(0,0,0);
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
        for (size_t i = 0; i < 8; ++i) {
          vertexIDs[i] = i;
          auto p = getNodeCoord(cell, vertexIDs[i]);
          grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
          hexa.XYZ(i, 0, p[0]);
          hexa.XYZ(i, 1, p[1]);
          hexa.XYZ(i, 2, p[2]);
          // std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }

        size_t iNrIps = hexa.IntegrationPoints();
        for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
          hexa.JacobianAtIntegrationPoint(iIp);
          double64 jacdet = hexa.JacobianDeterminant();
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

        for (size_t i = 0; i < 5; ++i) {
          auto p = getGlobalNodeCoord(vertexIDs[i]);
          grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
          pyra.XYZ(i, 0, p[0]);
          pyra.XYZ(i, 1, p[1]);
          pyra.XYZ(i, 2, p[2]);
          // std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }
        
        size_t iNrIps = pyra.IntegrationPoints();
        for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
          pyra.JacobianAtIntegrationPoint(iIp);
          double64 jacdet = pyra.JacobianDeterminant();
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

        for (size_t i = 0; i < 4; ++i) {
          auto p = getGlobalNodeCoord(vertexIDs[i]);
          grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
          tetra.XYZ(i, 0, p[0]);
          tetra.XYZ(i, 1, p[1]);
          tetra.XYZ(i, 2, p[2]);
          // std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }
        
        size_t iNrIps = tetra.IntegrationPoints();
        for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
          tetra.JacobianAtIntegrationPoint(iIp);
          double64 jacdet = tetra.JacobianDeterminant();
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

        for (size_t i = 0; i < 6; ++i) {
          auto p = getGlobalNodeCoord(vertexIDs[i]);
          grid.ConvertFromReservoirToCSMPcoordinateSystem(p);
          prism.XYZ(i, 0, p[0]);
          prism.XYZ(i, 1, p[1]);
          prism.XYZ(i, 2, p[2]);
          // std::cerr << "p" << i << " = " << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }
        
        size_t iNrIps = prism.IntegrationPoints();
        for (size_t iIp = 0; iIp < iNrIps; ++iIp) {
          prism.JacobianAtIntegrationPoint(iIp);
          double64 jacdet = prism.JacobianDeterminant();
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
            }
            else return FACE_TYPE::FULL_QUAD; 
          }
          case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_1: {
            if (cell.z[1][0] < cell.z[1][1]) {
              return FACE_TYPE::SPLIT_13_OR_57;
            }
            else return FACE_TYPE::FULL_QUAD; 
          }
          case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_2: {
            if (cell.z[2][0] < cell.z[2][1]) { 
              return FACE_TYPE::SPLIT_02_OR_46;
            }
            else return FACE_TYPE::FULL_QUAD; 
          }
          case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3: {
            if (cell.z[3][0] < cell.z[3][1]) {
              return FACE_TYPE::SPLIT_13_OR_57;
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

    /**
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

            // Sort and unique
            std::sort(zcoord.begin(), zcoord.end());
            zcoord.erase(std::unique(zcoord.begin(), zcoord.end()), zcoord.end());
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
            Pillar& p1 = (*this)(i + 0, j + 1);
            Pillar& p2 = (*this)(i + 1, j + 1);
            Pillar& p3 = (*this)(i + 1, j + 0);
            for (size_t k = 0; k < NZ_; ++k)
            {
              if (!CellActivity(i, j, k)) {
                continue;
              }

              ColumnCell cell;
              cell.k = k;
              cell.z[0][0] = p0.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 0]); // t_nw
              cell.z[0][1] = p0.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 4]); // b_nw
              cell.z[1][0] = p1.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 2]); // t_sw
              cell.z[1][1] = p1.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 6]); // b_sw
              cell.z[2][0] = p2.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 3]); // t_se
              cell.z[2][1] = p2.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 7]); // b_se
              cell.z[3][0] = p3.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 1]); // t_ne
              cell.z[3][1] = p3.FindPoint(zcorn[(i + j * NX_ + k * NXxNY) * 8 + 5]); // b_ne

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
      vset.HybridElementTypeMesh( true );

      //map<size_t, vector<long64>>  pfverts;

      // 1. Classify the cells

      size_t skewCells = 0;
      for (auto& index_column : columns_) {
        auto& column = index_column.second;
        size_t i = index_column.first.first;
        size_t j = index_column.first.second;

        const size_t iNrCells = column.cells_.size();
        for (size_t iCell = 0; iCell < iNrCells; ++iCell) {
          auto& cell = column.cells_[iCell];

          Pillar& p0 = (*this)(i + 0, j + 0);
          Pillar& p1 = (*this)(i + 0, j + 1);
          Pillar& p2 = (*this)(i + 1, j + 1);
          Pillar& p3 = (*this)(i + 1, j + 0);

          double64 z[4][2];
          z[0][0] = p0.GetZCoord(cell.z[0][0]);
          z[0][1] = p0.GetZCoord(cell.z[0][1]);
          z[1][0] = p1.GetZCoord(cell.z[1][0]);
          z[1][1] = p1.GetZCoord(cell.z[1][1]);
          z[2][0] = p2.GetZCoord(cell.z[2][0]);
          z[2][1] = p2.GetZCoord(cell.z[2][1]);
          z[3][0] = p3.GetZCoord(cell.z[3][0]);
          z[3][1] = p3.GetZCoord(cell.z[3][1]);

          // Classify the degeneracy
          uint8_t classification = 0;
          for (size_t v = 0; v < 4; ++v) {
            if (cell.z[v][0] == cell.z[v][1]) {
              classification |= (1 << (3-v));
            }
          }
          cell.classification = static_cast<ECLIPSE_CELL_CLASSIFICATION>(classification);
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

      for (auto column = columns_.begin(); column != columns_.end(); column++)  {
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
            cellAbove = &Col.cells_[k-1];
          }

          if (k < Col.cells_.size() - 1) {
            cellBeneath = &Col.cells_[k+1];
          }

          assert(generator.elementID == generator.plist.size());
          assert(generator.elementID == generator.fem_types.size());

          if (i == 34 && j == 31 && 99 <= k && k <= 103) {
              std::cerr << "Broken case\n";
          }

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

                    default: ;
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

                        default: ;
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

                        default: ;
                      }
                      break;
                    }
                    default: ;
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
#if 0
              // XXX BROKEN
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
#endif

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
              if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 0))) {
                addElementToMap(i, j, k, generator.EmitPyramid(cell));
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
              if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 1))) {
                addElementToMap(i, j, k, generator.EmitPyramid(cell));
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
              if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 2))) {
                addElementToMap(i, j, k, generator.EmitPyramid(cell));
              }
              else {
                ++badPyramids;
              }
              break;
            }

            case ECLIPSE_CELL_CLASSIFICATION::ECLIPSE_CELL_PYRAMID_3: {                    // 1110
              if (generator.ConstructPyramidOnFace(cell, 7, 6, 5, 4, generator.getNodeID(cell, 3))) {
                addElementToMap(i, j, k, generator.EmitPyramid(cell));
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
        for (auto& p : generator.extraNodes) {
          x.push_back(p[0]);
          y.push_back(p[1]);
          z.push_back(p[2]);
        }
        vset.AddXYZ(x, y, z);
      }

      // 4. Set up the rest of the vset
      vset.ResizePfverts(generator.fem_types.size());
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
      this->elementMap.emplace(ijk(i,j,k), elementID);
    };

    // CreateModel (MASTER - MONSTER METHOD)



#if 0
    /**
     Identifies the faces of the cells that make up the boundary.
     The cells are then grouped into boundary regions that are stored in the ModelTopology class.
     */
    void CornerPointGrid::EstablishBoundaries(csmp::VSet<3U>& vset,
                                              csmp::ModelTopology& model_topology,
                                              const std::vector<CornerPointCell>& poly)
    {
      if (!exclude_inactive_cells_) {
        throw csmp::Exception(ERROR, "CornerPointGrid::EstablishBoundaryRegions",
                              "Cannot create boundary regions if we are excluding inactive cells");
      }

      for (size_t j = 0; j < NY_; ++j) {
        for (size_t i = 0; i < NX_; ++i) {
          size_t kmin = NZ_;
          size_t kmax = 0;
          for (size_t k = 0; k < NZ_; ++k) {
            size_t hexa_cell_id = CellIndex(i, j, k);

            if (cell_activity_[hexa_cell_id] == 1) {
              kmin = std::min(kmin, k);
              kmax = std::max(kmax, k);
            }
          }

          for (size_t k = 0; k < NZ_; ++k) {
            /// hexa cell global id
            size_t hexa_cell_id = CellIndex(i, j, k);

            if (cell_activity_[hexa_cell_id] == 1) {
              const uint8_t TOP_FLAG = 1;
              const uint8_t BOTTOM_FLAG = 2;
              const uint8_t IRREGULAR_FLAG = 4;

              uint8_t node_flags[8];

              std::memset(node_flags, 0, sizeof(node_flags));

              if (i + 1 == NX_ || !cell_activity_[CellIndex(i + 1, j, k)]) {
                node_flags[1] |= IRREGULAR_FLAG;
                node_flags[3] |= IRREGULAR_FLAG;
                node_flags[5] |= IRREGULAR_FLAG;
                node_flags[7] |= IRREGULAR_FLAG;
              }
              if (i == 0 || !cell_activity_[CellIndex(i - 1, j, k)]) {
                node_flags[0] |= IRREGULAR_FLAG;
                node_flags[2] |= IRREGULAR_FLAG;
                node_flags[4] |= IRREGULAR_FLAG;
                node_flags[6] |= IRREGULAR_FLAG;
              }
              if (j + 1 == NY_ || !cell_activity_[CellIndex(i, j + 1, k)]) {
                node_flags[2] |= IRREGULAR_FLAG;
                node_flags[3] |= IRREGULAR_FLAG;
                node_flags[6] |= IRREGULAR_FLAG;
                node_flags[7] |= IRREGULAR_FLAG;
              }
              if (j == 0 || !cell_activity_[CellIndex(i, j - 1, k)]) {
                node_flags[0] |= IRREGULAR_FLAG;
                node_flags[2] |= IRREGULAR_FLAG;
                node_flags[4] |= IRREGULAR_FLAG;
                node_flags[6] |= IRREGULAR_FLAG;
              }
              if (k == kmax) {
                node_flags[0] |= BOTTOM_FLAG;
                node_flags[1] |= BOTTOM_FLAG;
                node_flags[2] |= BOTTOM_FLAG;
                node_flags[3] |= BOTTOM_FLAG;
              }
              if (k == kmin) {
                node_flags[4] |= TOP_FLAG;
                node_flags[5] |= TOP_FLAG;
                node_flags[6] |= TOP_FLAG;
                node_flags[7] |= TOP_FLAG;
              }

              auto& polygon = poly[hexa_cell_id];
              for (size_t n = 0; n < 8; ++n) {
                auto nid = polygon.GetPillarNodeGlobalIdOriginalOrder(n);
                switch (node_flags[n]) {
                  case 0:
                  {
                    break;
                  }

                  case 1: // TOP_FLAG
                  case 5: // IRREGULAR_FLAG | TOP_FLAG
                  {
                    vset.AddBFlag(nid, TOP_OUTSIDE);
                    break;
                  }

                  case 2: // BOTTOM_FLAG
                  case 6: // BOTTOM_FLAG | IRREGULAR_FLAG
                  {
                    vset.AddBFlag(nid, BOTTOM_OUTSIDE);
                    break;
                  }

                  case 4: // IRREGULAR_FLAG
                  {
                    vset.AddBFlag(nid, IRREGULAR_OUTSIDE);
                    break;
                  }

                  case 3: // TOP_FLAG | BOTTOM_FLAG
                  case 7: // TOP_FLAG | BOTTOM_FLAG | IRREGULAR_FLAG
                  {
                    std::cerr << "CornerPointGrid::EstablishBoundaries: "
                    << "Possibly erroneous node " << nid << " at "
                    << vset.Px(nid) << ',' << vset.Py(nid) << ',' << vset.Pz(nid) << '\n';
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }


    /**
     Identifies the faces of the cells that make up a fault.
     The cells are then grouped into separate fault regions that are stored in the ModelTopology class.

     The faces=side walls of the fault are added to the supplied VSet.

     this breaks the record in terms of Method size and nested loops ! (@@)

     TODO: this function never modifies the VSet although we tried several cases with faults present; why?
     */
    void CornerPointGrid::EstablishFaultRegions(csmp::VSet<3U>& vset,
                                                csmp::ModelTopology& model_topology,
                                                std::set<std::string>& faults,
                                                const std::vector<CornerPointCell>& poly)
    {
      /// cell data
      const size_t quad_fem_nodes(4U);
      size_t hexa_face_id(0);
      size_t hexa_cell_id(0);
      size_t hexa_cell_neighbor_id(0);
      size_t hexa_face_neighbor_id(0);
      std::pair<size_t, size_t> face, nface;

      /// opposite faces
      std::vector<size_t> opposite_face(6U, 0);
      opposite_face[CORNER_POINT_CELL_FACE_Xminus] = CORNER_POINT_CELL_FACE_Xplus;
      opposite_face[CORNER_POINT_CELL_FACE_Xplus] = CORNER_POINT_CELL_FACE_Xminus;
      opposite_face[CORNER_POINT_CELL_FACE_Yminus] = CORNER_POINT_CELL_FACE_Yplus;
      opposite_face[CORNER_POINT_CELL_FACE_Yplus] = CORNER_POINT_CELL_FACE_Yminus;
      opposite_face[CORNER_POINT_CELL_FACE_Zminus] = CORNER_POINT_CELL_FACE_Zplus;
      opposite_face[CORNER_POINT_CELL_FACE_Zplus] = CORNER_POINT_CELL_FACE_Zminus;

      /// neighbors from each side
      std::vector<long> opposite_cell_offset(6U, 0);
      opposite_cell_offset[CORNER_POINT_CELL_FACE_Xminus] = -1;
      opposite_cell_offset[CORNER_POINT_CELL_FACE_Xplus] = +1;
      opposite_cell_offset[CORNER_POINT_CELL_FACE_Yminus] = -NX_;
      opposite_cell_offset[CORNER_POINT_CELL_FACE_Yplus] = +NX_;
      opposite_cell_offset[CORNER_POINT_CELL_FACE_Zminus] = -NX_x_NY_;
      opposite_cell_offset[CORNER_POINT_CELL_FACE_Zplus] = +NX_x_NY_;

      std::vector<std::vector<size_t> > face_nodes(6U, std::vector<size_t>(5, 0));
      face_nodes[CORNER_POINT_CELL_FACE_Xminus][0] = 0;
      face_nodes[CORNER_POINT_CELL_FACE_Xminus][1] = 2;
      face_nodes[CORNER_POINT_CELL_FACE_Xminus][2] = 6;
      face_nodes[CORNER_POINT_CELL_FACE_Xminus][3] = 4;
      face_nodes[CORNER_POINT_CELL_FACE_Xminus][4] = 0;

      face_nodes[CORNER_POINT_CELL_FACE_Xplus][0] = 1;
      face_nodes[CORNER_POINT_CELL_FACE_Xplus][1] = 3;
      face_nodes[CORNER_POINT_CELL_FACE_Xplus][2] = 7;
      face_nodes[CORNER_POINT_CELL_FACE_Xplus][3] = 5;
      face_nodes[CORNER_POINT_CELL_FACE_Xplus][4] = 1;

      face_nodes[CORNER_POINT_CELL_FACE_Yminus][0] = 0;
      face_nodes[CORNER_POINT_CELL_FACE_Yminus][1] = 1;
      face_nodes[CORNER_POINT_CELL_FACE_Yminus][2] = 5;
      face_nodes[CORNER_POINT_CELL_FACE_Yminus][3] = 4;
      face_nodes[CORNER_POINT_CELL_FACE_Yminus][4] = 0;

      face_nodes[CORNER_POINT_CELL_FACE_Yplus][0] = 2;
      face_nodes[CORNER_POINT_CELL_FACE_Yplus][1] = 3;
      face_nodes[CORNER_POINT_CELL_FACE_Yplus][2] = 7;
      face_nodes[CORNER_POINT_CELL_FACE_Yplus][3] = 6;
      face_nodes[CORNER_POINT_CELL_FACE_Yplus][4] = 2;

      face_nodes[CORNER_POINT_CELL_FACE_Zminus][0] = 0;
      face_nodes[CORNER_POINT_CELL_FACE_Zminus][1] = 1;
      face_nodes[CORNER_POINT_CELL_FACE_Zminus][2] = 3;
      face_nodes[CORNER_POINT_CELL_FACE_Zminus][3] = 2;
      face_nodes[CORNER_POINT_CELL_FACE_Zminus][4] = 0;

      face_nodes[CORNER_POINT_CELL_FACE_Zplus][0] = 4;
      face_nodes[CORNER_POINT_CELL_FACE_Zplus][1] = 5;
      face_nodes[CORNER_POINT_CELL_FACE_Zplus][2] = 7;
      face_nodes[CORNER_POINT_CELL_FACE_Zplus][3] = 6;
      face_nodes[CORNER_POINT_CELL_FACE_Zplus][4] = 4;

      std::string fault_name, fault_name_a, fault_name_b;
      std::pair<std::map<std::string, std::vector<std::pair<size_t, size_t> > >::iterator, bool> faults_data_it;
      faults.clear();

      // loops over the faults assuming that their cells are adjacent to discontinuities = cell interfaces that are not shared between
      // adjacent cells.
      // To deal with these lower-dimensional surfaces, they are grouped in terms of sides.
      // These sides are logged in 2 new separate regions for furthert processing.
      if (!faults_data_.empty())
      {
        // DEBUGGING
        //for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
        //     it = faults_data_.begin(); it != faults_data_.end(); ++it ) {
        //    std::cerr <<"\n\n"<< (*it).first <<" (cell-ID/face-ID pairs): ";
        //    for ( std::vector<std::pair<size_t,size_t> >::const_iterator
        //          fit=(*it).second.begin(); fit!=(*it).second.end(); fit++ )
        //      std::cerr << (*fit).first <<","<< (*fit).second <<" ";
        // }
        //std::cerr <<"\n";

        // map of fault-name keys to vectors of cell-ID/face-ID pairs
        std::map<std::string, std::vector<std::pair<size_t, size_t> > > temp_faults_data;
        for (std::map<std::string, std::vector<std::pair<size_t, size_t> > >::const_iterator
             it = faults_data_.begin(); it != faults_data_.end(); ++it)
        {
          fault_name = (*it).first;

          /// in case if side is already defined
          if (fault_name.find("side") != std::string::npos)
          {
            temp_faults_data.insert(std::make_pair(fault_name, (*it).second));
          }
          else
          {
            /// add new fault region ( side A )
            fault_name_a = fault_name;
            fault_name_a += "_sideA";
            faults_data_it = temp_faults_data.insert(std::make_pair(fault_name_a, std::vector<std::pair<size_t, size_t> >()));
            std::vector<std::pair<size_t, size_t> >& fault_elmts_a((*faults_data_it.first).second);
            /// add new fault region ( side B )
            fault_name_b = fault_name;
            fault_name_b += "_sideB";
            faults_data_it = temp_faults_data.insert(std::make_pair(fault_name_b, std::vector<std::pair<size_t, size_t> >()));
            std::vector<std::pair<size_t, size_t> >& fault_elmts_b((*faults_data_it.first).second);

            // filling the cell-ID/face-ID vectors for the identified fault faces
            for (std::vector<std::pair<size_t, size_t> >::const_iterator
                 vit = (*it).second.begin(); vit != (*it).second.end(); ++vit)
            {
              /// assign data from provided face
              hexa_cell_id = (*vit).first;
              if (cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_)
              {
                hexa_face_id = (*vit).second;
                if (poly[hexa_cell_id].GetNumPolygonFaceNodes(hexa_face_id) > 2U)
                  fault_elmts_a.push_back(*vit);
              }

              /// assign data from opposite face
              hexa_cell_neighbor_id = hexa_cell_id + opposite_cell_offset[hexa_face_id];
              if (cell_activity_[hexa_cell_neighbor_id] == 1 || !exclude_inactive_cells_)
              {
                hexa_face_neighbor_id = opposite_face[hexa_face_id];
                if (poly[hexa_cell_neighbor_id].GetNumPolygonFaceNodes(hexa_face_neighbor_id) > 2U)
                  fault_elmts_b.push_back(std::make_pair(hexa_cell_neighbor_id, hexa_face_neighbor_id));
              }
            }
            if (fault_elmts_a.empty())
              temp_faults_data.erase(fault_name_a);
            if (fault_elmts_b.empty())
              temp_faults_data.erase(fault_name_b);
          }
        }
        // the old faults dataset is replaced by the new faults_data
        faults_data_.clear();
        faults_data_ = temp_faults_data;
        temp_faults_data.clear();

        // DEBUGGING - appears to recreate faults data exactly = duplicating it
        //std::cerr <<"\n\nrecreated fault data:\n";
        //for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
        //     it = faults_data_.begin(); it != faults_data_.end(); ++it ) {
        //    std::cerr <<"\n\n"<< (*it).first <<" (cell-ID/face-ID pairs): ";
        //    for ( std::vector<std::pair<size_t,size_t> >::const_iterator
        //          fit=(*it).second.begin(); fit!=(*it).second.end(); fit++ )
        //      std::cerr << (*fit).first <<","<< (*fit).second <<" ";
        // }
        //std::cerr <<"\n";

      }
      else // if no fault data are provided this information is created automatically (how?)
      {
        std::set<size_t> face_nidsA;
        std::set<size_t> face_nidsB;
        std::set<size_t> edge;
        std::set<std::pair<size_t, size_t> > fault_faces;

        /// face for each direction
        std::vector<std::vector<size_t> > direction_faces(3, std::vector<size_t>(2, 0));
        direction_faces[0][0] = CORNER_POINT_CELL_FACE_Xplus;
        direction_faces[0][1] = CORNER_POINT_CELL_FACE_Xminus;
        direction_faces[1][0] = CORNER_POINT_CELL_FACE_Yplus;
        direction_faces[1][1] = CORNER_POINT_CELL_FACE_Yminus;
        direction_faces[2][0] = CORNER_POINT_CELL_FACE_Zplus;
        direction_faces[2][1] = CORNER_POINT_CELL_FACE_Zminus;

        ///                  edge nodes         cell id , face id
        std::multimap<std::set<size_t>, std::pair<size_t, size_t> > nc;

        /// find fault faces in each direction
        for (size_t d = 0; d < 2U; ++d)
        {
          face.second = direction_faces[d][0];
          nface.second = direction_faces[d][1];
          for (size_t k = 0; k < NZ_; ++k)
            for (size_t j = 0; j < NY_; ++j)
              for (size_t i = 0; i < NX_; ++i)
              {
                hexa_cell_id = CellIndex(i, j, k);
                hexa_cell_neighbor_id = hexa_cell_id;
                hexa_cell_neighbor_id += opposite_cell_offset[direction_faces[d][0]];
                face.first = hexa_cell_id;
                nface.first = hexa_cell_neighbor_id;

                /// only internal faces
                if ((i == NX_ - 1 && d == 0) || (j == NY_ - 1 && d == 1) || (k == NZ_ - 1 && d == 2))
                  continue;

                /// only faces between active cells
                if ((exclude_inactive_cells_ && cell_activity_[hexa_cell_id] == 1 && cell_activity_[hexa_cell_neighbor_id] == 1)
                    || !exclude_inactive_cells_)
                {
                  /// comparing potentially neighboring faces
                  face_nidsA.clear();
                  for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
                    face_nidsA.insert(poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(face_nodes[face.second][nid]));

                  face_nidsB.clear();
                  for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
                    face_nidsB.insert(poly[hexa_cell_neighbor_id].GetPillarNodeGlobalIdOriginalOrder(face_nodes[nface.second][nid]));

                  if (face_nidsA != face_nidsB)
                  {
                    if (face_nidsA.size() > 2U)
                    {
                      fault_faces.insert(face);
                      for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
                      {
                        edge.clear();
                        edge.insert(poly[face.first].GetPillarNodeGlobalIdOriginalOrder(face_nodes[face.second][nid]));
                        edge.insert(poly[face.first].GetPillarNodeGlobalIdOriginalOrder(face_nodes[face.second][nid + 1]));
                        if (edge.size() == 2)
                          nc.insert(std::make_pair(edge, face));
                      }
                    }
                    if (face_nidsB.size() > 2U)
                    {
                      fault_faces.insert(nface);
                      for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
                      {
                        edge.clear();
                        edge.insert(poly[nface.first].GetPillarNodeGlobalIdOriginalOrder(face_nodes[nface.second][nid]));
                        edge.insert(poly[nface.first].GetPillarNodeGlobalIdOriginalOrder(face_nodes[nface.second][nid + 1]));
                        if (edge.size() == 2)
                          nc.insert(std::make_pair(edge, nface));
                      }
                    }
                  }
                }
              }
        }

        // if no fault data were provided, but they were recreated, potential subsets are identified as above for the faults data (repetition of code)
        if (!fault_faces.empty())
        {
          /// identify fault subsets
          std::string faultname_prefix = "FAULT";
          size_t      faultname_index = 1;
          std::pair<std::multimap<std::set<size_t>, std::pair<size_t, size_t> >::iterator, std::multimap<std::set<size_t>, std::pair<size_t, size_t> >::iterator> neighbor_faces;
          std::set<std::pair<size_t, size_t> >::iterator it;
          std::set<std::pair<size_t, size_t> > neighbors;
          bool found_face(false);
          while (!fault_faces.empty())
          {
            std::ostringstream oss;
            oss << faultname_prefix << faultname_index++;
            fault_name = oss.str();

            /// add new fault region ( side A )
            fault_name_a = fault_name;
            fault_name_a += "_sideA";
            faults_data_it = faults_data_.insert(std::make_pair(fault_name_a, std::vector<std::pair<size_t, size_t> >()));
            std::vector<std::pair<size_t, size_t> >& fault_elmts_a((*faults_data_it.first).second);

            /// add new fault region ( side B )
            fault_name_b = fault_name;
            fault_name_b += "_sideB";
            faults_data_it = faults_data_.insert(std::make_pair(fault_name_b, std::vector<std::pair<size_t, size_t> >()));
            std::vector<std::pair<size_t, size_t> >& fault_elmts_b((*faults_data_it.first).second);

            /// add first element
            it = fault_faces.begin();
            neighbors.insert(*it);
            fault_elmts_a.push_back(*it);
            fault_faces.erase(it);
            while (!neighbors.empty())
            {
              it = neighbors.begin();
              face = *it;
              neighbors.erase(it);

              for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
              {
                edge.clear();
                edge.insert(poly[face.first].GetPillarNodeGlobalIdOriginalOrder(face_nodes[face.second][nid]));
                edge.insert(poly[face.first].GetPillarNodeGlobalIdOriginalOrder(face_nodes[face.second][nid + 1]));
                if (edge.size() == 2)
                {
                  neighbor_faces = nc.equal_range(edge);
                  if (neighbor_faces.first != nc.end())
                  {
                    found_face = false;
                    /// first priority: faces of current cell
                    for (std::multimap<std::set<size_t>, std::pair<size_t, size_t> >::iterator
                         fit = neighbor_faces.first; fit != neighbor_faces.second; ++fit)
                    {
                      nface = (*fit).second;
                      if (nface.first == face.first && nface.second != face.second)
                      {
                        it = fault_faces.find(nface);
                        if (it != fault_faces.end())
                        {
                          fault_elmts_a.push_back(nface);
                          neighbors.insert(nface);
                          fault_faces.erase(it);
                          found_face = true;
                          break;
                        }
                      }
                    }
                    /// second priority: faces which does not belong to opposite cell
                    if (!found_face)
                    {
                      for (std::multimap<std::set<size_t>, std::pair<size_t, size_t> >::iterator
                           fit = neighbor_faces.first; fit != neighbor_faces.second; ++fit)
                      {
                        nface = (*fit).second;
                        if ((nface.first != face.first) && (nface.first != (face.first + opposite_cell_offset[face.second])))
                        {
                          it = fault_faces.find(nface);
                          if (it != fault_faces.end())
                          {
                            fault_elmts_a.push_back(nface);
                            neighbors.insert(nface);
                            fault_faces.erase(it);
                            break;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }

            /// add elements from opposite side
            for (std::vector<std::pair<size_t, size_t> >::const_iterator
                 fit = fault_elmts_a.begin(); fit != fault_elmts_a.end(); ++fit)
            {
              face = *fit;
              nface.first = face.first + opposite_cell_offset[face.second];
              nface.second = opposite_face[face.second];
              it = fault_faces.find(nface);
              if (it != fault_faces.end())
              {
                fault_elmts_b.push_back(*it);
                fault_faces.erase(it);
              }
            }
          }
        }

        // outputs "autogen" text file to record the fault information in Eclipse format
        if (!faults_data_.empty())
        {
          /// face map
          std::map<size_t, std::string> facemap;
          facemap.insert(std::make_pair(CORNER_POINT_CELL_FACE_Xminus, "X-"));
          facemap.insert(std::make_pair(CORNER_POINT_CELL_FACE_Yminus, "Y-"));
          facemap.insert(std::make_pair(CORNER_POINT_CELL_FACE_Zminus, "Z-"));
          facemap.insert(std::make_pair(CORNER_POINT_CELL_FACE_Xplus, "X+"));
          facemap.insert(std::make_pair(CORNER_POINT_CELL_FACE_Yplus, "Y+"));
          facemap.insert(std::make_pair(CORNER_POINT_CELL_FACE_Zplus, "Z+"));

          /// output
          std::string filename = model_topology.ModelName();
          filename += "-autogen-faults.txt";
          std::ofstream fout(filename, std::ofstream::out);
          fout << "-- Format      : ECLIPSE fault data (ASCII)" << std::endl;
          fout << "-- KEYWORD \"FAULTS\" HAS BEEN WRITTEN FROM CSMP++." << std::endl;
          fout << "FAULTS" << std::endl;
          fout << "-- Matrix Faults" << std::endl;
          for (std::map<std::string, std::vector<std::pair<size_t, size_t> > >::const_iterator
               fit = faults_data_.begin(); fit != faults_data_.end(); ++fit)
          {
            fout << "\n\n-- NAME\tIX1\tIX2\tIY1\tIY2\tIZ1\tIZ2\tFACE" << std::endl;
            fault_name = (*fit).first;
            for (std::vector<std::pair<size_t, size_t> >::const_iterator
                 ffit = (*fit).second.begin(); ffit != (*fit).second.end(); ++ffit)
            {
              hexa_cell_id = (*ffit).first;
              hexa_face_id = (*ffit).second;
              const size_t i = (hexa_cell_id%NX_);
              const size_t j = ((hexa_cell_id / NX_) % NY_);
              const size_t k = (hexa_cell_id / NX_x_NY_);
              fout << "'" << fault_name << "'"
              << "\t" << i << "\t" << i << "\t" << j << "\t" << j << "\t" << k << "\t" << k
              << "\t" << "'" << facemap[hexa_face_id] << "'" << std::endl;
            }
          }
          fout << '/' << std::endl;
          fout.close();
        }
      }

      /* process matching fault regions
       normally these shoud not match because the faults displace the grid.

       checks whether there is an offset in the faults.

       TODO: has this code been tested? - in all our examples it was not reached

       TODO: the algorithm can be extended to find out groups of mathing elements
       */
      if (!faults_data_.empty())
      {
        std::string fault_name_prefix;
#ifndef NEW
        std::map<std::string, std::vector<std::pair<size_t, size_t> > >::iterator fit;
        std::map<std::string, std::vector<std::pair<size_t, size_t> > > temp_faults_data;
        for (std::map<std::string, std::vector<std::pair<size_t, size_t> > >::iterator
#else /* NEW */
             std::map<std::string, vector<std::pair<size_t, size_t> > >::iterator fit;
             std::map<std::string, vector<std::pair<size_t, size_t> > > temp_faults_data;
             for (std::map<std::string, vector<std::pair<size_t, size_t> > >::iterator
#endif /* NEW */
                  it = faults_data_.begin(); it != faults_data_.end(); ++it)
             {
               fault_name = (*it).first;

               if (fault_name.find("sideA") != std::string::npos)
               {
                 // recovering the original name of the fault
                 fault_name_prefix = fault_name;
                 fault_name_prefix.erase(fault_name_prefix.find("_sideA"), 6);
                 fault_name_a = fault_name;

                 // creates a backup copy of the contents of fault sideA
                 std::set<std::set<size_t> > elmts_a;
                 std::vector<std::pair<size_t, size_t> >& fault_elmts_a((*it).second);
                 for (std::vector<std::pair<size_t, size_t> >::const_iterator
                      vit = fault_elmts_a.begin(); vit != fault_elmts_a.end(); ++vit)
                 {
                   /// assign data from provided face
                   hexa_cell_id = (*vit).first;
                   hexa_face_id = (*vit).second;

                   std::set<size_t> elmnt;
                   for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
                     elmnt.insert(poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(face_nodes[hexa_face_id][nid]));
                   elmts_a.insert(elmnt);
                 }

                 // finds the other side of the fault
                 fault_name_b = fault_name_prefix;
                 fault_name_b += "_sideB";
                 fit = faults_data_.find(fault_name_b);
                 std::set<std::set<size_t> > elmts_b;
                 if (fit != faults_data_.end())
                 {
                   std::vector<std::pair<size_t, size_t> >& fault_elmts_b((*fit).second);
                   for (std::vector<std::pair<size_t, size_t> >::const_iterator
                        vit = fault_elmts_b.begin(); vit != fault_elmts_b.end(); ++vit)
                   {
                     /// assign data from provided face
                     hexa_cell_id = (*vit).first;
                     hexa_face_id = (*vit).second;

                     std::set<size_t> elmnt;
                     for (size_t nid = 0; nid < quad_fem_nodes; ++nid)
                       elmnt.insert(poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(face_nodes[hexa_face_id][nid]));
                     elmts_b.insert(elmnt);
                   }
                 }

                 // if the 2 fault sides are the same, the separation in sides a and b is removed?
                 if (elmts_a == elmts_b)
                 {
                   temp_faults_data.insert(std::make_pair(fault_name_prefix, (*it).second));
                 }
                 else
                 {
                   temp_faults_data.insert(std::make_pair(fault_name_a, (*it).second));
                   temp_faults_data.insert(std::make_pair(fault_name_b, (*fit).second));
                 }
               }
             }

             //  ???? here we go again?
             /// reassign faults data
             faults_data_.clear();
             faults_data_ = temp_faults_data;
             temp_faults_data.clear();
             }


             // DEBUGGING - appears to recreate faults data exactly = duplicating it
             //std::cerr <<"\n\n twice recreated fault data:\n";
             //for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
             //     it = faults_data_.begin(); it != faults_data_.end(); ++it ) {
             //    std::cerr <<"\n\n"<< (*it).first <<" (cell-ID/face-ID pairs): ";
             //    for ( std::vector<std::pair<size_t,size_t> >::const_iterator
             //          fit=(*it).second.begin(); fit!=(*it).second.end(); fit++ )
             //      std::cerr << (*fit).first <<","<< (*fit).second <<" ";
             // }
             //std::cerr <<"\n";



             /* writing results to ModelTopology and VSet

              - for each fault two regions of surface elements will be created, being called _sideA and _sideB of the
              the given fault name

              - the lower-dimensional fault-side elements are being pushed back into the VSet
              */
             size_t                cell_id(vset.Elements());
             std::vector<size_t>   fault_elmts;
             std::set<std::string> fault_fem_types;
             size_t                num_sub_faces;

             for (std::map<std::string, std::vector<std::pair<size_t, size_t> > >::const_iterator
                  it = faults_data_.begin(); it != faults_data_.end(); ++it)
             {
               fault_name = (*it).first;
               fault_elmts.clear();
               fault_fem_types.clear();
               for (std::vector<std::pair<size_t, size_t> >::const_iterator
                    fit = (*it).second.begin(); fit != (*it).second.end(); ++fit)
               {
                 hexa_cell_id = (*fit).first;
                 if (cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_)
                 {
                   hexa_face_id = (*fit).second;
                   assert(hexa_face_id <= 5);
                   // DEBUGGING: NODE ORDER OF FACES IS DIFFERENT FROM DEFINITION FOR CORNER-POINT CELL
                   //poly[ hexa_cell_id ].Out();

                   // used to always evaluate to zero; now replaced by subclass method
                   num_sub_faces = poly[hexa_cell_id].GetNumSubFaces(hexa_face_id);

                   for (size_t fid = 0; fid < num_sub_faces; ++fid)
                   {
                     const size_t num_face_nodes(poly[hexa_cell_id].GetNumFaceNodes(hexa_face_id, fid));
                     if (num_face_nodes > 2U)
                     {
                       /// fill VSet
                       vset.ResizePlist(cell_id + 1, CSMP_ElementSpecifications::NodesPerElementOfType(
                                                                                                       poly[hexa_cell_id].GetFaceType(hexa_face_id, fid)));
                       vset.ResizeElementTypes(cell_id + 1);
                       for (size_t nid = 0U; nid < num_face_nodes; ++nid)
                         vset.Plist(cell_id, nid, poly[hexa_cell_id].GetFaceNodeGlobalId(hexa_face_id, fid, nid));
                       vset.ElementType(cell_id, poly[hexa_cell_id].GetFaceType(hexa_face_id, fid));
                       /// add fault elements
                       fault_elmts.push_back(cell_id);
                       fault_fem_types.insert(csmp::parseFiniteElementType(poly[hexa_cell_id].GetFaceType(hexa_face_id, fid)));
                       embedded_cells_[hexa_cell_id].push_back(std::make_pair(cell_id, poly[hexa_cell_id].GetFaceType(hexa_face_id, fid)));
                       ++cell_id;
                     }
                   }
                 }
               }

               // if a fault has been identified, its member elements are added as a region to the model topology
               if (!fault_elmts.empty())
               {
                 model_topology.AddRegion(fault_name.c_str(), fault_fem_types, fault_elmts);
                 faults.insert(fault_name.c_str());
               }
             }

             } // end







             void CornerPointGrid::EstablishWellRegions(csmp::VSet<3U>& vset,
                                                        csmp::ModelTopology& model_topology,
                                                        std::set<std::string>& wells,
                                                        const std::vector<CornerPointCell>& poly)
             {
               size_t                hexa_cell_id;
               size_t                nid1, nid2;
               size_t                cell_id(vset.Elements());
               std::string           well_name;
               std::vector<size_t>   well_elmts;
               std::set<std::string> well_fem_types;
               std::set<size_t>      edge_nodes;
               const csmp::CSMP_FEM_TYPE edge_fem_type(csmp::ISOPARAMETRIC_LINEAR_BAR);
               std::string           edge_fem_type_name(csmp::parseFiniteElementType(edge_fem_type));
               size_t                num_edge_nodes(2U);

               for (std::map<std::string, std::vector<std::pair<size_t, std::pair<size_t, size_t> > > >::const_iterator
                    it = well_face_path_.begin(); it != well_face_path_.end(); ++it)
               {
                 well_name = (*it).first;
                 well_elmts.clear();
                 well_fem_types.clear();
                 for (std::vector<std::pair<size_t, std::pair<size_t, size_t> > >::const_iterator
                      fit = (*it).second.begin(); fit != (*it).second.end(); ++fit)
                 {
                   hexa_cell_id = (*fit).first;
                   if (cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_)
                   {
                     /// fill VSet
                     const size_t num_well_sections(poly[hexa_cell_id].GetNumWellSections());
                     for (size_t eid = 0; eid < num_well_sections; ++eid)
                     {
                       vset.ResizePlist(cell_id + 1, num_edge_nodes);
                       vset.ResizeElementTypes(cell_id + 1);
                       for (size_t nid = 0U; nid < num_edge_nodes; ++nid)
                         vset.Plist(cell_id, nid, poly[hexa_cell_id].GetWellSectionNodeGlobalId(eid, nid));
                       vset.ElementType(cell_id, edge_fem_type);
                       /// add well elements
                       well_elmts.push_back(cell_id);
                       well_fem_types.insert(edge_fem_type_name);
                       embedded_cells_[hexa_cell_id].push_back(std::make_pair(cell_id, edge_fem_type));
                       ++cell_id;
                     }
                   }
                 }
                 if (!well_elmts.empty())
                 {
                   model_topology.AddRegion(well_name.c_str(), well_fem_types, well_elmts);
                   wells.insert(well_name.c_str());
                 }
               }

               for (std::map<std::string, std::vector<std::pair<size_t, std::pair<size_t, size_t> > > >::const_iterator
                    it = well_edge_path_.begin(); it != well_edge_path_.end(); ++it)
               {
                 well_name = (*it).first;
                 well_elmts.clear();
                 well_fem_types.clear();
                 for (std::vector<std::pair<size_t, std::pair<size_t, size_t> > >::const_iterator
                      fit = (*it).second.begin(); fit != (*it).second.end(); ++fit)
                 {
                   hexa_cell_id = (*fit).first;
                   if (cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_)
                   {
                     nid1 = (*fit).second.first;
                     nid2 = (*fit).second.second;
                     edge_nodes.clear();
                     edge_nodes.insert(poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(nid1));
                     edge_nodes.insert(poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(nid2));
                     num_edge_nodes = edge_nodes.size();
                     if (num_edge_nodes == 2)
                     {
                       /// fill VSet
                       vset.ResizePlist(cell_id + 1, num_edge_nodes);
                       vset.ResizeElementTypes(cell_id + 1);
                       vset.Plist(cell_id, 0, poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(nid1));
                       vset.Plist(cell_id, 1, poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(nid2));
                       vset.ElementType(cell_id, edge_fem_type);
                       /// add well elements
                       well_elmts.push_back(cell_id);
                       well_fem_types.insert(edge_fem_type_name);
                       embedded_cells_[hexa_cell_id].push_back(std::make_pair(cell_id, edge_fem_type));
                       ++cell_id;
                     }
                   }
                 }
                 if (!well_elmts.empty())
                 {
                   model_topology.AddRegion(well_name.c_str(), well_fem_types, well_elmts);
                   wells.insert(well_name.c_str());
                 }
               }
             }
#endif

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



             std::vector<ScalarVariable>& CellCenteredGrid::GetCellSizes( size_t i )
        {
          if( i == 0 )
            return dx_;
          if( i == 1 )
            return dy_;
          return dz_;
        }



             void CellCenteredGrid::AssignDimensionX( size_t NX )
        {
          NX_ = NX;
        }



             void CellCenteredGrid::AssignDimensionY( size_t NY )
        {
          NY_ = NY;
        }



             void CellCenteredGrid::AssignDimensionZ( size_t NZ )
        {
          NZ_ = NZ;
        }



#if 0
             void CellCenteredGrid::AssignCellCoordinatesToPillars( std::vector<std::vector<Pillar> >& pillars )
        {
          const size_t NXY( NX_*NY_ );
          size_t cell_id(0);
          double64 x_offset(0.0);
          double64 y_offset(0.0);
          double64 z_offset(0.0);
          double64 dx(0.0);
          double64 dy(0.0);

          csmp::Point<3U> pt;
          pillars.resize( (NY_ + 1), std::vector<Pillar>( NX_ + 1 ) );

          for( size_t k = 0; k < NZ_; k++ )
          {
            /// top level
            y_offset = 0.0;
            for( size_t j = 0; j < NY_; j++ )
            {
              x_offset = 0.0;
              for( size_t i = 0; i < NX_; i++ )
              {
                cell_id = i + j*NX_ + k*NXY;

                dx = dx_[cell_id]();
                dy = dy_[cell_id]();

                z_offset = tops_[ cell_id ]();

                // assign point coordinates
                ///(0,0,0) NW,top
                pt[0] = x_offset;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][ i ].AddPoint( pt );
                ///(1,0,0) NE,top
                pt[0] = x_offset + dx;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][i+1].AddPoint( pt );
                ///(0,1,0) SW,top
                pt[0] = x_offset;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][ i ].AddPoint( pt );
                ///(1,1,0) SE,top
                pt[0] = x_offset + dx;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][i+1].AddPoint( pt );

                // assign new x offset
                x_offset += dx;
              }
              // assign new y offset ( take dy value from the last cell in x row )
              y_offset += dy;
            }
            /// bottom level
            y_offset = 0.0;
            for( size_t j = 0; j < NY_; j++ )
            {
              x_offset = 0.0;
              for( size_t i = 0; i < NX_; i++ )
              {
                cell_id = i + j*NX_ + k*NXY;

                dx = dx_[cell_id]();
                dy = dy_[cell_id]();
                double64 dz = dz_[cell_id]();

                z_offset = ( ( k != NZ_-1 ) ? tops_[ cell_id + NXY ]() : tops_[ cell_id ]() + dz );

                // assign point coordinates
                ///(0,0,1) NW,btm
                pt[0] = x_offset;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][ i ].AddPoint( pt );
                ///(1,0,1) NE,btm
                pt[0] = x_offset + dx;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][i+1].AddPoint( pt );
                ///(0,1,1) SW,btm
                pt[0] = x_offset;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][ i ].AddPoint( pt );
                ///(1,1,1) SE,btm
                pt[0] = x_offset + dx;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][i+1].AddPoint( pt );

                // assign new x offset
                x_offset += dx;
              }
              // assign new y offset ( take dy value from the last cell in x row )
              y_offset += dy;
            }
          }

          for( size_t j = 0; j <= NY_; j++ )
            for( size_t i = 0; i <= NX_; i++ )
              pillars[ j ][ i ].AssignEnds( );
        }
#endif


             // WELLS

             /// add well path based on symmetry assumption ( neighbouring cell defines the direction )
             /// by default well is assumed to be vertical
             void addWellPath( size_t NX, size_t NY, size_t NZ,
                              const std::string& well_name,
                              const std::vector<ijk>& cell_ids,
                              std::map<std::string,EclipseWellPath>& well_path )
        {
          const size_t NX_x_NY( NX*NY );
          std::map<int64_t,CORNER_POINT_CELL_FACE_INDEX> face_map;
          face_map.insert( std::make_pair( -1,       CORNER_POINT_CELL_FACE_Xminus ) );
          face_map.insert( std::make_pair( +1,       CORNER_POINT_CELL_FACE_Xplus  ) );
          face_map.insert( std::make_pair( -NX,      CORNER_POINT_CELL_FACE_Yminus ) );
          face_map.insert( std::make_pair( +NX,      CORNER_POINT_CELL_FACE_Yplus  ) );
          face_map.insert( std::make_pair( -NX_x_NY, CORNER_POINT_CELL_FACE_Zminus ) );
          face_map.insert( std::make_pair( +NX_x_NY, CORNER_POINT_CELL_FACE_Zplus  ) );

          EclipseWellPath  wpath;
          size_t num_cells( cell_ids.size() );
          for( size_t cid=0; cid<num_cells; ++cid )
          {
            wpath.path.emplace_back( cell_ids[cid], CORNER_POINT_CELL_FACE_Zminus, CORNER_POINT_CELL_FACE_Zplus);
          }
          if( wpath.path.size() > 1 )
          {
            size_t nid( wpath.path.size() - num_cells + 1 ); /// neighbour is a next cell
            for( size_t cid = 0; cid <(num_cells-1); ++cid, ++nid )
            {
              //              neighbor id          cell id
              int64_t face_id  = ((int64_t)wpath.path[ nid ].cell.i - (int)cell_ids[ cid ].i)
              + ((int64_t)wpath.path[ nid ].cell.j - (int64_t)cell_ids[ cid ].j) * (int64_t)NX
              + ((int64_t)wpath.path[ nid ].cell.k - (int64_t)cell_ids[ cid ].k) * (int64_t)NX_x_NY;
              if( nid-1 == 0 ) wpath.path[ nid-1 ].from = face_map[ -face_id ];
              wpath.path[ nid-1 ].to = face_map[ face_id ];
              wpath.path[ nid ].from = face_map[ -face_id ];
              wpath.path[ nid ].to   = face_map[ face_id ];
            }
          }
          well_path.emplace(well_name, wpath);
        }





             /// add well path with explicitly specified faces
             void addWellPath( const std::string& well_name,
                              const std::vector<size_t>& cell_ids,
                              const std::vector<std::pair<size_t,size_t> >& face_ids,
                              std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
        {
          /// temp data
          std::pair<size_t,size_t> direction;
          std::pair<size_t,std::pair<size_t,size_t> > path;
          size_t cell_id;
          size_t num_cells( cell_ids.size() );
          std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
          well_path.insert( std::make_pair( well_name, empty_path ) );
          std::vector<std::pair<size_t,std::pair<size_t,size_t> > >& wpath( well_path[ well_name ] );
          for( size_t cid=0; cid<num_cells; ++cid )
          {
            cell_id = cell_ids[ cid ];
            direction.first  = face_ids[ cid ].first;
            direction.second = face_ids[ cid ].second;
            path.first  = cell_id;
            path.second = direction;
            wpath.push_back( path );
          }
        }




             /// add well path with explicitly specified faces ( same for all cells )
             void addWellPath( const std::string& well_name,
                              const std::vector<size_t>& cell_ids,
                              std::pair<size_t,size_t> face_id,
                              std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
        {
          /// temp data
          std::pair<size_t,size_t> direction;
          std::pair<size_t,std::pair<size_t,size_t> > path;
          size_t cell_id;
          size_t num_cells( cell_ids.size() );
          std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
          well_path.insert( std::make_pair( well_name, empty_path ) );
          std::vector<std::pair<size_t,std::pair<size_t,size_t> > >& wpath( well_path[ well_name ] );
          for( size_t i=0; i<num_cells; ++i )
          {
            cell_id = cell_ids[ i ];
            direction.first  = face_id.first;
            direction.second = face_id.second;
            path.first  = cell_id;
            path.second = direction;
            wpath.push_back( path );
          }
        }




             /// add well path with explicitly specified faces ( for single cell )
             void addWellPath( const std::string& well_name,
                              size_t cell_id,
                              std::pair<size_t,size_t> face_id,
                              std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
        {
          /// temp data
          std::pair<size_t,size_t> direction;
          std::pair<size_t,std::pair<size_t,size_t> > path;
          direction.first  = face_id.first;
          direction.second = face_id.second;
          path.first  = cell_id;
          path.second = direction;
          std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
          well_path.insert( std::make_pair( well_name, empty_path ) );
          well_path[ well_name ].push_back( path );
        }

  } // eclipse

} // end namespace csmp
