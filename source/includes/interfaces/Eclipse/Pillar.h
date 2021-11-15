#ifndef PILLAR_UOM_H
#define PILLAR_UOM_H

#include "Point.h"
#include <map>

namespace csmp {

/**

@class Pillar  Pillar "Pillar.h"
@author R. Manasipov
@author A.J. Bromage
@date 2015-2017

*/
class Pillar {
  public:

    Pillar();
    Pillar( const csmp::Point<3U>& pStart, const csmp::Point<3U>& pEnd );

    ~Pillar();
  
    /// reserve capacity
    void Reserve( size_t zpoints );

    /// beads on the pillar
    size_t GetNumPoints() const;
  
    /// returns the starting point of the pillar at the smallest Z
    const Point<3>& StartPoint() const { return start_; }

    /// returns the starting point of the pillar at the largest Z
    const Point<3>& EndPoint() const { return end_; }

    void AddZCoord( double z );
  
    void SetZCoords( const std::vector<double>& zcoords );

    double GetZCoord( size_t zid ) const;

    csmp::Point<3U> GetPoint( size_t zid ) const;

    size_t FindPoint( double z ) const;

    void SetFirstNodeNum(size_t node_num);
  
    size_t FirstNodeNum() const;

    /// prints the data stored in the pillar to screen
    void Out() const;

  private:
    size_t firstNodeNum_;

    /// line equations
    csmp::Point<3U> start_;  /// start of line
    csmp::Point<3U> end_;  /// end of line

    /// cells data
    std::vector<double>  points_;                    ///< points, called beads
};


enum class ECLIPSE_CELL_CLASSIFICATION : uint8_t {
  ECLIPSE_CELL_HEXAHEDRON,					// 0000
  ECLIPSE_CELL_PYRAMIDS_310_312,			// 0001
  ECLIPSE_CELL_PYRAMIDS_201_203,			// 0010
  ECLIPSE_CELL_PRISM_23,					// 0011
  ECLIPSE_CELL_PYRAMIDS_130_132,			// 0100

  ECLIPSE_CELL_TETRAHEDRONS_130_132,		// 0101
  ECLIPSE_CELL_PRISM_12,					// 0110
  ECLIPSE_CELL_PYRAMID_0,					// 0111
  ECLIPSE_CELL_PYRAMIDS_021_023,			// 1000
  ECLIPSE_CELL_PRISM_03,					// 1001

  ECLIPSE_CELL_TETRAHEDRONS_021_023,		// 1010
  ECLIPSE_CELL_PYRAMID_1,					// 1011
  ECLIPSE_CELL_PRISM_01,					// 1100
  ECLIPSE_CELL_PYRAMID_2,					// 1101
  ECLIPSE_CELL_PYRAMID_3,					// 1110

  ECLIPSE_CELL_DEGENERATE					// 1111
};

/// A column cell, in Eclipse convention ordering
struct ColumnCell {
  size_t k;
  bool skew = false;
  ECLIPSE_CELL_CLASSIFICATION classification;
  size_t z[4][2];
};

/// A column
struct Column {
  std::vector<ColumnCell> cells_;
  
  Column(const std::vector<ColumnCell>& cells)
    : cells_(cells.begin(), cells.end())
  {
  }
};

}// end namespace csmp

#endif

