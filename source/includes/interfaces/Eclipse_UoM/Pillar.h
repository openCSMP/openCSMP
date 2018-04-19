#ifndef PILLAR_H
#define PILLAR_H

#include "Point.h"
#include <map>

namespace csmp {

namespace eclipse {

/**

@class Pillar  Pillar "Pillar.h"
@author R. Manasipov
@date 2015

Pillar objects are the main building blocks of corner point geometry.
They store the x,y,z coordinates of the corner points for each layer
of the grid. There is an origin point at the top far corner and 
the direction which marks the (downward) path.

Pillars are described by:

    Pillar parametric equation in 3D: r=a+t*b, where
    'a' - origin,
    'b' - direction,
     t  - scalar parameter

Pillars stack up polygonal cells.
Pillars are managed independently of the shape of the cell.
They are essentially lines that connect sub-vertical cell edges imposing the 
semi-structured (unstructured in plane, but structured in cross-section)
shape of the pillar grid.

Pillar objects are defined via COORD keyword
Pillar objects are subdivided by ZCORN

*/
class Pillar {
  public:

    Pillar();
    Pillar( const csmp::Point<3U>& pStart, const csmp::Point<3U>& pEnd );

    Pillar( const Pillar& );
    Pillar& operator=( const Pillar& );

    ~Pillar();

    /// assignment
    void AssignAttachedCellId( size_t i, size_t j, size_t offset );
    void AssignEnds( const csmp::Point<3U>& pStart, const csmp::Point<3U>& pEnd );
    void AssignEnds();

    /// resolve missing coordinates
    void ResolveXYcoords();
    void ResolveXYcoords( csmp::Point<3U>& r );
    csmp::Point<3U> ResolveXYcoords( double z );

    /// beads on the pillar
    size_t GetNumPoints() const;
  
    /// returns the starting point of the pillar at the smallest Z
    Point<3> Origin() const { return a_; }
  
    /// Z-coordinate value of the origin
    double64 OriginZ() const { return a_[2]; }
  
    /// end point of the pillar
    Point<3> Destination() const { return b_ + a_; }
  
    /// Z-coordinate value of pillar end point
    double64 DestinationZ() const { return (b_ + a_)[2]; }
  
    /// returns the downward-pointing direction vector of the (straight) pillar
    Point<3> Direction() const { return b_; }
  
    /// determines the XY coordinates of the supplied point that must store the Z coordinate and stores the point as a pillar bead
    void AddPoint( csmp::Point<3U> );
  
    const csmp::Point<3U>& GetPoint( size_t zid ) const;

    /// cell points
    size_t GetNumAttachedCells() const;
    const csmp::Point<3U>& GetPoint( size_t xid, size_t yid, size_t zid, size_t nid ) const;
    void   AssignPoint( size_t xid, size_t yid, size_t zid, size_t nid, const csmp::Point<3U>& pt );
  
    /// prints the data stored in the pillar to screen
    void Out() const;

  protected:

    size_t NodeIndexIncrementI( size_t nid ) const;
    size_t NodeIndexIncrementJ( size_t nid ) const;
    size_t NodeIndexIncrementK( size_t nid ) const;

  private:

    /// line equations
    csmp::Point<3U> a_;  /// origin of line
    csmp::Point<3U> b_;  /// direction of line

    /// cells data
    std::map<std::pair<size_t,size_t>,size_t> attached_cells_; ///< indices of attached cells (i,j)
    std::vector<csmp::Point<3U> >  points_;                    ///< points, called beads
};

} // eclipse

}// end namespace csmp

#endif

