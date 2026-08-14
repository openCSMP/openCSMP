#ifndef CSMP_BOX_H
#define CSMP_BOX_H

#include <set>
#include "Point.h"

namespace csmp {

template<uint32_t> class VSet;
template<uint32_t> class Node;
template<uint32_t> class Model;

/**
@file Box.h

@brief Anything to do with definitions and operations on rectangular and box-shaped models(2D and 3D that is).
This also includes stand-alone functions for the processing of Boundary flags and values.

@author S.K. Matthai
@date 1997

*/

/**
@addtogroup CSMPglobalEnums
@{
*/

/// fixed boundary identifiers; @attention do not alter numbering or sequence because it is used in iterations
static constexpr std::int8_t
    IRREGULAR_OUTSIDE{-1},
    LEFT_OUTSIDE{-2},  ///< model boundary flags
    RIGHT_OUTSIDE{-3},  ///< ...
    BOTTOM_OUTSIDE{-4},
    TOP_OUTSIDE{-5},
    FRONT_OUTSIDE{-6},
    BACK_OUTSIDE{-7},
    CNR_MIN{-8},      ///< Min-x, min-y, min-z
    CNR_X{-9},        ///< Max x, min y, z
    CNR_XY{-10},      ///< Max x, y, min z
    CNR_Y{-11},       ///< Max y, min x, z
    CNR_Z{-12},       ///< Max z, min x, y
    CNR_XZ{-13},      ///< Max x, z, min y
    CNR_MAX{-14},     ///< Max-x, max-y, max-z
    CNR_YZ{-15},      ///< Max y, z, Min x
    BACK_BOTTOM{-16}, ///< model edges: BACK and BOTTOM
    BACK_RIGHT{-17},  ///< BACK and RIGHT
    BACK_TOP{-18},    ///< BACK and TOP
    BACK_LEFT{-19},   ///< BACK and LEFT
    BOTTOM_LEFT{-20}, ///< BOTTOM and LEFT
    BOTTOM_RIGHT{-21},///< BOTTOM and RIGHT
    TOP_RIGHT{-22},   ///< TOP and RIGHT
    TOP_LEFT{-23},    ///< TOP and LEFT
    FRONT_BOTTOM{-24},///< FRONT and BOTTOM
    FRONT_RIGHT{-25}, ///< FRONT and RIGHT
    FRONT_TOP{-26},   ///< FRONT and TOP
    FRONT_LEFT{-27},  ///< FRONT and LEFT
    REGION_BOUNDARY{-28}, 
    MULTIPLE_BOUNDARIES{-29};

/// @enum BOX_BOUNDARY uniquely identifies placement of nodes on the boundary of a box-shaped model
enum BOX_BOUNDARY : std::int8_t {
  NOT       = 0,                  ///< not located on a model boundary
  IRREGULAR = IRREGULAR_OUTSIDE,  ///< located on a not-specified outside boundary of model (usually in the bounding box)
  LEFT      = LEFT_OUTSIDE,
  RIGHT     = RIGHT_OUTSIDE,
  BOTTOM    = BOTTOM_OUTSIDE,
  TOP       = TOP_OUTSIDE,
  FRONT = FRONT_OUTSIDE,
  BACK  = BACK_OUTSIDE,
  CNR1  = CNR_MIN,
  CNR2  = CNR_X,
  CNR3  = CNR_XY,
  CNR4  = CNR_Y,
  CNR5  = CNR_Z,
  CNR6  = CNR_XZ,
  CNR7  = CNR_MAX,
  CNR8  = CNR_YZ,
  EDGE1 = BACK_BOTTOM,
  EDGE2 = BACK_RIGHT,
  EDGE3 = BACK_TOP,
  EDGE4 = BACK_LEFT,
  EDGE5 = BOTTOM_LEFT,
  EDGE6 = BOTTOM_RIGHT,
  EDGE7 = TOP_RIGHT,
  EDGE8 = TOP_LEFT,
  EDGE9 = FRONT_BOTTOM,
  EDGE10 = FRONT_RIGHT,
  EDGE11 = FRONT_TOP,
  EDGE12 = FRONT_LEFT,
  INTERNAL = REGION_BOUNDARY,    ///<  internal model boundary (usually inside bounding box, with neighbors on either side)
  MULTIPLE = MULTIPLE_BOUNDARIES ///<  can result when an element is at the front and back at the same time because model is only a single element thick or similar
};

/**
@}
*/

/**
@brief brings together all functionality needed to deal with Box shaped models whose boundaries
are identified by BOX_BOUNDARY flags.
*/
class Box {
public:
  void UnitNormalTo( BOX_BOUNDARY bdry, uint32_t dim,
                     std::vector<double>& nrml ) const;

};

constexpr BOX_BOUNDARY intToBOX_BOUNDARY( int8_t i ) noexcept;

/// turn enumeration into string
std::string  parseBoundary( BOX_BOUNDARY ) noexcept;

/// interprets enumeration from string
BOX_BOUNDARY  parseBoundary( const std::string& box_boundary_name ) noexcept;

/// returns whether supplied boundary flag refers to a model edge
bool isEdge( BOX_BOUNDARY ) noexcept;

/// returns whether flag is a corner
bool isCorner( BOX_BOUNDARY ) noexcept;

/// returns whether flag is a side of the model (LEFT, RIGHT, TOP...)
bool isSide( BOX_BOUNDARY ) noexcept;

/// can be located on the supplied side (any node on LEFT, RIGHT.., including edges and corners)
bool sharesSide(BOX_BOUNDARY b) noexcept;

/// returns whether supplied boundary flag is located at a model boundary edge
constexpr bool belongsToEdge( BOX_BOUNDARY edge, BOX_BOUNDARY )  noexcept;

/// returns whether supplied boundary flag is located at a model boundary side
constexpr bool belongsToSide( BOX_BOUNDARY side, BOX_BOUNDARY )  noexcept;

/// returns whether boundary flag belongs to boundary LEFT
bool isLEFT( BOX_BOUNDARY ) noexcept;
/// returns whether boundary flag belongs to boundary RIGHT
bool isRIGHT( BOX_BOUNDARY ) noexcept;
/// returns whether boundary flag belongs to boundary TOP
bool isTOP( BOX_BOUNDARY ) noexcept;
/// returns whether boundary flag belongs to boundary BOTTOM
bool isBOTTOM( BOX_BOUNDARY ) noexcept;
/// returns whether boundary flag belongs to boundary FRONT (3D only)
bool isFRONT( BOX_BOUNDARY ) noexcept;
/// returns whether boundary flag belongs to boundary BACK (3D only)
bool isBACK( BOX_BOUNDARY ) noexcept;

/// returns true also for edges and corners
bool canBeIRREGULAR( BOX_BOUNDARY ) noexcept;

std::array<BOX_BOUNDARY,2> facesOfEdge(BOX_BOUNDARY) noexcept;
std::array<BOX_BOUNDARY,3> facesOfCorner(BOX_BOUNDARY) noexcept;
std::set<BOX_BOUNDARY> impliedSides(const std::set<BOX_BOUNDARY>& ) noexcept;
BOX_BOUNDARY uniqueEdgeOrNot(const std::set<BOX_BOUNDARY>&) noexcept;
BOX_BOUNDARY uniqueImpliedSideOrNot(const std::set<BOX_BOUNDARY>& ) noexcept;

/// returns nodes that are either flagged  CNR1 or CNR2 from one-dimensional model
Node<1U>* const cornerFlaggedNode( Model<1U>&, BOX_BOUNDARY );

/// is the cell located at a model boundary? - for line or surface element in 3D models  this is inferred when all their nodes are on boundary; for volume elements all nodes of one face must be
template<uint32_t dim, template<uint32_t> class CELL>
bool atBoundary( const CELL<dim>* const );

/// infers from cell type and flags of the corner nodes of the face, on which BOX_BUNDARY the cell face lies on including INTERNAL boundaries
template<uint32_t dim, template<uint32_t> class CELL>
BOX_BOUNDARY atBoundary( const CELL<dim>* const, uint32_t boundary_face );

/// master function that does all the work to identify which boundary a lower dimensional element is located on based on its node flags
template<uint32_t dim>
BOX_BOUNDARY atBoundary( const std::set<BOX_BOUNDARY>& node_flags, uint32_t boundary_face_corner_nodes );

/// prints a summary of the current flags of the nodes and elements to screen.
template<uint32_t dim> void printBoxBoundaryFlags( const Model<dim>& );

/// gives the extreme coordinates of the bounding box of the BOUNDARY (not the model!) - 1D model
void boundaryMinMaxCoordinates( BOX_BOUNDARY,
                                csmp::Point<1U>& model_coord_min,
                                csmp::Point<1U>& model_coord_max );

/// gives the extreme coordinates of the bounding box of the BOUNDARY (not the model!) - 2D model
void boundaryMinMaxCoordinates( BOX_BOUNDARY, csmp::Point<2U>&, csmp::Point<2U>& ) noexcept;

/// gives the extreme coordinates of the bounding box of the BOUNDARY (not the model!) - 3D model
void boundaryMinMaxCoordinates( BOX_BOUNDARY, csmp::Point<3U>&, csmp::Point<3U>& ) noexcept;

/// compares the supplied string with valid BOX_BOUNDARY classifications; returns true if it is among them
bool isDiagnosticBoxBoundaryClassifier( const std::string& ) noexcept;

/// Applies box boundary flags assuming that the model stored in the VSet is box-shaped with the sides aligned with the coordinate axes.
template<uint32_t dim>
void establishBoundaryFlagsForBoxModel( VSet<dim>& vset, double tolerance );

/// using existing Boundary objects, establish node and element BOX_BOUNDARY flagging
void recreateBoxBoundaryFlags( Model<1U>& );
void recreateBoxBoundaryFlags( Model<2U>& );
void recreateBoxBoundaryFlags( Model<3U>& );

/// boundary flags for 2D models consisting of quadrilaterals only; since these models are regular method is rather fast
void recreateBoxBoundaryFlagsForQuadrilateralModel( Model<2U>& );

/// boundary flags for 3D models consisting of hexahedra only; since these models are regular method is rather fast
void recreateBoxBoundaryFlagsForHexahedralModel( Model<3U>& );

/// permits to create variables values from BOX_BOUNDARY flag enumeration values
template<uint32_t dim>
void boxFlagsToVariable( Model<dim>&, const char* node_variable, const char* elmt_variable );

/// using the coordinate locations of the vertices defining the bounding rectangle; node flagging is performed; error prone method!
template<uint32_t dim>
void flagCornerNodes( VSet<dim>&, double tol, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax );

/// Finds edges in boxed shaped model and gibes them a BOX_BOUDARY_FLAG
template<uint32_t dim>
void flagEdges( VSet<dim>&, double tol, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax );

/// returns true if the model contains some box boundary identifiers, but also irregular boundaries
bool hasAllSideBoundaries( const Model<2U>& );
bool hasAllSideBoundaries( const Model<3U>& );

/// returns true if the model contains the complete set of box boundary identifiers
bool isStrictlyBoxShaped( const Model<2U>& );
bool isStrictlyBoxShaped( const Model<3U>& );

/// returns the dim-2 edge  (if any) that lies between the sides of the box that are given as arguments
BOX_BOUNDARY  whichEdge( BOX_BOUNDARY side1, BOX_BOUNDARY side2 );

/// returns the corner of the BOX (1-6) that is shared by the side boundaries
BOX_BOUNDARY  whichCorner( BOX_BOUNDARY side1, BOX_BOUNDARY side2, BOX_BOUNDARY side3 );

/// return which boundary the edge or lower-dim face with the two end-node flags is on
BOX_BOUNDARY  whichBoundary( BOX_BOUNDARY node_flag1, BOX_BOUNDARY node_flag2 );


/// reports the range of property values on the nodes flagged with the BOX_BOUNDARY identifier
template<uint32_t dim>
void boxBoundaryPropertyRange( const Model<dim>& sg, BOX_BOUNDARY boundary,
                               const char* node_property, double& bmin, double& bmax );

/// linear interpolation of values to point between 2 points in 1D (use for assigning property gradients along straight boundaries)
double linearInterpolate( const std::pair<Point<1U>, double>&,
                          const std::pair<Point<1U>, double>&,
                          const Point<1U>& );

/// bilinear interpolation to point between 2 points in 2D
double linearInterpolate( const std::pair<Point<2U>, double>&,
                          const std::pair<Point<2U>, double>&,
                          const Point<2U>& );

/// 3D trilinear interpolation
double linearInterpolate( const std::pair<Point<3U>, double>&,
                          const std::pair<Point<3U>, double>&,
                          const Point<3U>& );

/// just a stub (ignores the bilinear)
double bilinearInterpolate( uint32_t idx_x, uint32_t idx_y,
                            const Point<1U>& xy1,
                            const Point<1U>& xy2,
                            const Point<1U>& coord,
                            double p1, double p2, double, double );

double bilinearInterpolate( uint32_t idx_x, uint32_t idx_y,
                            const Point<2U>& xy1,
                            const Point<2U>& xy2,
                            const Point<2U>& coord,
                            // val@x0,y0  val@x1,y0  val@x1,y1  val@x0,y1
                            double p1, double p2, double p3, double p4 );

double bilinearInterpolate( uint32_t idx_x, uint32_t idx_y,
                            const Point<3U>& xy1,
                            const Point<3U>& xy2,
                            const Point<3U>& coord,
                            // val@x0,y0  val@x1,y0  val@x1,y1  val@x0,y1
                            double p1, double p2, double p3, double p4 );





// ======================================================================
//
//          INLINE OF CONSTEXPR FUNCTIONS (LOOKUP TABLES
//
// ======================================================================



/**
   Parses the BOX_BOUNDARY identifier (see Box.h). If the boundary flag cannot be resolved a value of NOT is returned if it is positive and IRREGULAR if negative.
   
   @attention use this only to convert ANSYS or outside-of CSMP generated flags of the enlisted names into BOX_BOUNDARY enums
*/
inline constexpr BOX_BOUNDARY intToBOX_BOUNDARY(std::int8_t i) noexcept
 {
    // valid range: NOT (0) or any negative value down to MULTIPLE (-29)
    if (i == 0)
        return BOX_BOUNDARY::NOT;

    if (i <= -1 && i >= MULTIPLE_BOUNDARIES)
        return static_cast<BOX_BOUNDARY>(i);

    // fallback for out-of-range values
    return BOX_BOUNDARY::NOT;
 }


/**
Returns true if @param bd belongs to the corresponding edge.
Example: CNR1 and CNR2 belong to EDGE1.
*/
inline constexpr bool belongsToEdge(BOX_BOUNDARY edge, BOX_BOUNDARY bd) noexcept {
    // Mapping each edge to its associated flags
    constexpr std::array<std::pair<BOX_BOUNDARY, std::array<BOX_BOUNDARY,3>>, 12> edge_map{{
        {EDGE1,  {EDGE1, CNR1, CNR2}},
        {EDGE2,  {EDGE2, CNR2, CNR3}},
        {EDGE3,  {EDGE3, CNR3, CNR4}},
        {EDGE4,  {EDGE4, CNR1, CNR4}},
        {EDGE5,  {EDGE5, CNR1, CNR5}},
        {EDGE6,  {EDGE6, CNR2, CNR6}},
        {EDGE7,  {EDGE7, CNR3, CNR7}},
        {EDGE8,  {EDGE8, CNR4, CNR8}},
        {EDGE9,  {EDGE9, CNR5, CNR6}},
        {EDGE10, {EDGE10, CNR6, CNR7}},
        {EDGE11, {EDGE11, CNR7, CNR8}},
        {EDGE12, {EDGE12, CNR5, CNR8}}
    }};

    for (const auto& p : edge_map) {
        if (p.first == edge) {
            for (const auto f : p.second) {
                if (f == bd) return true;
            }
            return false;
        }
    }
    return false;
}


/**
Returns true if @param bd belongs to the corresponding side of the box-shaped model.
Example: CNR1 - CNR4 belong to BACK.
*/
inline constexpr bool belongsToSide(BOX_BOUNDARY side, BOX_BOUNDARY bd) noexcept {
    switch (side) {
        case LEFT:      return isLEFT(bd);
        case RIGHT:     return isRIGHT(bd);
        case TOP:       return isTOP(bd);
        case BOTTOM:    return isBOTTOM(bd);
        case FRONT:     return isFRONT(bd);
        case BACK:      return isBACK(bd);
        case IRREGULAR: return bd == IRREGULAR;
        default:        return false;
    }
}



} // end namespace csmp

#endif

