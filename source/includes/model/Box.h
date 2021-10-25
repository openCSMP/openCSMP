#ifndef CSMP_BOX_H
#define CSMP_BOX_H

#include <iostream>
#include <vector>
#include <deque>
#include "CSMP_number_types.h"
#include "Point.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Model;

/**
@file Box.h

@brief definitions for quadrilateral and box-shaped models(2D and 3D that is) captured by the class box.
@author S.K. Matthai
*/

/**
@addtogroup CSMPglobalEnums
@{
*/

/// fixed boundary identifiers; @attention do not alter numbering or sequence because it is used in iterations
enum : std::int8_t {
  IRREGULAR_OUTSIDE = -1,
  LEFT_OUTSIDE = -2,  ///< model boundary flags
  RIGHT_OUTSIDE = -3,  ///< ...
  BOTTOM_OUTSIDE = -4,
  TOP_OUTSIDE = -5,
  FRONT_OUTSIDE = -6,
  BACK_OUTSIDE = -7,
  CNR_MIN = -8,  ///< min-x, min-y, min-z
  CNR_MAX = -9,  ///< max-x, max-y, max-z
  CNR_MIN_MAXX = -10, ///< see users guide
  CNR_MIN_MAXXZ = -11,
  CNR_MIN_MAXZ = -12,
  CNR_MAX_MINXZ = -13,
  CNR_MAX_MAXX = -14,
  CNR_MAX_MAXZ = -15,
  BACK_BOTTOM = -16, ///< model edges: BACK and BOTTOM
  BACK_RIGHT = -17, ///< BACK and RIGHT
  BACK_TOP = -18, ///< BACK and TOP
  BACK_LEFT = -19, ///< BACK and LEFT
  BOTTOM_RIGHT = -20, ///< BOTTOM and RIGHT
  TOP_RIGHT = -21, ///< TOP and RIGHT
  TOP_LEFT = -22, ///< TOP and LEFT
  BOTTOM_LEFT = -23, ///< BOTTOM and LEFT
  FRONT_BOTTOM = -24, ///< FRONT and BOTTOM
  FRONT_RIGHT = -25, ///< FRONT and RIGHT
  FRONT_TOP = -26, ///< FRONT and TOP
  FRONT_LEFT = -27, ///< FRONT and LEFT
  REGION_BOUNDARY = -28, ///
  MULTIPLE_BOUNDARIES = -29
};

/// @enum BOX_BOUNDARY uniquely identifies placement of nodes on the boundary of a box-shaped model
enum BOX_BOUNDARY : std::int8_t {
  NOT       = 0,                  ///< not located on a model boundary
  IRREGULAR = IRREGULAR_OUTSIDE,  ///< located on a not-specified outside boundary of model (usually in the bounding box)
  TOP       = TOP_OUTSIDE,
  BOTTOM    = BOTTOM_OUTSIDE,
  LEFT   = LEFT_OUTSIDE,
  RIGHT = RIGHT_OUTSIDE,
  FRONT = FRONT_OUTSIDE,
  BACK  = BACK_OUTSIDE,
  CNR1  = CNR_MIN,
  CNR2  = CNR_MIN_MAXX,
  CNR3  = CNR_MAX_MAXX,
  CNR4  = CNR_MAX_MINXZ,
  CNR5  = CNR_MIN_MAXZ,
  CNR6  = CNR_MIN_MAXXZ,
  CNR7  = CNR_MAX,
  CNR8  = CNR_MAX_MAXZ,
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
  void UnitNormalTo( BOX_BOUNDARY bdry, size_t dim,
                     std::vector<double64>& nrml ) const;

};

BOX_BOUNDARY  intToBOX_BOUNDARY( int8_t i );

/// turn enumeration into string
std::string  parseBoundary( BOX_BOUNDARY );

/// interprets enumeration from string
BOX_BOUNDARY  parseBoundary( const std::string& box_boundary_name );

/// returns whether supplied boundary flag refers to a model edge
bool isEdge( BOX_BOUNDARY );

/// returns whether flag is a corner
bool isCorner( BOX_BOUNDARY );

/// returns whether flag is a side of the model
bool isSide( BOX_BOUNDARY );

/// returns whether supplied boundary flag is located at a model boundary edge
bool belongsToEdge( BOX_BOUNDARY edge, BOX_BOUNDARY );

/// returns whether supplied boundary flag is located at a model boundary side
bool belongsToSide( BOX_BOUNDARY side, BOX_BOUNDARY );

/// returns whether boundary flag belongs to boundary LEFT
bool isLEFT( BOX_BOUNDARY );
/// returns whether boundary flag belongs to boundary RIGHT
bool isRIGHT( BOX_BOUNDARY );
/// returns whether boundary flag belongs to boundary TOP
bool isTOP( BOX_BOUNDARY );
/// returns whether boundary flag belongs to boundary BOTTOM
bool isBOTTOM( BOX_BOUNDARY );
/// returns whether boundary flag belongs to boundary FRONT (3D only)
bool isFRONT( BOX_BOUNDARY );
/// returns whether boundary flag belongs to boundary BACK (3D only)
bool isBACK( BOX_BOUNDARY );

/// infers from node flags, and cell types, which boundary the element face lies on including INTERNAL ones
template<size_t dim, template<size_t> class CELL>
BOX_BOUNDARY atBoundary( const CELL<dim>* const, size_t boundary_face );

/// infers from node flags, missing neighbors, and cell types, which boundary the element lies on including INTERNAL ones; returns MULTIPLE if ambiguous; throws if inconsistent
template<size_t dim, template<size_t> class CELL>
BOX_BOUNDARY atBoundary( const CELL<dim>* const );

/// prints a summary of the current flags of the nodes and elements to screen.
template<size_t dim> void printBoxBoundaryFlags( const Model<dim>& );

/// gives the extreme coordinates of the bounding box of the BOUNDARY (not the model!) - 1D model
void boundaryMinMaxCoordinates( BOX_BOUNDARY,
                                csmp::Point<1U>& model_coord_min,
                                csmp::Point<1U>& model_coord_max );

/// gives the extreme coordinates of the bounding box of the BOUNDARY (not the model!) - 2D model
void boundaryMinMaxCoordinates( BOX_BOUNDARY, csmp::Point<2U>&, csmp::Point<2U>& );

/// gives the extreme coordinates of the bounding box of the BOUNDARY (not the model!) - 3D model
void boundaryMinMaxCoordinates( BOX_BOUNDARY, csmp::Point<3U>&, csmp::Point<3U>& );

/// compares the supplied string with valid BOX_BOUNDARY classifications; returns true if it is among them
bool isDiagnosticBoxBoundaryClassifier( const std::string& );

/// using existing Boundary objects, establish node and element BOX_BOUNDARY flagging
void recreateBoxBoundaryFlags( Model<1U>& );
void recreateBoxBoundaryFlags( Model<2U>& );
void recreateBoxBoundaryFlags( Model<3U>& );

/// boundary flags for 2D models consisting of quadrilaterals only; since these models are regular method is rather fast
void recreateBoxBoundaryFlagsForQuadrilateralModel( Model<2U>& );

/// boundary flags for 3D models consisting of hexahedra only; since these models are regular method is rather fast
void recreateBoxBoundaryFlagsForHexahedralModel( Model<3U>& );

/// permits to create variables values from BOX_BOUNDARY flag enumeration values
template<size_t dim>
void boxFlagsToVariable( Model<dim>&, const char* node_variable, const char* elmt_variable );

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
template<size_t dim>
void boxBoundaryPropertyRange( const Model<dim>& sg, BOX_BOUNDARY boundary,
                               const char* node_property, double64& bmin, double64& bmax );


} // end namespace csmp

#endif

