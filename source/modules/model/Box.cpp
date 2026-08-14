#include <numbers>
#include "Box.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "compareFloats.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// returns wether a node (or line element) lies on an edge of the model
bool isEdge(BOX_BOUNDARY b) noexcept {
  switch (b) {
    case EDGE1: case EDGE2: case EDGE3: case EDGE4:
    case EDGE5: case EDGE6: case EDGE7: case EDGE8:
    case EDGE9: case EDGE10: case EDGE11: case EDGE12:
      return true;
    default: return false;
  }
}


/// of rectangular (brick-shaped) model
bool isCorner(BOX_BOUNDARY b) noexcept {
  switch (b) {
    case CNR1: case CNR2: case CNR3: case CNR4:
    case CNR5: case CNR6: case CNR7: case CNR8:
      return true;
    default: return false;
  }
}


/// of rectangular (brick-shaped) model
bool isSide(BOX_BOUNDARY b) noexcept {
  switch (b) {
    case BOTTOM: case RIGHT: case TOP: case LEFT: case FRONT: case BACK:
      return true;
    default: return false;
  }
}


/**
   of rectangular (brick-shaped) model
   
    @attention checks whether the supplied node flag may be part of the box boundary b
*/
bool sharesSide(BOX_BOUNDARY b) noexcept {
    return isLEFT(b) ||
           isRIGHT(b) ||
           isBOTTOM(b) ||
           isTOP(b) ||
           isBACK(b) ||
           isFRONT(b);
}


bool isLEFT(BOX_BOUNDARY bd) noexcept {
    constexpr std::array<BOX_BOUNDARY, 9> left_flags{
        LEFT, EDGE4, EDGE5, EDGE8, EDGE12, CNR1, CNR4, CNR5, CNR8
    };
    return std::any_of(left_flags.begin(), left_flags.end(),
                       [bd](BOX_BOUNDARY f){ return f == bd; });
}

bool isRIGHT(BOX_BOUNDARY bd) noexcept {
    constexpr std::array<BOX_BOUNDARY, 9> right_flags{
        RIGHT, EDGE2, EDGE6, EDGE7, EDGE10, CNR2, CNR3, CNR6, CNR7
    };
    return std::any_of(right_flags.begin(), right_flags.end(),
                       [bd](BOX_BOUNDARY f){ return f == bd; });
}

bool isBOTTOM(BOX_BOUNDARY bd) noexcept {
    constexpr std::array<BOX_BOUNDARY, 9> bottom_flags{
        BOTTOM, EDGE1, EDGE5, EDGE6, EDGE9, CNR1, CNR2, CNR5, CNR6
    };
    return std::any_of(bottom_flags.begin(), bottom_flags.end(),
                       [bd](BOX_BOUNDARY f){ return f == bd; });
}

bool isTOP(BOX_BOUNDARY bd) noexcept {
    constexpr std::array<BOX_BOUNDARY, 9> top_flags{
        TOP, EDGE3, EDGE7, EDGE8, EDGE11, CNR3, CNR4, CNR7, CNR8
    };
    return std::any_of(top_flags.begin(), top_flags.end(),
                       [bd](BOX_BOUNDARY f){ return f == bd; });
}

bool isFRONT(BOX_BOUNDARY bd) noexcept {
    constexpr std::array<BOX_BOUNDARY, 9> front_flags{
        FRONT, EDGE9, EDGE10, EDGE11, EDGE12, CNR5, CNR6, CNR7, CNR8
    };
    return std::any_of(front_flags.begin(), front_flags.end(),
                       [bd](BOX_BOUNDARY f){ return f == bd; });
}

bool isBACK(BOX_BOUNDARY bd) noexcept {
    constexpr std::array<BOX_BOUNDARY, 9> back_flags{
        BACK, EDGE1, EDGE2, EDGE3, EDGE4, CNR1, CNR2, CNR3, CNR4
    };
    return std::any_of(back_flags.begin(), back_flags.end(),
                       [bd](BOX_BOUNDARY f){ return f == bd; });
}



/**
     For irregular shaped boundaries that were created from surfaces
     @todo figure out whether edges and corners should be included?
*/
bool canBeIRREGULAR( BOX_BOUNDARY bd ) noexcept
{
  if ( bd == IRREGULAR ) return true;
  if ( isEdge( bd ) )    return true;
  if ( isCorner( bd ) )  return true;
  return false;
}



// ---------- Tables: faces implied by each edge/corner (per your convention) ----------

/// checked: correct
constexpr std::array<BOX_BOUNDARY, 2> sidesOfEdge(BOX_BOUNDARY e) noexcept {
    switch(e) {
        case EDGE1: return {BACK, BOTTOM};
        case EDGE2: return {BACK, RIGHT};
        case EDGE3: return {BACK, TOP};
        case EDGE4: return {BACK, LEFT};
        case EDGE5: return {BOTTOM, LEFT};
        case EDGE6: return {BOTTOM, RIGHT};
        case EDGE7: return {TOP, RIGHT};
        case EDGE8: return {TOP, LEFT};
        case EDGE9: return {FRONT, BOTTOM};
        case EDGE10: return {FRONT, RIGHT};
        case EDGE11: return {FRONT, TOP};
        case EDGE12: return {FRONT, LEFT};
        default: return {IRREGULAR, IRREGULAR};
    }
}

/** Corners (three-face intersections):
 
    C1 (0,0,0): LEFT,  BOTTOM, BACK
    C2 (3,0,0): RIGHT, BOTTOM, BACK
    C3 (3,3,0): RIGHT, TOP,    BACK
    C4 (0,3,0): LEFT,  TOP,    BACK
    C5 (0,0,3): LEFT,  BOTTOM, FRONT
    C6 (3,0,3): RIGHT, BOTTOM, FRONT
    C7 (3,3,3): RIGHT, TOP,    FRONT
    C8 (0,3,3): LEFT,  TOP,    FRONT
*/
std::array<BOX_BOUNDARY,3> facesOfCorner(BOX_BOUNDARY c) noexcept {
  switch (c) {
    case CNR1: return { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) };
    case CNR2: return { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) };
    case CNR3: return { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(TOP_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) };
    case CNR4: return { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(TOP_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) };
    case CNR5: return { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) };
    case CNR6: return { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) };
    case CNR7: return { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(TOP_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) };
    case CNR8: return { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(TOP_OUTSIDE),
                        static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) };
    default:   return { static_cast<BOX_BOUNDARY>(IRREGULAR),
                        static_cast<BOX_BOUNDARY>(IRREGULAR),
                        static_cast<BOX_BOUNDARY>(IRREGULAR) };
  }
}



/** Box faces implied by a set of flags:
 
    - Add side flags directly
    - For edges, add their two faces
    - For corners, add their three faces
*/
std::set<BOX_BOUNDARY>  impliedSides(const std::set<BOX_BOUNDARY>& flags) noexcept
{
  std::set<BOX_BOUNDARY> faces;
  for (auto b : flags) {
        if (b == LEFT_OUTSIDE || b == RIGHT_OUTSIDE || b == BOTTOM_OUTSIDE ||
            b == TOP_OUTSIDE || b == FRONT_OUTSIDE || b == BACK_OUTSIDE) {
            faces.insert(b);
        } else if (isEdge(b)) {
            auto ff = sidesOfEdge(b);
            if (ff[0] != IRREGULAR) faces.insert(ff[0]);
            if (ff[1] != IRREGULAR) faces.insert(ff[1]);
        } else if (isCorner(b)) {
            auto ff = facesOfCorner(b);
            if (ff[0] != IRREGULAR) faces.insert(ff[0]);
            if (ff[1] != IRREGULAR) faces.insert(ff[1]);
            if (ff[2] != IRREGULAR) faces.insert(ff[2]);
        }
    }
    if (faces.empty() && flags.count(IRREGULAR)) {
        faces.insert(IRREGULAR);
    }
    return faces;
}

// If exactly one face is implied (possibly via edges/corners), return it.
BOX_BOUNDARY uniqueImpliedSideOrNot(const std::set<BOX_BOUNDARY>& flags) noexcept {
  auto sides = impliedSides(flags);
  // Remove IRREGULAR if it snuck in
  sides.erase(IRREGULAR);
  if (sides.size() == 1) return *sides.begin();
  return IRREGULAR;
}


// If exactly one unique edge flag exists in the set (ignoring corners), return it.
BOX_BOUNDARY uniqueEdgeOrNot(const std::set<BOX_BOUNDARY>& flags) noexcept {
  BOX_BOUNDARY found = IRREGULAR;
  for (auto b : flags) {
    if (!isEdge(b)) continue;
    if (found == IRREGULAR) {
      found = b;
    } else if (found != b) {
      return IRREGULAR; // more than one distinct edge present → ambiguous
    }
  }
  return found;
}






/**
    Finds the corner nodes of a boxed shaped model so that boundary conditions can be assigned.
    Function uses domain 'Model'  to avoid linear search.
    
     @return returns nodes that are either flagged  CNR1, CNR2 or INTERNAL from one-dimensional model
 */
Node<1U>* const cornerFlaggedNode( Model<1U>& model, BOX_BOUNDARY corner_flag )
 {
    // the corner nodes are the perimeter nodes of the model domain
    Region<1U>& model_domain = model.Region("Model");
    assert( model_domain.PerimeterNodes() >= 2 ); // at least the corners must be there
    
    auto nit = model_domain.PerimeterNodesBegin();
    while ( nit != model_domain.NodesEnd() ) {
         if ( (*nit)->AtBoundary() == corner_flag )
           return (*nit);
         nit++;
      }
      
    cerr <<"\n"<<"cornerFlaggedNode: corner node flagged "<< parseBoundary(corner_flag);
    cerr <<" was not found."<< endl;
    return nullptr;
 }





/**
   Parses the BOX_BOUNDARY identifier. If the boundary flag cannot be resolved a value of NOT is returned if it is positive and IRREGULAR if negative.
   
    SKM tested 2/12/2025
*/
std::string parseBoundary( BOX_BOUNDARY idx ) noexcept {
    // compile-time mapping of all BOX_BOUNDARY values to strings
    constexpr std::array<const char*,29> lookup{
        // 0 does not need to be part of this
        /* -1 */ "IRREGULAR",
        /* -2 */ "LEFT",
        /* -3 */ "RIGHT",
        /* -4 */ "BOTTOM",
        /* -5 */ "TOP",
        /* -6 */ "FRONT",
        /* -7 */ "BACK",
        /* -8 */ "CNR1",
        /* -9 */ "CNR2",
        /* -10 */ "CNR3",
        /* -11 */ "CNR4",
        /* -12 */ "CNR5",
        /* -13 */ "CNR6",
        /* -14 */ "CNR7",
        /* -15 */ "CNR8",
        /* -16 */ "EDGE1",
        /* -17 */ "EDGE2",
        /* -18 */ "EDGE3",
        /* -19 */ "EDGE4",
        /* -20 */ "EDGE5",
        /* -21 */ "EDGE6",
        /* -22 */ "EDGE7",
        /* -23 */ "EDGE8",
        /* -24 */ "EDGE9",
        /* -25 */ "EDGE10",
        /* -26 */ "EDGE11",
        /* -27 */ "EDGE12",
        /* -28 */ "INTERNAL",
        /* -29 */ "MULTIPLE",
    };

    // Rule: any non-negative → "NOT"
    if (idx >= 0) return "NOT";

    // Convert negative index to table offset
    const int i = -idx;

    // Out-of-range negative → "NOT" (or could choose another sentinel)
    if ( i > static_cast<int>(lookup.size()) )
        return "NOT";

    return lookup[ static_cast<size_t>(i-1) ];
}


/// text to boundary enum
BOX_BOUNDARY parseBoundary(const std::string& s) noexcept {
    using sv = std::string_view;

    // Compile-time table of string → enum
    constexpr std::array lookup{
        std::pair{sv{"NOT"}, NOT},
        std::pair{sv{"IRREGULAR"}, IRREGULAR},
        std::pair{sv{"LEFT"}, LEFT},
        std::pair{sv{"RIGHT"}, RIGHT},
        std::pair{sv{"BOTTOM"}, BOTTOM},
        std::pair{sv{"TOP"}, TOP},
        std::pair{sv{"FRONT"}, FRONT},
        std::pair{sv{"BACK"}, BACK},
        std::pair{sv{"CNR1"}, CNR1},
        std::pair{sv{"CNR2"}, CNR2},
        std::pair{sv{"CNR3"}, CNR3},
        std::pair{sv{"CNR4"}, CNR4},
        std::pair{sv{"CNR5"}, CNR5},
        std::pair{sv{"CNR6"}, CNR6},
        std::pair{sv{"CNR7"}, CNR7},
        std::pair{sv{"CNR8"}, CNR8},
        std::pair{sv{"EDGE1"}, EDGE1},
        std::pair{sv{"EDGE2"}, EDGE2},
        std::pair{sv{"EDGE3"}, EDGE3},
        std::pair{sv{"EDGE4"}, EDGE4},
        std::pair{sv{"EDGE5"}, EDGE5},
        std::pair{sv{"EDGE6"}, EDGE6},
        std::pair{sv{"EDGE7"}, EDGE7},
        std::pair{sv{"EDGE8"}, EDGE8},
        std::pair{sv{"EDGE9"}, EDGE9},
        std::pair{sv{"EDGE10"}, EDGE10},
        std::pair{sv{"EDGE11"}, EDGE11},
        std::pair{sv{"EDGE12"}, EDGE12},
        std::pair{sv{"INTERNAL"}, INTERNAL},
        std::pair{sv{"MULTIPLE"}, MULTIPLE}
    };

    // linear search (small table, compile-time friendly)
    for (const auto& p : lookup) {
        if (p.first == s)
            return p.second;
    }

    // fallback: NOT (your rule for invalid or non-existent flags)
    return NOT;
}





/**
Sets supplied model min, max coordinates to the coordinate value
range that applies to the specified model boundary.
*/
void boundaryMinMaxCoordinates( BOX_BOUNDARY boundary, csmp::Point<1U>&, csmp::Point<1U>& )
{ assert( !isCorner( boundary ) ); }


void boundaryMinMaxCoordinates( BOX_BOUNDARY boundary,
                                csmp::Point<2U>& model_coord_min,
                                csmp::Point<2U>& model_coord_max) noexcept
{
    assert(!isCorner(boundary));

    switch (boundary) {
        case LEFT:   model_coord_max[0] = model_coord_min[0]; return;
        case RIGHT:  model_coord_min[0] = model_coord_max[0]; return;
        case TOP:    model_coord_min[1] = model_coord_max[1]; return;
        case BOTTOM: model_coord_max[1] = model_coord_min[1]; return;
        default: 
            // Unrecognized 2D boundary; leave coords unchanged
            assert(false && "boundaryMinMaxCoordinates: invalid 2D boundary");
    }
}



void boundaryMinMaxCoordinates( BOX_BOUNDARY boundary,
                                csmp::Point<3U>& model_coord_min,
                                csmp::Point<3U>& model_coord_max) noexcept
{
    assert(!isCorner(boundary));

    switch (boundary) {
        // Sides
        case LEFT:   model_coord_max[0] = model_coord_min[0]; return;
        case RIGHT:  model_coord_min[0] = model_coord_max[0]; return;
        case BOTTOM: model_coord_max[1] = model_coord_min[1]; return;
        case TOP:    model_coord_min[1] = model_coord_max[1]; return;
        case BACK:   model_coord_max[2] = model_coord_min[2]; return;
        case FRONT:  model_coord_min[2] = model_coord_max[2]; return;

        // Edges in XY-plane
        case EDGE1:  model_coord_max[1] = model_coord_min[1]; model_coord_max[2] = model_coord_min[2]; return;
        case EDGE2:  model_coord_min[0] = model_coord_max[0]; model_coord_max[2] = model_coord_min[2]; return;
        case EDGE3:  model_coord_min[1] = model_coord_max[1]; model_coord_max[2] = model_coord_min[2]; return;
        case EDGE4:  model_coord_max[0] = model_coord_min[0]; model_coord_max[2] = model_coord_min[2]; return;

        // Edges along Z-axis
        case EDGE5:  model_coord_max[0] = model_coord_min[0]; model_coord_max[1] = model_coord_min[1]; return;
        case EDGE6:  model_coord_min[0] = model_coord_max[0]; model_coord_max[1] = model_coord_min[1]; return;
        case EDGE7:  model_coord_min[0] = model_coord_max[0]; model_coord_min[1] = model_coord_max[1]; return;
        case EDGE8:  model_coord_max[0] = model_coord_min[0]; model_coord_min[1] = model_coord_max[1]; return;

        // Edges in XY-plane at z=max
        case EDGE9:  model_coord_max[1] = model_coord_min[1]; model_coord_min[2] = model_coord_max[2]; return;
        case EDGE10: model_coord_min[0] = model_coord_max[0]; model_coord_min[2] = model_coord_max[2]; return;
        case EDGE11: model_coord_min[1] = model_coord_max[1]; model_coord_min[2] = model_coord_max[2]; return;
        case EDGE12: model_coord_max[0] = model_coord_min[0]; model_coord_min[2] = model_coord_max[2]; return;

        default:
            std::cerr << "\nboundaryMinMaxCoordinates(3D): boundary could not be parsed." << std::endl;
    }
} // end boundaryMinMaxCoordinates (3D)



/**
Returns a normal (vector) to the specified boundary of a box-shaped model.
*/
void Box::UnitNormalTo(BOX_BOUNDARY bdry, uint32_t dim, std::vector<double>& nrml) const
{
    nrml.resize(dim,0.);

    // C++23 const double s45 = sin(45.0 * numbers::pi / 180.0); // convert 45 deg to rad
    const double s45 = sin(45.0 * csmp::CSMP_PI / 180.0); // convert 45 deg to rad

    switch (dim) {
        case 1U:
            switch (bdry) {
                case LEFT:  nrml[0] = -1.; return;
                case RIGHT: nrml[0] =  1.; return;
                default:
                    throw csmp::Exception(ERROR, "Box::UnitNormalTo",
                                          "1D models only have LEFT and RIGHT boundaries.");
            }

        case 2U:
            switch (bdry) {
                case LEFT:   nrml = {-1.,  0.}; return;
                case RIGHT:  nrml = { 1.,  0.}; return;
                case TOP:    nrml = { 0.,  1.}; return;
                case BOTTOM: nrml = { 0., -1.}; return;
                default:
                    throw csmp::Exception(ERROR, "Box::UnitNormalTo",
                                          "2D models only have LEFT, RIGHT, TOP, and BOTTOM boundaries.");
            }

        case 3U:
            switch (bdry) {
                // surfaces
                case LEFT:   nrml = {-1.,  0.,  0.}; return;
                case RIGHT:  nrml = { 1.,  0.,  0.}; return;
                case BOTTOM: nrml = { 0., -1.,  0.}; return;
                case TOP:    nrml = { 0.,  1.,  0.}; return;
                case BACK:   nrml = { 0.,  0., -1.}; return;
                case FRONT:  nrml = { 0.,  0.,  1.}; return;

                // edges (back)
                case EDGE1:  nrml = {0., -s45, -s45}; return;
                case EDGE2:  nrml = {s45,  0., -s45}; return;
                case EDGE3:  nrml = {0.,  s45, -s45}; return;
                case EDGE4:  nrml = {-s45, 0., -s45}; return;

                // edges (center)
                case EDGE5:  nrml = { -s45, -s45, 0.}; return;
                case EDGE6:  nrml = {  s45, -s45, 0.}; return;
                case EDGE7:  nrml = {  s45,  s45, 0.}; return;
                case EDGE8:  nrml = { -s45,  s45, 0.}; return;

                // edges (front)
                case EDGE9:  nrml = {0., -s45, s45}; return;
                case EDGE10: nrml = { s45, 0., s45}; return;
                case EDGE11: nrml = {0., s45,  s45}; return;
                case EDGE12: nrml = {-s45, 0., s45}; return;

                default:
                    throw csmp::Exception(ERROR, "Box::UnitNormalTo",
                                          "Boundary flag could not be identified for 3D Box.");
            }

        default:
            throw csmp::Exception(ERROR, "Box::UnitNormalTo",
                                  "Invalid dimension. Only dim=1,2,3 supported.");
    }
} // end UnitNormalTo



/**
compares the supplied string with valid BOX_BOUNDARY classifications;
returns true if string is equivaled to identifiers, except for
NOT, INTERNAL, IRREGULAR

@author SKM (10/2/2016)
*/
bool isDiagnosticBoxBoundaryClassifier(const std::string& i) noexcept {
    static constexpr std::array<const char*, 26> valid_names = {
        "LEFT", "RIGHT", "BOTTOM", "TOP", "FRONT", "BACK",
        "CNR1", "CNR2", "CNR3", "CNR4", "CNR5", "CNR6", "CNR7", "CNR8",
        "EDGE1", "EDGE2", "EDGE3", "EDGE4", "EDGE5", "EDGE6",
        "EDGE7", "EDGE8", "EDGE9", "EDGE10", "EDGE11", "EDGE12"
    };
    
    return std::any_of(valid_names.begin(), valid_names.end(),
                       [&i](const char* name){ return i == name; });

} // end isBoxBoundaryClassifier



/**
prints a summary of the current flags of the nodes and elements to screen.
*/
template<uint32_t dim>
void printBoxBoundaryFlags( const Model<dim>& model )
 {
    const Region<dim>& modeldomain(model.Region("Model"));
    
    cout <<"\n\nprintBoxBoundaryFlags: Current flags and numbers assigned to node objects:\n";
    size_t equal_entries(0);
    multiset<BOX_BOUNDARY> node_flags;
    for ( auto nit=modeldomain.NodesBegin(); nit!=modeldomain.NodesEnd(); ++nit )
      node_flags.insert( (*nit)->AtBoundary() );

    multiset<BOX_BOUNDARY>::const_iterator it=node_flags.begin();

    if ( !node_flags.empty() ) {
          cout <<"\n\t"<<"node flags:"<< endl;
          while( it!=node_flags.end() ) {
               cout <<"\n\t\t"<< parseBoundary( (*it) ) <<": "<<  (equal_entries=node_flags.count( (*it) )) <<" nodes.";
               // avoiding printing of duplicates
               advance( it, equal_entries );
            }
          cout << endl;
      }
      
    multiset<BOX_BOUNDARY> elmt_flags;
    for ( auto eit=modeldomain.CellsBegin(); eit!=modeldomain.CellsEnd(); ++eit )
      for ( uint32_t i{0U}; i<(*eit)->Neighbors(); ++i )
        elmt_flags.insert( (*eit)->AtBoundary(i) );
     
    if ( !elmt_flags.empty() ) {
          cout <<"\n\t"<<"cell flags (not stored):"<< endl;
          it = elmt_flags.begin();
          while( it!=elmt_flags.end() ) {
               cout <<"\n\t\t"<< parseBoundary( (*it) ) <<": "<<  (equal_entries=elmt_flags.count( (*it) )) <<" elements.";
               advance( it, equal_entries );
            }
          cout << endl;
      }
 }

template void printBoxBoundaryFlags( const Model<1U>& );
template void printBoxBoundaryFlags( const Model<2U>& );
template void printBoxBoundaryFlags( const Model<3U>& );



/**
*/
void recreateBoxBoundaryFlags( Model<1U>& )
{
  throw csmp::Exception( ERROR, "recreateBoxBoundaryFlags(from Box.h)", "not implemented yet." );
}


/**
    Using already existing side csmp::Bboundary objects in the model,
    method recreates a corresponding, consistent box-boundary flagging.

    @attention the boundaries LEFT, RIGHT, TOP(or IRREGULAR), BOTTOM must be present
*/
void recreateBoxBoundaryFlags( Model<2U>& model )
{
  if ( distance( model.BoundariesBegin(), model.BoundariesEnd() ) == 0 )
    csmp::Exception( ERROR, "recreateBoxBoundaryFlags(2D):", "model contains no Boundary objects; nothing could be done." );

  // flagging the sides where they could be identified as box boundaries
  int counter( 0 );
  for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it )
  {
    const string boundary( (*it).first );
    if ( isDiagnosticBoxBoundaryClassifier( boundary ) ) {
      cout << "\nrecreateBoxBoundaryFlags(2D): found Boundary '" << boundary << "'";
      const BOX_BOUNDARY bflag( parseBoundary( boundary ) );
      for ( auto nit = (*it).second.NodesBegin(); nit != (*it).second.PerimeterNodesBegin(); ++nit )
        (*nit)->AtBoundary( bflag );
      counter++;
    }
  }
  cout << endl;

  if ( counter < 4 ) return;

  // if all (4) sides exist, the (4) edges can potentially be found by intersecting respective boundaries
  Boundary<2U>&  left( model.Boundary( "LEFT" ) );
  vector<Node<2U>*>  left_pnodes( left.PerimeterNodesBegin(), left.NodesEnd() );
  Boundary<2U>&  right( model.Boundary( "RIGHT" ) );
  vector<Node<2U>*>  right_pnodes( right.PerimeterNodesBegin(), right.NodesEnd() );
  Boundary<2U>&  top( model.Boundary( "TOP" ) );
  vector<Node<2U>*>  top_pnodes( top.PerimeterNodesBegin(), top.NodesEnd() );
  Boundary<2U>&  bottom( model.Boundary( "BOTTOM" ) );
  vector<Node<2U>*>  bottom_pnodes( bottom.PerimeterNodesBegin(), bottom.NodesEnd() );

  // finding the corners as the intersections between sides
  // CNR1
  vector<Node<2U>*> corner1;
  set_intersection( left_pnodes.begin(), left_pnodes.end(),
                    bottom_pnodes.begin(), bottom_pnodes.end(),
                    back_inserter( corner1 ) );
  if ( !corner1.empty() ) {
    (*corner1.begin())->AtBoundary( CNR1 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR1 (min_x,min_y)";
  }

  // CNR2
  vector<Node<2U>*> corner2;
  set_intersection( right_pnodes.begin(), right_pnodes.end(),
                    bottom_pnodes.begin(), bottom_pnodes.end(),
                    back_inserter( corner2 ) );
  if ( !corner2.empty() ) {
    (*corner2.begin())->AtBoundary( CNR2 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR2 (max_x,min_y)";
  }

  // CNR3
  vector<Node<2U>*> corner3;
  set_intersection( right_pnodes.begin(), right_pnodes.end(),
                    top_pnodes.begin(), top_pnodes.end(),
                    back_inserter( corner3 ) );
  if ( !corner3.empty() ) {
    (*corner3.begin())->AtBoundary( CNR3 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR3 (max_x,max_y)";
  }

  // CNR4
  vector<Node<2U>*> corner4;
  set_intersection( left_pnodes.begin(), left_pnodes.end(),
                    top_pnodes.begin(), top_pnodes.end(),
                    back_inserter( corner4 ) );
  if ( !corner4.empty() ) {
    (*corner4.begin())->AtBoundary( CNR4 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR4 (min_x,max_y)";
  }

} // end recreateBoxBoundaryFlags (2D version)




/**
re-establishes (in as much is possible) the Box boundary flagging of nodes and elements
using existing model Boundary objects, including sides, edges and corners.

@attention methods assumes that nodes and elements on the interior of the domain are already flagged as NOT.

@author SKM 10/3/2016
*/
void recreateBoxBoundaryFlags( Model<3>& model )
{
  if ( distance( model.BoundariesBegin(), model.BoundariesEnd() ) == 0 )
    csmp::Exception( ERROR, "recreateBoxBoundaryFlags:", "model contains no Boundary objects; nothing could be done." );

  // flagging the sides where box boundaries can be identified
  int counter( 0 );
  for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it )
  {
    const string boundary( (*it).first );
    if ( isDiagnosticBoxBoundaryClassifier( boundary ) ) {
      cout << "\nrecreateBoxBoundaryFlags: found Boundary '" << boundary << "'";
      const BOX_BOUNDARY bflag( parseBoundary( boundary ) );
      for ( auto nit = (*it).second.NodesBegin(); nit != (*it).second.PerimeterNodesBegin(); ++nit )
        (*nit)->AtBoundary( bflag );
      counter++;
    }
  }
  cout << endl;

  if ( counter < 6 ) return;

  // if all (6) sides exist, the (8) edges can potentially be found by intersecting respective boundaries
  Boundary<3>&  back( model.Boundary( "BACK" ) );
  vector<Node<3>*>  back_pnodes( back.PerimeterNodesBegin(), back.NodesEnd() ); // these are already sorted ranges
  assert( std::is_sorted( back_pnodes.begin(), back_pnodes.end() ) );

  Boundary<3>&  front( model.Boundary( "FRONT" ) );
  vector<Node<3>*>  front_pnodes( front.PerimeterNodesBegin(), front.NodesEnd() );
  assert( std::is_sorted( front_pnodes.begin(), front_pnodes.end() ) );

  Boundary<3>&  left( model.Boundary( "LEFT" ) );
  vector<Node<3>*>  left_pnodes( left.PerimeterNodesBegin(), left.NodesEnd() );
  assert( std::is_sorted( left_pnodes.begin(), left_pnodes.end() ) );

  Boundary<3>&  right( model.Boundary( "RIGHT" ) );
  vector<Node<3>*>  right_pnodes( right.PerimeterNodesBegin(), right.NodesEnd() );
  assert( std::is_sorted( right_pnodes.begin(), right_pnodes.end() ) );

  Boundary<3>&  top( model.Boundary( "TOP" ) );
  vector<Node<3>*>  top_pnodes( top.PerimeterNodesBegin(), top.NodesEnd() );
  assert( std::is_sorted( top_pnodes.begin(), top_pnodes.end() ) );

  Boundary<3>&  bottom( model.Boundary( "BOTTOM" ) );
  vector<Node<3>*>  bottom_pnodes( bottom.PerimeterNodesBegin(), bottom.NodesEnd() );
  assert( std::is_sorted( bottom_pnodes.begin(), bottom_pnodes.end() ) );

  // BACK_BOTTOM  = -16, ///< model edges: BACK and BOTTOM
  vector<Node<3>*> back_bottom;
  set_intersection( back_pnodes.begin(), back_pnodes.end(),
                    bottom_pnodes.begin(), bottom_pnodes.end(),
                    back_inserter( back_bottom ) );

  if ( !back_bottom.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE1' (BACK_BOTTOM)";
    for ( auto nit = back_bottom.begin(); nit != back_bottom.end(); ++nit )
      (*nit)->AtBoundary( EDGE1 );
  }

  // BACK_RIGHT   = -17, ///< BACK and RIGHT
  vector<Node<3>*> back_right;
  set_intersection( back_pnodes.begin(), back_pnodes.end(),
                    right_pnodes.begin(), right_pnodes.end(),
                    back_inserter( back_right ) );

  if ( !back_right.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE2' (BACK_RIGHT)";
    for ( auto nit = back_right.begin(); nit != back_right.end(); ++nit )
      (*nit)->AtBoundary( EDGE2 );
  }

  // BACK_TOP     = -18, ///< BACK and TOP
  vector<Node<3>*> back_top;
  set_intersection( back_pnodes.begin(), back_pnodes.end(),
                    top_pnodes.begin(), top_pnodes.end(),
                    back_inserter( back_top ) );

  if ( !back_top.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE3' (BACK_TOP)";
    for ( auto nit = back_top.begin(); nit != back_top.end(); ++nit )
      (*nit)->AtBoundary( EDGE3 );
  }

  // BACK_LEFT    = -19, ///< BACK and LEFT
  vector<Node<3>*> back_left;
  set_intersection( back_pnodes.begin(), back_pnodes.end(),
                    left_pnodes.begin(), left_pnodes.end(),
                    back_inserter( back_left ) );

  if ( !back_left.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE4' (BACK_LEFT)";
    for ( auto nit = back_left.begin(); nit != back_left.end(); ++nit )
      (*nit)->AtBoundary( EDGE4 );
  }

  // BOTTOM_LEFT  = -23, ///< BOTTOM and LEFT
  vector<Node<3>*> bottom_left;
  set_intersection( bottom_pnodes.begin(), bottom_pnodes.end(),
                    left_pnodes.begin(), left_pnodes.end(),
                    back_inserter( bottom_left ) );

  if ( !bottom_left.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE5' (BOTTOM_LEFT)";
    for ( auto nit = bottom_left.begin(); nit != bottom_left.end(); ++nit )
      (*nit)->AtBoundary( EDGE5 );
  }

  // BOTTOM_RIGHT = -20, ///< BOTTOM and RIGHT
  vector<Node<3>*> bottom_right;
  set_intersection( bottom_pnodes.begin(), bottom_pnodes.end(),
                    right_pnodes.begin(), right_pnodes.end(),
                    back_inserter( bottom_right ) );

  if ( !bottom_right.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE6' (BOTTOM_RIGHT)";
    for ( auto nit = bottom_right.begin(); nit != bottom_right.end(); ++nit )
      (*nit)->AtBoundary( EDGE6 );
  }

  // TOP_RIGHT    = -21, ///< TOP and RIGHT
  vector<Node<3>*> top_right;
  set_intersection( top_pnodes.begin(), top_pnodes.end(),
                    right_pnodes.begin(), right_pnodes.end(),
                    back_inserter( top_right ) );

  if ( !top_right.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE7' (TOP_RIGHT)";
    for ( auto nit = top_right.begin(); nit != top_right.end(); ++nit )
      (*nit)->AtBoundary( EDGE7 );
  }

  // TOP_LEFT     = -22, ///< TOP and LEFT
  vector<Node<3>*> top_left;
  set_intersection( top_pnodes.begin(), top_pnodes.end(),
                    left_pnodes.begin(), left_pnodes.end(),
                    back_inserter( top_left ) );

  if ( !top_left.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE8' (TOP_LEFT)";
    for ( auto nit = top_left.begin(); nit != top_left.end(); ++nit )
      (*nit)->AtBoundary( EDGE8 );
  }

  // FRONT_BOTTOM = -24, ///< FRONT and BOTTOM
  vector<Node<3>*> front_bottom;
  set_intersection( front_pnodes.begin(), front_pnodes.end(),
                    bottom_pnodes.begin(), bottom_pnodes.end(),
                    back_inserter( front_bottom ) );

  if ( !front_bottom.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE9' (FRONT_BOTTOM)";
    for ( auto nit = front_bottom.begin(); nit != front_bottom.end(); ++nit )
      (*nit)->AtBoundary( EDGE9 );
  }

  // FRONT_RIGHT  = -25, ///< FRONT and RIGHT
  vector<Node<3>*> front_right;
  set_intersection( front_pnodes.begin(), front_pnodes.end(),
                    right_pnodes.begin(), right_pnodes.end(),
                    back_inserter( front_right ) );

  if ( !front_right.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE10' (FRONT_RIGHT)";
    for ( auto nit = front_right.begin(); nit != front_right.end(); ++nit )
      (*nit)->AtBoundary( EDGE10 );
  }

  // FRONT_TOP    = -26, ///< FRONT and TOP
  vector<Node<3>*> front_top;
  set_intersection( front_pnodes.begin(), front_pnodes.end(),
                    top_pnodes.begin(), top_pnodes.end(),
                    back_inserter( front_top ) );

  if ( !front_top.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE11' (FRONT_TOP)";
    for ( auto nit = front_top.begin(); nit != front_top.end(); ++nit )
      (*nit)->AtBoundary( EDGE11 );
  }

  // FRONT_LEFT   = -27, ///< FRONT and LEFT
  vector<Node<3>*> front_left;
  set_intersection( front_pnodes.begin(), front_pnodes.end(),
                    left_pnodes.begin(), left_pnodes.end(),
                    back_inserter( front_left ) );

  if ( !front_left.empty() ) {
    cout << "\nrecreateBoxBoundaryFlags: found Boundary 'EDGE12'(FRONT_LEFT)";
    for ( auto nit = front_left.begin(); nit != front_left.end(); ++nit )
      (*nit)->AtBoundary( EDGE12 );
  }
  cout << endl;

  // and the model corners considering all possible permutations
  // CNR1
  vector<Node<3>*> corner1;
  set_intersection( back_left.begin(), back_left.end(),
                    back_bottom.begin(), back_bottom.end(),
                    back_inserter( corner1 ) );
  if ( corner1.empty() ) {
    set_intersection( back_left.begin(), back_left.end(),
                      bottom_left.begin(), bottom_left.end(),
                      back_inserter( corner1 ) );
  }
  if ( corner1.empty() ) {
    set_intersection( back_bottom.begin(), back_bottom.end(),
                      bottom_left.begin(), bottom_left.end(),
                      back_inserter( corner1 ) );
  }
  if ( !corner1.empty() ) {
    (*corner1.begin())->AtBoundary( CNR1 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR1 (min_x,min_y,min_z)";
  }

  // CNR2
  vector<Node<3>*> corner2;
  set_intersection( back_right.begin(), back_right.end(),
                    back_bottom.begin(), back_bottom.end(),
                    back_inserter( corner2 ) );
  if ( corner2.empty() ) {
    set_intersection( back_right.begin(), back_right.end(),
                      bottom_right.begin(), bottom_right.end(),
                      back_inserter( corner2 ) );
  }
  if ( corner2.empty() ) {
    set_intersection( back_bottom.begin(), back_bottom.end(),
                      bottom_right.begin(), bottom_right.end(),
                      back_inserter( corner2 ) );
  }
  if ( !corner2.empty() ) {
    (*corner2.begin())->AtBoundary( CNR2 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR2 (max_x,min_y,min_z)";
  }

  // CNR3
  vector<Node<3>*> corner3;
  set_intersection( back_right.begin(), back_right.end(),
                    back_top.begin(), back_top.end(),
                    back_inserter( corner3 ) );
  if ( corner3.empty() ) {
    set_intersection( back_right.begin(), back_right.end(),
                      top_right.begin(), top_right.end(),
                      back_inserter( corner3 ) );
  }
  if ( corner3.empty() ) {
    set_intersection( back_top.begin(), back_top.end(),
                      top_right.begin(), top_right.end(),
                      back_inserter( corner3 ) );
  }
  if ( !corner3.empty() ) {
    (*corner3.begin())->AtBoundary( CNR3 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR3 (max_x,max_y,min_z)";
  }

  // CNR4
  vector<Node<3>*> corner4;
  set_intersection( back_left.begin(), back_left.end(),
                    back_top.begin(), back_top.end(),
                    back_inserter( corner4 ) );
  if ( corner4.empty() ) {
    set_intersection( back_left.begin(), back_left.end(),
                      top_left.begin(), top_left.end(),
                      back_inserter( corner4 ) );
  }
  if ( corner4.empty() ) {
    set_intersection( top_left.begin(), top_left.end(),
                      back_top.begin(), back_top.end(),
                      back_inserter( corner4 ) );
  }
  if ( !corner4.empty() ) {
    (*corner4.begin())->AtBoundary( CNR4 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR4 (min_x,max_y,min_z)";
  }

  // corners on the FRONT of the model

  // CNR5
  vector<Node<3>*> corner5;
  set_intersection( bottom_left.begin(), bottom_left.end(),
                    front_left.begin(), front_left.end(),
                    back_inserter( corner5 ) );
  if ( corner5.empty() ) {
    set_intersection( front_left.begin(), front_left.end(),
                      bottom_left.begin(), bottom_left.end(),
                      back_inserter( corner4 ) );
  }
  if ( corner5.empty() ) {
    set_intersection( bottom_left.begin(), bottom_left.end(),
                      front_bottom.begin(), front_bottom.end(),
                      back_inserter( corner5 ) );
  }
  if ( !corner5.empty() ) {
    (*corner5.begin())->AtBoundary( CNR5 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR5 (min_x,min_y,max_z)";
  }

  // CNR6
  vector<Node<3>*> corner6;
  set_intersection( front_right.begin(), front_right.end(),
                    front_bottom.begin(), front_bottom.end(),
                    back_inserter( corner6 ) );
  if ( corner6.empty() ) {
    set_intersection( front_right.begin(), front_right.end(),
                      bottom_right.begin(), bottom_right.end(),
                      back_inserter( corner6 ) );
  }
  if ( corner6.empty() ) {
    set_intersection( front_bottom.begin(), front_bottom.end(),
                      bottom_right.begin(), bottom_right.end(),
                      back_inserter( corner6 ) );
  }
  if ( !corner6.empty() ) {
    (*corner6.begin())->AtBoundary( CNR6 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR6 (max_x,min_y,max_z)";
  }

  // CNR7
  vector<Node<3>*> corner7;
  set_intersection( front_right.begin(), front_right.end(),
                    front_top.begin(), front_top.end(),
                    back_inserter( corner7 ) );
  if ( corner7.empty() ) {
    set_intersection( front_right.begin(), front_right.end(),
                      top_right.begin(), top_right.end(),
                      back_inserter( corner7 ) );
  }
  if ( corner7.empty() ) {
    set_intersection( front_top.begin(), front_top.end(),
                      top_right.begin(), top_right.end(),
                      back_inserter( corner7 ) );
  }
  if ( !corner7.empty() ) {
    (*corner7.begin())->AtBoundary( CNR7 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR7 (max_x,max_y,max_z)";
  }

  // CNR8
  vector<Node<3>*> corner8;
  set_intersection( front_left.begin(), front_left.end(),
                    front_top.begin(), front_top.end(),
                    back_inserter( corner8 ) );
  if ( corner8.empty() ) {
    set_intersection( front_left.begin(), front_left.end(),
                      top_left.begin(), top_left.end(),
                      back_inserter( corner8 ) );
  }
  if ( corner8.empty() ) {
    set_intersection( top_left.begin(), top_left.end(),
                      front_top.begin(), front_top.end(),
                      back_inserter( corner8 ) );
  }
  if ( !corner8.empty() ) {
    (*corner8.begin())->AtBoundary( CNR8 );
    cout << "\nrecreateBoxBoundaryFlags: found CNR8 (min_x,max_y,max_z)";
  }
  cout << endl;

} // end recreateAtBoundaryFlags







/**
       Establishes the BOX_BOUNDARY flagging for a 2D model consisting out of quadrilateral elements. 
       
       @note this method is quick.
*/ 
void recreateBoxBoundaryFlagsForQuadrilateralModel( Model<2U>& model )
 { 
    Region<2U>&  modeldomain(model.Region("Model"));
    for ( vector<Element<2U>*>::const_iterator it=modeldomain.CellsBegin(); it!=modeldomain.CellsEnd(); ++it )
      {  // current version only works for linear quadrilaterals
         assert( (*it)->Nodes() <= 5 );
         
         if ( !isQuadrilateral( (*it)->FE_Type() ) )
           throw csmp::Exception( ERROR, "recreateBoxBoundaryFlagsForQuadrilateralModel", "this method only works for quadrilateral elements");
           
         for ( auto i{0U}; i<(*it)->Nodes(); ++i )
           (*it)->N(i)->AtBoundary( NOT );
        // BOTTOM 
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              (*it)->N(0)->AtBoundary( BOTTOM );
              (*it)->N(1)->AtBoundary( BOTTOM );
              continue;
           }
         // RIGHT 
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(2) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              (*it)->N(1)->AtBoundary( RIGHT );
              (*it)->N(2)->AtBoundary( RIGHT );
              continue;
           }
         // TOP 
         if ( (*it)->Neighbor(2) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              (*it)->N(2)->AtBoundary( TOP );
              (*it)->N(3)->AtBoundary( TOP );
              continue;
           }
         // LEFT 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              (*it)->N(0)->AtBoundary( LEFT );
              (*it)->N(3)->AtBoundary( LEFT );
              continue;
           }
         // CORNERS 
         // CNR1
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(3) == nullptr ) {
              (*it)->N(0)->AtBoundary( CNR1 );
              (*it)->N(1)->AtBoundary( BOTTOM );
              (*it)->N(3)->AtBoundary( LEFT );
              continue;
           } 
         // CNR2
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) == nullptr ) {
              (*it)->N(1)->AtBoundary( CNR2 );
              (*it)->N(0)->AtBoundary( BOTTOM );
              (*it)->N(2)->AtBoundary( RIGHT );
              continue;
           } 
         // CNR3
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(2) == nullptr ) {
              (*it)->N(2)->AtBoundary( CNR3 );
              (*it)->N(1)->AtBoundary( RIGHT );
              (*it)->N(3)->AtBoundary( TOP );
              continue;
           } 
         // CNR4
         if ( (*it)->Neighbor(2) == nullptr && (*it)->Neighbor(3) == nullptr ) {
              (*it)->N(3)->AtBoundary( CNR4 );
              (*it)->N(2)->AtBoundary( TOP );
              (*it)->N(0)->AtBoundary( LEFT );
              continue;
           } 
      }
      
    // testing nodes and elements
    printBoxBoundaryFlags( model ); 
       
 } // end recreateBoxBoundaryFlagsForQuadrilateralModel








/**
       Establishes the BOX_BOUNDARY flagging for a 3D model consisting out of hexahedral elements. 
       
       @note this method is quick.
*/ 
void recreateBoxBoundaryFlagsForHexahedralModel( Model<3U>& model )
 { 
    ErrorHandler&    csmp_error( ErrorHandler::Instance() );

    Region<3U>&      modeldomain(model.Region("Model"));
    vector<uint32_t> fnids;
    bool non_hex_cells_found{false};
    
    for ( auto it=modeldomain.CellsBegin(); it!=modeldomain.CellsEnd(); ++it )
      {
         // current version only works for linear hexahedra
         assert( (*it)->Nodes() <= 9 );

         if ( !isHexahedral( (*it)->FE_Type() ) ) {
              if ( !non_hex_cells_found ) {
                   csmp_error.Note( INFO, "recreateBoxBoundaryFlagsForHexahedralModel", "note that this method works only for hexahedra.");
                   non_hex_cells_found = true;
                }
              continue;
           }
           
         // DEFAULT (not at any boundary)
         for ( uint32_t i{0U}; i<(*it)->Nodes(); ++i )
           (*it)->N(i)->AtBoundary( NOT );
           
         // BOTTOM 
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              for ( const auto& j : (*it)->FE()->NodesOfFace(0U) )
                (*it)->N(j)->AtBoundary( BOTTOM );
              continue;
           }
         // RIGHT 
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(2) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              for ( const auto& j : (*it)->FE()->NodesOfFace(2U) )
                (*it)->N(j)->AtBoundary( RIGHT );
              continue;
           }
         // TOP 
         if ( (*it)->Neighbor(2) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              for ( const auto& j : (*it)->FE()->NodesOfFace(5U) )
                (*it)->N(j)->AtBoundary( TOP );
              continue;
           }
         // LEFT 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              for ( const auto& j : (*it)->FE()->NodesOfFace(4U) )
                (*it)->N(j)->AtBoundary( LEFT );
              continue;
           }
         // FRONT 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              for ( const auto& j : (*it)->FE()->NodesOfFace(1U) )
                (*it)->N(j)->AtBoundary( FRONT );
              continue;
           }
         // BACK 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              for ( const auto& j : (*it)->FE()->NodesOfFace(3U) )
                (*it)->N(j)->AtBoundary( BACK );
              continue;
           }
         // BACK CORNERS 
         // CNR1
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(3) == nullptr ) {
              (*it)->N(3)->AtBoundary( CNR1 );
              // sides
              (*it)->N(1)->AtBoundary( BOTTOM );
              (*it)->N(2)->AtBoundary( BACK );
              (*it)->N(4)->AtBoundary( LEFT );
              // edges
              (*it)->N(0)->AtBoundary( EDGE5 );
              (*it)->N(2)->AtBoundary( EDGE1 );
              (*it)->N(7)->AtBoundary( EDGE4 );
              continue;
           } 
         // CNR2
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) == nullptr ) {
              (*it)->N(2)->AtBoundary( CNR2 );
              // sides
              (*it)->N(0)->AtBoundary( BOTTOM );
              (*it)->N(5)->AtBoundary( RIGHT );
              (*it)->N(7)->AtBoundary( BACK );
              // edges
              (*it)->N(1)->AtBoundary( EDGE6 );
              (*it)->N(3)->AtBoundary( EDGE1 );
              (*it)->N(6)->AtBoundary( EDGE2 );
              continue;
           } 
         // CNR3
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(2) == nullptr ) {
              (*it)->N(6)->AtBoundary( CNR3 );
              // sides
              (*it)->N(1)->AtBoundary( RIGHT );
              (*it)->N(3)->AtBoundary( BACK );
              (*it)->N(4)->AtBoundary( TOP );
              // edges
              (*it)->N(2)->AtBoundary( EDGE2 );
              (*it)->N(5)->AtBoundary( EDGE7 );
              (*it)->N(7)->AtBoundary( EDGE3 );
              continue;
           } 
         // CNR4
         if ( (*it)->Neighbor(2) == nullptr && (*it)->Neighbor(3) == nullptr ) {
              (*it)->N(7)->AtBoundary( CNR4 );
              // sides
              (*it)->N(2)->AtBoundary( BACK );
              (*it)->N(0)->AtBoundary( LEFT );
              (*it)->N(5)->AtBoundary( TOP );
              // edges
              (*it)->N(3)->AtBoundary( EDGE4 );
              (*it)->N(6)->AtBoundary( EDGE3 );
              (*it)->N(4)->AtBoundary( EDGE8 );
              continue;
           } 
         // FRONT CORNERS 
         // CNR5
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(4) == nullptr ) {
              (*it)->N(0)->AtBoundary( CNR5 );
              // sides
              (*it)->N(2)->AtBoundary( BOTTOM );
              (*it)->N(5)->AtBoundary( FRONT );
              (*it)->N(7)->AtBoundary( LEFT );
              // edges
              (*it)->N(1)->AtBoundary( EDGE9 );
              (*it)->N(3)->AtBoundary( EDGE5 );
              (*it)->N(4)->AtBoundary( EDGE12 );
              continue;
           } 
         // CNR6
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(2) == nullptr ) {
              (*it)->N(1)->AtBoundary( CNR6 );
              // sides
              (*it)->N(3)->AtBoundary( BOTTOM );
              (*it)->N(6)->AtBoundary( RIGHT );
              (*it)->N(4)->AtBoundary( FRONT );
              // edges
              (*it)->N(0)->AtBoundary( EDGE9 );
              (*it)->N(2)->AtBoundary( EDGE6 );
              (*it)->N(5)->AtBoundary( EDGE10 );
              continue;
           } 
         // CNR7
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(2) == nullptr && (*it)->Neighbor(5) == nullptr ) {
              (*it)->N(5)->AtBoundary( CNR7 );
              // sides
              (*it)->N(1)->AtBoundary( RIGHT );
              (*it)->N(0)->AtBoundary( FRONT );
              (*it)->N(7)->AtBoundary( TOP );
              // edges
              (*it)->N(1)->AtBoundary( EDGE10 );
              (*it)->N(6)->AtBoundary( EDGE7 );
              (*it)->N(4)->AtBoundary( EDGE11 );
              continue;
           } 
         // CNR8
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(4) == nullptr && (*it)->Neighbor(5) == nullptr ) {
              (*it)->N(4)->AtBoundary( CNR8 );
              // sides
              (*it)->N(1)->AtBoundary( FRONT );
              (*it)->N(3)->AtBoundary( LEFT );
              (*it)->N(6)->AtBoundary( TOP );
              // edges
              (*it)->N(0)->AtBoundary( EDGE12 );
              (*it)->N(5)->AtBoundary( EDGE11 );
              (*it)->N(7)->AtBoundary( EDGE8 );
              continue;
           } 
      }
      
    // testing nodes and elements
    printBoxBoundaryFlags( model ); 
       
 } // end recreateBoxBoundaryFlagsForHexahedralModel




/// helper function for atBoundary() below (substitute for C++20 contains)
static bool isIn( const set<BOX_BOUNDARY>& bflags, initializer_list<BOX_BOUNDARY> flags ) {
     for ( auto it : flags )
       if ( bflags.find(it) != bflags.end() ) return true;
     return false;
  }






/**
    Returns whether the cell is at the model boundary
    
    @atttention For equidimensional elements, this assessment is  based on whether the cell is missing any of its neighbor elements
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool atBoundary( const CELL<dim>* const cptr )
 {
    // only in debug mode
    assert(cptr != nullptr);
    
    // Enforce noexcept on key CELL methods at compile time
    static_assert(noexcept(std::declval<const CELL<dim>&>().Nodes()),
                  "CELL::Nodes() must be noexcept");
    static_assert(noexcept(std::declval<const CELL<dim>&>().N(0)),
                  "CELL::N(size_t) must be noexcept");
    static_assert(noexcept(std::declval<const CELL<dim>&>().Neighbors()),
                  "CELL::Neighbors() must be noexcept");
    static_assert(noexcept(std::declval<const CELL<dim>&>().Neighbor(0)),
                  "CELL::Neighbor(size_t) must be noexcept");
    static_assert(noexcept(std::declval<const CELL<dim>&>().IsEquidimensional()),
                  "CELL::IsEquidimensional() must be noexcept");

    const uint32_t n_nbors = cptr->Neighbors();

    // Case 1: equidimensional elements (volumes in 3D, surfaces in 2D, lines in 1D)
    if ( cptr->IsEquidimensional() ) {
        for ( uint32_t i = 0u; i < n_nbors; ++i )
          if ( !cptr->Neighbor(i) )
            return true; // at least one node is interior -> not a boundary element
      }

    // Case 2: lower-dimensional elements (surfaces in 3D, lines in 2D, points in 1D)
    if constexpr( dim == 1U ) {
         vector<uint32_t> cnr_nodes;
         cptr->FE()->CornerNodes( cnr_nodes );
         if ( cptr->N(cnr_nodes[0])->AtBoundary() == NOT ) return true;
         if ( cptr->N(cnr_nodes[1])->AtBoundary() == NOT ) return true;
      }
      // 2D or 3D
    else {
        for ( uint32_t i = 0u; i < n_nbors; ++i ) {
            BOX_BOUNDARY b = atBoundary( cptr, i );
            if (  b != NOT )
              return true;
          }
      }
      
    // otherwise: all faces have neighbors -> interior
    return false;
       
 } // end atBoundary(bool)

template bool atBoundary( const Element<1U>* const );
template bool atBoundary( const Element<2U>* const );
template bool atBoundary( const Element<3U>* const );




/**
    Infers from the node flags and the boundary face (local) number, which boundary the cell lies on;
    returns INTERNAL if all flags are not but the cell has no neighbor
    
    @attention only cells which have a face on the BOX_BOUNDARY are considered boundary elements
    
    @attention lower dimensional cells may be flagged as INTERNAL boundary if they have no neighbor
    but occur inside of the model.
    
      @author SKM
      @date 14/10/21
      @test OK
*/
template<uint32_t dim, template<uint32_t> class CELL>
BOX_BOUNDARY atBoundary( const CELL<dim>* const eptr, uint32_t b_face )
{
    assert(eptr != nullptr);
    assert(b_face < eptr->Faces() && eptr->Faces() == eptr->Neighbors());

    // if the lower dimensional element or Face has an equidimensional neighbour element
    if (eptr->Neighbor(b_face) != nullptr) return NOT;

    // collecting boundary flags from the face nodes
    set<BOX_BOUNDARY> fflags;
    uint32_t n_face_nodes{0u};
    for (const auto* node_ptr : eptr->CornerNodesOfFace(b_face)) {
        //                            ^^^^^^^^^^^^^^^^^^^^^^^^^^
        const BOX_BOUNDARY b = node_ptr->AtBoundary();
        if (b != NOT) fflags.insert(b);
        n_face_nodes++;
    }

    // if the element is part of- and at the edge of a lower-dimensional region inside of a model
    // all edge nodes should be flagged INTERNAL (but not the geometry flags)
    if ( fflags.empty() ) return INTERNAL;
    
    // if there are nodes with boundary flags, further diagnostics are applied
    return atBoundary<dim>( fflags, n_face_nodes );
}
template BOX_BOUNDARY atBoundary( const Element<1U>* const, uint32_t );
template BOX_BOUNDARY atBoundary( const Element<2U>* const, uint32_t );
template BOX_BOUNDARY atBoundary( const Element<3U>* const, uint32_t );




/**
    The atBoundary function is designed to determine the boundary type (side, edge, or corner) of a boundary face
    based on the flags of its corner nodes and the number of corner nodes (boundary_face_corner_nodes).
    
    It thus aims at identifying the boundary that a lower dimensional element is located on
    from the BOX_BOUNDARY flags assigned to its nodes.
 */
template<uint32_t dim>
BOX_BOUNDARY atBoundary( const set<BOX_BOUNDARY>& flags, uint32_t boundary_face_corner_nodes )
{
    assert( flags.find(NOT) == flags.end() );

    ErrorHandler& csmp_error(ErrorHandler::Instance());

    if (flags.empty()) {
        csmp_error.Note(WARNING, "atBoundary(flag_set,n_cnr_odes)", "face contains no boundary-flagged nodes");
        return NOT;
    }
    if (flags.size() > 4 ) {
        csmp_error.Note(ERROR, "atBoundary(flag_set,n_cnr_odes)", "function does not yet consider cells with mid-side or bubble nodes");
        return NOT;
    }

    // if only a single node on the face has a boundary flag, this one is returned
    // (case of a line element)
    if (flags.size() == 1U)  return *flags.begin();

    // rectangular models
    // ==================
    // tested: OK
    if constexpr (dim == 2U) {
        if (boundary_face_corner_nodes == 2U) {
            // if a line-element face forms part of an internal boundary that classification takes precedence
            if (flags.count(INTERNAL))  return INTERNAL;
            // other cases
            auto it = flags.begin();
            BOX_BOUNDARY a = *it++;
            BOX_BOUNDARY b = *it;
            if (isSide(a) && !isSide(b)) return a;
            if (isSide(b) && !isSide(a)) return b;
            if (a == BOTTOM && (b == CNR1 || b == CNR2 || b == CNR3 || b == CNR4)) return BOTTOM;
            if (b == BOTTOM && (a == CNR1 || a == CNR2 || a == CNR3 || a == CNR4)) return BOTTOM;
            if (a == TOP && (b == CNR1 || b == CNR2 || b == CNR3 || b == CNR4)) return TOP;
            if (b == TOP && (a == CNR1 || a == CNR2 || a == CNR3 || a == CNR4)) return TOP;
            if (a == LEFT && (b == CNR1 || b == CNR2 || b == CNR3 || b == CNR4)) return LEFT;
            if (b == LEFT && (a == CNR1 || a == CNR2 || a == CNR3 || a == CNR4)) return LEFT;
            if (a == RIGHT && (b == CNR1 || b == CNR2 || b == CNR3 || b == CNR4)) return RIGHT;
            if (b == RIGHT && (a == CNR1 || a == CNR2 || a == CNR3 || a == CNR4)) return RIGHT;
            if ((a == CNR3 && b == CNR4) || (a == CNR4 && b == CNR3)) return BOTTOM;
            if ((a == CNR1 && b == CNR2) || (a == CNR2 && b == CNR1)) return BOTTOM;
            if ((a == CNR2 && b == CNR4) || (a == CNR4 && b == CNR2)) return RIGHT;
            if ((a == CNR1 && b == CNR3) || (a == CNR3 && b == CNR1)) return RIGHT;
            if ((a == CNR4 && b == CNR1) || (a == CNR1 && b == CNR4)) return LEFT;
            if ((a == CNR3 && b == CNR2) || (a == CNR2 && b == CNR3)) return TOP;
            if ((a == LEFT_OUTSIDE && b == BOTTOM_OUTSIDE) || (a == BOTTOM_OUTSIDE && b == LEFT_OUTSIDE)) return CNR1;
            if ((a == BOTTOM_OUTSIDE && b == RIGHT_OUTSIDE) || (a == RIGHT_OUTSIDE && b == BOTTOM_OUTSIDE)) return CNR2;
            if ((a == RIGHT_OUTSIDE && b == TOP_OUTSIDE) || (a == TOP_OUTSIDE && b == RIGHT_OUTSIDE)) return CNR3;
            if ((a == TOP_OUTSIDE && b == LEFT_OUTSIDE) || (a == LEFT_OUTSIDE && b == TOP_OUTSIDE)) return CNR4;
        }
        // model edges
        if (auto e = uniqueEdgeOrNot(flags); isEdge(e)) return e;
        if (auto f = uniqueImpliedSideOrNot(flags); isSide(f)) return f;
        csmp_error.Note(ERROR, "atBoundary(flag_set,n_cnr_odes,2D)", "could not deduce unique boundary from node flags.");
        return IRREGULAR;
    }

    // box-shaped models
    // =================
    if constexpr (dim == 3U) {
        // if the face is a lower-dimensional (surface) element
        if (boundary_face_corner_nodes == 2U) {
            // all flags are the same
            if ( flags.size() == 1U ) return (*flags.begin());
            else if ( flags.size() == 2U ) {
                 // sides and edges
                 return whichBoundary( (*flags.begin()), (*flags.rbegin()) );
              }
        }
        // if the face is equidimensional (a triangle or quadrilateral)
        if ( boundary_face_corner_nodes > 2U ) {
            // Prioritize side flags
            for ( const auto& bit : flags ) if (isSide(bit)) return bit;
            // When there are only corners
            unsigned int n_corners{0u};
            for ( const auto& bit : flags ) if ( isCorner(bit) ) n_corners++;
            if ( n_corners == flags.size() ) {
              const auto end = flags.end();
              // four corner-node sides
              if ( n_corners == 4 ) {
                   if ( flags.find(CNR1) != end && flags.find(CNR2) != end &&
                        flags.find(CNR5) != end && flags.find(CNR6) != end ) return BOTTOM;
                        
                   if ( flags.find(CNR5) != end && flags.find(CNR6) != end &&
                        flags.find(CNR7) != end && flags.find(CNR8) != end ) return FRONT;
                        
                   if ( flags.find(CNR2) != end && flags.find(CNR3) != end &&
                        flags.find(CNR7) != end && flags.find(CNR6) != end ) return RIGHT;
                        
                   if ( flags.find(CNR1) != end && flags.find(CNR2) != end &&
                        flags.find(CNR3) != end && flags.find(CNR4) != end ) return BACK;
                        
                   if ( flags.find(CNR1) != end && flags.find(CNR5) != end &&
                        flags.find(CNR8) != end && flags.find(CNR4) != end ) return LEFT;
                        
                   if ( flags.find(CNR3) != end && flags.find(CNR4) != end &&
                        flags.find(CNR7) != end && flags.find(CNR8) != end ) return TOP;
                }
              // 3-corner-node sides
              else if ( n_corners == 3 ) {
                   if ( (flags.find(CNR1) != end && flags.find(CNR5) != end && flags.find(CNR6) != end) ||
                        (flags.find(CNR1) != end && flags.find(CNR6) != end && flags.find(CNR2) != end) ||
                        (flags.find(CNR2) != end && flags.find(CNR5) != end && flags.find(CNR6) != end) ||
                        (flags.find(CNR1) != end && flags.find(CNR2) != end && flags.find(CNR5) != end) ) return BOTTOM;

                   if ( (flags.find(CNR5) != end && flags.find(CNR6) != end && flags.find(CNR8) != end) ||
                        (flags.find(CNR6) != end && flags.find(CNR7) != end && flags.find(CNR8) != end) ||
                        (flags.find(CNR5) != end && flags.find(CNR6) != end && flags.find(CNR7) != end) ||
                        (flags.find(CNR5) != end && flags.find(CNR7) != end && flags.find(CNR8) != end) ) return FRONT;

                   if ( (flags.find(CNR2) != end && flags.find(CNR7) != end && flags.find(CNR6) != end) ||
                        (flags.find(CNR2) != end && flags.find(CNR3) != end && flags.find(CNR7) != end) ||
                        (flags.find(CNR2) != end && flags.find(CNR3) != end && flags.find(CNR6) != end) ||
                        (flags.find(CNR3) != end && flags.find(CNR7) != end && flags.find(CNR6) != end) ) return RIGHT;

                   if ( (flags.find(CNR1) != end && flags.find(CNR2) != end && flags.find(CNR4) != end) ||
                        (flags.find(CNR2) != end && flags.find(CNR3) != end && flags.find(CNR4) != end) ||
                        (flags.find(CNR1) != end && flags.find(CNR2) != end && flags.find(CNR3) != end) ||
                        (flags.find(CNR1) != end && flags.find(CNR3) != end && flags.find(CNR4) != end) ) return BACK;

                   if ( (flags.find(CNR1) != end && flags.find(CNR5) != end && flags.find(CNR4) != end) ||
                        (flags.find(CNR4) != end && flags.find(CNR5) != end && flags.find(CNR8) != end) ||
                        (flags.find(CNR1) != end && flags.find(CNR5) != end && flags.find(CNR8) != end) ||
                        (flags.find(CNR1) != end && flags.find(CNR8) != end && flags.find(CNR4) != end) ) return LEFT;

                   if ( (flags.find(CNR3) != end && flags.find(CNR8) != end && flags.find(CNR7) != end) ||
                        (flags.find(CNR3) != end && flags.find(CNR4) != end && flags.find(CNR8) != end) ||
                        (flags.find(CNR3) != end && flags.find(CNR4) != end && flags.find(CNR7) != end) ||
                        (flags.find(CNR4) != end && flags.find(CNR8) != end && flags.find(CNR7) != end) ) return TOP;
                }
            }

            // oblique cuts through single-cell model or other
            return IRREGULAR;
        }
        csmp_error.Note(ERROR, "atBoundary(flag_set,n_cnr_odes,3D)", "invalid number of corner nodes: " + std::to_string(boundary_face_corner_nodes));
        return IRREGULAR;
    }

    return INTERNAL;
}

template BOX_BOUNDARY atBoundary<2>( const set<BOX_BOUNDARY>&, uint32_t );
template BOX_BOUNDARY atBoundary<3>( const set<BOX_BOUNDARY>&, uint32_t );




/**
returns true if the model contains the complete set of box boundary identifiers
all nodes are considered.
*/
bool isStrictlyBoxShaped( const Model<2U>& model )
{
  set<BOX_BOUNDARY> flags2d{ LEFT, RIGHT, TOP, BOTTOM, CNR1, CNR2, CNR3, CNR4, EDGE1, EDGE2, EDGE3, EDGE4 },
                    flags_of_model;

  const Region<2>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    if ( (*it)->AtBoundary() != NOT )
      flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags2d in flags_of_model
  auto it = flags_of_model.begin();
  for ( auto i = flags2d.begin(); i != flags2d.end() && it != flags_of_model.end(); ++i )
  {
    it = std::lower_bound( it, flags_of_model.end(), (*i) );
    // make sure the found item is a match
    if ( it != flags_of_model.end() && *i < *it )
      it = flags_of_model.end(); // break out early
  }
  if ( it != flags_of_model.end() ) return true;

  return false;

} // end isStrictlyBoxShaped






/**
tests (only) if all boundary identifiers are there
*/
bool isStrictlyBoxShaped( const Model<3U>& model )
{
  const set<BOX_BOUNDARY> flags3d{ LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK, CNR1, CNR2, CNR3, CNR4, CNR5, CNR6, CNR7, CNR8,
                                   EDGE1, EDGE2, EDGE3, EDGE4, EDGE5, EDGE6, EDGE7, EDGE8, EDGE9, EDGE10, EDGE11, EDGE12 };

  set<BOX_BOUNDARY>  flags_of_model;

  const Region<3>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags2d in flags_of_model
  size_t flags_not_found{0U};
  bool   first_error{true};
  for ( const auto& i : flags3d )
    if ( flags_of_model.count(i) == 0 ) {
         if ( first_error ) {
              cerr <<"\n"<<"isStrictlyBoxShaped: NO! - model misses nodes flagged";
              first_error = false;
           }
         cerr <<" "<< parseBoundary(i);
         flags_not_found++;
      }
  if ( !first_error ) cout << endl;

  if ( flags_not_found == 0U ) return true;
  return false;
  
} // end isStrictlyBoxShaped




/**
returns true if the model contains all the side boundary identifiers of the box
*/
bool hasAllSideBoundaries( const Model<3U>& model )
{
  set<BOX_BOUNDARY> flags3d{ LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK}, flags_of_model;

  const Region<3>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    if ( (*it)->AtBoundary() != NOT )
      flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags3d in flags_of_model
  for ( const auto& i : flags3d )
    if ( flags_of_model.count(i) < 1 ) return false;

  return true;
}




/// 2D side boudaries only
bool hasAllSideBoundaries( const Model<2U>& model )
{
  set<BOX_BOUNDARY> flags2d{ LEFT, RIGHT, TOP, BOTTOM }, flags_of_model;

  const Region<2>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    if ( (*it)->AtBoundary() != NOT )
      flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags2d in flags_of_model
  auto it = flags_of_model.begin();
  for ( auto i = flags2d.begin(); i != flags2d.end() && it != flags_of_model.end(); ++i )
  {
    it = std::lower_bound( it, flags_of_model.end(), (*i) );
    // make sure the found item is a match
    if ( it != flags_of_model.end() && *i < *it )
      it = flags_of_model.end(); // break out early
  }
  if ( it != flags_of_model.end() ) return true;

  return false;
}


/**
    Encodes boundary flag values into scalar variables which can be examined for correctness.

    @attention for element flags, this only works for elements that are not simultaneously sharing different boundaries.
    The latter just map to one of the multiple possible boundaries.

    @author SKM 7/3/16
*/
template<uint32_t dim>
void boxFlagsToVariable( Model<dim>& model, const char* node_variable, const char* elmt_variable )
{
  const csmp::Index nprop_key( model.Database().StorageKey( node_variable ) );
  assert( nprop_key.place == NODE );
  assert( nprop_key.type == SCALAR );
  const csmp::Index eprop_key( model.Database().StorageKey( elmt_variable ) );
  assert( eprop_key.place == ELEMENT );
  assert( eprop_key.type == SCALAR );
  Region<dim>&  mregion( model.Region( "Model" ) );

  // node flag conversion
  for ( auto nit = mregion.NodesBegin(); nit != mregion.NodesEnd(); nit++ )
    (*nit)->Store( nprop_key, makeScalar( ANY, static_cast<double>((*nit)->AtBoundary()) ) );

   // element flag to 'emt_variable' conversion
  for ( auto eit = mregion.CellsBegin(); eit != mregion.CellsEnd(); eit++ ) {
       uint32_t face = numeric_limits<uint32_t>::max();
       bool at_boundary{false};
       for ( uint32_t i{0U}; i<(*eit)->Neighbors(); ++i )
         if ( (*eit)->Neighbor(i) == nullptr ) {
              face        = i;
              at_boundary = true;
              break;
           }
       if ( at_boundary )
         (*eit)->Store( eprop_key, makeScalar( ANY, static_cast<double>(atBoundary(*eit,face)) ) );
       else //                                                          ^^^^^^^^^^^^^^^^^^^^^
         (*eit)->Store( eprop_key, makeScalar( ANY, static_cast<double>(NOT) ) );
    }

} // end boxFlagsToVariable

template void boxFlagsToVariable( Model<2>&, const char*, const char* );
template void boxFlagsToVariable( Model<3>&, const char*, const char* );





template<uint32_t dim>
void boxBoundaryPropertyRange( const Model<dim>& sg, BOX_BOUNDARY boundary,
                               const char* node_property, double& bmin, double& bmax )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const csmp::Index  prop_key = sg.Database().StorageKey( node_property );
  
  if ( prop_key.place != NODE || prop_key.type != SCALAR )
    csmp_error.Note( ERROR, "boxBoundaryPropertyRange",
                      "Property must a scalar placed on the nodes; nothing was done.");

  bool  first_value( true ), verbose( false );

  const csmp::Region<dim>& sgref = sg.Region( "Model" );

  // measuring the property ranges
  for ( auto nit = sgref.PerimeterNodesBegin(); nit != sgref.NodesEnd(); nit++ )
    if ( (*nit)->AtBoundary() == boundary )
      {
        if ( first_value ) {
          bmin = bmax = (*nit)->Read( prop_key );
          first_value = false;
        }
        else {
          bmin = std::min( bmin, (*nit)->Read( prop_key ) );
          bmax = std::max( bmax, (*nit)->Read( prop_key ) );
        }
        if ( verbose )
          cout << "\nboundary node " << (*nit)->Idx() << ": " << (*nit)->x() << " " << (*nit)->y() << " " << (*nit)->z();
      }

  if ( verbose ) {
      cout << "\nFiniteVolumeTransport::BoundaryPropertyRanges: of physical variable '" << node_property << "':";
      cout << "\nrange of '" << node_property << "' at boundary:  " << bmin << " to " << bmax << endl;
    }

} // end BoundaryPropertyRanges




BOX_BOUNDARY  whichEdge( BOX_BOUNDARY side1, BOX_BOUNDARY side2 )
{
  // dealing with the case where the 2 flags are the same so that this is not an EDGE
  if (side1 == side2) return side1;

  if (side1 > side2) swap( side1, side2 );

  if (side1 == BACK and side2 == BOTTOM) return EDGE1;
  if (side1 == BACK and side2 == RIGHT)  return EDGE2;
  if (side1 == BACK and side2 == TOP)    return EDGE3;
  if (side1 == BACK and side2 == LEFT)  return EDGE4;

  if (side1 == FRONT and side2 == BOTTOM) return EDGE9;
  if (side1 == FRONT and side2 == RIGHT) return EDGE10;
  if (side1 == FRONT and side2 == TOP) return EDGE11;
  if (side1 == FRONT and side2 == LEFT) return EDGE12;

  if (side1 == TOP and side2 == RIGHT) return EDGE7;
  if (side1 == TOP and side2 == LEFT) return EDGE8;

  if (side1 == BOTTOM and side2 == RIGHT) return EDGE6;
  if (side1 == BOTTOM and side2 == LEFT) return EDGE5;

  return NOT;
  
} // end whichEdge


BOX_BOUNDARY  whichCorner( BOX_BOUNDARY side1, BOX_BOUNDARY side2, BOX_BOUNDARY side3 )
{
  set<BOX_BOUNDARY> sides({ side1,side2,side3 });
  auto s1 = (*sides.begin());
  auto s2 = (*(next(sides.begin(), 1)));
  auto s3 = (*sides.rbegin());

  if (sides.empty()) return NOT;

  if (sides.size() == 1U) return (*sides.begin());

  // dealing with duplicates (there will only be 1 or 2 entries so this can only be a corner if the entries are edges)
  if (sides.size() == 2U) {
    // 2 edges EDGE12, EDGE11... EDGE1
    if (s1 == EDGE12 and s2 == EDGE11) return CNR8;
    if (s1 == EDGE12 and s2 == EDGE9) return CNR5;

    if (s1 == EDGE11 and s2 == EDGE10) return CNR7;

    if (s1 == EDGE10 and s2 == EDGE9) return CNR6;
    if (s1 == EDGE10 and s2 == EDGE6) return CNR6;

    if (s1 == EDGE9 and s2 == EDGE6) return CNR6;
    if (s1 == EDGE9 and s2 == EDGE5) return CNR5;

    if (s1 == EDGE8 and s2 == EDGE4) return CNR4;
    if (s1 == EDGE8 and s2 == EDGE3) return CNR4;

    if (s1 == EDGE7 and s2 == EDGE3) return CNR3;
    if (s1 == EDGE7 and s2 == EDGE2) return CNR3;

    if (s1 == EDGE6 and s2 == EDGE2) return CNR2;
    if (s1 == EDGE6 and s2 == EDGE1) return CNR2;

    if (s1 == EDGE5 and s2 == EDGE4) return CNR1;
    if (s1 == EDGE5 and s2 == EDGE1) return CNR1;

    if (s1 == EDGE4 and s2 == EDGE3) return CNR4;
    if (s1 == EDGE4 and s2 == EDGE1) return CNR1;

    if (s1 == EDGE3 and s2 == EDGE2) return CNR3;

    if (s1 == EDGE2 and s2 == EDGE1) return CNR2;

    // not a corner
    return whichEdge(side1, side2);
  }

  // obeying increasing value constraint: BACK, FRONT, TOP, BOTTOM, RIGHT, LEFT
  if (s1 == BACK and s2 == BOTTOM and s3 == LEFT)  return CNR1;
  if (s1 == BACK and s2 == BOTTOM and s3 == RIGHT) return CNR2;
  if (s1 == BACK and s2 == TOP    and s3 == RIGHT) return CNR3;
  if (s1 == BACK and s2 == TOP    and s3 == LEFT)  return CNR4;

  if (s1 == FRONT and s2 == BOTTOM and s3 == LEFT)  return CNR5;
  if (s1 == FRONT and s2 == BOTTOM and s3 == RIGHT) return CNR6;
  if (s1 == FRONT and s2 == TOP    and s3 == RIGHT) return CNR7;
  if (s1 == FRONT and s2 == TOP    and s3 == LEFT)  return CNR8;

  return NOT;
} // end whichCorner



/**
   Returns which boundary the edge or lower-dim face with  two end-node flags is on
   
   TODO: Does this method need to distinguish between 2D and 3D models?
*/
BOX_BOUNDARY  whichBoundary( BOX_BOUNDARY node_flag1, BOX_BOUNDARY node_flag2 )
 {
    if ( node_flag1 == NOT && node_flag2 == NOT ) return INTERNAL;
    if ( node_flag1 == node_flag2 ) return node_flag1;
    
    // sides
    if ( isLEFT(node_flag1)   && isLEFT(node_flag2) )   return LEFT;
    if ( isRIGHT(node_flag1)  && isRIGHT(node_flag2) )  return RIGHT;
    if ( isBOTTOM(node_flag1) && isBOTTOM(node_flag2) ) return BOTTOM;
    if ( isTOP(node_flag1)    && isTOP(node_flag2) )    return TOP;
    if ( isFRONT(node_flag1)  && isFRONT(node_flag2) )  return FRONT;
    if ( isBACK(node_flag1)   && isBACK(node_flag2) )   return BACK;
    
    // edges (using the fact that edges identifiers are all more negative than corners!)
    // back
    if ( node_flag1 == EDGE1 || node_flag2 == CNR1 ) return EDGE1;
    if ( node_flag1 == EDGE1 || node_flag2 == CNR2 ) return EDGE1;
    
    if ( node_flag1 == EDGE2 || node_flag2 == CNR2 ) return EDGE2;
    if ( node_flag1 == EDGE2 || node_flag2 == CNR3 ) return EDGE2;

    if ( node_flag1 == EDGE3 || node_flag2 == CNR3 ) return EDGE3;
    if ( node_flag1 == EDGE3 || node_flag2 == CNR4 ) return EDGE3;

    if ( node_flag1 == EDGE4 || node_flag2 == EDGE4 ) return EDGE4;
    if ( node_flag1 == EDGE4 || node_flag2 == EDGE1 ) return EDGE4;
    // neither back nor front
    if ( node_flag1 == EDGE5 || node_flag2 == CNR1 ) return EDGE5;
    if ( node_flag1 == EDGE5 || node_flag2 == CNR5 ) return EDGE5;

    if ( node_flag1 == EDGE6 || node_flag2 == CNR2 ) return EDGE6;
    if ( node_flag1 == EDGE6 || node_flag2 == CNR6 ) return EDGE6;

    if ( node_flag1 == EDGE7 || node_flag2 == CNR3 ) return EDGE7;
    if ( node_flag1 == EDGE7 || node_flag2 == CNR7 ) return EDGE7;

    if ( node_flag1 == EDGE8 || node_flag2 == CNR4 ) return EDGE8;
    if ( node_flag1 == EDGE8 || node_flag2 == CNR8 ) return EDGE8;
    // front
    if ( node_flag1 == EDGE9 || node_flag2 == CNR5 ) return EDGE9;
    if ( node_flag1 == EDGE9 || node_flag2 == CNR6 ) return EDGE9;

    if ( node_flag1 == EDGE10 || node_flag2 == CNR6 ) return EDGE10;
    if ( node_flag1 == EDGE10 || node_flag2 == CNR7 ) return EDGE10;

    if ( node_flag1 == EDGE11 || node_flag2 == CNR7 ) return EDGE11;
    if ( node_flag1 == EDGE11 || node_flag2 == CNR8 ) return EDGE11;

    if ( node_flag1 == EDGE12 || node_flag2 == CNR8 ) return EDGE12;
    if ( node_flag1 == EDGE12 || node_flag2 == CNR5 ) return EDGE12;
    
    // should not go here
    return MULTIPLE;
 
 } // end whichBoundary





double bilinearInterpolate( uint32_t idx_x, uint32_t,
                            const Point<1U>& xy1,
                            const Point<1U>& xy2,
                            const Point<1U>& coord,
                            double p1, double p2, double, double )
{
      // If min-coords. are equivalent to max-coords. the boundary-value average
      // is assigned.
      if ( xy1 == xy2 )
        {
           cout <<"\nbilinearInterpolate: min/max coordinates are identical:"<< endl;
           xy1.Out();
           xy2.Out();
           return (p1+p2) / 2.0;
        }

      if ( p1 == p2 ) return p1;
    
      // 4. computing interpolation function. Num. Recip. p. 105
      double t = (coord[idx_x] - xy1[idx_x] ) / (xy2[idx_x] - xy1[idx_x] );
      
      // 5 bi-linear interpolation
      return (1. - t) * p1 + t * p2;

} // end bilinearInterpolate





/**

interpolateXY() uses bilinear interpolation to find the value of a
property specified at the corner points of a rectangle, at the coordinates
of a point located inside of this rectangle.

Since interpolateXY() carries out an interpolation on a planar surface
in 3D space, it also needs the integer indices which give the axis
in the reference coordinate system.

@section arguments Input Arguments

The first two arguments of interpolateXY() specify the coordinate axis
indices of the following Point<dim> arguments which shall be used in
the interpolation. For instance, for i=0, j=1, the interpolation will be
carried out in the XY plane.

The three following VectorVariable<2U> arguments specify the lower left and
upper right right corners of the rectangle in which the variable value shall
be interpolated at a point given by the third VectorVariable<2U> argument.

The last four floating point arguments define the values of the variable
which is to be interpolated. They are the corner points of the rectangle
listed in counter-clockwise fashion (e.g., lower left, lower right, upper
right and upper left corners, respectively).

@return The result of the interpolation is returned into a double type variable.

@section application Application

Function is used by AssignBoundaryValues().

@section messages Messages

The function will report an error and return the average value of the
four cornerpoints if their coordinates are identical.
*/
double bilinearInterpolate( uint32_t idx_x, uint32_t idx_y,
                            const Point<2U>& xy1,
                            const Point<2U>& xy2,
                            const Point<2U>& coord,
                            double p1, double p2, double p3, double p4 )
{
      // If min-coords. are equivalent to max-coords. the boundary-value average
      // is assigned.
      if ( xy1 == xy2 ) {
           cout <<"\nbilinearInterpolate: min/max coordinates are identical:"<< endl;
           xy1.Out();
           xy2.Out();
           return (p1+p2+p3+p4) / 4.0;
        }

      if ( p1 == p2 && p2 == p3 && p3 == p4 ) return p1;
    
      // 4. computing interpolation functions. Num. Recip. p. 105
      double t = (coord[idx_x]-xy1[idx_x] ) / (xy2[idx_x] - xy1[idx_x] );
      double u = (coord[idx_y]-xy1[idx_y] ) / (xy2[idx_y] - xy1[idx_y] );
      
      // 5 bi-linear interpolation
      return (1.-t) * (1.-u) * p1 + t * (1.-u) * p2 + t * u * p3 + (1.-t) * u * p4;

} // end bilinearInterpolate




/**

interpolateXY() uses bilinear interpolation to find the value of a
property specified at the corner points of a rectangle, at the coordinates
of a point located inside of this rectangle.

Since interpolateXY() carries out an interpolation on a planar surface
in 3D space, it also needs the integer indices which give the axis
in the reference coordinate system.

@section arguments Input Arguments

The first two arguments of interpolateXY() specify the coordinate axis
indices of the following VectorVariable<dim> arguments which shall be used in
the interpolation. For instance, for i=0, j=1, the interpolation will be
carried out in the XY plane.

The three following VectorVariable<dim> arguments specify the lower left and
upper right right corners of the rectangle in which the variable value shall
be interpolated at a point given by the third VectorVariable<dim> argument.

The last four floating point arguments define the values of the variable
which is to be interpolated. They are the corner points of the rectangle
listed in counter-clockwise fashion (e.g., lower left, lower right, upper
right and upper left corners, respectively).

@return The result of the interpolation is returned into a double type variable.

@section application Application

Function is used by AssignBoundaryValues().

@section messages Messages

The function will report an error and return the average value of the
four cornerpoints if their coordinates are identical.
*/
double bilinearInterpolate( uint32_t idx_x, uint32_t idx_y,
                            const Point<3U>& xy1,
                            const Point<3U>& xy2,
                            const Point<3U>& coord,
                            double p1, double p2, double p3, double p4 )
{
      // If min-coords. are equivalent to max-coords. the boundary-value average
      // is assigned.
      if ( xy1 == xy2 ) {
           cout <<"\nbilinearInterpolate: min/max coordinates are identical:"<< endl;
           xy1.Out();
           xy2.Out();
           return (p1+p2+p3+p4) / 4.0;
        }

      if ( p1 == p2 && p2 == p3 && p3 == p4 ) return p1;
    
      // 4. computing interpolation functions. Num. Recip. p. 105
      double t = (coord[idx_x]-xy1[idx_x] ) / (xy2[idx_x] - xy1[idx_x] );
      double u = (coord[idx_y]-xy1[idx_y] ) / (xy2[idx_y] - xy1[idx_y] );
      
      // 5 bi-linear interpolation
      return (1.-t) * (1.-u) * p1 + t * (1.-u) * p2 + t * u * p3 + (1.-t) * u * p4;

} // end bilinearInterpolate


/// interpolate along boundaries of 2D rectangle-shaped model
double linearInterpolate( const pair<Point<1U>,double>& p1, // endpoint1, value1
                          const pair<Point<1U>,double>& p2, // endpoint2, value2
                          const Point<1U>& pt )                // current point x,y,z
{
    // endmember value range
    double dval = p2.second - p1.second;
    
    // distance between endpoints
    double dx = p2.first[0] - p1.first[0];
    
    // distance between current point and point 1
    double dist = pt[0] - p1.first[0];

    // computing the interpolated value (y = b + mx)
    //     min     normalized distance    gradient
    return p1.second + (dist * dval) / dx;

} // end



/// interpolate along boundaries of 2D rectangle-shaped model
double linearInterpolate( const pair<Point<2U>,double>& p1, // endpoint1, value1
                          const pair<Point<2U>,double>& p2, // endpoint2, value2
                          const Point<2U>& pt )             // current point x,y,z
{
    // endmember value range
    double dval = p2.second - p1.second;
    
    // distance between endpoints
    double dx   = distance( p1.first, p2.first );
    
    // distance between current point and point 1
    double dist = distance( p1.first, pt );

    // computing the interpolated value (y = b + mx)
    //     min     normalized distance    gradient
    return p1.second + (dist * dval) / dx;

} // end


/// interpolate along boundaries of 2D rectangle-shaped model
double linearInterpolate( const pair<Point<3U>,double>& p1, // endpoint1, value1
                          const pair<Point<3U>,double>& p2, // endpoint2, value2
                          const Point<3U>& pt )             // current point x,y,z
{
    // endmember value range
    double dval = p2.second - p1.second;
    
    // distance between endpoints
    double dx   = distance( p1.first, p2.first );
    
    // distance between current point and point 1
    double dist = distance( p1.first, pt );

    // computing the interpolated value (y = b + mx)
    //     min     normalized distance    gradient
    return p1.second + (dist * dval) / dx;

} // end






/**
     Finds edges in boxed shaped model and gibes them a BOX_BOUDARY_FLAG
     dependent on their location.
*/
template<uint32_t dim>
void flagEdges( VSet<dim>& vset, double tol,
                double xmin, double xmax,
                double ymin, double ymax,
                double zmin, double zmax )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     // finding min-max of x, and y coordinates of boundary nodes
     map<size_t,double>  bflags;

     // 3D CASE finding the edges on the basis of their coordinates
     // and flagging them accordingly. note there is no 2D case since in 2D we
     // are splitting up triangles that belong to two boundaries so that
     // they belong to only one side
     size_t n_node(0U);
     for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
       {
          // no flagging at all
          if ( (*bit) == 0 || (*bit) == IRREGULAR )
            {
               csmp_error.Note( ERROR, "flagEdges",
                               "Skipping node, the boundary flag of which could not be identified" );
               cout <<"\nNode "<< n_node <<", flagged: "<< (*bit) << endl;
            }
          
          // Edges along the back side (Edges 1 to 4)
          if ( approximatelyEqual(vset.Pz( n_node ),zmin,tol) )
             {
                // EDGE1
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = BACK_BOTTOM;
                // EDGE2
                if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = BACK_RIGHT;
                // EDGE3
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = BACK_TOP;
                // EDGE4
                if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = BACK_LEFT;
             }
             
          // Edges along the front side (Edges 9 to 12)
          if ( approximatelyEqual(vset.Pz( n_node ),zmax) )
             {
                // EDGE9
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = FRONT_BOTTOM;
                // EDGE10
                if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = FRONT_RIGHT;
                // EDGE11
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Px( n_node ) > xmin && vset.Px( n_node ) < xmax )
                     (*bit) = FRONT_TOP;
                // EDGE12
                if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) &&
                     vset.Py( n_node ) > ymin && vset.Py( n_node ) < ymax )
                     (*bit) = FRONT_LEFT;
             }
          
          // Edges along the right side (Edges 6 and 7)
          if ( approximatelyEqual(vset.Px( n_node ),xmax)  )
             {
                // EDGE6
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = BOTTOM_RIGHT;
                // EDGE7
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = TOP_RIGHT;
             }
          
          // Edges along the left side (Edges 5 and 8)
          if ( approximatelyEqual(vset.Px( n_node ),xmin)  )
             {
                // EDGE5
                if ( approximatelyEqual(vset.Py( n_node ),ymin,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = BOTTOM_LEFT;
                // EDGE8
                if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) &&
                     vset.Pz( n_node ) > zmin && vset.Pz( n_node ) < zmax )
                     (*bit) = TOP_LEFT;
             }
       }
       
 } // end flagEdges




/**
    Assigns corner-specific BOX_BOUNDARY flags (CNR1 to CNR8 in 3D, CNR1 to CNR4 in 2D)
    to nodes based on their coordinates at the corners of a bounding box.
 */
template<uint32_t dim>
void flagCornerNodes( VSet<dim>& vset, double tol,
                      double xmin, double xmax,
                      double ymin, double ymax,
                      double zmin, double zmax )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     // finding min-max of x, and y coordinates of boundary nodes
     size_t flagging_count(0U);
     
     // 2D CASE finding corners nodes on the basis of their coordinates
     // and flagging them accordingly
     if constexpr ( dim != 3 )
       {
         size_t n_node(0U);
         for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
           {
              if ( (*bit) < INTERNAL ) {
                   csmp_error.Note( WARNING, "VSetConverter<dim>::FlagCornerNodes (2D case)",
                                     "Skipping node, the BOX_boundary flag of which could not be identified" );
                   cerr <<"\nNode "<< n_node <<", flagged: "<< parseBoundary( static_cast<BOX_BOUNDARY>(*bit) ) << endl;
                }
              else if ( (*bit) != NOT )
                {
                   // CNR1
                   if ( approximatelyEqual( vset.Px( n_node ), xmin, tol ) && approximatelyEqual( vset.Py( n_node ), ymin ) )
                     { (*bit) = CNR_MIN; flagging_count++; }
                   // CNR2
                   else if ( approximatelyEqual( vset.Px( n_node ), xmax, tol ) && approximatelyEqual( vset.Py( n_node ), ymin ) )
                     { (*bit) = CNR_X; flagging_count++; }
                   // CNR3
                   else if ( approximatelyEqual( vset.Px( n_node ), xmax, tol ) && approximatelyEqual( vset.Py( n_node ), ymax ) )
                     { (*bit) = CNR_XY; flagging_count++; }
                   // CNR4
                   else if ( approximatelyEqual( vset.Px( n_node ), xmin, tol ) && approximatelyEqual( vset.Py( n_node ), ymax ) )
                     { (*bit) = CNR_Y; flagging_count++; }
                }
           }
         if ( flagging_count < 4U ) {
              csmp_error.Note( WARNING, "flagCornerNodes (2D case)",
                                          "Less than 4 corner nodes could be flagged" );
              cout <<"\nNumber of flagged nodes: "<< flagging_count << endl;
           }
         return;
       } // end 2D case

     // 3D CASE finding corners nodes on the basis of their coordinates
     // and flagging them accordingly
     if constexpr( dim == 3 )
       {
         size_t n_node(0U);
         for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
           {
              if ( (*bit) < INTERNAL ) {
                   csmp_error.Note( WARNING, "flagCornerNodes (3D case)",
                                             "Skipping node, the boundary flag of which could not be identified" );
                   cout <<"\nNode "<< n_node <<", flagged: "<< (*bit) << endl;
                }
              else if ( (*bit) != NOT )
                {
                   if ( approximatelyEqual( vset.Pz( n_node ), zmin, tol ) )
                     {
                        // CNR1
                        if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                          { (*bit) = CNR_MIN; flagging_count++; }
                        // CNR2
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                          { (*bit) = CNR_X; flagging_count++; }
                        // CNR3
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_XY; flagging_count++; }
                        // CNR4
                        else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_Y; flagging_count++; }
                     }
                   else if ( approximatelyEqual( vset.Pz( n_node ), zmax, tol ) )
                     {
                        // CNR5
                        if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                          { (*bit) = CNR_Z; flagging_count++; }
                        // CNR6
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                          { (*bit) = CNR_XZ; flagging_count++; }
                        // CNR7
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_MAX; flagging_count++; }
                        // CNR8
                        else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_Z; flagging_count++; }
                     }
                }
           }

         if ( flagging_count < 8U ) {
              csmp_error.Note( WARNING, "flagCornerNodes (3D case)",
                                         "Less than 8 nodes were identified as model corners" );
              cout <<"\nNumber of flagged nodes: "<< flagging_count << endl;
           }
           
      } // end 3D
        
 } // end FlagCornerNodes





/**
     Applies box boundary flags assuming that the model stored in the VSet is box-shaped.
     Brute force approach: based on the position the nodes are flagged as boundary nodes
     
     @param tol gives the precision with which the nodes must be matched.
     
// tested: O.K.
*/
template<uint32_t dim>
void establishBoundaryFlagsForBoxModel( VSet<dim>& vset, double tol )
 {
     // 0. getting rid of the existing boundary conditions
     // --------------------------------------------------
     vset.RemoveBflags();

     // 1. finding the dimensions of the box-shaped model
     //    assuming that the boundary have been identified
     //    but the flags are incorrect thus far
     // -------------------------------------------------
     double xmin, xmax, ymin, ymax, zmin, zmax;
     xmin = xmax = vset.Px(0);
     ymin = ymax = vset.Py(0);
     zmin = zmax = vset.Pz(0);

     double  x, y, z;

     for ( size_t i{0U}; i<vset.Vertices(); i++ )
       {
          x = vset.Px(i);
          y = vset.Py(i);
          z = vset.Pz(i);
          if ( x < xmin ) xmin = x;
          if ( x > xmax ) xmax = x;
          if ( y < ymin ) ymin = y;
          if ( y > ymax ) ymax = y;
          if ( z < zmin ) zmin = z;
          if ( z > zmax ) zmax = z;
       }

     // making new boundary conditions
     vector<std::int8_t>  new_bflags( vset.Vertices(), 0 );

     for ( size_t i{0U}; i<vset.Vertices(); i++ )  
       if ( approximatelyEqual(vset.Px(i),xmin,tol) ||
            approximatelyEqual(vset.Px(i),xmax,tol) ||
            approximatelyEqual(vset.Py(i),ymin,tol) ||
            approximatelyEqual(vset.Py(i),ymax,tol) ||
            approximatelyEqual(vset.Pz(i),zmin,tol) ||
            approximatelyEqual(vset.Pz(i),zmax,tol) )
         { 
            new_bflags[ i ] = UNSPECIFIED;
         }

     vset.AddBFlags( new_bflags.begin(), new_bflags.end() );

               
     // 2. flagging the sides of the model using the user-defined
     //    tolerances (the edges are dealt with later)
     // ---------------------------------------------------------
     size_t n_node(0U);
     for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++, n_node++ )
       {
          // bottom (y=ymin)
          if      ( approximatelyEqual(vset.Py( n_node ),ymin,tol) ) (*bit) = BOTTOM_OUTSIDE;
          // top    (y=ymax)
          else if ( approximatelyEqual(vset.Py( n_node ),ymax,tol) ) (*bit) = TOP_OUTSIDE;
          // left   (x=xmin)
          else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) ) (*bit) = LEFT_OUTSIDE;
          // right  (x=xmax)
          else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) ) (*bit) = RIGHT_OUTSIDE;
          // back   (z=zmin)
          else if ( approximatelyEqual(vset.Pz( n_node ),zmin,tol) ) (*bit) = BACK_OUTSIDE;
          // front  (z=zmax)
          else if ( approximatelyEqual(vset.Pz( n_node ),zmax,tol) ) (*bit) = FRONT_OUTSIDE;
       }
       
     // 3. now the corners & edges are flagged
     // --------------------------------------
     flagEdges( vset, tol, xmin, xmax, ymin, ymax, zmin, zmax );  
     flagCornerNodes( vset, tol, xmin, xmax, ymin, ymax, zmin, zmax );
       
 } // end

template void establishBoundaryFlagsForBoxModel( VSet<3U>&, double );
template void establishBoundaryFlagsForBoxModel( VSet<2U>&, double );
template void establishBoundaryFlagsForBoxModel( VSet<1U>&, double );




} // end namespace csmp
