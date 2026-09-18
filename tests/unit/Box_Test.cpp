#include "Box_Test.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "vsetMakers.h"
#include "CSMP_mathUtilities.h"
#include "CSMP_highLevelUtilities.h"
#include "CSMP_physical_constants.h"
#include "meshManagementUtilities.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

/**
     Tests the functionality of Box class and the flagging of related box-shaped models 
     
     @todo SKM add testing in 2D
*/
void Box_Test::run()
{
  // are the static functions parsing correctly ?
  TestParsingOfFlags();
  
  // testing method UnitNormal
  vector<double> unitNormal;
  Box().UnitNormalTo( LEFT, 1, unitNormal );
  _test( unitNormal.size() == 1 );
  _test( unitNormal.at(0) == -1 );
  Box().UnitNormalTo( RIGHT, 1, unitNormal );
  _test( unitNormal.size() == 1 );
  _test( unitNormal.at(0) == 1 );
  Box().UnitNormalTo( LEFT, 2, unitNormal );
  _test( unitNormal.size() == 2 );
  _test( unitNormal.at(0) == -1 );
  _test( unitNormal.at(1) == 0 );
  Box().UnitNormalTo( RIGHT, 2, unitNormal );
  _test( unitNormal.size() == 2 );
  _test( unitNormal.at(0) == 1 );
  _test( unitNormal.at(1) == 0 );
  Box().UnitNormalTo( TOP, 2, unitNormal );
  _test( unitNormal.size() == 2 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == 1 );
  Box().UnitNormalTo( BOTTOM, 2, unitNormal );
  _test( unitNormal.size() == 2 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -1 );
   Box().UnitNormalTo( LEFT, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -1 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == 0 );
   Box().UnitNormalTo( RIGHT, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 1 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == 0 );
   Box().UnitNormalTo( TOP, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == 1 );
  _test( unitNormal.at(2) == 0 );
   Box().UnitNormalTo( BOTTOM, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -1 );
  _test( unitNormal.at(2) == 0 );
   Box().UnitNormalTo( FRONT, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == 1 );
   Box().UnitNormalTo( BACK, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == -1 );
  // edges of 3D model
  const double sine45 = sin(45.0 * csmp::CSMP_PI / 180.0);
  Box().UnitNormalTo( EDGE1, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -sine45 );
  _test( unitNormal.at(2) == -sine45 );
  Box().UnitNormalTo( EDGE2, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sine45 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == -sine45 );
  Box().UnitNormalTo( EDGE3, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) ==  sine45 );
  _test( unitNormal.at(2) == -sine45 );
  Box().UnitNormalTo( EDGE4, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sine45 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == -sine45 );
  Box().UnitNormalTo( EDGE5, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sine45 );
  _test( unitNormal.at(1) == -sine45 );
  _test( unitNormal.at(2) == 0 );
  Box().UnitNormalTo( EDGE6, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sine45 );
  _test( unitNormal.at(1) == -sine45 );
  _test( unitNormal.at(2) == 0 );
  Box().UnitNormalTo( EDGE7, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sine45 );
  _test( unitNormal.at(1) == sine45 );
  _test( unitNormal.at(2) == 0 );
  Box().UnitNormalTo( EDGE8, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sine45 );
  _test( unitNormal.at(1) == sine45 );
  _test( unitNormal.at(2) == 0 );
   Box().UnitNormalTo( EDGE9, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -sine45 );
  _test( unitNormal.at(2) ==  sine45 );
  Box().UnitNormalTo( EDGE10, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sine45 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == sine45 );
  Box().UnitNormalTo( EDGE11, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == sine45 );
  _test( unitNormal.at(2) == sine45 );
  Box().UnitNormalTo( EDGE12, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sine45 );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == sine45 );
  
   // testing parse boundary
  _test( parseBoundary( NOT ) == "NOT" );
  _test( parseBoundary( IRREGULAR ) == "IRREGULAR" );
  _test( parseBoundary( LEFT ) == "LEFT" );
  _test( parseBoundary( RIGHT ) == "RIGHT" );
  _test( parseBoundary( BOTTOM ) == "BOTTOM" );
  _test( parseBoundary( TOP ) == "TOP" );
  _test( parseBoundary( FRONT ) == "FRONT" );
  _test( parseBoundary( BACK ) == "BACK" );
  _test( parseBoundary( CNR1 ) == "CNR1" );
  _test( parseBoundary( CNR2 ) == "CNR2" );
  _test( parseBoundary( CNR3 ) == "CNR3" );
  _test( parseBoundary( CNR4 ) == "CNR4" );
  _test( parseBoundary( CNR5 ) == "CNR5" );
  _test( parseBoundary( CNR6 ) == "CNR6" );
  _test( parseBoundary( CNR7 ) == "CNR7" );
  _test( parseBoundary( CNR8 ) == "CNR8" );
  _test( parseBoundary( EDGE1 ) == "EDGE1" );
  _test( parseBoundary( EDGE2 ) == "EDGE2" );
  _test( parseBoundary( EDGE3 ) == "EDGE3" );
  _test( parseBoundary( EDGE4 ) == "EDGE4" );
  _test( parseBoundary( EDGE5 ) == "EDGE5" );
  _test( parseBoundary( EDGE6 ) == "EDGE6" );
  _test( parseBoundary( EDGE7 ) == "EDGE7" );
  _test( parseBoundary( EDGE8 ) == "EDGE8" );
  _test( parseBoundary( EDGE9 ) == "EDGE9" );
  _test( parseBoundary( EDGE10 ) == "EDGE10" );
  _test( parseBoundary( EDGE11 ) == "EDGE11" );
  _test( parseBoundary( EDGE12 ) == "EDGE12" );
  _test( parseBoundary( INTERNAL ) == "INTERNAL" );
  
  // intToBOX_BOUNDARY - non automatic conversion
  _test( intToBOX_BOUNDARY( static_cast<int8_t>(LEFT_OUTSIDE) ) == LEFT );
  _test( intToBOX_BOUNDARY( -8 ) == CNR_MIN );
  _test( intToBOX_BOUNDARY( static_cast<int8_t>(-28) ) == REGION_BOUNDARY );

  // for vsetMaker model Pyra_Hexa with mixed element types (32elmts, 64 nodes): SKM tested bflags: OK (6/7/24)
  _test( TestWhetherSideBoundaryFlagsArePresent() );
  
  // _test( TestWhetherAllBoxFlagsArePresent() ); // fails because the corners are missing
  _test( TestBoundaryFlagging() ); // model Tetra
  _test( TestWhetherBoundaryFlagsArePreservedInBinaryFile1() ); // model Prism_Hexa
  _test( TestBoundaryVersusBOX_BOUNDARY_Flagging() ); // model FracBox
  
  // tests whether the function recreateBoxBoundaryFlags() manages to reconstruct edges and boundaries correctly
  _test( TestBoundaryFlagRecreation() );
  
  // tests whether the boundaries of a 2D box model and an ANSYS model are assigned correctly
  _test( TestBoundaryFlagAssigment2D() );

  // establishes model boundaries by using the outward pointing normals of the box model
  TestWhetherElementNormalsAreOutwardPointing();

} // end run



// Old functions pre-dating the compile-time versions used here to verify the new ones
/*
    std::string  parseBoundary( BOX_BOUNDARY i )
    BOX_BOUNDARY  parseBoundary( const string& i )

    bool isEdge( BOX_BOUNDARY bd )
    bool isCorner( BOX_BOUNDARY bd )
    isSide( BOX_BOUNDARY bd )

    isLEFT( BOX_BOUNDARY bd )
    isRIGHT( BOX_BOUNDARY bd )
    isTOP( BOX_BOUNDARY bd )
    bool isBOTTOM( BOX_BOUNDARY bd )
    bool isFRONT( BOX_BOUNDARY bd )
    bool isBACK( BOX_BOUNDARY bd )

    belongsToSide( BOX_BOUNDARY side, BOX_BOUNDARY bd )
    belongsToEdge( BOX_BOUNDARY edge, BOX_BOUNDARY bd )

    BOX_BOUNDARY intToBOX_BOUNDARY( int i )
*/
namespace box_test {

/// returns wether a node (or line element) lies on an edge of the model
static bool isEdge( BOX_BOUNDARY bd )
 {
    if ( bd == EDGE1 )  return true;
    if ( bd == EDGE2 )  return true;
    if ( bd == EDGE3 )  return true;
    if ( bd == EDGE4 )  return true;
    if ( bd == EDGE5 )  return true;
    if ( bd == EDGE6 )  return true;
    if ( bd == EDGE7 )  return true;
    if ( bd == EDGE8 )  return true;
    if ( bd == EDGE9 )  return true;
    if ( bd == EDGE10 )  return true;
    if ( bd == EDGE11 )  return true;
    if ( bd == EDGE12 )  return true;
    return false;
 }

/// of rectangular (brick-shaped) model; @test OK
static bool isCorner( BOX_BOUNDARY bd )
 {
    if ( bd == CNR1 )  return true;
    if ( bd == CNR2 )  return true;
    if ( bd == CNR3 )  return true;
    if ( bd == CNR4 )  return true;
    if ( bd == CNR5 )  return true;
    if ( bd == CNR6 )  return true;
    if ( bd == CNR7 )  return true;
    if ( bd == CNR8 )  return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isSide( BOX_BOUNDARY bd )
 {
    if ( bd == RIGHT )  return true;
    if ( bd == LEFT )   return true;
    if ( bd == TOP )    return true;
    if ( bd == BOTTOM ) return true;
    if ( bd == FRONT )  return true;
    if ( bd == BACK )   return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isLEFT( BOX_BOUNDARY bd )
 {
    if ( bd == LEFT )   return true;
    if ( bd == EDGE4 )  return true;
    if ( bd == EDGE8 )  return true;
    if ( bd == EDGE5 )  return true;
    if ( bd == EDGE12 ) return true;
    if ( bd == CNR1 )  return true;
    if ( bd == CNR4 )  return true;
    if ( bd == CNR5 )  return true;
    if ( bd == CNR8 )  return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isRIGHT( BOX_BOUNDARY bd )
 {
    if ( bd == RIGHT )  return true;
    if ( bd == EDGE2 )  return true;
    if ( bd == EDGE6 )  return true;
    if ( bd == EDGE7 )  return true;
    if ( bd == EDGE10 ) return true;
    if ( bd == CNR2 )  return true;
    if ( bd == CNR3 )  return true;
    if ( bd == CNR6 )  return true;
    if ( bd == CNR7 )  return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isTOP( BOX_BOUNDARY bd )
 {
    if ( bd == TOP )    return true;
    if ( bd == EDGE3 )  return true;
    if ( bd == EDGE8 )  return true;
    if ( bd == EDGE7 )  return true;
    if ( bd == EDGE11 ) return true;
    if ( bd == CNR3 )  return true;
    if ( bd == CNR4 )  return true;
    if ( bd == CNR7 )  return true;
    if ( bd == CNR8 )  return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isBOTTOM( BOX_BOUNDARY bd )
 {
    if ( bd == BOTTOM ) return true;
    if ( bd == EDGE1 )  return true;
    if ( bd == EDGE5 )  return true;
    if ( bd == EDGE6 )  return true;
    if ( bd == EDGE9 )  return true;
    if ( bd == CNR1 )  return true;
    if ( bd == CNR2 )  return true;
    if ( bd == CNR5 )  return true;
    if ( bd == CNR6 )  return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isFRONT( BOX_BOUNDARY bd )
 {
    if ( bd == FRONT )  return true;
    if ( bd == EDGE9 )  return true;
    if ( bd == EDGE10 ) return true;
    if ( bd == EDGE11 ) return true;
    if ( bd == EDGE12 ) return true;
    if ( bd == CNR5 )  return true;
    if ( bd == CNR6 )  return true;
    if ( bd == CNR7 )  return true;
    if ( bd == CNR8 )  return true;
    return false;
 }


/// of rectangular (brick-shaped) model; @test OK
static bool isBACK( BOX_BOUNDARY bd )
 {
    if ( bd == BACK )   return true;
    if ( bd == EDGE1 )  return true;
    if ( bd == EDGE2 )  return true;
    if ( bd == EDGE3 )  return true;
    if ( bd == EDGE4 )  return true;
    if ( bd == CNR1 )  return true;
    if ( bd == CNR2 )  return true;
    if ( bd == CNR3 )  return true;
    if ( bd == CNR4 )  return true;
    return false;
 }


static BOX_BOUNDARY intToBOX_BOUNDARY( int i )
  {
    if ( i == NOT )             return NOT;
    if ( i == IRREGULAR )       return IRREGULAR;
    if ( i == LEFT_OUTSIDE )    return LEFT;
    if ( i == RIGHT_OUTSIDE )   return RIGHT;   
    if ( i == BOTTOM_OUTSIDE )  return BOTTOM;  
    if ( i == TOP_OUTSIDE )     return TOP;       
    if ( i == FRONT_OUTSIDE )   return FRONT;   
    if ( i == BACK_OUTSIDE )    return BACK;
    if ( i == CNR_MIN )         return CNR1;
    if ( i == CNR_X )           return CNR2;
    if ( i == CNR_XY )          return CNR3;
    if ( i == CNR_Y )           return CNR4;
    if ( i == CNR_Z )           return CNR5;
    if ( i == CNR_XZ )          return CNR6;
    if ( i == CNR_MAX )         return CNR7;
    if ( i == CNR_YZ )          return CNR8;
    if ( i == BACK_BOTTOM )     return EDGE1;
    if ( i == BACK_RIGHT )      return EDGE2;
    if ( i == BACK_TOP )        return EDGE3;
    if ( i == BACK_LEFT )       return EDGE4;
    if ( i == BOTTOM_LEFT )     return EDGE5;
    if ( i == BOTTOM_RIGHT )    return EDGE6;
    if ( i == TOP_RIGHT )       return EDGE7;
    if ( i == TOP_LEFT )        return EDGE8;
    if ( i == FRONT_BOTTOM )    return EDGE9;
    if ( i == FRONT_RIGHT )     return EDGE10;
    if ( i == FRONT_TOP )       return EDGE11;
    if ( i == FRONT_LEFT )      return EDGE12;
    if ( i == REGION_BOUNDARY ) return INTERNAL;
    if ( i == MULTIPLE )        return MULTIPLE;

    //cout <<"\nintToSG_BOUNDARY(int): unable to parse integer: "<< i << endl;
    return NOT;   
  }


static std::string  parseBoundary( BOX_BOUNDARY i )
 {
    if ( i == NOT )      return string("NOT");
    if ( i == IRREGULAR )return string("IRREGULAR");    
    if ( i == LEFT )     return string("LEFT");   
    if ( i == RIGHT )    return string("RIGHT");   
    if ( i == BOTTOM )   return string("BOTTOM");  
    if ( i == TOP )      return string("TOP");       
    if ( i == FRONT )    return string("FRONT");   
    if ( i == BACK )     return string("BACK");
    if ( i == CNR1 )     return string("CNR1");
    if ( i == CNR2 )     return string("CNR2");
    if ( i == CNR3 )     return string("CNR3");
    if ( i == CNR4 )     return string("CNR4");
    if ( i == CNR5 )     return string("CNR5");
    if ( i == CNR6 )     return string("CNR6");
    if ( i == CNR7 )     return string("CNR7");
    if ( i == CNR8 )     return string("CNR8");
    if ( i == EDGE1 )    return string("EDGE1");
    if ( i == EDGE2 )    return string("EDGE2");
    if ( i == EDGE3 )    return string("EDGE3");
    if ( i == EDGE4 )    return string("EDGE4");
    if ( i == EDGE5 )    return string("EDGE5");
    if ( i == EDGE6 )    return string("EDGE6");
    if ( i == EDGE7 )    return string("EDGE7");
    if ( i == EDGE8 )    return string("EDGE8");
    if ( i == EDGE9 )    return string("EDGE9");
    if ( i == EDGE10 )   return string("EDGE10");
    if ( i == EDGE11 )   return string("EDGE11");
    if ( i == EDGE12 )   return string("EDGE12");
    if ( i == INTERNAL ) return string("INTERNAL");
    if ( i == MULTIPLE ) return string("MULTIPLE");

    // Unable to parse BOX_BOUNDARY: returning NOT
    return string("NOT");
 }


/// text to boundary enum
static BOX_BOUNDARY  parseBoundary( const string& i )
 {
    if ( i == "NOT" )      return NOT;
    if ( i == "IRREGULAR" )return IRREGULAR;
    if ( i == "LEFT" )     return LEFT;
    if ( i == "RIGHT" )    return RIGHT;
    if ( i == "BOTTOM" )   return BOTTOM;
    if ( i == "TOP" )      return TOP;
    if ( i == "FRONT" )    return FRONT;
    if ( i == "BACK" )     return BACK;
    if ( i == "CNR1" )     return CNR1;
    if ( i == "CNR2" )     return CNR2;
    if ( i == "CNR3" )     return CNR3;
    if ( i == "CNR4" )     return CNR4;
    if ( i == "CNR5" )     return CNR5;
    if ( i == "CNR6" )     return CNR6;
    if ( i == "CNR7" )     return CNR7;
    if ( i == "CNR8" )     return CNR8;
    if ( i == "EDGE1" )    return EDGE1;
    if ( i == "EDGE2" )    return EDGE2;
    if ( i == "EDGE3" )    return EDGE3;
    if ( i == "EDGE4" )    return EDGE4;
    if ( i == "EDGE5" )    return EDGE5;
    if ( i == "EDGE6" )    return EDGE6;
    if ( i == "EDGE7" )    return EDGE7;
    if ( i == "EDGE8" )    return EDGE8;
    if ( i == "EDGE9" )    return EDGE9;
    if ( i == "EDGE10" )   return EDGE10;
    if ( i == "EDGE11" )   return EDGE11;
    if ( i == "EDGE12" )   return EDGE12;
    if ( i == "INTERNAL" ) return INTERNAL;
    if ( i == "MULTIPLE" ) return MULTIPLE;

    // unable to parse BOX_BOUNDARY: returning NOT
    return NOT;   
 }

/**
    Returns true if @param bd belongs to the corresponding edge.
    Example: CNR1 and CNR2 belong to EDGE1. 
*/
static bool belongsToEdge( BOX_BOUNDARY edge, BOX_BOUNDARY bd )
 {
    if      ( edge == EDGE1 and (bd == EDGE1 or bd == CNR1 or bd == CNR2) ) return true;
    else if ( edge == EDGE2 and (bd == EDGE2 or bd == CNR2 or bd == CNR3) ) return true;
    else if ( edge == EDGE3 and (bd == EDGE3 or bd == CNR3 or bd == CNR4) ) return true;
    else if ( edge == EDGE4 and (bd == EDGE4 or bd == CNR1 or bd == CNR4) ) return true;
    else if ( edge == EDGE5 and (bd == EDGE5 or bd == CNR1 or bd == CNR5) ) return true;
    else if ( edge == EDGE6 and (bd == EDGE6 or bd == CNR2 or bd == CNR6) ) return true;
    else if ( edge == EDGE7 and (bd == EDGE7 or bd == CNR3 or bd == CNR7) ) return true;
    else if ( edge == EDGE8 and (bd == EDGE8 or bd == CNR4 or bd == CNR8) ) return true;
    else if ( edge == EDGE9 and (bd == EDGE9 or bd == CNR5 or bd == CNR6) ) return true;
    else if ( edge == EDGE10 and (bd == EDGE10 or bd == CNR6 or bd == CNR7) ) return true;
    else if ( edge == EDGE11 and (bd == EDGE11 or bd == CNR7 or bd == CNR8) ) return true;
    else if ( edge == EDGE12 and (bd == EDGE12 or bd == CNR5 or bd == CNR8) ) return true;
    return false;
 }


/**
    Returns true if @param bd belongs to the corresponding side of the box-shaped model.
    Example: CNR1 - CNR4 belong to BACK. 
*/
static bool belongsToSide( BOX_BOUNDARY side, BOX_BOUNDARY bd )
 {
    if      ( side == LEFT )   return box_test::isLEFT( bd );
    else if ( side == RIGHT )  return box_test::isRIGHT( bd );
    else if ( side == TOP )    return box_test::isTOP( bd );
    else if ( side == BOTTOM ) return box_test::isBOTTOM( bd );
    else if ( side == FRONT )  return box_test::isFRONT( bd );
    else if ( side == BACK )   return box_test::isBACK( bd );
    return false;
 }

} // end namespace box_test


void Box_Test::TestParsingOfFlags()
 {
   // dumb basics
   _test( isLEFT(LEFT) == true );
   _test( isRIGHT(RIGHT) == true );
   _test( isTOP(TOP) == true );
   _test( isBOTTOM(BOTTOM) == true );
   _test( isFRONT(FRONT) == true );
   _test( isBACK(BACK) == true );
   _test( canBeIRREGULAR(IRREGULAR) == true );
   
   // and in more detail: isBACK() EDGE1, EDGE2, EDGE3, EDGE4, CNR1, CNR2, CNR3, CNR4
   _test( isBACK(EDGE1) == true );
   _test( isBACK(EDGE2) == true );   
   _test( isBACK(EDGE3) == true );
   _test( isBACK(EDGE4) == true );
   _test( isBACK(CNR1) == true );
   _test( isBACK(CNR2) == true );
   _test( isBACK(CNR3) == true );
   _test( isBACK(CNR4) == true );
   
   // the main diagnostics tool
   _test( parseBoundary( BOTTOM ) == "BOTTOM" );
   _test( parseBoundary( RIGHT ) == "RIGHT" );
   _test( parseBoundary( TOP ) == "TOP" );
   _test( parseBoundary( LEFT ) == "LEFT" );
   _test( parseBoundary( BACK ) == "BACK" );
   _test( parseBoundary( FRONT ) == "FRONT" );
   _test( parseBoundary( IRREGULAR ) == "IRREGULAR" );
   _test( parseBoundary( CNR2 ) == "CNR2" );
   _test( parseBoundary( CNR5 ) == "CNR5" );
   _test( parseBoundary( EDGE3 ) == "EDGE3" );
   _test( parseBoundary( EDGE6 ) == "EDGE6" );

   _test( parseBoundary( "BOTTOM" ) == BOTTOM );
   _test( parseBoundary( "RIGHT" ) == RIGHT );
   _test( parseBoundary( "TOP" ) == TOP );
   _test( parseBoundary( "LEFT" ) == LEFT );
   _test( parseBoundary( "BACK" ) == BACK );
   _test( parseBoundary( "FRONT" ) == FRONT );
   _test( parseBoundary( "IRREGULAR" ) == IRREGULAR );
   _test( parseBoundary( "CNR2" ) == CNR2 );
   _test( parseBoundary( "CNR5" ) == CNR5 );
   _test( parseBoundary( "EDGE3" ) == EDGE3 );
   _test( parseBoundary( "EDGE6" ) == EDGE6 );

    // for all possible flag values
    for ( int8_t i{3}; i>MULTIPLE_BOUNDARIES-1; --i )
      {
         auto flag            = static_cast<BOX_BOUNDARY>(i);
         auto boundary_string = box_test::parseBoundary(flag);
         // enum argument
         _test( parseBoundary( flag ) == box_test::parseBoundary(flag) );
         if ( parseBoundary( flag ) != box_test::parseBoundary(flag) )
           cout <<" "<<"parseBoundary(BOX_BOUNDARY) used with: "<< boundary_string <<" returned "<< parseBoundary( flag ) << endl;
           
         // string argument
         _test( parseBoundary( boundary_string ) == box_test::parseBoundary( boundary_string ) );
         if ( parseBoundary( boundary_string ) != box_test::parseBoundary( boundary_string ) )
           cout <<" "<<"parseBoundary(string) used with: "<< boundary_string << endl;

         _test( isLEFT( flag ) == box_test::isLEFT( flag ) );
         _test( isRIGHT( flag ) == box_test::isRIGHT( flag ) );
         _test( isTOP( flag ) == box_test::isTOP( flag ) );
         _test( isBOTTOM( flag ) == box_test::isBOTTOM( flag ) );
         _test( isFRONT( flag ) == box_test::isFRONT( flag ) );
         _test( isBACK( flag ) == box_test::isBACK( flag ) );

         _test( isEdge( flag ) == box_test::isEdge( flag ) );
         _test( isCorner( flag ) == box_test::isCorner( flag ) );
          if ( isCorner( flag ) != box_test::isCorner( flag ) ) cout <<" "<< boundary_string << endl;
         
         // old function is more restrictive
         if ( box_test::isSide( flag ) )
         _test( isSide( flag ) == box_test::isSide( flag ) );

         _test( belongsToSide( BOTTOM, flag ) == box_test::belongsToSide( BOTTOM, flag ) );
         _test( belongsToSide( RIGHT, flag )  == box_test::belongsToSide( RIGHT, flag ) );
         _test( belongsToSide( TOP, flag )    == box_test::belongsToSide( TOP, flag ) );
         _test( belongsToSide( LEFT, flag )   == box_test::belongsToSide( LEFT, flag ) );
         _test( belongsToSide( FRONT, flag )  == box_test::belongsToSide( FRONT, flag ) );
         _test( belongsToSide( BACK, flag )   == box_test::belongsToSide( BACK, flag ) );

         _test( belongsToEdge( EDGE1, flag ) == box_test::belongsToEdge( EDGE1, flag ) );
         _test( belongsToEdge( EDGE2, flag ) == box_test::belongsToEdge( EDGE2, flag ) );
         _test( belongsToEdge( EDGE3, flag ) == box_test::belongsToEdge( EDGE3, flag ) );
         _test( belongsToEdge( EDGE4, flag ) == box_test::belongsToEdge( EDGE4, flag ) );
         _test( belongsToEdge( EDGE5, flag ) == box_test::belongsToEdge( EDGE5, flag ) );
         _test( belongsToEdge( EDGE6, flag ) == box_test::belongsToEdge( EDGE6, flag ) );
         _test( belongsToEdge( EDGE7, flag ) == box_test::belongsToEdge( EDGE7, flag ) );
         _test( belongsToEdge( EDGE8, flag ) == box_test::belongsToEdge( EDGE8, flag ) );
         _test( belongsToEdge( EDGE9, flag ) == box_test::belongsToEdge( EDGE9, flag ) );
         _test( belongsToEdge( EDGE10, flag ) == box_test::belongsToEdge( EDGE10, flag ) );
         _test( belongsToEdge( EDGE11, flag ) == box_test::belongsToEdge( EDGE11, flag ) );
         _test( belongsToEdge( EDGE12, flag ) == box_test::belongsToEdge( EDGE12, flag ) );

         _test( intToBOX_BOUNDARY(i) == box_test::intToBOX_BOUNDARY( static_cast<int>(i) ) );
         if ( intToBOX_BOUNDARY(i) != box_test::intToBOX_BOUNDARY( static_cast<int>(i) ) ) cout <<" "<< boundary_string << endl;
      }
          
 } // end TestParsingOfFlags
  





bool Box_Test::TestBoundaryFlagAssigment2D()
 {
    // 1. standard model construction
    // ------------------------------
    VSet<2U> vset;
    create_TrianglePatch_VSet( vset );
    Model<2U>          model( vset, "CSMP-variables.txt" );
    const Region<2U>&  mregion(model.Region("Model"));
   
    // storing the flags in node and element order in a list for comparison
    set<BOX_BOUNDARY>  rectangle_flags;
    //rectangle_flags.push_back(LEFT); // not present in this model
    rectangle_flags.insert(RIGHT);
    rectangle_flags.insert(TOP);
    //rectangle_flags.push_back(BOTTOM); // not present in this model
    rectangle_flags.insert(CNR1);
    rectangle_flags.insert(CNR2);
    rectangle_flags.insert(CNR3);
    rectangle_flags.insert(CNR4);
   
    set<BOX_BOUNDARY>  actual_node_flags;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      actual_node_flags.insert( (*nit)->AtBoundary() );
    // checking the flags
    bool flags_are_correct(true);
    for ( auto it=rectangle_flags.begin(); it!=rectangle_flags.end(); it++ ) {
         if ( actual_node_flags.find( (*it) ) == actual_node_flags.end() ) {
              flags_are_correct = false;
              break;
           }
      }
    _test( flags_are_correct == true );
   
    // 2. ANSYS_Model2D
    // ----------------
    ANSYS_Model2D  model2( "BoxHalfs2D", "CSMP-variables.txt" );
    Region<2>&     mregion2( model2.Region("Model") );
    // now we should also have top and bottom boundaries
    rectangle_flags.clear();
    rectangle_flags.insert(RIGHT);
    rectangle_flags.insert(TOP);
    rectangle_flags.insert(LEFT);
    rectangle_flags.insert(BOTTOM);
    actual_node_flags.clear();
    flags_are_correct = true;

    // visualising the flags that were created earlier
    /*
    const string node_variable("nodal box flag"), elmt_variable("element box flag");
    model2.CreateProperty( node_variable.c_str(), "flag", SCALAR, NODE );
    model2.CreateProperty( elmt_variable.c_str(), "flag", SCALAR, ELEMENT );
    boxFlagsToVariable( model2, node_variable.c_str(), elmt_variable.c_str() );
    VTU_Interface<2>  vtu2(model2);
    std::list<std::string> outputProps;
    outputProps.push_back(node_variable.c_str());
    outputProps.push_back(elmt_variable.c_str());
    vtu2.OutputDataToVTU( "box_flags", outputProps, "Model", 2 );
    */
     
    for ( auto nit=mregion2.NodesBegin(); nit!=mregion2.NodesEnd(); nit++ )
      actual_node_flags.insert( (*nit)->AtBoundary() );
    for ( auto it=rectangle_flags.begin(); it!=rectangle_flags.end(); it++ ) {
         if ( actual_node_flags.find( (*it) ) == actual_node_flags.end() ) {
              flags_are_correct = false;
              break;
           }
      }
    _test( flags_are_correct == true );

    return flags_are_correct;
 }


/**
     converts flags to scalar variables that are subsequently tested using node numbers.
*/
bool Box_Test::TestBoundaryFlagging()
 {
    VSet<3U>  vset;
    create_Tetra_VSet( vset );
    Model<3U> model( vset, "CSMP-variables.txt" );
   
    string node_variable("nodal box flag"), elmt_variable("element box flag");
    model.CreateProperty( node_variable.c_str(), "nbf", "flag", SCALAR, NODE );
//    model.UpdateIndices();
    model.CreateProperty( elmt_variable.c_str(), "ebf", "flag", SCALAR, ELEMENT );
//    model.UpdateIndices();
   
    boxFlagsToVariable( model, node_variable.c_str(), elmt_variable.c_str() );

    if ( verbose_ ) {
         VTU_Interface<3>  vtu(model);
         std::list<std::string> outputProps;
         outputProps.push_back(node_variable.c_str());
         outputProps.push_back(elmt_variable.c_str());
         vtu.OutputDataToVTU( "box_flags", outputProps, "Model", 1 );
      }

    return true;
  
 } // end TestBoundaryFlagging




static void printNeighboursOfElement( const VSet<3U>& vset, size_t elmt )
 {
    assert( elmt < vset.Elements() == true );
    assert( vset.HybridElementTypeMesh() == true ); // expecting 'pelmt' size > 1

    const auto n_nbors = distance( vset.PfvertsBegin(elmt), vset.PfvertsEnd(elmt) );
    const auto e_type  = vset.ElementType( elmt );

    cout <<"\nprintNeighboursOfElement: element "<< elmt << endl;
    cout <<"\t"<<"faces and their nodes:"<< endl;
    for ( uint32_t face{0u}; face<n_nbors; ++face ) {
          const uint32_t nodes_per_face = CSMP_ElementSpecifications::NodesPerFaceForElementOfType( e_type, face );
          vector<uint32_t> fnids( nodes_per_face );
          for ( uint32_t n{0U}; n<nodes_per_face; ++n )
            fnids[n] = CSMP_ElementSpecifications::FaceNodeForElementOfType( e_type, face, n );
          // printing the relevant information
          cout <<"\t\t"<< face <<": ";
          for ( const auto& i : fnids ) cout <<" "<< vset.Plist( elmt, i );
          cout << endl;
    }
    
 } // end printNeighboursOfElement



/**
    Comparing the unit normals with the those of the sides 
    of the box-shaped model. 
    The expectation is that are are pointing in the direction
    as the normals of the elements or faces on the outside boundary
    of the model.
*/
void Box_Test::TestWhetherElementNormalsAreOutwardPointing()
 {
    VSet<3U>  vset;
    const bool bSkewed{false};
    create_Prism_Hexa_VSet( vset, bSkewed ); // correct neighbors & node flagging
//    create_Pyramid_Hexa_VSet( vset, bSkewed ); // correct neighbors & node flagging
    vset.InitialiseNodeTopologyIdentifiers(); // needed later for error reporting

    // checking the Face nodes of element 1
    const size_t element{1ul};
    printNeighboursOfElement( vset, element );
    
    if ( verbose_ ) cout <<"\nBox_Test::TestWhetherElementNormalsAreOutwardPointing: building model 'Prism_Hexa'"<< endl;
    Model<3U>  model( vset, "CSMP-variables.txt" );
    // ANSYS_Model3D model( "prism_test", "CSMP-variables.txt", true ); 
    
    printModelDimensions( model );
    Region<3U>& model_domain(model.Region("Model"));
    // verifying that the perimeter of the Model region matches that of the overall model
    const double model_surface_area = model_domain.SurfaceArea();
    _equal( model_surface_area, 6. * 3. * 3., 10. ); // 6-faces with 9m2, tolerance=10 eps
    
    // testing that Bflags are correct in original model 'Prism_Hexa'
    model.CreateProperty( "node flag", "NF", "none" );
    model.CreateProperty( "element flag", "EF", "none", SCALAR, ELEMENT );
    boxFlagsToVariable( model, "node flag", "element flag" ); // element number initialised by VSetMaker
    list<string> out_vars{ "node flag", "element flag", "element number" };
    VTU_Interface<3U> vtu_output( model );
    string file_name = "PrismHexa";
    vtu_output.OutputDataToVTU( file_name, out_vars, "Model", 0 );
    _test( isStrictlyBoxShaped(model) == true );
    
    // testing that the normals of elements of the model are outward-pointing
    for ( const auto& eit : model_domain.CellVector() )
      _test( areUnitNormalsToFacesAreOutwardPointing(eit) == true );
      
    // visually
    if ( verbose_ ) {
       for ( const auto& eit : model_domain.CellVector() ) {
             if ( eit->FE_Type() == ISOPARAMETRIC_LINEAR_PRISM ) {
                  string filename{"prism"}; filename += numberToString( eit->Idx() );
                  writePrismWithNormalsToVTK( eit, filename );
               }
             else if ( eit->FE_Type() == ISOPARAMETRIC_LINEAR_PYRAMID ) {
                  string filename{"pyramid"}; filename += numberToString( eit->Idx() );
                  writePyramidWithNormalsToVTK( eit, filename );
               }
         }
    }
    
    
    // 0. Creating Boundary subdomains relying on Element-face normals
    // ---------------------------------------------------------------
    _test( model.EstablishBoxBoundariesFromOrientation() );
    _test( isStrictlyBoxShaped(model) == true );
    
    // 0.1 testing that the boundary Face normals are pointing in the same direction as inner element face normals
    for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit ) {
        // normal to boundary Face
        Point<3> Face_nrml  = (*fit).UnitNormal();
        // normal to boundary facing face of inner element
        vector<double> eface_normal;
        (*fit).InnerParent()->UnitNormalToFace( (*fit).InnerParentFaceID(), eface_normal );
        Point<3> eface_nrml( eface_normal );
        if ( approximatelyEqual( dotProduct(Face_nrml,eface_nrml), 1. ) == false ) {
            cout <<"\n\t"<<"Face "<< (*fit).Idx() <<": normal "<< Face_nrml <<" vs "<<"inner Element ";
            cout << (*fit).InnerParent()->Idx() <<": normal "<< eface_nrml << endl;
            // (*fit).Out();
          }
      }
      
     // 0.2 do all FiniteVolumePolicy unit normal calculations give the same result?
     for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit ) {
          // UnitNormal( VectorVariable<dim>& nrml )
          VectorVariable<3> VV_nrml;
          (*fit).UnitNormal( VV_nrml );
          // UnitNormal( vector<double>& nrml ) -> used by EstablishBoxBoundariesFromOrientation()
          vector<double> vec_nrml;
          (*fit).UnitNormal( vec_nrml );
          // Point<dim>  UnitNormal()
          Point<3> pt_nrml = (*fit).UnitNormal();
          for ( uint32_t i{0}; i<3; ++ i ) _test( approximatelyEqual(pt_nrml[i],VV_nrml[i]) == true );
          for ( uint32_t i{0}; i<3; ++ i ) _test( approximatelyEqual(pt_nrml[i],vec_nrml[i]) == true );
       }

    // 0.3 creating Face normal for testing
    {
      const csmp::Index nrml_key = model.Database().StorageKey("face vector");
      size_t n_face{0u};
      for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit, ++n_face ) {
           // vector (unit normal)
           auto nrml = (*fit).UnitNormal();
           (*fit).Store( nrml_key, makeVector(ANY,ANY,ANY,nrml[0],nrml[1],nrml[2]) );
       }
      // Output boundaries with normals to VTU
      vtu_output.OutputDataToVTU( "Face-normals", "face vector", model.Boundary("BACK"), 0 );
      vtu_output.OutputDataToVTU( "Face-normals", "face vector", model.Boundary("BOTTOM"), 0 );
      vtu_output.OutputDataToVTU( "Face-normals", "face vector", model.Boundary("RIGHT"), 0 );
      vtu_output.OutputDataToVTU( "Face-normals", "face vector", model.Boundary("TOP"), 0 );
      vtu_output.OutputDataToVTU( "Face-normals", "face vector", model.Boundary("LEFT"), 0 );
      vtu_output.OutputDataToVTU( "Face-normals", "face vector", model.Boundary("FRONT"), 0 );
   }

     // 0.4 visual checks: are nodal bflags OK for Prism_Hexa
     boxFlagsToVariable( model, "node flag", "element flag" ); // element number initialised by VSetMaker
     file_name = "PrismHexa-rebuilt";
     vtu_output.OutputDataToVTU( file_name, out_vars, "Model", 0 );
     // outputting at problematic boundaries BACK and FRONT where prism elements surface
     vtu_output.OutputDataToVTU( file_name, "face variable", model.Boundary("BACK"), 0 );
     vtu_output.OutputDataToVTU( file_name, "face variable", model.Boundary("FRONT"), 0 );


    // 1. Testing that the boundaries are correctly flagged
    // ----------------------------------------------------
    {
      const Boundary<3U>& back   = model.Boundary("BACK");
      const Boundary<3U>& bottom = model.Boundary("BOTTOM");
      const Boundary<3U>& right  = model.Boundary("RIGHT");
      const Boundary<3U>& top    = model.Boundary("TOP");
      const Boundary<3U>& left   = model.Boundary("LEFT");
      const Boundary<3U>& front  = model.Boundary("FRONT");
      bool wrong_flag{false};
      
      for ( auto nit=back.NodesBegin(); nit!=back.NodesEnd(); ++nit ) {
           _test( isBACK( (*nit)->AtBoundary() ) == true );
           if ( isBACK( (*nit)->AtBoundary() ) == false ) {
                cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
                wrong_flag = true;
             }
        }
      if ( wrong_flag ) {
           // checking nodes
           cout <<"\nflags of boundary BACK: "<< endl;
           for ( auto nit=back.NodesBegin(); nit!=back.NodesEnd(); ++nit )
             cout <<"  "<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() ) <<" z="<< (*nit)->z();
           cout << endl;
           // checking faces
           cout <<"\nFace finite-element types of boundary BACK: "<< endl;
           for ( const auto& eit : back.CellVector() ) {
                cout <<"  "<< eit->Idx() <<": "<< parseAbbreviated_FE_Type( eit->FE_Type() );
                cout <<", npe: "<< eit->Nodes();
             }
           cout << endl << endl;
           wrong_flag = false;
        }
        
      for ( auto nit=bottom.NodesBegin(); nit!=bottom.NodesEnd(); ++nit ) {
           _test( isBOTTOM( (*nit)->AtBoundary() ) == true );
           if ( !isBOTTOM( (*nit)->AtBoundary() ) ) {
                cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
             }
        }
      for ( auto nit=right.NodesBegin(); nit!=right.NodesEnd(); ++nit ) {
           _test( isRIGHT( (*nit)->AtBoundary() ) == true );
           if ( !isRIGHT( (*nit)->AtBoundary() ) ) {
                cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
             }
        }
      for ( auto nit=top.NodesBegin(); nit!=top.NodesEnd(); ++nit ) {
           _test( isTOP( (*nit)->AtBoundary() ) == true );
           if ( !isTOP( (*nit)->AtBoundary() ) ) {
                cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
             }
        }
      for ( auto nit=left.NodesBegin(); nit!=left.NodesEnd(); ++nit ) {
           _test( isLEFT( (*nit)->AtBoundary() ) == true );
           if ( !isLEFT( (*nit)->AtBoundary() ) ) {
                cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
             }
        }
      for ( auto nit=front.NodesBegin(); nit!=front.NodesEnd(); ++nit ) {
           _test( isFRONT( (*nit)->AtBoundary() ) == true );
           if ( !isFRONT( (*nit)->AtBoundary() ) ) {
                cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
             }
        }
        
      if ( model.ContainsBoundary("IRREGULAR") ) {
           const Boundary<3U>& irregular = model.Boundary("IRREGULAR");
           for ( auto nit=irregular.NodesBegin(); nit!=irregular.NodesEnd(); ++nit ) {
                _test( canBeIRREGULAR( (*nit)->AtBoundary() ) == true );
                 if ( !canBeIRREGULAR( (*nit)->AtBoundary() ) ) {
                      cout <<"\t"<< (*nit)->Idx() <<": "<< parseBoundary( (*nit)->AtBoundary() );
                      cout <<", geometry: "<< parseTopology( (*nit)->Attribute() ) << endl << endl;
                   }
               }
        }
    } // end testing boundary flags
    
    
    // 1. testing the unit normals of the (volumetric elements)
    // --------------------------------------------------------
    vector<double> leftNormal, rightNormal, topNormal, bottomNormal, frontNormal, backNormal, eUnitNormal;
    Box().UnitNormalTo( LEFT,   3, leftNormal );
    Box().UnitNormalTo( RIGHT,  3, rightNormal );
    Box().UnitNormalTo( TOP,    3, topNormal );
    Box().UnitNormalTo( BOTTOM, 3, bottomNormal );
    Box().UnitNormalTo( FRONT,  3, frontNormal );
    Box().UnitNormalTo( BACK,   3, backNormal );

    double accumulated_area{ 0. };
    for ( size_t i=model_domain.InteriorCells(); i<model_domain.Cells(); ++i )
       {
         for ( uint32_t j=0U; j<model_domain.PerimeterFaces(i); ++j ) {
                const BOX_BOUNDARY flag = model_domain.E(i)->AtBoundary( model_domain.PerimeterFace(i,j) );
                _test( flag != NOT );
                if ( flag == NOT ) {
                     cout <<"\n"<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                     cout <<": Element face boundary identification is not correct.";
                     model_domain.E(i)->Out();
                  }
                // verifying that area of perimeter faces adds up to model surface area
                const auto          fnids = model_domain.E(i)->FE()->NodesOfFace( model_domain.PerimeterFace(i,j) );
                const CSMP_FEM_TYPE etype = model_domain.E(i)->FE()->ElementTypeOfFace( model_domain.PerimeterFace(i,j) );
                if ( isTriangular(etype) )
                  accumulated_area += triangleArea( model_domain.E(i)->N( fnids[2] )->Coordinate(), // counter-clockwise nodes
                                                    model_domain.E(i)->N( fnids[1] )->Coordinate(),
                                                    model_domain.E(i)->N( fnids[0] )->Coordinate() );
                else if ( isQuadrilateral(etype) )
                  accumulated_area += facetArea4( model_domain.E(i)->N( fnids[3] )->Coordinate(),
                                                  model_domain.E(i)->N( fnids[2] )->Coordinate(),
                                                  model_domain.E(i)->N( fnids[1] )->Coordinate(),
                                                  model_domain.E(i)->N( fnids[0] )->Coordinate() );
                  
                // verifying alignment of the element's unit normal with that of the model boundary
                model_domain.E(i)->UnitNormalToFace( model_domain.PerimeterFace(i,j), eUnitNormal );
                // checking whether the normals are aligned and of of same unit magnitude
                if ( flag == LEFT )  {
                     const double dotProduct(vector_product<3U,double>(leftNormal,eUnitNormal));
                     // testing for alignment
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at LEFT boundary");
                     // testing for unit length
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nLEFT: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
                else if ( flag == RIGHT )  {
                     const double dotProduct(vector_product<3U,double>(rightNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at RIGHT boundary");
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nRIGHT: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
                else if ( flag == BOTTOM )  {
                     const double dotProduct(vector_product<3U,double>(bottomNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at BOTTOM boundary");
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nBOTTOM: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
                else if ( flag == TOP )  {
                     const double dotProduct(vector_product<3U,double>(topNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at TOP boundary");
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nTOP: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
                else if ( flag == BACK )  {
                     const double dotProduct(vector_product<3U,double>(backNormal,eUnitNormal));
                     assert( backNormal[0]*eUnitNormal[0] + backNormal[1]*eUnitNormal[1] + backNormal[2]*eUnitNormal[2] == dotProduct );
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at BACK boundary");
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nBACK: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
                else if ( flag == FRONT )  {
                     const double dotProduct(vector_product<3U,double>(frontNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at FRONT boundary");
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nFRONT: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
                else if ( flag == IRREGULAR )  {
                     const double dotProduct(vector_product<3U,double>(frontNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     //_fail("negative dot product for normal at FRONT boundary");
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 10. );
                      // debugging
                      if ( dotProduct < 0. ) {
                           cout <<"\nIRREGULAR: "<< parseFiniteElementType( model_domain.E(i)->FE_Type() );
                           cout <<", face: "<< j <<", dotproduct: "<< dotProduct << endl;
                           out( eUnitNormal );
                        }
                  }
            }
      }

    _equal( model_surface_area, accumulated_area, 10. ); // 6-faces with 9m2, tolerance=10 eps

   
    // 2. testing whether the unit normals of the faces making up the outside boundaries of the model
    //    are outward pointing and aligned
    // -----------------------------------
    Boundary<3U>  left(model.Boundary("LEFT"));
    for ( auto it=left.CellsBegin(); it!=left.CellsEnd(); ++it ) {
         if ( (*it)->IsSurface() ) {
             (*it)->UnitNormal( eUnitNormal );
             const double dotProduct(vector_product<3U,double>(leftNormal,eUnitNormal));
             // testing for alignment
             _test( dotProduct > 0. );
             // testing for unit length
             _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
          }
      }
   
    Boundary<3U>  right(model.Boundary("RIGHT"));
    for ( auto it=right.CellsBegin(); it!=right.CellsEnd(); ++it ) {
         if ( (*it)->IsSurface() ) {
             (*it)->UnitNormal( eUnitNormal );
             const double dotProduct(vector_product<3U,double>(rightNormal,eUnitNormal));
             // testing for alignment
             _test( dotProduct > 0. );
             // testing for unit length
             _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
          }
      }
   
    Boundary<3U>  bottom(model.Boundary("BOTTOM"));
    for ( auto it=bottom.CellsBegin(); it!=bottom.CellsEnd(); ++it ) {
         if ( (*it)->IsSurface() ) {
             (*it)->UnitNormal( eUnitNormal );
             const double dotProduct(vector_product<3U,double>(bottomNormal,eUnitNormal));
             // testing for alignment
             _test( dotProduct > 0. );
             // testing for unit length
             _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
          }
      }
   
    Boundary<3U>  top(model.Boundary("TOP"));
    for ( auto it=top.CellsBegin(); it!=top.CellsEnd(); ++it ) {
         if ( (*it)->IsSurface() ) {
             (*it)->UnitNormal( eUnitNormal );
             const double dotProduct(vector_product<3U,double>(topNormal,eUnitNormal));
             // testing for alignment
             _test( dotProduct > 0. );
             // testing for unit length
             _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
          }
      }
   
    Boundary<3U>  back(model.Boundary("BACK"));
    for ( auto it=back.CellsBegin(); it!=back.CellsEnd(); ++it ) {
         if ( (*it)->IsSurface() ) {
             (*it)->UnitNormal( eUnitNormal );
             const double dotProduct(vector_product<3U,double>(backNormal,eUnitNormal));
             // testing for alignment
             _test( dotProduct > 0. );
             // testing for unit length
             _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
          }
      }
   
    Boundary<3U>  front(model.Boundary("FRONT"));
    for ( auto it=front.CellsBegin(); it!=front.CellsEnd(); ++it ) {
         if ( (*it)->IsSurface() ) {
             (*it)->UnitNormal( eUnitNormal );
             const double dotProduct(vector_product<3U,double>(frontNormal,eUnitNormal));
             // testing for alignment
             _test( dotProduct > 0. );
             // testing for unit length
             _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
          }
      }
   
 } // end TestWhetherSimplexNormalsAreOutwardPointing

  


bool Box_Test::TestWhetherSideBoundaryFlagsArePresent()
 {
    VSet<3U>  vset;
    // 64 elements and including nodes on the side of the model
    create_Pyramid_Hexa_VSet( vset );
    Model<3U> model( vset, "CSMP-variables.txt" );

    return hasAllSideBoundaries( model );
}



/** 
    All flags for 3D model including edges and corners
*/
bool Box_Test::TestWhetherAllBoxFlagsArePresent()
 {
    // const bool irregular_mesh(true), binary_file(true);
    // ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file );
    
    VSet<3U>  vset;
    const bool bSkewed{false};
    create_Prism_Hexa_VSet( vset, bSkewed );
    Model<3U>  model( vset, "CSMP-variables.txt" );

    return isStrictlyBoxShaped( model );
}
  
  

// HELPER
/**
 * @brief Writes a vector of 3D points to a legacy VTK file format.
 * * The legacy VTK format is human-readable and ideal for simple unstructured data.
 * * @param filename The name of the output VTK file (e.g., "output.vtk").
 * @param points The vector containing the 3D Point structures.
 * @return true if the file was successfully written, false otherwise.
 */
static bool writePointVectorToVTK(const std::string& filename, const std::vector<Point<3>>& points )
 {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    // 1. VTK File Header
    // The header is mandatory for all legacy VTK files.
    outfile << "# vtk DataFile Version 3.0\n";
    outfile << "Point Set Example\n";
    outfile << "ASCII\n"; // Specify ASCII for human-readable format

    // 2. DATASET Structure
    // Defines the type of data structure, here UNSTRUCTURED_GRID is suitable
    // for a collection of disconnected points, but POLYDATA is more common
    // and simpler for just points. We will use POLYDATA for simplicity.
    outfile << "DATASET POLYDATA\n";

    // 3. POINTS Section
    // Define the number of points and the data type (float or double)
    outfile << "POINTS " << points.size() << " double\n";

    // Write all point coordinates
    for (const auto& pt : points) {
        // Use a space-separated format for coordinates
        outfile << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    // 4. VERTICES Section (Optional, but recommended for visualization)
    // This section tells visualization software (like Paraview) that the points
    // should be rendered as individual vertices (dots).
    
    // Format: VERTICES N_cells Cell_list_size
    // For N points, N_cells is N. Cell_list_size is N + N (index + size_of_list)
    outfile << "VERTICES " << points.size() << " " << points.size() * 2 << "\n";

    // Write the connectivity list: each "cell" is just one vertex
    for (size_t i = 0; i < points.size(); ++i) {
        // Format: <Number of points in cell> <Index of point 1>
        // Here, 1 is the point count, and i is the zero-based index of the point
        outfile << "1 " << i << "\n";
    }

    // 5. POINT DATA (Optional attributes associated with each point)
    // If you had scalar or vector data (e.g., temperature, velocity) for each point, 
    // you would write it here. Example:
    /*
    outfile << "POINT_DATA " << points.size() << "\n";
    outfile << "SCALARS PointTemperatures float 1\n";
    outfile << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < points.size(); ++i) {
        outfile << (10.0 + i) << "\n"; // Example scalar data
    }
    */

    outfile.close();
    return true;
    
} // end writePointVectorToVTK



/**
       Returns map of points representing nodes in the Model with the chosen boundary flags.
*/
map<string,vector<Point<3>>> Box_Test::FormPointCloudsFromBOX_BOUNDARY_Flags( Model<3>& model )
  {
     map<BOX_BOUNDARY,vector<Point<3>>>  box_boundary_node_points;
     const Region<3>&                    model_domain = model.Region("Model");
     model_domain.RenumberNodes();
     
     // organising nodes with the same flag into vectors
     for ( const auto& nit : model_domain.NodeVector() ) {
        auto bbit = box_boundary_node_points.insert( make_pair( nit->AtBoundary(), vector<Point<3>>{nit->Coordinate()} ) );
        // if this flag already exists no new insertion was made
        if ( !bbit.second ) (*bbit.first).second.push_back( nit->Coordinate() );
     }
      
     // creating the VTK files from the point clouds
     map<string,vector<Point<3>>>  point_clouds_created;
     for ( auto& node_vec : box_boundary_node_points ) {
          // writing VTK files
          string new_region = parseBoundary( node_vec.first );
          bool success = writePointVectorToVTK( new_region + "_point_cloud.vtk", node_vec.second );
          if ( !success || node_vec.second.size() < 1 ) {
               cerr <<"\n"<<"Box_Test::PointCloudFromBOX_BOUNDARY_Flags: no cells in region '"<< new_region <<"'";
            }
          else point_clouds_created.insert( make_pair(new_region,node_vec.second) );
       }

     if ( verbose_ && !point_clouds_created.empty() ) {
        cout <<"\n"<<"Box_Test::PointCloudFromBOX_BOUNDARY_Flags: formed new model regions: "<< endl;
        for ( const auto& region : point_clouds_created )
          cout <<"  "<< region.first;
        cout << endl << endl;
     }

     return point_clouds_created;

  } // end FormRegionsFromBOX_BOUNDARY_Flags
  
  
  
  

bool Box_Test::TestBoundaryVersusBOX_BOUNDARY_Flagging()
 {
    VSet<3U>      vset;
    ModelTopology topo = create_FracBox( vset ); // without boundary flags
    
    // adding 'node number' as a variable
    PropertyData node_nums( NODE, SCALAR, 3U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );

    // renaming the boundaries according to a box-shaped model
    // (see assignments by normal orientation shown further below
    topo.ChangeDomainName( "BOUNDARY1", "LEFT");
    topo.ChangeDomainName( "BOUNDARY2", "RIGHT" );
    topo.ChangeDomainName( "BOUNDARY3", "FRONT" );
    topo.ChangeDomainName( "BOUNDARY4", "BACK" );
    topo.ChangeDomainName( "BOUNDARY5", "TOP" );
    topo.ChangeDomainName( "BOUNDARY6", "BOTTOM" );
    // reassigning BOX_BOUNDARY flags
    topo.AssignBoxShapedModelFlags( vset );
    
    // change the name so that the regions file is not found (and all regions are used)
    topo.ModelName("FracBox_without_regions_file");
    
    // building model with boundaries, converting surface elements to faces
    const bool create_boundaries_from_surf_elmts{ true };
    Model<3U> model( topo, vset, "CSMP-variables.txt", create_boundaries_from_surf_elmts );
    if ( verbose_ ) printBoxBoundaryFlags( model );

    _test( model.Mesh().Elements() + model.Mesh().Faces() == vset.Cells() );
    _test( model.Mesh().Nodes() == vset.Vertices() );
    
    // creating Face numbers and normals
    {
      const csmp::Index face_key = model.Database().StorageKey("face number");
      const csmp::Index nrml_key = model.Database().StorageKey("face vector");
      size_t n_face{0u};
      for ( auto fit=model.Mesh().FacesBegin(); fit!=model.Mesh().FacesEnd(); ++fit, ++n_face ) {
           // scalar (face number)
           (*fit).Store( face_key, makeScalar(ANY,n_face) );
           // vector (unit normal)
           auto unrml = (*fit).UnitNormal();
           VectorVariable<3> vc;
           for ( uint32_t i{0}; i<3; ++i ) vc(i) = unrml[i];
           (*fit).Store( nrml_key, vc );
       }
    }
    // saving node flags to 'nodal variable'
    {
      const csmp::Index nvar_key = model.Database().StorageKey("nodal variable");
      Region<3>& model_domain = model.Region("Model");
      for ( auto& nit : model_domain.NodeVector() )
        nit->Store( nvar_key, makeScalar(ANY, static_cast<double>(nit->AtBoundary()) ) );
    }

    if ( verbose_ ) {
         const bool vtk_output{ true };
         if ( vtk_output ) { // tested: 5/12/2025: FracBox normals are all outward-pointing
             VTU_Interface<3U> vtu_output( model ); // checked: SKM 6/7/24 (these are the correct boundaries)
             vtu_output.OutputDataToVTU( "BoxTest_BACK_fn",   "face vector", model.Boundary("BACK"), 0 );
             vtu_output.OutputDataToVTU( "BoxTest_BOTTOM_fn", "face vector", model.Boundary("BOTTOM"), 0 );
             vtu_output.OutputDataToVTU( "BoxTest_RIGHT_fn",  "face vector", model.Boundary("RIGHT"), 0 );
             vtu_output.OutputDataToVTU( "BoxTest_TOP_fn",    "face vector", model.Boundary("TOP"), 0 );
             vtu_output.OutputDataToVTU( "BoxTest_LEFT_fn",   "face vector", model.Boundary("LEFT"), 0 );
             vtu_output.OutputDataToVTU( "BoxTest_FRONT_fn",  "face vector", model.Boundary("FRONT"), 0 );
             // boundary flags turned into regions
             FormPointCloudsFromBOX_BOUNDARY_Flags( model );
          }
         // printing the normal to first Face
         cout <<"\nBoxTest: boundary flag recreation:\n";
         cout <<"\nBACK "<< model.Boundary("BACK").E(0)->UnitNormal();
         cout <<"\nBOTTOM "<< model.Boundary("BOTTOM").E(0)->UnitNormal();
         cout <<"\nRIGHT "<< model.Boundary("RIGHT").E(0)->UnitNormal();
         cout <<"\nTOP "<< model.Boundary("TOP").E(0)->UnitNormal();
         cout <<"\nLEFT "<< model.Boundary("LEFT").E(0)->UnitNormal();
         cout <<"\nFRONT "<< model.Boundary("FRONT").E(0)->UnitNormal();
         cout << endl;
      }
   
    // testing that the normals are correct
    {
      auto ba_n = model.Boundary("BACK").E(0)->UnitNormal();
      auto bo_n = model.Boundary("BOTTOM").E(0)->UnitNormal();
      auto ri_n = model.Boundary("RIGHT").E(0)->UnitNormal();
      auto to_n = model.Boundary("TOP").E(0)->UnitNormal();
      auto le_n = model.Boundary("LEFT").E(0)->UnitNormal();
      auto fr_n = model.Boundary("FRONT").E(0)->UnitNormal();
      _test( approximatelyEqual(le_n[0],-1) && approximatelyEqual(le_n[1],0) && approximatelyEqual(le_n[2],0) );
      _test( approximatelyEqual(ri_n[0],1) && approximatelyEqual(ri_n[1],0) && approximatelyEqual(ri_n[2],0) );
      _test( approximatelyEqual(bo_n[0],0) && approximatelyEqual(bo_n[1],-1) && approximatelyEqual(bo_n[2],0) );
      _test( approximatelyEqual(to_n[0],0) && approximatelyEqual(to_n[1],1) && approximatelyEqual(to_n[2],0) );
      _test( approximatelyEqual(ba_n[0],0) && approximatelyEqual(ba_n[1],0) && approximatelyEqual(ba_n[2],-1) );
      _test( approximatelyEqual(fr_n[0],0) && approximatelyEqual(fr_n[1],0) && approximatelyEqual(fr_n[2],1) );
    }

   // checking that the BoxBoundary flags match the boundary names
   {
      Boundary<3U>& left   = model.Boundary("LEFT");
      for ( const auto& nit : left.NodeVector() )
       _test( isLEFT( nit->AtBoundary() ) && isLEFT( left.AtBoundary() ) );
      Boundary<3U>& right  = model.Boundary("RIGHT");
      for ( const auto& nit : right.NodeVector() )
       _test( isRIGHT( nit->AtBoundary() ) && isRIGHT( right.AtBoundary() ) );
      Boundary<3U>& bottom = model.Boundary("BOTTOM");
      for ( const auto& nit : bottom.NodeVector() )
       _test( isBOTTOM( nit->AtBoundary() ) && isBOTTOM( bottom.AtBoundary() ) );
      Boundary<3U>& top    = model.Boundary("TOP");
      for ( const auto& nit : top.NodeVector() )
       _test( isTOP( nit->AtBoundary() ) && isTOP( top.AtBoundary() ) );
      Boundary<3U>& back   = model.Boundary("BACK");
      for ( const auto& nit : back.NodeVector() )
       _test( isBACK( nit->AtBoundary() ) && isBACK( back.AtBoundary() ) );
      Boundary<3U>& front  = model.Boundary("FRONT");
      for ( const auto& nit : front.NodeVector() )
       _test( isFRONT( nit->AtBoundary() ) && isFRONT( front.AtBoundary() ) );
   }
   
   return true;
   
} // end Boundary vs BOX_BOUNDARY flagging
  
  
  
/** Edges (intersections of box-side pairs):
 
    E1:  x-axis, y=0, z=0  -> BOTTOM & BACK
    E2:  y-axis, x=3, z=0  -> RIGHT & BACK
    E3:  x-axis, y=3, z=0  -> TOP   & BACK
    E4:  y-axis, x=0, z=0  -> LEFT  & BACK
    E5:  z-axis, x=0, y=0  -> LEFT  & BOTTOM
    E6:  z-axis, x=3, y=0  -> RIGHT & BOTTOM
    E7:  z-axis, x=3, y=3  -> RIGHT & TOP
    E8:  z-axis, x=0, y=3  -> LEFT  & TOP
    E9:  x-axis, y=0, z=3  -> BOTTOM & FRONT
    E10: y-axis, x=3, z=3  -> RIGHT  & FRONT
    E11: x-axis, y=3, z=3  -> TOP    & FRONT
    E12: y-axis, x=0, z=3  -> LEFT   & FRONT
    
     SKM checked
*/
constexpr std::array<std::array<BOX_BOUNDARY, 2>, 12> sidesOfEdge() noexcept {
    return {{
        { static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE), static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) },  // EDGE1
        { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),  static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) },  // EDGE2
        { static_cast<BOX_BOUNDARY>(TOP_OUTSIDE),    static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) },  // EDGE3
        { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),   static_cast<BOX_BOUNDARY>(BACK_OUTSIDE) },  // EDGE4
        
        { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),   static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE) }, // EDGE5
        { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),  static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE) }, // EDGE6
        { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),  static_cast<BOX_BOUNDARY>(TOP_OUTSIDE) },    // EDGE7
        { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),   static_cast<BOX_BOUNDARY>(TOP_OUTSIDE) },    // EDGE8
        
        { static_cast<BOX_BOUNDARY>(BOTTOM_OUTSIDE), static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) },  // EDGE9
        { static_cast<BOX_BOUNDARY>(RIGHT_OUTSIDE),  static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) },  // EDGE10
        { static_cast<BOX_BOUNDARY>(TOP_OUTSIDE),    static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) },  // EDGE11
        { static_cast<BOX_BOUNDARY>(LEFT_OUTSIDE),   static_cast<BOX_BOUNDARY>(FRONT_OUTSIDE) }   // EDGE12
    }};
}
 
  
  

/**
   tests whether method recreateBoxBoundaryFlags() which  allows to recreate boundary flags works correctly
   
   @note also tests EstablishBoxBoundariesFromNodeFlags()
*/
bool Box_Test::TestBoundaryFlagRecreation()
 {
    VSet<3U>   vset;
    const bool bSkewed{false};
    create_Pyramid_Hexa_VSet( vset, bSkewed ); // used because boundary flags have been verified
    
    // adding 'node number' as a variable
    PropertyData node_nums( NODE, SCALAR, 3U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );

    // building model (usin constructor that does not attempt to recreate boundary flags)
    Model<3U> model( vset, "CSMP-variables.txt" );
    model.Name("Pyramid_Hexa_without_regions_file");
    _test( model.Mesh().Elements() + model.Mesh().Faces() == vset.Cells() );
    _test( model.Mesh().Nodes() == vset.Vertices() );
    const bool recreate_box_boundary_flags_before{false}; // do not change any node flags
    model.EstablishBoxBoundariesFromNodeFlags( recreate_box_boundary_flags_before );
    _test( isStrictlyBoxShaped( model ) );

    // testing whether the reflagging works correctly
    recreateBoxBoundaryFlags( model );
    
    const Region<3U>& model_domain = model.Region("Model");
    model_domain.RenumberNodes();
    const csmp::Index nn_key = model.Database().StorageKey("node number");
    // testing that the boundary flags match
    for ( size_t i{0U}; i<vset.Vertices(); i++ ) { // converting double value into integer
         size_t       node_number = static_cast<size_t>(model_domain.N(i)->Read(nn_key));
         BOX_BOUNDARY nodal_bflag = model_domain.N(i)->AtBoundary();
         // IMPORTANT: because the nodes in regions are sorted by pointer values, 'node number' is needed to find the corresponding bflag in VSet
         BOX_BOUNDARY vset_bflag  = static_cast<BOX_BOUNDARY>(vset.BFlag(node_number));
         _test( vset_bflag == nodal_bflag );
         if ( verbose_ && vset_bflag != nodal_bflag ) {
              cout <<"\nnode "<< node_number <<", idx "<< model_domain.N(node_number)->Idx() <<": "<< parseBoundary( vset_bflag );
              cout <<"(vset) vs. "<< parseBoundary( nodal_bflag ) <<"(rebuilt)"<< endl;
           }
      }

    // repeat on model with recreated boundary conditions
    return isStrictlyBoxShaped( model );
    
 } // end TestBoundaryFlagRecreation




/** 
    writes model to CSMP binary and then re-reads it to see whether the box flags survive
*/
bool Box_Test::TestWhetherBoundaryFlagsArePreservedInBinaryFile()
 {
    VSet<3U>  vset;
    const bool bSkewed{false};
    create_Prism_Hexa_VSet( vset, bSkewed );
    Model<3U>  model( vset, "CSMP-variables.txt" );
    
    const Region<3>&  mregion(model.Region("Model"));
   
    // storing the flags in node and element order in a list for comparison
    list<BOX_BOUNDARY>  node_flags_before;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      node_flags_before.push_back( (*nit)->AtBoundary() );
    list<BOX_BOUNDARY>  elmt_flags_before;
    for ( auto eit=mregion.CellsBegin(); eit!=mregion.CellsEnd(); eit++ )
      for ( uint32_t i{0u}; i<(*eit)->Neighbors(); ++i )
        elmt_flags_before.push_back( (*eit)->AtBoundary(i) );
   
    // Saving the model to disk
    const string output_model("box_test_temporary");
    model.OutputToBinaryFile( output_model.c_str() );
   
    // Recovering the model from file (all variables are read as subset is empty)
    Model<3U>* mptr = new Model<3>( output_model, set<string>({}) );

    // recovering the flags again for comparison
    list<BOX_BOUNDARY>  node_flags_after;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      node_flags_after.push_back( (*nit)->AtBoundary() );
    list<BOX_BOUNDARY>  elmt_flags_after;
    for ( auto eit=mregion.CellsBegin(); eit!=mregion.CellsEnd(); eit++ )
      for ( uint32_t i{0u}; i<(*eit)->Neighbors(); ++i )
        elmt_flags_after.push_back( (*eit)->AtBoundary(i) );

    delete mptr;
   
    // comparing node and element flags of the original and the restored model
    if ( node_flags_after != node_flags_before ) return false;
    if ( elmt_flags_after != elmt_flags_before ) return false;
    return true;
 }


/**
     Tests new output functionality.
*/
bool Box_Test::TestWhetherBoundaryFlagsArePreservedInBinaryFile1()
 {
    VSet<3U>  vset;
    const bool bSkewed{false};
    create_Prism_Hexa_VSet( vset, bSkewed );
    Model<3U>  model( vset, "CSMP-variables.txt" );
    const Region<3>&  mregion(model.Region("Model"));

    // storing the flags in node and element order in a list for comparison
    list<BOX_BOUNDARY>  node_flags_before;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      node_flags_before.push_back( (*nit)->AtBoundary() );
    list<BOX_BOUNDARY>  elmt_flags_before;
    for ( auto eit=mregion.CellsBegin(); eit!=mregion.CellsEnd(); eit++ )
      for ( uint32_t i{0u}; i<(*eit)->Neighbors(); ++i )
        elmt_flags_before.push_back( (*eit)->AtBoundary(i) );
   
    // Saving the model to disk
    const string output_model("box_test_temporary");
    // NEW model.OutputToDisk( output_model.c_str() );
    model.OutputToBinaryFile("box_test_temporary");
   
    // Recovering the model from file (read all variables)
    Model<3U>* mptr = new Model<3>( output_model, set<string>({}) );

    // recovering the flags again for comparison
    list<BOX_BOUNDARY>  node_flags_after;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      node_flags_after.push_back( (*nit)->AtBoundary() );
    list<BOX_BOUNDARY>  elmt_flags_after;
    for ( auto eit=mregion.CellsBegin(); eit!=mregion.CellsEnd(); eit++ )
      for ( uint32_t i{0u}; i<(*eit)->Neighbors(); ++i )
        elmt_flags_after.push_back( (*eit)->AtBoundary(i) );

    delete mptr;
   
    // comparing node and element flags of the original and the restored model
    if ( node_flags_after != node_flags_before ) return false;
    if ( elmt_flags_after != elmt_flags_before ) return false;
    return true;
 }



/** How to use the output
 
 This method generates a legacy VTK unstructured grid output file for visualisation with Paraview or other.
 
 @code
 int main() {
    const std::string filename = "prism_face_normals.vtk";
    if (write_prism_normals_to_vtk(filename)) {
        std::cout << "Successfully wrote the prism faces and normals to " << filename << std::endl;
        std::cout << "\nTo visualize the normals in Paraview:\n";
        std::cout << "1. Load the file.\n";
        std::cout << "2. Apply the 'Glyph' filter.\n";
        std::cout << "3. In Glyph Properties:\n";
        std::cout << "   - Set 'Orientation Array' to 'FaceNormals'.\n";
        std::cout << "   - Set 'Source' to 'Cell Centers'.\n";
        std::cout << "   - Adjust 'Scale Factor' (e.g., 0.1) for visibility.\n";
        std::cout << "4. Apply the Glyph filter.\n";
    } else {
        std::cerr << "Failed to write VTK file." << std::endl;
    }
#endcode

*/
bool writePrismWithNormalsToVTK( const csmp::Element<3>* eptr, const std::string& filename )
{
    assert( eptr != nullptr );
    if ( eptr->FE_Type() != ISOPARAMETRIC_LINEAR_PRISM ) return false;
    
    std::ofstream outfile(filename + ".vtk");

    if (!outfile.is_open()) {
        std::cerr << "writePrismWithNormalsToVTK: Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    // Set precision for floating-point output
    outfile << std::fixed << std::setprecision(6);

    // --- 1. Define the Prism's 6 Vertices (Points) ---
    // A standard prism aligned with the Z-axis (height 1, side 1).
    const std::vector<Point<3>> points = { eptr->N(0)->Coordinate(), eptr->N(1)->Coordinate(), eptr->N(2)->Coordinate(), // Base Triangle (0, 1, 2)
                                           eptr->N(3)->Coordinate(), eptr->N(4)->Coordinate(), eptr->N(5)->Coordinate() }; // Top Triangle (3, 4, 5)

    const size_t num_points = points.size();

    // --- 2. Define the 5 Faces (Decomposed Cells) using Point Indices ---
    // The faces are defined in counter-clockwise order for correct normal visualization.
    // Face 0: Quad (bottom)
    // Face 1: Quad (front)
    // Face 2: Quad (side)
    // Face 3: Triangle (back-left)
    // Face 4: Triangle (top)
    
    // Connectivity list: (Count, Index0, Index1, ...)
    std::vector<int> connectivity; 
    std::vector<int> face_counts;
    
    // F0: Bottom Triangle (P0, P2, P1)
    connectivity.insert(connectivity.end(), {3, 0, 2, 1});
    face_counts.push_back(3);

    // F1: Quad (P0, P1, P4, P3)
    connectivity.insert(connectivity.end(), {4, 0, 1, 4, 3});
    face_counts.push_back(4);
    
    // F2: Quad (P1, P2, P5, P4)
    connectivity.insert(connectivity.end(), {4, 1, 2, 5, 4});
    face_counts.push_back(4);

    // F3: Quad (P0, P3, P5, P2)
    connectivity.insert(connectivity.end(), {4, 0, 3, 5, 2});
    face_counts.push_back(4);

    // F4: Top Triangle (P3, P4, P5) - Note: Normal direction is opposite of F3, hence the index swap.
    connectivity.insert(connectivity.end(), {3, 3, 4, 5});
    face_counts.push_back(3);

    const size_t num_faces = face_counts.size();

    // --- 3. Define the Unit Normal Vector for Each Face ---
    eptr->CoordinateMatrix();
    const std::vector<Point<3>> face_normals = {
        eptr->UnitNormalToFace(0),
        eptr->UnitNormalToFace(1), eptr->UnitNormalToFace(2), eptr->UnitNormalToFace(3),
        eptr->UnitNormalToFace(4)
    };

    // --- 4. Write VTK File Header ---
    outfile << "# vtk DataFile Version 3.0\n";
    outfile << "VTK Prism Faces with Normals\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n"; // POLYDATA is suitable for surface meshes

    // --- 5. POINTS Block ---
    outfile << "POINTS " << num_points << " double\n";
    for (const auto& pt : points) {
        outfile << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    // --- 6. POLYGONS (Faces/Cells) Block ---
    // Total size is sum of all face counts + number of faces (to store the count value).
    int list_size = 0;
    for( int count : face_counts) list_size += (count + 1);

    outfile << "POLYGONS " << num_faces << " " << list_size << "\n";
    
    int conn_index = 0;
    for ( unsigned int i = 0; i < num_faces; ++i) {
        // First value is the number of points in the polygon (3 for triangle, 4 for quad)
        outfile << face_counts[i]; 
        for ( int j = 0; j < face_counts[i]; ++j) {
            outfile << " " << connectivity[ static_cast<size_t>(conn_index + j + 1) ];
        }
        outfile << "\n";
        conn_index += face_counts[i] + 1; // Move to start of next face
    }

    // --- 7. CELL DATA Block (Normals) ---
    // Attach the normal vector as VECTORS to each face (cell)
    outfile << "CELL_DATA " << num_faces << "\n";
    outfile << "VECTORS FaceNormals double\n";

    for (const auto& normal : face_normals) {
        outfile << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    }

    outfile.close();
    return true;
    
} // end writePrismWithNormalsToVTK



bool writePyramidWithNormalsToVTK( const csmp::Element<3>* eptr, const std::string& filename )
{
    assert( eptr != nullptr );
    if ( eptr->FE_Type() != ISOPARAMETRIC_LINEAR_PYRAMID ) return false;
    
    std::ofstream outfile(filename + ".vtk");

    if (!outfile.is_open()) {
        std::cerr << "writePrismWithNormalsToVTK: Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    // Set precision for floating-point output
    outfile << std::fixed << std::setprecision(6);

    // --- 1. Define the Pyramids 5 Vertices (Points) ---
    // A standard pyramid aligned with the Z-axis (height 1, side 1).
    const std::vector<Point<3>> points = { eptr->N(0)->Coordinate(), eptr->N(1)->Coordinate(), eptr->N(2)->Coordinate(),
                                           eptr->N(3)->Coordinate(), eptr->N(4)->Coordinate() };

    const size_t num_points = points.size();

    // --- 2. Define the 5 Faces (Decomposed Cells) using Point Indices ---
    // The faces are defined in counter-clockwise order for correct normal visualization.
    // Face 0: Triangle (bottom)
    // Face 1: Triangle (front)
    // Face 2: Triangle (side)
    // Face 3: Triangle (back-left)
    // Face 4: Quadrilateral (bottom)
    
    // Connectivity list: (Count, Index0, Index1, ...)
    std::vector<int> connectivity; 
    std::vector<int> face_counts;
    
    // F0: Bottom Triangle (P0, P2, P1)
    connectivity.insert(connectivity.end(), {3, 1,0,4});
    face_counts.push_back(3);

    // F1: Quad (P1, P2, P5, P4)
    connectivity.insert(connectivity.end(), {3, 1,2,4});
    face_counts.push_back(3);

    // F2: Quad (P0, P3, P5, P2)
    connectivity.insert(connectivity.end(), {3, 2,3,4});
    face_counts.push_back(3);

    // F3: Quad (P0, P3, P5, P2)
    connectivity.insert(connectivity.end(), {3, 0,4,3});
    face_counts.push_back(3);

    // F4: Base
    connectivity.insert(connectivity.end(), {4, 0,3,2,1 });
    face_counts.push_back(4);

    const size_t num_faces = face_counts.size();

    // --- 3. Define the Unit Normal Vector for Each Face ---
    eptr->CoordinateMatrix();
    const std::vector<Point<3>> face_normals = {
        eptr->UnitNormalToFace(0),
        eptr->UnitNormalToFace(1),
        eptr->UnitNormalToFace(2),
        eptr->UnitNormalToFace(3),
        eptr->UnitNormalToFace(4)
    };

    // --- 4. Write VTK File Header ---
    outfile << "# vtk DataFile Version 3.0\n";
    outfile << "VTK Pyramid Faces with Normals\n";
    outfile << "ASCII\n";
    outfile << "DATASET POLYDATA\n"; // POLYDATA is suitable for surface meshes

    // --- 5. POINTS Block ---
    outfile << "POINTS " << num_points << " double\n";
    for (const auto& pt : points) {
        outfile << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    // --- 6. POLYGONS (Faces/Cells) Block ---
    // Total size is sum of all face counts + number of faces (to store the count value).
    int list_size = 0;
    for( int count : face_counts) list_size += (count + 1);

    outfile << "POLYGONS " << num_faces << " " << list_size << "\n";
    
    int conn_index = 0;
    for ( unsigned int i = 0; i < num_faces; ++i) {
        // First value is the number of points in the polygon (3 for triangle, 4 for quad)
        outfile << face_counts[i]; 
        for ( int j = 0; j < face_counts[i]; ++j) {
            outfile << " " << connectivity[ static_cast<size_t>(conn_index + j + 1) ];
        }
        outfile << "\n";
        conn_index += face_counts[i] + 1; // Move to start of next face
    }

    // --- 7. CELL DATA Block (Normals) ---
    // Attach the normal vector as VECTORS to each face (cell)
    outfile << "CELL_DATA " << num_faces << "\n";
    outfile << "VECTORS FaceNormals double\n";

    for (const auto& normal : face_normals) {
        outfile << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    }

    outfile.close();
    return true;
    
} // end writePyramidithNormalsToVTK








} // csmp










