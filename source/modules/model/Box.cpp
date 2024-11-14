#include "Box.h"
#include "Region.h"
#include "Boundary.h"
#include "Element.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

/**
@addtogroup CSMPglobalFunctions
@{
*/

/// returns wether a node (or line element) lies on an edge of the model
bool isEdge( BOX_BOUNDARY bd )
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
bool isCorner( BOX_BOUNDARY bd )
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
bool isSide( BOX_BOUNDARY bd )
{
  if ( bd == RIGHT )  return true;
  if ( bd == LEFT )   return true;
  if ( bd == TOP )    return true;
  if ( bd == BOTTOM ) return true;
  if ( bd == FRONT )  return true;
  if ( bd == BACK )   return true;
  if ( bd == IRREGULAR ) return true;
  return false;
}


/// of rectangular (brick-shaped) model; @test OK
bool isLEFT( BOX_BOUNDARY bd )
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
bool isRIGHT( BOX_BOUNDARY bd )
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
bool isTOP( BOX_BOUNDARY bd )
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
bool isBOTTOM( BOX_BOUNDARY bd )
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
bool isFRONT( BOX_BOUNDARY bd )
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
bool isBACK( BOX_BOUNDARY bd )
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



/**
     For irregular shaped boundaries that were created from surfaces
     @todo figure out whether edges and corners should be included?
*/
bool canBeIRREGULAR( BOX_BOUNDARY bd )
{
  if ( bd == IRREGULAR ) return true;
  if ( isEdge( bd ) )  return true;
  if ( isCorner( bd ) )  return true;
  return false;
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
   Parses the BOX_BOUNDARY identifier (see Box.h). If the boundary flag cannot be resolved a value of NOT is returned if it is positive and IRREGULAR if negative.
   
   @attention use this only to convert ANSYS or outside-of CSMP generated flags of the enlisted names into BOX_BOUNDARY enums
*/
BOX_BOUNDARY intToBOX_BOUNDARY( int8_t i )
{
  if ( i == 0 )                 return NOT;
  if ( i == IRREGULAR_OUTSIDE ) return IRREGULAR;
  if ( i == LEFT_OUTSIDE )    return LEFT;
  if ( i == RIGHT_OUTSIDE )   return RIGHT;
  if ( i == BOTTOM_OUTSIDE )  return BOTTOM;
  if ( i == TOP_OUTSIDE )     return TOP;
  if ( i == FRONT_OUTSIDE )   return FRONT;
  if ( i == BACK_OUTSIDE )    return BACK;
  if ( i == CNR_MIN )         return CNR1;
  if ( i == CNR_MIN_MAXX )    return CNR2;
  if ( i == CNR_MAX_MAXX )    return CNR3;
  if ( i == CNR_MAX_MINXZ )   return CNR4;
  if ( i == CNR_MIN_MAXZ )    return CNR5;
  if ( i == CNR_MIN_MAXXZ )   return CNR6;
  if ( i == CNR_MAX )         return CNR7;
  if ( i == CNR_MAX_MAXZ )    return CNR8;
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
  if ( i == MULTIPLE_BOUNDARIES ) return MULTIPLE;
  if ( i < 0 ) return IRREGULAR;

  //cout <<"\nintToSG_BOUNDARY(int): unable to parse integer: "<< i << endl;
  return NOT;
}




/**
   Parses the BOX_BOUNDARY identifier (see Box.h). If the boundary flag cannot be resolved a value of NOT is returned if it is positive and IRREGULAR if negative.
*/
std::string  parseBoundary( BOX_BOUNDARY i )
{
  if ( i == NOT )      return string( "NOT" );
  if ( i == IRREGULAR )return string( "IRREGULAR" );
  if ( i == LEFT )     return string( "LEFT" );
  if ( i == RIGHT )    return string( "RIGHT" );
  if ( i == BOTTOM )   return string( "BOTTOM" );
  if ( i == TOP )      return string( "TOP" );
  if ( i == FRONT )    return string( "FRONT" );
  if ( i == BACK )     return string( "BACK" );
  if ( i == CNR1 )     return string( "CNR1" );
  if ( i == CNR2 )     return string( "CNR2" );
  if ( i == CNR3 )     return string( "CNR3" );
  if ( i == CNR4 )     return string( "CNR4" );
  if ( i == CNR5 )     return string( "CNR5" );
  if ( i == CNR6 )     return string( "CNR6" );
  if ( i == CNR7 )     return string( "CNR7" );
  if ( i == CNR8 )     return string( "CNR8" );
  if ( i == EDGE1 )    return string( "EDGE1" );
  if ( i == EDGE2 )    return string( "EDGE2" );
  if ( i == EDGE3 )    return string( "EDGE3" );
  if ( i == EDGE4 )    return string( "EDGE4" );
  if ( i == EDGE5 )    return string( "EDGE5" );
  if ( i == EDGE6 )    return string( "EDGE6" );
  if ( i == EDGE7 )    return string( "EDGE7" );
  if ( i == EDGE8 )    return string( "EDGE8" );
  if ( i == EDGE9 )    return string( "EDGE9" );
  if ( i == EDGE10 )   return string( "EDGE10" );
  if ( i == EDGE11 )   return string( "EDGE11" );
  if ( i == EDGE12 )   return string( "EDGE12" );
  if ( i == INTERNAL ) return string( "INTERNAL" );
  if ( i == MULTIPLE ) return string( "MULTIPLE" );
  if ( i < 0 ) return string( "IRREGULAR" );

//  cerr << "\nparseBoundary(BOX_BOUNDARY): unable to parse BOX_BOUNDARY: " << i << endl;
  return string( "NOT" );
}



/// text to boundary enum
BOX_BOUNDARY  parseBoundary( const string& i )
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
  
  // if the name is not recognized it is not a BOX_BOUNDARY
  return MULTIPLE;
}





/**
Returns true if @param bd belongs to the corresponding edge.
Example: CNR1 and CNR2 belong to EDGE1.
*/
bool belongsToEdge( BOX_BOUNDARY edge, BOX_BOUNDARY bd )
{
  if ( edge == EDGE1 and (bd == EDGE1 or bd == CNR1 or bd == CNR2) ) return true;
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
bool belongsToSide( BOX_BOUNDARY side, BOX_BOUNDARY bd )
{
  if ( side == LEFT )   return isLEFT( bd );
  else if ( side == RIGHT )  return isRIGHT( bd );
  else if ( side == TOP )    return isTOP( bd );
  else if ( side == BOTTOM ) return isBOTTOM( bd );
  else if ( side == FRONT )  return isFRONT( bd );
  else if ( side == BACK )   return isBACK( bd );
  else if ( side == IRREGULAR && bd == IRREGULAR ) return true;
  return false;
}


/**
Sets supplied model min, max coordinates to the coordinate value
range that applies to the specified model boundary.
*/
void boundaryMinMaxCoordinates( BOX_BOUNDARY boundary, csmp::Point<1U>&, csmp::Point<1U>& )
{ assert( !isCorner( boundary ) ); }


void boundaryMinMaxCoordinates( BOX_BOUNDARY boundary,
                                csmp::Point<2U>& model_coord_min,
                                csmp::Point<2U>& model_coord_max )
{
  assert( !isCorner( boundary ) );
  if ( boundary == LEFT ) {
    model_coord_max[0] = model_coord_min[0];
    return;
  }
  if ( boundary == RIGHT ) {
    model_coord_min[0] = model_coord_max[0];
    return;
  }
  if ( boundary == TOP ) {
    model_coord_min[1] = model_coord_max[1];
    return;
  }
  if ( boundary == BOTTOM ) {
    model_coord_max[1] = model_coord_min[1];
    return;
  }
  std::cout << "\nboundaryMinMaxCoordinates(2D): boundary could not be parsed." << std::endl;

} // end boundaryMinMaxCoordinates




void boundaryMinMaxCoordinates( BOX_BOUNDARY boundary,
                                csmp::Point<3U>& model_coord_min,
                                csmp::Point<3U>& model_coord_max )
{
  assert( !isCorner( boundary ) );
  if ( boundary == LEFT ) {
    model_coord_max[0] = model_coord_min[0];
    return;
  }
  if ( boundary == RIGHT ) {
    model_coord_min[0] = model_coord_max[0];
    return;
  }
  if ( boundary == TOP ) {
    model_coord_min[1] = model_coord_max[1];
    return;
  }
  if ( boundary == BOTTOM ) {
    model_coord_max[1] = model_coord_min[1];
    return;
  }
  if ( boundary == FRONT ) {
    model_coord_min[2] = model_coord_max[2];
    return;
  }
  if ( boundary == BACK ) {
    model_coord_max[2] = model_coord_min[2];
    return;
  }

  // edges
  // 1-4 (in XY-plane)
  if ( boundary == EDGE1 ) { // x is variable
    model_coord_max[2] = model_coord_min[2];
    model_coord_max[1] = model_coord_min[1];
    return;
  }
  if ( boundary == EDGE2 ) { // y is variable
    model_coord_max[2] = model_coord_min[2];
    model_coord_min[0] = model_coord_max[0];
    return;
  }
  if ( boundary == EDGE3 ) { // x is variable
    model_coord_max[2] = model_coord_min[2];
    model_coord_min[1] = model_coord_max[1];
    return;
  }
  if ( boundary == EDGE4 ) { // y is variable
    model_coord_max[2] = model_coord_min[2];
    model_coord_max[0] = model_coord_min[0];
    return;
  }

  // 5-8 aligned with the Z axis
  if ( boundary == EDGE5 ) { // z is variable
    model_coord_max[0] = model_coord_min[0]; // min XY
    model_coord_max[1] = model_coord_min[1];
    return;
  }
  if ( boundary == EDGE6 ) { // z is variable
    model_coord_min[0] = model_coord_max[0]; // max X min Y
    model_coord_max[1] = model_coord_min[1];
    return;
  }
  if ( boundary == EDGE7 ) { // z is variable
    model_coord_min[0] = model_coord_max[0]; // max XY
    model_coord_min[1] = model_coord_max[1];
    return;
  }
  if ( boundary == EDGE8 ) { // z is variable
    model_coord_max[0] = model_coord_min[0]; // min X max Y
    model_coord_min[1] = model_coord_max[1];
    return;
  }

  // 9-12 (in XY-plane)
  if ( boundary == EDGE9 ) { // x is variable
    model_coord_min[2] = model_coord_max[2];
    model_coord_max[1] = model_coord_min[1];
    return;
  }
  if ( boundary == EDGE10 ) { // y is variable
    model_coord_min[2] = model_coord_max[2];
    model_coord_min[0] = model_coord_max[0];
    return;
  }
  if ( boundary == EDGE11 ) { // x is variable
    model_coord_min[2] = model_coord_max[2];
    model_coord_min[1] = model_coord_max[1];
    return;
  }
  if ( boundary == EDGE12 ) { // y is variable
    model_coord_min[2] = model_coord_max[2];
    model_coord_max[0] = model_coord_min[0];
    return;
  }
  std::cout << "\nboundaryMinMaxCoordinates(3D): boundary could not be parsed." << std::endl;

} // end boundaryMinMaxCoordinates (3D)


/**
Returns a normal (vector) to the specified boundary of a box-shaped model.
*/
void Box::UnitNormalTo( BOX_BOUNDARY bdry, uint32_t dim, vector<double>& nrml ) const
{
  nrml.resize( dim );

  // only two cases are possible
  if ( dim == 1U ) {
    if ( bdry == LEFT )  nrml[0] = -1.;
    else if ( bdry == RIGHT ) nrml[0] = 1.;
    else
      // 1D models only have x-coordinate and therefore only left and right
      throw csmp::Exception( ERROR, "Box::UnitNormalTo",
                             "1D models only have a LEFT and RIGHT boundary." );
  }

  else if ( dim == 2U ) {
    if ( bdry == LEFT ) {
      nrml[0] = -1.; nrml[1] = 0.;
    }
    else if ( bdry == RIGHT ) {
      nrml[0] = 1.; nrml[1] = 0.;
    }
    else if ( bdry == TOP ) {
      nrml[0] = 0.; nrml[1] = 1.;
    }
    else if ( bdry == BOTTOM ) {
      nrml[0] = 0.; nrml[1] = -1.;
    }
    else
      // 2D models only have x,y-coordinate and therefore only 4 boundaries
      throw csmp::Exception( ERROR, "Box::UnitNormalTo",
                             "2D models only have a LEFT, RIGHT, TOP & BOTTOM boundary." );
  }


  else if ( dim == 3U ) {
    if ( bdry == LEFT ) {
      nrml[0] = -1.; nrml[1] = 0.; nrml[2] = 0.;
    }
    else if ( bdry == RIGHT ) {
      nrml[0] = 1.; nrml[1] = 0.; nrml[2] = 0.;
    }
    else if ( bdry == TOP ) {
      nrml[0] = 0.; nrml[1] = 1.; nrml[2] = 0.;
    }
    else if ( bdry == BOTTOM ) {
      nrml[0] = 0.; nrml[1] = -1.; nrml[2] = 0.;
    }
    else if ( bdry == FRONT ) {
      nrml[0] = 0.; nrml[1] = 0.; nrml[2] = 1.;
    }
    else if ( bdry == BACK ) {
      nrml[0] = 0.; nrml[1] = 0.; nrml[2] = -1.;
    }
    // 12 edges are treated with 45o normals
    // back
    else if ( bdry == EDGE1 ) {
      nrml[0] = 0.; nrml[1] = sin( 45. ); nrml[2] = sin( 45. );
    }
    else if ( bdry == EDGE2 ) {
      nrml[0] = -sin( 45. ); nrml[1] = 0.; nrml[2] = sin( 45. );
    }
    else if ( bdry == EDGE3 ) {
      nrml[0] = 0.; nrml[1] = -sin( 45. ); nrml[2] = sin( 45. );
    }
    else if ( bdry == EDGE4 ) {
      nrml[0] = sin( 45. ); nrml[1] = 0.; nrml[2] = sin( 45. );
    }
    // center
    else if ( bdry == EDGE5 ) {
      nrml[0] = sin( 45. ); nrml[1] = sin( 45. ); nrml[2] = 0.;
    }
    else if ( bdry == EDGE6 ) {
      nrml[0] = -sin( 45. ); nrml[1] = sin( 45. ); nrml[2] = 0.;
    }
    else if ( bdry == EDGE7 ) {
      nrml[0] = -sin( 45. ); nrml[1] = -sin( 45. ); nrml[2] = 0.;
    }
    else if ( bdry == EDGE8 ) {
      nrml[0] = sin( 45. ); nrml[1] = -sin( 45. ); nrml[2] = 0.;
    }
    // front  
    else if ( bdry == EDGE9 ) {
      nrml[0] = 0.; nrml[1] = -sin( 45. ); nrml[2] = sin( 45. );
    }
    else if ( bdry == EDGE10 ) {
      nrml[0] = -sin( 45. ); nrml[1] = 0.; nrml[2] = -sin( 45. );
    }
    else if ( bdry == EDGE11 ) {
      nrml[0] = 0.; nrml[1] = -sin( 45. ); nrml[2] = -sin( 45. );
    }
    else if ( bdry == EDGE12 ) {
      nrml[0] = sin( 45. ); nrml[1] = 0.; nrml[2] = -sin( 45. );
    }
    else
      // 3D models only have surface and edge boundaries
      throw csmp::Exception( ERROR, "Box::UnitNormalTo",
                             "Boundary flag could not be identified." );

  }

} // end UnitNormalTo



/**
compares the supplied string with valid BOX_BOUNDARY classifications;
returns true if string is equivaled to identifiers, except for
NOT, INTERNAL, IRREGULAR

@author SKM (10/2/2016)
*/
bool isDiagnosticBoxBoundaryClassifier( const string& i )
{
  if ( i == "LEFT" )     return true;
  if ( i == "RIGHT" )    return true;
  if ( i == "BOTTOM" )   return true;
  if ( i == "TOP" )      return true;
  if ( i == "FRONT" )    return true;
  if ( i == "BACK" )     return true;
  if ( i == "CNR1" )     return true;
  if ( i == "CNR2" )     return true;
  if ( i == "CNR3" )     return true;
  if ( i == "CNR4" )     return true;
  if ( i == "CNR5" )     return true;
  if ( i == "CNR6" )     return true;
  if ( i == "CNR7" )     return true;
  if ( i == "CNR8" )     return true;
  if ( i == "EDGE1" )    return true;
  if ( i == "EDGE2" )    return true;
  if ( i == "EDGE3" )    return true;
  if ( i == "EDGE4" )    return true;
  if ( i == "EDGE5" )    return true;
  if ( i == "EDGE6" )    return true;
  if ( i == "EDGE7" )    return true;
  if ( i == "EDGE8" )    return true;
  if ( i == "EDGE9" )    return true;
  if ( i == "EDGE10" )   return true;
  if ( i == "EDGE11" )   return true;
  if ( i == "EDGE12" )   return true;
  return false;

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
    while( it!=node_flags.end() ) {
         cout <<"\n\t"<< parseBoundary( (*it) ) <<": "<<  (equal_entries=node_flags.count( (*it) )) <<" nodes.";
         // avoiding printing of duplicates
         advance( it, equal_entries );
      }
    cout << endl;  
      
    multiset<BOX_BOUNDARY> elmt_flags;
    for ( auto eit=modeldomain.CellsBegin(); eit!=modeldomain.CellsEnd(); ++eit )
      for ( auto i{0U}; i<(*eit)->Neighbors(); ++i )
        elmt_flags.insert( (*eit)->AtBoundary(i) );
      
    it = elmt_flags.begin();
    while( it!=elmt_flags.end() ) {
         cout <<"\n\t"<< parseBoundary( (*it) ) <<": "<<  (equal_entries=elmt_flags.count( (*it) )) <<" elements.";
         advance( it, equal_entries );
      }
    cout << endl;  
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
    assert( eptr != nullptr );
    assert( b_face < eptr->Faces() );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( b_face >= eptr->Neighbors() ) {
         csmp_error.Note( ERROR, "atBoundary:", "element face number is out of range:", to_string(b_face) );
         return NOT;
      }

    // storing the box-boundary flags of the boundary face in a set
    set<BOX_BOUNDARY>  eflags;
    for ( const auto& fnit : eptr->CornerNodesOfFace(b_face) )
      if ( fnit->AtBoundary() != NOT ) eflags.insert( fnit->AtBoundary() );
    if ( eflags.empty() ) return NOT;

    // Case 1: one-dimensional model
    // -----------------------------
    if constexpr( dim == 1U )
      if ( eflags.size() == 1U ) return (*eflags.begin());
    
    // Case 2: two-dimensional model
    // -----------------------------
    if constexpr ( dim == 2U ) {
         // elements are considered boundary elements only if they have a face on the model boundary
         if ( eflags.size() == 1U ) return (*eflags.begin());
         // if there are two different flags
         else if ( eflags.size() == 2U ) {
              // if there is a corner involved, the other flag is chosen because an element must not span a corner
              // and its boundary face will lie on one of the sides of the model
              auto flag_it = eflags.begin();
              const BOX_BOUNDARY flag1 = (*flag_it);
              flag_it++;
              const BOX_BOUNDARY flag2 = (*flag_it);
              // intersections between internal and external boundaries
              if ( flag1 != NOT && flag2 == INTERNAL ) return flag1;
              if ( flag2 != NOT && flag1 == INTERNAL ) return flag2;
              // corner cases
              if ( isCorner(flag1) || flag1 == MULTIPLE ) return flag2;
              if ( isCorner(flag2) || flag2 == MULTIPLE ) return flag1;
              // degenerate cases where the corner was not flagged correctly
              if ( (flag1 == LEFT && flag2 == BOTTOM)  || (flag2 == LEFT   && flag1 == BOTTOM) ) return CNR1;
              if ( (flag1 == BOTTOM && flag2 == RIGHT) || (flag2 == BOTTOM && flag1 == RIGHT) ) return CNR2;
              if ( (flag1 == RIGHT && flag2 == TOP)    || (flag2 == RIGHT && flag1 == TOP) ) return CNR3;
              if ( (flag1 == TOP && flag2 == LEFT)     || (flag2 == TOP && flag1 == LEFT) ) return CNR4;
              // there should be no other cases because the 2D model has no edges
              cerr <<"\n\tmissed case: ";
              for ( const auto& boundary : eflags )
                cerr << parseBoundary( boundary ) << " ";
              csmp_error.Note( ERROR, "atBoundary(2D):", "did not succeed in finding a unique box boundary flag for cell.");
              return flag1;
           }
         // A face can only have 2D corner nodes
         else {
                 cerr <<"\n\n\tcell: "<< eptr->Idx() <<" ("<< parseFiniteElementType(eptr->FE_Type()) <<"), boundary flags:\n\t\t\t";
                 for ( auto i{0U}; i<eptr->Nodes(); ++i )
                   cerr <<" "<< eptr->N(i)->Idx() <<": "<< parseBoundary( eptr->N(i)->AtBoundary() );
                 cerr << endl;
                 csmp_error.Note( ERROR, "atBoundary(2D):", "unable to identify meaningful BOX_BOUNDARY flag for cell face.");
            }
          return MULTIPLE;
            
      } // if dim=2


    // Case 3: three-dimensional model
    // -------------------------------
    if constexpr ( dim == 3U )
      {
         // all nodes of the boundary face have the same flag (covers also line-element case)
         if ( eflags.size() == 1U ) return (*eflags.begin());
         
         // only for surface elements two different flags if they belong to the same face indicate a boundary position
         if ( eflags.size() >= 2U ) {
              // for faces of volume elements only the sides of a model are valid options
              if ( eptr->IsVolume() ) {
                   for ( auto bit : eflags )
                     if ( isSide(bit) )
                      return bit;
                }
              // if this is a face of a surface element it has only 2 corner nodes
              else if ( eptr->IsSurface() ) {
                   assert( eflags.size() == 2U );
                   const BOX_BOUNDARY flag1 = (*eflags.begin());
                   const BOX_BOUNDARY flag2 = (*next(eflags.begin(),1));
                   // if there is a corner involved, the other flag is chosen because the face of a surface cell cannot span a corner
                   if ( isCorner(flag1) ) return flag2;
                   if ( isCorner(flag2) ) return flag1;
                   // if this is a surface element with three nodes on the outside of the model
                   if ( isBOTTOM(flag1) && isBOTTOM(flag2) ) return BOTTOM;
                   if ( isRIGHT(flag1) && isRIGHT(flag2) ) return RIGHT;
                   if ( isTOP(flag1) && isTOP(flag2) ) return TOP;
                   if ( isLEFT(flag1) && isLEFT(flag2) ) return LEFT;
                   if ( isBACK(flag1) && isBACK(flag2) ) return BACK;
                   if ( isFRONT(flag1) && isFRONT(flag2) ) return BOTTOM;
                   assert( !isEdge(flag1) && !isEdge(flag2) );
                   // intersections between internal and external boundaries
                   if ( flag1 != NOT && flag2 == INTERNAL ) return flag1;
                   if ( flag2 != NOT && flag1 == INTERNAL ) return flag2;
                   return whichBoundary( flag1, flag2 );
                }
              // line element case was already covered
              else csmp_error.Note( ERROR, "atBoundary(3D):", "Line element face should only have a single flag.");
           }
         
      } // end dim=3

   return INTERNAL;
    
 } // end atBoundary

template BOX_BOUNDARY atBoundary( const Element<1U>* const, uint32_t );
template BOX_BOUNDARY atBoundary( const Element<2U>* const, uint32_t  );
template BOX_BOUNDARY atBoundary( const Element<3U>* const, uint32_t  );

// TESTING
/*
if ( eptr->FE()->IsLine() ) {
     cerr <<"\nelement "<< eptr->Idx() <<": "<< parseFiniteElementType(eptr->FE_Type()) <<", nodes:\n";
     for ( auto i{0U}; i<eptr->Nodes(); ++i )
       cerr <<" "<< eptr->N(i)->Idx() <<": "<< parseBoundary( eptr->N(i)->AtBoundary() );
     cerr << endl;
  }
*/



/// returns whether the cell is at the model boundary; this is so if all nodes of a line or surface element are flagged boundary or one face of a volume element
template<uint32_t dim, template<uint32_t> class CELL>
bool atBoundary( const CELL<dim>* const cptr )
 {
     // for surface and line elements all nodes must be at the boundary
     if ( !cptr->IsEquidimensional() ) {
          const size_t n_nodes{ cptr->Nodes() };
          for( auto i{0U}; i<n_nodes; ++i )
            if ( cptr->N(i)->AtBoundary() == NOT )
              return false;
       }
     else { // equidimensional elements
          const auto n_nbors{ cptr->Neighbors() };
          for ( auto i{0U}; i<n_nbors; ++i )
            if ( cptr->Neighbor(i) == nullptr )
              return true;
          return false;
       }
       
    return true;
       
 } // end atBoundary(bool)

template bool atBoundary( const Element<1U>* const );
template bool atBoundary( const Element<2U>* const );
template bool atBoundary( const Element<3U>* const );



/* FAILS AT TIMES

template<uint32_t dim, template<uint32_t> class CELL>
BOX_BOUNDARY atBoundary( const CELL<dim>* const eptr )
 {
    assert( eptr != nullptr );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    // nodal bflags are stored in a set
    set<BOX_BOUNDARY>  eflags;
    const size_t nodes( eptr->Nodes() );
    for ( auto i{0U}; i<nodes; ++i ) {
         assert( eptr->N(i) != nullptr );
         if ( eptr->N(i)->AtBoundary() != NOT )
           eflags.insert( eptr->N( i )->AtBoundary() );
      }
      
    // even if all flags=NOT, the element may lie on an internal boundary
    // in which case it is missing on of its neighbors
    if ( eflags.empty() ) {
         size_t n_nbors{eptr->Neighbors()};
         for ( auto i{0U}; i<n_nbors; ++i )
           if ( eptr->Neighbor(i) == nullptr )
             return INTERNAL;
         return NOT;
      }

    // if there were boundary flags
    // ----------------------------
    
    // Case 1: one-dimensional model
    // -----------------------------
    if constexpr( dim == 1 ) {
         if ( eflags.size() == 1U ) {
              if ( eptr->N(0)->AtBoundary() != NOT && eptr->N(1)->AtBoundary() != NOT )
                csmp_error.Note( ERROR, "atBoundary(1D):", "1D line element has the same 2 bflags.");
              return (*eflags.begin());
           }
         // there may only be 2 different flags when this is a single element model
         csmp_error.Note( ERROR, "atBoundary(1D):", "1D line element has the two different 2 bflags.");
      }
    
    // Case 2: two-dimensional model
    // -----------------------------
    if constexpr ( dim == 2 ) {
         // elements are considered boundary elements only if they have a face on the model boundary
         if ( eflags.size() == 1U ) {
              if ( eptr->FE()->IsLine() ) return (*eflags.begin());
              else {
                   // checking that there is indeed a face on the model boundary
                   const size_t n_nbors{eptr->Neighbors()};
                   for ( auto i{0U}; i<n_nbors; ++i )
                     if ( eptr->Neighbor(i) == nullptr )
                       return (*eflags.begin());
                }
           }
         // if there are two different flags
         if ( eflags.size() == 2U ) {
              // if there is a corner involved, the other flag is chosen because an element must not span a corner
              // and its boundary face will lie on one of the sides of the model
              auto flag_it = eflags.begin();
              const BOX_BOUNDARY flag1 = (*flag_it);
              flag_it++;
              const BOX_BOUNDARY flag2 = (*flag_it);
              if ( isCorner(flag1) ) return flag2;
              if ( isCorner(flag2) ) return flag1;
              // there should be no other cases because the 2D model has no edges
              cerr <<"\n\tmissed case: ";
              for ( auto boundary : eflags )
                cerr << parseBoundary( boundary ) << " ";
              csmp_error.Note( ERROR, "atBoundary(2D):", "did not succeed in finding a unique box boundary flag for cell.");
              return flag1;
           }
         // corner cases
         else if ( eflags.size() == 3U ) {
              // 3 flags are legitimate only for quadrilaterals at the corners of a rectangular model
              if ( !isQuadrilateral( eptr->FE_Type() ) ) {
                   cerr <<"\n\n\telement: "<< eptr->Idx() <<" ("<< parseFiniteElementType(eptr->FE_Type()) <<"), boundary flags:\n\t\t\t";
                   for ( auto i{0U}; i<eptr->Nodes(); ++i )
                     cerr <<" "<< eptr->N(i)->Idx() <<": "<< parseBoundary( eptr->N(i)->AtBoundary() );
                   cerr << endl;
                   csmp_error.Note( ERROR, "atBoundary(2D):",
                    "triangular element with two faces at boundary should be removed because it may cause problems when applying boundary conditions.");
                }
              else {
                   // from C++20 onwards use eflags.contains( BOTTOM ) );
                   if ( isIn( eflags, {BOTTOM,LEFT} ) ) return CNR1;
                   if ( isIn( eflags, {BOTTOM,RIGHT} ) ) return CNR2;
                   if ( isIn( eflags, {TOP,RIGHT} ) ) return CNR3;
                   if ( isIn( eflags, {LEFT,TOP} ) ) return CNR4;
                }
           }
         else {
                   cerr <<"\n\n\telement: "<< eptr->Idx() <<" ("<< parseFiniteElementType(eptr->FE_Type()) <<"), boundary flags:\n\t\t\t";
                   for ( auto i{0U}; i<eptr->Nodes(); ++i )
                     cerr <<" "<< eptr->N(i)->Idx() <<": "<< parseBoundary( eptr->N(i)->AtBoundary() );
                   cerr << endl;
                   csmp_error.Note( ERROR, "atBoundary(2D):", "all nodes of 2D triangular element appear to be located on boundary.");
              }
          return MULTIPLE;
            
      } // if dim=2


    // Case 3: three-dimensional model
    // -------------------------------
    if constexpr ( dim == 3U )
      {
         // only line and surface elements may be at boundary if there is only one boundary flag
         if ( eflags.size() == 1U && !eptr->IsVolume() )
           return (*eflags.begin());
         
         // only for surface elements two different flags if they belong to the same face indicate a boundary position
         if ( eflags.size() == 2U ) {
              if ( eptr->FE()->IsSurface() ) {
                   // if there is a corner involved, the other flag is chosen because an element must not span a corner
                   // and its boundary face will lie on one of the sides of the model
                   auto flag_it = eflags.begin();
                   const BOX_BOUNDARY flag1 = (*flag_it);
                   flag_it++;
                   const BOX_BOUNDARY flag2 = (*flag_it);
                   if ( isCorner(flag1) ) return flag2;
                   if ( isCorner(flag2) ) return flag1;
                   // the element face must be located on a model edge
                   return whichBoundary( flag1, flag2 );
                }
              // corner quadrilaterals that should not exist unless the element is a non simplex element
              if ( !isQuadrilateral( eptr->FE_Type() ) ) {
                   cerr <<"\n\n\telement: "<< eptr->Idx() <<" ("<< parseFiniteElementType(eptr->FE_Type()) <<"), boundary flags:\n\t\t\t";
                   for ( auto i{0U}; i<eptr->Nodes(); ++i )
                     cerr <<" "<< eptr->N(i)->Idx() <<": "<< parseBoundary( eptr->N(i)->AtBoundary() );
                   cerr << endl;
                   csmp_error.Note( ERROR, "atBoundary(3D):", "simplex element with two faces at boundary.");
                }
              else {
                   if ( isIn( eflags, {BOTTOM,LEFT} ) ) return CNR1;
                   if ( isIn( eflags, {BOTTOM,RIGHT} ) ) return CNR2;
                   if ( isIn( eflags, {TOP,RIGHT} ) ) return CNR3;
                   if ( isIn( eflags, {LEFT,TOP} ) ) return CNR4;
                }
           }
         // three flags may imply that a tetrahedral, pyramid, or prism element is located at boundary
         else if ( eflags.size() == 3U ) {
            const size_t n_nbors{eptr->Neighbors()};
            for ( auto i{0U}; i<n_nbors; ++i )
              if ( eptr->Neighbor(i) == nullptr )
                if ( isTriangular( eptr->FE()->ElementTypeOfFace(i) ) )
                  {
                     // diagnosing which side of the model this boundary face with 3 different flags is on
                     // (one of the flags must denote a side)
                     if ( isIn( eflags, {BOTTOM} ) ) return BOTTOM;
                     else if ( isIn( eflags, {RIGHT} ) ) return RIGHT;
                     else if ( isIn( eflags, {TOP} ) ) return TOP;
                     else if ( isIn( eflags, {LEFT} ) ) return LEFT;
                     else if ( isIn( eflags, {BACK} ) ) return BACK;
                     else if ( isIn( eflags, {FRONT} ) ) return FRONT;
                     else if ( isIn( eflags, {IRREGULAR} ) ) return IRREGULAR;
                     else {
                          // something went wrong
                          cerr <<"\n\n\ttriangular element face: "<< i <<", boundary flags:\n\t\t\t";
                          vector<uint32_t> fnids;
                          eptr->FE()->NodesOfFace( i, fnids );
                          for ( size_t j{0}; j<fnids.size(); ++j )
                            cerr <<" "<< eptr->N( fnids[j] )->Idx() <<": "<< parseBoundary( eptr->N( fnids[j] )->AtBoundary() );
                          cerr << endl;
                          csmp_error.Note( ERROR, "atBoundary(3D):", "could not resolve placement of triangular boundary face.");
                       }
                  }
           }
         // 4 flags may imply that a pyramid, prism or hexahedral element is located at boundary
         else if ( eflags.size() == 4U ) {
            const size_t n_nbors{eptr->Neighbors()};
            for ( auto i{0U}; i<n_nbors; ++i )
              if ( eptr->Neighbor(i) == nullptr )
                if ( isQuadrilateral( eptr->FE()->ElementTypeOfFace(i) ) )
                  {
                     // diagnosing which side of the model this boundary face with 3 different flags is on
                     // (one of the flags must denote a side)
                     if ( isIn( eflags, {BOTTOM} ) ) return BOTTOM;
                     else if ( isIn( eflags, {RIGHT} ) ) return RIGHT;
                     else if ( isIn( eflags, {TOP} ) ) return TOP;
                     else if ( isIn( eflags, {LEFT} ) ) return LEFT;
                     else if ( isIn( eflags, {BACK} ) ) return BACK;
                     else if ( isIn( eflags, {FRONT} ) ) return FRONT;
                     else if ( isIn( eflags, {IRREGULAR} ) ) return IRREGULAR;
                     else {
                          // something went wrong
                          cerr <<"\n\n\tquadrilateral element face: "<< i <<", boundary flags:\n\t\t\t";
                          vector<uint32_t> fnids;
                          eptr->FE()->NodesOfFace( i, fnids );
                          for ( size_t j{0}; j<fnids.size(); ++j )
                            cerr <<" "<< eptr->N( fnids[j] )->Idx() <<": "<< parseBoundary( eptr->N( fnids[j] )->AtBoundary() );
                          cerr << endl;
                          csmp_error.Note( ERROR, "atBoundary(3D):", "could not resolve placement of quadrilateral boundary face.");
                       }
                  }
           }
         // more than 4 flags
         else {
              // edge and corner cases
              auto sit = eflags.begin();
              const BOX_BOUNDARY flag1 = (*sit); sit++;
              const BOX_BOUNDARY flag2 = (*sit); sit++;
              const BOX_BOUNDARY flag3 = (*sit);
              // corner combinations
              // CNR1 (order from most negative integer to highest value)
              if ( flag1 == BACK and flag2 == BOTTOM  and flag3 == LEFT ) return CNR1;
              // CNR2
              if ( flag1 == BACK and flag2 == BOTTOM  and flag3 == RIGHT ) return CNR2;
              // CNR3
              if ( flag1 == BACK and flag2 == TOP     and flag3 == RIGHT ) return CNR3;
              // CNR4
              if ( flag1 == BACK and flag2 == TOP     and flag3 == LEFT ) return CNR4;
              // CNR5
              if ( flag1 == FRONT and flag2 == BOTTOM and flag3 == LEFT ) return CNR5;
              // CNR6
              if ( flag1 == FRONT and flag2 == BOTTOM and flag3 == RIGHT ) return CNR6;
              // CNR7
              if ( flag1 == FRONT and flag2 == TOP    and flag3 == RIGHT ) return CNR7;
              // CRN8
              if ( flag1 == FRONT and flag2 == TOP    and flag3 == LEFT ) return CNR8;
              
              // erase basic boundary options to focus on edges
              eflags.erase( LEFT );
              eflags.erase( RIGHT );
              eflags.erase( BOTTOM );
              eflags.erase( TOP );
              eflags.erase( FRONT );
              eflags.erase( BACK );

              if ( !eflags.empty() )
                {
                  const BOX_BOUNDARY flag_min = (*min_element( eflags.begin(), eflags.end() ));
                  const BOX_BOUNDARY flag_max = (*max_element( eflags.begin(), eflags.end() ));
                  // if a corner is contained that corner flag is choosen
                  if ( flag_max <= CNR1 and flag_max >= CNR8 ) return flag_max;
                  if ( flag_min <= CNR1 and flag_min >= CNR8 ) return flag_min;
                  // if the first integer entry in the set is an edge, that flag is chosen
                  if ( flag_max <= EDGE1 and flag_max > INTERNAL ) return flag_max;
                  if ( flag_min <= EDGE1 and flag_min > INTERNAL ) return flag_min;
                  // if a corner is contained that corner flag is choosen
                  if ( (*eflags.begin()) <= CNR1 and
                       (*eflags.begin()) >= CNR8 ) return (*eflags.begin());
                  // if the first integer entry in the set is an edge, that flag is chosen
                  if ( (*eflags.begin()) >= EDGE1 and
                       (*eflags.begin()) <  INTERNAL ) return (*eflags.begin());
                  if ( (*eflags.begin()) == IRREGULAR ) return IRREGULAR;
                  // no success
                  cerr << "\n\n\natBoundary(3D): unable to determine box boundary flag for element: ";
                  cerr << eptr->Idx() <<" ("<< parseFiniteElementType(eptr->FE_Type()) <<") with the boundary flags:\n\t\t\t";
                  for ( auto boundary : eflags )
                    cerr << parseBoundary( boundary ) << " ";
                  cerr << endl;
                }
              csmp_error.Note( ERROR, "atBoundary(3D):", "implausible case where more than 4 element nodes are located on model boundary.");
           }
         
        return MULTIPLE;
         
      } // end dim=3

   return MULTIPLE;
    
 } // end atBoundary(element pointer)

template BOX_BOUNDARY atBoundary( const Element<1U>* const );
template BOX_BOUNDARY atBoundary( const Element<2U>* const );
template BOX_BOUNDARY atBoundary( const Element<3U>* const );
*/



/// identifying the boundary that a lower dimensional element is located on from the BOX_BOUNDARY flags assigned to its nodes
template<uint32_t dim>
BOX_BOUNDARY atBoundary( const set<BOX_BOUNDARY>& node_flags, uint32_t boundary_face_corner_nodes )
 {
    static_assert( dim != 1U, "use atBoundary( const set<BOX_BOUNDARY>&, uint32_t ): only for 2 or 3D models");
    assert( boundary_face_corner_nodes >= 2U );
    
    if ( node_flags.empty() ) return NOT;
    
    const auto n_flags = node_flags.size();
   
    // TWO-DIMENSIONAL MODELS - lower-dim elements are line elements
    if constexpr ( dim == 2U ) {
         // one flag
         if ( n_flags == 1U) {
              switch( (*node_flags.begin())  ) {
                   case LEFT:
                     return LEFT;
                   case RIGHT:
                     return RIGHT;
                   case BOTTOM:
                     return BOTTOM;
                   case TOP:
                     return TOP;
                   // edges
                   case EDGE1:
                     return EDGE1;
                   case EDGE2:
                     return EDGE2;
                   case EDGE3:
                     return EDGE3;
                   case EDGE4:
                     return EDGE4;
                   default:
                     return IRREGULAR;
                }
           }
         else if ( n_flags == 2U ) { // set with 2 flags must contain a corner flag
             const BOX_BOUNDARY b1{ (*node_flags.begin()) }, b2{ *next(node_flags.begin(),1) };
             // the higher bflag value would be the corner value
             if ( isCorner(b2) ) {
                  switch( b1 ) {
                       // edges have lower negative int values than corners
                       case EDGE1:
                         return EDGE1;
                       case EDGE2:
                         return EDGE2;
                       case EDGE3:
                         return EDGE3;
                       case EDGE4:
                         return EDGE4;
                       default:
                         return IRREGULAR;
                    }
               }
             else if ( isCorner(b1) ) {
                  switch( b2 ) {
                       case LEFT:
                         return LEFT;
                       case RIGHT:
                         return RIGHT;
                       case BOTTOM:
                         return BOTTOM;
                       case TOP:
                         return TOP;
                       default:
                         return IRREGULAR;
                    }
              }
          }
         else { // error
            cerr <<"\n"<<"BOX_BOUNDARY atBoundary<2>( const set<BOX_BOUNDARY>&, face_cnr_nodes ): Error: could not resolve boundary status."<< endl;
            cerr <<"\t"<<"input data: corner node flags of element face: "<< boundary_face_corner_nodes << endl;
            cerr <<"\t"<<"boundary node flags:\n";
            for ( const auto& fit : node_flags ) cerr <<" "<< parseBoundary(fit);
            cerr << endl;
          }
          
      } // end 2D
      
      
    // THREE-DIMENSIONAL MODELS - lower dimensional elements are line or surface elements
    if constexpr ( dim == 3U ) {
         if ( boundary_face_corner_nodes == 2U ) { // line element case
             // all node flags are the same (size 1)
             if ( n_flags == 1U ) {
                    switch( (*node_flags.begin())  ) {
                         case EDGE1:
                           return EDGE1;
                         case EDGE2:
                           return EDGE2;
                         case EDGE3:
                           return EDGE3;
                         case EDGE4:
                           return EDGE4;
                         case EDGE5:
                           return EDGE5;
                         case EDGE6:
                           return EDGE6;
                         case EDGE7:
                           return EDGE7;
                         case EDGE8:
                           return EDGE8;
                         case EDGE9:
                           return EDGE9;
                         case EDGE10:
                           return EDGE10;
                         case EDGE11:
                           return EDGE11;
                         case EDGE12:
                           return EDGE12;
                         case LEFT:
                           return LEFT;
                         case RIGHT:
                           return RIGHT;
                         case BOTTOM:
                           return BOTTOM;
                         case TOP:
                           return TOP;
                         case FRONT:
                           return FRONT;
                         case BACK:
                           return BACK;
                         default:
                           return IRREGULAR;
                      }
               }
             // line segments with two node flags (line segments)
             else if ( n_flags == 2U ) {
                 const BOX_BOUNDARY b1{ (*node_flags.begin()) }, b2{ *next(node_flags.begin(),1) };
                 // when the higher bflag value is the corner flag value
                 if ( isCorner(b2) ) {
                      switch( b1 ) {
                           // edges have lower negative int values than corners
                           case EDGE1:
                             return EDGE1;
                           case EDGE2:
                             return EDGE2;
                           case EDGE3:
                             return EDGE3;
                           case EDGE4:
                             return EDGE4;
                           case EDGE5:
                             return EDGE5;
                           case EDGE6:
                             return EDGE6;
                           case EDGE7:
                             return EDGE7;
                           case EDGE8:
                             return EDGE8;
                           case EDGE9:
                             return EDGE9;
                           case EDGE10:
                             return EDGE10;
                           case EDGE11:
                             return EDGE11;
                           case EDGE12:
                             return EDGE12;
                           default:
                             return IRREGULAR;
                        }
                   }
                 // 2 different flags and the corner flag value has a lower value than these boundary flags
                 else if ( isCorner(b1) ) {
                      switch( b2 ) {
                           // edges have lower negative int values than corners
                           case EDGE1:
                             return EDGE1;
                           case EDGE2:
                             return EDGE2;
                           case EDGE3:
                             return EDGE3;
                           case EDGE4:
                             return EDGE4;
                           case EDGE5:
                             return EDGE5;
                           case EDGE6:
                             return EDGE6;
                           case EDGE7:
                             return EDGE7;
                           case EDGE8:
                             return EDGE8;
                           case EDGE9:
                             return EDGE9;
                           case EDGE10:
                             return EDGE10;
                           case EDGE11:
                             return EDGE11;
                           case EDGE12:
                             return EDGE12;
                           default:
                             return IRREGULAR;
                        }
                  }
               }
           }
           
         // when the element face at the boundary is a triangle or quadrilateral
         else if ( boundary_face_corner_nodes > 2U )
           {
             if ( n_flags == 1U ) {
                    switch( (*node_flags.begin())  ) {
                         case LEFT:
                           return LEFT;
                         case RIGHT:
                           return RIGHT;
                         case BOTTOM:
                           return BOTTOM;
                         case TOP:
                           return TOP;
                         case FRONT:
                           return FRONT;
                         case BACK:
                           return BACK;
                         default:
                           return IRREGULAR;
                      }
               }
             // two node flags
             else if ( n_flags == 2U ) {
                 const BOX_BOUNDARY b1{ (*node_flags.begin()) }, b2{ *next(node_flags.begin(),1) };
                  // corners
                  if ( isCorner(b1) || isCorner(b2) ) {
                      cerr <<"\n"<<"BOX_BOUNDARY atBoundary<3>(): Error: corner face must have a mininum of 3 different flags."<< endl;
                      cerr <<"\t"<<"input data: corner node flags of element face: "<< boundary_face_corner_nodes << endl;
                      cerr <<"\t"<<"boundary node flags:\n";
                      for ( const auto& fit : node_flags ) cerr <<" "<< parseBoundary(fit);
                      cerr << endl;
                    }
                 // edge flags should be discerned by lower values than side flags
                 if ( isEdge(b1) ) {
                      switch( b2 ) {
                           // edges have lower negative int values than corners
                           case BACK:
                             return BACK;
                           case FRONT:
                             return FRONT;
                           case LEFT:
                             return LEFT;
                           case RIGHT:
                             return RIGHT;
                           case TOP:
                             return TOP;
                           case BOTTOM:
                             return BOTTOM;
                           default:
                             return IRREGULAR;
                        }
                   }
                 if ( isEdge(b2) ) {
                      switch( b1 ) {
                           // edges have lower negative int values than corners
                           case BACK:
                             return BACK;
                           case FRONT:
                             return FRONT;
                           case LEFT:
                             return LEFT;
                           case RIGHT:
                             return RIGHT;
                           case TOP:
                             return TOP;
                           case BOTTOM:
                             return BOTTOM;
                           default:
                             return IRREGULAR;
                        }
                   }
               }
             // three node flags (surfaces)
             else if ( n_flags == 3U ) {
                  const BOX_BOUNDARY b1{ (*node_flags.begin()) }, b2{ *next(node_flags.begin(),1) }, b3{ *next(node_flags.begin(),2) };
                  if ( isLEFT(b1) && isLEFT(b2) && isLEFT(b3) ) return LEFT;
                  if ( isRIGHT(b1) && isRIGHT(b2) && isRIGHT(b3) ) return RIGHT;
                  if ( isBOTTOM(b1) && isBOTTOM(b2) && isBOTTOM(b3) ) return BOTTOM;
                  if ( isTOP(b1) && isTOP(b2) && isTOP(b3) ) return TOP;
                  if ( isFRONT(b1) && isFRONT(b2) && isFRONT(b3) ) return FRONT;
                  if ( isBACK(b1) && isBACK(b2) && isBACK(b3) ) return BACK;
                  return IRREGULAR;
               }
             // boundary face is a quadrilateral (surface) containing a corner node
             else if ( n_flags == 4U ) {
                  const BOX_BOUNDARY b1{ (*node_flags.begin()) }, b2{ *next(node_flags.begin(),1) },
                                     b3{ *next(node_flags.begin(),2) }, b4{ *next(node_flags.begin(),3) };
                  if ( isLEFT(b1) && isLEFT(b2) && isLEFT(b3) && isLEFT(b4) ) return LEFT;
                  if ( isRIGHT(b1) && isRIGHT(b2) && isRIGHT(b3) && isRIGHT(b4) ) return RIGHT;
                  if ( isBOTTOM(b1) && isBOTTOM(b2) && isBOTTOM(b3) && isBOTTOM(b4) ) return BOTTOM;
                  if ( isTOP(b1) && isTOP(b2) && isTOP(b3) && isTOP(b4) ) return TOP;
                  if ( isFRONT(b1) && isFRONT(b2) && isFRONT(b3) && isFRONT(b4) ) return FRONT;
                  if ( isBACK(b1) && isBACK(b2) && isBACK(b3) && isBACK(b4) ) return BACK;
                  return IRREGULAR;
               }
             else { // error
                  cerr <<"\n"<<"BOX_BOUNDARY atBoundary<3>(): Error: face can have a maximum of 4 flags."<< endl;
                  cerr <<"\t"<<"input data: corner node flags of element face: "<< boundary_face_corner_nodes << endl;
                  cerr <<"\t"<<"boundary node flags:\n";
                  for ( const auto& fit : node_flags ) cerr <<" "<< parseBoundary(fit);
                  cerr << endl;
                  return IRREGULAR;
               }
             
           } // end for surface faces

      } // end 3D

  // error
  cerr <<"\n"<<"BOX_BOUNDARY atBoundary<3>( const set<BOX_BOUNDARY>& ): Error: could not resolve boundary status."<< endl;
  cerr <<"\t"<<"input data: corner nodes of element face: "<< boundary_face_corner_nodes << endl;
  cerr <<"\t"<<"boundary node flags:\n";
  for ( const auto& fit : node_flags ) cerr <<" "<< parseBoundary(fit);
  cerr << endl;

  return IRREGULAR;

} // end atBoundary(set<BOX_BOUNDARY>)

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

  for ( auto nit = mregion.NodesBegin(); nit != mregion.NodesEnd(); nit++ )
    (*nit)->Store( nprop_key, makeScalar( ANY, static_cast<double>((*nit)->AtBoundary()) ) );

  for ( auto eit = mregion.CellsBegin(); eit != mregion.CellsEnd(); eit++ ) {
       uint32_t face = numeric_limits<uint32_t>::max();
       bool at_boundary{false};
       for ( auto i{0U}; i<(*eit)->Neighbors(); ++i )
         if ( (*eit)->Neighbor(i) == nullptr ) {
              face        = i;
              at_boundary = true;
              break;
           }
       if ( at_boundary )
         (*eit)->Store( eprop_key, makeScalar( ANY, static_cast<double>(atBoundary(*eit,face)) ) );
       else
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



/// return which boundary the edge or lower-dim face with the two end-node flags is on
BOX_BOUNDARY  whichBoundary( BOX_BOUNDARY node_flag1, BOX_BOUNDARY node_flag2 )
 {
    if ( node_flag1 == node_flag2 ) return node_flag1;
    // sides
    if ( isLEFT(node_flag1)   && isLEFT(node_flag2) )   return LEFT;
    if ( isRIGHT(node_flag1)  && isRIGHT(node_flag2) )  return RIGHT;
    if ( isBOTTOM(node_flag1) && isBOTTOM(node_flag2) ) return BOTTOM;
    if ( isTOP(node_flag1)    && isTOP(node_flag2) )    return TOP;
    if ( isFRONT(node_flag1)  && isFRONT(node_flag2) )  return FRONT;
    if ( isBACK(node_flag1)   && isBACK(node_flag2) )   return BACK;
    // edges
    if ( node_flag1 == EDGE1 || node_flag2 == EDGE1 ) return EDGE1;
    if ( node_flag1 == EDGE2 || node_flag2 == EDGE2 ) return EDGE2;
    if ( node_flag1 == EDGE3 || node_flag2 == EDGE3 ) return EDGE3;
    if ( node_flag1 == EDGE4 || node_flag2 == EDGE4 ) return EDGE4;
    if ( node_flag1 == EDGE5 || node_flag2 == EDGE5 ) return EDGE5;
    if ( node_flag1 == EDGE6 || node_flag2 == EDGE6 ) return EDGE6;
    if ( node_flag1 == EDGE7 || node_flag2 == EDGE7 ) return EDGE7;
    if ( node_flag1 == EDGE8 || node_flag2 == EDGE8 ) return EDGE8;
    if ( node_flag1 == EDGE9 || node_flag2 == EDGE9 ) return EDGE9;
    if ( node_flag1 == EDGE10 || node_flag2 == EDGE10 ) return EDGE10;
    if ( node_flag1 == EDGE11 || node_flag2 == EDGE11 ) return EDGE11;
    if ( node_flag1 == EDGE12 || node_flag2 == EDGE12 ) return EDGE12;
    
    // corner cases are handled in first block
    
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
void flagEdges( VSet<dim>& vset, double tol, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax )
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





template<uint32_t dim>
void flagCornerNodes( VSet<dim>& vset, double tol, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     // finding min-max of x, and y coordinates of boundary nodes
     map<size_t,int8_t>  bflags;
     size_t              flagging_count(0U);
     
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
                     { (*bit) = CNR_MIN_MAXX; flagging_count++; }
                   // CNR3
                   else if ( approximatelyEqual( vset.Px( n_node ), xmax, tol ) && approximatelyEqual( vset.Py( n_node ), ymax ) )
                     { (*bit) = CNR_MAX_MAXX; flagging_count++; }
                   // CNR4
                   else if ( approximatelyEqual( vset.Px( n_node ), xmin, tol ) && approximatelyEqual( vset.Py( n_node ), ymax ) )
                     { (*bit) = CNR_MAX_MINXZ; flagging_count++; }
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
                          { (*bit) = CNR_MIN_MAXX; flagging_count++; }
                        // CNR3
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_MAX_MAXX; flagging_count++; }
                        // CNR4
                        else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_MAX_MINXZ; flagging_count++; }
                     }
                   else if ( approximatelyEqual( vset.Pz( n_node ), zmax, tol ) )
                     {
                        // CNR5
                        if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                          { (*bit) = CNR_MIN_MAXZ; flagging_count++; }
                        // CNR6
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymin,tol) )
                          { (*bit) = CNR_MIN_MAXXZ; flagging_count++; }
                        // CNR7
                        else if ( approximatelyEqual(vset.Px( n_node ),xmax,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_MAX; flagging_count++; }
                        // CNR8
                        else if ( approximatelyEqual(vset.Px( n_node ),xmin,tol) && approximatelyEqual(vset.Py( n_node ),ymax,tol) )
                          { (*bit) = CNR_MAX_MAXZ; flagging_count++; }
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
