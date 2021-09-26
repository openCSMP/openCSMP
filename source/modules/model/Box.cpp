#include "Box.h"
#include "Boundary.h"
#include "Element.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "Exception.h"

#include "SmallSet.h"

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



BOX_BOUNDARY intToBOX_BOUNDARY( long64 i )
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

  //cout <<"\nintToSG_BOUNDARY(int): unable to parse integer: "<< i << endl;
  return NOT;
}


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

  cerr << "\nparseBoundary(BOX_BOUNDARY): unable to parse BOX_BOUNDARY: " << i << endl;
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
  return NOT;
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
void Box::UnitNormalTo( BOX_BOUNDARY bdry, size_t dim, vector<double64>& nrml ) const
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
template<size_t dim>
void printBoxBoundaryFlags( const Model<dim>& model )
 {
    const Region<dim>& modeldomain(model.Region("Model"));
    
    cout <<"\n\nprintBoxBoundaryFlags: Current flags and numbers assigned to node and element objects:\n";
    size_t equal_entries(0);
    multiset<BOX_BOUNDARY> node_flags;
    for ( typename vector<Node<dim>*>::const_iterator 
          nit=modeldomain.NodesBegin(); nit!=modeldomain.NodesEnd(); ++nit ) 
      node_flags.insert( (*nit)->AtBoundary() );
    multiset<BOX_BOUNDARY>::const_iterator it=node_flags.begin();
    while( it!=node_flags.end() ) {
         cout <<"\n\t"<< parseBoundary( (*it) ) <<": "<<  (equal_entries=node_flags.count( (*it) )) <<" nodes.";
         // avoiding printing of duplicates
         advance( it, equal_entries );
      }
    cout << endl;  
      
    multiset<BOX_BOUNDARY> elmt_flags;
    for ( typename vector<Element<dim>*>::const_iterator 
          eit=modeldomain.ElementsBegin(); eit!=modeldomain.ElementsEnd(); ++eit ) 
      elmt_flags.insert( (*eit)->AtBoundary() );
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
Using the side boundaries of the model,
the method recreateBoxBoundaryFlags recreates the corresponding box-boundary flagging.

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
    for ( vector<Element<2U>*>::const_iterator it=modeldomain.ElementsBegin(); it!=modeldomain.ElementsEnd(); ++it )
      {  // current version only works for linear quadrilaterals
         assert( (*it)->Nodes() <= 5 );
         
         if ( !isQuadrilateral( (*it)->FE_Type() ) )
           throw csmp::Exception( ERROR, "recreateBoxBoundaryFlagsForQuadrilateralModel", "this method only works for quadrilateral elements");
           
         for ( size_t i=0; i<(*it)->Nodes(); ++i )
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
    Region<3U>&    modeldomain(model.Region("Model"));
    vector<size_t> fnids;
    
    for ( vector<Element<3U>*>::const_iterator it=modeldomain.ElementsBegin(); it!=modeldomain.ElementsEnd(); ++it )
      {
         // current version only works for linear hexahedra
         assert( (*it)->Nodes() <= 9 );

         if ( !isHexahedral( (*it)->FE_Type() ) )
           throw csmp::Exception( ERROR, "recreateBoxBoundaryFlagsForHexahedralModel", "this method only works for quadrilateral elements");
           
         // DEFAULT (not at any boundary)
         for ( size_t i=0; i<(*it)->Nodes(); ++i )
           (*it)->N(i)->AtBoundary( NOT );
           
         // BOTTOM 
         if ( (*it)->Neighbor(0) == nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              (*it)->FE()->NodesOfFace( 0U, fnids );
              for ( size_t j=0U; j<fnids.size(); ++j )
                (*it)->N( fnids[j] )->AtBoundary( BOTTOM );
              continue;
           }
         // RIGHT 
         if ( (*it)->Neighbor(1) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(2) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              (*it)->FE()->NodesOfFace( 2U, fnids );
              for ( size_t j=0U; j<fnids.size(); ++j )
                (*it)->N( fnids[j] )->AtBoundary( RIGHT );
              continue;
           }
         // TOP 
         if ( (*it)->Neighbor(2) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(3) != nullptr ) {
              (*it)->FE()->NodesOfFace( 5U, fnids );
              for ( size_t j=0U; j<fnids.size(); ++j )
                (*it)->N( fnids[j] )->AtBoundary( TOP );
              continue;
           }
         // LEFT 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              (*it)->FE()->NodesOfFace( 4U, fnids );
              for ( size_t j=0U; j<fnids.size(); ++j )
                (*it)->N( fnids[j] )->AtBoundary( LEFT );
              continue;
           }
         // FRONT 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              (*it)->FE()->NodesOfFace( 1U, fnids );
              for ( size_t j=0U; j<fnids.size(); ++j )
                (*it)->N( fnids[j] )->AtBoundary( FRONT );
              continue;
           }
         // BACK 
         if ( (*it)->Neighbor(3) == nullptr && (*it)->Neighbor(0) != nullptr && (*it)->Neighbor(1) != nullptr && (*it)->Neighbor(2) != nullptr ) {
              (*it)->FE()->NodesOfFace( 3U, fnids );
              for ( size_t j=0U; j<fnids.size(); ++j )
                (*it)->N( fnids[j] )->AtBoundary( BACK );
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




/**
    infers from the node flags which boundary the element lies on; returns MULTIPLE is ambiguous; throws exception if inconsistent
    
      @author SKM
      @date 4/7/20
      @test OK
*/
template<size_t dim, template<size_t> class CELL>
BOX_BOUNDARY atBoundary( const CELL<dim>* const eptr ) 
 {
    assert( eptr != nullptr );
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    // nodal bflags are stored in a set
    SmallSet<BOX_BOUNDARY>  eflags; 
    const size_t nodes( eptr->Nodes() );
    for ( size_t i = 0U; i<nodes; ++i ) {
         assert( eptr->N(i) != nullptr );                
         if ( eptr->N(i)->AtBoundary() != NOT )
           eflags.insert( eptr->N( i )->AtBoundary() );
      }
      
    // diagnosis
    if ( eflags.empty() ) return NOT;

    // if there were boundary flags
    // ----------------------------
    // if only a single unique flag was contained the decision is easy
    if ( eflags.size() == 1U ) return (*eflags.begin());
    
    // if this is a 2D model so that there are only CNR1 to CRN4 and BOTTOM, TOP, RIGHT, LEFT, INTERNAL
    if ( dim == 2U ) {
         if ( eflags.size() == 2U ) {
              // if there is a corner involved, note: corners will always have more negative values than sides
              const BOX_BOUNDARY flag1 = (*eflags.begin());
              if ( isCorner(flag1) ) return flag1;
              if ( flag1 == INTERNAL ) return MULTIPLE;
              // there should be no other cases because the 2D model has no edges
              cerr <<"\n\tmissed case: ";
              for ( auto boundary : eflags )
                cerr << parseBoundary( boundary ) << " ";
              csmp_error.notice( ERROR, "atBoundary(2D):", "did not succeed in finding a unique box boundary flag for cell.");
              return flag1;
           }
         else {
            // if there is a corner involved, note: corners will always have more negative values than sides
            const BOX_BOUNDARY flag1 = (*eflags.begin());
            if ( isCorner(flag1) ) return flag1;
            if ( flag1 == INTERNAL ) return MULTIPLE;
            if ( isTriangular( eptr->FE_Type() ) ) {
                 cerr <<"\n\tmissed case: ";
                 for ( auto boundary : eflags )
                   cerr << parseBoundary( boundary ) << " ";
                 csmp_error.notice( ERROR, "atBoundary:", "all nodes of 2D triangular element appear to be located on boundary.");
              } 
         }
      }
      
    // if there are 2 flags and one of them is IRREGULAR, it is removed
    if ( eflags.size() == 2U )
      {
        // 1. case handled already: if all not-NOT boundary flags are IRREGULAR so should be the element 
        // 2. if one of the 2 flags is IRREGULAR, it is removed
        if ( eflags.count( IRREGULAR ) ) {
            eflags.erase( IRREGULAR );
            BOX_BOUNDARY flag = (*eflags.begin());
            if ( isEdge(flag) || isCorner(flag) ) return flag;
            // the edge is determined by the intersection of the side boundary with the irregular boundary assumed to have been on the top
            // NOTE: this follows the convention that box-shaped models can have an IRREGULAR top boundary
            if ( flag == BACK )  return EDGE3;
            if ( flag == FRONT ) return EDGE11;
            if ( flag == LEFT )  return EDGE8;
            if ( flag == RIGHT ) return EDGE7;
            return  flag; // choosing the flag if a more specific one cannot be found
          }
        else {
            // 12 possibilities relating to edge elements
            // NOTE: the order of the flags is important here !
            auto sit = eflags.begin();
            const BOX_BOUNDARY flag1 = (*sit); sit++;
            const BOX_BOUNDARY flag2 = (*sit);
            // BACK_BOTTOM
            if ( flag1 == BACK and flag2 == BOTTOM ) return EDGE1;
            // BACK_RIGHT
            if ( flag1 == BACK and flag2 == RIGHT ) return EDGE2;
            // BACK_TOP
            if ( flag1 == BACK and flag2 == TOP ) return EDGE3;
            // BACK_LEFT
            if ( flag1 == BACK and flag2 == LEFT ) return EDGE4;

            // BOTTOM_LEFT
            if ( flag1 == BOTTOM and flag2 == LEFT ) return EDGE5;
            // BOTTOM_RIGHT
            if ( flag1 == BOTTOM and flag2 == RIGHT ) return EDGE6;
            // TOP_RIGHT
            if ( flag1 == TOP and flag2 == RIGHT ) return EDGE7;
            // TOP_LEFT
            if ( flag1 == TOP and flag2 == LEFT ) return EDGE8;
            
            // FRONT_BOTTOM
            if ( flag1 == FRONT and flag2 == BOTTOM ) return EDGE9;
            // FRONT_RIGHT
            if ( flag1 == FRONT and flag2 == RIGHT ) return EDGE10;
            // FRONT_TOP
            if ( flag1 == FRONT and flag2 == TOP ) return EDGE11;
            // FRONT_LEFT
            if ( flag1 == FRONT and flag2 == LEFT ) return EDGE12;
            
            // if a corner is contained that corner flag is choosen
            if ( flag1 <= CNR1 and flag1 >= CNR8 ) return flag1;
            if ( flag2 <= CNR1 and flag2 >= CNR8 ) return flag2;
            // if the first integer entry in the set is an edge, that flag is chosen
            if ( flag1 <= EDGE1 and flag1 > INTERNAL ) return flag1;
            if ( flag2 <= EDGE1 and flag2 > INTERNAL ) return flag2;
            cerr << "\n\tmissed case: " << parseBoundary( flag1 ) << ", " << parseBoundary( flag2 ) << endl;
         }
      }
    else if ( eflags.size() > 2U ) 
      {
        // corner cases
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
            cerr << "\n\n\natBoundary: unable to determine box boundary flag for element:\n";
            for ( auto boundary : eflags )
              cerr << parseBoundary( boundary ) << " ";
            cerr << endl;
          }
        else {
            csmp_error.notice( WARNING, "atBoundary:", "could not discern which boundary elemnt is located on.");
          }
     }

   return MULTIPLE;
    
 } // end atBoundary

template BOX_BOUNDARY atBoundary( const Element<1U>* const ); 
template BOX_BOUNDARY atBoundary( const Element<2U>* const ); 
template BOX_BOUNDARY atBoundary( const Element<3U>* const ); 




/**
Assuming that the nodes of the box-shaped model are flagged correctly,
this method considers the existing combinations of node flags and assign
the appropriate BOX_BOUNDARY flags to all elements of the model.

Furthermore the methods assumes that provided iterator range encompasses
all elements.
*/
/*
template<size_t dim>
void flagElementsUsingNodal_BOX_BOUNDARY_Flags( typename std::deque<csmp::Element<dim>* >::iterator it,
                                                typename std::deque<csmp::Element<dim>* >::iterator last_elmt )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  assert( it != last_elmt );
  
  size_t nodes_at_boundary(0U);
  
  while ( it != last_elmt ) 
    {
      assert( (*it) != nullptr );
      // the bflags of each element are stored in a set
      SmallSet<BOX_BOUNDARY>  eflags; // this actually is a multiset!
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0U; i<nodes; ++i ) {
           assert( (*it)->N(i) != nullptr );                
           if ( (*it)->N(i)->AtBoundary() != NOT && (*it)->N(i)->AtBoundary() != INTERNAL ) {
                eflags.insert( (*it)->N( i )->AtBoundary() );
                nodes_at_boundary++;
             }
        }

// testing
//if ( !eflags.empty() ) {
//     cerr <<"\n"<< (*it)->Idx() <<": ";
//     for ( auto i=eflags.begin(); i!=eflags.end(); ++ i )
//       cerr << parseBoundary( (*i) ) <<" ";
//  }       
        
      // if no identifier could be found the boundary flag is set to NOT
      if ( eflags.empty() )
        (*it)->AtBoundary( NOT );
        
      // if only a single flag is contained the decision is easy
      else if ( eflags.size() == 1U )
        (*it)->AtBoundary( (*eflags.begin()) );
        
      // if there are 2 flags and one of them is IRREGULAR, it is removed
      else if ( eflags.size() == 2U )
        {
          // 1. case handled already: if all not-NOT boundary flags are IRREGULAR so should be the element 
          // 2. if one of the 2 flags is IRREGULAR, it is removed
          if ( eflags.count( IRREGULAR ) ) {
              eflags.erase( IRREGULAR );
              (*it)->AtBoundary( (*eflags.begin()) );
            }
          else {
              // 12 possibilities relating to edge elements
              // NOTE: the order of the flags is important here !
              auto sit = eflags.begin();
              const BOX_BOUNDARY flag1 = (*sit); sit++;
              const BOX_BOUNDARY flag2 = (*sit);
              // BACK_BOTTOM
              if ( flag1 == BACK and flag2 == BOTTOM ) (*it)->AtBoundary( EDGE1 );
              // BACK_RIGHT
              else if ( flag1 == BACK and flag2 == RIGHT ) (*it)->AtBoundary( EDGE2 );
              // BACK_TOP
              else if ( flag1 == BACK and flag2 == TOP ) (*it)->AtBoundary( EDGE3 );
              // BACK_LEFT
              else if ( flag1 == BACK and flag2 == LEFT ) (*it)->AtBoundary( EDGE4 );
  
              // BOTTOM_LEFT
              else if ( flag1 == LEFT  and flag2 == BOTTOM ) (*it)->AtBoundary( EDGE5 );
              // BOTTOM_RIGHT
              else if ( flag1 == RIGHT and flag2 == BOTTOM ) (*it)->AtBoundary( EDGE6 );
              // TOP_RIGHT
              else if ( flag1 == RIGHT and flag2 == TOP ) (*it)->AtBoundary( EDGE7 );
              // TOP_LEFT
              else if ( flag1 == LEFT  and flag2 == TOP ) (*it)->AtBoundary( EDGE8 );
              
              // FRONT_BOTTOM
              else if ( flag1 == FRONT and flag2 == BOTTOM ) (*it)->AtBoundary( EDGE9 );
              // FRONT_RIGHT
              else if ( flag1 == FRONT and flag2 == RIGHT ) (*it)->AtBoundary( EDGE10 );
              // FRONT_TOP
              else if ( flag1 == FRONT and flag2 == TOP ) (*it)->AtBoundary( EDGE11 );
              // FRONT_LEFT
              else if ( flag1 == FRONT and flag2 == LEFT ) (*it)->AtBoundary( EDGE12 );
              
              // if a corner is contained that corner flag is choosen
              else if ( flag1 <= CNR1 and flag1 >= CNR8 ) (*it)->AtBoundary( flag1 );
              else if ( flag2 <= CNR1 and flag2 >= CNR8 ) (*it)->AtBoundary( flag2 );
              // if the first integer entry in the set is an edge, that flag is chosen
              else if ( flag1 <= EDGE1 and flag1 > INTERNAL ) (*it)->AtBoundary( flag1 );
              else if ( flag2 <= EDGE1 and flag2 > INTERNAL ) (*it)->AtBoundary( flag2 );
              else  cerr << "\n\tmissed case: " << parseBoundary( flag1 ) << ", " << parseBoundary( flag2 ) << endl;
          }
        }
      else if ( eflags.size() > 2U ) 
        {
          if ( dim == 2 && isTriangular( (*it)->FE_Type() ) ) {
               (*it)->Out();
               csmp_error.notice( ERROR, "flagElementsUsingNodal_BOX_BOUNDARY_Flags:", "all nodes of 2D triangular element appear to be located on boundary.");
            }
          // case when all boundary flags are the same was already considered
            
          // erase the basic boundary options
          eflags.erase( LEFT );
          eflags.erase( RIGHT );
          eflags.erase( TOP );
          eflags.erase( BOTTOM );
          eflags.erase( BACK );
          eflags.erase( FRONT );

          if ( !eflags.empty() ) 
            {
              const BOX_BOUNDARY flag_min = (*min_element( eflags.begin(), eflags.end() ));
              const BOX_BOUNDARY flag_max = (*max_element( eflags.begin(), eflags.end() ));

              // if a corner is contained that corner flag is choosen
              if ( flag_max <= CNR1 and flag_max >= CNR8 ) (*it)->AtBoundary( flag_max );
              else if ( flag_min <= CNR1 and flag_min >= CNR8 ) (*it)->AtBoundary( flag_min );
              // if the first integer entry in the set is an edge, that flag is chosen
              else if ( flag_max <= EDGE1 and flag_max > INTERNAL ) (*it)->AtBoundary( flag_max );
              else if ( flag_min <= EDGE1 and flag_min > INTERNAL ) (*it)->AtBoundary( flag_min );
              // if a corner is contained that corner flag is choosen
              else if ( (*eflags.begin()) <= CNR1 and
                        (*eflags.begin()) >= CNR8 ) (*it)->AtBoundary( (*eflags.begin()) );
              // if the first integer entry in the set is an edge, that flag is chosen
              else if ( (*eflags.begin()) >= EDGE1 and
                        (*eflags.begin()) <  INTERNAL ) (*it)->AtBoundary( (*eflags.begin()) );
              else if ( (*eflags.begin()) == IRREGULAR ) (*it)->AtBoundary( IRREGULAR );
              else {
                  cerr << "\n\n\nflagElementsUsingNodal_BOX_BOUNDARY_Flags: unable to determine box boundary flag for element:\n";
                  for ( auto boundary : eflags )
                    cerr << parseBoundary( boundary ) << " ";
                  cerr << endl;
                  (*it)->Out();
                }
             }
           else  {
                cerr << "\n\tmissed case: ";
                (*it)->Out();
                (*it)->AtBoundary( MULTIPLE );
                csmp_error.notice( WARNING, "flagElementsUsingNodal_BOX_BOUNDARY_Flags:", "could not discern which boundary elemnt is located on.");
             }
        }
      ++it;
    }
    
   if ( nodes_at_boundary == 0 ) 
     csmp_error.notice( ERROR, "flagElementsUsingNodal_BOX_BOUNDARY_Flags:", "none of the nodes in the model had BOX_BOUNDARY flags.");

} // end flagElementUsingNodalAtBoundaryFlags


template void flagElementsUsingNodal_BOX_BOUNDARY_Flags<1U>( typename std::deque<csmp::Element<1U>*>::iterator, typename std::deque<csmp::Element<1U>*>::iterator );
template void flagElementsUsingNodal_BOX_BOUNDARY_Flags<2U>( typename std::deque<csmp::Element<2U>*>::iterator, typename std::deque<csmp::Element<2U>*>::iterator );
template void flagElementsUsingNodal_BOX_BOUNDARY_Flags<3U>( typename std::deque<csmp::Element<3U>*>::iterator, typename std::deque<csmp::Element<3U>*>::iterator );
*/





/**
returns true if the model contains the complete set of box boundary identifiers
all nodes are considered.
*/
bool isStrictlyBoxShaped( const Model<2U>& model )
{
  SmallSet<BOX_BOUNDARY> flags2d, flags_of_model;
  flags2d.insert( LEFT );
  flags2d.insert( RIGHT );
  flags2d.insert( TOP );
  flags2d.insert( BOTTOM );
  flags2d.insert( CNR1 );
  flags2d.insert( CNR2 );
  flags2d.insert( CNR3 );
  flags2d.insert( CNR4 );

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
  SmallSet<BOX_BOUNDARY> flags3d, flags_of_model;
  flags3d.insert( LEFT );
  flags3d.insert( RIGHT );
  flags3d.insert( TOP );
  flags3d.insert( BOTTOM );
  flags3d.insert( FRONT );
  flags3d.insert( BACK );

  flags3d.insert( CNR1 );
  flags3d.insert( CNR2 );
  flags3d.insert( CNR3 );
  flags3d.insert( CNR4 );
  flags3d.insert( CNR5 );
  flags3d.insert( CNR6 );
  flags3d.insert( CNR7 );
  flags3d.insert( CNR8 );

  flags3d.insert( EDGE1 );
  flags3d.insert( EDGE2 );
  flags3d.insert( EDGE3 );
  flags3d.insert( EDGE4 );
  flags3d.insert( EDGE5 );
  flags3d.insert( EDGE6 );
  flags3d.insert( EDGE7 );
  flags3d.insert( EDGE8 );
  flags3d.insert( EDGE9 );
  flags3d.insert( EDGE10 );
  flags3d.insert( EDGE11 );
  flags3d.insert( EDGE12 );

  const Region<3>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    if ( (*it)->AtBoundary() != NOT )
      flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags2d in flags_of_model
  auto it = flags_of_model.begin();
  for ( auto i = flags3d.begin(); i != flags3d.end() && it != flags_of_model.end(); ++i )
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
returns true if the model contains all the side boundary identifiers of the box
*/
bool hasAllSideBoundaries( const Model<3U>& model )
{
  SmallSet<BOX_BOUNDARY> flags3d, flags_of_model;
  flags3d.insert( LEFT );
  flags3d.insert( RIGHT );
  flags3d.insert( TOP );
  flags3d.insert( BOTTOM );
  flags3d.insert( FRONT );
  flags3d.insert( BACK );

  const Region<3>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    if ( (*it)->AtBoundary() != NOT )
      flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags2d in flags_of_model
  auto it = flags_of_model.begin();
  for ( auto i = flags3d.begin(); i != flags3d.end() && it != flags_of_model.end(); ++i )
  {
    it = std::lower_bound( it, flags_of_model.end(), (*i) );
    // make sure the found item is a match
    if ( it != flags_of_model.end() && *i < *it )
      it = flags_of_model.end(); // break out early
  }
  if ( it != flags_of_model.end() ) return true;

  return false;
}

/// 2D side boudaries only
bool hasAllSideBoundaries( const Model<2U>& model )
{
  SmallSet<BOX_BOUNDARY> flags3d, flags_of_model;
  flags3d.insert( LEFT );
  flags3d.insert( RIGHT );
  flags3d.insert( TOP );
  flags3d.insert( BOTTOM );

  const Region<2>&  model_domain( model.Region( "Model" ) );
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    if ( (*it)->AtBoundary() != NOT )
      flags_of_model.insert( (*it)->AtBoundary() );

  // searching for flags2d in flags_of_model
  auto it = flags_of_model.begin();
  for ( auto i = flags3d.begin(); i != flags3d.end() && it != flags_of_model.end(); ++i )
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
template<size_t dim>
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
    (*nit)->Store( nprop_key, makeScalar( ANY, static_cast<double64>((*nit)->AtBoundary()) ) );

  for ( auto eit = mregion.ElementsBegin(); eit != mregion.ElementsEnd(); eit++ )
    (*eit)->Store( eprop_key, makeScalar( ANY, static_cast<double64>((*eit)->AtBoundary()) ) );

} // end boxFlagsToVariable

template void boxFlagsToVariable( Model<2>&, const char*, const char* );
template void boxFlagsToVariable( Model<3>&, const char*, const char* );





template<size_t dim>
void boxBoundaryPropertyRange( const Model<dim>& sg, BOX_BOUNDARY boundary,
                               const char* node_property, double64& bmin, double64& bmax )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  const csmp::Index  prop_key = sg.Database().StorageKey( node_property );
  
  if ( prop_key.place != NODE || prop_key.type != SCALAR )
    csmp_error.notice( ERROR, "boxBoundaryPropertyRange",
                      "Property must a scalar placed on the nodes; nothing was done.");

  bool  first_value( true ), verbose( false );

  const csmp::Region<dim>& sgref = sg.Region( "Model" );

  // measuring the property ranges
  for ( typename vector<Node<dim>*>::const_iterator
        nit = sgref.PerimeterNodesBegin(); nit != sgref.NodesEnd(); nit++ )
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

  if (side1 == LEFT and side2 == BOTTOM) return EDGE1;
  if (side1 == LEFT and side2 == RIGHT)  return EDGE2;
  if (side1 == LEFT and side2 == TOP)    return EDGE3;
  if (side1 == LEFT and side2 == FRONT)  return EDGE4;

  if (side1 == BOTTOM and side2 == FRONT) return EDGE5;
  if (side1 == BOTTOM and side2 == LEFT) return EDGE1;
  if (side1 == BOTTOM and side2 == BACK) return EDGE6;
  if (side1 == BOTTOM and side2 == TOP) return EDGE7;

  if (side1 == TOP and side2 == FRONT) return EDGE8;
  if (side1 == TOP and side2 == LEFT) return EDGE3;
  if (side1 == TOP and side2 == BACK) return EDGE7;
  if (side1 == TOP and side2 == RIGHT) return EDGE12;

  if (side1 == FRONT and side2 == BOTTOM) return EDGE5;
  if (side1 == FRONT and side2 == LEFT) return EDGE4;
  if (side1 == FRONT and side2 == TOP) return EDGE8;
  if (side1 == FRONT and side2 == RIGHT) return EDGE12;

  if (side1 == RIGHT and side2 == BOTTOM) return EDGE9;
  if (side1 == RIGHT and side2 == BACK) return EDGE10;
  if (side1 == RIGHT and side2 == TOP) return EDGE11;
  if (side1 == RIGHT and side2 == FRONT) return EDGE12;
  return NOT;
}


BOX_BOUNDARY  whichCorner( BOX_BOUNDARY side1, BOX_BOUNDARY side2, BOX_BOUNDARY side3 )
{
  set<BOX_BOUNDARY> sides({ side1,side2,side3 });
  auto s3 = (*sides.begin());
  auto s2 = (*(next(sides.begin(), 1)));
  auto s1 = (*sides.rbegin());

  if (sides.empty()) return NOT;

  if (sides.size() == 1U) return (*sides.begin());

  // dealing with duplicates (there will only be 1 or 2 entries so this can only be a corner if the entries are edges)
  if (sides.size() == 2U) {
    // 2 edges EDGE12, EDGE11... EDGE1
    if (side1 == EDGE12 and side2 == EDGE11) return CNR8;
    if (side1 == EDGE12 and side2 == EDGE9) return CNR5;

    if (side1 == EDGE11 and side2 == EDGE10) return CNR7;

    if (side1 == EDGE10 and side2 == EDGE9) return CNR6;
    if (side1 == EDGE10 and side2 == EDGE6) return CNR6;

    if (side1 == EDGE9 and side2 == EDGE6) return CNR6;
    if (side1 == EDGE9 and side2 == EDGE5) return CNR5;

    if (side1 == EDGE8 and side2 == EDGE4) return CNR4;
    if (side1 == EDGE8 and side2 == EDGE3) return CNR4;

    if (side1 == EDGE7 and side2 == EDGE3) return CNR3;
    if (side1 == EDGE7 and side2 == EDGE2) return CNR3;

    if (side1 == EDGE6 and side2 == EDGE2) return CNR2;
    if (side1 == EDGE6 and side2 == EDGE1) return CNR2;

    if (side1 == EDGE5 and side2 == EDGE4) return CNR1;
    if (side1 == EDGE5 and side2 == EDGE1) return CNR1;

    if (side1 == EDGE4 and side2 == EDGE3) return CNR4;
    if (side1 == EDGE4 and side2 == EDGE1) return CNR1;

    if (side1 == EDGE3 and side2 == EDGE2) return CNR3;

    if (side1 == EDGE2 and side2 == EDGE1) return CNR2;

    // not a corner
    if (side1 == BACK and side2 == TOP) return EDGE7;
    if (side1 == BACK and side2 == BOTTOM) return EDGE6;
    if (side1 == BACK and side2 == RIGHT) return EDGE10;
    if (side1 == BACK and side2 == LEFT) return EDGE2;

    if (side1 == FRONT and side2 == TOP) return EDGE8;
    if (side1 == FRONT and side2 == BOTTOM) return EDGE5;
    if (side1 == FRONT and side2 == RIGHT) return EDGE12;
    if (side1 == FRONT and side2 == LEFT) return EDGE4;

    if (side1 == TOP and side2 == RIGHT) return EDGE11;
    if (side1 == TOP and side2 == LEFT) return EDGE3;

    if (side1 == BOTTOM and side2 == LEFT) return EDGE1;
    if (side1 == BOTTOM and side2 == RIGHT) return EDGE9;
  }

  // obeying increasing value constraint: BACK, FRONT, TOP, BOTTOM, RIGHT, LEFT
  if (s1 == FRONT and s2 == BOTTOM and s3 == LEFT) return CNR1;
  if (s1 == BACK  and s2 == BOTTOM and s3 == LEFT) return CNR2;
  if (s1 == BACK  and s2 == TOP    and s3 == BOTTOM) return CNR3;
  if (s1 == FRONT and s2 == TOP    and s3 == LEFT) return CNR4;

  if (s1 == FRONT and s2 == BOTTOM and s3 == RIGHT) return CNR5;
  if (s1 == BACK  and s2 == BOTTOM and s3 == RIGHT) return CNR6;
  if (s1 == BACK  and s2 == TOP    and s3 == RIGHT) return CNR7;
  if (s1 == FRONT and s2 == TOP    and s3 == RIGHT) return CNR5;

  return NOT;
}


template void boxBoundaryPropertyRange( const Model<1U>&, BOX_BOUNDARY, const char*, double64&, double64& );
template void boxBoundaryPropertyRange( const Model<2U>&, BOX_BOUNDARY, const char*, double64&, double64& );
template void boxBoundaryPropertyRange( const Model<3U>&, BOX_BOUNDARY, const char*, double64&, double64& );





} // end namespace csmp
