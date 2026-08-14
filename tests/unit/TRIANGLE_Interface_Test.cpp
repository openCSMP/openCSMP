//
//  TRIANGLE_Interface_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 14/11/2024.
//

#include "TRIANGLE_Interface_Test.h"
#include "TRIANGLE_Interface.h"
#include "Model.h"
#include "Boundary.h"
#include "VSet.h"
#include "Box.h"

using namespace std;

namespace csmp {

void TRIANGLE_Interface_Test::run()
 {
     Test_BOX_BOUNDARY_NodeFlags();
 }
 
void TRIANGLE_Interface_Test::Test_BOX_BOUNDARY_NodeFlags()
 {
    // ---------------------------------------------------------------------------------------
    // 0. Import model from Shewchuk's Triangle mesher
    // ---------------------------------------------------------------------------------------
    TRIANGLE_Interface  mesh_interface;
    VSet<2U>            mesh_container;
//    string              file_name = "irish1x.1";
    string              file_name = "blunt30deg.1";
    const bool isoparametric{true};
    mesh_interface.ReadTriangle2DMesh( file_name.c_str(), mesh_container, isoparametric );
//    establishBoundaryFlagsForBoxModel( mesh_container, 0.1 );

    Model<2U>  model( mesh_container, "Minimum-variables.txt" );
    mesh_container.Erase();

    // create csmp::Boundary objects to test whether the right corners are contained each of the boundaries
    // (EstablishBoxBoundariesFromOrientation will not change flags that were set earlier)
    model.EstablishBoxBoundariesFromOrientation();
    
    if ( verbose_ ) cout <<"\n"<<"TRIANGLE_Interface_Test(verbose): checking validity of boundary flags:"<< endl;
    for ( auto bit=model.BoundariesBegin(); bit!=model.BoundariesEnd(); ++bit ) {
         if ( verbose_ ) {
             cout <<"\n"<<"'"<< (*bit).first <<"', with boudary flags: "<< endl;
             for ( auto b : (*bit).second.NodeVector() )
               cout <<"  "<< parseBoundary(b->AtBoundary());
           }
         if ( (*bit).first == "BOTTOM" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              _test( flags.find(CNR1) != flags.end() );
              _test( flags.find(CNR2) != flags.end() );
           }
         else if ( (*bit).first == "RIGHT" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              _test( flags.find(CNR2) != flags.end() );
              _test( flags.find(CNR3) != flags.end() );
           }
         else if ( (*bit).first == "TOP" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              _test( flags.find(CNR3) != flags.end() );
              _test( flags.find(CNR4) != flags.end() );
           }
         else if ( (*bit).first == "LEFT" ) {
              set<BOX_BOUNDARY> flags;
              for ( auto b : (*bit).second.NodeVector() )
                if ( b->AtBoundary() != NOT )
                  flags.insert( b->AtBoundary() );
              _test( flags.find(CNR1) != flags.end() );
              _test( flags.find(CNR4) != flags.end() );
           }
       }
    cout << endl;

//    model.InputPropertyValue( "fluid pressure", makeScalar(ANY,0.) );
//    model.InputPropertyValue( "permeability", makeScalar(ANY,0.) );
//    boxFlagsToVariable( model, "fluid pressure", "permeability" );
    
} // end Test_BOX_BOUNDARY_NodeFlags



} // end csmp
