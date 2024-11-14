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
    string              file_name = "irish1x.1";
//    cout <<"\nmain: Enter name of 'Triangle' input file set: ";
//    cin >> file_name;
    const bool isoparametric{true};
    mesh_interface.ReadTriangle2DMesh( file_name.c_str(), mesh_container, isoparametric );
//    auto x_range = mesh_container.X_Range();
//    auto y_range = mesh_container.Y_Range();
//    const double tolerance = distance( Point<2>(x_range.second,y_range.second), Point<2>(x_range.first,y_range.first) );
//    flagCornerNodes( mesh_container, tolerance, x_range.first, x_range.second, y_range.first, y_range.second, 0., 0. );

    // 'UG4_ProMeshOutput_Example-variables.txt' is the text file that defines the variables used in this example
    Model<2U>  model( mesh_container, "CSMP-variables.txt" );
    mesh_container.Erase();
//    printModelDimensions( model, true );
//    printRangeOfVariable( model, "permeability", true );


    // create csmp::Boundary objects to test whether the right corners are contained each of the boundaries
    // (EstablishBoxBoundariesFromOrientation will not change flags that were set earlier)
    model.EstablishBoxBoundariesFromOrientation();
    
    if ( verbose_ ) cout <<"\nmain: checking validity of boundary flags:"<< endl;
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
