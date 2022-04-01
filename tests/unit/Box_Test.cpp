#include "Box_Test.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "vsetMakers.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

/**
     Tests the functionality of Box class and the flagging of related box-shaped models 
     
     @todo SKM add testing in 2D
*/
void Box_Test::run()
{
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
  Box().UnitNormalTo( EDGE1, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == sin(45.) );
  _test( unitNormal.at(2) == sin(45.) );
  Box().UnitNormalTo( EDGE2, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sin(45.) );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == sin(45.) );
  Box().UnitNormalTo( EDGE3, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -sin(45.) );
  _test( unitNormal.at(2) == sin(45.) );
  Box().UnitNormalTo( EDGE4, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sin(45.) );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == sin(45.) );
  Box().UnitNormalTo( EDGE5, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sin(45.) );
  _test( unitNormal.at(1) == sin(45.) );
  _test( unitNormal.at(2) == 0 );
  Box().UnitNormalTo( EDGE6, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sin(45.) );
  _test( unitNormal.at(1) == sin(45.) );
  _test( unitNormal.at(2) == 0 );
  Box().UnitNormalTo( EDGE7, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sin(45.) );
  _test( unitNormal.at(1) == -sin(45.) );
  _test( unitNormal.at(2) == 0 );
  Box().UnitNormalTo( EDGE8, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sin(45.) );
  _test( unitNormal.at(1) == -sin(45.) );
  _test( unitNormal.at(2) == 0 );
   Box().UnitNormalTo( EDGE9, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -sin(45.) );
  _test( unitNormal.at(2) == sin(45.) );
  Box().UnitNormalTo( EDGE10, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == -sin(45.) );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == -sin(45.) );
  Box().UnitNormalTo( EDGE11, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == 0 );
  _test( unitNormal.at(1) == -sin(45.) );
  _test( unitNormal.at(2) == -sin(45.) );
  Box().UnitNormalTo( EDGE12, 3, unitNormal );
  _test( unitNormal.size() == 3 );
  _test( unitNormal.at(0) == sin(45.) );
  _test( unitNormal.at(1) == 0 );
  _test( unitNormal.at(2) == -sin(45.) );
  
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

  _test( TestWhetherSideBoundaryFlagsArePresent() );
  
  // _test( TestWhetherAllBoxFlagsArePresent() ); // fails because the corners are missing
  _test( TestBoundaryFlagging() );
  _test( TestWhetherBoundaryFlagsArePreservedInBinaryFile1() );
  
  // tests whether the function recreateBoxBoundaryFlags() manages to reconstruct edges and boundaries correctly
  _test( TestBoundaryFlagRecreation() );
  
  // tests whether the boundaries of a 2D box model are assigned correctly
  _test( TestBoundaryFlagAssigment2D() );

  TestWhetherSimplexNormalsAreOutwardPointing();

} // end run



bool Box_Test::TestBoundaryFlagAssigment2D()
 {
    // 1. standard model construction
    // ------------------------------
    VSet<2U> vset;
    test_Create_TrianglePatch_VSet( vset );
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
    const bool irregular_mesh(true);
    const bool binary_file(true);
    const bool use_regions_file(true);
    const bool debug(true);
   
    ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );
   
    string node_variable("nodal box flag"), elmt_variable("element box flag");
    model.CreateProperty( node_variable.c_str(), "flag", SCALAR, NODE );
//    model.UpdateIndices();
    model.CreateProperty( elmt_variable.c_str(), "flag", SCALAR, ELEMENT );
//    model.UpdateIndices();
   
    boxFlagsToVariable( model, node_variable.c_str(), elmt_variable.c_str() );

    if ( debug ) {
         VTU_Interface<3>  vtu(model);
         std::list<std::string> outputProps;
         outputProps.push_back(node_variable.c_str());
         outputProps.push_back(elmt_variable.c_str());
         vtu.OutputDataToVTU( "box_flags", outputProps, "Model", 1 );
      }

    return true;
  
 } // end TestBoundaryFlagging




/**
    Comparing the unit normals with the those of the sides 
    of the box-shaped model. 
    The expectation is that are are pointing in the direction
    as the normals of the elements or faces on the outside boundary
    of the model.
*/
void Box_Test::TestWhetherSimplexNormalsAreOutwardPointing()
 {
    const bool irregular_mesh(true);
    const bool binary_file(true);
    const bool use_regions_file(true);
   
    ANSYS_Model3D model( "prism_test", "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );
    Region<3U>    model_domain(model.Region("Model"));

    // 1. testing the unit normals of the (volumetric elements)
    // --------------------------------------------------------
    vector<double> leftNormal, rightNormal, topNormal, bottomNormal, frontNormal, backNormal, eUnitNormal;
    Box().UnitNormalTo( LEFT,   3, leftNormal );
    Box().UnitNormalTo( RIGHT,  3, rightNormal );
    Box().UnitNormalTo( TOP,    3, topNormal );
    Box().UnitNormalTo( BOTTOM, 3, bottomNormal );
    Box().UnitNormalTo( FRONT,  3, frontNormal );
    Box().UnitNormalTo( BACK,   3, backNormal );

    for ( size_t i=model_domain.InteriorCells(); i<model_domain.Cells(); ++i ) {
         for ( auto j=0U; j<model_domain.PerimeterFaces(i); ++j ) {
                const BOX_BOUNDARY flag = model_domain.E(i)->AtBoundary(j);
                assert( flag != NOT );
                // verifying alignment of the element's unit normal with that of the model boundary
                model_domain.E(i)->UnitNormalToFace( model_domain.PerimeterFace(i,j), eUnitNormal );
                // checking whether the normals are aligned and of of same unit magnitude
                if ( flag == LEFT )  {
                     const double dotProduct(vector_product<3U,double>(leftNormal,eUnitNormal));
                     // testing for alignment
                     _test( dotProduct > 0. );
                     // testing for unit length
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
                  }
                else if ( flag == RIGHT )  {
                     const double dotProduct(vector_product<3U,double>(rightNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
                  }
                else if ( flag == BOTTOM )  {
                     const double dotProduct(vector_product<3U,double>(bottomNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
                  }
                else if ( flag == TOP )  {
                     const double dotProduct(vector_product<3U,double>(topNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
                  }
                else if ( flag == BACK )  {
                     const double dotProduct(vector_product<3U,double>(backNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
                  }
                else if ( flag == FRONT )  {
                     const double dotProduct(vector_product<3U,double>(frontNormal,eUnitNormal));
                     _test( dotProduct > 0. );
                     _equal( dotProduct, 1., numeric_limits<double>::epsilon() * 5. );
                  }
            }
      }
   
   
    // 2. testing whether the unit normals of the faces making up the outside boundaries of the model
    //    are outward pointing and aligned
    // -----------------------------------
    Boundary<3U>  left(model.Boundary("LEFT"));
    for ( auto it=left.CellsBegin(); it!=left.CellsEnd(); ++it ) {
         if ( (*it)->IsSurfaceElement() ) {
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
         if ( (*it)->IsSurfaceElement() ) {
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
         if ( (*it)->IsSurfaceElement() ) {
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
         if ( (*it)->IsSurfaceElement() ) {
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
         if ( (*it)->IsSurfaceElement() ) {
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
         if ( (*it)->IsSurfaceElement() ) {
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
    const bool irregular_mesh(false);
    const bool binary_file(true);
    const bool use_regions_file(true);
    ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );
    return hasAllSideBoundaries( model );
}



/** 
    All flags for 3D model including edges and corners
*/
bool Box_Test::TestWhetherAllBoxFlagsArePresent()
 {
    const bool irregular_mesh(true);
    const bool binary_file(true);
    const bool use_regions_file(true);
    ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );
    return isStrictlyBoxShaped( model );
}
  

/**
   tests method which recreates boundary flags
*/
bool Box_Test::TestBoundaryFlagRecreation()
 {
    const bool irregular_mesh(true);
    const bool binary_file(true);
    const bool use_regions_file(true);
   
    // build a model with valid boundaries
    ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );

    // testing whether the reflagging works correctly
    recreateBoxBoundaryFlags( model );

    return isStrictlyBoxShaped( model );
 }



/** 
    writes model to CSMP binary and then re-reads it to see whether the box flags survive
*/
bool Box_Test::TestWhetherBoundaryFlagsArePreservedInBinaryFile()
 {
    const bool irregular_mesh(true);
    const bool binary_file(true);
    const bool use_regions_file(true);

    ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );
    const Region<3>&  mregion(model.Region("Model"));
   
    // storing the flags in node and element order in a list for comparison
    list<BOX_BOUNDARY>  node_flags_before;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      node_flags_before.push_back( (*nit)->AtBoundary() );
    list<BOX_BOUNDARY>  elmt_flags_before;
    for ( auto eit=mregion.CellsBegin(); eit!=mregion.CellsEnd(); eit++ )
      for ( auto i{0}; i<(*eit)->Neighbors(); ++i )
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
      for ( auto i{0}; i<(*eit)->Neighbors(); ++i )
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
    const bool irregular_mesh(true);
    const bool binary_file(true);
    const bool use_regions_file(true);

    ANSYS_Model3D model( model_name_.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file );
    const Region<3>&  mregion(model.Region("Model"));
   
    // storing the flags in node and element order in a list for comparison
    list<BOX_BOUNDARY>  node_flags_before;
    for ( auto nit=mregion.NodesBegin(); nit!=mregion.NodesEnd(); nit++ )
      node_flags_before.push_back( (*nit)->AtBoundary() );
    list<BOX_BOUNDARY>  elmt_flags_before;
    for ( auto eit=mregion.CellsBegin(); eit!=mregion.CellsEnd(); eit++ )
      for ( auto i{0}; i<(*eit)->Neighbors(); ++i )
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
      for ( auto i{0}; i<(*eit)->Neighbors(); ++i )
        elmt_flags_after.push_back( (*eit)->AtBoundary(i) );

    delete mptr;
   
    // comparing node and element flags of the original and the restored model
    if ( node_flags_after != node_flags_before ) return false;
    if ( elmt_flags_after != elmt_flags_before ) return false;
    return true;
 }




} // csmp


















