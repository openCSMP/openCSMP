#include "SplitBoundaryInterface_Test.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "VTU_Interface.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "EclipseModel.h"
#include "VTK_Interface.h"
#include "compareFloats.h"

using namespace std;

namespace csmp
{





//E.P Added tests which rigorously check the nodes, elements and interfaces have been correctly split and updated on inside and outside of the splitboundary.
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::Test_splitboundary_from_lower_dim_region() {

    if constexpr (dim == 2){
        Test_NodeCorrespondance_2D("InternalBoundary_test");
    } else if constexpr (dim == 3){
        Test_NodeAndElementsCorrespondance_3D("InternalBoundary3D_test");
        Test_NodeAndElementsCorrespondance_3D_X_Intersection("InternalBoundary3D_Intersect_Test");
        Test_NodeAndElementsCorrespondance_3D_X_Intersection_Reverse("InternalBoundary3D_Intersect_Test");
    }
}




/**
 * A Semi Rigorous test: Test that nodes match on either side of splitboundary (Basic test, more rigorous to follow).
 *
 * 0) Tests nodes are matching when using Interface::MatchingN(i,side)
 * 1) Test Unit Normal
**/
template<uint32_t dim>
bool SplitBoundaryInterface_Test<dim>::Test_NodeCorrespondance_2D( const char* mesh_file ) {

  const uint32_t DIM{2U};
  int32_t material_id = 1;
  //model construction
  const char* variables_file("SplitBoundary_Test-variables.txt");
  const char* regions_file("InternalBoundary_Test");

// TODO: replace these ANSYS models with CSMP native models that we construct once and then hand over to the tests
  ANSYS_Model2D model1(mesh_file, regions_file, variables_file, false, true, true);

  ///Testing Differnt model creations
  //Region -> SplitBoundary
  string sb_name1  = *((model1.CreateSplitBoundaryFrom( "FRACTURE" ) ).first.begin()) ;

  ///Inserting lower dimensional region in each
  model1.InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id );

  //Getting splitboundaries
  SplitBoundary<DIM>& sb1 = model1.SplitBoundary( sb_name1 );

  std::vector<SplitBoundary<DIM>> sb_vec{sb1};
  for ( auto& sb : sb_vec){
    for (auto& ifp : sb.CellVector() ){
      uint32_t n_nodes = ifp->FE()->Nodes();

      std::vector<uint32_t> nids_in = ifp->InnerParent()->FE()->NodesOfFace( ifp->InnerParentFaceID() );
      std::vector<uint32_t> nids_out = ifp->OuterParent()->FE()->NodesOfFace( ifp->OuterParentFaceID() );

      for ( uint32_t n{0U}; n<n_nodes;++n){
        //Nodes should match that of face
        _test( ifp->N(n,INSIDE)  == ifp->InnerParent()->N(nids_in[n] ));
        //_test( ifp->N(n,OUTSIDE) == ifp->OuterParent()->N(nids_out[n])) ; The nodes of face of OUTSIDE Volumetric Element are not collocated with the Outside of Interface node_connector.

        //Testing matching assignment
        //inside remains unchanged
        _test( ifp->MatchingN(n,INSIDE) == ifp->N(n,INSIDE) );

        //Outside coordinates must match with Inside and middle
        _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,OUTSIDE)->Coordinate() );
        _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,MIDDLE)->Coordinate() );
      }
    }

  }

  sb1.PullApartSplitBoundary(0.1);

  //Uncomment to visualise
  /*
  list<string> outputProps;
  outputProps.push_back( "nodal id" );
  model1.InputPropertyValue("nodal id", ScalarVariable(PLAIN,1.0));
  VTU_Interface<DIM> vtu1( model1 );
  vtu1.OmitZeroInFileName(true);
  vtu1.OutputDataToVTU( "../Output/InternalBoundary_TestA", outputProps, "Model", static_cast<int>(0) );
  */

  return true;

}






/**
 * A VERY RIGOROUS TEST: Test that splitboundary is properly created from Lowerdim Region. Intersections are not handled.
 *
 * 0) Tests nodes are matching when using Interface::MatchingN(i,side)
 * 1) Test Unit normals
 * 2) Test Inside nodes of Interface match exactly with the face of the Inner Parent nodes (and not slightly rotated, which is not good enough).
 * 3) Tests all OUTSIDE elements parents have the CORRECT new nodes on OUTSIDE, including elements that do not share a face with the interface
 *
*/
template<uint32_t dim>
bool SplitBoundaryInterface_Test<dim>::Test_NodeAndElementsCorrespondance_3D( const char* mesh_file){

  const uint32_t DIM{3U};       //we do not use template parameter - thats stupid in a test...
  int32_t material_id = 1;
  //model construction
  const char* variables_file("SplitBoundary_Test-variables.txt");
  const char* regions_file("InternalBoundary3D_Test");

  bool irregular = false;
  bool reduce_to_regions = true;

  ANSYS_Model3D model1(mesh_file, regions_file, variables_file, irregular, reduce_to_regions);

  //Getting volumetric regions for later use
  Region<DIM>& BottomUnit1 = model1.Region("BOTTOMUNIT");
  Region<DIM>& TopUnit1 = model1.Region("TOPUNIT");

  //Getting perimeter nodes of region, before splitboundary is created and information is lost
  std::set<Node<DIM>*> perim_nodes1, perim_nodes2;
  Region<DIM>& frac_region1 = model1.Region( "FRACTURE" );
  for ( auto nit = frac_region1.PerimeterNodesBegin(); nit !=frac_region1.NodesEnd(); nit++){
    perim_nodes1.insert(*nit);
  }
  _test( perim_nodes1.size() == frac_region1.PerimeterNodes() );


  ///Testing splitboundary creation

  //Region -> SplitBoundary
  string sb_name1  = *((model1.CreateSplitBoundaryFrom( "FRACTURE" ) ).first.begin()) ;

  ///Inserting lower dimensional region in each
  model1.InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id );

  //Getting splitboundaries
  SplitBoundary<DIM>& sb1 = model1.SplitBoundary( sb_name1 );


  std::vector<SplitBoundary<DIM>> sb_vec{sb1};
  for ( auto& sb : sb_vec){
    for (auto& ifp : sb.CellVector() ){

      //Check unit normals are 0 1 0 on INSIDE and MIDDLE and the opposite on OUTSIDE
      Point<DIM> nrml_in  = ifp->UnitNormal(INSIDE);
      Point<DIM> nrml_out = ifp->UnitNormal(OUTSIDE);
      Point<DIM> nrml_mid = ifp->UnitNormal(MIDDLE);
      //Inside opposite to outside
      _equal( nrml_in[0] ,  nrml_out[0] , std::numeric_limits<double>::epsilon()  );
      _equal( nrml_in[1] , -nrml_out[1] , std::numeric_limits<double>::epsilon()  );
      _equal( nrml_in[2] ,  nrml_out[2] , std::numeric_limits<double>::epsilon()  );
      //inside is 0 1 0
      _equal( nrml_in[1] , 1.0 , 0.0001  );
      _test(  std::fabs(nrml_in[0]) < 0.0001 );
      _test(  std::fabs(nrml_in[2]) < 0.0001 );
      //middle is 0 1 0
      _equal( std::fabs(nrml_mid[1]) , 1.0 , 0.0001  );
      _test(  std::fabs(nrml_mid[0]) < 0.0001 );
      _test(  std::fabs(nrml_mid[2]) < 0.0001 );

      //Retrieve the nodes of face that match with INSIDE OUTSIDE
      std::vector<uint32_t> nids_in  = ifp->InnerParent()->FE()->NodesOfFace( ifp->InnerParentFaceID()),
                            nids_out = ifp->OuterParent()->FE()->NodesOfFace( ifp->OuterParentFaceID());


      std::set<Node<DIM>*> outside_nds_interface;
      std::set<Node<DIM>*> outside_nds_face;
      uint32_t n_nodes = ifp->FE()->Nodes();
      for ( uint32_t n{0U}; n<n_nodes;++n){
        //Nodes should match that of INSIDE face
        _test( ifp->N(n,INSIDE)  == ifp->InnerParent()->N(nids_in[n] ));

        //Nodes should contain same nodes of OUTSIDE face, but may be rotated, so may not match
        outside_nds_interface.insert( ifp->N(n,OUTSIDE));
        outside_nds_face.insert( ifp->OuterParent()->N(nids_out[n]));

        //Testing matching assignment
        //inside remains unchanged
        _test( ifp->MatchingN(n,INSIDE) == ifp->N(n,INSIDE) );
        //Outside coordinates must match with Inside and middle
        _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,OUTSIDE)->Coordinate() );
        _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,MIDDLE)->Coordinate() );

        //Testing nodes are duplicated in interior and the same on perimeter
        //if node is not on perimeter
        if ( perim_nodes1.find( ifp->N(n,INSIDE) ) == perim_nodes1.end() ){
          //test node is duplicated
          _test( ifp->MatchingN(n,INSIDE) != ifp->MatchingN(n,OUTSIDE) );

          //Inside Node is in bottom half
          _test( BottomUnit1.Contains( ifp->N(n,INSIDE) ) );
          //matching Outside Node is in top half
          _test( TopUnit1.Contains( ifp->MatchingN(n,OUTSIDE) ) );
          //matching outside node not in bottom half
          _test( !BottomUnit1.Contains( ifp->MatchingN(n,OUTSIDE) ) );

          //Testing inside parents assignmnet by checking they are part of BottomUnit region
          for (uint32_t p{0U} ; p < ifp->N(n,INSIDE)->Parents() ; ++p){
            _test( BottomUnit1.Contains( ifp->N(n,INSIDE)->Parent(p) )); //bottom unit has parent of inside node
            _test( !TopUnit1.Contains( ifp->N(n,INSIDE)->Parent(p)));    //top unit doenst have parent of inside node
          }
          for (uint32_t p{0U} ; p < ifp->MatchingN(n,OUTSIDE)->Parents() ; ++p){
            _test( !BottomUnit1.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) )); //bottom unit Does Not have parent of matching outside node
            _test( TopUnit1.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)));    //top unit does have parent of matching outside node
          }

        }  else _test( ifp->MatchingN(n,INSIDE) == ifp->MatchingN(n,OUTSIDE) ) ; //test nodes match if on perimeter


      }

      //outside nodes must match in content but may be rotated, so we test just the content
      _test(outside_nds_interface == outside_nds_face);

    }
  }

  //For visualisation purposes
  double gap = 1.0;
  sb1.PullApartSplitBoundary(gap);

  //uncomment to visualise
  /*
  list<string> outputProps;
  outputProps.push_back( "nodal id" );
  model1.InputPropertyValue("nodal id", ScalarVariable(PLAIN,1.0));
  VTU_Interface<DIM> vtu1( model1 );
  vtu1.OmitZeroInFileName(true);
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_TestA", outputProps, "Model", static_cast<int>(0) );
  */

  return true;

}






/**
 * A VERY RIGOROUS TEST: Test that splitboundary is properly created from Lowerdim Region WITH an X Intersection being handled and incrementally split.
 *
 *  This test requires the split boundaries to be partly seperating four Regions (bottom left, bottom right, Top left, Top right) for it to work
 *
 * 0) Tests nodes are matching when using Interface::MatchingN(i,side)
 * 1) Test Unit normals
 * 2) Test Inside nodes of Interface match exactly with the face of the Inner Parent nodes (and not slightly rotated, which is not good enough).
 * 3) Tests all OUTSIDE elements parents have the CORRECT new nodes on OUTSIDE, including elements that do not share a face with the interface
 *
*/
template<uint32_t dim>
bool SplitBoundaryInterface_Test<dim>::Test_NodeAndElementsCorrespondance_3D_X_Intersection( const char* mesh_file){

  const uint32_t DIM{3U};
  int32_t material_id = 1;
  //model construction
  const char* variables_file("SplitBoundary_Test-variables.txt");
  const char* regions_file("InternalBoundary3D_Intersect_Test");

  bool irregular = false;
  bool reduce_to_regions = true;

  ANSYS_Model3D model1(mesh_file, regions_file, variables_file, irregular, reduce_to_regions);

  //Getting volumetric regions for later use
  Region<DIM>& left_back_bottom     = model1.Region("LEFT_BACK_BOTTOM");
  Region<DIM>& left_back_top        = model1.Region("LEFT_BACK_TOP");
  Region<DIM>& right_front_bottom   = model1.Region("RIGHT_FRONT_BOTTOM");
  Region<DIM>& right_front_top      = model1.Region("RIGHT_FRONT_TOP");

  //Getting perimeter nodes of region, before splitboundary is created and information is lost
  std::set<Node<DIM>*> perim_nodes_planar, perim_nodes_diagonal;
  Region<DIM>& frac_region_planar = model1.Region( "FRACTURE_PLANAR" );
  for ( auto nit = frac_region_planar.PerimeterNodesBegin(); nit !=frac_region_planar.NodesEnd(); nit++){
    perim_nodes_planar.insert(*nit);
  }
  _test( perim_nodes_planar.size() == frac_region_planar.PerimeterNodes() );

  Region<DIM>& frac_region_diagonal = model1.Region( "FRACTURE_DIAGONAL" );
  for ( auto nit = frac_region_diagonal.PerimeterNodesBegin(); nit !=frac_region_diagonal.NodesEnd(); nit++){
    perim_nodes_diagonal.insert(*nit);
  }
  _test( perim_nodes_diagonal.size() == frac_region_diagonal.PerimeterNodes() );


  Node<DIM>* split_perimter_node = nullptr;
  for ( auto n : perim_nodes_planar ){
    if ( approximatelyEqual(n->Coordinate()[0] , 2.5, 0.01) &&
         approximatelyEqual(n->Coordinate()[1] , 5.0, 0.01) &&
         approximatelyEqual(n->Coordinate()[2] , 7.5, 0.01)  ){
      split_perimter_node = n;                          //this node is on the perimter of fracture planar, but split by fracture diagonal
    }
  }
  assert(split_perimter_node!=nullptr);

  ///Testing splitboundary creation


  //Region -> SplitBoundary
  std::set<string> sb_names    = (model1.CreateSplitBoundaryFrom( "FRACTURE_PLANAR" ) ).first ;
  std::set<string> sb_names_2  = (model1.CreateSplitBoundaryFrom( "FRACTURE_DIAGONAL" ) ).first ;

  //Need to find node on perimeter that was split by diagonal fracture
  Node<DIM>* second_split_perimeter_node = nullptr;
  assert(split_perimter_node->Manifold()->Branches()==2);
  split_perimter_node->Manifold()->N(0) == split_perimter_node ?  second_split_perimeter_node = split_perimter_node->Manifold()->N(1) : second_split_perimeter_node->Manifold()->N(0);


  std::string sb1 = model1.MergeSplitBoundaries( "FRACTURE_PLANAR",   sb_names );
  std::string sb2 = model1.MergeSplitBoundaries( "FRACTURE_DIAGONAL", sb_names_2);

  //For Planar SplitBoundary
  for (auto ifp : model1.SplitBoundary(sb1).CellVector() ){

    //Check unit normals are 0 1 0 on INSIDE and MIDDLE and the opposite on OUTSIDE
    Point<DIM> nrml_in  = ifp->UnitNormal(INSIDE);
    Point<DIM> nrml_out = ifp->UnitNormal(OUTSIDE);
    //Point<dim> nrml_mid = ifp->UnitNormal(MIDDLE);
    //Inside opposite to outside
    _equal( nrml_in[0] ,  nrml_out[0] , std::numeric_limits<double>::epsilon()  );
    _equal( nrml_in[1] , -nrml_out[1] , std::numeric_limits<double>::epsilon()  );
    _equal( nrml_in[2] ,  nrml_out[2] , std::numeric_limits<double>::epsilon()  );
    //inside is 0 1 0
    _equal( nrml_in[1] , 1.0 , 0.0001  );
    _test(  std::fabs(nrml_in[0]) < 0.0001 );
    _test(  std::fabs(nrml_in[2]) < 0.0001 );

    //Retrieve the nodes of face that match with INSIDE OUTSIDE
    std::vector<uint32_t> nids_in  = ifp->InnerParent()->FE()->NodesOfFace( ifp->InnerParentFaceID()),
                          nids_out = ifp->OuterParent()->FE()->NodesOfFace( ifp->OuterParentFaceID());


    std::set<Node<DIM>*> outside_nds_interface;
    std::set<Node<DIM>*> outside_nds_face;
    uint32_t n_nodes = ifp->FE()->Nodes();
    for ( uint32_t n{0U}; n<n_nodes;++n){
      //Nodes should match that of INSIDE face
      _test( ifp->N(n,INSIDE)  == ifp->InnerParent()->N(nids_in[n] ));

      //Nodes should contain same nodes of OUTSIDE face, but may be rotated, so may not match
      outside_nds_interface.insert( ifp->N(n,OUTSIDE));
      outside_nds_face.insert( ifp->OuterParent()->N(nids_out[n]));

      //Testing matching assignment
      //inside remains unchanged
      _test( ifp->MatchingN(n,INSIDE) == ifp->N(n,INSIDE) );
      //Outside coordinates must match with Inside and middle
      _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,OUTSIDE)->Coordinate() );

      //Testing nodes are duplicated in interior and the same on perimeter
      //if node is not on perimeter
      if ( perim_nodes_planar.find( ifp->N(n,INSIDE) ) == perim_nodes_planar.end() && ifp->N(n,INSIDE) != second_split_perimeter_node  ){
        //test node is duplicated
        _test( ifp->MatchingN(n,INSIDE) != ifp->MatchingN(n,OUTSIDE) );

        //Inside Node is in bottom half
        _test( left_back_bottom.Contains( ifp->N(n,INSIDE) ) || right_front_bottom.Contains( ifp->N(n,INSIDE)) );
        //matching Outside Node is in top half
        _test( left_back_top.Contains( ifp->MatchingN(n,OUTSIDE) ) || right_front_top.Contains(ifp->MatchingN(n,OUTSIDE)) );
        //matching outside node not in bottom half
        _test( !left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE) ) && !right_front_bottom.Contains(ifp->MatchingN(n,OUTSIDE)) );

        //Testing inside parents assignmnet by checking they are part of BottomUnit region
        for (uint32_t p{0U} ; p < ifp->N(n,INSIDE)->Parents() ; ++p){
          if ( ifp->N(n,INSIDE)->Parent(p)->IsEquidimensional()){
            _test( left_back_bottom.Contains( ifp->N(n,INSIDE)->Parent(p)) || right_front_bottom.Contains( ifp->N(n,INSIDE)->Parent(p) )); //bottom unit has parent of inside node
            _test( !left_back_top.Contains( ifp->N(n,INSIDE)->Parent(p)) && !right_front_top.Contains( ifp->N(n,INSIDE)->Parent(p)) );    //top unit doenst have parent of inside node
          }
        }
        for (uint32_t p{0U} ; p < ifp->MatchingN(n,OUTSIDE)->Parents() ; ++p){
          if ( ifp->MatchingN(n,OUTSIDE)->Parent(p)->IsEquidimensional()){
            _test( !left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) ) && !right_front_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) )); //bottom unit Does Not have parent of matching outside node
            _test( left_back_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)) || right_front_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)));    //top unit does have parent of matching outside node
          }
        }

      }  else _test( ifp->MatchingN(n,INSIDE) == ifp->MatchingN(n,OUTSIDE) ) ; //test nodes match if on perimeter


    }

    //outside nodes must match in content but may be rotated, so we test just the content
    _test(outside_nds_interface == outside_nds_face);
  }


  //For Diagonal SplitBoundary
  for (auto ifp : model1.SplitBoundary(sb2).CellVector() ){

    //Check unit normals are 0 1 0 on INSIDE and MIDDLE and the opposite on OUTSIDE
    Point<DIM> nrml_in  = ifp->UnitNormal(INSIDE);
    Point<DIM> nrml_out = ifp->UnitNormal(OUTSIDE);
    //Point<dim> nrml_mid = ifp->UnitNormal(MIDDLE);
    //Inside opposite to outside
    _equal( nrml_in[0] , -nrml_out[0] , 10.0*std::numeric_limits<double>::epsilon()  );
    _test(  std::fabs(nrml_in[1]) < 0.001 );     //normal should be along x-z plane
    _test(  std::fabs(nrml_out[1]) < 0.001 );    //normal should be along x-z plane
    _equal( nrml_in[2] , -nrml_out[2] , 10.0*std::numeric_limits<double>::epsilon()  );
    //inside is -0.707 0 -0.707
    _equal( nrml_in[0] , -1.0/std::sqrt(2) , 0.001  );
    _equal( nrml_in[2] , -1.0/std::sqrt(2) , 0.001  );
    _test(  std::fabs(nrml_in[1]) < 0.001 );

    //Retrieve the nodes of face that match with INSIDE OUTSIDE
    std::vector<uint32_t> nids_in  = ifp->InnerParent()->FE()->NodesOfFace( ifp->InnerParentFaceID()),
                          nids_out = ifp->OuterParent()->FE()->NodesOfFace( ifp->OuterParentFaceID());


    std::set<Node<DIM>*> outside_nds_interface;
    std::set<Node<DIM>*> outside_nds_face;
    uint32_t n_nodes = ifp->FE()->Nodes();
    for ( uint32_t n{0U}; n<n_nodes;++n){
      //Nodes should match that of INSIDE face
      _test( ifp->N(n,INSIDE)  == ifp->InnerParent()->N(nids_in[n] ));

      //Nodes should contain same nodes of OUTSIDE face, but may be rotated, so may not match
      outside_nds_interface.insert( ifp->N(n,OUTSIDE));
      outside_nds_face.insert( ifp->OuterParent()->N(nids_out[n]));

      //Testing matching assignment
      //inside remains unchanged
      _test( ifp->MatchingN(n,INSIDE) == ifp->N(n,INSIDE) );
      //Outside coordinates must match with Inside and middle
      _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,OUTSIDE)->Coordinate() );

      //Testing nodes are duplicated in interior and the same on perimeter
      //if node is not on perimeter
      if ( perim_nodes_diagonal.find( ifp->N(n,INSIDE) ) == perim_nodes_diagonal.end()  ){
        //test node is duplicated
        _test( ifp->MatchingN(n,INSIDE) != ifp->MatchingN(n,OUTSIDE) );

        //Inside Node is in bottom half
        _test( right_front_bottom.Contains( ifp->N(n,INSIDE) ) || right_front_top.Contains( ifp->N(n,INSIDE)) );
        //matching Outside Node is in top half
        _test( left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE) ) || left_back_top.Contains(ifp->MatchingN(n,OUTSIDE)) );
        //matching outside node not in inside diagonal
        _test( !right_front_bottom.Contains( ifp->MatchingN(n,OUTSIDE) ) && !right_front_top.Contains(ifp->MatchingN(n,OUTSIDE)) );
        //inside node not in outside diagonal
        _test( !left_back_bottom.Contains( ifp->MatchingN(n,INSIDE) ) && !left_back_top.Contains(ifp->MatchingN(n,INSIDE)) );

        //Testing inside parents assignmnet by checking they are part of BottomUnit region
        for (uint32_t p{0U} ; p < ifp->N(n,INSIDE)->Parents() ; ++p){
          if ( ifp->N(n,INSIDE)->Parent(p)->IsEquidimensional()){
            _test( right_front_bottom.Contains( ifp->N(n,INSIDE)->Parent(p)) || right_front_top.Contains( ifp->N(n,INSIDE)->Parent(p) )); //bottom unit has parent of inside node
            _test( !left_back_bottom.Contains( ifp->N(n,INSIDE)->Parent(p)) && !left_back_top.Contains( ifp->N(n,INSIDE)->Parent(p)) );    //top unit doenst have parent of inside node
          }
        }
        for (uint32_t p{0U} ; p < ifp->MatchingN(n,OUTSIDE)->Parents() ; ++p){
          if ( ifp->MatchingN(n,OUTSIDE)->Parent(p)->IsEquidimensional()){
            _test( !right_front_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) ) && !right_front_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) )); //bottom unit Does Not have parent of matching outside node
            _test( left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)) || left_back_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)));    //top unit does have parent of matching outside node
          }
        }

      }  else _test( ifp->MatchingN(n,INSIDE) == ifp->MatchingN(n,OUTSIDE) ) ; //test nodes match if on perimeter


    }

    //outside nodes must match in content but may be rotated, so we test just the content
    _test(outside_nds_interface == outside_nds_face);

  }


  //For visualisation purposes
  double gap = 0.5;
  model1.SplitBoundary(sb1).PullApartSplitBoundary(gap);
  model1.SplitBoundary(sb2).PullApartSplitBoundary(0.5);

  //uncommment to visualise
  /*
  list<string> outputProps;
  outputProps.push_back( "nodal id" );
  model1.InputPropertyValue("nodal id", ScalarVariable(PLAIN,1.0));
  VTU_Interface<DIM> vtu1( model1 );
  vtu1.OmitZeroInFileName(true);
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "LEFT_BACK_BOTTOM", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "LEFT_BACK_TOP", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "RIGHT_FRONT_BOTTOM", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "RIGHT_FRONT_TOP", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "Model", static_cast<int>(0) );
  */

  return true;

}















/**
 * Very Rigorous Test:
 *
 * Duplicate of previous test Test_NodeCorrespondance_Intersection_3D but with splitboundary creation in reverse order
 *
*/

template<uint32_t DIM>
bool SplitBoundaryInterface_Test<DIM>::Test_NodeAndElementsCorrespondance_3D_X_Intersection_Reverse( const char* mesh_file){

  const uint32_t dim{3U};
  //model construction
  const char* variables_file("SplitBoundary_Test-variables.txt");
  const char* regions_file("InternalBoundary3D_Intersect_Test");

  bool irregular = false;
  bool reduce_to_regions = true;

  ANSYS_Model3D model1(mesh_file, regions_file, variables_file, irregular, reduce_to_regions);

  //Getting volumetric regions for later use
  Region<dim>& left_back_bottom     = model1.Region("LEFT_BACK_BOTTOM");
  Region<dim>& left_back_top        = model1.Region("LEFT_BACK_TOP");
  Region<dim>& right_front_bottom   = model1.Region("RIGHT_FRONT_BOTTOM");
  Region<dim>& right_front_top      = model1.Region("RIGHT_FRONT_TOP");

  //Getting perimeter nodes of region, before splitboundary is created and information is lost
  std::set<Node<dim>*> perim_nodes_planar, perim_nodes_diagonal;
  Region<dim>& frac_region_planar = model1.Region( "FRACTURE_PLANAR" );
  for ( auto nit = frac_region_planar.PerimeterNodesBegin(); nit !=frac_region_planar.NodesEnd(); nit++){
    perim_nodes_planar.insert(*nit);
  }
  _test( perim_nodes_planar.size() == frac_region_planar.PerimeterNodes() );

  Region<dim>& frac_region_diagonal = model1.Region( "FRACTURE_DIAGONAL" );
  for ( auto nit = frac_region_diagonal.PerimeterNodesBegin(); nit !=frac_region_diagonal.NodesEnd(); nit++){
    perim_nodes_diagonal.insert(*nit);
  }
  _test( perim_nodes_diagonal.size() == frac_region_diagonal.PerimeterNodes() );


  Node<dim>* split_perimter_node = nullptr;
  for ( auto n : perim_nodes_planar ){
    if ( approximatelyEqual(n->Coordinate()[0] , 2.5, 0.01) &&
         approximatelyEqual(n->Coordinate()[1] , 5.0, 0.01) &&
         approximatelyEqual(n->Coordinate()[2] , 7.5, 0.01)  ){
      split_perimter_node = n;                          //this node is on the perimter of fracture planar, but split by fracture diagonal
    }
  }
  assert(split_perimter_node!=nullptr);

  ///Testing splitboundary creation
  //Region -> SplitBoundary
  std::set<string> sb_names_2  = (model1.CreateSplitBoundaryFrom( "FRACTURE_DIAGONAL" ) ).first ;
  std::set<string> sb_names    = (model1.CreateSplitBoundaryFrom( "FRACTURE_PLANAR" ) ).first ;

  //Need to find node on perimeter that was split by diagonal fracture
  Node<dim>* second_split_perimeter_node = nullptr;
  assert(split_perimter_node->Manifold()->Branches()==2);
  split_perimter_node->Manifold()->N(0) == split_perimter_node ?  second_split_perimeter_node = split_perimter_node->Manifold()->N(1) : second_split_perimeter_node->Manifold()->N(0);


  std::string sb1 = model1.MergeSplitBoundaries( "FRACTURE_PLANAR",   sb_names );
  std::string sb2 = model1.MergeSplitBoundaries( "FRACTURE_DIAGONAL", sb_names_2);

  //For Planar SplitBoundary
  for (auto ifp : model1.SplitBoundary(sb1).CellVector() ){

    //Check unit normals are 0 1 0 on INSIDE and MIDDLE and the opposite on OUTSIDE
    Point<dim> nrml_in  = ifp->UnitNormal(INSIDE);
    Point<dim> nrml_out = ifp->UnitNormal(OUTSIDE);
    //Point<dim> nrml_mid = ifp->UnitNormal(MIDDLE);
    //Inside opposite to outside
    _equal( nrml_in[0] ,  nrml_out[0] , std::numeric_limits<double>::epsilon()  );
    _equal( nrml_in[1] , -nrml_out[1] , std::numeric_limits<double>::epsilon()  );
    _equal( nrml_in[2] ,  nrml_out[2] , std::numeric_limits<double>::epsilon()  );
    //inside is 0 1 0
    _equal( nrml_in[1] , 1.0 , 0.0001  );
    _test(  std::fabs(nrml_in[0]) < 0.0001 );
    _test(  std::fabs(nrml_in[2]) < 0.0001 );

    //Retrieve the nodes of face that match with INSIDE OUTSIDE
    std::vector<uint32_t> nids_in  = ifp->InnerParent()->FE()->NodesOfFace( ifp->InnerParentFaceID()),
                          nids_out = ifp->OuterParent()->FE()->NodesOfFace( ifp->OuterParentFaceID());


    std::set<Node<dim>*> outside_nds_interface;
    std::set<Node<dim>*> outside_nds_face;
    uint32_t n_nodes = ifp->FE()->Nodes();
    for ( uint32_t n{0U}; n<n_nodes;++n){
      //Nodes should match that of INSIDE face
      _test( ifp->N(n,INSIDE)  == ifp->InnerParent()->N(nids_in[n] ));

      //Nodes should contain same nodes of OUTSIDE face, but may be rotated, so may not match
      outside_nds_interface.insert( ifp->N(n,OUTSIDE));
      outside_nds_face.insert( ifp->OuterParent()->N(nids_out[n]));

      //Testing matching assignment
      //inside remains unchanged
      _test( ifp->MatchingN(n,INSIDE) == ifp->N(n,INSIDE) );
      //Outside coordinates must match with Inside and middle
      _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,OUTSIDE)->Coordinate() );

      //Testing nodes are duplicated in interior and the same on perimeter
      //if node is not on perimeter
      if ( perim_nodes_planar.find( ifp->N(n,INSIDE) ) == perim_nodes_planar.end() && ifp->N(n,INSIDE) != second_split_perimeter_node  ){
        //test node is duplicated
        _test( ifp->MatchingN(n,INSIDE) != ifp->MatchingN(n,OUTSIDE) );

        //Inside Node is in bottom half
        _test( left_back_bottom.Contains( ifp->N(n,INSIDE) ) || right_front_bottom.Contains( ifp->N(n,INSIDE)) );
        //matching Outside Node is in top half
        _test( left_back_top.Contains( ifp->MatchingN(n,OUTSIDE) ) || right_front_top.Contains(ifp->MatchingN(n,OUTSIDE)) );
        //matching outside node not in bottom half
        _test( !left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE) ) && !right_front_bottom.Contains(ifp->MatchingN(n,OUTSIDE)) );

        //Testing inside parents assignmnet by checking they are part of BottomUnit region
        for (uint32_t p{0U} ; p < ifp->N(n,INSIDE)->Parents() ; ++p){
          if ( ifp->N(n,INSIDE)->Parent(p)->IsEquidimensional()){
            _test( left_back_bottom.Contains( ifp->N(n,INSIDE)->Parent(p)) || right_front_bottom.Contains( ifp->N(n,INSIDE)->Parent(p) )); //bottom unit has parent of inside node
            _test( !left_back_top.Contains( ifp->N(n,INSIDE)->Parent(p)) && !right_front_top.Contains( ifp->N(n,INSIDE)->Parent(p)) );    //top unit doenst have parent of inside node
          }
        }
        for (uint32_t p{0U} ; p < ifp->MatchingN(n,OUTSIDE)->Parents() ; ++p){
          if ( ifp->MatchingN(n,OUTSIDE)->Parent(p)->IsEquidimensional()){
            _test( !left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) ) && !right_front_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) )); //bottom unit Does Not have parent of matching outside node
            _test( left_back_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)) || right_front_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)));    //top unit does have parent of matching outside node
          }
        }

      }  else _test( ifp->MatchingN(n,INSIDE) == ifp->MatchingN(n,OUTSIDE) ) ; //test nodes match if on perimeter


    }

    //outside nodes must match in content but may be rotated, so we test just the content
    _test(outside_nds_interface == outside_nds_face);
  }


  //For Diagonal SplitBoundary
  for (auto ifp : model1.SplitBoundary(sb2).CellVector() ){

    //Check unit normals are 0 1 0 on INSIDE and MIDDLE and the opposite on OUTSIDE
    Point<dim> nrml_in  = ifp->UnitNormal(INSIDE);
    Point<dim> nrml_out = ifp->UnitNormal(OUTSIDE);
    //Point<dim> nrml_mid = ifp->UnitNormal(MIDDLE);
    //Inside opposite to outside
    _equal( nrml_in[0] , -nrml_out[0] , 10.0*std::numeric_limits<double>::epsilon()  );
    _test(  std::fabs(nrml_in[1]) < 0.001 );     //normal should be along x-z plane
    _test(  std::fabs(nrml_out[1]) < 0.001 );    //normal should be along x-z plane
    _equal( nrml_in[2] , -nrml_out[2] , 10.0*std::numeric_limits<double>::epsilon()  );
    //inside is -0.707 0 -0.707
    _equal( nrml_in[0] , -1.0/std::sqrt(2) , 0.001  );
    _equal( nrml_in[2] , -1.0/std::sqrt(2) , 0.001  );
    _test(  std::fabs(nrml_in[1]) < 0.001 );

    //Retrieve the nodes of face that match with INSIDE OUTSIDE
    std::vector<uint32_t> nids_in  = ifp->InnerParent()->FE()->NodesOfFace( ifp->InnerParentFaceID()),
                          nids_out = ifp->OuterParent()->FE()->NodesOfFace( ifp->OuterParentFaceID());


    std::set<Node<dim>*> outside_nds_interface;
    std::set<Node<dim>*> outside_nds_face;
    uint32_t n_nodes = ifp->FE()->Nodes();
    for ( uint32_t n{0U}; n<n_nodes;++n){
      //Nodes should match that of INSIDE face
      _test( ifp->N(n,INSIDE)  == ifp->InnerParent()->N(nids_in[n] ));

      //Nodes should contain same nodes of OUTSIDE face, but may be rotated, so may not match
      outside_nds_interface.insert( ifp->N(n,OUTSIDE));
      outside_nds_face.insert( ifp->OuterParent()->N(nids_out[n]));

      //Testing matching assignment
      //inside remains unchanged
      _test( ifp->MatchingN(n,INSIDE) == ifp->N(n,INSIDE) );
      //Outside coordinates must match with Inside and middle
      _test( ifp->MatchingN(n,INSIDE)->Coordinate() == ifp->MatchingN(n,OUTSIDE)->Coordinate() );

      //Testing nodes are duplicated in interior and the same on perimeter
      //if node is not on perimeter
      if ( perim_nodes_diagonal.find( ifp->N(n,INSIDE) ) == perim_nodes_diagonal.end()  ){
        //test node is duplicated
        _test( ifp->MatchingN(n,INSIDE) != ifp->MatchingN(n,OUTSIDE) );

        //Inside Node is in bottom half
        _test( right_front_bottom.Contains( ifp->N(n,INSIDE) ) || right_front_top.Contains( ifp->N(n,INSIDE)) );
        //matching Outside Node is in top half
        _test( left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE) ) || left_back_top.Contains(ifp->MatchingN(n,OUTSIDE)) );
        //matching outside node not in inside diagonal
        _test( !right_front_bottom.Contains( ifp->MatchingN(n,OUTSIDE) ) && !right_front_top.Contains(ifp->MatchingN(n,OUTSIDE)) );
        //inside node not in outside diagonal
        _test( !left_back_bottom.Contains( ifp->MatchingN(n,INSIDE) ) && !left_back_top.Contains(ifp->MatchingN(n,INSIDE)) );

        //Testing inside parents assignmnet by checking they are part of BottomUnit region
        for (uint32_t p{0U} ; p < ifp->N(n,INSIDE)->Parents() ; ++p){
          if ( ifp->N(n,INSIDE)->Parent(p)->IsEquidimensional()){
            _test( right_front_bottom.Contains( ifp->N(n,INSIDE)->Parent(p)) || right_front_top.Contains( ifp->N(n,INSIDE)->Parent(p) )); //bottom unit has parent of inside node
            _test( !left_back_bottom.Contains( ifp->N(n,INSIDE)->Parent(p)) && !left_back_top.Contains( ifp->N(n,INSIDE)->Parent(p)) );    //top unit doenst have parent of inside node
          }
        }
        for (uint32_t p{0U} ; p < ifp->MatchingN(n,OUTSIDE)->Parents() ; ++p){
          if ( ifp->MatchingN(n,OUTSIDE)->Parent(p)->IsEquidimensional()){
            _test( !right_front_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) ) && !right_front_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p) )); //bottom unit Does Not have parent of matching outside node
            _test( left_back_bottom.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)) || left_back_top.Contains( ifp->MatchingN(n,OUTSIDE)->Parent(p)));    //top unit does have parent of matching outside node
          }
        }

      }  else _test( ifp->MatchingN(n,INSIDE) == ifp->MatchingN(n,OUTSIDE) ) ; //test nodes match if on perimeter


    }

    //outside nodes must match in content but may be rotated, so we test just the content
    _test(outside_nds_interface == outside_nds_face);

  }


  //For visualisation purposes
  double gap = 0.5;
  model1.SplitBoundary(sb1).PullApartSplitBoundary(gap);
  model1.SplitBoundary(sb2).PullApartSplitBoundary(0.5);

  // Uncomment to visualise
  /*
  list<string> outputProps;
  outputProps.push_back( "nodal id" );
  model1.InputPropertyValue("nodal id", ScalarVariable(PLAIN,1.0));
  VTU_Interface<dim> vtu1( model1 );
  vtu1.OmitZeroInFileName(true);
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "LEFT_BACK_BOTTOM", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "LEFT_BACK_TOP", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "RIGHT_FRONT_BOTTOM", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "RIGHT_FRONT_TOP", static_cast<int>(0) );
  vtu1.OutputDataToVTU( "../Output/InternalBoundary3D_Intersect_TestA", outputProps, "Model", static_cast<int>(0) );
  */

  return true;

}

























































template<uint32_t dim>
void shiftSplitBoundary( SplitBoundary<dim>& splitboundary, INTERFACE_SIDE side, double xShift, double yShift, double zShift )
{
  const typename vector<InterFace<dim>*>::const_iterator facesEnd( splitboundary.CellsEnd() );
  if constexpr ( dim == 2 )
    for ( auto it = splitboundary.CellsBegin(); it != facesEnd; ++it )
    {
      const size_t nodes( (*it)->Nodes() );
      for ( auto i = 0; i < nodes; ++i )
      {
        if ( (*it)->N( i, INSIDE )->Idx() != (*it)->N( i, OUTSIDE )->Idx() )
        {
          (*it)->N( i, side )->x( (*it)->N( i, side )->x() + xShift );
          (*it)->N( i, side )->y( (*it)->N( i, side )->y() + yShift );
        }
      }
    }
  else if constexpr ( dim == 3 )
    for ( auto it = splitboundary.CellsBegin(); it != facesEnd; ++it )
    {
      const size_t nodes( (*it)->FE()->Nodes() );
      for ( auto i = 0; i < nodes; ++i )
      {
        if ( (*it)->N( i, INSIDE )->Idx() != (*it)->N( i, OUTSIDE )->Idx() )
        {
          (*it)->N( i, side )->x( (*it)->N( i, side )->x() + xShift );
          (*it)->N( i, side )->y( (*it)->N( i, side )->y() + yShift );
          (*it)->N( i, side )->z( (*it)->N( i, side )->z() + zShift );
        }
      }
    }
}


template<uint32_t dim>
void shiftSplitBoundary( SplitBoundary<dim>& splitboundary, double shift )
{
  Point<dim> displacementPerpedicularToInterface;
  Point<dim> displacementToBaryCenter;
  Point<dim> baryCenter;
  size_t null_neighbors( 0 );
  Element<dim>* parentElement( nullptr );

  const typename vector<InterFace<dim>*>::const_iterator facesEnd( splitboundary.CellsEnd() );
  if constexpr ( dim == 2U )
    for ( auto it = splitboundary.CellsBegin(); it != facesEnd; ++it )
    {
      if ( (*it)->InterveningElement() != nullptr )
        displacementPerpedicularToInterface = (*it)->InterveningElement()->UnitNormal();
      else {
           (*it)->CurrentSide( INSIDE );
           displacementPerpedicularToInterface = (*it)->UnitNormal();
        }
      displacementPerpedicularToInterface *= shift;
      const size_t nodes( (*it)->FE()->Nodes() );
      for ( auto i{0U}; i < nodes; ++i )
      {
        if ( (*it)->N( i, INSIDE ) != (*it)->N( i, OUTSIDE ) )
        {
          null_neighbors = 0;
          for ( uint32_t k{0U}; k<(*it)->N( i, OUTSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, OUTSIDE )->Parent( k );
            for ( uint32_t j{0U}; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == nullptr )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( OUTSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, OUTSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementToBaryCenter[1] );
          }
          (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementPerpedicularToInterface[0] );
          (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementPerpedicularToInterface[1] );

          null_neighbors = 0;
          for ( uint32_t k{0U}; k<(*it)->N( i, INSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, INSIDE )->Parent( k );
            for ( uint32_t j{0U}; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == nullptr )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( INSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, INSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() + displacementToBaryCenter[1] );
          }
          (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() - displacementPerpedicularToInterface[0] );
          (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() - displacementPerpedicularToInterface[1] );
        }
      }
    }
  else if constexpr ( dim == 3 )
    for ( auto it = splitboundary.CellsBegin(); it != facesEnd; ++it )
    {
      const auto nodes( (*it)->FE()->Nodes() );
      for ( auto i{0U}; i < nodes; ++i )
      {
        if ( (*it)->InterveningElement() != nullptr )
          displacementPerpedicularToInterface = (*it)->InterveningElement()->UnitNormal();
        else {
            (*it)->CurrentSide( INSIDE );
            displacementPerpedicularToInterface = (*it)->UnitNormal();
          }
        displacementPerpedicularToInterface *= shift;
        if ( (*it)->N( i, INSIDE ) != (*it)->N( i, OUTSIDE ) )
        {
          null_neighbors = 0;
          for ( uint32_t k{0U}; k<(*it)->N( i, OUTSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, OUTSIDE )->Parent( k );
            for ( uint32_t j{0U}; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == nullptr )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( OUTSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, OUTSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementToBaryCenter[1] );
            (*it)->N( i, OUTSIDE )->z( (*it)->N( i, OUTSIDE )->z() + displacementToBaryCenter[2] );
          }
          (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementPerpedicularToInterface[0] );
          (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementPerpedicularToInterface[1] );
          (*it)->N( i, OUTSIDE )->z( (*it)->N( i, OUTSIDE )->z() + displacementPerpedicularToInterface[2] );

          null_neighbors = 0;
          for ( auto k{0U}; k<(*it)->N( i, INSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, INSIDE )->Parent( k );
            for ( auto j{0U}; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == nullptr )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( INSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, INSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() + displacementToBaryCenter[1] );
            (*it)->N( i, INSIDE )->z( (*it)->N( i, INSIDE )->z() + displacementToBaryCenter[2] );
          }
          (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() - displacementPerpedicularToInterface[0] );
          (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() - displacementPerpedicularToInterface[1] );
          (*it)->N( i, INSIDE )->z( (*it)->N( i, INSIDE )->z() - displacementPerpedicularToInterface[2] );
        }
      }
    }
}


template<uint32_t dim>
void shiftInterfaceTips( Region<dim>& region, double shift )
{
  Point<dim> displacement;
  Point<dim> baryCenter;

  if constexpr ( dim == 2U )
    for ( auto it = region.PerimeterCellsBegin(); it != region.CellsEnd(); ++it )
    {
      baryCenter = (*it)->BaryCenter();
      const auto nodes( (*it)->FE()->Nodes() );
      for ( auto i{0U}; i < nodes; ++i )
      {
        displacement = baryCenter - (*it)->N( i )->Coordinate();
        displacement.NormalizeLengthTo( 1.0 );
        displacement *= shift;
        (*it)->N( i )->x( (*it)->N( i )->x() + displacement[0] );
        (*it)->N( i )->y( (*it)->N( i )->y() + displacement[1] );
      }
    }
  else if constexpr ( dim == 3U )
    for ( auto it = region.PerimeterCellsBegin(); it != region.CellsEnd(); ++it )
    {
      baryCenter = (*it)->BaryCenter();
      const size_t nodes( (*it)->Nodes() );
      for ( auto i{0U}; i < nodes; ++i )
      {
        displacement = baryCenter - (*it)->N( i )->Coordinate();
        displacement.NormalizeLengthTo( 1.0 );
        displacement *= shift;
        (*it)->N( i )->x( (*it)->N( i )->x() + displacement[0] );
        (*it)->N( i )->y( (*it)->N( i )->y() + displacement[1] );
        (*it)->N( i )->z( (*it)->N( i )->z() + displacement[2] );
      }
    }
}


template<uint32_t dim>
void shiftRegion( Region<dim>& region, double xShift, double yShift, double zShift )
{
  if constexpr ( dim == 2U )
    for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x() + xShift );
      (*it)->y( (*it)->y() + yShift );
    }
  else if constexpr ( dim == 3U )
    for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x() + xShift );
      (*it)->y( (*it)->y() + yShift );
      (*it)->z( (*it)->z() + zShift );
    }
}


template<uint32_t dim>
void shiftRegionAboveLine( Region<dim>& region, size_t x_or_y_or_z, double line_coordinate, double shift, double eps )
{
  if constexpr ( dim == 2U )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() > (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() > (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
  }
  else if constexpr ( dim == 3U )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() > (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() > (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
    else {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->z() > (line_coordinate + eps) )
          (*it)->z( (*it)->z() + shift );
    }
  }
}


template<uint32_t dim>
void shiftRegionBelowLine( Region<dim>& region, size_t x_or_y_or_z, double line_coordinate, double shift, double eps )
{
  if constexpr ( dim == 2U )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() < (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() < (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
  }
  else if constexpr ( dim == 3U )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() < (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() < (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
    else {
      for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->z() < (line_coordinate + eps) )
          (*it)->z( (*it)->z() + shift );
    }
  }
}


template<uint32_t dim>
void scaleRegionSymmetricOverZero( Region<dim>& region, double xScale, double yScale, double zScale )
{
  if constexpr ( dim == 2U )
    for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x()*xScale );
      (*it)->y( (*it)->y()*yScale );
    }
  else if constexpr ( dim == 3U )
    for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x()*xScale );
      (*it)->y( (*it)->y()*yScale );
      (*it)->z( (*it)->z()*zScale );
    }
}


template<uint32_t dim>
void scaleRegion( Region<dim>& region, double xScale, double yScale, double zScale )
{
  Point<dim> min_point;
  Point<dim> max_point;
  Point<dim> mid_point;

  region.MinMaxCoordinates( min_point, max_point );
  mid_point = (min_point + max_point) / 2.0;

  if constexpr ( dim == 2U )
    for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( mid_point[0] + ((*it)->x() - mid_point[0])*xScale );
      (*it)->y( mid_point[1] + ((*it)->y() - mid_point[1])*yScale );
    }
  else if constexpr ( dim == 3U )
    for ( auto it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( mid_point[0] + ((*it)->x() - mid_point[0])*xScale );
      (*it)->y( mid_point[1] + ((*it)->y() - mid_point[1])*yScale );
      (*it)->z( mid_point[2] + ((*it)->z() - mid_point[2])*zScale );
    }
}


template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::NodeParents( const Region<dim>& region, size_t minParentCount /* = 2 */ )
{
  const auto nodesEnd( region.NodesEnd() );
  for ( auto node( region.NodesBegin() ); node != nodesEnd; ++node )
    _test( (*node)->Parents() >= minParentCount );
}


template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::ElementNodes( const Region<dim>& region )
{
  const auto elementsEnd( region.CellsEnd() );
  for ( auto element( region.CellsBegin() ); element != elementsEnd; ++element )
    for ( auto i{0U}; i<(*element)->Nodes(); i++ )
    _test( (*element)->N(i) != nullptr );
}


/// could be considered a BoundaryInterface/Element/Node test
template<uint32_t dim>
void csmp::SplitBoundaryInterface_Test<dim>::CheckRemovedLowDimParents( Boundary<dim>& boundary )
{
  const auto nodesEnd( boundary.NodesEnd() );
  for ( auto node( boundary.NodesBegin() ); node != nodesEnd; ++node )
    for ( auto parent{0U}; parent < (*node)->Parents(); ++parent )
      _test( !( (*node)->Parent( parent )->IsEquidimensional() ) );
}


// OK
template<uint32_t dim>
bool SplitBoundaryInterface_Test<dim>::NoNeighborNull( const csmp::InterFace<dim>& interFace )
{
  for ( auto n{0U}; n < interFace.Neighbors(); ++n )
    if ( interFace.Neighbor( n ) == nullptr )
      return false;
  return true;
}


/// Test whether the splitted nodes share parent elements, located on different sides of interfaces
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::TestSplitNodeAssignment( const Model<dim>& model )
{
  set<size_t> outerParentsNodes;
  set<size_t> innerParentsNodes;
  set<size_t> not_duplicated_nodes;
  not_duplicated_nodes.clear();

  // loop over SplitBoundaries
  for ( auto spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
  {
    // loop over InterFaces
    for ( auto ifit = spbit->second.CellsBegin(); ifit != spbit->second.CellsEnd(); ++ifit )
    {
      if ( NoNeighborNull( *(*ifit) ) )
      {
        // interior interfaces

        const Element<dim>* const outerParent( (*ifit)->Parent( OUTSIDE ) );
        for ( auto n{0U}; n < outerParent->Nodes(); ++n )
          outerParentsNodes.insert( outerParent->N( n )->Idx() );

        const Element<dim>* const innerParent( (*ifit)->Parent( INSIDE ) );
        for ( auto n{0U}; n < innerParent->Nodes(); ++n )
          innerParentsNodes.insert( innerParent->N( n )->Idx() );

        // if these test fails, the boundary was not split correctly
        for ( auto it( outerParentsNodes.begin() ); it != outerParentsNodes.end(); ++it )
          _test( innerParentsNodes.find( *it ) == innerParentsNodes.end() );

        for ( auto it( innerParentsNodes.begin() ); it != innerParentsNodes.end(); ++it )
          _test( outerParentsNodes.find( *it ) == outerParentsNodes.end() );

        outerParentsNodes.clear();
        innerParentsNodes.clear();

        for ( auto n{0U}; n < (*ifit)->Nodes(); ++n )
        {

          _test( (*ifit)->N( n, INSIDE )->Idx() != (*ifit)->N( n, OUTSIDE )->Idx() );

          set<size_t> parents;
          for ( auto p{0U}; p < (*ifit)->N( n, OUTSIDE )->Parents(); ++p )
            parents.insert( (*ifit)->N( n, OUTSIDE )->Parent( p )->Idx() );

          for ( auto p{0U}; p < (*ifit)->N( n, INSIDE )->Parents(); ++p )
            _test( parents.find( (*ifit)->N( n, INSIDE )->Parent( p )->Idx() ) == parents.end() );
          parents.clear();

          for ( auto p{0U}; p < (*ifit)->N( n, INSIDE )->Parents(); ++p )
            parents.insert( (*ifit)->N( n, INSIDE )->Parent( p )->Idx() );
          for ( auto p{0U}; p < (*ifit)->N( n, OUTSIDE )->Parents(); ++p )
            _test( parents.find( (*ifit)->N( n, OUTSIDE )->Parent( p )->Idx() ) == parents.end() );
          parents.clear();

        }

      }
      else {
        // perimeter interfaces
        for ( auto i{0U}; i<(*ifit)->Nodes(); ++i )
          if ( (*ifit)->N( i, INSIDE )->Idx() == (*ifit)->N( i, OUTSIDE )->Idx() )
            not_duplicated_nodes.insert( (*ifit)->N( i )->Idx() );
      }
    }
  }
  
} // TestSplitNodeAssignment





template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::TestUnitNormals( Model<dim>& model, const string& test_name )
{
  // -----------------------------
  // SplitBoundary Normals
  // -----------------------------

  vector<double> zero_vec( dim, 0.0 );

  VectorVariable<dim> nrml_in( zero_vec );
  VectorVariable<dim> nrml_out( zero_vec );
  VectorVariable<dim> nrml_sum( zero_vec );
  VectorVariable<dim> nrml_zero( zero_vec );

  model.InputPropertyValue( "element vector", nrml_zero );

  Index normal_idx( model.Database().StorageKey( "element vector" ) );

  // loop over SplitBoundaries
  for ( auto spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
  {
    // loop over InterFaces
    for ( auto ifit = spbit->second.CellsBegin(); ifit != spbit->second.CellsEnd(); ++ifit )
    {
      nrml_in  = (*ifit)->UnitNormal( INSIDE );
      nrml_out = (*ifit)->UnitNormal( OUTSIDE );

      nrml_sum = nrml_in + nrml_out;

      if(dim < 3) _test( nrml_sum == nrml_zero );
      (*ifit)->Parent( INSIDE )->Store( normal_idx, nrml_in );
      (*ifit)->Parent( OUTSIDE )->Store( normal_idx, nrml_out );
    }
  }

} // end TestUnitNormals






template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::VisualiseSplitBoundaries( Model<dim>& model, const string& test_name )
{
  // visualization
  VTU_Interface<dim> vtu( model );
  vtu.OmitZeroInFileName( false );
  list<string> outputProps;
  outputProps.push_back( "nodal id" );
  outputProps.push_back( "nodal variable" );
  outputProps.push_back( "element variable" );

  // -----------------------------
  // SplitBoundary Normals
  // -----------------------------

  const ScalarVariable matrixValue( ANY, 0.0 );
  ScalarVariable       interfaceValue( ANY, 0.0 );
  ScalarVariable       interfaceValueWrite( ANY, 0.0 );
  ScalarVariable       interfaceValueRead( ANY, 0.0 );

  model.InputPropertyValue( "nodal id", matrixValue );
  model.InputPropertyValue( "nodal variable", matrixValue );
  model.InputPropertyValue( "element variable", matrixValue );

  Index node_idx( model.Database().StorageKey( "nodal id" ) );
  Index node_prop_idx( model.Database().StorageKey( "nodal variable" ) );
  Index element_prop_idx( model.Database().StorageKey( "element variable" ) );

  // loop over SplitBoundaries
  for ( auto spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
  {
    interfaceValue += 1.0;

    // loop over InterFaces
    for ( auto ifit = spbit->second.CellsBegin(); ifit != spbit->second.CellsEnd(); ++ifit )
    {
      (*ifit)->Parent( INSIDE )->Read( element_prop_idx, interfaceValueRead );
      interfaceValueWrite = +1;
      (*ifit)->Parent( INSIDE )->Store( element_prop_idx, interfaceValueWrite );
      (*ifit)->Parent( OUTSIDE )->Read( element_prop_idx, interfaceValueRead );
      interfaceValueWrite = -1;
      (*ifit)->Parent( OUTSIDE )->Store( element_prop_idx, interfaceValueWrite );

      for ( auto i = 0; i<(*ifit)->FE()->Nodes(); ++i )
      {
        if ( (*ifit)->N( i, OUTSIDE )->Idx() != (*ifit)->N( i, INSIDE )->Idx() )
        {
          interfaceValueWrite = -1.0;
          (*ifit)->N( i, INSIDE )->Store( node_prop_idx, interfaceValueWrite );
          interfaceValueWrite = (*ifit)->N( i, INSIDE )->Idx();
          (*ifit)->N( i, INSIDE )->Store( node_idx, interfaceValueWrite );

          interfaceValueWrite = 1.0;
          (*ifit)->N( i, OUTSIDE )->Store( node_prop_idx, interfaceValueWrite );
          interfaceValueWrite = (*ifit)->N( i, OUTSIDE )->Idx();
          (*ifit)->N( i, OUTSIDE )->Store( node_idx, interfaceValueWrite );
        }
      }
    }
  }

  if ( verbose_ ) vtu.OutputDataToVTU( test_name.c_str(), outputProps, "Model", 0.0 );

  for ( auto it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); it++ )
    vtu.OutputDataToVTU( (*it).first.c_str(), outputProps, (*it).first.c_str(), 0.0 );
  return;
}


template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::PullApartSplitboundaries( Model<dim>& model, vector<string>& interfaces, double displacement )
{
  // pull apart SplitBoundaries
  string split_boundary_name;

  size_t split_boundaries( model.SplitBoundaries() );
  for ( auto i = 0; i<split_boundaries; i++ )
  {
    split_boundary_name = "SPLITBOUNDARY_" + interfaces[i];
    shiftSplitBoundary( model.SplitBoundary( split_boundary_name.c_str() ), displacement );
    shiftInterfaceTips( model.Region( interfaces[i].c_str() ), displacement );
  }
}


template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::PullApartSplitboundaries( Model<dim>& model, double displacement )
{
  //pPull apart SplitBoundaries
  for ( auto spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
    shiftSplitBoundary( (*spbit).second, displacement );
}


template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::LoadModel( const string& model_name )
{
  // Model initialization
  const string variables_file("SplitBoundary_Test-variables.txt");
  Model<dim>* model = NULL;
  if constexpr ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  else if constexpr ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));

  // visualization
  VTU_Interface<dim> vtu( *model );
  vtu.OmitZeroInFileName( false );
  list<string> outputProps;
  outputProps.push_back( "element variable" );
  ScalarVariable regionValue( ANY, 0.0 );

  if ( verbose_ ) cout << "\n\n\nSplitBoundary_Test::PrepareModel: the following interface / interface(s) sets will be considered:\n\n";
  model->InputPropertyValue( "element variable", regionValue );
  for ( auto it = model->RegionsBegin(); it != model->RegionsEnd(); it++ )
  {
    regionValue += 1.0;
    (*it).second.InputPropertyValue( "element variable", regionValue );
    vtu.OutputDataToVTU( model_name.c_str(), outputProps, (*it).first.c_str(), 0.0 );
  }
  if ( verbose_ ) vtu.OutputDataToVTU( model_name.c_str(), outputProps, "Model", 0.0 );

  model->OutputToBinaryFile( model_name.c_str() );
}


template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::LoadModel( const string& model_name,
                                             vector<string>& regions )
{
  // Model initialization
  const string variables_file("SplitBoundary_Test-variables.txt");
  Model<dim>* model = NULL;
  if constexpr ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  else if constexpr ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));

  // visualization
  VTU_Interface<dim> vtu( *model );
  vtu.OmitZeroInFileName( false );
  list<string> outputProps;
  outputProps.push_back( "element variable" );
  ScalarVariable regionValue( ANY, 0.0 );

  string  spliboundary_regions_file( model_name );

  // load regions from the input file
  InputFromFile( string( spliboundary_regions_file + "-disconnected-interface-regions.txt" ).c_str(), regions );

  if ( verbose_ ) cout << "\n\n\nSplitBoundary_Test::PrepareModel: the following interface / interface(s) sets will be considered:\n\n";
  model->InputPropertyValue( "element variable", regionValue );
  for ( auto it = regions.begin(); it != regions.end(); it++ )
    {
      regionValue += 1.0;
      model->Region( (*it).c_str() ).InputPropertyValue( "element variable", regionValue );
      vtu.OutputDataToVTU( model_name.c_str(), outputProps, (*it).c_str(), 0.0 );
    }
  if ( verbose_ ) vtu.OutputDataToVTU( model_name.c_str(), outputProps, "Model", 0.0 );

  model->OutputToBinaryFile( model_name.c_str() );
}





template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::LoadContiguousModel( const string& model_name,
                                                       vector<string>& interfaces )
{
  // 0. Model initialization
  const string variables_file("SplitBoundary_Test-variables.txt");
  Model<dim>* model(nullptr);
  if constexpr ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  else if constexpr ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));

  string  spliboundary_regions_file( model_name );

  // 1. InterFace sets
  set<string> interface_basic_sets;
  set<string> interface_sets;
  InputFromFile( string( spliboundary_regions_file + "-disconnected-interface-regions.txt" ).c_str(), interface_basic_sets );

  // 2. splitting input regions if they are discontigouos
  bool discontiguous_regions( false );
  for ( auto it = interface_basic_sets.begin(); it != interface_basic_sets.end(); ++it )
    if ( !model->Region( (*it) ).IsContiguous() ) {
      if ( verbose_ ) cerr << "\n\tSplitBoundary_Test::PrepareModel: discovered discontiguous region: " << (*it);
      discontiguous_regions = true;
    }
  if ( discontiguous_regions ) {
    set<string>  original_region_names;
    for ( auto it = model->UniqueRegionsBegin(); it != model->UniqueRegionsEnd(); ++it )
      if ( (*it).first != "Model" )
        original_region_names.insert( (*it).first.c_str() );

    // partitioning regions without revisiting new partitions that can inserted into region map
    for ( set<string>::const_iterator it = original_region_names.begin(); it != original_region_names.end(); ++it )
      model->PartitionRegionIntoContiguousSubRegions( (*it).c_str() );
  }

  // 2. preparing lower dimensional regions for making SplitBoundaries around
  EstablishContiguousRegionsList( *model, interface_basic_sets, interface_sets );

  model->MergeRegions( interface_sets, "interfaces" );

  for ( auto& name : interface_sets )
    model->RemoveRegion( name.c_str(), false );

  interfaces.clear();
  interfaces.push_back( "interfaces" );
  OutputToFile( string( spliboundary_regions_file + "-connected-interface-regions.txt" ).c_str(), interfaces );

  if ( verbose_ ) cout << "\n\n\nSplitBoundary_Test::PrepareModel: the following interface(s) / interface sets will be considered:\n\n";

  // 3. creating SplitBoundaries
  for ( auto it = interfaces.begin(); it != interfaces.end(); it++ )
    model->CreateSplitBoundaryFrom( (*it).c_str() );

  // 4. creating lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  //    insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  const int32_t material_id_for_new_elements{2};
  set<string> new_regions = model->InsertLowerDimensionalRegionsIntoSplitBoundaries(material_id_for_new_elements);

  // 5. Writing it 
  model->OutputToBinaryFile( model_name.c_str() );
}




/// Input name of regions to split
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::InputFromFile( const char* file_name,
                                                      set<string>& interface_basic_set )
{
  assert( file_name != NULL );

  ifstream ifs( file_name );
  if ( !ifs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile",
                           file_name,
                           "file specifying regions is missing" );
  if ( !interface_basic_set.empty() )
    interface_basic_set.clear();

  // reading header line printing it to screen and swallowing empty line thereafter
  char text[256];
  ifs.getline( text, 256 );
  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: file header: " << text << endl;
  ifs.getline( text, 256 );

  int n_interfaces( 0 );
  ifs >> n_interfaces;
  assert( n_interfaces > 0 );
  assert( n_interfaces < 10000 );

  string  interface_name;
  for ( int n = 0; n<n_interfaces; n++ ) {
    ifs >> interface_name;
    if ( !interface_name.empty() )
      interface_basic_set.insert( interface_name );
    else
      throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile:", "encountered empty region name." );
  }

  ifs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: region names stored in '" << file_name << "' read successfully." << endl;

} // end inputFromFile


/// Input name of regions to split
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::EstablishContiguousRegionsList( Model<dim>& model,
                                                                  const set<string>& interface_basic_set,
                                                                  set<string>& interface_sets )
{
  if ( verbose_ ) cout << "\nSplitBoundary_Test::EstablishContiguosRegionsList:" << endl;

  if ( !interface_sets.empty() )
    interface_sets.clear();
  bool first_call( true );
  for ( auto git = model.UniqueRegionsBegin(); git != model.UniqueRegionsEnd(); git++ )
  {
    for ( auto it = interface_basic_set.begin(); it != interface_basic_set.end(); ++it )
      if ( (*git).first.find( *it ) != string::npos )
      {
        interface_sets.insert( (*git).first );
        if ( first_call ) {
          cout << "\n\tregion(s) incorporated into the input list: ";
          first_call = false;
        }
        cout << (*git).first << " ";
        break;
      }
  }
  if ( verbose_ ) cout << "\n\n";

} // end inputFromFile


/// Input name of regions to split
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::InputFromFile( const char* file_name,
                                                 vector<string>& interfaces )
{
  assert( file_name != NULL );

  ifstream ifs( file_name );
  if ( !ifs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile",
                           file_name,
                           "file specifying regions is missing" );
  if ( !interfaces.empty() )
    interfaces.clear();

  // reading header line printing it to screen and swallowing empty line thereafter
  char text[256];
  ifs.getline( text, 256 );
  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: file header: " << text << endl;
  ifs.getline( text, 256 );

  int n_interfaces( 0 );
  ifs >> n_interfaces;
  assert( n_interfaces > 0 );
  assert( n_interfaces < 10000 );

  string  interface_name;
  for ( int n = 0; n<n_interfaces; n++ ) {
    ifs >> interface_name;
    if ( !interface_name.empty() )
      interfaces.push_back( interface_name );
    else
      throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile:", "encountered empty region name." );
  }

  ifs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: region names stored in '" << file_name << "' read successfully." << endl;

} // end inputFromFile





/// Write contiguous regions
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::OutputToFile( const char* file_name,
                                                const vector<string>& interfaces )
{
  assert( file_name != NULL );

  ofstream ofs( file_name );
  if ( !ofs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::outputToFile",
                           file_name, "file specifying permeability-model input variables is missing" );

  ofs << "'" << file_name << "' interface regions to be included.\n\n";

  ofs << interfaces.size() << " ";
  for ( vector<string>::const_iterator it = interfaces.begin(); it != interfaces.end(); ++it )
    ofs << (*it) << " ";

  ofs << "\n";
  ofs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::outputToFile: '" << file_name << "' written successfully." << endl;

} // end outputToFile




/// Write interface regions
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::OutputToFile( const char* file_name,
                                                const set<string>& interfaces )
{
  assert( file_name != NULL );

  ofstream ofs( file_name );
  if ( !ofs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::outputToFile",
                           file_name, "file specifying permeability-model input variables is missing" );

  ofs << "'" << file_name << "' interface regions to be included.\n\n";

  ofs << interfaces.size() << " ";
  for ( auto it = interfaces.begin(); it != interfaces.end(); ++it )
    ofs << (*it) << " ";

  ofs << "\n";
  ofs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::outputToFile: '" << file_name << "' written successfully." << endl;

} // end outputToFile


/// TESTS
/// SPLITBOUNDARY BETWEEN REGIONS
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::Test_splitboundary_between_regions( const string& model_name )
{
  ostringstream ostr;
  string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) cerr << "\nStart " << dimension << " SplitBoundary Test: SplitBoundary between Regions\n";

  string test_name( "SPLITBOUNDARY_TEST_BETWEEN_REGIONS_" );
  test_name += dimension;
  test_name += "_";
  test_name += model_name;

  // load Model
  const string variables_file( "SplitBoundary_Test-variables.txt" );

  Model<dim>* modelIN(nullptr);

  if constexpr ( dim == 2U )
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  else if constexpr ( dim == 3U )
    // SKM FIX: irregular = true
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), true, true ));

  // validating the model
  if ( verbose_ ) cout << "\nSplitBoundary_Test<" << dim << ">::test_splitboundary_between_regions: model contains the regions:";
  for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ )
    if ( verbose_ ) cout << "\n\tboundary: " << (*it).first;
  if ( verbose_ ) cout << endl;
  if ( verbose_ ) cout << "\nSplitBoundary_Test<" << dim << ">::test_splitboundary_between_regions: model contains the boundaries:";
  for ( auto it = modelIN->BoundariesBegin(); it != modelIN->BoundariesEnd(); it++ )
    if ( verbose_ ) cout << "\n\tboundary: " << (*it).first;
  if ( verbose_ ) cout << endl;

  // create SplitBoundaries
  vector<string> regions;
  regions.reserve( modelIN->UniqueRegions() );

  const pair<int32_t, int32_t>  model_dim = modelIN->Region( "Model" ).SpatialDimensions();
  for ( auto it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); ++it ) {
    const pair<int32_t, int32_t>  sub_dim = (*it).second.SpatialDimensions();
    if ( sub_dim.second == model_dim.second ) // check whether the highest dimension of the region is equal to the highest dimension of the model
      regions.push_back( (*it).second.Name() );
  }

  // search the existing unique sub-regions
  set<pair<string, string>>	discovered;
  deque<string>	current_regions;
  vector<pair<string, string>>  region_final_pairs;
  for ( auto i = 0U; i < regions.size(); i++ ) {
    string root = regions[i];
    // starting at the first region
    current_regions.push_back( root );
    while ( !current_regions.empty() ) {
      string current_region( *current_regions.begin() );
      set<string> neighbors;
      // for all neighbor sub-regions of the current region
      const csmp::Region<dim>&  gref1( modelIN->Region( current_region ) );
      for ( size_t j = 0; j < regions.size(); j++ ) {
        if ( current_region.compare( regions[j] ) == 0 ) continue;
        const csmp::Region<dim>&  gref2( modelIN->Region( regions[j] ) );
        const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
        if ( shared_nodes > 1 )
          neighbors.insert( regions[j] );
      }

      for ( auto neighbour_region : neighbors ) {
        // if this neighbor is new one        
        auto new_region = discovered.insert( make_pair( current_region, neighbour_region ) );
        if ( new_region.second ) {
          discovered.insert( make_pair( neighbour_region, current_region ) );
          current_regions.push_back( neighbour_region );
          region_final_pairs.push_back( make_pair( current_region, neighbour_region ) );
        }
      }
      // removing the sub-region from the discovered (but not yet explored) deque
      current_regions.pop_front();
    }
  }

  // CREATION OF SPLITBOUNDARY BETWEEN 2 REGIONS
  // -------------------------------------------
  // do this sequentially according to neighbours, otherwise boundaries are not assgiend properly
  for ( auto it : region_final_pairs ) // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
    modelIN->CreateSplitBoundaryBetween( it.first.c_str(), it.second.c_str() );
  
  // create lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  // insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  const int32_t rocktype(8);
  const int32_t material_id_for_new_elements(rocktype);
  set<string> new_regions = modelIN->InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id_for_new_elements );

  // read Model from Binary
  modelIN->OutputToBinaryFile( model_name.c_str() );
  /* //added
  set<string> subset_variables;
  Model<dim> model( model_name ); //added
  cout << "\nNodes: " << model.Mesh().Nodes() << "\n";
  cout << "\nNode Groups: " << model.Mesh().NodeGroups() << "\n";
  cout << "\nElements: " << model.Mesh().Elements() << "\n";
  cout << "\nElement Groups: " << model.Mesh().ElementGroups() << "\n";
  cout << "\nFaces: " << model.Mesh().Faces() << "\n";
  cout << "\nFace Groups: " << model.Mesh().FaceGroups() << "\n";
  cout << "\nInterfaces: " << model.Mesh().InterFaces() << "\n";
  cout << "\nInterface Groups: " << model.Mesh().InterFaceGroups() << "\n";

  // tests
  TestSplitNodeAssignment( model ); 
  TestUnitNormals( model, test_name.c_str() );
  double displacement( 0.001 );
  PullApartSplitboundaries( model, displacement );
  VisualiseSplitBoundaries( model, test_name );
  */ //added

  //added
  // tests
  double displacement( 0.001 );
  PullApartSplitboundaries( *modelIN, displacement );
  VisualiseSplitBoundaries( *modelIN, test_name );  


  
  // screen output
  modelIN->RegionsOut();
  modelIN->BoundariesOut();
  modelIN->SplitBoundariesOut();

  return;
}


/// SPLITBOUNDARY AROUND REGIONS
template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::Test_splitboundary_around_regions( const string& model_name )
{
  ostringstream ostr;
  string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) cerr << "\nStart " << dimension << " SplitBoundary Test: SplitBoundary around Regions\n";
    
  // load Model
  vector<string> interfaces;
  LoadContiguousModel( model_name, interfaces );

  // read from binary
  string spliboundary_regions_file( model_name );
  spliboundary_regions_file += "-connected-interface-regions.txt";
  InputFromFile( spliboundary_regions_file.c_str(), interfaces );

// TODO: where is the test?

  Model<dim>  model_out( model_name );
  if ( verbose_ ) {
      cout << "\nNodes: " << model_out.Mesh().Nodes();
      cout << "\nElements: " << model_out.Mesh().Elements();
      cout << "\nFaces: " << model_out.Mesh().Faces();
      cout << "\nInterfaces: " << model_out.Mesh().InterFaces();
      cout << "\nFinish " << dimension << " SplitBoundary Test: SplitBoundary around Regions\n";
    }
}




template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::Detect_and_create_splitboundaries( const string& model_name )
{
  const string variables_file( "SplitBoundary_Test-variables.txt" );

  // 1. convert ansys model into CSMP model
  Model<dim>* model = NULL;
  if constexpr ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  else if constexpr ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));

  // 2. build CSMP SplitBoundary
  // ---------------------------------------------------------------------------------------
  model->DetectAndCreateSplitBoundaries();

  // 3. see whether the split boundary survives being writting to and recovered from file
  // ------------------------------------------------------------------------------------------
  model->OutputToBinaryFile( model_name.c_str() );

  csmp::Model<dim> model_out( model_name );
  cout << "\nNodes: " << model_out.Mesh().Nodes() << "\n";
  cout << "\nElements: " << model_out.Mesh().Elements() << "\n";
  cout << "\nFaces: " << model_out.Mesh().Faces() << "\n";
  cout << "\nInterfaces: " << model_out.Mesh().InterFaces() << "\n";

  // 4. visualising
  string test_name( "DETECTED_SPLITBOUNDARY_TEST_FROM_" );
  test_name += model_name;
  VisualiseSplitBoundaries( model_out, test_name );

  return;
}



template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::Detect_and_create_splitboundaries_from_constructor( const string& model_name )
{
  const string variables_file( "SplitBoundary_Test-variables.txt" );
  
  // 1. convert ansys model into CSMP model
  Model<dim>* model = NULL;
  if constexpr ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  else if constexpr ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true ));
  
  // 2. create lower - dimensional stand - alone meshes from SplitBoundary objects, and
  //    insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  int32_t material_id{1};
  auto new_regions = model->InsertLowerDimensionalRegionsIntoSplitBoundaries( material_id );

  // 3. see whether the split boundary survives being writting to and recovered from file
  // ------------------------------------------------------------------------------------------
  model->OutputToBinaryFile( model_name.c_str() );

  csmp::Model<dim> model_out( model_name );
  cout << "\nNodes: " << model_out.Mesh().Nodes() << "\n";
  cout << "\nElements: " << model_out.Mesh().Elements() << "\n";
  cout << "\nFaces: " << model_out.Mesh().Faces() << "\n";
  cout << "\nInterfaces: " << model_out.Mesh().InterFaces() << "\n";

} // endf Detect_and_create_splitboundaries_from_constructor




template<uint32_t dim>
void SplitBoundaryInterface_Test<dim>::run()
{

  //THE ONLY RIGOROUS TEST
  Test_splitboundary_from_lower_dim_region();

  //TESTS WHICH RUN BUT HAVE NO DIAGNOSTICS
  if constexpr( dim == 2U ) {
      Test_splitboundary_between_regions( "BoxHalfs2D" );
      Test_splitboundary_between_regions( "ThreeZones2D" );

      ///E.P These remaining tests are not working yet
      /*
      // JC: working on the QC process which is requried for the following models
      Test_splitboundary_around_regions( "UnitSquareFracs_xline" );
      Test_splitboundary_around_regions( "UnitSquareFracs_yline" );
      Test_splitboundary_around_regions( "UnitSquareFracs_orthogonal" );
      Test_splitboundary_around_regions( "UnitSquareFracs_irregular" );
      Test_splitboundary_between_regions( "kueper_one_interface" );
      Test_splitboundary_between_regions( "lens2D" ); //added
      */
 }
    
  // test splitboundary between 3D regions
  if constexpr( dim == 3U ) {
      Test_splitboundary_between_regions( "BoxHalfs3D" );

      ///E.P These remaining tests are not working yet
      /*
      Test_splitboundary_between_regions( "ThreeZones3D" );
      // test splitboundary from constructor of ansys model
      // JC: working on the QC process which is requried for the following models
      //detect_and_create_splitboundaries_from_constructor<2U>( "Jura-slope1" );
      Detect_and_create_splitboundaries_from_constructor( "Dyke_Split" );
      // test splitboundary for complex ansys models
      Test_splitboundary_between_regions( "lamination" ); VERY SLOWWWW and negative jacobians warning ...
      */
    }

}

template class SplitBoundaryInterface_Test<2U>;
template class SplitBoundaryInterface_Test<3U>;


} // csmp
