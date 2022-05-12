#include "ANSYS_Model2D_Test.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "MeshManagementUtilities.h"
#include "vsetMakers.h"

#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "DenseMatrix.h"
#include "plf_colony.h"

using namespace std;

namespace csmp {

// for model building from ANSYS
 static void create_ANSYS2D_Model( bool reconstruct_from_file )
 {
    string model2d_name_ = "box2d_fault";
    string varFileName = "CSMP-variables.txt";
    Model<2>* model2d_ = new ANSYS_Model2D(model2d_name_.c_str(), varFileName.c_str());

    // ansys 2d model - contiguous
    cout << "\n------------------------------------------";
    cout << "\nMeshManager_Test: ANSYS model 'box2d_fault'";
    cout << "\n------------------------------------------";
    cout << "\nNodes: " << model2d_->Mesh().Nodes() << "\n";
    set<Element<2>*> elements;
    cout << "\nInterconnected elements: " << findContiguousMeshPatch<2,Element>( &(*model2d_->Mesh().ElementsBegin()), elements ) << "\n";
    cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
    std::map<std::string,std::vector<Element<2U>*> > patch_map;
    cout << "\nElement Groups: " << findContiguousMeshPatches( model2d_->Mesh().ElementsBegin(), model2d_->Mesh().ElementsEnd(), patch_map ) << "\n";
    cout << "\nFaces: " << model2d_->Mesh().Faces() << "\n";
    std::map<std::string,std::vector<Face<2U>*> >  face_map;
    cout << "\nFace Groups: " << findContiguousMeshPatches( model2d_->Mesh().FacesBegin(), model2d_->Mesh().FacesEnd(), face_map ) << "\n";
    cout << "\nInterfaces: " << model2d_->Mesh().InterFaces() << "\n";
    std::map<std::string,std::vector<InterFace<2U>*> >  iface_map;
    cout << "\nInterface Groups: " << findContiguousMeshPatches( model2d_->Mesh().InterFacesBegin(), model2d_->Mesh().InterFacesEnd(), iface_map ) << "\n";
    
    if (reconstruct_from_file) {
        model2d_->OutputToBinaryFile(model2d_name_.c_str());
        delete model2d_;
        model2d_ = new Model<2U>(model2d_name_);
        MeshManager<2>& mesh(model2d_->Mesh());
        cout << "\nNodes: " << mesh.Nodes() << "\n";
        cout << "\nNode Groups: " << findContiguousMeshPatch<2,Element>( &(*mesh.ElementsBegin()), elements ) << "\n";
        cout << "\nElements: " << model2d_->Mesh().Elements() << "\n";
        cout << "\nElement Groups: " << findContiguousMeshPatches( mesh.ElementsBegin(), mesh.ElementsEnd(), patch_map ) << "\n";
        cout << "\nFaces: " << model2d_->Mesh().Faces() << "\n";
        cout << "\nFace Groups: " << findContiguousMeshPatches( mesh.FacesBegin(), mesh.FacesEnd(), face_map ) << "\n";
        cout << "\nInterfaces: " << model2d_->Mesh().InterFaces() << "\n";
        cout << "\nInterface Groups: " << findContiguousMeshPatches( mesh.InterFacesBegin(), mesh.InterFacesEnd(), iface_map ) << "\n";
      }
      
 } // end create_ANSYS2D_Model
 
 
 

/**
    Using a whole suite of 2D Ansys models with boundaries and even split  boundaries being created :
    
    - BoxHalfs2D
    - Fluid_Flower
        
*/
void ANSYS_Model2D_Test::run()
  {
    const bool verbose(false);
    
    ANSYS_Model2D model( "BoxHalfs2D", "CSMP-variables.txt" );
    Region<2>& rref( model.Region( "Model" ) );

    DenseMatrix<DM_MIN> dm;
    const auto elementsEnd( rref.CellsEnd() );
    try
      {    
        for( auto it = rref.CellsBegin(); it != elementsEnd; ++it )
          (*it)->CoordinateMatrix();
      }
    catch(...)
      {
        _test(false);
      }

    if ( verbose ) cout << "\nModel Node Count: " << rref.Nodes() << endl;
    _test( rref.Nodes() == 106 );

    const Index nodalKey( model.Database().StorageKey( "nodal variable" ) );
    rref.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0. ) );
    size_t boundaryNodeCount( 0 );
    const vector<Node<2>*>::const_iterator nodesEnd( rref.NodesEnd() );
    for( vector<Node<2>*>::const_iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
      {
        const BOX_BOUNDARY boxBoundary( (*it)->AtBoundary() );
        if( boxBoundary != NOT )
          (*it)->Store( nodalKey, makeScalar( PLAIN, 1.0 ) );
        else
          continue;   
        ++boundaryNodeCount;
      }
    if ( verbose ) cout << "\nBoundary Node Count: " << boundaryNodeCount << endl;
    _test ( boundaryNodeCount == 40 );
    if ( verbose ) {
         VTU_Interface<2> vtu( model );
         vtu.OmitZeroInFileName( true );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-BoxBoundary", "nodal variable", "Model", static_cast<int>(0) );
      }
    _test( model.ContainsBoundary( "LEFT" ) );
    _test( model.ContainsBoundary( "RIGHT" ) );
    _test( model.ContainsBoundary( "BOTTOM" ) );
    _test( model.ContainsBoundary( "TOP" ) );

    _test( model.ContainsRegion( "MATRIX_LEFT" ) );
    _test( model.ContainsRegion( "MATRIX_RIGHT" ) );


    Boundary<2>& left( model.Boundary( std::string("LEFT") ) );
    Boundary<2>& right( model.Boundary( std::string("RIGHT") ) );
    Boundary<2>& bottom( model.Boundary( std::string("BOTTOM") ) );
    Boundary<2>& top( model.Boundary( std::string("TOP") ) );
    
    // testing that the corner nodes have been flagged correctly
    const Region<2>& mref(model.Region("Model"));
    size_t  counter(0U);
    for ( auto it=mref.PerimeterNodesBegin(); it!=mref.NodesEnd(); it++ ) {
         if ( (*it)->AtBoundary() == CNR1 ) counter++;
         else if ( (*it)->AtBoundary() == CNR2 ) counter++;
         else if ( (*it)->AtBoundary() == CNR3 ) counter++;
         else if ( (*it)->AtBoundary() == CNR4 ) counter++;
      }
    _test( counter == 4U );

    right.InputPropertyValue( "face variable", makeScalar( PLAIN, 2. ) );
    ArrayVariable av( "nodal array", model.Database() );
    // SKM FIX: if av is not initialized its values are NAN and a comparison with another array variable will always evaluate as false
    av = 5.;
    model.InputPropertyValue( "nodal array", av );
      
    const size_t leftNodes( left.Nodes() );
    const size_t rightNodes( right.Nodes() );
    const size_t bottomNodes( bottom.Nodes() );
    const size_t topNodes( top.Nodes() );
    const size_t leftFaces( left.Cells() );
    const size_t rightFaces( right.Cells() );
    const size_t bottomFaces( bottom.Cells() );
    const size_t topFaces( top.Cells() );
    
    if ( verbose ) cout << "\nModel Node Count: " << rref.Nodes() << endl;
    _test( rref.Nodes() == 106 );

    string bin1name("ANSYS2D_bin");
    model.OutputToBinaryFile(bin1name.c_str());
    // model.Out();

    Model<2U> modelBinIn0(bin1name);
    Index nodalArrayKey0( modelBinIn0.Database().StorageKey("nodal array") );
    Index faceVariableKey0( modelBinIn0.Database().StorageKey("face variable") );
    _test( (*modelBinIn0.Boundary("RIGHT").CellsBegin())->Read(faceVariableKey0) == 2. );
    ArrayVariable avBin0( "nodal array", modelBinIn0.Database() );
    (*modelBinIn0.Region("Model").NodesBegin())->Read( nodalArrayKey0, avBin0 );
    _test( avBin0 == av );
    
    Model<2U> modelBinIn1(bin1name);
    Index nodalArrayKey( modelBinIn1.Database().StorageKey("nodal array") );
    Index faceVariableKey( modelBinIn1.Database().StorageKey("face variable") );
    _test( leftNodes == modelBinIn1.Boundary("LEFT").Nodes() );
    _test( rightNodes == modelBinIn1.Boundary("RIGHT").Nodes() );
    _test( topNodes == modelBinIn1.Boundary("TOP").Nodes() );
    _test( bottomNodes == modelBinIn1.Boundary("BOTTOM").Nodes() );
    _test( leftFaces == modelBinIn1.Boundary("LEFT").Cells() );
    _test( rightFaces == modelBinIn1.Boundary("RIGHT").Cells() );
    _test( topFaces == modelBinIn1.Boundary("TOP").Cells() );
    _test( bottomFaces == modelBinIn1.Boundary("BOTTOM").Cells() );

    _test( (*modelBinIn1.Boundary("RIGHT").CellsBegin())->Read(faceVariableKey) == 2. );
    ArrayVariable avBin( "nodal array", modelBinIn1.Database() );
    (*modelBinIn1.Region("Model").NodesBegin())->Read( nodalArrayKey, avBin );
    _test( avBin == av );
    
    // ADVANCED FUNCTIONALITY TESTS in relation to ANSYS models
    // --------------------------------------------------------
    // SKM (27/3/22)
    // assunming that the line-element connectivity was rebuilt when the model was created
    Test_printLineElementRegion();
    
    Test_CreateConsistentLineElementOrientations2D();
    
    Test_CreateInternalBoundary();

    Test_CreateSplitBoundaries();

    Test_CreateSplitBoundariesBetweenUniqueRegions();

  } // end run
  

/*  for test development
 
void ANSYS_Model2D_Test::run()
  {
    // Test_printLineElementRegion();
    // Test_CreateConsistentLineElementOrientations2D();
    // Test_CreateSplitBoundaries();
    Test_CreateSplitBoundariesBetweenUniqueRegions();
  }
*/
 
 
 
 
       // using line element VSet
void ANSYS_Model2D_Test::Test_printLineElementRegion()
 {
    if ( verbose_ ) {
         cout <<"\n\n"<<"ANSYS_Model2D_Test::Test_printLineElementRegion: running test on model '"<< endl;
         cout << "MeshPatchWithLineElements" <<"'"<< endl;
         cout.flush();
      }
    
    // verifying function with predefined correct dataset
    {
      VSet<2U>      vset;
      ModelTopology topo = test_Create_MeshPatchWithLineElements_VSet( vset );
      Model<2U>     model( topo, vset, "CSMP-1phase-variables.txt", true );
       
      // works fine
      const bool renumber_nodes{false};
      _test( printLineElementRegion( model, "FRAC1", renumber_nodes ) == 3 );
      _test( printLineElementRegion( model, "FRAC2", renumber_nodes ) == 3 );
      _test( printLineElementRegion( model, "FRAC3", renumber_nodes ) == 1 ); // only one element
    }
 
    // rebuilding neighbor connectivty and line element connectivity
    {
      VSet<2U>      vset;
      ModelTopology topo = test_Create_MeshPatchWithLineElements_VSet( vset );
      vset.EstablishElementConnectivity2D(); // also deals with line-element orientations
      Model<2U>     model( topo, vset, "CSMP-1phase-variables.txt", true );
       
      const bool renumber_nodes{false};
      _test( printLineElementRegion( model, "FRAC1", renumber_nodes ) == 3 );
      _test( printLineElementRegion( model, "FRAC2", renumber_nodes ) == 3 );
      _test( printLineElementRegion( model, "FRAC3", renumber_nodes ) == 1 ); // only one element
    }
    
} // end Test_printLineElementRegion

 
 
 


void  ANSYS_Model2D_Test::Test_CreateConsistentLineElementOrientations2D()
 {
    string model2d_name_ = "three_layers"; // TODO: use model that is already in the testing fixtures
    if ( verbose_ ) {
         cout <<"\n\n"<<"ANSYS_Model2D_Test::Test_CreateConsistentLineElementOrientations2D: running test on model '"<< endl;
         cout << model2d_name_ <<"'"<< endl;
         cout.flush();
      }

    string varFileName = "CSMP-variables.txt";
    ANSYS_Model2D  model( model2d_name_.c_str(), varFileName.c_str() );
     
    // 1. accessing the line-element regions representing the boundaries between layers
    Region<2U>&  interface1{ model.Region("INTERFACE1") }, interface2{ model.Region("INTERFACE2") };
    
    // 2. checking whether the elements are forming a chain that can be traversed neighbor to neighbors
    //    (traversing until end is discovered, but terminating after at most as many steps that the region has elements)
    if ( verbose_ ) {
         interface1.Out();
         const bool renumber_elmts{true};
         _test( printLineElementRegion( model, interface1.Name().c_str(), renumber_elmts ) == interface1.Cells() );
         interface2.Out();
         _test( printLineElementRegion( model, interface2.Name().c_str(), renumber_elmts ) == interface2.Cells() );
      }
    
    // 3. The endpoints of these regions must be at the vertical model boundaries
    _test( interface1.PerimeterCells() == 2 );
    _test( interface2.PerimeterCells() == 2 );
    _test( interface1.PerimeterNodes() == 2 );
    _test( interface2.PerimeterNodes() == 2 );
    
    // 4. Are the regions contiguous (all elements are interconnected except for those at the end
    int contiguous{2}, n_missing_nbors{0};
    // interface 1
    for ( auto it=interface1.CellsBegin(); it!=interface1.CellsEnd(); ++it )
      for ( int i{0}; i<(*it)->Neighbors(); ++i )
        if ( (*it)->Neighbor(i) == nullptr )
          n_missing_nbors++;
    _test( n_missing_nbors <= contiguous );
    // interface 2
    n_missing_nbors = 0;
    for ( auto it=interface2.CellsBegin(); it!=interface2.CellsEnd(); ++it )
      for ( int i{0}; i<(*it)->Neighbors(); ++i )
        if ( (*it)->Neighbor(i) == nullptr )
          n_missing_nbors++;
    _test( n_missing_nbors <= contiguous );

    // 4. getting pointers to the elements at the opposite ends of the line-element sequence and checking these elements
    // interface 1
    {
      auto eit1{ *interface1.PerimeterCellsBegin() };
      auto eit2{ (*next(interface1.CellsEnd(),-1)) };
      _test( eit1->ConnectedNeighbors() == 1 );
      _test( eit2->ConnectedNeighbors() == 1 );
      // establishing the direction in which the line-element sequence is to be traversed
      // should be consistent for beginnging and end
      // (perimeter element 1 is expected to be at beginning)
      bool forward_at_beginning = ( eit1->Neighbor(1) == nullptr ) ? true : false;
      bool forward_at_end       = ( eit2->Neighbor(0) == nullptr ) ? true : false;
      _test( forward_at_beginning == forward_at_end );
    }
    // interface 2
    {
      auto eit1{ *interface2.PerimeterCellsBegin() };
      auto eit2{ (*next(interface2.CellsEnd(),-1)) };
      _test( eit1->ConnectedNeighbors() == 1 );
      _test( eit2->ConnectedNeighbors() == 1 );
      // establishing the direction in which the line-element sequence is to be traversed
      // (perimeter element 1 is expected to be at beginning)
      bool forward_at_beginning = ( eit1->Neighbor(1) == nullptr ) ? true : false;
      bool forward_at_end       = ( eit2->Neighbor(0) == nullptr ) ? true : false;
      _test( forward_at_beginning == forward_at_end );
    }
        
 } // end Test_CreateConsistentLineElementOrientations2D





/**
    For model   'three_layers'   converts the line-element regions INTERFACE1 and INTERFACE2 into internal boundaries
*/
void  ANSYS_Model2D_Test::Test_CreateInternalBoundary()
{
    string model2d_name_ = "three_layers"; // TODO: use model that is already in the testing fixtures
    if ( verbose_ ) {
         cout <<"\n\n"<<"ANSYS_Model2D_Test::Test_CreateInternalBoundary: running test on model '";
         cout << model2d_name_ <<"'"<< endl;
         cout.flush();
      }
    string varFileName = "CSMP-variables.txt";
    ANSYS_Model2D model( model2d_name_.c_str(), varFileName.c_str() );
     
    // 1. accessing the line-element regions representing the boundaries between layers
    const Region<2U>& interface1{ model.Region("INTERFACE1") }, interface2{ model.Region("INTERFACE2") };
    size_t n_elmts_region1{ interface1.Cells() };
    size_t n_elmts_region2{ interface2.Cells() };
    
    // removes input region
    pair<set<string>,bool> boundaryName1 = model.CreateInternalBoundaryFrom( "INTERFACE1" );
    const Boundary<2U>& boundary1(model.Boundary( (*(boundaryName1.first).begin()) ) );
    _test( boundary1.Cells() == n_elmts_region1 );
    
    pair<set<string>,bool> boundaryName2 = model.CreateInternalBoundaryFrom( "INTERFACE2" );
    const Boundary<2U>& boundary2(model.Boundary( (*(boundaryName2.first).begin()) ) );
    _test( boundary2.Cells() == n_elmts_region2 );

} // end Test_CreateInternalBoundary





/**
    For model   'three_layers'   converts the line-element regions INTERFACE1 and INTERFACE2 into internal boundaries
*/
void  ANSYS_Model2D_Test::Test_CreateSplitBoundaries()
{
    if ( verbose_ ) {
         cout <<"\n\n"<<"ANSYS_Model2D_Test::Test_CreateSplitBoundaries: running test..."<< endl;
         cout.flush();
      }
    string model2d_name_ = "three_layers"; // TODO: use model that is already in the testing fixtures
    string varFileName = "CSMP-variables.txt";
    ANSYS_Model2D model( model2d_name_.c_str(), varFileName.c_str() );
    const size_t n_original_cells = model.Mesh().Elements() + model.Mesh().Faces();  // only those
    const size_t n_original_faces = model.Mesh().Faces();
     
    // 1. creating the SplitBoundary from lower-dimensional region
    size_t n_elmts_region1 = model.Region( "INTERFACE1" ).Cells();
    size_t n_elmts_region2 = model.Region( "INTERFACE2" ).Cells();

    pair<set<string>,bool> splitBoundaryName1 = model.CreateSplitBoundaryFrom( "INTERFACE1" );
    pair<set<string>,bool> splitBoundaryName2 = model.CreateSplitBoundaryFrom( "INTERFACE2" );
    assert( splitBoundaryName1.first.size() == 1 );
    assert( splitBoundaryName2.first.size() == 1 );
    _test( model.SplitBoundary( (*splitBoundaryName1.first.begin()) ).Cells() == n_elmts_region1 );
    _test( model.SplitBoundary( (*splitBoundaryName2.first.begin()) ).Cells() == n_elmts_region2 );
    // the number of total cells in the model should not have changed
    size_t volume_elmts{0U}, surface_elmts{0U}, line_elmts{0U};
    _test( currentCellTypes( model.Mesh(), ELEMENT, volume_elmts, surface_elmts, line_elmts ) == n_original_cells );
    // now the line elements should be gone and Interfaces there instead
    _test( line_elmts == 0 );
    // the number of Faces should have stayed constant because the once created new must have been removed again
    currentCellTypes( model.Mesh(), FACE, volume_elmts, surface_elmts, line_elmts );
    _test( line_elmts == n_original_faces );
    // there should be x interfaces
    currentCellTypes( model.Mesh(), INTER_FACE, volume_elmts, surface_elmts, line_elmts );
    _test( line_elmts == n_elmts_region1 + n_elmts_region2 );
    
    model.SplitBoundariesOut();

    // 2. remove split boundary 1 here before creating new ones in the same place
    const bool erase_interfaces{ true };
    model.RemoveSplitBoundary( (*splitBoundaryName1.first.begin()).c_str(), erase_interfaces ); // INTERFACE1
        
    model.RegionsOut();

    // 3. Directly creating a split boundary between regions (should not work leading to detection that there already is a split boundary)
    model.CreateSplitBoundaryBetween( "LOWER_REGION", "MIDDLE_REGION" ); // INTERFACE 2

    // 4. Creating a split boundary between regions (detecting the disjointed nodes?)
    ANSYS_Model2D model2( model2d_name_.c_str(), varFileName.c_str() );
    // now the nodes are shared so this should work
    model2.CreateSplitBoundaryBetween( "MIDDLE_REGION", "UPPER_REGION" ); // INTERFACE 1

} // end Test_CreateSplitBoundaries



/**
       Direct creation of SplitBoundaries between all the unique regions of the model.
*/
void  ANSYS_Model2D_Test::Test_CreateSplitBoundariesBetweenUniqueRegions()
{
    if ( verbose_ ) {
         cout <<"\n\n"<<"ANSYS_Model2D_Test::Test_CreateSplitBoundaries: running test..."<< endl;
         cout.flush();
      }
    string model2d_name_ = "Fluid_Flower"; // TODO: use model that is already in the testing fixtures
    string varFileName = "CSMP-variables.txt";
    ANSYS_Model2D model( model2d_name_.c_str(), varFileName.c_str() );
    model.InputPropertyValue( "nodal variable", makeScalar(ANY,0.) );
    VTU_Interface<2> vtu( model );
    vtu.OmitZeroInFileName( true );
     
    // some initial tests on the model
    vector<pair<pair<Element<2U>*,uint32_t>,pair<Element<2U>*,uint32_t> > >  matching_cells;
    vector<Node<2U>*>  matching_nodes;
    size_t n_nodes = sharedPerimeterNodes( model.Region("F_T"), model.Region("E_BB_T"), matching_nodes );
    size_t n_cells = sharedPerimeterCells( model.Region("F_T"), model.Region("E_BB_T"), matching_cells );
    // in this case the regions touch at single nodes do not share a parent element
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_F_T", "nodal variable", "F_T", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_E_BB_T", "nodal variable", "E_BB_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("D_BB"), model.Region("C_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("D_BB"), model.Region("C_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_D_BB", "nodal variable", "D_BB", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_C_T", "nodal variable", "C_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("D_BB"), model.Region("E_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("D_BB"), model.Region("E_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_D_BB", "nodal variable", "D_BB", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_E_T", "nodal variable", "E_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("F_BB"), model.Region("E_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("F_BB"), model.Region("E_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_F_BB", "nodal variable", "F_BB", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_E_T", "nodal variable", "E_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("E_BB_T"), model.Region("D_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("E_BB_T"), model.Region("D_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_E_BB_T", "nodal variable", "E_BB_T", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_D_T", "nodal variable", "D_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("C_BB_T"), model.Region("D_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("C_BB_T"), model.Region("D_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_C_BB_T", "nodal variable", "C_BB_T", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_D_T", "nodal variable", "D_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("C_BB_T"), model.Region("ESF_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("C_BB_T"), model.Region("ESF_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_C_BB_T", "nodal variable", "C_BB_T", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_ESF_T", "nodal variable", "ESF_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("C_T"), model.Region("ESF_BB_T"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("C_T"), model.Region("ESF_BB_T"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_ESF_BB_T", "nodal variable", "ESF_BB_T", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_C_T", "nodal variable", "C_T", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("ESF"), model.Region("F_BA"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("ESF"), model.Region("F_BA"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_ESF", "nodal variable", "ESF", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_F_BA", "nodal variable", "F_BA", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("ESF_BA"), model.Region("F"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("ESF_BA"), model.Region("F"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_ESF_BA", "nodal variable", "ESF_BA", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_F", "nodal variable", "F", static_cast<int>(0) );
      }
    n_nodes = sharedPerimeterNodes( model.Region("F"), model.Region("G_BA"), matching_nodes );
    n_cells = sharedPerimeterCells( model.Region("F"), model.Region("G_BA"), matching_cells );
    if ( verbose_ && ( n_nodes==2U && n_cells==0U ) ) {
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region1_F", "nodal variable", "F", static_cast<int>(0) );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-region2_G_BA", "nodal variable", "G_BA", static_cast<int>(0) );
      }
     
    // 1. creating the SplitBoundary objects everywhere
    /* Observations
       - super slow!
       - most computational effort goes into UpdateConnectivity and PartitionCellVector
    */
    size_t n_split_boundaries = model.SeparateUniqueRegionsBySplitBoundaries();
    _test( n_split_boundaries < 93 );
    model.SplitBoundariesOut();
    
    // using 'nodal variable' to visualise which nodes are manifolds
    Region<2U>  model_domain = model.Region("Model");
    const csmp::Index var_key = model.Database().StorageKey("nodal variable");
    
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit )
      if ( (*nit)->IsManifold() )
        // there should be as many branches as materials come together
        (*nit)->Store( var_key, makeScalar(PLAIN,(*nit)->Manifold()->Branches()));
     
    // visualisation
    if ( verbose_ ) {
         vtu.OmitZeroInFileName( true );
         vtu.OutputDataToVTU( "ANSYS_Model2D_Test-node_manifolds", "nodal variable", "Model", static_cast<int>(0) );
      }

    // TODO: perform some testing
    cout <<"\n"<<"ANSYS_Model2D_Test::Test_CreateSplitBoundariesBetweenUniqueRegions: completed successfully"<< endl;

} // end Test_CreateSplitBoundaries



} // csmp
