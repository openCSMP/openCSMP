#include "Boundary_Test.h"

#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"
#include "Boundary.h"
#include "Region.h"
#include "PL_Utilities.h"
#include "VTU_Interface.h"
#include "ANSYS_Interface.h"
#include "variableOperations.h"


using namespace std;

namespace csmp{

template<uint32_t dim>
void innerOuterParents( Model<dim>& model, VTU_Interface<dim>& vtu, Boundary<dim>& boundary, string fileName )
{
  const ScalarVariable zero( PLAIN, 0. );
  const ScalarVariable one( PLAIN, 1. );
  const ScalarVariable two( PLAIN, 2. );
  Index elVarKey( model.Database().StorageKey("element variable") );
  model.InputPropertyValue( "element variable", zero);
  for( typename vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != boundary.CellsEnd(); ++it )
    {
      (*it)->Parent(OUTSIDE)->Store( elVarKey, two );
      (*it)->Parent(INSIDE)->Store( elVarKey, one );
    }
  vtu.OutputDataToVTU( fileName.data(), "element variable", "Model", static_cast<int>(0) );
}



  /*
  face var in var file does not crash - is it in faces storage?
  accessing an element var in face does not work - that's ok
  accessing an face var in face does not work - not ok
  face local variable storage is not correctly initialized if we want to access face vars
  face lvs empty - fix initialization. - DONE ok (for CreateBetween() only!!!)
   - initialization if DB from var file seems fine - ok
   - initialization in Boundary::Create...() missing ! - DONE ok
   - do we want CPVS in face?? do we even want FACE vars?? - yep, hence ctors. - DONE ok
  implement AddFaces()
  Face inherits Read() (as does InterFace for that matter), so asserts for ELEMENT placement.
    do we want FACE vars only? or eELEMENT vars only? or extra LVS and keep both? before there was empty LVS whatsoever
    FACE placement has to be asserted external!!  - ok
  UpdateMemberIndices() was missing in CreateBetween() - DONE ok
  CreateBetween() tested ok P. Lang 28/08/2011
  CreateFrom only sets BOX_BOUNDARY flags for elements and INNER nodes, so none for nodes at EDGES !!
  CreateNodeVector() in Face::Assign( Element, bool ) was missing - DONE ok
  CreateNodeVector() did not account for single parent element, hence new Assign implementation for single parent
  CreateBetweenAllElements() does not work, for the time being it is removed as an interface and became a private method
  No boundary flags assigned to edges of boundaries yet. maybe as new parameter 'bool assignEdgeFlags' in CreateFrom and CreateFor
  CreateFor used size_t i as running variable in two nested for loops - fixed ok
  Test CreateFor !!
  Test CreateFrom for internal surface
  Test 2D
  Remove Face::Assign(elmt,bool) - will never work, only used in CreateBetweenAllElements()
  In 2D, ANSYS interface does not create regions for LEFT, RIGHT etc... 3D??
  should the InsertBoundary methods return a reference to the boundary??
  the boundary and region, split boundary should not themselves be responsible for their creation. pure container class - ops/mods - outsource
  2D CreateBetween Test
  CreateBetweenAllElements
  ANSYS_Interface (at least in 2D) removes boundary regions for box shape, so calling RectShape returns false consecutively
  */


// this also includes Face<dim> tests
void Boundary_Test::run()
{
  runLegacy();
  runCurrent();
  UnitNormalTest3D();
}


template<uint32_t dim>
void Boundary_Test::ElementNodes( const Region<dim>& region )
  {
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.CellsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator element( region.CellsBegin() ); element != elementsEnd; ++element )
    _test( (*element)->Nodes() > 1 );
  }


template<uint32_t dim>
void Boundary_Test::NoSurfaceElementsAsNodeParents( const Region<dim>& region )
{
  const auto domainNodesEnd( region.NodesEnd() );
  for( auto it = region.NodesBegin(); it != domainNodesEnd; ++it )
    for( auto i{0U}; i < (*it)->Parents(); ++i )
      _test( !(*it)->Parent(i)->IsSurface() );
}


template <uint32_t dim>
void Boundary_Test::TestBoxBoundary( Model<dim>& model, const string& boxBoundary, VTU_Interface<dim>& vtu )
{
  std::string boundary_name( boxBoundary );
  _test( model.ContainsBoundary( boundary_name ) );
  Boundary<dim>& boundary( model.Boundary( boundary_name ) ); 
  _test( InputElementAreaAsVolumeVariable<dim>( model, boundary, "face variable" ) > 0 );
  if( dim == 3 )
    vtu.OutputDataToVTU( "TestFaceVariable3D", "face variable", boundary, static_cast<int>(0) );
  else
    vtu.OutputDataToVTU( "TestFaceVariable2D", "face variable", boundary, static_cast<int>(0) );
  CheckFaceNeighbors(boundary);
  BOX_BOUNDARY boxBoundaryFlag;
  if( boxBoundary == "TOP" )
    boxBoundaryFlag = TOP;
  else if( boxBoundary == "BOTTOM" )
    boxBoundaryFlag = BOTTOM;
  else if( boxBoundary == "LEFT" )
    boxBoundaryFlag = LEFT;
  else if( boxBoundary == "RIGHT" )
    boxBoundaryFlag = RIGHT;
  else if( boxBoundary == "FRONT" )
    boxBoundaryFlag = FRONT;
  else if( boxBoundary == "BACK" )
    boxBoundaryFlag = BACK;
  else
    throw csmp::Exception(ERROR, "Boundary_Test::TestBoxBoundary", "Box boundary not recognized!" );
  CheckNodeFlags( boundary, boxBoundaryFlag, true );
  CheckFaceUnitNormalOrientation( boundary );
  CheckNodeParents( boundary );

  return;
}


template <uint32_t dim>
size_t Boundary_Test::InputElementAreaAsVolumeVariable( Model<dim>& model, Boundary<dim>& boundary, const char* variableName )
{
  Index areaKey( model.Database().StorageKey( variableName ) );
  assert( areaKey.place == FACE ); // this needs to be done for Face::Read bc Read is inherited from Element and can
                                   // either be Face, Element or InterFace
  ScalarVariable area( PLAIN, 0. );
  size_t surfaceElementCount( 0 );
  const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd( boundary.CellsEnd() );
  for( typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it )
  {
    area = (*it)->Area();
    (*it)->Store( areaKey, area );
    ++surfaceElementCount;
  }
  

  return surfaceElementCount;
  }


/// checks that at least one neighbor is present
template <uint32_t dim>
void Boundary_Test::CheckFaceNeighbors( const Boundary<dim>& boundary )
 {
    const auto domainElementsEnd( boundary.CellsEnd() );
    for( auto it = boundary.CellsBegin(); it != domainElementsEnd; ++it )
      {
        int notNullNeighbors{0};
        const auto neighbors( (*it)->Neighbors() );
        for( auto i = 0; i < neighbors; ++i )
          {
            if( (*it)->Neighbor(i) )
              ++notNullNeighbors;
          }
        _test( notNullNeighbors > 0 );
        // diagnostics
        if ( notNullNeighbors == 0 ) {
             cerr <<"\n"<< boundary.Name() <<": Face "<< (*it)->Idx() <<": ";
             cerr << parseFiniteElementType( (*it)->FE_Type() ) <<"  ";
          }
      }
  }


/// checks for that face unit normal points toward outer parent element
template <uint32_t dim>
void Boundary_Test::CheckFaceUnitNormalOrientation( const Boundary<dim>& boundary )
  {
	// IMPORTANT: this is ignored since the unit normal vector of that face such as the 2D-line face on a boundary
	//			  is not available. That means the function UnitNormal is not available for this case.
	return;

    VectorVariable<dim> unFace( PLAIN, 9999999. ), faceToInner( PLAIN, 9999999. );
    size_t              inward_pointing_normals(0U);
    const auto domainElementsEnd( boundary.CellsEnd() );
    for( auto it = boundary.CellsBegin(); it != domainElementsEnd; ++it )
      {
        // IMPORTANT - this is the method that is tested (CoordinateMatrix() is called inside)
        (*it)->UnitNormal( unFace );
        // creating a vector that points from face barycenter to higher-dimensional parent element barycenter
        Point<dim> bcFace = (*it)->BaryCenter();
        Point<dim> bcInner = (*it)->Parent(INSIDE)->BaryCenter();
        // constructing vector that points to inner element
        for( size_t d(0); d < dim; ++d )
            faceToInner(d) = bcInner[d] - bcFace[d];

        // the dotproduct of the inward pointing vector and the unit normal should be negative		
        _test( dotProduct(faceToInner,unFace) < 0. );
        // debugging diagnostics
        if ( dotProduct(faceToInner,unFace) > 0. ) {
             cerr <<"\nBoundary: '"<< boundary.Name() <<"', inward-pointing normal detected: "<< unFace <<"\n";
             (*it)->Out();
             inward_pointing_normals++;
        }
      }
    _test( inward_pointing_normals == 0 );
  }



/**
    Tests FiniteVolumePolicy:
       - UnitNormalToFace( size_t face, std::vector<double>& );
 
    Tests prism_test model because it contains elements of all
    types.
 
    @todo 2D model has to be tested as well
*/
void Boundary_Test::UnitNormalTest3D()
 {
// TODO: create a FiniteVolumePolicy_Test
     // ------------------------------------------------------------
     // 1. building model from ANSYS data files
     // ------------------------------------------------------------
      string  model_name("prism_test");
      // TODO: UnitNormalTest does not require any variables; remove property file
      ANSYS_Model3D  model( model_name.c_str(), "example25.txt");

     // ------------------------------------------------------------
     // 2. looping over all highest-dimensional elements
     //    testing whether normals are aligned with vectors
     //    between barycenter and face barycenters
     // ------------------------------------------------------------
     std::vector<double> unrml;
     const Region<3U>& model_domain(model.Region("Model"));
     if ( verbose_ ) cout <<"\nElement_Test::UnitNormalTest: testing normal directions...\n";
     for ( auto it=model_domain.CellsBegin(); it!=model_domain.CellsEnd(); ++it )
       {
          Point<3U> bctr((*it)->BaryCenter());
          // for all the faces of the element
          for ( auto face{0U}; face<(*it)->Faces(); ++face ) {
               // constructing a vector from element to face barycenter
               Point<3U> fbctr((*it)->FaceBaryCenter( face ));
               Point<3U> outward_vec(fbctr - bctr);
               // testing that the face unit normal is aligned with the outward
               // pointing vector
               (*it)->UnitNormalToFace( face, unrml );
               Point<3U> unitnormal(unrml);
               // the normals are aligned if dotproduct is positive
               double dotproduct = dotProduct( outward_vec, unitnormal );
               _test( dotproduct > 0. );
               if ( verbose_ and dotproduct < 0. ) {
                    cerr <<"\nunit normal to face "<< face <<" is inward pointing:";
                    (*it)->Out();
                 }
            }
       }
   
 } // end UnitNormalTest







/// checks for flags of the boundary nodes
template <uint32_t dim>
void Boundary_Test::CheckNodeFlags( const Boundary<dim>& boundary, BOX_BOUNDARY flag, bool interiorOnly )
  {
    const auto interiorDomainNodesEnd( boundary.PerimeterNodesBegin() );
    for( auto it = boundary.NodesBegin(); it != interiorDomainNodesEnd; ++it )
      {
        if( interiorOnly && it < interiorDomainNodesEnd )
          _test( (*it)->AtBoundary() == flag );
        else if( interiorOnly && it >= interiorDomainNodesEnd )
          _test( (*it)->AtBoundary() != NOT );
        else
          _test( (*it)->AtBoundary() != NOT );
      }
  }


/// checks for flags of the boundary nodes
template <uint32_t dim>
void Boundary_Test::CheckNodeParents( const Boundary<dim>& boundary )
  {
      const auto domainNodesEnd( boundary.NodesEnd() );
      for( auto it = boundary.NodesBegin(); it != domainNodesEnd; ++it )
      _test( (*it)->Parents() > 1 );
  }



void Boundary_Test::runLegacy()
  {
    const size_t DIM3 = 3U;
    // const scalars
    const ScalarVariable zero( PLAIN, 0. );
    const ScalarVariable one( PLAIN, 1. );
    const ScalarVariable two( PLAIN, 2. );

    // 2D TESTS I
    // ==========
    ANSYS_Model2D model2D( "BoxHalfs2D", "CSMP-variables.txt" );
    CheckFaceUnitNormalOrientation( model2D.Boundary("LEFT") );
    CheckFaceUnitNormalOrientation( model2D.Boundary("RIGHT") );
    CheckFaceUnitNormalOrientation( model2D.Boundary("TOP") );
    CheckFaceUnitNormalOrientation( model2D.Boundary("BOTTOM") );
    Region<2>& mref2D( model2D.Region( "Model" ) );
    mref2D.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 1.0 ) );
    mref2D.InputPropertyValue( "element variable", makeScalar( PLAIN, 2.0 ) );
    VTU_Interface<2> vtu2D( model2D ); vtu2D.OmitZeroInFileName( true );
    vtu2D.OutputDataToVTU( "TestNodeVariable2D", "nodal variable", "Model", static_cast<int>(0) );
    // RIGHT
    TestBoxBoundary<2>( model2D, "RIGHT", vtu2D );
    // TOP
    TestBoxBoundary<2>( model2D, "TOP", vtu2D );
    // BOTTOM
    TestBoxBoundary<2>( model2D, "BOTTOM", vtu2D );
    // LEFT
    TestBoxBoundary<2>( model2D, "LEFT", vtu2D );
    // BETWEEN
    pair<string,bool> result = model2D.CreateBoundaryBetween( "MATRIX_LEFT", "MATRIX_RIGHT" );
    _test( result.second == true );
    // this is the name according to the CSMP convention
    const string boundary2D_BETWEENname( string("BOUNDARY_MATRIX_LEFT") + string("_") + string("MATRIX_RIGHT") );
    _test( model2D.ContainsBoundary( boundary2D_BETWEENname ) );
    Boundary<2>& boundary2D_BETWEEN( model2D.Boundary( result.first ) ); 
    _test( InputElementAreaAsVolumeVariable<2>( model2D, boundary2D_BETWEEN, "face variable" ) > 0 );
    vtu2D.OutputDataToVTU( "TestFaceVariable2D_BETWEEN", "face variable", boundary2D_BETWEEN, static_cast<int>(0) );
    CheckFaceNeighbors(boundary2D_BETWEEN);
    CheckNodeFlags( boundary2D_BETWEEN, IRREGULAR );
    CheckFaceUnitNormalOrientation( boundary2D_BETWEEN );
    CheckNodeParents( boundary2D_BETWEEN );
    innerOuterParents( model2D, vtu2D, boundary2D_BETWEEN, "Parents3" );

    ///  getting rid of the boundary so that a new one can be created
    const bool erase_faces{true};
    model2D. RemoveBoundary( boundary2D_BETWEENname.c_str(), erase_faces );
    
    pair<set<string>,bool>  result2 = model2D.CreateInternalBoundaryFrom( "STANDARD" );
    _test( result2.first.size() == 1 );
    string boundary2D_FROMname( (*result2.first.begin()) );
    _test( model2D.ContainsBoundary( boundary2D_FROMname ) );
    Boundary<2>& boundary2D_FROM( model2D.Boundary( boundary2D_FROMname ) );
    boundary2D_FROM.Out();
    _test( InputElementAreaAsVolumeVariable<2>( model2D, boundary2D_FROM, "face variable" ) > 0 );
    vtu2D.OutputDataToVTU( "TestFaceVariable2D_INTERNAL", "face variable", boundary2D_FROM, static_cast<int>(0) );
    CheckFaceNeighbors(boundary2D_FROM);
    CheckFaceUnitNormalOrientation( boundary2D_FROM );
    CheckNodeFlags( boundary2D_FROM, IRREGULAR );
    CheckNodeParents( boundary2D_FROM );
    
    
    // 2D TESTS II
    // ===========

    // intersecting boundaries
    ANSYS_Model2D model2( "HorFracs2D", "CSMP-variables.txt" );
    CheckFaceUnitNormalOrientation( model2.Boundary("LEFT") );
    CheckFaceUnitNormalOrientation( model2.Boundary("RIGHT") );
    CheckFaceUnitNormalOrientation( model2.Boundary("TOP") );
    CheckFaceUnitNormalOrientation( model2.Boundary("BOTTOM") );
    Index eKey2( model2.Database().StorageKey( "element variable" ) );
    ElementNodes( model2.Region( "Model" ) );
    VTU_Interface<2> vtu2( model2 ); 

    pair<set<string>,bool>  result3 = model2.CreateInternalBoundaryFrom( "FRACTURE1" ); //, IRREGULAR, true );
    _test( result3.first.size() == 1 );
    Boundary<2>& fracture1( model2.Boundary(*result3.first.begin()) );
    CheckFaceNeighbors(fracture1);
    CheckFaceUnitNormalOrientation( fracture1 );
    CheckNodeFlags( fracture1, IRREGULAR );
    CheckNodeParents( fracture1 );
    pair<set<string>,bool>  result4 = model2.CreateInternalBoundaryFrom( "FRACTURE2" );
    Boundary<2>& fracture2( model2.Boundary(*result4.first.begin()) );
    CheckFaceNeighbors(fracture2);
    CheckFaceUnitNormalOrientation( fracture2 );
    CheckNodeFlags( fracture2, IRREGULAR );
    CheckNodeParents( fracture2 );

    //model2.InsertBoundary( "FRACTURE3", IRREGULAR, true );
    //Boundary<2>& fracture3( model2.Boundary("FRACTURE3") );
    //CheckFaceNeighbors(fracture3);
    //CheckFaceUnitNormalOrientation( fracture3 );
    //CheckNodeFlags( fracture3, IRREGULAR );
    //CheckNodeParents( fracture3 );
    //model2.InsertBoundary( "FRACTURE4", IRREGULAR, true );
    //Boundary<2>& fracture4( model2.Boundary("FRACTURE4") );
    //CheckFaceNeighbors(fracture4);
    //CheckFaceUnitNormalOrientation( fracture4 );
    //CheckNodeFlags( fracture4, IRREGULAR );
    //CheckNodeParents( fracture4 );    
    //model2.InsertBoundary( "FRACTURE5", IRREGULAR, true );
    //Boundary<2>& fracture5( model2.Boundary("FRACTURE5") );
    //CheckFaceNeighbors(fracture5);
    //CheckFaceUnitNormalOrientation( fracture5 );
    //CheckNodeFlags( fracture5, IRREGULAR );
    //CheckNodeParents( fracture5 );
    //model2.InsertBoundary( "FRACTURE6", IRREGULAR, true );
    //Boundary<2>& fracture6( model2.Boundary("FRACTURE6") );
    //CheckFaceNeighbors(fracture6);
    //CheckFaceUnitNormalOrientation( fracture6 );
    //CheckNodeFlags( fracture6, IRREGULAR );
    //CheckNodeParents( fracture6 );
    
    
    // 3D TESTS I
    // ==========

    ANSYS_Model3D model3D( "BoxHalfs3D", "CSMP-variables.txt" );
    Region<DIM3>& mref3D( model3D.Region( "Model" ) );
    mref3D.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 1.0 ) );
    mref3D.InputPropertyValue( "element variable", makeScalar( PLAIN, 2.0 ) );
    VTU_Interface<3> vtu3D( model3D ); vtu3D.OmitZeroInFileName( true );
    vtu3D.OutputDataToVTU( "TestNodeVariable3D", "nodal variable", "Model", static_cast<int>(0) );
    // RIGHT
    TestBoxBoundary<3>( model3D, "RIGHT", vtu3D );
    // TOP
    TestBoxBoundary<3>( model3D, "TOP", vtu3D );
    // BOTTOM
    TestBoxBoundary<3>( model3D, "BOTTOM", vtu3D );
    // LEFT
    TestBoxBoundary<3>( model3D, "LEFT", vtu3D );
    // FRONT
    TestBoxBoundary<3>( model3D, "FRONT", vtu3D );
    // BACK
    TestBoxBoundary<3>( model3D, "BACK", vtu3D );


    // 3D TESTS II
    // ===========
    
    ANSYS_Model3D model( "BoxHalfs", "CSMP-variables.txt", true );
    VTU_Interface<DIM3> vtu( model ); vtu.OmitZeroInFileName( true );
    Region<DIM3>& mref( model.Region( "Model" ) );
    mref.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 1.0 ) );
    mref.InputPropertyValue( "element variable", makeScalar( PLAIN, 2.0 ) );
    string region1Name( "MATRIX_LEFT" ), region2Name( "MATRIX_RIGHT" );

    // VARIABLE PLACED ON FACE
    Index faceKey( model.Database().StorageKey( "face variable" ) );
    if ( verbose_ ) {
        cout << "\n\nFace variable:\n";
        faceKey.Out();
      }

    // CREATE BETWEEN
    // --------------
    if ( verbose_ ) cout << "\nAttempting to insert csmp:: Boundary for region1  "<< region1Name << " and region2 " << region2Name << " ...\n";
    auto betweenBoundary = model.CreateBoundaryBetween( region1Name.data(), region2Name.data() );
    std::string boundary12Name("BOUNDARY_");
    boundary12Name += region1Name;
    boundary12Name += "_";
    boundary12Name += region2Name;
    bool boundary12Test(  model.ContainsBoundary( boundary12Name ) );
    if( boundary12Test && verbose_ )
      cout << "\nBoundary set up successful." << endl;
    else {
      if ( verbose_ ) cout << "\nBoundary set up MATRIX_LEFT-MATRIX_RIGHT failed." << endl;
      }
    _test( boundary12Test );
    // testing proper parent assignment
    Boundary<DIM3>& boundaryOne( model.Boundary( boundary12Name ) );
    CheckFaceNeighbors( boundaryOne );
    CheckFaceUnitNormalOrientation( boundaryOne );
    CheckNodeFlags( boundaryOne, IRREGULAR );
    CheckNodeParents( boundaryOne );
    model.Region( "Model" ).UpdateMemberIndexes();
    vector<size_t> outerParentsIDs;
    for( vector<Face<3>*>::const_iterator it = boundaryOne.CellsBegin(); it != boundaryOne.CellsEnd(); ++it )
      outerParentsIDs.push_back( (*it)->Parent(OUTSIDE)->Idx() );
    model.FormRegionFrom( "BOUNDARY OUTER PARENTS", outerParentsIDs );
    vector<size_t> innerParentsIDs;
    for( vector<Face<3>*>::const_iterator it = boundaryOne.CellsBegin(); it != boundaryOne.CellsEnd(); ++it )
      innerParentsIDs.push_back( (*it)->Parent(INSIDE)->Idx() );
    model.FormRegionFrom( "BOUNDARY INNER PARENTS", innerParentsIDs );
    if ( verbose_ ) {
         vtu.OutputDataToVTU( "BoundaryOuterParents", "element variable", "BOUNDARY OUTER PARENTS", static_cast<int>(0) );
         vtu.OutputDataToVTU( "BoundaryInnerParents", "element variable", "BOUNDARY INNER PARENTS", static_cast<int>(0) );
      }
    innerOuterParents( model, vtu, boundaryOne, "Parents1" );
    // BOUNDARY FACE COUNT
    Boundary<DIM3>& boundary12( model.Boundary( boundary12Name ) );
    size_t boundary12FaceCount( boundary12.Cells() );
    if ( verbose_ ) cout << "\nBoundary element count: " << boundary12FaceCount << endl;
    _test( boundary12FaceCount != 0 );
    CheckFaceNeighbors( boundary12 );
    CheckFaceUnitNormalOrientation( boundary12 );
    CheckNodeFlags( boundary12, IRREGULAR );
    CheckNodeParents( boundary12 );
    // BOUNDARY AREA
    double boundary12Area( boundary12.Area() );
    bool boundary12AreaNotZero( !withinTolerance( 0., boundary12Area, 0.1 ) );
    if ( verbose_ ) cout << "\nBoundary area: " << boundary12Area << endl;
    _test( boundary12AreaNotZero );
    if ( verbose_ ) vtu.OutputDataToVTU( "ElementVariable", "element variable", region1Name.data(), static_cast<int>(0) );
    model.Region( region2Name.data() ).InputPropertyValue( "element variable", makeScalar( PLAIN, 2.5 ) );
    if ( verbose_ )  {
         vtu.OutputDataToVTU( "ElementVariable", "element variable", region2Name.data(), static_cast<int>(0) );
         vtu.OutputDataToVTU( "NodalVariable", "nodal variable", boundary12, static_cast<int>(0) );
      }
    boundary12.InputPropertyValue( "face variable", makeScalar( PLAIN, 9999.0 ) );
    if ( verbose_ ) {
         vtu.OutputDataToVTU( "FaceVariable_A", "face variable", boundary12, static_cast<int>(0) );
         cout << "\nHULL_LEFT area: " << model.Region( "HULL_LEFT" ).Volume() << endl;
      }
    _test( InputElementAreaAsVolumeVariable<DIM3>( model, boundary12, "face variable" ) > 0 );
    if ( verbose_ ) vtu.OutputDataToVTU( "FaceVariable_B", "face variable", boundary12, static_cast<int>(0) );
    // FACE AREA
    const std::vector<Face<3>*>::const_iterator boundaryElementsEnd( boundary12.CellsEnd() );
    for( std::vector<Face<3>*>::const_iterator it = boundary12.CellsBegin(); it != boundaryElementsEnd; ++it )
    {
      Face<DIM3> currentFace = (*(*it));
      _test( !withinTolerance( currentFace.Area(), 0., 1.0E-5 ) );
    }


    // CREATE AROUND
    // -------------
    // removing pre-existing boundary so that new one can be created in its place
    Boundary<DIM3>& bdry = model.Boundary( betweenBoundary.first );
    model.Mesh().DeleteAndRepairConnnectivity( bdry.CellVector().begin(), bdry.CellVector().end() ); // getting rid of the faces
    model.RemoveBoundary( (betweenBoundary.first).c_str(), true ); // so that a new boundary can be created
    // creating and checking the new boundary
    _test( model.CreateBoundaryAround( "MATRIX_LEFT" ) );
    Boundary<DIM3>& boundaryHullLeft = model.Boundary( "BOUNDARY_MATRIX_LEFT_HULL" );
    _test( InputElementAreaAsVolumeVariable<DIM3>( model, boundaryHullLeft, "face variable" ) > 0 );
    vtu.OutputDataToVTU( "AddFaces", "face variable", boundaryHullLeft, static_cast<int>(0) );
    CheckFaceNeighbors(boundaryHullLeft);
    CheckFaceUnitNormalOrientation( boundaryHullLeft );
    CheckNodeFlags( boundaryHullLeft, IRREGULAR );
    CheckNodeParents( boundaryHullLeft );
    
  } // runLegacy




/**
       refactored tests from Philip Lang
*/
// TODO: move to ANSYS_Model3D_Test
void Boundary_Test::runCurrent()
  {
      size_t n_nodes_model{0U};
      
      // test that model contains no surface elements after boundary construction
      {
         ANSYS_Model3D m0( "BoxHalfs3D", "BoxHalfs3DirregularNoHalf", "CSMP-variables.txt", true, true );
         NoSurfaceElementsAsNodeParents( m0.Region("Model") );
         n_nodes_model = m0.Mesh().Nodes();
      }
      
      // test standard model with box boundaries
      {
        ANSYS_Model3D m01( "BoxHalfs3D", "CSMP-variables.txt", false, true );
        _test( m01.Boundaries() == 6 ); // six surfaces of the box-shaped model
        _test( m01.Regions() == 4 );
        _test( n_nodes_model == m01.Region("Model").Nodes() );
        if ( verbose_ ) {
            cout << "m03.Boundaries: " << m01.Boundaries() << "\n";
            cout << "m03.Regions: " << m01.Regions() << "\n";
            cout << "nodeCount: " << m01.Mesh().Nodes() << "\n";
            cout << "m03.Region(Model).Nodes(): " << m01.Region("Model").Nodes() << "\n";
          }
        // creating boundaries from the edges
        m01.EstablishEdgeBoundariesOfBoxShapedModel();
        _test( m01.Boundaries() == 18 );
      }
 
      // testing case where no boundaries are created because there are no boundary-regions listed in regions file
      {
        const bool irregular_mesh(true);   /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
        const bool binary_file(true);      /* true = binary, false = ascii */
        ANSYS_Model3D m02( "BoxHalfs3D", "BoxHalfs3DirregularNoBoundaries", "CSMP-variables.txt",
                            irregular_mesh, binary_file );
                            
        _test( m02.Boundaries() == 0 );
        _test( n_nodes_model == m02.Region("Model").Nodes() );

        size_t nullNeighborCount(0);
        Region<3>& rref( m02.Region("Model") );
        for ( auto eit = rref.CellsBegin(); eit != rref.CellsEnd(); ++eit )
          {
            for ( auto i{0U}; i < (*eit)->Neighbors(); ++i )
              if ( !(*eit)->Neighbor(i) )
                ++nullNeighborCount;
          }
        _test( nullNeighborCount != 0 );

        if ( verbose_ ) cout << "\nNull neighbor count: " << nullNeighborCount << endl;
      }
      
      // irregular model again
      {
        ANSYS_Model3D m04( "BoxHalfs3D", "BoxHalfs3Dirregular", "CSMP-1phase-variables.txt", true );
        cout << "m04.Boundaries(): " << m04.Boundaries() << "\n";
        cout << "m04.Regions(): " << m04.Regions() << "\n";
        cout << "nodeCount: " << m04.Mesh().Nodes() << "\n";
        cout << "m04.Region(Model).Nodes(): " << m04.Region("Model").Nodes() << "\n";
        _test( m04.Boundaries() == 6 );
        _test( m04.Regions() == 4 );
        _test( n_nodes_model == m04.Region("Model").Nodes() );
      }

  } // end runCurrent

} // csmp
