#include "Boundary_Test.h"

#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"
#include "Boundary.h"
#include "PL_Utilities.h"
#include "VTU_Interface.h"
#include "ANSYS_Interface.h"
#include "variableOperations.h"


using namespace std;

namespace csmp{

template<size_t dim>
void innerOuterParents( Model<dim>& model, VTU_Interface<dim>& vtu, Boundary<dim>& boundary, string fileName )
{
  const ScalarVariable zero( PLAIN, 0. );
  const ScalarVariable one( PLAIN, 1. );
  const ScalarVariable two( PLAIN, 2. );
  Index elVarKey( model.Database().StorageKey("element variable") );
  model.InputPropertyValue( "element variable", zero);
  for( typename vector<Face<dim>*>::const_iterator it = boundary.ElementsBegin(); it != boundary.ElementsEnd(); ++it )
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
}


template<size_t dim>
void Boundary_Test::ElementNodes( const Region<dim>& region )
  {
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator element( region.ElementsBegin() ); element != elementsEnd; ++element )
    _test( (*element)->Nodes() > 1 );
  }

template<size_t dim>
void Boundary_Test::NoSurfaceElementsAsNodeParents( const Region<dim>& region )
{
  const typename std::vector<Node<dim>*>::const_iterator domainNodesEnd( region.NodesEnd() );
  for( typename std::vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != domainNodesEnd; ++it )
    for( size_t i(0); i < (*it)->Parents(); ++i )
      _test( !(*it)->Parent(i)->IsSurfaceElement() );
}


template <size_t dim>
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


template <size_t dim>
size_t Boundary_Test::InputElementAreaAsVolumeVariable( Model<dim>& model, Boundary<dim>& boundary, const char* variableName )
{
  Index areaKey( model.Database().StorageKey( variableName ) );
  assert( areaKey.place == FACE ); // this needs to be done for Face::Read bc Read is inherited from Element and can
                                   // either be Face, Element or InterFace
  ScalarVariable area( PLAIN, 0. );
  size_t surfaceElementCount( 0 );
  const typename std::vector<Face<dim>*>::iterator domainElementsEnd( boundary.ElementsEnd() );
  for( typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it )
  {
    area = (*it)->Area();
    (*it)->Store( areaKey, area );
    ++surfaceElementCount;
  }

  return surfaceElementCount;
  }


/// checks for that at least one neigbor is present
template <size_t dim>
void Boundary_Test::CheckFaceNeighbors( const Boundary<dim>& boundary )
  {

  const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd( boundary.ElementsEnd() );
  for( typename std::vector<Face<dim>*>::const_iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it )
    {
      size_t notNullNeighbors(0);
      const size_t neighbors( (*it)->Neighbors() );
      for( size_t i = 0; i < neighbors; ++i )
        {
          if( (*it)->Neighbor(i) )
            ++notNullNeighbors;
        }
      _test( notNullNeighbors > 0 );
    }
  }


/// checks for that face unit normal points toward outer parent element
template <size_t dim>
void Boundary_Test::CheckFaceUnitNormalOrientation( const Boundary<dim>& boundary )
  {
    VectorVariable<dim> unFace( PLAIN, 9999999. ), faceToInner( PLAIN, 9999999. );
    Point<dim> bcFace, bcInner;
    const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd( boundary.ElementsEnd() );
    for( typename std::vector<Face<dim>*>::const_iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it )
      {
        (*it)->UnitNormal( unFace );          
        bcFace = (*it)->BaryCenter();
        bcInner = (*it)->Parent(INSIDE)->BaryCenter();
        for( size_t d(0); d < dim; ++d )
            faceToInner(d) = bcInner[d] - bcFace[d];
        _test( dotProduct(faceToInner,unFace) < 0. );
      }
  }


/// checks for flags of the boundary nodes
template <size_t dim>
void Boundary_Test::CheckNodeFlags( const Boundary<dim>& boundary, BOX_BOUNDARY flag, bool interiorOnly )
  {

  const typename std::vector<Node<dim>*>::const_iterator domainNodesEnd( boundary.NodesEnd() );
  const typename std::vector<Node<dim>*>::const_iterator interiorDomainNodesEnd( boundary.PerimeterNodesBegin() );
  for( typename std::vector<Node<dim>*>::const_iterator it = boundary.NodesBegin(); it != interiorDomainNodesEnd; ++it )
    {
      if( interiorOnly && it < interiorDomainNodesEnd )
        _test( (*it)->AtBoundary() == flag );
      else if( interiorOnly && it >= interiorDomainNodesEnd )
        _test( (*it)->AtBoundary() != NOT );
      else
        _test( (*it)->AtBoundary() == flag );
    }
  }


/// checks for flags of the boundary nodes
template <size_t dim>
void Boundary_Test::CheckNodeParents( const Boundary<dim>& boundary )
  {
  const typename std::vector<Node<dim>*>::const_iterator domainNodesEnd( boundary.NodesEnd() );
  for( typename std::vector<Node<dim>*>::const_iterator it = boundary.NodesBegin(); it != domainNodesEnd; ++it )
      _test( (*it)->Parents() > 1 );
  }

void Boundary_Test::runLegacy()
  {
    const size_t SPACE = 3U;
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
    _test( model2D.InsertBoundary( "MATRIX_LEFT", "MATRIX_RIGHT" ) );
    std::string boundary2D_BETWEENname( std::string("MATRIX_LEFT") + std::string("_") + std::string("MATRIX_RIGHT") );
    _test( model2D.ContainsBoundary( boundary2D_BETWEENname ) );
    Boundary<2>& boundary2D_BETWEEN( model2D.Boundary( boundary2D_BETWEENname ) ); 
    _test( InputElementAreaAsVolumeVariable<2>( model2D, boundary2D_BETWEEN, "face variable" ) > 0 );
    vtu2D.OutputDataToVTU( "2DTestFaceVariable", "face variable", boundary2D_BETWEEN, static_cast<int>(0) );
    CheckFaceNeighbors(boundary2D_BETWEEN);
    CheckNodeFlags( boundary2D_BETWEEN, IRREGULAR );
    CheckFaceUnitNormalOrientation( boundary2D_BETWEEN );
    CheckNodeParents( boundary2D_BETWEEN );
    innerOuterParents( model2D, vtu2D, boundary2D_BETWEEN, "Parents3" );
    // FROM
    _test( model2D.InsertBoundary( "STANDARD" ) );
    std::string boundary2D_FROMname( std::string("STANDARD") );
    _test( model2D.ContainsBoundary( boundary2D_FROMname ) );
    Boundary<2>& boundary2D_FROM( model2D.Boundary( boundary2D_FROMname ) ); 
    _test( InputElementAreaAsVolumeVariable<2>( model2D, boundary2D_FROM, "face variable" ) > 0 );
    vtu2D.OutputDataToVTU( "TestFaceVariable2D", "face variable", boundary2D_FROM, static_cast<int>(0) );
    CheckFaceNeighbors(boundary2D_FROM);
    CheckFaceUnitNormalOrientation( boundary2D_BETWEEN );
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

    model2.InsertBoundary( "FRACTURE1", IRREGULAR, true );
    Boundary<2>& fracture1( model2.Boundary("FRACTURE1") );
    CheckFaceNeighbors(fracture1);
    CheckFaceUnitNormalOrientation( fracture1 );
    CheckNodeFlags( fracture1, IRREGULAR );
    CheckNodeParents( fracture1 );
    model2.InsertBoundary( "FRACTURE2", IRREGULAR, true );
    Boundary<2>& fracture2( model2.Boundary("FRACTURE2") );
    CheckFaceNeighbors(fracture2);
    CheckFaceUnitNormalOrientation( fracture2 );
    CheckNodeFlags( fracture2, IRREGULAR );
    CheckNodeParents( fracture2 );
    model2.InsertBoundary( "FRACTURE3", IRREGULAR, true );
    Boundary<2>& fracture3( model2.Boundary("FRACTURE3") );
    CheckFaceNeighbors(fracture3);
    CheckFaceUnitNormalOrientation( fracture3 );
    CheckNodeFlags( fracture3, IRREGULAR );
    CheckNodeParents( fracture3 );
    model2.InsertBoundary( "FRACTURE4", IRREGULAR, true );
    Boundary<2>& fracture4( model2.Boundary("FRACTURE4") );
    CheckFaceNeighbors(fracture4);
    CheckFaceUnitNormalOrientation( fracture4 );
    CheckNodeFlags( fracture4, IRREGULAR );
    CheckNodeParents( fracture4 );    
    model2.InsertBoundary( "FRACTURE5", IRREGULAR, true );
    Boundary<2>& fracture5( model2.Boundary("FRACTURE5") );
    CheckFaceNeighbors(fracture5);
    CheckFaceUnitNormalOrientation( fracture5 );
    CheckNodeFlags( fracture5, IRREGULAR );
    CheckNodeParents( fracture5 );
    model2.InsertBoundary( "FRACTURE6", IRREGULAR, true );
    Boundary<2>& fracture6( model2.Boundary("FRACTURE6") );
    CheckFaceNeighbors(fracture6);
    CheckFaceUnitNormalOrientation( fracture6 );
    CheckNodeFlags( fracture6, IRREGULAR );
    CheckNodeParents( fracture6 );
    
    // 3D TESTS I
    // ==========

    ANSYS_Model3D model3D( "BoxHalfs3D", "CSMP-variables.txt" );
    Region<SPACE>& mref3D( model3D.Region( "Model" ) );
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
    VTU_Interface<SPACE> vtu( model ); vtu.OmitZeroInFileName( true );
    Region<SPACE>& mref( model.Region( "Model" ) );
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
    if ( verbose_ ) cout << "\nAttempting to insert csmp:: Boundary for region1  "<< region1Name << " and region2 " << region2Name << " ...\n";
    model.InsertBoundary( region1Name.data(), region2Name.data() );
    std::string boundary12Name;
    boundary12Name = region1Name;
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
    Boundary<SPACE>& boundaryOne( model.Boundary( boundary12Name ) );
    CheckFaceNeighbors( boundaryOne );
    CheckFaceUnitNormalOrientation( boundaryOne );
    CheckNodeFlags( boundaryOne, IRREGULAR );
    CheckNodeParents( boundaryOne );
    model.Region( "Model" ).UpdateMemberIndexes();
    vector<size_t> outerParentsIDs;
    for( vector<Face<3>*>::const_iterator it = boundaryOne.ElementsBegin(); it != boundaryOne.ElementsEnd(); ++it )
      outerParentsIDs.push_back( (*it)->Parent(OUTSIDE)->Idx() );
    model.FormRegionFrom( "BOUNDARY OUTER PARENTS", outerParentsIDs );
    vector<size_t> innerParentsIDs;
    for( vector<Face<3>*>::const_iterator it = boundaryOne.ElementsBegin(); it != boundaryOne.ElementsEnd(); ++it )
      innerParentsIDs.push_back( (*it)->Parent(INSIDE)->Idx() );
    model.FormRegionFrom( "BOUNDARY INNER PARENTS", innerParentsIDs );
    if ( verbose_ ) {
         vtu.OutputDataToVTU( "BoundaryOuterParents", "element variable", "BOUNDARY OUTER PARENTS", static_cast<int>(0) );
         vtu.OutputDataToVTU( "BoundaryInnerParents", "element variable", "BOUNDARY INNER PARENTS", static_cast<int>(0) );
      }
    innerOuterParents( model, vtu, boundaryOne, "Parents1" );
    // BOUNDARY FACE COUNT
    Boundary<SPACE>& boundary12( model.Boundary( boundary12Name ) );
    size_t boundary12FaceCount( boundary12.Elements() );
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
    _test( InputElementAreaAsVolumeVariable<SPACE>( model, boundary12, "face variable" ) > 0 );
    if ( verbose_ ) vtu.OutputDataToVTU( "FaceVariable_B", "face variable", boundary12, static_cast<int>(0) );
    // FACE AREA
    const std::vector<Face<3>*>::iterator boundaryElementsEnd( boundary12.ElementsEnd() );
    for( std::vector<Face<3>*>::iterator it = boundary12.ElementsBegin(); it != boundaryElementsEnd; ++it )
    {
      Face<SPACE> currentFace = (*(*it));
      _test( !withinTolerance( currentFace.Area(), 0., 1.0E-5 ) );
    }
    // CREATE FROM
    _test( model.InsertBoundary( "HULL_RIGHT" ) );
    Boundary<SPACE>& boundaryHullRight( model.Boundary( std::string("HULL_RIGHT") ) );
    _test( InputElementAreaAsVolumeVariable<SPACE>( model, boundaryHullRight, "face variable" ) > 0 );
    vtu.OutputDataToVTU( "FaceVariable", "face variable", boundaryHullRight, static_cast<int>(0) );
    double boundaryHullRightArea( boundaryHullRight.Area() );
    bool boundaryHullRightAreaNotZero( !withinTolerance( 0., boundaryHullRightArea, 0.1 ) );
    cout << "\nBoundary area: " << boundaryHullRightArea << endl;
    _test( boundary12AreaNotZero );
    CheckFaceNeighbors(boundaryHullRight);
    CheckFaceUnitNormalOrientation( boundaryHullRight );
    CheckNodeFlags( boundaryHullRight, IRREGULAR );
    CheckNodeParents( boundaryHullRight );
    // CREATE AROUND
    _test( model.AddFaces( "MATRIX_LEFT" ) );
    Boundary<SPACE>& boundaryHullLeft( model.Boundary( std::string("MATRIX_LEFT") ) );
    _test( InputElementAreaAsVolumeVariable<SPACE>( model, boundaryHullLeft, "face variable" ) > 0 );
    vtu.OutputDataToVTU( "AddFaces", "face variable", boundaryHullLeft, static_cast<int>(0) );
    CheckFaceNeighbors(boundaryHullLeft);
    CheckFaceUnitNormalOrientation( boundaryHullLeft );
    CheckNodeFlags( boundaryHullLeft, IRREGULAR );
    CheckNodeParents( boundaryHullLeft );

    // 3D TESTS III
    // ============
    /*
    ANSYS_Model3D model3Dsurface( "CircleHalfs2D", "CSMP-variables.txt", true );
    Region<3>& mref3Dsurface( model3Dsurface.Region( "Model" ) );
    Region<3>& upper3Dsurface( model3Dsurface.Region( "UPPER" ) );
    VTU_Interface<3> vtu3Dsurface( model3Dsurface ); vtu3Dsurface.OmitZeroInFileName( true );
    // AUTO INSERTED BY ANSYS INTERFACE
    std::string boundary3Dsurface_AUTOname( std::pair<string,string>( "Face", "BOUNDARY" ) );
    _test( model3Dsurface.ContainsBoundary( boundary3Dsurface_AUTOname ) );
    Boundary<3>& boundary3Dsurface_AUTO( model3Dsurface.Boundary( boundary3Dsurface_AUTOname ) ); 
    _test( InputElementAreaAsVolumeVariable<3>( model3Dsurface, boundary3Dsurface_AUTO, "face variable" ) > 0 );
    vtu3Dsurface.OutputDataToVTU( "TestFaceVariable3Dsurface", "face variable", boundary3Dsurface_AUTO ); 
    CheckFaceNeighbors(boundary3Dsurface_AUTO);
    CheckFaceUnitNormalOrientation( boundary3Dsurface_AUTO );
    CheckNodeFlags( boundary3Dsurface_AUTO, IRREGULAR );
    CheckNodeParents( boundary3Dsurface_AUTO );
    */

    
  }

  void Boundary_Test::runCurrent()
    {      
      ANSYS_Model3D m00( "BoxHalfs3D", "BoxHalfs3DirregularNoBoundaries", "CSMP-variables.txt", true, true, false );
      m00.InsertBoundary("HALF");
      NoSurfaceElementsAsNodeParents( m00.Region("Model") );

      ANSYS_Model3D m0( "BoxHalfs3D", "BoxHalfs3DirregularNoHalf", "CSMP-variables.txt", true, true, true );
      NoSurfaceElementsAsNodeParents( m0.Region("Model") );
      
      // this are a redundant checks to make sure the ANSYS_Model no csmp::Boundary constructor works
      // for both legacy box and irregular models (legacy functionality)    
      ANSYS_Model3D m01( "BoxHalfs3D", "CSMP-variables.txt", false, true, true, false );
        _test( m01.Boundaries() == 0 );
      ANSYS_Model3D m02( "BoxHalfs3D", "BoxHalfs3Dirregular", "CSMP-variables.txt", true, true, false );
        _test( m02.Boundaries() == 0 );  
      const size_t nodeCount( m01.Region("Model").Nodes() );
      _test( nodeCount == m02.Region("Model").Nodes() );

      size_t nullNeighborCount(0);
      Region<3>& rref( m02.Region("Model") );
      for ( vector<Element<3>*>::const_iterator eit = rref.ElementsBegin(); eit != rref.ElementsEnd(); ++eit ) 
        {
          for ( size_t i(0); i < (*eit)->Neighbors(); ++i )
            if ( !(*eit)->Neighbor(i) )
              ++nullNeighborCount; 
        }
      _test( nullNeighborCount != 0 );

      if ( verbose_ ) cout << "\nNull neighbor count: " << nullNeighborCount << endl;

      // for the case of an irregular model, check CreateAround and Split
      m02.AddFaces("Model");
      _test( m02.Boundaries() == 1 );

      _test( nodeCount == m02.Region("Model").Nodes() );

      Boundary<3>& b0102( m02.Boundary( "Model" ) );

      VTU_Interface<3> v02(m02);
      b0102.InputPropertyValue( "face variable", makeScalar( PLAIN, 1. ) );
      v02.OutputDataToVTU( "ModelBoundary", "face variable", b0102, static_cast<int>(0) );

      Model<3>::boundaryIterator boundary( m02.Boundary(b0102) );
      deque<string> regionsToRemove;
      for( Model<3>::regionIterator it( m02.UniqueRegionsBegin() ); it != m02.UniqueRegionsEnd(); ++it )
        {
          if( m02.IsBoundaryName( it->first ) )
            {
              m02.DivideBoundary( boundary, it );
              regionsToRemove.push_back( it->first );
            }
          if( boundary->second.Elements() == 0 )
            m02.RemoveBoundary( boundary->second );
        }

      for( size_t i(0); i < regionsToRemove.size(); ++i ) {
          // SKM FIX m02.RemoveRegion( regionsToRemove.at(i).c_str() );
          m02.RemoveFromRegion( "Model", regionsToRemove.at(i).c_str() );
          m02.MoveToNonUniqueRegions( regionsToRemove.at(i).c_str() );
        }

      _test( nodeCount == m02.Region("Model").Nodes() );

      v02.DeleteConnectivity();

      for( Model<3>::boundaryIterator it( m02.BoundariesBegin() ); it != m02.BoundariesEnd(); ++it )
        {
          it->second.InputPropertyValue( "face variable", makeScalar( PLAIN, 1. ) );
          if ( verbose_ ) v02.OutputDataToVTU( "ModelBoundarySplit", "face variable", it->second, static_cast<int>(0) );
        }

      m02.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 1. ) );
      
      m02.Boundary("TOP").InputPropertyValue( "nodal variable", makeScalar( PLAIN, 2. ) );
      if ( verbose_ ) v02.OutputDataToVTU( "BoundaryValue", "nodal variable", "Model", static_cast<int>(0) );

      ANSYS_Model3D m03( "BoxHalfs3D", "CSMP-variables.txt" );
      _test( m03.Boundaries() == 6 );
      _test( m03.Regions() == 4 );
      _test( nodeCount == m03.Region("Model").Nodes() );

      ANSYS_Model3D m04( "BoxHalfs3D", "BoxHalfs3Dirregular", "CSMP-1phase-variables.txt", true );
      _test( m04.Boundaries() == 6 );  
      _test( m04.Regions() == 4 );
      _test( nodeCount == m04.Region("Model").Nodes() );

      return;
    }



} // csmp
