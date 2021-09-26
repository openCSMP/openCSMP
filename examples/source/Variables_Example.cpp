#include "Variables_Example.h"

#include "ANSYS_Model3D.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "VTU_Interface.h"
#include "Boundary.h"
#include "SplitBoundary.h"

using namespace std;

namespace csmp {

void Variables_Example::Specifications()
  {
    SetTitle( "Variable Operations" );
    SetDifficulty( 1 );
    SetCategory( "Software Functionality" );
    AddAuthor( "P. Lang" );
    AddDescription( "source in: Variables_Example.cpp" );
    AddDescription( "basic operations with csmp variables" );
    AddRequirement( "input model: 'FracBox'");
    AddRequirement( "VariablesTutorial.txt");
  }


void Variables_Example::Run()
  {
    // creating a model from an ANSYS mesh and inserting a  split boundary
    const size_t D(3);
    ANSYS_Model3D model( "FracBox", "VariablesTutorial.txt", true ); 
    
    const bool remove_dim_minus1_region(true);
    pair<set<string>,bool> boundary_patches = model.CreateInternalBoundaryFrom( "FRACTURE", remove_dim_minus1_region );  
    assert( boundary_patches.second == true );
    assert( boundary_patches.first.size() == 1 );
    const string boundary_name = (*boundary_patches.first.begin());
    Boundary<D>& fractureBoundary = model.Boundary( boundary_name.c_str() ); 
    pair<string,bool> split_boundary = model.CreateSplitBoundaryFrom( fractureBoundary );
    assert( split_boundary.second == true );
    
    // having a look at which regions and boundaries we have at the moment
    model.RegionsOut();
    model.BoundariesOut();
    model.SplitBoundariesOut();


    // getting keys to variables from property database (inside of model)
    const csmp::Index elementScalarKey = model.Database().StorageKey("element scalar");
    const csmp::Index elementIpVectorKey = model.Database().StorageKey("element ip vector");
    const csmp::Index elementFaipScalarKey = model.Database().StorageKey("element faip scalar");
    const csmp::Index elementSeipArrayKey = model.Database().StorageKey("element seip array");
    const csmp::Index faceScalarKey = model.Database().StorageKey("face scalar");
    const csmp::Index interfaceScalarKey = model.Database().StorageKey("interface scalar");
    const csmp::Index nodalTensorKey = model.Database().StorageKey("nodal tensor");
    const csmp::Index nodalArrayKey = model.Database().StorageKey("nodal array");
    const csmp::Index regionScalarKey = model.Database().StorageKey("region scalar");
    const csmp::Index boundaryVectorKey = model.Database().StorageKey("boundary vector");
    const csmp::Index splitBoundaryVectorKey = model.Database().StorageKey( "split boundary vector" );
    const csmp::Index modelArrayKey = model.Database().StorageKey("model array");


    // working variables
    ScalarVariable scalarVariablePlain;
    ScalarVariable scalarVariable( PLAIN, 0. );
    VectorVariable<D> vectorVariablePlain;
    VectorVariable<D> vectorVariable( PLAIN, PLAIN, DIRICH, 0., 1., 2. );
    TensorVariable<D> tensorVariablePlain;
    TensorVariable<D> tensorVariable( PLAIN, PLAIN, DIRICH, 0., 1., 2., 3., 4., 5., 6., 7., 8. ); // new: as many flags as dim
    ArrayVariable arrayVariablePlain;
    ArrayVariable arrayVariable( 100, 2., ANY );
    ArrayVariable nodalArrayVariable( "nodal array", model.Database() );
    ArrayVariable modelArrayVariable( "model array", model.Database() );


    // a single array variable stored on model
    model.Store( modelArrayKey, modelArrayVariable );
    arrayVariablePlain.Resize( modelArrayVariable.Size() );
    model.Read( modelArrayKey, arrayVariablePlain );
    
    for( size_t i(0); i < arrayVariablePlain.Size(); ++i )
      arrayVariablePlain(i) = double(i)*1.2;

    model.Store( modelArrayKey, arrayVariablePlain );


    // model subdomains (regions, boundaries, splitboundaries...)
    Region<D>& fracture = model.Region("FRACTURE");
    fracture.Store( regionScalarKey, scalarVariable );
    fracture.Read( regionScalarKey, scalarVariablePlain );
    cout << "\nRegion scalar: " << scalarVariablePlain << endl;

    Boundary<D>& boundary = model.Boundary("BOUNDARY1");
    boundary.Store( boundaryVectorKey, vectorVariable );
    boundary.Read( boundaryVectorKey, vectorVariablePlain );
    cout << "\nBoundary vector: " << vectorVariablePlain << endl;

    SplitBoundary<D>& splitboundary = model.SplitBoundary( split_boundary.first.c_str() );
    splitboundary.Store( splitBoundaryVectorKey, vectorVariable );
    splitboundary.Read( splitBoundaryVectorKey, vectorVariablePlain );
    cout << "\nSplitBoundary scalar: " << vectorVariablePlain << endl;

    // nodal operation
    const vector<Node<D>*>::const_iterator modelNodesEnd( model.Region("Model").NodesEnd() );
    for( vector<Node<D>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != modelNodesEnd; ++it )
      (*it)->Store( nodalTensorKey, tensorVariable );

    // is equivalent to
    model.InputPropertyValue( "nodal tensor", tensorVariable );

    // finite-element integration (Gauss quadrature) points
    const vector<Element<D>*>::const_iterator modelElementsEnd( model.Region("Model").ElementsEnd() );
    for( vector<Element<D>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != modelElementsEnd; ++it )
      {
        (*it)->Store( elementScalarKey, scalarVariable );

        for( size_t ip(0); ip < (*it)->IntegrationPoints(); ++ip )
          (*it)->Store( ip, elementIpVectorKey, vectorVariable );
      }

    // instatiate a finite volume scheme to establish fv ip variables
    NodeCenteredFiniteVolumeTransport<3> fvModule( "Model", model, "element scalar", "nodal scalar", "element vector", "nodal scalar", false, false );

    // finite volume integration points
    ArrayVariable seipArray( "element seip array", model.Database(), 0.2, ROBIN );
    for( vector<Element<D>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != modelElementsEnd; ++it )
      {
      // element
      (*it)->Store( elementScalarKey, scalarVariable );

      // element ip
      for( size_t ip(0); ip < (*it)->IntegrationPoints(); ++ip )
        (*it)->Store( ip, elementIpVectorKey, vectorVariable );

      // fv facets ip
      for( size_t fc(0); fc < (*it)->Facets(); ++fc )
        for( size_t fcip(0); fcip < (*it)->IntegrationPointsPerFacet(); ++fcip )
          (*it)->Store( fc, fcip, elementFaipScalarKey, scalarVariable );

      // fv sector ip
      for( size_t se(0); se < (*it)->Sectors(); ++se )
        for( size_t seip(0); seip < (*it)->IntegrationPointsPerSector(); ++seip )
          (*it)->Store( se, seip, elementSeipArrayKey, seipArray );
      }

    // adding a property at runtime
    model.CreateProperty( "new element scalar", "X", SCALAR, ELEMENT );
    model.InputPropertyValue( "new element scalar",  makeScalar( PLAIN, 1. ) );
    fracture.InputPropertyValue( "new element scalar",  makeScalar( PLAIN, 2. ) );

    for( Model<D>::boundaryConstIterator bit( model.BoundariesBegin() ); bit != model.BoundariesEnd(); ++bit )
      for( vector<Face<D>*>::const_iterator fit( bit->second.ElementsBegin() ); fit != bit->second.ElementsEnd(); ++fit )
        (*fit)->Store( faceScalarKey, makeScalar( PLAIN, (*fit)->Volume() ) );


    // XML (ASCII) output for the Visualisation Toolkit (VTK) / Paraview
    VTU_Interface<D> vtu(model);
    list<string> outputProps;
    outputProps.push_back("element scalar");
    outputProps.push_back("new element scalar");
    vtu.OutputDataToVTU( "VTU1", outputProps, "Model", static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", outputProps, fracture, static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", "face scalar", model.Boundary("BOUNDARY1"), static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", "face scalar", model.Boundary("BOUNDARY2"), static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", "face scalar", model.Boundary("BOUNDARY3"), static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", "face scalar", model.Boundary("BOUNDARY4"), static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", "face scalar", model.Boundary("BOUNDARY5"), static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU1", "face scalar", model.Boundary("BOUNDARY6"), static_cast<int>(0) );


    // displace nodes of splitboudnaries
    SplitBoundary<D>& fractureSplitBoundary = model.SplitBoundary( "SPLITBOUNDARY_FRACTURE" );
    const vector<InterFace<D>*>::const_iterator elementsEnd( fractureSplitBoundary.ElementsEnd() );
    for( vector<InterFace<D>*>::const_iterator it( fractureSplitBoundary.ElementsBegin() ); it != elementsEnd; ++it )
      {
        ScalarVariable volume( PLAIN, (*it)->InnerParent()->Volume() );
        (*it)->Store( interfaceScalarKey, volume );
        (*it)->N(0)->z( (*it)->N(0)->z() + 0.3 );
      }


    // output
    vtu.DeleteConnectivity();
    vtu.OutputDataToVTU( "VTU2", outputProps, "Model", static_cast<int>(0) );
    vtu.OutputDataToVTU( "VTU2", "interface scalar", fractureSplitBoundary, static_cast<int>(0) );

    cout << "\nDone...\n";
 } // Run()

} // csmp
