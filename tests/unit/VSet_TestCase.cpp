#include "VSet_TestCase.h"
#include "Boundary.h"
#include "Region.h"
#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "CSMP_highLevelUtilities.h"
#include "ScalarVariable.h"
#include "ArrayVariable.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "vsetMakers.h"

// File I/O and Initialization
#include "ANSYS_Model3D.h"

using namespace std;

namespace csmp {
VSet_TestCase::VSet_TestCase( bool verbose )
 : verbose_(verbose)
{
    this->setName("VSet_TestCase");
}


VSet_TestCase::~VSet_TestCase()
{
}


/**
      Testing whether processing of the VSet produces valid model with correct line-element normal
      orientations etc.
*/
void VSet_TestCase::run()
{
   _test( Test_EstablishElementConnectivity2D() );
   _test( Test_ModelConstructionAndSaving2D() );
   Test_ANSYS_ModelConstructionAndSaving2D( "HorFracs2D" );
    
} // end VSet_TestCase



bool VSet_TestCase::Test_ModelConstructionAndSaving2D()
  {
    enum{DIM=2U};
    if ( verbose_ ) cout <<"\nStart  of - "<<this->getName()<<endl<<endl;
    
    VSet<DIM> vset, vset2;
    test_Create_MeshPatchWithLineElements_VSet( vset );
    
    // creating a matching model topology
    ModelTopology mesh_topology( "VSet_TestCase", true );
    // all surface elements are "MATRIX"
    mesh_topology.AddRegion( "MATRIX", set<string>{"ISOPARAMETRIC_LINEAR_TRIANGLE", "ISOPARAMETRIC_LINEAR_QUADRILATERAL"},
                              vector<size_t>{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24} );

    // fracture line-element regions
    mesh_topology.AddRegion( "FRAC1", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{25,29,31} );
    mesh_topology.AddRegion( "FRAC2", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{26,27,28} );
    mesh_topology.AddRegion( "FRAC3", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{30} );
    mesh_topology.AddRegion( "FRAC4", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{32} );
    // boundaries
    mesh_topology.AddRegion( "BOTTOM", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{33,34,35} );
    mesh_topology.AddRegion( "RIGHT",  set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{36,37,38} );
    mesh_topology.AddRegion( "TOP",    set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{39,40,41} );
    mesh_topology.AddRegion( "LEFT",   set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{42,43,44} );

    _test( mesh_topology.Elements() == vset.Elements() );

    // adding corresponding materials to VSet
    vector<int32> pmtrl(45,1); // matrix
    fill( next(pmtrl.begin(),25), next(pmtrl.begin(),31), 2 ); // fine because wrong values will be overwritten next
    fill( next(pmtrl.begin(),26), next(pmtrl.begin(),28), 3 );
    fill( next(pmtrl.begin(),33), next(pmtrl.begin(),35), 4 );
    fill( next(pmtrl.begin(),36), next(pmtrl.begin(),38), 5 );
    fill( next(pmtrl.begin(),39), next(pmtrl.begin(),41), 6 );
    fill( next(pmtrl.begin(),42), pmtrl.end(), 7 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    // adding node and element numbers for comparisons
    PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
    elmt_nums.Reserve( vset.Elements() );
    for ( size_t i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );
    // node numbers
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );


    // build model from mesh
    Model<DIM>  model( mesh_topology, vset, "Vset_TestCase-variables.txt" );
    printModelDimensions( model, true );
    _test( printRangeOfVariable( model, "element number" ) <= vset.Elements() );
    _test( printRangeOfVariable( model, "node number" ) <= vset.Vertices() );
    const bool   get_indices_from_stored_variables{true};
    const size_t zero_errors{0};
    _test( model.Mesh().CheckElementConnectivity() == zero_errors );
    model.OutputMeshTo( vset2, get_indices_from_stored_variables );
    _test( vset2 == vset );
    vset2.Out();
    
    // saving model to binary
    const string test_model_name( string(model.Name()) + "Vset_TestCase" );
    model.OutputToBinaryFile( test_model_name.c_str() );
    
    // bringing the model back (calling reconstructor)
    Model<DIM>  model2( test_model_name );
    printModelDimensions( model, true );
    model2.OutputMeshTo( vset2 );
    
    // comparing it to original VSet
    if ( vset2 == vset ) return true;
    return false;
    
  } // end Test_ANSYS_ModelConstructionAndSaving2D
 



void VSet_TestCase::Test_ANSYS_ModelConstructionAndSaving2D( const std::string& input_file_name )
  {
    enum{DIM=2U};
    if ( verbose_ ) cout <<"\nStart  of - "<<this->getName()<<endl<<endl;
    
    // read ANSYS model data and build model
    ANSYS_Interface mesh_interface(true); // true = isoparametric elements
    ModelTopology   mesh_topology(true);
    VSet<DIM>       vset;

    const bool binary_file( true );
    const bool irregular_mesh( false );
    mesh_interface.Read_ANSYS_Mesh( input_file_name.c_str(), vset, mesh_topology, binary_file, irregular_mesh );

    // keep all mesh regions from topology and vset
    mesh_topology.ReduceToRegions( input_file_name.c_str() );
    map<size_t,size_t>  old_and_new_elmtids;
    mesh_topology.CreateNewElementNumbers( old_and_new_elmtids );
    vset.ReduceTo( old_and_new_elmtids );
    old_and_new_elmtids.clear();
    
    // processing the (deliberately) inconsistent VSet
    const size_t rotated_elements = vset.RenumberElementsCounterClockwise2D();
    if ( rotated_elements == 0U )
      ErrorHandler::Instance().notice( WARNING, "VSet_TestCase::TestModelConstructionAndSaving2D",
                                                "non-diagnostic test: element orientations are already correct.");
      
    // computes connectivity between equidimensional elements, faces and interfaces and replaces existing connectivity with it
    vset.RemovePfverts();
    vset.EstablishElementConnectivity2D();

    // build model from mesh
    Model<DIM>  model( mesh_topology, vset, "Vset_TestCase-variables.txt" );
    printModelDimensions( model, true );
    
    // saving model to binary
    model.OutputToBinaryFile( string( string(model.Name()) + "Vset_TestCase" ).c_str() );
    
  } // end Test_ANSYS_ModelConstructionAndSaving2D
 
 
 
 
  
  
void VSet_TestCase::Test_ANSYS_ModelConstructionAndSaving3D( const std::string& input_file_name )
  {
    enum{DIM=3U};
  
    //end Parameters section
    //------------------------------------
    if ( verbose_ ) cout <<"\nStart  of - "<<this->getName()<<endl<<endl;

    //------------------------------------
    // Model Output Regions only test
    //------------------------------------
    
    ScalarVariable diff( ANY, 1. );
    VectorVariable<3> vv( DIRICH, 2. );
    VectorVariable<3> vvPlain;
    TensorVariable<3> tv( ANY, 3. );
    TensorVariable<3> tvPlain;
    
    // Testing model without boundaries, variable&topology tests
    if ( verbose_ ) cout <<"Building ModelOutput..."<<endl;
    ANSYS_Model3D modelOutput1( input_file_name.data(),(this->getName()+".txt").c_str(),true,true,true,false);
    ArrayVariable na( "nodal array", modelOutput1.Database(), 2., ROBIN );
    const size_t elementCount1( modelOutput1.Region("Model").Elements() );
    const size_t nodeCount1( modelOutput1.Region("Model").Nodes() );
    const size_t regionCount1( modelOutput1.Regions() );
    modelOutput1.Region("Model").InputPropertyValue( "diffusivity", diff );
    modelOutput1.Region("Model").InputPropertyValue( "nodal array", na );
    modelOutput1.OutputToBinaryFile("model1");
    if ( verbose_ ) {
        cout<<"The Output of Model without boundaries done..."<<endl;
        cout <<"Building ModelInput..."<<endl;
      }
    Model<3U> modelInput1( string("model1") );
    Index dKey1( modelOutput1.Database().StorageKey("diffusivity") );
    Index naKey1( modelOutput1.Database().StorageKey("nodal array") );
    _test( elementCount1 == modelInput1.Region("Model").Elements() );
    _test( nodeCount1 == modelInput1.Region("Model").Nodes() );
    _test( regionCount1 == modelInput1.Regions() );
    ScalarVariable scalVal;
    ArrayVariable aVal( "nodal array", modelInput1.Database() );
    for( vector<Element<3>*>::const_iterator it( modelInput1.Region("Model").ElementsBegin() ); it != modelInput1.Region("Model").ElementsEnd(); ++it )
      {
        (*it)->Read( dKey1, scalVal );
        _test( scalVal == diff );
      }
    for( vector<Node<3>*>::const_iterator it( modelInput1.Region("Model").NodesBegin() ); it != modelInput1.Region("Model").NodesEnd(); ++it )
      {
      (*it)->Read( naKey1, aVal );
      _test( aVal == na );
      }
    
    // Testing model without boundaries, variable&topology tests
    if ( verbose_ ) cout <<"Building ModelOutput..."<<endl;
    ANSYS_Model3D modelOutput2( input_file_name.data(),(this->getName()+".txt").c_str(),true,true,true,true);
    Index boundaryScalarKey = modelOutput2.Database().StorageKey("boundary scalar");
    Index boundaryArrayKey = modelOutput2.Database().StorageKey("boundary array");
    Index regionVectorKey = modelOutput2.Database().StorageKey("region vector");
    Index modelTensorKey = modelOutput2.Database().StorageKey("model tensor");
    ArrayVariable ba( "boundary array",  modelOutput2.Database() );
    ArrayVariable baPlain( "boundary array",  modelOutput2.Database() );
    ba = 99.;
    const size_t nodeCount2( modelOutput2.Region("Model").Nodes() );
    const size_t regionCount2( modelOutput2.Regions() );
    const size_t boundaryCount2( modelOutput2.Boundaries() );
    modelOutput2.Region("Model").InputPropertyValue( "diffusivity", diff );
    modelOutput2.Region("Model").Store( regionVectorKey, vv );
    modelOutput2.Region("Model").InputPropertyValue( "nodal array", na );
    modelOutput2.Store( modelTensorKey, tv );
    modelOutput2.Boundary("BOUNDARY1").InputPropertyValue("boundary scalar", makeScalar( PLAIN, 1. ) );
    modelOutput2.Boundary("BOUNDARY2").InputPropertyValue("boundary array", ba );
    modelOutput2.OutputToBinaryFile("model2");
    if ( verbose_ ) {
        cout<<"The Output of Model with boundaries done..."<<endl;
        cout <<"Building ModelInput with boundaries..."<<endl;
      }
    Model<3U> modelInput2( string("model2") );
    _test( modelInput2.Boundary("BOUNDARY1").Read(boundaryScalarKey) == 1. );
    modelInput2.Region("Model").Read( regionVectorKey, vvPlain );
    _test( vvPlain == vv );
    modelInput2.Read( modelTensorKey, tvPlain );
    modelInput2.Boundary("BOUNDARY2").Read( boundaryArrayKey, baPlain );
    _test( tvPlain == tv );
    _test( baPlain == ba );
    _test( boundaryCount2 == modelInput2.Boundaries() );
    _test( regionCount2 == modelInput2.Regions() );    
    _test( nodeCount2 == modelInput2.Region("Model").Nodes() );  
    

    // testing model with finite volume variables
    ANSYS_Model3D modelOutput3( input_file_name.data(),(this->getName()+".txt").c_str(),true,true,true,true);
    NodeCenteredFiniteVolumeTransport<3> fvModule1( "Model", modelOutput3, "diffusivity", "nodal variable", "element vector", "nodal variable", false, false );
    Index faipVectorKey( modelOutput3.Database().StorageKey("faip vector") );
    Index seipTensorKey( modelOutput3.Database().StorageKey("seip tensor") );
    
    modelOutput3.InputPropertyValue( "faip vector", vv );
    modelOutput3.InputPropertyValue( "seip tensor", tv );
    Element<3>* ePtr = *modelOutput3.Region("Model").ElementsBegin();
    for( size_t f(0); f < ePtr->Facets(); ++f )
      for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
        }
    modelOutput3.OutputToBinaryFile("VSet_TestCase_modelOutput3");
    
    Model<3> modelInput3( string("VSet_TestCase_modelOutput3") );
    
    ePtr = *modelInput3.Region("Model").ElementsBegin();
    size_t ctr(0);
    for( size_t f(0); f < ePtr->Facets(); ++f )
      for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
      {
        ePtr->Read( f, fip, faipVectorKey, vvPlain );
        _test( vvPlain == vv );
        ++ctr;
      }
    for( size_t s(0); s < ePtr->Sectors(); ++s )
      for( size_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
      {
        ePtr->Read( s, sip, seipTensorKey, tvPlain );
        _test( tvPlain == tv );
        ++ctr;
      }
      _test( ctr == 10 );

      NodeCenteredFiniteVolumeTransport<3> fvModule2( "Model", modelInput3, "diffusivity", "nodal variable", "element vector", "nodal variable", false, false );
      for( size_t f(0); f < ePtr->Facets(); ++f )
        for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
          ++ctr;
        }
    _test( ctr == 16 );
    for( size_t s(0); s < ePtr->Sectors(); ++s )
      for( size_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
      {
        ePtr->Read( s, sip, seipTensorKey, tvPlain );
        _test( tvPlain == tv );
        ++ctr;
      }
    
    if ( verbose_ ) cout <<"\n\n"<<this->getName()<<" FINISHED!!!"<<endl;

  } // end 




/**
     uses specific line element model from vset_makers
     SKM 1/10/2022
*/
bool VSet_TestCase::Test_EstablishElementConnectivity2D()
 {
    VSet<2U> vset;
    test_Create_MeshPatchWithLineElements_VSet( vset );
    // making a backup copy
    VSet<2U> backup_vset( vset );
    
    // removing the neighbor information and recreating it
    vset.RemovePfverts();
    vset.EstablishElementConnectivity2D();

    // comparison
    cout <<"\nVSet_TestCase::Test_EstablishElementConnectivity2D: errors if any: ";
    auto itb=backup_vset.PfvertsBegin();
    size_t elmt{0U}, vec_mismatches{0U};
    for ( auto it=vset.PfvertsBegin(); it!=vset.PfvertsEnd(); ++it, ++itb ) {
        for ( size_t i=0U; i<(*it).size(); ++i )
          if ( (*it) != (*itb) ) {
               cerr <<"\n\t"<< elmt <<":";
               for ( auto i : (*itb) ) cerr <<" "<< i;
               cerr <<" vs. ";
               for ( auto i : (*it) ) cerr <<" "<< i;
               vec_mismatches++;
            }
        elmt++;
      }
    if ( vec_mismatches == 0 ) cout <<"NONE\n";
    cout << endl;
    
    if ( vec_mismatches > 0 ) return false;
    return true;
    
 } // end Test_CreateConsistentLineElementOrientations2D




} // namespace csmp
