#include "VSet_TestCase.h"
#include "Boundary.h"
#include "Region.h"
#include "ModelTopology.h"
#include "ScalarVariable.h"
#include "ArrayVariable.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "vsetMakers.h"
#include "VTK_Interface.h"
#include "ANSYS_Interface.h"

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
   _test( Test_EstablishElementConnectivity3D() );
   _test( Test_ModelConstructionAndSaving2D() );
   Test_ModelConstructionAndSaving3D(); // uses FracBox and other models
    
} // end VSet_TestCase



bool VSet_TestCase::Test_ModelConstructionAndSaving2D()
  {
    enum{DIM=2U};
    if ( verbose_ ) cout <<"\nStart  of - "<<this->getName()<<endl<<endl;
    
    // building and testing first 2D model from mesh & topology (already including faces and interfaces)
    // -------------------------------------------------------------------------------------------------
    VSet<DIM> vset, vset2;
    ModelTopology mesh_topology = test_Create_BoundarySplitBoundaryPatch( vset );
    _test( mesh_topology.Cells() == vset.Cells() );
    {
      // adding the original element numbers to VSet, assigning the same numbers
      // as face numbers as these will be converted later
      PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
      elmt_nums.Reserve( vset.Elements() );
      for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
      vset.AddData( "element number", elmt_nums );
      // face numbers
      PropertyData face_nums( FACE, SCALAR, 2U );
      face_nums.Reserve( vset.Faces() );
      for ( size_t i{ vset.Elements() }; i<vset.Elements() + vset.Faces(); ++i ) pushBack( face_nums, makeScalar( ANY, i ) );
      vset.AddData( "face number", face_nums );
      // interface numbers
      PropertyData iface_nums( INTER_FACE, SCALAR, 2U );
      iface_nums.Reserve( vset.Interfaces() );
      for ( size_t i{ vset.Elements() + vset.Faces() }; i<vset.Cells(); ++i ) pushBack( iface_nums, makeScalar( ANY, i ) );
      vset.AddData( "interface number", iface_nums );
      // node numbers
      PropertyData node_nums( NODE, SCALAR, 2U );
      node_nums.Reserve( vset.Vertices() );
      for ( size_t i{0U}; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
      vset.AddData( "node number", node_nums );
    }
    // ----------------------------------------------------
    // creating and testing model with faces and interfaces
    // ----------------------------------------------------
    {
      const bool vset_only_contains_elements{ false };
      Model<DIM>  model( mesh_topology, vset, "VSet_TestCase-variables.txt", vset_only_contains_elements );
      printModelDimensions( model, true );
      _test( printRangeOfVariable( model, "element number" ) <= vset.Elements() );
      _test( printRangeOfVariable( model, "node number" ) <= vset.Vertices() );
      // correct reproduction of what is in the VSet, see  model.Mesh().Out();
      
      // copying "element number" to "face number" for the faces created from lower-dimensional elements
      const csmp::Index fn_key = model.Database().StorageKey("face number");
      _test( mesh_topology.CellsWithinDomain("BOTTOM") == model.Boundary("BOTTOM").Cells() );
      _test( mesh_topology.CellsWithinDomain("RIGHT")  == model.Boundary("RIGHT").Cells() );
      _test( mesh_topology.CellsWithinDomain("TOP")    == model.Boundary("TOP").Cells() );
      _test( mesh_topology.CellsWithinDomain("LEFT")   == model.Boundary("LEFT").Cells() );
      Boundary<2U>& bottom{ model.Boundary("BOTTOM") }, right{ model.Boundary("RIGHT") },
                    top{ model.Boundary("TOP") }, left{ model.Boundary("LEFT") };
      // BOTTOM
      size_t n_face{0U};
      const auto end1 = mesh_topology.CellsOfDomainEnd("BOTTOM");
      for ( auto it=mesh_topology.CellsOfDomainBegin("BOTTOM"); it!=end1; ++it )
        bottom.E( n_face++ )->Store( fn_key, makeScalar(FIELD_DATA,*it) );
      // RIGHT
      n_face = 0U;
      const auto end2 = mesh_topology.CellsOfDomainEnd("RIGHT");
      for ( auto it=mesh_topology.CellsOfDomainBegin("RIGHT"); it!=end2; ++it )
        right.E( n_face++ )->Store( fn_key, makeScalar(FIELD_DATA,*it) );
      // TOP
      n_face = 0U;
      const auto end3 = mesh_topology.CellsOfDomainEnd("TOP");
      for ( auto it=mesh_topology.CellsOfDomainBegin("TOP"); it!=end3; ++it )
        top.E( n_face++ )->Store( fn_key, makeScalar(FIELD_DATA,*it) );
      // LEFT
      n_face = 0U;
      const auto end4 = mesh_topology.CellsOfDomainEnd("LEFT");
      for ( auto it=mesh_topology.CellsOfDomainBegin("LEFT"); it!=end4; ++it )
        left.E( n_face++ )->Store( fn_key, makeScalar(FIELD_DATA,*it) );
      
      // testing whether original Face and Interface numbers are preserved in output
      const bool   get_indices_from_stored_variables{true};
      const auto zero_errors{0};
      _test( model.Mesh().CheckElementConnectivity() == zero_errors );
      model.OutputMeshTo( vset2, get_indices_from_stored_variables );
      _test( vset2 == vset );
    }
    
    // -------------------------------------------------------------------
    // Building second model with faces created from elements
    // -------------------------------------------------------------------
    {
      // model consists only of elements
      mesh_topology = test_Create_MeshPatchWithLineElements_VSet( vset );
      const bool vset_only_contains_elements{ true };
      Model<DIM>  model( mesh_topology, vset, "VSet_TestCase-variables.txt", vset_only_contains_elements );
      printModelDimensions( model, true );
      // checking single element regions
      //model.RegionsOut();
      const Region<2U>& fracs(model.Region("FRAC3"));
      //cerr <<"\ninterior vs perimeter: "<< fracs.InteriorCells() <<" "<< fracs.PerimeterCells();
      _test( fracs.InteriorCells()  == 0 );
      _test( fracs.PerimeterCells() == 1 );
      // saving model to binary
      const string test_model_name( string(model.Name()) + "Vset_TestCase" );
      model.OutputToBinaryFile( test_model_name.c_str() );
      // bringing the model back (calling reconstructor)
      Model<DIM>  model2( test_model_name );
      // getting this a distinct name for the destruction process
      model2.Name( (string(model2.Name()) + "_reconstructed").c_str() );
      printModelDimensions( model, true );
      model2.OutputMeshTo( vset2 );
    }
    // comparing it to original VSet
    if ( vset2 == vset ) return true;
    return false;
    
  } // end Test_ANSYS_ModelConstructionAndSaving2D
 



 
 
 
 
  
// various input models from vsetMakers.h
void VSet_TestCase::Test_ModelConstructionAndSaving3D()
  {
    enum{DIM=3U};
  
    //end Parameters section
    //------------------------------------
    if ( verbose_ ) cout <<"\nStart  of - "<<this->getName()<<endl<<endl;

    const ScalarVariable    diff( ANY, 1. );
    const VectorVariable<3> vv( DIRICH, 2. );
    VectorVariable<3>       vvPlain;
    const TensorVariable<3> tv( ANY, 3. );
    TensorVariable<3>       tvPlain;
    
    //------------------------------------
    // 1. Model Output Regions only test
    //------------------------------------
    {
    // Testing model without boundaries, variable&topology tests
    if ( verbose_ ) cout <<"Building ModelOutput (only regions)..."<<endl;
    VSet<3U> vset;
    testCreateTetra_VSet( vset );
    Model<3U>    modelOutput1( vset, "VSet_TestCase-variables.txt" );
    Region<3U>&  model_domain{ modelOutput1.Region("Model") };
    const size_t elementCount1( model_domain.Cells() );
    const size_t nodeCount1( model_domain.Nodes() );
    const size_t regionCount1( modelOutput1.Regions() );
    model_domain.InputPropertyValue( "element variable", diff );
    modelOutput1.OutputToBinaryFile("model1");
    if ( verbose_ ) {
        cout<<"The Output of Model without boundaries done..."<<endl;
        cout <<"Building ModelInput..."<<endl;
      }
    // bringing the model back from disk
    Model<3U> modelInput1( string("model1") );
    const Index dKey1( modelInput1.Database().StorageKey("element variable") );
    const Index naKey1( modelInput1.Database().StorageKey("nodal variable") );
    _test( elementCount1 == model_domain.Cells() );
    _test( nodeCount1 == model_domain.Nodes() );
    _test( regionCount1 == modelInput1.Regions() );
    ScalarVariable scalVal;
    ArrayVariable aVal( "nodal array", modelInput1.Database() );
    for( auto it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
      {
        (*it)->Read( dKey1, scalVal );
        _test( scalVal == diff );
      }
      
    } // Model Output Regions only test
    
    
    
    // ------------------------------------------------------------------------------------
    // 2. Testing model with boundaries, variable&topology tests (requires 'FracBox' model)
    // ------------------------------------------------------------------------------------
    {
        if ( verbose_ ) cout <<"Building ModelOutput (regions & irregular boundaries)..."<<endl;
        ModelTopology topology;
        VSet<3U>      vset;
        test_Create_FracBox( topology, vset );
        Model<3U>    modelOutput2( topology, vset, "VSet_TestCase-variables.txt", true );
        
        const Index boundaryScalarKey = modelOutput2.Database().StorageKey("boundary scalar");
        const Index boundaryArrayKey = modelOutput2.Database().StorageKey("boundary array");
        const Index regionVectorKey = modelOutput2.Database().StorageKey("region vector");
        const Index modelTensorKey = modelOutput2.Database().StorageKey("model tensor");
        
        ArrayVariable na( "nodal array", modelOutput2.Database() );
        ArrayVariable ba( "boundary array",  modelOutput2.Database() );
        ArrayVariable baPlain( "boundary array",  modelOutput2.Database() );
        ba = 99.;
        
        Region<3U>& model_domain2{ modelOutput2.Region("Model") };
        model_domain2.InputPropertyValue( "element variable", diff );
        model_domain2.Store( regionVectorKey, vv );
        model_domain2.InputPropertyValue( "nodal array", na );
        modelOutput2.Store( modelTensorKey, tv );
        modelOutput2.Boundary("BOUNDARY1").InputPropertyValue("boundary scalar", makeScalar( PLAIN, 1. ) );
        modelOutput2.Boundary("BOUNDARY2").InputPropertyValue("boundary array", ba );
        modelOutput2.OutputToBinaryFile("model2");
        if ( verbose_ ) {
            cout<<"The Output of Model with boundaries done..."<<endl;
            cout <<"Building ModelInput with boundaries..."<<endl;
          }
    }
    // bringing the model back from disk
    {
        Model<3U> modelInput2( string("model2") );
        Region<3U>& model_domain2{ modelInput2.Region("Model") };
        
        const Index boundaryScalarKey = modelInput2.Database().StorageKey("boundary scalar");
        const Index boundaryArrayKey = modelInput2.Database().StorageKey("boundary array");
        const Index regionVectorKey = modelInput2.Database().StorageKey("region vector");
        const Index modelTensorKey = modelInput2.Database().StorageKey("model tensor");

        _test( modelInput2.Boundary("BOUNDARY1").Read(boundaryScalarKey) == 1. );
        model_domain2.Read( regionVectorKey, vvPlain );
        _test( vvPlain == vv );
        modelInput2.Read( modelTensorKey, tvPlain );
        ArrayVariable baPlain( "boundary array",  modelInput2.Database() );
        modelInput2.Boundary("BOUNDARY2").Read( boundaryArrayKey, baPlain );
        _test( tvPlain == tv );
        _test( baPlain == baPlain );
    
    } // end test with irregular bundaries



    // ------------------------------------------------------------------------------------
    // 3. testing model creation with finite volume variables
    // ------------------------------------------------------------------------------------
    {
      VSet<3U> vset;
      testCreateTetra_VSet( vset );
      Model<3U> modelOutput3( vset, "VSet_TestCase-variables.txt" );
      
      // 'diffusity' for FV scheme
      csmp::Index diff_key = modelOutput3.CreateProperty( "diffusivity", "D", "m2/s", SCALAR, ELEMENT );
      _test( diff_key.place == ELEMENT );
      _test( diff_key.type  == SCALAR );
      
      NodeCenteredFiniteVolumeTransport<3> fvModule1( "Model", modelOutput3,
                                                      "element variable", // porosity var
                                                      "diffusivity",      // diffusivity var
                                                      "nodal variable",   // advected var
                                                      "element vector",   // transport var
                                                      "nodal variable",   // source var
                                                      false, false );
                                                      
      const Index faipVectorKey( modelOutput3.Database().StorageKey("faip vector 1") );
      const Index seipTensorKey( modelOutput3.Database().StorageKey("seip tensor 1") );
      
      modelOutput3.InputPropertyValue( "faip vector 1", vv );
      modelOutput3.InputPropertyValue( "seip tensor 1", tv );
      Element<3>* ePtr = (*modelOutput3.Region("Model").CellsBegin());
      for( auto f(0); f < ePtr->Facets(); ++f )
        for( auto fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
          }
      for( auto s(0); s < ePtr->Sectors(); ++s )
        for( auto sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
        {
          ePtr->Read( s, sip, seipTensorKey, tvPlain );
          _test( tvPlain == tv );
        }
      modelOutput3.OutputToBinaryFile("VSet_TestCase_modelOutput3");
    }
    
    // bringing the model back from disk
    // ---------------------------------
    {
      Model<3> modelInput3( string("VSet_TestCase_modelOutput3") );

      const Index faipVectorKey( modelInput3.Database().StorageKey("faip vector 1") );
      const Index seipTensorKey( modelInput3.Database().StorageKey("seip tensor 1") );

      Element<3>*  ePtr = *modelInput3.Region("Model").CellsBegin();
      auto ctr(0);
      for( auto f{0U}; f < ePtr->Facets(); ++f )
        for( auto fip{0U}; fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
          ++ctr;
        }
      for( auto s{0U}; s < ePtr->Sectors(); ++s )
        for( auto sip{0U}; sip < ePtr->IntegrationPointsPerSector(); ++sip )
        {
          ePtr->Read( s, sip, seipTensorKey, tvPlain );
          _test( tvPlain == tv );
          ++ctr;
        }
        _test( ctr == 10 );
      
    } // end test case with FV functionality
    
    if ( verbose_ )
      cout <<"\n\n"<<this->getName()<<" FINISHED!!!"<<endl;

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
    auto elmt{0U}, vec_mismatches{0U};
    for ( auto it=vset.PfvertsBegin(); it!=vset.PfvertsEnd(); ++it, ++itb ) {
        for ( uint32_t i{0u}; i<(*it).size(); ++i )
          if ( (*it) != (*itb) ) {
               cerr <<"\n\t"<< elmt <<":";
               for ( auto j : (*itb) ) cerr <<" "<< j;
               cerr <<" vs. ";
               for ( auto j : (*it) ) cerr <<" "<< j;
               vec_mismatches++;
            }
        elmt++;
      }
    if ( vec_mismatches == 0 ) cout <<"NONE\n";
    cout << endl;
    
    if ( vec_mismatches > 0 ) return false;
    return true;
    
 } // end Test_CreateConsistentLineElementOrientations2D







bool VSet_TestCase::Test_EstablishElementConnectivity3D()
 {
    VSet<3U> vset;
    test_Create_Prism_Hexa_VSet( vset);
    BoundaryFlagsToVTK( vset );

    // making a backup copy
    VSet<3U> backup_vset( vset );
    
    // removing the neighbor information and recreating it
    vset.RemovePfverts();
    vset.EstablishElementConnectivity3D();

    // comparison - but only checking for the existing element neighbors since atBoundary(elmt) is not used further
    cout <<"\nVSet_TestCase::Test_EstablishElementConnectivity3D: errors if any: ";
    auto itb=backup_vset.PfvertsBegin();
    auto elmt{0U}, vec_mismatches{0U};
    for ( auto it=vset.PfvertsBegin(); it!=vset.PfvertsEnd(); ++it, ++itb ) {
        for ( auto i{0}; i<(*it).size(); ++i )
          if ( (*it)[i] >= 0 && (*it)[i] != (*itb)[i] ) {
               cerr <<"\n\t"<< elmt <<":";
               for ( auto j : (*itb) ) cerr <<" "<< j;
               cerr <<" vs. ";
               for ( auto j : (*it) ) cerr <<" "<< j;
               vec_mismatches++;
            }
        elmt++;
      }
    if ( vec_mismatches == 0 ) cout <<"NONE\n";
    cout << endl;
    
    if ( vec_mismatches > 0 ) return false;
    return true;
    
 } // end Test_CreateConsistentLineElementOrientations3D





/**
   Creates model and writes boundary flags to 'nodal variable' variable
*/
void VSet_TestCase::BoundaryFlagsToVTK( VSet<3>& vset )
 {
    const string variable_file{"CSMP-variables.txt"};
    Model<3> model( vset, variable_file.c_str() );
    
    VTK_Interface<3>  vtk_out;
    
    boxFlagsToVariable( model, "nodal variable", "element number" );
    
    vtk_out.OutputDataToVTK( model, "BoundaryFlagsToVTK_output", "nodal variable", static_cast<int>(0) );
    
 } // BoundaryFlagsToVTK





} // namespace csmp
