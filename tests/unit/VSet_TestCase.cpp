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
#include "VTK_Interface.h"

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
   _test( Test_EstablishElementConnectivity3D() );
   _test( Test_ModelConstructionAndSaving2D() );
   Test_ANSYS_ModelConstructionAndSaving2D( "HorFracs2D" );
   Test_ANSYS_ModelConstructionAndSaving3D( "FracBox" ); 
    
} // end VSet_TestCase



bool VSet_TestCase::Test_ModelConstructionAndSaving2D()
  {
    enum{DIM=2U};
    if ( verbose_ ) cout <<"\nStart  of - "<<this->getName()<<endl<<endl;
    
    VSet<DIM> vset, vset2;
    ModelTopology mesh_topology = test_Create_MeshPatchWithLineElements_VSet( vset );
    
    _test( mesh_topology.Cells() == vset.Elements() );

    // build model from mesh
    const bool vset_only_contains_elements{ true };
    Model<DIM>  model( mesh_topology, vset, "VSet_TestCase-variables.txt", vset_only_contains_elements );
    printModelDimensions( model, true );
    _test( printRangeOfVariable( model, "element number" ) <= vset.Elements() );
    _test( printRangeOfVariable( model, "node number" ) <= vset.Vertices() );
    const bool   get_indices_from_stored_variables{true};
    const auto zero_errors{0};
    _test( model.Mesh().CheckElementConnectivity() == zero_errors );
    model.OutputMeshTo( vset2, get_indices_from_stored_variables );
    _test( vset2 == vset );
    
    // checking the single element regions
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

    const bool binary_file( true ), recreate_bflags(true);
    mesh_interface.Read_ANSYS_Mesh( input_file_name.c_str(), vset, mesh_topology, binary_file, recreate_bflags );
    // for ( auto i{0}; i<vset.BFlags(); ++i )
    //  cout <<" "<< static_cast<int>(vset.BoundaryFlag(i) );
    //cout << endl;

    // keep all mesh regions from topology and vset
    // calls CheckTopology and re-numbers nodes counter-clockwise if necessary
    mesh_topology.ReduceToDomains( input_file_name.c_str() );
    map<size_t,size_t>  old_and_new_elmtids;
    mesh_topology.CreateNewCellNumbers( old_and_new_elmtids );
    vset.ReduceTo( old_and_new_elmtids );
    old_and_new_elmtids.clear();
    _test( mesh_topology.Cells() == vset.Elements() );
    
    // computes connectivity between equidimensional elements, faces and interfaces and replaces existing connectivity with it
    vset.RemovePfverts();
    vset.EstablishElementConnectivity2D();
    
    // testing whether connectivity of the boundary faces has been achieved
    // looping over element faces that have a neighbor, reporting those where both nodes are at the boundary
    auto dodgy_neighbors{0};
    for ( size_t eidx{0}; eidx < vset.Elements(); ++eidx ) {
        CSMP_FEM_TYPE etype = parseFiniteElementTypeEnum( vset.ElementType(eidx) );
        // faces=neighbors
        auto face{0};
        for ( auto j=vset.PfvertsBegin(eidx); j!=vset.PfvertsEnd(eidx); ++j, ++face )
          if ( isTriangularElement(etype) && (*j) >= 0 ) {
             // face 0
             if ( face == 0 && vset.BoundaryFlag(vset.Plist(eidx,1)) != NOT && vset.BoundaryFlag(vset.Plist(eidx,2)) != NOT ) {
                  cerr <<"\nelement "<< eidx <<": face "<< face << " is at boundary but has neighbor: "<< *j;
                  cerr <<", node flags: "<< parseBoundary(intToBOX_BOUNDARY(vset.BoundaryFlag(vset.Plist(eidx,1))));
                  cerr <<" "<<              parseBoundary(intToBOX_BOUNDARY(vset.BoundaryFlag(vset.Plist(eidx,2))));
                  dodgy_neighbors++;
               }
             if ( face == 1 && vset.BoundaryFlag(vset.Plist(eidx,2)) != NOT && vset.BoundaryFlag(vset.Plist(eidx,0)) != NOT ) {
                  cerr <<"\nelement "<< eidx <<": face "<< face << " is at boundary but has neighbor: "<< *j;
                  cerr <<", node flags: "<< parseBoundary(intToBOX_BOUNDARY(vset.BoundaryFlag(vset.Plist(eidx,2))));
                  cerr <<" "<<              parseBoundary(intToBOX_BOUNDARY(vset.BoundaryFlag(vset.Plist(eidx,0))));
                  dodgy_neighbors++;
               }
             if ( face == 2 && vset.BoundaryFlag(vset.Plist(eidx,0)) != NOT && vset.BoundaryFlag(vset.Plist(eidx,1)) != NOT ) {
                  cerr <<"\nelement "<< eidx <<": face "<< face << " is at boundary but has neighbor: "<< *j;
                  cerr <<", node flags: "<< parseBoundary(intToBOX_BOUNDARY(vset.BoundaryFlag(vset.Plist(eidx,0))));
                  cerr <<" "<<              parseBoundary(intToBOX_BOUNDARY(vset.BoundaryFlag(vset.Plist(eidx,1))));
                  dodgy_neighbors++;
               }
          }
      }
    _test( dodgy_neighbors == 0 );

    // build model from mesh
    const bool get_domain_info_from_regions_file{true};
    Model<DIM>  model( mesh_topology, vset, "VSet_TestCase-variables.txt", get_domain_info_from_regions_file );
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

    ScalarVariable    diff( ANY, 1. );
    VectorVariable<3> vv( DIRICH, 2. );
    VectorVariable<3> vvPlain;
    TensorVariable<3> tv( ANY, 3. );
    TensorVariable<3> tvPlain;
    
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
      csmp::Index diff_key = modelOutput3.CreateProperty( "diffusivity", "m2/s", SCALAR, ELEMENT );
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
      modelOutput3.InputPropertyValue( "seip tensor 2", tv );
      Element<3>* ePtr = (*modelOutput3.Region("Model").CellsBegin());
      for( auto f(0); f < ePtr->Facets(); ++f )
        for( auto fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
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
      for( auto f(0); f < ePtr->Facets(); ++f )
        for( auto fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
          ++ctr;
        }
      for( auto s(0); s < ePtr->Sectors(); ++s )
        for( auto sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
        {
          ePtr->Read( s, sip, seipTensorKey, tvPlain );
          _test( tvPlain == tv );
          ++ctr;
        }
        _test( ctr == 10 );

       NodeCenteredFiniteVolumeTransport<3> fvModule2( "Model", modelInput3,
                                                      "element variable", // porosity var
                                                      "diffusivity",      // diffusivity var
                                                      "nodal variable",   // advected var
                                                      "element vector",   // transport var
                                                      "nodal variable",   // source var
                                                      false, false );

        for( auto f{0U}; f < ePtr->Facets(); ++f )
          for( auto fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
            ++ctr;
          }
      _test( ctr == 16 );
      for( auto s{0U}; s < ePtr->Sectors(); ++s )
        for( auto sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
        {
          ePtr->Read( s, sip, seipTensorKey, tvPlain );
          _test( tvPlain == tv );
          ++ctr;
        }
      
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
        for ( auto i{0}; i<(*it).size(); ++i )
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
