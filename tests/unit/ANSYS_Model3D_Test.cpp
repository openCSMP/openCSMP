#include "ANSYS_Model3D_Test.h"
#include "ANSYS_Model3D.h"
#include "ModelSubDomain_Test.h" // for the use of CompareModelSubDomains
#include "Boundary.h"
#include "Region.h"
#include "VTU_Interface.h"
//#include "NodeCenteredFiniteVolumeTransport.h"
#include "Timer.hpp"

using namespace std;

namespace csmp {
  
void create_ANSYS3D_Model( bool contiguous, bool reconstruct_from_file )
 {
   Model<3>* model3d_(0);
   string    model3d_name_;
   
    // ansys 3d model - discontiguous
   if ( !contiguous ) {
        cout << "\n-------------------------------------------------------";
        cout << "\nMeshManager_Test: ANSYS model 'ModelDykeAllLayersSplit'";
        cout << "\n-------------------------------------------------------";
        string varFileName = "ANSYS_SplitBoundaryMatch_Test-variables.txt";
        model3d_name_ = "ModelDykeAllLayersSplit";
        model3d_ = new ANSYS_Model3D(model3d_name_.c_str(), varFileName.c_str(), true );
 
        //writing ansys model to file deleting it and then recreating a csmp native model from the file
        if ( reconstruct_from_file ) {
            model3d_->OutputToBinaryFile(model3d_name_.c_str());
            delete model3d_;
            model3d_ = new Model<3U>(model3d_name_);
           }
        delete model3d_;
        model3d_ = nullptr;
        return;
     }

    // ansys 3d model - contiguous
    cout << "\n-------------------------------------------------------";
    cout << "\nMeshManager_Test: ANSYS model 'prism_test'";
    cout << "\n-------------------------------------------------------";
    string varFileName = "CSMP-variables.txt";
    model3d_name_ = "prism_test";
    model3d_ = new ANSYS_Model3D(model3d_name_.c_str(), varFileName.c_str());

    if ( reconstruct_from_file ) {
        // writing ansys model to file deleting it and then recreating a csmp native model from the file
        model3d_->OutputToBinaryFile(model3d_name_.c_str());
        delete model3d_;
        model3d_ = new Model<3U>(model3d_name_);
      }
  
 } // end create_ANSYS3D_Model






  void ANSYS_Model3D_Test::run()
    {
      const bool verbose(false);
      // BINARY IO
      // ---------
      if ( verbose ) {
           cout <<"\nStart simulation of - "<<this->getName()<<endl<<endl;
           cout <<"Building ModelOutput..."<<endl;
        }
      const string variablesFile("UIVariables.txt");
      Timer timer;
      timer.Start();
      //                          fileset       regions-file
      ANSYS_Model3D modelOutput1( "BoxHalfs3D", "BoxHalfs3D", variablesFile.c_str(), true );
      const double icemModelTime( timer.Stop() );

      size_t nullNeighborsOut(0);
      const Region<3U>& model_domain1(modelOutput1.Region("Model"));
      for ( auto it = model_domain1.CellsBegin(); it != model_domain1.CellsEnd(); ++it )
        for ( uint32_t n{0}; n < (*it)->Neighbors(); ++n )
          if( (*it)->Neighbor(n) == nullptr )
            ++nullNeighborsOut;

      _test( nullNeighborsOut != 0 );

      const double matrixLeftValue(2.);
      const double matrixRightValue(3.);
      const double tolerance(1.0e-5);
      modelOutput1.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 1.0 ) );
      modelOutput1.InputPropertyValue( "permeability", makeScalar( PLAIN, 1.0 ) ); // to avoid NaN values
      modelOutput1.Region("MATRIX_LEFT").InputPropertyValue( "permeability", makeScalar( PLAIN, matrixLeftValue ) );
      modelOutput1.Region("MATRIX_RIGHT").InputPropertyValue( "permeability", makeScalar( PLAIN, matrixRightValue ) );

      const size_t elementCount( modelOutput1.Region("Model").Cells() );
      const size_t nodeCount( modelOutput1.Region("Model").Nodes() );
      const size_t regionCount( modelOutput1.Regions() );
      const size_t matrixLeftNodeCount(  modelOutput1.Region("MATRIX_LEFT").Nodes() );
      const size_t matrixRightNodeCount(  modelOutput1.Region("MATRIX_RIGHT").Nodes() );
      const size_t matrixLeftElementCount(  modelOutput1.Region("MATRIX_LEFT").Cells() );
      const size_t matrixRightElementCount(  modelOutput1.Region("MATRIX_RIGHT").Cells() );

      _test( modelOutput1.ContainsBoundary("TOP") );
      _test( modelOutput1.ContainsBoundary("BOTTOM") );
      _test( modelOutput1.ContainsBoundary("RIGHT") );
      _test( modelOutput1.ContainsBoundary("FRONT") );
      _test( modelOutput1.ContainsBoundary("BACK") );
      _test( modelOutput1.ContainsBoundary("LEFT") );

      const size_t topFaceCount( modelOutput1.Boundary("TOP").Cells() );
      const size_t bottomFaceCount( modelOutput1.Boundary("BOTTOM").Cells() );
      const size_t leftFaceCount( modelOutput1.Boundary("LEFT").Cells() );
      const size_t rightFaceCount( modelOutput1.Boundary("RIGHT").Cells() );
      const size_t frontFaceCount( modelOutput1.Boundary("FRONT").Cells() );
      const size_t backFaceCount( modelOutput1.Boundary("BACK").Cells() );

      const size_t topNodeCount( modelOutput1.Boundary("TOP").Nodes() );
      const size_t bottomNodeCount( modelOutput1.Boundary("BOTTOM").Nodes() );
      const size_t leftNodeCount( modelOutput1.Boundary("LEFT").Nodes() );
      const size_t rightNodeCount( modelOutput1.Boundary("RIGHT").Nodes() );
      const size_t frontNodeCount( modelOutput1.Boundary("FRONT").Nodes() );
      const size_t backNodeCount( modelOutput1.Boundary("BACK").Nodes() );
      
      modelOutput1.OutputToBinaryFile("ANSYS_Model3D_Test_modelOutput1");

      //------------------------------------
      // Verification Binary
      //------------------------------------
      timer.Start();
      const string file_name("ANSYS_Model3D_Test_modelOutput1");
      Model<3> modelInput1( file_name, set<string>({}) );
      const Region<3U>& model_domain2(modelInput1.Region("Model"));
      const double binaryModelTime( timer.Stop() );
      _test( elementCount == modelInput1.Region("Model").Cells() );
      _test( nodeCount    == modelInput1.Region("Model").Nodes() );
      
      if ( verbose ) cout <<"\n\nrun: Model reconstructed from file:\n";
      size_t nullNeighbors(0);
      for ( auto it = model_domain2.CellsBegin(); it != model_domain2.CellsEnd(); ++it )
        for ( uint32_t n{0u}; n < (*it)->Neighbors(); ++n )
          if( (*it)->Neighbor(n) == nullptr )
            ++nullNeighbors;

      _test( nullNeighbors == nullNeighborsOut );

      _test( modelInput1.ContainsRegion("Model") );
      _test( modelInput1.ContainsBoundary("LEFT") );
      _test( modelInput1.ContainsBoundary("RIGHT") );
      _test( modelInput1.ContainsBoundary("TOP") );
      _test( modelInput1.ContainsBoundary("FRONT") );
      _test( modelInput1.ContainsBoundary("BACK") );
      _test( modelInput1.ContainsBoundary("BOTTOM") );
      _test( elementCount == modelInput1.Region("Model").Cells() );
      _test( nodeCount == modelInput1.Region("Model").Nodes() );
      _test( regionCount == modelInput1.Regions() );
      _test( matrixLeftNodeCount == modelInput1.Region("MATRIX_LEFT").Nodes() );
      _test( matrixRightNodeCount == modelInput1.Region("MATRIX_RIGHT").Nodes() );
      _test( matrixLeftElementCount == modelInput1.Region("MATRIX_LEFT").Cells() );
      _test( matrixRightElementCount == modelInput1.Region("MATRIX_RIGHT").Cells() );
      _test( topFaceCount == modelInput1.Boundary("TOP").Cells() );
      _test( bottomFaceCount == modelInput1.Boundary("BOTTOM").Cells() );
      _test( leftFaceCount == modelInput1.Boundary("LEFT").Cells() );
      _test( rightFaceCount == modelInput1.Boundary("RIGHT").Cells() );
      _test( frontFaceCount == modelInput1.Boundary("FRONT").Cells() );
      _test( backFaceCount == modelInput1.Boundary("BACK").Cells() );

      const size_t topNodeCountIn( modelInput1.Boundary("TOP").Nodes() );
      const size_t bottomNodeCountIn( modelInput1.Boundary("BOTTOM").Nodes() );
      const size_t leftNodeCountIn( modelInput1.Boundary("LEFT").Nodes() );
      const size_t rightNodeCountIn( modelInput1.Boundary("RIGHT").Nodes() );
      const size_t frontNodeCountIn( modelInput1.Boundary("FRONT").Nodes() );
      const size_t backNodeCountIn( modelInput1.Boundary("BACK").Nodes() );
      
      _test( topNodeCount == topNodeCountIn );
      _test( bottomNodeCount == bottomNodeCountIn );
      _test( leftNodeCount == leftNodeCountIn );
      _test( rightNodeCount == rightNodeCountIn );
      _test( frontNodeCount == frontNodeCountIn );
      _test( backNodeCount == backNodeCountIn );    
      
      _equal( modelInput1.Region("MATRIX_RIGHT").Average("permeability"), matrixRightValue, tolerance );
      _equal( modelInput1.Region("MATRIX_LEFT").Average("permeability"), matrixLeftValue, tolerance );
      modelInput1.InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 1. ) );
      if ( verbose ) cout << "\nICEM time: " << icemModelTime << " vs BINARY time: " << binaryModelTime << endl;
      
      ScalarVariable diff( ANY, 1. );
      VectorVariable<3> vv( DIRICH, 2. );
      VectorVariable<3> vvPlain;
      TensorVariable<3> tv( ANY, 3. );
      TensorVariable<3> tvPlain;
      
      if ( verbose ) cout <<"Building ModelOutput..."<<endl;
      ANSYS_Model3D modelOutput2( "FracBox", "Vset_TestCase.txt",true );
      ArrayVariable na( "nodal array", modelOutput2.Database(), 2., ROBIN );
      Index boundaryScalarKey = modelOutput2.Database().StorageKey("boundary scalar");
      Index boundaryArrayKey = modelOutput2.Database().StorageKey("boundary array");
      Index regionVectorKey = modelOutput2.Database().StorageKey("region vector");
      Index modelTensorKey = modelOutput2.Database().StorageKey("model tensor");
      ArrayVariable zero_array( "boundary array",  modelOutput2.Database() );
      ArrayVariable ba( "boundary array",  modelOutput2.Database() );
      ArrayVariable baPlain( "boundary array",  modelOutput2.Database() );
      zero_array = 0.;
      ba         = 99.;

      //Geometry
      const size_t nodeCount2( modelOutput2.Region("Model").Nodes() );
      const size_t regionCount2( modelOutput2.Regions() );
      const size_t boundaryCount2( modelOutput2.Boundaries() );

      //Variables
      modelOutput2.InputPropertyValue( "diffusivity", diff );
      modelOutput2.Region("Model").Store( regionVectorKey, vv );
      modelOutput2.InputPropertyValue( "nodal array", na );
      modelOutput2.Store( modelTensorKey, tv );
      modelOutput2.InputPropertyValue( "boundary scalar", makeScalar( PLAIN, 0. ) );
      modelOutput2.Boundary("BOUNDARY1").InputPropertyValue("boundary scalar", makeScalar( PLAIN, 1. ) );
      modelOutput2.InputPropertyValue( "boundary array", zero_array );
      modelOutput2.Boundary("BOUNDARY2").InputPropertyValue("boundary array", ba );

      modelOutput2.OutputToBinaryFile("ANSYS_Model3D_Test_modelOutput2");
      if ( verbose ) cout<<"The Output of Model with boundaries done..."<<endl;

      cout <<"Re-building ModelInput with boundaries..."<<endl;
      const string file_name2("ANSYS_Model3D_Test_modelOutput2");
      Model<3U> modelInput2( file_name2, set<string>({}) );
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

      // gui consturctor
      ANSYS_Model3D guiModel1( "FracBox", true );
      _test( guiModel1.UniqueRegions() == 2 ); // FRACTURES and MATRIX
      _test( guiModel1.Boundaries() == 6 );

      // irregular shaped
      ANSYS_Model3D modelI( "FracBox", "CSMP-variables.txt", true );
      _test( modelI.ContainsBoundary( "BOUNDARY1" ) );
      _test( modelI.ContainsBoundary( "BOUNDARY2" ) );
      _test( modelI.ContainsBoundary( "BOUNDARY3" ) );
      _test( modelI.ContainsBoundary( "BOUNDARY4" ) );
      _test( modelI.ContainsBoundary( "BOUNDARY5" ) );
      _test( modelI.ContainsBoundary( "BOUNDARY6" ) );

      _test( modelI.ContainsRegion( "MATRIX" ) );
      _test( modelI.ContainsRegion( "FRACTURE" ) );

      // box shaped
      ANSYS_Model3D model( "BoxHalfs3D", "CSMP-variables.txt" );
      Index nodalKey( model.Database().StorageKey( "nodal variable" ) );
      Region<3>& rref( model.Region( "Model" ) );
      rref.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 0. ) );
      size_t boundaryNodeCount( 0 );
      const vector<Node<3>*>::const_iterator nodesEnd( rref.NodesEnd() );
      for( vector<Node<3>*>::const_iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
        {
        const BOX_BOUNDARY boxBoundary( (*it)->AtBoundary() );
        if( boxBoundary == LEFT or boxBoundary == RIGHT or boxBoundary == TOP or
            boxBoundary == BOTTOM  or boxBoundary == FRONT or boxBoundary == BACK )
          (*it)->Store( nodalKey, makeScalar( PLAIN, 1.0 ) );
        else
          continue;   
        ++boundaryNodeCount;
        }
      if ( verbose ) cout << "\nBoundary Node Count: " << boundaryNodeCount << endl;
      _test ( boundaryNodeCount > 0 );
      VTU_Interface<3> vtu( model ); vtu.OmitZeroInFileName( true );
      //vtu.OutputDataToVTU( "ANSYS_Model3D_Test-PerimeterNodes", "nodal variable" );

      _test( model.ContainsBoundary( "LEFT" ) );
      _test( model.ContainsBoundary( "RIGHT" ) );
      _test( model.ContainsBoundary( "BOTTOM" ) );
      _test( model.ContainsBoundary( "TOP" ) );
      _test( model.ContainsBoundary( "FRONT" ) );
      _test( model.ContainsBoundary( "BACK" ) );

      _test( model.ContainsRegion( "MATRIX_LEFT" ) );
      _test( model.ContainsRegion( "MATRIX_RIGHT" ) );
      _test( model.ContainsRegion( "HALF" ) );

      Boundary<3>& left( model.Boundary( std::string("LEFT") ) );
      Boundary<3>& right( model.Boundary( std::string("RIGHT") ) );
      Boundary<3>& bottom( model.Boundary( std::string("BOTTOM") ) );
      Boundary<3>& top( model.Boundary( std::string("TOP") ) );
      Boundary<3>& front( model.Boundary( std::string("FRONT") ) );
      Boundary<3>& back( model.Boundary( std::string("BACK") ) );

      _test( left.Cells() != 0 );
      _test( right.Cells() != 0 );
      _test( bottom.Cells() != 0 );
      _test( top.Cells() != 0 );
      _test( front.Cells() != 0 );
      _test( back.Cells() != 0 );
      
      // testing model with finite volume variables
      ANSYS_Model3D modelOutput3( "FracBox", "FracBoxNoFrac", "Vset_TestCase.txt",true);

      Index faipVectorKey( modelOutput3.Database().StorageKey("faip vector") );
      Index seipTensorKey( modelOutput3.Database().StorageKey("seip tensor") );
      const size_t modelOutput3_boundar1NodesOut( modelOutput3.Boundary("BOUNDARY1").Nodes() );
      const size_t modelOutput3_boundar2NodesOut( modelOutput3.Boundary("BOUNDARY2").Nodes() );
      const size_t modelOutput3_boundar3NodesOut( modelOutput3.Boundary("BOUNDARY3").Nodes() );
      const size_t modelOutput3_boundar4NodesOut( modelOutput3.Boundary("BOUNDARY4").Nodes() );
      const size_t modelOutput3_boundar5NodesOut( modelOutput3.Boundary("BOUNDARY5").Nodes() );
      const size_t modelOutput3_boundar6NodesOut( modelOutput3.Boundary("BOUNDARY6").Nodes() );

      modelOutput3.InputPropertyValue( "faip vector", vv );
      modelOutput3.InputPropertyValue( "seip tensor", tv );
      Element<3>* ePtr = *modelOutput3.Region("Model").CellsBegin();
      for( uint32_t f(0u); f < ePtr->Facets(); ++f )
        for( uint32_t fip(0u); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
        }
        modelOutput3.OutputToBinaryFile("ANSYS_Model3D_Test_modelOutput3");

        const string file_name3("ANSYS_Model3D_Test_modelOutput3");
        Model<3> modelInput3( file_name3, set<string>({}) );
        _test( modelInput3.Boundary("BOUNDARY1").Nodes() == modelOutput3_boundar1NodesOut );
        _test( modelInput3.Boundary("BOUNDARY2").Nodes() == modelOutput3_boundar2NodesOut );
        _test( modelInput3.Boundary("BOUNDARY3").Nodes() == modelOutput3_boundar3NodesOut );
        _test( modelInput3.Boundary("BOUNDARY4").Nodes() == modelOutput3_boundar4NodesOut );
        _test( modelInput3.Boundary("BOUNDARY5").Nodes() == modelOutput3_boundar5NodesOut );
        _test( modelInput3.Boundary("BOUNDARY6").Nodes() == modelOutput3_boundar6NodesOut );

        ePtr = *modelInput3.Region("Model").CellsBegin();
        size_t ctrFaIps(0);
        for( uint32_t f(0u); f < ePtr->Facets(); ++f )
          for( uint32_t fip(0u); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
            ++ctrFaIps;
          }

          _test( ctrFaIps == 6 );
          
          size_t ctrSeIps(0);
          for( uint32_t s(0); s < ePtr->Sectors(); ++s )
            for( uint32_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
            {
              ePtr->Read( s, sip, seipTensorKey, tvPlain );
              _test( tvPlain == tv );
              ++ctrSeIps;
            }
         _test( ctrSeIps == 4 );

        if ( verbose ) {
            cout << "\nElement Type: " <<parseFiniteElementType( ePtr->FE_Type() ) << endl;
            cout << "\nSectors: " << ePtr->Sectors() << endl;
            cout << "\nFacets: " << ePtr->Facets() << endl;
            cout << "\nFacet IPs accessed: " << ctrFaIps << endl;
            cout << "\nSector IPs accessed: " << ctrSeIps << endl;
          }

        ctrFaIps = 0;
        ctrSeIps = 0;

        for( uint32_t f(0); f < ePtr->Facets(); ++f )
          for( uint32_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
            ++ctrFaIps;
          }

          for( uint32_t s(0); s < ePtr->Sectors(); ++s )
            for( uint32_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
            {
              ePtr->Read( s, sip, seipTensorKey, tvPlain );
              _test( tvPlain == tv );
               ++ctrSeIps;
            }
            
     _test( ctrSeIps == 4 );
     _test( ctrFaIps == 6 );
     
      // Testing input / output of the model with quadratic FEM basis functions (midside nodes)
      Test_ReadWriteQuadraticFEM_Model();
     
       if ( verbose ) cout <<"\n\n"<<this->getName()<<" FINISHED!!!"<<endl;
    }
    
    
    
    
    
    


void ANSYS_Model3D_Test::ModelRecoveryFromFileTest()
 {
    // Test: box-shaped model with boundary information (Boundary->ModelSubDomain)
    // -----------------------------------------------------------------------------
         ANSYS_Model3D model1( "cube_flag", "CSMP-variables.txt",
                                 true   /* binary_file */
                             );
                             
        _test( model1.EstablishEdgeBoundariesOfBoxShapedModel() );
        
         // loop over model boundary verifying consistency between AtBoundary() and region flags
         bool left_fail(false), right_fail(false), bottom_fail(false), top_fail(false), front_fail(false), back_fail(false);
	       for ( auto it = model1.BoundariesBegin(); it != model1.BoundariesEnd(); ++it ) {
               const string bname((*it).first);
               for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.PerimeterNodesBegin(); ++nit ) {
                     if ( bname == "LEFT"   and (*nit)->AtBoundary() != LEFT )   left_fail = true;
                     if ( bname == "RIGHT"  and (*nit)->AtBoundary() != RIGHT )  right_fail = true;
                     if ( bname == "TOP"    and (*nit)->AtBoundary() != TOP )    top_fail = true;
                     if ( bname == "BOTTOM" and (*nit)->AtBoundary() != BOTTOM ) bottom_fail = true;
                     if ( bname == "FRONT"  and (*nit)->AtBoundary() != FRONT )  front_fail = true;
                     if ( bname == "BACK"   and (*nit)->AtBoundary() != BACK )   back_fail = true;
                 }
            }
        
         // testing the side boundaries
         _test( left_fail   == false );
         _test( right_fail  == false );
         _test( bottom_fail == false );
         _test( top_fail    == false );
         _test( front_fail  == false );
         _test( back_fail   == false );
        
         // testing the edges and their perimeter nodes=model corners as well
         bool edge1_fail(false), edge2_fail(false), edge3_fail(false), edge4_fail(false),
              edge5_fail(false), edge6_fail(false), edge7_fail(false), edge8_fail(false),
              edge9_fail(false), edge10_fail(false), edge11_fail(false), edge12_fail(false);

         // E.P Need to create EDGE boundaries in the first place for this to Work!!! NEED TO CALL METHOD BELOW FIRST
         // model1.EstablishEdgeBoundariesOfBoxShapedModel();


	       for ( auto it = model1.BoundariesBegin(); it != model1.BoundariesEnd(); ++it ) {
               const string bname((*it).first);
               // the interior part of the edges comes first
               for ( auto nit=(*it).second.NodesBegin(); nit!=(*it).second.PerimeterNodesBegin(); ++nit ) {
                     if ( bname == "EDGE1" and (*nit)->AtBoundary() != EDGE1 ) edge1_fail = true;
                     if ( bname == "EDGE2" and (*nit)->AtBoundary() != EDGE2 ) edge2_fail = true;
                     if ( bname == "EDGE3" and (*nit)->AtBoundary() != EDGE3 ) edge3_fail = true;
                     if ( bname == "EDGE4" and (*nit)->AtBoundary() != EDGE4 ) edge4_fail = true;
                     if ( bname == "EDGE5" and (*nit)->AtBoundary() != EDGE5 ) edge5_fail = true;
                     if ( bname == "EDGE6" and (*nit)->AtBoundary() != EDGE6 ) edge6_fail = true;
                     if ( bname == "EDGE7" and (*nit)->AtBoundary() != EDGE7 ) edge7_fail = true;
                     if ( bname == "EDGE8" and (*nit)->AtBoundary() != EDGE8 ) edge8_fail = true;
                     if ( bname == "EDGE9" and (*nit)->AtBoundary() != EDGE9 ) edge9_fail = true;
                     if ( bname == "EDGE10" and (*nit)->AtBoundary() != EDGE10 ) edge10_fail = true;
                     if ( bname == "EDGE11" and (*nit)->AtBoundary() != EDGE11 ) edge11_fail = true;
                     if ( bname == "EDGE12" and (*nit)->AtBoundary() != EDGE12 ) edge12_fail = true;
                 }
            }

         // testing the side boundaries
         _test( edge1_fail == false );
         _test( edge2_fail == false );
         _test( edge3_fail == false );
         _test( edge4_fail == false );
         _test( edge5_fail == false );
         _test( edge6_fail == false );
         _test( edge7_fail == false );
         _test( edge8_fail == false );
         _test( edge9_fail == false );
         _test( edge10_fail == false );
         _test( edge11_fail == false );
         _test( edge12_fail == false );
        
         // verifying that the perimeter nodes of the edges contain the correct corner nodes
         set<BOX_BOUNDARY> bnodes;
         // testing edge1
         const Boundary<3U> edge1(model1.Boundary("EDGE1"));
         _test( edge1.PerimeterNodes() == 2 );
         for ( auto nit=edge1.PerimeterNodesBegin(); nit!=edge1.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR1) > 0 );
         _test( bnodes.count(CNR2) > 0 );
         bnodes.clear();
         // testing edge2
         const Boundary<3U> edge2(model1.Boundary("EDGE2"));
         _test( edge2.PerimeterNodes() == 2 );
         for ( auto nit=edge2.PerimeterNodesBegin(); nit!=edge2.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR3) > 0 );
         _test( bnodes.count(CNR2) > 0 );
         bnodes.clear();
         // testing edge3
         const Boundary<3U> edge3(model1.Boundary("EDGE3"));
         _test( edge3.PerimeterNodes() == 2 );
         for ( auto nit=edge3.PerimeterNodesBegin(); nit!=edge3.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR3) > 0 );
         _test( bnodes.count(CNR4) > 0 );
         bnodes.clear();
         // testing edge4
         const Boundary<3U> edge4(model1.Boundary("EDGE4"));
         _test( edge4.PerimeterNodes() == 2 );
         for ( auto nit=edge4.PerimeterNodesBegin(); nit!=edge4.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR1) > 0 );
         _test( bnodes.count(CNR4) > 0 );
         bnodes.clear();
         // testing edge5
         const Boundary<3U> edge5(model1.Boundary("EDGE5"));
         _test( edge5.PerimeterNodes() == 2 );
         for ( auto nit=edge5.PerimeterNodesBegin(); nit!=edge5.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR1) > 0 );
         _test( bnodes.count(CNR5) > 0 );
         bnodes.clear();
         // testing edge6
         const Boundary<3U> edge6(model1.Boundary("EDGE6"));
         _test( edge6.PerimeterNodes() == 2 );
         for ( auto nit=edge6.PerimeterNodesBegin(); nit!=edge6.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR2) > 0 );
         _test( bnodes.count(CNR6) > 0 );
         bnodes.clear();
         // testing edge7
         const Boundary<3U> edge7(model1.Boundary("EDGE7"));
         _test( edge7.PerimeterNodes() == 2 );
         for ( auto nit=edge7.PerimeterNodesBegin(); nit!=edge7.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR3) > 0 );
         _test( bnodes.count(CNR7) > 0 );
         bnodes.clear();
         // testing edge1
         const Boundary<3U> edge8(model1.Boundary("EDGE8"));
         _test( edge8.PerimeterNodes() == 2 );
         for ( auto nit=edge8.PerimeterNodesBegin(); nit!=edge8.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR4) > 0 );
         _test( bnodes.count(CNR8) > 0 );
         bnodes.clear();
         // testing edge9
         const Boundary<3U> edge9(model1.Boundary("EDGE9"));
         _test( edge9.PerimeterNodes() == 2 );
         for ( auto nit=edge9.PerimeterNodesBegin(); nit!=edge9.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR5) > 0 );
         _test( bnodes.count(CNR6) > 0 );
         bnodes.clear();
         // testing edge10
         const Boundary<3U> edge10(model1.Boundary("EDGE10"));
         _test( edge10.PerimeterNodes() == 2 );
         for ( auto nit=edge10.PerimeterNodesBegin(); nit!=edge10.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR6) > 0 );
         _test( bnodes.count(CNR7) > 0 );
         bnodes.clear();
         // testing edge11
         const Boundary<3U> edge11(model1.Boundary("EDGE11"));
         _test( edge11.PerimeterNodes() == 2 );
         for ( auto nit=edge11.PerimeterNodesBegin(); nit!=edge11.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR7) > 0 );
         _test( bnodes.count(CNR8) > 0 );
         bnodes.clear();
         // testing edge12
         const Boundary<3U> edge12(model1.Boundary("EDGE12"));
         _test( edge1.PerimeterNodes() == 2 );
         for ( auto nit=edge12.PerimeterNodesBegin(); nit!=edge12.NodesEnd(); ++nit )
           bnodes.insert( (*nit)->AtBoundary() );
         _test( bnodes.count(CNR8) > 0 );
         _test( bnodes.count(CNR5) > 0 );
         bnodes.clear();

         // mapping flags to values to test assignments
// TODO: adding properties upsets indices for the access of the flagged array
// call  UpdateParametersAndDatabase(); or  UpdateIndexReferences(); but they are private?

/* CREATE PROPERTY UPSETS flagged array variable storage offset
         model1.CreateProperty( "box flag", "none", SCALAR, NODE );
         model1.CreateProperty( "box flag element", "none", SCALAR, ELEMENT );
         boxFlagsToVariable( model1, "box flag", "box flag element" );

         VTK_Interface<3U>  vtk_output;
         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
        
         // checking whether the side boundary interior and boundary flags are identified correctly
         //const csmp::Index  prop_key(model1.Database().StorageKey("box flag"));
         Boundary<3U>&      bref(model1.Boundary("FRONT"));
         bref.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(FRONT)), INTERIOR );
         bref.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(REGION_BOUNDARY)), PERIMETER );

         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 1 );

         Boundary<3U>&      brefl(model1.Boundary("LEFT"));
         brefl.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(LEFT)), INTERIOR );
         brefl.InputPropertyValue( "box flag", makeScalar(ANY,static_cast<double>(REGION_BOUNDARY)), PERIMETER );

         vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 2 );
        
         // testing by comparison with BOX_BOUNDARY flags
*/

         // saving the model to csmp native binary file format
         model1.OutputToBinaryFile("ModelSubDomain_Test2");
        
         // bringing the model back from binary file
         Model<3U>  model2( string("ModelSubDomain_Test2") );

         ModelSubDomain_Test comparitor;

         _test( comparitor.CompareModelSubdomains( model1.Region("Model"), model2.Region("Model"), true ) );


   // TESTING THE SAME WITH  A MORE COMPLEX MODEL
      {
         ANSYS_Model3D modelOut( "prism_test", "CSMP-variables.txt",
                                 true   /* binary_file */
                             );
        
         model1.OutputToBinaryFile("ModelSubDomain_Test3");
         Model<3U>  modelIn( string("ModelSubDomain_Test3") );

         _test( comparitor.CompareModelSubdomains( modelOut.Region("FRAC_VOLUMES"),
                                                   modelIn.Region("FRAC_VOLUMES"), true ) );
      }

} // end ModelRecoveryFromFileTest




void ANSYS_Model3D_Test::Test_ReadWriteQuadraticFEM_Model()
 {
   // box with a vertical well in the middle that is meshed volumetrically
   // 4617 nodes, 3878 elements, 8 families                 :
   ANSYS_Model3D model( "box_with_hole2", "CSMP-variables.txt",  true /* binary_file */ );
   model.OutputToBinaryFile("Test_ReadWriteQuadraticFEM_Model");

   Model<3U>  modelIn( string("Test_ReadWriteQuadraticFEM_Model") );

   ModelSubDomain_Test comparitor;
   _test( comparitor.CompareModelSubdomains( model.Region("WELL"),
                                             modelIn.Region("WELL"), true ) );
 }




} // csmp
