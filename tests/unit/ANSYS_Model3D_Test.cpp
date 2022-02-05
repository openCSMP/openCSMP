#include "ANSYS_Model3D_Test.h"
#include "ANSYS_Model3D.h"
#include "Boundary.h"
#include "Region.h"
#include "VTU_Interface.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "Timer.hpp"

using namespace std;

namespace csmp
  {

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
      ANSYS_Model3D modelOutput1( "BoxHalfs3D", "BoxHalfs3Dirregular", variablesFile.c_str(), true );
      const double icemModelTime( timer.Stop() );

      size_t nullNeighborsOut(0);
      const Region<3U> model_domain1(modelOutput1.Region("Model"));
      for ( auto it = model_domain1.ElementsBegin(); it != model_domain1.ElementsEnd(); ++it )
        for ( size_t n(0); n < (*it)->Neighbors(); ++n )
          if( (*it)->Neighbor(n) == nullptr )
            ++nullNeighborsOut;

      _test( nullNeighborsOut != 0 );

      const double matrixLeftValue(2.);
      const double matrixRightValue(3.);
      const double tolerance(1.0e-5);
      modelOutput1.Region("Model").InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 1.0 ) );
      modelOutput1.Region("MATRIX_LEFT").InputPropertyValue( "permeability", makeScalar( PLAIN, matrixLeftValue ) );
      modelOutput1.Region("MATRIX_RIGHT").InputPropertyValue( "permeability", makeScalar( PLAIN, matrixRightValue ) );

      const size_t elementCount( modelOutput1.Region("Model").Elements() );
      const size_t nodeCount( modelOutput1.Region("Model").Nodes() );
      const size_t regionCount( modelOutput1.Regions() );
      const size_t matrixLeftNodeCount(  modelOutput1.Region("MATRIX_LEFT").Nodes() );
      const size_t matrixRightNodeCount(  modelOutput1.Region("MATRIX_RIGHT").Nodes() );
      const size_t matrixLeftElementCount(  modelOutput1.Region("MATRIX_LEFT").Elements() );
      const size_t matrixRightElementCount(  modelOutput1.Region("MATRIX_RIGHT").Elements() );

      _test( modelOutput1.ContainsBoundary("TOP") );
      _test( modelOutput1.ContainsBoundary("BOTTOM") );
      _test( modelOutput1.ContainsBoundary("RIGHT") );
      _test( modelOutput1.ContainsBoundary("FRONT") );
      _test( modelOutput1.ContainsBoundary("BACK") );
      _test( modelOutput1.ContainsBoundary("LEFT") );

      const size_t topFaceCount( modelOutput1.Boundary("TOP").Elements() );
      const size_t bottomFaceCount( modelOutput1.Boundary("BOTTOM").Elements() );
      const size_t leftFaceCount( modelOutput1.Boundary("LEFT").Elements() );
      const size_t rightFaceCount( modelOutput1.Boundary("RIGHT").Elements() );
      const size_t frontFaceCount( modelOutput1.Boundary("FRONT").Elements() );
      const size_t backFaceCount( modelOutput1.Boundary("BACK").Elements() );

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
      const Region<3U> model_domain2(modelInput1.Region("Model"));
      const double binaryModelTime( timer.Stop() );
      _test( elementCount == modelInput1.Region("Model").Elements() );
      _test( nodeCount    == modelInput1.Region("Model").Nodes() );
      
      if ( verbose ) cout <<"\n\nrun: Model reconstructed from file:\n";
      size_t nullNeighbors(0);
      for ( auto it = model_domain2.ElementsBegin(); it != model_domain2.ElementsEnd(); ++it )
        for ( size_t n(0); n < (*it)->Neighbors(); ++n )
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
      _test( elementCount == modelInput1.Region("Model").Elements() );
      _test( nodeCount == modelInput1.Region("Model").Nodes() );
      _test( regionCount == modelInput1.Regions() );
      _test( matrixLeftNodeCount == modelInput1.Region("MATRIX_LEFT").Nodes() );
      _test( matrixRightNodeCount == modelInput1.Region("MATRIX_RIGHT").Nodes() );
      _test( matrixLeftElementCount == modelInput1.Region("MATRIX_LEFT").Elements() );
      _test( matrixRightElementCount == modelInput1.Region("MATRIX_RIGHT").Elements() );
      _test( topFaceCount == modelInput1.Boundary("TOP").Elements() );
      _test( bottomFaceCount == modelInput1.Boundary("BOTTOM").Elements() );
      _test( leftFaceCount == modelInput1.Boundary("LEFT").Elements() );
      _test( rightFaceCount == modelInput1.Boundary("RIGHT").Elements() );
      _test( frontFaceCount == modelInput1.Boundary("FRONT").Elements() );
      _test( backFaceCount == modelInput1.Boundary("BACK").Elements() );

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
      ANSYS_Model3D modelOutput2( "FracBox", "Vset_TestCase.txt",true,true,true );
      ArrayVariable na( "nodal array", modelOutput2.Database(), 2., ROBIN );
      Index boundaryScalarKey = modelOutput2.Database().StorageKey("boundary scalar");
      Index boundaryArrayKey = modelOutput2.Database().StorageKey("boundary array");
      Index regionVectorKey = modelOutput2.Database().StorageKey("region vector");
      Index modelTensorKey = modelOutput2.Database().StorageKey("model tensor");
      ArrayVariable ba( "boundary array",  modelOutput2.Database() );
      ArrayVariable baPlain( "boundary array",  modelOutput2.Database() );
      ba = 99.;

      //Geometry
      const size_t nodeCount2( modelOutput2.Region("Model").Nodes() );
      const size_t regionCount2( modelOutput2.Regions() );
      const size_t boundaryCount2( modelOutput2.Boundaries() );

      //Variables
      modelOutput2.Region("Model").InputPropertyValue( "diffusivity", diff );
      modelOutput2.Region("Model").Store( regionVectorKey, vv );
      modelOutput2.Region("Model").InputPropertyValue( "nodal array", na );
      modelOutput2.Store( modelTensorKey, tv );
      modelOutput2.Boundary("BOUNDARY1").InputPropertyValue("boundary scalar", makeScalar( PLAIN, 1. ) );
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
      ANSYS_Model3D guiModel1( "FracBox", true, true );
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
      const vector<Node<3>*>::iterator nodesEnd( rref.NodesEnd() );
      for( vector<Node<3>*>::iterator it = rref.NodesBegin(); it != nodesEnd; ++it )
        {
        const BOX_BOUNDARY boxBoundary( (*it)->AtBoundary() );
        if( boxBoundary == LEFT or boxBoundary == RIGHT or boxBoundary == TOP or boxBoundary == BOTTOM  or boxBoundary == FRONT or boxBoundary == BACK )
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

      _test( left.Elements() != 0 );
      _test( right.Elements() != 0 );
      _test( bottom.Elements() != 0 );
      _test( top.Elements() != 0 );
      _test( front.Elements() != 0 );
      _test( back.Elements() != 0 );
      
      // testing model with finite volume variables
      ANSYS_Model3D modelOutput3( "FracBox", "FracBoxNoFrac", "Vset_TestCase.txt",true,true,true);

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
      Element<3>* ePtr = *modelOutput3.Region("Model").ElementsBegin();
      for( size_t f(0); f < ePtr->Facets(); ++f )
        for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
        }
        modelOutput3.OutputToBinaryFile("ANSYS_Model3D_Test_modelOutput3");

        const string file_name3("ANSYS_Model3D_Test_modelOutput3");
        Model<3> modelInput3( file_name, set<string>({}) );
        _test( modelInput3.Boundary("BOUNDARY1").Nodes() == modelOutput3_boundar1NodesOut );
        _test( modelInput3.Boundary("BOUNDARY2").Nodes() == modelOutput3_boundar2NodesOut );
        _test( modelInput3.Boundary("BOUNDARY3").Nodes() == modelOutput3_boundar3NodesOut );
        _test( modelInput3.Boundary("BOUNDARY4").Nodes() == modelOutput3_boundar4NodesOut );
        _test( modelInput3.Boundary("BOUNDARY5").Nodes() == modelOutput3_boundar5NodesOut );
        _test( modelInput3.Boundary("BOUNDARY6").Nodes() == modelOutput3_boundar6NodesOut );

        ePtr = *modelInput3.Region("Model").ElementsBegin();
        size_t ctrFaIps(0);
        for( size_t f(0); f < ePtr->Facets(); ++f )
          for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
            ++ctrFaIps;
          }

          _test( ctrFaIps == 6 );
          
          size_t ctrSeIps(0);
          for( size_t s(0); s < ePtr->Sectors(); ++s )
            for( size_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
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

        for( size_t f(0); f < ePtr->Facets(); ++f )
          for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
          {
            ePtr->Read( f, fip, faipVectorKey, vvPlain );
            _test( vvPlain == vv );
            ++ctrFaIps;
          }

          for( size_t s(0); s < ePtr->Sectors(); ++s )
            for( size_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
            {
              ePtr->Read( s, sip, seipTensorKey, tvPlain );
              _test( tvPlain == tv );
               ++ctrSeIps;
            }
            
     _test( ctrSeIps == 4 );
     _test( ctrFaIps == 6 );

       if ( verbose ) cout <<"\n\n"<<this->getName()<<" FINISHED!!!"<<endl;
    }

  } // csmp
