#include "Variables_Test.h"
#include "Boundary.h"

#include "ANSYS_Model3D.h"
#include "NodeCenteredFiniteVolumeTransport.h"

using namespace std;

namespace csmp {

/**
     @param prefix is the model name. 
*/
Variables_Test::Variables_Test( const char* prefix )
  :  prefix_(prefix)
  {
    this->setName("Variables_Test");
  }
  
  
struct IndexTrackerTestStruct
  {
    csmp::Index key1_, key2_;
  };




void Variables_Test::run()
  {
    // Run Test for 3D Model constructed by ANSYS mesh reader
	  const string variables_filename = (string)(this->getName() + ".txt");
    ANSYS_Model3D m0(prefix_, variables_filename.c_str(), true, true, true, true );
    m0.OutputToBinaryFile("Variables_Test_BinaryModel");
    runModel(m0);

    // Run Test for 3D Model loaded from CSMP++ binary format
    Model<3> m1( string("Variables_Test_BinaryModel") );
    runModel(m1);
  }
  
  
  
/**
 @fn  void Variables_Test::run()

 @todo (2-F) Test for InterFace and SplitBoundary properties

 @author  P. Lang
 @date  9/24/2012
 */
void Variables_Test::runModel( Model<3>& model )
    {
      // Indices
      IndexTrackerTestStruct indexCache;
      // ========================================================================
      // PLACEMENT: Model:
      // TYPE: Scalar
      Index mKey1       ( model.Database().StorageKey("model variable 1") );
      // TYPE: Array
      Index maKey1      ( model.Database().StorageKey("model array 1") );
      // TYPE: Flagged Array
      Index mfaKey1     ( model.Database().StorageKey("model flagged array 1") );

      // ========================================================================
      // PLACEMENT: Subdomain ( Region, Boundary ):
      // TYPE: Vector property
      Index rvKey1      ( model.Database().StorageKey("region vector 1") );
      // Tensor property
      Index btKey1      ( model.Database().StorageKey("boundary tensor 1") );

      // ========================================================================
      // PLACEMENT: Node:
      // TYPE: Scalar
      Index nKey1       ( model.Database().StorageKey("nodal variable 1") );
      Index nKey2       ( model.Database().StorageKey("nodal variable 2") );
      // TYPE: Vector
      Index nvKey1      ( model.Database().StorageKey("nodal vector 1") );
      Index nvKey2      ( model.Database().StorageKey("nodal vector 2") );
      // TYPE: Tensor
      Index ntKey1      ( model.Database().StorageKey("nodal tensor 1") );
      Index ntKey2      ( model.Database().StorageKey("nodal tensor 2") );
      // TYPE: Array
      Index naKey1      ( model.Database().StorageKey("nodal array 1") );
      Index naKey2      ( model.Database().StorageKey("nodal array 2") );
      Index naKey3      ( model.Database().StorageKey("nodal array 3") );
      // TYPE: Flagged Array
      Index nfaKey1     ( model.Database().StorageKey("nodal flagged array 1") );
      Index nfaKey2     ( model.Database().StorageKey("nodal flagged array 2") );
      Index nfaKey3     ( model.Database().StorageKey("nodal flagged array 3") );

      // ========================================================================
      // PLACEMENt: Element
      // TYPE: Scalar
      Index eKey1       ( model.Database().StorageKey("element variable 1") );
      Index eKey2       ( model.Database().StorageKey("element variable 2") );
      // TYPE: Vector
      Index evKey1      ( model.Database().StorageKey("element vector 1") );
      Index evKey2      ( model.Database().StorageKey("element vector 2") );
      // TYPE: Tensor
      Index etKey1      ( model.Database().StorageKey("element tensor 1") );
      Index etKey2      ( model.Database().StorageKey("element tensor 2") );
      // TYPE: Array
      Index eaKey1      ( model.Database().StorageKey("element array 1") );
      Index eaKey2      ( model.Database().StorageKey("element array 2") );
      Index eaKey3      ( model.Database().StorageKey("element array 3") );
      // TYPE: Flagged Array
      Index efaKey1     ( model.Database().StorageKey("element flagged array 1") );
      Index efaKey2     ( model.Database().StorageKey("element flagged array 2") );
      Index efaKey3     ( model.Database().StorageKey("element flagged array 3") );
      indexCache.key1_ = eaKey1;

      // ========================================================================
      // PLACEMENT: FACE
      // TYPE: Scalar
      Index fKey1       ( model.Database().StorageKey("face variable 1") );
      Index fKey2       ( model.Database().StorageKey("face variable 2") );
      // TYPE: Vector
      Index fvKey1      ( model.Database().StorageKey("face vector 1") );
      Index fvKey2      ( model.Database().StorageKey("face vector 2") );
      // TYPE: Tensor
      Index ftKey1      ( model.Database().StorageKey("face tensor 1") );
      Index ftKey2      ( model.Database().StorageKey("face tensor 2") );
      // TYPE: Array
      Index faKey1      ( model.Database().StorageKey("face array 1") );
      Index faKey2      ( model.Database().StorageKey("face array 2") );
      Index faKey3      ( model.Database().StorageKey("face array 3") );
      // TYPE: Flagged Array
      Index ffaKey1     ( model.Database().StorageKey("face flagged array 1") );
      Index ffaKey2     ( model.Database().StorageKey("face flagged array 2") );
      Index ffaKey3     ( model.Database().StorageKey("face flagged array 3") );

      // ========================================================================
      // PLACEMENT: Element Integration Point
      // TYPE: Scalar
      Index eipKey1     ( model.Database().StorageKey("eip scalar 1") );
      Index eipKey2     ( model.Database().StorageKey("eip scalar 2") );
      // TYPE: Vector
      Index eipvKey1    ( model.Database().StorageKey("eip vector 1") );
      Index eipvKey2    ( model.Database().StorageKey("eip vector 2") );
      // TYPE: Tensor
      Index eiptKey1    ( model.Database().StorageKey("eip tensor 1") );
      Index eiptKey2    ( model.Database().StorageKey("eip tensor 2") );
      // TYPE: Array
      Index eipaKey1    ( model.Database().StorageKey("eip array 1") );
      Index eipaKey2    ( model.Database().StorageKey("eip array 2") );
      Index eipaKey3    ( model.Database().StorageKey("eip array 3") );
      // TYPE: Flagged Array
      Index eipfaKey1   ( model.Database().StorageKey("eip flagged array 1") );
      Index eipfaKey2   ( model.Database().StorageKey("eip flagged array 2") );
      Index eipfaKey3   ( model.Database().StorageKey("eip flagged array 3") );
      indexCache.key2_ = model.Database().StorageKey("eip array 1");

      // ========================================================================
      // PLACEMENT: Face Integration Point
      // TYPE: Scalar
      Index fipKey1     ( model.Database().StorageKey("fip scalar 1") );
      Index fipKey2     ( model.Database().StorageKey("fip scalar 2") );
      // TYPE: Vector
      Index fipvKey1    ( model.Database().StorageKey("fip vector 1") );
      Index fipvKey2    ( model.Database().StorageKey("fip vector 2") );
      // TYPE: Tensor
      Index fiptKey1    ( model.Database().StorageKey("fip tensor 1") );
      Index fiptKey2    ( model.Database().StorageKey("fip tensor 2") );
      // TYPE: Array
      Index fipaKey1    ( model.Database().StorageKey("fip array 1") );
      Index fipaKey2    ( model.Database().StorageKey("fip array 2") );
      Index fipaKey3    ( model.Database().StorageKey("fip array 3") );
      // TYPE: Flagged Array
      Index fipfaKey1   ( model.Database().StorageKey("fip flagged array 1") );
      Index fipfaKey2   ( model.Database().StorageKey("fip flagged array 2") );
      Index fipfaKey3   ( model.Database().StorageKey("fip flagged array 3") );

      // ========================================================================
      // PLACEMENT: Sector Integration Point
      // TYPE: Scalar
      Index seipKey1    ( model.Database().StorageKey("seip scalar 1") );
      Index seipKey2    ( model.Database().StorageKey("seip scalar 2") );
      // TYPE: Vector
      Index seipvKey1   ( model.Database().StorageKey("seip vector 1") );
      Index seipvKey2   ( model.Database().StorageKey("seip vector 2") );
      // TYPE: Tensor
      Index seiptKey1   ( model.Database().StorageKey("seip tensor 1") );
      Index seiptKey2   ( model.Database().StorageKey("seip tensor 2") );
      // TYPE: Array
      Index seipaKey1   ( model.Database().StorageKey("seip array 1") );
      Index seipaKey2   ( model.Database().StorageKey("seip array 2") );
      Index seipaKey3   ( model.Database().StorageKey("seip array 3") );
      // TYPE: Flagged Array
      Index seipfaKey1  ( model.Database().StorageKey("seip flagged array 1") );
      Index seipfaKey2  ( model.Database().StorageKey("seip flagged array 2") );
      Index seipfaKey3  ( model.Database().StorageKey("seip flagged array 3") );

      // ========================================================================
      // PLACEMENT: Facet Integration Point
      // TYPE: Scalar
      Index faipKey1    ( model.Database().StorageKey("faip scalar 1") );
      Index faipKey2    ( model.Database().StorageKey("faip scalar 2") );
      // TYPE: Vector
      Index faipvKey1   ( model.Database().StorageKey("faip vector 1") );
      Index faipvKey2   ( model.Database().StorageKey("faip vector 2") );
      // TYPE: Tensor
      Index faiptKey1   ( model.Database().StorageKey("faip tensor 1") );
      Index faiptKey2   ( model.Database().StorageKey("faip tensor 2") );
      // TYPE: Array
      Index faipaKey1   ( model.Database().StorageKey("faip array 1") );
      Index faipaKey2   ( model.Database().StorageKey("faip array 2") );
      Index faipaKey3   ( model.Database().StorageKey("faip array 3") );
      // TYPE: Flagged Array
      Index faipfaKey1  ( model.Database().StorageKey("faip flagged array 1") );
      Index faipfaKey2  ( model.Database().StorageKey("faip flagged array 2") );
      Index faipfaKey3  ( model.Database().StorageKey("faip flagged array 3") );

      // ========================================================================
      // Some working variables
      // Scalar's
      ScalarVariable scalarV;
      ScalarVariable oneS               ( ANY, 1.       );
      ScalarVariable twoS               ( DIRICH, 2.    );
      ScalarVariable threeS             ( ROBIN, 3.     );
      // Vector's
      VectorVariable<3> vectorV;
      VectorVariable<3> fourV           ( NEUMANN, 4.   );
      VectorVariable<3> fiveV           ( ANY, 5.       );
      VectorVariable<3> sixV            ( ANY, 6.       );
      // Tensor's
      TensorVariable<3> tensorV;
      TensorVariable<3> sevenT          ( DIRICH, 7.    );
      TensorVariable<3> eightT          ( ROBIN, 8.     );
      TensorVariable<3> nineT           ( NEUMANN, 9.   );
      // Array's
      ArrayVariable ArrayV(28); // size of elemnt array variable 3
      ArrayVariable modelArray1         ( "model array 1",  model.Database(), 1.,   DIRICH  );
      ArrayVariable modelArray1_plain   ( "model array 1",  model.Database()                );
      ArrayVariable nodalArray1         ( "nodal array 1",  model.Database(), 1.,   DIRICH  );
      ArrayVariable nodalArray2         ( "nodal array 2",  model.Database(), 2.,   NEUMANN );
      ArrayVariable nodalArray3         ( "nodal array 3",  model.Database(), 3.,   ROBIN   );
      ArrayVariable nodalArray1_plain   ( "nodal array 1",  model.Database() );
      ArrayVariable nodalArray2_plain   ( "nodal array 2",  model.Database() );
      ArrayVariable nodalArray3_plain   ( "nodal array 3",  model.Database() );
      ArrayVariable elementArray1       ( "element array 1", model.Database(), 1.,  DIRICH  );
      ArrayVariable elementArray2       ( "element array 2", model.Database(), 2.,  NEUMANN );
      ArrayVariable elementArray3       ( "element array 3", model.Database(), 3.,  ROBIN   );
      ArrayVariable elementArray1_plain ( "element array 1", model.Database() );
      ArrayVariable elementArray2_plain ( "element array 2", model.Database() );
      ArrayVariable elementArray3_plain ( "element array 3", model.Database() );
      ArrayVariable eipArray1           ( "eip array 1",    model.Database(), 1.,   DIRICH  );
      ArrayVariable eipArray2           ( "eip array 2",    model.Database(), 2.,   NEUMANN );
      ArrayVariable eipArray3           ( "eip array 3",    model.Database(), 3.,   ROBIN   );
      ArrayVariable eipArray1_plain     ( "eip array 1",    model.Database() );
      ArrayVariable eipArray2_plain     ( "eip array 2",    model.Database() );
      ArrayVariable eipArray3_plain     ( "eip array 3",    model.Database() );
      ArrayVariable fvipArray1          ( "seip array 1",   model.Database(), 1.,   DIRICH  );
      ArrayVariable fvipArray2          ( "seip array 2",   model.Database(), 2.,   NEUMANN );
      ArrayVariable fvipArray3          ( "seip array 3",   model.Database(), 3.,   ROBIN   );
      ArrayVariable fvipArray1_plain    ( "seip array 1",   model.Database() );
      ArrayVariable fvipArray2_plain    ( "seip array 2",   model.Database() );
      ArrayVariable fvipArray3_plain    ( "seip array 3",   model.Database() );
      // Flagged Array's
      FlaggedArrayVariable FlaggedArrayV(3); // size of flagged array variable 3
      FlaggedArrayVariable modelFlaggedArray1         ( "model flagged array 1", model.Database(), 1.,    DIRICH  );
      FlaggedArrayVariable modelFlaggedArray1_plain   ( "model flagged array 1", model.Database() );
      FlaggedArrayVariable nodalFlaggedArray1         ( "nodal flagged array 1", model.Database(), 1.,    DIRICH  );
      FlaggedArrayVariable nodalFlaggedArray2         ( "nodal flagged array 2", model.Database(), 2.,    NEUMANN );
      FlaggedArrayVariable nodalFlaggedArray3         ( "nodal flagged array 3", model.Database(), 3.,    ROBIN   );
      FlaggedArrayVariable nodalFlaggedArray1_plain   ( "nodal flagged array 1", model.Database() );
      FlaggedArrayVariable nodalFlaggedArray2_plain   ( "nodal flagged array 2", model.Database() );
      FlaggedArrayVariable nodalFlaggedArray3_plain   ( "nodal flagged array 3", model.Database() );
      FlaggedArrayVariable elementFlaggedArray1       ( "element flagged array 1", model.Database(), 1.,  DIRICH  );
      FlaggedArrayVariable elementFlaggedArray2       ( "element flagged array 2", model.Database(), 2.,  NEUMANN );
      FlaggedArrayVariable elementFlaggedArray3       ( "element flagged array 3", model.Database(), 3.,  ROBIN   );
      FlaggedArrayVariable elementFlaggedArray1_plain ( "element flagged array 1", model.Database() );
      FlaggedArrayVariable elementFlaggedArray2_plain ( "element flagged array 2", model.Database() );
      FlaggedArrayVariable elementFlaggedArray3_plain ( "element flagged array 3", model.Database() );
      FlaggedArrayVariable eipFlaggedArray1           ( "eip flagged array 1", model.Database(), 1.,      DIRICH  );
      FlaggedArrayVariable eipFlaggedArray2           ( "eip flagged array 2", model.Database(), 2.,      NEUMANN );
      FlaggedArrayVariable eipFlaggedArray3           ( "eip flagged array 3", model.Database(), 3.,      ROBIN   );
      FlaggedArrayVariable eipFlaggedArray1_plain     ( "eip flagged array 1", model.Database() );
      FlaggedArrayVariable eipFlaggedArray2_plain     ( "eip flagged array 2", model.Database() );
      FlaggedArrayVariable eipFlaggedArray3_plain     ( "eip flagged array 3", model.Database() );
      FlaggedArrayVariable fvipFlaggedArray1          ( "seip flagged array 1", model.Database(), 1.,     DIRICH  );
      FlaggedArrayVariable fvipFlaggedArray2          ( "seip flagged array 2", model.Database(), 2.,     NEUMANN );
      FlaggedArrayVariable fvipFlaggedArray3          ( "seip flagged array 3", model.Database(), 3.,     ROBIN   );
      FlaggedArrayVariable fvipFlaggedArray1_plain    ( "seip flagged array 1", model.Database() );
      FlaggedArrayVariable fvipFlaggedArray2_plain    ( "seip flagged array 2", model.Database() );
      FlaggedArrayVariable fvipFlaggedArray3_plain    ( "seip flagged array 3", model.Database() );


      // ========================================================================
      // Region Nodal Scalar Average
      // TYPE: Scalar
      model.InputPropertyValue( "nodal variable 1", makeScalar( PLAIN, 1. ) );
      _test( model.Region("Model").Average("nodal variable 1") == 1. );
      _test( ( *model.Region("Model").NodesBegin() )->Read(nKey1) == 1. );

      // single face variable test
      // TYPE: Scalar
      model.Boundary("BOUNDARY1").InputPropertyValue( "face variable 1", makeScalar( PLAIN, 4. ) );
      _test( ( *model.Boundary("BOUNDARY1").ElementsBegin() )->Read(fKey1) == 4. );

      model.Database().WriteVariablesFile("CSMP-variables-output.txt");


      // ========================================================================
      // Nodal Scalar Variable Test Thorough
      ScalarVariable sv( PLAIN, 99999. );
      ScalarVariable sv1( PLAIN, 1. );
      model.InputPropertyValue( "nodal variable 1", sv1 );
      model.InputPropertyValue( "nodal variable 2", makeScalar( DIRICH, 99. ) );
      for( vector<Node<3>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != model.Region("Model").NodesEnd(); ++it )
        {
          sv = makeScalar( DIRICH, 9999. );
          (*it)->Read( nKey1, sv );
          _test( sv == sv1 );
        }
      
      // ========================================================================
      // Nodal Array Variable Test Thorough
      ArrayVariable av0(4, 10., DIRICH );
      for( vector<Node<3>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != model.Region("Model").NodesEnd(); ++it )
        (*it)->Store( naKey1, av0 );
      ArrayVariable av00( "nodal array 1", model.Database() );
      for( vector<Node<3>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != model.Region("Model").NodesEnd(); ++it )
        {
          // TYPE: Array
          (*it)->Read( naKey1, av00 );
          _test( av00 == av0 );
        }

      // ========================================================================
      // Element Array Variable Test Thorough
      av00 = ArrayVariable( "element array 2", model.Database(), 0. );
      model.InputPropertyValue( "element array 2", av00 );
      ArrayVariable av1( 22, 77., DIRICH );
      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
          // TYPE: Array
          (*it)->Read( eaKey2, av1 );
          _test( av1 == av00 );
        }


      // ========================================================================
      // Nodal Flagged Array Variable Test Thorough
      FlaggedArrayVariable fav0(1, 10., DIRICH );
      for( vector<Node<3>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != model.Region("Model").NodesEnd(); ++it )
        (*it)->Store( nfaKey1, fav0 );
      FlaggedArrayVariable fav00( "nodal flagged array 1", model.Database() );
      for( vector<Node<3>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != model.Region("Model").NodesEnd(); ++it )
        {
          // TYPE: Flagged Array
          (*it)->Read( nfaKey1, fav00 );
          _test( fav00 == fav0 );
        }

      // ========================================================================
      // Element Flagged Array Variable Test Thorough
      fav00 = FlaggedArrayVariable( "element flagged array 2", model.Database(), 0. );
      model.InputPropertyValue( "element flagged array 2", fav00 );
      FlaggedArrayVariable fav1( 2, 77., DIRICH );
      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
          // TYPE: Flagged Array
          (*it)->Read( efaKey2, fav1 );
          _test( fav1 == fav00 );
        }


      // ========================================================================
      //      Tests For all avaliable PLACEMENT's and TYPE's
      // ========================================================================

      // Initialization:

      // PLACEMENT: Model:
      // TYPE: Scalar
      model.InputPropertyValue( "model variable 1",         oneS );
      // TYPE: Array
      model.InputPropertyValue( "model array 1",            modelArray1 );
      // TYPE: Flagged Array
      model.InputPropertyValue( "model flagged array 1",    modelFlaggedArray1 );

      // PLACEMENT: Subdomain ( Region, Boundary )
      // TYPE: Vector
      model.InputPropertyValue( "region vector 1",          fourV );
      // TYPE: Tensor (placed on boundary)
      model.Boundary("BOUNDARY2").InputPropertyValue( "boundary tensor 1", sevenT );

      // PLACEMENT: Node
      // TYPE: Scalar
      model.InputPropertyValue( "nodal variable 1",         oneS );
      model.InputPropertyValue( "nodal variable 2",         twoS );
      // TYPE: Vector
      model.InputPropertyValue( "nodal vector 1",           fourV );
      model.InputPropertyValue( "nodal vector 2",           fiveV );
      // TYPE: Tensor
      model.InputPropertyValue( "nodal tensor 1",           sevenT );
      model.InputPropertyValue( "nodal tensor 2",           eightT );
      // TYPE: Array
      model.InputPropertyValue( "nodal array 1",            nodalArray1 );
      model.InputPropertyValue( "nodal array 2",            nodalArray2 );
      model.InputPropertyValue( "nodal array 3",            nodalArray3 );
      // TYPE: Flagged Array
      model.InputPropertyValue( "nodal flagged array 1",    nodalFlaggedArray1 );
      model.InputPropertyValue( "nodal flagged array 2",    nodalFlaggedArray2 );
      model.InputPropertyValue( "nodal flagged array 3",    nodalFlaggedArray3 );

      // PLACEMENT: Element
      // TYPE: Scalar
      model.InputPropertyValue( "element variable 1",       oneS );
      model.InputPropertyValue( "element variable 2",       twoS );
      // TYPE: Vector
      model.InputPropertyValue( "element vector 1",         fourV );
      model.InputPropertyValue( "element vector 2",         fiveV );
      // TYPE: Tensor
      model.InputPropertyValue( "element tensor 1",         sevenT );
      model.InputPropertyValue( "element tensor 2",         eightT );
      // TYPE: Array
      model.InputPropertyValue( "element array 1",          elementArray1 );
      model.InputPropertyValue( "element array 2",          elementArray2 );
      model.InputPropertyValue( "element array 3",          elementArray3 );
      // TYPE: Flagged Array
      model.InputPropertyValue( "element flagged array 1",  elementFlaggedArray1 );
      model.InputPropertyValue( "element flagged array 2",  elementFlaggedArray2 );
      model.InputPropertyValue( "element flagged array 3",  elementFlaggedArray3 );

      // PLACEMENT: Face
      // TYPE: Scalar
      model.InputPropertyValue( "face variable 1",          oneS );
      model.InputPropertyValue( "face variable 2",          twoS );
      // TYPE: Vector
      model.InputPropertyValue( "face vector 1",            fourV );
      model.InputPropertyValue( "face vector 2",            fiveV );
      // TYPE: Tensor
      model.InputPropertyValue( "face tensor 1",            sevenT );
      model.InputPropertyValue( "face tensor 2",            eightT );
      // TYPE: Array
      model.InputPropertyValue( "face array 1",             elementArray1 );
      model.InputPropertyValue( "face array 2",             elementArray2 );
      model.InputPropertyValue( "face array 3",             elementArray3 );
      // TYPE: Flagged Array
      model.InputPropertyValue( "face flagged array 1",     elementFlaggedArray1 );
      model.InputPropertyValue( "face flagged array 2",     elementFlaggedArray2 );
      model.InputPropertyValue( "face flagged array 3",     elementFlaggedArray3 );

      // PLACEMENT: Element Integration Point
      // TYPE: Scalar
      model.InputPropertyValue( "eip scalar 1",             oneS );
      model.InputPropertyValue( "eip scalar 2",             twoS );
      // TYPE: Vector
      model.InputPropertyValue( "eip vector 1",             fourV );
      model.InputPropertyValue( "eip vector 2",             fiveV );
      // TYPE: Tensor
      model.InputPropertyValue( "eip tensor 1",             sevenT );
      model.InputPropertyValue( "eip tensor 2",             eightT );
      // TYPE: Array
      model.InputPropertyValue( "eip array 1",              eipArray1 );
      model.InputPropertyValue( "eip array 2",              eipArray2 );
      model.InputPropertyValue( "eip array 3",              eipArray3 );
      // TYPE: Flagged Array
      model.InputPropertyValue( "eip flagged array 1",      eipFlaggedArray1 );
      model.InputPropertyValue( "eip flagged array 2",      eipFlaggedArray2 );
      model.InputPropertyValue( "eip flagged array 3",      eipFlaggedArray3 );
      
      // PLACEMENT: Face Integration Point
      // TYPE: Scalar
      model.InputPropertyValue( "fip scalar 1",             oneS );
      model.InputPropertyValue( "fip scalar 2",             twoS );
      // TYPE: Vector
      model.InputPropertyValue( "fip vector 1",             fourV );
      model.InputPropertyValue( "fip vector 2",             fiveV );
      // TYPE: Tensor
      model.InputPropertyValue( "fip tensor 1",             sevenT );
      model.InputPropertyValue( "fip tensor 2",             eightT );
      // TYPE: Array
      model.InputPropertyValue( "fip array 1",              eipArray1 );
      model.InputPropertyValue( "fip array 2",              eipArray2 );
      model.InputPropertyValue( "fip array 3",              eipArray3 );
      // TYPE: Flagged Array
      model.InputPropertyValue( "fip flagged array 1",      eipFlaggedArray1 );
      model.InputPropertyValue( "fip flagged array 2",      eipFlaggedArray2 );
      model.InputPropertyValue( "fip flagged array 3",      eipFlaggedArray3 );

      // ========================================================================
      // PLACEMENT: Model
      // TYPE: Scalar
      model.Read( mKey1, scalarV );
      _test( scalarV == oneS );
      // TYPE: Array
      model.Read( maKey1, modelArray1_plain );
      _test( modelArray1_plain == modelArray1 );
      // TYPE: Flagged Array
      model.Read( mfaKey1, modelFlaggedArray1_plain );
      _test( modelFlaggedArray1_plain == modelFlaggedArray1 );

      // ========================================================================
      // PLACEMENT: Subdomain ( Region, Boundary )
      // TYPE: Vector
      model.Region("Model").Read( rvKey1, vectorV );
      _test( vectorV == fourV );
      // TYPE: Tensor
      model.Boundary("BOUNDARY2").Read( btKey1, tensorV );
      _test( tensorV == sevenT );

      // ========================================================================
      // PLACEMENT: Node
      for( vector<Node<3>*>::const_iterator it( model.Region("Model").NodesBegin() ); it != model.Region("Model").NodesEnd(); ++it )
        {
          // TYPE: Scalar
          _test( (*it)->Read( nKey1 ) == oneS() );
          (*it)->Read( nKey1, scalarV );
          _test( scalarV == oneS );
          (*it)->Read( nKey2, scalarV );
          _test( scalarV == twoS );
          // TYPE: Vector
          (*it)->Read( nvKey1, vectorV );
          _test( vectorV == fourV );
          (*it)->Read( nvKey2, vectorV );
          _test( vectorV == fiveV );
          // TYPE: Tensor
          (*it)->Read( ntKey1, tensorV );
          _test( tensorV == sevenT );
          (*it)->Read( ntKey2, tensorV );
          _test( tensorV == eightT );
          // TYPE: Array
          (*it)->Read( naKey1, nodalArray1_plain );
          _test( nodalArray1_plain == nodalArray1 );
          (*it)->Read( naKey2, nodalArray2_plain );
          _test( nodalArray2_plain == nodalArray2 );
          (*it)->Read( naKey3, nodalArray3_plain );
          _test( nodalArray3_plain == nodalArray3 );
          // TYPE: Flagged Array
          (*it)->Read( nfaKey1, nodalFlaggedArray1_plain );
          _test( nodalFlaggedArray1_plain == nodalFlaggedArray1 );
          (*it)->Read( nfaKey2, nodalFlaggedArray2_plain );
          _test( nodalFlaggedArray2_plain == nodalFlaggedArray2 );
          (*it)->Read( nfaKey3, nodalFlaggedArray3_plain );
          _test( nodalFlaggedArray3_plain == nodalFlaggedArray3 );

        }

      // ========================================================================
      // PLACEMENT: Element & Element Integration Points

      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
          // ========================================================================
          // PLACEMENT: Element

          // TYPE: Scalar
          _test( (*it)->Read( eKey1 ) == oneS() );
          (*it)->Read( eKey1, scalarV );
          _test( scalarV == oneS );
          (*it)->Read( eKey2, scalarV );
          _test( scalarV == twoS );
          // TYPE: Vector
          (*it)->Read( evKey1, vectorV );
          _test( vectorV == fourV );
          (*it)->Read( evKey2, vectorV );
          _test( vectorV == fiveV );
          // TYPE: Tensor
          (*it)->Read( etKey1, tensorV );
          _test( tensorV == sevenT );
          (*it)->Read( etKey2, tensorV );
          _test( tensorV == eightT );
          // TYPE: Array
          (*it)->Read( indexCache.key1_, elementArray1_plain );
          _test( elementArray1_plain == elementArray1 );
          (*it)->Read( eaKey1, elementArray1_plain );
          _test( elementArray1_plain == elementArray1 );
          (*it)->Read( eaKey2, elementArray2_plain );
          _test( elementArray2_plain == elementArray2 );
          (*it)->Read( eaKey3, elementArray3_plain );
          _test( elementArray3_plain == elementArray3 );
          // TYPE: Flagged Array
          (*it)->Read( efaKey1, elementFlaggedArray1_plain );
          _test( elementFlaggedArray1_plain == elementFlaggedArray1 );
          (*it)->Read( efaKey2, elementFlaggedArray2_plain );
          _test( elementFlaggedArray2_plain == elementFlaggedArray2 );
          (*it)->Read( efaKey3, elementFlaggedArray3_plain );
          _test( elementFlaggedArray3_plain == elementFlaggedArray3 );

          // ========================================================================
          // PLACEMENT: Element Integration Point
          for( size_t i(0); i < (*it)->IntegrationPoints(); ++i )
            {
              // TYPE: Scalar
              (*it)->Read( i, eipKey1, scalarV );
              _test( scalarV == oneS );
              (*it)->Read( i, eipKey2, scalarV );
              _test( scalarV == twoS );
              // TYPE: Vecto
              (*it)->Read( i, eipvKey1, vectorV );
              _test( vectorV == fourV );
              (*it)->Read( i, eipvKey2, vectorV );
              _test( vectorV == fiveV );
              // TYPE: Tensor
              (*it)->Read( i, eiptKey1, tensorV );
              _test( tensorV == sevenT );
              (*it)->Read( i, eiptKey2, tensorV );
              _test( tensorV == eightT );
              // TYPE: Array
              (*it)->Read( i, indexCache.key2_, eipArray1_plain );
              _test( eipArray1_plain == eipArray1 );
              (*it)->Read( i, eipaKey1, eipArray1_plain );
              _test( eipArray1_plain == eipArray1 );
              (*it)->Read( i, eipaKey2, eipArray2_plain );
              _test( eipArray2_plain == eipArray2 );
              (*it)->Read( i, eipaKey3, eipArray3_plain );
              _test( eipArray3_plain == eipArray3 );
              // TYPE: Flagged Array
              (*it)->Read( i, eipfaKey1, eipFlaggedArray1_plain );
              _test( eipFlaggedArray1_plain == eipFlaggedArray1 );
              (*it)->Read( i, eipfaKey2, eipFlaggedArray2_plain );
              _test( eipFlaggedArray2_plain == eipFlaggedArray2 );
              (*it)->Read( i, eipfaKey3, eipFlaggedArray3_plain );
              _test( eipFlaggedArray3_plain == eipFlaggedArray3 );
            }

        }

      // ========================================================================
      // PLACEMENT: Face
      for( vector<Face<3>*>::const_iterator it( model.Boundary("BOUNDARY3").ElementsBegin() ); it != model.Boundary("BOUNDARY3").ElementsEnd(); ++it )
        {
        // TYPE: Scalar
        (*it)->Read( fKey1, scalarV );
        _test( scalarV == oneS );
        (*it)->Read( fKey2, scalarV );
        _test( scalarV == twoS );
        // TYPE: Vector
        (*it)->Read( fvKey1, vectorV );
        _test( vectorV == fourV );
        (*it)->Read( fvKey2, vectorV );
        _test( vectorV == fiveV );
        // TYPE: Tensor
        (*it)->Read( ftKey1, tensorV );
        _test( tensorV == sevenT );
        (*it)->Read( ftKey2, tensorV );
        _test( tensorV == eightT );
        // TYPE: Array
        (*it)->Read( faKey1, elementArray1_plain );
        _test( elementArray1_plain == elementArray1 );
        (*it)->Read( faKey2, elementArray2_plain );
        _test( elementArray2_plain == elementArray2 );
        (*it)->Read( faKey3, elementArray3_plain );
        _test( elementArray3_plain == elementArray3 );
        // TYPE: Flagged Array
        (*it)->Read( ffaKey1, elementFlaggedArray1_plain );
        _test( elementFlaggedArray1_plain == elementFlaggedArray1 );
        (*it)->Read( ffaKey2, elementFlaggedArray2_plain );
        _test( elementFlaggedArray2_plain == elementFlaggedArray2 );
        (*it)->Read( ffaKey3, elementFlaggedArray3_plain );
        _test( elementFlaggedArray3_plain == elementFlaggedArray3 );

        // ========================================================================
        // PLACEMENT: Face Integration Points
        for( size_t i(0); i < (*it)->IntegrationPoints(); ++i )
          {
          // TYPE: Scalar
          (*it)->Read( i, fipKey1, scalarV );
          _test( scalarV == oneS );
          (*it)->Read( i, fipKey2, scalarV );
          _test( scalarV == twoS );
          // TYPE: Vector
          (*it)->Read( i, fipvKey1, vectorV );
          _test( vectorV == fourV );
          (*it)->Read( i, fipvKey2, vectorV );
          _test( vectorV == fiveV );
          // TYPE: Tensor
          (*it)->Read( i, fiptKey1, tensorV );
          _test( tensorV == sevenT );
          (*it)->Read( i, fiptKey2, tensorV );
          _test( tensorV == eightT );
          // TYPE: Array
          (*it)->Read( i, fipaKey1, eipArray1_plain );
          _test( eipArray1_plain == eipArray1 );
          (*it)->Read( i, fipaKey2, eipArray2_plain );
          _test( eipArray2_plain == eipArray2 );
          (*it)->Read( i, fipaKey3, eipArray3_plain );
          _test( eipArray3_plain == eipArray3 );
          // TYPE: Flagged Array
          (*it)->Read( i, fipfaKey1, eipFlaggedArray1_plain );
          _test( eipFlaggedArray1_plain == eipFlaggedArray1 );
          (*it)->Read( i, fipfaKey2, eipFlaggedArray2_plain );
          _test( eipFlaggedArray2_plain == eipFlaggedArray2 );
          (*it)->Read( i, fipfaKey3, eipFlaggedArray3_plain );
          _test( eipFlaggedArray3_plain == eipFlaggedArray3 );
          }

        }


      // ========================================================================
      // Note:  Repeat previous tests for element's in order to check that
      //        creation of FV stencil's and variables placed on it's integration points
      //        didn't affect the other variables placed on element and it's integration point
      // Instantiation of Finite Volumes
      model.InstantiateFiniteVolumes();

      // Initialization of Variables Placed on FV Integration Points
      // PLACEMENT: Sector Integration Point
      // TYPE: Scalar Properties
      model.InputPropertyValue( "seip scalar 1", oneS );
      model.InputPropertyValue( "seip scalar 2", twoS );
      // TYPE: Vector Properties
      model.InputPropertyValue( "seip vector 1", fourV );
      model.InputPropertyValue( "seip vector 2", fiveV );
      model.InputPropertyValue( "seip tensor 1", sevenT );
      model.InputPropertyValue( "seip tensor 2", eightT );
      // TYPE: Array Properties
      model.InputPropertyValue( "seip array 1", fvipArray1 );
      model.InputPropertyValue( "seip array 2", fvipArray2 );
      model.InputPropertyValue( "seip array 3", fvipArray3 );
      // TYPE: Flagged Array Properties
      model.InputPropertyValue( "seip flagged array 1", fvipFlaggedArray1 );
      model.InputPropertyValue( "seip flagged array 2", fvipFlaggedArray2 );
      model.InputPropertyValue( "seip flagged array 3", fvipFlaggedArray3 );

      // PLACEMENT: Facet  Integration Point
      // TYPE: Scalar Properties
      model.InputPropertyValue( "faip scalar 1", oneS );
      model.InputPropertyValue( "faip scalar 2", twoS );
      // TYPE: Vector Properties
      model.InputPropertyValue( "faip vector 1", fourV );
      model.InputPropertyValue( "faip vector 2", fiveV );
      // TYPE: Tensor Properties
      model.InputPropertyValue( "faip tensor 1", sevenT );
      model.InputPropertyValue( "faip tensor 2", eightT );
      // TYPE: Array Properties
      model.InputPropertyValue( "faip array 1", fvipArray1 );
      model.InputPropertyValue( "faip array 2", fvipArray2 );
      model.InputPropertyValue( "faip array 3", fvipArray3 );
      // TYPE: Flagged Array Properties
      model.InputPropertyValue( "faip flagged array 1", fvipFlaggedArray1 );
      model.InputPropertyValue( "faip flagged array 2", fvipFlaggedArray2 );
      model.InputPropertyValue( "faip flagged array 3", fvipFlaggedArray3 );

      // ========================================================================
      // PLACEMENT: Element & Integration Points

      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
          // ========================================================================
          // PLACEMENT: Element

          // TYPE: Scalar
          _test( (*it)->Read( eKey1 ) == oneS() );
          (*it)->Read( eKey1, scalarV );
          _test( scalarV == oneS );
          (*it)->Read( eKey2, scalarV );
          _test( scalarV == twoS );
          // TYPE: Vector
          (*it)->Read( evKey1, vectorV );
          _test( vectorV == fourV );
          (*it)->Read( evKey2, vectorV );
          _test( vectorV == fiveV );
          // TYPE: Tensor
          (*it)->Read( etKey1, tensorV );
          _test( tensorV == sevenT );
          (*it)->Read( etKey2, tensorV );
          _test( tensorV == eightT );
          // TYPE: Array
          (*it)->Read( indexCache.key1_, elementArray1_plain );
          _test( elementArray1_plain == elementArray1 );
          (*it)->Read( eaKey1, elementArray1_plain );
          _test( elementArray1_plain == elementArray1 );
          (*it)->Read( eaKey2, elementArray2_plain );
          _test( elementArray2_plain == elementArray2 );
          (*it)->Read( eaKey3, elementArray3_plain );
          _test( elementArray3_plain == elementArray3 );
          // TYPE: Flagged Array
          (*it)->Read( efaKey1, elementFlaggedArray1_plain );
          _test( elementFlaggedArray1_plain == elementFlaggedArray1 );
          (*it)->Read( efaKey2, elementFlaggedArray2_plain );
          _test( elementFlaggedArray2_plain == elementFlaggedArray2 );
          (*it)->Read( efaKey3, elementFlaggedArray3_plain );
          _test( elementFlaggedArray3_plain == elementFlaggedArray3 );

          // ========================================================================
          // PLACEMENT: Element Integration Point
          for( size_t i(0); i < (*it)->IntegrationPoints(); ++i )
            {
              // TYPE: Scalar
              (*it)->Read( i, eipKey1, scalarV );
              _test( scalarV == oneS );
              (*it)->Read( i, eipKey2, scalarV );
              _test( scalarV == twoS );
              // TYPE: Vecto
              (*it)->Read( i, eipvKey1, vectorV );
              _test( vectorV == fourV );
              (*it)->Read( i, eipvKey2, vectorV );
              _test( vectorV == fiveV );
              // TYPE: Tensor
              (*it)->Read( i, eiptKey1, tensorV );
              _test( tensorV == sevenT );
              (*it)->Read( i, eiptKey2, tensorV );
              _test( tensorV == eightT );
              // TYPE: Array
              (*it)->Read( i, indexCache.key2_, eipArray1_plain );
              _test( eipArray1_plain == eipArray1 );
              (*it)->Read( i, eipaKey1, eipArray1_plain );
              _test( eipArray1_plain == eipArray1 );
              (*it)->Read( i, eipaKey2, eipArray2_plain );
              _test( eipArray2_plain == eipArray2 );
              (*it)->Read( i, eipaKey3, eipArray3_plain );
              _test( eipArray3_plain == eipArray3 );
              // TYPE: Flagged Array
              (*it)->Read( i, eipfaKey1, eipFlaggedArray1_plain );
              _test( eipFlaggedArray1_plain == eipFlaggedArray1 );
              (*it)->Read( i, eipfaKey2, eipFlaggedArray2_plain );
              _test( eipFlaggedArray2_plain == eipFlaggedArray2 );
              (*it)->Read( i, eipfaKey3, eipFlaggedArray3_plain );
              _test( eipFlaggedArray3_plain == eipFlaggedArray3 );
            }
            
          
          // ========================================================================
          // PLACEMENT: Sector Integration Point
          const size_t sectors( (*it)->Sectors() );
          assert( sectors != 0 );
          for( size_t j(0); j < sectors; ++j )
            for( size_t i(0); i < (*it)->IntegrationPointsPerSector(); ++i )
              {
                // TYPE: Scalar
                (*it)->Read( j, i, seipKey1, scalarV );
                _test( scalarV == oneS );
                (*it)->Read( j, i, seipKey2, scalarV );
                _test( scalarV == twoS );
                // TYPE: Vector
                (*it)->Read( j, i, seipvKey1, vectorV );
                _test( vectorV == fourV );
                (*it)->Read( j, i, seipvKey2, vectorV );
                _test( vectorV == fiveV );
                // TYPE: Tensor
                (*it)->Read( j, i, seiptKey1, tensorV );
                _test( tensorV == sevenT );
                (*it)->Read( j, i, seiptKey2, tensorV );
                _test( tensorV == eightT );
                // TYPE: Array
                (*it)->Read( j, i, seipaKey1, fvipArray1_plain );
                _test( fvipArray1_plain == fvipArray1 );
                (*it)->Read( j, i, seipaKey2, fvipArray2_plain );
                _test( fvipArray2_plain == fvipArray2 );
                (*it)->Read( j, i, seipaKey3, fvipArray3_plain );
                _test( fvipArray3_plain == fvipArray3 );
                // TYPE: Flagged Array
                (*it)->Read( j, i, seipfaKey1, fvipFlaggedArray1_plain );
                _test( fvipFlaggedArray1_plain == fvipFlaggedArray1 );
                (*it)->Read( j, i, seipfaKey2, fvipFlaggedArray2_plain );
                _test( fvipFlaggedArray2_plain == fvipFlaggedArray2 );
                (*it)->Read( j, i, seipfaKey3, fvipFlaggedArray3_plain );
                _test( fvipFlaggedArray3_plain == fvipFlaggedArray3 );
              }

            // ========================================================================
            // PLACEMENT: Facet Integration Point
            const size_t facets( (*it)->Facets() );
            assert( facets != 0 );
            for( size_t j(0); j < facets; ++j )
              for( size_t i(0); i < (*it)->IntegrationPointsPerFacet(); ++i )
                {
                  // TYPE: Scalar
                  (*it)->Read( j, i, faipKey1, scalarV );
                  _test( scalarV == oneS );
                  (*it)->Read( j, i, faipKey2, scalarV );
                  _test( scalarV == twoS );
                  // TYPE: Vector
                  (*it)->Read( j, i, faipvKey1, vectorV );
                  _test( vectorV == fourV );
                  (*it)->Read( j, i, faipvKey2, vectorV );
                  _test( vectorV == fiveV );
                  // TYPE: Tensor
                  (*it)->Read( j, i, faiptKey1, tensorV );
                  _test( tensorV == sevenT );
                  (*it)->Read( j, i, faiptKey2, tensorV );
                  _test( tensorV == eightT );
                  // TYPE: Array
                  (*it)->Read( j, i, faipaKey1, fvipArray1_plain );
                  _test( fvipArray1_plain == fvipArray1 );
                  (*it)->Read( j, i, faipaKey2, fvipArray2_plain );
                  _test( fvipArray2_plain == fvipArray2 );
                  (*it)->Read( j, i, faipaKey3, fvipArray3_plain );
                  _test( fvipArray3_plain == fvipArray3 );
                  // TYPE: Flagged Array
                  (*it)->Read( j, i, faipfaKey1, fvipFlaggedArray1_plain );
                  _test( fvipFlaggedArray1_plain == fvipFlaggedArray1 );
                  (*it)->Read( j, i, faipfaKey2, fvipFlaggedArray2_plain );
                  _test( fvipFlaggedArray2_plain == fvipFlaggedArray2 );
                  (*it)->Read( j, i, faipfaKey3, fvipFlaggedArray3_plain );
                  _test( fvipFlaggedArray3_plain == fvipFlaggedArray3 );
                }
              
        }
       
      IndexTrackerTestStruct* scopedIndices = new IndexTrackerTestStruct;
      *scopedIndices = indexCache;
      
      // ========================================================================
      // Deleting Property
      model.DeleteProperty("element variable 1");
      model.DeleteProperty("element vector 1");
      model.DeleteProperty("element tensor 1");
      model.DeleteProperty("element array 1");
      model.DeleteProperty("element flagged array 1");
      model.DeleteProperty("eip scalar 1");
      model.DeleteProperty("eip vector 1");
      model.DeleteProperty("eip tensor 1");
      model.DeleteProperty("eip array 1");
      model.DeleteProperty("eip flagged array 1");
      model.DeleteProperty("seip scalar 1");
      model.DeleteProperty("seip vector 1");
      model.DeleteProperty("seip tensor 1");
      model.DeleteProperty("seip array 1");
      model.DeleteProperty("seip flagged array 1");
      model.DeleteProperty("faip scalar 1");
      model.DeleteProperty("faip vector 1");
      model.DeleteProperty("faip tensor 1");
      model.DeleteProperty("faip array 1");
      model.DeleteProperty("faip flagged array 1");

      // ========================================================================
      // Note:  Repeat previous tests for element's in order to check that
      //        deleting some of the properties in element
      //        didn't affect the other variables placed on element and it's integration point

      // ========================================================================
      // PLACEMENT: Element & Integration Points
      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
        // ========================================================================
        // PLACEMENT: Element
        // TYPE: Scalar
        (*it)->Read( eKey2, scalarV );
        _test( scalarV == twoS );
        // TYPE: Vector
        (*it)->Read( evKey2, vectorV );
        _test( vectorV == fiveV );
        // TYPE: Tensor
        (*it)->Read( etKey2, tensorV );
        _test( tensorV == eightT );
        // TYPE: Array
        (*it)->Read( eaKey2, elementArray2_plain );
        _test( elementArray2_plain == elementArray2 );
        // TYPE: Flagged Array
        (*it)->Read( efaKey2, elementFlaggedArray2_plain );
        _test( elementFlaggedArray2_plain == elementFlaggedArray2 );

        // ========================================================================
        // PLACEMENT: Element Integration Point
        for( size_t i(0); i < (*it)->IntegrationPoints(); ++i )
          {
          // TYPE: Scalar
          (*it)->Read( i, eipKey2, scalarV );
          _test( scalarV == twoS );
          // TYPE: Vector
          (*it)->Read( i, eipvKey2, vectorV );
          _test( vectorV == fiveV );
          // TYPE: Tensor
          (*it)->Read( i, eiptKey2, tensorV );
          _test( tensorV == eightT );
          // TYPE: Array
          (*it)->Read( i, eipaKey2, eipArray2_plain );
          _test( eipArray2_plain == eipArray2 );
          (*it)->Read( i, eipaKey3, eipArray3_plain );
          _test( eipArray3_plain == eipArray3 );
          // TYPE: Flagged Array
          (*it)->Read( i, eipfaKey2, eipFlaggedArray2_plain );
          _test( eipFlaggedArray2_plain == eipFlaggedArray2 );
          (*it)->Read( i, eipfaKey3, eipFlaggedArray3_plain );
          _test( eipFlaggedArray3_plain == eipFlaggedArray3 );
          }


        // ========================================================================
        // PLACEMENT: Sector Integration Point
        const size_t sectors( (*it)->Sectors() );
        assert( sectors != 0 );
        for( size_t j(0); j < sectors; ++j )
          for( size_t i(0); i < (*it)->IntegrationPointsPerSector(); ++i )
            {
            // TYPE: Scalar
            (*it)->Read( j, i, seipKey2, scalarV );
            _test( scalarV == twoS );
            // TYPE: Vector
            (*it)->Read( j, i, seipvKey2, vectorV );
            _test( vectorV == fiveV );
            // TYPE: Tensor
            (*it)->Read( j, i, seiptKey2, tensorV );
            _test( tensorV == eightT );
            // TYPE: Array
            (*it)->Read( j, i, seipaKey2, fvipArray2_plain );
            _test( fvipArray2_plain == fvipArray2 );
            (*it)->Read( j, i, seipaKey3, fvipArray3_plain );
            _test( fvipArray3_plain == fvipArray3 );
            // TYPE: Flagged Array
            (*it)->Read( j, i, seipfaKey2, fvipFlaggedArray2_plain );
            _test( fvipFlaggedArray2_plain == fvipFlaggedArray2 );
            (*it)->Read( j, i, seipfaKey3, fvipFlaggedArray3_plain );
            _test( fvipFlaggedArray3_plain == fvipFlaggedArray3 );
            }

          // ========================================================================
          // PLACEMENT: Facet Integration Points
          const size_t facets( (*it)->Facets() );
          assert( facets != 0 );
          for( size_t j(0); j < facets; ++j )
            for( size_t i(0); i < (*it)->IntegrationPointsPerFacet(); ++i )
              {
              // TYPE: Scalar
              (*it)->Read( j, i, faipKey2, scalarV );
              _test( scalarV == twoS );
              // TYPE: Vector
              (*it)->Read( j, i, faipvKey2, vectorV );
              _test( vectorV == fiveV );
              // TYPE: Tensor
              (*it)->Read( j, i, faiptKey2, tensorV );
              _test( tensorV == eightT );
              // TYPE: Array
              (*it)->Read( j, i, faipaKey2, fvipArray2_plain );
              _test( fvipArray2_plain == fvipArray2 );
              (*it)->Read( j, i, faipaKey3, fvipArray3_plain );
              _test( fvipArray3_plain == fvipArray3 );
              // TYPE: Flagged Array
              (*it)->Read( j, i, faipfaKey2, fvipFlaggedArray2_plain );
              _test( fvipFlaggedArray2_plain == fvipFlaggedArray2 );
              (*it)->Read( j, i, faipfaKey3, fvipFlaggedArray3_plain );
              _test( fvipFlaggedArray3_plain == fvipFlaggedArray3 );
              }

        }

      delete scopedIndices;

      // ========================================================================
      // Note:  Repeat previous tests for element's in order to check that
      //        creating a new variables
      //        didn't affect the other variables placed on element and it's integration point

      // ========================================================================
      // Adding a Property

      // PLACEMENT: Element
      // TYPE: Scalar
      Index neweKey3  = model.CreateProperty( "new element scalar 3", "X", SCALAR, ELEMENT );
      model.InputPropertyValue( "new element scalar 3", threeS );
      // TYPE: Vector
      Index newevKey3 = model.CreateProperty( "new element vector 3", "X", VECTOR, ELEMENT );
      model.InputPropertyValue( "new element vector 3", sixV );
      // TYPE: Tensor
      Index newetKey3 = model.CreateProperty( "new element tensor 3", "X", TENSOR, ELEMENT );
      model.InputPropertyValue( "new element tensor 3", nineT );
      // TYPE: Array
      Index neweaKey3 = model.CreateProperty( "new element array 3", "X", ARRAY, ELEMENT, 28 );
      model.InputPropertyValue( "new element array 3", elementArray3 );
      // TYPE: Flagged Array
      Index newefaKey3 = model.CreateProperty( "new element flagged array 3", "X", FLAGGEDARRAY, ELEMENT, 3 );
      model.InputPropertyValue( "new element flagged array 3", elementFlaggedArray3 );

      // PLACEMENT: Element Integration Point
      // TYPE: Scalar
      Index neweipKey3  = model.CreateProperty( "new eip scalar 3", "X", SCALAR, ELEMENT_INTEGRATION_POINT );
      model.InputPropertyValue( "new eip scalar 3", threeS );
      // TYPE: Vector
      Index neweipvKey3 = model.CreateProperty( "new eip vector 3", "X", VECTOR, ELEMENT_INTEGRATION_POINT );
      model.InputPropertyValue( "new eip vector 3", sixV );
      // TYPE: Tensor
      Index neweiptKey3 = model.CreateProperty( "new eip tensor 3", "X", TENSOR, ELEMENT_INTEGRATION_POINT );
      model.InputPropertyValue( "new eip tensor 3", nineT );
      // TYPE: Array
      Index neweipaKey3 = model.CreateProperty( "new eip array 3", "X", ARRAY, ELEMENT_INTEGRATION_POINT, 28 );
      model.InputPropertyValue( "new eip array 3", elementArray3 );
      // TYPE: Flagged Array
      Index neweipfaKey3 = model.CreateProperty( "new eip flagged array 3", "X", FLAGGEDARRAY, ELEMENT_INTEGRATION_POINT, 3 );
      model.InputPropertyValue( "new eip flagged array 3", elementFlaggedArray3 );

      // PLACEMENT: Sector Integration Point
      // TYPE: Scalar
      Index newseipKey3  = model.CreateProperty( "new seip scalar 3", "X", SCALAR, SECTOR_INTEGRATION_POINT );
      model.InputPropertyValue( "new seip scalar 3", threeS );
      // TYPE: Vector
      Index newseipvKey3 = model.CreateProperty( "new seip vector 3", "X", VECTOR, SECTOR_INTEGRATION_POINT );
      model.InputPropertyValue( "new seip vector 3", sixV );
      // TYPE: Tensor
      Index newseiptKey3 = model.CreateProperty( "new seip tensor 3", "X", TENSOR, SECTOR_INTEGRATION_POINT );
      model.InputPropertyValue( "new seip tensor 3", nineT );
      // TYPE: Array
      Index newseipaKey3 = model.CreateProperty( "new seip array 3", "X", ARRAY, SECTOR_INTEGRATION_POINT, 28 );
      model.InputPropertyValue( "new seip array 3", elementArray3 );
      // TYPE: Flagged Array
      Index newseipfaKey3 = model.CreateProperty( "new seip flagged array 3", "X", FLAGGEDARRAY, SECTOR_INTEGRATION_POINT, 3 );
      model.InputPropertyValue( "new seip flagged array 3", elementFlaggedArray3 );

      // PLACEMENT: Facet Integration Point
      // TYPE: Scalar
      Index newfaipKey3  = model.CreateProperty( "new faip scalar 3", "X", SCALAR, FACET_INTEGRATION_POINT );
      model.InputPropertyValue( "new faip scalar 3", threeS );
      // TYPE: Vector
      Index newfaipvKey3 = model.CreateProperty( "new faip vector 3", "X", VECTOR, FACET_INTEGRATION_POINT );
      model.InputPropertyValue( "new faip vector 3", sixV );
      // TYPE: Tensor
      Index newfaiptKey3 = model.CreateProperty( "new faip tensor 3", "X", TENSOR, FACET_INTEGRATION_POINT );
      model.InputPropertyValue( "new faip tensor 3", nineT );
      // TYPE: Array
      Index newfaipaKey3 = model.CreateProperty( "new faip array 3", "X", ARRAY, FACET_INTEGRATION_POINT, 28 );
      model.InputPropertyValue( "new faip array 3", elementArray3 );
      // TYPE: Flagged Array
      Index newfaipfaKey3 = model.CreateProperty( "new faip flagged array 3", "X", FLAGGEDARRAY, FACET_INTEGRATION_POINT, 3 );
      model.InputPropertyValue( "new faip flagged array 3", elementFlaggedArray3 );

      // ========================================================================
      // PLACEMENT: Element & Integration Points
      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
        // TYPE: Scalar
        (*it)->Read( eKey2, scalarV );
        _test( scalarV == twoS );
        (*it)->Read( neweKey3, scalarV );
        _test( scalarV == threeS );
        // TYPE: Vector
        (*it)->Read( evKey2, vectorV );
        _test( vectorV == fiveV );
        (*it)->Read( newevKey3, vectorV );
        _test( vectorV == sixV );
        // TYPE: Tensor
        (*it)->Read( etKey2, tensorV );
        _test( tensorV == eightT );
        (*it)->Read( newetKey3, tensorV );
        _test( tensorV == nineT );
        // TYPE: Array
        (*it)->Read( eaKey2, elementArray2_plain );
        _test( elementArray2_plain == elementArray2 );
        (*it)->Read( neweaKey3, ArrayV );
        _test( ArrayV == elementArray3 );
        // TYPE: Flagged Array
        (*it)->Read( efaKey2, elementFlaggedArray2_plain );
        _test( elementFlaggedArray2_plain == elementFlaggedArray2 );
        (*it)->Read( newefaKey3, FlaggedArrayV );
        _test( FlaggedArrayV == elementFlaggedArray3 );

        // ========================================================================
        // PLACEMENT: Element Integration Point
        for( size_t i(0); i < (*it)->IntegrationPoints(); ++i )
          {
          // TYPE: Scalar
          (*it)->Read( i, eipKey2, scalarV );
          _test( scalarV == twoS );
          (*it)->Read( i, neweipKey3, scalarV );
          _test( scalarV == threeS );
          // TYPE: Vector
          (*it)->Read( i, eipvKey2, vectorV );
          _test( vectorV == fiveV );
          (*it)->Read( i, neweipvKey3, vectorV );
          _test( vectorV == sixV );
          // TYPE: Tensor
          (*it)->Read( i, eiptKey2, tensorV );
          _test( tensorV == eightT );
          (*it)->Read( i, neweiptKey3, tensorV );
          _test( tensorV == nineT );
          // TYPE: Array
          (*it)->Read( i, eipaKey2, eipArray2_plain );
          _test( eipArray2_plain == eipArray2 );
          (*it)->Read( i, eipaKey3, eipArray3_plain );
          _test( eipArray3_plain == eipArray3 );
          (*it)->Read( i, neweipaKey3, ArrayV );
          _test( ArrayV == elementArray3 );
          // TYPE: Flagged Array
          (*it)->Read( i, eipfaKey2, eipFlaggedArray2_plain );
          _test( eipFlaggedArray2_plain == eipFlaggedArray2 );
          (*it)->Read( i, eipfaKey3, eipFlaggedArray3_plain );
          _test( eipFlaggedArray3_plain == eipFlaggedArray3 );
          (*it)->Read( i, neweipfaKey3, FlaggedArrayV );
          _test( FlaggedArrayV == elementFlaggedArray3 );
          }


        // ========================================================================
        // PLACEMENT: Sector Integration Points
        const size_t sectors( (*it)->Sectors() );
        assert( sectors != 0 );
        for( size_t j(0); j < sectors; ++j )
          for( size_t i(0); i < (*it)->IntegrationPointsPerSector(); ++i )
            {
            // TYPE: Scalar
            (*it)->Read( j, i, seipKey2, scalarV );
            _test( scalarV == twoS );
            (*it)->Read( j, i, newseipKey3, scalarV );
            _test( scalarV == threeS );
            // TYPE: Vector
            (*it)->Read( j, i, seipvKey2, vectorV );
            _test( vectorV == fiveV );
            (*it)->Read( j, i, newseipvKey3, vectorV );
            _test( vectorV == sixV );
            // TYPE: Tensor
            (*it)->Read( j, i, seiptKey2, tensorV );
            _test( tensorV == eightT );
            (*it)->Read( j, i, newseiptKey3, tensorV );
            _test( tensorV == nineT );
            // TYPE: Array
            (*it)->Read( j, i, seipaKey2, fvipArray2_plain );
            _test( fvipArray2_plain == fvipArray2 );
            (*it)->Read( j, i, seipaKey3, fvipArray3_plain );
            _test( fvipArray3_plain == fvipArray3 );
            (*it)->Read( j, i, newseipaKey3, ArrayV );
            _test( ArrayV == elementArray3 );
            // TYPE: Flagged Array
            (*it)->Read( j, i, seipfaKey2, fvipFlaggedArray2_plain );
            _test( fvipFlaggedArray2_plain == fvipFlaggedArray2 );
            (*it)->Read( j, i, seipfaKey3, fvipFlaggedArray3_plain );
            _test( fvipFlaggedArray3_plain == fvipFlaggedArray3 );
            (*it)->Read( j, i, newseipfaKey3, FlaggedArrayV );
            _test( FlaggedArrayV == elementFlaggedArray3 );
            }

          // ========================================================================
          // PLACEMENT: Facet Integration Point
          const size_t facets( (*it)->Facets() );
          assert( facets != 0 );
          for( size_t j(0); j < facets; ++j )
            for( size_t i(0); i < (*it)->IntegrationPointsPerFacet(); ++i )
              {
              // TYPE: Scalar
              (*it)->Read( j, i, faipKey2, scalarV );
              _test( scalarV == twoS );
              (*it)->Read( j, i, newfaipKey3, scalarV );
              _test( scalarV == threeS );
              // TYPE: Vector
              (*it)->Read( j, i, faipvKey2, vectorV );
              _test( vectorV == fiveV );
              (*it)->Read( j, i, newfaipvKey3, vectorV );
              _test( vectorV == sixV );
              // TYPE: Tensor
              (*it)->Read( j, i, faiptKey2, tensorV );
              _test( tensorV == eightT );
              (*it)->Read( j, i, newfaiptKey3, tensorV );
              _test( tensorV == nineT );
              // TYPE: Array
              (*it)->Read( j, i, faipaKey2, fvipArray2_plain );
              _test( fvipArray2_plain == fvipArray2 );
              (*it)->Read( j, i, faipaKey3, fvipArray3_plain );
              _test( fvipArray3_plain == fvipArray3 );
              (*it)->Read( j, i, newfaipaKey3, ArrayV );
              _test( ArrayV == elementArray3 );
              // TYPE: Flagged Array
              (*it)->Read( j, i, faipfaKey2, fvipFlaggedArray2_plain );
              _test( fvipFlaggedArray2_plain == fvipFlaggedArray2 );
              (*it)->Read( j, i, faipfaKey3, fvipFlaggedArray3_plain );
              _test( fvipFlaggedArray3_plain == fvipFlaggedArray3 );
              (*it)->Read( j, i, newfaipfaKey3, FlaggedArrayV );
              _test( FlaggedArrayV == elementFlaggedArray3 );
              }
         }



      // ========================================================================
      // Status tests

      // PLACEMENT: Element
      for( vector<Element<3>*>::const_iterator it( model.Region("Model").ElementsBegin() ); it != model.Region("Model").ElementsEnd(); ++it )
        {
        // TYPE: Scalar
        _test( (*it)->Status(eKey2)     == twoS.Flag() );
        _test( (*it)->Status(neweKey3)  == threeS.Flag() );
        for( size_t d(0); d<3; ++d )
          {
          // TYPE: Vector
          _test( (*it)->Status(evKey2,d)    == fiveV.Flag(d) );
          _test( (*it)->Status(newevKey3,d) == sixV.Flag(d) );
          // TYPE: Tensor
          _test( (*it)->Status(etKey2,d)    == eightT.Flag(d) );
          _test( (*it)->Status(newetKey3,d) == nineT.Flag(d) );
          }
        // TYPE: Array
        _test( (*it)->Status(eaKey2) == elementArray2.Flag() );
        _test( (*it)->Status(neweaKey3) == elementArray3.Flag() );
        // TYPE: Flagged Array
        for( size_t d(0); d<elementFlaggedArray2.Size(); ++d )
            _test( (*it)->Status(efaKey2,d) == elementFlaggedArray2.Flag(d) );
        for( size_t d(0); d<elementFlaggedArray3.Size(); ++d )
            _test( (*it)->Status(newefaKey3,d) == elementFlaggedArray3.Flag(d) );

        // ========================================================================
        // PLACEMENT: Element Integration Point
        for( size_t i(0); i < (*it)->IntegrationPoints(); ++i )
          {
          // TYPE: Scalar
          _test( (*it)->Status(i,eipKey2)       == twoS.Flag() );
          _test( (*it)->Status(i,neweipKey3)    == threeS.Flag() );
          for( size_t d(0); d<3; ++d )
            {
            // TYPE: Vector
            _test( (*it)->Status(i,eipvKey2,d)      == fiveV.Flag(d) );
            _test( (*it)->Status(i,neweipvKey3,d)   == sixV.Flag(d) );
            // TYPE: Tensor
            _test( (*it)->Status(i,eiptKey2,d)      == eightT.Flag(d) );
            _test( (*it)->Status(i,neweiptKey3,d)   == nineT.Flag(d) );
            }
          // TYPE: Array
          _test( (*it)->Status(i,eipaKey2)      == eipArray2.Flag() );
          _test( (*it)->Status(i,eipaKey3)      == eipArray3.Flag() );
          _test( (*it)->Status(i,neweipaKey3)   == elementArray3.Flag() );
          // TYPE: Flagged Array
          for( size_t d(0); d<eipFlaggedArray3.Size(); ++d )
              _test( (*it)->Status(i,eipfaKey3,d) == eipFlaggedArray3.Flag(d) );
          for( size_t d(0); d<elementFlaggedArray3.Size(); ++d )
              _test( (*it)->Status(i,neweipfaKey3,d) == elementFlaggedArray3.Flag(d) );
          }


        // ========================================================================
        // PLACEMENT: Sector Integration Point
        const size_t sectors( (*it)->Sectors() );
        assert( sectors != 0 );
        for( size_t j(0); j < sectors; ++j )
          for( size_t i(0); i < (*it)->IntegrationPointsPerSector(); ++i )
            {
            // TYPE: Scalar
            _test( (*it)->Status(j,i,seipKey2)      == twoS.Flag() );
            _test( (*it)->Status(j,i,newseipKey3)   == threeS.Flag() );
            for( size_t d(0); d<3; ++d )
              {
              // TYPE: Vector
              _test( (*it)->Status(j,i,seipvKey2,d)     == fiveV.Flag(d) );
              _test( (*it)->Status(j,i,newseipvKey3,d)  == sixV.Flag(d) );
              // TYPE: Tensor
              _test( (*it)->Status(j,i,seiptKey2,d)     == eightT.Flag(d) );
              _test( (*it)->Status(j,i,newseiptKey3,d)  == nineT.Flag(d) );
              }
            // TYPE: Array
            _test( (*it)->Status(j,i,seipaKey2)     == fvipArray2.Flag() );
            _test( (*it)->Status(j,i,seipaKey3)     == fvipArray3.Flag() );
            _test( (*it)->Status(j,i,newseipaKey3)  == elementArray3.Flag() );
            // TYPE: Flagged Array
            for( size_t d(0); d<fvipFlaggedArray2.Size(); ++d )
                _test( (*it)->Status(j,i,seipfaKey2,d) == fvipFlaggedArray2.Flag(d) );
            for( size_t d(0); d<fvipFlaggedArray3.Size(); ++d )
                _test( (*it)->Status(j,i,seipfaKey3,d) == fvipFlaggedArray3.Flag(d) );
            for( size_t d(0); d<elementFlaggedArray3.Size(); ++d )
                _test( (*it)->Status(j,i,newseipfaKey3,d) == elementFlaggedArray3.Flag(d) );
            }

          // ========================================================================
          // PLACEMENT: Facet Integration Point
          const size_t facets( (*it)->Facets() );
          assert( facets != 0 );
          for( size_t j(0); j < facets; ++j )
            for( size_t i(0); i < (*it)->IntegrationPointsPerFacet(); ++i )
              {
              // TYPE: Scalar
              _test( (*it)->Status(j,i,faipKey2)    == twoS.Flag() );
              _test( (*it)->Status(j,i,newfaipKey3) == threeS.Flag() );
              for( size_t d(0); d<3; ++d )
                {
                // TYPE: Vector
                _test( (*it)->Status(j,i,faipvKey2,d)       == fiveV.Flag(d) );
                _test( (*it)->Status(j,i,newfaipvKey3,d)    == sixV.Flag(d) );
                // TYPE: Tensor
                _test( (*it)->Status(j,i,faiptKey2,d)       == eightT.Flag(d) );
                _test( (*it)->Status(j,i,newfaiptKey3,d)    == nineT.Flag(d) );
                }
              // TYPE: Array
              _test( (*it)->Status(j,i,faipaKey2) == fvipArray2.Flag() );
              _test( (*it)->Status(j,i,faipaKey3) == fvipArray3.Flag() );
              _test( (*it)->Status(j,i,newfaipaKey3)  == elementArray3.Flag() );
              // TYPE: Flagged Array
              for( size_t d(0); d<fvipFlaggedArray2.Size(); ++d )
                  _test( (*it)->Status(j,i,faipfaKey2,d) == fvipFlaggedArray2.Flag(d) );
              for( size_t d(0); d<fvipFlaggedArray3.Size(); ++d )
                  _test( (*it)->Status(j,i,faipfaKey3,d) == fvipFlaggedArray3.Flag(d) );
              for( size_t d(0); d<elementFlaggedArray3.Size(); ++d )
                  _test( (*it)->Status(j,i,newfaipfaKey3,d) == elementFlaggedArray3.Flag(d) );

              }
        }
      
    } // end runModel

} // csmp
