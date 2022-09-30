#include "INDEXandVariables_Test.h"
#include "Region.h"
#include "Boundary.h"
#include "ModelTopology.h"
#include "vsetMakers.h"

#include "NodeCenteredFiniteVolumeTransport.h"

using namespace std;

namespace csmp {

/**
     SPLIT22_BASIC   is the test model model name.
     
     The test model should contain Boundary and SplitBoundary objects
*/
INDEXandVariables_Test::INDEXandVariables_Test()
  {
    this->setName("INDEXandVariables_Test");
  }
  
  
struct IndexTrackerTestStruct
  {
    csmp::Index key1_, key2_;
  };




void INDEXandVariables_Test::run()
  {
    // Run Test for 3D Model constructed by ANSYS mesh reader
	  string variables_filename = (string)(this->getName() + ".txt");
    VSet<2U>       vset;
    ModelTopology  model_topo = test_Create_BoundarySplitBoundaryPatch( vset );
    const bool treat_domains_as_regions_and_use_regions_file_if_any{ false };
    Model m0( model_topo, vset, variables_filename.c_str(), treat_domains_as_regions_and_use_regions_file_if_any );
    m0.OutputToBinaryFile("INDEXandVariables_Test_BinaryModel");
    runModel(m0);

    // Run Test for 3D Model loaded from CSMP++ binary format
    Model<2U> m1( string("INDEXandVariables_Test_BinaryModel") );
    runModel(m1);
  }
  
  
  
/**
 @fn  void INDEXandVariables_Test::run()

 @todo (2-F) Test for InterFace and SplitBoundary properties

 @author  P. Lang
 @date  9/24/2012
 */
void INDEXandVariables_Test::runModel( Model<2>& model )
  {
      const int dim{ 2U };
      
      // Indices
      IndexTrackerTestStruct indexCache;
      // ========================================================================
      // PLACEMENT: Model:
      // TYPE: Scalar
      INDEX<SCALAR,MODEL>  mKey1( model.Database().StorageKey("model variable 1") );
      // TYPE: Array
      INDEX<ARRAY,MODEL>  maKey1( model.Database().StorageKey("model array 1") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,MODEL>  mfaKey1( model.Database().StorageKey("model flagged array 1") );

      // ========================================================================
      // PLACEMENT: Subdomain ( Region, Boundary ):
      // TYPE: Vector property
      INDEX<VECTOR,REGION>  rvKey1( model.Database().StorageKey("region vector 1") );
      // Tensor property
      INDEX<TENSOR,BOUNDARY>  btKey1( model.Database().StorageKey("boundary tensor 1") );

      // ========================================================================
      // PLACEMENT: Node:
      // TYPE: Scalar
      INDEX<SCALAR,NODE>  nKey1       ( model.Database().StorageKey("nodal variable 1") );
      INDEX<SCALAR,NODE>  nKey2       ( model.Database().StorageKey("nodal variable 2") );
      // TYPE: Vector
      INDEX<VECTOR,NODE>  nvKey1      ( model.Database().StorageKey("nodal vector 1") );
      INDEX<VECTOR,NODE>  nvKey2      ( model.Database().StorageKey("nodal vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,NODE>  ntKey1      ( model.Database().StorageKey("nodal tensor 1") );
      INDEX<TENSOR,NODE>  ntKey2      ( model.Database().StorageKey("nodal tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,NODE>  naKey1      ( model.Database().StorageKey("nodal array 1") );
      INDEX<ARRAY,NODE>  naKey2      ( model.Database().StorageKey("nodal array 2") );
      INDEX<ARRAY,NODE>  naKey3      ( model.Database().StorageKey("nodal array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,NODE>  nfaKey1     ( model.Database().StorageKey("nodal flagged array 1") );
      INDEX<FLAGGEDARRAY,NODE>  nfaKey2     ( model.Database().StorageKey("nodal flagged array 2") );
      INDEX<FLAGGEDARRAY,NODE>  nfaKey3     ( model.Database().StorageKey("nodal flagged array 3") );

      // ========================================================================
      // PLACEMENt: Element
      // TYPE: Scalar
      INDEX<SCALAR,ELEMENT>  eKey1       ( model.Database().StorageKey("element variable 1") );
      INDEX<SCALAR,ELEMENT>  eKey2       ( model.Database().StorageKey("element variable 2") );
      // TYPE: Vector
      INDEX<VECTOR,ELEMENT>  evKey1      ( model.Database().StorageKey("element vector 1") );
      INDEX<VECTOR,ELEMENT>  evKey2      ( model.Database().StorageKey("element vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,ELEMENT>  etKey1      ( model.Database().StorageKey("element tensor 1") );
      INDEX<TENSOR,ELEMENT>  etKey2      ( model.Database().StorageKey("element tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,ELEMENT>  eaKey1      ( model.Database().StorageKey("element array 1") );
      INDEX<ARRAY,ELEMENT>  eaKey2      ( model.Database().StorageKey("element array 2") );
      INDEX<ARRAY,ELEMENT>  eaKey3      ( model.Database().StorageKey("element array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,ELEMENT>  efaKey1     ( model.Database().StorageKey("element flagged array 1") );
      INDEX<FLAGGEDARRAY,ELEMENT>  efaKey2     ( model.Database().StorageKey("element flagged array 2") );
      INDEX<FLAGGEDARRAY,ELEMENT>  efaKey3     ( model.Database().StorageKey("element flagged array 3") );
      indexCache.key1_ = eaKey1;

      // ========================================================================
      // PLACEMENT: FACE
      // TYPE: Scalar
      INDEX<SCALAR,FACE>  fKey1       ( model.Database().StorageKey("face variable 1") );
      INDEX<SCALAR,FACE>  fKey2       ( model.Database().StorageKey("face variable 2") );
      // TYPE: Vector
      INDEX<VECTOR,FACE>  fvKey1      ( model.Database().StorageKey("face vector 1") );
      INDEX<VECTOR,FACE>  fvKey2      ( model.Database().StorageKey("face vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,FACE>  ftKey1      ( model.Database().StorageKey("face tensor 1") );
      INDEX<TENSOR,FACE>  ftKey2      ( model.Database().StorageKey("face tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,FACE>  faKey1      ( model.Database().StorageKey("face array 1") );
      INDEX<ARRAY,FACE>  faKey2      ( model.Database().StorageKey("face array 2") );
      INDEX<ARRAY,FACE>  faKey3      ( model.Database().StorageKey("face array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,FACE>  ffaKey1     ( model.Database().StorageKey("face flagged array 1") );
      INDEX<FLAGGEDARRAY,FACE>  ffaKey2     ( model.Database().StorageKey("face flagged array 2") );
      INDEX<FLAGGEDARRAY,FACE>  ffaKey3     ( model.Database().StorageKey("face flagged array 3") );

      // ========================================================================
      // PLACEMENT: Element Integration Point
      // TYPE: Scalar
      INDEX<SCALAR,ELEMENT_INTEGRATION_POINT>  eipKey1( model.Database().StorageKey("eip scalar 1") );
      INDEX<SCALAR,ELEMENT_INTEGRATION_POINT>  eipKey2( model.Database().StorageKey("eip scalar 2") );
      // TYPE: Vector
      INDEX<VECTOR,ELEMENT_INTEGRATION_POINT>  eipvKey1( model.Database().StorageKey("eip vector 1") );
      INDEX<VECTOR,ELEMENT_INTEGRATION_POINT>  eipvKey2( model.Database().StorageKey("eip vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,ELEMENT_INTEGRATION_POINT>  eiptKey1( model.Database().StorageKey("eip tensor 1") );
      INDEX<TENSOR,ELEMENT_INTEGRATION_POINT>  eiptKey2( model.Database().StorageKey("eip tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,ELEMENT_INTEGRATION_POINT>  eipaKey1( model.Database().StorageKey("eip array 1") );
      INDEX<ARRAY,ELEMENT_INTEGRATION_POINT>  eipaKey2( model.Database().StorageKey("eip array 2") );
      INDEX<ARRAY,ELEMENT_INTEGRATION_POINT>  eipaKey3( model.Database().StorageKey("eip array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,ELEMENT_INTEGRATION_POINT>  eipfaKey1( model.Database().StorageKey("eip flagged array 1") );
      INDEX<FLAGGEDARRAY,ELEMENT_INTEGRATION_POINT>  eipfaKey2( model.Database().StorageKey("eip flagged array 2") );
      INDEX<FLAGGEDARRAY,ELEMENT_INTEGRATION_POINT>  eipfaKey3( model.Database().StorageKey("eip flagged array 3") );
      indexCache.key2_ = model.Database().StorageKey("eip array 1");

      // ========================================================================
      // PLACEMENT: Face Integration Point
      // TYPE: Scalar
      INDEX<SCALAR,FACE_INTEGRATION_POINT>  fipKey1     ( model.Database().StorageKey("fip scalar 1") );
      INDEX<SCALAR,FACE_INTEGRATION_POINT>  fipKey2     ( model.Database().StorageKey("fip scalar 2") );
      // TYPE: Vector
      INDEX<VECTOR,FACE_INTEGRATION_POINT>  fipvKey1    ( model.Database().StorageKey("fip vector 1") );
      INDEX<VECTOR,FACE_INTEGRATION_POINT>  fipvKey2    ( model.Database().StorageKey("fip vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,FACE_INTEGRATION_POINT>  fiptKey1    ( model.Database().StorageKey("fip tensor 1") );
      INDEX<TENSOR,FACE_INTEGRATION_POINT>  fiptKey2    ( model.Database().StorageKey("fip tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,FACE_INTEGRATION_POINT>  fipaKey1    ( model.Database().StorageKey("fip array 1") );
      INDEX<ARRAY,FACE_INTEGRATION_POINT>  fipaKey2    ( model.Database().StorageKey("fip array 2") );
      INDEX<ARRAY,FACE_INTEGRATION_POINT>  fipaKey3    ( model.Database().StorageKey("fip array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,FACE_INTEGRATION_POINT>  fipfaKey1( model.Database().StorageKey("fip flagged array 1") );
      INDEX<FLAGGEDARRAY,FACE_INTEGRATION_POINT>  fipfaKey2( model.Database().StorageKey("fip flagged array 2") );
      INDEX<FLAGGEDARRAY,FACE_INTEGRATION_POINT>  fipfaKey3( model.Database().StorageKey("fip flagged array 3") );

      // ========================================================================
      // PLACEMENT: Sector Integration Point
      // TYPE: Scalar
      INDEX<SCALAR,SECTOR_INTEGRATION_POINT>  seipKey1    ( model.Database().StorageKey("seip scalar 1") );
      INDEX<SCALAR,SECTOR_INTEGRATION_POINT>  seipKey2    ( model.Database().StorageKey("seip scalar 2") );
      // TYPE: Vector
      INDEX<VECTOR,SECTOR_INTEGRATION_POINT>  seipvKey1   ( model.Database().StorageKey("seip vector 1") );
      INDEX<VECTOR,SECTOR_INTEGRATION_POINT>  seipvKey2   ( model.Database().StorageKey("seip vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,SECTOR_INTEGRATION_POINT>  seiptKey1   ( model.Database().StorageKey("seip tensor 1") );
      INDEX<TENSOR,SECTOR_INTEGRATION_POINT>  seiptKey2   ( model.Database().StorageKey("seip tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,SECTOR_INTEGRATION_POINT>  seipaKey1   ( model.Database().StorageKey("seip array 1") );
      INDEX<ARRAY,SECTOR_INTEGRATION_POINT>  seipaKey2   ( model.Database().StorageKey("seip array 2") );
      INDEX<ARRAY,SECTOR_INTEGRATION_POINT>  seipaKey3   ( model.Database().StorageKey("seip array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,SECTOR_INTEGRATION_POINT>  seipfaKey1( model.Database().StorageKey("seip flagged array 1") );
      INDEX<FLAGGEDARRAY,SECTOR_INTEGRATION_POINT>  seipfaKey2( model.Database().StorageKey("seip flagged array 2") );
      INDEX<FLAGGEDARRAY,SECTOR_INTEGRATION_POINT>  seipfaKey3( model.Database().StorageKey("seip flagged array 3") );

      // ========================================================================
      // PLACEMENT: Facet Integration Point
      // TYPE: Scalar
      INDEX<SCALAR,FACET_INTEGRATION_POINT>  faipKey1    ( model.Database().StorageKey("faip scalar 1") );
      INDEX<SCALAR,FACET_INTEGRATION_POINT>  faipKey2    ( model.Database().StorageKey("faip scalar 2") );
      // TYPE: Vector
      INDEX<VECTOR,FACET_INTEGRATION_POINT>  faipvKey1   ( model.Database().StorageKey("faip vector 1") );
      INDEX<VECTOR,FACET_INTEGRATION_POINT>  faipvKey2   ( model.Database().StorageKey("faip vector 2") );
      // TYPE: Tensor
      INDEX<TENSOR,FACET_INTEGRATION_POINT>  faiptKey1   ( model.Database().StorageKey("faip tensor 1") );
      INDEX<TENSOR,FACET_INTEGRATION_POINT>  faiptKey2   ( model.Database().StorageKey("faip tensor 2") );
      // TYPE: Array
      INDEX<ARRAY,FACET_INTEGRATION_POINT>  faipaKey1   ( model.Database().StorageKey("faip array 1") );
      INDEX<ARRAY,FACET_INTEGRATION_POINT>  faipaKey2   ( model.Database().StorageKey("faip array 2") );
      INDEX<ARRAY,FACET_INTEGRATION_POINT>  faipaKey3   ( model.Database().StorageKey("faip array 3") );
      // TYPE: Flagged Array
      INDEX<FLAGGEDARRAY,FACET_INTEGRATION_POINT>  faipfaKey1( model.Database().StorageKey("faip flagged array 1") );
      INDEX<FLAGGEDARRAY,FACET_INTEGRATION_POINT>  faipfaKey2( model.Database().StorageKey("faip flagged array 2") );
      INDEX<FLAGGEDARRAY,FACET_INTEGRATION_POINT>  faipfaKey3( model.Database().StorageKey("faip flagged array 3") );

      // ========================================================================
      // Some working variables
      // Scalar's
      ScalarVariable scalarV;
      ScalarVariable oneS               ( ANY, 1.       );
      ScalarVariable twoS               ( DIRICH, 2.    );
      ScalarVariable threeS             ( ROBIN, 3.     );
      // Vector's
      VectorVariable<dim> vectorV;
      VectorVariable<dim> fourV           ( NEUMANN, 4.   );
      VectorVariable<dim> fiveV           ( ANY, 5.       );
      VectorVariable<dim> sixV            ( ANY, 6.       );
      // Tensor's
      TensorVariable<dim> tensorV;
      TensorVariable<dim> sevenT          ( DIRICH, 7.    );
      TensorVariable<dim> eightT          ( ROBIN, 8.     );
      TensorVariable<dim> nineT           ( NEUMANN, 9.   );
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
      model.Boundary("BOTTOM").InputPropertyValue( "face variable 1", makeScalar( PLAIN, 4. ) );
      _test( ( *model.Boundary("BOTTOM").CellsBegin() )->Read(fKey1) == 4. );

      model.Database().WriteVariablesFile("CSMP-variables-output.txt");


      // ========================================================================
      // Nodal Scalar Variable Test Thorough
      ScalarVariable sv( PLAIN, 99999. );
      ScalarVariable sv1( PLAIN, 1. );
      model.InputPropertyValue( "nodal variable 1", sv1 );
      model.InputPropertyValue( "nodal variable 2", makeScalar( DIRICH, 99. ) );
      Region<dim>&  model_domain = model.Region("Model");
      for( vector<Node<dim>*>::const_iterator it( model_domain.NodesBegin() ); it != model_domain.NodesEnd(); ++it )
        {
          sv = makeScalar( DIRICH, 9999. );
          (*it)->Read( nKey1, sv );
          _test( sv == sv1 );
        }
      
      // ========================================================================
      // Nodal Array Variable Test Thorough
      ArrayVariable av0(4, 10., DIRICH );
      for( vector<Node<dim>*>::const_iterator it( model_domain.NodesBegin() ); it != model_domain.NodesEnd(); ++it )
        (*it)->Store( naKey1, av0 );
      ArrayVariable av00( "nodal array 1", model.Database() );
      for( vector<Node<dim>*>::const_iterator it( model_domain.NodesBegin() ); it != model_domain.NodesEnd(); ++it )
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
      for( vector<Element<dim>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
        {
          // TYPE: Array
          (*it)->Read( eaKey2, av1 );
          _test( av1 == av00 );
        }


      // ========================================================================
      // Nodal Flagged Array Variable Test Thorough
      FlaggedArrayVariable fav0(1, 10., DIRICH );
      for( vector<Node<dim>*>::const_iterator it( model_domain.NodesBegin() ); it != model_domain.NodesEnd(); ++it )
        (*it)->Store( nfaKey1, fav0 );
      FlaggedArrayVariable fav00( "nodal flagged array 1", model.Database() );
      for( vector<Node<2>*>::const_iterator it( model_domain.NodesBegin() ); it != model_domain.NodesEnd(); ++it )
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
      for( vector<Element<2>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
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
      model.Boundary("TOP").InputPropertyValue( "boundary tensor 1", sevenT );

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
      model_domain.Read( rvKey1, vectorV );
      _test( vectorV == fourV );
      // TYPE: Tensor
      Boundary<2U>& boundary2 = model.Boundary("TOP");
      boundary2.Read( btKey1, tensorV );
      _test( tensorV == sevenT );

      // ========================================================================
      // PLACEMENT: Node
      for( vector<Node<dim>*>::const_iterator it( model_domain.NodesBegin() ); it != model_domain.NodesEnd(); ++it )
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

      for( vector<Element<dim>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
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
          for( auto i{0U}; i < (*it)->IntegrationPoints(); ++i )
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
      Boundary<dim>& left{ model.Boundary("LEFT") };
      for( vector<Face<dim>*>::const_iterator it( left.CellsBegin() ); it !=left.CellsEnd(); ++it )
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
        for( auto i{0U}; i < (*it)->IntegrationPoints(); ++i )
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

      for( vector<Element<2>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
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
          for( auto i{0U}; i < (*it)->IntegrationPoints(); ++i )
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
          const auto sectors( (*it)->Sectors() );
          assert( sectors != 0 );
          for( auto j(0U); j < sectors; ++j )
            for( auto i{0U}; i < (*it)->IntegrationPointsPerSector(); ++i )
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
            const auto facets( (*it)->Facets() );
            assert( facets != 0 );
            for( auto j(0U); j < facets; ++j )
              for( auto i{0U}; i < (*it)->IntegrationPointsPerFacet(); ++i )
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
      for( vector<Element<2>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
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
        for( auto i{0U}; i < (*it)->IntegrationPoints(); ++i )
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
        const auto sectors( (*it)->Sectors() );
        assert( sectors != 0 );
        for( auto j(0U); j < sectors; ++j )
          for( auto i{0U}; i < (*it)->IntegrationPointsPerSector(); ++i )
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
          const auto facets( (*it)->Facets() );
          assert( facets != 0 );
          for( auto j(0U); j < facets; ++j )
            for( auto i{0U}; i < (*it)->IntegrationPointsPerFacet(); ++i )
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
      csmp::Index  neweKey3  = model.CreateProperty( "new element scalar 3", "nes3", "X", SCALAR, ELEMENT );
      model.InputPropertyValue( "new element scalar 3", threeS );
      // TYPE: Vector
      csmp::Index  newevKey3 = model.CreateProperty( "new element vector 3", "nev3", "X", VECTOR, ELEMENT );
      model.InputPropertyValue( "new element vector 3", sixV );
      // TYPE: Tensor
      csmp::Index  newetKey3 = model.CreateProperty( "new element tensor 3", "net3", "X", TENSOR, ELEMENT );
      model.InputPropertyValue( "new element tensor 3", nineT );
      // TYPE: Array
      csmp::Index  neweaKey3 = model.CreateProperty( "new element array 3", "nea3", "X", ARRAY, ELEMENT, 28 );
      model.InputPropertyValue( "new element array 3", elementArray3 );
      // TYPE: Flagged Array
      csmp::Index  newefaKey3 = model.CreateProperty( "new element flagged array 3", "nefa3", "X", FLAGGEDARRAY, ELEMENT, 3 );
      model.InputPropertyValue( "new element flagged array 3", elementFlaggedArray3 );

      // PLACEMENT: Element Integration Point
      // TYPE: Scalar
      csmp::Index  neweipKey3  = model.CreateProperty( "new eip scalar 3", "nsei3", "X", SCALAR, ELEMENT_INTEGRATION_POINT );
      model.InputPropertyValue( "new eip scalar 3", threeS );
      // TYPE: Vector
      csmp::Index  neweipvKey3 = model.CreateProperty( "new eip vector 3", "neiv3", "X", VECTOR, ELEMENT_INTEGRATION_POINT );
      model.InputPropertyValue( "new eip vector 3", sixV );
      // TYPE: Tensor
      csmp::Index  neweiptKey3 = model.CreateProperty( "new eip tensor 3", "neit3", "X", TENSOR, ELEMENT_INTEGRATION_POINT );
      model.InputPropertyValue( "new eip tensor 3", nineT );
      // TYPE: Array
      csmp::Index  neweipaKey3 = model.CreateProperty( "new eip array 3", "neia3", "X", ARRAY, ELEMENT_INTEGRATION_POINT, 28 );
      model.InputPropertyValue( "new eip array 3", elementArray3 );
      // TYPE: Flagged Array
      csmp::Index  neweipfaKey3 = model.CreateProperty( "new eip flagged array 3", "neifa3", "X", FLAGGEDARRAY, ELEMENT_INTEGRATION_POINT, 3 );
      model.InputPropertyValue( "new eip flagged array 3", elementFlaggedArray3 );

      // PLACEMENT: Sector Integration Point
      // TYPE: Scalar
      csmp::Index  newseipKey3  = model.CreateProperty( "new seip scalar 3", "nseis3", "X", SCALAR, SECTOR_INTEGRATION_POINT );
      model.InputPropertyValue( "new seip scalar 3", threeS );
      // TYPE: Vector
      csmp::Index  newseipvKey3 = model.CreateProperty( "new seip vector 3", "nseiv3", "X", VECTOR, SECTOR_INTEGRATION_POINT );
      model.InputPropertyValue( "new seip vector 3", sixV );
      // TYPE: Tensor
      csmp::Index  newseiptKey3 = model.CreateProperty( "new seip tensor 3", "nseit3", "X", TENSOR, SECTOR_INTEGRATION_POINT );
      model.InputPropertyValue( "new seip tensor 3", nineT );
      // TYPE: Array
      csmp::Index  newseipaKey3 = model.CreateProperty( "new seip array 3", "nseia3", "X", ARRAY, SECTOR_INTEGRATION_POINT, 28 );
      model.InputPropertyValue( "new seip array 3", elementArray3 );
      // TYPE: Flagged Array
      csmp::Index  newseipfaKey3 = model.CreateProperty( "new seip flagged array 3", "nseifa3", "X", FLAGGEDARRAY, SECTOR_INTEGRATION_POINT, 3 );
      model.InputPropertyValue( "new seip flagged array 3", elementFlaggedArray3 );

      // PLACEMENT: Facet Integration Point
      // TYPE: Scalar
      csmp::Index  newfaipKey3  = model.CreateProperty( "new faip scalar 3", "nfais3", "X", SCALAR, FACET_INTEGRATION_POINT );
      model.InputPropertyValue( "new faip scalar 3", threeS );
      // TYPE: Vector
      csmp::Index  newfaipvKey3 = model.CreateProperty( "new faip vector 3", "nfaiv3", "X", VECTOR, FACET_INTEGRATION_POINT );
      model.InputPropertyValue( "new faip vector 3", sixV );
      // TYPE: Tensor
      csmp::Index  newfaiptKey3 = model.CreateProperty( "new faip tensor 3", "nfait3", "X", TENSOR, FACET_INTEGRATION_POINT );
      model.InputPropertyValue( "new faip tensor 3", nineT );
      // TYPE: Array
      csmp::Index  newfaipaKey3 = model.CreateProperty( "new faip array 3", "nfaia3", "X", ARRAY, FACET_INTEGRATION_POINT, 28 );
      model.InputPropertyValue( "new faip array 3", elementArray3 );
      // TYPE: Flagged Array
      csmp::Index  newfaipfaKey3 = model.CreateProperty( "new faip flagged array 3", "nfaifa3", "X", FLAGGEDARRAY, FACET_INTEGRATION_POINT, 3 );
      model.InputPropertyValue( "new faip flagged array 3", elementFlaggedArray3 );

      // ========================================================================
      // PLACEMENT: Element & Integration Points
      for( vector<Element<dim>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
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
        for( auto i{0U}; i < (*it)->IntegrationPoints(); ++i )
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
        const auto sectors( (*it)->Sectors() );
        assert( sectors != 0 );
        for( auto j(0U); j < sectors; ++j )
          for( auto i{0U}; i < (*it)->IntegrationPointsPerSector(); ++i )
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
          const auto facets( (*it)->Facets() );
          assert( facets != 0 );
          for( auto j(0U); j < facets; ++j )
            for( auto i{0U}; i < (*it)->IntegrationPointsPerFacet(); ++i )
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
      for( vector<Element<dim>*>::const_iterator it( model_domain.CellsBegin() ); it != model_domain.CellsEnd(); ++it )
        {
        // TYPE: Scalar
        _test( (*it)->Status(eKey2)     == twoS.Flag() );
        _test( (*it)->Status(neweKey3)  == threeS.Flag() );
        for( auto d{0U}; d<dim; ++d )
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
        for( auto d{0U}; d<elementFlaggedArray2.Size(); ++d )
            _test( (*it)->Status(efaKey2,d) == elementFlaggedArray2.Flag(d) );
        for( auto d{0U}; d<elementFlaggedArray3.Size(); ++d )
            _test( (*it)->Status(newefaKey3,d) == elementFlaggedArray3.Flag(d) );

        // ========================================================================
        // PLACEMENT: Element Integration Point
        for( auto i{0U}; i < (*it)->IntegrationPoints(); ++i )
          {
          // TYPE: Scalar
          _test( (*it)->Status(i,eipKey2)       == twoS.Flag() );
          _test( (*it)->Status(i,neweipKey3)    == threeS.Flag() );
          for( auto d{0U}; d<dim; ++d )
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
          for( auto d{0U}; d<eipFlaggedArray3.Size(); ++d )
              _test( (*it)->Status(i,eipfaKey3,d) == eipFlaggedArray3.Flag(d) );
          for( auto d{0U}; d<elementFlaggedArray3.Size(); ++d )
              _test( (*it)->Status(i,neweipfaKey3,d) == elementFlaggedArray3.Flag(d) );
          }


        // ========================================================================
        // PLACEMENT: Sector Integration Point
        const auto sectors( (*it)->Sectors() );
        assert( sectors != 0 );
        for( auto j(0U); j < sectors; ++j )
          for( auto i{0U}; i < (*it)->IntegrationPointsPerSector(); ++i )
            {
            // TYPE: Scalar
            _test( (*it)->Status(j,i,seipKey2)      == twoS.Flag() );
            _test( (*it)->Status(j,i,newseipKey3)   == threeS.Flag() );
            for( auto d{0U}; d<dim; ++d )
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
            for( auto d{0U}; d<fvipFlaggedArray2.Size(); ++d )
                _test( (*it)->Status(j,i,seipfaKey2,d) == fvipFlaggedArray2.Flag(d) );
            for( auto d{0U}; d<fvipFlaggedArray3.Size(); ++d )
                _test( (*it)->Status(j,i,seipfaKey3,d) == fvipFlaggedArray3.Flag(d) );
            for( auto d{0U}; d<elementFlaggedArray3.Size(); ++d )
                _test( (*it)->Status(j,i,newseipfaKey3,d) == elementFlaggedArray3.Flag(d) );
            }

          // ========================================================================
          // PLACEMENT: Facet Integration Point
          const auto facets( (*it)->Facets() );
          assert( facets != 0 );
          for( auto j(0U); j < facets; ++j )
            for( auto i{0U}; i < (*it)->IntegrationPointsPerFacet(); ++i )
              {
              // TYPE: Scalar
              _test( (*it)->Status(j,i,faipKey2)    == twoS.Flag() );
              _test( (*it)->Status(j,i,newfaipKey3) == threeS.Flag() );
              for( auto d{0U}; d<dim; ++d )
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
              for( auto d{0U}; d<fvipFlaggedArray2.Size(); ++d )
                  _test( (*it)->Status(j,i,faipfaKey2,d) == fvipFlaggedArray2.Flag(d) );
              for( auto d{0U}; d<fvipFlaggedArray3.Size(); ++d )
                  _test( (*it)->Status(j,i,faipfaKey3,d) == fvipFlaggedArray3.Flag(d) );
              for( auto d{0U}; d<elementFlaggedArray3.Size(); ++d )
                  _test( (*it)->Status(j,i,newfaipfaKey3,d) == elementFlaggedArray3.Flag(d) );

              }
        }
      
    } // end runModel

} // csmp
