// PropertyHandle_Test.cpp
// Unit test for PropertyHandle.
// 2024

#include "PropertyHandle_Test.h"
#include "ANSYS_Model3D.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "PropertyHandle.h"

using namespace std;

namespace csmp {

// ============================================================================
//  Constructor / destructor
// ============================================================================

/** Creates the split boundary on the model and returns the model pointer.
    Called during member initialisation so the split boundary exists
    before any SplitBoundary PropertyHandle objects are constructed.
*/
static csmp::Model<3>* initialiseSplitBoundary( csmp::Model<3>* model )
{
    const bool retain_elmts_as_intervening_elements = true;
    model->CreateSplitBoundaryFrom( "HALF",
                                    retain_elmts_as_intervening_elements );
    model->SplitBoundariesOut();
    return model;
}


PropertyHandle_Test::PropertyHandle_Test( double tolerance, bool verbose )
    : model_( initialiseSplitBoundary(
                  new ANSYS_Model3D( "BoxHalfs3D",
                                     "PropertyHandle_Test-variables.txt" ) ) ),
      TOLERANCE( tolerance ),
      verbose_( verbose ),

      // --- Region SCALAR ---
      reg_sc_elem1( *model_, "reg sc elem 1", SCALAR, ELEMENT ),
      reg_sc_elem2( *model_, "reg sc elem 2", SCALAR, ELEMENT ),
      reg_sc_node1( *model_, "reg sc node 1", SCALAR, NODE ),
      reg_sc_node2( *model_, "reg sc node 2", SCALAR, NODE ),
      reg_sc_ip1  ( *model_, "reg sc ip 1",   SCALAR, ELEMENT_INTEGRATION_POINT ),
      reg_sc_ip2  ( *model_, "reg sc ip 2",   SCALAR, ELEMENT_INTEGRATION_POINT ),
      reg_sc_sec1 ( *model_, "reg sc sec 1",  SCALAR, SECTOR_INTEGRATION_POINT ),
      reg_sc_sec2 ( *model_, "reg sc sec 2",  SCALAR, SECTOR_INTEGRATION_POINT ),
      reg_sc_fac1 ( *model_, "reg sc fac 1",  SCALAR, FACET_INTEGRATION_POINT ),
      reg_sc_fac2 ( *model_, "reg sc fac 2",  SCALAR, FACET_INTEGRATION_POINT ),

      // --- Region VECTOR ---
      reg_vc_elem1( *model_, "reg vc elem 1", VECTOR, ELEMENT ),
      reg_vc_elem2( *model_, "reg vc elem 2", VECTOR, ELEMENT ),
      reg_vc_node1( *model_, "reg vc node 1", VECTOR, NODE ),
      reg_vc_node2( *model_, "reg vc node 2", VECTOR, NODE ),
      reg_vc_ip1  ( *model_, "reg vc ip 1",   VECTOR, ELEMENT_INTEGRATION_POINT ),
      reg_vc_ip2  ( *model_, "reg vc ip 2",   VECTOR, ELEMENT_INTEGRATION_POINT ),

      // --- Region TENSOR ---
      reg_ts_elem1( *model_, "reg ts elem 1", TENSOR, ELEMENT ),
      reg_ts_elem2( *model_, "reg ts elem 2", TENSOR, ELEMENT ),
      reg_ts_node1( *model_, "reg ts node 1", TENSOR, NODE ),
      reg_ts_node2( *model_, "reg ts node 2", TENSOR, NODE ),
      reg_ts_ip1  ( *model_, "reg ts ip 1",   TENSOR, ELEMENT_INTEGRATION_POINT ),
      reg_ts_ip2  ( *model_, "reg ts ip 2",   TENSOR, ELEMENT_INTEGRATION_POINT ),

      // --- Boundary SCALAR ---
      bnd_sc_node1( *model_, "LEFT", "bnd sc node 1", SCALAR, NODE ),
      bnd_sc_node2( *model_, "LEFT", "bnd sc node 2", SCALAR, NODE ),
      bnd_sc_face1( *model_, "LEFT", "bnd sc face 1", SCALAR, FACE ),
      bnd_sc_face2( *model_, "LEFT", "bnd sc face 2", SCALAR, FACE ),
      bnd_sc_fip1 ( *model_, "LEFT", "bnd sc fip 1",  SCALAR, FACE_INTEGRATION_POINT ),
      bnd_sc_fip2 ( *model_, "LEFT", "bnd sc fip 2",  SCALAR, FACE_INTEGRATION_POINT ),
      bnd_sc_fsec1( *model_, "LEFT", "bnd sc fsec 1", SCALAR, FACE_SECTOR_INTEGRATION_POINT ),
      bnd_sc_fsec2( *model_, "LEFT", "bnd sc fsec 2", SCALAR, FACE_SECTOR_INTEGRATION_POINT ),
      bnd_sc_ffac1( *model_, "LEFT", "bnd sc ffac 1", SCALAR, FACE_FACET_INTEGRATION_POINT ),
      bnd_sc_ffac2( *model_, "LEFT", "bnd sc ffac 2", SCALAR, FACE_FACET_INTEGRATION_POINT ),

      // --- Boundary VECTOR ---
      bnd_vc_face1( *model_, "LEFT", "bnd vc face 1", VECTOR, FACE ),
      bnd_vc_face2( *model_, "LEFT", "bnd vc face 2", VECTOR, FACE ),

      // --- Boundary TENSOR ---
      bnd_ts_face1( *model_, "LEFT", "bnd ts face 1", TENSOR, FACE ),
      bnd_ts_face2( *model_, "LEFT", "bnd ts face 2", TENSOR, FACE ),

      // --- SplitBoundary SCALAR ---
      // Safe to construct now: CreateSplitBoundaryFrom was called
      // inside initialiseSplitBoundary before any member initialisation.
      spb_sc_node1 ( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc node 1",  SCALAR, NODE ),
      spb_sc_node2 ( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc node 2",  SCALAR, NODE ),
      spb_sc_iface1( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc iface 1", SCALAR, INTER_FACE ),
      spb_sc_iface2( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc iface 2", SCALAR, INTER_FACE ),
      spb_sc_ifip1 ( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc ifip 1",  SCALAR, INTER_FACE_INTEGRATION_POINT ),
      spb_sc_ifip2 ( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc ifip 2",  SCALAR, INTER_FACE_INTEGRATION_POINT ),
      spb_sc_ifsec1( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc ifsec 1", SCALAR, INTER_FACE_SECTOR_INTEGRATION_POINT ),
      spb_sc_ifsec2( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc ifsec 2", SCALAR, INTER_FACE_SECTOR_INTEGRATION_POINT ),
      spb_sc_iffac1( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc iffac 1", SCALAR, INTER_FACE_FACET_INTEGRATION_POINT ),
      spb_sc_iffac2( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb sc iffac 2", SCALAR, INTER_FACE_FACET_INTEGRATION_POINT ),

      // --- SplitBoundary VECTOR ---
      spb_vc_iface1( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb vc iface 1", VECTOR, INTER_FACE ),
      spb_vc_iface2( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb vc iface 2", VECTOR, INTER_FACE ),

      // --- SplitBoundary TENSOR ---
      spb_ts_iface1( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb ts iface 1", TENSOR, INTER_FACE ),
      spb_ts_iface2( *model_, "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT", "spb ts iface 2", TENSOR, INTER_FACE ),

      // --- whole-model sentinels ---
      whole_sc_elem( *model_, "whole sc elem", SCALAR, ELEMENT ),
      whole_sc_node( *model_, "whole sc node", SCALAR, NODE )
{
}



PropertyHandle_Test::~PropertyHandle_Test()
{
}

// ============================================================================
//  ResetVariables
// ============================================================================

void PropertyHandle_Test::ResetVariables()
{
    // Region SCALAR
    reg_sc_elem1 = 1.5;   reg_sc_elem2 = 2.0;
    reg_sc_node1 = 3.0;   reg_sc_node2 = 4.0;
    reg_sc_ip1   = 5.0;   reg_sc_ip2   = 6.0;
    reg_sc_sec1  = 7.0;   reg_sc_sec2  = 8.0;
    reg_sc_fac1  = 9.0;   reg_sc_fac2  = 10.0;

    // Region VECTOR
    reg_vc_elem1 = 1.5;   reg_vc_elem2 = 2.0;
    reg_vc_node1 = 3.0;   reg_vc_node2 = 4.0;
    reg_vc_ip1   = 5.0;   reg_vc_ip2   = 6.0;

    // Region TENSOR
    reg_ts_elem1 = 1.5;   reg_ts_elem2 = 2.0;
    reg_ts_node1 = 3.0;   reg_ts_node2 = 4.0;
    reg_ts_ip1   = 5.0;   reg_ts_ip2   = 6.0;

    // Boundary SCALAR
    bnd_sc_node1 = 1.0;   bnd_sc_node2 = 2.0;
    bnd_sc_face1 = 3.0;   bnd_sc_face2 = 4.0;
    bnd_sc_fip1  = 5.0;   bnd_sc_fip2  = 6.0;
    bnd_sc_fsec1 = 7.0;   bnd_sc_fsec2 = 8.0;
    bnd_sc_ffac1 = 9.0;   bnd_sc_ffac2 = 10.0;

    // Boundary VECTOR / TENSOR
    bnd_vc_face1 = 3.0;   bnd_vc_face2 = 4.0;
    bnd_ts_face1 = 3.0;   bnd_ts_face2 = 4.0;

    // SplitBoundary SCALAR
    spb_sc_node1  = 1.0;  spb_sc_node2  = 2.0;
    spb_sc_iface1 = 3.0;  spb_sc_iface2 = 4.0;
    spb_sc_ifip1  = 5.0;  spb_sc_ifip2  = 6.0;
    spb_sc_ifsec1 = 7.0;  spb_sc_ifsec2 = 8.0;
    spb_sc_iffac1 = 9.0;  spb_sc_iffac2 = 10.0;

    // SplitBoundary VECTOR / TENSOR
    spb_vc_iface1 = 3.0;  spb_vc_iface2 = 4.0;
    spb_ts_iface1 = 3.0;  spb_ts_iface2 = 4.0;
}

// ============================================================================
//  run
// ============================================================================

void PropertyHandle_Test::run()
{
    if ( verbose_ )
    {
        cout << "\n=========================" << endl;
        cout << "Testing PropertyHandle"    << endl;
        cout << "========================="  << endl;
    }

    TestRegionPlacements();
    TestBoundaryPlacements();
    TestSplitBoundaryPlacements();
    TestRegionCrossPlacement();
    TestBoundaryCrossPlacement();
    TestSplitBoundaryCrossPlacement();
    TestSubdomainIsolation();
    TestModelToSubdomainAssignment();
    TestApplyMethod();
    TestMathMethods();
    TestArraySupport();
    TestRegionCrossPlacementFull();
    TestFlaggedArrayFacetIP();
}



// ============================================================================
//  TestRegionPlacements
// ============================================================================

void PropertyHandle_Test::TestRegionPlacements()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- Region placements ---" << endl;

    ScalarVariable    sc;
    VectorVariable<3> vc;
    TensorVariable<3> ts;

    // ------------------------------------------------------------------
    //  SCALAR ELEMENT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR ELEMENT" << endl;
    {
        Index key = model_->Database().StorageKey( "reg sc elem 2" );

        // += double
        reg_sc_elem2 += 1.0;
        // expected: 2.0 + 1.0 = 3.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );

        // -= double
        reg_sc_elem2 -= 1.0;
        // expected: 3.0 - 1.0 = 2.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 2.0, numeric_limits<double>::epsilon() * 2.0 );

        // *= double
        reg_sc_elem2 *= 3.0;
        // expected: 2.0 * 3.0 = 6.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 6.0, numeric_limits<double>::epsilon() * 6.0 );

        // /= double
        reg_sc_elem2 /= 2.0;
        // expected: 6.0 / 2.0 = 3.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR NODE
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR NODE" << endl;
    {
        Index key = model_->Database().StorageKey( "reg sc node 2" );

        reg_sc_node2 += 1.0;
        // expected: 4.0 + 1.0 = 5.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );

        reg_sc_node2 -= 2.0;
        // expected: 5.0 - 2.0 = 3.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );

        reg_sc_node2 *= 4.0;
        // expected: 3.0 * 4.0 = 12.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 12.0, numeric_limits<double>::epsilon() * 12.0 );

        reg_sc_node2 /= 3.0;
        // expected: 12.0 / 3.0 = 4.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR ELEMENT_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR ELEMENT_INTEGRATION_POINT" << endl;
    {
        Index key = model_->Database().StorageKey( "reg sc ip 2" );

        reg_sc_ip2 += 1.0;
        // expected: 6.0 + 1.0 = 7.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 7.0, numeric_limits<double>::epsilon() * 7.0 );

        reg_sc_ip2 -= 2.0;
        // expected: 7.0 - 2.0 = 5.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );

        reg_sc_ip2 *= 2.0;
        // expected: 5.0 * 2.0 = 10.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 10.0, numeric_limits<double>::epsilon() * 10.0 );

        reg_sc_ip2 /= 5.0;
        // expected: 10.0 / 5.0 = 2.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 2.0, numeric_limits<double>::epsilon() * 2.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR SECTOR_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR SECTOR_INTEGRATION_POINT" << endl;
    {
        reg_sc_sec2 += 1.0;
        // expected: 8.0 + 1.0 = 9.0
        // Verify via Apply reading back the value.
        double readback = 0.0;
        reg_sc_sec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 9.0, numeric_limits<double>::epsilon() * 9.0 );

        reg_sc_sec2 -= 1.0;
        // expected: 9.0 - 1.0 = 8.0
        reg_sc_sec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 8.0, numeric_limits<double>::epsilon() * 8.0 );

        reg_sc_sec2 *= 2.0;
        // expected: 8.0 * 2.0 = 16.0
        reg_sc_sec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 16.0, numeric_limits<double>::epsilon() * 16.0 );

        reg_sc_sec2 /= 4.0;
        // expected: 16.0 / 4.0 = 4.0
        reg_sc_sec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR FACET_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR FACET_INTEGRATION_POINT" << endl;
    {
        reg_sc_fac2 += 1.0;
        // expected: 10.0 + 1.0 = 11.0
        double readback = 0.0;
        reg_sc_fac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 11.0, numeric_limits<double>::epsilon() * 11.0 );

        reg_sc_fac2 -= 1.0;
        reg_sc_fac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 10.0, numeric_limits<double>::epsilon() * 10.0 );

        reg_sc_fac2 *= 3.0;
        reg_sc_fac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 30.0, numeric_limits<double>::epsilon() * 30.0 );

        reg_sc_fac2 /= 3.0;
        reg_sc_fac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 10.0, numeric_limits<double>::epsilon() * 10.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  VECTOR ELEMENT — same-placement PropertyHandle operators
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  VECTOR ELEMENT (PropertyHandle operators)" << endl;
    {
        Index key = model_->Database().StorageKey( "reg vc elem 2" );

        reg_vc_elem2 += reg_vc_elem1;
        // expected: 2.0 + 1.5 = 3.5
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 3.5, TOLERANCE );

        reg_vc_elem2 -= reg_vc_elem1;
        // expected: 3.5 - 1.5 = 2.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 2.0, TOLERANCE );

        reg_vc_elem2 *= reg_vc_elem1;
        // expected: 2.0 * 1.5 = 3.0 (component-wise)
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 3.0, TOLERANCE );

        reg_vc_elem2 /= reg_vc_elem1;
        // expected: 3.0 / 1.5 = 2.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 2.0, TOLERANCE );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  TENSOR ELEMENT — same-placement PropertyHandle operators
    //  Tensor *= tensor is matrix product: result(0,0) = 3 * a * b
    //  for uniform 3x3 matrices with all entries a and b.
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  TENSOR ELEMENT (PropertyHandle operators)" << endl;
    {
        Index key = model_->Database().StorageKey( "reg ts elem 2" );

        reg_ts_elem2 += reg_ts_elem1;
        // expected: 2.0 + 1.5 = 3.5
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 3.5, TOLERANCE );

        reg_ts_elem2 -= reg_ts_elem1;
        // expected: 3.5 - 1.5 = 2.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 2.0, TOLERANCE );

        reg_ts_elem2 *= reg_ts_elem1;
        // expected: 3 * 2.0 * 1.5 = 9.0 (matrix product)
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 3.0 * 2.0 * 1.5, TOLERANCE );

        ResetVariables();
        reg_ts_elem2 /= reg_ts_elem1;
        // For uniform matrices this is not well-defined (singular);
        // test /= double instead.
        reg_ts_elem2 = 6.0;
        reg_ts_elem2 /= 2.0;
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 3.0, TOLERANCE );
    }
}

// ============================================================================
//  TestBoundaryPlacements
// ============================================================================

void PropertyHandle_Test::TestBoundaryPlacements()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- Boundary placements ---" << endl;

    ScalarVariable    sc;
    VectorVariable<3> vc;
    TensorVariable<3> ts;

    const auto& bnd = model_->Boundary( "LEFT" );

    // ------------------------------------------------------------------
    //  SCALAR FACE
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR FACE" << endl;
    {
        Index key = model_->Database().StorageKey( "bnd sc face 2" );

        bnd_sc_face2 += 1.0;
        // expected: 4.0 + 1.0 = 5.0
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );

        bnd_sc_face2 -= 2.0;
        // expected: 5.0 - 2.0 = 3.0
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );

        bnd_sc_face2 *= 4.0;
        // expected: 3.0 * 4.0 = 12.0
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 12.0, numeric_limits<double>::epsilon() * 12.0 );

        bnd_sc_face2 /= 3.0;
        // expected: 12.0 / 3.0 = 4.0
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR NODE on Boundary
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR NODE (Boundary)" << endl;
    {
        Index key = model_->Database().StorageKey( "bnd sc node 2" );

        bnd_sc_node2 += 1.0;
        // expected: 2.0 + 1.0 = 3.0
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );

        bnd_sc_node2 -= 1.0;
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 2.0, numeric_limits<double>::epsilon() * 2.0 );

        bnd_sc_node2 *= 5.0;
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 10.0, numeric_limits<double>::epsilon() * 10.0 );

        bnd_sc_node2 /= 2.0;
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR FACE_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR FACE_INTEGRATION_POINT" << endl;
    {
        Index key = model_->Database().StorageKey( "bnd sc fip 2" );

        bnd_sc_fip2 += 1.0;
        // expected: 6.0 + 1.0 = 7.0
        ( *(bnd.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 7.0, numeric_limits<double>::epsilon() * 7.0 );

        bnd_sc_fip2 -= 2.0;
        ( *(bnd.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );

        bnd_sc_fip2 *= 2.0;
        ( *(bnd.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 10.0, numeric_limits<double>::epsilon() * 10.0 );

        bnd_sc_fip2 /= 5.0;
        ( *(bnd.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 2.0, numeric_limits<double>::epsilon() * 2.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR FACE_SECTOR_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR FACE_SECTOR_INTEGRATION_POINT" << endl;
    {
        double readback = 0.0;

        bnd_sc_fsec2 += 1.0;
        // expected: 8.0 + 1.0 = 9.0
        bnd_sc_fsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 9.0, numeric_limits<double>::epsilon() * 9.0 );

        bnd_sc_fsec2 -= 1.0;
        bnd_sc_fsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 8.0, numeric_limits<double>::epsilon() * 8.0 );

        bnd_sc_fsec2 *= 2.0;
        bnd_sc_fsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 16.0, numeric_limits<double>::epsilon() * 16.0 );

        bnd_sc_fsec2 /= 4.0;
        bnd_sc_fsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR FACE_FACET_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR FACE_FACET_INTEGRATION_POINT" << endl;
    {
        double readback = 0.0;

        bnd_sc_ffac2 += 1.0;
        // expected: 10.0 + 1.0 = 11.0
        bnd_sc_ffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 11.0, numeric_limits<double>::epsilon() * 11.0 );

        bnd_sc_ffac2 -= 1.0;
        bnd_sc_ffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 10.0, numeric_limits<double>::epsilon() * 10.0 );

        bnd_sc_ffac2 *= 3.0;
        bnd_sc_ffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 30.0, numeric_limits<double>::epsilon() * 30.0 );

        bnd_sc_ffac2 /= 3.0;
        bnd_sc_ffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 10.0, numeric_limits<double>::epsilon() * 10.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  VECTOR FACE — same-placement operators
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  VECTOR FACE (PropertyHandle operators)" << endl;
    {
        Index key = model_->Database().StorageKey( "bnd vc face 2" );

        bnd_vc_face2 += bnd_vc_face1;
        // expected: 4.0 + 3.0 = 7.0
        ( *(bnd.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 7.0, TOLERANCE );

        bnd_vc_face2 -= bnd_vc_face1;
        // expected: 7.0 - 3.0 = 4.0
        ( *(bnd.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 4.0, TOLERANCE );

        bnd_vc_face2 *= 2.0;
        // expected: 4.0 * 2.0 = 8.0
        ( *(bnd.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 8.0, TOLERANCE );

        bnd_vc_face2 /= 2.0;
        // expected: 8.0 / 2.0 = 4.0
        ( *(bnd.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 4.0, TOLERANCE );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  TENSOR FACE — same-placement operators
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  TENSOR FACE (PropertyHandle operators)" << endl;
    {
        Index key = model_->Database().StorageKey( "bnd ts face 2" );

        bnd_ts_face2 += bnd_ts_face1;
        // expected: 4.0 + 3.0 = 7.0
        ( *(bnd.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 7.0, TOLERANCE );

        bnd_ts_face2 -= bnd_ts_face1;
        // expected: 7.0 - 3.0 = 4.0
        ( *(bnd.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 4.0, TOLERANCE );

        bnd_ts_face2 *= 2.0;
        // expected: 4.0 * 2.0 = 8.0
        ( *(bnd.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 8.0, TOLERANCE );

        bnd_ts_face2 /= 2.0;
        // expected: 8.0 / 2.0 = 4.0
        ( *(bnd.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 4.0, TOLERANCE );
    }
}

// ============================================================================
//  TestSplitBoundaryPlacements
// ============================================================================

void PropertyHandle_Test::TestSplitBoundaryPlacements()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- SplitBoundary placements ---" << endl;

    ScalarVariable    sc;
    VectorVariable<3> vc;
    TensorVariable<3> ts;

    const auto& spb = model_->SplitBoundary( "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT" );

    // ------------------------------------------------------------------
    //  SCALAR INTER_FACE
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR INTER_FACE" << endl;
    {
        Index key = model_->Database().StorageKey( "spb sc iface 2" );

        spb_sc_iface2 += 1.0;
        // expected: 4.0 + 1.0 = 5.0
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );

        spb_sc_iface2 -= 2.0;
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );

        spb_sc_iface2 *= 4.0;
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 12.0, numeric_limits<double>::epsilon() * 12.0 );

        spb_sc_iface2 /= 3.0;
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR NODE on SplitBoundary (InterFace FE nodes, deduplicated)
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR NODE (SplitBoundary, deduplicated)" << endl;
    {
        // Set to known value and verify via Apply readback.
        spb_sc_node2 = 5.0;
        double readback = 0.0;
        spb_sc_node2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 5.0, numeric_limits<double>::epsilon() * 5.0 );

        spb_sc_node2 += 2.0;
        // expected: 5.0 + 2.0 = 7.0
        spb_sc_node2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 7.0, numeric_limits<double>::epsilon() * 7.0 );

        spb_sc_node2 -= 3.0;
        // expected: 7.0 - 3.0 = 4.0
        spb_sc_node2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 4.0, numeric_limits<double>::epsilon() * 4.0 );

        spb_sc_node2 *= 3.0;
        // expected: 4.0 * 3.0 = 12.0
        spb_sc_node2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 12.0, numeric_limits<double>::epsilon() * 12.0 );

        spb_sc_node2 /= 4.0;
        // expected: 12.0 / 4.0 = 3.0
        spb_sc_node2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 3.0, numeric_limits<double>::epsilon() * 3.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR INTER_FACE_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR INTER_FACE_INTEGRATION_POINT" << endl;
    {
        Index key = model_->Database().StorageKey( "spb sc ifip 2" );

        spb_sc_ifip2 += 1.0;
        // expected: 6.0 + 1.0 = 7.0
        ( *(spb.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 7.0, numeric_limits<double>::epsilon() * 7.0 );

        spb_sc_ifip2 -= 2.0;
        ( *(spb.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 5.0, numeric_limits<double>::epsilon() * 5.0 );

        spb_sc_ifip2 *= 2.0;
        ( *(spb.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 10.0, numeric_limits<double>::epsilon() * 10.0 );

        spb_sc_ifip2 /= 5.0;
        ( *(spb.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 2.0, numeric_limits<double>::epsilon() * 2.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR INTER_FACE_SECTOR_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR INTER_FACE_SECTOR_INTEGRATION_POINT" << endl;
    {
        double readback = 0.0;

        spb_sc_ifsec2 += 1.0;
        // expected: 8.0 + 1.0 = 9.0
        spb_sc_ifsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 9.0, numeric_limits<double>::epsilon() * 9.0 );

        spb_sc_ifsec2 -= 1.0;
        spb_sc_ifsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 8.0, numeric_limits<double>::epsilon() * 8.0 );

        spb_sc_ifsec2 *= 2.0;
        spb_sc_ifsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 16.0, numeric_limits<double>::epsilon() * 16.0 );

        spb_sc_ifsec2 /= 4.0;
        spb_sc_ifsec2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  SCALAR INTER_FACE_FACET_INTEGRATION_POINT
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR INTER_FACE_FACET_INTEGRATION_POINT" << endl;
    {
        double readback = 0.0;

        spb_sc_iffac2 += 1.0;
        // expected: 10.0 + 1.0 = 11.0
        spb_sc_iffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 11.0, numeric_limits<double>::epsilon() * 11.0 );

        spb_sc_iffac2 -= 1.0;
        spb_sc_iffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 10.0, numeric_limits<double>::epsilon() * 10.0 );

        spb_sc_iffac2 *= 3.0;
        spb_sc_iffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 30.0, numeric_limits<double>::epsilon() * 30.0 );

        spb_sc_iffac2 /= 3.0;
        spb_sc_iffac2.Apply( [&readback]( double& x ) { readback = x; } );
        _equal( readback, 10.0, numeric_limits<double>::epsilon() * 10.0 );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  VECTOR INTER_FACE — same-placement operators
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  VECTOR INTER_FACE (PropertyHandle operators)" << endl;
    {
        Index key = model_->Database().StorageKey( "spb vc iface 2" );

        spb_vc_iface2 += spb_vc_iface1;
        // expected: 4.0 + 3.0 = 7.0
        ( *(spb.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 7.0, TOLERANCE );

        spb_vc_iface2 -= spb_vc_iface1;
        ( *(spb.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 4.0, TOLERANCE );

        spb_vc_iface2 *= 2.0;
        ( *(spb.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 8.0, TOLERANCE );

        spb_vc_iface2 /= 2.0;
        ( *(spb.CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 4.0, TOLERANCE );
    }

    ResetVariables();

    // ------------------------------------------------------------------
    //  TENSOR INTER_FACE — same-placement operators
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  TENSOR INTER_FACE (PropertyHandle operators)" << endl;
    {
        Index key = model_->Database().StorageKey( "spb ts iface 2" );

        spb_ts_iface2 += spb_ts_iface1;
        // expected: 4.0 + 3.0 = 7.0
        ( *(spb.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 7.0, TOLERANCE );

        spb_ts_iface2 -= spb_ts_iface1;
        ( *(spb.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 4.0, TOLERANCE );

        spb_ts_iface2 *= 2.0;
        ( *(spb.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 8.0, TOLERANCE );

        spb_ts_iface2 /= 2.0;
        ( *(spb.CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 4.0, TOLERANCE );
    }
}

// ============================================================================
//  TestRegionCrossPlacement
// ============================================================================

void PropertyHandle_Test::TestRegionCrossPlacement()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- Region cross-placement ---" << endl;

    ScalarVariable sc;

    // NODE += ELEMENT  (ELEMENT interpolated to NODE)
    {
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 += reg_sc_elem1;
        // expected: 4.0 + 1.5 = 5.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 5.5, TOLERANCE );
    }

    ResetVariables();

    // NODE += ELEMENT_INTEGRATION_POINT
    {
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 += reg_sc_ip1;
        // expected: 4.0 + 5.0 = 9.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    ResetVariables();

    // ELEMENT += NODE
    {
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 += reg_sc_node1;
        // expected: 2.0 + 3.0 = 5.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, TOLERANCE );
    }

    ResetVariables();

    // ELEMENT += ELEMENT_INTEGRATION_POINT
    {
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 += reg_sc_ip1;
        // expected: 2.0 + 5.0 = 7.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 7.0, TOLERANCE );
    }

    ResetVariables();

    // ELEMENT_INTEGRATION_POINT += NODE
    {
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 += reg_sc_node1;
        // expected: 6.0 + 3.0 = 9.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    ResetVariables();

    // ELEMENT_INTEGRATION_POINT += ELEMENT
    {
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 += reg_sc_elem1;
        // expected: 6.0 + 1.5 = 7.5
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 7.5, TOLERANCE );
    }

    if ( verbose_ ) cout << "  Region cross-placement: passed" << endl;
}

// ============================================================================
//  TestBoundaryCrossPlacement
// ============================================================================

void PropertyHandle_Test::TestBoundaryCrossPlacement()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- Boundary cross-placement ---" << endl;

    ScalarVariable sc;
    const auto& bnd = model_->Boundary( "LEFT" );

    // NODE += FACE
    {
        Index key = model_->Database().StorageKey( "bnd sc node 2" );
        bnd_sc_node2 += bnd_sc_face1;
        // expected: 2.0 + 3.0 = 5.0
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, TOLERANCE );
    }

    ResetVariables();

    // FACE += NODE
    {
        Index key = model_->Database().StorageKey( "bnd sc face 2" );
        bnd_sc_face2 += bnd_sc_node1;
        // expected: 4.0 + 1.0 = 5.0
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, TOLERANCE );
    }

    ResetVariables();

    // FACE += FACE_INTEGRATION_POINT
    {
        Index key = model_->Database().StorageKey( "bnd sc face 2" );
        bnd_sc_face2 += bnd_sc_fip1;
        // expected: 4.0 + 5.0 = 9.0
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    ResetVariables();

    // FACE_INTEGRATION_POINT += FACE
    {
        Index key = model_->Database().StorageKey( "bnd sc fip 2" );
        bnd_sc_fip2 += bnd_sc_face1;
        // expected: 6.0 + 3.0 = 9.0
        ( *(bnd.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    if ( verbose_ ) cout << "  Boundary cross-placement: passed" << endl;
}

// ============================================================================
//  TestSplitBoundaryCrossPlacement
// ============================================================================

void PropertyHandle_Test::TestSplitBoundaryCrossPlacement()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- SplitBoundary cross-placement ---" << endl;

    ScalarVariable sc;
    const auto& spb = model_->SplitBoundary( "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT" );

    // INTER_FACE += INTER_FACE_INTEGRATION_POINT
    {
        Index key = model_->Database().StorageKey( "spb sc iface 2" );
        spb_sc_iface2 += spb_sc_ifip1;
        // expected: 4.0 + 5.0 = 9.0
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    ResetVariables();

    // INTER_FACE_INTEGRATION_POINT += INTER_FACE
    {
        Index key = model_->Database().StorageKey( "spb sc ifip 2" );
        spb_sc_ifip2 += spb_sc_iface1;
        // expected: 6.0 + 3.0 = 9.0
        ( *(spb.CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    if ( verbose_ ) cout << "  SplitBoundary cross-placement: passed" << endl;
}

// ============================================================================
//  TestSubdomainIsolation
// ============================================================================

void PropertyHandle_Test::TestSubdomainIsolation()
{
    if ( verbose_ ) cout << "\n--- Subdomain isolation ---" << endl;

    ScalarVariable sc;

    // Set whole-model element variable to sentinel -999.
    whole_sc_elem = -999.0;

    // Set boundary face variable to 5.0.
    bnd_sc_face1 = 5.0;

    // Operate on boundary only.
    bnd_sc_face1 += 1.0;
    // expected on boundary: 6.0

    // Verify boundary face was modified.
    {
        Index key = model_->Database().StorageKey( "bnd sc face 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 6.0, numeric_limits<double>::epsilon() * 6.0 );
    }

    // Verify whole-model element variable was NOT modified.
    {
        Index key = model_->Database().StorageKey( "whole sc elem" );
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), -999.0, numeric_limits<double>::epsilon() * 999.0 );
    }

    // Same test for SplitBoundary.
    whole_sc_elem = -999.0;
    spb_sc_iface1 = 5.0;
    spb_sc_iface1 += 1.0;

    {
        Index key = model_->Database().StorageKey( "spb sc iface 1" );
        const auto& spb = model_->SplitBoundary( "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT" );
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 6.0, numeric_limits<double>::epsilon() * 6.0 );
    }

    {
        Index key = model_->Database().StorageKey( "whole sc elem" );
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), -999.0, numeric_limits<double>::epsilon() * 999.0 );
    }

    if ( verbose_ ) cout << "  Subdomain isolation: passed" << endl;
}

// ============================================================================
//  TestModelToSubdomainAssignment
// ============================================================================

void PropertyHandle_Test::TestModelToSubdomainAssignment()
{
    if ( verbose_ ) cout << "\n--- Model-to-subdomain assignment ---" << endl;

    ScalarVariable sc;

    // ------------------------------------------------------------------
    //  Test 1: whole-model NODE → Boundary NODE.
    //  Boundary nodes are a subset of the model nodes.
    // ------------------------------------------------------------------
    whole_sc_node = 7.0;
    bnd_sc_node1  = whole_sc_node;

    {
        Index key = model_->Database().StorageKey( "bnd sc node 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 7.0, TOLERANCE );
    }

    // ------------------------------------------------------------------
    //  Test 2: whole-model NODE → SplitBoundary NODE.
    // ------------------------------------------------------------------
    whole_sc_node = 11.0;
    spb_sc_node1  = whole_sc_node;

    double readback = 0.0;
    spb_sc_node1.Apply( [&readback]( double& x ) { readback = x; } );
    _equal( readback, 11.0, TOLERANCE );

    // ------------------------------------------------------------------
    //  Test 3: whole-model NODE += Boundary NODE.
    // ------------------------------------------------------------------
    whole_sc_node = 3.0;
    bnd_sc_node1  = 4.0;
    bnd_sc_node1 += whole_sc_node;

    {
        Index key = model_->Database().StorageKey( "bnd sc node 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.NodesBegin()) )->Read( key, sc );
        _equal( sc(), 7.0, TOLERANCE );
    }

    // ------------------------------------------------------------------
    //  Tests 4-6: rejected assignments verified via sentinel values.
    //
    //  In debug builds, csmp::Exception prints and calls getchar()
    //  in its constructor before being thrown, which interferes with
    //  automated try/catch testing. We therefore verify rejection by
    //  checking that the destination variable was not modified.
    //
    //  The try/catch with catch(...) suppresses the exception so the
    //  test continues; the _equal check on the sentinel value confirms
    //  the assignment was rejected.
    // ------------------------------------------------------------------

    // Test 4: whole-model ELEMENT → Boundary FACE must be rejected.
    // FACE placement exists exclusively on Boundary subdomains.
    // ELEMENT placement exists exclusively on Region subdomains.
    {
        bnd_sc_face1  = 99.0;   // sentinel
        whole_sc_elem = 7.0;
        try { bnd_sc_face1 = whole_sc_elem; } catch ( ... ) {}
        Index key = model_->Database().StorageKey( "bnd sc face 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 99.0,
                numeric_limits<double>::epsilon() * 99.0 );
        if ( verbose_ )
            cout << "  Test 4 (ELEMENT->FACE rejected): passed" << endl;
    }

    // Test 5: Boundary NODE → whole-model NODE must be rejected.
    // Assigning from a sub-domain into the whole model would only
    // partially overwrite the model.
    {
        whole_sc_node = 55.0;   // sentinel
        bnd_sc_node1  = 3.0;
        try { whole_sc_node = bnd_sc_node1; } catch ( ... ) {}
        Index key = model_->Database().StorageKey( "whole sc node" );
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 55.0,
                numeric_limits<double>::epsilon() * 55.0 );
        if ( verbose_ )
            cout << "  Test 5 (subdomain->Model rejected): passed" << endl;
    }

    // Test 6: whole-model NODE → Boundary FACE must be rejected.
    // Incompatible placements even though source is whole-model.
    {
        bnd_sc_face1  = 99.0;   // sentinel
        whole_sc_node = 7.0;
        try { bnd_sc_face1 = whole_sc_node; } catch ( ... ) {}
        Index key = model_->Database().StorageKey( "bnd sc face 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 99.0,
                numeric_limits<double>::epsilon() * 99.0 );
        if ( verbose_ )
            cout << "  Test 6 (NODE->FACE rejected): passed" << endl;
    }

    if ( verbose_ ) cout << "  Model-to-subdomain assignment: passed" << endl;
}



// ============================================================================
//  TestApplyMethod
// ============================================================================

void PropertyHandle_Test::TestApplyMethod()
{
    ResetVariables();

    if ( verbose_ ) cout << "\n--- Apply() method ---" << endl;

    ScalarVariable sc;

    // Region ELEMENT: apply x = x * 2 + 1
    reg_sc_elem1.Apply( []( double& x ) { x = x * 2.0 + 1.0; } );
    // expected: 1.5 * 2 + 1 = 4.0
    {
        Index key = model_->Database().StorageKey( "reg sc elem 1" );
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 4.0, numeric_limits<double>::epsilon() * 4.0 );
    }

    // Boundary FACE: apply x = x + 10
    bnd_sc_face1.Apply( []( double& x ) { x += 10.0; } );
    // expected: 3.0 + 10 = 13.0
    {
        Index key = model_->Database().StorageKey( "bnd sc face 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 13.0, numeric_limits<double>::epsilon() * 13.0 );
    }

    // SplitBoundary INTER_FACE: apply x = x * 3
    spb_sc_iface1.Apply( []( double& x ) { x *= 3.0; } );
    // expected: 3.0 * 3 = 9.0
    {
        Index key = model_->Database().StorageKey( "spb sc iface 1" );
        const auto& spb = model_->SplitBoundary( "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT" );
        ( *(spb.CellsBegin()) )->Read( key, sc );
        _equal( sc(), 9.0, numeric_limits<double>::epsilon() * 9.0 );
    }

    // SplitBoundary NODE (deduplicated): apply x = x + 5
    spb_sc_node1.Apply( []( double& x ) { x += 5.0; } );
    // expected: 1.0 + 5 = 6.0
    double readback = 0.0;
    spb_sc_node1.Apply( [&readback]( double& x ) { readback = x; } );
    _equal( readback, 6.0, numeric_limits<double>::epsilon() * 6.0 );

    if ( verbose_ ) cout << "  Apply() method: passed" << endl;
}

// ============================================================================
//  TestMathMethods
// ============================================================================

void PropertyHandle_Test::TestMathMethods()
{
    if ( verbose_ ) cout << "\n--- Math methods (Region ELEMENT) ---" << endl;

    ScalarVariable sc;
    Index key = model_->Database().StorageKey( "reg sc elem 1" );

    // Squared: 1.5^2 = 2.25
    reg_sc_elem1 = 1.5;
    reg_sc_elem1.Squared();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 2.25, numeric_limits<double>::epsilon() * 2.25 );

    // Sqrt: sqrt(2.25) = 1.5
    reg_sc_elem1.Sqrt();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 1.5, numeric_limits<double>::epsilon() * 1.5 );

    // Ln: ln(1.5)
    reg_sc_elem1.Ln();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), std::log(1.5),
            numeric_limits<double>::epsilon() * std::log(1.5) );

    // Exp: exp(ln(1.5)) = 1.5
    reg_sc_elem1.Exp();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 1.5, numeric_limits<double>::epsilon() * 1.5 );

    // Log10: log10(1.5)
    reg_sc_elem1.Log10();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), std::log10(1.5),
            numeric_limits<double>::epsilon() * std::log10(1.5) );

    // Pow: 2^3 = 8
    reg_sc_elem1 = 2.0;
    reg_sc_elem1.Pow( 3.0 );
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 8.0, numeric_limits<double>::epsilon() * 8.0 );

    // Abs: abs(-8) = 8
    reg_sc_elem1 = -8.0;
    reg_sc_elem1.Abs();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 8.0, numeric_limits<double>::epsilon() * 8.0 );

    // ZapNAN: NaN -> 99
    reg_sc_elem1 = numeric_limits<double>::quiet_NaN();
    reg_sc_elem1.ZapNAN( 99.0 );
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 99.0, numeric_limits<double>::epsilon() * 99.0 );

    // Clamp: 99 clamped to [0, 10] = 10
    reg_sc_elem1.Clamp( 0.0, 10.0 );
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 10.0, numeric_limits<double>::epsilon() * 10.0 );

    // Clamp: -5 clamped to [0, 10] = 0
    reg_sc_elem1 = -5.0;
    reg_sc_elem1.Clamp( 0.0, 10.0 );
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 0.0, numeric_limits<double>::epsilon() );

    // Sin: sin(30 deg) = 0.5
    reg_sc_elem1 = 30.0;
    reg_sc_elem1.Sin();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 0.5, numeric_limits<double>::epsilon() * 10.0 );

    // Cos: cos(60 deg) = 0.5
    reg_sc_elem1 = 60.0;
    reg_sc_elem1.Cos();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 0.5, numeric_limits<double>::epsilon() * 10.0 );

    // Tan: tan(45 deg) = 1.0
    reg_sc_elem1 = 45.0;
    reg_sc_elem1.Tan();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 1.0, numeric_limits<double>::epsilon() * 10.0 );

    // Acos: acos(0.5) = 60 deg
    reg_sc_elem1 = 0.5;
    reg_sc_elem1.Acos();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 60.0, numeric_limits<double>::epsilon() * 60.0 );

    // Asin: asin(0.5) = 30 deg
    reg_sc_elem1 = 0.5;
    reg_sc_elem1.Asin();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 30.0, numeric_limits<double>::epsilon() * 30.0 );

    // Atan: atan(1.0) = 45 deg
    reg_sc_elem1 = 1.0;
    reg_sc_elem1.Atan();
    ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
    _equal( sc(), 45.0, numeric_limits<double>::epsilon() * 45.0 );

    if ( verbose_ ) cout << "  Math methods: passed" << endl;

    // ------------------------------------------------------------------
    //  Verify math methods also work on Boundary and SplitBoundary.
    //  Use a single representative method (Sqrt) on each.
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  Math methods on Boundary and SplitBoundary" << endl;

    // Boundary FACE: sqrt(4.0) = 2.0
    bnd_sc_face1 = 4.0;
    bnd_sc_face1.Sqrt();
    {
        Index bkey = model_->Database().StorageKey( "bnd sc face 1" );
        const auto& bnd = model_->Boundary( "LEFT" );
        ( *(bnd.CellsBegin()) )->Read( bkey, sc );
        _equal( sc(), 2.0, numeric_limits<double>::epsilon() * 2.0 );
    }

    // SplitBoundary INTER_FACE: sqrt(9.0) = 3.0
    spb_sc_iface1 = 9.0;
    spb_sc_iface1.Sqrt();
    {
        Index skey = model_->Database().StorageKey( "spb sc iface 1" );
        const auto& spb = model_->SplitBoundary( "HALF_SPLITBOUNDARY0_MATRIX_RIGHT_MATRIX_LEFT" );
        ( *(spb.CellsBegin()) )->Read( skey, sc );
        _equal( sc(), 3.0, numeric_limits<double>::epsilon() * 3.0 );
    }

    // SplitBoundary SECTOR_IP: sqrt(16.0) = 4.0
    spb_sc_ifsec1 = 16.0;
    spb_sc_ifsec1.Sqrt();
    double readback = 0.0;
    spb_sc_ifsec1.Apply( [&readback]( double& x ) { readback = x; } );
    _equal( readback, 4.0, numeric_limits<double>::epsilon() * 4.0 );

    if ( verbose_ ) cout << "  Math methods on Boundary/SplitBoundary: passed" << endl;
}

// ============================================================================
//  TestArray treatment
// ============================================================================

void PropertyHandle_Test::TestArraySupport()
{
    if ( verbose_ ) cout << "\n--- Array variable support ---" << endl;

    const uint32_t array_size = 4U;

    // -----------------------------------------------------------------------
    //  ARRAY — construction, assignment, compound scalar operators,
    //  and custom ApplyArray lambda.
    //  All verification uses direct Read from the first cell to avoid
    //  accumulation errors from iterating all mesh points.
    // -----------------------------------------------------------------------
    if ( verbose_ ) cout << "  ARRAY" << endl;
    {
        PropertyHandle<3> ph_array( *model_, "test array var",
                                     ARRAY, ELEMENT, array_size );

        const Index key =
            model_->Database().StorageKey( "test array var" );
        ArrayVariable av;
        auto firstCell = model_->Region("Model").CellsBegin();

        // --- operator=( double ) ---
        ph_array = 2.0;
        ( *firstCell )->Read( key, av );
        _equal( av[0], 2.0, numeric_limits<double>::epsilon() * 2.0 );
        _equal( av[1], 2.0, numeric_limits<double>::epsilon() * 2.0 );
        _equal( av[2], 2.0, numeric_limits<double>::epsilon() * 2.0 );
        _equal( av[3], 2.0, numeric_limits<double>::epsilon() * 2.0 );

        // --- operator=( const ArrayVariable& ) ---
        ArrayVariable init( array_size, 5.0 );
        ph_array = init;
        ( *firstCell )->Read( key, av );
        _equal( av[0], 5.0, numeric_limits<double>::epsilon() * 5.0 );

        // --- operator+=( double ) ---
        ph_array += 1.0;
        ( *firstCell )->Read( key, av );
        // expected: 5.0 + 1.0 = 6.0
        _equal( av[0], 6.0, numeric_limits<double>::epsilon() * 6.0 );

        // --- operator-=( double ) ---
        ph_array -= 2.0;
        ( *firstCell )->Read( key, av );
        // expected: 6.0 - 2.0 = 4.0
        _equal( av[0], 4.0, numeric_limits<double>::epsilon() * 4.0 );

        // --- operator*=( double ) ---
        ph_array *= 3.0;
        ( *firstCell )->Read( key, av );
        // expected: 4.0 * 3.0 = 12.0
        _equal( av[0], 12.0, numeric_limits<double>::epsilon() * 12.0 );

        // --- operator/=( double ) ---
        ph_array /= 4.0;
        ( *firstCell )->Read( key, av );
        // expected: 12.0 / 4.0 = 3.0
        _equal( av[0], 3.0, numeric_limits<double>::epsilon() * 3.0 );

        // --- ApplyArray with custom lambda ---
        // Set element i to i+1: [1, 2, 3, 4]
        ph_array.ApplyArray( []( ArrayVariable& a )
        {
            for ( uint32_t i = 0; i < a.Size(); ++i )
                a(i) = static_cast<double>(i) + 1.0;
        });
        ( *firstCell )->Read( key, av );
        _equal( av[0], 1.0, numeric_limits<double>::epsilon() );
        _equal( av[1], 2.0, numeric_limits<double>::epsilon() * 2.0 );
        _equal( av[2], 3.0, numeric_limits<double>::epsilon() * 3.0 );
        _equal( av[3], 4.0, numeric_limits<double>::epsilon() * 4.0 );

        // --- operator+=( const PropertyHandle3& ) same placement ---
        // ph_array = [1,2,3,4]; add itself: expected [2,4,6,8]
        ph_array += ph_array;
        ( *firstCell )->Read( key, av );
        _equal( av[0], 2.0, TOLERANCE );
        _equal( av[1], 4.0, TOLERANCE );
        _equal( av[2], 6.0, TOLERANCE );
        _equal( av[3], 8.0, TOLERANCE );

        // --- operator-=( const PropertyHandle3& ) ---
        // [2,4,6,8] - [2,4,6,8] = [0,0,0,0]
        ph_array -= ph_array;
        ( *firstCell )->Read( key, av );
        _equal( av[0], 0.0, TOLERANCE );

    // --- Named math methods on ARRAY ---
    if ( verbose_ ) cout << "  Named math methods on ARRAY" << endl;
    {
        PropertyHandle<3> ph( *model_, "test array math",
                               ARRAY, ELEMENT, array_size );
        const Index keya = model_->Database().StorageKey( "test array math" );
        ArrayVariable av;
        auto firstCell = model_->Region("Model").CellsBegin();

        // Squared: 2^2 = 4
        ph = 2.0;
        ph.Squared();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 4.0, numeric_limits<double>::epsilon() * 4.0 );

        // Sqrt: sqrt(4) = 2
        ph.Sqrt();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 2.0, numeric_limits<double>::epsilon() * 2.0 );

        // Ln: ln(2)
        ph.Ln();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], std::log(2.0),
                numeric_limits<double>::epsilon() * std::log(2.0) );

        // Exp: exp(ln(2)) = 2
        ph.Exp();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 2.0, numeric_limits<double>::epsilon() * 2.0 );

        // Log10: log10(2)
        ph.Log10();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], std::log10(2.0),
                numeric_limits<double>::epsilon() * std::log10(2.0) );

        // Pow: 2^3 = 8
        ph = 2.0;
        ph.Pow( 3.0 );
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 8.0, numeric_limits<double>::epsilon() * 8.0 );

        // Abs: abs(-8) = 8
        ph = -8.0;
        ph.Abs();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 8.0, numeric_limits<double>::epsilon() * 8.0 );

        // ZapNAN: NaN -> 99
        ph = numeric_limits<double>::quiet_NaN();
        ph.ZapNAN( 99.0 );
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 99.0, numeric_limits<double>::epsilon() * 99.0 );

        // Clamp: 99 -> [0,10] = 10
        ph.Clamp( 0.0, 10.0 );
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 10.0, numeric_limits<double>::epsilon() * 10.0 );

        // Sin: sin(30 deg) = 0.5
        ph = 30.0;
        ph.Sin();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 0.5, numeric_limits<double>::epsilon() * 10.0 );

        // Cos: cos(60 deg) = 0.5
        ph = 60.0;
        ph.Cos();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 0.5, numeric_limits<double>::epsilon() * 10.0 );

        // Tan: tan(45 deg) = 1
        ph = 45.0;
        ph.Tan();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 1.0, numeric_limits<double>::epsilon() * 10.0 );

        // Acos: acos(0.5) = 60 deg
        ph = 0.5;
        ph.Acos();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 60.0, numeric_limits<double>::epsilon() * 60.0 );

        // Asin: asin(0.5) = 30 deg
        ph = 0.5;
        ph.Asin();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 30.0, numeric_limits<double>::epsilon() * 30.0 );

        // Atan: atan(1) = 45 deg
        ph = 1.0;
        ph.Atan();
        ( *firstCell )->Read( keya, av );
        _equal( av[0], 45.0, numeric_limits<double>::epsilon() * 45.0 );

        if ( verbose_ )
            cout << "    Named math methods on ARRAY: passed" << endl;
    }

    // --- Named math methods on FLAGGEDARRAY ---
    // Verify that DIRICH-flagged elements are not modified.
    if ( verbose_ ) cout << "  Named math methods on FLAGGEDARRAY" << endl;
    {
        PropertyHandle<3> ph( *model_, "test farray math",
                               FLAGGEDARRAY, ELEMENT, array_size );
        const Index key1 = model_->Database().StorageKey( "test farray math" );
        FlaggedArrayVariable fav;
        auto firstCell1 = model_->Region("Model").CellsBegin();

        // Initialise: all elements = 4.0, element 0 = DIRICH.
        ph = 4.0;
        ph.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
        {
            f.Flag( 0, DIRICH );
        });

        // Sqrt: sqrt(4) = 2 for ANY elements; DIRICH element unchanged.
        ph.Sqrt();
        ( *firstCell1 )->Read( key1, fav );
        // element 0 (DIRICH): must remain 4.0
        _equal( fav[0], 4.0, numeric_limits<double>::epsilon() * 4.0 );
        // element 1 (ANY): sqrt(4) = 2.0
        _equal( fav[1], 2.0, numeric_limits<double>::epsilon() * 2.0 );

        // Squared: 2^2 = 4 for ANY; DIRICH unchanged.
        ph.Squared();
        ( *firstCell )->Read( key1, fav );
        _equal( fav[0], 4.0, numeric_limits<double>::epsilon() * 4.0 );
        _equal( fav[1], 4.0, numeric_limits<double>::epsilon() * 4.0 );

        // Pow: 4^2 = 16 for ANY; DIRICH unchanged.
        ph.Pow( 2.0 );
        ( *firstCell1 )->Read( key1, fav );
        _equal( fav[0], 4.0,  numeric_limits<double>::epsilon() * 4.0 );
        _equal( fav[1], 16.0, numeric_limits<double>::epsilon() * 16.0 );

        if ( verbose_ )
            cout << "    Named math methods on FLAGGEDARRAY: passed" << endl;
    }

        if ( verbose_ )
            cout << "    ARRAY: passed" << endl;
    }

    // -----------------------------------------------------------------------
    //  FLAGGEDARRAY — construction, assignment, per-element flag checking,
    //  compound scalar operators, and ApplyToFlaggedElements helper.
    // -----------------------------------------------------------------------
    if ( verbose_ ) cout << "  FLAGGEDARRAY" << endl;
    {
        PropertyHandle<3> ph_farray( *model_, "test flagged array var",
                                      FLAGGEDARRAY, ELEMENT, array_size );

        const Index key2 =
            model_->Database().StorageKey( "test flagged array var" );
        FlaggedArrayVariable fav;
        auto firstCell = model_->Region("Model").CellsBegin();

        // --- operator=( double ) ---
        ph_farray = 3.0;
        ( *firstCell )->Read( key2, fav );
        _equal( fav[0], 3.0, numeric_limits<double>::epsilon() * 3.0 );
        _equal( fav[1], 3.0, numeric_limits<double>::epsilon() * 3.0 );

        // --- operator=( const FlaggedArrayVariable& ) ---
        FlaggedArrayVariable init( array_size, 7.0, ANY );
        ph_farray = init;
        ( *firstCell )->Read( key2, fav );
        _equal( fav[0], 7.0, numeric_limits<double>::epsilon() * 7.0 );

        // --- Per-element flag: fix element 0, leave others as ANY ---
        ph_farray.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
        {
            f.Flag( 0, DIRICH );
        });

        // --- operator+=( double ) respects per-element flags ---
        ph_farray += 1.0;
        ( *firstCell )->Read( key2, fav );
        // element 0 is FIXED — must remain 7.0
        _equal( fav[0], 7.0, numeric_limits<double>::epsilon() * 7.0 );
        // element 1 is ANY — expected: 7.0 + 1.0 = 8.0
        _equal( fav[1], 8.0, numeric_limits<double>::epsilon() * 8.0 );

        // --- operator-=( double ) ---
        ph_farray -= 1.0;
        ( *firstCell )->Read( key2, fav );
        // element 0 (FIXED): unchanged at 7.0
        _equal( fav[0], 7.0, numeric_limits<double>::epsilon() * 7.0 );
        // element 1 (ANY): 8.0 - 1.0 = 7.0
        _equal( fav[1], 7.0, numeric_limits<double>::epsilon() * 7.0 );

        // --- operator*=( double ) ---
        ph_farray *= 2.0;
        ( *firstCell )->Read( key2, fav );
        // element 0 (FIXED): unchanged at 7.0
        _equal( fav[0], 7.0, numeric_limits<double>::epsilon() * 7.0 );
        // element 1 (ANY): 7.0 * 2.0 = 14.0
        _equal( fav[1], 14.0, numeric_limits<double>::epsilon() * 14.0 );

        // --- operator/=( double ) ---
        ph_farray /= 2.0;
        ( *firstCell )->Read( key2, fav );
        // element 0 (FIXED): unchanged at 7.0
        _equal( fav[0], 7.0, numeric_limits<double>::epsilon() * 7.0 );
        // element 1 (ANY): 14.0 / 2.0 = 7.0
        _equal( fav[1], 7.0, numeric_limits<double>::epsilon() * 7.0 );

        // --- ApplyToFlaggedElements helper ---
        // Multiply only ANY elements by 3.
        ph_farray.ApplyFlaggedArray(
            [this]( FlaggedArrayVariable& f )
            {
                PropertyHandle<3>::ApplyToFlaggedElements(
                    f, ANY,
                    []( double& x ) { x *= 3.0; } );
            });
        ( *firstCell )->Read( key2, fav );
        // element 0 (FIXED): unchanged at 7.0
        _equal( fav[0], 7.0, numeric_limits<double>::epsilon() * 7.0 );
        // element 1 (ANY): 7.0 * 3.0 = 21.0
        _equal( fav[1], 21.0, numeric_limits<double>::epsilon() * 21.0 );

        // --- ApplyFlaggedArray with custom lambda ---
        ph_farray.ApplyFlaggedArray( []( FlaggedArrayVariable& f )
        {
            // Set all elements to their index value, ignoring flags.
            for ( uint32_t i = 0; i < f.Size(); ++i )
                f(i) = static_cast<double>(i);
            // elements: [0, 1, 2, 3]
        });
        ( *firstCell )->Read( key2, fav );
        _equal( fav[0], 0.0, numeric_limits<double>::epsilon() );
        _equal( fav[1], 1.0, numeric_limits<double>::epsilon() );
        _equal( fav[2], 2.0, numeric_limits<double>::epsilon() * 2.0 );
        _equal( fav[3], 3.0, numeric_limits<double>::epsilon() * 3.0 );

        if ( verbose_ )
            cout << "    FLAGGEDARRAY: passed" << endl;
    }

    // -----------------------------------------------------------------------
    //  Cross-type rejection: ARRAY <-> SCALAR must not modify destination.
    //  Verified via sentinel values since csmp::Exception constructor
    //  may pause in debug builds.
    // -----------------------------------------------------------------------
    if ( verbose_ ) cout << "  Cross-type rejection" << endl;
    {
        PropertyHandle<3> ph_array ( *model_, "test array reject",
                                      ARRAY,  ELEMENT, array_size );
        PropertyHandle<3> ph_scalar( *model_, "test scalar reject",
                                      SCALAR, ELEMENT );

        const Index akey =
            model_->Database().StorageKey( "test array reject" );
        const Index skey =
            model_->Database().StorageKey( "test scalar reject" );

        ph_array  = 1.0;   // sentinel
        ph_scalar = 99.0;  // sentinel

        // SCALAR = ARRAY must be rejected — ph_scalar must remain 99.0.
        try { ph_scalar = ph_array; } catch ( ... ) {}
        {
            ScalarVariable sc;
            ( *(model_->Region("Model").CellsBegin()) )->Read( skey, sc );
            _equal( sc(), 99.0, numeric_limits<double>::epsilon() * 99.0 );
        }

        // ARRAY = SCALAR must be rejected — ph_array element 0 must remain 1.0.
        try { ph_array = ph_scalar; } catch ( ... ) {}
        {
            ArrayVariable av;
            ( *(model_->Region("Model").CellsBegin()) )->Read( akey, av );
            _equal( av[0], 1.0, numeric_limits<double>::epsilon() );
        }

        if ( verbose_ )
            cout << "    Cross-type rejection: passed" << endl;
    }

    // -----------------------------------------------------------------------
    //  BOUNDARY placement rejection still applies.
    // -----------------------------------------------------------------------
    if ( verbose_ ) cout << "  BOUNDARY placement rejection" << endl;
    {
        bool threw = false;
        try
        {
            PropertyHandle<3> bad( *model_, "bad placement",
                                    SCALAR, BOUNDARY );
        }
        catch ( const csmp::Exception& ) { threw = true; }
        catch ( ... )                    { threw = true; }
        _equal( threw, true, 0.5 );

        if ( verbose_ )
            cout << "    BOUNDARY placement rejection: passed" << endl;
    }

    if ( verbose_ ) cout << "  Array variable support: passed" << endl;
}





void PropertyHandle_Test::TestRegionCrossPlacementFull()
{
    if ( verbose_ ) cout << "\n--- Region cross-placement full ---" << endl;

    ScalarVariable    sc;
    VectorVariable<3> vc;
    TensorVariable<3> ts;

    // ------------------------------------------------------------------
    //  SCALAR — all four operators, all six placement combinations
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  SCALAR cross-placement" << endl;

    // += 
    ResetVariables();
    { // NODE += ELEMENT
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 += reg_sc_elem1;
        // expected: 4.0 + 1.5 = 5.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 5.5, TOLERANCE );
    }
    ResetVariables();
    { // NODE += IP
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 += reg_sc_ip1;
        // expected: 4.0 + 5.0 = 9.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT += NODE
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 += reg_sc_node1;
        // expected: 2.0 + 3.0 = 5.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 5.0, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT += IP
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 += reg_sc_ip1;
        // expected: 2.0 + 5.0 = 7.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 7.0, TOLERANCE );
    }
    ResetVariables();
    { // IP += NODE
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 += reg_sc_node1;
        // expected: 6.0 + 3.0 = 9.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }
    ResetVariables();
    { // IP += ELEMENT
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 += reg_sc_elem1;
        // expected: 6.0 + 1.5 = 7.5
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 7.5, TOLERANCE );
    }

    // -=
    ResetVariables();
    { // NODE -= ELEMENT
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 -= reg_sc_elem1;
        // expected: 4.0 - 1.5 = 2.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 2.5, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT -= NODE
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 -= reg_sc_node1;
        // expected: 2.0 - 3.0 = -1.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), -1.0, TOLERANCE );
    }
    ResetVariables();
    { // IP -= NODE
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 -= reg_sc_node1;
        // expected: 6.0 - 3.0 = 3.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 3.0, TOLERANCE );
    }

    // *=
    ResetVariables();
    { // NODE *= ELEMENT
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 *= reg_sc_elem1;
        // expected: 4.0 * 1.5 = 6.0
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 6.0, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT *= NODE
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 *= reg_sc_node1;
        // expected: 2.0 * 3.0 = 6.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 6.0, TOLERANCE );
    }
    ResetVariables();
    { // IP *= ELEMENT
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 *= reg_sc_elem1;
        // expected: 6.0 * 1.5 = 9.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 9.0, TOLERANCE );
    }

    // /=
    ResetVariables();
    { // NODE /= ELEMENT
        Index key = model_->Database().StorageKey( "reg sc node 2" );
        reg_sc_node2 /= reg_sc_elem1;
        // expected: 4.0 / 1.5 = 8/3
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, sc );
        _equal( sc(), 8.0/3.0, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT /= NODE
        Index key = model_->Database().StorageKey( "reg sc elem 2" );
        reg_sc_elem2 /= reg_sc_node1;
        // expected: 2.0 / 3.0 = 2/3
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, sc );
        _equal( sc(), 2.0/3.0, TOLERANCE );
    }
    ResetVariables();
    { // IP /= ELEMENT
        Index key = model_->Database().StorageKey( "reg sc ip 2" );
        reg_sc_ip2 /= reg_sc_elem1;
        // expected: 6.0 / 1.5 = 4.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( 0U, key, sc );
        _equal( sc(), 4.0, TOLERANCE );
    }

    // ------------------------------------------------------------------
    //  VECTOR — representative cross-placement combinations
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  VECTOR cross-placement" << endl;

    ResetVariables();
    { // NODE += ELEMENT
        Index key = model_->Database().StorageKey( "reg vc node 2" );
        reg_vc_node2 += reg_vc_elem1;
        // expected: 4.0 + 1.5 = 5.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, vc );
        _equal( vc(0), 5.5, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT += NODE
        Index key = model_->Database().StorageKey( "reg vc elem 2" );
        reg_vc_elem2 += reg_vc_node1;
        // expected: 2.0 + 3.0 = 5.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 5.0, TOLERANCE );
    }
    ResetVariables();
    { // NODE -= ELEMENT
        Index key = model_->Database().StorageKey( "reg vc node 2" );
        reg_vc_node2 -= reg_vc_elem1;
        // expected: 4.0 - 1.5 = 2.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, vc );
        _equal( vc(0), 2.5, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT *= NODE (component-wise)
        Index key = model_->Database().StorageKey( "reg vc elem 2" );
        reg_vc_elem2 *= reg_vc_node1;
        // expected: 2.0 * 3.0 = 6.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, vc );
        _equal( vc(0), 6.0, TOLERANCE );
    }
    ResetVariables();
    { // NODE /= ELEMENT
        Index key = model_->Database().StorageKey( "reg vc node 2" );
        reg_vc_node2 /= reg_vc_elem1;
        // expected: 4.0 / 1.5 = 8/3
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, vc );
        _equal( vc(0), 8.0/3.0, TOLERANCE );
    }

    // ------------------------------------------------------------------
    //  TENSOR — representative cross-placement combinations
    //  Tensor *= tensor is matrix product: result(0,0) = 3 * a * b
    //  for uniform 3x3 matrices with all entries a and b.
    // ------------------------------------------------------------------
    if ( verbose_ ) cout << "  TENSOR cross-placement" << endl;

    ResetVariables();
    { // NODE += ELEMENT
        Index key = model_->Database().StorageKey( "reg ts node 2" );
        reg_ts_node2 += reg_ts_elem1;
        // expected: 4.0 + 1.5 = 5.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, ts );
        _equal( ts(0,0), 5.5, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT += NODE
        Index key = model_->Database().StorageKey( "reg ts elem 2" );
        reg_ts_elem2 += reg_ts_node1;
        // expected: 2.0 + 3.0 = 5.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 5.0, TOLERANCE );
    }
    ResetVariables();
    { // NODE -= ELEMENT
        Index key = model_->Database().StorageKey( "reg ts node 2" );
        reg_ts_node2 -= reg_ts_elem1;
        // expected: 4.0 - 1.5 = 2.5
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, ts );
        _equal( ts(0,0), 2.5, TOLERANCE );
    }
    ResetVariables();
    { // ELEMENT *= NODE (matrix product)
        Index key = model_->Database().StorageKey( "reg ts elem 2" );
        reg_ts_elem2 *= reg_ts_node1;
        // expected: 3 * 2.0 * 3.0 = 18.0
        ( *(model_->Region("Model").CellsBegin()) )->Read( key, ts );
        _equal( ts(0,0), 3.0 * 2.0 * 3.0, TOLERANCE );
    }
    ResetVariables();
    { // NODE /= double (singular matrix issue — use scalar divisor)
        reg_ts_node2 = 6.0;
        reg_ts_node2 /= 2.0;
        // expected: 6.0 / 2.0 = 3.0
        Index key = model_->Database().StorageKey( "reg ts node 2" );
        ( *(model_->Region("Model").NodesBegin()) )->Read( key, ts );
        _equal( ts(0,0), 3.0, TOLERANCE );
    }

    if ( verbose_ ) cout << "  Region cross-placement full: passed" << endl;
}





void PropertyHandle_Test::TestFlaggedArrayFacetIP()
{
    if ( verbose_ )
        cout << "\n--- FlaggedArray on FACET_INTEGRATION_POINT ---" << endl;

    const uint32_t array_size = 10U;

    // ------------------------------------------------------------------
    //  Create a FLAGGEDARRAY on FACET_INTEGRATION_POINT and a matching
    //  ARRAY handle to receive the output.
    // ------------------------------------------------------------------
    PropertyHandle<3> ph_farray( *model_,
                                   "facet farray test",
                                   FLAGGEDARRAY,
                                   FACET_INTEGRATION_POINT,
                                   array_size );

    PropertyHandle<3> ph_array( *model_,
                                  "facet array test",
                                  ARRAY,
                                  FACET_INTEGRATION_POINT,
                                  array_size );

    const Index farray_key =
        model_->Database().StorageKey( "facet farray test" );
    const Index array_key =
        model_->Database().StorageKey( "facet array test" );

    auto firstCell = model_->Region("Model").CellsBegin();
    const uint32_t nFacets = ( *firstCell )->Facets();

    // ------------------------------------------------------------------
    //  Step 1: Initialise all elements to 0.0 with flag ANY.
    // ------------------------------------------------------------------
    ph_farray = 0.0;

    // Verify initialisation on first cell, all facets.
    {
        FlaggedArrayVariable fav;
        for ( uint32_t fac = 0; fac < nFacets; ++fac )
        {
            ( *firstCell )->Read( fac, 0U, farray_key, fav );
            for ( uint32_t i = 0; i < array_size; ++i )
            {
                _equal( fav[i], 0.0,
                        numeric_limits<double>::epsilon() );
                _equal( static_cast<int>( fav.Flag(i) ),
                        static_cast<int>( ANY ), 0.5 );
            }
        }
        if ( verbose_ )
            cout << "  Step 1 (initialisation): passed" << endl;
    }

    // ------------------------------------------------------------------
    //  Step 2: Set elements at prime indices (1, 3, 5, 7) to DIRICH
    //  on every facet of every element.
    // ------------------------------------------------------------------
    ph_farray.ApplyFlaggedArray( []( FlaggedArrayVariable& fav )
    {
        fav.Flag( 1, DIRICH );
        fav.Flag( 3, DIRICH );
        fav.Flag( 5, DIRICH );
        fav.Flag( 7, DIRICH );
    });

    // Verify flags on first cell, all facets.
    {
        FlaggedArrayVariable fav;
        for ( uint32_t fac = 0; fac < nFacets; ++fac )
        {
            ( *firstCell )->Read( fac, 0U, farray_key, fav );

            // Non-prime indices: ANY
            _equal( static_cast<int>( fav.Flag(0) ),
                    static_cast<int>( ANY ), 0.5 );
            _equal( static_cast<int>( fav.Flag(2) ),
                    static_cast<int>( ANY ), 0.5 );
            _equal( static_cast<int>( fav.Flag(4) ),
                    static_cast<int>( ANY ), 0.5 );
            _equal( static_cast<int>( fav.Flag(6) ),
                    static_cast<int>( ANY ), 0.5 );
            _equal( static_cast<int>( fav.Flag(8) ),
                    static_cast<int>( ANY ), 0.5 );
            _equal( static_cast<int>( fav.Flag(9) ),
                    static_cast<int>( ANY ), 0.5 );

            // Prime indices: DIRICH
            _equal( static_cast<int>( fav.Flag(1) ),
                    static_cast<int>( DIRICH ), 0.5 );
            _equal( static_cast<int>( fav.Flag(3) ),
                    static_cast<int>( DIRICH ), 0.5 );
            _equal( static_cast<int>( fav.Flag(5) ),
                    static_cast<int>( DIRICH ), 0.5 );
            _equal( static_cast<int>( fav.Flag(7) ),
                    static_cast<int>( DIRICH ), 0.5 );
        }
        if ( verbose_ )
            cout << "  Step 2 (flag assignment): passed" << endl;
    }

    // ------------------------------------------------------------------
    //  Step 3: Overwrite with -1.0 using OutputCondition = ANY.
    //  Only elements with flag ANY should be modified.
    //  Elements with flag DIRICH must remain 0.0.
    // ------------------------------------------------------------------
    ph_farray.OutputCondition( ANY );
    ph_farray = -1.0;

    // Verify on every facet of the first cell.
    {
        FlaggedArrayVariable fav;
        for ( uint32_t fac = 0; fac < nFacets; ++fac )
        {
            ( *firstCell )->Read( fac, 0U, farray_key, fav );

            // Non-prime indices (ANY): must be -1.0
            _equal( fav[0], -1.0, numeric_limits<double>::epsilon() );
            _equal( fav[2], -1.0, numeric_limits<double>::epsilon() );
            _equal( fav[4], -1.0, numeric_limits<double>::epsilon() );
            _equal( fav[6], -1.0, numeric_limits<double>::epsilon() );
            _equal( fav[8], -1.0, numeric_limits<double>::epsilon() );
            _equal( fav[9], -1.0, numeric_limits<double>::epsilon() );

            // Prime indices (DIRICH): must remain 0.0
            _equal( fav[1], 0.0, numeric_limits<double>::epsilon() );
            _equal( fav[3], 0.0, numeric_limits<double>::epsilon() );
            _equal( fav[5], 0.0, numeric_limits<double>::epsilon() );
            _equal( fav[7], 0.0, numeric_limits<double>::epsilon() );
        }
        if ( verbose_ )
            cout << "  Step 3 (selective overwrite with -1.0): passed"
                 << endl;
    }

    // ------------------------------------------------------------------
    //  Step 4: Cross-type assignment ARRAY = FLAGGEDARRAY.
    //  Uses CopyValuesOnly at each point — all values transferred
    //  regardless of per-element flags. Flags are not present in ARRAY.
    //  Expected: non-prime positions = -1.0, prime positions = 0.0.
    // ------------------------------------------------------------------
    ph_array = ph_farray;

    // Verify on every facet of the first cell.
    {
        ArrayVariable        av;
        FlaggedArrayVariable fav;

        for ( uint32_t fac = 0; fac < nFacets; ++fac )
        {
            ( *firstCell )->Read( fac, 0U, farray_key, fav );
            ( *firstCell )->Read( fac, 0U, array_key,  av  );

            // Every element of the ARRAY must match the FLAGGEDARRAY.
            for ( uint32_t i = 0; i < array_size; ++i )
                _equal( av[i], fav[i],
                        numeric_limits<double>::epsilon()
                        * std::max( 1.0, std::abs( fav[i] ) ) );

            // Explicitly verify key values.
            // Non-prime positions: -1.0
            _equal( av[0], -1.0, numeric_limits<double>::epsilon() );
            _equal( av[2], -1.0, numeric_limits<double>::epsilon() );
            _equal( av[4], -1.0, numeric_limits<double>::epsilon() );
            _equal( av[6], -1.0, numeric_limits<double>::epsilon() );
            _equal( av[8], -1.0, numeric_limits<double>::epsilon() );
            _equal( av[9], -1.0, numeric_limits<double>::epsilon() );

            // Prime positions: 0.0 (preserved by DIRICH protection)
            _equal( av[1], 0.0, numeric_limits<double>::epsilon() );
            _equal( av[3], 0.0, numeric_limits<double>::epsilon() );
            _equal( av[5], 0.0, numeric_limits<double>::epsilon() );
            _equal( av[7], 0.0, numeric_limits<double>::epsilon() );
        }
        if ( verbose_ )
            cout << "  Step 4 (ARRAY = FLAGGEDARRAY CopyValuesOnly): passed"
                 << endl;
    }

    // ------------------------------------------------------------------
    //  Step 5: Verify that the ARRAY has a single ANY flag — ArrayVariable
    //  has one flag for the whole array, not per-element flags.
    // ------------------------------------------------------------------
    {
        ArrayVariable av;
        ( *firstCell )->Read( 0U, 0U, array_key, av );
        _equal( static_cast<int>( av.Flag() ),
                static_cast<int>( ANY ), 0.5 );

        if ( verbose_ )
            cout << "  Step 5 (ARRAY has single ANY flag): passed" << endl;
    }

    if ( verbose_ )
        cout << "  FlaggedArray on FACET_INTEGRATION_POINT: passed" << endl;
}                





} // namespace csmp

