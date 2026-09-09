#ifndef CSMP_PROPERTY_HANDLE4_TEST_H
#define CSMP_PROPERTY_HANDLE4_TEST_H

#include "Test.h"
#include "PropertyHandle.h"
#include "Model.h"

namespace csmp {

/**
@brief Unit test for PropertyHandle.

@author 2024

@section overview Overview

Tests PropertyHandle across all supported subdomain types
(Region, Boundary, SplitBoundary) and all supported placements:

  Region:        NODE, ELEMENT, ELEMENT_INTEGRATION_POINT,
                 SECTOR_INTEGRATION_POINT, FACET_INTEGRATION_POINT

  Boundary:      NODE, FACE, FACE_INTEGRATION_POINT,
                 FACE_SECTOR_INTEGRATION_POINT,
                 FACE_FACET_INTEGRATION_POINT

  SplitBoundary: NODE (via InterFace FE nodes, deduplicated),
                 INTER_FACE, INTER_FACE_INTEGRATION_POINT,
                 INTER_FACE_SECTOR_INTEGRATION_POINT,
                 INTER_FACE_FACET_INTEGRATION_POINT

For each placement, all four compound assignment operators
(+=, -=, *=, /=) are tested for SCALAR, VECTOR, and TENSOR
variable types. Cross-placement interpolation is tested for
the primary placement pairs on each subdomain type.

Sub-domain isolation is verified: operations on a Boundary or
SplitBoundary handle must not modify variables on the Region.

@section mesh Test mesh

Uses the ANSYS mesh "BoxHalfs3D" which provides:
  - Region "Model"
  - Region "BoxHalf1", "BoxHalf2"
  - Boundary "LEFT_BOUNDARY", "RIGHT_BOUNDARY"
  - SplitBoundary "INTERNAL_SPLITBOUNDARY"

Variables are declared in "PropertyHandle_Test-variables.txt".

@section tolerance Tolerance

Default tolerance 5e-11 is appropriate for values in the range
[-3000, 3000] that arise from cross-placement interpolation.
*/
class PropertyHandle_Test : public Test
{
public:

    /**
    Constructs the test and initialises all PropertyHandle objects.

    @param tolerance  Absolute tolerance for floating-point comparisons.
    @param verbose    If true, progress messages are printed to stdout.
    */
    explicit PropertyHandle_Test( double tolerance = 5.0e-11,
                                   bool   verbose   = false );

    /**
    Destructor. Deletes the dynamically allocated model.
    All PropertyHandle members that created variables delete them
    via their own destructors before the model is deleted.
    */
    ~PropertyHandle_Test();

    /**
    Runs all sub-tests in sequence.
    */
    void run();

private:

    // --- sub-test methods ---

    /**
    Resets all variables to their canonical initial values.
    Called at the start of each operator test block.
    */
    void ResetVariables();

    /**
    Tests all operators on Region placements:
    NODE, ELEMENT, ELEMENT_INTEGRATION_POINT,
    SECTOR_INTEGRATION_POINT, FACET_INTEGRATION_POINT.
    */
    void TestRegionPlacements();

    /**
    Tests all operators on Boundary placements:
    NODE, FACE, FACE_INTEGRATION_POINT,
    FACE_SECTOR_INTEGRATION_POINT, FACE_FACET_INTEGRATION_POINT.
    */
    void TestBoundaryPlacements();

    /**
    Tests all operators on SplitBoundary placements:
    NODE (deduplicated InterFace FE nodes),
    INTER_FACE, INTER_FACE_INTEGRATION_POINT,
    INTER_FACE_SECTOR_INTEGRATION_POINT,
    INTER_FACE_FACET_INTEGRATION_POINT.
    */
    void TestSplitBoundaryPlacements();

    /**
    Tests cross-placement interpolation for Region handles.
    Verifies NODE <-> ELEMENT <-> ELEMENT_INTEGRATION_POINT paths.
    */
    void TestRegionCrossPlacement();

    /**
    Tests cross-placement interpolation for Boundary handles.
    Verifies NODE <-> FACE <-> FACE_INTEGRATION_POINT paths.
    */
    void TestBoundaryCrossPlacement();

    /**
    Tests cross-placement interpolation for SplitBoundary handles.
    Verifies NODE <-> INTER_FACE <-> INTER_FACE_INTEGRATION_POINT paths.
    */
    void TestSplitBoundaryCrossPlacement();

    /**
    Tests that operations on a Boundary handle do not modify
    variables on the Region, and vice versa.
    */
    void TestSubdomainIsolation();

    /**
    Tests that assignment from a whole-Model handle into a
    sub-domain handle is permitted, and that the reverse is
    correctly rejected.
    */
    void TestModelToSubdomainAssignment();

    /**
    Tests the Apply() method with a custom lambda on each
    subdomain type and placement.
    */
    void TestApplyMethod();

    /**
    Tests all named math methods (Ln, Exp, Sqrt, Pow, ZapNAN,
    Abs, Clamp, Sin, Cos, Tan, Acos, Asin, Atan) on a Region
    ELEMENT handle.
    */
    void TestMathMethods();

    /**
    Tests full support for ArrayVariable and FlaggedArrayVariable in
    PropertyHandle3, including:
      - Construction with ARRAY and FLAGGEDARRAY types
      - operator=( double ), operator=( ArrayVariable ),
        operator=( FlaggedArrayVariable )
      - Compound scalar operators +=, -=, *=, /= with per-element
        flag checking for FlaggedArrayVariable
      - ApplyArray with custom lambda
      - ApplyFlaggedArray with custom lambda
      - ApplyToFlaggedElements helper
      - Cross-type assignment rejection (ARRAY <-> SCALAR)
      - BOUNDARY placement rejection

    All verification uses direct Read from the first mesh cell to avoid
    accumulation errors from iterating all mesh points.
    */
    void TestArraySupport();

    /**
    Tests FlaggedArrayVariable on FACET_INTEGRATION_POINT placement with
    selective per-element flag protection.

    Test sequence:
      1. Initialise 10-element FLAGGEDARRAY to 0.0 with ALL flags = ANY.
      2. Set elements at prime indices (1,3,5,7) to DIRICH on every facet.
      3. Overwrite with -1.0 using OutputCondition = ANY:
           - ANY elements (0,2,4,6,8,9) become -1.0
           - DIRICH elements (1,3,5,7) remain 0.0
      4. Assign FLAGGEDARRAY handle into ARRAY handle via CopyValuesOnly:
           - All values transferred regardless of flags
           - ARRAY receives -1.0 at non-prime positions, 0.0 at prime positions
      5. Verify ARRAY has a single ANY flag (no per-element flags).

    Verified on every facet of the first mesh element.
    */
    void TestFlaggedArrayFacetIP();

    /**
    Tests all four compound operators (+=, -=, *=, /=) with cross-placement
    interpolation for SCALAR, VECTOR, and TENSOR variables on Region,
    mirroring the full coverage of PropertyHandle2_Test.
    */
    void TestRegionCrossPlacementFull();


  private:

    // -----------------------------------------------------------------------
    //  Model — declared FIRST so it is destroyed LAST
    //  PropertyHandle members are destroyed before model_ because
    //  members are destroyed in reverse declaration order.
    // -----------------------------------------------------------------------
    std::unique_ptr<csmp::Model<3>>  model_;  ///< owns the model

    const double TOLERANCE;
    const bool   verbose_;

    // -----------------------------------------------------------------------
    //  Region handles — SCALAR
    // -----------------------------------------------------------------------
    PropertyHandle<3> reg_sc_elem1;   ///< SCALAR ELEMENT,                   value 1.5
    PropertyHandle<3> reg_sc_elem2;   ///< SCALAR ELEMENT,                   value 2.0
    PropertyHandle<3> reg_sc_node1;   ///< SCALAR NODE,                      value 3.0
    PropertyHandle<3> reg_sc_node2;   ///< SCALAR NODE,                      value 4.0
    PropertyHandle<3> reg_sc_ip1;     ///< SCALAR ELEMENT_INTEGRATION_POINT, value 5.0
    PropertyHandle<3> reg_sc_ip2;     ///< SCALAR ELEMENT_INTEGRATION_POINT, value 6.0
    PropertyHandle<3> reg_sc_sec1;    ///< SCALAR SECTOR_INTEGRATION_POINT,  value 7.0
    PropertyHandle<3> reg_sc_sec2;    ///< SCALAR SECTOR_INTEGRATION_POINT,  value 8.0
    PropertyHandle<3> reg_sc_fac1;    ///< SCALAR FACET_INTEGRATION_POINT,   value 9.0
    PropertyHandle<3> reg_sc_fac2;    ///< SCALAR FACET_INTEGRATION_POINT,   value 10.0

    // -----------------------------------------------------------------------
    //  Region handles — VECTOR
    // -----------------------------------------------------------------------
    PropertyHandle<3> reg_vc_elem1;   ///< VECTOR ELEMENT,                   value 1.5
    PropertyHandle<3> reg_vc_elem2;   ///< VECTOR ELEMENT,                   value 2.0
    PropertyHandle<3> reg_vc_node1;   ///< VECTOR NODE,                      value 3.0
    PropertyHandle<3> reg_vc_node2;   ///< VECTOR NODE,                      value 4.0
    PropertyHandle<3> reg_vc_ip1;     ///< VECTOR ELEMENT_INTEGRATION_POINT, value 5.0
    PropertyHandle<3> reg_vc_ip2;     ///< VECTOR ELEMENT_INTEGRATION_POINT, value 6.0

    // -----------------------------------------------------------------------
    //  Region handles — TENSOR
    // -----------------------------------------------------------------------
    PropertyHandle<3> reg_ts_elem1;   ///< TENSOR ELEMENT,                   value 1.5
    PropertyHandle<3> reg_ts_elem2;   ///< TENSOR ELEMENT,                   value 2.0
    PropertyHandle<3> reg_ts_node1;   ///< TENSOR NODE,                      value 3.0
    PropertyHandle<3> reg_ts_node2;   ///< TENSOR NODE,                      value 4.0
    PropertyHandle<3> reg_ts_ip1;     ///< TENSOR ELEMENT_INTEGRATION_POINT, value 5.0
    PropertyHandle<3> reg_ts_ip2;     ///< TENSOR ELEMENT_INTEGRATION_POINT, value 6.0

    // -----------------------------------------------------------------------
    //  Boundary handles — SCALAR
    // -----------------------------------------------------------------------
    PropertyHandle<3> bnd_sc_node1;   ///< SCALAR NODE,                      value 1.0
    PropertyHandle<3> bnd_sc_node2;   ///< SCALAR NODE,                      value 2.0
    PropertyHandle<3> bnd_sc_face1;   ///< SCALAR FACE,                      value 3.0
    PropertyHandle<3> bnd_sc_face2;   ///< SCALAR FACE,                      value 4.0
    PropertyHandle<3> bnd_sc_fip1;    ///< SCALAR FACE_INTEGRATION_POINT,    value 5.0
    PropertyHandle<3> bnd_sc_fip2;    ///< SCALAR FACE_INTEGRATION_POINT,    value 6.0
    PropertyHandle<3> bnd_sc_fsec1;   ///< SCALAR FACE_SECTOR_IP,            value 7.0
    PropertyHandle<3> bnd_sc_fsec2;   ///< SCALAR FACE_SECTOR_IP,            value 8.0
    PropertyHandle<3> bnd_sc_ffac1;   ///< SCALAR FACE_FACET_IP,             value 9.0
    PropertyHandle<3> bnd_sc_ffac2;   ///< SCALAR FACE_FACET_IP,             value 10.0

    // -----------------------------------------------------------------------
    //  Boundary handles — VECTOR
    // -----------------------------------------------------------------------
    PropertyHandle<3> bnd_vc_face1;   ///< VECTOR FACE,                      value 3.0
    PropertyHandle<3> bnd_vc_face2;   ///< VECTOR FACE,                      value 4.0

    // -----------------------------------------------------------------------
    //  Boundary handles — TENSOR
    // -----------------------------------------------------------------------
    PropertyHandle<3> bnd_ts_face1;   ///< TENSOR FACE,                      value 3.0
    PropertyHandle<3> bnd_ts_face2;   ///< TENSOR FACE,                      value 4.0

    // -----------------------------------------------------------------------
    //  SplitBoundary handles — SCALAR
    // -----------------------------------------------------------------------
    PropertyHandle<3> spb_sc_node1;   ///< SCALAR NODE (InterFace FE),       value 1.0
    PropertyHandle<3> spb_sc_node2;   ///< SCALAR NODE (InterFace FE),       value 2.0
    PropertyHandle<3> spb_sc_iface1;  ///< SCALAR INTER_FACE,                value 3.0
    PropertyHandle<3> spb_sc_iface2;  ///< SCALAR INTER_FACE,                value 4.0
    PropertyHandle<3> spb_sc_ifip1;   ///< SCALAR INTER_FACE_IP,             value 5.0
    PropertyHandle<3> spb_sc_ifip2;   ///< SCALAR INTER_FACE_IP,             value 6.0
    PropertyHandle<3> spb_sc_ifsec1;  ///< SCALAR INTER_FACE_SECTOR_IP,      value 7.0
    PropertyHandle<3> spb_sc_ifsec2;  ///< SCALAR INTER_FACE_SECTOR_IP,      value 8.0
    PropertyHandle<3> spb_sc_iffac1;  ///< SCALAR INTER_FACE_FACET_IP,       value 9.0
    PropertyHandle<3> spb_sc_iffac2;  ///< SCALAR INTER_FACE_FACET_IP,       value 10.0

    // -----------------------------------------------------------------------
    //  SplitBoundary handles — VECTOR
    // -----------------------------------------------------------------------
    PropertyHandle<3> spb_vc_iface1;  ///< VECTOR INTER_FACE,                value 3.0
    PropertyHandle<3> spb_vc_iface2;  ///< VECTOR INTER_FACE,                value 4.0

    // -----------------------------------------------------------------------
    //  SplitBoundary handles — TENSOR
    // -----------------------------------------------------------------------
    PropertyHandle<3> spb_ts_iface1;  ///< TENSOR INTER_FACE,                value 3.0
    PropertyHandle<3> spb_ts_iface2;  ///< TENSOR INTER_FACE,                value 4.0

    // -----------------------------------------------------------------------
    //  Isolation sentinels — whole-model handles used to set known
    //  values across the entire mesh before sub-domain tests.
    // -----------------------------------------------------------------------
    PropertyHandle<3> whole_sc_elem;  ///< SCALAR ELEMENT whole-model sentinel
    PropertyHandle<3> whole_sc_node;  ///< SCALAR NODE    whole-model sentinel
};

} // namespace csmp

#endif // CSMP_PROPERTY_HANDLE3_TEST_H
