// PropertyHandle_MathTest.h

#ifndef CSMP_PROPERTY_HANDLE4_MATH_TEST_H
#define CSMP_PROPERTY_HANDLE4_MATH_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "PropertyHandle.h"

namespace csmp {

template<uint32_t> class Model;

/**
@brief Tests all named math methods of PropertyHandle for all five
       variable types (SCALAR, VECTOR, TENSOR, ARRAY, FLAGGEDARRAY)
       and representative placements.

@author 2024

@section overview Overview

PropertyHandle exposes the following named math methods:

  Squared, Sqrt, Ln, Log10, Exp, Pow, ZapNAN, Abs, Clamp,
  Sin, Cos, Tan, Acos, Asin, Atan

Each method is tested for:

  SCALAR  — ELEMENT placement, verified via direct Read from first cell.
  VECTOR  — ELEMENT placement, verified component-wise.
  TENSOR  — ELEMENT placement, verified component-wise.
  ARRAY   — ELEMENT placement, verified element-wise.
  FLAGGEDARRAY — ELEMENT placement, verified with DIRICH-protected
                 elements to confirm per-element flag checking.

@section angles Angles

All trigonometric methods expect and produce values in degrees.

@section flags FlaggedArrayVariable flag checking

For FLAGGEDARRAY handles, each math method must only modify elements
whose flag matches OutputCondition() (default ANY). Elements flagged
DIRICH must remain unchanged. This is verified explicitly for a
representative subset of methods.

@section mesh Test mesh

Uses "BoxHalfs3D" with variable file "PropertyHandle_MathTest-variables.txt".
*/
class PropertyHandle_MathTest : public Test
{
public:

    /**
    Constructs the test and initialises all PropertyHandle objects.

    @param tolerance  Absolute tolerance for floating-point comparisons.
    @param verbose    If true, progress messages are printed to stdout.
    */
    explicit PropertyHandle_MathTest( double tolerance = 1.0e-10,
                                       bool   verbose   = false );

    /** Destructor. Deletes the dynamically allocated model. */
    ~PropertyHandle_MathTest();

    /** Runs all math method sub-tests. */
    void run();

private:

    /** Tests all math methods for SCALAR ELEMENT handles. */
    void TestScalarMath();

    /** Tests all math methods for VECTOR ELEMENT handles. */
    void TestVectorMath();

    /** Tests all math methods for TENSOR ELEMENT handles. */
    void TestTensorMath();

    /** Tests all math methods for ARRAY ELEMENT handles. */
    void TestArrayMath();

    /**
    Tests all math methods for FLAGGEDARRAY ELEMENT handles.
    Verifies that DIRICH-flagged elements are not modified.
    */
    void TestFlaggedArrayMath();

    // -----------------------------------------------------------------------
    //  Model — declared FIRST so it is destroyed LAST
    //  PropertyHandle members are destroyed before model_ because
    //  members are destroyed in reverse declaration order.
    // -----------------------------------------------------------------------
    std::unique_ptr<csmp::Model<3>>  model_;  ///< owns the model

    const double    TOLERANCE;
    const bool      verbose_;

    // --- SCALAR ---
    PropertyHandle<3> sc_elem;   ///< SCALAR ELEMENT

    // --- VECTOR ---
    PropertyHandle<3> vc_elem;   ///< VECTOR ELEMENT

    // --- TENSOR ---
    PropertyHandle<3> ts_elem;   ///< TENSOR ELEMENT

    // --- ARRAY ---
    PropertyHandle<3> ar_elem;   ///< ARRAY ELEMENT, size 4

    // --- FLAGGEDARRAY ---
    PropertyHandle<3> fa_elem;   ///< FLAGGEDARRAY ELEMENT, size 4
};

} // namespace csmp

#endif // CSMP_PROPERTY_HANDLE3_MATH_TEST_H

