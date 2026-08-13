/*
 *  PropertyConstraints_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *  Revised and extended for full interface coverage.
 */

#ifndef CSMP_PROPERTY_CONSTRAINTS_TEST_H
#define CSMP_PROPERTY_CONSTRAINTS_TEST_H

#include "Test.h"
#include "Model.h"

namespace csmp {

/**
    Unit test for PropertyConstraints object.

    Covers:
    - Default and parameterised construction
    - Copy construction and copy assignment
    - Move construction and move assignment
    - AddConstraint / ChangeConstraint / DeleteConstraint
    - InitializePropertyIndices
    - WithIndexes / Constraints accessors
    - CheckConstraints (all three node-matching modes)
    - CheckConstraints with Index& (failed_upon) overload
    - CheckLengthOfVectorVariables
    - Erase
    - Building regions via Model::FormRegionFrom
    - PointInVolumeElement utility

    @note Requires a FracBox model and associated variable file.
    @note Requires a prism_test model for PointInVolumeElementTest.
*/
class PropertyConstraints_Test : public Test {
  public:
    explicit PropertyConstraints_Test( bool verbose = false );
    ~PropertyConstraints_Test();

    virtual void run();

  private:
    /// Default construction produces empty, usable object
    bool TestDefaultConstruction();

    /// Parameterised construction with a database initialises indices immediately
    bool TestParameterisedConstruction();

    /// Copy construction produces independent object with identical state
    bool TestCopyConstruction();

    /// Copy assignment produces independent object with identical state
    bool TestCopyAssignment();

    /// Move construction transfers state; source is left empty
    bool TestMoveConstruction();

    /// Move assignment transfers state; source is left empty
    bool TestMoveAssignment();

    /// AddConstraint returns true on first add, false on duplicate
    bool TestAddConstraint();

    /// ChangeConstraint updates both maps consistently
    bool TestChangeConstraint();

    /// DeleteConstraint removes from both maps
    bool TestDeleteConstraint();

    /// InitializePropertyIndices populates check_list_ from criteria_
    bool TestInitializePropertyIndices();

    /// WithIndexes and Constraints reflect current state correctly
    bool TestAccessors();

    /// Erase clears all state
    bool TestErase();

    /// CheckConstraints default mode: all nodes must satisfy constraints
    bool TestCheckConstraintsAllNodes();

    /// CheckConstraints with SatisfyConstraintsForAtLeastOneNode
    bool TestCheckConstraintsSingleNode();

    /// CheckConstraints with SatisfyConstraintsForNodalAverage
    bool TestCheckConstraintsNodalAverage();

    /// CheckConstraints with Index& overload reports failing index
    bool TestCheckConstraintsFailedUpon();

    /// CheckLengthOfVectorVariables switches vector checking mode
    bool TestCheckLengthOfVectorVariables();

    /// Mutual exclusivity of one_node_only_ and nodal_average_
    bool TestMutualExclusivity();

    /// Builds a region from property constraints via Model::FormRegionFrom
    bool TestBuildRegionsFromPropertyConstraints();

    /// Tests that pointInVolumeElement works using prism_test model
    bool PointInVolumeElementTest();

  private:
  
    const bool verbose_;
    Model<3U>* model_ = nullptr;
};

} // end csmp

#endif /* CSMP_PROPERTY_CONSTRAINTS_TEST_H */
