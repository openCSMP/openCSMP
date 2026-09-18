#ifndef FINITE_ELEMENT_TEST_H
#define FINITE_ELEMENT_TEST_H

#include "Test.h"
#include "FiniteElement.h"
#include "FiniteElement_TestData.h"

namespace csmp{

/// Generic Test for any FEM subclass
class FiniteElement_Test : public Test
{
  public:
    FiniteElement_Test( FiniteElement* testee, const char* results_file, bool verbose );
    
    virtual void run();
    
    ~FiniteElement_Test();
    
    /// testing declerative programming extensions
    void TestLocalCoordinateFunctionsReturnedAs_vector( int interpolation_order );
    
  private:
    std::string fileName_;
    FiniteElement* femPtr_ = nullptr;
    FiniteElement_TestData femData_;
    const bool verbose_;
};

/**

@class FiniteElement_Test
@author Philipp Lang

Uses the FiniteElement_TestData class to test FEMs. Please see there for instructions
on the text file that is used for testing.

@todo (3) ConsecutiveNodesAtBoundary Test(not implemented for most FEMs)
@todo (3) Think about _equal vs _test in some cases
@todo (3) Add more feedback when ordering matters

@section prerequisites FEM Prerequisites

Coordinate matrix has to be at least of size 2x2, meaning no less
than a 2D line element may be tested.

*/

}

#endif // FINITE_ELEMENT_TEST_H
