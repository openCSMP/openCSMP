#ifndef TDDUNITTEST_EXAMPLE_H
#define TDDUNITTEST_EXAMPLE_H

#include "Example.h"

namespace csmp {

/// Introduces to Test Driven Development in general and attempts to provide a CSMP++ specific approach.
class  TDDUnitTest_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
};

/**
@class  TDDUnitTest_Example TDDUnitTest_Example "example_suite/TDDUnitTest_Example.h"

UnitTests and Test Driven Development increase confidence in the code you write and
drastically reduce debugging needs. In an environment like CSMP that is not tested
bottom up, a different approach to testing new classes may be chosen.
Also for classes that act on a csmp::Model, the user might skip test
driven development/unit testing due to the high effort required. This document tries
to provide an approach to CSMP specific TDD.

@author P Lang
@date   September 2010
@note  WORK IN PROGRESS

@section outline Outline
- Unit Tests and Functional Tests
 - in general XP
 - interpreted for CSMP
- Introduction to TDD
- Implementation

@section tests Unit Tests and Functional Tests

@subsection general In general XP

A Unit Test is designed to test the interface of a single class through proper pass/fail
criteria. A Functional Test tests the result of class interactions within a bigger context.
Hence, pass/fail criteria are more complex. This also includes performance test and the
validation of applications in their whole.

@subsection csmp Interpreted for CSMP

To test classes that require a model to act on, multiple Unit Tests may be run within
one test function and work with the same model. Classes are tested indepently from another,
and the TestSuite will spot failures individually. To sum it up, we put multiple tests
inside a single test function that establishes a model for all of them. This is
done to reduce the effort that goes into creating a test.
Functional Tests, or Test Cases, measure performance and compare results of complex applications.

@section tdd Introduction to TDD

Write your Unit Test first, let it fail(because there is no class yet, afterall) then
implement your class piecewise to make the tests pass. Once the tests pass, improve your code.
Benefits:
- improves design
 - think it through
 - then implement
- reduces amount of code
 - by sticking to what you really planned to do
 - by having code that does the job quickly
  - and improve it when there's time
- used together with SVN, often replaces debugger
 - revert to revision that works, go form there
- YOU CAN TRUST IN YOUR CODE AN SUBSEQUENT USE OF IT
- IT IS TESTED(what's the alternative?)

How it's done:
-# Write your test
-# See it fail
-# Implement your code
-# Make the test pass

@section implementation Implementation

Basic Unit Test and TDD with CSMP's csmp::Test and csmp::TestSuite class.
Available Macros:
-_test( cond ) tests for logical statement
-_fail( stri ) fails
-_equal( arg1, arg2, tolerance ) tests if arguments equal within tolerance

@code
  // implementation follos
  class Unit
  {
  public:
    int foo1( int a ) { return a/a; }
    int foo2( int a ) { return a*a; }
  };

  // this comes first(failing)
  class  UnitTest : public Test
  {
  public:
    virtual void run() { _test( Unit().foo1( 5 ) == 1 );
                         _test( Unit().foo2( 5 ) > 5 ); }
  };

int main( int argc, char** argv)
{
  TestSuite s("TDD SUITE", &cout);

  s.addTest( new UnitTest() );

  s.run();
  int failures = s.report();

  return failures;
 }
@endcode

Multiple Unit Tests within one test function for a group of classes to be tested,
e.g. Interrelations.

@code

// global function to reload model from binary
void resetModel( Model<DIM>* model )
{
  if( model != 0 )
    delete model;
  model = new ANSYS_Model3D( "small_prism_bin",
                             "CSMP-2phase-variables.txt" );
}

class  MultipleUnitTest : public Test
{
public:
  virtual void run() {

    // establishing a test model
    Model<DIM>* model;
    resetModel( model );


    // getting indices to evaluate test criterion
    const Region<dim>& rref = model->Region( "Model" );
    Index porosityKey				= model->Database().StorageKey( "porosity" );
    Index frapKey					= model->Database().StorageKey( "fracture aperture" );

    // backing up property that will be manipulated
    PropertyHandle<DIM>* cachePorosity = new PropertyHandle<DIM>( (*model), "cache porosity", SCALAR, ELEMENT );
    Index cachePorosityKey = model->Database().StorageKey( "cache porosity" );
    (*model).CopyReplace( "porosity", "cache porosity" );

    // Applying an Interrelation to the model
    FracturePorosityFromAperture<DIM> fracturePorosityInterrelation( (*model).Database() );
    (*model).Apply( fracturePorosityInterrelation );

    // Running a trivial test criterion
    bool passed( true );
    const std::vector<Element<DIM>*>::const_iterator elementsEnd = rref.ElementsEnd();
    for( std::vector<Element<DIM>*>::const_iterator it = rref.ElementsBegin();
         it != elementsEnd; ++it )
      {
        if( (*it)->Read( porosityKey ) != (*it)->Read( cachePorosityKey ) *
                                          (*it)->Read( frapKey ) ) 
		  {
          passed = false;
		  break;
		  }

		// if we wish to test more than one interrelation within one loop, we do
		// _equal( (*it)->Read( cacheValue1Key ), (*it)->Read( value1Key ) ), tolerance_ );
		// _equal( (*it)->Read( cacheValue2Key ), (*it)->Read( value2Key ) ), tolerance_ );
		// ....
		// as long as the tests still qualify as independent
      }

      // evaluating test criterion
      _test( passed );

      // restoring model
      (*model).CopyReplace( "cache porosity", "porosity" );
      delete cachePorosity;

      //...tests

      // resetting model if required
      resetModel( model );

      //...tests


  } // run

};

int main( int argc, char** argv)
{
  TestSuite s("TDD SUITE", &cout);

  s.addTest( new MultipleUnitTest() );

  s.run();
  int failures = s.report();

  return failures;
 }
@endcode

Points to consider:
- load model from a CSMP binary
- use model.CopyReplace(...) if possible

For developers, IMHO there's no alternativ to writing a single
test per core class:
@code
int main( int argc, char** argv)
{
  TestSuite s("TDD SUITE", &cout);

  // establishing model
  enum{ DIM = 2 };
  Model<DIM>* model;
  resetModel( model )

  // the tests themselve have the responsability to restore the model
  // to its default
  s.addTest( new UnitTest1( model ) );
  s.addTest( new UnitTest2( model ) );
  s.addTest( new UnitTest3( model ) );
  s.addTest( new UnitTest4( model ) );

  s.run();
  int failures = s.report();

  return failures;
 }
@endcode

@note 10 MIN BUILD
@todo (3) Implement Suite features for functioal tests(csmp binary comparison, speed, vtk comparison(DONE)...)
*/

} // csmp

#endif // TDDUNITTEST_EXAMPLE_H
