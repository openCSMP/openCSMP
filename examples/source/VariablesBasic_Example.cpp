#include "VariablesBasic_Example.h"
#include "CSMP_definitions.h"
#include "ScalarVariable.h"
#include "TensorVariable.h"
#include "VectorVariable.h"

#define DIM 3U

using namespace std;

namespace csmp{

void VariablesBasic_Example::Specifications()
{
  SetTitle( "Basic Variable Operations" );
  SetDifficulty( 1 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "P. Lang" );
  SetCategory( "Software Functionality" );
  AddDescription( "source in: VariablesBasic_Example.cpp" );
  AddDescription( "basic operations with csmp variable types Scalar, Vector, Tensor" );
}


void VariablesBasic_Example::Run()
{
  // a scalar variable, flagged PLAIN, instantiated with value of 5.
  ScalarVariable       sc( PLAIN, 5. );
  //                       ^^^^^ enumeration for later use in computations
  // may be used like built-in double
  sc *= 5.;
  // accessing the flag
  sc.Flag() = ANY;
  // offers Out() functionality
  sc.Out();

  // creating filled vector variables
  VectorVariable<DIM>  vc( ANY, ANY, ANY, 1, 2, 3 );
  //                       ^^^  ^^^  ^^^ component flags
  VectorVariable<DIM>  vc1( ANY, ANY, ANY, 1, 0, 0 );
  VectorVariable<DIM>  vc2( ANY, ANY, ANY, 0, 1, 0 );
  TensorVariable<DIM>  ts( ANY, ANY, ANY,
                           1, 2, 3,
                           4, 5, 6,
                           7, 8, 9 );
  // copying a variable
  TensorVariable<DIM>  ts3( ts );

  TensorVariable<DIM>  ts4 = ts * ts3;

  VectorVariable<DIM> result;

  result = vc1 % vc2 ;
  result.Out();

  ts4 = ts - ts3;
  ts4.Out();

  VectorVariable<DIM> res = vc * ts;
  res.Out();

} // Run()

} // csmp
