#ifndef FEM_DATA_TEST_H
#define FEM_DATA_TEST_H

#include "Test.h"
#include "FEM_Data.h"
#include "Index.h"

namespace csmp
{

/// PL Nov 2010
class FEM_Data_Test : public Test
{
public:
  virtual void run();
};

void FEM_Data_Test::run()
{
  const size_t simplexCount( 100 );

  // .) CONSTRUCTORS
  // . .) with placement and count
  FEM_Data<ScalarVariable> femDataScalarElement( ELEMENT, simplexCount );
  FEM_Data<VectorVariable<1U> > femDataVectorElement1( ELEMENT, simplexCount );
  FEM_Data<VectorVariable<2U> > femDataVectorElement2( ELEMENT, simplexCount );
  FEM_Data<VectorVariable<3U> > femDataVectorElement3( ELEMENT, simplexCount );
  FEM_Data<TensorVariable<2U> > femDataTensorElement2( ELEMENT, simplexCount );
  FEM_Data<TensorVariable<3U> > femDataTensorElement3( ELEMENT, simplexCount );

  FEM_Data<ScalarVariable> femDataScalarNode( NODE, simplexCount );
  FEM_Data<VectorVariable<1U> > femDataVectorNode1( NODE, simplexCount );
  FEM_Data<VectorVariable<2U> > femDataVectorNode2( NODE, simplexCount );
  FEM_Data<VectorVariable<3U> > femDataVectorNode3( NODE, simplexCount );
  FEM_Data<VectorVariable<3U> > femDataVectorNode3copy1( NODE, simplexCount );
  FEM_Data<TensorVariable<2U> > femDataTensorNode2( NODE, simplexCount );
  FEM_Data<TensorVariable<3U> > femDataTensorNode3( NODE, simplexCount );

  // . .) copy constructor
  FEM_Data<TensorVariable<3U> > femDataTensorNode3copy1( femDataTensorNode3 );
  FEM_Data<VectorVariable<3U> > femDataVectorElement3copy1( femDataVectorElement3 );

  // . .) with input vector
  std::vector<VectorVariable<3U> > inputVector;
  std::vector<TensorVariable<3U> > inputVectorTensor;
  VectorVariable<3U> vectorVariable( PLAIN, 1. );
  TensorVariable<3U> tensorVariable( PLAIN, 1. );
  for( size_t i = 0; i < simplexCount; ++i )
  {
    inputVector.push_back( vectorVariable );
    inputVectorTensor.push_back( tensorVariable );
  }
  FEM_Data<VectorVariable<3U> > femDataVectorNode3inputVector( NODE, inputVector );
  FEM_Data<VectorVariable<3U> > femDataVectorNode3inputVectorCopy( femDataVectorNode3inputVector );
  FEM_Data<TensorVariable<3U> > femDataVectorNode3inputTensor( NODE, inputVectorTensor );
  FEM_Data<TensorVariable<3U> > femDataVectorNode3inputTensorCopy( femDataVectorNode3inputTensor );


  // . .) default constructors
  FEM_Data<ScalarVariable> femDataScalarDefault;
  FEM_Data<VectorVariable<3U> > femDataVectorDefault;
  FEM_Data<TensorVariable<2U> > femDataTensorDefault;


  // .) PLACEMENT
  _test( femDataVectorNode2.Placement() == NODE );
  _test( femDataVectorElement2.Placement() == ELEMENT );
  _test( femDataScalarElement.Placement() == ELEMENT );
  _test( femDataTensorNode3.Placement() == NODE );

  // .) SIZE
  _test( femDataTensorNode2.Size() == simplexCount );
  _test( femDataVectorElement2.Size() == simplexCount );
  _test( femDataScalarElement.Size() == simplexCount );
  _test( femDataTensorNode3.Size() == simplexCount );

  // .) COMPARISON OP
  _test( femDataVectorNode3inputVector == femDataVectorNode3inputVectorCopy );
  _test( femDataVectorNode3inputTensorCopy == femDataVectorNode3inputTensor );

  // .) ACCESS OPS
  bool passed( true );
  for( size_t i = 0; i < femDataVectorNode3inputVector.Size(); ++i )
  {
    if( femDataVectorNode3inputVector( i ) != vectorVariable )
      passed = false;
  }
  _test( passed );
  passed = true;
  vectorVariable( 1 ) = 5.;
  femDataVectorNode3inputVector[ 25 ] = vectorVariable;
  _test( femDataVectorNode3inputVector( 25 ) == vectorVariable );
  _test( femDataVectorNode3inputVector[ 25 ]( 1 ) == 5. );

  // .) RESET
  const size_t simplexCount2( 50 );
  Index key( TENSOR, NODE, 99999 );
  femDataTensorElement3.Reset( key, simplexCount2, tensorVariable );
  _test( femDataTensorElement3.Placement() == NODE );
  _test( femDataTensorElement3.Size() == simplexCount2 );

  // .) EXTEND
  const size_t extendeCount( 200 );
  femDataScalarNode.Extend( extendeCount );
  _test( femDataScalarNode.Placement() == NODE );
  _test( femDataScalarNode.Size() == extendeCount );

  // .) REDUCETO
  const size_t simplexCount3( 10 );
  ScalarVariable scalarVariable( PLAIN, 1. );
  std::map<size_t,size_t> newAndOld;
  std::vector<ScalarVariable> scalarInputVector;
  for( size_t i = 0; i < simplexCount3; ++i )
  {
    scalarVariable() = (double)(i+1);
    scalarInputVector.push_back( scalarVariable );
    newAndOld.insert( std::make_pair( i, (simplexCount3-1)-i ) );
  }
  FEM_Data<ScalarVariable> reduceToTestFEM_Data( NODE, scalarInputVector );
  reduceToTestFEM_Data.ReduceTo( newAndOld );
  for( size_t i = 0; i < simplexCount3; ++i )
    if( reduceToTestFEM_Data[ i ]() != (double)(simplexCount3-i) )
      passed = false;
  _test( passed );
  passed = true;

  // .) LOG
/*
  _test( femDataTensorElement2.Logarithmitized() == false );
  ScalarVariable scalarVariableLog( PLAIN, 100. );
  std::vector<ScalarVariable> inputVectorPreLog( 3, scalarVariableLog );
  FEM_Data<ScalarVariable> logTest( NODE, inputVectorPreLog );
  logTest.LogarithmOfValues();
  for( size_t i = 0; i < logTest.Size(); ++i )
    _equal( logTest[i](), 4.60517, 0.01 );
  FEM_Data<ScalarVariable> logTest2( NODE, inputVectorPreLog );
  logTest2.DecadicLogarithmOfValues();
  for( size_t i = 0; i < logTest2.Size(); ++i )
    _equal( logTest2[i](), 2., 0.01 );

  // .) SQUARE ROOT
  FEM_Data<ScalarVariable> sqrtTest( NODE, inputVectorPreLog );
  sqrtTest.SquareRootOfValues();
  for( size_t i = 0; i < logTest.Size(); ++i )
    _equal( sqrtTest[i](), 10., 0.01 );
*/
  // .) RANGE OPERATIONS, GIVEMINMAX
  VectorVariable<3U> vectorVariableRange1( DIRICH, 5. );
  VectorVariable<3U> vectorVariableRange2( DIRICH, 15. );
  VectorVariable<3U> vectorVariableRange3( DIRICH, 25. );
  VectorVariable<3U> vectorVariableRange4( DIRICH, 35. );
  std::vector<VectorVariable<3U> > vectorVariableRange;
  vectorVariableRange.push_back( vectorVariableRange1 );
  vectorVariableRange.push_back( vectorVariableRange2 );
  vectorVariableRange.push_back( vectorVariableRange3 );
  vectorVariableRange.push_back( vectorVariableRange4 );
  FEM_Data<VectorVariable<3U> > femDataRange( ELEMENT, vectorVariableRange );
  VectorVariable<3U> vectorVariableRangeTest1, vectorVariableRangeTest2;
  VectorVariable<3U> vectorVariableRangeTest3, vectorVariableRangeTest4;
  femDataRange.MinMaxOf( vectorVariableRangeTest1, vectorVariableRangeTest2 );
  femDataRange.MinMaxOf( vectorVariableRangeTest3, vectorVariableRangeTest4 );
  _test( vectorVariableRangeTest1 == vectorVariableRangeTest3 );
  _test( vectorVariableRangeTest2 == vectorVariableRangeTest4 );
  _test( vectorVariableRangeTest1 == vectorVariableRange1 );
  _test( vectorVariableRangeTest2 == vectorVariableRange4 );
  femDataRange.ScaleRangeTo( vectorVariableRange2, vectorVariableRange3 );
  for( size_t i = 0; i < femDataRange.Size(); ++i )
  {
    _test( femDataRange[ i ] !=  vectorVariableRange1 );
    _test( femDataRange[ i ] !=  vectorVariableRange4 );
  }
  femDataRange.OffsetRangeBy( vectorVariable );
  _test( femDataRange[ 0 ] != vectorVariableRange2 );

  // .) BINARY OPS
  std::FILE * femDataBinary;
  femDataBinary = std::fopen( "femDataBinary", "w" );
  femDataVectorNode3inputVector.OutBinary( femDataBinary );
  std::fclose( femDataBinary );
  femDataBinary = std::fopen( "femDataBinary", "r" );
  femDataVectorNode3inputVectorCopy.InBinary( femDataBinary );
  std::fclose( femDataBinary );
  _test( femDataVectorNode3inputVector == femDataVectorNode3inputVectorCopy );


} // run()


} // csmp

#endif // FEM_DATA_TEST_H
