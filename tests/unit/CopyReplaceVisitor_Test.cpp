#include "CopyReplaceVisitor_Test.h"
#include "PropertyHandle.h"
#include "CSMP_definitions.h"
#include "CopyReplaceVisitor.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"

namespace csmp {

  
void CopyReplaceVisitor_Test::run()
{
    // Test for 3D only, creating auxilliaries
    enum{DIM=3U};
    ScalarVariable scalarValue( ROBIN, 2. );
    VectorVariable<DIM> vectorValue( ROBIN, 2. );
    TensorVariable<DIM> tensorValue( ROBIN, 2. );

    // nodal scalar
    PropertyHandle<DIM> nodeScalar1( *model_, "node scalar 1", SCALAR, NODE );
    PropertyHandle<DIM> nodeScalar2( *model_, "node scalar 2", SCALAR, NODE );
    model_->InputPropertyValue( "node scalar 1", scalarValue );
    CopyReplaceVisitor<ScalarVariable,DIM> cpyNodeScalar1to2( model_->Database(), "node scalar 1", "node scalar 2");
    model_->Accept(cpyNodeScalar1to2);
    testNodes( model_, scalarValue, "node scalar 2" );  

    // nodal vector
    PropertyHandle<DIM> nodeVector1( *model_, "node vector 1", VECTOR, NODE );
    PropertyHandle<DIM> nodeVector2( *model_, "node vector 2", VECTOR, NODE );
    model_->InputPropertyValue( "node vector 1", vectorValue );
    CopyReplaceVisitor<VectorVariable<DIM>,DIM> cpyNodeVector1to2( model_->Database(), "node vector 1", "node vector 2" );
    model_->Accept(cpyNodeVector1to2);
    testNodes( model_, vectorValue, "node vector 2" );  

    // nodal tensor
    PropertyHandle<DIM> nodeTensor1( *model_, "node tensor 1", TENSOR, NODE );
    PropertyHandle<DIM> nodeTensor2( *model_, "node tensor 2", TENSOR, NODE );
    model_->InputPropertyValue( "node tensor 1", tensorValue );
    CopyReplaceVisitor<TensorVariable<DIM>,DIM> cpyNodeTensor1to2( model_->Database(), "node tensor 1", "node tensor 2" );
    model_->Accept(cpyNodeTensor1to2);
    testNodes( model_, tensorValue, "node tensor 2" );

    // element scalar
    PropertyHandle<DIM> elementScalar1( *model_, "element scalar 1", SCALAR, ELEMENT );
    PropertyHandle<DIM> elementScalar2( *model_, "element scalar 2", SCALAR, ELEMENT );
    model_->InputPropertyValue( "element scalar 1", scalarValue );
    CopyReplaceVisitor<ScalarVariable,DIM> cpyElementScalar1to2( model_->Database(), "element scalar 1", "element scalar 2" );
    model_->Accept(cpyElementScalar1to2);
    testElements( model_, scalarValue, "element scalar 2" );  

    // element vector
    PropertyHandle<DIM> elementVector1( *model_, "element vector 1", VECTOR, ELEMENT );
    PropertyHandle<DIM> elementVector2( *model_, "element vector 2", VECTOR, ELEMENT );
    model_->InputPropertyValue( "element vector 1", vectorValue );
    CopyReplaceVisitor<VectorVariable<DIM>,DIM> cpyElementVector1to2( model_->Database(), "element vector 1", "element vector 2" );
    model_->Accept(cpyElementVector1to2);
    testElements( model_, vectorValue, "element vector 2" );  

    // element tensor
    PropertyHandle<DIM> elementTensor1( *model_, "element tensor 1", TENSOR, ELEMENT );
    PropertyHandle<DIM> elementTensor2( *model_, "element tensor 2", TENSOR, ELEMENT );
    model_->InputPropertyValue( "element tensor 1", tensorValue );
    CopyReplaceVisitor<TensorVariable<DIM>,DIM> cpyElementTensor1to2( model_->Database(), "element tensor 1", "element tensor 2" );
    model_->Accept(cpyElementTensor1to2);
    testElements( model_, tensorValue, "element tensor 2" );
} // run

CopyReplaceVisitor_Test::~CopyReplaceVisitor_Test()
{
}

template<class V>
void CopyReplaceVisitor_Test::testNodes( Model<3>* model, const V& value, const char* propertyName )
{
  V cache;
  Index key( model->Database().StorageKey(propertyName) );
  Region<3U>& domain{ model->Region("Model") };
  
  for( auto it( domain.NodesBegin() ); it != domain.NodesEnd(); ++it )
  {
    (*it)->Read( key, cache );
    _test( cache == value );
  }  
}


template<class V>
void CopyReplaceVisitor_Test::testElements( Model<3>* model, const V& value, const char* propertyName )
{
  V cache;
  Index key( model->Database().StorageKey(propertyName) );
  Region<3U>& domain{ model->Region("Model") };

  for( auto it( domain.CellsBegin() ); it != domain.CellsEnd(); ++it )
  {
    (*it)->Read( key, cache );
    _test( cache == value );
  }  
}

template<class V>
void CopyReplaceVisitor_Test::testElementIntegrationPoints( Model<3>* model, const V& value, const char* propertyName )
{
  V cache;
  Index key( model->Database().StorageKey(propertyName) );
  Region<3U>& domain{ model->Region("Model") };

  for( auto it( domain.CellsBegin() ); it != domain.CellsEnd(); ++it )
    for( auto ip{0U}; ip < (*it)->IntegrationPoints(); ++ip )
      {
        (*it)->Read( key, cache );
        _test( cache == value );
      }  
}


} // csmp
